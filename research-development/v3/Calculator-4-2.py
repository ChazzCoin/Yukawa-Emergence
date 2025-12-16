#!/usr/bin/env python3
import json
from numpy.linalg import eigh
from scipy.linalg import expm
import numpy as np
from dataclasses import dataclass, asdict
from typing import Optional, Dict, Any, List, Tuple
import csv

@dataclass
class ScanResult:
    ok_3d: bool
    blocks: Any
    g_diag: float
    g_off: float
    eps: float
    phase_pi: float  # phase expressed as multiple of pi (e.g. 0.075 means 0.075π)
    beta: float
    rel_cut: float
    tol_rel_blocks: float
    seesaw_M: float

    theta12: Optional[float] = None
    theta13: Optional[float] = None
    theta23: Optional[float] = None

    ue3: Optional[float] = None
    um3: Optional[float] = None
    ut3: Optional[float] = None

    m_light: Optional[List[float]] = None
    dm21: Optional[float] = None
    dm31: Optional[float] = None

    score: Optional[float] = None
    note: str = ""


def set_geometry_weights_pmns_two_phase(
    cfg,
    g_diag=0.894,
    g_off=1.00,
    eps=0.012,
    phase12=0.075*np.pi,   # gen1 ↔ plane
    phase23=0.010*np.pi    # weak plane-breaking phase
):
    base = {
        (0,0): 1.00,
        (0,1): 0.92 * np.exp(1j * phase12),
        (0,2): 0.85,
        (1,0): 0.95 * np.exp(-1j * phase12),
        (1,1): 1.10,
        (1,2): 0.88 * np.exp(1j * phase23),   # NEW
        (2,0): 0.80,
        (2,1): 0.90 * np.exp(-1j * phase23),  # NEW
        (2,2): 1.05,
    }

    diag = {(0,0), (1,1), (2,2)}
    gen1_mix = {(0,1), (1,0)}

    for g in base:
        if g in diag:
            base[g] *= g_diag
        else:
            base[g] *= g_off
        if g in gen1_mix:
            base[g] *= (1.0 + eps)

    cfg.geometry_weights = base

def angles_from_absU(absU):
    s13 = absU[0,2]
    c13 = np.sqrt(max(0.0, 1.0 - s13**2))

    s12 = absU[0,1] / c13 if c13 > 0 else np.nan
    s23 = absU[1,2] / c13 if c13 > 0 else np.nan

    th12 = np.degrees(np.arcsin(np.clip(s12, -1, 1)))
    th13 = np.degrees(np.arcsin(np.clip(s13, -1, 1)))
    th23 = np.degrees(np.arcsin(np.clip(s23, -1, 1)))

    return th12, th13, th23


# ============================================================
#  CONFIG (same as yours, unchanged)
# ============================================================

class AlignmentV33Config:
    group_elements = [(i, j) for i in range(3) for j in range(3)]
    subgroup_H = [(0, 0), (1, 1), (2, 2)]
    triad_shifts = [(0, 0), (1, 0), (0, 1)]

    kernel_characters = [
        (1, 0,  1.0),
        (0, 1,  0.6),
        (1, 1,  0.35),
    ]

    geometry_weights = {
        (0,0): 1.00, (0,1): 0.92, (0,2): 0.85,
        (1,0): 0.95, (1,1): 1.10, (1,2): 0.88,
        (2,0): 0.80, (2,1): 0.90, (2,2): 1.05,
    }

    damping_strength = 0.35

    compression_characters = [
        (0, 0),
        (1, 0),
        (1, 1),
    ]

    higgs_vev = 174.0


# ============================================================
#  PIPELINE
# ============================================================

class AlignmentPipeline:
    """
    AlignmentPipeline encapsulates the operator sequence:

      K --(misalignment flow)--> K_flow
        --(emergent C360 projector)--> K_proj
        --(compression)--> Y
        --(block diagnosis)--> blocks on Y
        --(conditional seesaw)--> neutrino-like light eigenvalues

    This makes quark/lepton emergence automatic:
      - 1D harmonic blocks => Dirac-like (no seesaw)
      - >=2D harmonic blocks => closure failure => seesaw extension
    """

    def __init__(
        self,
        cfg: AlignmentV33Config,
        beta: float = 1.5,
        rel_cut: float = 0.15,
        tol_rel_blocks: float = 0.03,
        seesaw_M: float = 1e6,
        compression_shear: float = 0.0
    ):
        self.cfg = cfg
        self.beta = float(beta)
        self.rel_cut = float(rel_cut)
        self.tol_rel_blocks = float(tol_rel_blocks)
        self.seesaw_M = float(seesaw_M)
        self.compression_shear = float(compression_shear)
        # caches (filled by run())
        self.triads = None
        self.K = None
        self.S = None
        self.K_flow = None
        self.P_C360 = None
        self.kept_indices = None
        self.K_proj = None
        self.Y = None

    # ---------- group + characters ----------

    @staticmethod
    def add_g(a, b):
        return ((a[0]+b[0]) % 3, (a[1]+b[1]) % 3)

    @staticmethod
    def sub_g(a, b):
        return ((a[0]-b[0]) % 3, (a[1]-b[1]) % 3)

    @staticmethod
    def chi(g, p, q):
        i, j = g
        return np.exp(2j*np.pi*(p*i+q*j)/3.0)

    def build_S(self, triads, shear=0.0):
        """
        Build compression S with an optional tiny phase shear applied only
        to the gen-1 row (i=0) on its non-anchor sites inside the triad.
        This gently breaks residual symmetry that pins |Ue3|.
        """
        cfg = self.cfg
        G = cfg.group_elements
        S = np.zeros((3, 9), dtype=complex)

        for i, triad in enumerate(triads):
            p, q = cfg.compression_characters[i]
            for k, idx in enumerate(triad):
                g = G[idx]
                ph = self.chi(g, p, q)

                # shear only on generation-1 row, excluding the first triad element
                if (i == 0) and (k > 0) and (shear != 0.0):
                    ph *= np.exp(1j * shear)

                S[i, idx] = ph / np.sqrt(3)

        return S

    # ---------- charged-lepton (Dirac) rotation ----------

    @staticmethod
    def _rotation_matrix_13(theta_rad):
        """
        Small unitary rotation in the 1–3 plane.
        """
        c = np.cos(theta_rad)
        s = np.sin(theta_rad)
        return np.array([
            [ c, 0.0,  s],
            [0.0, 1.0, 0.0],
            [-s, 0.0,  c],
        ], dtype=complex)

    def extract_charged_lepton_rotation(self, Y_e):
        """
        Extract a *small* charged-lepton left-handed rotation Ue in the 1–3 plane.

        We do NOT use full Ue from diagonalizing Y_e as the PMNS factor because that
        permutes/flips columns unpredictably across scan points.

        Instead we extract a robust small angle from the charged-lepton eigenvectors.
        """
        evals, U_full = eigh(Y_e)

        # Heaviest charged lepton eigenstate ~ tau is the largest |eigenvalue|
        # (works even if signs flip; we only need a stable direction)
        idx_tau = int(np.argmax(np.abs(np.real(evals))))
        v_tau = U_full[:, idx_tau]

        # Extract a small 1–3 mixing angle: electron component inside tau eigenvector
        # Use magnitude (stable); sign from real parts if possible.
        s = float(np.clip(np.abs(v_tau[0]), 0.0, 1.0))
        theta = float(np.arcsin(s))

        # Give theta a sign using the relative real parts (optional but helps continuity)
        sign = 1.0
        if np.real(v_tau[0]) * np.real(v_tau[2]) < 0:
            sign = -1.0
        theta *= sign

        # Clamp to "small charged-lepton correction" regime (safety)
        theta = float(np.clip(theta, -0.10, 0.10))  # ~ ±5.7°

        Ue = self._rotation_matrix_13(theta)
        return Ue, theta

    # ---------- build K (geometry-selected, non-convolution) ----------

    def build_kernel(self):
        cfg = self.cfg
        G = cfg.group_elements
        n = len(G)
        K = np.zeros((n, n), dtype=complex)

        for a, g in enumerate(G):
            for b, h in enumerate(G):
                # spectral part
                F = sum(
                    w * self.chi(self.sub_g(g, h), p, q)
                    for (p, q, w) in cfg.kernel_characters
                )

                # geometry weights (B)
                alpha_g = cfg.geometry_weights[g]
                alpha_h = cfg.geometry_weights[h]

                # misalignment damping W(g,h)
                dist = min(abs(g[0]-h[0]), 3-abs(g[0]-h[0])) + \
                       min(abs(g[1]-h[1]), 3-abs(g[1]-h[1]))
                W = np.exp(-cfg.damping_strength * dist)

                K[a, b] = alpha_g * F * np.conj(alpha_h) * W

        # enforce Hermitian
        return 0.5 * (K + K.conj().T)

    # ---------- triads + compression S ----------

    def build_triads(self):
        cfg = self.cfg
        index = {g: i for i, g in enumerate(cfg.group_elements)}
        triads = []
        for s in cfg.triad_shifts:
            triads.append([index[self.add_g(h, s)] for h in cfg.subgroup_H])
        return triads


    # ---------- operators: flow, C360 projector, compression ----------

    def misalignment_flow(self, K):
        return expm(-self.beta * K)

    def emergent_C360_projector(self, K):
        evals, evecs = eigh(K)
        flowed = np.exp(-self.beta * evals)

        max_val = flowed.max()
        keep = flowed >= self.rel_cut * max_val
        kept = np.where(keep)[0]

        P = np.zeros_like(K, dtype=complex)
        for i in kept:
            v = evecs[:, i:i+1]
            P += v @ v.conj().T

        return P, kept.tolist()

    @staticmethod
    def effective_yukawa(K_like, S):
        return S @ K_like @ S.conj().T

    # ---------- block diagnosis on Y ----------

    def harmonic_blocks_on_Y(self, Y):
        evals, _ = eigh(Y)
        evals = np.sort(np.real(evals))

        blocks = []
        block = [0]
        for i in range(1, len(evals)):
            scale = max(1.0, abs(evals[i-1]), abs(evals[i]))
            if abs(evals[i] - evals[i-1]) <= self.tol_rel_blocks * scale:
                block.append(i)
            else:
                blocks.append(block)
                block = [i]
        blocks.append(block)
        return blocks, evals

    # ---------- conditional seesaw on >=2D blocks ----------

    @staticmethod
    def _seesaw_light_eigs_for_block(Y, block, M):
        Yb = Y[np.ix_(block, block)]
        n = len(block)

        zero = np.zeros_like(Yb)
        MR = M * np.eye(n)

        big = np.block([
            [zero, Yb],
            [Yb.conj().T, MR]
        ])

        eigvals, _ = eigh(big)
        eigvals = np.sort(np.abs(eigvals))
        return eigvals[:n]

    def _build_Y_with_geometry(self, geometry_weights_override=None):
        """
        Build Y = S K_proj S† but optionally using a different geometry_weights
        (charged leptons vs neutrinos), while keeping kernel_characters/damping/triads/S consistent.
        """
        cfg = self.cfg

        old_geom = cfg.geometry_weights
        if geometry_weights_override is not None:
            cfg.geometry_weights = geometry_weights_override

        try:
            K = self.build_kernel()
            K_flow = self.misalignment_flow(K)
            P, _kept = self.emergent_C360_projector(K)
            K_proj = P @ K_flow @ P

            triads = self.build_triads()

            # IMPORTANT: keep shear consistent across sectors too
            S = self.build_S(triads, shear=self.compression_shear)

            Y = self.effective_yukawa(K_proj, S)
            return Y
        finally:
            cfg.geometry_weights = old_geom

    def apply_conditional_seesaw(self, Y, blocks):
        """
        Returns:
          - dirac_blocks: list of (block, eigenvalues) for 1D blocks
          - majorana_blocks: list of (block, light_eigs_after_seesaw) for >=2D blocks
        """
        dirac = []
        majorana = []

        # we use eigenvalues of the sub-block as "Dirac-like" for 1D
        # (you can later map these to charged lepton vs quark readouts)
        for block in blocks:
            if len(block) == 1:
                i = block[0]
                dirac.append((block, np.array([np.real(eigh(Y)[0][i])])))
            else:
                light = self._seesaw_light_eigs_for_block(Y, block, self.seesaw_M)
                majorana.append((block, light))

        return dirac, majorana

    def _canonicalize_pmns(self, U):
        """
        Make PMNS-like matrices stable across scans:
          - Put the '3' column as the one with smallest |U_ei|  (i=3 => smallest electron component)
          - Order remaining two by |U_ei| descending (so col0 has largest |U_e|)
          - Fix phases so U[0,k] is real and >= 0 for each column
        """
        U = np.array(U, dtype=complex, copy=True)

        ue = np.abs(U[0, :])

        # column with smallest electron component -> column 2 (mass state 3)
        k3 = int(np.argmin(ue))
        rest = [k for k in range(3) if k != k3]

        # among the remaining two, order by decreasing |Ue|
        rest = sorted(rest, key=lambda k: ue[k], reverse=True)

        perm = rest + [k3]
        U = U[:, perm]

        # phase-fix each column so U[0,k] is real positive
        for k in range(3):
            ph = np.angle(U[0, k])
            U[:, k] *= np.exp(-1j * ph)
            if np.real(U[0, k]) < 0:
                U[:, k] *= -1.0

        return U

    def apply_full_seesaw_3x3(self, Y, include_charged_lepton=True, geometry_weights_charged=None):
        """
        Full seesaw on 3x3 Y -> Unu.
        Optionally apply a *small* charged-lepton rotation Ue† (NOT full diagonalization).
        Returns:
            light_masses, U_pmns_or_Unu, diagnostics
        """
        assert Y.shape == (3, 3)

        # --- build 6x6 seesaw operator ---
        zero = np.zeros_like(Y)
        MR = self.seesaw_M * np.eye(3)

        big = np.block([
            [zero, Y],
            [Y.conj().T, MR]
        ])

        eigvals, eigvecs = eigh(big)
        eigvals = np.real(eigvals)

        # sort by absolute value => light first
        order = np.argsort(np.abs(eigvals))
        eigvals = eigvals[order]
        eigvecs = eigvecs[:, order]

        light_masses = np.abs(eigvals[:3])
        Unu = eigvecs[:3, :3].copy()

        # normalize columns
        for k in range(3):
            nrm = np.linalg.norm(Unu[:, k])
            if nrm > 0:
                Unu[:, k] /= nrm

        # --- canonicalize neutrino mixing (column order + phase fix) ---
        Unu = self._canonicalize_pmns(Unu)

        diagnostics = {}

        if not include_charged_lepton:
            return light_masses, Unu, diagnostics

        # --- build charged-lepton operator Y_e (default: strip phases only) ---
        if geometry_weights_charged is None:
            geom_e = {k: np.abs(v) for k, v in self.cfg.geometry_weights.items()}
            diagnostics["Y_e_geometry"] = "abs(phased_geometry)"
        else:
            geom_e = geometry_weights_charged
            diagnostics["Y_e_geometry"] = "override"

        Y_e = self._build_Y_with_geometry(geometry_weights_override=geom_e)

        # extract small Ue (1–3 plane) and apply
        Ue, theta_e = self.extract_charged_lepton_rotation(Y_e)
        diagnostics["theta_e_rad"] = float(theta_e)
        diagnostics["theta_e_deg"] = float(np.degrees(theta_e))
        diagnostics["Ue_from"] = "small_13_rotation_from_Ye"

        U_pmns = Ue.conj().T @ Unu

        # canonicalize again (small Ue can swap phases slightly)
        U_pmns = self._canonicalize_pmns(U_pmns)

        return light_masses, U_pmns, diagnostics

    # ---------- run the whole pipeline ----------

    def run(self):
        cfg = self.cfg

        self.triads = self.build_triads()
        self.K = self.build_kernel()

        # IMPORTANT: actually use compression_shear here
        self.S = self.build_S(self.triads, shear=self.compression_shear)

        self.K_flow = self.misalignment_flow(self.K)
        self.P_C360, self.kept_indices = self.emergent_C360_projector(self.K)
        self.K_proj = self.P_C360 @ self.K_flow @ self.P_C360

        self.Y = self.effective_yukawa(self.K_proj, self.S)

        # Y spectrum
        y_vals, U = eigh(self.Y)
        masses = np.abs(y_vals) * cfg.higgs_vev

        # block structure
        blocks, y_sorted = self.harmonic_blocks_on_Y(self.Y)

        # conditional seesaw where needed
        dirac_blocks, majorana_blocks = self.apply_conditional_seesaw(self.Y, blocks)

        return {
            "Y": self.Y,
            "Y_eigvals": y_vals,
            "masses_GeV": masses,
            "U": U,
            "absU": np.abs(U),
            "blocks": blocks,
            "Y_sorted_eigvals": y_sorted,
            "dirac_blocks": dirac_blocks,
            "majorana_blocks": majorana_blocks,
            "kept_indices_C360": self.kept_indices,
        }



import numpy as np
from numpy.linalg import eigh

# ============================================================
#  Physics-grade PMNS utilities (no pipeline changes)
# ============================================================

def _phase_fix_columns(U: np.ndarray) -> np.ndarray:
    """Rephase columns so U[0,k] is real and >= 0 (column-wise)."""
    U = np.array(U, dtype=complex, copy=True)
    for k in range(U.shape[1]):
        ph = np.angle(U[0, k])
        U[:, k] *= np.exp(-1j * ph)
        if np.real(U[0, k]) < 0:
            U[:, k] *= -1.0
    return U

def _jarlskog(U: np.ndarray) -> float:
    """J = Im(Ue1 Uμ2 Ue2* Uμ1*)."""
    return float(np.imag(U[0,0] * U[1,1] * np.conj(U[0,1]) * np.conj(U[1,0])))

def _delta_cp_from_U(U: np.ndarray) -> float:
    """
    Extract δ_CP (degrees) using:
      J = c12 c13^2 c23 s12 s13 s23 sin δ
    Returns principal value in [-90, 90] degrees (arcsin branch).
    """
    absU = np.abs(U)
    th12, th13, th23 = angles_from_absU(absU)

    s12 = np.sin(np.radians(th12)); c12 = np.cos(np.radians(th12))
    s13 = np.sin(np.radians(th13)); c13 = np.cos(np.radians(th13))
    s23 = np.sin(np.radians(th23)); c23 = np.cos(np.radians(th23))

    J = _jarlskog(U)
    denom = (c12 * (c13**2) * c23 * s12 * s13 * s23)

    if denom == 0.0:
        return float("nan")

    sdelta = np.clip(J / denom, -1.0, 1.0)
    delta = np.degrees(np.arcsin(sdelta))
    return float(delta)

def _label_solar_pair(m: np.ndarray) -> tuple:
    """
    Given 3 positive light masses (any order), choose:
      (i1,i2) = solar pair = indices with smallest |Δm^2|
      i3 = remaining index
    Returns (i1,i2,i3) with m[i2] >= m[i1].
    """
    m = np.array(m, dtype=float)
    m2 = m**2

    pairs = [(0,1), (0,2), (1,2)]
    dms = [abs(m2[a] - m2[b]) for (a,b) in pairs]
    a,b = pairs[int(np.argmin(dms))]
    k = [i for i in [0,1,2] if i not in (a,b)][0]

    # order solar so dm21 positive
    if m2[b] < m2[a]:
        a,b = b,a

    return a,b,k

def _order_by_physical_labels(m_light: np.ndarray, U: np.ndarray) -> tuple:
    """
    Reorder columns of U to [ν1, ν2, ν3] where (1,2) is the solar pair.
    Returns (m_ordered, U_ordered, ordering_tag, dm21, dm31, dm32).
    """
    m_light = np.array(m_light, dtype=float)
    U = np.array(U, dtype=complex)

    i1, i2, i3 = _label_solar_pair(m_light)

    perm = [i1, i2, i3]
    m_ord = m_light[perm]
    U_ord = U[:, perm]

    # phase-fix columns after permutation
    U_ord = _phase_fix_columns(U_ord)

    m2 = m_ord**2
    dm21 = float(m2[1] - m2[0])
    dm31 = float(m2[2] - m2[0])
    dm32 = float(m2[2] - m2[1])

    ordering = "NO" if dm31 > 0 else "IO"
    return m_ord, U_ord, ordering, dm21, dm31, dm32

def compute_pmns_physics_grade(
    pipe,
    Y_nu: np.ndarray,
    *,
    use_higgs_vev_in_seesaw: bool = True,
    charged_geom_override=None
):
    """
    Build physics-grade PMNS without touching pipe internals.

    - Neutrino: do 6x6 seesaw from Y_nu
      If use_higgs_vev_in_seesaw: m_D = v * Y_nu (v in GeV from cfg.higgs_vev)
    - Charged leptons: Y_e built from same pipeline but with phases stripped (default)
      Apply small 1–3 Ue extracted by your method.

    Returns dict with angles, deltaCP, splittings, sum_mnu, ordering, plus matrices.
    """
    cfg = pipe.cfg

    # --- charged lepton operator Y_e (strip phases by default) ---
    if charged_geom_override is None:
        geom_e = {k: np.abs(v) for k, v in cfg.geometry_weights.items()}
    else:
        geom_e = charged_geom_override

    Y_e = pipe._build_Y_with_geometry(geometry_weights_override=geom_e)
    Ue, theta_e = pipe.extract_charged_lepton_rotation(Y_e)

    # --- seesaw build ---
    Y = np.array(Y_nu, dtype=complex)
    if use_higgs_vev_in_seesaw:
        mD = cfg.higgs_vev * Y  # GeV if Y is Yukawa-like
        MR = pipe.seesaw_M * np.eye(3)  # treat seesaw_M as GeV
    else:
        mD = Y
        MR = pipe.seesaw_M * np.eye(3)

    zero = np.zeros_like(mD)
    big = np.block([[zero, mD],
                    [mD.conj().T, MR]])

    eigvals, eigvecs = eigh(big)
    eigvals = np.real(eigvals)

    # light first by |eig|
    order = np.argsort(np.abs(eigvals))
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]

    m_light = np.abs(eigvals[:3])
    Unu = eigvecs[:3, :3].copy()

    # normalize columns
    for k in range(3):
        nrm = np.linalg.norm(Unu[:, k])
        if nrm > 0:
            Unu[:, k] /= nrm

    # --- eigenvalue-based labeling ν1,ν2,ν3 + phase fixing ---
    m_ord, Unu_ord, ordering, dm21, dm31, dm32 = _order_by_physical_labels(m_light, Unu)

    # --- PMNS ---
    U_pmns = Ue.conj().T @ Unu_ord
    U_pmns = _phase_fix_columns(U_pmns)

    absU = np.abs(U_pmns)
    th12, th13, th23 = angles_from_absU(absU)
    J = _jarlskog(U_pmns)
    delta = _delta_cp_from_U(U_pmns)

    out = {
        "theta12": float(th12),
        "theta13": float(th13),
        "theta23": float(th23),
        "deltaCP_deg": float(delta),
        "J": float(J),
        "m_light": m_ord.tolist(),
        "sum_mnu": float(np.sum(m_ord)),
        "dm21": float(dm21),
        "dm31": float(dm31),
        "dm32": float(dm32),
        "ordering": ordering,
        "U_pmns": U_pmns,
        "Ue_theta_deg": float(np.degrees(theta_e)),
    }
    return out

def pmns_global_score(
    *,
    theta12, theta13, theta23,
    dm21, dm31,
    # supply your preferred targets (units consistent with your seesaw convention)
    target_angles=(33.4, 8.6, 45.0),
    target_dm=(7.4e-5, 2.5e-3),
    w_angles=(1.0, 6.0, 1.0),
    w_dm=(1.0, 1.0),
):
    """
    A simple combined score:
      score = angle_score + dm_score
    Use units consistently:
      - If you use Higgs VEV and seesaw_M in GeV, dm's will be in GeV^2.
      - If you want eV^2, scale masses before computing dm targets.
    """
    t12, t13, t23 = target_angles
    w12, w13, w23 = w_angles
    s_ang = (w12*(theta12 - t12)**2 + w13*(theta13 - t13)**2 + w23*(theta23 - t23)**2)

    tdm21, tdm31 = target_dm
    wdm21, wdm31 = w_dm
    s_dm = (wdm21*(dm21 - tdm21)**2 + wdm31*(dm31 - tdm31)**2)

    return float(s_ang + s_dm)

# ============================================================
#  Physics-grade 2D scan wrapper
# ============================================================

def run_pmns_2d_scan_physics_grade(
    *,
    phase12_vals_pi,
    phase23_vals_pi,
    eps_vals,
    g_diag=0.894,
    g_off=1.0,
    beta=1.5,
    rel_cut=0.15,
    tol_rel_blocks=0.03,
    seesaw_M=1e6,
    require_3d=True,
    use_higgs_vev_in_seesaw=True,
    top_k=15,
    csv_path="pmns_2d_scan_physics.csv",
    verbose=True
):
    rows = []
    total = len(eps_vals) * len(phase12_vals_pi) * len(phase23_vals_pi)
    n = 0

    for eps in eps_vals:
        for ph12 in phase12_vals_pi:
            for ph23 in phase23_vals_pi:
                n += 1
                cfg = AlignmentV33Config()

                set_geometry_weights_pmns_two_phase(
                    cfg,
                    g_diag=g_diag,
                    g_off=g_off,
                    eps=float(eps),
                    phase12=float(ph12) * np.pi,
                    phase23=float(ph23) * np.pi
                )

                pipe = AlignmentPipeline(
                    cfg,
                    beta=beta,
                    rel_cut=rel_cut,
                    tol_rel_blocks=tol_rel_blocks,
                    seesaw_M=seesaw_M
                )

                out = pipe.run()
                blocks = out["blocks"]
                ok_3d = (blocks == [[0, 1, 2]])

                if require_3d and not ok_3d:
                    rows.append({
                        "eps": float(eps),
                        "phase12_pi": float(ph12),
                        "phase23_pi": float(ph23),
                        "ok_3d": False,
                        "blocks": str(blocks),
                        "theta12": None, "theta13": None, "theta23": None,
                        "deltaCP_deg": None,
                        "dm21": None, "dm31": None, "sum_mnu": None,
                        "ordering": None,
                        "score": None
                    })
                    continue

                phys = compute_pmns_physics_grade(
                    pipe,
                    out["Y"],
                    use_higgs_vev_in_seesaw=use_higgs_vev_in_seesaw
                )

                # score: angles + (optional) splittings (targets must match your units)
                score = pmns_global_score(
                    theta12=phys["theta12"],
                    theta13=phys["theta13"],
                    theta23=phys["theta23"],
                    dm21=phys["dm21"],
                    dm31=phys["dm31"],
                    # Set dm targets to something consistent with your chosen mass units.
                    # If you want to only score angles, set w_dm=(0,0) in pmns_global_score.
                    target_dm=(0.0, 0.0),
                    w_dm=(0.0, 0.0),
                )

                row = {
                    "eps": float(eps),
                    "phase12_pi": float(ph12),
                    "phase23_pi": float(ph23),
                    "ok_3d": True,
                    "blocks": str(blocks),

                    "theta12": phys["theta12"],
                    "theta13": phys["theta13"],
                    "theta23": phys["theta23"],
                    "deltaCP_deg": phys["deltaCP_deg"],
                    "J": phys["J"],

                    "dm21": phys["dm21"],
                    "dm31": phys["dm31"],
                    "dm32": phys["dm32"],
                    "sum_mnu": phys["sum_mnu"],
                    "ordering": phys["ordering"],

                    "Ue_theta_deg": phys["Ue_theta_deg"],
                    "score": float(score),
                }
                rows.append(row)

                if verbose and (n % 25 == 0 or n == total):
                    print(f"[{n}/{total}] eps={eps:.3f} ph12={ph12:.4f}π ph23={ph23:.4f}π "
                          f"θ12={phys['theta12']:.2f} θ13={phys['theta13']:.2f} θ23={phys['theta23']:.2f} "
                          f"δ={phys['deltaCP_deg']:.1f}°  sum={phys['sum_mnu']:.3e}  score={score:.2f}")

    # rank best
    valid = [r for r in rows if r["score"] is not None]
    valid.sort(key=lambda r: r["score"])

    print("\n===== TOP PHYSICS-GRADE CANDIDATES =====")
    for i, r in enumerate(valid[:top_k], start=1):
        print(f"{i:2d}) score={r['score']:.2f} eps={r['eps']:.3f} "
              f"ph12={r['phase12_pi']:.4f}π ph23={r['phase23_pi']:.4f}π "
              f"| θ12={r['theta12']:.2f} θ13={r['theta13']:.2f} θ23={r['theta23']:.2f} "
              f"δ={r['deltaCP_deg']:.1f}° "
              f"dm21={r['dm21']:.3e} dm31={r['dm31']:.3e} sum={r['sum_mnu']:.3e} {r['ordering']}")

    # save csv
    if csv_path and rows:
        import csv
        keys = list(rows[0].keys())
        with open(csv_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            for r in rows:
                w.writerow(r)
        print(f"\nWrote CSV: {csv_path}")

    return rows


# Example call (mirrors your current __main__)
if __name__ == "__main__":
    run_pmns_2d_scan_physics_grade(
        eps_vals=(0.010, 0.011, 0.012, 0.013),
        phase12_vals_pi=np.round(np.linspace(0.072, 0.088, 9), 5),
        phase23_vals_pi=np.round(np.linspace(0.012, 0.032, 11), 5),
        seesaw_M=1e6,
        use_higgs_vev_in_seesaw=True,   # switch to False if you insist Y already has mass units
        verbose=True
    )




"""
RESULTS:
===== TOP 2D CANDIDATES =====
 1) score=0.26  eps=0.014  ph12=0.0895π  ph23=0.0124π  | θ12=33.29 θ13=8.41 θ23=44.85  Ue3=0.1462
 2) score=0.32  eps=0.014  ph12=0.0895π  ph23=0.0120π  | θ12=33.02 θ13=8.46 θ23=44.76  Ue3=0.1471
 3) score=0.44  eps=0.013  ph12=0.0900π  ph23=0.0132π  | θ12=33.84 θ13=8.43 θ23=44.73  Ue3=0.1466
 4) score=0.45  eps=0.013  ph12=0.0910π  ph23=0.0128π  | θ12=33.49 θ13=8.35 θ23=44.77  Ue3=0.1452
 5) score=0.46  eps=0.013  ph12=0.0885π  ph23=0.0128π  | θ12=33.29 θ13=8.33 θ23=44.91  Ue3=0.1449
 6) score=0.49  eps=0.013  ph12=0.0905π  ph23=0.0124π  | θ12=32.90 θ13=8.40 θ23=44.96  Ue3=0.1461
 7) score=0.52  eps=0.014  ph12=0.0905π  ph23=0.0140π  | θ12=33.31 θ13=8.31 θ23=45.03  Ue3=0.1445
 8) score=0.52  eps=0.013  ph12=0.0915π  ph23=0.0156π  | θ12=33.84 θ13=8.38 θ23=44.79  Ue3=0.1458
 9) score=0.55  eps=0.013  ph12=0.0910π  ph23=0.0124π  | θ12=33.40 θ13=8.34 θ23=45.36  Ue3=0.1450
10) score=0.57  eps=0.013  ph12=0.0905π  ph23=0.0124π  | θ12=33.46 θ13=8.29 θ23=45.02  Ue3=0.1443
11) score=0.61  eps=0.014  ph12=0.0915π  ph23=0.0136π  | θ12=33.84 θ13=8.35 θ23=44.79  Ue3=0.1452
12) score=0.61  eps=0.013  ph12=0.0910π  ph23=0.0144π  | θ12=33.10 θ13=8.34 θ23=45.33  Ue3=0.1450
13) score=0.64  eps=0.014  ph12=0.0885π  ph23=0.0132π  | θ12=33.46 θ13=8.30 θ23=44.69  Ue3=0.1443
14) score=0.67  eps=0.014  ph12=0.0905π  ph23=0.0136π  | θ12=33.19 θ13=8.28 θ23=45.15  Ue3=0.1441
15) score=0.72  eps=0.014  ph12=0.0910π  ph23=0.0144π  | θ12=33.28 θ13=8.26 θ23=44.94  Ue3=0.1436
"""