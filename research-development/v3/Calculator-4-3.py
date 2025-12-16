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

def wrap_deg(x):
    return float((x + 360.0) % 360.0)

def circ_diff_deg(a, b):
    # returns signed difference in (-180, 180]
    d = (a - b + 180.0) % 360.0 - 180.0
    return float(d)

def delta_cp_from_U(U_pmns):
    """
    Extract δCP in degrees using:
      J = Im(Ue1 Uμ2 Ue2* Uμ1*)
      J = c12 c23 c13^2 s12 s23 s13 sinδ
    and cosδ from |U_{μ1}|^2 formula.
    Works from moduli + one rephasing-invariant (J).
    """
    U = np.array(U_pmns, dtype=complex)
    absU = np.abs(U)

    th12, th13, th23 = angles_from_absU(absU)
    t12, t13, t23 = np.radians([th12, th13, th23])

    s12, c12 = np.sin(t12), np.cos(t12)
    s13, c13 = np.sin(t13), np.cos(t13)
    s23, c23 = np.sin(t23), np.cos(t23)

    # Jarlskog invariant (rephasing invariant)
    J = float(np.imag(U[0,0] * U[1,1] * np.conj(U[0,1]) * np.conj(U[1,0])))

    denom = c12 * c23 * (c13**2) * s12 * s23 * s13
    if abs(denom) < 1e-15:
        return np.nan, J, th12, th13, th23

    sin_delta = np.clip(J / denom, -1.0, 1.0)

    # cosδ from |U_{μ1}|^2 = s12^2 c23^2 + c12^2 s23^2 s13^2 + 2 s12 c12 c23 s23 s13 cosδ
    Um1_sq = float(absU[1,0]**2)
    denom_c = 2.0 * s12 * c12 * c23 * s23 * s13
    if abs(denom_c) < 1e-15:
        cos_delta = 0.0
    else:
        cos_delta = (Um1_sq - (s12**2)*(c23**2) - (c12**2)*(s23**2)*(s13**2)) / denom_c
        cos_delta = float(np.clip(cos_delta, -1.0, 1.0))

    delta = np.degrees(np.arctan2(sin_delta, cos_delta))
    return wrap_deg(delta), J, th12, th13, th23

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

def dm2_from_masses(light_masses, ordering="NO"):
    """
    Returns (dm21, dm32) in *raw* units (whatever your light_masses carry).
    ordering:
      NO: m1<m2<m3  (assume sorted ascending)
      IO: m3<m1<m2  (assume sorted ascending, then relabel)
    """
    m = np.sort(np.array(light_masses, dtype=float))

    if ordering.upper() == "NO":
        m1, m2, m3 = m[0], m[1], m[2]
    else:  # IO: lightest is m3
        m3, m1, m2 = m[0], m[1], m[2]

    dm21 = m2**2 - m1**2
    dm32 = m3**2 - m2**2
    return float(dm21), float(dm32)

def scale_to_match_dm32(dm21_raw, dm32_raw, dm32_target):
    """
    Choose scale s so that (s^2)*dm32_raw matches dm32_target in magnitude/sign.
    You can enforce sign if you want; here we match magnitude and keep raw sign.
    """
    if abs(dm32_raw) < 1e-30:
        return np.nan, np.nan, np.nan
    s2 = abs(dm32_target) / abs(dm32_raw)  # scale^2
    dm21 = s2 * dm21_raw
    dm32 = s2 * dm32_raw
    s = np.sqrt(s2)
    return float(dm21), float(dm32), float(s)

def unified_lepton_objective(U_pmns, light_masses,
                             targets=None,
                             sigmas=None,
                             ordering="NO",
                             w_angles=1.0, w_delta=1.0, w_dm=1.0):
    """
    Returns: total_score, diagnostics_dict
    """
    if targets is None:
        # PDG Table 14.7 (one of the displayed global fits), Normal Ordering :contentReference[oaicite:2]{index=2}
        targets = dict(
            th12=33.41, th13=8.54, th23=49.1,
            delta=197.0,
            dm21=7.41e-5, dm32=2.437e-3
        )

    if sigmas is None:
        # tune these weights as you like (units: degrees, eV^2)
        sigmas = dict(
            th12=0.8, th13=0.15, th23=1.5,
            delta=25.0,
            dm21=0.25e-5, dm32=0.06e-3
        )

    delta, J, th12, th13, th23 = delta_cp_from_U(U_pmns)

    # dm^2 (scale-fixed by matching |dm32|)
    dm21_raw, dm32_raw = dm2_from_masses(light_masses, ordering=ordering)
    dm21, dm32, s = scale_to_match_dm32(dm21_raw, dm32_raw, targets["dm32"])

    # angle penalties
    ang = (
        ((th12 - targets["th12"]) / sigmas["th12"])**2 +
        ((th13 - targets["th13"]) / sigmas["th13"])**2 +
        ((th23 - targets["th23"]) / sigmas["th23"])**2
    )

    # circular δ penalty
    if np.isnan(delta):
        dcp = 1e6
    else:
        d = circ_diff_deg(delta, targets["delta"])
        dcp = (d / sigmas["delta"])**2

    # dm penalties
    if np.isnan(dm21) or np.isnan(dm32):
        dm = 1e6
    else:
        dm = (
            ((dm21 - targets["dm21"]) / sigmas["dm21"])**2 +
            ((dm32 - targets["dm32"]) / sigmas["dm32"])**2
        )

    total = float(w_angles*ang + w_delta*dcp + w_dm*dm)

    diag = dict(
        ordering=ordering,
        th12=float(th12), th13=float(th13), th23=float(th23),
        delta=float(delta) if not np.isnan(delta) else None,
        J=float(J),
        dm21=float(dm21) if not np.isnan(dm21) else None,
        dm32=float(dm32) if not np.isnan(dm32) else None,
        scale_s=float(s) if not np.isnan(s) else None,
        ang_score=float(ang), dcp_score=float(dcp), dm_score=float(dm),
        total=float(total)
    )
    return total, diag

def mass_splittings(m_light):
    m = np.array(m_light, dtype=float)
    dm21 = m[1]**2 - m[0]**2
    dm31 = m[2]**2 - m[0]**2
    return float(dm21), float(dm31)

def delta_cp_deg(U):
    U = np.array(U, dtype=complex)
    A = np.abs(U)

    s13 = A[0,2]
    c13 = np.sqrt(max(0.0, 1.0 - s13**2))
    s12 = A[0,1] / c13
    c12 = A[0,0] / c13
    s23 = A[1,2] / c13
    c23 = A[2,2] / c13

    # Jarlskog invariant
    J = np.imag(U[0,0]*U[1,1]*np.conj(U[0,1])*np.conj(U[1,0]))
    pref = c12*c23*(c13**2)*s12*s23*s13
    sin_delta = J/pref if abs(pref) > 1e-15 else 0.0

    # cos(delta) from |U_mu1|^2 identity
    Umu1 = A[1,0]
    num = Umu1**2 - (s12**2 * c23**2) - (c12**2 * s23**2 * s13**2)
    den = 2*s12*c12*c23*s23*s13
    cos_delta = num/den if abs(den) > 1e-15 else 1.0

    sin_delta = float(np.clip(sin_delta, -1.0, 1.0))
    cos_delta = float(np.clip(cos_delta, -1.0, 1.0))

    delta = np.arctan2(sin_delta, cos_delta)   # (-pi, pi]
    return float(np.degrees(delta))

def wrap_deg(x):
    return ((x + 180.0) % 360.0) - 180.0

def unified_objective(th12, th13, th23, dm21, dm31, delta_deg,
                      target_angles=(33.4, 8.6, 45.0),
                      r_target=0.03,
                      delta_target_deg=177.0,
                      w_angles=1.0, w_r=5.0, w_delta=1.0, w_hier=50.0):
    t12, t13, t23 = target_angles

    # scale-free ratio (works even if your overall neutrino scale is arbitrary)
    r = dm21 / abs(dm31) if dm31 != 0 else 1e9

    # hierarchy penalty (want NO => dm31 > 0)
    hier_pen = 0.0 if dm31 > 0 else 1.0

    ddelta = wrap_deg(delta_deg - delta_target_deg)

    L_angles = (th12 - t12)**2 + 6.0*(th13 - t13)**2 + (th23 - t23)**2
    L = (w_angles*L_angles
         + w_r*(r - r_target)**2
         + w_delta*(ddelta/30.0)**2
         + w_hier*hier_pen)
    return float(L), float(r)

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
        U = np.array(U, dtype=complex, copy=True)

        ue = np.abs(U[0, :])
        k3 = int(np.argmin(ue))
        rest = [k for k in range(3) if k != k3]
        rest = sorted(rest, key=lambda k: ue[k], reverse=True)

        perm = rest + [k3]
        U = U[:, perm]

        # phase-fix columns so U[0,k] real and >=0
        for k in range(3):
            ph = np.angle(U[0, k])
            U[:, k] *= np.exp(-1j * ph)
            if np.real(U[0, k]) < 0:
                U[:, k] *= -1.0

        return U, perm

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
        Unu, perm = self._canonicalize_pmns(Unu)
        light_masses = light_masses[perm]

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
        U_pmns, perm2 = self._canonicalize_pmns(U_pmns)
        light_masses = light_masses[perm2]

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


# ============================================================
#  Example usage
# ============================================================

def pmns_score(theta12, theta13, theta23,
               target=(33.4, 8.6, 45.0),
               weights=(1.0, 2.0, 1.0)):
    """
    Lower is better. Weighted squared error in degrees.
    By default, θ13 gets heavier weight because it's been the hardest to land.
    """
    t12, t13, t23 = target
    w12, w13, w23 = weights
    return (w12*(theta12 - t12)**2 +
            w13*(theta13 - t13)**2 +
            w23*(theta23 - t23)**2)

def run_pmns_2d_scan(
    *,
    phase12_vals_pi=None,
    phase23_vals_pi=None,
    eps_vals=(0.013,0.010),          # keep fixed unless you want the optional slice
    g_diag=0.894,
    g_off=1.0,
    beta=1.5,
    rel_cut=0.15,
    tol_rel_blocks=0.03,
    seesaw_M=1e6,
    require_3d=True,
    top_k=15,
    csv_path="pmns_2d_scan.csv",
    verbose=True
):
    if phase12_vals_pi is None:
        phase12_vals_pi = np.round(np.linspace(0.070, 0.085, 7), 5)  # 0.070π..0.085π
    if phase23_vals_pi is None:
        phase23_vals_pi = np.round(np.linspace(0.010, 0.030, 9), 5)  # 0.010π..0.030π

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
                    eps=eps,
                    phase12=ph12 * np.pi,
                    phase23=ph23 * np.pi
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
                        "eps": eps, "phase12_pi": ph12, "phase23_pi": ph23,
                        "ok_3d": False, "blocks": str(blocks),
                        "theta12": None, "theta13": None, "theta23": None,
                        "Ue3": None, "score": None
                    })
                    continue

                m_light, U_pmns, diag = pipe.apply_full_seesaw_3x3(out["Y"], include_charged_lepton=True)
                score_total, diag2 = unified_lepton_objective(U_pmns, m_light, ordering="NO")
                # if you want to allow both orderings and pick best:
                score_NO, dNO = unified_lepton_objective(U_pmns, m_light, ordering="NO")
                score_IO, dIO = unified_lepton_objective(U_pmns, m_light, ordering="IO")
                if score_IO < score_NO:
                    score_total, diag2 = score_IO, dIO
                else:
                    score_total, diag2 = score_NO, dNO

                print(f"θe (charged lepton) = {diag.get('theta_e_deg', None)}°")

                absU = np.abs(U_pmns)
                th12, th13, th23 = angles_from_absU(absU)
                Ue3 = float(absU[0, 2])
                dm21, dm31 = mass_splittings(m_light)
                dcp = delta_cp_deg(U_pmns)
                obj, r = unified_objective(th12, th13, th23, dm21, dm31, dcp)
                # scoring: emphasize θ13 without ignoring θ12, θ23
                score = float(pmns_score(th12, th13, th23,
                                         target=(33.4, 8.6, 45.0),
                                         weights=(1.0, 6.0, 1.0)))

                rows.append({
                    "eps": float(eps),
                    "phase12_pi": float(ph12),
                    "phase23_pi": float(ph23),
                    "ok_3d": True,
                    "blocks": str(blocks),
                    "theta12": float(th12),
                    "theta13": float(th13),
                    "theta23": float(th23),
                    "Ue3": float(Ue3),
                    "score": float(score),
                    "deltaCP_deg": diag2["delta"],
                    "dm21_eV2": diag2["dm21"],
                    "dm32_eV2": diag2["dm32"],
                    "score_total": score_total,
                    "dm21": dm21, "dm31": dm31, "r": r,
                    "delta_cp_deg": dcp,
                    "objective": obj
                })

                if verbose and (n % 25 == 0 or n == total):
                    print(f"[{n}/{total}] eps={eps:.3f} ph12={ph12:.4f}π ph23={ph23:.4f}π "
                          f"θ12={th12:.2f} θ13={th13:.2f} θ23={th23:.2f} score={score:.2f}")

    # rank best
    valid = [r for r in rows if r["score"] is not None]
    valid.sort(key=lambda r: r["score"])

    print("\n===== TOP 2D CANDIDATES =====")
    for i, r in enumerate(valid[:top_k], start=1):
        print(f"{i:2d}) score={r['score']:.2f}  eps={r['eps']:.3f}  "
              f"ph12={r['phase12_pi']:.4f}π  ph23={r['phase23_pi']:.4f}π  "
              f"| θ12={r['theta12']:.2f} θ13={r['theta13']:.2f} θ23={r['theta23']:.2f}  Ue3={r['Ue3']:.4f}")

    # save csv
    if csv_path:
        import csv
        keys = list(rows[0].keys()) if rows else []
        with open(csv_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            for r in rows:
                w.writerow(r)
        print(f"\nWrote CSV: {csv_path}")

    return rows


if __name__ == "__main__":
    run_pmns_2d_scan(
        eps_vals=(0.0125, 0.0130, 0.0135, 0.0140),
        phase12_vals_pi=np.round(np.linspace(0.088, 0.092, 9), 5),
        phase23_vals_pi=np.round(np.linspace(0.012, 0.016, 11), 5),
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