#!/usr/bin/env python3
import json
import csv
import argparse
from dataclasses import dataclass
from typing import Optional, Any, List, Dict, Tuple

import numpy as np
from numpy.linalg import eigh
from scipy.linalg import expm, svd


# ============================================================
#  DATA
# ============================================================

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
    dm32: Optional[float] = None

    delta_cp_deg: Optional[float] = None
    jarlskog: Optional[float] = None

    score: Optional[float] = None
    note: str = ""

GEV_TO_EV = 1e9

def eV(m_GeV: np.ndarray) -> np.ndarray:
    return np.asarray(m_GeV, dtype=float) * GEV_TO_EV
def MR_texture(M: float, r12: float = 0.20, phi12: float = 0.5*np.pi):
    MR = M * np.eye(3, dtype=complex)
    z = M * r12 * np.exp(1j * phi12)
    MR[0, 1] = z
    MR[1, 0] = z
    return MR

# ============================================================
#  Takagi + invariants
# ============================================================

def takagi(M: np.ndarray, atol: float = 1e-12) -> Tuple[np.ndarray, np.ndarray]:
    """
    Takagi factorization for complex symmetric M:
        U^T M U = diag(m_i),  m_i >= 0, U unitary
    Returns: (m, U) with m >= 0
    """
    M = np.asarray(M, dtype=complex)
    if not np.allclose(M, M.T, atol=atol):
        max_asym = np.max(np.abs(M - M.T))
        raise ValueError(f"Takagi requires M == M.T (max |M-M^T| = {max_asym:.3e}).")

    U, s, Vh = svd(M)

    # phase-fix columns of U so diag(U^T M U) becomes real nonnegative
    T = U.T @ M @ U
    ph = np.angle(np.diag(T))
    U = U @ np.diag(np.exp(-0.5j * ph))

    D = np.array(np.real(np.diag(U.T @ M @ U)), copy=True)
    D[D < 0] = 0.0  # clamp tiny negatives (numerical)

    return D, U


def dm2_from_masses(m: np.ndarray) -> Tuple[float, float, float]:
    m = np.asarray(m, dtype=float)
    m2 = m**2
    dm21 = float(m2[1] - m2[0])
    dm31 = float(m2[2] - m2[0])
    dm32 = float(m2[2] - m2[1])
    return dm21, dm31, dm32


def delta_cp_from_U(U: np.ndarray) -> Tuple[Optional[float], Optional[float]]:
    """
    Extract delta_CP using a rephasing-invariant quartet + PDG relations.
    Returns (delta_rad, J). delta in [-pi, pi].
    """
    U = np.asarray(U, dtype=complex)

    s13 = float(np.clip(np.abs(U[0, 2]), 0.0, 1.0))
    c13 = float(np.sqrt(max(0.0, 1.0 - s13*s13)))
    if c13 < 1e-12:
        return None, None

    s12 = float(np.clip(np.abs(U[0, 1]) / c13, 0.0, 1.0))
    s23 = float(np.clip(np.abs(U[1, 2]) / c13, 0.0, 1.0))

    c12 = float(np.sqrt(max(0.0, 1.0 - s12*s12)))
    c23 = float(np.sqrt(max(0.0, 1.0 - s23*s23)))

    P = U[0, 0] * U[1, 1] * np.conj(U[0, 1]) * np.conj(U[1, 0])
    J = float(np.imag(P))
    R = float(np.real(P))

    denom = (c12*s12*c23*s23*(c13**2)*s13)
    if abs(denom) < 1e-15:
        return None, J

    sin_delta = float(np.clip(J / denom, -1.0, 1.0))

    # cosδ from Re(P) = A + B cosδ
    cos2_12 = c12*c12 - s12*s12
    A = -(c12*c12*s12*s12*(c13**2)) * (c23*c23 - s23*s23*s13*s13)
    B = -(c12*s12*(c13**2)) * (cos2_12 * c23*s23*s13)

    if abs(B) < 1e-15:
        return float(np.arcsin(sin_delta)), J

    cos_delta = float(np.clip((R - A) / B, -1.0, 1.0))
    delta = float(np.arctan2(sin_delta, cos_delta))
    return delta, J


# ============================================================
#  Phase label helpers
# ============================================================

def block_sizes(blocks) -> List[int]:
    return sorted([len(b) for b in blocks])

def regime_label_from_blocks(blocks) -> str:
    s = block_sizes(blocks)
    if s == [1, 1, 1]: return "split"
    if s == [1, 2]:     return "1⊕2"
    if s == [3]:        return "3D"
    return f"other:{s}"


# ============================================================
#  Geometry weights helper
# ============================================================

def set_geometry_weights_pmns_two_phase(
    cfg,
    g_diag=0.894,
    g_off=1.00,
    eps=0.012,
    phase12=0.075*np.pi,
    phase23=0.010*np.pi
):
    base = {
        (0,0): 1.00,
        (0,1): 0.92 * np.exp(1j * phase12),
        (0,2): 0.85,
        (1,0): 0.95 * np.exp(-1j * phase12),
        (1,1): 1.10,
        (1,2): 0.88 * np.exp(1j * phase23),
        (2,0): 0.80,
        (2,1): 0.90 * np.exp(-1j * phase23),
        (2,2): 1.05,
    }

    diag = {(0,0), (1,1), (2,2)}
    gen1_mix = {(0,1), (1,0)}

    for g in base:
        base[g] *= (g_diag if g in diag else g_off)
        if g in gen1_mix:
            base[g] *= (1.0 + eps)

    cfg.geometry_weights = base


# ============================================================
#  Angles helper
# ============================================================

def angles_from_absU(absU: np.ndarray) -> Tuple[float, float, float]:
    s13 = float(absU[0, 2])
    c13 = float(np.sqrt(max(0.0, 1.0 - s13**2)))

    s12 = float(absU[0, 1] / c13) if c13 > 0 else float("nan")
    s23 = float(absU[1, 2] / c13) if c13 > 0 else float("nan")

    th12 = float(np.degrees(np.arcsin(np.clip(s12, -1, 1))))
    th13 = float(np.degrees(np.arcsin(np.clip(s13, -1, 1))))
    th23 = float(np.degrees(np.arcsin(np.clip(s23, -1, 1))))
    return th12, th13, th23


# ============================================================
#  CONFIG
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
    def add_g(a, b): return ((a[0]+b[0]) % 3, (a[1]+b[1]) % 3)
    @staticmethod
    def sub_g(a, b): return ((a[0]-b[0]) % 3, (a[1]-b[1]) % 3)
    @staticmethod
    def chi(g, p, q):
        i, j = g
        return np.exp(2j*np.pi*(p*i+q*j)/3.0)

    def build_S(self, triads, shear=0.0):
        cfg = self.cfg
        G = cfg.group_elements
        S = np.zeros((3, 9), dtype=complex)

        for i, triad in enumerate(triads):
            p, q = cfg.compression_characters[i]
            for k, idx in enumerate(triad):
                g = G[idx]
                ph = self.chi(g, p, q)
                if (i == 0) and (k > 0) and (shear != 0.0):
                    ph *= np.exp(1j * shear)
                S[i, idx] = ph / np.sqrt(3)
        return S

    # ---------- charged-lepton rotation ----------
    @staticmethod
    def _rotation_matrix_13(theta_rad):
        c = np.cos(theta_rad)
        s = np.sin(theta_rad)
        return np.array([
            [ c, 0.0,  s],
            [0.0, 1.0, 0.0],
            [-s, 0.0,  c],
        ], dtype=complex)

    def charged_lepton_left_rotation(self, Y_e: np.ndarray):
        """
        Correct LH rotation: diagonalize H_e = m_e m_e^\dagger.
        Returns Ue (LH), plus diagnostic masses (GeV) for ordering sanity.
        """
        me = self.cfg.higgs_vev * np.asarray(Y_e, dtype=complex)  # GeV
        He = me @ me.conj().T  # Hermitian

        evals, Ue = eigh(He)  # columns: eigenvectors
        order = np.argsort(np.real(evals))  # e, mu, tau (ascending)
        evals = np.real(evals[order])
        Ue = Ue[:, order]

        # Unphysical column phase-fix (optional but stabilizes scans)
        Ue = self._phase_fix_columns(Ue)

        m_e = np.sqrt(np.clip(evals, 0.0, None))  # GeV-ish (model units)
        diag = {"m_e_like_GeV": [float(x) for x in m_e.tolist()]}
        return Ue, diag

    def extract_charged_lepton_rotation(self, Y_e):
        evals, U_full = eigh(Y_e)
        idx_tau = int(np.argmax(np.abs(np.real(evals))))
        v_tau = U_full[:, idx_tau]

        s = float(np.clip(np.abs(v_tau[0]), 0.0, 1.0))
        theta = float(np.arcsin(s))

        sign = 1.0
        if np.real(v_tau[0]) * np.real(v_tau[2]) < 0:
            sign = -1.0
        theta *= sign

        theta = float(np.clip(theta, -0.10, 0.10))
        Ue = self._rotation_matrix_13(theta)
        return Ue, theta

    # ---------- build K ----------
    def build_kernel(self):
        cfg = self.cfg
        G = cfg.group_elements
        n = len(G)
        K = np.zeros((n, n), dtype=complex)

        for a, g in enumerate(G):
            for b, h in enumerate(G):
                F = sum(
                    w * self.chi(self.sub_g(g, h), p, q)
                    for (p, q, w) in cfg.kernel_characters
                )
                alpha_g = cfg.geometry_weights[g]
                alpha_h = cfg.geometry_weights[h]

                dist = min(abs(g[0]-h[0]), 3-abs(g[0]-h[0])) + \
                       min(abs(g[1]-h[1]), 3-abs(g[1]-h[1]))
                W = np.exp(-cfg.damping_strength * dist)

                K[a, b] = alpha_g * F * np.conj(alpha_h) * W

        return 0.5 * (K + K.conj().T)

    # ---------- triads ----------
    def build_triads(self):
        cfg = self.cfg
        index = {g: i for i, g in enumerate(cfg.group_elements)}
        triads = []
        for s in cfg.triad_shifts:
            triads.append([index[self.add_g(h, s)] for h in cfg.subgroup_H])
        return triads

    # ---------- operators ----------
    def misalignment_flow(self, K): return expm(-self.beta * K)

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
    def effective_yukawa(K_like, S): return S @ K_like @ S.conj().T

    # ---------- block diagnosis ----------
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

    def _phase_fix_columns(self, U):
        """
        Keep column order (mass order), only fix unphysical column phases so U[0,k] is real >= 0.
        """
        U = np.array(U, dtype=complex, copy=True)
        for k in range(U.shape[1]):
            ph = np.angle(U[0, k])
            U[:, k] *= np.exp(-1j * ph)
            if np.real(U[0, k]) < 0:
                U[:, k] *= -1.0
        return U

    # ---------- conditional seesaw (block-level) ----------
    @staticmethod
    def _seesaw_light_eigs_for_block_majorana(Y, block, M, v=174.0):
        """
        Effective seesaw on a block (dimensionless Yukawa):
          mnu_b = - (v Yb) (M I)^{-1} (v Yb)^T
        Takagi -> light masses (GeV).
        """
        Yb = np.asarray(Y[np.ix_(block, block)], dtype=complex)
        mDb = v * Yb
        invMR = np.linalg.inv(M * np.eye(len(block), dtype=complex))
        mnu_b = - mDb @ invMR @ mDb.T
        mnu_b = 0.5 * (mnu_b + mnu_b.T)
        m, _U = takagi(mnu_b)
        m = np.sort(m)
        return m  # GeV

    def apply_conditional_seesaw(self, Y, blocks):
        dirac = []
        majorana = []

        y_eigs = eigh(Y)[0]
        for block in blocks:
            if len(block) == 1:
                i = block[0]
                dirac.append((block, np.array([float(np.real(y_eigs[i]))])))
            else:
                light = self._seesaw_light_eigs_for_block_majorana(Y, block, self.seesaw_M)
                majorana.append((block, light))
        return dirac, majorana

    # ---------- canonicalize ----------
    def _canonicalize_pmns(self, U, return_perm=False):
        U = np.array(U, dtype=complex, copy=True)
        ue = np.abs(U[0, :])

        k3 = int(np.argmin(ue))
        rest = [k for k in range(3) if k != k3]
        rest = sorted(rest, key=lambda k: ue[k], reverse=True)

        perm = rest + [k3]
        U = U[:, perm]

        for k in range(3):
            ph = np.angle(U[0, k])
            U[:, k] *= np.exp(-1j * ph)
            if np.real(U[0, k]) < 0:
                U[:, k] *= -1.0

        return (U, perm) if return_perm else U

    def _build_Y_with_geometry(self, geometry_weights_override=None):
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
            S = self.build_S(triads, shear=self.compression_shear)
            Y = self.effective_yukawa(K_proj, S)
            return Y
        finally:
            cfg.geometry_weights = old_geom

    # ---------- full 3x3 Majorana seesaw via Takagi ----------
    def apply_full_seesaw_3x3(
            self,
            Y: np.ndarray,
            include_charged_lepton: bool = True,
            geometry_weights_charged=None,
            MR_matrix: Optional[np.ndarray] = None,
            output_units: str = "eV",  # "eV" or "GeV"
    ):
        """
        Physically-scaled type-I seesaw (effective 3x3):
          m_D = v Y   (GeV)
          m_nu = - m_D M_R^{-1} m_D^T  (GeV), complex symmetric
        Then Takagi(m_nu) -> (m_i >= 0, U_nu).

        PMNS = Ue^\dagger U_nu where Ue is LH charged-lepton rotation from m_e m_e^\dagger.
        """
        Y = np.asarray(Y, dtype=complex)
        assert Y.shape == (3, 3)

        diagnostics: Dict[str, Any] = {}
        v = float(self.cfg.higgs_vev)

        # Build MR
        if MR_matrix is None:
            MR = self.seesaw_M * np.eye(3, dtype=complex)  # GeV
            diagnostics["MR"] = "M * I"
        else:
            MR = np.asarray(MR_matrix, dtype=complex)
            diagnostics["MR"] = "custom"
        diagnostics["seesaw_M_GeV"] = float(self.seesaw_M)

        # Effective light Majorana mass
        mD = v * Y  # GeV
        invMR = np.linalg.inv(MR)
        mnu = - mD @ invMR @ mD.T  # GeV, symmetric (should be)
        # enforce symmetry numerically
        mnu = 0.5 * (mnu + mnu.T)

        # Takagi diagonalization
        m_light_GeV, Unu = takagi(mnu)
        order = np.argsort(m_light_GeV)  # m1<=m2<=m3
        m_light_GeV = m_light_GeV[order]
        Unu = Unu[:, order]
        Unu = self._phase_fix_columns(Unu)

        # Charged-lepton LH rotation
        if include_charged_lepton:
            if geometry_weights_charged is None:
                # stable choice: strip phases from geometry for Ye
                geom_e = {k: np.abs(v) for k, v in self.cfg.geometry_weights.items()}
                diagnostics["Y_e_geometry"] = "abs(phased_geometry)"
            else:
                geom_e = geometry_weights_charged
                diagnostics["Y_e_geometry"] = "override"

            Y_e = self._build_Y_with_geometry(geometry_weights_override=geom_e)
            Ue, diag_e = self.charged_lepton_left_rotation(Y_e)
            diagnostics.update(diag_e)

            U_pmns = Ue.conj().T @ Unu
            U_pmns = self._phase_fix_columns(U_pmns)
        else:
            U_pmns = Unu

        # Output units
        if output_units.lower() == "ev":
            m_light = eV(m_light_GeV)
            diagnostics["m_units"] = "eV"
        else:
            m_light = m_light_GeV
            diagnostics["m_units"] = "GeV"

        return m_light, U_pmns, diagnostics

    # ---------- run ----------
    def run(self):
        cfg = self.cfg

        self.triads = self.build_triads()
        self.K = self.build_kernel()
        self.S = self.build_S(self.triads, shear=self.compression_shear)

        self.K_flow = self.misalignment_flow(self.K)
        self.P_C360, self.kept_indices = self.emergent_C360_projector(self.K)
        self.K_proj = self.P_C360 @ self.K_flow @ self.P_C360

        self.Y = self.effective_yukawa(self.K_proj, self.S)

        y_vals, U = eigh(self.Y)
        masses = np.abs(y_vals) * cfg.higgs_vev

        blocks, y_sorted = self.harmonic_blocks_on_Y(self.Y)
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
#  Scoring + observable packing
# ============================================================

def pmns_score(theta12, theta13, theta23,
               target=(33.4, 8.6, 45.0),
               weights=(1.0, 2.0, 1.0)) -> float:
    t12, t13, t23 = target
    w12, w13, w23 = weights
    return float(w12*(theta12 - t12)**2 +
                 w13*(theta13 - t13)**2 +
                 w23*(theta23 - t23)**2)


def compute_obs_from_pmns(m_light: np.ndarray, U_pmns: np.ndarray) -> Dict[str, Any]:
    absU = np.abs(U_pmns)
    th12, th13, th23 = angles_from_absU(absU)

    delta_rad, J = delta_cp_from_U(U_pmns)
    delta_deg = None if delta_rad is None else float(np.degrees(delta_rad))

    dm21, dm31, dm32 = dm2_from_masses(m_light)

    return {
        "theta12": float(th12),
        "theta13": float(th13),
        "theta23": float(th23),
        "delta_cp_deg": delta_deg,
        "jarlskog": None if J is None else float(J),
        "ue3": float(absU[0, 2]),
        "um3": float(absU[1, 2]),
        "ut3": float(absU[2, 2]),
        "m_light": [float(x) for x in m_light.tolist()],
        "dm21": float(dm21),
        "dm31": float(dm31),
        "dm32": float(dm32),
    }


def write_csv(rows: List[Dict[str, Any]], path: str):
    if not rows:
        print(f"[write_csv] no rows to write: {path}")
        return
    keys = sorted({k for r in rows for k in r.keys()})
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"Wrote CSV: {path}")


def write_json(obj: Any, path: str):
    with open(path, "w") as f:
        json.dump(obj, f, indent=2)
    print(f"Wrote JSON: {path}")


# ============================================================
#  Phase-map + robustness (drop-in upgrades)
# ============================================================

def phase_scan_2d(
    make_pipe,
    xs: np.ndarray,
    ys: np.ndarray,
    set_knobs,
    *,
    tol_rel_blocks: Optional[float] = None,
    rel_cut: Optional[float] = None,
    include_obs: bool = True
) -> Tuple[List[List[str]], List[Dict[str, Any]]]:
    labels = [["" for _ in xs] for _ in ys]
    rows: List[Dict[str, Any]] = []

    for iy, y in enumerate(ys):
        for ix, x in enumerate(xs):
            pipe: AlignmentPipeline = make_pipe()
            if tol_rel_blocks is not None:
                pipe.tol_rel_blocks = float(tol_rel_blocks)
            if rel_cut is not None:
                pipe.rel_cut = float(rel_cut)

            set_knobs(pipe, float(x), float(y))
            out = pipe.run()

            blocks = out["blocks"]
            label = regime_label_from_blocks(blocks)
            labels[iy][ix] = label

            row = {
                "x": float(x),
                "y": float(y),
                "label": label,
                "blocks": str(blocks),
                "beta": float(pipe.beta),
                "rel_cut": float(pipe.rel_cut),
                "tol_rel_blocks": float(pipe.tol_rel_blocks),
                "seesaw_M": float(pipe.seesaw_M),
            }

            if include_obs:
                m_light, U_pmns, diag = pipe.apply_full_seesaw_3x3(out["Y"], include_charged_lepton=True)
                obs = compute_obs_from_pmns(m_light, U_pmns)
                row.update(obs)
                row["theta_e_deg"] = float(diag.get("theta_e_deg", 0.0))

            rows.append(row)

    return labels, rows


def robustness_suite(
    make_pipe,
    xs: np.ndarray,
    ys: np.ndarray,
    set_knobs,
    base_tol: float,
    base_relcut: float,
    *,
    include_obs: bool = False
) -> List[Dict[str, Any]]:
    sweeps: List[Dict[str, Any]] = []
    for tol in [0.5*base_tol, base_tol, 1.5*base_tol]:
        for rc in [0.7*base_relcut, base_relcut, 1.3*base_relcut]:
            labels, rows = phase_scan_2d(
                make_pipe, xs, ys, set_knobs,
                tol_rel_blocks=float(tol),
                rel_cut=float(rc),
                include_obs=include_obs
            )
            sweeps.append({
                "tol_rel_blocks": float(tol),
                "rel_cut": float(rc),
                "labels": labels,
                "rows": rows
            })
    return sweeps


# ============================================================
#  Original PMNS scan, upgraded with observables + fixes
# ============================================================

def run_pmns_2d_scan(
    *,
    phase12_vals_pi=None,
    phase23_vals_pi=None,
    eps_vals=(0.013, 0.010),
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
        phase12_vals_pi = np.round(np.linspace(0.070, 0.085, 7), 5)
    if phase23_vals_pi is None:
        phase23_vals_pi = np.round(np.linspace(0.010, 0.030, 9), 5)

    rows: List[Dict[str, Any]] = []
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
                label = regime_label_from_blocks(blocks)
                ok_3d = (label == "3D")

                row: Dict[str, Any] = {
                    "eps": float(eps),
                    "phase12_pi": float(ph12),
                    "phase23_pi": float(ph23),
                    "label": label,
                    "ok_3d": bool(ok_3d),
                    "blocks": str(blocks),
                    "beta": float(beta),
                    "rel_cut": float(rel_cut),
                    "tol_rel_blocks": float(tol_rel_blocks),
                    "seesaw_M": float(seesaw_M),
                }

                if require_3d and not ok_3d:
                    rows.append(row)
                    continue

                m_light, U_pmns, diag = pipe.apply_full_seesaw_3x3(out["Y"], include_charged_lepton=True)
                obs = compute_obs_from_pmns(m_light, U_pmns)
                row.update(obs)
                row["theta_e_deg"] = float(diag.get("theta_e_deg", 0.0))

                # score (θ13 emphasized)
                row["score"] = float(pmns_score(
                    row["theta12"], row["theta13"], row["theta23"],
                    target=(33.4, 8.6, 45.0),
                    weights=(1.0, 6.0, 1.0)
                ))

                rows.append(row)

                if verbose and (n % 25 == 0 or n == total):
                    print(f"[{n}/{total}] eps={eps:.4f} ph12={ph12:.5f}π ph23={ph23:.5f}π "
                          f"label={label} θ12={row['theta12']:.2f} θ13={row['theta13']:.2f} θ23={row['theta23']:.2f} "
                          f"δCP={row['delta_cp_deg'] if row['delta_cp_deg'] is not None else 'NA'} "
                          f"score={row.get('score', None)}")

    valid = [r for r in rows if r.get("score") is not None]
    valid.sort(key=lambda r: r["score"])

    print("\n===== TOP CANDIDATES =====")
    for i, r in enumerate(valid[:top_k], start=1):
        print(f"{i:2d}) score={r['score']:.2f}  eps={r['eps']:.4f}  "
              f"ph12={r['phase12_pi']:.5f}π  ph23={r['phase23_pi']:.5f}π  "
              f"θ12={r['theta12']:.2f} θ13={r['theta13']:.2f} θ23={r['theta23']:.2f}  "
              f"|Ue3|={r['ue3']:.4f}  δCP={r['delta_cp_deg'] if r['delta_cp_deg'] is not None else 'NA'}  "
              f"dm21={r['dm21']:.3e} dm31={r['dm31']:.3e}")

    if csv_path:
        write_csv(rows, csv_path)

    return rows

def run_eps_gdiag_scan(
    *,
    eps_vals=np.round(np.linspace(0.0, 0.08, 17), 5),
    g_diag_vals=np.round(np.linspace(0.55, 1.15, 25), 5),
    phase12_pi=0.0880,
    phase23_pi=0.0160,
    g_off=1.0,
    beta=1.5,
    rel_cut=0.15,
    tol_rel_blocks=0.03,
    seesaw_M=1e6,
    require_3d=False,
    csv_path="eps_gdiag_scan.csv",
    verbose=True
):
    rows = []
    total = len(eps_vals) * len(g_diag_vals)
    n = 0

    for eps in eps_vals:
        for g_diag in g_diag_vals:
            n += 1
            cfg = AlignmentV33Config()
            set_geometry_weights_pmns_two_phase(
                cfg,
                g_diag=float(g_diag),
                g_off=float(g_off),
                eps=float(eps),
                phase12=float(phase12_pi) * np.pi,
                phase23=float(phase23_pi) * np.pi,
            )

            pipe = AlignmentPipeline(cfg, beta=beta, rel_cut=rel_cut, tol_rel_blocks=tol_rel_blocks, seesaw_M=seesaw_M)
            out = pipe.run()

            label = regime_label_from_blocks(out["blocks"])
            ok_3d = (label == "3D")

            row = {
                "eps": float(eps),
                "g_diag": float(g_diag),
                "g_off": float(g_off),
                "phase12_pi": float(phase12_pi),
                "phase23_pi": float(phase23_pi),
                "label": label,
                "blocks": str(out["blocks"]),
                "ok_3d": bool(ok_3d),
            }

            if require_3d and not ok_3d:
                rows.append(row)
                continue

            m_light, U_pmns, diag = pipe.apply_full_seesaw_3x3(out["Y"], include_charged_lepton=True)
            obs = compute_obs_from_pmns(m_light, U_pmns)

            row.update(obs)
            row["theta_e_deg"] = float(diag.get("theta_e_deg", 0.0))

            row["score"] = float(pmns_score(
                row["theta12"], row["theta13"], row["theta23"],
                target=(33.4, 8.6, 45.0),
                weights=(1.0, 6.0, 1.0)
            ))

            rows.append(row)

            if verbose and (n % 50 == 0 or n == total):
                print(f"[{n}/{total}] eps={eps:.4f} g_diag={g_diag:.4f} label={label} "
                      f"θ12={row['theta12']:.2f} θ13={row['theta13']:.2f} θ23={row['theta23']:.2f} score={row['score']:.2f}")

    rows_sorted = [r for r in rows if r.get("score") is not None]
    rows_sorted.sort(key=lambda r: r["score"])

    print("\n===== TOP eps × g_diag CANDIDATES =====")
    for i, r in enumerate(rows_sorted[:15], start=1):
        print(f"{i:2d}) score={r['score']:.2f} label={r['label']} eps={r['eps']:.4f} g_diag={r['g_diag']:.4f} "
              f"θ12={r['theta12']:.2f} θ13={r['theta13']:.2f} θ23={r['theta23']:.2f} |Ue3|={r['ue3']:.4f}")

    if csv_path:
        write_csv(rows, csv_path)

    return rows

def run_phase_map_gdiag_beta(
    *,
    g_diag_vals=np.round(np.linspace(0.6, 1.2, 31), 5),
    beta_vals=np.round(np.linspace(0.5, 3.0, 26), 5),
    eps=0.013,
    g_off=1.0,
    phase12_pi=0.0880,
    phase23_pi=0.0160,
    rel_cut=0.15,
    tol_rel_blocks=0.03,
    seesaw_M=1e6,
    out_prefix="phase_gdiag_beta",
):
    xs = np.array(g_diag_vals, dtype=float)
    ys = np.array(beta_vals, dtype=float)

    def make_pipe():
        cfg = AlignmentV33Config()
        return AlignmentPipeline(cfg, beta=1.5, rel_cut=rel_cut, tol_rel_blocks=tol_rel_blocks, seesaw_M=seesaw_M)

    def set_knobs(pipe: AlignmentPipeline, x: float, y: float):
        pipe.beta = float(y)
        set_geometry_weights_pmns_two_phase(
            pipe.cfg,
            g_diag=float(x),
            g_off=float(g_off),
            eps=float(eps),
            phase12=float(phase12_pi) * np.pi,
            phase23=float(phase23_pi) * np.pi
        )

    labels, rows = phase_scan_2d(make_pipe, xs, ys, set_knobs, include_obs=False)
    write_csv(rows, f"{out_prefix}.csv")
    write_json({
        "x_axis": {"name": "g_diag", "values": [float(v) for v in xs.tolist()]},
        "y_axis": {"name": "beta", "values": [float(v) for v in ys.tolist()]},
        "labels": labels,
        "params": {
            "eps": float(eps),
            "g_off": float(g_off),
            "phase12_pi": float(phase12_pi),
            "phase23_pi": float(phase23_pi),
            "rel_cut": float(rel_cut),
            "tol_rel_blocks": float(tol_rel_blocks),
            "seesaw_M": float(seesaw_M),
        }
    }, f"{out_prefix}.json")

    return labels, rows

# ============================================================
#  New: Phase map driver for (phase12_pi, phase23_pi)
# ============================================================

def run_phase_map_over_phases(
    *,
    phase12_vals_pi: np.ndarray,
    phase23_vals_pi: np.ndarray,
    eps: float,
    g_diag: float,
    g_off: float,
    beta: float,
    rel_cut: float,
    tol_rel_blocks: float,
    seesaw_M: float,
    out_prefix: str = "phase_map",
    include_obs: bool = True,
    do_robustness: bool = True
):
    xs = np.array(phase12_vals_pi, dtype=float)
    ys = np.array(phase23_vals_pi, dtype=float)

    def make_pipe():
        cfg = AlignmentV33Config()
        pipe = AlignmentPipeline(cfg, beta=beta, rel_cut=rel_cut, tol_rel_blocks=tol_rel_blocks, seesaw_M=seesaw_M)
        return pipe

    def set_knobs(pipe: AlignmentPipeline, x: float, y: float):
        set_geometry_weights_pmns_two_phase(
            pipe.cfg,
            g_diag=g_diag, g_off=g_off, eps=eps,
            phase12=x * np.pi,
            phase23=y * np.pi
        )

    labels, rows = phase_scan_2d(
        make_pipe, xs, ys, set_knobs,
        tol_rel_blocks=tol_rel_blocks,
        rel_cut=rel_cut,
        include_obs=include_obs
    )

    write_csv(rows, f"{out_prefix}.csv")
    write_json({
        "x_axis": {"name": "phase12_pi", "values": [float(v) for v in xs.tolist()]},
        "y_axis": {"name": "phase23_pi", "values": [float(v) for v in ys.tolist()]},
        "labels": labels,
        "params": {
            "eps": float(eps),
            "g_diag": float(g_diag),
            "g_off": float(g_off),
            "beta": float(beta),
            "rel_cut": float(rel_cut),
            "tol_rel_blocks": float(tol_rel_blocks),
            "seesaw_M": float(seesaw_M),
        }
    }, f"{out_prefix}.json")

    if do_robustness:
        sweeps = robustness_suite(
            make_pipe, xs, ys, set_knobs,
            base_tol=tol_rel_blocks,
            base_relcut=rel_cut,
            include_obs=False  # keep it light; turn on if you want heavy robustness
        )
        write_json(sweeps, f"{out_prefix}.robustness.json")

    return labels, rows


# ============================================================
#  MAIN
# ============================================================

def main():
    ap = argparse.ArgumentParser(description="Alignment v3.3 phase scans + PMNS scans (Takagi seesaw upgraded).")
    ap.add_argument("--mode", choices=["pmns_scan", "phase_map"], default="pmns_scan")

    # shared
    ap.add_argument("--beta", type=float, default=1.5)
    ap.add_argument("--rel_cut", type=float, default=0.15)
    ap.add_argument("--tol_rel_blocks", type=float, default=0.03)
    ap.add_argument("--seesaw_M", type=float, default=1e6)
    ap.add_argument("--g_diag", type=float, default=0.894)
    ap.add_argument("--g_off", type=float, default=1.0)
    ap.add_argument("--eps", type=float, default=0.013)

    # pmns scan ranges
    ap.add_argument("--ph12_lo", type=float, default=0.088)
    ap.add_argument("--ph12_hi", type=float, default=0.092)
    ap.add_argument("--ph12_n", type=int, default=9)

    ap.add_argument("--ph23_lo", type=float, default=0.012)
    ap.add_argument("--ph23_hi", type=float, default=0.016)
    ap.add_argument("--ph23_n", type=int, default=11)

    ap.add_argument("--csv", type=str, default=None)
    ap.add_argument("--out_prefix", type=str, default="phase_map")
    ap.add_argument("--require_3d", action="store_true", default=False)
    ap.add_argument("--no_robustness", action="store_true", default=False)
    ap.add_argument("--no_obs", action="store_true", default=False)

    args = ap.parse_args()

    ph12_vals = np.round(np.linspace(args.ph12_lo, args.ph12_hi, args.ph12_n), 6)
    ph23_vals = np.round(np.linspace(args.ph23_lo, args.ph23_hi, args.ph23_n), 6)

    if args.mode == "pmns_scan":
        csv_path = args.csv if args.csv is not None else "pmns_2d_scan.csv"
        run_pmns_2d_scan(
            eps_vals=(args.eps,),
            phase12_vals_pi=ph12_vals,
            phase23_vals_pi=ph23_vals,
            g_diag=args.g_diag,
            g_off=args.g_off,
            beta=args.beta,
            rel_cut=args.rel_cut,
            tol_rel_blocks=args.tol_rel_blocks,
            seesaw_M=args.seesaw_M,
            require_3d=args.require_3d,
            csv_path=csv_path,
            verbose=True
        )

    else:
        run_phase_map_over_phases(
            phase12_vals_pi=ph12_vals,
            phase23_vals_pi=ph23_vals,
            eps=args.eps,
            g_diag=args.g_diag,
            g_off=args.g_off,
            beta=args.beta,
            rel_cut=args.rel_cut,
            tol_rel_blocks=args.tol_rel_blocks,
            seesaw_M=args.seesaw_M,
            out_prefix=args.out_prefix,
            include_obs=(not args.no_obs),
            do_robustness=(not args.no_robustness)
        )


if __name__ == "__main__":
    main()


