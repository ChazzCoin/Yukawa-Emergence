#!/usr/bin/env python3
"""
DROP-IN UPGRADE (Majorana-symmetric seesaw + Takagi, delta_CP + Δm^2 extraction, scan scoring)

What changed vs your version:
  1) Seesaw is now the *symmetric* Majorana form:
        [ 0   mD ]
        [ mD^T MR ]
     and we diagonalize it with a Takagi factorization (not Hermitian eigh).

  2) We extract:
        - delta_CP (in radians + degrees) via rephasing-invariant quartets
        - Jarlskog J
        - dm21, dm31, dm32 from light masses

  3) Canonicalization returns the permutation so masses stay aligned with PMNS columns.

This is a complete file you can replace yours with.
"""

import json
import numpy as np
from numpy.linalg import eigh
from scipy.linalg import expm, svd
from dataclasses import dataclass, asdict
from typing import Optional, Dict, Any, List, Tuple
import csv


# ============================================================
#  RESULTS STRUCT
# ============================================================

@dataclass
class ScanResult:
    ok_3d: bool
    blocks: Any
    g_diag: float
    g_off: float
    eps: float
    phase_pi: float  # kept for compatibility; you now have phase12/phase23
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

    delta_cp_rad: Optional[float] = None
    delta_cp_deg: Optional[float] = None
    j_cp: Optional[float] = None

    score: Optional[float] = None
    note: str = ""


# ============================================================
#  NEW: Takagi + invariants
# ============================================================

def takagi(M: np.ndarray, atol: float = 1e-12):
    """
    Takagi factorization for complex symmetric M:
        U^T M U = diag(m_i),  m_i >= 0, U unitary
    Returns: (m, U) with m >= 0
    """
    M = np.asarray(M, dtype=complex)
    if not np.allclose(M, M.T, atol=atol):
        # helpful diagnostic if it’s only *nearly* symmetric
        max_asym = np.max(np.abs(M - M.T))
        raise ValueError(f"Takagi requires M == M.T (max |M-M^T| = {max_asym:.3e}).")

    U, s, Vh = svd(M)

    # phase-fix columns of U so diag(U^T M U) becomes real nonnegative
    T = U.T @ M @ U
    ph = np.angle(np.diag(T))
    U = U @ np.diag(np.exp(-0.5j * ph))

    # IMPORTANT: make it writable
    D = np.array(np.real(np.diag(U.T @ M @ U)), copy=True)

    # clamp tiny negatives from numerical noise
    D[D < 0] = 0.0

    return D, U


def dm2_from_masses(m):
    m = np.asarray(m, dtype=float)
    m2 = m**2
    dm21 = float(m2[1] - m2[0])
    dm31 = float(m2[2] - m2[0])
    dm32 = float(m2[2] - m2[1])
    return dm21, dm31, dm32


def delta_cp_from_U(U):
    """
    Extract delta_CP using rephasing-invariant quartet + PDG relations.
    Returns (delta_rad, J). delta in [-pi, pi].

    NOTE: Majorana phases are not extracted here (they require different invariants).
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

    # rephasing-invariant quartet
    P = U[0, 0] * U[1, 1] * np.conj(U[0, 1]) * np.conj(U[1, 0])
    J = float(np.imag(P))
    R = float(np.real(P))

    denom = (c12*s12*c23*s23*(c13**2)*s13)
    if abs(denom) < 1e-15:
        return None, J

    sin_delta = float(np.clip(J / denom, -1.0, 1.0))

    # Solve cosδ from Re(P) = A + B cosδ  (PDG parameterization identity)
    cos2_12 = c12*c12 - s12*s12
    A = -(c12*c12*s12*s12*(c13**2)) * (c23*c23 - s23*s23*s13*s13)
    B = -(c12*s12*(c13**2)) * (cos2_12 * c23*s23*s13)

    if abs(B) < 1e-15:
        # fallback: principal arcsin
        return float(np.arcsin(sin_delta)), J

    cos_delta = float(np.clip((R - A) / B, -1.0, 1.0))
    delta = float(np.arctan2(sin_delta, cos_delta))
    return delta, J


# ============================================================
#  GEOMETRY WEIGHTS (yours)
# ============================================================

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
        (1,2): 0.88 * np.exp(1j * phase23),
        (2,0): 0.80,
        (2,1): 0.90 * np.exp(-1j * phase23),
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
#  CONFIG (yours)
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
#  PIPELINE (upgraded seesaw + canonicalization perm)
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

    @staticmethod
    def _rotation_matrix_13(theta_rad):
        c = np.cos(theta_rad)
        s = np.sin(theta_rad)
        return np.array([
            [ c, 0.0,  s],
            [0.0, 1.0, 0.0],
            [-s, 0.0,  c],
        ], dtype=complex)

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

        theta = float(np.clip(theta, -0.10, 0.10))  # ~ ±5.7°
        Ue = self._rotation_matrix_13(theta)
        return Ue, theta

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

    def build_triads(self):
        cfg = self.cfg
        index = {g: i for i, g in enumerate(cfg.group_elements)}
        triads = []
        for s in cfg.triad_shifts:
            triads.append([index[self.add_g(h, s)] for h in cfg.subgroup_H])
        return triads

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
    def effective_yukawa(K_like, S_left, S_right=None):
        if S_right is None:
            S_right = S_left
        return S_left @ K_like @ S_right.conj().T

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

    # ---------- upgraded: symmetric majorana seesaw for blocks ----------
    @staticmethod
    def _seesaw_light_eigs_for_block_majorana(Y, block, M):
        Yb = Y[np.ix_(block, block)]
        n = len(block)

        MR = M * np.eye(n, dtype=complex)
        big = np.block([
            [np.zeros((n, n), dtype=complex), Yb],
            [Yb.T,                           MR]
        ])

        m, _U = takagi(big)
        m = np.sort(m)   # light first
        return m[:n]

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
            S_L = self.build_S(triads, shear=self.compression_shear)
            S_R = self.build_S(triads, shear=-self.compression_shear)
            Y = self.effective_yukawa(K_proj, S_L, S_R)

            return Y
        finally:
            cfg.geometry_weights = old_geom

    def apply_conditional_seesaw(self, Y, blocks):
        dirac = []
        majorana = []

        # 1D blocks treated "Dirac-like" using the (Hermitian) Y eigenvalue slice
        y_eigs = eigh(Y)[0]
        for block in blocks:
            if len(block) == 1:
                i = block[0]
                dirac.append((block, np.array([np.real(y_eigs[i])])))
            else:
                light = self._seesaw_light_eigs_for_block_majorana(Y, block, self.seesaw_M)
                majorana.append((block, light))

        return dirac, majorana

    def _canonicalize_pmns(self, U, return_perm=False):
        """
        Stable ordering + column phase-fix.
          - put "3" column as smallest |U_ei|
          - order remaining two by |U_ei| descending
          - phase-fix each column so U[0,k] real and >= 0
        Returns U, and optionally the permutation used.
        """
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

    # ---------- upgraded: full Takagi Majorana seesaw ----------
    def apply_full_seesaw_3x3(self, Y, include_charged_lepton=True, geometry_weights_charged=None):
        assert Y.shape == (3, 3)

        MR = self.seesaw_M * np.eye(3, dtype=complex)
        big = np.block([
            [np.zeros((3,3), dtype=complex), Y],
            [Y.T,                            MR]
        ])

        m6, U6 = takagi(big)  # masses >= 0, U6 unitary
        order = np.argsort(m6)  # light first
        m6 = m6[order]
        U6 = U6[:, order]

        light_masses = m6[:3].copy()
        Unu = U6[:3, :3].copy()

        for k in range(3):
            nrm = np.linalg.norm(Unu[:, k])
            if nrm > 0:
                Unu[:, k] /= nrm

        Unu, perm = self._canonicalize_pmns(Unu, return_perm=True)
        light_masses = light_masses[perm]

        diagnostics = {}

        if not include_charged_lepton:
            return light_masses, Unu, diagnostics

        if geometry_weights_charged is None:
            geom_e = {k: np.abs(v) for k, v in self.cfg.geometry_weights.items()}
            diagnostics["Y_e_geometry"] = "abs(phased_geometry)"
        else:
            geom_e = geometry_weights_charged
            diagnostics["Y_e_geometry"] = "override"

        Y_e = self._build_Y_with_geometry(geometry_weights_override=geom_e)

        Ue, theta_e = self.extract_charged_lepton_rotation(Y_e)
        diagnostics["theta_e_rad"] = float(theta_e)
        diagnostics["theta_e_deg"] = float(np.degrees(theta_e))
        diagnostics["Ue_from"] = "small_13_rotation_from_Ye"

        U_pmns = Ue.conj().T @ Unu

        # canonicalize again and re-order masses to stay aligned
        U_pmns, perm2 = self._canonicalize_pmns(U_pmns, return_perm=True)
        light_masses = light_masses[perm2]

        return light_masses, U_pmns, diagnostics

    def run(self):
        cfg = self.cfg

        self.triads = self.build_triads()
        self.K = self.build_kernel()
        self.K_flow = self.misalignment_flow(self.K)
        self.P_C360, self.kept_indices = self.emergent_C360_projector(self.K)
        self.K_proj = self.P_C360 @ self.K_flow @ self.P_C360

        S_L = self.build_S(self.triads, shear=self.compression_shear)
        S_R = self.build_S(self.triads, shear=-self.compression_shear)  # simplest asymmetry knob

        self.S = S_L
        self.Y = self.effective_yukawa(self.K_proj, S_L, S_R)

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
#  SCORING
# ============================================================

def pmns_score(theta12, theta13, theta23,
               target=(33.4, 8.6, 45.0),
               weights=(1.0, 2.0, 1.0)):
    t12, t13, t23 = target
    w12, w13, w23 = weights
    return (w12*(theta12 - t12)**2 +
            w13*(theta13 - t13)**2 +
            w23*(theta23 - t23)**2)


def phys_score(th12, th13, th23, dm21, dm31, delta_deg,
               ang_target=(33.4, 8.6, 45.0),
               r_target=0.0296,
               delta_target_deg=-90.0,
               w_ang=(1.0, 6.0, 1.0),
               w_r=2e4,
               w_delta=0.02):
    a = pmns_score(th12, th13, th23, target=ang_target, weights=w_ang)

    if dm31 == 0 or np.isnan(dm21) or np.isnan(dm31):
        b = 10.0
    else:
        r = float(dm21 / abs(dm31))
        b = w_r * (r - r_target)**2

    if delta_deg is None:
        c = 10.0
    else:
        d = ((delta_deg - delta_target_deg + 180) % 360) - 180
        c = w_delta*(d**2)

    return float(a + b + c)



# ============================================================
#  2D SCAN (upgraded outputs + optional physical score)
# ============================================================

def run_pmns_2d_scan(
    *,
    phase12_vals_pi=None,
    phase23_vals_pi=None,
    eps_vals=(0.013,0.010),
    g_diag=0.894,
    g_off=1.0,
    beta=1.5,
    rel_cut=0.15,
    tol_rel_blocks=0.03,
    seesaw_M=1e6,
    require_3d=True,
    top_k=15,
    csv_path="pmns_2d_scan.csv",
    verbose=True,
    optimize_physical=True
):
    if phase12_vals_pi is None:
        phase12_vals_pi = np.round(np.linspace(0.070, 0.085, 7), 5)
    if phase23_vals_pi is None:
        phase23_vals_pi = np.round(np.linspace(0.010, 0.030, 9), 5)

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
                        "eps": float(eps), "phase12_pi": float(ph12), "phase23_pi": float(ph23),
                        "ok_3d": False, "blocks": str(blocks),
                        "theta12": None, "theta13": None, "theta23": None,
                        "Ue3": None,
                        "dm21": None, "dm31": None, "dm32": None,
                        "delta_cp_deg": None, "J_cp": None,
                        "score": None
                    })
                    continue

                m_light, U_pmns, diag = pipe.apply_full_seesaw_3x3(out["Y"], include_charged_lepton=True)

                # fixed: print the angle, not the dict
                if verbose and (n % 25 == 0 or n == total):
                    print(f"θe (charged lepton) = {diag['theta_e_deg']:.3f}°")

                absU = np.abs(U_pmns)
                th12, th13, th23 = angles_from_absU(absU)
                Ue3 = float(absU[0, 2])

                dm21, dm31, dm32 = dm2_from_masses(m_light)
                delta_rad, J = delta_cp_from_U(U_pmns)
                delta_deg = None if delta_rad is None else float(np.degrees(delta_rad))
                if verbose and (n % 25 == 0 or n == total):
                    print(f"J_CP={0.0 if J is None else J:.3e}  δ={('NA' if delta_deg is None else f'{delta_deg:.1f}°')}")

                if optimize_physical:
                    score = phys_score(th12, th13, th23, dm21, dm31, delta_deg)
                else:
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
                    "m1": float(m_light[0]),
                    "m2": float(m_light[1]),
                    "m3": float(m_light[2]),
                    "dm21": float(dm21),
                    "dm31": float(dm31),
                    "dm32": float(dm32),
                    "delta_cp_deg": (None if delta_deg is None else float(delta_deg)),
                    "J_cp": (None if J is None else float(J)),
                    "score": float(score),
                })

                if verbose and (n % 25 == 0 or n == total):
                    print(f"[{n}/{total}] eps={eps:.4f} ph12={ph12:.5f}π ph23={ph23:.5f}π "
                          f"θ12={th12:.2f} θ13={th13:.2f} θ23={th23:.2f} "
                          f"dm21={dm21:.3e} dm31={dm31:.3e} δ={('NA' if delta_deg is None else f'{delta_deg:.1f}°')} "
                          f"score={score:.3f}")

    valid = [r for r in rows if r["score"] is not None]
    valid.sort(key=lambda r: r["score"])

    print("\n===== TOP 2D CANDIDATES =====")
    for i, r in enumerate(valid[:top_k], start=1):
        dstr = "NA" if r["delta_cp_deg"] is None else f"{r['delta_cp_deg']:.1f}°"
        print(f"{i:2d}) score={r['score']:.3f}  eps={r['eps']:.4f}  "
              f"ph12={r['phase12_pi']:.5f}π  ph23={r['phase23_pi']:.5f}π  "
              f"| θ12={r['theta12']:.2f} θ13={r['theta13']:.2f} θ23={r['theta23']:.2f}  "
              f"Ue3={r['Ue3']:.4f}  dm21={r['dm21']:.3e} dm31={r['dm31']:.3e}  δ={dstr}")

    if csv_path and rows:
        keys = list(rows[0].keys())
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
        verbose=True,
        optimize_physical=True
    )
