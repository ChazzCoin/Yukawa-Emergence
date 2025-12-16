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


def set_geometry_weights_pmns_safe(
    cfg,
    g_diag=0.90,   # back to known 3D-safe value
    g_off=1.00,
    eps=0.02,      # smaller gen-1 boost
):
    """
    PMNS-safe refinement that preserves the 3D harmonic block.
    """

    base = {
        (0,0): 1.00, (0,1): 0.92, (0,2): 0.85,
        (1,0): 0.95, (1,1): 1.10, (1,2): 0.88,
        (2,0): 0.80, (2,1): 0.90, (2,2): 1.05,
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
def set_geometry_weights_pmns_with_phase(
    cfg,
    g_diag=0.90,
    g_off=1.00,
    eps=0.02,
    phase=0.15*np.pi
):
    """
    PMNS refinement with a single internal phase.
    Preserves the 3D harmonic block while increasing θ12 and θ13.
    """

    base = {
        (0,0): 1.00,
        (0,1): 0.92 * np.exp(1j * phase),   # phased coupling
        (0,2): 0.85,
        (1,0): 0.95 * np.exp(-1j * phase),  # conjugate phase
        (1,1): 1.10,
        (1,2): 0.88,
        (2,0): 0.80,
        (2,1): 0.90,
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

    def extract_charged_lepton_rotation(self, Y, blocks):
        """
        Extract a small left-handed rotation Ue from 1D Dirac blocks.

        Strategy:
        - Use the eigenvector of the dominant 1D block
        - Map its misalignment relative to flavor basis into a 1–3 rotation
        """
        # find first 1D block
        one_d_blocks = [b for b in blocks if len(b) == 1]
        if not one_d_blocks:
            return np.eye(3, dtype=complex), 0.0

        # diagonalize Y to get eigenvectors
        _, U = eigh(Y)

        # choose the eigenvector corresponding to the 1D block
        i = one_d_blocks[0][0]
        v = U[:, i]

        # project misalignment onto 1–3 plane
        # (this is the physically relevant suppression direction)
        theta_e = np.arctan2(
            np.real(v[2]),
            np.real(v[0]) if abs(v[0]) > 1e-12 else 1e-12
        )

        # clamp to small angles (safety)
        theta_e = float(np.clip(theta_e, -0.05, 0.05))  # ~ ±3°

        Ue = self._rotation_matrix_13(theta_e)
        return Ue, theta_e

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
        to represent a different sector (charged leptons vs neutrinos),
        while keeping kernel_characters / damping / triads / S the same.
        """
        cfg = self.cfg

        # temporarily override geometry
        old_geom = cfg.geometry_weights
        if geometry_weights_override is not None:
            cfg.geometry_weights = geometry_weights_override

        try:
            K = self.build_kernel()
            K_flow = self.misalignment_flow(K)
            P, _kept = self.emergent_C360_projector(K)
            K_proj = P @ K_flow @ P

            triads = self.build_triads()
            S = self.build_S(triads)

            Y = self.effective_yukawa(K_proj, S)
            return Y
        finally:
            cfg.geometry_weights = old_geom

    def build_lepton_sector_PMNS(self, geometry_weights_charged=None):
        """
        Build:
          - neutrino operator Y_nu from current cfg.geometry_weights (phased)
          - charged lepton operator Y_e from a (usually unphased) geometry override

        Returns:
          Ue, Unu, Upmns
        """
        # neutrino operator (uses current geometry)
        Y_nu = self.Y  # assumes self.run() already built this

        # charged lepton operator
        if geometry_weights_charged is None:
            # default: "same magnitudes, phases removed" (safe minimal choice)
            geom_e = {}
            for k, v in self.cfg.geometry_weights.items():
                geom_e[k] = np.abs(v)
        else:
            geom_e = geometry_weights_charged

        Y_e = self._build_Y_with_geometry(geometry_weights_override=geom_e)

        # diagonalize charged leptons: Y_e = Ue diag Ue†
        _, Ue = eigh(Y_e)

        # neutrinos: Unu comes from seesaw
        m_light, Unu, _diag = self.apply_full_seesaw_3x3(self.Y, include_charged_lepton=False)

        # physical PMNS
        U_pmns = Ue.conj().T @ Unu
        return Y_e, Y_nu, Ue, Unu, U_pmns

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

    def apply_full_seesaw_3x3(self, Y, include_charged_lepton=True, geometry_weights_charged=None):
        """
        Full seesaw on Y -> Unu, and optionally apply Ue† from a charged-lepton Y_e.
        """
        zero = np.zeros_like(Y)
        MR = self.seesaw_M * np.eye(3)

        big = np.block([
            [zero, Y],
            [Y.conj().T, MR]
        ])

        eigvals, eigvecs = eigh(big)
        eigvals = np.real(eigvals)

        order = np.argsort(np.abs(eigvals))
        eigvals = eigvals[order]
        eigvecs = eigvecs[:, order]

        light_masses = np.abs(eigvals[:3])
        Unu = eigvecs[:3, :3]

        # normalize columns
        for k in range(3):
            nrm = np.linalg.norm(Unu[:, k])
            if nrm > 0:
                Unu[:, k] /= nrm

        diagnostics = {}

        if not include_charged_lepton:
            return light_masses, Unu, diagnostics

        # build Ue from charged-lepton operator
        if geometry_weights_charged is None:
            geom_e = {k: np.abs(v) for k, v in self.cfg.geometry_weights.items()}
        else:
            geom_e = geometry_weights_charged

        Y_e = self._build_Y_with_geometry(geometry_weights_override=geom_e)
        _, Ue = eigh(Y_e)

        U_pmns = Ue.conj().T @ Unu
        diagnostics["Ue_from"] = "charged_geometry_override"
        return light_masses, U_pmns, diagnostics

    # ---------- run the whole pipeline ----------

    def run(self):
        cfg = self.cfg

        self.triads = self.build_triads()
        self.K = self.build_kernel()
        self.S = self.build_S(self.triads)

        self.K_flow = self.misalignment_flow(self.K)
        self.P_C360, self.kept_indices = self.emergent_C360_projector(self.K)
        self.K_proj = self.P_C360 @ self.K_flow @ self.P_C360

        self.Y = self.effective_yukawa(self.K_proj, self.S)

        # Y spectrum
        y_vals, U = eigh(self.Y)
        masses = np.abs(y_vals) * cfg.higgs_vev

        # block structure (this is the “emergence diagnosis”)
        blocks, y_sorted = self.harmonic_blocks_on_Y(self.Y)

        # seesaw only where needed
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

def run_lepton_sector_refined():
    cfg = AlignmentV33Config()

    set_geometry_weights_pmns_safe(
        cfg,
        g_diag=0.88,
        g_off=1.00,
        eps=0.02
    )

    pipe = AlignmentPipeline(
        cfg,
        beta=1.5,
        rel_cut=0.15,
        tol_rel_blocks=0.03,
        seesaw_M=1e6
    )

    out = pipe.run()
    print("\nLepton-sector harmonic blocks:", out["blocks"])

    if out["blocks"] != [[0, 1, 2]]:
        print("⚠️  Warning: not in full 3D harmonic regime")
        return

    light_masses, U_pmns = pipe.apply_full_seesaw_3x3(out["Y"])

    print("\nLight neutrino masses:", light_masses)
    print("\n|U_PMNS|:\n", np.abs(U_pmns))

    th12, th13, th23 = angles_from_absU(np.abs(U_pmns))
    print(f"\nPMNS angles: θ12={th12:.1f}°, θ13={th13:.1f}°, θ23={th23:.1f}°")

    return light_masses, U_pmns

def run_phased_scan():
    cfg = AlignmentV33Config()
    for ph in [0.08, 0.075, 0.07, 0.065, 0.06]:
        print(f"\n-> Running phased scan on [ {ph} ] <-")
        run_lepton_sector_phased(cfg=cfg, g_diag=0.8945, g_off=1.00, eps=0.012, phase=ph)
        print(f"\n XX Finished phased scan on [ {ph} ] XX \n")

def run_lepton_sector_phased(cfg, g_diag, g_off, eps, phase):
    if not cfg:
        print("\nNo configuration provided, creating a new one...\n")
        cfg = AlignmentV33Config()

    set_geometry_weights_pmns_with_phase(
        cfg,
        g_diag=g_diag,
        g_off=g_off,
        eps=eps,
        phase=phase*np.pi
    )

    pipe = AlignmentPipeline(
        cfg,
        beta=1.5,
        rel_cut=0.15,
        tol_rel_blocks=0.03,
        seesaw_M=1e6
    )

    out = pipe.run()
    print("\nLepton-sector harmonic blocks:", out["blocks"])

    if out["blocks"] != [[0,1,2]]:
        print("⚠️  Lost 3D harmonic block — reduce phase or eps")
        return

    light_masses, U_pmns = pipe.apply_full_seesaw_3x3(out["Y"])

    print("\nLight neutrino masses:", light_masses)
    print("\n|U_PMNS|:\n", np.abs(U_pmns))

    th12, th13, th23 = angles_from_absU(np.abs(U_pmns))
    print(f"\nPMNS angles: θ12={th12:.1f}°, θ13={th13:.1f}°, θ23={th23:.1f}°")

    return light_masses, U_pmns

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


def evaluate_point(
    *,
    g_diag: float,
    g_off: float,
    eps: float,
    phase_pi: float,
    beta: float = 1.5,
    rel_cut: float = 0.15,
    tol_rel_blocks: float = 0.03,
    seesaw_M: float = 1e6,
    require_3d: bool = True,
    cfg_factory= None,
    r1=None, r2=None, r3=None,
    theta13_deg=None, phase12=None, phase13=None
) -> ScanResult:
    """
    Runs one lepton-sector evaluation and returns a ScanResult.
    Creates a fresh cfg each time by default to avoid state bleed.
    """
    if cfg_factory is None:
        cfg_factory = AlignmentV33Config

    cfg = cfg_factory()

    # Apply geometry (phase_pi is multiple of pi)
    set_geometry_weights_pmns_two_phase(
        cfg,
        g_diag=g_diag,
        g_off=g_off,
        eps=eps,
        phase12=phase12,
        phase23=phase13
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

    res = ScanResult(
        ok_3d=ok_3d,
        blocks=blocks,
        g_diag=g_diag,
        g_off=g_off,
        eps=eps,
        phase_pi=phase_pi,
        beta=beta,
        rel_cut=rel_cut,
        tol_rel_blocks=tol_rel_blocks,
        seesaw_M=seesaw_M,
    )

    if require_3d and not ok_3d:
        res.note = "lost_3d_block"
        return res

    # Full 3x3 seesaw
    theta13_deg = 0.0 if theta13_deg is None else float(theta13_deg)

    m_light, U_pmns, diag = pipe.apply_full_seesaw_3x3(out["Y"], include_charged_lepton=True)

    absU = np.abs(U_pmns)

    th12, th13, th23 = angles_from_absU(absU)

    res.theta12, res.theta13, res.theta23 = float(th12), float(th13), float(th23)
    res.ue3, res.um3, res.ut3 = float(absU[0, 2]), float(absU[1, 2]), float(absU[2, 2])
    res.m_light = [float(x) for x in m_light]

    # Simple mass-splitting diagnostics (dimensionless, squared)
    m = np.array(m_light, dtype=float)
    # Ensure ascending
    m = np.sort(m)
    res.dm21 = float(m[1]**2 - m[0]**2)
    res.dm31 = float(m[2]**2 - m[0]**2)

    return res

def run_full_phased_scan(
    *,
    phases_pi=(0.08, 0.075, 0.07, 0.065, 0.06),
    eps_vals=(0.012,),
    g_diag_vals=(0.8945,),
    g_off_vals=(1.00,),
    beta=1.5,
    rel_cut=0.15,
    tol_rel_blocks=0.03,
    seesaw_M=1e6,
    require_3d=True,
    # Optional angle windows to filter:
    angle_window=None,  # e.g. {"t12":(25,40), "t13":(5,12), "t23":(35,60)}
    # Scoring:
    target_angles=(33.4, 8.6, 45.0),
    weights=(1.0, 2.0, 1.0),
    top_k=10,
    csv_path="pmns_scan_results.csv",
    verbose=True,
):
    results: List[ScanResult] = []

    total = len(phases_pi)*len(eps_vals)*len(g_diag_vals)*len(g_off_vals)
    done = 0

    for g_off in g_off_vals:
        for g_diag in g_diag_vals:
            for eps in eps_vals:
                for ph in phases_pi:
                    done += 1
                    if verbose:
                        print(f"\n[{done}/{total}] g_diag={g_diag} g_off={g_off} eps={eps} phase={ph}π")

                    res = evaluate_point(
                        g_diag=g_diag,
                        g_off=g_off,
                        eps=eps,
                        phase_pi=ph,
                        beta=beta,
                        rel_cut=rel_cut,
                        tol_rel_blocks=tol_rel_blocks,
                        seesaw_M=seesaw_M,
                        require_3d=require_3d
                    )

                    if not res.ok_3d and require_3d:
                        if verbose:
                            print(f"  blocks={res.blocks}  -> rejected (lost 3D)")
                        results.append(res)
                        continue

                    # Apply angle window filter if requested
                    if angle_window is not None:
                        t12_lo, t12_hi = angle_window.get("t12", (-np.inf, np.inf))
                        t13_lo, t13_hi = angle_window.get("t13", (-np.inf, np.inf))
                        t23_lo, t23_hi = angle_window.get("t23", (-np.inf, np.inf))

                        if not (t12_lo <= res.theta12 <= t12_hi and
                                t13_lo <= res.theta13 <= t13_hi and
                                t23_lo <= res.theta23 <= t23_hi):
                            res.note = "outside_angle_window"
                            if verbose:
                                print(f"  θ12={res.theta12:.1f} θ13={res.theta13:.1f} θ23={res.theta23:.1f} -> window reject")
                            results.append(res)
                            continue

                    # Score
                    res.score = float(pmns_score(res.theta12, res.theta13, res.theta23,
                                                 target=target_angles, weights=weights))

                    if verbose:
                        print(f"  blocks={res.blocks}")
                        print(f"  θ12={res.theta12:.2f}°, θ13={res.theta13:.2f}°, θ23={res.theta23:.2f}°  score={res.score:.2f}")
                        print(f"  |Ue3|={res.ue3:.4f}  dm21={res.dm21:.3e}  dm31={res.dm31:.3e}")

                    results.append(res)

    # Write CSV
    if csv_path:
        fieldnames = list(asdict(results[0]).keys()) if results else []
        with open(csv_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            for r in results:
                d = asdict(r)
                # flatten list for CSV readability
                if d["m_light"] is not None:
                    d["m_light"] = ";".join(f"{x:.12e}" for x in d["m_light"])
                w.writerow(d)

        if verbose:
            print(f"\nWrote CSV: {csv_path}")

    # Rank top results that are valid + scored
    valid = [r for r in results if r.score is not None]
    valid.sort(key=lambda r: r.score)

    if verbose:
        print("\n===== TOP CANDIDATES =====")
        for i, r in enumerate(valid[:top_k], start=1):
            print(f"{i:2d}) score={r.score:.2f}  phase={r.phase_pi}π eps={r.eps} g_diag={r.g_diag} g_off={r.g_off} "
                  f"| θ12={r.theta12:.2f} θ13={r.theta13:.2f} θ23={r.theta23:.2f}  blocks={r.blocks}")

    return results

def run_pmns_micro_scan():
    # Centered around your best solution
    phases_pi = np.round(np.linspace(0.06, 0.085, 11), 5)   # 0.060π ... 0.085π
    eps_vals  = [0.010, 0.011, 0.012, 0.013]
    g_diag_vals = [0.892, 0.8935, 0.8945, 0.8955]
    g_off_vals  = [1.00]

    angle_window = {
        "t12": (28, 38),
        "t13": (6, 12),
        "t23": (38, 55),
    }

    return run_full_phased_scan(
        phases_pi=phases_pi,
        eps_vals=eps_vals,
        g_diag_vals=g_diag_vals,
        g_off_vals=g_off_vals,
        beta=1.5,
        rel_cut=0.15,
        tol_rel_blocks=0.03,
        seesaw_M=1e6,
        require_3d=True,
        angle_window=angle_window,
        target_angles=(33.4, 8.6, 45.0),
        weights=(1.0, 3.0, 1.0),   # emphasize θ13
        top_k=12,
        csv_path="pmns_micro_scan.csv",
        verbose=True
    )
def run_MR_ratio_scan():
    base = dict(
        g_diag=0.894,
        g_off=1.00,
        eps=0.012,
        phase_pi=0.075,
        beta=1.5,
        rel_cut=0.15,
        tol_rel_blocks=0.03,
        seesaw_M=1e6,
        require_3d=True
    )
    print("- Base Dict - ")
    print(json.dumps(base, indent=2))
    ratios = [0.5, 0.7, 1.0, 1.5, 2.0]
    thetas = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0]  # degrees

    best = None
    best_pack = None

    for r1 in ratios:
        for r2 in ratios:
            for r3 in ratios:
                for th in thetas:
                    res = evaluate_point(**base, r1=r1, r2=r2, r3=r3, theta13_deg=th)

                    if not res.ok_3d:
                        continue

                    # emphasize theta13
                    res.score = float(pmns_score(res.theta12, res.theta13, res.theta23,
                                                 target=(33.4, 8.6, 45.0),
                                                 weights=(1.0, 6.0, 1.0)))

                    if best is None or res.score < best:
                        best = res.score
                        best_pack = (r1, r2, r3, th, res)

    print("\nBEST ROTATED-MR POINT:")
    if best_pack:
        r1, r2, r3, th, res = best_pack
        print(f"score={res.score:.2f}  MR ratios=(r1={r1}, r2={r2}, r3={r3})  theta13_R={th}°")
        print(f"θ12={res.theta12:.2f}  θ13={res.theta13:.2f}  θ23={res.theta23:.2f}")
        print(f"|Ue3|={res.ue3:.4f}  dm21={res.dm21:.3e}  dm31={res.dm31:.3e}")
        print(f"m_light={res.m_light}")
    else:
        print("No valid 3D-block points found.")

def run_pmns_scan(
    *,
    phase23_vals=(0.0, 0.005, 0.01, 0.015, 0.02),
    g_diag=0.894,
    g_off=1.0,
    eps=0.012,
    phase12_pi=0.075,   # interpreted as multiple of pi
    beta=1.5,
    rel_cut=0.15,
    tol_rel_blocks=0.03,
    seesaw_M=1e6,
    require_3d=True,
):
    results = []

    print("\nPMNS two-phase scan (varying phase23)")
    print("---------------------------------------------------------------")
    print("ph23π    blocks        θ12     θ13     θ23     |Ue3|     score")
    print("---------------------------------------------------------------")

    for ph23 in phase23_vals:
        # IMPORTANT: fresh cfg each loop to avoid state bleed
        cfg = AlignmentV33Config()

        set_geometry_weights_pmns_two_phase(
            cfg,
            g_diag=g_diag,
            g_off=g_off,
            eps=eps,
            phase12=phase12_pi * np.pi,
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

        row = {
            "phase23_pi": float(ph23),
            "blocks": blocks,
            "theta12": None,
            "theta13": None,
            "theta23": None,
            "Ue3": None,
            "score": None,
            "m_light": None,
        }

        ok_3d = (blocks == [[0, 1, 2]])

        if require_3d and not ok_3d:
            print(f"{ph23:0.4f}  {str(blocks):<12}   --      --      --      --       --  (lost 3D)")
            results.append(row)
            continue

        # Full 3x3 seesaw on the effective Y
        m_light, U_pmns, diag = pipe.apply_full_seesaw_3x3(out["Y"], include_charged_lepton=True)

        absU = np.abs(U_pmns)

        th12, th13, th23 = angles_from_absU(absU)
        Ue3 = float(absU[0, 2])

        # same score convention as your scans
        score = float(pmns_score(th12, th13, th23, target=(33.4, 8.6, 45.0), weights=(1.0, 4.0, 1.0)))

        row.update({
            "theta12": float(th12),
            "theta13": float(th13),
            "theta23": float(th23),
            "Ue3": Ue3,
            "score": score,
            "m_light": [float(x) for x in m_light],
        })
        results.append(row)

        print(f"{ph23:0.4f}  {str(blocks):<12}  {th12:6.2f}  {th13:6.2f}  {th23:6.2f}  {Ue3:7.4f}  {score:7.2f}")

    print("---------------------------------------------------------------")

    # Print best (among valid)
    valid = [r for r in results if r["score"] is not None]
    if valid:
        best = min(valid, key=lambda r: r["score"])
        print("\nBEST phase23 point:")
        print(f"phase23 = {best['phase23_pi']}π | blocks={best['blocks']}")
        print(f"θ12={best['theta12']:.2f}  θ13={best['theta13']:.2f}  θ23={best['theta23']:.2f}  |Ue3|={best['Ue3']:.4f}")
        print(f"score={best['score']:.2f}")
        print(f"m_light={best['m_light']}")

    return results

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

                print(f"θe (charged lepton) = {diag}°")

                absU = np.abs(U_pmns)
                th12, th13, th23 = angles_from_absU(absU)
                Ue3 = float(absU[0, 2])


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
        eps_vals=(0.010, 0.011, 0.012, 0.013),
        phase12_vals_pi=np.round(np.linspace(0.072, 0.088, 9), 5),
        phase23_vals_pi=np.round(np.linspace(0.012, 0.032, 11), 5),
        verbose=True
    )