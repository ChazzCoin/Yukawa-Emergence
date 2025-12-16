#!/usr/bin/env python3
"""
Alignment Spectral Triple v4.0 — Engine Outline (Python)

Goal:
  Build a production-ready v4.0 engine that:
    (1) defines the internal SM basis (32 states per generation),
    (2) defines lifted projectors, gamma_SM, and real structure J,
    (3) builds channel maps V_u,V_d,V_e,V_nu,W_R,
    (4) assembles D_int(Y[K]) as a 32N×32N sparse block matrix,
    (5) assembles the even product Dirac D = D_geom⊗1 + γ_geom⊗D_int,
    (6) builds one-forms and inner fluctuations D_A = D + A + JAJ^{-1},
    (7) runs commutator gates (order-zero, first-order, fluctuation stability),
    (8) provides sector projectors + extraction utilities for Yukawas/PMNS,
    (9) provides serialization + reproducibility and scan harness hooks.

This is an OUTLINE: signatures + responsibilities + minimal stubs.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple, Callable

import numpy as np
from numpy.linalg import eigh, norm
from scipy.linalg import expm


# ============================================================
#  Core numerics helpers
# ============================================================

def op_norm(A: np.ndarray) -> float:
    """Operator norm (spectral norm): ||A||_2."""
    # robust: largest singular value
    s = np.linalg.svd(A, compute_uv=False)
    return float(s[0]) if s.size else 0.0

def hermitian_part(A: np.ndarray) -> np.ndarray:
    return 0.5 * (A + A.conj().T)

def comm(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    return A @ B - B @ A

def anti_comm(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    return A @ B + B @ A

def is_hermitian(A: np.ndarray, atol: float = 1e-12) -> bool:
    return bool(np.allclose(A, A.conj().T, atol=atol))

def is_projector(P: np.ndarray, atol: float = 1e-12) -> bool:
    return bool(np.allclose(P @ P, P, atol=atol) and is_hermitian(P, atol=atol))

def put_block(D: np.ndarray, i: int, j: int, B: np.ndarray) -> None:
    """Place an NxN block into SM slot (i,j) of a (32N)x(32N) matrix."""
    N = B.shape[0]
    r = slice(i * N, (i + 1) * N)
    c = slice(j * N, (j + 1) * N)
    D[r, c] += B


# ============================================================
#  v4.0 Configuration and parameter surfaces
# ============================================================

@dataclass(frozen=True)
class V40Params:
    """
    Parameters that vary scan-to-scan (textures, scales, cutoffs).
    Treat these as immutable per run for reproducibility.
    """
    N_flav: int

    # Higgs vev (for physical scaling, optional)
    higgs_vev_GeV: float = 174.0

    # Majorana scale (can be a matrix via builder hook)
    seesaw_M_GeV: float = 2.0e14

    # Gate tolerances
    eps_order0: float = 1e-10
    eps_first: float = 1e-10
    eps_sa: float = 1e-10
    eps_fluct: float = 1e-10

    # Optional scan knobs
    beta: float = 1.5
    rel_cut: float = 0.15
    tol_rel_blocks: float = 0.03


@dataclass(frozen=True)
class TexturePack:
    """
    Flavor textures (operators on H_flav), one per fermion sector.
    Each Y_* is NxN complex.
    """
    Yu: np.ndarray
    Yd: np.ndarray
    Ye: np.ndarray
    Ynu: np.ndarray
    MR: Optional[np.ndarray] = None  # NxN complex symmetric preferred


# ============================================================
#  SM basis (32 states) + projectors + grading
# ============================================================

class SMBasis32:
    """
    Owns the canonical 32-state ordering (16 particle + 16 conjugate),
    and provides projector matrices P_X in M_32(C) and grading gamma_SM.
    """

    def __init__(self) -> None:
        self.dim = 32
        self._index_map = self._build_index_map()  # names -> indices

        # cached projectors (32x32)
        self.P = self._build_projectors()

        # grading gamma_SM (32x32)
        self.gamma_SM = self._build_gamma_SM()

    # -------- basis + indices --------

    def _build_index_map(self) -> Dict[str, List[int]]:
        """
        Return dict mapping semantic labels -> index lists in the 32 ordering.
        Must match the fixed ordering used in the v4.0 paper.
        """
        # Particle: 0..15
        # (uL r,g,b, dL r,g,b, nuL, eL, uR r,g,b, dR r,g,b, eR, nuR)
        # Conjugates: +16
        uL = [0, 1, 2]
        dL = [3, 4, 5]
        nuL = [6]
        eL = [7]
        uR = [8, 9, 10]
        dR = [11, 12, 13]
        eR = [14]
        nuR = [15]

        # conjugates
        uL_c = [i + 16 for i in uL]
        dL_c = [i + 16 for i in dL]
        nuL_c = [i + 16 for i in nuL]
        eL_c = [i + 16 for i in eL]
        uR_c = [i + 16 for i in uR]
        dR_c = [i + 16 for i in dR]
        eR_c = [i + 16 for i in eR]
        nuR_c = [i + 16 for i in nuR]

        return {
            "uL": uL, "dL": dL, "nuL": nuL, "eL": eL,
            "uR": uR, "dR": dR, "eR": eR, "nuR": nuR,
            "uL_c": uL_c, "dL_c": dL_c, "nuL_c": nuL_c, "eL_c": eL_c,
            "uR_c": uR_c, "dR_c": dR_c, "eR_c": eR_c, "nuR_c": nuR_c,
            "part": list(range(0, 16)),
            "conj": list(range(16, 32)),
            "L_part": list(range(0, 8)),
            "R_part": list(range(8, 16)),
            "L_conj": list(range(24, 32)),
            "R_conj": list(range(16, 24)),
        }

    def indices(self, key: str) -> List[int]:
        return list(self._index_map[key])

    # -------- projectors --------

    def _proj_from_indices(self, idx: Sequence[int]) -> np.ndarray:
        P = np.zeros((self.dim, self.dim), dtype=complex)
        for i in idx:
            P[i, i] = 1.0
        return P

    def _build_projectors(self) -> Dict[str, np.ndarray]:
        """
        Build the basic sector projectors needed for v4.0 extraction and assembly.
        """
        P: Dict[str, np.ndarray] = {}
        for k, idx in self._index_map.items():
            P[k] = self._proj_from_indices(idx)
        return P

    # -------- grading gamma_SM --------

    def _build_gamma_SM(self) -> np.ndarray:
        """
        v4.0 grading on the SM finite space (32x32).

        Convention used earlier:
          +1 on: L particles (0..7) and L conjugates (24..31)
          -1 on: R particles (8..15) and R conjugates (16..23)
        """
        g = np.zeros((self.dim,), dtype=float)
        for i in self._index_map["L_part"] + self._index_map["L_conj"]:
            g[i] = +1.0
        for i in self._index_map["R_part"] + self._index_map["R_conj"]:
            g[i] = -1.0
        return np.diag(g).astype(complex)

    # -------- lifted operators --------

    def lift_projector(self, P32: np.ndarray, N: int) -> np.ndarray:
        """Lift P (32x32) to (32N)x(32N): P ⊗ I_N."""
        return np.kron(P32, np.eye(N, dtype=complex))

    def lift_gamma(self, N: int) -> np.ndarray:
        return self.lift_projector(self.gamma_SM, N)


# ============================================================
#  Real structures (antiunitary) on finite spaces
# ============================================================

@dataclass(frozen=True)
class AntiUnitary:
    """
    Antiunitary J implemented as J = U K where K is entrywise conjugation.
    We implement conjugation action on operators: J X J^{-1} = U \bar{X} U^\dagger.
    """
    U: np.ndarray  # unitary matrix

    def conj_op(self, X: np.ndarray) -> np.ndarray:
        """Return J X J^{-1}."""
        return self.U @ np.conj(X) @ self.U.conj().T

    def square_sign(self, atol: float = 1e-12) -> Optional[int]:
        """
        Return sign of J^2 if inferable in this basis (optional).
        For J=U K: J^2 = U \bar{U}.
        """
        J2 = self.U @ np.conj(self.U)
        # if J2 ≈ ±I, return that sign
        I = np.eye(J2.shape[0], dtype=complex)
        if np.allclose(J2, I, atol=atol):
            return +1
        if np.allclose(J2, -I, atol=atol):
            return -1
        return None


class RealStructureFactory:
    """
    Factory for building J_SM, J_flav, and lifted product J_int = J_SM ⊗ J_flav.
    """

    @staticmethod
    def build_J_flav(N: int, U_flav: Optional[np.ndarray] = None) -> AntiUnitary:
        """
        Default: J_flav = K (so U_flav = I).
        """
        if U_flav is None:
            U_flav = np.eye(N, dtype=complex)
        return AntiUnitary(U=U_flav)

    @staticmethod
    def build_J_SM(basis: SMBasis32) -> AntiUnitary:
        """
        Placeholder: in production you pin the exact U_SM implementing charge conjugation.
        For now, provide identity (i.e. J=K) as a scaffold.
        """
        U_SM = np.eye(basis.dim, dtype=complex)
        return AntiUnitary(U=U_SM)

    @staticmethod
    def kron(J1: AntiUnitary, J2: AntiUnitary) -> AntiUnitary:
        """Product antiunitary: (U1⊗U2) K."""
        return AntiUnitary(U=np.kron(J1.U, J2.U))


# ============================================================
#  Channel maps V_u,V_d,V_e,V_nu and W_R
# ============================================================

class ChannelMaps:
    """
    Builds partial isometries on H_SM (32x32):
      V_u,V_d,V_e,V_nu  and W_R (nu_R -> nu_R^c).

    The minimal concrete implementation uses the fixed index ordering
    and builds each map as a sum of matrix units E_ij (color-preserving).
    """

    def __init__(self, basis: SMBasis32) -> None:
        self.basis = basis
        self.Vu = self._build_Vu()
        self.Vd = self._build_Vd()
        self.Ve = self._build_Ve()
        self.Vnu = self._build_Vnu()
        self.WR = self._build_WR()

    def _E(self, i: int, j: int) -> np.ndarray:
        E = np.zeros((self.basis.dim, self.basis.dim), dtype=complex)
        E[i, j] = 1.0
        return E

    def _sum_E(self, pairs: Sequence[Tuple[int, int]]) -> np.ndarray:
        V = np.zeros((self.basis.dim, self.basis.dim), dtype=complex)
        for (i, j) in pairs:
            V += self._E(i, j)
        return V

    def _build_Vu(self) -> np.ndarray:
        uL = self.basis.indices("uL")
        uR = self.basis.indices("uR")
        pairs = list(zip(uR, uL))  # uR <- uL (color-preserving)
        return self._sum_E(pairs)

    def _build_Vd(self) -> np.ndarray:
        dL = self.basis.indices("dL")
        dR = self.basis.indices("dR")
        pairs = list(zip(dR, dL))
        return self._sum_E(pairs)

    def _build_Ve(self) -> np.ndarray:
        eL = self.basis.indices("eL")[0]
        eR = self.basis.indices("eR")[0]
        return self._sum_E([(eR, eL)])

    def _build_Vnu(self) -> np.ndarray:
        nuL = self.basis.indices("nuL")[0]
        nuR = self.basis.indices("nuR")[0]
        return self._sum_E([(nuR, nuL)])

    def _build_WR(self) -> np.ndarray:
        nuR = self.basis.indices("nuR")[0]
        nuR_c = self.basis.indices("nuR_c")[0]
        return self._sum_E([(nuR_c, nuR)])

    # ---- conjugate maps via J_SM ----

    def conjugate_map(self, V: np.ndarray, J_SM: AntiUnitary) -> np.ndarray:
        """V^c = J_SM V J_SM^{-1}."""
        return J_SM.conj_op(V)


# ============================================================
#  D_int builder (32N x 32N sparse block matrix)
# ============================================================

class DIntBuilder:
    """
    Assemble D_int(Y[K]) on H_SM ⊗ H_flav:
      D_int = Dirac(part) + Dirac(conj) + Majorana
    using channel maps and lifted textures.
    """

    def __init__(self, basis: SMBasis32, maps: ChannelMaps) -> None:
        self.basis = basis
        self.maps = maps

    def build_D_int(self, textures: TexturePack, J_SM: AntiUnitary, N: int) -> np.ndarray:
        Yu, Yd, Ye, Ynu = textures.Yu, textures.Yd, textures.Ye, textures.Ynu
        MR = textures.MR

        D = np.zeros((32 * N, 32 * N), dtype=complex)

        # --- helpers: place block in SM slot (i,j) ---
        def slot(i: int, j: int, B: np.ndarray) -> None:
            put_block(D, i, j, B)

        # --- particle Dirac blocks (indices explicit, robust and sparse) ---
        uL = self.basis.indices("uL")
        uR = self.basis.indices("uR")
        dL = self.basis.indices("dL")
        dR = self.basis.indices("dR")
        eL = self.basis.indices("eL")[0]
        eR = self.basis.indices("eR")[0]
        nuL = self.basis.indices("nuL")[0]
        nuR = self.basis.indices("nuR")[0]

        for k in range(3):
            slot(uR[k], uL[k], Yu)
            slot(uL[k], uR[k], Yu.conj().T)

            slot(dR[k], dL[k], Yd)
            slot(dL[k], dR[k], Yd.conj().T)

        slot(eR, eL, Ye)
        slot(eL, eR, Ye.conj().T)

        slot(nuR, nuL, Ynu)
        slot(nuL, nuR, Ynu.conj().T)

        # --- conjugate Dirac blocks (+16 shift), code matches paper lemma ---
        for k in range(3):
            slot(uL[k] + 16, uR[k] + 16, Yu.conj())
            slot(uR[k] + 16, uL[k] + 16, Yu.T)

            slot(dL[k] + 16, dR[k] + 16, Yd.conj())
            slot(dR[k] + 16, dL[k] + 16, Yd.T)

        slot(eL + 16, eR + 16, Ye.conj())
        slot(eR + 16, eL + 16, Ye.T)

        slot(nuL + 16, nuR + 16, Ynu.conj())
        slot(nuR + 16, nuL + 16, Ynu.T)

        # --- Majorana block nu_R <-> nu_R^c ---
        if MR is not None:
            slot(nuR + 16, nuR, MR)
            slot(nuR, nuR + 16, MR.conj().T)

        return D

    def check_evenness(self, D_int: np.ndarray, gamma_SM_lifted: np.ndarray, atol: float = 1e-10) -> Dict[str, float]:
        """
        v4.0 hard checks:
          (i) self-adjointness
          (ii) oddness: {gamma_SM⊗I, D_int} = 0
        """
        sa = op_norm(D_int - D_int.conj().T)
        odd = op_norm(anti_comm(gamma_SM_lifted, D_int))
        return {"sa": float(sa), "odd": float(odd)}


# ============================================================
#  Geometry factor hooks (D_geom, gamma_geom, J_geom)
# ============================================================

class GeometryFactor:
    """
    Minimal hooks for the geometric triple.
    In production, D_geom is an operator on H_geom with compact resolvent, etc.
    For numerical experiments, you can use a finite truncation.
    """

    def __init__(self, modes: np.ndarray) -> None:
        self.modes = np.asarray(modes, dtype=float)

        # finite truncation diagonal Φ|n> = n|n>
        self.D_geom = np.diag(self.modes).astype(complex)

        # simplest chirality on this truncation; replace with the true γ_geom in your model
        self.gamma_geom = np.diag(np.sign(self.modes + 1e-18)).astype(complex)

        # placeholder J_geom
        self.J_geom = AntiUnitary(U=np.eye(self.D_geom.shape[0], dtype=complex))

    def dim(self) -> int:
        return int(self.D_geom.shape[0])


# ============================================================
#  Even product Dirac operator D = D_geom⊗1 + γ_geom⊗D_int
# ============================================================

class EvenProductDirac:
    """
    Build the v4.0 even product Dirac on:
      H = H_geom ⊗ (H_SM ⊗ H_flav)
    """

    def __init__(self, geom: GeometryFactor, basis: SMBasis32, params: V40Params) -> None:
        self.geom = geom
        self.basis = basis
        self.params = params

    def build_D(self, D_int: np.ndarray) -> np.ndarray:
        dG = self.geom.D_geom
        gG = self.geom.gamma_geom
        I_int = np.eye(D_int.shape[0], dtype=complex)
        return np.kron(dG, I_int) + np.kron(gG, D_int)

    def build_Gamma_total(self, N: int) -> np.ndarray:
        """
        Total grading Γ = γ_geom ⊗ γ_SM ⊗ 1_N
        """
        gamma_SM_lift = self.basis.lift_gamma(N)
        return np.kron(self.geom.gamma_geom, gamma_SM_lift)


# ============================================================
#  Representations π(a) and generators for A_SM (minimal scaffold)
# ============================================================

class SMRepresentation:
    """
    Minimal representation layer needed to run commutator gates.
    In production, π is the actual NCG SM representation on H_SM.
    Here we provide an interface and placeholders.
    """

    def __init__(self, basis: SMBasis32) -> None:
        self.basis = basis

    def _build_index_map(self) -> Dict[str, List[int]]:
        uL = [0, 1, 2]
        dL = [3, 4, 5]
        nuL = [6]
        eL = [7]
        uR = [8, 9, 10]
        dR = [11, 12, 13]
        eR = [14]
        nuR = [15]

        uL_c = [i + 16 for i in uL]
        dL_c = [i + 16 for i in dL]
        nuL_c = [i + 16 for i in nuL]
        eL_c = [i + 16 for i in eL]
        uR_c = [i + 16 for i in uR]
        dR_c = [i + 16 for i in dR]
        eR_c = [i + 16 for i in eR]
        nuR_c = [i + 16 for i in nuR]

        # multiplet aliases (particle)
        Q_L = uL + dL
        L_L = nuL + eL
        Q_R = uR + dR
        L_R = eR + nuR  # (singlets)

        # multiplet aliases (conjugate)
        Q_L_c = uL_c + dL_c
        L_L_c = nuL_c + eL_c
        Q_R_c = uR_c + dR_c
        L_R_c = eR_c + nuR_c

        return {
            # fine-grained
            "uL": uL, "dL": dL, "nuL": nuL, "eL": eL,
            "uR": uR, "dR": dR, "eR": eR, "nuR": nuR,
            "uL_c": uL_c, "dL_c": dL_c, "nuL_c": nuL_c, "eL_c": eL_c,
            "uR_c": uR_c, "dR_c": dR_c, "eR_c": eR_c, "nuR_c": nuR_c,

            # aliases used by SMRepresentation / paper language
            "Q_L": Q_L,
            "L_L": L_L,
            "Q_R": Q_R,
            "L_R": L_R,
            "Q_L_c": Q_L_c,
            "L_L_c": L_L_c,
            "Q_R_c": Q_R_c,
            "L_R_c": L_R_c,

            # bookkeeping
            "part": list(range(0, 16)),
            "conj": list(range(16, 32)),
            "L_part": list(range(0, 8)),
            "R_part": list(range(8, 16)),
            "L_conj": list(range(24, 32)),
            "R_conj": list(range(16, 24)),
        }

    def generators(self) -> List[np.ndarray]:
        """
        Return a finite generating/spanning set G_A ⊂ π(A_SM).
        MUST be replaced by the true SM algebra reps for real gates.
        """
        # placeholder: projectors onto major summands
        return [
            self.basis.P["uL"] + self.basis.P["dL"],  # Q_L
            self.basis.P["nuL"] + self.basis.P["eL"],  # L_L
        ]

    # Example placeholders; in v4.0 you implement actual SM algebra action
    @property
    def P_Q_L(self) -> np.ndarray:
        # You would add these keys if you expose those subspaces explicitly.
        raise NotImplementedError


# ============================================================
#  One-forms Ω^1_D(A) and inner fluctuations
# ============================================================

class InnerFluctuations:
    """
    Build one-forms A = Σ π(a_i)[D,π(b_i)] and D_A = D + A + JAJ^{-1}.
    """

    def __init__(self, J_total: AntiUnitary) -> None:
        self.J_total = J_total

    def one_form(self, D: np.ndarray, pi_a: np.ndarray, pi_b: np.ndarray) -> np.ndarray:
        return pi_a @ comm(D, pi_b)

    def build_A(self, D: np.ndarray, pairs: Sequence[Tuple[np.ndarray, np.ndarray]]) -> np.ndarray:
        A = np.zeros_like(D, dtype=complex)
        for (pa, pb) in pairs:
            A += self.one_form(D, pa, pb)
        return A

    def fluctuate(self, D: np.ndarray, A: np.ndarray, hermitize: bool = True) -> np.ndarray:
        A_use = hermitian_part(A) if hermitize else A
        return D + A_use + self.J_total.conj_op(A_use)


# ============================================================
#  v4.0 production gates (order-zero, first-order, stability)
# ============================================================

@dataclass
class GateReport:
    order0: float
    first: float
    self_adj: float
    odd_int: float
    note: str = ""


class V4Gates:
    """
    Implements Section 10 production commutator gates.
    """

    def __init__(self, eps0: float, eps1: float, eps_sa: float, epsA: float) -> None:
        self.eps0 = eps0
        self.eps1 = eps1
        self.eps_sa = eps_sa
        self.epsA = epsA

    def gate_order0(self, pis: List[np.ndarray], J_total: AntiUnitary) -> float:
        mx = 0.0
        for a in pis:
            for b in pis:
                mx = max(mx, op_norm(comm(a, J_total.conj_op(b))))
        return float(mx)

    def gate_first(self, D: np.ndarray, pis: List[np.ndarray], J_total: AntiUnitary) -> float:
        mx = 0.0
        for a in pis:
            Da = comm(D, a)
            for b in pis:
                mx = max(mx, op_norm(comm(Da, J_total.conj_op(b))))
        return float(mx)

    def gate_selfadj(self, X: np.ndarray) -> float:
        return float(op_norm(X - X.conj().T))

    def pass_fail(self, report: GateReport) -> bool:
        return (report.order0 <= self.eps0 and
                report.first <= self.eps1 and
                report.self_adj <= self.eps_sa and
                report.odd_int <= 1e-10)  # oddness is a hard axiom check


# ============================================================
#  Texture generation hooks (your current v3.3 kernel pipeline)
# ============================================================

class TextureGenerator:
    """
    Wrap your existing K→flow→P_C360→S compression machinery as a texture generator:
      returns (Yu,Yd,Ye,Ynu,MR).

    In the short term you may set Yu=Yd=0 and use Ye,Ynu from different geometry weights.
    """

    def __init__(self, cfg: Any) -> None:
        self.cfg = cfg

    def build_textures(self, params: V40Params) -> TexturePack:
        """
        MUST return NxN matrices for each sector.
        Stub: return identity textures.
        """
        N = params.N_flav
        Yu = np.zeros((N, N), dtype=complex)
        Yd = np.zeros((N, N), dtype=complex)
        Ye = np.eye(N, dtype=complex)
        Ynu = np.eye(N, dtype=complex)
        MR = params.seesaw_M_GeV * np.eye(N, dtype=complex)
        return TexturePack(Yu=Yu, Yd=Yd, Ye=Ye, Ynu=Ynu, MR=MR)


# ============================================================
#  v4.0 Engine: orchestrates everything
# ============================================================

class AlignmentV40Engine:
    """
    One-stop v4.0 engine:
      - builds textures
      - builds D_int
      - builds even-product D
      - (optionally) builds fluctuations D_A
      - runs gates
      - exposes sector extraction utilities
    """

    def __init__(
        self,
        geom: GeometryFactor,
        params: V40Params,
        texture_gen: TextureGenerator,
    ) -> None:
        self.params = params
        self.geom = geom
        self.texture_gen = texture_gen

        # fixed structures
        self.basis = SMBasis32()
        self.maps = ChannelMaps(self.basis)

        self.J_SM = RealStructureFactory.build_J_SM(self.basis)
        self.J_flav = RealStructureFactory.build_J_flav(params.N_flav)
        self.J_int = RealStructureFactory.kron(self.J_SM, self.J_flav)

        # total J on full Hilbert space H_geom ⊗ H_int
        self.J_total = RealStructureFactory.kron(self.geom.J_geom, self.J_int)

        self.dint_builder = DIntBuilder(self.basis, self.maps)
        self.product_dirac = EvenProductDirac(self.geom, self.basis, params)

        # representation layer (must be replaced with true π(A_SM))
        self.pi_SM = SMRepresentation(self.basis)

        # gates
        self.gates = V4Gates(
            eps0=params.eps_order0,
            eps1=params.eps_first,
            eps_sa=params.eps_sa,
            epsA=params.eps_fluct
        )

    # -------- build core objects --------

    def build_textures(self) -> TexturePack:
        return self.texture_gen.build_textures(self.params)

    def build_D_int(self, textures: TexturePack) -> np.ndarray:
        return self.dint_builder.build_D_int(textures, J_SM=self.J_SM, N=self.params.N_flav)

    def build_D(self, D_int: np.ndarray) -> np.ndarray:
        return self.product_dirac.build_D(D_int)

    def build_fluctuated(self, D: np.ndarray, one_form_pairs: Sequence[Tuple[np.ndarray, np.ndarray]]) -> np.ndarray:
        fl = InnerFluctuations(self.J_total)
        A = fl.build_A(D, one_form_pairs)
        return fl.fluctuate(D, A, hermitize=True)

    # -------- diagnostics / gates --------

    def run_axiom_checks(self, D_int: np.ndarray) -> Dict[str, float]:
        gamma_SM_lift = self.basis.lift_gamma(self.params.N_flav)
        return self.dint_builder.check_evenness(D_int, gamma_SM_lifted=gamma_SM_lift)

    def run_commutator_gates(self, D: np.ndarray) -> GateReport:
        pis = self._lift_pi_generators_to_full()
        order0 = self.gates.gate_order0(pis, self.J_total)
        first = self.gates.gate_first(D, pis, self.J_total)
        sa = self.gates.gate_selfadj(D)
        # oddness is checked at D_int level; include placeholder here
        odd_int = 0.0
        return GateReport(order0=order0, first=first, self_adj=sa, odd_int=odd_int)

    def _lift_pi_generators_to_full(self) -> List[np.ndarray]:
        """
        Lift π(A) generators to the full product Hilbert space:
          H_geom ⊗ H_SM ⊗ H_flav
        In v4.0: π(a) acts trivially on flavor.
        """
        N = self.params.N_flav
        gens_SM = self.pi_SM.generators()  # 32x32
        gens_int = [np.kron(g, np.eye(N, dtype=complex)) for g in gens_SM]  # 32N x 32N

        # lift to full space by kron with I_geom
        I_geom = np.eye(self.geom.dim(), dtype=complex)
        gens_full = [np.kron(I_geom, g_int) for g_int in gens_int]
        return gens_full

    # -------- sector extraction utilities (projector-defined) --------

    def sector_projector(self, key: str) -> np.ndarray:
        """Return lifted projector (32N x 32N) onto a sector key in basis.P."""
        P32 = self.basis.P[key]
        return self.basis.lift_projector(P32, self.params.N_flav)

    def extract_lr_block(self, D_int: np.ndarray, left_key: str, right_key: str) -> np.ndarray:
        """Projector-defined block extraction: P_L D_int P_R."""
        PL = self.sector_projector(left_key)
        PR = self.sector_projector(right_key)
        return PL @ D_int @ PR

    # -------- end-to-end run --------

    def run(self) -> Dict[str, Any]:
        textures = self.build_textures()
        D_int = self.build_D_int(textures)
        D = self.build_D(D_int)

        ax = self.run_axiom_checks(D_int)

        gates = self.run_commutator_gates(D)

        return {
            "textures": textures,
            "D_int": D_int,
            "D": D,
            "axioms": ax,
            "gates": gates,
        }


# ============================================================
#  Optional: Spectral action evaluators (conceptual / numerical)
# ============================================================

class SpectralAction:
    """
    Bosonic spectral action approximations on finite truncations.
    In production you may:
      - compute Tr f(D_A^2/Λ^2) via eigenvalues
      - use heat-kernel coefficients if available
    """

    def __init__(self, cutoff_Lambda: float, f: Callable[[np.ndarray], np.ndarray]) -> None:
        self.Lambda = float(cutoff_Lambda)
        self.f = f

    def trace_f(self, D: np.ndarray) -> float:
        """Compute Tr f(D^2/Λ^2) by diagonalizing D (finite truncation)."""
        evals = eigh(hermitian_part(D))[0]
        x = (np.real(evals) ** 2) / (self.Lambda ** 2)
        return float(np.sum(self.f(x)))

    def fermionic_term(self, psi: np.ndarray, D: np.ndarray) -> complex:
        """<psi, D psi>."""
        return complex(np.vdot(psi, D @ psi))


# ============================================================
#  Example wiring (scaffold)
# ============================================================

def example_build_engine() -> AlignmentV40Engine:
    params = V40Params(N_flav=3)
    geom = GeometryFactor(modes=np.array([1, 2, 3, 4, 5], dtype=float))
    texgen = TextureGenerator(cfg={})
    return AlignmentV40Engine(geom=geom, params=params, texture_gen=texgen)

def main() -> None:
    eng = example_build_engine()
    out = eng.run()

    ax = out["axioms"]
    print("Axiom checks (D_int):", ax)

    gates: GateReport = out["gates"]
    print("Gate report:", gates)

if __name__ == "__main__":
    main()
