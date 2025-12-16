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
        self.state_labels = self._build_state_labels()  # 32 labels in canonical order
        self._label_to_index = {lbl: i for i, lbl in enumerate(self.state_labels)}

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


    def label(self, i: int) -> str:
        """Human-readable label for basis index i (0..31)."""
        return self.state_labels[i]

    def index_of(self, label: str) -> int:
        """Inverse lookup: label -> index."""
        return int(self._label_to_index[label])

    def _build_state_labels(self) -> List[str]:
        """
        Canonical 32-state labels matching the fixed ordering:
          0..15  : particles
          16..31 : charge-conjugates (suffix _c)
        """
        # particles
        labels: List[str] = [
            "uL_r","uL_g","uL_b",
            "dL_r","dL_g","dL_b",
            "nuL","eL",
            "uR_r","uR_g","uR_b",
            "dR_r","dR_g","dR_b",
            "eR","nuR",
        ]
        # conjugates
        labels += [f"{x}_c" for x in labels]
        return labels

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
        Concrete finite real structure for the 32-state basis.

        We choose the standard "swap" charge-conjugation that maps each particle
        basis vector to its conjugate partner and vice-versa:
            J |i> = |i+16>   for i=0..15
            J |i> = |i-16>   for i=16..31

        Implemented as J = U K with:
            U = [[0, I_16],
                 [I_16, 0]]

        This choice is unitary, involutive (J^2 = +1 in this basis), and matches
        the fixed ordering used by SMBasis32.
        """
        if basis.dim != 32:
            raise ValueError(f"Expected SMBasis32 dim=32, got {basis.dim}")

        I16 = np.eye(16, dtype=complex)
        Z16 = np.zeros((16, 16), dtype=complex)
        U_SM = np.block([[Z16, I16],
                         [I16, Z16]])
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



    def validate(self, atol: float = 1e-12) -> Dict[str, float]:
        """
        Sanity checks for partial isometries:
          Vu†Vu = P_uL, VuVu† = P_uR, etc.
        Returns operator-norm deviations for each identity.
        """
        b = self.basis
        out: Dict[str, float] = {}

        def dev(A: np.ndarray, B: np.ndarray) -> float:
            return float(op_norm(A - B))

        Vu, Vd, Ve, Vnu, WR = self.Vu, self.Vd, self.Ve, self.Vnu, self.WR

        out["Vu_dagVu"] = dev(Vu.conj().T @ Vu, b.P["uL"])
        out["VuVu_dag"] = dev(Vu @ Vu.conj().T, b.P["uR"])

        out["Vd_dagVd"] = dev(Vd.conj().T @ Vd, b.P["dL"])
        out["VdVd_dag"] = dev(Vd @ Vd.conj().T, b.P["dR"])

        out["Ve_dagVe"] = dev(Ve.conj().T @ Ve, b.P["eL"])
        out["VeVe_dag"] = dev(Ve @ Ve.conj().T, b.P["eR"])

        out["Vnu_dagVnu"] = dev(Vnu.conj().T @ Vnu, b.P["nuL"])
        out["VnuVnu_dag"] = dev(Vnu @ Vnu.conj().T, b.P["nuR"])

        # WR is a partial isometry nuR -> nuR_c
        out["WR_dagWR"] = dev(WR.conj().T @ WR, b.P["nuR"])
        out["WRWR_dag"] = dev(WR @ WR.conj().T, b.P["nuR_c"])

        # also assert they're close to partial isometries (V†V and VV† are projectors)
        # (return numeric diagnostics; raising is left to the caller)
        return out

# ============================================================
#  D_int builder (32N x 32N sparse block matrix)
# ============================================================

class DIntBuilder:
    """
    Assemble D_int(Y[K]) on H_SM ⊗ H_flav:

        D_int = Dirac(particles) + Dirac(conjugates) + Majorana

    Conventions are tied to SMBasis32.gamma_SM:
      - "left"  := gamma_SM = +1
      - "right" := gamma_SM = -1

    In the conjugate sector chirality is flipped (as encoded in SMBasis32),
    so the conjugate Dirac blocks use (V^c)† ⊗ \bar{Y} and V^c ⊗ Y^T.
    """

    def __init__(self, basis: SMBasis32, maps: ChannelMaps) -> None:
        self.basis = basis
        self.maps = maps

    # ---------------- internal validators ----------------

    @staticmethod
    def _require_square_NxN(name: str, A: np.ndarray, N: int) -> None:
        if A is None:
            raise ValueError(f"{name} is None (expected {N}x{N} complex array).")
        A = np.asarray(A)
        if A.shape != (N, N):
            raise ValueError(f"{name} must be shape ({N},{N}), got {A.shape}.")

    def _validate_textures(self, textures: TexturePack, N: int) -> None:
        self._require_square_NxN("Yu", textures.Yu, N)
        self._require_square_NxN("Yd", textures.Yd, N)
        self._require_square_NxN("Ye", textures.Ye, N)
        self._require_square_NxN("Ynu", textures.Ynu, N)
        if textures.MR is not None:
            self._require_square_NxN("MR", textures.MR, N)

    # ---------------- build ----------------

    def build_D_int(self, textures: TexturePack, J_SM: AntiUnitary, N: int) -> np.ndarray:
        """
        Return the internal Dirac operator as a dense (32N)x(32N) matrix.
        """
        self._validate_textures(textures, N)

        Yu, Yd, Ye, Ynu = (np.asarray(textures.Yu, dtype=complex),
                          np.asarray(textures.Yd, dtype=complex),
                          np.asarray(textures.Ye, dtype=complex),
                          np.asarray(textures.Ynu, dtype=complex))
        MR = None if textures.MR is None else np.asarray(textures.MR, dtype=complex)

        Vu, Vd, Ve, Vnu, WR = self.maps.Vu, self.maps.Vd, self.maps.Ve, self.maps.Vnu, self.maps.WR

        # Conjugate-channel maps on H_SM: V^c = J_SM V J_SM^{-1}.
        Vu_c = self.maps.conjugate_map(Vu, J_SM)
        Vd_c = self.maps.conjugate_map(Vd, J_SM)
        Ve_c = self.maps.conjugate_map(Ve, J_SM)
        Vnu_c = self.maps.conjugate_map(Vnu, J_SM)

        # --- particle Dirac: (R <- L)⊗Y  +  (L <- R)⊗Y†
        D_part = (
            np.kron(Vu, Yu) + np.kron(Vu.conj().T, Yu.conj().T) +
            np.kron(Vd, Yd) + np.kron(Vd.conj().T, Yd.conj().T) +
            np.kron(Ve, Ye) + np.kron(Ve.conj().T, Ye.conj().T) +
            np.kron(Vnu, Ynu) + np.kron(Vnu.conj().T, Ynu.conj().T)
        )

        # --- conjugate Dirac (lemma-consistent with your explicit +16 version)
        D_conj = (
            np.kron(Vu_c.conj().T, np.conj(Yu)) + np.kron(Vu_c, Yu.T) +
            np.kron(Vd_c.conj().T, np.conj(Yd)) + np.kron(Vd_c, Yd.T) +
            np.kron(Ve_c.conj().T, np.conj(Ye)) + np.kron(Ve_c, Ye.T) +
            np.kron(Vnu_c.conj().T, np.conj(Ynu)) + np.kron(Vnu_c, Ynu.T)
        )

        # --- Majorana: nu_R <-> nu_R^c (WR: nuR_c <- nuR)
        D_M = np.zeros_like(D_part)
        if MR is not None:
            D_M = np.kron(WR, MR) + np.kron(WR.conj().T, MR.conj().T)

        D = np.asarray(D_part + D_conj + D_M, dtype=complex)
        if D.shape != (32 * N, 32 * N):
            raise RuntimeError(f"D_int has wrong shape {D.shape}, expected {(32*N, 32*N)}.")
        return D

    # ---------------- axioms ----------------

    def check_evenness(self, D_int: np.ndarray, gamma_SM_lifted: np.ndarray, atol: float = 1e-10) -> Dict[str, float]:
        """
        v4.0 hard checks:
          (i) self-adjointness
          (ii) oddness: {gamma_SM⊗I, D_int} = 0
        """
        D_int = np.asarray(D_int, dtype=complex)
        gamma_SM_lifted = np.asarray(gamma_SM_lifted, dtype=complex)

        sa = op_norm(D_int - D_int.conj().T)
        odd = op_norm(anti_comm(gamma_SM_lifted, D_int))

        return {"sa": float(sa), "odd": float(odd), "sa_ok": sa <= atol, "odd_ok": odd <= atol}


# ============================================================
#  Geometry factor hooks (D_geom, gamma_geom, J_geom)
# ============================================================

class GeometryFactor:
    """
    Minimal hooks for the geometric triple (H_geom, D_geom, gamma_geom, J_geom).

    Defaults (backward-compatible with the template):
      - D_geom: diagonal truncation diag(modes)
      - gamma_geom: diagonal sign(modes) (zeros treated as +1)
      - J_geom: complex conjugation (AntiUnitary with U = I)

    If you want an even truncation satisfying {D_geom, gamma_geom} = 0, use
    `GeometryFactor.even_truncation(lambdas)`, which builds a chiral doubling.
    """

    def __init__(
        self,
        modes: np.ndarray,
        *,
        gamma: np.ndarray | None = None,
        J: AntiUnitary | None = None,
        zero_sign: float = +1.0,
        atol: float = 1e-12,
    ) -> None:
        self.modes = np.asarray(modes, dtype=float).reshape(-1)

        # Finite truncation diagonal Φ|n> = n|n>
        self.D_geom = np.diag(self.modes).astype(complex)

        # Chirality: default sign(modes) with a stable convention at 0.
        if gamma is None:
            s = np.sign(self.modes)
            if zero_sign not in (+1.0, -1.0):
                raise ValueError("zero_sign must be +1.0 or -1.0")
            s = np.where(s == 0.0, float(zero_sign), s)
            self.gamma_geom = np.diag(s).astype(complex)
        else:
            g = np.asarray(gamma, dtype=complex)
            if g.shape != self.D_geom.shape:
                raise ValueError(f"gamma must have shape {self.D_geom.shape}, got {g.shape}")
            self.gamma_geom = g

        # Real structure (antiunitary): default is complex conjugation in this basis.
        self.J_geom = J if J is not None else AntiUnitary(U=np.eye(self.D_geom.shape[0], dtype=complex))

        self._validate_basics(atol=atol)

    @classmethod
    def even_truncation(
        cls,
        lambdas: np.ndarray,
        *,
        J: AntiUnitary | None = None,
        atol: float = 1e-12,
    ) -> "GeometryFactor":
        """
        Build a chiral-doubled even finite triple:

            H_geom = C^n ⊕ C^n
            D_geom = [[0, Λ],
                      [Λ, 0]]   where Λ = diag(lambdas)
            gamma  = [[ I, 0],
                      [ 0,-I]]

        This satisfies {D_geom, gamma} = 0 (up to atol).
        """
        lam = np.asarray(lambdas, dtype=float).reshape(-1)
        n = lam.shape[0]
        Z = np.zeros((n, n), dtype=complex)
        Lam = np.diag(lam).astype(complex)

        obj = cls.__new__(cls)
        obj.modes = lam.copy()
        obj.D_geom = np.block([[Z, Lam], [Lam, Z]]).astype(complex)
        obj.gamma_geom = np.block([[np.eye(n, dtype=complex), Z], [Z, -np.eye(n, dtype=complex)]]).astype(complex)
        obj.J_geom = J if J is not None else AntiUnitary(U=np.eye(2 * n, dtype=complex))

        obj._validate_basics(atol=atol)

        anti = anti_comm(obj.gamma_geom, obj.D_geom)
        if op_norm(anti) > atol:
            raise ValueError(f"even_truncation: {{gamma_geom, D_geom}} not ~0 (norm={op_norm(anti)})")
        return obj

    def dim(self) -> int:
        return int(self.D_geom.shape[0])

    def lift_D(self, fin_dim: int) -> np.ndarray:
        """Return D_geom ⊗ I_fin."""
        return np.kron(self.D_geom, np.eye(fin_dim, dtype=complex))

    def lift_gamma(self, fin_dim: int) -> np.ndarray:
        """Return gamma_geom ⊗ I_fin."""
        return np.kron(self.gamma_geom, np.eye(fin_dim, dtype=complex))

    def lift_J(self, J_fin: AntiUnitary) -> AntiUnitary:
        """
        Return J_geom ⊗ J_fin as AntiUnitary with U_total = U_geom ⊗ U_fin.
        """
        return AntiUnitary(U=np.kron(self.J_geom.U, J_fin.U))

    def _validate_basics(self, *, atol: float = 1e-12) -> None:
        # D_geom should be self-adjoint.
        if op_norm(self.D_geom - self.D_geom.conj().T) > atol:
            raise ValueError("D_geom must be self-adjoint (within atol).")

        # gamma_geom should be self-adjoint and an involution.
        if op_norm(self.gamma_geom - self.gamma_geom.conj().T) > atol:
            raise ValueError("gamma_geom must be self-adjoint (within atol).")
        I = np.eye(self.gamma_geom.shape[0], dtype=complex)
        if op_norm(self.gamma_geom @ self.gamma_geom - I) > 1e-9:
            raise ValueError("gamma_geom must satisfy gamma_geom^2 = I (within tolerance).")

        # J_geom's unitary part should be unitary.
        U = self.J_geom.U
        if U.shape != self.D_geom.shape:
            raise ValueError(f"J_geom.U must have shape {self.D_geom.shape}, got {U.shape}")
        if op_norm(U.conj().T @ U - I) > 1e-9:
            raise ValueError("J_geom.U must be unitary (within tolerance).")

# ============================================================
#  Even product Dirac operator D = D_geom⊗1 + γ_geom⊗D_int
# ============================================================

class EvenProductDirac:
    """
    Build the v4.0 even product Dirac on:
        H = H_geom ⊗ (H_SM ⊗ H_flav)

    Core constructions:
        D_total  = (D_geom ⊗ I_int) + (gamma_geom ⊗ D_int)
        Gamma    = gamma_geom ⊗ (gamma_SM ⊗ I_N)

    Notes:
      - Oddness of D_total w.r.t Gamma holds if both
            {D_geom, gamma_geom} = 0   and   {D_int, gamma_SM⊗I_N} = 0.
      - This class provides hard shape checks and basic axioms checks.
    """

    def __init__(self, geom: GeometryFactor, basis: SMBasis32, params: V40Params) -> None:
        self.geom = geom
        self.basis = basis
        self.params = params

    # ---------------- builders ----------------

    def build_D(
        self,
        D_int: np.ndarray,
        *,
        scale_geom: float = 1.0,
        scale_int: float = 1.0,
    ) -> np.ndarray:
        """
        Return D_total = scale_geom*(D_geom⊗I) + scale_int*(gamma_geom⊗D_int).
        """
        D_int = np.asarray(D_int, dtype=complex)
        if D_int.ndim != 2 or D_int.shape[0] != D_int.shape[1]:
            raise ValueError(f"D_int must be square, got {D_int.shape}.")

        dG = np.asarray(self.geom.D_geom, dtype=complex)
        gG = np.asarray(self.geom.gamma_geom, dtype=complex)

        I_int = np.eye(D_int.shape[0], dtype=complex)
        D_total = (float(scale_geom) * np.kron(dG, I_int)) + (float(scale_int) * np.kron(gG, D_int))

        expected = (self.geom.dim() * D_int.shape[0], self.geom.dim() * D_int.shape[0])
        if D_total.shape != expected:
            raise RuntimeError(f"D_total has shape {D_total.shape}, expected {expected}.")
        return D_total

    def build_Gamma_total(self, N: int) -> np.ndarray:
        """
        Total grading Γ = γ_geom ⊗ (γ_SM ⊗ 1_N).
        """
        gamma_SM_lift = self.basis.lift_gamma(N)  # (32N)x(32N)
        Gamma = np.kron(np.asarray(self.geom.gamma_geom, dtype=complex), gamma_SM_lift)
        return Gamma

    def build_J_internal(self, J_SM: AntiUnitary, N: int, J_flav: AntiUnitary | None = None) -> AntiUnitary:
        """
        Internal real structure:
            J_int = J_SM ⊗ J_flav
        Default J_flav is complex conjugation on C^N (U = I_N).
        """
        if J_flav is None:
            J_flav = AntiUnitary(U=np.eye(N, dtype=complex))
        return AntiUnitary(U=np.kron(J_SM.U, J_flav.U))

    def build_J_total(self, J_int: AntiUnitary) -> AntiUnitary:
        """
        Total real structure:
            J_total = J_geom ⊗ J_int.
        """
        return self.geom.lift_J(J_int)

    # ---------------- checks ----------------

    def check_axioms(
        self,
        D_total: np.ndarray,
        D_int: np.ndarray,
        N: int,
        *,
        atol: float | None = None,
    ) -> Dict[str, float]:
        """
        Basic even-product hard checks:
          (i)  self-adjointness of D_total
          (ii) oddness: {Gamma_total, D_total} = 0
          (iii) diagnostics:
               {gamma_geom, D_geom} and {gamma_int, D_int}

        Returns numeric norms and boolean flags.
        """
        eps = float(self.params.eps_sa if atol is None else atol)

        D_total = np.asarray(D_total, dtype=complex)
        D_int = np.asarray(D_int, dtype=complex)

        Gamma = self.build_Gamma_total(N)

        sa = op_norm(D_total - D_total.conj().T)
        odd_total = op_norm(anti_comm(Gamma, D_total))

        # diagnostics: the two sufficient conditions for odd_total≈0
        dG = np.asarray(self.geom.D_geom, dtype=complex)
        gG = np.asarray(self.geom.gamma_geom, dtype=complex)
        odd_geom = op_norm(anti_comm(gG, dG))

        gamma_int = self.basis.lift_gamma(N)
        odd_int = op_norm(anti_comm(gamma_int, D_int))

        return {
            "sa_total": float(sa),
            "odd_total": float(odd_total),
            "odd_geom": float(odd_geom),
            "odd_int": float(odd_int),
            "sa_ok": sa <= eps,
            "odd_ok": odd_total <= eps,
        }


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

    # ----------------- full algebra representation π(λ,q,m) -----------------

    @staticmethod
    def _as2(name: str, A: np.ndarray) -> np.ndarray:
        A = np.asarray(A, dtype=complex)
        if A.shape != (2, 2):
            raise ValueError(f"{name} must be 2x2, got {A.shape}")
        return A

    @staticmethod
    def _as3(name: str, A: np.ndarray) -> np.ndarray:
        A = np.asarray(A, dtype=complex)
        if A.shape != (3, 3):
            raise ValueError(f"{name} must be 3x3, got {A.shape}")
        return A

    def pi(self, lam: complex, q: np.ndarray, m: np.ndarray, *, enforce_quaternion: bool = False) -> np.ndarray:
        """
        Algebra representation for a=(lam,q,m) ∈ C ⊕ H ⊕ M3(C).

        Minimal NCG-consistent block action:
          - On Q_L (uL,dL):         q ⊗ m
          - On L_L (nuL,eL):        q
          - On uR,dR:              lam * m
          - On nuR,eR:             lam
          - On conjugates:         complex-conjugate parameters (lam̄, q̄, m̄)

        Returns 32x32 complex matrix in the SMBasis32 ordering.
        """
        lam = complex(lam)
        q = self._as2("q", q)
        m = self._as3("m", m)

        if enforce_quaternion:
            # Optional light check: q should be in the standard complex 2x2 image of H:
            # q = [[a, b], [-b̄, ā]]
            a = q[0, 0]
            b = q[0, 1]
            target = np.array([[a, b], [-np.conj(b), np.conj(a)]], dtype=complex)
            if op_norm(q - target) > 1e-10:
                raise ValueError("q does not look like a quaternion-embedded 2x2 matrix.")

        X = np.zeros((32, 32), dtype=complex)

        # --- particles (0..15) ---
        # Q_L block: (uL_r,uL_g,uL_b,dL_r,dL_g,dL_b)
        idx_QL = self.uL + self.dL
        X[np.ix_(idx_QL, idx_QL)] = np.kron(q, m)

        # L_L block: (nuL,eL)
        idx_LL = [self.nuL[0], self.eL[0]]
        X[np.ix_(idx_LL, idx_LL)] = q

        # uR, dR: lam*m on each color triplet
        X[np.ix_(self.uR, self.uR)] = lam * m
        X[np.ix_(self.dR, self.dR)] = lam * m

        # nuR, eR: lam (singlets)
        X[self.nuR[0], self.nuR[0]] = lam
        X[self.eR[0], self.eR[0]] = lam

        # --- conjugates (16..31): use conjugate parameters ---
        lamc = np.conj(lam)
        qc = np.conj(q)
        mc = np.conj(m)

        idx_QLc = self.uL_c + self.dL_c
        X[np.ix_(idx_QLc, idx_QLc)] = np.kron(qc, mc)

        idx_LLc = [self.nuL_c[0], self.eL_c[0]]
        X[np.ix_(idx_LLc, idx_LLc)] = qc

        X[np.ix_(self.uR_c, self.uR_c)] = lamc * mc
        X[np.ix_(self.dR_c, self.dR_c)] = lamc * mc

        X[self.nuR_c[0], self.nuR_c[0]] = lamc
        X[self.eR_c[0], self.eR_c[0]] = lamc

        return X

    def pi_lifted(self, lam: complex, q: np.ndarray, m: np.ndarray, N: int, *, enforce_quaternion: bool = False) -> np.ndarray:
        """
        Lift π(a) to H_SM ⊗ H_flav: π(a) ⊗ I_N.
        """
        Pi32 = self.pi(lam, q, m, enforce_quaternion=enforce_quaternion)
        return np.kron(Pi32, np.eye(N, dtype=complex))

# ============================================================
#  One-forms Ω^1_D(A) and inner fluctuations
# ============================================================

class InnerFluctuations:
    """
    Build one-forms A = Σ π(a_i)[D, π(b_i)] and D_A = D + A + JAJ^{-1}.

    Hermitization options:
      - hermitize="A":        replace A by (A+A†)/2 before adding JAJ^{-1}
      - hermitize="total":    hermitize the full fluctuation term (A + JAJ^{-1})  [recommended]
      - hermitize=False:      do nothing (debug only)
    """

    def __init__(self, J_total: AntiUnitary) -> None:
        self.J_total = J_total

    # ---------------- core primitives ----------------

    @staticmethod
    def _require_same_shape(name: str, X: np.ndarray, Y: np.ndarray) -> None:
        if X.shape != Y.shape:
            raise ValueError(f"{name}: shape mismatch {X.shape} vs {Y.shape}")

    @staticmethod
    def _require_square(name: str, X: np.ndarray) -> None:
        if X.ndim != 2 or X.shape[0] != X.shape[1]:
            raise ValueError(f"{name} must be square, got {X.shape}")

    def one_form(
        self,
        D: np.ndarray,
        pi_a: np.ndarray,
        pi_b: np.ndarray,
        *,
        remove_trace: bool = False,
    ) -> np.ndarray:
        """
        A(a,b) = π(a)[D,π(b)].

        remove_trace=True subtracts (tr(A)/n) I to keep a traceless contribution (often useful
        when you want to drop U(1) pieces for diagnostics).
        """
        D = np.asarray(D, dtype=complex)
        pa = np.asarray(pi_a, dtype=complex)
        pb = np.asarray(pi_b, dtype=complex)

        self._require_square("D", D)
        self._require_square("pi_a", pa)
        self._require_square("pi_b", pb)
        self._require_same_shape("one_form inputs", D, pa)
        self._require_same_shape("one_form inputs", D, pb)

        A = pa @ comm(D, pb)

        if remove_trace:
            n = A.shape[0]
            tr = np.trace(A) / complex(n)
            A = A - tr * np.eye(n, dtype=complex)

        return A

    def build_A(
        self,
        D: np.ndarray,
        pairs: Sequence[Tuple[np.ndarray, np.ndarray]],
        *,
        weights: Sequence[complex] | None = None,
        remove_trace: bool = False,
    ) -> np.ndarray:
        """
        Sum A = Σ w_i π(a_i)[D,π(b_i)].
        """
        D = np.asarray(D, dtype=complex)
        self._require_square("D", D)

        if weights is not None and len(weights) != len(pairs):
            raise ValueError(f"weights must have same length as pairs: {len(weights)} vs {len(pairs)}")

        A = np.zeros_like(D, dtype=complex)
        for idx, (pa, pb) in enumerate(pairs):
            w = 1.0 if weights is None else complex(weights[idx])
            A += w * self.one_form(D, pa, pb, remove_trace=remove_trace)
        return A

    def fluctuate(
        self,
        D: np.ndarray,
        A: np.ndarray,
        *,
        hermitize: str | bool = "total",
        stabilize_output: bool = True,
    ) -> np.ndarray:
        """
        D_A = D + A + JAJ^{-1}, with optional hermitization policy.

        stabilize_output=True forces final D_A -> (D_A + D_A†)/2 to suppress tiny numerical skew.
        """
        D = np.asarray(D, dtype=complex)
        A = np.asarray(A, dtype=complex)

        self._require_square("D", D)
        self._require_square("A", A)
        self._require_same_shape("fluctuate inputs", D, A)

        if hermitize is True:
            hermitize = "A"
        if hermitize not in (False, "A", "total"):
            raise ValueError("hermitize must be False, 'A', or 'total'.")

        if hermitize == "A":
            A_use = hermitian_part(A)
            DA = D + A_use + self.J_total.conj_op(A_use)
        else:
            JA = self.J_total.conj_op(A)
            total_fluct = A + JA
            if hermitize == "total":
                total_fluct = hermitian_part(total_fluct)
            DA = D + total_fluct

        return hermitian_part(DA) if stabilize_output else DA

    # ---------------- convenience ----------------

    def fluctuate_from_pairs(
        self,
        D: np.ndarray,
        pairs: Sequence[Tuple[np.ndarray, np.ndarray]],
        *,
        weights: Sequence[complex] | None = None,
        remove_trace: bool = False,
        hermitize: str | bool = "total",
        stabilize_output: bool = True,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Convenience: build A from pairs then return (D_A, A).
        """
        A = self.build_A(D, pairs, weights=weights, remove_trace=remove_trace)
        D_A = self.fluctuate(D, A, hermitize=hermitize, stabilize_output=stabilize_output)
        return D_A, A

    # ---------------- diagnostics ----------------

    def check_self_adjoint(self, X: np.ndarray) -> float:
        X = np.asarray(X, dtype=complex)
        self._require_square("X", X)
        return float(op_norm(X - X.conj().T))

    def check_oddness(self, X: np.ndarray, Gamma_total: np.ndarray) -> float:
        X = np.asarray(X, dtype=complex)
        G = np.asarray(Gamma_total, dtype=complex)
        self._require_square("X", X)
        self._require_square("Gamma_total", G)
        self._require_same_shape("oddness inputs", X, G)
        return float(op_norm(anti_comm(G, X)))

    def check_commutator(self, D: np.ndarray, pi: np.ndarray) -> float:
        D = np.asarray(D, dtype=complex)
        pi = np.asarray(pi, dtype=complex)
        self._require_square("D", D)
        self._require_square("pi", pi)
        self._require_same_shape("commutator inputs", D, pi)
        return float(op_norm(comm(D, pi)))

# ============================================================
#  v4.0 production gates (order-zero, first-order, stability)
# ============================================================

@dataclass(frozen=True)
class GateReport:
    order0: float
    first: float
    self_adj: float
    odd_int: float
    # optional extras
    order0_pair: Optional[tuple[int, int]] = None
    first_pair: Optional[tuple[int, int]] = None




class V4Gates:
    """
    Implements Section 10 production commutator gates.

    Gate definitions:
      - order-0:  max_{a,b in Π} || [π(a), J π(b) J^{-1}] ||
      - first-order: max_{a,b in Π} || [[D, π(a)], J π(b) J^{-1}] ||
      - self-adjoint: ||X - X†|| (for whatever X you're checking)
      - oddness: ||{Γ, D}|| or ||{γ_int, D_int}|| depending on what you pass in

    Production details:
      * Precomputes Jπ(b)J^{-1} and [D,π(a)] once.
      * Optional early-exit if a gate exceeds its epsilon.
      * Optional sampling over Π for quick runs.
    """

    def __init__(self, eps0: float, eps1: float, eps_sa: float, epsA: float) -> None:
        self.eps0 = float(eps0)
        self.eps1 = float(eps1)
        self.eps_sa = float(eps_sa)
        self.epsA = float(epsA)

    # ---------- basic checks ----------

    @staticmethod
    def _require_square(name: str, X: np.ndarray) -> None:
        if X.ndim != 2 or X.shape[0] != X.shape[1]:
            raise ValueError(f"{name} must be square, got {X.shape}")

    @staticmethod
    def _require_same_dim(name: str, A: np.ndarray, B: np.ndarray) -> None:
        if A.shape != B.shape:
            raise ValueError(f"{name}: shape mismatch {A.shape} vs {B.shape}")

    def gate_selfadj(self, X: np.ndarray) -> float:
        X = np.asarray(X, dtype=complex)
        self._require_square("X", X)
        return float(op_norm(X - X.conj().T))

    def _sanitize_sample(self, n: int, sample: Optional[Sequence[int]]) -> List[int]:
        """
        Normalize and validate sample indices for a list of length n.
        - Supports negative indices like Python (-1 = last).
        - Drops out-of-range entries.
        - Preserves order and removes duplicates.
        """
        if n <= 0:
            raise ValueError("pis must be non-empty")

        if sample is None:
            return list(range(n))

        out: List[int] = []
        seen = set()
        for s in sample:
            i = int(s)
            if -n <= i < 0:
                i = i + n
            if 0 <= i < n and i not in seen:
                out.append(i)
                seen.add(i)

        if not out:
            raise ValueError(f"sample has no valid indices for pis length {n}. "
                             f"Valid range is 0..{n - 1} (or -{n}..-1). Got: {list(sample)}")
        return out

    # ---------- Section 10 gates ----------

    def gate_order0(
        self,
        pis: List[np.ndarray],
        J_total: AntiUnitary,
        *,
        early_exit: bool = False,
        eps: Optional[float] = None,
        sample: Optional[Sequence[int]] = None,
        return_witness: bool = False,
    ) -> float | tuple[float, tuple[int, int]]:
        """
        max_{a,b} || [a, J(b)] || where J(b)=J_total.conj_op(b)

        sample: optional list of indices into pis to restrict the sweep.
        """
        if len(pis) == 0:
            raise ValueError("pis must be non-empty")

        n = len(pis)
        idxs = self._sanitize_sample(n, sample)

        mats = [np.asarray(pis[i], dtype=complex) for i in idxs]
        for k, A in enumerate(mats):
            self._require_square(f"pis[{idxs[k]}]", A)
            if k > 0:
                self._require_same_dim("pis", mats[0], A)

        # Precompute J(b)
        Jpis = [J_total.conj_op(B) for B in mats]

        thresh = self.eps0 if eps is None else float(eps)
        mx = 0.0
        w = (idxs[0], idxs[0])

        for i, A in enumerate(mats):
            for j, JB in enumerate(Jpis):
                v = float(op_norm(comm(A, JB)))
                if v > mx:
                    mx = v
                    w = (idxs[i], idxs[j])
                    if early_exit and mx > thresh:
                        return (mx, w) if return_witness else mx

        return (mx, w) if return_witness else mx

    def gate_first(
        self,
        D: np.ndarray,
        pis: List[np.ndarray],
        J_total: AntiUnitary,
        *,
        early_exit: bool = False,
        eps: Optional[float] = None,
        sample: Optional[Sequence[int]] = None,
        return_witness: bool = False,
    ) -> float | tuple[float, tuple[int, int]]:
        """
        max_{a,b} || [[D,a], J(b)] ||.
        """
        D = np.asarray(D, dtype=complex)
        self._require_square("D", D)

        if len(pis) == 0:
            raise ValueError("pis must be non-empty")

        n = len(pis)
        idxs = self._sanitize_sample(n, sample)

        mats = [np.asarray(pis[i], dtype=complex) for i in idxs]
        for k, A in enumerate(mats):
            self._require_square(f"pis[{idxs[k]}]", A)
            self._require_same_dim("D vs pi", D, A)

        # Precompute [D, a] and J(b)
        Das = [comm(D, A) for A in mats]
        Jpis = [J_total.conj_op(B) for B in mats]

        thresh = self.eps1 if eps is None else float(eps)
        mx = 0.0
        w = (idxs[0], idxs[0])

        for i, Da in enumerate(Das):
            for j, JB in enumerate(Jpis):
                v = float(op_norm(comm(Da, JB)))
                if v > mx:
                    mx = v
                    w = (idxs[i], idxs[j])
                    if early_exit and mx > thresh:
                        return (mx, w) if return_witness else mx

        return (mx, w) if return_witness else mx

    # ---------- report helpers ----------

    def build_report(
        self,
        *,
        order0: float,
        first: float,
        self_adj: float,
        odd_int: float,
        order0_pair: Optional[tuple[int, int]] = None,
        first_pair: Optional[tuple[int, int]] = None,
    ) -> GateReport:
        return GateReport(
            order0=float(order0),
            first=float(first),
            self_adj=float(self_adj),
            odd_int=float(odd_int),
            order0_pair=order0_pair,
            first_pair=first_pair,
        )

    def pass_fail(self, report: GateReport) -> bool:
        # oddness is a hard axiom check: use epsA as the “hard” tolerance unless you prefer eps_sa
        return (
            report.order0 <= self.eps0 and
            report.first <= self.eps1 and
            report.self_adj <= self.eps_sa and
            report.odd_int <= self.epsA
        )
# ============================================================
#  Texture generation hooks (your current v3.3 kernel pipeline)
# ============================================================

class TextureGenerator:
    """
    Wrap your existing K→flow→P_C360→S compression machinery as a texture generator.

    Supported cfg styles (pick one):

    1) cfg.builder: callable(params) -> TexturePack | dict with keys Yu,Yd,Ye,Ynu,MR
    2) cfg.npz_path: path to .npz containing arrays named Yu,Yd,Ye,Ynu,MR
    3) cfg.mode: one of:
         - "identity_leptons" (default): Yu=Yd=0, Ye=I, Ynu=I, MR=seesaw*I
         - "identity_all": Yu=Yd=Ye=Ynu=I, MR=seesaw*I
         - "zeros_all": Yu=Yd=Ye=Ynu=0, MR=None
         - "random": random complex textures (seedable), MR symmetric

    Common cfg knobs:
      - disable_quarks: bool (forces Yu=Yd=0)
      - scales: dict with optional keys {"Yu","Yd","Ye","Ynu","MR"}
      - enforce_MR_symmetric: bool (default True)
      - seed: int (for mode="random")
      - normalize: bool (rescale each Yukawa so ||Y||_2 ~= 1 before scales)
    """

    def __init__(self, cfg: Any) -> None:
        self.cfg = cfg

    # ----------------- helpers -----------------

    @staticmethod
    def _as_mat(name: str, X: Any, N: int, allow_none: bool = False) -> np.ndarray | None:
        if X is None:
            if allow_none:
                return None
            raise ValueError(f"{name} is None, expected {N}x{N} matrix.")
        A = np.asarray(X, dtype=complex)
        if A.shape != (N, N):
            raise ValueError(f"{name} must be shape ({N},{N}), got {A.shape}.")
        return A

    @staticmethod
    def _symmetrize_majorana(M: np.ndarray) -> np.ndarray:
        # Majorana mass matrix is typically complex symmetric: M = M^T.
        return 0.5 * (M + M.T)

    @staticmethod
    def _opnorm_rescale(Y: np.ndarray, target: float = 1.0, eps: float = 1e-14) -> np.ndarray:
        n = op_norm(Y)
        if n < eps:
            return Y
        return (target / n) * Y

    def _get_cfg(self, key: str, default: Any = None) -> Any:
        # Supports dict-like cfg or object-with-attrs cfg
        if isinstance(self.cfg, dict):
            return self.cfg.get(key, default)
        return getattr(self.cfg, key, default)

    # ----------------- main entry -----------------

    def build_textures(self, params: V40Params) -> TexturePack:
        N = int(params.N_flav)

        # Highest priority: a callable builder (your v3.3 pipeline hook)
        builder = self._get_cfg("builder", None)
        if builder is not None:
            out = builder(params)
            if isinstance(out, TexturePack):
                tp = out
            elif isinstance(out, dict):
                tp = TexturePack(
                    Yu=out.get("Yu"), Yd=out.get("Yd"), Ye=out.get("Ye"),
                    Ynu=out.get("Ynu"), MR=out.get("MR")
                )
            else:
                raise ValueError("cfg.builder must return TexturePack or dict.")
        else:
            # Next: load from npz if provided
            npz_path = self._get_cfg("npz_path", None)
            if npz_path is not None:
                data = np.load(npz_path, allow_pickle=False)
                tp = TexturePack(
                    Yu=data.get("Yu"), Yd=data.get("Yd"), Ye=data.get("Ye"),
                    Ynu=data.get("Ynu"), MR=data.get("MR", None)
                )
            else:
                # Fallback: built-in modes
                mode = self._get_cfg("mode", "identity_leptons")

                if mode == "identity_leptons":
                    Yu = np.zeros((N, N), dtype=complex)
                    Yd = np.zeros((N, N), dtype=complex)
                    Ye = np.eye(N, dtype=complex)
                    Ynu = np.eye(N, dtype=complex)
                    MR = complex(params.seesaw_M_GeV) * np.eye(N, dtype=complex)
                elif mode == "identity_all":
                    Yu = np.eye(N, dtype=complex)
                    Yd = np.eye(N, dtype=complex)
                    Ye = np.eye(N, dtype=complex)
                    Ynu = np.eye(N, dtype=complex)
                    MR = complex(params.seesaw_M_GeV) * np.eye(N, dtype=complex)
                elif mode == "zeros_all":
                    Yu = np.zeros((N, N), dtype=complex)
                    Yd = np.zeros((N, N), dtype=complex)
                    Ye = np.zeros((N, N), dtype=complex)
                    Ynu = np.zeros((N, N), dtype=complex)
                    MR = None
                elif mode == "random":
                    seed = self._get_cfg("seed", 0)
                    rng = np.random.default_rng(seed)
                    def rc():
                        return rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N))
                    Yu, Yd, Ye, Ynu = rc(), rc(), rc(), rc()
                    MR = rc()
                    MR = self._symmetrize_majorana(MR)
                    MR *= complex(params.seesaw_M_GeV) / max(op_norm(MR), 1e-14)
                else:
                    raise ValueError(f"Unknown cfg.mode={mode!r}")

                tp = TexturePack(Yu=Yu, Yd=Yd, Ye=Ye, Ynu=Ynu, MR=MR)

        # Validate shapes + dtype
        Yu = self._as_mat("Yu", tp.Yu, N)
        Yd = self._as_mat("Yd", tp.Yd, N)
        Ye = self._as_mat("Ye", tp.Ye, N)
        Ynu = self._as_mat("Ynu", tp.Ynu, N)
        MR = self._as_mat("MR", tp.MR, N, allow_none=True)

        # Optional toggles / post-processing
        if bool(self._get_cfg("disable_quarks", False)):
            Yu = np.zeros_like(Yu)
            Yd = np.zeros_like(Yd)

        if bool(self._get_cfg("normalize", False)):
            Yu = self._opnorm_rescale(Yu)
            Yd = self._opnorm_rescale(Yd)
            Ye = self._opnorm_rescale(Ye)
            Ynu = self._opnorm_rescale(Ynu)

        scales = self._get_cfg("scales", {}) or {}
        Yu *= complex(scales.get("Yu", 1.0))
        Yd *= complex(scales.get("Yd", 1.0))
        Ye *= complex(scales.get("Ye", 1.0))
        Ynu *= complex(scales.get("Ynu", 1.0))
        if MR is not None:
            MR *= complex(scales.get("MR", 1.0))

        if MR is not None and bool(self._get_cfg("enforce_MR_symmetric", True)):
            MR = self._symmetrize_majorana(MR)

        return TexturePack(Yu=Yu, Yd=Yd, Ye=Ye, Ynu=Ynu, MR=MR)


# ============================================================
#  v4.0 Engine: orchestrates everything
# ============================================================

class AlignmentV40Engine:
    """
    One-stop v4.0 engine:
      - builds textures
      - builds D_int
      - builds even-product D_total
      - optionally builds fluctuations D_A
      - runs axioms + commutator gates
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

        # internal real structures
        self.J_SM = RealStructureFactory.build_J_SM(self.basis)
        self.J_flav = RealStructureFactory.build_J_flav(params.N_flav)
        self.J_int = RealStructureFactory.kron(self.J_SM, self.J_flav)

        # total J on H_geom ⊗ H_int
        self.J_total = RealStructureFactory.kron(self.geom.J_geom, self.J_int)

        # builders
        self.dint_builder = DIntBuilder(self.basis, self.maps)
        self.product_dirac = EvenProductDirac(self.geom, self.basis, params)

        # representation layer
        self.pi_SM = SMRepresentation(self.basis)

        # gates
        self.gates = V4Gates(
            eps0=params.eps_order0,
            eps1=params.eps_first,
            eps_sa=params.eps_sa,
            epsA=params.eps_fluct,   # use as hard oddness tolerance too
        )

    # -------- build core objects --------

    def build_textures(self) -> TexturePack:
        return self.texture_gen.build_textures(self.params)

    def build_D_int(self, textures: TexturePack) -> np.ndarray:
        return self.dint_builder.build_D_int(
            textures,
            J_SM=self.J_SM,
            N=self.params.N_flav,
        )

    def build_D(self, D_int: np.ndarray) -> np.ndarray:
        return self.product_dirac.build_D(D_int)

    def build_Gamma_total(self) -> np.ndarray:
        return self.product_dirac.build_Gamma_total(self.params.N_flav)

    # -------- fluctuations --------

    def build_fluctuated(
        self,
        D: np.ndarray,
        one_form_pairs: Sequence[Tuple[np.ndarray, np.ndarray]],
        *,
        weights: Sequence[complex] | None = None,
        hermitize: str | bool = "total",
        remove_trace: bool = False,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Returns (D_A, A).
        """
        fl = InnerFluctuations(self.J_total)
        A = fl.build_A(D, one_form_pairs, weights=weights, remove_trace=remove_trace)
        D_A = fl.fluctuate(D, A, hermitize=hermitize, stabilize_output=True)
        return D_A, A

    # -------- diagnostics / axioms --------

    def run_axiom_checks_internal(self, D_int: np.ndarray) -> Dict[str, float]:
        gamma_int = self.basis.lift_gamma(self.params.N_flav)
        sa = op_norm(D_int - D_int.conj().T)
        odd_int = op_norm(anti_comm(gamma_int, D_int))
        return {"sa": float(sa), "odd": float(odd_int), "sa_ok": sa <= self.params.eps_sa, "odd_ok": odd_int <= self.params.eps_fluct}

    def run_axiom_checks_total(self, D_total: np.ndarray) -> Dict[str, float]:
        Gamma = self.build_Gamma_total()
        sa = op_norm(D_total - D_total.conj().T)
        odd_total = op_norm(anti_comm(Gamma, D_total))
        return {"sa_total": float(sa), "odd_total": float(odd_total), "sa_ok": sa <= self.params.eps_sa, "odd_ok": odd_total <= self.params.eps_fluct}

    # -------- commutator gates --------

    def run_commutator_gates(
        self,
        D_total: np.ndarray,
        *,
        return_witness: bool = True,
        early_exit: bool = False,
        sample: Sequence[int] | None = None,
        odd_int_value: float | None = None,
    ) -> GateReport:
        pis = self._lift_pi_generators_to_full()

        if return_witness:
            order0, w0 = self.gates.gate_order0(pis, self.J_total, return_witness=True, early_exit=early_exit, sample=sample)
            first, w1 = self.gates.gate_first(D_total, pis, self.J_total, return_witness=True, early_exit=early_exit, sample=sample)
        else:
            order0 = self.gates.gate_order0(pis, self.J_total, early_exit=early_exit, sample=sample)
            first = self.gates.gate_first(D_total, pis, self.J_total, early_exit=early_exit, sample=sample)
            w0 = None
            w1 = None

        sa = self.gates.gate_selfadj(D_total)

        # oddness is fundamentally internal; accept caller value or compute from D_int elsewhere.
        odd_int = float(odd_int_value) if odd_int_value is not None else 0.0

        return GateReport(
            order0=float(order0),
            first=float(first),
            self_adj=float(sa),
            odd_int=float(odd_int),
            order0_pair=w0 if return_witness else None,
            first_pair=w1 if return_witness else None,
        )

    def _lift_pi_generators_to_full(self) -> List[np.ndarray]:
        """
        Lift π(A) generators to full H_geom ⊗ H_SM ⊗ H_flav
        (π acts trivially on flavor and geometry).
        """
        N = self.params.N_flav
        gens_SM = self.pi_SM.generators()  # 32x32
        gens_int = [np.kron(g, np.eye(N, dtype=complex)) for g in gens_SM]  # 32N x 32N

        I_geom = np.eye(self.geom.dim(), dtype=complex)
        return [np.kron(I_geom, g_int) for g_int in gens_int]

    # -------- sector extraction utilities --------

    def sector_projector_int(self, key: str) -> np.ndarray:
        """(32N x 32N) projector onto an internal sector key."""
        P32 = self.basis.P[key]
        return self.basis.lift_projector(P32, self.params.N_flav)

    def sector_projector_full(self, key: str) -> np.ndarray:
        """(dim_geom*32N x dim_geom*32N) projector onto a sector in the full product space."""
        I_geom = np.eye(self.geom.dim(), dtype=complex)
        return np.kron(I_geom, self.sector_projector_int(key))

    def extract_lr_block_int(self, D_int: np.ndarray, left_key: str, right_key: str) -> np.ndarray:
        """P_L D_int P_R on H_SM ⊗ H_flav."""
        PL = self.sector_projector_int(left_key)
        PR = self.sector_projector_int(right_key)
        return PL @ D_int @ PR

    def extract_lr_block_full(self, D_total: np.ndarray, left_key: str, right_key: str) -> np.ndarray:
        """(I⊗P_L) D_total (I⊗P_R) on H_geom ⊗ H_int."""
        PL = self.sector_projector_full(left_key)
        PR = self.sector_projector_full(right_key)
        return PL @ D_total @ PR

    # -------- end-to-end run --------

    def run(
        self,
        *,
        do_fluctuate: bool = False,
        one_form_pairs: Sequence[Tuple[np.ndarray, np.ndarray]] | None = None,
        fluct_weights: Sequence[complex] | None = None,
        hermitize: str | bool = "total",
        gate_sample: Sequence[int] | None = None,
        early_exit_gates: bool = False,
    ) -> Dict[str, Any]:
        textures = self.build_textures()
        D_int = self.build_D_int(textures)
        D_total = self.build_D(D_int)

        ax_int = self.run_axiom_checks_internal(D_int)
        ax_tot = self.run_axiom_checks_total(D_total)

        D_used = D_total
        A = None
        if do_fluctuate:
            if one_form_pairs is None:
                raise ValueError("do_fluctuate=True requires one_form_pairs.")
            D_used, A = self.build_fluctuated(
                D_total,
                one_form_pairs,
                weights=fluct_weights,
                hermitize=hermitize,
            )

        gates = self.run_commutator_gates(
            D_used,
            return_witness=True,
            early_exit=early_exit_gates,
            sample=gate_sample,
            odd_int_value=ax_int["odd"],
        )

        return {
            "textures": textures,
            "D_int": D_int,
            "D_total": D_total,
            "A": A,
            "D_used": D_used,
            "axioms_internal": ax_int,
            "axioms_total": ax_tot,
            "gates": gates,
            "pass": self.gates.pass_fail(gates),
        }


# ============================================================
#  Optional: Spectral action evaluators (conceptual / numerical)
# ============================================================

class SpectralAction:
    """
    Bosonic spectral action approximations on finite truncations.

    Supported evaluation modes:
      1) exact diagonalization (dense):     Tr f(D^2/Λ^2)
      2) partial eigenvalues (k largest):   good if f decays fast and you only need top modes
      3) stochastic trace (Hutchinson):    Tr f(D^2/Λ^2) ≈ (1/M) Σ z^T f(D^2/Λ^2) z
         using Chebyshev polynomial approximation for f on [0, R].

    Notes:
      - For stability we always hermitize D first.
      - We prefer working with D^2 as a PSD operator when possible.
    """

    def __init__(
        self,
        cutoff_Lambda: float,
        f: Callable[[np.ndarray], np.ndarray],
        *,
        atol_sa: float = 1e-10,
    ) -> None:
        self.Lambda = float(cutoff_Lambda)
        self.f = f
        self.atol_sa = float(atol_sa)

    # ---------------- basics ----------------

    @staticmethod
    def _require_square(name: str, X: np.ndarray) -> None:
        if X.ndim != 2 or X.shape[0] != X.shape[1]:
            raise ValueError(f"{name} must be square, got {X.shape}")

    def _herm(self, D: np.ndarray) -> np.ndarray:
        D = np.asarray(D, dtype=complex)
        self._require_square("D", D)
        return hermitian_part(D)

    # ---------------- exact trace ----------------

    def trace_f(self, D: np.ndarray) -> float:
        """
        Exact: Tr f(D^2/Λ^2) by full diagonalization of hermitian(D).
        """
        Dh = self._herm(D)
        evals = np.linalg.eigvalsh(Dh)  # real
        x = (evals.astype(float) ** 2) / (self.Lambda ** 2)
        fx = np.asarray(self.f(x))
        return float(np.sum(fx))

    def trace_f_from_eigs(self, eigs: np.ndarray) -> float:
        """
        If you already have eigenvalues of hermitian(D), reuse them.
        """
        eigs = np.asarray(eigs, dtype=float).reshape(-1)
        x = (eigs ** 2) / (self.Lambda ** 2)
        return float(np.sum(self.f(x)))

    # ---------------- partial spectrum (optional) ----------------

    def trace_f_topk(self, D: np.ndarray, k: int) -> float:
        """
        Approximate Tr f(D^2/Λ^2) using only the k largest-magnitude eigenvalues of D.
        Useful only if f(x) decays quickly and truncation is justified.
        Dense fallback implementation (still O(n^3) worst-case), but isolates the API.

        If you want true sparse partial eigensolve, plug scipy.sparse.linalg.eigsh in here.
        """
        Dh = self._herm(D)
        n = Dh.shape[0]
        k = int(k)
        if k <= 0 or k > n:
            raise ValueError(f"k must be in 1..{n}, got {k}")

        # Dense "partial" via full eigvals then select; replace with sparse eigsh in production.
        evals = np.linalg.eigvalsh(Dh)
        # select k largest |λ|
        idx = np.argsort(np.abs(evals))[::-1][:k]
        sel = evals[idx]
        x = (sel.astype(float) ** 2) / (self.Lambda ** 2)
        return float(np.sum(self.f(x)))

    # ---------------- fermionic term ----------------

    def fermionic_term(self, psi: np.ndarray, D: np.ndarray, *, normalize: bool = False) -> complex:
        """
        <psi, D psi>. If normalize=True, uses psi/||psi||.
        """
        Dh = self._herm(D)
        psi = np.asarray(psi, dtype=complex).reshape(-1)
        if psi.shape[0] != Dh.shape[0]:
            raise ValueError(f"psi length {psi.shape[0]} does not match D dim {Dh.shape[0]}")
        if normalize:
            nrm = np.linalg.norm(psi)
            if nrm > 0:
                psi = psi / nrm
        return complex(np.vdot(psi, Dh @ psi))

    # ---------------- stochastic trace (Hutchinson) ----------------

    def trace_f_hutchinson(
        self,
        D: np.ndarray,
        *,
        num_samples: int = 64,
        rng_seed: int = 0,
    ) -> float:
        """
        Hutchinson estimator for Tr f(D^2/Λ^2) using *exact* f(D^2/Λ^2) application
        via diagonalization of D (so this is still expensive but shows the interface).

        In production, you’d replace the body with:
          - build a LinearOperator for A = D^2/Λ^2
          - approximate f(A)z via Chebyshev / Lanczos
        """
        Dh = self._herm(D)
        evals, U = np.linalg.eigh(Dh)  # Dh = U diag(evals) U†
        x = (evals.astype(float) ** 2) / (self.Lambda ** 2)
        fx = np.asarray(self.f(x), dtype=float)

        rng = np.random.default_rng(rng_seed)
        n = Dh.shape[0]
        M = int(num_samples)
        if M <= 0:
            raise ValueError("num_samples must be > 0")

        # z are Rademacher vectors for Hutchinson
        acc = 0.0
        for _ in range(M):
            z = rng.choice([-1.0, 1.0], size=(n,))
            # compute z^T f(A) z where f(A)=U diag(fx) U†
            Uz = U.conj().T @ z
            acc += float(np.sum(fx * (np.abs(Uz) ** 2)))
        return float(acc / M)

# ============================================================
#  Minimal self-tests (fast numerical sanity)
# ============================================================

def _self_test_sm_structures() -> None:
    b = SMBasis32()

    # gamma_SM sanity
    assert is_hermitian(b.gamma_SM)
    assert np.allclose(b.gamma_SM @ b.gamma_SM, np.eye(32), atol=1e-12)

    # J_SM swap sanity
    J = RealStructureFactory.build_J_SM(b)
    U = J.U
    assert np.allclose(U.conj().T @ U, np.eye(32), atol=1e-12)
    # swap property on basis vectors (up to exact equality in this basis)
    e0 = np.zeros((32,), dtype=complex); e0[0] = 1.0
    e16 = np.zeros((32,), dtype=complex); e16[16] = 1.0
    assert np.allclose(U @ e0, e16)
    assert np.allclose(U @ e16, e0)

    # ChannelMaps partial-isometry diagnostics
    maps = ChannelMaps(b)
    diags = maps.validate()
    # they should be exactly 0 in floating arithmetic for these integer matrices
    for k, v in diags.items():
        assert v < 1e-12, f"{k} deviation too large: {v}"

    # Representation generator sanity
    rep = SMRepresentation(b)
    gens = rep.generators()
    assert all(g.shape == (32, 32) for g in gens)

def _self_test_dirac_builders() -> None:
    # tiny sizes for speed
    N = 2

    # SM core
    b = SMBasis32()
    maps = ChannelMaps(b)
    J_SM = RealStructureFactory.build_J_SM(b)

    # textures: simple deterministic (nontrivial but stable)
    Yu = np.zeros((N, N), dtype=complex)
    Yd = np.zeros((N, N), dtype=complex)
    Ye = np.eye(N, dtype=complex)
    Ynu = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)  # swap
    MR = 1e9 * np.eye(N, dtype=complex)
    textures = TexturePack(Yu=Yu, Yd=Yd, Ye=Ye, Ynu=Ynu, MR=MR)

    # Build D_int
    dint = DIntBuilder(b, maps).build_D_int(textures, J_SM=J_SM, N=N)
    assert dint.shape == (32 * N, 32 * N)
    assert is_hermitian(dint)

    # Internal oddness: {gamma_SM⊗I, D_int} = 0
    gamma_int = b.lift_gamma(N)
    odd_int = op_norm(anti_comm(gamma_int, dint))
    assert odd_int < 1e-10, f"odd_int too large: {odd_int}"

    # Geometry: use even truncation so {gamma_geom, D_geom}=0 by construction
    geom = GeometryFactor.even_truncation(lambdas=np.array([1.0, 2.0]))

    # Product Dirac
    class _P:  # minimal params shim
        eps_sa = 1e-10
        eps_fluct = 1e-10
    prod = EvenProductDirac(geom, b, _P())

    Dtot = prod.build_D(dint)
    assert Dtot.shape == (geom.dim() * 32 * N, geom.dim() * 32 * N)
    assert is_hermitian(Dtot)

    # Total oddness: {Gamma_total, D_total} = 0
    Gamma = prod.build_Gamma_total(N)
    odd_tot = op_norm(anti_comm(Gamma, Dtot))
    assert odd_tot < 1e-10, f"odd_total too large: {odd_tot}"

def _self_test_engine_smoke() -> None:
    # small even geometry
    geom = GeometryFactor.even_truncation(lambdas=np.array([1.0, 3.0]))

    # minimal params object (match the fields your engine uses)
    class _Params:
        N_flav = 2
        seesaw_M_GeV = 1e9

        eps_order0 = 1e-10
        eps_first = 1e-10
        eps_sa = 1e-10
        eps_fluct = 1e-10

    params = _Params()

    # textures: deterministic identity leptons (fast)
    texgen = TextureGenerator(cfg={"mode": "identity_leptons"})

    eng = AlignmentV40Engine(geom=geom, params=params, texture_gen=texgen)

    # build core
    textures = eng.build_textures()
    D_int = eng.build_D_int(textures)
    D_tot = eng.build_D(D_int)

    # internal axioms
    ax_int = eng.run_axiom_checks_internal(D_int)
    assert ax_int["sa_ok"] and ax_int["odd_ok"], f"internal axioms failed: {ax_int}"

    # total axioms
    ax_tot = eng.run_axiom_checks_total(D_tot)
    assert ax_tot["sa_ok"] and ax_tot["odd_ok"], f"total axioms failed: {ax_tot}"

    # commutator gates (sample a few generators for speed)
    rep = eng.run_commutator_gates(
        D_tot,
        return_witness=True,
        early_exit=True,
        sample=[0, 1, 2, 3],          # I, Y, first SU2 gens...
        odd_int_value=ax_int["odd"],
    )
    assert np.isfinite(rep.order0) and np.isfinite(rep.first) and np.isfinite(rep.self_adj)

    # optional: do one tiny fluctuation using a couple of lifted generators
    pis = eng._lift_pi_generators_to_full()
    pairs = [(pis[1], pis[2]), (pis[2], pis[1])]  # simple symmetric-ish pair set
    D_A, A = eng.build_fluctuated(D_tot, pairs, hermitize="total")
    assert D_A.shape == D_tot.shape
    assert is_hermitian(D_A)

def run_self_tests() -> None:
    _self_test_sm_structures()
    _self_test_dirac_builders()
    _self_test_engine_smoke()

# ============================================================
#  Example wiring (scaffold)
# ============================================================

def example_build_engine() -> AlignmentV40Engine:
    # Params: assumes your V40Params has these defaults; override eps if desired.
    params = V40Params(N_flav=3)

    # Use an even truncation so {D_geom, gamma_geom}=0 holds.
    # lambdas are the geometric mode eigenvalues in the doubled construction.
    geom = GeometryFactor.even_truncation(lambdas=np.array([1, 2, 3, 4, 5], dtype=float))

    # Default stub textures (identity leptons) unless cfg says otherwise.
    texgen = TextureGenerator(cfg={"mode": "identity_leptons"})

    return AlignmentV40Engine(geom=geom, params=params, texture_gen=texgen)


def main() -> None:
    eng = example_build_engine()

    out = eng.run(
        do_fluctuate=False,
        gate_sample=[0, 1, 2, 3],     # quick gate smoke run (optional)
        early_exit_gates=True,
    )

    print("Axioms (internal):", out["axioms_internal"])
    print("Axioms (total):   ", out["axioms_total"])
    print("Gate report:", out["gates"])
    print("PASS:", out["pass"])


if __name__ == "__main__":
    main()

"""
/Users/chazzromeo/rAI/Flavor/.venv/bin/python /Users/chazzromeo/rAI/Flavor/alignment_v40_engine_step1.py 
Axioms (internal): {'sa': 0.0, 'odd': 0.0, 'sa_ok': True, 'odd_ok': True}
Axioms (total):    {'sa_total': 0.0, 'odd_total': 0.0, 'sa_ok': True, 'odd_ok': True}
Gate report: GateReport(order0=0.0, first=0.0, self_adj=0.0, odd_int=0.0, order0_pair=(0, 0), first_pair=(0, 0))
PASS: True
"""