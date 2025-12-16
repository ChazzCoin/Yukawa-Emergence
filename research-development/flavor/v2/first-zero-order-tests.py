#!/usr/bin/env python3
"""
ncg_tests.py

Reusable Noncommutative Geometry test harness for finite spectral triples.

Provides:
  - AlgebraElement dataclass
  - First-order condition checker
  - Zero-order condition checker
  - Grading & reality checker
  - A convenience function to run the full test suite and pretty-print results
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional
import numpy as np

# ============================================================
# 1. Basis for 1-generation internal Hilbert space
# ============================================================

@dataclass
class BasisState:
    name: str          # e.g. "nu_L", "e_L", "u_L_r", "d_R_b", ...
    chirality: str     # "L" or "R"
    particle: bool     # True = particle, False = antiparticle
    is_quark: bool     # True for quarks, False for leptons
    color: Optional[str]  # "r","g","b" or None
    generation: int    # 1 (for now)


def build_sm_basis_1gen(include_nu_R: bool = True) -> Tuple[List[BasisState], Dict[str, int]]:
    """
    One generation of SM fermions + antiparticles.

    Particle sector:
      Leptons:
        L_L: (nu_L, e_L)
        R :  (nu_R (optional), e_R)

      Quarks (color = r,g,b):
        Q_L: (u_L^r, d_L^r, u_L^g, d_L^g, u_L^b, d_L^b)
        R :  (u_R^r, u_R^g, u_R^b, d_R^r, d_R^g, d_R^b)

    Antiparticle sector: charge conjugates of all above, with same chirality label
    (J will swap particle ↔ antiparticle sectors).
    """
    basis: List[BasisState] = []

    def add(name, chirality, particle, is_quark, color=None):
        basis.append(BasisState(
            name=name, chirality=chirality, particle=particle,
            is_quark=is_quark, color=color, generation=1
        ))

    # --- Particle Leptons ---
    add("nu_L", "L", True, False)
    add("e_L",  "L", True, False)
    if include_nu_R:
        add("nu_R", "R", True, False)
    add("e_R",  "R", True, False)

    # --- Particle Quarks (3 colors) ---
    colors = ["r", "g", "b"]
    for col in colors:
        add(f"u_L_{col}", "L", True, True, color=col)
        add(f"d_L_{col}", "L", True, True, color=col)
    for col in colors:
        add(f"u_R_{col}", "R", True, True, color=col)
    for col in colors:
        add(f"d_R_{col}", "R", True, True, color=col)

    # At this point, count particle states:
    # leptons: 2L + (1 or 2)R = 3 or 4
    # quarks   6L + 3R + 3R = 12
    # total particle states = 15 or 16
    n_particle = len(basis)

    # --- Antiparticles: one conjugate state per particle state ---
    for bs in list(basis):
        add(bs.name + "_c", bs.chirality, False, bs.is_quark, bs.color)

    idx: Dict[str, int] = {bs.name: i for i, bs in enumerate(basis)}
    return basis, idx


# ============================================================
# 2. Gamma_F and J_F (swap matrix)
# ============================================================

def build_gamma_F_SM(basis: List[BasisState]) -> np.ndarray:
    """
    gamma_F = -1 on left-handed states, +1 on right-handed states,
    for both particles and antiparticles.
    """
    dimH = len(basis)
    gamma = np.zeros((dimH, dimH), dtype=complex)
    for i, bs in enumerate(basis):
        sgn = -1.0 if bs.chirality == "L" else +1.0
        gamma[i, i] = sgn
    return gamma


def build_swap_particle_antiparticle(basis: List[BasisState], idx: Dict[str, int]) -> np.ndarray:
    """
    Build S such that:
      S |particle> = |particle_c>
      S |particle_c> = |particle>
    i.e. S^2 = I.
    """
    dimH = len(basis)
    S = np.zeros((dimH, dimH), dtype=complex)

    for bs in basis:
        if bs.particle:
            i = idx[bs.name]
            j = idx[bs.name + "_c"]
            S[i, j] = 1.0
            S[j, i] = 1.0

    return S


def J_action_SM(M: np.ndarray, S: np.ndarray, phase: complex = 1.0) -> np.ndarray:
    """
    Real structure on matrices:
        J M J^{-1} = phase * S M^* S^T
    """
    return phase * (S @ M.conj() @ S.T)


# ============================================================
# 3. Representation of A_F = C ⊕ H ⊕ M_3(C)
# ============================================================

@dataclass
class SMAlgebraElement:
    label: str
    op: np.ndarray


def quaternion_to_2x2(a: complex, b: complex) -> np.ndarray:
    """
    Represent q = a + b j as 2x2 complex matrix:
      [ a   b ]
      [-b* a* ]
    For our purposes, we only need a toy faithful representation.
    """
    a = complex(a)
    b = complex(b)
    return np.array([[a, b], [-b.conjugate(), a.conjugate()]], dtype=complex)


def rep_A_SM(
    lam: complex,
    q: np.ndarray,
    m3: np.ndarray,
    basis: List[BasisState],
    idx: Dict[str, int],
) -> np.ndarray:
    """
    Represent (lam ∈ C, q ∈ H≈2x2, m3 ∈ M3(C)) on H_F.

    Simplified rules:
      - lambda (C-part) acts as scalar on all states.
      - q (H-part) acts non-trivially on SU(2)_L doublets:
          (nu_L, e_L) and (u_L^c, d_L^c) for each color,
        and acts trivially on SU(2) singlets (all R states).
      - m3 acts on color indices of quarks (3-dim rep), trivially on leptons.

    Antiparticle sector: use same representation (J takes care of conjugation).
    """
    dimH = len(basis)
    A = np.zeros((dimH, dimH), dtype=complex)

    # C-part: global scalar
    A += lam * np.eye(dimH, dtype=complex)

    # H-part: SU(2)_L doublets
    # (nu_L, e_L)
    if "nu_L" in idx and "e_L" in idx:
        i_nu = idx["nu_L"]
        i_e  = idx["e_L"]
        # insert q on the (nu_L, e_L) subspace
        A[np.ix_([i_nu, i_e], [i_nu, i_e])] += q

    # Quark doublets Q_L: (u_L_col, d_L_col) for each color
    colors = ["r", "g", "b"]
    for col in colors:
        u_name = f"u_L_{col}"
        d_name = f"d_L_{col}"
        if u_name in idx and d_name in idx:
            i_u = idx[u_name]
            i_d = idx[d_name]
            A[np.ix_([i_u, i_d], [i_u, i_d])] += q

    # H-part acts trivially on R states and we also keep it trivial on antiparticles
    # for this first-pass implementation (can be refined later).

    # M3-part: color action on quarks
    # For each chirality, quark multiplet (u, d) share the same color rep.
    # We treat leptons as color singlets (no action).
    for bs in basis:
        if bs.is_quark and bs.color is not None:
            i = idx[bs.name]
            # Build an ordering of colors: r,g,b → 0,1,2
            col_index = {"r": 0, "g": 1, "b": 2}[bs.color]

            # For simplicity, we let m3 act as diag(m3[col,col]) in this basis;
            # a more faithful representation would mix colors explicitly.
            A[i, i] += m3[col_index, col_index]

    return A


def build_SM_algebra_generators(
    basis: List[BasisState],
    idx: Dict[str, int],
) -> List[SMAlgebraElement]:
    """
    Build a small generating set of A_F elements:
      - I (identity)
      - a couple of quaternion directions
      - a couple of color generators
    """
    dimH = len(basis)
    I = np.eye(dimH, dtype=complex)

    # Basic quaternion directions
    q_id = quaternion_to_2x2(1.0, 0.0)
    q_j  = quaternion_to_2x2(0.0, 1.0)

    # Basic color matrices: identity and λ3-like diagonal
    m3_id = np.eye(3, dtype=complex)
    m3_diag = np.diag([1.0, -1.0, 0.0])

    ops: List[SMAlgebraElement] = []

    # Identity in A_F: (lam=1, q=I_2, m3=I_3)
    A_I = rep_A_SM(1.0 + 0j, q_id, m3_id, basis, idx)
    ops.append(SMAlgebraElement("I", A_I))

    # Pure quaternion j on SU(2)_L
    A_qj = rep_A_SM(0.0 + 0j, q_j, m3_id, basis, idx)
    ops.append(SMAlgebraElement("H_j", A_qj))

    # Pure color diagonal
    A_color_diag = rep_A_SM(0.0 + 0j, q_id, m3_diag, basis, idx)
    ops.append(SMAlgebraElement("color_diag", A_color_diag))

    return ops


# ============================================================
# 4. Dirac operator D_F for 1 generation
# ============================================================

def build_DF_SM_1gen(
    Y_e: complex,
    Y_nu: complex,
    Y_u: complex,
    Y_d: complex,
    basis: List[BasisState],
    idx: Dict[str, int],
    include_nu_R: bool = True,
) -> np.ndarray:
    """
    Very minimal 1-generation Dirac operator:

    - Acts only between particle L and R in each sector using Yukawa couplings:
        D_F |e_L>  ~ Y_e |e_R>
        D_F |nu_L> ~ Y_nu|nu_R> (if present)
        D_F |u_L_c> ~ Y_u |u_R_c>
        D_F |d_L_c> ~ Y_d |d_R_c>
    - Antiparticle block mirrors the same structure.

    This is just enough structure to let you:
      - plug in canonical SM Yukawas,
      - plug in aligned Yukawas from your pipeline,
      - and run order tests against the SM-like algebra above.
    """
    dimH = len(basis)
    D = np.zeros((dimH, dimH), dtype=complex)

    def couple(L_name: str, R_name: str, Y: complex):
        if L_name in idx and R_name in idx:
            iL = idx[L_name]
            iR = idx[R_name]
            D[iL, iR] = Y.conjugate()
            D[iR, iL] = Y

    # --- Particle sector couplings ---
    # Leptons
    if include_nu_R and "nu_L" in idx and "nu_R" in idx:
        couple("nu_L", "nu_R", Y_nu)
    couple("e_L", "e_R", Y_e)

    # Quarks (3 colors)
    for col in ["r", "g", "b"]:
        couple(f"u_L_{col}", f"u_R_{col}", Y_u)
        couple(f"d_L_{col}", f"d_R_{col}", Y_d)

    # --- Antiparticle sector couplings ---
    # mirror the same pattern for the conjugate states
    def conj_name(name: str) -> str:
        return name + "_c"

    if include_nu_R and "nu_L_c" in idx and "nu_R_c" in idx:
        couple("nu_L_c", "nu_R_c", Y_nu)
    couple("e_L_c", "e_R_c", Y_e)

    for col in ["r", "g", "b"]:
        couple(f"u_L_{col}_c", f"u_R_{col}_c", Y_u)
        couple(f"d_L_{col}_c", f"d_R_{col}_c", Y_d)

    return D


def main():
    basis, idx = build_sm_basis_1gen(include_nu_R=True)

    # SM-like Yukawas (1 generation, rough magnitudes)
    Y_e  = 2.94e-6     # me / v
    Y_nu = 1.0e-12     # tiny
    Y_u  = 1.3e-5      # up-type
    Y_d  = 2.8e-5      # down-type

    D_F_SM = build_DF_SM_1gen(Y_e, Y_nu, Y_u, Y_d, basis, idx)

    # Build SM algebra generators
    sm_ops = build_SM_algebra_generators(basis, idx)

    # If your ncg_tests expects AlgebraElement, you can just wrap:
    algebra = [AlgebraElement(op.label, op.op) for op in sm_ops]

    # Run tests
    run_ncg_test_suite(D_F_SM, algebra, name="SM-like 1gen finite triple")


# ===========================================
# Basic helpers
# ===========================================

def assert_square(M: np.ndarray, name: str = "matrix") -> None:
    if M.shape[0] != M.shape[1]:
        raise ValueError(f"{name} must be square, got {M.shape}.")


def assert_even_dim(M: np.ndarray, name: str = "matrix") -> None:
    assert_square(M, name)
    n = M.shape[0]
    if n % 2 != 0:
        raise ValueError(f"{name} must have even dimension, got {n}.")


def frob_norm(M: np.ndarray) -> float:
    return float(np.linalg.norm(M, ord="fro"))


# ===========================================
# Data structures
# ===========================================

@dataclass
class AlgebraElement:
    label: str
    op: np.ndarray


@dataclass
class FirstOrderResult:
    max_norm: float
    worst_pair: Optional[Tuple[str, str]]
    good_pairs: List[Tuple[str, str, float]]


@dataclass
class ZeroOrderResult:
    max_norm: float
    worst_pair: Optional[Tuple[str, str]]
    bad_pairs: List[Tuple[str, str, float]]


@dataclass
class GradingRealityResult:
    anticom_norm: float          # || {gamma_F, D_F} ||_F
    max_comm_gamma: float        # max ||[gamma_F, a]||_F over a in A
    J2_deviation: float          # ||J^2 - I||_F
    norm_plus: float             # ||J D_F J^-1 - D_F||_F
    norm_minus: float            # ||J D_F J^-1 + D_F||_F
    ko_sign: Optional[int]       # +1, -1, or None


# ===========================================
# Real structure & grading builders
# ===========================================

def build_swap_LR(dim_left: int) -> np.ndarray:
    """
    Swap matrix S on H = H_L ⊕ H_R, dim(H_L) = dim(H_R) = dim_left.
    """
    S = np.zeros((2 * dim_left, 2 * dim_left), dtype=complex)
    S[:dim_left, dim_left:] = np.eye(dim_left)
    S[dim_left:, :dim_left] = np.eye(dim_left)
    return S


def build_gamma_F(dim_left: int) -> np.ndarray:
    """
    Grading operator gamma_F with eigenvalue -1 on H_L and +1 on H_R.
    """
    g = np.zeros((2 * dim_left, 2 * dim_left), dtype=complex)
    g[:dim_left, :dim_left] = -np.eye(dim_left)
    g[dim_left:, dim_left:] = +np.eye(dim_left)
    return g


def J_action(M: np.ndarray, S: np.ndarray, phase: complex = 1.0) -> np.ndarray:
    """
    Implement J M J^{-1} = phase * S M^* S^T, where S encodes the L/R swap
    (and potentially more structure in the future).
    """
    return phase * (S @ M.conj() @ S.T)


# ===========================================
# First-order condition
# ===========================================

def check_first_order_condition(
    D_F: np.ndarray,
    algebra: List[AlgebraElement],
    eps: float = 1e-12,
    J_phase: complex = 1.0,
) -> FirstOrderResult:
    """
    First-order condition:
        [[D_F, a], J b J^{-1}] = 0
    for all a,b in the algebra.

    Returns FirstOrderResult with max norm and list of "good pairs"
    whose violation is below eps.
    """
    assert_even_dim(D_F, "D_F")
    n = D_F.shape[0]
    dim_left = n // 2

    S = build_swap_LR(dim_left)

    max_norm = 0.0
    worst_pair = None
    good_pairs: List[Tuple[str, str, float]] = []

    for a in algebra:
        Da = D_F @ a.op - a.op @ D_F
        for b in algebra:
            b_tilde = J_action(b.op, S, phase=J_phase)
            comm2 = Da @ b_tilde - b_tilde @ Da
            norm = frob_norm(comm2)

            if norm > max_norm:
                max_norm = norm
                worst_pair = (a.label, b.label)

            if norm < eps:
                good_pairs.append((a.label, b.label, norm))

    return FirstOrderResult(
        max_norm=max_norm,
        worst_pair=worst_pair,
        good_pairs=good_pairs,
    )


def print_first_order_result(
    result: FirstOrderResult,
    eps: float = 1e-12,
    title: str = "First-order condition",
) -> None:
    print(f"=== {title} ===")
    print(f"Max Frobenius norm over all pairs (a,b): {result.max_norm:.3e}")
    if result.worst_pair is not None:
        la, lb = result.worst_pair
        print(f"Worst offender: (a={la}, b={lb})")
    if result.good_pairs:
        print(f"Pairs with norm < {eps:.1e}:")
        for la, lb, nrm in result.good_pairs:
            print(f"  (a={la:>12s}, b={lb:>12s}) → ||[[D,a],J b J^-1]||_F = {nrm:.3e}")
    else:
        print(f"No pairs with norm < {eps:.1e}")
    print()


# ===========================================
# Zero-order condition
# ===========================================

def check_zero_order_condition(
    algebra: List[AlgebraElement],
    eps: float = 1e-12,
    J_phase: complex = 1.0,
) -> ZeroOrderResult:
    """
    Zero-order condition:
        [a, J b J^{-1}] = 0
    for all a,b in the algebra.
    """
    if not algebra:
        raise ValueError("Algebra must contain at least one element.")

    n = algebra[0].op.shape[0]
    assert_even_dim(algebra[0].op, "algebra[0].op")
    dim_left = n // 2

    S = build_swap_LR(dim_left)

    max_norm = 0.0
    worst_pair = None
    bad_pairs: List[Tuple[str, str, float]] = []

    for a in algebra:
        for b in algebra:
            b_tilde = J_action(b.op, S, phase=J_phase)
            comm = a.op @ b_tilde - b_tilde @ a.op
            norm = frob_norm(comm)

            if norm > max_norm:
                max_norm = norm
                worst_pair = (a.label, b.label)

            if norm > eps:
                bad_pairs.append((a.label, b.label, norm))

    return ZeroOrderResult(
        max_norm=max_norm,
        worst_pair=worst_pair,
        bad_pairs=bad_pairs,
    )


def print_zero_order_result(
    result: ZeroOrderResult,
    eps: float = 1e-12,
    title: str = "Zero-order condition",
) -> None:
    print(f"=== {title} ===")
    print(f"Max Frobenius norm over all pairs (a,b): {result.max_norm:.3e}")
    if result.worst_pair is not None:
        la, lb = result.worst_pair
        print(f"Worst offender: (a={la}, b={lb})")
    if result.bad_pairs:
        print("Pairs with significant violation:")
        for la, lb, nrm in result.bad_pairs:
            print(f"  (a={la:>12s}, b={lb:>12s}) → ||[a, J b J^-1]||_F = {nrm:.3e}")
    else:
        print(f"All pairs satisfy [a, J b J^-1] ≈ 0 within eps={eps:.1e}")
    print()


# ===========================================
# Grading & reality
# ===========================================

def check_grading_and_reality(
    D_F: np.ndarray,
    algebra: List[AlgebraElement],
    J_phase: complex = 1.0,
) -> GradingRealityResult:
    """
    - Check gamma_F anti-commutes with D_F and commutes with A_F.
    - Check J^2 = 1 (via S^2).
    - Estimate KO-dimension sign via J D_F J^{-1} = ± D_F.
    """
    assert_even_dim(D_F, "D_F")
    n = D_F.shape[0]
    dim_left = n // 2

    gamma_F = build_gamma_F(dim_left)
    S = build_swap_LR(dim_left)

    # {gamma_F, D_F}
    anti = gamma_F @ D_F + D_F @ gamma_F
    anticom_norm = frob_norm(anti)

    # [gamma_F, a]
    max_comm_gamma = 0.0
    for a in algebra:
        comm = gamma_F @ a.op - a.op @ gamma_F
        max_comm_gamma = max(max_comm_gamma, frob_norm(comm))

    # J^2 - I
    S2 = S @ S
    J2_deviation = frob_norm(S2 - np.eye(n))

    # KO-sign
    JDJ = J_action(D_F, S, phase=J_phase)
    norm_plus = frob_norm(JDJ - D_F)
    norm_minus = frob_norm(JDJ + D_F)

    scale = frob_norm(D_F) + 1e-16
    rel_plus = norm_plus / scale
    rel_minus = norm_minus / scale

    ko_sign: Optional[int]
    if rel_plus < 1e-12 and rel_minus > rel_plus:
        ko_sign = +1
    elif rel_minus < 1e-12 and rel_plus > rel_minus:
        ko_sign = -1
    else:
        ko_sign = None

    return GradingRealityResult(
        anticom_norm=anticom_norm,
        max_comm_gamma=max_comm_gamma,
        J2_deviation=J2_deviation,
        norm_plus=norm_plus,
        norm_minus=norm_minus,
        ko_sign=ko_sign,
    )


def print_grading_and_reality_result(
    result: GradingRealityResult,
    title: str = "Grading & reality tests",
) -> None:
    print(f"=== {title} ===")
    print(f"||{{gamma_F, D_F}}||_F       = {result.anticom_norm:.3e}")
    print(f"max ||[gamma_F, a]||_F       = {result.max_comm_gamma:.3e}")
    print(f"||J^2 - I||_F                = {result.J2_deviation:.3e}")
    print(f"||J D_F J^-1 - D_F||_F       = {result.norm_plus:.3e}")
    print(f"||J D_F J^-1 + D_F||_F       = {result.norm_minus:.3e}")
    if result.ko_sign == +1:
        print("→ KO-sign: J D_F J^-1 = + D_F (J-even Dirac operator)")
    elif result.ko_sign == -1:
        print("→ KO-sign: J D_F J^-1 = - D_F (J-odd Dirac operator)")
    else:
        print("→ KO-sign: ambiguous or not clean at numerical precision.")
    print()


# ===========================================
# Master convenience function
# ===========================================

def run_ncg_test_suite(
    D_F: np.ndarray,
    algebra: List[AlgebraElement],
    eps_first: float = 1e-12,
    eps_zero: float = 1e-12,
    J_phase: complex = 1.0,
    name: str = "",
) -> Tuple[FirstOrderResult, ZeroOrderResult, GradingRealityResult]:
    """
    Run first-order, zero-order, and grading/reality tests and pretty-print a summary.

    Returns (FirstOrderResult, ZeroOrderResult, GradingRealityResult).
    """
    if name:
        print(f"=== NCG test suite for {name} ===")

    fo_res = check_first_order_condition(D_F, algebra, eps=eps_first, J_phase=J_phase)
    print_first_order_result(fo_res, eps=eps_first)

    zo_res = check_zero_order_condition(algebra, eps=eps_zero, J_phase=J_phase)
    print_zero_order_result(zo_res, eps=eps_zero)

    gr_res = check_grading_and_reality(D_F, algebra, J_phase=J_phase)
    print_grading_and_reality_result(gr_res)

    return fo_res, zo_res, gr_res

if __name__ == "__main__":
    main()