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
from scipy.linalg import kron

from v3.spectral import SMState, build_internal_algebra_ops

# =========================
# v4.0 EVEN PRODUCT TOGGLES
# =========================
USE_V4_EVEN_PRODUCT = True   # set False to keep your current "odd product" behavior
EPS_GAMMA = 1e-12            # tolerance for grading anti/commutation sanity
TOPK_VIOLATIONS = 10         # show worst offenders for debug


def hermitian_basis(n: int):
    """
    Full Hermitian basis {E_k} for n x n Hermitian matrices.
    """
    basis = []
    # diagonal elements
    for i in range(n):
        M = np.zeros((n, n), dtype=complex)
        M[i, i] = 1.0
        basis.append(M)
    # off-diagonal: real-symmetric and imaginary-antisymmetric parts
    for i in range(n):
        for j in range(i + 1, n):
            M_re = np.zeros((n, n), dtype=complex)
            M_im = np.zeros((n, n), dtype=complex)
            M_re[i, j] = 1.0
            M_re[j, i] = 1.0
            M_im[i, j] = 1.0j
            M_im[j, i] = -1.0j
            basis.append(M_re)
            basis.append(M_im)
    return basis


def filter_basis_gamma_J(basis_D, gamma: np.ndarray, S: np.ndarray,
                         tol: float = 1e-12):
    """
    Filter Hermitian basis elements to those that are:
      - gamma-odd:  {gamma, E} = 0
      - J-even:     S E* S^T - E = 0
    Returns filtered_basis.
    """
    filtered = []
    for E in basis_D:
        Cg = gamma @ E + E @ gamma
        CJ = S @ E.conj() @ S.T - E
        if (np.linalg.norm(Cg, ord="fro") < tol and
            np.linalg.norm(CJ, ord="fro") < tol):
            filtered.append(E)
    return filtered

# ============================================================
# 1. Basis for 1-generation internal Hilbert space
# ============================================================

@dataclass
class BasisState:
    name: str          # e.g. "nu_L", "e_L", "u_L_r", "d_R_b", ...
    chirality: str     # "L" or "R"
    particle: bool     # True = particle, False = antiparticle
    is_quark: bool     # True for quarks, False for leptons
    # before: color: Optional[str]
    color: Optional[Tuple[str, str]]  # (left_color, right_color) or None
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

    # --- Particle Quarks as color bimodules: 3x3 = 9 states each species ---
    colors = ["r", "g", "b"]

    def add_quark(species: str, chir: str, particle: bool):
        for cl in colors:
            for cr in colors:
                add(f"{species}_{chir}_{cl}{cr}", chir, particle, True, color=(cl, cr))

    # Q_L doublets split by species; we keep species labels explicit
    add_quark("u", "L", True)
    add_quark("d", "L", True)
    add_quark("u", "R", True)
    add_quark("d", "R", True)


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

def gell_mann():
    """Return a small SU(3) generator set (not normalized)."""
    zero = 0.0 + 0.0j
    one  = 1.0 + 0.0j
    i    = 1.0j

    lam1 = np.array([[0,1,0],[1,0,0],[0,0,0]], dtype=complex)
    lam2 = np.array([[0,-i,0],[i,0,0],[0,0,0]], dtype=complex)
    lam3 = np.array([[1,0,0],[0,-1,0],[0,0,0]], dtype=complex)
    lam4 = np.array([[0,0,1],[0,0,0],[1,0,0]], dtype=complex)
    lam5 = np.array([[0,0,-i],[0,0,0],[i,0,0]], dtype=complex)
    lam6 = np.array([[0,0,0],[0,0,1],[0,1,0]], dtype=complex)
    lam7 = np.array([[0,0,0],[0,0,-i],[0,i,0]], dtype=complex)
    lam8 = np.array([[1,0,0],[0,1,0],[0,0,-2]], dtype=complex)

    return {"lam1": lam1, "lam2": lam2, "lam3": lam3, "lam6": lam6, "lam8": lam8}


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

    # M3-part: faithful bimodule color action on quarks
    # Particle quarks:  L(m3) = m3 ⊗ I3  on (cl, cr)
    # Antiparticle quarks: R(m3) = I3 ⊗ m3^T on (cl, cr)
    I3 = np.eye(3, dtype=complex)
    cpos = {"r": 0, "g": 1, "b": 2}

    def block_indices(species: str, chir: str, conj: bool) -> Optional[List[int]]:
        suffix = "_c" if conj else ""
        names = []
        for cl in ["r", "g", "b"]:
            for cr in ["r", "g", "b"]:
                names.append(f"{species}_{chir}_{cl}{cr}{suffix}")
        if all(n in idx for n in names):
            return [idx[n] for n in names]
        return None

    def add_color_block(species: str, chir: str, conj: bool):
        ii = block_indices(species, chir, conj)
        if ii is None:
            return
        if not conj:
            B = np.kron(m3, I3)  # left action on particles
        else:
            B = np.kron(I3, m3.T)  # right action on antiparticles
        A[np.ix_(ii, ii)] += B

    for conj in (False, True):
        for chir in ("L", "R"):
            add_color_block("u", chir, conj)
            add_color_block("d", chir, conj)

    return A


def build_SM_algebra_generators(
    basis: List[BasisState],
    idx: Dict[str, int],
) -> List[SMAlgebraElement]:
    """
    Upgraded generating set for A_F = C ⊕ H ⊕ M_3(C), represented on H_F.

    What’s improved vs your draft:
      - Uses *true* direct-sum generators (separates C, H, M3 identities)
      - Adds a full Hermitian SU(2) generator set (Pauli σ1,σ2,σ3) on doublets
      - Adds several non-diagonal SU(3) (Gell-Mann) generators to force color mixing
      - Ensures every returned operator is Hermitian numerically
    """
    # Zero elements for the summands
    q0  = np.zeros((2, 2), dtype=complex)
    m30 = np.zeros((3, 3), dtype=complex)

    # SU(2) (Hermitian) generators acting on L-doublets
    qI = np.eye(2, dtype=complex)
    sigma1 = np.array([[0, 1], [1, 0]], dtype=complex)
    sigma2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sigma3 = np.array([[1, 0], [0, -1]], dtype=complex)

    # SU(3) color generators (Hermitian) acting on quark color triplets
    gm = gell_mann()
    m3I = np.eye(3, dtype=complex)

    # Pick a small but genuinely noncommuting / non-diagonal subset
    color_gens: list[tuple[str, np.ndarray]] = []
    for key in ("lam1", "lam2", "lam3", "lam6", "lam8"):
        if key not in gm:
            raise KeyError(f"gell_mann() did not provide '{key}'. Available: {list(gm.keys())}")
        color_gens.append((f"color_{key}", gm[key]))

    # Build generators as (label, lambda, q, m3)
    # Note: these are *direct-sum* generators; do NOT collapse them into one “identity triple”.
    gen_specs: list[tuple[str, complex, np.ndarray, np.ndarray]] = [
        # C-summand (scalar on all states)
        ("C_1", 1.0 + 0j, q0,  m30),

        # H-summand (acts on SU(2)_L doublets)
        ("H_1",        0.0 + 0j, qI,     m30),
        ("H_sigma1",   0.0 + 0j, sigma1, m30),
        ("H_sigma2",   0.0 + 0j, sigma2, m30),
        ("H_sigma3",   0.0 + 0j, sigma3, m30),

        # M3-summand (acts on color triplets)
        ("M3_1", 0.0 + 0j, q0,  m3I),
    ]

    # Add color generators
    for lab, m3 in color_gens:
        gen_specs.append((lab, 0.0 + 0j, q0, m3))

    ops: List[SMAlgebraElement] = []

    for label, lam, q, m3 in gen_specs:
        A = rep_A_SM(lam, q, m3, basis, idx)

        # Enforce Hermiticity numerically (your tests assume *-structure sanity)
        A = 0.5 * (A + A.conj().T)

        ops.append(SMAlgebraElement(label, A))

    # Convenience: include the literal identity on H_F (useful for debugging)
    # (This is redundant with C_1 in your simplified rep, but harmless and explicit.)
    dimH = len(basis)
    ops.append(SMAlgebraElement("I_HF", np.eye(dimH, dtype=complex)))

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

    # Quarks (3x3 color bimodule indices)
    for cl in ["r", "g", "b"]:
        for cr in ["r", "g", "b"]:
            couple(f"u_L_{cl}{cr}", f"u_R_{cl}{cr}", Y_u)
            couple(f"d_L_{cl}{cr}", f"d_R_{cl}{cr}", Y_d)

    # Antiparticles mirror
    for cl in ["r", "g", "b"]:
        for cr in ["r", "g", "b"]:
            couple(f"u_L_{cl}{cr}_c", f"u_R_{cl}{cr}_c", Y_u)
            couple(f"d_L_{cl}{cr}_c", f"d_R_{cl}{cr}_c", Y_d)

    return D





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
    df_norm: float               # ||D_F||_F
    rel_plus: float              # norm_plus / ||D_F||_F
    rel_minus: float             # norm_minus / ||D_F||_F
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

def evenize_geom_dirac_and_grading(D_geom: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Standard evenization of an odd geometric triple:
      H -> H ⊕ H
      D_even = [[0, D],[D, 0]]
      gamma  = [[I, 0],[0,-I]]
      lift(A)= [[A, 0],[0, A]]

    Returns: (D_even, gamma_geom, lift_matrix)
      lift_matrix is a 2x2 block operator that can lift any A via:
        A_even = lift(A) = L @ kron(A, I) ??? (we implement explicit lift below)
    """
    dim = D_geom.shape[0]
    I = np.eye(dim, dtype=complex)
    Z = np.zeros((dim, dim), dtype=complex)

    D_even = np.block([[Z, D_geom.astype(complex)],
                       [D_geom.astype(complex), Z]])
    gamma = np.block([[ I, Z],
                      [ Z,-I]])

    # We return gamma; lifting is done by lift_even(A) below.
    return D_even, gamma, I  # I returned only for convenience


def lift_even(A: np.ndarray) -> np.ndarray:
    """Lift an operator A on H to diag(A, A) on H ⊕ H."""
    dim = A.shape[0]
    Z = np.zeros((dim, dim), dtype=complex)
    return np.block([[A.astype(complex), Z],
                     [Z, A.astype(complex)]])
def build_product_dirac_v4(D_geom_even: np.ndarray, gamma_geom: np.ndarray, D_F: np.ndarray) -> np.ndarray:
    """
    v4.0 even product:
      D = D_geom_even ⊗ I_F + gamma_geom ⊗ D_F
    """
    dim_geom = D_geom_even.shape[0]
    dimF = D_F.shape[0]
    I_F = np.eye(dimF, dtype=complex)
    return kron(D_geom_even, I_F) + kron(gamma_geom.astype(complex), D_F.astype(complex))

def build_geom_algebra_generators(N: int) -> dict:
    dim = 2 * N + 1
    n_vals = np.arange(-N, N + 1, dtype=int)
    I_geom = np.eye(dim, dtype=complex)

    gens = {"I_geom": I_geom}

    def proj_div(d: int) -> np.ndarray:
        mask = (n_vals % d == 0)
        return np.diag(mask.astype(float))

    for d in [2, 3, 5]:
        gens[f"P_div_{d}"] = proj_div(d)

    return gens

def build_product_algebra(N: int, *, even_geom: bool) -> tuple[list[np.ndarray], list[str], list[str]]:
    """
    Returns:
      ops_prod, labels_prod, tags_prod
    where tags_prod ∈ {"geom", "finite"} for category diagnostics.
    """
    geom_gens = build_geom_algebra_generators(N)  # acts on H_geom (undoubled)
    I_geom = geom_gens["I_geom"]

    ops_F, labels_F = build_internal_algebra_ops()
    dimF = ops_F[0].shape[0]
    I_F = np.eye(dimF, dtype=complex)

    ops_prod: list[np.ndarray] = []
    labels_prod: list[str] = []
    tags_prod: list[str] = []

    # Lift geom generators if we are in v4 (evenized geom)
    for name, A_geom in geom_gens.items():
        A = lift_even(A_geom) if even_geom else A_geom
        ops_prod.append(kron(A, I_F))
        labels_prod.append(f"{name}⊗I_F")
        tags_prod.append("geom")

    # Finite part always acts as identity on geom factor (even or odd)
    I_geom_used = lift_even(I_geom) if even_geom else I_geom

    for A_F, lab in zip(ops_F, labels_F):
        ops_prod.append(kron(I_geom_used, A_F))
        labels_prod.append(f"I_geom⊗{lab}")
        tags_prod.append("finite")

    return ops_prod, labels_prod, tags_prod

def check_v4_grading_axioms(
    D_geom_even: np.ndarray,
    gamma_geom: np.ndarray,
    ops_prod: list[np.ndarray],
    tags_prod: list[str],
    dimF: int,
    eps: float = 1e-12,
) -> None:
    """
    v4.0 core sanity:
      {gamma_geom, D_geom_even} = 0
      [gamma_geom, pi(a_geom)] = 0   (geom algebra commutes with gamma)
      [gamma_geom, I_geom] = 0       (trivial)
    We check these after tensoring to the full Hilbert space where needed.
    """
    print("=== v4.0 grading sanity checks ===")

    # 1) geometric anti-commutation
    anti = gamma_geom @ D_geom_even + D_geom_even @ gamma_geom
    n1 = np.linalg.norm(anti, ord="fro")
    print(f"||{{γ_geom, D_geom_even}}||_F = {n1:.3e}  (target ~0)")

    # 2) gamma commutes with geom algebra (lifted), checked on full product reps that are tagged "geom"
    # Build γ_full = γ_geom ⊗ I_F
    I_F = np.eye(dimF, dtype=complex)
    gamma_full = kron(gamma_geom.astype(complex), I_F)

    max_comm = 0.0
    for A, tag in zip(ops_prod, tags_prod):
        if tag != "geom":
            continue
        comm = gamma_full @ A - A @ gamma_full
        max_comm = max(max_comm, np.linalg.norm(comm, ord="fro"))

    print(f"max ||[γ_full, π(a_geom)]||_F = {max_comm:.3e}  (target ~0)")
    if n1 > eps or max_comm > eps:
        print(f"WARNING: grading sanity exceeds eps={eps:.1e}")
    print()

# ===========================================
# Order-zero and first-order checkers (finite triple)
# ===========================================

def check_first_order_condition(
    D_F: np.ndarray,
    algebra: List[AlgebraElement],
    eps: float = 1e-12,
    J_phase: complex = 1.0,
    S: np.ndarray | None = None,
    topk: int = 10,
) -> FirstOrderResult:
    """
    First-order: [[D, a], J b J^{-1}] = 0 for all a,b in algebra.

    If S is None, default to the simple H_L ⊕ H_R swap.
    """
    assert_square(D_F, "D_F")
    n = D_F.shape[0]

    if S is None:
        assert_even_dim(D_F, "D_F")
        S = build_swap_LR(n // 2)

    max_norm = 0.0
    worst_pair: Optional[Tuple[str, str]] = None

    # We store violating pairs (label_a, label_b, norm) up to topk largest.
    viols: List[Tuple[str, str, float]] = []

    # Precompute J b J^{-1}
    Jb_list = [J_action(b.op, S, phase=J_phase) for b in algebra]

    for a in algebra:
        Da = D_F @ a.op - a.op @ D_F
        for b, Jb in zip(algebra, Jb_list):
            comm2 = Da @ Jb - Jb @ Da
            nrm = frob_norm(comm2)

            if nrm > max_norm:
                max_norm = nrm
                worst_pair = (a.label, b.label)

            if nrm > eps:
                viols.append((a.label, b.label, nrm))

    # keep only topk by norm
    viols.sort(key=lambda x: x[2], reverse=True)
    viols = viols[:topk]

    return FirstOrderResult(
        max_norm=max_norm,
        worst_pair=worst_pair,
        good_pairs=viols,   # NOTE: field name kept for backward compatibility; these are VIOLATIONS.
    )


def check_zero_order_condition(
    algebra: List[AlgebraElement],
    eps: float = 1e-12,
    J_phase: complex = 1.0,
    S: np.ndarray | None = None,
    topk: int = 10,
) -> ZeroOrderResult:
    """
    Zero-order: [a, J b J^{-1}] = 0 for all a,b in algebra.

    If S is None, default to the simple H_L ⊕ H_R swap.
    """
    # infer n from first algebra element
    if not algebra:
        raise ValueError("algebra must be non-empty")

    n = algebra[0].op.shape[0]
    for a in algebra:
        assert_square(a.op, f"algebra element {a.label}")

    if S is None:
        if n % 2 != 0:
            raise ValueError("Default LR swap requires even-dimensional representation. Provide S explicitly.")
        S = build_swap_LR(n // 2)

    max_norm = 0.0
    worst_pair: Optional[Tuple[str, str]] = None

    bad_pairs: List[Tuple[str, str, float]] = []

    # Precompute J b J^{-1}
    Jb_list = [J_action(b.op, S, phase=J_phase) for b in algebra]

    for a in algebra:
        for b, Jb in zip(algebra, Jb_list):
            comm = a.op @ Jb - Jb @ a.op
            nrm = frob_norm(comm)

            if nrm > max_norm:
                max_norm = nrm
                worst_pair = (a.label, b.label)

            if nrm > eps:
                bad_pairs.append((a.label, b.label, nrm))

    bad_pairs.sort(key=lambda x: x[2], reverse=True)
    bad_pairs = bad_pairs[:topk]

    return ZeroOrderResult(
        max_norm=max_norm,
        worst_pair=worst_pair,
        bad_pairs=bad_pairs,
    )


def print_first_order_result(res: FirstOrderResult, eps: float = 1e-12) -> None:
    print("=== First-order condition ===")
    print(f"max ||[[D,a], JbJ^-1]||_F = {res.max_norm:.3e}")
    if res.worst_pair is not None:
        print(f"worst pair: a={res.worst_pair[0]}, b={res.worst_pair[1]}")
    if res.good_pairs:  # these are actually violations (see note in checker)
        print(f"Top violations (> eps={eps:.1e}):")
        for la, lb, nrm in res.good_pairs:
            print(f"  a={la:>12s}, b={lb:>12s}  -> {nrm:.3e}")
    else:
        print(f"All pairs satisfy first-order within eps={eps:.1e}")
    print()


def print_zero_order_result(res: ZeroOrderResult, eps: float = 1e-12) -> None:
    print("=== Zero-order condition ===")
    print(f"max ||[a, JbJ^-1]||_F = {res.max_norm:.3e}")
    if res.worst_pair is not None:
        print(f"worst pair: a={res.worst_pair[0]}, b={res.worst_pair[1]}")
    if res.bad_pairs:
        print(f"Top violations (> eps={eps:.1e}):")
        for la, lb, nrm in res.bad_pairs:
            print(f"  a={la:>12s}, b={lb:>12s}  -> {nrm:.3e}")
    else:
        print(f"All pairs satisfy zero-order within eps={eps:.1e}")
    print()

# ===========================================
# First-order condition
# ===========================================

def _topk_insert(lst, item, k):
    # lst is list[(norm, info)] sorted descending by norm
    lst.append(item)
    lst.sort(key=lambda x: x[0], reverse=True)
    del lst[k:]


def test_first_order_condition_product(
    D: np.ndarray,
    ops: list[np.ndarray],
    labels: list[str],
    S_prod: np.ndarray,
    eps: float = 1e-12,
    tags: list[str] | None = None,
    topk: int = 10,
) -> None:
    print("=== First-order condition test (v4 gate) ===")
    max_norm = 0.0
    worst: list[tuple[float, str]] = []

    # Optional category maxima
    cat_max = {
        ("geom","geom"): 0.0,
        ("geom","finite"): 0.0,
        ("finite","geom"): 0.0,
        ("finite","finite"): 0.0,
    } if tags is not None else None

    for i, a in enumerate(ops):
        Da = D @ a - a @ D
        for j, b in enumerate(ops):
            b_tilde = J_action(S_prod, b)
            comm2 = Da @ b_tilde - b_tilde @ Da
            norm = np.linalg.norm(comm2, ord="fro")

            if norm > max_norm:
                max_norm = norm

            if norm > eps:
                info = f"(a={labels[i]}, b={labels[j]})  ||[[D,a], JbJ^-1]||_F={norm:.3e}"
                _topk_insert(worst, (norm, info), topk)

            if cat_max is not None:
                key = (tags[i], tags[j])
                cat_max[key] = max(cat_max[key], norm)

    print(f"Max Frobenius norm over all pairs (a,b): {max_norm:.3e}")
    if cat_max is not None:
        print("Category maxima (a-tag, b-tag):")
        for k, v in cat_max.items():
            print(f"  {k}: {v:.3e}")

    if worst:
        print(f"Top violations > eps={eps:.1e}:")
        for _, info in worst:
            print(" ", info)
    else:
        print(f"All pairs satisfy first-order within eps={eps:.1e}")
    print()


def test_zero_order_condition_product(
    ops: list[np.ndarray],
    labels: list[str],
    S_prod: np.ndarray,
    eps: float = 1e-12,
    tags: list[str] | None = None,
    topk: int = 10,
) -> None:
    print("=== Zero-order condition test (v4 gate) ===")
    max_norm = 0.0
    worst: list[tuple[float, str]] = []

    cat_max = {
        ("geom","geom"): 0.0,
        ("geom","finite"): 0.0,
        ("finite","geom"): 0.0,
        ("finite","finite"): 0.0,
    } if tags is not None else None

    for i, a in enumerate(ops):
        for j, b in enumerate(ops):
            b_tilde = J_action(S_prod, b)
            comm = a @ b_tilde - b_tilde @ a
            norm = np.linalg.norm(comm, ord="fro")

            max_norm = max(max_norm, norm)

            if norm > eps:
                info = f"(a={labels[i]}, b={labels[j]})  ||[a, JbJ^-1]||_F={norm:.3e}"
                _topk_insert(worst, (norm, info), topk)

            if cat_max is not None:
                key = (tags[i], tags[j])
                cat_max[key] = max(cat_max[key], norm)

    print(f"Max Frobenius norm over all pairs (a,b): {max_norm:.3e}")
    if cat_max is not None:
        print("Category maxima (a-tag, b-tag):")
        for k, v in cat_max.items():
            print(f"  {k}: {v:.3e}")

    if worst:
        print(f"Top violations > eps={eps:.1e}:")
        for _, info in worst:
            print(" ", info)
    else:
        print(f"All pairs satisfy zero-order within eps={eps:.1e}")
    print()


# ===========================================
# Grading & reality
# ===========================================

# in ncg_tests.py

def check_grading_and_reality(
    D_F: np.ndarray,
    algebra: List[AlgebraElement],
    J_phase: complex = 1.0,
    gamma_F: np.ndarray | None = None,
    S: np.ndarray | None = None,
) -> GradingRealityResult:
    """
    - Check gamma_F anticommutes with D_F and commutes with A_F.
    - Check J^2 = 1 (via S^2).
    - Estimate KO-dimension sign via J D_F J^{-1} = ± D_F.

    If gamma_F or S are None, fall back to the simple H_L ⊕ H_R split.
    """
    assert_even_dim(D_F, "D_F")
    n = D_F.shape[0]
    dim_left = n // 2

    if gamma_F is None:
        gamma_F = build_gamma_F(dim_left)      # old behavior
    if S is None:
        S = build_swap_LR(dim_left)            # old behavior

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
    norm_plus  = frob_norm(JDJ - D_F)
    norm_minus = frob_norm(JDJ + D_F)

    df_norm = frob_norm(D_F)
    scale = df_norm + 1e-16

    rel_plus = norm_plus / scale
    rel_minus = norm_minus / scale


    if rel_plus < 1e-12 and rel_minus > rel_plus:
        ko_sign: int | None = +1
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
        df_norm=df_norm,
        rel_plus=rel_plus,
        rel_minus=rel_minus,
        ko_sign=ko_sign,
    )


def print_grading_and_reality_result(
    result: GradingRealityResult,
    title: str = "Grading & reality tests",
) -> None:
    print(f"=== {title} ===")
    print(f"||{{gamma_F, D_F}}||_F         = {result.anticom_norm:.3e}")
    print(f"max ||[gamma_F, a]||_F         = {result.max_comm_gamma:.3e}")
    print(f"||J^2 - I||_F                  = {result.J2_deviation:.3e}")
    print(f"||D_F||_F                      = {result.df_norm:.3e}")
    print(f"||J D_F J^-1 - D_F||_F         = {result.norm_plus:.3e}   (rel {result.rel_plus:.3e})")
    print(f"||J D_F J^-1 + D_F||_F         = {result.norm_minus:.3e}  (rel {result.rel_minus:.3e})")

    if result.ko_sign == +1:
        print("→ KO-sign: J D_F J^-1 = + D_F (J-even). Expect rel_minus ≈ 2.")
    elif result.ko_sign == -1:
        print("→ KO-sign: J D_F J^-1 = - D_F (J-odd). Expect rel_plus ≈ 2.")
    else:
        print("→ KO-sign: ambiguous at numerical precision.")
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
    gamma_F: np.ndarray | None = None,
    S: np.ndarray | None = None,
    name: str = "",
) -> Tuple[FirstOrderResult, ZeroOrderResult, GradingRealityResult]:
    """
    Run first-order, zero-order, and grading/reality tests and pretty-print a summary.

    If gamma_F and/or S are provided, they are used for grading/reality;
    and S is used for order tests (so your SM particle↔antiparticle J is honored).
    """
    if name:
        print(f"=== NCG test suite for {name} ===")

    fo_res = check_first_order_condition(D_F, algebra, eps=eps_first, J_phase=J_phase, S=S, topk=TOPK_VIOLATIONS)
    print_first_order_result(fo_res, eps=eps_first)

    zo_res = check_zero_order_condition(algebra, eps=eps_zero, J_phase=J_phase, S=S, topk=TOPK_VIOLATIONS)
    print_zero_order_result(zo_res, eps=eps_zero)

    gr_res = check_grading_and_reality(D_F, algebra, J_phase=J_phase, gamma_F=gamma_F, S=S)
    print_grading_and_reality_result(gr_res)

    return fo_res, zo_res, gr_res



def hermitian_parametrization(n: int):
    """
    Return:
      - basis: list of n x n Hermitian matrices E_k
      so that any Hermitian D can be written as D = sum_k x_k E_k
      with real coefficients x_k.
    """
    basis = []
    # diagonal basis
    for i in range(n):
        M = np.zeros((n, n), dtype=complex)
        M[i, i] = 1.0
        basis.append(M)
    # off-diagonal (i<j): real-symmetric and imaginary-antisymmetric parts
    for i in range(n):
        for j in range(i+1, n):
            M_re = np.zeros((n, n), dtype=complex)
            M_im = np.zeros((n, n), dtype=complex)
            M_re[i, j] = 1.0
            M_re[j, i] = 1.0
            M_im[i, j] = 1.0j
            M_im[j, i] = -1.0j
            basis.append(M_re)
            basis.append(M_im)
    return basis  # length = n + 2 * n*(n-1)/2 = n^2
def build_first_order_constraint_matrix(basis_D_reduced,
                                        ops: list[np.ndarray],
                                        S: np.ndarray) -> np.ndarray:
    """
    Build constraint matrix V for first-order condition only:

       [[D, a_i], J a_j J^{-1}] = 0   for all i,j,

    where D = sum_k x_k E_k, with E_k in basis_D_reduced (already gamma-odd, J-even).

    Returns:
      V: complex matrix of shape (num_constraints, num_unknowns)
         such that V @ x = 0 encodes all constraints.
    """
    n = ops[0].shape[0]
    num_unknowns = len(basis_D_reduced)
    num_ops = len(ops)
    block_size = n * n
    num_blocks = num_ops * num_ops
    vec_len = num_blocks * block_size

    V = np.zeros((vec_len, num_unknowns), dtype=complex)

    # Precompute J a_j J^{-1}
    Jops = [S @ a.conj() @ S.T for a in ops]

    for k, E in enumerate(basis_D_reduced):
        offset = 0
        for a in ops:
            Da = E @ a - a @ E
            for btilde in Jops:
                C_first = Da @ btilde - btilde @ Da
                V[offset:offset + block_size, k] = C_first.reshape(-1)
                offset += block_size

    return V
def find_DF_solution_basis(basis_D_reduced, V: np.ndarray,
                           tol: float = 1e-12) -> list[np.ndarray]:
    """
    Given reduced basis {E_k} and constraint matrix V (M x K),
    find a basis of Hermitian matrices D_alpha = sum_k x_k^{(alpha)} E_k
    spanning the nullspace of V @ x = 0.

    Returns:
      DF_basis: list of n x n Hermitian matrices D_alpha.
    """
    # SVD: V = U diag(s) Vh, right-singular vectors in rows of Vh
    U, s, Vh = np.linalg.svd(V, full_matrices=False)

    null_mask = (s < tol)
    null_vectors = Vh[null_mask, :]   # shape: (num_null, K)

    DF_basis = []
    for vec in null_vectors:
        D = np.zeros_like(basis_D_reduced[0])
        for coeff, E in zip(vec, basis_D_reduced):
            D += coeff * E
        # ensure Hermitian numerically
        D = 0.5 * (D + D.conj().T)
        DF_basis.append(D)

    return DF_basis
def print_coupling_pattern(D: np.ndarray,
                           basis: List[SMState],
                           name_to_index: Dict[str, int],
                           thresh: float = 1e-6):
    """
    Print which L↔R pairs (for each species) are significantly coupled by D.
    """
    def show_pair(a, b):
        i = name_to_index[a]
        j = name_to_index[b]
        z = D[i, j]
        if abs(z) > thresh:
            print(f"  {a:>10s} <-> {b:<10s}  |D_ij| = {abs(z):.3e}")

    print("Lepton couplings:")
    if "nu_R" in name_to_index:
        show_pair("nu_L", "nu_R")
    show_pair("e_L", "e_R")

    print("Quark couplings (per color):")
    colors = ["r", "g", "b"]
    for c in colors:
        show_pair(f"u_L_{c}", f"u_R_{c}")
    for c in colors:
        show_pair(f"d_L_{c}", f"d_R_{c}")

    print("Antiparticle couplings (sanity check):")
    if "bar_nu_R" in name_to_index:
        show_pair("bar_nu_L", "bar_nu_R")
    show_pair("bar_e_L", "bar_e_R")
    for c in colors:
        show_pair(f"bar_u_L_{c}", f"bar_u_R_{c}")
    for c in colors:
        show_pair(f"bar_d_L_{c}", f"bar_d_R_{c}")

def main():
    basis, idx = build_sm_basis_1gen(include_nu_R=True)

    # SM-like Yukawas (1 generation, rough magnitudes)
    Y_e  = 2.94e-6     # me / v
    Y_nu = 1.0e-12     # tiny
    Y_u  = 1.3e-5      # up-type
    Y_d  = 2.8e-5      # down-type

    D_F_SM = build_DF_SM_1gen(Y_e, Y_nu, Y_u, Y_d, basis, idx)
    sm_generators = build_SM_algebra_generators(basis, idx)
    algebra_SM = [AlgebraElement(op.label, op.op) for op in sm_generators]

    # SM-specific gamma and J swap
    gamma_SM = build_gamma_F_SM(basis)
    S_SM = build_swap_particle_antiparticle(basis, idx)

    # Full suite using SM gamma + SM J
    run_ncg_test_suite(
        D_F_SM,
        algebra_SM,
        eps_first=1e-12,
        eps_zero=1e-12,
        gamma_F=gamma_SM,
        S=S_SM,
        name="SM-like 1gen finite triple (SM gamma,J)",
    )



if __name__ == "__main__":
    main()

