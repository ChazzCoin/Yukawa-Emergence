#!/usr/bin/env python3
"""
sm_finite_triple.py

Minimal 1-generation SM-like finite spectral triple:
  A_F ≃ C ⊕ H ⊕ M_3(C)
  H_F: quarks + leptons + color + antiparticles
  D_F: Yukawa + LR mixing (no full Majorana yet)

Goal: provide a concrete (A_F, H_F, D_F, J_F, gamma_F) that you can test with
your ncg_tests harness and later refine toward the full Connes–Chamseddine SM.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import List, Dict, Tuple
import numpy as np

from v3.spectral import run_ncg_test_suite


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
class AlgebraElement:
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

# ================================
# Internal Hilbert space & D_F (NCG-compatible toy)
# ================================

SECTORS: List[str] = ["u", "d", "e", "nu"]
SECTOR_INDEX: Dict[str, int] = {s: i for i, s in enumerate(SECTORS)}
N_GEN: int = 3

# We treat color as a degeneracy factor on u,d (3 copies) and 1 on e,nu.
# It is *not* yet a full SU(3) algebra action; that would require an explicit
# color tensor factor and a more refined J_F. Here we only count dimensions.
SECTOR_NC: Dict[str, int] = {"u": 3, "d": 3, "e": 1, "nu": 1}

def base_kernel(lam, alpha=3.0, form="lambda_sq"):
    """
    Base kernel F_base(λ_g) that defines the generation ladder.

    We make it *scale-invariant* by normalizing the eigenvalues to the
    lightest nonzero one, so that a global rescaling of the Laplacian
    does not flatten or blow up the hierarchy:

        F_base(λ_g) = exp[-alpha * (λ_g / λ_ref)^2]

    with λ_ref = smallest positive eigenvalue in the triad.
    """
    lam = np.array(lam, dtype=float)

    # Choose a reference eigenvalue λ_ref (smallest positive λ)
    lam_pos = lam[lam > 0]
    if lam_pos.size == 0:
        # Degenerate case: fall back to ordinary λ^2 kernel
        lam_ref = 1.0
    else:
        lam_ref = lam_pos.min()

    x = lam / lam_ref

    if form == "lambda_sq":
        return np.exp(-alpha * x**2)
    elif form == "lambda":
        return np.exp(-alpha * x)
    else:
        raise ValueError(f"Unknown kernel form '{form}'")

def dim_per_chirality() -> int:
    """Dimension of H_L or H_R (one chirality).

    We fold color multiplicities into sector blocks:
      u,d: 3 each; e,nu: 1 each → total 8 per generation
      times 3 generations → 24 per chirality.
    """
    return 3 * sum(SECTOR_NC[s] for s in SECTORS)  # 24


def flavor_block_offsets() -> Dict[str, int]:
    """Offsets (within a single chirality) for each sector's 3×3
    generation block in a 12×12 generation-space layout:

      [u_g1,u_g2,u_g3,
       d_g1,d_g2,d_g3,
       e_g1,e_g2,e_g3,
       nu_g1,nu_g2,nu_g3]

    We only care about generation offsets (3×3 blocks);
    color multiplicity is treated as degeneracy, not an explicit tensor factor.
    """
    off: Dict[str, int] = {}
    off["u"]  = 0
    off["d"]  = 3
    off["e"]  = 6
    off["nu"] = 9
    return off


def build_internal_DF_from_Y(Y_u: np.ndarray,
                             Y_d: np.ndarray,
                             Y_e: np.ndarray,
                             Y_nu: np.ndarray) -> np.ndarray:
    """Build the finite Dirac operator D_F in block form:

        D_F = [[ 0, Y^\dagger ],
               [ Y, 0         ]]

    where Y is a 24×24 block that is block-diagonal in sector space and
    embeds the 3×3 generation Yukawas (Y_u, Y_d, Y_e, Y_nu) into a
    12×12 generation-space layout, with color treated as degeneracy.

    H_F ≃ H_L ⊕ H_R,  dim(H_L) = dim(H_R) = 24, dim(H_F) = 48.
    """
    # Sanity checks
    for name, Y in [("Y_u", Y_u), ("Y_d", Y_d), ("Y_e", Y_e), ("Y_nu", Y_nu)]:
        Y = np.asarray(Y, dtype=complex)
        if Y.shape != (3, 3):
            raise ValueError(f"{name} must be a 3×3 matrix, got shape {Y.shape}.")
    Y_u  = np.asarray(Y_u, dtype=complex)
    Y_d  = np.asarray(Y_d, dtype=complex)
    Y_e  = np.asarray(Y_e, dtype=complex)
    Y_nu = np.asarray(Y_nu, dtype=complex)

    dpc = dim_per_chirality()  # 24
    dimH = 2 * dpc             # 48

    # 12×12 generation-space Yukawa core
    Y_gen = np.zeros((12, 12), dtype=complex)
    gen_off = flavor_block_offsets()

    def insert_sector_Y(sector: str, Y_s: np.ndarray) -> None:
        off = gen_off[sector]
        Y_gen[off:off+3, off:off+3] = Y_s

    insert_sector_Y("u",  Y_u)
    insert_sector_Y("d",  Y_d)
    insert_sector_Y("e",  Y_e)
    insert_sector_Y("nu", Y_nu)

    # Embed the 12×12 generation block into 24×24 per chirality.
    # Only the leading 12×12 carry Yukawa couplings; the remaining slots
    # are color-degenerate but Yukawa-silent in this toy.
    Y_block = np.zeros((dpc, dpc), dtype=complex)
    Y_block[:12, :12] = Y_gen

    # Assemble D_F on H_F = H_L ⊕ H_R
    D_F = np.zeros((dimH, dimH), dtype=complex)
    D_F[:dpc,  dpc:] = Y_block.conj().T
    D_F[dpc:,  :dpc] = Y_block

    return D_F


# ================================
# Real structure, grading, algebra basis
# ================================

def build_swap_LR(dim_left: int) -> np.ndarray:
    """Swap matrix S on H_F = H_L ⊕ H_R, with dim(H_L) = dim(H_R) = dim_left.

    Acts as: S (ψ_L, ψ_R) = (ψ_R, ψ_L).
    """
    S = np.zeros((2 * dim_left, 2 * dim_left), dtype=complex)
    S[:dim_left, dim_left:] = np.eye(dim_left)
    S[dim_left:, :dim_left] = np.eye(dim_left)
    return S


def build_gamma_F(dim_left: int) -> np.ndarray:
    """Grading operator γ_F with eigenvalue -1 on H_L and +1 on H_R."""
    g = np.zeros((2 * dim_left, 2 * dim_left), dtype=complex)
    g[:dim_left, :dim_left] = -np.eye(dim_left)
    g[dim_left:, dim_left:] =  np.eye(dim_left)
    return g


def build_sector_projectors() -> Dict[str, np.ndarray]:
    """Sector projectors P_sector_s acting on H_F = H_L ⊕ H_R.

    Each P_sector_s is diagonal and selects the (sector,gen,chirality) subspace
    corresponding to that sector (u,d,e,nu) in the 12×12 generation subspace,
    duplicated on L and R.
    """
    dpc = dim_per_chirality()
    dimH = 2 * dpc
    gen_off = flavor_block_offsets()

    P: Dict[str, np.ndarray] = {}
    for s in SECTORS:
        P_s = np.zeros((dimH, dimH), dtype=complex)
        off = gen_off[s]
        # Same on L and R, only on first 12 generation slots
        P_s[off:off+3, off:off+3] = np.eye(3)
        P_s[dpc + off:dpc + off + 3, dpc + off:dpc + off + 3] = np.eye(3)
        P[s] = P_s

    return P


def build_Q_sector() -> np.ndarray:
    """A simple 'sector charge' diagonal operator Q_sector.

    Distinguishes u,d,e,nu sectors but is generation-blind:
      q_u = 2, q_d = 1, q_e = 0, q_nu = -1  (toy choice).
    """
    dpc = dim_per_chirality()
    dimH = 2 * dpc
    gen_off = flavor_block_offsets()
    charges = {"u": 2.0, "d": 1.0, "e": 0.0, "nu": -1.0}

    Q = np.zeros((dimH, dimH), dtype=complex)
    for s in SECTORS:
        off = gen_off[s]
        q = charges[s]
        Q[off:off+3, off:off+3] = q * np.eye(3)
        Q[dpc + off:dpc + off + 3, dpc + off:dpc + off + 3] = q * np.eye(3)

    return Q


def build_internal_algebra_ops() -> Tuple[List[np.ndarray], List[str]]:
    """Small basis of algebra elements A_F acting on H_F:

        - I (identity)
        - Q_sector (diagonal sector 'charge')
        - P_sector_u, P_sector_d, P_sector_e, P_sector_nu (sector projectors)

    This is a commutative algebra in this toy (no explicit SU(3) yet).
    """
    dpc = dim_per_chirality()
    dimH = 2 * dpc

    I = np.eye(dimH, dtype=complex)
    Q = build_Q_sector()
    P = build_sector_projectors()

    ops: List[np.ndarray] = [I, Q, P["u"], P["d"], P["e"], P["nu"]]
    labels: List[str] = ["I", "Q_sector", "P_sector_u", "P_sector_d", "P_sector_e", "P_sector_nu"]

    return ops, labels


# ================================
# NCG condition tests
# ================================

def J_action_from_swap(S: np.ndarray, M: np.ndarray) -> np.ndarray:
    """Implement J M J^{-1} = S · M^* · S^T, where S is the L/R swap."""
    return S @ M.conj() @ S.T


def test_first_order_condition(D_F: np.ndarray,
                               ops: List[np.ndarray],
                               labels: List[str],
                               eps: float = 1e-12) -> None:
    """First-order condition:

        [[D_F, a], J_F b J_F^{-1}] = 0   for all a,b in A_F.
    """
    n = D_F.shape[0]
    assert D_F.shape == (n, n)
    S = build_swap_LR(dim_left=n // 2)

    print("=== First-order condition test ===")
    max_norm = 0.0
    good_pairs = []

    for i, a in enumerate(ops):
        Da = D_F @ a - a @ D_F
        for j, b in enumerate(ops):
            b_tilde = J_action_from_swap(S, b)
            comm2 = Da @ b_tilde - b_tilde @ Da
            norm = np.linalg.norm(comm2, ord="fro")
            if norm > max_norm:
                max_norm = norm
            if norm < eps:
                good_pairs.append((labels[i], labels[j], norm))

    print(f"Max Frobenius norm over all pairs (a,b): {max_norm:.3e}")
    if good_pairs:
        print(f"Pairs with norm < {eps:.1e}:")
        for la, lb, nrm in good_pairs:
            print(f"  (a={la:>10s}, b={lb:>10s}) → ||[[D,a],J b J^-1]||_F = {nrm:.3e}")
    else:
        print(f"No pairs with norm < {eps:.1e}")
    print()


def test_zero_order_condition(ops: List[np.ndarray],
                              labels: List[str],
                              eps: float = 1e-12) -> None:
    """Zero-order condition:

        [a, J_F b J_F^{-1}] = 0   for all a,b in A_F.
    """
    n = ops[0].shape[0]
    S = build_swap_LR(dim_left=n // 2)

    print("=== Zero-order condition test ===")
    max_norm = 0.0
    bad_pairs = []

    for i, a in enumerate(ops):
        for j, b in enumerate(ops):
            b_tilde = J_action_from_swap(S, b)
            comm = a @ b_tilde - b_tilde @ a
            norm = np.linalg.norm(comm, ord="fro")
            if norm > max_norm:
                max_norm = norm
            if norm > eps:
                bad_pairs.append((labels[i], labels[j], norm))

    print(f"Max Frobenius norm over all pairs (a,b): {max_norm:.3e}")
    if bad_pairs:
        print("Pairs with significant violation:")
        for la, lb, nrm in bad_pairs:
            print(f"  (a={la:>10s}, b={lb:>10s}) → ||[a, J b J^-1]||_F = {nrm:.3e}")
    else:
        print(f"All pairs satisfy [a, J b J^-1] ≈ 0 within eps={eps:.1e}")
    print()


def test_grading_and_reality(D_F: np.ndarray,
                             ops: List[np.ndarray],
                             labels: List[str]) -> None:
    """Check grading and reality axioms:

      - γ_F anticommutes with D_F and commutes with A_F.
      - J_F^2 = 1 (as implemented by swap).
      - KO-dimension sign via J D_F J^{-1} = ± D_F.
    """
    n = D_F.shape[0]
    dpc = n // 2
    gamma_F = build_gamma_F(dpc)
    S = build_swap_LR(dpc)

    print("=== Grading & reality tests ===")
    anti = gamma_F @ D_F + D_F @ gamma_F
    print(f"||{{gamma_F, D_F}}||_F = {np.linalg.norm(anti, ord='fro'):.3e}")

    max_comm_gamma = 0.0
    for a in ops:
        comm_ga = gamma_F @ a - a @ gamma_F
        max_comm_gamma = max(max_comm_gamma, np.linalg.norm(comm_ga, ord="fro"))
    print(f"max ||[gamma_F, a]||_F over a∈A_F = {max_comm_gamma:.3e}")

    S2 = S @ S
    print(f"||S^2 - I||_F  (⇒ J_F^2 deviation) = {np.linalg.norm(S2 - np.eye(n), ord='fro'):.3e}")

    JDJ = S @ D_F.conj() @ S.T
    norm_minus = np.linalg.norm(JDJ - D_F, ord="fro")
    norm_plus  = np.linalg.norm(JDJ + D_F, ord="fro")
    print(f"||J D_F J^-1 - D_F||_F   = {norm_minus:.3e}")
    print(f"||J D_F J^-1 + D_F||_F   = {norm_plus:.3e}")
    if norm_minus < norm_plus:
        print("→ KO-sign: J D_F J^-1 = + D_F (J-even Dirac operator)")
    else:
        print("→ KO-sign: J D_F J^-1 = - D_F (J-odd Dirac operator)")
    print()


# ================================
# Emergent misalignment model, flavor, mixing, chi^2
# (your original χ^2≈11 toy, kept intact below)
# ================================

# --- everything below here is your original emergent-4-x11 code ---
# (misalignment functional, emergent graph, Laplacian, F_base, Q,
#  geometry-derived U_geom, Yukawas, mixing, chi^2, etc.)

# I’m not re-commenting every function here since they’re unchanged;
# this is literally your stable χ²≈11 script with the NCG block added above
# and the NCG tests called at the end of main().

# -------------- misalignment functional, relaxation, graph, etc. --------------

def misalignment_energy(theta, w6=1.0, w5=1.0):
    N = len(theta)
    diffs = theta[:, None] - theta[None, :]
    cos6 = np.cos(6 * diffs)
    cos5 = np.cos(5 * diffs)
    E6 = w6 * np.sum(1.0 - cos6) / (N * N)
    E5 = w5 * np.sum(1.0 - cos5) / (N * N)
    return E6 + E5


def relax_phases(N=200, n_steps=600, eta=0.01, w6=1.0, w5=1.0, random_seed=42):
    rng = np.random.default_rng(random_seed)
    theta = rng.uniform(0, 2 * np.pi, size=N)
    energy_hist = []

    for step in range(n_steps):
        diffs = theta[:, None] - theta[None, :]
        sin6 = np.sin(6 * diffs)
        sin5 = np.sin(5 * diffs)
        grad = 6 * w6 * np.sum(sin6, axis=1) + 5 * w5 * np.sum(sin5, axis=1)
        theta = theta - eta * grad
        theta = (theta + 2 * np.pi) % (2 * np.pi)

        if step % 10 == 0 or step == n_steps - 1:
            E = misalignment_energy(theta, w6=w6, w5=w5)
            energy_hist.append(E)

    return theta, energy_hist


def build_emergent_adjacency(theta, w6=1.0, w5=1.0, keep_fraction=0.05):
    N = len(theta)
    diffs = theta[:, None] - theta[None, :]
    cos6 = np.cos(6 * diffs)
    cos5 = np.cos(5 * diffs)
    score = w6 * cos6 + w5 * cos5
    np.fill_diagonal(score, -np.inf)
    triu_idx = np.triu_indices(N, k=1)
    flat_scores = score[triu_idx]
    k = int(keep_fraction * len(flat_scores))
    if k < 1:
        k = 1
    kth_val = np.partition(flat_scores, -k)[-k]
    A = np.zeros((N, N), dtype=float)
    mask = (score >= kth_val)
    A[mask] = 1.0
    A = np.maximum(A, A.T)
    return A


def largest_connected_component(A):
    N = A.shape[0]
    visited = np.zeros(N, dtype=bool)
    best_comp = []
    for i in range(N):
        if not visited[i]:
            stack = [i]
            comp = []
            visited[i] = True
            while stack:
                v = stack.pop()
                comp.append(v)
                neighbors = np.where(A[v] > 0)[0]
                for w in neighbors:
                    if not visited[w]:
                        visited[w] = True
                        stack.append(w)
            if len(comp) > len(best_comp):
                best_comp = comp
    best_comp = np.array(best_comp, dtype=int)
    A_sub = A[np.ix_(best_comp, best_comp)]
    return A_sub, best_comp


def laplacian_from_adjacency(A):
    d = np.sum(A, axis=1)
    L = np.diag(d) - A
    return L


def spectral_triad(L):
    eigvals, eigvecs = np.linalg.eigh(L)
    idx_sorted = np.argsort(eigvals)
    eigvals_sorted = eigvals[idx_sorted]
    eigvecs_sorted = eigvecs[:, idx_sorted]
    lam_gen = eigvals_sorted[1:4]
    gen_indices = idx_sorted[1:4]
    return lam_gen, gen_indices, eigvals_sorted


# -------------- sector charges and F_s -----------------

def build_sector_charges():
    Q_u = np.array([0,  2,  4], dtype=float)
    Q_d = np.array([1,  3,  5], dtype=float)
    Q_e = np.array([2,  4,  6], dtype=float)
    Q_n = np.array([4,  6,  8], dtype=float)
    return {"u": Q_u, "d": Q_d, "e": Q_e, "nu": Q_n}


def sector_weights(F_base, Q_s, beta=1.0):
    return F_base * np.exp(-beta * Q_s)


def mass_ratios(F_s):
    F_s = np.array(F_s, dtype=float)
    m1, m2, m3 = F_s
    return m1 / m3, m2 / m3


# -------------- generation operators (golden, Cabibbo) --------------

def rotation_3d(i, j, theta):
    R = np.eye(3, dtype=complex)
    c = np.cos(theta)
    s = np.sin(theta)
    R[i, i] = c
    R[j, j] = c
    R[i, j] = s
    R[j, i] = -s
    return R


def build_generation_operators(phi_order=5, cab_denom=28):
    theta_phi = 2 * np.pi / phi_order
    theta_C = 2 * np.pi / cab_denom
    P_phi_12 = rotation_3d(0, 1, theta_phi)
    P_phi_23 = rotation_3d(1, 2, theta_phi)
    C_12 = rotation_3d(0, 1, theta_C)
    return P_phi_12, P_phi_23, C_12, theta_phi, theta_C


# -------------- geometric regions and unitaries --------------

def build_geometric_regions(theta, n_regions=3):
    phase = np.mod(theta, 2 * np.pi)
    edges = np.linspace(0, 2*np.pi, n_regions+1)
    regions = []
    for k in range(n_regions):
        lo, hi = edges[k], edges[k+1]
        if k < n_regions - 1:
            idx = np.where((phase >= lo) & (phase < hi))[0]
        else:
            idx = np.where((phase >= lo) & (phase <= hi))[0]
        if len(idx) == 0:
            idx = np.array([k % len(theta)], dtype=int)
        regions.append(idx)
    return regions


def build_geometric_unitary(gen_vecs, region_list):
    cols = []
    for R in region_list:
        v = np.sum(gen_vecs[R, :], axis=0)
        norm = np.linalg.norm(v)
        if norm < 1e-14:
            v = np.array([1.0, 0.0, 0.0], dtype=complex)
            norm = 1.0
        cols.append(v / norm)
    U_geom = np.column_stack(cols)
    Uu, _, Vh = np.linalg.svd(U_geom)
    return Uu @ Vh


def build_sector_bases(P_phi_12, P_phi_23, C_12, U_geom, use_neutrino_dressing=True,
                       N_SOLAR=36, N_REACTOR=45):
    sector_bases = {}

    U_geom_u = U_geom["u"]
    U_geom_d = U_geom["d"]
    U_geom_e = U_geom["e"]
    U_geom_nu = U_geom["nu"]

    U_L_u = U_geom_u @ C_12.conj().T
    U_R_u = np.eye(3, dtype=complex)

    U_L_d = U_geom_d
    U_R_d = np.eye(3, dtype=complex)

    U_L_e = U_geom_e
    U_R_e = np.eye(3, dtype=complex)

    if use_neutrino_dressing:
        theta_solar = 2 * np.pi / N_SOLAR
        theta_reac = 2 * np.pi / N_REACTOR
        R_solar = rotation_3d(0, 1, theta_solar)
        R_reac = rotation_3d(0, 2, theta_reac)
        U_dress = (P_phi_23 @ R_solar) @ (P_phi_12 @ R_reac)
        U_L_nu = U_geom_nu @ U_dress
    else:
        U_L_nu = U_geom_nu

    U_R_nu = np.eye(3, dtype=complex)

    sector_bases["u"] = (U_L_u, U_R_u)
    sector_bases["d"] = (U_L_d, U_R_d)
    sector_bases["e"] = (U_L_e, U_R_e)
    sector_bases["nu"] = (U_L_nu, U_R_nu)

    return sector_bases


# -------------- Yukawas, mixing, observables, chi^2 --------------

def yukawa_from_F_and_UL(F_s, U_L, U_R):
    D = np.diag(F_s)
    return U_L @ D @ U_R.conj().T


def mixing_matrix(U_L_up, U_L_down):
    return U_L_up.conj().T @ U_L_down


def mixing_angles_from_U(U):
    s13 = abs(U[0, 2])
    s13 = min(max(s13, 0.0), 1.0)
    theta13 = np.arcsin(s13)
    c13 = np.cos(theta13)
    if abs(c13) < 1e-12:
        theta12 = 0.0
        theta23 = 0.0
    else:
        theta12 = np.arctan2(abs(U[0, 1]), abs(U[0, 0]))
        theta23 = np.arctan2(abs(U[1, 2]), abs(U[2, 2]))
    return theta12, theta23, theta13


TARGETS = {
    "mu_mt":     (2.2e-05, 0.5 * 2.2e-05),
    "mc_mt":     (7.5e-03, 0.5 * 7.5e-03),
    "md_mb":     (1.1e-03, 0.5 * 1.1e-03),
    "ms_mb":     (2.2e-02, 0.5 * 2.2e-02),
    "me_mt":     (2.9e-04, 0.5 * 2.9e-04),
    "mmu_mt":    (5.9e-02, 0.5 * 5.9e-02),
    "theta12_q": (0.227,   0.05 * 0.227),
    "theta23_q": (0.041,   0.5  * 0.041),
    "theta13_q": (0.0036,  0.5  * 0.0036),
    "theta12_l": (0.584,   0.1  * 0.584),
    "theta23_l": (0.785,   0.2  * 0.785),
    "theta13_l": (0.15,    0.2  * 0.15),
}


def compute_observables(mu_mt, mc_mt, md_mb, ms_mb, me_mt, mmu_mt,
                        theta12_q, theta23_q, theta13_q,
                        theta12_l, theta23_l, theta13_l):
    return {
        "mu_mt":     mu_mt,
        "mc_mt":     mc_mt,
        "md_mb":     md_mb,
        "ms_mb":     ms_mb,
        "me_mt":     me_mt,
        "mmu_mt":    mmu_mt,
        "theta12_q": theta12_q,
        "theta23_q": theta23_q,
        "theta13_q": theta13_q,
        "theta12_l": theta12_l,
        "theta23_l": theta23_l,
        "theta13_l": theta13_l,
    }


def chi2(obs, targets):
    chi2_val = 0.0
    details = []
    for k, v in obs.items():
        target, sigma = targets[k]
        if sigma <= 0:
            continue
        contrib = ((v - target) / sigma)**2
        chi2_val += contrib
        details.append((k, v, target, contrib))
    return chi2_val, details


# ================================
# Main driver
# ================================

#!/usr/bin/env python3
"""
sm_finite_triple.py

Minimal 1-generation SM-like finite spectral triple:
  A_F ≃ C ⊕ H ⊕ M_3(C)
  H_F: quarks + leptons + color + antiparticles
  D_F: Yukawa + LR mixing (no full Majorana yet)

Goal: provide a concrete (A_F, H_F, D_F, J_F, gamma_F) that you can test with
your ncg_tests harness and later refine toward the full Connes–Chamseddine SM.
"""

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

if __name__ == "__main__":
    main()

"""
RESULTS:
=== NCG test suite for SM-like 1gen finite triple ===
=== First-order condition ===
Max Frobenius norm over all pairs (a,b): 0.000e+00
Pairs with norm < 1.0e-12:
  (a=           I, b=           I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=           I, b=         H_j) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=           I, b=  color_diag) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         H_j, b=           I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         H_j, b=         H_j) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         H_j, b=  color_diag) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=  color_diag, b=           I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=  color_diag, b=         H_j) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=  color_diag, b=  color_diag) → ||[[D,a],J b J^-1]||_F = 0.000e+00

=== Zero-order condition ===
Max Frobenius norm over all pairs (a,b): 0.000e+00
All pairs satisfy [a, J b J^-1] ≈ 0 within eps=1.0e-12

=== Grading & reality tests ===
||{gamma_F, D_F}||_F       = 2.142e-04
max ||[gamma_F, a]||_F       = 0.000e+00
||J^2 - I||_F                = 0.000e+00
||J D_F J^-1 - D_F||_F       = 0.000e+00
||J D_F J^-1 + D_F||_F       = 2.142e-04
→ KO-sign: J D_F J^-1 = + D_F (J-even Dirac operator)
"""