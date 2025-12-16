#!/usr/bin/env python3
"""
check_product_triple.py

Combine:
  - Harmonic geometric Dirac D_geom on modes n = -N,...,N
  - Finite Dirac D_F from emergent-9.py (3-gen internal toy triple)

into a product Dirac:
    D = D_geom ⊗ I_F + I_geom ⊗ D_F

and run first-order / zero-order tests for the full (odd) product triple.

Dependencies:
    numpy
    emergent-9.py sitting on the PYTHONPATH (same directory is fine)

This is an ODD spectral triple test: no grading γ is used at the product level.
"""

import numpy as np

# Import your internal triple machinery
import _emergent_9 as em  # rename the file emergent-9.py -> emergent_9.py or adjust import


# =========================
# CONFIGURATION
# =========================

N_MODES = 20        # modes n = -N,...,N  → dim(H_geom) = 2N+1
ZETA_EPS = 1e-8     # cutoff for zeta (if you want it later)
EPS_FIRST = 1e-12   # tolerance for first-order condition
EPS_ZERO = 1e-12    # tolerance for zero-order condition


# =========================
# GEOMETRIC PART
# =========================

def build_geom_dirac(N: int) -> np.ndarray:
    """
    Truncated geometric Dirac on modes n = -N,...,N:
        D_geom |n> = n |n>
    """
    n_vals = np.arange(-N, N + 1, dtype=float)
    return np.diag(n_vals)


def build_geom_algebra_generators(N: int) -> dict:
    """
    Simple geometric algebra generators on H_geom:
      - I_geom
      - P_div_d: projectors onto modes divisible by d
    These all commute with D_geom by construction (diagonal in same basis).
    """
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


# =========================
# PRODUCT TRIPLE
# =========================

def kron(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.kron(a, b)


def build_product_dirac(D_geom: np.ndarray, D_F: np.ndarray) -> np.ndarray:
    """
    D = D_geom ⊗ I_F + I_geom ⊗ D_F
    """
    dim_geom = D_geom.shape[0]
    dimF = D_F.shape[0]

    I_geom = np.eye(dim_geom, dtype=complex)
    I_F = np.eye(dimF, dtype=complex)

    D1 = kron(D_geom, I_F)
    D2 = kron(I_geom, D_F)
    return D1 + D2


def build_product_algebra(N: int) -> tuple[list[np.ndarray], list[str]]:
    """
    Build a small basis of product algebra generators acting on H_geom ⊗ H_F:

        - a_geom ⊗ I_F   for a_geom in A_geom
        - I_geom ⊗ a_F   for a_F in A_F

    where A_F is the internal algebra from emergent_9.build_internal_algebra_ops().
    """
    # Geometric generators
    geom_gens = build_geom_algebra_generators(N)
    I_geom = geom_gens["I_geom"]

    # Internal generators
    ops_F, labels_F = em.build_internal_algebra_ops()
    dimF = ops_F[0].shape[0]
    I_F = np.eye(dimF, dtype=complex)

    ops_prod: list[np.ndarray] = []
    labels_prod: list[str] = []

    # Geometric part tensored with identity on finite
    for name, A_geom in geom_gens.items():
        ops_prod.append(kron(A_geom, I_F))
        labels_prod.append(f"{name}⊗I_F")

    # Finite part tensored with identity on geometric
    for A_F, lab in zip(ops_F, labels_F):
        ops_prod.append(kron(I_geom, A_F))
        labels_prod.append(f"I_geom⊗{lab}")

    return ops_prod, labels_prod


# =========================
# PRODUCT J (REAL STRUCTURE)
# =========================

def build_swap_LR_full(dim_geom: int, dim_left_F: int) -> np.ndarray:
    """
    Build the product swap S_prod implementing J on H_geom ⊗ H_F,
    assuming J acts trivially on geometry and as LR-swap on internal space:

        J M J^{-1} = S_prod · M^* · S_prod^T

    where S_prod = I_geom ⊗ S_F, and S_F swaps L/R blocks of size dim_left_F.
    """
    S_F = em.build_swap_LR(dim_left_F)
    I_geom = np.eye(dim_geom, dtype=complex)
    return kron(I_geom, S_F)


def J_action(S_prod: np.ndarray, M: np.ndarray) -> np.ndarray:
    """
    J M J^{-1} implemented as S_prod · M^* · S_prod^T.
    """
    return S_prod @ M.conj() @ S_prod.T


# =========================
# FIRST- & ZERO-ORDER TESTS (PRODUCT)
# =========================

def test_first_order_condition_product(
    D: np.ndarray,
    ops: list[np.ndarray],
    labels: list[str],
    S_prod: np.ndarray,
    eps: float = 1e-12,
) -> None:
    """
    First-order condition for the product triple:

        [[D, a], J b J^{-1}] = 0   for all a,b in A.

    Here J is implemented via S_prod and complex conjugation.
    """
    n = D.shape[0]
    assert D.shape == (n, n)
    print("=== First-order condition test (product triple) ===")

    max_norm = 0.0
    good_pairs: list[tuple[str, str, float]] = []

    for i, a in enumerate(ops):
        Da = D @ a - a @ D
        for j, b in enumerate(ops):
            b_tilde = J_action(S_prod, b)
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
            print(f"  (a={la:>20s}, b={lb:>20s}) → ||[[D,a],J b J^-1]||_F = {nrm:.3e}")
    else:
        print(f"No pairs with norm < {eps:.1e}")
    print()


def test_zero_order_condition_product(
    ops: list[np.ndarray],
    labels: list[str],
    S_prod: np.ndarray,
    eps: float = 1e-12,
) -> None:
    """
    Zero-order condition for the product triple:

        [a, J b J^{-1}] = 0   for all a,b in A.
    """
    n = ops[0].shape[0]
    print("=== Zero-order condition test (product triple) ===")
    max_norm = 0.0
    bad_pairs: list[tuple[str, str, float]] = []

    for i, a in enumerate(ops):
        for j, b in enumerate(ops):
            b_tilde = J_action(S_prod, b)
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
            print(f"  (a={la:>20s}, b={lb:>20s}) → ||[a, J b J^-1]||_F = {nrm:.3e}")
    else:
        print(f"All pairs satisfy [a, J b J^-1] ≈ 0 within eps={eps:.1e}")
    print()


# =========================
# MAIN DRIVER
# =========================

def main() -> None:
    N = N_MODES

    print("=== Product Spectral Triple Diagnostics ===")
    print(f"Truncation N          = {N}   (geom dimension = {2*N+1})")

    # 1) Geometric Dirac
    D_geom = build_geom_dirac(N)
    dim_geom = D_geom.shape[0]

    # 2) Finite Dirac from emergent-9: build_internal_DF_from_Y(run_emergent_alignment())
    print("\n--- Running emergent alignment to get Yukawas for D_F ---")
    align = em.run_emergent_alignment()
    Y_u = align["Y_u"]
    Y_d = align["Y_d"]
    Y_e = align["Y_e"]
    Y_nu = align["Y_nu"]

    D_F = em.build_internal_DF_from_Y(Y_u, Y_d, Y_e, Y_nu)
    dimF = D_F.shape[0]
    print(f"Finite internal dim(H_F) = {dimF}")

    # 3) Build product Dirac
    D = build_product_dirac(D_geom, D_F)
    dimH = D.shape[0]
    print(f"Total Hilbert space dim(H) = dim_geom * dim_F = {dimH}")
    print()

    # Basic Hermiticity check
    herm_norm = np.linalg.norm(D - D.T.conj(), ord=2)
    print("=== Basic operator check ===")
    print(f"||D - D^†||_2 = {herm_norm:.3e}")
    print()

    # 4) Product algebra
    ops_prod, labels_prod = build_product_algebra(N)

    # 5) Product J (via swap on finite sector)
    dpc = em.dim_per_chirality()  # size of one chirality in finite space
    S_prod = build_swap_LR_full(dim_geom, dpc)

    # 6) First-order and zero-order tests for the product triple
    test_first_order_condition_product(D, ops_prod, labels_prod, S_prod, eps=EPS_FIRST)
    test_zero_order_condition_product(ops_prod, labels_prod, S_prod, eps=EPS_ZERO)

    # (Optionally, you could also add zeta / spectral-action diagnostics here
    #  for D, but the primary goal is first-/zero-order.)

    print("Product triple tests complete.")


if __name__ == "__main__":
    main()

"""
=== Product Spectral Triple Diagnostics ===
Truncation N          = 20   (geom dimension = 41)

--- Running emergent alignment to get Yukawas for D_F ---
Finite internal dim(H_F) = 48
Total Hilbert space dim(H) = dim_geom * dim_F = 1968

=== Basic operator check ===
||D - D^†||_2 = 0.000e+00

=== First-order condition test (product triple) ===
Max Frobenius norm over all pairs (a,b): 0.000e+00
Pairs with norm < 1.0e-12:
  (a=          I_geom⊗I_F, b=          I_geom⊗I_F) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=          I_geom⊗I_F, b=         P_div_2⊗I_F) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=          I_geom⊗I_F, b=         P_div_3⊗I_F) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=          I_geom⊗I_F, b=         P_div_5⊗I_F) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=          I_geom⊗I_F, b=            I_geom⊗I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=          I_geom⊗I_F, b=     I_geom⊗Q_sector) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=          I_geom⊗I_F, b=   I_geom⊗P_sector_u) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=          I_geom⊗I_F, b=   I_geom⊗P_sector_d) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=          I_geom⊗I_F, b=   I_geom⊗P_sector_e) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=          I_geom⊗I_F, b=  I_geom⊗P_sector_nu) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         P_div_2⊗I_F, b=          I_geom⊗I_F) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         P_div_2⊗I_F, b=         P_div_2⊗I_F) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         P_div_2⊗I_F, b=         P_div_3⊗I_F) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         P_div_2⊗I_F, b=         P_div_5⊗I_F) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         P_div_2⊗I_F, b=            I_geom⊗I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         P_div_2⊗I_F, b=     I_geom⊗Q_sector) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         P_div_2⊗I_F, b=   I_geom⊗P_sector_u) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         P_div_2⊗I_F, b=   I_geom⊗P_sector_d) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         P_div_2⊗I_F, b=   I_geom⊗P_sector_e) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         P_div_2⊗I_F, b=  I_geom⊗P_sector_nu) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         P_div_3⊗I_F, b=          I_geom⊗I_F) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         P_div_3⊗I_F, b=         P_div_2⊗I_F) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         P_div_3⊗I_F, b=         P_div_3⊗I_F) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         P_div_3⊗I_F, b=         P_div_5⊗I_F) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         P_div_3⊗I_F, b=            I_geom⊗I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         P_div_3⊗I_F, b=     I_geom⊗Q_sector) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         P_div_3⊗I_F, b=   I_geom⊗P_sector_u) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         P_div_3⊗I_F, b=   I_geom⊗P_sector_d) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         P_div_3⊗I_F, b=   I_geom⊗P_sector_e) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         P_div_3⊗I_F, b=  I_geom⊗P_sector_nu) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         P_div_5⊗I_F, b=          I_geom⊗I_F) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         P_div_5⊗I_F, b=         P_div_2⊗I_F) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         P_div_5⊗I_F, b=         P_div_3⊗I_F) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         P_div_5⊗I_F, b=         P_div_5⊗I_F) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         P_div_5⊗I_F, b=            I_geom⊗I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         P_div_5⊗I_F, b=     I_geom⊗Q_sector) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         P_div_5⊗I_F, b=   I_geom⊗P_sector_u) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         P_div_5⊗I_F, b=   I_geom⊗P_sector_d) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         P_div_5⊗I_F, b=   I_geom⊗P_sector_e) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         P_div_5⊗I_F, b=  I_geom⊗P_sector_nu) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=            I_geom⊗I, b=          I_geom⊗I_F) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=            I_geom⊗I, b=         P_div_2⊗I_F) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=            I_geom⊗I, b=         P_div_3⊗I_F) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=            I_geom⊗I, b=         P_div_5⊗I_F) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=            I_geom⊗I, b=            I_geom⊗I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=            I_geom⊗I, b=     I_geom⊗Q_sector) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=            I_geom⊗I, b=   I_geom⊗P_sector_u) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=            I_geom⊗I, b=   I_geom⊗P_sector_d) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=            I_geom⊗I, b=   I_geom⊗P_sector_e) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=            I_geom⊗I, b=  I_geom⊗P_sector_nu) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=     I_geom⊗Q_sector, b=          I_geom⊗I_F) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=     I_geom⊗Q_sector, b=         P_div_2⊗I_F) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=     I_geom⊗Q_sector, b=         P_div_3⊗I_F) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=     I_geom⊗Q_sector, b=         P_div_5⊗I_F) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=     I_geom⊗Q_sector, b=            I_geom⊗I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=     I_geom⊗Q_sector, b=     I_geom⊗Q_sector) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=     I_geom⊗Q_sector, b=   I_geom⊗P_sector_u) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=     I_geom⊗Q_sector, b=   I_geom⊗P_sector_d) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=     I_geom⊗Q_sector, b=   I_geom⊗P_sector_e) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=     I_geom⊗Q_sector, b=  I_geom⊗P_sector_nu) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=   I_geom⊗P_sector_u, b=          I_geom⊗I_F) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=   I_geom⊗P_sector_u, b=         P_div_2⊗I_F) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=   I_geom⊗P_sector_u, b=         P_div_3⊗I_F) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=   I_geom⊗P_sector_u, b=         P_div_5⊗I_F) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=   I_geom⊗P_sector_u, b=            I_geom⊗I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=   I_geom⊗P_sector_u, b=     I_geom⊗Q_sector) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=   I_geom⊗P_sector_u, b=   I_geom⊗P_sector_u) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=   I_geom⊗P_sector_u, b=   I_geom⊗P_sector_d) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=   I_geom⊗P_sector_u, b=   I_geom⊗P_sector_e) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=   I_geom⊗P_sector_u, b=  I_geom⊗P_sector_nu) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=   I_geom⊗P_sector_d, b=          I_geom⊗I_F) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=   I_geom⊗P_sector_d, b=         P_div_2⊗I_F) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=   I_geom⊗P_sector_d, b=         P_div_3⊗I_F) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=   I_geom⊗P_sector_d, b=         P_div_5⊗I_F) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=   I_geom⊗P_sector_d, b=            I_geom⊗I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=   I_geom⊗P_sector_d, b=     I_geom⊗Q_sector) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=   I_geom⊗P_sector_d, b=   I_geom⊗P_sector_u) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=   I_geom⊗P_sector_d, b=   I_geom⊗P_sector_d) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=   I_geom⊗P_sector_d, b=   I_geom⊗P_sector_e) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=   I_geom⊗P_sector_d, b=  I_geom⊗P_sector_nu) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=   I_geom⊗P_sector_e, b=          I_geom⊗I_F) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=   I_geom⊗P_sector_e, b=         P_div_2⊗I_F) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=   I_geom⊗P_sector_e, b=         P_div_3⊗I_F) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=   I_geom⊗P_sector_e, b=         P_div_5⊗I_F) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=   I_geom⊗P_sector_e, b=            I_geom⊗I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=   I_geom⊗P_sector_e, b=     I_geom⊗Q_sector) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=   I_geom⊗P_sector_e, b=   I_geom⊗P_sector_u) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=   I_geom⊗P_sector_e, b=   I_geom⊗P_sector_d) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=   I_geom⊗P_sector_e, b=   I_geom⊗P_sector_e) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=   I_geom⊗P_sector_e, b=  I_geom⊗P_sector_nu) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=  I_geom⊗P_sector_nu, b=          I_geom⊗I_F) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=  I_geom⊗P_sector_nu, b=         P_div_2⊗I_F) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=  I_geom⊗P_sector_nu, b=         P_div_3⊗I_F) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=  I_geom⊗P_sector_nu, b=         P_div_5⊗I_F) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=  I_geom⊗P_sector_nu, b=            I_geom⊗I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=  I_geom⊗P_sector_nu, b=     I_geom⊗Q_sector) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=  I_geom⊗P_sector_nu, b=   I_geom⊗P_sector_u) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=  I_geom⊗P_sector_nu, b=   I_geom⊗P_sector_d) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=  I_geom⊗P_sector_nu, b=   I_geom⊗P_sector_e) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=  I_geom⊗P_sector_nu, b=  I_geom⊗P_sector_nu) → ||[[D,a],J b J^-1]||_F = 0.000e+00

=== Zero-order condition test (product triple) ===
Max Frobenius norm over all pairs (a,b): 0.000e+00
All pairs satisfy [a, J b J^-1] ≈ 0 within eps=1.0e-12

Product triple tests complete.

"""