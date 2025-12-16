#!/usr/bin/env python3
"""
check_product_triple.py

Product triple:
  - Geometric: harmonic Dirac D_geom on modes n = -N,...,N
  - Finite:    D_F from emergent_9 (emergent alignment finite triple)

D = D_geom ⊗ I_F + I_geom ⊗ D_F

We check:
  - Hermiticity
  - First-order condition
  - Zero-order condition
  - Zeta-function approximation
  - Spectral action scaling

Triple is treated as ODD (no grading γ on the product).
"""

import numpy as np
import _emergent_9 as em  # ensure file is named emergent_9.py


# =========================
# CONFIGURATION
# =========================

N_MODES   = 20       # n = -N,...,N  → dim(H_geom) = 2N + 1
EPS_FIRST = 1e-12    # tolerance for first-order
EPS_ZERO  = 1e-12    # tolerance for zero-order

# Zeta / spectral action settings
ZETA_S_LIST      = [2.0, 3.0]        # Re(s) > 1
ZETA_EPS_CUTOFFS = [1e-1, 1e-2, 1e-3]  # IR cutoffs to probe UV vs IR
LAMBDA_LIST      = [5.0, 10.0, 20.0]  # spectral-action scales


# =========================
# GEOMETRIC / PRODUCT PART
# =========================

def build_geom_dirac(N: int) -> np.ndarray:
    n_vals = np.arange(-N, N + 1, dtype=float)
    return np.diag(n_vals)


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


def kron(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.kron(a, b)


def build_product_dirac(D_geom: np.ndarray, D_F: np.ndarray) -> np.ndarray:
    dim_geom = D_geom.shape[0]
    dimF     = D_F.shape[0]

    I_geom = np.eye(dim_geom, dtype=complex)
    I_F    = np.eye(dimF, dtype=complex)

    return kron(D_geom, I_F) + kron(I_geom, D_F)


def build_product_algebra(N: int) -> tuple[list[np.ndarray], list[str]]:
    geom_gens = build_geom_algebra_generators(N)
    I_geom    = geom_gens["I_geom"]

    ops_F, labels_F = em.build_internal_algebra_ops()
    dimF = ops_F[0].shape[0]
    I_F  = np.eye(dimF, dtype=complex)

    ops_prod:   list[np.ndarray] = []
    labels_prod: list[str]       = []

    for name, A_geom in geom_gens.items():
        ops_prod.append(kron(A_geom, I_F))
        labels_prod.append(f"{name}⊗I_F")

    for A_F, lab in zip(ops_F, labels_F):
        ops_prod.append(kron(I_geom, A_F))
        labels_prod.append(f"I_geom⊗{lab}")

    return ops_prod, labels_prod


# =========================
# REAL STRUCTURE J
# =========================

def build_swap_LR_full(dim_geom: int, dim_left_F: int) -> np.ndarray:
    S_F   = em.build_swap_LR(dim_left_F)
    I_geom = np.eye(dim_geom, dtype=complex)
    return kron(I_geom, S_F)


def J_action(S_prod: np.ndarray, M: np.ndarray) -> np.ndarray:
    return S_prod @ M.conj() @ S_prod.T


# =========================
# FIRST- & ZERO-ORDER TESTS
# =========================

def test_first_order_condition_product(
    D: np.ndarray,
    ops: list[np.ndarray],
    labels: list[str],
    S_prod: np.ndarray,
    eps: float = 1e-12,
) -> None:
    print("=== First-order condition test (product triple) ===")
    max_norm = 0.0

    for i, a in enumerate(ops):
        Da = D @ a - a @ D
        for j, b in enumerate(ops):
            b_tilde = J_action(S_prod, b)
            comm2   = Da @ b_tilde - b_tilde @ Da
            norm    = np.linalg.norm(comm2, ord="fro")
            if norm > max_norm:
                max_norm = norm

    print(f"Max Frobenius norm over all pairs (a,b): {max_norm:.3e}\n")


def test_zero_order_condition_product(
    ops: list[np.ndarray],
    labels: list[str],
    S_prod: np.ndarray,
    eps: float = 1e-12,
) -> None:
    print("=== Zero-order condition test (product triple) ===")
    max_norm  = 0.0
    bad_pairs: list[tuple[str, str, float]] = []

    for i, a in enumerate(ops):
        for j, b in enumerate(ops):
            b_tilde = J_action(S_prod, b)
            comm    = a @ b_tilde - b_tilde @ a
            norm    = np.linalg.norm(comm, ord="fro")
            if norm > max_norm:
                max_norm = norm
            if norm > eps:
                bad_pairs.append((labels[i], labels[j], norm))

    print(f"Max Frobenius norm over all pairs (a,b): {max_norm:.3e}")
    if bad_pairs:
        print("Pairs with significant violation (> eps):")
        for la, lb, nrm in bad_pairs:
            print(f"  (a={la:>20s}, b={lb:>20s}) → ||[a, J b J^-1]||_F = {nrm:.3e}")
    else:
        print(f"All pairs satisfy [a, J b J^-1] ≈ 0 within eps={eps:.1e}")
    print()


# =========================
# ZETA & SPECTRAL ACTION
# =========================

def eigenvalues(D: np.ndarray) -> np.ndarray:
    vals, _ = np.linalg.eigh(D)
    return vals


def zeta_approx(D: np.ndarray, s: float, eps_cutoff: float) -> float:
    """
    zeta_D(s) ~ sum_{|λ|>eps} |λ|^{-s} on the truncated spectrum.
    """
    vals = eigenvalues(D)
    mask = np.abs(vals) > eps_cutoff
    vals = np.abs(vals[mask])
    return float(np.sum(vals ** (-s)))


def spectral_action(D: np.ndarray, Lambda: float) -> float:
    """
    S(Λ) = Tr exp(-(D/Λ)^2).
    """
    vals = eigenvalues(D)
    x = vals / Lambda
    return float(np.sum(np.exp(-x**2)))


# =========================
# MAIN
# =========================

def main() -> None:
    N = N_MODES

    print("=== Product Spectral Triple Diagnostics ===")
    print(f"Truncation N          = {N}   (geom dimension = {2*N+1})")

    # 1) Geometric Dirac
    D_geom   = build_geom_dirac(N)
    dim_geom = D_geom.shape[0]

    # 2) Finite Dirac from emergent_9
    print("\n--- Running emergent alignment to get Yukawas for D_F ---")
    align = em.run_emergent_alignment()
    Y_u, Y_d = align["Y_u"], align["Y_d"]
    Y_e, Y_nu = align["Y_e"], align["Y_nu"]

    D_F = em.build_internal_DF_from_Y(Y_u, Y_d, Y_e, Y_nu)
    dimF = D_F.shape[0]
    print(f"Finite internal dim(H_F) = {dimF}")

    # 3) Product Dirac
    D    = build_product_dirac(D_geom, D_F)
    dimH = D.shape[0]
    print(f"Total Hilbert space dim(H) = {dimH}\n")

    # Basic Hermiticity
    herm_norm = np.linalg.norm(D - D.T.conj(), ord=2)
    print("=== Basic operator check ===")
    print(f"||D - D^†||_2 = {herm_norm:.3e}\n")

    # 4) Product algebra
    ops_prod, labels_prod = build_product_algebra(N)

    # 5) Product J
    dpc    = em.dim_per_chirality()
    S_prod = build_swap_LR_full(dim_geom, dpc)

    # 6) First- and zero-order tests
    test_first_order_condition_product(D, ops_prod, labels_prod, S_prod, eps=EPS_FIRST)
    test_zero_order_condition_product(ops_prod, labels_prod, S_prod, eps=EPS_ZERO)

    # 7) Zeta & spectral action diagnostics
    print("=== Zeta-function approximation for full D ===")
    for eps_cut in ZETA_EPS_CUTOFFS:
        for s in ZETA_S_LIST:
            z = zeta_approx(D, s, eps_cutoff=eps_cut)
            print(f"eps={eps_cut:>5.0e}, s={s:.1f}: zeta_D(s) ≈ {z:.6e}")
    print()

    print("=== Spectral action S(Λ) = Tr exp(-(D/Λ)^2) ===")
    for Lam in LAMBDA_LIST:
        S_L = spectral_action(D, Lam)
        print(f"Λ={Lam:>5.1f} : S(Λ) ≈ {S_L:.6f}")
    print()

    # Optional: show smallest eigenvalues to see IR structure
    vals = eigenvalues(D)
    abs_vals = np.abs(vals)
    print("=== Smallest |λ| for full D ===")
    print("Min |λ| =", abs_vals.min())
    print("10 smallest |λ|:", np.sort(abs_vals)[:10])


if __name__ == "__main__":
    main()

"""
=== Product Spectral Triple Diagnostics ===
Truncation N          = 20   (geom dimension = 41)

--- Running emergent alignment to get Yukawas for D_F ---
Finite internal dim(H_F) = 48
Total Hilbert space dim(H) = 1968

=== Basic operator check ===
||D - D^†||_2 = 0.000e+00

=== First-order condition test (product triple) ===
Max Frobenius norm over all pairs (a,b): 0.000e+00

=== Zero-order condition test (product triple) ===
Max Frobenius norm over all pairs (a,b): 0.000e+00
All pairs satisfy [a, J b J^-1] ≈ 0 within eps=1.0e-12

=== Zeta-function approximation for full D ===
eps=1e-01, s=2.0: zeta_D(s) ≈ 1.532689e+02
eps=1e-01, s=3.0: zeta_D(s) ≈ 1.153549e+02
eps=1e-02, s=2.0: zeta_D(s) ≈ 6.922043e+03
eps=1e-02, s=3.0: zeta_D(s) ≈ 3.418311e+05
eps=1e-03, s=2.0: zeta_D(s) ≈ 5.097497e+04
eps=1e-03, s=3.0: zeta_D(s) ≈ 6.879866e+06

=== Spectral action S(Λ) = Tr exp(-(D/Λ)^2) ===
Λ=  5.0 : S(Λ) ≈ 425.388922
Λ= 10.0 : S(Λ) ≈ 847.618739
Λ= 20.0 : S(Λ) ≈ 1451.266015

=== Smallest |λ| for full D ===
Min |λ| = 0.0
10 smallest |λ|: [0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]

"""