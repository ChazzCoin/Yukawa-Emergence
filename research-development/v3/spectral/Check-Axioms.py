#!/usr/bin/env python3
"""
check_axioms.py

Diagnostic script for a truncated harmonic / divisor-based spectral triple.

What it does:
- Builds truncated geometric Dirac D_geom on modes n = -N,...,N
- Builds a simple 3x3 finite Dirac D_F (alignment-style)
- Forms the product Dirac: D = D_geom ⊗ I_F + gamma_geom ⊗ D_F
- Constructs some sample algebra generators A (geometric + finite)
- Checks operator norms of commutators [D, a]
- Computes eigenvalues of D, a zeta-function approximation, and spectral action

Dependencies:
    numpy
"""

import numpy as np


# =========================
# CONFIGURATION
# =========================

# Truncation parameter: geometric modes n = -N,...,N (dimension = 2N+1)
N_MODES = 20

# Internal (finite) dimension
N_FINITE = 3

# Alignment parameter for finite Dirac
KAPPA = 0.24

# Zeta-function parameters
ZETA_S_LIST = [2.0, 3.0]   # Re(s) > 1 for 1D-like spectrum
ZETA_EPS_CUTOFF = 1e-12    # Ignore eigenvalues with |lambda| < eps

# Spectral action parameters
LAMBDA_LIST = [5.0, 10.0]
# f(x) = exp(-x^2) is used below as the test function


# =========================
# GEOMETRIC PART
# =========================

def build_geom_dirac(N: int) -> np.ndarray:
    """
    Build truncated geometric Dirac D_geom on modes n = -N,...,N:
        D_geom |n> = n |n>
    Returned as a diagonal matrix of shape (2N+1, 2N+1).
    """
    n_vals = np.arange(-N, N + 1, dtype=float)
    return np.diag(n_vals)


def build_geom_gamma(N: int) -> np.ndarray:
    n_vals = np.arange(-N, N + 1, dtype=float)
    gamma_vals = np.sign(n_vals)
    # Ensure gamma^2 = I by setting gamma(0) = +1
    gamma_vals[N] = 1.0
    return np.diag(gamma_vals)



# =========================
# FINITE / INTERNAL PART
# =========================

def build_finite_dirac(kappa: float, n_finite: int = 3) -> np.ndarray:
    """
    Simple alignment-style finite Dirac D_F with entries kappa^{|i-j|}.
    Indices i,j = 0,...,n_finite-1.
    """
    D_F = np.zeros((n_finite, n_finite), dtype=float)
    for i in range(n_finite):
        for j in range(n_finite):
            D_F[i, j] = kappa ** abs(i - j)
    return D_F


def build_finite_gamma(n_finite: int = 3) -> np.ndarray:
    """
    Simple finite grading gamma_F.
    For demonstration, use diag(+1, -1, +1), but you can adjust as needed.
    """
    vals = np.ones(n_finite, dtype=float)
    if n_finite >= 2:
        vals[1] = -1.0
    return np.diag(vals)


def build_finite_J(n_finite: int = 3) -> np.ndarray:
    """
    Finite part of the real structure J_F.
    Here we take it to be the identity matrix; the anti-linearity will be
    "complex conjugation" outside of this script.
    """
    return np.eye(n_finite, dtype=complex)


# =========================
# PRODUCT TRIPLE
# =========================

def kron(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Convenience wrapper for Kronecker product."""
    return np.kron(a, b)


def build_product_dirac(N: int, kappa: float, n_finite: int) -> np.ndarray:
    """
    Build full truncated Dirac:
        D = D_geom ⊗ I_F + gamma_geom ⊗ D_F
    """
    D_geom = build_geom_dirac(N)
    gamma_geom = build_geom_gamma(N)
    D_F = build_finite_dirac(kappa, n_finite)

    I_F = np.eye(n_finite, dtype=float)
    D1 = kron(D_geom, I_F)
    D2 = kron(gamma_geom, D_F)
    return D1 + D2


def build_product_gamma(N: int, n_finite: int) -> np.ndarray:
    """Build full grading gamma = gamma_geom ⊗ gamma_F."""
    gamma_geom = build_geom_gamma(N)
    gamma_F = build_finite_gamma(n_finite)
    return kron(gamma_geom, gamma_F)


# =========================
# ALGEBRA GENERATORS
# =========================

def build_divisor_projector(N: int, d: int) -> np.ndarray:
    """
    Geometric algebra generator: projector onto modes with n divisible by d.
    On modes n = -N,...,N.

    Note: 0 is treated as divisible by any d; you can change this if desired.
    """
    n_vals = np.arange(-N, N + 1, dtype=int)
    mask = (n_vals % d == 0)
    proj = np.diag(mask.astype(float))
    return proj


def build_geom_algebra_generators(N: int) -> dict:
    """
    Build a small set of sample geometric algebra generators.

    Returns a dict mapping names -> matrices on H_geom.
    """
    I_geom = np.eye(2 * N + 1, dtype=float)
    generators = {"I_geom": I_geom}

    # Example: projectors onto modes divisible by 2, 3, 5
    for d in [2, 3, 5]:
        generators[f"P_div_{d}"] = build_divisor_projector(N, d)

    return generators


def build_finite_algebra_generators(n_finite: int) -> dict:
    """
    Build simple finite algebra generators: diagonal idempotents e_ii.
    A_F ~ C^n, represented as diagonal matrices.
    """
    gens = {}
    I_F = np.eye(n_finite, dtype=float)
    gens["I_F"] = I_F
    for i in range(n_finite):
        e = np.zeros((n_finite, n_finite), dtype=float)
        e[i, i] = 1.0
        gens[f"e_{i}"] = e
    return gens


def build_product_algebra_generators(N: int, n_finite: int) -> dict:
    """
    Build product algebra generators acting on H_geom ⊗ H_F:
        a_geom ⊗ I_F and I_geom ⊗ a_F
    """
    geom_gens = build_geom_algebra_generators(N)
    finite_gens = build_finite_algebra_generators(n_finite)

    I_geom = geom_gens["I_geom"]
    I_F = finite_gens["I_F"]

    product_gens = {}

    # Geometric generators tensored with identity on finite part
    for name, A_geom in geom_gens.items():
        product_gens[f"{name}⊗I_F"] = kron(A_geom, I_F)

    # Finite generators tensored with identity on geom part
    for name, A_F in finite_gens.items():
        product_gens[f"I_geom⊗{name}"] = kron(I_geom, A_F)

    return product_gens


# =========================
# COMMUTATOR NORMS
# =========================

def commutator(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """[A, B] = AB - BA"""
    return A @ B - B @ A


def op_norm(A: np.ndarray, ord: int = 2) -> float:
    """
    Operator norm of a matrix.
    ord=2 gives the spectral norm (largest singular value).
    """
    return np.linalg.norm(A, ord)


def check_commutator_norms(D: np.ndarray, A_gens: dict) -> dict:
    """
    For each generator a in A_gens, compute ||[D, a]||_2.
    Returns a dict name -> norm.
    """
    norms = {}
    for name, A in A_gens.items():
        C = commutator(D, A)
        norms[name] = op_norm(C, ord=2)
    return norms


# =========================
# SPECTRUM, ZETA, SPECTRAL ACTION
# =========================

def eigenvalues(D: np.ndarray) -> np.ndarray:
    """
    Compute eigenvalues of D (Hermitian assumed).
    """
    # eigh is for Hermitian, which D is (in this construction)
    vals, _ = np.linalg.eigh(D)
    return vals


def zeta_approx(D: np.ndarray, s: float, eps_cutoff: float = 1e-12) -> float:
    """
    Approximate zeta_D(s) = sum_{lambda != 0} |lambda|^{-s}
    on the truncated spectrum.
    """
    vals = eigenvalues(D)
    # Exclude very small |lambda| to avoid divergence / numerical issues
    mask = np.abs(vals) > eps_cutoff
    vals = np.abs(vals[mask])
    return np.sum(vals ** (-s))


def spectral_action(D: np.ndarray, Lambda: float) -> float:
    """
    Simple spectral action S(Λ) = Tr f(D/Λ) with f(x) = exp(-x^2).

    You can change f if desired.
    """
    vals = eigenvalues(D)
    x = vals / Lambda
    f_vals = np.exp(-x ** 2)
    return float(np.sum(f_vals))


# =========================
# MAIN CHECK ROUTINE
# =========================

def main() -> None:
    N = N_MODES
    n_finite = N_FINITE
    kappa = KAPPA

    print("=== Harmonic / Divisor-Based Spectral Triple Diagnostics ===")
    print(f"Truncation N          = {N}   (geom dimension = {2 * N + 1})")
    print(f"Finite dim            = {n_finite}")
    print(f"Kappa (alignment)     = {kappa}")
    print()

    # Build full Dirac and algebra generators
    D = build_product_dirac(N, kappa, n_finite)
    gamma = build_product_gamma(N, n_finite)
    A_gens = build_product_algebra_generators(N, n_finite)

    dimH = D.shape[0]
    print(f"Hilbert space dimension = {dimH}")
    print()

    # Basic Hermiticity checks
    herm_diff = np.linalg.norm(D - D.T.conj(), ord=2)
    gamma_sq_diff = np.linalg.norm(gamma @ gamma - np.eye(dimH), ord=2)
    anti_comm_diff = np.linalg.norm(gamma @ D + D @ gamma, ord=2)

    print("=== Basic operator checks ===")
    print(f"||D - D^†||_2                = {herm_diff:.3e} (should ~ 0)")
    print(f"||gamma^2 - I||_2           = {gamma_sq_diff:.3e} (should ~ 0)")
    print(f"||gamma D + D gamma||_2     = {anti_comm_diff:.3e} (should ~ 0 for even triple)")
    print()

    # Commutator norms
    print("=== Commutator norms ||[D, a]||_2 for sample generators ===")
    comm_norms = check_commutator_norms(D, A_gens)
    for name, norm in comm_norms.items():
        print(f"{name:20s} : {norm:.3e}")
    print()

    # Spectral analysis
    print("=== Spectrum of D (truncated) ===")
    vals = eigenvalues(D)
    print(f"Number of eigenvalues  = {len(vals)}")
    print("First 10 eigenvalues (sorted):")
    print(np.round(vals[:10], 6))
    print("Last 10 eigenvalues (sorted):")
    print(np.round(vals[-10:], 6))
    print()

    # Zeta approximation
    print("=== Zeta-function approximation (truncated) ===")
    for s in ZETA_S_LIST:
        zeta_val = zeta_approx(D, s, eps_cutoff=ZETA_EPS_CUTOFF)
        print(f"zeta_D({s}) ≈ {zeta_val:.6e} (with |lambda| > {ZETA_EPS_CUTOFF})")
    print()

    # Spectral action
    print("=== Spectral action S(Λ) = Tr exp(-(D/Λ)^2) ===")
    for L in LAMBDA_LIST:
        S_L = spectral_action(D, L)
        print(f"S({L}) ≈ {S_L:.6f}")
    print()

    print("Diagnostics complete.")


if __name__ == "__main__":
    main()

"""
=== Harmonic / Divisor-Based Spectral Triple Diagnostics ===
Truncation N          = 20   (geom dimension = 41)
Finite dim            = 3
Kappa (alignment)     = 0.24

Hilbert space dimension = 123

=== Basic operator checks ===
||D - D^†||_2                = 0.000e+00 (should ~ 0)
||gamma^2 - I||_2           = 0.000e+00 (should ~ 0)
||gamma D + D gamma||_2     = 4.212e+01 (should ~ 0 for even triple)

=== Commutator norms ||[D, a]||_2 for sample generators ===
I_geom⊗I_F           : 0.000e+00
P_div_2⊗I_F          : 0.000e+00
P_div_3⊗I_F          : 0.000e+00
P_div_5⊗I_F          : 0.000e+00
I_geom⊗e_0           : 2.468e-01
I_geom⊗e_1           : 3.394e-01
I_geom⊗e_2           : 2.468e-01

=== Spectrum of D (truncated) ===
Number of eigenvalues  = 123
First 10 eigenvalues (sorted):
[-21.369431 -20.9424   -20.688169 -20.369431 -19.9424   -19.688169
 -19.369431 -18.9424   -18.688169 -18.369431]
Last 10 eigenvalues (sorted):
[18.369431 18.688169 18.9424   19.369431 19.688169 19.9424   20.369431
 20.688169 20.9424   21.369431]

=== Zeta-function approximation (truncated) ===
zeta_D(2.0) ≈ 7.477203e+00 (with |lambda| > 1e-12)
zeta_D(3.0) ≈ 5.963309e+00 (with |lambda| > 1e-12)

=== Spectral action S(Λ) = Tr exp(-(D/Λ)^2) ===
S(5.0) ≈ 20.722721
S(10.0) ≈ 47.082763

Diagnostics complete.
"""