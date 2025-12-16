#!/usr/bin/env python3
"""
Universally honest spectral-resonance toy model.

Core idea:
- There is ONE parent space: a ring of N = 360 sites (C_360).
- Its Laplacian has eigenmodes labelled by k, with eigenvalues λ_k.
- We pick the triad of periods (60, 120, 180), which correspond to
  wavenumbers k = 360/60 = 6, 360/120 = 3, 360/180 = 2.
- A single universal kernel K is defined as a function of the spectrum: f(λ) = exp(-λ).
- Yukawa matrices are overlaps between:
    - a "left-handed" triad basis (3 modes),
    - and "right-handed" triad bases rotated by fixed, parameter-free 3×3 unitaries.
- All structure comes from:
    - the eigenvalues λ_k of the ring Laplacian (geometry),
    - the universal function f(λ) = exp(-λ),
    - and simple internal Z3 rotations (parameterless matrices).

No:
- No random Majorana matrices.
- No hand-chosen exponent tables like [4,2,0], [3,2,0], etc.
- No per-sector α scaling parameters.
- No RGEs here: this is a *static* spectral snapshot.

This script is a clean "spectrum → kernel → Yukawas → mixing" pipeline.
"""

import numpy as np

# ---------------------------------------------------------------------------
# 1. Geometry and spectrum of the parent space: C_360 ring
# ---------------------------------------------------------------------------

def ring_laplacian_eigenvalues(N: int) -> np.ndarray:
    """
    Eigenvalues of the Laplacian on a 1D ring (cycle graph) of N sites.

    For the ring (C_N) with nearest-neighbor coupling,
    the Laplacian eigenvalues are:

        λ_k = 2 - 2 cos(2π k / N),  k = 0,1,...,N-1

    This is a standard result in spectral graph theory and condensed matter.
    """
    k = np.arange(N)
    lambdas = 2.0 - 2.0 * np.cos(2.0 * np.pi * k / N)
    return lambdas


def dft_mode(N: int, k: int) -> np.ndarray:
    """
    Normalized complex plane-wave eigenmode on the ring:

        φ_k(n) = exp(2π i k n / N) / sqrt(N),  n = 0,...,N-1

    This satisfies L φ_k = λ_k φ_k, where L is the ring Laplacian.
    """
    n = np.arange(N)
    mode = np.exp(2j * np.pi * k * n / N) / np.sqrt(N)
    return mode


# ---------------------------------------------------------------------------
# 2. Triad selection: periods 60, 120, 180 ⇒ wavenumbers k = 6, 3, 2
# ---------------------------------------------------------------------------

def triad_wavenumbers(N: int):
    """
    For a ring of length N=360, a mode with *spatial period* P
    repeats every P sites, so its wavenumber is k = N / P.

    Our triad of periods: 60, 120, 180.

    That gives k = 6, 3, 2. These are divisors of 360 in wavenumber space.
    """
    periods = np.array([60, 120, 180])
    k_vals = N // periods  # integer division, works because 60,120,180 divide 360
    return periods, k_vals


def triad_spectral_data(N: int):
    """
    Return:
        periods  - [60, 120, 180]
        k_vals   - corresponding wavenumbers [6, 3, 2]
        lambdas  - corresponding Laplacian eigenvalues λ_k
    """
    lambdas = ring_laplacian_eigenvalues(N)
    periods, k_vals = triad_wavenumbers(N)
    triad_lambdas = lambdas[k_vals]
    return periods, k_vals, triad_lambdas


# ---------------------------------------------------------------------------
# 3. Universal kernel f(λ) = exp(-λ)
# ---------------------------------------------------------------------------

def universal_kernel_values(lambdas: np.ndarray) -> np.ndarray:
    """
    Universal scalar kernel f(λ) applied to eigenvalues.

    We choose here the simplest nontrivial positive function:

        f(λ) = exp(-λ)

    - No per-sector scale β.
    - No extra parameters.

    This is the discrete heat kernel at "time" t=1 in units where the Laplacian is dimensionless.
    """
    return np.exp(-lambdas)


# ---------------------------------------------------------------------------
# 4. Internal Z3 rotations: parameter-free mixing matrices
# ---------------------------------------------------------------------------

def unitary_F3() -> np.ndarray:
    """
    3×3 discrete Fourier transform matrix on Z3:

        (F3)_{jk} = 1/sqrt(3) * exp(2π i j k / 3), j,k=0,1,2.

    This is a canonical unitary; no continuous parameters.
    """
    j = np.arange(3)[:, None]  # 3×1
    k = np.arange(3)[None, :]  # 1×3
    omega = np.exp(2j * np.pi / 3.0)
    F = omega ** (j * k)
    F /= np.sqrt(3.0)
    return F


def real_Z3_rotation() -> np.ndarray:
    """
    A simple real orthogonal 3×3 matrix with a Z3-like structure.
    For example, a 'Hadamard-like' mixing in the (2,3) sector:

        R = [[1,       0,        0      ],
             [0,  1/sqrt(2), 1/sqrt(2)],
             [0,  1/sqrt(2), -1/sqrt(2)]]

    This is parameter-free and orthonormal.
    """
    R = np.array([
        [1.0, 0.0, 0.0],
        [0.0, 1.0 / np.sqrt(2.0), 1.0 / np.sqrt(2.0)],
        [0.0, 1.0 / np.sqrt(2.0), -1.0 / np.sqrt(2.0)]
    ])
    return R


# ---------------------------------------------------------------------------
# 5. Build Yukawa matrices in the triad subspace
# ---------------------------------------------------------------------------

def build_yukawa_matrices_from_triad(N: int):
    """
    Construct 3×3 Yukawa matrices for four sectors (u, d, e, ν)
    purely from:

    - The triad eigenvalues λ_k for periods 60,120,180.
    - The universal kernel f(λ) = exp(-λ).
    - Parameter-free internal unitaries on the 3-dimensional triad space.

    We work entirely in the triad eigenbasis, so we never need the full 360×360 matrices
    explicitly: everything is spectral.

    Setup:
    - Left-handed triad basis = identity in this 3D space (we label states by (k=6,3,2)).
    - Up-right basis      = identity (no extra mixing).
    - Down-right basis    = F3 (discrete Fourier on Z3).
    - Electron-right      = real rotation R.
    - Neutrino-right      = F3^† (conjugate transpose of F3).

    Yukawas:
        Yu  = diag(f(λ_k)) @ U_uR
        Yd  = diag(f(λ_k)) @ U_dR
        Ye  = diag(f(λ_k)) @ U_eR
        Ynu = diag(f(λ_k)) @ U_nuR

    Left-handed flavor/mass misalignment (CKM / PMNS analogues)
    come from diagonalizing these Yukawas.
    """
    # Triad spectral data
    periods, k_vals, triad_lambdas = triad_spectral_data(N)
    f_vals = universal_kernel_values(triad_lambdas)  # f(λ_k)

    # Put them in a diagonal matrix in the triad basis
    K_tri = np.diag(f_vals)  # 3×3

    # Internal unitaries (parameter-free)
    F3 = unitary_F3()
    R  = real_Z3_rotation()

    # Right-handed "flavor" bases in triad space
    U_uR  = np.eye(3, dtype=complex)   # up: same as triad basis
    U_dR  = F3                         # down: Z3 Fourier rotation
    U_eR  = R.astype(complex)          # charged leptons: real mixing
    U_nR  = F3.conj().T                # neutrinos: conjugate of F3

    # Yukawa matrices in the left-triad basis (3×3)
    # (Left mixing = identity; all nontrivial structure on the right)
    Yu  = K_tri @ U_uR
    Yd  = K_tri @ U_dR
    Ye  = K_tri @ U_eR
    Ynu = K_tri @ U_nR

    return {
        "periods": periods,
        "k_vals": k_vals,
        "lambdas": triad_lambdas,
        "kernel_values": f_vals,
        "Yu": Yu,
        "Yd": Yd,
        "Ye": Ye,
        "Ynu": Ynu,
    }


# ---------------------------------------------------------------------------
# 6. Diagonalization and mixing matrices (CKM / PMNS analogues)
# ---------------------------------------------------------------------------

def diagonalize_dirac(Y: np.ndarray):
    """
    Dirac-like Yukawa diagonalization via SVD:

        Y = U_L diag(s) U_R^†

    Returns:
        U_L, s, U_R
    where s are the singular values (≥ 0).
    """
    U_L, s_vals, U_Rh = np.linalg.svd(Y)
    U_R = U_Rh.conj().T
    return U_L, s_vals, U_R


def mixing_matrix(U_L_up: np.ndarray, U_L_down: np.ndarray) -> np.ndarray:
    """
    CKM/PMNS-like mixing matrix:

        V = U_L_up^† U_L_down
    """
    return U_L_up.conj().T @ U_L_down


def mixing_angles_from_U(U: np.ndarray):
    """
    Extract approximate (θ12, θ23, θ13) mixing angles from a unitary 3×3 matrix U
    using the standard parametrization (ignoring CP phase).

    Uses absolute values of matrix elements:

        s13 = |U_13|
        s12 = |U_12| / sqrt(1 - |U_13|^2)
        s23 = |U_23| / sqrt(1 - |U_13|^2)
    """
    U_abs = np.abs(U)
    s13 = U_abs[0, 2]
    c13 = np.sqrt(max(0.0, 1.0 - s13**2))
    if c13 < 1e-12:
        # pathological corner; just return zeros
        return 0.0, 0.0, np.pi / 2.0

    s12 = U_abs[0, 1] / c13
    s23 = U_abs[1, 2] / c13

    # Clamp numerical drift
    s12 = np.clip(s12, -1.0, 1.0)
    s23 = np.clip(s23, -1.0, 1.0)

    theta12 = np.arcsin(s12)
    theta23 = np.arcsin(s23)
    theta13 = np.arcsin(s13)

    return theta12, theta23, theta13


# ---------------------------------------------------------------------------
# 7. Main routine to show the emergent structure
# ---------------------------------------------------------------------------

def main():
    N = 360  # size of the parent ring (D_360 environment)

    data = build_yukawa_matrices_from_triad(N)
    periods = data["periods"]
    k_vals  = data["k_vals"]
    lambdas = data["lambdas"]
    f_vals  = data["kernel_values"]

    Yu  = data["Yu"]
    Yd  = data["Yd"]
    Ye  = data["Ye"]
    Ynu = data["Ynu"]

    print("=== Parent space: ring C_{} ===".format(N))
    print("Triad periods (in sites):", periods)
    print("Corresponding wavenumbers k:", k_vals)
    print("Laplacian eigenvalues λ_k for triad modes:", lambdas)
    print("Universal kernel values f(λ_k) = exp(-λ_k):", f_vals)
    print()

    # Diagonalize Yukawas
    Uu_L, su, Uu_R = diagonalize_dirac(Yu)
    Ud_L, sd, Ud_R = diagonalize_dirac(Yd)
    Ue_L, se, Ue_R = diagonalize_dirac(Ye)
    Un_L, sn, Un_R = diagonalize_dirac(Ynu)

    # Compute "CKM" (quark mixing) and "PMNS" (lepton mixing) analogues
    V_ckm  = mixing_matrix(Uu_L, Ud_L)
    U_pmns = mixing_matrix(Ue_L, Un_L)

    theta12_q, theta23_q, theta13_q = mixing_angles_from_U(V_ckm)
    theta12_l, theta23_l, theta13_l = mixing_angles_from_U(U_pmns)

    print("=== Yukawa singular values (up to overall scale) ===")
    print("Up-type (su):", su)
    print("Down-type (sd):", sd)
    print("Charged leptons (se):", se)
    print("Neutrino Dirac (sn):", sn)
    print()

    print("=== CKM-like mixing matrix (quarks) ===")
    print(V_ckm)
    print("Approx mixing angles (radians):")
    print("theta12_q ≈ {:.3f}, theta23_q ≈ {:.3f}, theta13_q ≈ {:.3f}".format(
        theta12_q, theta23_q, theta13_q))
    print()

    print("=== PMNS-like mixing matrix (leptons) ===")
    print(U_pmns)
    print("Approx mixing angles (radians):")
    print("theta12_l ≈ {:.3f}, theta23_l ≈ {:.3f}, theta13_l ≈ {:.3f}".format(
        theta12_l, theta23_l, theta13_l))
    print()

    print("NOTE:")
    print("- There were NO random matrices.")
    print("- There were NO sector-specific scaling parameters (α_u, α_d, etc.).")
    print("- There were NO hand-tuned exponent tables.")
    print("- EVERYTHING came from:")
    print("    * The geometry of the ring (C_{}) → λ_k spectrum.".format(N))
    print("    * The universal kernel f(λ) = exp(-λ).")
    print("    * A small set of parameter-free 3×3 unitaries (Z3 rotations).")
    print("- This is NOT tuned to match SM data; it's a clean demonstration of")
    print("  'kernel/Hilbert space → triad modes → Yukawas → mixing' with minimal hand input.")


if __name__ == "__main__":
    main()