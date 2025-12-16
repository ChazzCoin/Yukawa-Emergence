import numpy as np
import math

"""
Resonant Spectral Triple Flavor Toy
-----------------------------------

This script is a "production-level" toy model that combines:

- Spectral-triple style thinking (à la Connes) for the finite/internal space.
- A cyclic "base-360" flavor symmetry implemented as a unitary operator R of order 360.
- A universal spectral kernel F acting on the eigenvalues of R.
- Minimal, discrete left/right flavor rotations to generate nontrivial Yukawa textures
  and mixing matrices, without arbitrary random matrices or sector-specific exponents.

It is NOT a realistic SM fit. The goal is to have a clean, axioms-first,
operator-driven construction where flavor hierarchies and mixing are emergent
from a small algebra of operators, not hand-tuned numbers.
"""

# ============================================================
# 1. Base-360 cyclic operator on generation space
# ============================================================

# We work in a 3-dimensional "generation space" H_gen spanned by |1>, |2>, |3>.
# We introduce a unitary "rotation" operator R such that:
#   R^360 = I
# and its eigenvalues are 360-th roots of unity e^{ 2πi k_j / 360 }.
#
# These three eigenvalues encode the three generations via characters of Z_360.
# We pick k_j corresponding to periods (360 / gcd(360, k_j)) = {60, 120, 180}
# as a nod to the 60–120–180 triad:
#
#   k_1 = 360 / 60  = 6
#   k_2 = 360 / 120 = 3
#   k_3 = 360 / 180 = 2

BASE_ANGLE = 2.0 * math.pi / 360.0
k_vals = np.array([6, 3, 2], dtype=int)  # triadic "frequencies" in base-360

# Eigenvalues (characters) of R on the 3 generations
chi = np.exp(1j * BASE_ANGLE * k_vals)

# In the eigenbasis of R, R is diagonal:
R_diag = np.diag(chi)


# ============================================================
# 2. Universal spectral kernel F(R)
# ============================================================

# We now define a universal kernel F, which is a function of R (and thus of its spectrum).
# In spectral-triple language, this is analogous to using a function of the Dirac/Laplacian
# operator to build a "propagator" or "heat kernel".
#
# Here we want F to:
#   - be the same for all sectors (up, down, leptons, neutrinos),
#   - respect the base-360 rotation symmetry,
#   - introduce a mild hierarchy between the different characters chi_j.
#
# A simple choice is:
#   F(chi_j) = exp( -lambda * ( 1 - Re(chi_j) ) )
#
# This is a discrete analogue of a heat kernel on the circle: it damps modes based
# on how far their phase is from 1.
#
# We fix lambda = 1.0 here to avoid any continuous tuning; a user can change it if desired.

def spectral_kernel(chi_vals, lam=1.0):
    chi_real = np.real(chi_vals)
    return np.exp(-lam * (1.0 - chi_real))

F_vals = spectral_kernel(chi)            # shape (3,)
F_diag = np.diag(F_vals)                 # 3×3 diagonal matrix of spectral weights


# ============================================================
# 3. Left/right flavor bases and Yukawa construction
# ============================================================

# In Connes-style noncommutative geometry, the finite Dirac operator D_F encodes
# Yukawa couplings and mass matrices. Here we mimic this structure as follows:
#
# For each sector s in {u, d, e, nu} we define:
#   - A left-handed generation basis (3×3 unitary matrix U_L^s).
#   - A right-handed generation basis (3×3 unitary matrix U_R^s).
#
# We then define the sector's Yukawa matrix in the R-eigenbasis as:
#   Y_s = U_L^s†  F_diag  U_R^s
#
# where F_diag is the same universal spectral kernel for all sectors.
#
# All flavor structure (hierarchy + mixing) thus comes from the choice of the
# LEFT/RIGHT unitaries, which we restrict to simple, parameter-free, group-like
# matrices (discrete Fourier transforms, permutations, and fixed real rotations).


def unitary_F3():
    """3×3 discrete Fourier transform on Z_3."""
    omega = np.exp(2j * math.pi / 3.0)
    j = np.arange(3)[:, None]
    k = np.arange(3)[None, :]
    F = omega ** (j * k)
    F /= math.sqrt(3.0)
    return F


def real_rotation_23(theta):
    """
    3×3 real rotation in the (2,3) subspace by angle theta.
    The first generation is left invariant.
    """
    c = math.cos(theta)
    s = math.sin(theta)
    R = np.array([
        [1.0, 0.0, 0.0],
        [0.0, c,   s  ],
        [0.0, -s,  c  ],
    ], dtype=float)
    return R


# Build a small library of "axiomatically allowed" unitaries: F3, its powers,
# and fixed rotations/permutations. In a more elaborate version these would be
# tied to discrete flavor groups (A_4, S_4, etc.). Here we keep it minimal.

F3 = unitary_F3()
I3 = np.eye(3, dtype=complex)

# A fixed Cabibbo-like twist in (2,3), chosen at a nice angle (π/6) without fitting.
R23_30deg = real_rotation_23(math.pi / 6.0).astype(complex)

# A simple permutation matrix that swaps generations 2 and 3.
P_23 = np.array([
    [1.0, 0.0, 0.0],
    [0.0, 0.0, 1.0],
    [0.0, 1.0, 0.0],
], dtype=complex)


# Define sector-specific left/right unitaries from these building blocks.
# This is where one can encode discrete choices reminiscent of different
# residual flavor symmetries for quarks vs leptons.

def sector_unitaries():
    """
    Return a dictionary mapping sector names to (U_L, U_R).
    Choices are discrete and parameter-free except for the fixed 30° rotation,
    which is a symbolic analogue of a small mixing angle.
    """
    sectors = {}

    # Up-type quarks: nearly aligned with the R-eigenbasis.
    U_L_u = I3
    U_R_u = I3

    # Down-type quarks: left-handed fields twisted by a small 2–3 rotation;
    # right-handed fields rotated by a discrete Fourier transform.
    U_L_d = R23_30deg @ I3
    U_R_d = F3

    # Charged leptons: left-handed basis given by F3, right-handed almost diagonal.
    U_L_e = F3
    U_R_e = I3

    # Neutrinos: left-handed aligned to a rotated F3, right-handed in a permuted basis.
    U_L_n = R23_30deg @ F3
    U_R_n = P_23 @ F3

    sectors["u"]  = (U_L_u, U_R_u)
    sectors["d"]  = (U_L_d, U_R_d)
    sectors["e"]  = (U_L_e, U_R_e)
    sectors["nu"] = (U_L_n, U_R_n)

    return sectors


def build_yukawas(F_diag, sectors):
    """
    Build Yukawa matrices for all sectors:
        Y_s = U_L^s†  F_diag  U_R^s
    """
    Y = {}
    for name, (U_L, U_R) in sectors.items():
        Y[name] = U_L.conj().T @ F_diag @ U_R
    return Y


# ============================================================
# 4. Diagonalization and mixing (CKM/PMNS analogues)
# ============================================================

def diagonalize_dirac(Y):
    """
    Dirac-like Yukawa diagonalization via SVD:
        Y = U_L diag(s) U_R†
    Returns:
        U_L, s_vals, U_R
    """
    U_L, s_vals, U_Rh = np.linalg.svd(Y)
    U_R = U_Rh.conj().T
    return U_L, s_vals, U_R


def mixing_matrix(U_L_up, U_L_down):
    """CKM/PMNS-like mixing matrix: V = U_L_up† U_L_down."""
    return U_L_up.conj().T @ U_L_down


def mixing_angles_from_U(U):
    """
    Extract approximate (θ12, θ23, θ13) mixing angles from a unitary 3×3 matrix U
    using PDG-like conventions on |U|, ignoring CP phases:
        s13 = |U_13|
        c13 = sqrt(1 - s13^2)
        s12 = |U_12| / c13
        s23 = |U_23| / c13
    """
    U_abs = np.abs(U)
    s13 = U_abs[0, 2]
    c13 = math.sqrt(max(0.0, 1.0 - s13**2))
    if c13 < 1e-12:
        return 0.0, 0.0, math.pi / 2.0
    s12 = U_abs[0, 1] / c13
    s23 = U_abs[1, 2] / c13
    s12 = max(-1.0, min(1.0, s12))
    s23 = max(-1.0, min(1.0, s23))
    theta12 = math.asin(s12)
    theta23 = math.asin(s23)
    theta13 = math.asin(s13)
    return theta12, theta23, theta13


# ============================================================
# 5. Main routine: build everything and print emergent structure
# ============================================================

def main():
    print("=== Base-360 spectral data on generation space ===")
    print("k-values (triad indices):", k_vals)
    print("Eigenvalues of R (chi_j):", chi)
    print("Universal kernel values F(chi_j):", F_vals)
    print()

    sectors = sector_unitaries()
    Y = build_yukawas(F_diag, sectors)

    # Diagonalize Yukawas
    Uu_L, su, Uu_R = diagonalize_dirac(Y["u"])
    Ud_L, sd, Ud_R = diagonalize_dirac(Y["d"])
    Ue_L, se, Ue_R = diagonalize_dirac(Y["e"])
    Un_L, sn, Un_R = diagonalize_dirac(Y["nu"])

    print("=== Yukawa singular values (up to overall scale) ===")
    print("Up-type (su):        ", su)
    print("Down-type (sd):      ", sd)
    print("Charged leptons (se):", se)
    print("Neutrino Dirac (sn): ", sn)
    print()

    # Mixing matrices
    V_ckm  = mixing_matrix(Uu_L, Ud_L)
    U_pmns = mixing_matrix(Ue_L, Un_L)

    theta12_q, theta23_q, theta13_q = mixing_angles_from_U(V_ckm)
    theta12_l, theta23_l, theta13_l = mixing_angles_from_U(U_pmns)

    print("=== CKM-like mixing matrix (quarks) ===")
    print(V_ckm)
    print("Approx CKM mixing angles (radians):")
    print(f"theta12_q ≈ {theta12_q:.3f}, theta23_q ≈ {theta23_q:.3f}, theta13_q ≈ {theta13_q:.3f}")
    print()

    print("=== PMNS-like mixing matrix (leptons) ===")
    print(U_pmns)
    print("Approx PMNS mixing angles (radians):")
    print(f"theta12_l ≈ {theta12_l:.3f}, theta23_l ≈ {theta23_l:.3f}, theta13_l ≈ {theta13_l:.3f}")
    print()

    print("NOTES:")
    print("- We did NOT choose a concrete geometric graph (360-ring, 24-cell, etc.).")
    print("  Instead we postulated:")
    print("    * A base-360 unitary R with three characters chi_j = e^{2πi k_j / 360}.")
    print("    * A universal spectral kernel F(chi_j) = exp[-(1 - Re chi_j)].")
    print("    * Simple, parameter-free left/right unitaries built from discrete")
    print("      Fourier transforms, fixed rotations, and permutations.")
    print("- Yukawas in each sector are Y_s = U_L^s† F_diag U_R^s, i.e. a finite")
    print("  Dirac-operator-like object constrained by the base-360 flavor symmetry.")
    print("- All hierarchies and mixing arise from:")
    print("    * The spectrum of R (via F), and")
    print("    * The discrete choice of left/right flavor bases per sector.")
    print("- This is the finite/internal part of a spectral triple toy model, ready")
    print("  to be tensored with a spacetime Dirac operator in a full NCG setup.")
    print("- To connect with your previous ring/24-cell geometries, one would next:")
    print("    * Realize R as an actual rotation/translation on a 360-node graph,")
    print("    * Or embed this Z_360 action into a larger root lattice (e.g. F4), and")
    print("    * Use that geometry to justify specific choices of U_L^s, U_R^s.")

if __name__ == "__main__":
    main()