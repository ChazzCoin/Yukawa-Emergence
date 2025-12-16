import numpy as np
import math

"""
Resonant Spectral Triple Flavor Toy (with sector charge operator Q)
-------------------------------------------------------------------

This extends the previous operator-first flavor toy by adding:

- A second internal operator Q, represented here by discrete "sector charges" q_s.
- For each sector s in {u, d, e, nu}, we define a real charge q_s ∈ Z.
- The universal spectral kernel is modified per sector as:

    F_s(chi_j) = exp( -(1 - Re chi_j) ) * exp( -beta * q_s )

  where beta is a fixed constant (no continuous fitting), and q_s are
  discrete integers. This yields different overall weights ("mass scales")
  for quarks vs leptons vs neutrinos, while keeping the same triadic
  pattern across generations.

We still do NOT introduce random matrices or sector-specific exponent tables.
All structure is from:

- The base-360 unitary R (order 360),
- The spectral function F(chi),
- The discrete sector charges q_s (second internal operator Q),
- Simple, group-like left/right unitaries per sector.
"""

# ============================================================
# 1. Base-360 cyclic operator on generation space
# ============================================================

BASE_ANGLE = 2.0 * math.pi / 360.0
k_vals = np.array([6, 3, 2], dtype=int)  # triadic "frequencies" in base-360

# Eigenvalues (characters) of R on the 3 generations
chi = np.exp(1j * BASE_ANGLE * k_vals)
R_diag = np.diag(chi)


# ============================================================
# 2. Universal spectral kernel F(chi) (before sector charges)
# ============================================================

def spectral_kernel_base(chi_vals, lam=1.0):
    """
    Base spectral kernel f(chi) = exp( -lam * (1 - Re chi) ),
    same for all sectors BEFORE including Q.
    """
    chi_real = np.real(chi_vals)
    return np.exp(-lam * (1.0 - chi_real))

F_vals_base = spectral_kernel_base(chi)   # shape (3,)
F_diag_base = np.diag(F_vals_base)        # 3×3 base diagonal kernel


# ============================================================
# 3. Sector charge operator Q: discrete charges q_s
# ============================================================

# We now introduce a second internal "operator" Q, implemented as a sector label
# with a discrete charge q_s ∈ Z. Conceptually, Q is a diagonal operator in
# sector space, but here we just use its eigenvalues q_s as a multiplicative
# factor in the kernel.
#
# The sector-dependent kernel is:
#
#    F_s = exp( - (1 - Re chi) ) * exp( -beta * q_s )
#
# beta is fixed; q_s are discrete integers. This yields:
# - a common triadic pattern across generations (from chi),
# - different overall scales for different sectors (from q_s).
#
# Choice (example):
#   q_u  = 0   (up quarks: reference scale, heaviest sector)
#   q_d  = 1   (down quarks: somewhat suppressed)
#   q_e  = 2   (charged leptons: more suppressed)
#   q_nu = 3   (neutrinos: most suppressed)

beta = 1.0  # fixed scale for Q's contribution (no fitting)

sector_charges = {
    "u":  0,  # up-type quarks
    "d":  1,  # down-type quarks
    "e":  2,  # charged leptons
    "nu": 3,  # neutrinos
}


def sector_kernel_diag(q_s):
    """
    Given a sector charge q_s, return the corresponding diagonal kernel:

        F_s_diag = exp( -beta * q_s ) * F_diag_base.

    This is equivalent to applying the operator exp( -beta Q_s ) with
    Q_s = q_s * I_3 on top of the base kernel.
    """
    scale = math.exp(-beta * q_s)
    return scale * F_diag_base


# ============================================================
# 4. Left/right flavor bases and Yukawa construction
# ============================================================

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


F3 = unitary_F3()
I3 = np.eye(3, dtype=complex)
R23_30deg = real_rotation_23(math.pi / 6.0).astype(complex)

P_23 = np.array([
    [1.0, 0.0, 0.0],
    [0.0, 0.0, 1.0],
    [0.0, 1.0, 0.0],
], dtype=complex)


def sector_unitaries():
    """
    Return a dictionary mapping sector names to (U_L, U_R).
    Choices are discrete and parameter-free except for the fixed 30° rotation.
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


def build_yukawas(F_diag_base, sector_charges, sectors):
    """
    Build Yukawa matrices for all sectors:

        F_s_diag = exp( -beta * q_s ) * F_diag_base
        Y_s      = U_L^s†  F_s_diag  U_R^s
    """
    Y = {}
    for name, (U_L, U_R) in sectors.items():
        q_s = sector_charges[name]
        F_s_diag = sector_kernel_diag(q_s)
        Y[name] = U_L.conj().T @ F_s_diag @ U_R
    return Y


# ============================================================
# 5. Diagonalization and mixing (CKM/PMNS analogues)
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
    Extract approximate (θ12, θ23, θ13) from a unitary 3×3 matrix U
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
# 6. Main: build everything and print emergent structure
# ============================================================

def main():
    print("=== Base-360 spectral data on generation space ===")
    print("k-values (triad indices):", k_vals)
    print("Eigenvalues of R (chi_j):", chi)
    print("Base kernel values F(chi_j):", F_vals_base)
    print()

    print("Sector charges q_s (for Q):")
    for name, q in sector_charges.items():
        print(f"  {name}: q_{name} = {q}")
    print(f"beta (fixed) = {beta}")
    print()

    sectors = sector_unitaries()
    Y = build_yukawas(F_diag_base, sector_charges, sectors)

    # Diagonalize Yukawas
    Uu_L, su, Uu_R = diagonalize_dirac(Y["u"])
    Ud_L, sd, Ud_R = diagonalize_dirac(Y["d"])
    Ue_L, se, Ue_R = diagonalize_dirac(Y["e"])
    Un_L, sn, Un_R = diagonalize_dirac(Y["nu"])

    print("=== Yukawa singular values (sector spectra, up to overall scale) ===")
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
    print("- We introduced a second internal operator Q via sector charges q_s.")
    print("- The sector-dependent kernel is F_s = exp(-(1 - Re chi)) * exp(-beta * q_s).")
    print("- This modifies the overall scale of Yukawa eigenvalues per sector,")
    print("  while preserving the triadic pattern from the base-360 spectrum.")
    print("- No randomness, no sector-specific exponent tables; all differences")
    print("  between sectors arise from discrete charges q_s and discrete")
    print("  flavor rotations (U_L^s, U_R^s).")
    print("- Next refinements would:")
    print("    * Allow Q to act nontrivially in generation space (not just scalar),")
    print("    * Or derive q_s and U_L^s, U_R^s from a discrete flavor group (A4, etc.),")
    print("    * Or embed this internal Z_360 × Q structure into a full NCG spectral triple.")

if __name__ == "__main__":
    main()