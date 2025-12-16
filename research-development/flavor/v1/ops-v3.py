import numpy as np
import math

"""
Resonant Spectral Triple Flavor Toy (Q generation-dependent)
------------------------------------------------------------

Upgrades the previous model by letting the second internal operator Q
act nontrivially in generation space:

- For each sector s ∈ {u, d, e, nu}, we assign a 3-component integer
  charge vector q_s = (q_{s,1}, q_{s,2}, q_{s,3}).

- The sector-dependent kernel in the R-eigenbasis becomes:

    F_s(chi_j) = exp( -(1 - Re chi_j) ) * exp( -beta * q_{s,j} )

  so each generation j in sector s gets its own discrete weight.

This generates *intra-sector* 3-level hierarchies (generations) on top of
*inter-sector* hierarchies (quarks vs leptons vs neutrinos), still with:

- One base-360 unitary R (order 360),
- One universal spectral function of R,
- Q as discrete charges (no continuous fitting),
- Simple, group-like left/right unitaries.
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
# 2. Universal spectral kernel F_base(chi) (before Q)
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
# 3. Generation-dependent sector charges Q
# ============================================================

# beta is fixed: overall strength of Q's contribution (no continuous tuning).
beta = 1.0

# Generation-dependent integer charges q_{s,j} (3-vector per sector).
# Larger q => stronger suppression => lighter generation.
# Chosen qualitatively to mimic:
#   - 3rd gen heaviest in each sector,
#   - quarks heavier than leptons,
#   - neutrinos lightest.

sector_charges_gen = {
    # [q_1, q_2, q_3] in some fixed generation ordering
    "u":  np.array([2, 1, 0], dtype=float),  # u, c, t
    "d":  np.array([3, 2, 1], dtype=float),  # d, s, b
    "e":  np.array([4, 3, 2], dtype=float),  # e, mu, tau
    "nu": np.array([6, 5, 4], dtype=float),  # nu1, nu2, nu3
}


def sector_kernel_diag(q_vec):
    """
    Given a 3-component charge vector q_vec for a sector, return the
    corresponding diagonal kernel in the R-eigenbasis:

        F_s_diag(j,j) = F_base(j) * exp( -beta * q_vec[j] ).

    This implements a generation-dependent internal operator Q_s = diag(q_vec)
    acting on top of the base kernel.
    """
    # Elementwise: F_base * exp(-beta * q_j)
    weights = F_vals_base * np.exp(-beta * q_vec)
    return np.diag(weights)


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


def build_yukawas(sector_charges_gen, sectors):
    """
    Build Yukawa matrices for all sectors:

        F_s_diag = diag( F_base(j) * exp( -beta * q_{s,j} ) )
        Y_s      = U_L^s†  F_s_diag  U_R^s
    """
    Y = {}
    for name, (U_L, U_R) in sectors.items():
        q_vec = sector_charges_gen[name]
        F_s_diag = sector_kernel_diag(q_vec)
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
    print("Base kernel values F_base(chi_j):", F_vals_base)
    print()

    print("Generation-dependent sector charges q_{s,j}:")
    for name, q_vec in sector_charges_gen.items():
        print(f"  {name}: q_{name} =", q_vec.tolist())
    print(f"beta (fixed) = {beta}")
    print()

    sectors = sector_unitaries()
    Y = build_yukawas(sector_charges_gen, sectors)

    # Diagonalize Yukawas
    Uu_L, su, Uu_R = diagonalize_dirac(Y["u"])
    Ud_L, sd, Ud_R = diagonalize_dirac(Y["d"])
    Ue_L, se, Ue_R = diagonalize_dirac(Y["e"])
    Un_L, sn, Un_R = diagonalize_dirac(Y["nu"])

    print("=== Yukawa singular values (sector + generation hierarchies) ===")
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
    print("- Q now acts nontrivially in generation space via integer charge vectors q_{s,j}.")
    print("- Sector kernels are F_s(j) = F_base(chi_j) * exp(-beta * q_{s,j}).")
    print("- This generates 3-level hierarchies within each sector and between sectors,")
    print("  all from discrete charges, one base-360 operator R, and one spectral function.")
    print("- Mixing angles are still controlled by the discrete choices of U_L^s, U_R^s.")
    print("  To differentiate CKM vs PMNS more, the next step is to tie these unitaries")
    print("  to a discrete flavor group or additional internal operators.")

if __name__ == "__main__":
    main()