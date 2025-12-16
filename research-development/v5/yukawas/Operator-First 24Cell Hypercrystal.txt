import numpy as np
import math

"""
Operator-first hypercrystal flavor toy
======================================

Internal space:
  - 4D regular 24-cell (a hypercrystal with 24 vertices).
  - Vertices at all permutations of (±1, ±1, 0, 0).
  - Edges connect vertices at Euclidean distance sqrt(2).
  - Graph Laplacian L encodes the elastic backbone of the internal medium.

Flavor construction:
  - Take three distinct nonzero Laplacian eigenvalues as a "generation triad".
  - Define a spectral kernel F(lambda) = exp(-alpha * lambda).
  - Add a discrete charge operator Q_s per sector & generation, with weights
        F_s(g) = F_base(g) * exp(-beta * q_{s,g}).
  - Build Yukawas as
        Y_s = U_L^s† diag(F_s) U_R^s
    using simple discrete unitaries (I, F3, 30° rotations, permutations).
  - Diagonalize Y_s via SVD to get singular values (mass hierarchies)
    and left-handed unitaries (mixing matrices).
  - Compute a simple chi^2 against rough SM-like targets.

This is NOT a realistic SM model; it's an honest operator-first hypercrystal
prototype with no tuned continuous parameters beyond fixed alpha, beta.
"""

# ----------------------------------------------------------------------
# 1. Build the 24-cell hypercrystal graph and Laplacian
# ----------------------------------------------------------------------

def vertices_24cell():
    """
    Return the 24 vertices of the 4D regular 24-cell as 4D vectors.
    One standard representation: all permutations of (±1, ±1, 0, 0).
    """
    verts = []
    coords = [-1.0, 1.0]
    for i in range(4):
        for j in range(i + 1, 4):
            for s1 in coords:
                for s2 in coords:
                    v = [0.0, 0.0, 0.0, 0.0]
                    v[i] = s1
                    v[j] = s2
                    verts.append(v)
    return np.array(verts, dtype=float)  # shape (24, 4)


def adjacency_24cell(verts):
    """
    Build adjacency matrix for the 24-cell:
      - Connect vertices whose Euclidean distance is sqrt(2)
        (i.e. squared distance == 2).
      - This yields a regular graph where each vertex has degree 8.
    """
    N = len(verts)
    A = np.zeros((N, N), dtype=float)
    for i in range(N):
        for j in range(i + 1, N):
            d2 = np.sum((verts[i] - verts[j]) ** 2)
            if abs(d2 - 2.0) < 1e-8:
                A[i, j] = 1.0
                A[j, i] = 1.0
    return A


def laplacian_from_adjacency(A):
    D = np.diag(A.sum(axis=1))
    return D - A


# Build hypercrystal
verts_24 = vertices_24cell()
A_hyper = adjacency_24cell(verts_24)
L_hyper = laplacian_from_adjacency(A_hyper)

# Diagonalize Laplacian
eigvals, eigvecs = np.linalg.eigh(L_hyper)
# Unique eigenvalues (rounded to avoid tiny numerical noise)
lam_unique = sorted(set(round(v, 6) for v in eigvals))

# The 24-cell Laplacian eigenvalues are:
# 0 (once), 4 (deg 4), 8 (deg 9), 10 (deg 8), 12 (deg 2).
# Use three distinct nonzero ones as generation triad: 4, 8, 10.
lam_gen = np.array(lam_unique[1:4], dtype=float)

# ----------------------------------------------------------------------
# 2. Spectral kernel F(lambda) from the hypercrystal
# ----------------------------------------------------------------------

def base_kernel(lams, alpha=0.5):
    """
    Spectral kernel from the hypercrystal:
        F(lambda) = exp(-alpha * lambda),
    where alpha > 0 sets how strongly higher Laplacian modes are suppressed.
    """
    return np.exp(-alpha * lams)

alpha = 0.5  # fixed, not fitted
F_base = base_kernel(lam_gen, alpha=alpha)  # shape (3,)

# ----------------------------------------------------------------------
# 3. Discrete charge operator Q_s: sector + generation hierarchy
# ----------------------------------------------------------------------

beta = 1.0  # fixed, not fitted

# Generation-dependent integer charges q_{s,g} (3-vector per sector).
# Larger q => stronger suppression => lighter generation.
# Pattern qualitatively mimics:
#   - up-type: 3rd >> 2nd >> 1st
#   - quarks heavier than leptons
#   - neutrinos lightest
sector_charges_gen = {
    "u":  np.array([2.0, 1.0, 0.0], dtype=float),  # u, c, t
    "d":  np.array([3.0, 2.0, 1.0], dtype=float),  # d, s, b
    "e":  np.array([4.0, 3.0, 2.0], dtype=float),  # e, mu, tau
    "nu": np.array([6.0, 5.0, 4.0], dtype=float),  # nu1, nu2, nu3
}

def sector_weights(F_base, q_vec, beta=1.0):
    """
    For a given sector with charge vector q_vec, build its diagonal weights:
        F_s(g) = F_base(g) * exp(-beta * q_vec[g]).
    """
    return F_base * np.exp(-beta * q_vec)

# ----------------------------------------------------------------------
# 4. Flavor bases: discrete unitaries per sector
# ----------------------------------------------------------------------

def unitary_F3():
    """3x3 discrete Fourier transform on Z_3."""
    omega = np.exp(2j * math.pi / 3.0)
    j = np.arange(3)[:, None]
    k = np.arange(3)[None, :]
    F = omega ** (j * k)
    F /= math.sqrt(3.0)
    return F

def real_rotation_23(theta):
    """
    3x3 real rotation in the (2,3) subspace by angle theta.
    The first generation is left invariant.
    """
    c = math.cos(theta)
    s = math.sin(theta)
    return np.array([
        [1.0, 0.0, 0.0],
        [0.0, c,   s  ],
        [0.0, -s,  c  ],
    ], dtype=float)

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
    Sector-dependent left/right unitaries (U_L, U_R).
    These are discrete, group-like choices, not fitted:
      - up: nearly aligned with the eigenbasis
      - down: small 2-3 rotation on the left, F3 on the right
      - charged leptons: F3 on the left, identity on the right
      - neutrinos: rotated F3 on the left, permuted F3 on the right
    """
    sectors = {}
    # Up-type quarks
    U_L_u = I3
    U_R_u = I3

    # Down-type quarks
    U_L_d = R23_30deg @ I3
    U_R_d = F3

    # Charged leptons
    U_L_e = F3
    U_R_e = I3

    # Neutrinos
    U_L_n = R23_30deg @ F3
    U_R_n = P_23 @ F3

    sectors["u"]  = (U_L_u, U_R_u)
    sectors["d"]  = (U_L_d, U_R_d)
    sectors["e"]  = (U_L_e, U_R_e)
    sectors["nu"] = (U_L_n, U_R_n)
    return sectors

# ----------------------------------------------------------------------
# 5. Build Yukawa matrices from operators
# ----------------------------------------------------------------------

def build_yukawas(F_base, sector_charges_gen, beta=1.0):
    """
    Build Yukawa matrices for all sectors:
        F_s_diag = diag(F_base(g) * exp(-beta * q_{s,g}))
        Y_s      = U_L^s†  F_s_diag  U_R^s
    """
    sectors = sector_unitaries()
    Y = {}
    for name, (U_L, U_R) in sectors.items():
        q_vec = sector_charges_gen[name]
        weights = sector_weights(F_base, q_vec, beta=beta)
        F_s_diag = np.diag(weights.astype(complex))
        Y[name] = U_L.conj().T @ F_s_diag @ U_R
    return Y

Yukawas = build_yukawas(F_base, sector_charges_gen, beta=beta)

# ----------------------------------------------------------------------
# 6. Diagonalization and mixing
# ----------------------------------------------------------------------

def diagonalize_dirac(Y):
    """
    Dirac-like Yukawa diagonalization via SVD:
        Y = U_L diag(s) U_R†
    Returns U_L, singular values s, U_R.
    """
    U_L, s_vals, U_Rh = np.linalg.svd(Y)
    U_R = U_Rh.conj().T
    return U_L, s_vals, U_R

Uu_L, su, Uu_R = diagonalize_dirac(Yukawas["u"])
Ud_L, sd, Ud_R = diagonalize_dirac(Yukawas["d"])
Ue_L, se, Ue_R = diagonalize_dirac(Yukawas["e"])
Un_L, sn, Un_R = diagonalize_dirac(Yukawas["nu"])

def mixing_matrix(U_L_up, U_L_down):
    """CKM/PMNS-like mixing matrix: V = U_L_up† U_L_down."""
    return U_L_up.conj().T @ U_L_down

def mixing_angles_from_U(U):
    """
    Extract approximate (theta12, theta23, theta13) from a unitary 3x3 matrix U
    using PDG-like conventions on |U|, ignoring CP phases.
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

V_ckm  = mixing_matrix(Uu_L, Ud_L)
U_pmns = mixing_matrix(Ue_L, Un_L)

theta12_q, theta23_q, theta13_q = mixing_angles_from_U(V_ckm)
theta12_l, theta23_l, theta13_l = mixing_angles_from_U(U_pmns)

# ----------------------------------------------------------------------
# 7. Compare to rough SM targets via chi^2
# ----------------------------------------------------------------------

def compute_observables(su, sd, se, sn,
                        theta12_q, theta23_q, theta13_q,
                        theta12_l, theta23_l, theta13_l):
    """
    Build a dictionary of dimensionless observables:
      - mass ratios m1/m3, m2/m3 per sector (using singular values)
      - mixing angles (radians) for CKM- and PMNS-like matrices
    """
    def ratios(s):
        s_sorted = np.sort(s)  # ascending: [light, mid, heavy]
        m1, m2, m3 = s_sorted
        return m1 / m3, m2 / m3

    mu_mt, mc_mt   = ratios(su)
    md_mb, ms_mb   = ratios(sd)
    me_mt, mmu_mt  = ratios(se)

    return {
        "mu/mt":      mu_mt,
        "mc/mt":      mc_mt,
        "md/mb":      md_mb,
        "ms/mb":      ms_mb,
        "me/mtau":    me_mt,
        "mmu/mtau":   mmu_mt,
        "theta12_q":  theta12_q,
        "theta23_q":  theta23_q,
        "theta13_q":  theta13_q,
        "theta12_l":  theta12_l,
        "theta23_l":  theta23_l,
        "theta13_l":  theta13_l,
    }

# Rough target values (order-of-magnitude SM-like, not precise PDG)
targets = {
    # up-type mass ratios
    "mu/mt":     2.2e-5,
    "mc/mt":     7.5e-3,
    # down-type mass ratios
    "md/mb":     1.1e-3,
    "ms/mb":     2.2e-2,
    # charged lepton ratios
    "me/mtau":   2.9e-4,
    "mmu/mtau":  0.059,
    # CKM angles (radians)
    "theta12_q": 0.227,
    "theta23_q": 0.041,
    "theta13_q": 0.0036,
    # PMNS angles (radians)
    "theta12_l": 0.584,
    "theta23_l": 0.785,
    "theta13_l": 0.150,
}

def chi2(observables, targets):
    """
    Simple chi^2:
      - for ratios, use log10 error with sigma_log = 0.3 dex (~ factor 2)
      - for angles, use sigma = 0.2 rad
    """
    chi2_total = 0.0
    details = []

    ratio_keys = ["mu/mt", "mc/mt", "md/mb", "ms/mb", "me/mtau", "mmu/mtau"]
    angle_keys = ["theta12_q", "theta23_q", "theta13_q",
                  "theta12_l", "theta23_l", "theta13_l"]

    # ratios
    for k in ratio_keys:
        m = observables[k]
        t = targets[k]
        if m <= 0 or t <= 0:
            continue
        logm = math.log10(m)
        logt = math.log10(t)
        sigma_log = 0.3
        contrib = ((logm - logt) / sigma_log)**2
        chi2_total += contrib
        details.append((k, m, t, contrib))

    # angles
    for k in angle_keys:
        m = observables[k]
        t = targets[k]
        sigma = 0.2
        contrib = ((m - t) / sigma)**2
        chi2_total += contrib
        details.append((k, m, t, contrib))

    return chi2_total, details

obs = compute_observables(su, sd, se, sn,
                          theta12_q, theta23_q, theta13_q,
                          theta12_l, theta23_l, theta13_l)
chi2_value, chi2_details = chi2(obs, targets)

# ----------------------------------------------------------------------
# 8. Print summary
# ----------------------------------------------------------------------

def main():
    print("=== 24-cell hypercrystal Laplacian spectrum ===")
    print("Eigenvalues:", eigvals)
    print("Distinct eigenvalues (rounded):", lam_unique)
    print("Chosen generation triad (lam_gen):", lam_gen)
    print("Base kernel values F_base(lam_gen):", F_base)
    print()

    print("=== Yukawa singular values (up to overall scale) ===")
    print("Up-type (su):        ", su)
    print("Down-type (sd):      ", sd)
    print("Charged leptons (se):", se)
    print("Neutrino Dirac (sn): ", sn)
    print()

    print("=== CKM-like mixing matrix ===")
    print(V_ckm)
    print("Mixing angles (radians):")
    print(f"theta12_q ≈ {theta12_q:.3f}, theta23_q ≈ {theta23_q:.3f}, theta13_q ≈ {theta13_q:.3e}")
    print()

    print("=== PMNS-like mixing matrix ===")
    print(U_pmns)
    print("Mixing angles (radians):")
    print(f"theta12_l ≈ {theta12_l:.3f}, theta23_l ≈ {theta23_l:.3f}, theta13_l ≈ {theta13_l:.3e}")
    print()

    print("=== Observables vs rough targets ===")
    for k, m, t, contrib in chi2_details:
        print(f"{k:12s}: model={m:.3e}, target={t:.3e}, chi2_contrib={contrib:.2f}")
    print()
    print(f"Total chi^2 ≈ {chi2_value:.2f}")
    print()

    print("NOTES:")
    print("- Internal space is the 4D regular 24-cell hypercrystal.")
    print("- Generation triad comes from distinct nonzero eigenvalues {4,8,10}.")
    print("- F_base(lambda) = exp(-alpha * lambda) with alpha =", alpha)
    print("- Sector + generation hierarchies arise from discrete charges q_{s,g}")
    print("  via exp(-beta * q_{s,g}) with beta =", beta)
    print("- Left/right flavor bases are discrete unitaries (I, F3, 30° rotations, permutations).")
    print("- No continuous parameters were tuned to fit data; alpha, beta are fixed,")
    print("  and q_{s,g} are small integers chosen by hand.")
    print("- The resulting chi^2 is large (~O(10^2)), so this does NOT quantitatively")
    print("  reproduce SM flavor data, but it is a clean operator-first hypercrystal toy.")
    print("- Next steps to improve realism would be to:")
    print("    * derive q_{s,g} and U_L^s, U_R^s from the 24-cell's symmetry group,")
    print("    * or embed a genuinely aperiodic (quasicrystal) internal graph with")
    print("      similar operator structure, then re-run the same chi^2 test.")

if __name__ == "__main__":
    main()

"""
RESULTS:
=== 24-cell hypercrystal Laplacian spectrum ===
Eigenvalues: [7.70207911e-16 4.00000000e+00 4.00000000e+00 4.00000000e+00
 4.00000000e+00 8.00000000e+00 8.00000000e+00 8.00000000e+00
 8.00000000e+00 8.00000000e+00 8.00000000e+00 8.00000000e+00
 8.00000000e+00 8.00000000e+00 1.00000000e+01 1.00000000e+01
 1.00000000e+01 1.00000000e+01 1.00000000e+01 1.00000000e+01
 1.00000000e+01 1.00000000e+01 1.20000000e+01 1.20000000e+01]
Distinct eigenvalues (rounded): [np.float64(0.0), np.float64(4.0), np.float64(8.0), np.float64(10.0), np.float64(12.0)]
Chosen generation triad (lam_gen): [ 4.  8. 10.]
Base kernel values F_base(lam_gen): [0.13533528 0.01831564 0.00673795]

=== Yukawa singular values (up to overall scale) ===
Up-type (su):         [0.01831564 0.00673795 0.00673795]
Down-type (sd):       [0.00673795 0.00247875 0.00247875]
Charged leptons (se): [0.00247875 0.00091188 0.00091188]
Neutrino Dirac (sn):  [0.00033546 0.00012341 0.00012341]

=== CKM-like mixing matrix ===
[[-1.00000000e+00+2.84988656e-17j -5.55111512e-17-5.47781116e-17j
   6.14231114e-17+1.04756333e-16j]
 [-1.44603563e-16-6.88090638e-17j -9.65925826e-01+1.32258969e-16j
  -2.13198569e-16-2.58819045e-01j]
 [ 1.01880477e-17+4.51525510e-17j -2.58819045e-01-8.67884188e-17j
   3.39510135e-16+9.65925826e-01j]]
Mixing angles (radians):
theta12_q ≈ 0.000, theta23_q ≈ 0.262, theta13_q ≈ 1.214e-16

=== PMNS-like mixing matrix ===
[[ 1.00000000e+00-1.90800966e-16j  4.50590348e-17-5.42823935e-17j
   1.15626596e-16-8.29562354e-17j]
 [-2.11580703e-16+1.43219951e-16j  1.27578753e-01-4.30755418e-01j
   4.95572220e-01-7.43358330e-01j]
 [ 4.11584032e-17+1.15303861e-16j  2.68819898e-01-8.52003107e-01j
  -2.55771376e-01+3.69333957e-01j]]
Mixing angles (radians):
theta12_l ≈ 0.000, theta23_l ≈ 1.105, theta13_l ≈ 1.423e-16

=== Observables vs rough targets ===
mu/mt       : model=3.679e-01, target=2.200e-05, chi2_contrib=198.18
mc/mt       : model=3.679e-01, target=7.500e-03, chi2_contrib=31.76
md/mb       : model=3.679e-01, target=1.100e-03, chi2_contrib=70.80
ms/mb       : model=3.679e-01, target=2.200e-02, chi2_contrib=16.63
me/mtau     : model=3.679e-01, target=2.900e-04, chi2_contrib=107.01
mmu/mtau    : model=3.679e-01, target=5.900e-02, chi2_contrib=7.02
theta12_q   : model=7.799e-17, target=2.270e-01, chi2_contrib=1.29
theta23_q   : model=2.618e-01, target=4.100e-02, chi2_contrib=1.22
theta13_q   : model=1.214e-16, target=3.600e-03, chi2_contrib=0.00
theta12_l   : model=7.055e-17, target=5.840e-01, chi2_contrib=8.53
theta23_l   : model=1.105e+00, target=7.850e-01, chi2_contrib=2.56
theta13_l   : model=1.423e-16, target=1.500e-01, chi2_contrib=0.56

Total chi^2 ≈ 445.55
"""