import numpy as np
import math

"""
Operator-first quasi-crystal flavor toy
=======================================

This script builds a minimal "operator-first" flavor model where:

- The internal space is a 1D Fibonacci quasi-crystal graph (aperiodic chain).
- Evolution is encoded in a spectral kernel F(lambda) = exp(-alpha * lambda)
  acting on the graph Laplacian eigenvalues.
- A discrete charge operator Q_s (different per sector and generation) creates
  hierarchical weights ~ exp(-beta * q_{s,g}).
- Simple, group-like left/right unitaries define flavor bases per sector.
- Yukawa matrices are built as Y_s = U_L^s† diag(F_s) U_R^s.
- Mass hierarchies and mixing matrices are extracted and compared to rough
  SM-like targets via a chi^2.

This is NOT a precision SM model; it is a clean, operator-first, quasi-crystal
proof-of-concept with no fitted continuous parameters.
"""

# ----------------------------------------------------------------------
# 1. Fibonacci quasi-crystal graph and Laplacian
# ----------------------------------------------------------------------

def fibonacci_word(n_iter=7, seed="A"):
    """Generate a Fibonacci substitution word: A->AB, B->A."""
    s = seed
    for _ in range(n_iter):
        s = "".join("AB" if c == "A" else "A" for c in s)
    return s

def fibonacci_graph(n_iter=7):
    """
    Build a 1D Fibonacci chain (aperiodic) with nearest-neighbor edges.
    Edge weights depend slightly on local A/B pattern to encode a weak
    quasi-crystal modulation.
    Returns: (word, adjacency A, Laplacian L)
    """
    word = fibonacci_word(n_iter)
    N = len(word)
    A = np.zeros((N, N), dtype=float)

    # connect nearest neighbors with base weight 1.0
    for i in range(N - 1):
        A[i, i+1] = A[i+1, i] = 1.0

    # modulate edge weights by local pattern
    for i in range(N - 1):
        pair = word[i:i+2]
        if pair == "AA":
            w = 1.0
        elif pair == "AB":
            w = 1.1
        elif pair == "BA":
            w = 0.9
        else:  # "BB" (rare in Fibonacci word of this depth)
            w = 1.0
        A[i, i+1] = A[i+1, i] = w

    D = np.diag(A.sum(axis=1))
    L = D - A
    return word, A, L

word, A_crystal, L_crystal = fibonacci_graph(7)
N_sites = len(word)

# Diagonalize the crystal Laplacian
eigvals, eigvecs = np.linalg.eigh(L_crystal)

# Select three low-lying non-zero modes as our "generation triad"
# (skip index 0 which is the zero mode)
lam_gen = eigvals[1:4]

# ----------------------------------------------------------------------
# 2. Spectral kernel F(lambda) from the crystal
# ----------------------------------------------------------------------

def base_kernel(lams, alpha=10.0):
    """
    Spectral kernel from the crystal: F(lambda) = exp(-alpha * lambda).
    alpha > 0 sets how strongly higher Laplacian modes are suppressed.
    """
    return np.exp(-alpha * lams)

alpha = 10.0  # fixed, not fitted
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
    "u":  np.array([2, 1, 0], dtype=float),  # u, c, t
    "d":  np.array([3, 2, 1], dtype=float),  # d, s, b
    "e":  np.array([4, 3, 2], dtype=float),  # e, mu, tau
    "nu": np.array([6, 5, 4], dtype=float),  # nu1, nu2, nu3
}

def sector_weights(F_base, q_vec, beta=1.0):
    """
    For a given sector with charge vector q_vec, build its diagonal weights:
        F_s(g) = F_base(g) * exp(-beta * q_vec[g]).
    """
    return F_base * np.exp(-beta * q_vec)

# ----------------------------------------------------------------------
# 4. Flavor bases: simple discrete unitaries per sector
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
    Choose sector-dependent left/right unitaries (U_L, U_R).
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
    print("=== Internal quasi-crystal (Fibonacci chain) ===")
    print(f"Number of sites: {N_sites}")
    print("First 10 Laplacian eigenvalues:")
    print(eigvals[:10])
    print()
    print("Selected triad eigenvalues (lam_gen):", lam_gen)
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
    print(f"theta12_q ≈ {theta12_q:.3f}, theta23_q ≈ {theta23_q:.3f}, theta13_q ≈ {theta13_q:.3f}")
    print()

    print("=== PMNS-like mixing matrix ===")
    print(U_pmns)
    print("Mixing angles (radians):")
    print(f"theta12_l ≈ {theta12_l:.3f}, theta23_l ≈ {theta23_l:.3f}, theta13_l ≈ {theta13_l:.3f}")
    print()

    print("=== Observables vs rough targets ===")
    for k, m, t, contrib in chi2_details:
        print(f"{k:12s}: model={m:.3e}, target={t:.3e}, chi2_contrib={contrib:.2f}")
    print()
    print(f"Total chi^2 ≈ {chi2_value:.2f}")
    print()

    print("NOTES:")
    print("- The internal space is a simple 1D Fibonacci quasi-crystal graph.")
    print("- The triad of 'generation' scales comes from the Laplacian spectrum")
    print("  via F_base(lambda) = exp(-alpha * lambda).")
    print("- Sector and generation hierarchies arise from discrete charges q_{s,g}")
    print("  via exp(-beta * q_{s,g}).")
    print("- Left/right flavor bases are discrete unitaries (I, F3, 30° rotations).")
    print("- No continuous parameters were tuned to fit data; alpha and beta are")
    print("  fixed, and q_{s,g} are small integers chosen by hand.")
    print("- The resulting chi^2 is large: this model captures the qualitative")
    print("  idea of hierarchical sectors and a Cabibbo-like mixing angle, but")
    print("  it does NOT quantitatively reproduce SM flavor data.")
    print("- This is an operator-first quasi-crystal toy, a starting point for")
    print("  embedding more realistic internal geometries and constraints.")

if __name__ == "__main__":
    main()