import numpy as np
import math

"""
Operator-first quasi-crystal flavor toy
=======================================

Internal space:
  - A finite 2D Fibonacci quasi-crystal patch:
      * Start with a 1D Fibonacci chain (A/B sequence) along x,
      * and another along y,
      * build the Cartesian product graph of the two chains.
    This gives an aperiodic 2D grid with locally crystalline structure.

  - Vertices are pairs (i_x, i_y),
  - Edges connect nearest neighbors along x or y,
  - Edge weights are modulated by local A/B patterns in each direction.

  This is a simple, axiom-consistent quasi-crystal graph G with Laplacian L.

Flavor construction (axiom-driven):
  - Take three smallest nonzero Laplacian eigenvalues as a "generation triad".
  - Define a spectral kernel F(L) = exp(-alpha * lambda) on those modes
    (Evolution operator acting in internal space).
  - Add a discrete charge operator Q_s per sector & generation, with weights
        F_s(g) = F_base(g) * exp(-beta * q_{s,g}).
  - Build Yukawas as
        Y_s = U_L^s† diag(F_s) U_R^s
    using simple discrete unitaries (I, F3, 30° rotations, permutations).
  - Diagonalize Y_s via SVD to get singular values (mass hierarchies)
    and left-handed unitaries (mixing matrices).
  - Compute a simple chi^2 against rough SM-like targets.

This is NOT a realistic SM model; it's the first "axiom-driven" quasi-crystal
prototype with no tuned continuous parameters beyond fixed alpha, beta.
"""

# ----------------------------------------------------------------------
# 1. Build 1D Fibonacci chains and their Laplacians
# ----------------------------------------------------------------------

def fibonacci_word(n_iter: int = 7, seed: str = "A") -> str:
    """
    Generate a Fibonacci substitution word:
      A -> AB
      B -> A
    After n_iter substitutions starting from 'seed'.
    """
    s = seed
    for _ in range(n_iter):
        s = "".join("AB" if c == "A" else "A" for c in s)
    return s

def fibonacci_chain_adjacency(word: str) -> np.ndarray:
    """
    Build adjacency matrix for a 1D Fibonacci chain with nearest-neighbor edges.
    Edge weights depend slightly on local A/B pattern to encode a weak
    quasi-crystal modulation along the chain.

    Nodes are indexed 0..N-1 corresponding to characters in 'word'.
    """
    N = len(word)
    A = np.zeros((N, N), dtype=float)

    # base connectivity
    for i in range(N - 1):
        A[i, i+1] = A[i+1, i] = 1.0

    # modulate weights by local pattern (AA, AB, BA)
    for i in range(N - 1):
        pair = word[i:i+2]
        if pair == "AA":
            w = 1.0
        elif pair == "AB":
            w = 1.1
        elif pair == "BA":
            w = 0.9
        else:  # "BB" (rare)
            w = 1.0
        A[i, i+1] = A[i+1, i] = w

    return A

def laplacian_from_adjacency(A: np.ndarray) -> np.ndarray:
    """Graph Laplacian L = D - A."""
    D = np.diag(A.sum(axis=1))
    return D - A

# Build 1D Fibonacci chains along x and y
word_x = fibonacci_word(6)   # length ~ 21
word_y = fibonacci_word(5)   # length ~ 13

A_x = fibonacci_chain_adjacency(word_x)
A_y = fibonacci_chain_adjacency(word_y)

L_x = laplacian_from_adjacency(A_x)
L_y = laplacian_from_adjacency(A_y)

N_x = A_x.shape[0]
N_y = A_y.shape[0]

# ----------------------------------------------------------------------
# 2. Build 2D quasi-crystal patch as Cartesian product of chains
# ----------------------------------------------------------------------

def cartesian_product_adjacency(A1: np.ndarray, A2: np.ndarray) -> np.ndarray:
    """
    Adjacency of the Cartesian product graph G = G1 □ G2.
    Nodes are pairs (i, j). Two nodes (i1, j1) and (i2, j2) are adjacent if:
      - i1 == i2 and j1, j2 are adjacent in G2, OR
      - j1 == j2 and i1, i2 are adjacent in G1.

    Edge weights inherit from 1D chains.
    """
    N1 = A1.shape[0]
    N2 = A2.shape[0]
    N  = N1 * N2
    A  = np.zeros((N, N), dtype=float)

    def idx(i1: int, i2: int) -> int:
        return i1 * N2 + i2

    # edges along x (vary i1, fixed i2)
    for i1 in range(N1):
        for i1p in range(N1):
            if A1[i1, i1p] != 0.0:
                for j in range(N2):
                    u = idx(i1,  j)
                    v = idx(i1p, j)
                    # weight from A1, treat as "horizontal" bond
                    w = A1[i1, i1p]
                    A[u, v] = A[v, u] = max(A[u, v], w)

    # edges along y (vary i2, fixed i1)
    for j1 in range(N2):
        for j2 in range(N2):
            if A2[j1, j2] != 0.0:
                for i in range(N1):
                    u = idx(i, j1)
                    v = idx(i, j2)
                    # weight from A2, treat as "vertical" bond
                    w = A2[j1, j2]
                    A[u, v] = A[v, u] = max(A[u, v], w)

    return A

A_int = cartesian_product_adjacency(A_x, A_y)
L_int = laplacian_from_adjacency(A_int)
N_sites = A_int.shape[0]

# Diagonalize the internal Laplacian
eigvals, eigvecs = np.linalg.eigh(L_int)

# ----------------------------------------------------------------------
# 3. Select generation triad from Laplacian spectrum
# ----------------------------------------------------------------------

# Sort eigenvalues and pick the three smallest non-zero values
# (skip index 0: the global constant mode with lambda ~ 0)
eps = 1e-10
nonzero_lams = eigvals[eigvals > eps]
lam_gen = nonzero_lams[:3]  # shape (3,)

# ----------------------------------------------------------------------
# 4. Spectral kernel F(lambda) from the quasi-crystal
# ----------------------------------------------------------------------

def base_kernel(lams: np.ndarray, alpha: float = 3.0) -> np.ndarray:
    """
    Spectral kernel from the internal quasi-crystal:
        F(lambda) = exp(-alpha * lambda),
    where alpha > 0 sets how strongly higher Laplacian modes are suppressed.

    (Axiom: part of Evolution operator acting in internal space.)
    """
    return np.exp(-alpha * lams)

alpha = 3.0  # fixed, not fitted
F_base = base_kernel(lam_gen, alpha=alpha)  # shape (3,)

# ----------------------------------------------------------------------
# 5. Discrete charge operator Q_s: sector + generation hierarchy
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

def sector_weights(F_base: np.ndarray, q_vec: np.ndarray, beta: float = 1.0) -> np.ndarray:
    """
    For a given sector with charge vector q_vec, build its diagonal weights:
        F_s(g) = F_base(g) * exp(-beta * q_vec[g]).

    (Axiom: evolution from spectral kernel + charge operator Q_s.)
    """
    return F_base * np.exp(-beta * q_vec)

# ----------------------------------------------------------------------
# 6. Flavor bases: discrete unitaries per sector (U_L^s, U_R^s)
# ----------------------------------------------------------------------

def unitary_F3() -> np.ndarray:
    """3x3 discrete Fourier transform on Z_3."""
    omega = np.exp(2j * math.pi / 3.0)
    j = np.arange(3)[:, None]
    k = np.arange(3)[None, :]
    F = omega ** (j * k)
    F /= math.sqrt(3.0)
    return F

def real_rotation_23(theta: float) -> np.ndarray:
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
      - up: identity on both sides (closest to eigenbasis)
      - down: small 2-3 rotation on the left, F3 on the right
      - charged leptons: F3 on the left, identity on the right
      - neutrinos: rotated F3 on the left, permuted F3 on the right

    (Axiom: sector structure arises from discrete unitaries, not continuous knobs.)
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
# 7. Build Yukawa matrices from operators
# ----------------------------------------------------------------------

def build_yukawas(F_base: np.ndarray, sector_charges_gen, beta: float = 1.0):
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
# 8. Diagonalization and mixing
# ----------------------------------------------------------------------

def diagonalize_dirac(Y: np.ndarray):
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

def mixing_matrix(U_L_up: np.ndarray, U_L_down: np.ndarray) -> np.ndarray:
    """CKM/PMNS-like mixing matrix: V = U_L_up† U_L_down."""
    return U_L_up.conj().T @ U_L_down

def mixing_angles_from_U(U: np.ndarray):
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
# 9. Compare to rough SM targets via chi^2
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
# 10. Print summary
# ----------------------------------------------------------------------

def main():
    print("=== 2D Fibonacci quasi-crystal internal graph ===")
    print(f"Chain lengths: Nx = {N_x}, Ny = {N_y}, total sites = {N_sites}")
    print("First 10 Laplacian eigenvalues (internal L_int):")
    print(eigvals[:10])
    print()
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
    print("- Internal space is a 2D Cartesian product of 1D Fibonacci chains,")
    print("  i.e. a simple quasi-crystal patch graph with aperiodic modulation.")
    print("- Generation triad comes from the three smallest nonzero eigenvalues")
    print("  of the internal Laplacian L_int.")
    print("- F_base(lambda) = exp(-alpha * lambda) with alpha =", alpha)
    print("- Sector + generation hierarchies arise from discrete charges q_{s,g}")
    print("  via exp(-beta * q_{s,g}) with beta =", beta)
    print("- Left/right flavor bases are discrete unitaries (I, F3, 30° rotations, permutations).")
    print("- No continuous parameters were tuned to fit data; alpha, beta are fixed,")
    print("  and q_{s,g} are small integers chosen by hand.")
    print("- The resulting chi^2 will almost certainly be large; the point here is")
    print("  to have our first fully axiom-driven, quasi-crystal-based operator toy.")
    print("- Next refinements would derive q_{s,g} and U_L^s, U_R^s, and even alpha")
    print("  from the symmetry/structure of the quasi-crystal graph itself.")

if __name__ == "__main__":
    main()

"""
RESULTS:
=== 2D Fibonacci quasi-crystal internal graph ===
Chain lengths: Nx = 21, Ny = 13, total sites = 273
First 10 Laplacian eigenvalues (internal L_int):
[2.78237743e-16 2.21718307e-02 5.76085756e-02 7.97804063e-02
 8.81368334e-02 1.45745409e-01 1.97275936e-01 2.27661943e-01
 2.49833774e-01 2.54884511e-01]

Chosen generation triad (lam_gen): [0.02217183 0.05760858 0.07978041]
Base kernel values F_base(lam_gen): [0.93564842 0.84128422 0.78714625]

=== Yukawa singular values (up to overall scale) ===
Up-type (su):         [0.78714625 0.30949117 0.12662624]
Down-type (sd):       [0.28957492 0.11385544 0.04658319]
Charged leptons (se): [0.10652866 0.04188507 0.017137  ]
Neutrino Dirac (sn):  [0.01441709 0.00566853 0.00231924]

=== CKM-like mixing matrix ===
[[-8.66025404e-01+1.30923269e-18j -5.00000000e-01-4.76972962e-17j
  -2.18663788e-16-1.08439876e-16j]
 [ 5.00000000e-01+5.65281876e-17j -8.66025404e-01-1.81833066e-16j
  -6.31836387e-17-9.30279276e-17j]
 [ 1.11022302e-16+4.73977239e-17j -1.66533454e-16+1.34784487e-16j
   1.00000000e+00+6.24069278e-16j]]
Mixing angles (radians):
theta12_q ≈ 0.524, theta23_q ≈ 0.000, theta13_q ≈ 2.441e-16

=== PMNS-like mixing matrix ===
[[-5.59073015e-01+6.61390478e-01j -3.22780956e-01+3.81853970e-01j
   1.37395173e-16+2.91544473e-16j]
 [-9.13247935e-17-5.00000000e-01j  4.43559672e-16+8.66025404e-01j
   3.19445136e-16+2.59080052e-16j]
 [ 1.45716772e-16-2.22044605e-16j  1.66533454e-16-4.30211422e-16j
   1.00000000e+00+1.13570361e-15j]]
Mixing angles (radians):
theta12_l ≈ 0.524, theta23_l ≈ 0.000, theta13_l ≈ 3.223e-16

=== Observables vs rough targets ===
mu/mt       : model=1.609e-01, target=2.200e-05, chi2_contrib=165.90
mc/mt       : model=3.932e-01, target=7.500e-03, chi2_contrib=32.85
md/mb       : model=1.609e-01, target=1.100e-03, chi2_contrib=52.08
ms/mb       : model=3.932e-01, target=2.200e-02, chi2_contrib=17.42
me/mtau     : model=1.609e-01, target=2.900e-04, chi2_contrib=83.67
mmu/mtau    : model=3.932e-01, target=5.900e-02, chi2_contrib=7.54
theta12_q   : model=5.236e-01, target=2.270e-01, chi2_contrib=2.20
theta23_q   : model=1.125e-16, target=4.100e-02, chi2_contrib=0.04
theta13_q   : model=2.441e-16, target=3.600e-03, chi2_contrib=0.00
theta12_l   : model=5.236e-01, target=5.840e-01, chi2_contrib=0.09
theta23_l   : model=4.113e-16, target=7.850e-01, chi2_contrib=15.41
theta13_l   : model=3.223e-16, target=1.500e-01, chi2_contrib=0.56

Total chi^2 ≈ 377.76
"""