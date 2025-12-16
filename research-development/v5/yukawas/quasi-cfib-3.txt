import numpy as np
import math

"""
Operator-first flavor toy with generic internal graph
=====================================================

This script is structured so that:

  - The *internal graph* (and its 'dimension') is just a pluggable choice
    that produces a Laplacian L_int and spectrum {lambda_k}.
  - The *operators* (kernel F(lambda), integer charges Q, golden P_phi,
    base-360 Cabibbo rotation, etc.) are specified once, up front, and
    then applied to whatever internal graph we choose.
  - There are NO sector-dependent continuous parameters:
      * alpha, beta are universal,
      * angles are discrete fractions of 2π (2π/5, 2π/28),
      * charges q_{s,g} are small integers.

The current internal graph model is "fib2d": a finite product of two
Fibonacci chains (one along x, one along y). This is a *test graph*,
not a claim about the true internal dimension; the operator pipeline
works for any other choice of internal graph that provides a Laplacian.
"""

# ----------------------------------------------------------------------
# 0. Global operator / model configuration (discrete, explicit)
# ----------------------------------------------------------------------

INTERNAL_MODEL = "fib2d"   # label for the internal graph model under test

ALPHA = 3.0                # spectral kernel exponent (universal)
BETA  = 1.0                # charge exponent (universal)

PHI_ORDER  = 5             # golden operator: angle = 2π / 5
CAB_DENOM  = 28            # Cabibbo-like angle: 2π / 28

KERNEL_FORM = "lambda_sq"  # "lambda_sq" → F = exp(-alpha * lambda^2)

# Neutrino dressing: extra base-360 rotations around golden core
USE_NEUTRINO_DRESSING = True
N_SOLAR   = 36   # 2π/36 ≈ 10° (1-2)
N_REACTOR = 45   # 2π/45 ≈ 8°  (1-3)

# ----------------------------------------------------------------------
# 1. Generic internal graph interface
# ----------------------------------------------------------------------
# Only requirement: build_internal_graph(model) returns a Laplacian L_int.
# The operator pipeline below never references dimension or geometry
# directly; it only uses the eigenvalues {lambda_k} of L_int.
# ----------------------------------------------------------------------

def laplacian_from_adjacency(A: np.ndarray) -> np.ndarray:
    """Graph Laplacian L = D - A."""
    D = np.diag(A.sum(axis=1))
    return D - A

def build_internal_graph(model: str = "fib2d"):
    """
    Build an internal graph Laplacian L_int according to the chosen model.
    This is the ONLY place where geometry/dimension is specified.

    Currently implemented:
      - "fib2d": product of two finite Fibonacci chains (quasi-crystal-like)
    """
    if model == "fib2d":
        L_int, meta = build_internal_graph_fib2d()
        return L_int, meta
    else:
        raise ValueError(f"Unknown internal model: {model}")

# ----------------------------------------------------------------------
# 1a. Example internal graph: "fib2d" (Fibonacci × Fibonacci)
# ----------------------------------------------------------------------

def fibonacci_word(n_iter: int = 6, seed: str = "A") -> str:
    """
    Generate a Fibonacci substitution word:
      A -> AB
      B -> A
    """
    s = seed
    for _ in range(n_iter):
        s = "".join("AB" if c == "A" else "A" for c in s)
    return s

def fibonacci_chain_adjacency(word: str) -> np.ndarray:
    """
    Adjacency for a 1D Fibonacci chain with nearest-neighbor edges.
    Edge weights are slightly modulated by local A/B pattern to encode
    quasi-crystal-like structure.
    """
    N = len(word)
    A = np.zeros((N, N), dtype=float)

    # base connectivity
    for i in range(N - 1):
        A[i, i+1] = A[i+1, i] = 1.0

    # pattern modulation
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

def cartesian_product_adjacency(A1: np.ndarray, A2: np.ndarray) -> np.ndarray:
    """
    Adjacency of the Cartesian product graph G = G1 □ G2.
    Nodes are pairs (i, j).
    """
    N1 = A1.shape[0]
    N2 = A2.shape[0]
    N  = N1 * N2
    A  = np.zeros((N, N), dtype=float)

    def idx(i1: int, i2: int) -> int:
        return i1 * N2 + i2

    # edges along x
    for i1 in range(N1):
        for i1p in range(N1):
            if A1[i1, i1p] != 0.0:
                w = A1[i1, i1p]
                for j in range(N2):
                    u = idx(i1,  j)
                    v = idx(i1p, j)
                    A[u, v] = A[v, u] = max(A[u, v], w)

    # edges along y
    for j1 in range(N2):
        for j2 in range(N2):
            if A2[j1, j2] != 0.0:
                w = A2[j1, j2]
                for i in range(N1):
                    u = idx(i, j1)
                    v = idx(i, j2)
                    A[u, v] = A[v, u] = max(A[u, v], w)

    return A

def build_internal_graph_fib2d():
    """
    Example internal graph: product of two Fibonacci chains.
    Returns L_int and some metadata for printing.
    """
    # 1D chains
    word_x = fibonacci_word(6)   # length ~ 21
    word_y = fibonacci_word(5)   # length ~ 13

    A_x = fibonacci_chain_adjacency(word_x)
    A_y = fibonacci_chain_adjacency(word_y)

    L_x = laplacian_from_adjacency(A_x)
    L_y = laplacian_from_adjacency(A_y)

    N_x = A_x.shape[0]
    N_y = A_y.shape[0]

    # 2D Cartesian product
    A_int = cartesian_product_adjacency(A_x, A_y)
    L_int = laplacian_from_adjacency(A_int)

    meta = {
        "model": "fib2d",
        "Nx": N_x,
        "Ny": N_y,
        "N_sites": N_x * N_y,
    }
    return L_int, meta


# ----------------------------------------------------------------------
# 2. Operator-level: spectrum → F_base(λ)
# ----------------------------------------------------------------------

def spectral_triad(L_int: np.ndarray, triad_rule: str = "lowest3_nonzero"):
    """
    Given L_int, select a 3-eigenvalue 'generation triad' according to a rule.
    Currently: three lowest nonzero eigenvalues.
    """
    eigvals, eigvecs = np.linalg.eigh(L_int)
    eps = 1e-10
    nonzero_indices = np.where(eigvals > eps)[0]
    gen_indices = nonzero_indices[:3]
    lam_gen = eigvals[gen_indices]
    return lam_gen, gen_indices, eigvals

def base_kernel(lams: np.ndarray, alpha: float, form: str = "lambda_sq") -> np.ndarray:
    """
    Spectral kernel F(lambda) from the internal Laplacian spectrum.

    Currently:
      - "lambda_sq": F = exp(-alpha * lambda^2)

    'alpha' is universal and specified in config.
    """
    if form == "lambda_sq":
        return np.exp(-alpha * (lams ** 2))
    elif form == "lambda":
        return np.exp(-alpha * lams)
    else:
        raise ValueError(f"Unknown kernel form: {form}")


# ----------------------------------------------------------------------
# 3. Integer charges Q_{s,g}: sector + generation hierarchy
# ----------------------------------------------------------------------

def build_sector_charges():
    """
    Discrete integer charges q_{s,g}. These are not tuned continuously;
    they encode how strongly each sector/generation is suppressed.

    Pattern is universal in shape across sectors (same offsets), with
    quarks heavier than leptons, and neutrinos lightest.
    """
    sector_charges_gen = {
        "u":  np.array([2.0, 1.0, 0.0]),  # u, c, t
        "d":  np.array([3.0, 2.0, 1.0]),  # d, s, b
        "e":  np.array([4.0, 3.0, 2.0]),  # e, mu, tau
        "nu": np.array([6.0, 5.0, 4.0]),  # nu1, nu2, nu3
    }
    return sector_charges_gen

def sector_weights(F_base: np.ndarray, q_vec: np.ndarray, beta: float) -> np.ndarray:
    """
    Sector + generation weights:
        F_{s,g} = F_base(g) * exp(-beta * q_{s,g})
    """
    return F_base * np.exp(-beta * q_vec)


# ----------------------------------------------------------------------
# 4. Discrete generation-space operators: P_phi, C_360
# ----------------------------------------------------------------------

def rot12(theta: float) -> np.ndarray:
    """Rotation in 1-2 plane by angle theta."""
    c = math.cos(theta)
    s = math.sin(theta)
    return np.array([
        [ c,  s, 0.0],
        [-s,  c, 0.0],
        [0.0, 0.0, 1.0]
    ], dtype=complex)

def rot23(theta: float) -> np.ndarray:
    """Rotation in 2-3 plane by angle theta."""
    c = math.cos(theta)
    s = math.sin(theta)
    return np.array([
        [1.0, 0.0, 0.0],
        [0.0,  c,  s],
        [0.0, -s,  c]
    ], dtype=complex)

def rot13(theta: float) -> np.ndarray:
    """Rotation in 1-3 plane by angle theta."""
    c = math.cos(theta)
    s = math.sin(theta)
    return np.array([
        [ c, 0.0,  s],
        [0.0, 1.0, 0.0],
        [-s, 0.0,  c]
    ], dtype=complex)

def build_generation_operators(phi_order: int, cab_denom: int):
    """
    Build:
      - golden operator P_phi with angle 2π / phi_order (default: 2π/5),
      - Cabibbo-like base-360 rotation with angle 2π / cab_denom (default: 2π/28).
    """
    theta_phi = 2.0 * math.pi / float(phi_order)
    theta_C   = 2.0 * math.pi / float(cab_denom)

    P_phi_12 = rot12(theta_phi)
    P_phi_23 = rot23(theta_phi)
    C_12     = rot12(theta_C)

    return P_phi_12, P_phi_23, C_12, theta_phi, theta_C


# ----------------------------------------------------------------------
# 5. Sector left/right flavor bases and Yukawas
# ----------------------------------------------------------------------

def build_sector_bases(P_phi_12, P_phi_23, C_12):
    """
    Define left-handed flavor bases U_L^s in terms of P_phi and C_12.
    Right-handed bases are identity in this operator toy.

    Up-type quarks:
      U_L^u = P_phi^(12)
    Down-type quarks:
      U_L^d = P_phi^(12) @ C_12
      -> CKM = U_L^u† U_L^d = C_12 (golden cancels)

    Charged leptons:
      U_L^e = I

    Neutrinos:
      If USE_NEUTRINO_DRESSING:
        U_L^nu = R_12(2π/N_SOLAR) @ P_phi^(23) @ R_13(2π/N_REACTOR)
      Else:
        U_L^nu = P_phi^(23)
    """
    I3 = np.eye(3, dtype=complex)

    # Quarks
    U_L_u  = P_phi_12
    U_L_d  = P_phi_12 @ C_12

    # Charged leptons
    U_L_e  = I3

    # Neutrinos
    if USE_NEUTRINO_DRESSING:
        theta12_nu = 2.0 * math.pi / float(N_SOLAR)
        theta13_nu = 2.0 * math.pi / float(N_REACTOR)
        U_L_nu = rot12(theta12_nu) @ P_phi_23 @ rot13(theta13_nu)
    else:
        U_L_nu = P_phi_23

    U_R_u  = I3
    U_R_d  = I3
    U_R_e  = I3
    U_R_nu = I3

    sector_bases = {
        "u":  (U_L_u,  U_R_u),
        "d":  (U_L_d,  U_R_d),
        "e":  (U_L_e,  U_R_e),
        "nu": (U_L_nu, U_R_nu),
    }
    return sector_bases

def yukawa_from_F_and_UL(F_s: np.ndarray, U_L: np.ndarray, U_R: np.ndarray) -> np.ndarray:
    """
    Build Yukawa-like operator:
        Y_s = U_L^† diag(F_s) U_R

    In this toy, we treat F_s as unnormalized 'mass' scales; mixing comes
    solely from the U_L^s definitions.
    """
    F_diag = np.diag(F_s.astype(complex))
    return U_L.conj().T @ F_diag @ U_R


# ----------------------------------------------------------------------
# 6. Mixing and observables
# ----------------------------------------------------------------------

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

def mass_ratios(F_s: np.ndarray):
    """
    From a 3-component F_s, compute ratios m1/m3 and m2/m3
    (sorted ascending).
    """
    s_sorted = np.sort(F_s)
    m1, m2, m3 = s_sorted
    return m1 / m3, m2 / m3


# ----------------------------------------------------------------------
# 7. Chi^2 comparison vs rough SM targets
# ----------------------------------------------------------------------

def compute_observables(mu_mt, mc_mt, md_mb, ms_mb, me_mt, mmu_mt,
                        theta12_q, theta23_q, theta13_q,
                        theta12_l, theta23_l, theta13_l):
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

TARGETS = {
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

    for k in angle_keys:
        m = observables[k]
        t = targets[k]
        sigma = 0.2
        contrib = ((m - t) / sigma)**2
        chi2_total += contrib
        details.append((k, m, t, contrib))

    return chi2_total, details


# ----------------------------------------------------------------------
# 8. Main: glue the operator pipeline together
# ----------------------------------------------------------------------

def main():
    # Internal graph + spectrum
    L_int, meta = build_internal_graph(INTERNAL_MODEL)
    lam_gen, gen_indices, eigvals = spectral_triad(L_int)

    F_base = base_kernel(lam_gen, alpha=ALPHA, form=KERNEL_FORM)

    # Integer charges
    sector_charges_gen = build_sector_charges()

    F_u = sector_weights(F_base, sector_charges_gen["u"],  BETA)
    F_d = sector_weights(F_base, sector_charges_gen["d"],  BETA)
    F_e = sector_weights(F_base, sector_charges_gen["e"],  BETA)
    F_n = sector_weights(F_base, sector_charges_gen["nu"], BETA)

    # Generation-space operators
    P_phi_12, P_phi_23, C_12, theta_phi, theta_C = build_generation_operators(
        PHI_ORDER, CAB_DENOM
    )
    sector_bases = build_sector_bases(P_phi_12, P_phi_23, C_12)

    # Yukawa-like operators (masses from F_s, mixing from U_L^s)
    U_L_u, U_R_u   = sector_bases["u"]
    U_L_d, U_R_d   = sector_bases["d"]
    U_L_e, U_R_e   = sector_bases["e"]
    U_L_n, U_R_n   = sector_bases["nu"]

    Y_u  = yukawa_from_F_and_UL(F_u, U_L_u,  U_R_u)
    Y_d  = yukawa_from_F_and_UL(F_d, U_L_d,  U_R_d)
    Y_e  = yukawa_from_F_and_UL(F_e, U_L_e,  U_R_e)
    Y_nu = yukawa_from_F_and_UL(F_n, U_L_n,  U_R_n)

    # Mass ratios (from F_s)
    mu_mt, mc_mt   = mass_ratios(F_u)
    md_mb, ms_mb   = mass_ratios(F_d)
    me_mt, mmu_mt  = mass_ratios(F_e)
    # (neutrino absolute scale/ratios omitted here)

    # Mixing matrices
    V_ckm  = mixing_matrix(U_L_u, U_L_d)
    U_pmns = mixing_matrix(U_L_e, U_L_n)

    theta12_q, theta23_q, theta13_q = mixing_angles_from_U(V_ckm)
    theta12_l, theta23_l, theta13_l = mixing_angles_from_U(U_pmns)

    # Observables + chi^2
    obs = compute_observables(mu_mt, mc_mt, md_mb, ms_mb, me_mt, mmu_mt,
                              theta12_q, theta23_q, theta13_q,
                              theta12_l, theta23_l, theta13_l)
    chi2_value, chi2_details = chi2(obs, TARGETS)

    # ------------------------------------------------------------------
    # Print summary (operator-first narrative, geometry as meta)
    # ------------------------------------------------------------------

    print("=== Internal graph model ===")
    print(f"Model label: {meta['model']}")
    print(f"Sites: {meta['N_sites']}")
    print("First 10 Laplacian eigenvalues (L_int):")
    print(eigvals[:10])
    print()

    print("Generation eigenvalue indices:", gen_indices)
    print("Generation triad lam_gen:", lam_gen)
    print("Base kernel F_base(lam_gen):", F_base)
    print()

    print("=== Yukawa-like mass scales F_s ===")
    print("Up-type (F_u):        ", F_u)
    print("Down-type (F_d):      ", F_d)
    print("Charged leptons (F_e):", F_e)
    print("Neutrino (F_n):       ", F_n)
    print()

    print("Mass ratios (m1/m3, m2/m3) from F_s:")
    print(f"mu/mt:     {mu_mt:.3e}, mc/mt:     {mc_mt:.3e}")
    print(f"md/mb:     {md_mb:.3e}, ms/mb:     {ms_mb:.3e}")
    print(f"me/mtau:   {me_mt:.3e}, mmu/mtau:  {mmu_mt:.3e}")
    print()

    print("=== CKM-like mixing matrix (operator-level) ===")
    print(V_ckm)
    print(f"theta12_q ≈ {theta12_q:.3f} rad, theta23_q ≈ {theta23_q:.3f}, theta13_q ≈ {theta13_q:.3e}")
    print(f"(Cabibbo operator angle = 2π/{CAB_DENOM} ≈ {theta_C:.3f} rad)")
    print()

    print("=== PMNS-like mixing matrix (operator-level) ===")
    print(U_pmns)
    print(f"theta12_l ≈ {theta12_l:.3f} rad, theta23_l ≈ {theta23_l:.3f}, theta13_l ≈ {theta13_l:.3e}")
    print(f"(Golden operator angle = 2π/{PHI_ORDER} ≈ {theta_phi:.3f} rad)")
    print()

    print("=== Observables vs rough targets ===")
    for k, m, t, contrib in chi2_details:
        print(f"{k:12s}: model={m:.3e}, target={t:.3e}, chi2_contrib={contrib:.2f}")
    print()
    print(f"Total chi^2 ≈ {chi2_value:.2f}")
    print()

    print("NOTES:")
    print("- This is an operator-first toy:")
    print("    * internal Laplacian L_int comes from a chosen test graph (label 'fib2d'),")
    print("      but the operator pipeline does not assume a particular dimension.")
    print("    * F_base(lambda) = exp(-alpha * lambda^2) with alpha =", ALPHA)
    print("    * integer charges q_{s,g} define sector+generation suppression via exp(-beta q),")
    print("      with beta =", BETA)
    print("    * P_phi encodes a golden 72° rotation on generation space (2π/5),")
    print("    * C_12 encodes a base-360 Cabibbo rotation with angle 2π/{}.".format(CAB_DENOM))
    print("- CKM ≈ C_12 (small Cabibbo-like 1–2 mixing).")
    print("- PMNS ≈ P_phi^(23) (large golden 2–3 mixing).")
    print("- Mass hierarchies are triadic and sector-ordered but too shallow compared to SM.")
    print("- No sector-specific continuous tuning is used; any improvement from here would")
    print("  come from changing the internal graph model (L_int), discrete choices of")
    print("  kernel form/alpha, or integer charge patterns, all treated explicitly.")

if __name__ == "__main__":
    main()

"""
=== Internal graph model ===
Model label: fib2d
Sites: 273
First 10 Laplacian eigenvalues (L_int):
[2.78237743e-16 2.21718307e-02 5.76085756e-02 7.97804063e-02
 8.81368334e-02 1.45745409e-01 1.97275936e-01 2.27661943e-01
 2.49833774e-01 2.54884511e-01]

Generation eigenvalue indices: [1 2 3]
Generation triad lam_gen: [0.02217183 0.05760858 0.07978041]
Base kernel F_base(lam_gen): [0.99852632 0.99009316 0.98108641]

=== Yukawa-like mass scales F_s ===
Up-type (F_u):         [0.13513584 0.36423492 0.98108641]
Down-type (F_d):       [0.0497137  0.13399454 0.36092152]
Charged leptons (F_e): [0.01828865 0.04929384 0.13277561]
Neutrino (F_n):        [0.0024751  0.0066712  0.01796922]

Mass ratios (m1/m3, m2/m3) from F_s:
mu/mt:     1.377e-01, mc/mt:     3.713e-01
md/mb:     1.377e-01, ms/mb:     3.713e-01
me/mtau:   1.377e-01, mmu/mtau:  3.713e-01

=== CKM-like mixing matrix (operator-level) ===
[[ 0.97492791+0.j  0.22252093+0.j  0.        +0.j]
 [-0.22252093+0.j  0.97492791+0.j  0.        +0.j]
 [ 0.        +0.j  0.        +0.j  1.        +0.j]]
theta12_q ≈ 0.224 rad, theta23_q ≈ 0.000, theta13_q ≈ 0.000e+00
(Cabibbo operator angle = 2π/28 ≈ 0.224 rad)

=== PMNS-like mixing matrix (operator-level) ===
[[ 0.95223934+0.j  0.05366024+0.j  0.30060076+0.j]
 [-0.30230886+0.j  0.30432233+0.j  0.90332567+0.j]
 [-0.04300685+0.j -0.95105652+0.j  0.30600966+0.j]]
theta12_l ≈ 0.056 rad, theta23_l ≈ 1.244, theta13_l ≈ 3.053e-01
(Golden operator angle = 2π/5 ≈ 1.257 rad)

=== Observables vs rough targets ===
mu/mt       : model=1.377e-01, target=2.200e-05, chi2_contrib=160.16
mc/mt       : model=3.713e-01, target=7.500e-03, chi2_contrib=31.91
md/mb       : model=1.377e-01, target=1.100e-03, chi2_contrib=48.89
ms/mb       : model=3.713e-01, target=2.200e-02, chi2_contrib=16.73
me/mtau     : model=1.377e-01, target=2.900e-04, chi2_contrib=79.61
mmu/mtau    : model=3.713e-01, target=5.900e-02, chi2_contrib=7.09
theta12_q   : model=2.244e-01, target=2.270e-01, chi2_contrib=0.00
theta23_q   : model=0.000e+00, target=4.100e-02, chi2_contrib=0.04
theta13_q   : model=0.000e+00, target=3.600e-03, chi2_contrib=0.00
theta12_l   : model=5.629e-02, target=5.840e-01, chi2_contrib=6.96
theta23_l   : model=1.244e+00, target=7.850e-01, chi2_contrib=5.27
theta13_l   : model=3.053e-01, target=1.500e-01, chi2_contrib=0.60

Total chi^2 ≈ 357.27

NOTES:
- This is an operator-first toy:
    * internal Laplacian L_int comes from a chosen test graph (label 'fib2d'),
      but the operator pipeline does not assume a particular dimension.
    * F_base(lambda) = exp(-alpha * lambda^2) with alpha = 3.0
    * integer charges q_{s,g} define sector+generation suppression via exp(-beta q),
      with beta = 1.0
    * P_phi encodes a golden 72° rotation on generation space (2π/5),
    * C_12 encodes a base-360 Cabibbo rotation with angle 2π/28.
- CKM ≈ C_12 (small Cabibbo-like 1–2 mixing).
- PMNS ≈ P_phi^(23) (large golden 2–3 mixing).
- Mass hierarchies are triadic and sector-ordered but too shallow compared to SM.
- No sector-specific continuous tuning is used; any improvement from here would
  come from changing the internal graph model (L_int), discrete choices of
  kernel form/alpha, or integer charge patterns, all treated explicitly.
"""