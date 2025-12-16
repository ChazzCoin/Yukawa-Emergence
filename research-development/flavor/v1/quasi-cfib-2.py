import numpy as np
import math

"""
Golden P_phi + base-360 + Q on a quasi-crystal internal space
=============================================================

- Internal space: 2D Fibonacci quasi-crystal patch (Cartesian product of
  two 1D Fibonacci chains).

- Spectrum (L_int) -> 3 low modes -> F_base(lambda) = exp(-alpha * lambda).

- Integer charges Q_{s,g} give sector+generation hierarchies.

- P_phi operator: golden 72° rotation in generation space.

- Base-360-inspired small rotation: C_12 with angle 15° (2π/24).

- Left-handed flavor bases:
    U_L^u  = P_phi^(12)
    U_L^d  = P_phi^(12) @ C_12
    U_L^e  = I
    U_L^nu = P_phi^(23)

- CKM  = U_L^u† U_L^d  = C_12  (small quark mixing)
- PMNS = U_L^e† U_L^nu = P_phi^(23) (large 2–3 neutrino mixing)

- Mass hierarchies from F_base(lambda) and integer charges as before.
"""

# ----------------------------------------------------------------------
# 1. Build 1D Fibonacci chains and 2D quasi-crystal Laplacian
# ----------------------------------------------------------------------

def fibonacci_word(n_iter: int = 6, seed: str = "A") -> str:
    """Generate a Fibonacci substitution word: A->AB, B->A."""
    s = seed
    for _ in range(n_iter):
        s = "".join("AB" if c == "A" else "A" for c in s)
    return s

def fibonacci_chain_adjacency(word: str) -> np.ndarray:
    """
    Adjacency matrix for a 1D Fibonacci chain with nearest-neighbor edges.
    Edge weights modulated by local A/B pattern.
    """
    N = len(word)
    A = np.zeros((N, N), dtype=float)

    for i in range(N - 1):
        A[i, i+1] = A[i+1, i] = 1.0

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
    D = np.diag(A.sum(axis=1))
    return D - A

# 1D chains
word_x = fibonacci_word(6)   # ~21
word_y = fibonacci_word(5)   # ~13

A_x = fibonacci_chain_adjacency(word_x)
A_y = fibonacci_chain_adjacency(word_y)

L_x = laplacian_from_adjacency(A_x)
L_y = laplacian_from_adjacency(A_y)

N_x = A_x.shape[0]
N_y = A_y.shape[0]

def cartesian_product_adjacency(A1: np.ndarray, A2: np.ndarray) -> np.ndarray:
    """
    Adjacency of Cartesian product G = G1 □ G2.
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

A_int = cartesian_product_adjacency(A_x, A_y)
L_int = laplacian_from_adjacency(A_int)
N_sites = A_int.shape[0]

# Diagonalize L_int
eigvals, eigvecs = np.linalg.eigh(L_int)

# ----------------------------------------------------------------------
# 2. Generation triad and spectral kernel F_base
# ----------------------------------------------------------------------

eps = 1e-10
nonzero_indices = np.where(eigvals > eps)[0]
gen_indices = nonzero_indices[:3]   # three lowest nonzero eigenmodes
lam_gen = eigvals[gen_indices]

def base_kernel(lams, alpha=5.0):
    return np.exp(-alpha * (lams**2))  # instead of linear in λ

alpha = 3.0
F_base = base_kernel(lam_gen, alpha=alpha)  # shape (3,)

# ----------------------------------------------------------------------
# 3. Integer charges Q_s: sector + generation hierarchy
# ----------------------------------------------------------------------

beta = 1.0  # universal

sector_charges_gen = {
    "u":  np.array([2.0, 1.0, 0.0]),  # u, c, t
    "d":  np.array([3.0, 2.0, 1.0]),  # d, s, b
    "e":  np.array([4.0, 3.0, 2.0]),  # e, mu, tau
    "nu": np.array([6.0, 5.0, 4.0]),  # nu1, nu2, nu3
}

def sector_weights(F_base: np.ndarray, q_vec: np.ndarray, beta: float = 1.0) -> np.ndarray:
    """F_s(g) = F_base(g) * exp(-beta * q_vec[g])."""
    return F_base * np.exp(-beta * q_vec)

F_u = sector_weights(F_base, sector_charges_gen["u"],  beta)
F_d = sector_weights(F_base, sector_charges_gen["d"],  beta)
F_e = sector_weights(F_base, sector_charges_gen["e"],  beta)
F_n = sector_weights(F_base, sector_charges_gen["nu"], beta)

# ----------------------------------------------------------------------
# 4. Golden P_phi and base-360 Cabibbo operator on generation space
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
    c = math.cos(theta)
    s = math.sin(theta)
    return np.array([
        [ c, 0.0,  s],
        [0.0, 1.0, 0.0],
        [-s, 0.0,  c]
    ], dtype=complex)

theta12_nu = 2.0 * math.pi / 36.0  # ~10°
theta13_nu = 2.0 * math.pi / 45.0  # ~8°


# Golden angle: 72° = 2π/5
theta_phi = 2.0 * math.pi / 5.0
P_phi_12 = rot12(theta_phi)
P_phi_23 = rot23(theta_phi)
U_L_nu = rot12(theta12_nu) @ P_phi_23 @ rot13(theta13_nu)
U_L_e  = np.eye(3, dtype=complex)
# Base-360-inspired small Cabibbo angle: 15° = 2π/24
# theta_C = 2.0 * math.pi / 24.0
# C_12 = rot12(theta_C)
# Base-360-inspired Cabibbo angle: 360 / 28 degrees
theta_C = 2.0 * math.pi / 28.0  # ≈ 12.857°
C_12 = rot12(theta_C)
# ----------------------------------------------------------------------
# 5. Left-handed flavor bases U_L^s (operator-level ansatz)
# ----------------------------------------------------------------------

U_L_u  = P_phi_12
U_L_d  = P_phi_12 @ C_12
U_L_e  = np.eye(3, dtype=complex)
U_L_nu = P_phi_23

# We keep right-handed bases trivial in this toy
U_R_u  = np.eye(3, dtype=complex)
U_R_d  = np.eye(3, dtype=complex)
U_R_e  = np.eye(3, dtype=complex)
U_R_nu = np.eye(3, dtype=complex)

# ----------------------------------------------------------------------
# 6. Yukawa "mass" matrices from F_s and U_L, U_R
# ----------------------------------------------------------------------

def yukawa_from_F_and_UL(F_s: np.ndarray, U_L: np.ndarray, U_R: np.ndarray) -> np.ndarray:
    """
    Build Yukawa as Y_s = U_L^† diag(F_s) U_R.
    We will treat F_s as the (unnormalized) mass eigenvalues;
    mixing comes from the U_L^s we specified.
    """
    F_diag = np.diag(F_s.astype(complex))
    return U_L.conj().T @ F_diag @ U_R

Y_u  = yukawa_from_F_and_UL(F_u, U_L_u,  U_R_u)
Y_d  = yukawa_from_F_and_UL(F_d, U_L_d,  U_R_d)
Y_e  = yukawa_from_F_and_UL(F_e, U_L_e,  U_R_e)
Y_nu = yukawa_from_F_and_UL(F_n, U_L_nu, U_R_nu)

# For this operator-level toy, we'll:
# - take |F_s| as "mass" singular values
# - use the chosen U_L^s to define mixing directly

def mass_ratios(F_s: np.ndarray):
    s_sorted = np.sort(F_s)  # ascending: [light, mid, heavy]
    m1, m2, m3 = s_sorted
    return m1 / m3, m2 / m3

mu_mt, mc_mt   = mass_ratios(F_u)
md_mb, ms_mb   = mass_ratios(F_d)
me_mt, mmu_mt  = mass_ratios(F_e)
# neutrino absolute scale irrelevant for this toy; skip ratios

# ----------------------------------------------------------------------
# 7. Mixing matrices and angles from U_L^s
# ----------------------------------------------------------------------

def mixing_matrix(U_L_up: np.ndarray, U_L_down: np.ndarray) -> np.ndarray:
    return U_L_up.conj().T @ U_L_down

V_ckm  = mixing_matrix(U_L_u, U_L_d)   # should be ~C_12
U_pmns = mixing_matrix(U_L_e, U_L_nu)  # should be ~P_phi_23

def mixing_angles_from_U(U: np.ndarray):
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

theta12_q, theta23_q, theta13_q = mixing_angles_from_U(V_ckm)
theta12_l, theta23_l, theta13_l = mixing_angles_from_U(U_pmns)

# ----------------------------------------------------------------------
# 8. Chi^2 vs rough SM targets
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

targets = {
    "mu/mt":     2.2e-5,
    "mc/mt":     7.5e-3,
    "md/mb":     1.1e-3,
    "ms/mb":     2.2e-2,
    "me/mtau":   2.9e-4,
    "mmu/mtau":  0.059,
    "theta12_q": 0.227,
    "theta23_q": 0.041,
    "theta13_q": 0.0036,
    "theta12_l": 0.584,
    "theta23_l": 0.785,
    "theta13_l": 0.150,
}

def chi2(observables, targets):
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

obs = compute_observables(mu_mt, mc_mt, md_mb, ms_mb, me_mt, mmu_mt,
                          theta12_q, theta23_q, theta13_q,
                          theta12_l, theta23_l, theta13_l)
chi2_value, chi2_details = chi2(obs, targets)

# ----------------------------------------------------------------------
# 9. Print summary
# ----------------------------------------------------------------------

def main():
    print("=== 2D Fibonacci quasi-crystal internal graph ===")
    print(f"Chain lengths: Nx = {N_x}, Ny = {N_y}, total sites = {N_sites}")
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

    print("=== CKM-like mixing matrix (P_phi-cancelled, base-360 Cabibbo) ===")
    print(V_ckm)
    print(f"theta12_q ≈ {theta12_q:.3f} rad, theta23_q ≈ {theta23_q:.3f}, theta13_q ≈ {theta13_q:.3e}")
    print()

    print("=== PMNS-like mixing matrix (golden 2-3 from P_phi^(23)) ===")
    print(U_pmns)
    print(f"theta12_l ≈ {theta12_l:.3f} rad, theta23_l ≈ {theta23_l:.3f}, theta13_l ≈ {theta13_l:.3e}")
    print()

    print("=== Observables vs rough targets ===")
    for k, m, t, contrib in chi2_details:
        print(f"{k:12s}: model={m:.3e}, target={t:.3e}, chi2_contrib={contrib:.2f}")
    print()
    print(f"Total chi^2 ≈ {chi2_value:.2f}")
    print()
    print("NOTES:")
    print("- P_phi^(12) encodes a golden 72° rotation on generations 1-2.")
    print("- For quarks, P_phi cancels between up and down sectors; CKM is just a")
    print("  small base-360 rotation C_12 by 15°, giving a Cabibbo-like angle.")
    print("- For leptons, U_L^nu = P_phi^(23) exposes the golden angle as a large")
    print("  2-3 mixing in the PMNS matrix, while charged leptons remain aligned.")
    print("- Mass hierarchies still come from the quasi-crystal spectrum via")
    print("  F_base(lambda_gen) and integer charges q_{s,g}.")
    print("- No continuous sector-dependent tuning: only universal alpha, beta, and")
    print("  discrete choices of angles (2π/5, 2π/24) tied to golden ratio and base-360.")

if __name__ == "__main__":
    main()
