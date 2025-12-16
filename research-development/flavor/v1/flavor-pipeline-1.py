#!/usr/bin/env python3
"""
Alignment + Seesaw + 1-loop SM RGE pipeline.

Implements the steps:
0. Global config
1. Alignment kernel K (9x9)
2. Proto-matrices (Gaussian, normalized)
3. Apply alignment (Schur product)
4. 9->3 Schur complement
5. Heavy Majorana M_R
6. Type-I seesaw
7. 1-loop RGEs (SM + Y_nu)
8. Diagonalization, CKM/PMNS
9. Angles, phases, Δm^2
10. Optional χ^2
"""

import numpy as np
from numpy.linalg import svd, eig, solve, cond

# If you want, you can later swap to scipy for Takagi factorization, etc.

# ============================================================
# 0. Global config
# ============================================================

class Config:
    v = 246.0              # Electroweak vev in GeV
    mu0 = 1.0e12           # High scale
    mu_EW = 91.1876        # m_Z
    Lambda_Maj = 1.0e14    # Heavy Majorana scale

    # Alignment: we use suppression convention: epsilon < 1
    kappa = 4.10
    eps = 1.0 / kappa      # ≈ 0.244

    seed = 12345           # Random seed

    # RGE integration
    t0 = np.log(mu0)
    t1 = np.log(mu_EW)
    dt = -0.001             # Run down in energy (t decreases)


# SM one-loop beta coefficients for gauge couplings
B1, B2, B3 = 41/6, -19/6, -7


# ============================================================
# 1. Alignment kernel K (9x9)
# ============================================================

def build_alignment_kernel(eps, N=9):
    """
    Build 9x9 kernel K_ij = eps^{|i-j|} for |i-j| in D_360^{(9)} = {1,2,3,4,5,6,8},
    K=0 for |i-j|=7, K_ii=1.
    We use 0-based indices in code, but |i-j| is the same.
    """
    D_360_9 = {1, 2, 3, 4, 5, 6, 8}
    K = np.zeros((N, N), dtype=float)
    for i in range(N):
        for j in range(N):
            d = abs(i - j)
            if i == j:
                K[i, j] = 1.0
            elif d in D_360_9:
                K[i, j] = eps**d
            elif d == 7:
                K[i, j] = 0.0
            else:
                # For N=9, d can only be 0..8; all allowed distances are covered
                K[i, j] = 0.0
    return K


# ============================================================
# 2. Proto-matrices (Gaussian, normalized)
# ============================================================

def random_complex_matrix(shape, rng):
    """
    Entries ~ complex Gaussian with variance 1:
    X_ij = (N1 + i N2) / sqrt(2), N1,N2 ~ N(0,1).
    """
    real = rng.normal(0.0, 1.0, size=shape)
    imag = rng.normal(0.0, 1.0, size=shape)
    return (real + 1j * imag) / np.sqrt(2.0)


def normalize_by_largest_singular_value(X):
    """
    Normalize X so that its largest singular value is 1.
    """
    s = svd(X, compute_uv=False)
    s_max = np.max(s)
    if s_max == 0:
        return X
    return X / s_max

def check_rge_sanity(Yu, Yd, Ye, Ynu, step_index=None):
    Yu_norm = np.linalg.norm(Yu)
    Yd_norm = np.linalg.norm(Yd)
    Ye_norm = np.linalg.norm(Ye)
    Ynu_norm = np.linalg.norm(Ynu)
    Y_max = max(Yu_norm, Yd_norm, Ye_norm, Ynu_norm)
    if Y_max > 10.0:
        raise RuntimeError(f"Yukawa matrix blew up (||Y|| ~ {Y_max}) at step {step_index}")

def check_rge_sanity2(Yu, Yd, Ye, Ynu, g1, g2, g3, step_index=None):
    """
    Simple sanity check to avoid numerical blow-ups.
    Raises RuntimeError if couplings get too large.
    """
    # thresholds are arbitrary but safe for this toy setup
    g_max = max(abs(g1), abs(g2), abs(g3))
    if g_max > 10.0:
        raise RuntimeError(f"Gauge coupling blew up (|g| ~ {g_max}) at step {step_index}")

    # use Frobenius norms for Yukawas
    Yu_norm = np.linalg.norm(Yu)
    Yd_norm = np.linalg.norm(Yd)
    Ye_norm = np.linalg.norm(Ye)
    Ynu_norm = np.linalg.norm(Ynu)

    Y_max = max(Yu_norm, Yd_norm, Ye_norm, Ynu_norm)
    if Y_max > 10.0:
        raise RuntimeError(f"Yukawa matrix blew up (||Y|| ~ {Y_max}) at step {step_index}")

def generate_proto_matrices(config):
    """
    Generate normalized 9x9 proto-matrices for:
    Y_u, Y_d, Y_l, Y_nu, M (Majorana).
    """
    rng = np.random.default_rng(config.seed)

    Yu0 = random_complex_matrix((9, 9), rng)
    Yd0 = random_complex_matrix((9, 9), rng)
    Ye0 = random_complex_matrix((9, 9), rng)
    Ynu0 = random_complex_matrix((9, 9), rng)
    M0  = random_complex_matrix((9, 9), rng)

    Yu0 = normalize_by_largest_singular_value(Yu0)
    Yd0 = normalize_by_largest_singular_value(Yd0)
    Ye0 = normalize_by_largest_singular_value(Ye0)
    Ynu0 = normalize_by_largest_singular_value(Ynu0)
    M0  = normalize_by_largest_singular_value(M0)

    # Enforce Majorana symmetry
    M0 = 0.5 * (M0 + M0.T)

    # --- NEW: globally rescale Yukawas to keep them below 1 ---
    yukawa_scale = 0.5   # you can try 0.3 or 0.2 if it still blows up
    Yu0  *= yukawa_scale
    Yd0  *= yukawa_scale
    Ye0  *= yukawa_scale
    Ynu0 *= yukawa_scale
    # M0 we leave as-is (it’s dimensionless, scale is attached via Λ_Maj)

    return Yu0, Yd0, Ye0, Ynu0, M0


# ============================================================
# 3. Apply alignment Φ: K ⊙ X
# ============================================================

def apply_alignment(K, X):
    """
    Schur (Hadamard) product: Φ(X) = K ⊙ X.
    K is real 9x9, X is complex 9x9.
    """
    return K * X


def align_all(K, Yu0, Yd0, Ye0, Ynu0, M0):
    Yu9  = apply_alignment(K, Yu0)
    Yd9  = apply_alignment(K, Yd0)
    Ye9  = apply_alignment(K, Ye0)
    Ynu9 = apply_alignment(K, Ynu0)
    M9   = apply_alignment(K, M0)
    return Yu9, Yd9, Ye9, Ynu9, M9


# ============================================================
# 4. Seesaw reduction 9→3 for Dirac Yukawas (Schur complement)
# ============================================================

def schur_9_to_3(Y9):
    """
    Y9 is 9x9. Light sites: 0,1,2; heavy: 3..8.
    Y_eff = A - B D^{-1} B†.
    """
    A = Y9[0:3, 0:3]
    B = Y9[0:3, 3:9]
    D = Y9[3:9, 3:9]

    # Check conditioning; if bad, you might want to resample
    if cond(D) > 1e12:
        print("Warning: D block is ill-conditioned (cond ~ {:.2e})".format(cond(D)))

    # Solve D X = B†  → X = D^{-1} B†
    X = solve(D, B.conj().T)       # shape (6,3)
    # Then B D^{-1} B† = B X
    BDinvBdag = B @ X             # shape (3,3)
    Y_eff = A - BDinvBdag
    return Y_eff


def reduce_all_dirac(Yu9, Yd9, Ye9, Ynu9):
    Yu_eff  = schur_9_to_3(Yu9)
    Yd_eff  = schur_9_to_3(Yd9)
    Ye_eff  = schur_9_to_3(Ye9)
    Ynu_eff = schur_9_to_3(Ynu9)
    return Yu_eff, Yd_eff, Ye_eff, Ynu_eff


# ============================================================
# 5. Heavy Majorana matrix M_R (3x3)
# ============================================================

def build_M_R(M9_aligned, Lambda_Maj):
    """
    Project 9x9 aligned Majorana M9 to 3x3 heavy M_R.

    Here we choose sites 3,4,5 as heavy sites 4,5,6 in the paper
    (0-based indices 3,4,5).
    P_H: 9x3 selector. Then M_R ~ Λ_Maj * (P_H^† M9 P_H)_sym.
    """
    PH = np.zeros((9, 3), dtype=float)
    # map (i,α) = (3,0),(4,1),(5,2)
    PH[3, 0] = 1.0
    PH[4, 1] = 1.0
    PH[5, 2] = 1.0

    M3 = PH.T @ M9_aligned @ PH       # 3x3
    M3 = 0.5 * (M3 + M3.T)            # enforce symmetry
    M_R = Lambda_Maj * M3
    return M_R


# ============================================================
# 6. Type-I seesaw: m_ν
# ============================================================

def seesaw_light_neutrinos(Ynu_eff, M_R, v):
    """
    m_D = v/√2 * Ynu_eff
    m_ν = - m_D M_R^{-1} m_D^T
    """
    m_D = (v / np.sqrt(2.0)) * Ynu_eff
    # Solve M_R X = m_D^T -> X = M_R^{-1} m_D^T
    X = solve(M_R, m_D.T)
    m_nu = - m_D @ X
    # Should be complex symmetric; enforce symmetry numerically
    m_nu = 0.5 * (m_nu + m_nu.T)
    return m_nu


# ============================================================
# 7. One-loop RGEs and integrator
# ============================================================

def beta_g(g1, g2, g3):
    """
    Freeze gauge couplings: no running.
    If you later want proper gauge running, replace this with the usual 1-loop betas.
    """
    return 0.0, 0.0, 0.0
# def beta_g(g1, g2, g3):
#     """
#     One-loop gauge beta functions.
#     """
#     dg1 = (B1 * g1**3) / (16*np.pi**2)
#     dg2 = (B2 * g2**3) / (16*np.pi**2)
#     dg3 = (B3 * g3**3) / (16*np.pi**2)
#     return dg1, dg2, dg3


def beta_Yukawas(Yu, Yd, Ye, Ynu, g1, g2, g3):
    """
    One-loop SM+Yν Yukawa RGEs.
    T = Tr(3Yu†Yu + 3Yd†Yd + Ye†Ye)  (we can add Tr(Yν†Yν) if desired).
    """
    Yu_dagYu = Yu.conj().T @ Yu
    Yd_dagYd = Yd.conj().T @ Yd
    Ye_dagYe = Ye.conj().T @ Ye
    Ynu_dagYnu = Ynu.conj().T @ Ynu

    # For simplicity we keep T as in v5 (no Ynu trace), but you can add it.
    T = np.trace(3*Yu_dagYu + 3*Yd_dagYd + Ye_dagYe)

    # Beta(Yu)
    factor_u = T - (17/20*g1**2 + 9/4*g2**2 + 8*g3**2)
    dYu = Yu * factor_u + (3/2)*(Yu @ Yu_dagYu - Yd @ (Yd_dagYd @ Yu))

    # Beta(Yd)
    factor_d = T - (1/4*g1**2 + 9/4*g2**2 + 8*g3**2)
    dYd = Yd * factor_d + (3/2)*(Yd @ Yd_dagYd - Yu @ (Yu_dagYu @ Yd))

    # Beta(Ye)
    factor_e = T - (9/4*g1**2 + 9/4*g2**2)
    dYe = Ye * factor_e + (3/2)*(Ye @ Ye_dagYe)

    # Beta(Ynu)
    factor_nu = T - (9/20*g1**2 + 9/4*g2**2)
    dYnu = Ynu * factor_nu + (3/2)*(Ynu @ Ynu_dagYnu - Ye @ (Ye_dagYe @ Ynu))

    # Overall 1/(16π²)
    dYu /= (16*np.pi**2)
    dYd /= (16*np.pi**2)
    dYe /= (16*np.pi**2)
    dYnu /= (16*np.pi**2)

    return dYu, dYd, dYe, dYnu

def rk4_step_yukawas_only(Yu, Yd, Ye, Ynu, g1, g2, g3, dt):
    """
    RK4 step evolving only Yukawas, with gauge couplings held fixed.
    """
    # k1
    dYu1, dYd1, dYe1, dYnu1 = beta_Yukawas(Yu, Yd, Ye, Ynu, g1, g2, g3)

    # k2
    Yu2 = Yu + 0.5 * dt * dYu1
    Yd2 = Yd + 0.5 * dt * dYd1
    Ye2 = Ye + 0.5 * dt * dYe1
    Ynu2 = Ynu + 0.5 * dt * dYnu1
    dYu2, dYd2, dYe2, dYnu2 = beta_Yukawas(Yu2, Yd2, Ye2, Ynu2, g1, g2, g3)

    # k3
    Yu3 = Yu + 0.5 * dt * dYu2
    Yd3 = Yd + 0.5 * dt * dYd2
    Ye3 = Ye + 0.5 * dt * dYe2
    Ynu3 = Ynu + 0.5 * dt * dYnu2
    dYu3, dYd3, dYe3, dYnu3 = beta_Yukawas(Yu3, Yd3, Ye3, Ynu3, g1, g2, g3)

    # k4
    Yu4 = Yu + dt * dYu3
    Yd4 = Yd + dt * dYd3
    Ye4 = Ye + dt * dYe3
    Ynu4 = Ynu + dt * dYnu3
    dYu4, dYd4, dYe4, dYnu4 = beta_Yukawas(Yu4, Yd4, Ye4, Ynu4, g1, g2, g3)

    Yu_next  = Yu  + (dt / 6.0) * (dYu1  + 2*dYu2  + 2*dYu3  + dYu4)
    Yd_next  = Yd  + (dt / 6.0) * (dYd1  + 2*dYd2  + 2*dYd3  + dYd4)
    Ye_next  = Ye  + (dt / 6.0) * (dYe1  + 2*dYe2  + 2*dYe3  + dYe4)
    Ynu_next = Ynu + (dt / 6.0) * (dYnu1 + 2*dYnu2 + 2*dYnu3 + dYnu4)

    return Yu_next, Yd_next, Ye_next, Ynu_next

def rk4_step(Yu, Yd, Ye, Ynu, g1, g2, g3, dt):
    """
    Single RK4 step for (Yu,Yd,Ye,Ynu,g1,g2,g3).
    """
    # k1
    dg1_1, dg2_1, dg3_1 = beta_g(g1, g2, g3)
    dYu1, dYd1, dYe1, dYnu1 = beta_Yukawas(Yu, Yd, Ye, Ynu, g1, g2, g3)

    # k2
    g1_2 = g1 + 0.5*dt*dg1_1
    g2_2 = g2 + 0.5*dt*dg2_1
    g3_2 = g3 + 0.5*dt*dg3_1
    Yu2 = Yu + 0.5*dt*dYu1
    Yd2 = Yd + 0.5*dt*dYd1
    Ye2 = Ye + 0.5*dt*dYe1
    Ynu2 = Ynu + 0.5*dt*dYnu1

    dg1_2, dg2_2, dg3_2 = beta_g(g1_2, g2_2, g3_2)
    dYu2, dYd2, dYe2, dYnu2 = beta_Yukawas(Yu2, Yd2, Ye2, Ynu2, g1_2, g2_2, g3_2)

    # k3
    g1_3 = g1 + 0.5*dt*dg1_2
    g2_3 = g2 + 0.5*dt*dg2_2
    g3_3 = g3 + 0.5*dt*dg3_2
    Yu3 = Yu + 0.5*dt*dYu2
    Yd3 = Yd + 0.5*dt*dYd2
    Ye3 = Ye + 0.5*dt*dYe2
    Ynu3 = Ynu + 0.5*dt*dYnu2

    dg1_3, dg2_3, dg3_3 = beta_g(g1_3, g2_3, g3_3)
    dYu3, dYd3, dYe3, dYnu3 = beta_Yukawas(Yu3, Yd3, Ye3, Ynu3, g1_3, g2_3, g3_3)

    # k4
    g1_4 = g1 + dt*dg1_3
    g2_4 = g2 + dt*dg2_3
    g3_4 = g3 + dt*dg3_3
    Yu4 = Yu + dt*dYu3
    Yd4 = Yd + dt*dYd3
    Ye4 = Ye + dt*dYe3
    Ynu4 = Ynu + dt*dYnu3

    dg1_4, dg2_4, dg3_4 = beta_g(g1_4, g2_4, g3_4)
    dYu4, dYd4, dYe4, dYnu4 = beta_Yukawas(Yu4, Yd4, Ye4, Ynu4, g1_4, g2_4, g3_4)

    # Combine
    g1_next = g1 + (dt/6.0)*(dg1_1 + 2*dg1_2 + 2*dg1_3 + dg1_4)
    g2_next = g2 + (dt/6.0)*(dg2_1 + 2*dg2_2 + 2*dg2_3 + dg2_4)
    g3_next = g3 + (dt/6.0)*(dg3_1 + 2*dg3_2 + 2*dg3_3 + dg3_4)

    Yu_next = Yu + (dt/6.0)*(dYu1 + 2*dYu2 + 2*dYu3 + dYu4)
    Yd_next = Yd + (dt/6.0)*(dYd1 + 2*dYd2 + 2*dYd3 + dYd4)
    Ye_next = Ye + (dt/6.0)*(dYe1 + 2*dYe2 + 2*dYe3 + dYe4)
    Ynu_next = Ynu + (dt/6.0)*(dYnu1 + 2*dYnu2 + 2*dYnu3 + dYnu4)

    return Yu_next, Yd_next, Ye_next, Ynu_next, g1_next, g2_next, g3_next

def run_RGE(Yu0, Yd0, Ye0, Ynu0, g1_const, g2_const, g3_const, config):
    """
    Integrate Yukawas from t0 to t1 with fixed gauge couplings g1_const,g2_const,g3_const.
    """
    Yu, Yd, Ye, Ynu = Yu0.copy(), Yd0.copy(), Ye0.copy(), Ynu0.copy()

    g1, g2, g3 = g1_const, g2_const, g3_const  # just for returning, not evolved

    t = config.t0
    step = 0
    while (config.dt < 0 and t > config.t1) or (config.dt > 0 and t < config.t1):
        # If you like, keep a *Yukawa-only* sanity check:
        # check_rge_sanity(Yu, Yd, Ye, Ynu, g1, g2, g3, step_index=step)
        Yu, Yd, Ye, Ynu = rk4_step_yukawas_only(Yu, Yd, Ye, Ynu, g1, g2, g3, config.dt)
        t += config.dt
        step += 1

    return Yu, Yd, Ye, Ynu, g1, g2, g3

def run_RGE2(Yu0, Yd0, Ye0, Ynu0, g1_0, g2_0, g3_0, config):
    """
    Integrate from t0 to t1 using RK4 with step dt.
    """
    Yu, Yd, Ye, Ynu = Yu0.copy(), Yd0.copy(), Ye0.copy(), Ynu0.copy()
    g1, g2, g3 = g1_0, g2_0, g3_0

    t = config.t0
    step = 0
    while (config.dt < 0 and t > config.t1) or (config.dt > 0 and t < config.t1):
        # --- NEW: sanity check before taking a step ---
        check_rge_sanity(Yu, Yd, Ye, Ynu, g1, g2, g3, step_index=step)

        Yu, Yd, Ye, Ynu, g1, g2, g3 = rk4_step(Yu, Yd, Ye, Ynu, g1, g2, g3, config.dt)
        t += config.dt
        step += 1

        # Optional: print norms every N steps to inspect behavior
        # if step % 500 == 0:
        #     print(f"step {step}, t={t:.2f}, ||Yu||={np.linalg.norm(Yu):.3g}, g3={g3:.3g}")

    return Yu, Yd, Ye, Ynu, g1, g2, g3



# ============================================================
# 8. Diagonalization, masses, CKM, PMNS
# ============================================================

def diag_dirac_Y(Y, v):
    """
    Bi-unitary diag: Y = U_L diag(y_i) U_R†.
    Here we use SVD: Y = U Σ V†  => U_L=U, U_R=V.
    """
    U_L, s, U_Rh = svd(Y)
    masses = (v / np.sqrt(2.0)) * s
    return U_L, np.diag(s), U_Rh, masses


def takagi_symmetric(m):
    """
    Takagi decomposition for complex symmetric m:
    m = U diag(m_i) U^T, m_i >= 0.
    Implementation via SVD: m = U Σ V†, for symmetric m we take U=V*.
    """
    U, s, Vh = svd(m)
    # For a perfectly symmetric matrix, Vh ≈ U.T
    # We enforce symmetry by using U only.
    # m ≈ U diag(s) U^T
    return U, s


def diagonalize_all(Yu, Yd, Ye, mnu, v):
    # Dirac sectors
    UuL, Yu_diag, UuR, mu = diag_dirac_Y(Yu, v)
    UdL, Yd_diag, UdR, md = diag_dirac_Y(Yd, v)
    UlL, Ye_diag, UlR, me = diag_dirac_Y(Ye, v)

    # Neutrino Takagi
    U_nu, mnu_vals = takagi_symmetric(mnu)
    mnu_masses = mnu_vals  # already positive

    # CKM and PMNS
    Vckm = UuL.conj().T @ UdL
    Vpmns = UlL.conj().T @ U_nu

    return (mu, md, me, mnu_masses, Vckm, Vpmns)


# ============================================================
# 9. Angles, phases, Δm^2
# ============================================================

def extract_angles_and_phase(V):
    """
    Approximate PDG-style extraction of (θ12, θ23, θ13, δ) in radians.
    Assumes V is ~unitary and phases reasonably PDG-like.
    """
    # Indices: (e,μ,τ) → (0,1,2)
    s13 = np.abs(V[0, 2])
    theta13 = np.arcsin(np.clip(s13, 0.0, 1.0))

    s12 = np.abs(V[0, 1])
    c12 = np.abs(V[0, 0])
    theta12 = np.arctan2(s12, c12)

    s23 = np.abs(V[1, 2])
    c23 = np.abs(V[2, 2])
    theta23 = np.arctan2(s23, c23)

    # Jarlskog
    J = np.imag(V[0, 0] * V[1, 1] * np.conj(V[0, 1]) * np.conj(V[1, 0]))
    denom = (np.sin(2*theta12) * np.sin(2*theta23) *
             np.sin(2*theta13) * np.cos(theta13))
    if np.abs(denom) < 1e-12:
        delta = 0.0
    else:
        x = 8 * J / denom
        x = np.clip(x, -1.0, 1.0)
        delta = np.arcsin(x)

    return theta12, theta23, theta13, delta

def rescale_yukawa_sector(Y, v, m_target_heaviest):
    """
    Rescale a 3x3 Yukawa matrix Y so that its largest mass eigenvalue
    matches m_target_heaviest.
    """
    # Diagonalize once to get current heaviest mass
    U_L, s, U_Rh = np.linalg.svd(Y)
    m_current = (v / np.sqrt(2.0)) * np.max(s)
    if m_current == 0:
        return Y, 1.0
    alpha = m_target_heaviest / m_current
    return alpha * Y, alpha

def neutrino_splittings(mnu_masses):
    """
    Sort masses ascending and then compute
    Δm^2_21 = m2^2 - m1^2,
    Δm^2_31 = m3^2 - m1^2.
    """
    m_sorted = np.sort(mnu_masses)   # [m1,m2,m3] with m1 <= m2 <= m3
    m1, m2, m3 = m_sorted
    dm2_21 = m2**2 - m1**2
    dm2_31 = m3**2 - m1**2
    return dm2_21, dm2_31


# ============================================================
# 10. χ^2 (placeholder)
# ============================================================

def chi2(observed, expected, sigma):
    """
    Simple χ^2 = sum((obs-exp)^2 / sigma^2).
    All args are 1D numpy arrays of same length.
    """
    return np.sum(((observed - expected) / sigma)**2)


# ============================================================
# Main driver: one full run
# ============================================================

def main():
    cfg = Config()

    # 1. Alignment kernel
    K = build_alignment_kernel(cfg.eps, N=9)

    # 2. Proto-matrices
    Yu0, Yd0, Ye0, Ynu0, M0 = generate_proto_matrices(cfg)

    # 3. Apply alignment
    Yu9, Yd9, Ye9, Ynu9, M9 = align_all(K, Yu0, Yd0, Ye0, Ynu0, M0)

    # 4. 9->3 Schur complement
    Yu_eff, Yd_eff, Ye_eff, Ynu_eff = reduce_all_dirac(Yu9, Yd9, Ye9, Ynu9)

    # 5. M_R
    M_R = build_M_R(M9, cfg.Lambda_Maj)

    # 6. Light neutrino mass matrix (at μ0)
    m_nu = seesaw_light_neutrinos(Ynu_eff, M_R, cfg.v)

    # 7. RGE: need initial gauge couplings at μ0.
    # For a first pass you can just pick approximate values or run them up from m_Z.
    # Here we choose some placeholders (they won't be physically accurate).
    g1_0, g2_0, g3_0 = 0.46, 0.63, 0.88  # replace with something better if you like

    Yu_EW, Yd_EW, Ye_EW, Ynu_EW, g1_EW, g2_EW, g3_EW = run_RGE(
        Yu_eff, Yd_eff, Ye_eff, Ynu_eff, g1_0, g2_0, g3_0, cfg
    )

    m_t_target   = 173.0
    m_b_target   = 4.18
    m_tau_target = 1.777

    Yu_EW, alpha_u  = rescale_yukawa_sector(Yu_EW, cfg.v, m_t_target)
    Yd_EW, alpha_d  = rescale_yukawa_sector(Yd_EW, cfg.v, m_b_target)
    Ye_EW, alpha_e  = rescale_yukawa_sector(Ye_EW, cfg.v, m_tau_target)

    # 8. Diagonalize at μ_EW
    mu, md, me, mnu_masses, Vckm, Vpmns = diagonalize_all(
        Yu_EW, Yd_EW, Ye_EW, m_nu, cfg.v
    )

    # 9. Angles, phases, splittings
    th12_q, th23_q, th13_q, delta_q = extract_angles_and_phase(Vckm)
    th12_l, th23_l, th13_l, delta_l = extract_angles_and_phase(Vpmns)
    dm2_21, dm2_31 = neutrino_splittings(mnu_masses)

    # Print a quick summary
    print("Up masses (GeV):   ", mu)
    print("Down masses (GeV): ", md)
    print("Lepton masses (GeV):", me)
    print("Neutrino masses (eV-ish units if v normalized):", mnu_masses)
    print()
    print("CKM angles (rad): θ12, θ23, θ13 =", th12_q, th23_q, th13_q)
    print("CKM δ (rad):", delta_q)
    print()
    print("PMNS angles (rad): θ12, θ23, θ13 =", th12_l, th23_l, th13_l)
    print("PMNS δ (rad):", delta_l)
    print("Δm^2_21, Δm^2_31:", dm2_21, dm2_31)



if __name__ == "__main__":
    main()

"""
Up masses (GeV):    [173.         144.89384508  94.23343557]
Down masses (GeV):  [4.18       2.93929816 2.05044361]
Lepton masses (GeV): [1.777      1.14296264 0.72240297]
Neutrino masses (eV-ish units if v normalized): [5.38710365e-11 1.39753165e-11 1.24948729e-11]

CKM angles (rad): θ12, θ23, θ13 = 0.982820232389264 0.6316849274502765 1.0903920791347799
CKM δ (rad): -0.03853724571445584

PMNS angles (rad): θ12, θ23, θ13 = 1.5063804531846523 0.897614572395708 0.7101432500646634
PMNS δ (rad): -0.4742350962605663
Δm^2_21, Δm^2_31: 3.918762256647952e-23 2.7459667281136067e-21
"""