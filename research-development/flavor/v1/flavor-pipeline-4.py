#!/usr/bin/env python3
import numpy as np
from numpy.linalg import svd, solve, cond

# =========================
# Global config
# =========================

class Config:
    v = 246.0                 # GeV
    mu0 = 1.0e12              # GeV
    mu_EW = 91.1876           # GeV
    Lambda_Maj = 1.0e14       # GeV

    # Alignment scale: Fibonacci / 360
    # eps = 89/360, kappa = 360/89
    kappa = 360.0 / 89.0
    eps   = 89.0  / 360.0

    seed = 12345              # overwritten per run

    t0 = np.log(mu0)
    t1 = np.log(mu_EW)
    dt = -0.001               # RG step in log μ (downwards)

    # Higgs quartic (approx EW value, treated constant here)
    lam = 0.13


# =========================
# Alignment kernel K (9x9)
# =========================

def build_alignment_kernel(eps, N=9):
    """
    K_ij = eps^{|i-j|} for |i-j| in D_360^{(9)} = {1,2,3,4,5,6,8},
    K_ij = 0 for |i-j| = 7, K_ii = 1.
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
            else:
                K[i, j] = 0.0
    return K


# =========================
# Proto-matrices
# =========================

def random_complex_matrix(shape, rng):
    real = rng.normal(0.0, 1.0, size=shape)
    imag = rng.normal(0.0, 1.0, size=shape)
    return (real + 1j * imag) / np.sqrt(2.0)

def random_weighted_proto(shape, rng, site_scales):
    """
    Gaussian proto-matrix with site-dependent variances:
      Var[X_ij] ~ site_scales[i] * site_scales[j].
    Then normalized to largest singular value = 1.
    """
    N = shape[0]
    X = np.zeros(shape, dtype=complex)
    for i in range(N):
        for j in range(N):
            s = site_scales[i] * site_scales[j]
            real = rng.normal(0.0, s)
            imag = rng.normal(0.0, s)
            X[i, j] = (real + 1j * imag) / np.sqrt(2.0)
    return normalize_by_largest_singular_value(X)

def build_site_scales_from_generations(gen_pattern):
    """
    gen_pattern: length-3 array [s1, s2, s3] for 'generation' 1,2,3.
    We assign site_scales[i] = gen_pattern[i % 3] for a 9-site chain.
    This enforces triadic repetition: (0,3,6)->s1, (1,4,7)->s2, (2,5,8)->s3.
    """
    gen_pattern = np.array(gen_pattern, dtype=float)
    scales = np.zeros(9, dtype=float)
    for i in range(9):
        g = i % 3
        scales[i] = gen_pattern[g]
    return scales


def normalize_by_largest_singular_value(X):
    s = svd(X, compute_uv=False)
    s_max = np.max(s)
    if s_max == 0:
        return X
    return X / s_max


def generate_proto_matrices(cfg: Config):
    rng = np.random.default_rng(cfg.seed)

    eps = cfg.eps  # 89/360

    # --- sector-dependent generation patterns (gen1, gen2, gen3) ---
    # up-type: strong hierarchy
    gen_u = [eps ** 4, eps ** 2, 1.0]

    # down-type: moderate hierarchy
    gen_d = [eps ** 3, eps, 1.0]

    # charged leptons: similar to down
    gen_e = [eps ** 3, eps, 1.0]

    # neutrino Dirac: weak hierarchy
    gen_nu = [eps, 1.0, 1.0]

    # build 9-site scales with triadic pattern
    site_scales_u  = build_site_scales_from_generations(gen_u)
    site_scales_d  = build_site_scales_from_generations(gen_d)
    site_scales_e  = build_site_scales_from_generations(gen_e)
    site_scales_nu = build_site_scales_from_generations(gen_nu)

    # --- draw weighted proto-matrices ---
    Yu0  = random_weighted_proto((9, 9), rng, site_scales_u)
    Yd0  = random_weighted_proto((9, 9), rng, site_scales_d)
    Ye0  = random_weighted_proto((9, 9), rng, site_scales_e)
    Ynu0 = random_weighted_proto((9, 9), rng, site_scales_nu)

    # Majorana proto: keep it O(1) and symmetric, no extra site hierarchy yet
    M0  = random_complex_matrix((9, 9), rng)
    M0  = normalize_by_largest_singular_value(M0)
    M0  = 0.5 * (M0 + M0.T)

    # optional overall Yukawa scale if you still want it (can set to 1.0 now)
    yukawa_scale = 0.5
    Yu0  *= yukawa_scale
    Yd0  *= yukawa_scale
    Ye0  *= yukawa_scale
    Ynu0 *= yukawa_scale

    return Yu0, Yd0, Ye0, Ynu0, M0

# =========================
# Alignment Φ: K ⊙ X
# =========================

def apply_alignment(K, X):
    return K * X


def align_all(K, Yu0, Yd0, Ye0, Ynu0, M0):
    Yu9  = apply_alignment(K, Yu0)
    Yd9  = apply_alignment(K, Yd0)
    Ye9  = apply_alignment(K, Ye0)
    Ynu9 = apply_alignment(K, Ynu0)
    M9   = apply_alignment(K, M0)
    return Yu9, Yd9, Ye9, Ynu9, M9


# =========================
# Schur complement 9→3
# =========================

def schur_9_to_3(Y9):
    """
    Y9 is 9x9. Light sites: 0,1,2; heavy: 3..8.
    Y_eff = A - B D^{-1} B†.
    """
    A = Y9[0:3, 0:3]
    B = Y9[0:3, 3:9]
    D = Y9[3:9, 3:9]

    if cond(D) > 1e12:
        # could resample here if you want
        pass

    X = solve(D, B.conj().T)      # D X = B†
    BDinvBdag = B @ X
    Y_eff = A - BDinvBdag
    return Y_eff


def reduce_all_dirac(Yu9, Yd9, Ye9, Ynu9):
    Yu_eff  = schur_9_to_3(Yu9)
    Yd_eff  = schur_9_to_3(Yd9)
    Ye_eff  = schur_9_to_3(Ye9)
    Ynu_eff = schur_9_to_3(Ynu9)
    return Yu_eff, Yd_eff, Ye_eff, Ynu_eff


# =========================
# Majorana sector: triadic projection 6→3
# =========================

def heavy_block(M9):
    """
    Extract 6x6 heavy block (sites 3..8, 0-based).
    """
    return M9[3:9, 3:9]


def triad_heavy_basis(Nh=6):
    """
    Build a 6x3 triadic basis in heavy space using DFT modes k = 1,2,3.
    Columns are normalized.
    """
    ks = np.array([1, 2, 3])
    i = np.arange(Nh)
    basis = []
    for k in ks:
        vec = np.exp(2j * np.pi * k * i / Nh)
        vec /= np.linalg.norm(vec)
        basis.append(vec)
    # shape (Nh, 3)
    return np.stack(basis, axis=1)


def build_M_R_triadic(M9_aligned, Lambda_Maj):
    """
    9x9 aligned Majorana → 6x6 heavy block → triadic 3x3 projection.
    """
    M_H = heavy_block(M9_aligned)   # 6x6
    B_H = triad_heavy_basis(6)      # 6x3
    # Project to 3x3 heavy triadic block
    M3 = B_H.conj().T @ M_H @ B_H   # 3x3
    M3 = 0.5 * (M3 + M3.T)          # enforce symmetry
    M_R = Lambda_Maj * M3
    return M_R


def seesaw_light_neutrinos(Ynu_eff, M_R, v):
    """
    m_D = v/√2 Ynu_eff,
    m_ν = - m_D M_R^{-1} m_D^T.
    """
    m_D = (v / np.sqrt(2.0)) * Ynu_eff
    X = solve(M_R, m_D.T)
    m_nu = - m_D @ X
    m_nu = 0.5 * (m_nu + m_nu.T)
    return m_nu


# =========================
# 1-loop Yukawa RGEs (g frozen)
# + Weinberg operator RGE
# =========================

def beta_Yukawas(Yu, Yd, Ye, Ynu, g1, g2, g3):
    Yu_dagYu   = Yu.conj().T  @ Yu
    Yd_dagYd   = Yd.conj().T  @ Yd
    Ye_dagYe   = Ye.conj().T  @ Ye
    Ynu_dagYnu = Ynu.conj().T @ Ynu

    T = np.trace(3*Yu_dagYu + 3*Yd_dagYd + Ye_dagYe)

    factor_u  = T - (17/20*g1**2 + 9/4*g2**2 + 8*g3**2)
    factor_d  = T - ( 1/4*g1**2 + 9/4*g2**2 + 8*g3**2)
    factor_e  = T - ( 9/4*g1**2 + 9/4*g2**2)
    factor_nu = T - (9/20*g1**2 + 9/4*g2**2)

    dYu  = Yu  * factor_u  + (3/2)*(Yu  @ Yu_dagYu  - Yd  @ (Yd_dagYd  @ Yu))
    dYd  = Yd  * factor_d  + (3/2)*(Yd  @ Yd_dagYd  - Yu  @ (Yu_dagYu  @ Yd))
    dYe  = Ye  * factor_e  + (3/2)*(Ye  @ Ye_dagYe)
    dYnu = Ynu * factor_nu + (3/2)*(Ynu @ Ynu_dagYnu - Ye @ (Ye_dagYe @ Ynu))

    dYu  /= (16*np.pi**2)
    dYd  /= (16*np.pi**2)
    dYe  /= (16*np.pi**2)
    dYnu /= (16*np.pi**2)

    return dYu, dYd, dYe, dYnu


def beta_kappa_L(kappa_L, Yu, Ye, g2, lam):
    """
    16π² dκ_L/dt = (-3 g2² + 2λ + 6 Tr(Yu†Yu)) κ_L
                   - 3/2 (Ye†Ye κ_L + κ_L (Ye†Ye)^T).
    We treat λ as constant, and ignore g1,g3 in this operator.
    """
    Yu_dagYu = Yu.conj().T @ Yu
    Ye_dagYe = Ye.conj().T @ Ye
    T_u = np.trace(Yu_dagYu)

    pref = (-3*g2**2 + 2*lam + 6*T_u)
    term1 = pref * kappa_L
    term2 = -1.5 * (Ye_dagYe @ kappa_L + kappa_L @ Ye_dagYe.T.conj())

    dkappa = (term1 + term2) / (16*np.pi**2)
    return dkappa


def rk4_step_full(Yu, Yd, Ye, Ynu, kappa_L, g1, g2, g3, lam, dt):
    """
    RK4 step evolving Yukawas + κ_L with fixed (g1,g2,g3,lam).
    """

    # k1
    dYu1, dYd1, dYe1, dYnu1 = beta_Yukawas(Yu, Yd, Ye, Ynu, g1, g2, g3)
    dkL1 = beta_kappa_L(kappa_L, Yu, Ye, g2, lam)

    # k2
    Yu2  = Yu  + 0.5*dt*dYu1
    Yd2  = Yd  + 0.5*dt*dYd1
    Ye2  = Ye  + 0.5*dt*dYe1
    Ynu2 = Ynu + 0.5*dt*dYnu1
    kL2  = kappa_L + 0.5*dt*dkL1

    dYu2, dYd2, dYe2, dYnu2 = beta_Yukawas(Yu2, Yd2, Ye2, Ynu2, g1, g2, g3)
    dkL2 = beta_kappa_L(kL2, Yu2, Ye2, g2, lam)

    # k3
    Yu3  = Yu  + 0.5*dt*dYu2
    Yd3  = Yd  + 0.5*dt*dYd2
    Ye3  = Ye  + 0.5*dt*dYe2
    Ynu3 = Ynu + 0.5*dt*dYnu2
    kL3  = kappa_L + 0.5*dt*dkL2

    dYu3, dYd3, dYe3, dYnu3 = beta_Yukawas(Yu3, Yd3, Ye3, Ynu3, g1, g2, g3)
    dkL3 = beta_kappa_L(kL3, Yu3, Ye3, g2, lam)

    # k4
    Yu4  = Yu  + dt*dYu3
    Yd4  = Yd  + dt*dYd3
    Ye4  = Ye  + dt*dYe3
    Ynu4 = Ynu + dt*dYnu3
    kL4  = kappa_L + dt*dkL3

    dYu4, dYd4, dYe4, dYnu4 = beta_Yukawas(Yu4, Yd4, Ye4, Ynu4, g1, g2, g3)
    dkL4 = beta_kappa_L(kL4, Yu4, Ye4, g2, lam)

    Yu_next  = Yu  + (dt/6.0)*(dYu1  + 2*dYu2  + 2*dYu3  + dYu4)
    Yd_next  = Yd  + (dt/6.0)*(dYd1  + 2*dYd2  + 2*dYd3  + dYd4)
    Ye_next  = Ye  + (dt/6.0)*(dYe1  + 2*dYe2  + 2*dYe3  + dYe4)
    Ynu_next = Ynu + (dt/6.0)*(dYnu1 + 2*dYnu2 + 2*dYnu3 + dYnu4)
    kL_next  = kappa_L + (dt/6.0)*(dkL1 + 2*dkL2 + 2*dkL3 + dkL4)

    return Yu_next, Yd_next, Ye_next, Ynu_next, kL_next


def run_RGE_full(Yu0, Yd0, Ye0, Ynu0, kappa_L0, g1_const, g2_const, g3_const, cfg: Config):
    Yu, Yd, Ye, Ynu = Yu0.copy(), Yd0.copy(), Ye0.copy(), Ynu0.copy()
    kappa_L = kappa_L0.copy()
    g1, g2, g3 = g1_const, g2_const, g3_const
    lam = cfg.lam

    t = cfg.t0
    while (cfg.dt < 0 and t > cfg.t1) or (cfg.dt > 0 and t < cfg.t1):
        Yu, Yd, Ye, Ynu, kappa_L = rk4_step_full(
            Yu, Yd, Ye, Ynu, kappa_L, g1, g2, g3, lam, cfg.dt
        )
        t += cfg.dt

    return Yu, Yd, Ye, Ynu, kappa_L, g1, g2, g3


# =========================
# Diagonalization and angles
# =========================

def diag_dirac_Y(Y, v):
    U_L, s, U_Rh = svd(Y)
    masses = (v / np.sqrt(2.0)) * s
    return U_L, np.diag(s), U_Rh, masses


def takagi_symmetric(m):
    U, s, Vh = svd(m)
    return U, s


def diagonalize_all(Yu, Yd, Ye, mnu, v):
    UuL, Yu_diag, UuR, mu = diag_dirac_Y(Yu, v)
    UdL, Yd_diag, UdR, md = diag_dirac_Y(Yd, v)
    UlL, Ye_diag, UlR, me = diag_dirac_Y(Ye, v)

    U_nu, mnu_vals = takagi_symmetric(mnu)
    mnu_masses = mnu_vals

    Vckm  = UuL.conj().T @ UdL
    Vpmns = UlL.conj().T @ U_nu

    return mu, md, me, mnu_masses, Vckm, Vpmns


def extract_angles_and_phase(V):
    s13 = np.abs(V[0, 2])
    theta13 = np.arcsin(np.clip(s13, 0.0, 1.0))

    s12 = np.abs(V[0, 1])
    c12 = np.abs(V[0, 0])
    theta12 = np.arctan2(s12, c12)

    s23 = np.abs(V[1, 2])
    c23 = np.abs(V[2, 2])
    theta23 = np.arctan2(s23, c23)

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


def neutrino_splittings(mnu_masses):
    m_sorted = np.sort(mnu_masses)
    m1, m2, m3 = m_sorted
    dm2_21 = m2**2 - m1**2
    dm2_31 = m3**2 - m1**2
    return dm2_21, dm2_31   # GeV^2


# =========================
# χ² and observables
# =========================

def chi2(observed, expected, sigma):
    return np.sum(((observed - expected) / sigma)**2)


def rescale_yukawa_sector(Y, v, m_target_heaviest):
    U_L, s, U_Rh = svd(Y)
    m_current = (v / np.sqrt(2.0)) * np.max(s)
    if m_current == 0:
        return Y, 1.0
    alpha = m_target_heaviest / m_current
    return alpha * Y, alpha


# experimental targets (rough)
x_exp = np.array([
    # mass ratios
    0.007,    # m_c/m_t
    1e-5,     # m_u/m_t
    0.02,     # m_s/m_b
    0.001,    # m_d/m_b
    0.06,     # m_mu/m_tau
    0.0003,   # m_e/m_tau
    # CKM angles (rad)
    0.226, 0.041, 0.0035,
    # PMNS angles (rad)
    0.59, 0.84, 0.15,
    # Δm² (eV²)
    7.4e-5, 2.5e-3
])

sigma = np.array([
    0.5*x_exp[0], 0.5*x_exp[1], 0.5*x_exp[2], 0.5*x_exp[3],
    0.5*x_exp[4], 0.5*x_exp[5],
    0.1*x_exp[6], 0.1*x_exp[7], 0.1*x_exp[8],
    0.1*x_exp[9], 0.1*x_exp[10], 0.1*x_exp[11],
    0.3*x_exp[12], 0.3*x_exp[13]
])

def make_observables(res):
    mu, md, me, mnu = res["mu"], res["md"], res["me"], res["mnu"]
    th12_q, th23_q, th13_q = res["th_q"]
    th12_l, th23_l, th13_l = res["th_l"]
    dm2_21, dm2_31 = res["dm2_eV2"]

    # sort ascending so index 2 is heaviest
    mu_sorted = np.sort(mu)
    md_sorted = np.sort(md)
    me_sorted = np.sort(me)

    obs = []

    # mass ratios
    obs.append(mu_sorted[1] / mu_sorted[2])   # m_c/m_t
    obs.append(mu_sorted[0] / mu_sorted[2])   # m_u/m_t
    obs.append(md_sorted[1] / md_sorted[2])   # m_s/m_b
    obs.append(md_sorted[0] / md_sorted[2])   # m_d/m_b
    obs.append(me_sorted[1] / me_sorted[2])   # m_mu/m_tau
    obs.append(me_sorted[0] / me_sorted[2])   # m_e/m_tau

    # CKM
    obs.append(th12_q)
    obs.append(th23_q)
    obs.append(th13_q)

    # PMNS
    obs.append(th12_l)
    obs.append(th23_l)
    obs.append(th13_l)

    # neutrino splittings (eV²)
    obs.append(dm2_21)
    obs.append(dm2_31)

    return np.array(obs)


def chi2_from_res(res):
    x_th = make_observables(res)
    return chi2(x_th, x_exp, sigma)


# =========================
# run_pipeline
# =========================

def run_pipeline(seed, cfg: Config, use_RGE=True):
    cfg.seed = seed

    # 1. kernel
    K = build_alignment_kernel(cfg.eps, N=9)

    # 2. proto
    Yu0, Yd0, Ye0, Ynu0, M0 = generate_proto_matrices(cfg)

    # 3. alignment
    Yu9, Yd9, Ye9, Ynu9, M9 = align_all(K, Yu0, Yd0, Ye0, Ynu0, M0)

    # 4. Schur
    Yu_eff, Yd_eff, Ye_eff, Ynu_eff = reduce_all_dirac(Yu9, Yd9, Ye9, Ynu9)

    # 5. M_R (triadic heavy projection)
    M_R = build_M_R_triadic(M9, cfg.Lambda_Maj)

    # 6. seesaw at μ0 → mν(μ0)
    m_nu_0 = seesaw_light_neutrinos(Ynu_eff, M_R, cfg.v)

    # 6b. Weinberg operator κ_L(μ0)
    kappa_L_0 = (2.0 / cfg.v**2) * m_nu_0

    # 7. RG
    g1_0, g2_0, g3_0 = 0.46, 0.63, 0.88
    if use_RGE:
        Yu_EW, Yd_EW, Ye_EW, Ynu_EW, kappa_L_EW, g1_EW, g2_EW, g3_EW = run_RGE_full(
            Yu_eff, Yd_eff, Ye_eff, Ynu_eff, kappa_L_0, g1_0, g2_0, g3_0, cfg
        )
        m_nu_EW = 0.5 * cfg.v**2 * kappa_L_EW
    else:
        Yu_EW, Yd_EW, Ye_EW, Ynu_EW = Yu_eff, Yd_eff, Ye_eff, Ynu_eff
        m_nu_EW = m_nu_0
        g1_EW, g2_EW, g3_EW = g1_0, g2_0, g3_0

    # 7b. rescale sectors to fix heavy masses
    m_t_target   = 173.0
    m_b_target   = 4.18
    m_tau_target = 1.777

    Yu_EW, alpha_u = rescale_yukawa_sector(Yu_EW, cfg.v, m_t_target)
    Yd_EW, alpha_d = rescale_yukawa_sector(Yd_EW, cfg.v, m_b_target)
    Ye_EW, alpha_e = rescale_yukawa_sector(Ye_EW, cfg.v, m_tau_target)

    # 8. diag at μ_EW
    mu, md, me, mnu_masses, Vckm, Vpmns = diagonalize_all(
        Yu_EW, Yd_EW, Ye_EW, m_nu_EW, cfg.v
    )

    # 9. angles, Δm²
    th12_q, th23_q, th13_q, delta_q = extract_angles_and_phase(Vckm)
    th12_l, th23_l, th13_l, delta_l = extract_angles_and_phase(Vpmns)
    dm2_21_GeV2, dm2_31_GeV2 = neutrino_splittings(mnu_masses)
    dm2_21_eV2 = dm2_21_GeV2 * 1e18
    dm2_31_eV2 = dm2_31_GeV2 * 1e18

    res = {
        "mu": mu,
        "md": md,
        "me": me,
        "mnu": mnu_masses,
        "Vckm": Vckm,
        "Vpmns": Vpmns,
        "th_q": (th12_q, th23_q, th13_q),
        "delta_q": delta_q,
        "th_l": (th12_l, th23_l, th13_l),
        "delta_l": delta_l,
        "dm2_GeV2": (dm2_21_GeV2, dm2_31_GeV2),
        "dm2_eV2": (dm2_21_eV2, dm2_31_eV2),
        "alphas": (alpha_u, alpha_d, alpha_e),
        "g_EW": (g1_EW, g2_EW, g3_EW),
    }

    res["chi2"] = chi2_from_res(res)
    return res


# =========================
# Scan driver
# =========================

if __name__ == "__main__":
    cfg = Config()
    N_seeds = 10

    all_results = []
    chi2_vals = []

    for seed in range(N_seeds):
        r = run_pipeline(seed, cfg, use_RGE=True)
        all_results.append(r)
        chi2_vals.append(r["chi2"])
        print(f"seed {seed}: chi2 = {r['chi2']:.3g}")

    best_idx = int(np.argmin(chi2_vals))
    best = all_results[best_idx]

    print("\nBest seed:", best_idx)
    print("chi2 =", best["chi2"])
    print("Up masses (GeV):   ", best["mu"])
    print("Down masses (GeV): ", best["md"])
    print("Lepton masses (GeV):", best["me"])
    print("Neutrino masses (GeV):", best["mnu"])
    print("Δm² (eV²):", best["dm2_eV2"])
    print("CKM angles (rad):", best["th_q"], "δq:", best["delta_q"])
    print("PMNS angles (rad):", best["th_l"], "δℓ:", best["delta_l"])

"""
seed 0: chi2 = 738
seed 1: chi2 = 1.54e+04
seed 2: chi2 = 1.58e+03
seed 3: chi2 = 517
seed 4: chi2 = 4.95e+05
seed 5: chi2 = 1.08e+03
seed 6: chi2 = 446
seed 7: chi2 = 1.13e+03
seed 8: chi2 = 1.14e+03
seed 9: chi2 = 1.49e+03

Best seed: 6
chi2 = 445.65575358004185
Up masses (GeV):    [1.73000000e+02 3.84363311e-01 2.19678249e-03]
Down masses (GeV):  [4.18000000e+00 5.13651859e-01 6.72225340e-04]
Lepton masses (GeV): [1.77700000e+00 4.09707294e-02 3.63836742e-04]
Neutrino masses (GeV): [1.07549491e-10 3.77417927e-12 4.75913293e-13]
Δm² (eV²): (np.float64(1.4017935683436117e-05), np.float64(0.011566666624029913))
CKM angles (rad): (np.float64(0.11779934551675526), np.float64(0.02453724461291187), np.float64(0.0017602388273065277)) δq: -0.12974536390414831
PMNS angles (rad): (np.float64(0.4118290886093317), np.float64(0.1787909573681276), np.float64(0.049780221388775726)) δℓ: -1.4416353033899119

"""