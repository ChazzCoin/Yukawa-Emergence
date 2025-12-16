import numpy as np
import math


def generation_index(i):
    """Generation index g(i) in {0,1,2} for site i in {0..8}."""
    return i % 3

def build_phase_profile_gen(n0_deg, delta_deg):
    """
    Returns φ_gen[g] for g=0,1,2 in radians:
        φ_gen[g] = (n0 + g*delta) * 2π/360
    """
    phi_gen = []
    for g in range(3):
        angle_deg = n0_deg + g * delta_deg
        phi_gen.append(2.0 * math.pi * angle_deg / 360.0)
    return np.array(phi_gen, dtype=float)

def build_site_phases(phi_gen):
    """
    Given φ_gen[g], build φ_i for i=0..8 via g(i)=i mod 3.
    """
    phi_site = np.zeros(9, dtype=float)
    for i in range(9):
        g = generation_index(i)
        phi_site[i] = phi_gen[g]
    return phi_site

def build_phase_matrix(phi_site):
    """
    P_ij = exp(i(φ_i - φ_j)) on 9x9.
    """
    N = len(phi_site)
    P = np.zeros((N, N), dtype=complex)
    for i in range(N):
        for j in range(N):
            P[i, j] = np.exp(1j * (phi_site[i] - phi_site[j]))
    return P

# =========================
# Basic constants
# =========================

v_HIGGS = 246.0          # GeV
Lambda_Maj = 7.0e13      # overall Majorana scale (GeV)
kappa = 360.0 / 89.0
eps = 1.0 / kappa
LIGHT = slice(0, 3)
HEAVY = slice(3, 9)

# EW-scale gauge couplings (at ~m_Z)
g1_EW, g2_EW, g3_EW = 0.357, 0.652, 1.221  # PDG/Wikipedia values
mu_EW = 173.0                              # reference EW scale (GeV)


# =========================
# Utility functions
# =========================

def random_complex_matrix(shape, rng):
    X = rng.normal(size=shape)
    Y = rng.normal(size=shape)
    return X + 1j * Y

def normalize_by_largest_singular_value(M):
    s = np.linalg.svd(M, compute_uv=False)
    s_max = s[0]
    return M if s_max == 0 else M / s_max

def generation_pattern(eps, exponents):
    a, b, c = exponents
    return np.array([eps**a, eps**b, eps**c], float)

def build_site_scales_from_generations(gen3):
    s = np.zeros(9, float)
    s[[0, 3, 6]] = gen3[0]
    s[[1, 4, 7]] = gen3[1]
    s[[2, 5, 8]] = gen3[2]
    return s

def random_weighted_proto(shape, rng, site_scales):
    N = shape[0]
    var = np.zeros((N, N), float)
    for i in range(N):
        for j in range(N):
            var[i, j] = site_scales[i] * site_scales[j]
    X = rng.normal(scale=np.sqrt(var), size=shape)
    Y = rng.normal(scale=np.sqrt(var), size=shape)
    M = X + 1j * Y
    return normalize_by_largest_singular_value(M)


# =========================
# Alignment kernel and proto matrices
# =========================

def build_alignment_kernel(eps, N=9, d_star=7):
    """
    K_ij = eps^{|i-j|} for 0 < |i-j| != d_star
         = 1 for i=j
         = 0 for |i-j| = d_star
    """
    K = np.zeros((N, N), float)
    for i in range(N):
        for j in range(N):
            d = abs(i - j)
            if d == 0:
                K[i, j] = 1.0
            elif d == d_star:
                K[i, j] = 0.0
            else:
                K[i, j] = eps**d
    return K

def apply_alignment(K, X):
    """Schur (Hadamard) alignment: Φ(X) = K ⊙ X."""
    return K * X
def generate_aligned_proto_matrices(
    seed,
    use_site_hierarchy=True,
    exponents_u=(4, 2, 0),
    exponents_d=(3, 1, 0),
    exponents_e=(4, 2, 0),
    exponents_nu=(1, 0, 0),
    # phase patterns: (n0_deg, delta_deg) for each sector
    phase_u=(0, 2),
    phase_d=(0, 3),
    phase_e=(0, 10),
    phase_nu=(0, 25),
    noise_level=0.05,
):
    """
    Build 'aligned' proto Yukawas and Majorana matrix:

    Y_f^(0) = A_f * (s_i s_j) * P_f_ij * (1 + noise_level * ξ_ij),

    where s_i encode generation exponents and P_f encodes
    a triadic phase gradient pattern fixed by (n0,delta) for each sector.
    """
    rng = np.random.default_rng(seed)

    # --- site-scale magnitudes from exponents (unchanged core idea) ---
    if use_site_hierarchy:
        gen_u = generation_pattern(eps, exponents_u)
        gen_d = generation_pattern(eps, exponents_d)
        gen_e = generation_pattern(eps, exponents_e)
        gen_nu = generation_pattern(eps, exponents_nu)
    else:
        gen_u = gen_d = gen_e = gen_nu = np.array([1.0, 1.0, 1.0])

    s_u = build_site_scales_from_generations(gen_u)
    s_d = build_site_scales_from_generations(gen_d)
    s_e = build_site_scales_from_generations(gen_e)
    s_nu = build_site_scales_from_generations(gen_nu)

    # Outer products of magnitudes
    Mag_u = np.outer(s_u, s_u)
    Mag_d = np.outer(s_d, s_d)
    Mag_e = np.outer(s_e, s_e)
    Mag_nu = np.outer(s_nu, s_nu)

    # --- phase patterns per sector ---
    # up
    phi_gen_u = build_phase_profile_gen(*phase_u)
    phi_site_u = build_site_phases(phi_gen_u)
    P_u = build_phase_matrix(phi_site_u)
    # down
    phi_gen_d = build_phase_profile_gen(*phase_d)
    phi_site_d = build_site_phases(phi_gen_d)
    P_d = build_phase_matrix(phi_site_d)
    # charged lepton
    phi_gen_e = build_phase_profile_gen(*phase_e)
    phi_site_e = build_site_phases(phi_gen_e)
    P_e = build_phase_matrix(phi_site_e)
    # neutrino (Dirac)
    phi_gen_nu = build_phase_profile_gen(*phase_nu)
    phi_site_nu = build_site_phases(phi_gen_nu)
    P_nu = build_phase_matrix(phi_site_nu)

    # --- small complex noise matrices ---
    def small_noise_matrix():
        # complex noise of order 1
        A = rng.normal(size=(9, 9))
        B = rng.normal(size=(9, 9))
        return A + 1j * B

    N_u = small_noise_matrix()
    N_d = small_noise_matrix()
    N_e = small_noise_matrix()
    N_nu = small_noise_matrix()

    # --- build proto Yukawas with alignment-dominated structure ---
    Yu0 = Mag_u * P_u * (1.0 + noise_level * N_u)
    Yd0 = Mag_d * P_d * (1.0 + noise_level * N_d)
    Ye0 = Mag_e * P_e * (1.0 + noise_level * N_e)
    Ynu0 = Mag_nu * P_nu * (1.0 + noise_level * N_nu)

    # Normalize each sector so largest singular value ≈ 1
    Yu0 = normalize_by_largest_singular_value(Yu0)
    Yd0 = normalize_by_largest_singular_value(Yd0)
    Ye0 = normalize_by_largest_singular_value(Ye0)
    Ynu0 = normalize_by_largest_singular_value(Ynu0)

    # --- Majorana proto: can also be aligned, but keep simple for now ---
    M0 = random_complex_matrix((9, 9), rng)
    M0 = normalize_by_largest_singular_value(M0)
    M0 = 0.5 * (M0 + M0.T)  # symmetric Majorana

    return Yu0, Yd0, Ye0, Ynu0, M0

def generate_proto_matrices(
    seed,
    use_site_hierarchy=True,
    exponents_u=(4, 2, 0),
    exponents_d=(3, 1, 0),
    exponents_e=(4, 2, 0),
    exponents_nu=(1, 0, 0),
):
    rng = np.random.default_rng(seed)

    if use_site_hierarchy:
        gen_u = generation_pattern(eps, exponents_u)
        gen_d = generation_pattern(eps, exponents_d)
        gen_e = generation_pattern(eps, exponents_e)
        gen_nu = generation_pattern(eps, exponents_nu)
    else:
        gen_u = gen_d = gen_e = gen_nu = np.array([1.0, 1.0, 1.0])

    site_scales_u = build_site_scales_from_generations(gen_u)
    site_scales_d = build_site_scales_from_generations(gen_d)
    site_scales_e = build_site_scales_from_generations(gen_e)
    site_scales_nu = build_site_scales_from_generations(gen_nu)

    Yu0 = random_weighted_proto((9, 9), rng, site_scales_u)
    Yd0 = random_weighted_proto((9, 9), rng, site_scales_d)
    Ye0 = random_weighted_proto((9, 9), rng, site_scales_e)
    Ynu0 = random_weighted_proto((9, 9), rng, site_scales_nu)

    M0 = random_complex_matrix((9, 9), rng)
    M0 = normalize_by_largest_singular_value(M0)
    M0 = 0.5 * (M0 + M0.T)   # symmetric Majorana proto

    return Yu0, Yd0, Ye0, Ynu0, M0

def align_all(K, Yu0, Yd0, Ye0, Ynu0, M0):
    Yu9 = apply_alignment(K, Yu0)
    Yd9 = apply_alignment(K, Yd0)
    Ye9 = apply_alignment(K, Ye0)
    Ynu9 = apply_alignment(K, Ynu0)
    M9 = apply_alignment(K, M0)
    return Yu9, Yd9, Ye9, Ynu9, M9


# =========================
# Schur complement 9→3
# =========================

def schur_9_to_3(Y9, cond_tol=1e12):
    A = Y9[LIGHT, LIGHT]
    B = Y9[LIGHT, HEAVY]
    D = Y9[HEAVY, HEAVY]

    s = np.linalg.svd(D, compute_uv=False)
    cond = s[0] / max(s[-1], 1e-18)
    if cond > cond_tol:
        D_inv = np.linalg.pinv(D)
        Y_eff = A - B @ D_inv @ B.conj().T
    else:
        X = np.linalg.solve(D, B.conj().T)
        Y_eff = A - B @ X
    return Y_eff

def reduce_all_dirac(Yu9, Yd9, Ye9, Ynu9):
    Yu_eff = schur_9_to_3(Yu9)
    Yd_eff = schur_9_to_3(Yd9)
    Ye_eff = schur_9_to_3(Ye9)
    Ynu_eff = schur_9_to_3(Ynu9)
    return Yu_eff, Yd_eff, Ye_eff, Ynu_eff


# =========================
# Majorana triadic heavy sector and seesaw
# =========================

def heavy_block(M9):
    return M9[HEAVY, HEAVY]

def triad_heavy_basis(Nh=6, ks=(1, 2, 3)):
    i = np.arange(Nh)
    cols = []
    for k in ks:
        v = np.exp(2j * np.pi * k * i / Nh)
        v = v / np.linalg.norm(v)
        cols.append(v)
    return np.stack(cols, axis=1)  # 6 x 3

def build_M_R_triadic(M9_aligned, Lambda_Maj, ks=(1, 2, 3)):
    M_H = heavy_block(M9_aligned)      # 6x6
    B_H = triad_heavy_basis(6, ks)     # 6x3
    M3 = B_H.conj().T @ M_H @ B_H
    M3 = 0.5 * (M3 + M3.T)             # enforce symmetry
    return Lambda_Maj * M3

def seesaw_light_neutrinos(Ynu_eff, M_R, v=v_HIGGS, cond_tol=1e12):
    m_D = (v / math.sqrt(2.0)) * Ynu_eff
    s = np.linalg.svd(M_R, compute_uv=False)
    cond = s[0] / max(s[-1], 1e-18)
    if cond > cond_tol:
        M_R_inv = np.linalg.pinv(M_R)
    else:
        M_R_inv = np.linalg.inv(M_R)
    m_nu = - m_D @ M_R_inv @ m_D.T
    m_nu = 0.5 * (m_nu + m_nu.T)
    return m_nu


# =========================
# Diagonalization and mixing
# =========================

def diag_dirac_Y(Y, v=v_HIGGS):
    U_L, s, U_Rh = np.linalg.svd(Y)
    masses = (v / math.sqrt(2.0)) * s
    return U_L, np.diag(s), U_Rh, masses

def takagi_symmetric(m):
    U, s, Vh = np.linalg.svd(m)
    return U, s

def diagonalize_all(Yu, Yd, Ye, mnu, v=v_HIGGS):
    UuL, Yu_diag, UuR, mu = diag_dirac_Y(Yu, v)
    UdL, Yd_diag, UdR, md = diag_dirac_Y(Yd, v)
    UlL, Ye_diag, UlR, me = diag_dirac_Y(Ye, v)
    U_nu, mnu_vals = takagi_symmetric(mnu)

    Vckm = UuL.conj().T @ UdL
    Vpmns = UlL.conj().T @ U_nu
    return mu, md, me, mnu_vals, Vckm, Vpmns

def extract_angles_and_phase(V):
    s13 = abs(V[0, 2])
    theta13 = math.asin(max(0.0, min(1.0, s13)))

    s12 = abs(V[0, 1])
    c12 = abs(V[0, 0])
    theta12 = math.atan2(s12, c12)

    s23 = abs(V[1, 2])
    c23 = abs(V[2, 2])
    theta23 = math.atan2(s23, c23)

    J = np.imag(V[0, 0] * V[1, 1] * np.conj(V[0, 1]) * np.conj(V[1, 0]))
    denom = (math.sin(2 * theta12) * math.sin(2 * theta23) *
             math.sin(2 * theta13) * math.cos(theta13))

    if abs(denom) < 1e-12:
        delta = 0.0
    else:
        x = 8 * J / denom
        x = max(-1.0, min(1.0, x))
        delta = math.asin(x)

    return theta12, theta23, theta13, delta


# =========================
# Alignment at high scale
# =========================

def run_alignment_high_scale(
    seed=0,
    exponents_u=(4, 2, 0),
    exponents_d=(3, 1, 0),
    exponents_e=(4, 2, 0),
    exponents_nu=(1, 0, 0),
    triad_ks=(1, 2, 3),
    use_site_hierarchy=True,
    # phase patterns can be forwarded from here if you like
):
    K = build_alignment_kernel(eps, N=9, d_star=7)

    Yu0, Yd0, Ye0, Ynu0, M0 = generate_aligned_proto_matrices(
        seed,
        use_site_hierarchy=use_site_hierarchy,
        exponents_u=exponents_u,
        exponents_d=exponents_d,
        exponents_e=exponents_e,
        exponents_nu=exponents_nu,
        # you can pass custom phase_u, phase_d, phase_e, phase_nu here
    )

    Yu9, Yd9, Ye9, Ynu9, M9 = align_all(K, Yu0, Yd0, Ye0, Ynu0, M0)
    Yu_eff, Yd_eff, Ye_eff, Ynu_eff = reduce_all_dirac(Yu9, Yd9, Ye9, Ynu9)
    M_R = build_M_R_triadic(M9, Lambda_Maj, ks=triad_ks)
    mnu = seesaw_light_neutrinos(Ynu_eff, M_R, v_HIGGS)
    return Yu_eff, Yd_eff, Ye_eff, mnu



# =========================
# 1-loop SM RGEs
# =========================

def beta_gauge(g1, g2, g3):
    # 1-loop SM beta for gauge couplings (SU(5)-normalized g1)
    b1, b2, b3 = 41.0 / 6.0, -19.0 / 6.0, -7.0
    factor = 1.0 / (16 * math.pi**2)
    dg1 = factor * b1 * g1**3
    dg2 = factor * b2 * g2**3
    dg3 = factor * b3 * g3**3
    return dg1, dg2, dg3

def beta_yukawas(Yu, Yd, Ye, g1, g2, g3):
    """
    1-loop SM matrix RGEs for Yu,Yd,Ye (Ramond).
    16π² dYu/dt = Yu βu, etc.
    """
    factor = 1.0 / (16 * math.pi**2)

    Hu = Yu.conj().T @ Yu
    Hd = Yd.conj().T @ Yd
    He = Ye.conj().T @ Ye

    T = np.trace(3 * Hu + 3 * Hd + He).real
    I = np.eye(3, dtype=complex)

    cu = (17.0 / 20.0) * g1**2 + (9.0 / 4.0) * g2**2 + 8.0 * g3**2
    cd = (1.0 / 4.0) * g1**2 + (9.0 / 4.0) * g2**2 + 8.0 * g3**2
    ce = (9.0 / 4.0) * (g1**2 + g2**2)

    beta_u_mat = 1.5 * (Hu - Hd) + T * I - cu * I
    beta_d_mat = 1.5 * (Hd - Hu) + T * I - cd * I
    beta_e_mat = 1.5 * He + T * I - ce * I

    dYu = factor * (Yu @ beta_u_mat)
    dYd = factor * (Yd @ beta_d_mat)
    dYe = factor * (Ye @ beta_e_mat)
    return dYu, dYd, dYe

def rge_run(Yu0, Yd0, Ye0, g1_0, g2_0, g3_0, mu_high, mu_low, steps=4000):
    """
    Run Yukawas + gauge couplings from mu_high down to mu_low in t = ln μ.
    Simple RK2 integrator.
    """
    t_high = math.log(mu_high)
    t_low = math.log(mu_low)
    dt = (t_low - t_high) / steps

    Yu, Yd, Ye = Yu0.copy(), Yd0.copy(), Ye0.copy()
    g1, g2, g3 = g1_0, g2_0, g3_0

    for _ in range(steps):
        # First stage
        dYu1, dYd1, dYe1 = beta_yukawas(Yu, Yd, Ye, g1, g2, g3)
        dg1_1, dg2_1, dg3_1 = beta_gauge(g1, g2, g3)

        Yu_mid = Yu + 0.5 * dYu1 * dt
        Yd_mid = Yd + 0.5 * dYd1 * dt
        Ye_mid = Ye + 0.5 * dYe1 * dt
        g1_mid = g1 + 0.5 * dg1_1 * dt
        g2_mid = g2 + 0.5 * dg2_1 * dt
        g3_mid = g3 + 0.5 * dg3_1 * dt

        # Second stage
        dYu2, dYd2, dYe2 = beta_yukawas(Yu_mid, Yd_mid, Ye_mid, g1_mid, g2_mid, g3_mid)
        dg1_2, dg2_2, dg3_2 = beta_gauge(g1_mid, g2_mid, g3_mid)

        Yu += dYu2 * dt
        Yd += dYd2 * dt
        Ye += dYe2 * dt
        g1 += dg1_2 * dt
        g2 += dg2_2 * dt
        g3 += dg3_2 * dt

    return Yu, Yd, Ye, g1, g2, g3

def gauge_run_analytic(g1_EW, g2_EW, g3_EW, mu_EW, mu_high):
    """
    Analytic 1-loop gauge running: 1/g^2(μ) = 1/g^2(μ0) - (2b/16π²) ln(μ/μ0)
    """
    b1, b2, b3 = 41.0 / 6.0, -19.0 / 6.0, -7.0

    def run_one(g0, b):
        L = math.log(mu_high / mu_EW)
        denom = 1.0 / g0**2 - (2 * b / (16 * math.pi**2)) * L
        return math.sqrt(1.0 / denom)

    return (run_one(g1_EW, b1),
            run_one(g2_EW, b2),
            run_one(g3_EW, b3))


# =========================
# Sector-wise rescaling
# =========================

def rescale_yukawa_to_heaviest_mass(Y, target_mass, v=v_HIGGS):
    _, _, _, masses = diag_dirac_Y(Y, v)
    m_max = max(masses)
    if m_max == 0:
        return Y, 1.0
    alpha = target_mass / m_max
    return alpha * Y, alpha

def rescale_neutrino_masses(mnu_matrix, target_m3):
    U, vals = takagi_symmetric(mnu_matrix)
    m3 = max(vals)
    if m3 == 0:
        return mnu_matrix, 1.0
    beta = target_m3 / m3
    return beta * mnu_matrix, beta


# =========================
# Full pipeline: alignment + RGE + rescaling
# =========================

def run_full_pipeline_with_RGE_and_rescaling(
    seed=0,
    mu_high=1.0e14,
    mu_low=mu_EW,
    triad_ks=(1, 2, 3),
    m_t_target=173.0,
    m_b_target=4.18,
    m_tau_target=1.77686,
    m3_nu_target_eV=0.058,
):
    # 1. Alignment at high scale
    Yu_high, Yd_high, Ye_high, mnu_high = run_alignment_high_scale(
        seed=seed,
        triad_ks=triad_ks,
    )

    # 2. Gauge couplings at high scale (from EW → high analytic run)
    g1_high, g2_high, g3_high = gauge_run_analytic(
        g1_EW, g2_EW, g3_EW, mu_EW, mu_high
    )

    # 3. 1-loop RGE down to EW
    Yu_low, Yd_low, Ye_low, g1_low, g2_low, g3_low = rge_run(
        Yu_high, Yd_high, Ye_high,
        g1_high, g2_high, g3_high,
        mu_high, mu_low,
        steps=4000,
    )

    # 4. Sector-wise rescaling of Yukawas
    Yu_res, alpha_u = rescale_yukawa_to_heaviest_mass(Yu_low, m_t_target, v_HIGGS)
    Yd_res, alpha_d = rescale_yukawa_to_heaviest_mass(Yd_low, m_b_target, v_HIGGS)
    Ye_res, alpha_e = rescale_yukawa_to_heaviest_mass(Ye_low, m_tau_target, v_HIGGS)

    # Neutrino mass rescaling: match heaviest eigenvalue to 0.058 eV
    m3_target_GeV = m3_nu_target_eV * 1e-9
    mnu_res, beta_nu = rescale_neutrino_masses(mnu_high, m3_target_GeV)

    # 5. Diagonalize at EW scale
    mu, md, me, mnu_vals, Vckm, Vpmns = diagonalize_all(
        Yu_res, Yd_res, Ye_res, mnu_res, v_HIGGS
    )

    mu_sorted = np.sort(mu)
    md_sorted = np.sort(md)
    me_sorted = np.sort(me)
    mnu_sorted = np.sort(mnu_vals)

    th12_q, th23_q, th13_q, delta_q = extract_angles_and_phase(Vckm)
    th12_l, th23_l, th13_l, delta_l = extract_angles_and_phase(Vpmns)

    return {
        "mu": mu_sorted,
        "md": md_sorted,
        "me": me_sorted,
        "mnu": mnu_sorted,
        "Vckm": Vckm,
        "Vpmns": Vpmns,
        "angles_quark": (th12_q, th23_q, th13_q, delta_q),
        "angles_lepton": (th12_l, th23_l, th13_l, delta_l),
        "alphas": (alpha_u, alpha_d, alpha_e),
        "beta_nu": beta_nu,
        "gauges_low": (g1_low, g2_low, g3_low),
    }


# =========================
# Example run (seed = 0)
# =========================

if __name__ == "__main__":
    res = run_full_pipeline_with_RGE_and_rescaling(seed=0)

    mu, md, me, mnu = res["mu"], res["md"], res["me"], res["mnu"]
    print("=== Masses at EW scale (GeV) ===")
    print("up:",   mu)
    print("down:", md)
    print("lep:",  me)
    print("nu (GeV):", mnu)

    print("\n=== Mass ratios (normalized to heaviest) ===")
    print("up   :", mu / mu[-1])
    print("down :", md / md[-1])
    print("lep  :", me / me[-1])
    print("nu   :", mnu / mnu[-1])

    thq = [math.degrees(x) for x in res["angles_quark"]]
    thl = [math.degrees(x) for x in res["angles_lepton"]]

    print("\n=== Quark mixing angles (deg) ===")
    print("theta12 =", thq[0])
    print("theta23 =", thq[1])
    print("theta13 =", thq[2])
    print("delta_CP (q) =", thq[3])

    print("\n=== Lepton mixing angles (deg) ===")
    print("theta12 =", thl[0])
    print("theta23 =", thl[1])
    print("theta13 =", thl[2])
    print("delta_CP (ℓ) =", thl[3])
