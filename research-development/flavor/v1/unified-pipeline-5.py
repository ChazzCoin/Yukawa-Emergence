import numpy as np
import math

# ======================================================
# Basic constants and configuration
# ======================================================
# Default context cycle: universal = Earth
DEFAULT_N_EFF = 360
v_HIGGS = 246.0          # GeV
Lambda_Maj = 7.0e13      # heavy RH neutrino scale (GeV)
kappa = 360.0 / 89.0
eps = 1.0 / kappa
LIGHT = slice(0, 3)
HEAVY = slice(3, 9)

# EW-scale gauge couplings (approx at m_Z)
g1_EW, g2_EW, g3_EW = 0.357, 0.652, 1.221
mu_EW = 173.0  # reference "EW" scale (GeV)

# Default phase patterns (deg): (n0, delta)
DEFAULT_PHASE_U = (0, 2)
DEFAULT_PHASE_D = (0, 3)
DEFAULT_PHASE_E = (0, 10)
DEFAULT_PHASE_NU = (0, 25)
DEFAULT_NOISE_LEVEL = 0.05
# ================================================
# Quark-only misalignment (X_q^2)
# ================================================

# Target quark-sector observables (can refine from PDG later)
x_exp_quark = np.array([
    # mass ratios
    0.007,    # m_c/m_t
    1e-5,     # m_u/m_t
    0.02,     # m_s/m_b
    0.001,    # m_d/m_b
    # CKM angles (rad)
    0.226, 0.041, 0.0035,
])

sigma_quark = np.array([
    0.5 * x_exp_quark[0],
    0.5 * x_exp_quark[1],
    0.5 * x_exp_quark[2],
    0.5 * x_exp_quark[3],
    0.1 * x_exp_quark[4],
    0.1 * x_exp_quark[5],
    0.1 * x_exp_quark[6],
])


def make_observables_quark(res):
    """
    Extract quark-only observables from a full pipeline result:
      - 4 mass ratios (up/down)
      - 3 CKM angles (rad)
    """
    mu_vals = res["mu"]
    md_vals = res["md"]
    th12_q, th23_q, th13_q, _ = res["th_q"]

    mu_sorted = np.sort(mu_vals)
    md_sorted = np.sort(md_vals)

    obs = []
    # mass ratios
    obs.append(mu_sorted[1] / mu_sorted[2])  # m_c/m_t
    obs.append(mu_sorted[0] / mu_sorted[2])  # m_u/m_t
    obs.append(md_sorted[1] / md_sorted[2])  # m_s/m_b
    obs.append(md_sorted[0] / md_sorted[2])  # m_d/m_b
    # CKM angles (rad)
    obs.append(th12_q)
    obs.append(th23_q)
    obs.append(th13_q)

    return np.array(obs)


def chi2_from_res_quark(res):
    """
    Quark-only misalignment energy X_q^2 = χ_q^2.
    """
    obs = make_observables_quark(res)
    return chi2(obs, x_exp_quark, sigma_quark)

# ================================================
# Minimal inverse problem: up+down sectors only
# ================================================

def build_target_Yukawas_from_masses(
    m_u, m_c, m_t,
    m_d, m_s, m_b,
    v=v_HIGGS
):
    """
    Build 3x3 target Yukawa matrices Y_u^target, Y_d^target
    from given quark masses at the EW scale.
    For this minimal version, we ignore CKM and take them diagonal.

        Y_f^target = diag(√2 m_i / v)
    """
    yu = np.array([m_u, m_c, m_t], dtype=float)
    yd = np.array([m_d, m_s, m_b], dtype=float)

    Yu_target = np.diag(math.sqrt(2.0) * yu / v)
    Yd_target = np.diag(math.sqrt(2.0) * yd / v)
    return Yu_target, Yd_target
def embed_target_into_9x9_trivial(Y_target_3x3, alpha_heavy=1.0):
    """
    Embed a 3x3 target Yukawa matrix into a 9x9 block structure
    such that the Schur complement (over heavy indices 3..8)
    returns exactly Y_target_3x3.

    Construction:
      A = Y_target_3x3   (3x3)
      B = 0_(3x6)
      C = 0_(6x3)
      D = alpha_heavy * I_6

    Then Y_eff = A - B D^{-1} C = A.
    """
    A = Y_target_3x3
    B = np.zeros((3, 6), dtype=complex)
    C = np.zeros((6, 3), dtype=complex)
    D = alpha_heavy * np.eye(6, dtype=complex)

    # assemble 9x9
    Y9 = np.zeros((9, 9), dtype=complex)
    Y9[0:3, 0:3] = A
    Y9[0:3, 3:9] = B
    Y9[3:9, 0:3] = C
    Y9[3:9, 3:9] = D
    return Y9
def test_trivial_inverse_up_down(
    Yu_target_3, Yd_target_3,
    alpha_heavy=1.0,
    verbose=True
):
    """
    Minimal consistency check:
    - Embed Yu_target_3, Yd_target_3 into 9x9 trivial Yukawas.
    - Apply schur_9_to_3.
    - Compare recovered Y_eff to the target.

    This verifies that, at least in a trivial embedding, the 9x9
    Schur machinery can reproduce a known 3x3 flavor structure exactly.
    """
    Yu_9 = embed_target_into_9x9_trivial(Yu_target_3, alpha_heavy=alpha_heavy)
    Yd_9 = embed_target_into_9x9_trivial(Yd_target_3, alpha_heavy=alpha_heavy)

    Yu_eff = schur_9_to_3(Yu_9)
    Yd_eff = schur_9_to_3(Yd_9)

    diff_u = Yu_eff - Yu_target_3
    diff_d = Yd_eff - Yd_target_3

    norm_u = np.linalg.norm(diff_u)
    norm_d = np.linalg.norm(diff_d)

    if verbose:
        print("=== Minimal inverse test: up/down ===")
        print("Yu_target_3 =\n", Yu_target_3)
        print("Yu_eff (from 9x9) =\n", Yu_eff)
        print("||Yu_eff - Yu_target_3||_F =", norm_u)
        print()
        print("Yd_target_3 =\n", Yd_target_3)
        print("Yd_eff (from 9x9) =\n", Yd_eff)
        print("||Yd_eff - Yd_target_3||_F =", norm_d)

    return norm_u, norm_d

# ======================================================
# Utility: random, normalization, generation exponents
# ======================================================

def random_complex_matrix(shape, rng):
    X = rng.normal(size=shape)
    Y = rng.normal(size=shape)
    return X + 1j * Y

def normalize_by_largest_singular_value(M):
    s = np.linalg.svd(M, compute_uv=False)
    s_max = s[0]
    return M if s_max == 0 else M / s_max

def generation_pattern(eps_local, exponents):
    a, b, c = exponents
    return np.array([eps_local**a, eps_local**b, eps_local**c], float)

def build_site_scales_from_generations(gen3):
    """
    Map generation scales (gen3[0], gen3[1], gen3[2])
    to 9 sites: (0,3,6)->0, (1,4,7)->1, (2,5,8)->2.
    """
    s = np.zeros(9, float)
    s[[0, 3, 6]] = gen3[0]
    s[[1, 4, 7]] = gen3[1]
    s[[2, 5, 8]] = gen3[2]
    return s
# --------------------------------------------------
# Sub-projector: D_360 -> D_{N_eff} for phase sector
# --------------------------------------------------

def project_phase_to_subcycle(n0_base, delta_base, N_eff):
    """
    Implement the sub-projector C_{N_eff} ⊂ C_{360} on phase indices.

    n0_base, delta_base are 'base' integers.
    We enforce that effective phase indices are multiples of q = 360 / N_eff,
    so that Earth phases live on a D_{N_eff} sub-lattice.

        n0_eff    = q * n0_base
        delta_eff = q * delta_base

    Assumes N_eff divides 360.
    """
    if 360 % N_eff != 0:
        raise ValueError(f"N_eff={N_eff} must divide 360.")
    q = 360 // N_eff
    n0_eff    = q * n0_base
    delta_eff = q * delta_base
    return n0_eff, delta_eff
def build_phase_profile_gen_contextual(n0_base, delta_base, N_eff):
    """
    Contextual phase generator:

        1. Project (n0_base, delta_base) through the sub-projector C_{N_eff},
        2. Build φ_gen[g] = (n0_eff + g * delta_eff) * 2π / 360.

    This means universal 360 is still the parent, but we restrict to a
    D_{N_eff} sub-lattice of allowed phases.
    """
    n0_eff, delta_eff = project_phase_to_subcycle(n0_base, delta_base, N_eff)
    phi_gen = []
    for g in range(3):
        angle_deg = n0_eff + g * delta_eff
        phi_gen.append(2.0 * math.pi * angle_deg / 360.0)
    return np.array(phi_gen, dtype=float)

# ======================================================
# Alignment kernel K and Schur alignment
# ======================================================

def build_alignment_kernel(eps_local, N=9, d_star=7):
    """
    K_ij = eps^{|i-j|} for 0<|i-j|!=d_star
         = 1 for i=j
         = 0 for |i-j|=d_star
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
                K[i, j] = eps_local**d
    return K

def apply_alignment(K, X):
    """Schur (Hadamard) alignment: Φ(X) = K ⊙ X."""
    return K * X

def align_all(K, Yu0, Yd0, Ye0, Ynu0, M0):
    Yu9 = apply_alignment(K, Yu0)
    Yd9 = apply_alignment(K, Yd0)
    Ye9 = apply_alignment(K, Ye0)
    Ynu9 = apply_alignment(K, Ynu0)
    M9 = apply_alignment(K, M0)
    return Yu9, Yd9, Ye9, Ynu9, M9

# ======================================================
# Schur complement 9 -> 3 (Dirac sectors)
# ======================================================

def schur_9_to_3(Y9, cond_tol=1e12):
    """
    Light sites: 0..2, heavy: 3..8.
    Y_eff = A - B D^{-1} B† (or pseudo-inverse if ill-conditioned).
    """
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

# ======================================================
# Majorana sector: triadic heavy projection + seesaw
# ======================================================

def heavy_block(M9):
    return M9[HEAVY, HEAVY]   # 6x6

def triad_heavy_basis(Nh=6, ks=(1, 2, 3)):
    """
    Simple DFT-based 6x3 triadic heavy basis.
    """
    i = np.arange(Nh)
    cols = []
    for k in ks:
        v = np.exp(2j * np.pi * k * i / Nh)
        v = v / np.linalg.norm(v)
        cols.append(v)
    return np.stack(cols, axis=1)  # 6 x 3

def build_M_R_triadic(M9_aligned, Lambda_Maj_local, ks=(1, 2, 3)):
    M_H = heavy_block(M9_aligned)      # 6x6
    B_H = triad_heavy_basis(6, ks)     # 6x3
    M3 = B_H.conj().T @ M_H @ B_H
    M3 = 0.5 * (M3 + M3.T)             # enforce symmetry
    return Lambda_Maj_local * M3

def seesaw_light_neutrinos(Ynu_eff, M_R, v=v_HIGGS, cond_tol=1e12):
    """
    m_D = v/sqrt(2) * Ynu_eff
    m_ν = - m_D M_R^{-1} m_D^T
    """
    m_D = (v / math.sqrt(2.0)) * Ynu_eff
    s = np.linalg.svd(M_R, compute_uv=False)
    cond = s[0] / max(s[-1], 1e-18)
    if cond > cond_tol:
        M_R_inv = np.linalg.pinv(M_R)
    else:
        M_R_inv = np.linalg.inv(M_R)
    m_nu = - m_D @ M_R_inv @ m_D.T
    m_nu = 0.5 * (m_nu + m_nu.T)   # enforce symmetry
    return m_nu

# ======================================================
# Diagonalization and mixing
# ======================================================

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
    """
    Approx PDG-like extraction of (θ12,θ23,θ13,δ) from a unitary matrix V.
    """
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

# ======================================================
# Phase gradients and aligned proto Yukawas
# ======================================================

def generation_index(i):
    return i % 3

def build_phase_profile_gen(n0_deg, delta_deg):
    """
    φ_gen[g] = (n0 + g*delta) * 2π/360,  g=0,1,2
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

def generate_aligned_proto_matrices(
    seed,
    use_site_hierarchy=True,
    exponents_u=(4, 2, 0),
    exponents_d=(3, 1, 0),
    exponents_e=(4, 2, 0),
    exponents_nu=(1, 0, 0),
    # 'base' phase patterns: (n0_base, delta_base) for each sector
    phase_u=(0, 2),
    phase_d=(0, 3),
    phase_e=(0, 10),
    phase_nu=(0, 25),
    noise_level=0.05,
    N_eff=DEFAULT_N_EFF,
):
    """
    Build aligned proto Yukawas and Majorana matrix.

    Universal parent: D_360.
    Context sub-projector: D_{N_eff} ⊂ D_360 applied to phase gradients.

    Y_f^(0) = (s_i s_j) * P_f_ij * (1 + noise_level * ξ_ij),
    with P_f built from phases living on the D_{N_eff} sub-lattice.
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

    Mag_u = np.outer(s_u, s_u)
    Mag_d = np.outer(s_d, s_d)
    Mag_e = np.outer(s_e, s_e)
    Mag_nu = np.outer(s_nu, s_nu)

    # --- phase patterns per sector, with sub-projector C_{N_eff} ---
    phi_gen_u  = build_phase_profile_gen_contextual(*phase_u,  N_eff)
    phi_site_u = build_site_phases(phi_gen_u)
    P_u        = build_phase_matrix(phi_site_u)

    phi_gen_d  = build_phase_profile_gen_contextual(*phase_d,  N_eff)
    phi_site_d = build_site_phases(phi_gen_d)
    P_d        = build_phase_matrix(phi_site_d)

    phi_gen_e  = build_phase_profile_gen_contextual(*phase_e,  N_eff)
    phi_site_e = build_site_phases(phi_gen_e)
    P_e        = build_phase_matrix(phi_site_e)

    phi_gen_nu  = build_phase_profile_gen_contextual(*phase_nu,  N_eff)
    phi_site_nu = build_site_phases(phi_gen_nu)
    P_nu        = build_phase_matrix(phi_site_nu)

    # --- small complex noise matrices (unchanged) ---
    def small_noise_matrix():
        A = rng.normal(size=(9, 9))
        B = rng.normal(size=(9, 9))
        return A + 1j * B

    N_u  = small_noise_matrix()
    N_d  = small_noise_matrix()
    N_e  = small_noise_matrix()
    N_nu = small_noise_matrix()

    # --- build proto Yukawas ---
    Yu0  = Mag_u  * P_u  * (1.0 + noise_level * N_u)
    Yd0  = Mag_d  * P_d  * (1.0 + noise_level * N_d)
    Ye0  = Mag_e  * P_e  * (1.0 + noise_level * N_e)
    Ynu0 = Mag_nu * P_nu * (1.0 + noise_level * N_nu)

    Yu0  = normalize_by_largest_singular_value(Yu0)
    Yd0  = normalize_by_largest_singular_value(Yd0)
    Ye0  = normalize_by_largest_singular_value(Ye0)
    Ynu0 = normalize_by_largest_singular_value(Ynu0)

    # Majorana proto (unchanged for now; can also be aligned/sub-projected later)
    M0 = random_complex_matrix((9, 9), rng)
    M0 = normalize_by_largest_singular_value(M0)
    M0 = 0.5 * (M0 + M0.T)

    return Yu0, Yd0, Ye0, Ynu0, M0
# ======================================================
# Alignment at high scale (with explicit phase control)
# ======================================================

def run_alignment_high_scale(
    seed=0,
    exponents_u=(4, 2, 0),
    exponents_d=(3, 1, 0),
    exponents_e=(4, 2, 0),
    exponents_nu=(1, 0, 0),
    triad_ks=(1, 2, 3),
    use_site_hierarchy=True,
    phase_u=(0, 2),
    phase_d=(0, 3),
    phase_e=(0, 10),
    phase_nu=(0, 25),
    noise_level=0.05,
    N_eff=DEFAULT_N_EFF,
):
    K = build_alignment_kernel(eps, N=9, d_star=7)
    Yu0, Yd0, Ye0, Ynu0, M0 = generate_aligned_proto_matrices(
        seed=seed,
        use_site_hierarchy=use_site_hierarchy,
        exponents_u=exponents_u,
        exponents_d=exponents_d,
        exponents_e=exponents_e,
        exponents_nu=exponents_nu,
        phase_u=phase_u,
        phase_d=phase_d,
        phase_e=phase_e,
        phase_nu=phase_nu,
        noise_level=noise_level,
        N_eff=N_eff,
    )
    Yu9, Yd9, Ye9, Ynu9, M9 = align_all(K, Yu0, Yd0, Ye0, Ynu0, M0)
    Yu_eff, Yd_eff, Ye_eff, Ynu_eff = reduce_all_dirac(Yu9, Yd9, Ye9, Ynu9)
    M_R = build_M_R_triadic(M9, Lambda_Maj, ks=triad_ks)
    mnu = seesaw_light_neutrinos(Ynu_eff, M_R, v_HIGGS)
    return Yu_eff, Yd_eff, Ye_eff, mnu


# ======================================================
# 1-loop SM RGEs (gauge + Yukawas)
# ======================================================

def beta_gauge(g1, g2, g3):
    b1, b2, b3 = 41.0 / 6.0, -19.0 / 6.0, -7.0
    factor = 1.0 / (16 * math.pi**2)
    dg1 = factor * b1 * g1**3
    dg2 = factor * b2 * g2**3
    dg3 = factor * b3 * g3**3
    return dg1, dg2, dg3

def beta_yukawas(Yu, Yd, Ye, g1, g2, g3):
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
    t_high = math.log(mu_high)
    t_low = math.log(mu_low)
    dt = (t_low - t_high) / steps

    Yu, Yd, Ye = Yu0.copy(), Yd0.copy(), Ye0.copy()
    g1, g2, g3 = g1_0, g2_0, g3_0

    for _ in range(steps):
        dYu1, dYd1, dYe1 = beta_yukawas(Yu, Yd, Ye, g1, g2, g3)
        dg1_1, dg2_1, dg3_1 = beta_gauge(g1, g2, g3)

        Yu_mid = Yu + 0.5 * dYu1 * dt
        Yd_mid = Yd + 0.5 * dYd1 * dt
        Ye_mid = Ye + 0.5 * dYe1 * dt
        g1_mid = g1 + 0.5 * dg1_1 * dt
        g2_mid = g2 + 0.5 * dg2_1 * dt
        g3_mid = g3 + 0.5 * dg3_1 * dt

        dYu2, dYd2, dYe2 = beta_yukawas(Yu_mid, Yd_mid, Ye_mid, g1_mid, g2_mid, g3_mid)
        dg1_2, dg2_2, dg3_2 = beta_gauge(g1_mid, g2_mid, g3_mid)

        Yu += dYu2 * dt
        Yd += dYd2 * dt
        Ye += dYe2 * dt
        g1 += dg1_2 * dt
        g2 += dg2_2 * dt
        g3 += dg3_2 * dt

    return Yu, Yd, Ye, g1, g2, g3

def gauge_run_analytic(g1_EW_local, g2_EW_local, g3_EW_local, mu_EW_local, mu_high):
    b1, b2, b3 = 41.0 / 6.0, -19.0 / 6.0, -7.0

    def run_one(g0, b):
        L = math.log(mu_high / mu_EW_local)
        denom = 1.0 / g0**2 - (2 * b / (16 * math.pi**2)) * L
        return math.sqrt(1.0 / denom)

    return (run_one(g1_EW_local, b1),
            run_one(g2_EW_local, b2),
            run_one(g3_EW_local, b3))

# ======================================================
# Sector-wise rescaling
# ======================================================

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

# ======================================================
# Full pipeline: alignment + RGE + rescaling
# ======================================================

def run_full_pipeline_with_RGE_and_rescaling(
    seed=0,
    mu_high=1.0e14,
    triad_ks=(1, 2, 3),
    m_t_target=173.0,
    m_b_target=4.18,
    m_tau_target=1.77686,
    m3_nu_target_eV=0.058,
    phase_u=(0, 2),
    phase_d=(0, 3),
    phase_e=(0, 10),
    phase_nu=(0, 25),
    noise_level=0.05,
    N_eff=DEFAULT_N_EFF,
):
    # 1. Alignment at high scale
    Yu_high, Yd_high, Ye_high, mnu_high = run_alignment_high_scale(
        seed=seed,
        triad_ks=triad_ks,
        phase_u=phase_u,
        phase_d=phase_d,
        phase_e=phase_e,
        phase_nu=phase_nu,
        noise_level=noise_level,
        N_eff=N_eff,
    )
    # 2. Gauge couplings at high scale
    g1_high, g2_high, g3_high = gauge_run_analytic(
        g1_EW, g2_EW, g3_EW, mu_EW, mu_high
    )

    # 3. RGE down to EW
    Yu_low, Yd_low, Ye_low, g1_low, g2_low, g3_low = rge_run(
        Yu_high, Yd_high, Ye_high,
        g1_high, g2_high, g3_high,
        mu_high, mu_EW,
        steps=4000,
    )

    # 4. Rescaling
    Yu_res, alpha_u = rescale_yukawa_to_heaviest_mass(Yu_low, m_t_target, v_HIGGS)
    Yd_res, alpha_d = rescale_yukawa_to_heaviest_mass(Yd_low, m_b_target, v_HIGGS)
    Ye_res, alpha_e = rescale_yukawa_to_heaviest_mass(Ye_low, m_tau_target, v_HIGGS)

    m3_target_GeV = m3_nu_target_eV * 1e-9
    mnu_res, beta_nu = rescale_neutrino_masses(mnu_high, m3_target_GeV)

    # 5. Diagonalization at EW scale
    mu_vals, md_vals, me_vals, mnu_vals, Vckm, Vpmns = diagonalize_all(
        Yu_res, Yd_res, Ye_res, mnu_res, v_HIGGS
    )

    return {
        "mu": mu_vals,
        "md": md_vals,
        "me": me_vals,
        "mnu": mnu_vals,
        "Vckm": Vckm,
        "Vpmns": Vpmns,
        "alphas": (alpha_u, alpha_d, alpha_e),
        "beta_nu": beta_nu,
        "gauges_low": (g1_low, g2_low, g3_low),
    }

# ======================================================
# Observables and χ²
# ======================================================

def neutrino_splittings(mnu_masses):
    m_sorted = np.sort(mnu_masses)
    m1, m2, m3 = m_sorted
    dm2_21 = m2**2 - m1**2
    dm2_31 = m3**2 - m1**2
    return dm2_21, dm2_31

def chi2(observed, expected, sigma):
    return np.sum(((observed - expected) / sigma) ** 2)

# Rough experimental targets (same structure as earlier)
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
    0.5 * x_exp[0], 0.5 * x_exp[1], 0.5 * x_exp[2], 0.5 * x_exp[3],
    0.5 * x_exp[4], 0.5 * x_exp[5],
    0.1 * x_exp[6], 0.1 * x_exp[7], 0.1 * x_exp[8],
    0.1 * x_exp[9], 0.1 * x_exp[10], 0.1 * x_exp[11],
    0.3 * x_exp[12], 0.3 * x_exp[13]
])

def make_observables(res):
    mu_vals, md_vals, me_vals, mnu_vals = res["mu"], res["md"], res["me"], res["mnu"]
    th12_q, th23_q, th13_q, _ = res["th_q"]
    th12_l, th23_l, th13_l, _ = res["th_l"]
    dm2_21_eV2, dm2_31_eV2 = res["dm2_eV2"]

    mu_sorted = np.sort(mu_vals)
    md_sorted = np.sort(md_vals)
    me_sorted = np.sort(me_vals)

    obs = []
    # mass ratios
    obs.append(mu_sorted[1] / mu_sorted[2])  # m_c/m_t
    obs.append(mu_sorted[0] / mu_sorted[2])  # m_u/m_t
    obs.append(md_sorted[1] / md_sorted[2])  # m_s/m_b
    obs.append(md_sorted[0] / md_sorted[2])  # m_d/m_b
    obs.append(me_sorted[1] / me_sorted[2])  # m_mu/m_tau
    obs.append(me_sorted[0] / me_sorted[2])  # m_e/m_tau
    # CKM angles (rad)
    obs.append(th12_q)
    obs.append(th23_q)
    obs.append(th13_q)
    # PMNS angles (rad)
    obs.append(th12_l)
    obs.append(th23_l)
    obs.append(th13_l)
    # Δm² (eV²)
    obs.append(dm2_21_eV2)
    obs.append(dm2_31_eV2)

    return np.array(obs)

def chi2_from_res(res):
    obs = make_observables(res)
    return chi2(obs, x_exp, sigma)
def scan_seeds(N_seeds=11,
               mu_high=1e14,
               phase_u=DEFAULT_PHASE_U,
               phase_d=DEFAULT_PHASE_D,
               phase_e=DEFAULT_PHASE_E,
               phase_nu=DEFAULT_PHASE_NU,
               noise_level=DEFAULT_NOISE_LEVEL,
               N_eff=DEFAULT_N_EFF):
    chi2_vals = []
    results = []
    for seed in range(N_seeds):
        res = run_pipeline_for_seed(
            seed,
            mu_high=mu_high,
            phase_u=phase_u,
            phase_d=phase_d,
            phase_e=phase_e,
            phase_nu=phase_nu,
            noise_level=noise_level,
            N_eff=N_eff,           # <-- here
        )
        chi2_vals.append(res["chi2"])
        results.append(res)
        print(f"seed {seed:2d}: chi2 = {res['chi2']:.3g}")
    chi2_vals = np.array(chi2_vals)
    best_idx = int(np.argmin(chi2_vals))
    return best_idx, results[best_idx], chi2_vals, results

def run_pipeline_for_seed(seed,
                          mu_high=1e14,
                          phase_u=DEFAULT_PHASE_U,
                          phase_d=DEFAULT_PHASE_D,
                          phase_e=DEFAULT_PHASE_E,
                          phase_nu=DEFAULT_PHASE_NU,
                          noise_level=DEFAULT_NOISE_LEVEL,
                          N_eff=DEFAULT_N_EFF):
    """
    Wrapper: run the full alignment+RGE+rescaling pipeline for a given seed
    and context cycle N_eff, then extract observables and chi2.
    """
    base = run_full_pipeline_with_RGE_and_rescaling(
        seed=seed,
        mu_high=mu_high,
        phase_u=phase_u,
        phase_d=phase_d,
        phase_e=phase_e,
        phase_nu=phase_nu,
        noise_level=noise_level,
        N_eff=N_eff,          # <-- pass sub-projector context down
    )

    mu_vals  = base["mu"]
    md_vals  = base["md"]
    me_vals  = base["me"]
    mnu_vals = base["mnu"]
    Vckm     = base["Vckm"]
    Vpmns    = base["Vpmns"]

    th12_q, th23_q, th13_q, delta_q = extract_angles_and_phase(Vckm)
    th12_l, th23_l, th13_l, delta_l = extract_angles_and_phase(Vpmns)

    dm2_21_GeV2, dm2_31_GeV2 = neutrino_splittings(mnu_vals)
    dm2_21_eV2 = dm2_21_GeV2 * 1e18
    dm2_31_eV2 = dm2_31_GeV2 * 1e18

    res = {
        "mu": mu_vals,
        "md": md_vals,
        "me": me_vals,
        "mnu": mnu_vals,
        "th_q": (th12_q, th23_q, th13_q, delta_q),
        "th_l": (th12_l, th23_l, th13_l, delta_l),
        "dm2_eV2": (dm2_21_eV2, dm2_31_eV2),
    }
    res["chi2"] = chi2_from_res(res)
    return res
# ================================================
# CKM-based target Yukawas (minimal inverse, up+down)
# ================================================

def ckm_from_angles(theta12, theta23, theta13, delta):
    """
    Build a CKM matrix from (θ12, θ23, θ13, δ) in the standard PDG-like parametrization.
    Angles in radians.
    """
    s12, c12 = math.sin(theta12), math.cos(theta12)
    s23, c23 = math.sin(theta23), math.cos(theta23)
    s13, c13 = math.sin(theta13), math.cos(theta13)

    e_minus_i_delta = math.cos(delta) - 1j * math.sin(delta)

    V = np.zeros((3, 3), dtype=complex)

    # First row
    V[0, 0] = c12 * c13
    V[0, 1] = s12 * c13
    V[0, 2] = s13 * np.conjugate(e_minus_i_delta)

    # Second row
    V[1, 0] = -s12 * c23 - c12 * s23 * s13 * e_minus_i_delta
    V[1, 1] =  c12 * c23 - s12 * s23 * s13 * e_minus_i_delta
    V[1, 2] =  s23 * c13

    # Third row
    V[2, 0] =  s12 * s23 - c12 * c23 * s13 * e_minus_i_delta
    V[2, 1] = -c12 * s23 - s12 * c23 * s13 * e_minus_i_delta
    V[2, 2] =  c23 * c13

    return V

def build_target_Yukawas_with_CKM(
    m_u, m_c, m_t,
    m_d, m_s, m_b,
    theta12, theta23, theta13, delta,
    v=v_HIGGS
):
    """
    Build 3x3 target Yukawas Y_u^target, Y_d^target consistent with
    given quark masses and CKM angles (up to unphysical phases).

    Gauge choice:
      U_L^u = I
      U_L^d = V_CKM
      U_R^u = U_R^d = I

    So:
      Y_u^target = diag(y_u, y_c, y_t)
      Y_d^target = V_CKM @ diag(y_d, y_s, y_b)
    """
    yu = np.array([m_u, m_c, m_t], dtype=float)
    yd = np.array([m_d, m_s, m_b], dtype=float)

    # Yukawa eigenvalues
    y_u_vals = math.sqrt(2.0) * yu / v
    y_d_vals = math.sqrt(2.0) * yd / v

    Yu_target = np.diag(y_u_vals)
    Vckm_target = ckm_from_angles(theta12, theta23, theta13, delta)
    Yd_target = Vckm_target @ np.diag(y_d_vals)

    return Yu_target, Yd_target, Vckm_target
def run_pipeline_for_seed_quark(seed,
                                mu_high=1e14,
                                phase_u=DEFAULT_PHASE_U,
                                phase_d=DEFAULT_PHASE_D,
                                phase_e=DEFAULT_PHASE_E,
                                phase_nu=DEFAULT_PHASE_NU,
                                noise_level=DEFAULT_NOISE_LEVEL,
                                N_eff=DEFAULT_N_EFF):
    """
    Run the full pipeline but return a result object + χ² only for the quark sector.
    """
    base = run_full_pipeline_with_RGE_and_rescaling(
        seed=seed,
        mu_high=mu_high,
        phase_u=phase_u,
        phase_d=phase_d,
        phase_e=phase_e,
        phase_nu=phase_nu,
        noise_level=noise_level,
        N_eff=N_eff,
    )

    mu_vals  = base["mu"]
    md_vals  = base["md"]
    me_vals  = base["me"]
    mnu_vals = base["mnu"]
    Vckm     = base["Vckm"]
    Vpmns    = base["Vpmns"]

    th12_q, th23_q, th13_q, delta_q = extract_angles_and_phase(Vckm)
    th12_l, th23_l, th13_l, delta_l = extract_angles_and_phase(Vpmns)

    dm2_21_GeV2, dm2_31_GeV2 = neutrino_splittings(mnu_vals)
    dm2_21_eV2 = dm2_21_GeV2 * 1e18
    dm2_31_eV2 = dm2_31_GeV2 * 1e18

    res = {
        "mu": mu_vals,
        "md": md_vals,
        "me": me_vals,
        "mnu": mnu_vals,
        "th_q": (th12_q, th23_q, th13_q, delta_q),
        "th_l": (th12_l, th23_l, th13_l, delta_l),
        "dm2_eV2": (dm2_21_eV2, dm2_31_eV2),
    }
    res["chi2_quark"] = chi2_from_res_quark(res)
    return res
def scan_quark_phases_and_seeds(
    phase_u_grid,
    phase_d_grid,
    seeds=range(5),
    mu_high=1e14,
    phase_e=DEFAULT_PHASE_E,
    phase_nu=DEFAULT_PHASE_NU,
    noise_level=DEFAULT_NOISE_LEVEL,
    N_eff=DEFAULT_N_EFF,
):
    """
    Scan over (phase_u, phase_d) pairs and seeds, recording quark-only χ².
    phase_u_grid, phase_d_grid: lists of (n0, delta) tuples.
    """
    best = None
    best_record = None

    for phase_u in phase_u_grid:
        for phase_d in phase_d_grid:
            print(f"\n=== phase_u={phase_u}, phase_d={phase_d}, N_eff={N_eff} ===")
            for seed in seeds:
                res = run_pipeline_for_seed_quark(
                    seed=seed,
                    mu_high=mu_high,
                    phase_u=phase_u,
                    phase_d=phase_d,
                    phase_e=phase_e,
                    phase_nu=phase_nu,
                    noise_level=noise_level,
                    N_eff=N_eff,
                )
                chi2_q = res["chi2_quark"]
                print(f"  seed {seed:2d}: chi2_quark = {chi2_q:.3g}")

                if best is None or chi2_q < best:
                    best = chi2_q
                    best_record = {
                        "seed": seed,
                        "phase_u": phase_u,
                        "phase_d": phase_d,
                        "res": res,
                    }

    print("\n=== Best quark-only configuration found ===")
    print(f"  seed      = {best_record['seed']}")
    print(f"  phase_u   = {best_record['phase_u']}")
    print(f"  phase_d   = {best_record['phase_d']}")
    print(f"  chi2_quark= {best:.3g}")
    return best_record

# ======================================================
# Simple seed scan helper
# ======================================================

def scan_seeds(N_seeds=11,
               mu_high=1e14,
               phase_u=DEFAULT_PHASE_U,
               phase_d=DEFAULT_PHASE_D,
               phase_e=DEFAULT_PHASE_E,
               phase_nu=DEFAULT_PHASE_NU,
               noise_level=DEFAULT_NOISE_LEVEL):
    chi2_vals = []
    results = []
    for seed in range(N_seeds):
        res = run_pipeline_for_seed(seed,
                                    mu_high=mu_high,
                                    phase_u=phase_u,
                                    phase_d=phase_d,
                                    phase_e=phase_e,
                                    phase_nu=phase_nu,
                                    noise_level=noise_level)
        chi2_vals.append(res["chi2"])
        results.append(res)
        print(f"seed {seed:2d}: chi2 = {res['chi2']:.3g}")
    chi2_vals = np.array(chi2_vals)
    best_idx = int(np.argmin(chi2_vals))
    return best_idx, results[best_idx], chi2_vals, results
def test_inverse_with_CKM(
    Yu_target_3, Yd_target_3, Vckm_target,
    alpha_heavy=1.0,
    v=v_HIGGS,
    verbose=True
):
    """
    Minimal inverse test with nontrivial CKM:

    - Embed Yu_target_3, Yd_target_3 into 9x9 trivial Yukawas.
    - Apply schur_9_to_3 to recover Yu_eff, Yd_eff.
    - Diagonalize to get Uu, Ud, and reconstruct Vckm_eff = Uu† Ud.
    - Compare CKM angles of Vckm_eff to those of Vckm_target.
    """
    Yu_9 = embed_target_into_9x9_trivial(Yu_target_3, alpha_heavy=alpha_heavy)
    Yd_9 = embed_target_into_9x9_trivial(Yd_target_3, alpha_heavy=alpha_heavy)

    Yu_eff = schur_9_to_3(Yu_9)
    Yd_eff = schur_9_to_3(Yd_9)

    # Diagonalize effective Yukawas
    UuL, _, _, mu_vals = diag_dirac_Y(Yu_eff, v=v)
    UdL, _, _, md_vals = diag_dirac_Y(Yd_eff, v=v)

    Vckm_eff = UuL.conj().T @ UdL

    # Extract angles
    th12_eff, th23_eff, th13_eff, delta_eff = extract_angles_and_phase(Vckm_eff)
    th12_target, th23_target, th13_target, delta_target = extract_angles_and_phase(Vckm_target)

    if verbose:
        print("=== Inverse CKM test: up/down ===")
        print("Yu_target_3 =\n", Yu_target_3)
        print("Yd_target_3 =\n", Yd_target_3)
        print("\nCKM target:")
        print("  theta12 =", math.degrees(th12_target), "deg")
        print("  theta23 =", math.degrees(th23_target), "deg")
        print("  theta13 =", math.degrees(th13_target), "deg")
        print("  delta   =", math.degrees(delta_target), "deg")

        print("\nCKM from 9x9 Schur + diagonalization:")
        print("  theta12 =", math.degrees(th12_eff), "deg")
        print("  theta23 =", math.degrees(th23_eff), "deg")
        print("  theta13 =", math.degrees(th13_eff), "deg")
        print("  delta   =", math.degrees(delta_eff), "deg")

        print("\nMass eigenvalues (up, down) from Yu_eff, Yd_eff:")
        print("  up   =", np.sort(mu_vals))
        print("  down =", np.sort(md_vals))

    return (th12_eff, th23_eff, th13_eff, delta_eff), (mu_vals, md_vals)

# ======================================================
# Example usage
# ======================================================
N_eff_candidates = [360, 180, 120, 90, 72, 60]

if __name__ == "__main__":
    # ... your existing main code ...

    # Example: run the minimal inverse test for up+down sectors
    # Rough CKM angles (PDG-ish) in radians
    theta12_CKM = 0.226  # ~13°
    theta23_CKM = 0.041  # ~2.35°
    theta13_CKM = 0.0035  # ~0.2°
    delta_CKM = 1.2  # ~69° as a placeholder

    # Quark masses at EW scale (GeV, rough)
    m_u_EW = 0.0022
    m_c_EW = 1.27
    m_t_EW = 173.0

    m_d_EW = 0.0047
    m_s_EW = 0.096
    m_b_EW = 4.18

    print("\n\n=== Minimal inverse CKM test ===")

    # Build target Yukawas with CKM
    Yu_target_3, Yd_target_3, Vckm_target = build_target_Yukawas_with_CKM(
        m_u_EW, m_c_EW, m_t_EW,
        m_d_EW, m_s_EW, m_b_EW,
        theta12_CKM, theta23_CKM, theta13_CKM, delta_CKM,
        v_HIGGS
    )

    (th12_eff, th23_eff, th13_eff, delta_eff), (mu_vals_eff, md_vals_eff) = test_inverse_with_CKM(
        Yu_target_3, Yd_target_3, Vckm_target,
        alpha_heavy=1.0,
        v=v_HIGGS,
        verbose=True
    )

    # Quark-only phase scan (most unified move for CKM + quark masses)
    phase_u_grid = [(0, 1), (0, 2), (0, 3), (5, 2)]
    phase_d_grid = [(0, 3), (0, 4), (0, 5), (5, 3)]

    best_q = scan_quark_phases_and_seeds(
        phase_u_grid=phase_u_grid,
        phase_d_grid=phase_d_grid,
        seeds=range(5),
        mu_high=1e14,
        phase_e=DEFAULT_PHASE_E,  # keep leptons fixed for now
        phase_nu=DEFAULT_PHASE_NU,
        noise_level=0.05,
        N_eff=360,  # full parent cycle
    )
