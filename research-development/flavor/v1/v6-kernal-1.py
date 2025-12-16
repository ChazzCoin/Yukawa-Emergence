import numpy as np
import math

# ==================================
# GEOMETRIC AXIOMS (FIXED)
# ==================================

KAPPA = 360.0 / 89.0
EPS_ALIGN = 1.0 / KAPPA

N_SITES = 9
LIGHT_SITES = [0, 1, 2]
HEAVY_SITES = [3, 4, 5, 6, 7, 8]

def generation_index(i: int) -> int:
    return i % 3

# Fixed triadic exponents per sector
EXP_UP   = (4, 2, 0)
EXP_DOWN = (3, 1, 0)
EXP_LEP  = (4, 2, 0)
EXP_NU   = (1, 0, 0)

def build_site_scales(exponents):
    a, b, c = exponents
    gen_scales = np.array([EPS_ALIGN**a, EPS_ALIGN**b, EPS_ALIGN**c], dtype=float)
    s = np.zeros(N_SITES, dtype=float)
    for i in range(N_SITES):
        s[i] = gen_scales[generation_index(i)]
    return s

def build_kernel():
    K = np.zeros((N_SITES, N_SITES), dtype=float)
    for i in range(N_SITES):
        for j in range(N_SITES):
            if i == j:
                K[i, j] = 1.0
            else:
                d = min(abs(i - j), N_SITES - abs(i - j))
                # For N_SITES = 9, d ∈ {1,2,3,4}; we do NOT impose a forbidden distance
                K[i, j] = EPS_ALIGN ** d
    return K

KERNEL = build_kernel()

# ==================================
# HARMONIC PHASES (D_360, N_eff)
# ==================================

def build_phase_profile_gen(n0_tilde: int, delta_tilde: int, N_eff: int = 360) -> np.ndarray:
    """
    Generation phases implementing D_360 -> D_Neff:
      q = 360 / N_eff,
      n0_eff = q * n0_tilde,
      delta_eff = q * delta_tilde,
      φ_g = (n0_eff + g * delta_eff) * 2π / 360.
    """
    q = 360 // N_eff
    n0_eff = q * n0_tilde
    delta_eff = q * delta_tilde
    phi_gen = np.zeros(3, dtype=float)
    for g in range(3):
        angle_deg = n0_eff + g * delta_eff
        phi_gen[g] = 2.0 * math.pi * angle_deg / 360.0
    return phi_gen

def build_site_phases(phi_gen: np.ndarray) -> np.ndarray:
    phi_site = np.zeros(N_SITES, dtype=float)
    for i in range(N_SITES):
        phi_site[i] = phi_gen[generation_index(i)]
    return phi_site

def build_phase_matrix(phi_site: np.ndarray) -> np.ndarray:
    N = len(phi_site)
    P = np.zeros((N, N), dtype=complex)
    for i in range(N):
        for j in range(N):
            P[i, j] = np.exp(1j * (phi_site[i] - phi_site[j]))
    return P

# ==================================
# Proto Majorana, alignment, Schur
# ==================================

def generate_proto_Majorana(rng):
    M0 = rng.normal(size=(N_SITES, N_SITES)) + 1j * rng.normal(size=(N_SITES, N_SITES))
    M0 = 0.5 * (M0 + M0.T)
    _u, sing, _vh = np.linalg.svd(M0)
    max_sing = np.max(np.abs(sing))
    if max_sing > 0:
        M0 /= max_sing
    return M0

def align_9x9(M: np.ndarray) -> np.ndarray:
    return KERNEL * M

def schur_9_to_3(Y9: np.ndarray) -> np.ndarray:
    ls = LIGHT_SITES
    hs = HEAVY_SITES
    A = Y9[np.ix_(ls, ls)]
    B = Y9[np.ix_(ls, hs)]
    C = Y9[np.ix_(hs, ls)]
    D = Y9[np.ix_(hs, hs)]
    D_inv = np.linalg.pinv(D)
    Y_eff = A - B @ D_inv @ C
    Y_eff = Y_eff + 1e-9 * np.eye(3)
    return Y_eff

def triadic_Majorana_seesaw(M9_aligned: np.ndarray,
                            Ynu_eff: np.ndarray,
                            v: float = 174.0,
                            Lambda_Maj: float = 7e13) -> np.ndarray:
    M_H = M9_aligned[3:9, 3:9]
    h_indices = np.arange(6)
    B_H = np.zeros((6, 3), dtype=complex)
    for col, k in enumerate([1, 2, 3]):
        B_H[:, col] = np.exp(2j * math.pi * k * h_indices / 6.0) / math.sqrt(6.0)
    M_R_dimless = B_H.conj().T @ M_H @ B_H
    M_R_dimless = 0.5 * (M_R_dimless + M_R_dimless.T)
    M_R_dimless += 1e-9 * np.eye(3)
    M_R = Lambda_Maj * M_R_dimless

    m_D = (v / math.sqrt(2.0)) * Ynu_eff
    M_R_inv = np.linalg.inv(M_R)
    m_nu = - m_D @ M_R_inv @ m_D.T
    return m_nu

# ==================================
# RGE, diagonalization, observables
# ==================================

def stub_rge_run(Y_eff: np.ndarray,
                 alpha: float = 0.1,
                 mu_high: float = 1e14,
                 mu_EW: float = 173.0) -> np.ndarray:
    factor = math.log(mu_EW / mu_high) * alpha
    return Y_eff * math.exp(factor)

def diagonalize_dirac(Y: np.ndarray):
    U_L, S, U_Rh = np.linalg.svd(Y)
    return U_L, np.diag(S), U_Rh.conj().T

def diagonalize_majorana(M: np.ndarray):
    H = 0.5 * (M + M.conj().T)
    vals, U = np.linalg.eigh(H)
    idx = np.argsort(np.abs(vals))[::-1]
    vals = vals[idx]
    U = U[:, idx]
    return U, vals

def extract_angles_from_U(U: np.ndarray):
    U = np.array(U, dtype=complex)
    s13 = abs(U[0, 2])
    s13 = min(max(s13, 0.0), 1.0)
    theta13 = math.asin(s13)
    c13 = math.cos(theta13)
    if c13 == 0:
        return 0.0, 0.0, theta13
    s12 = abs(U[0, 1]) / c13
    s12 = min(max(s12, 0.0), 1.0)
    theta12 = math.asin(s12)
    s23 = abs(U[1, 2]) / c13
    s23 = min(max(s23, 0.0), 1.0)
    theta23 = math.asin(s23)
    return theta12, theta23, theta13

exp_targets = {
    'm_c/m_t': 7e-3,
    'm_u/m_t': 1e-5,
    'm_s/m_b': 2e-2,
    'm_d/m_b': 1e-3,
    'm_mu/m_tau': 6e-2,
    'm_e/m_tau': 3e-4,
    'theta12_q': 0.226,
    'theta23_q': 0.041,
    'theta13_q': 0.0035,
    'theta12_l': 0.59,
    'theta23_l': 0.84,
    'theta13_l': 0.15,
    'Delta m2_21': 7.4e-5,
    'Delta m2_31': 2.5e-3,
}

sigma_targets = {
    'm_c/m_t': 0.5 * exp_targets['m_c/m_t'],
    'm_u/m_t': 0.5 * exp_targets['m_u/m_t'],
    'm_s/m_b': 0.5 * exp_targets['m_s/m_b'],
    'm_d/m_b': 0.5 * exp_targets['m_d/m_b'],
    'm_mu/m_tau': 0.5 * exp_targets['m_mu/m_tau'],
    'm_e/m_tau': 0.5 * exp_targets['m_e/m_tau'],
    'theta12_q': 0.1 * exp_targets['theta12_q'],
    'theta23_q': 0.1 * exp_targets['theta23_q'],
    'theta13_q': 0.1 * exp_targets['theta13_q'],
    'theta12_l': 0.1 * exp_targets['theta12_l'],
    'theta23_l': 0.1 * exp_targets['theta23_l'],
    'theta13_l': 0.1 * exp_targets['theta13_l'],
    'Delta m2_21': 0.3 * exp_targets['Delta m2_21'],
    'Delta m2_31': 0.3 * exp_targets['Delta m2_31'],
}

def compute_observables_from_matrices(Yu_EW, Yd_EW, Ye_EW, Mnu_EW):
    Uu_L, Su, _ = diagonalize_dirac(Yu_EW)
    Ud_L, Sd, _ = diagonalize_dirac(Yd_EW)
    Ue_L, Se, _ = diagonalize_dirac(Ye_EW)

    mu_vals = np.diag(Su)
    md_vals = np.diag(Sd)
    me_vals = np.diag(Se)

    mu_sorted = np.sort(np.abs(mu_vals))[::-1]
    md_sorted = np.sort(np.abs(md_vals))[::-1]
    me_sorted = np.sort(np.abs(me_vals))[::-1]

    Vckm = Uu_L.conj().T @ Ud_L
    th12_q, th23_q, th13_q = extract_angles_from_U(Vckm)

    U_nu, mnu_eig = diagonalize_majorana(Mnu_EW)
    mnu_sorted = np.sort(np.abs(mnu_eig))[::-1]

    U_pmns = Ue_L.conj().T @ U_nu
    th12_l, th23_l, th13_l = extract_angles_from_U(U_pmns)

    if mu_sorted[0] != 0:
        mu_sorted *= 173.0 / mu_sorted[0]
    if md_sorted[0] != 0:
        md_sorted *= 4.18 / md_sorted[0]
    if me_sorted[0] != 0:
        me_sorted *= 1.77686 / me_sorted[0]
    if mnu_sorted[0] != 0:
        mnu_sorted *= 0.058 / mnu_sorted[0]

    obs = {}
    obs['m_c/m_t']      = mu_sorted[1] / mu_sorted[0] if mu_sorted[0] != 0 else 0.0
    obs['m_u/m_t']      = mu_sorted[2] / mu_sorted[0] if mu_sorted[0] != 0 else 0.0
    obs['m_s/m_b']      = md_sorted[1] / md_sorted[0] if md_sorted[0] != 0 else 0.0
    obs['m_d/m_b']      = md_sorted[2] / md_sorted[0] if md_sorted[0] != 0 else 0.0
    obs['m_mu/m_tau']   = me_sorted[1] / me_sorted[0] if me_sorted[0] != 0 else 0.0
    obs['m_e/m_tau']    = me_sorted[2] / me_sorted[0] if me_sorted[0] != 0 else 0.0

    obs['theta12_q']    = th12_q
    obs['theta23_q']    = th23_q
    obs['theta13_q']    = th13_q
    obs['theta12_l']    = th12_l
    obs['theta23_l']    = th23_l
    obs['theta13_l']    = th13_l

    mnu_asc = np.sort(mnu_sorted)
    dm21 = mnu_asc[1]**2 - mnu_asc[0]**2
    dm31 = mnu_asc[2]**2 - mnu_asc[0]**2
    obs['Delta m2_21'] = dm21
    obs['Delta m2_31'] = dm31

    return obs, Vckm, U_pmns

def chi2_from_obs(obs):
    chi2 = 0.0
    pulls = {}
    for key, xexp in exp_targets.items():
        xth = obs.get(key, np.nan)
        if not np.isfinite(xth):
            continue
        sig = sigma_targets[key]
        pull = (xth - xexp) / sig
        chi2 += pull**2
        pulls[key] = pull
    return chi2, pulls

# ==================================
# FN-left dressing (best config)
# ==================================

def apply_left_FN_3x3(Y, QL, eps_L):
    QL = np.array(QL, dtype=float)
    F_L = np.diag(eps_L ** QL)
    return F_L @ Y

def apply_left_FN_Majorana_3x3(M, QL, eps_L):
    QL = np.array(QL, dtype=float)
    F_L = np.diag(eps_L ** QL)
    return F_L @ M @ F_L.T

# ==================================
# Generation-phase shift Yukawas (Δ-gen)
# ==================================

def generate_proto_Y_from_gen_shifts(exp_triple,
                                     n0_tilde, delta_tilde,
                                     delta_gen,
                                     N_eff=180):
    delta_gen = np.array(delta_gen, dtype=float)
    phi_gen_target = build_phase_profile_gen(n0_tilde, delta_tilde, N_eff)
    phi_gen = phi_gen_target + delta_gen
    phi_site = build_site_phases(phi_gen)
    P = build_phase_matrix(phi_site)
    s = build_site_scales(exp_triple)
    Mag = np.outer(s, s)
    Y0 = Mag * P
    _u, sing, _vh = np.linalg.svd(Y0)
    max_sing = np.max(np.abs(sing))
    if max_sing > 0:
        Y0 /= max_sing
    return Y0

def build_underlying_eff_from_gen_shifts(delta_gen_u,
                                         delta_gen_d,
                                         delta_gen_e,
                                         delta_gen_nu,
                                         M0,
                                         N_eff=180,
                                         phases_u=(0,6),
                                         phases_d=(0,3),
                                         phases_e=(0,8),
                                         phases_nu=(0,24)):
    n0u, delu = phases_u
    n0d, deld = phases_d
    n0e, dele = phases_e
    n0n, deln = phases_nu

    Yu0  = generate_proto_Y_from_gen_shifts(EXP_UP,   n0u, delu, delta_gen_u,  N_eff)
    Yd0  = generate_proto_Y_from_gen_shifts(EXP_DOWN, n0d, deld, delta_gen_d,  N_eff)
    Ye0  = generate_proto_Y_from_gen_shifts(EXP_LEP,  n0e, dele, delta_gen_e,  N_eff)
    Ynu0 = generate_proto_Y_from_gen_shifts(EXP_NU,   n0n, deln, delta_gen_nu, N_eff)

    Yu9  = align_9x9(Yu0)
    Yd9  = align_9x9(Yd0)
    Ye9  = align_9x9(Ye0)
    Ynu9 = align_9x9(Ynu0)
    M9   = align_9x9(M0)

    Yu_eff  = schur_9_to_3(Yu9)
    Yd_eff  = schur_9_to_3(Yd9)
    Ye_eff  = schur_9_to_3(Ye9)
    Ynu_eff = schur_9_to_3(Ynu9)

    Mnu_eff = triadic_Majorana_seesaw(M9, Ynu_eff)
    return Yu_eff, Yd_eff, Ye_eff, Mnu_eff

def evaluate_phase_shift_config(delta_gen_u,
                                delta_gen_d,
                                delta_gen_e,
                                delta_gen_nu,
                                M0,
                                N_eff=180,
                                lambda_geom=0.0):
    Yu_eff, Yd_eff, Ye_eff, Mnu_eff = build_underlying_eff_from_gen_shifts(
        delta_gen_u, delta_gen_d, delta_gen_e, delta_gen_nu, M0, N_eff
    )

    QL_u   = (0, 0, 0)
    QL_d   = (1, 0, 0)
    QL_e   = (1, 0, 0)
    QL_nu  = (0, 0, 0)
    eps_L  = 0.3

    Yu_FN   = apply_left_FN_3x3(Yu_eff,  QL_u,  eps_L)
    Yd_FN   = apply_left_FN_3x3(Yd_eff,  QL_d,  eps_L)
    Ye_FN   = apply_left_FN_3x3(Ye_eff,  QL_e,  eps_L)
    Mnu_FN  = apply_left_FN_Majorana_3x3(Mnu_eff, QL_nu, eps_L)

    Yu_EW  = stub_rge_run(Yu_FN)
    Yd_EW  = stub_rge_run(Yd_FN)
    Ye_EW  = stub_rge_run(Ye_FN)
    Mnu_EW = stub_rge_run(Mnu_FN)

    obs, Vckm, Upmns = compute_observables_from_matrices(Yu_EW, Yd_EW, Ye_EW, Mnu_EW)
    chi2, pulls = chi2_from_obs(obs)

    geom_penalty = (np.sum(np.array(delta_gen_u)**2) +
                    np.sum(np.array(delta_gen_d)**2) +
                    np.sum(np.array(delta_gen_e)**2) +
                    np.sum(np.array(delta_gen_nu)**2))
    cost = chi2 + lambda_geom * geom_penalty

    return cost, chi2, geom_penalty, obs, pulls, Vckm, Upmns

# ==================================
# Best Δ-gen and test run
# ==================================

def run_best_with_geom(lambda_geom=0.05):
    # Best Δ-gen from χ²-only optimization
    delta_gen_u  = np.array([-0.02188686,  0.34068281, -2.26590014])
    delta_gen_d  = np.array([-4.03359769, -0.28489373,  0.50558494])
    delta_gen_e  = np.array([ 0.38241537,  1.68425423, -0.74604266])
    delta_gen_nu = np.array([ 2.33108129, -0.20603719, -0.02546716])

    rng_M = np.random.default_rng(9)
    M0 = generate_proto_Majorana(rng_M)

    cost, chi2, geom_pen, obs, pulls, Vckm, Upmns = evaluate_phase_shift_config(
        delta_gen_u, delta_gen_d, delta_gen_e, delta_gen_nu, M0,
        N_eff=180, lambda_geom=lambda_geom
    )
    return cost, chi2, geom_pen, obs, pulls, Vckm, Upmns

if __name__ == "__main__":
    lambda_geom = 0.05
    cost, chi2_val, geom_pen, obs_val, pulls_val, Vckm_val, Upmns_val = run_best_with_geom(lambda_geom)
    print(f"Best (fixed Δ-gen) with λ_geom={lambda_geom}:")
    print(f"  χ²        ≈ {chi2_val:.3f}")
    print(f"  geom_pen  ≈ {geom_pen:.3f}")
    print(f"  cost      ≈ {cost:.3f}\n")

    print("Observables (model vs target, pull in σ):")
    for key in exp_targets.keys():
        print(f"  {key:12s}  model = {obs_val[key]:.6g}   "
              f"target = {exp_targets[key]:.6g}   pull = {pulls_val[key]: .3f}")

    print("\n|V_CKM| ≈")
    print(np.abs(Vckm_val))

    print("\n|U_PMNS| ≈")
    print(np.abs(Upmns_val))
