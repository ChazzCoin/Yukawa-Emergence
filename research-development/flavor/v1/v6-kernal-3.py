import numpy as np
import math
import cma

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

# Fixed triadic exponents per sector (BASELINES)
# These are the *alignment priors*; CMA-ES will optimize small shifts on top of them.
EXP_UP_BASE   = np.array((4.0, 2.0, 0.0), dtype=float)
EXP_DOWN_BASE = np.array((3.0, 2.0, 0.0), dtype=float)
EXP_LEP_BASE  = np.array((3.0, 2.0, 0.0), dtype=float)
EXP_NU_BASE   = np.array((1.0, 0.0, 0.0), dtype=float)

# Keep these names for backward compatibility where needed (if you still use them)
EXP_UP   = tuple(EXP_UP_BASE)
EXP_DOWN = tuple(EXP_DOWN_BASE)
EXP_LEP  = tuple(EXP_LEP_BASE)
EXP_NU   = tuple(EXP_NU_BASE)


# Forbidden distance on the 9-site ring (harmonic frustration point)
FORBIDDEN_D = 2

def build_site_scales(exponents):
    a, b, c = exponents
    gen_scales = np.array([EPS_ALIGN**a, EPS_ALIGN**b, EPS_ALIGN**c], dtype=float)
    s = np.zeros(N_SITES, dtype=float)
    for i in range(N_SITES):
        s[i] = gen_scales[generation_index(i)]
    return s

# ==================================
# GEOMETRIC KERNELS (SECTOR-DEPENDENT)
# ==================================

def build_kernel_sector(gamma: float, d_forbid: int = 2) -> np.ndarray:
    """
    Sector-dependent geometric kernel:
      - Toeplitz in distance d
      - forbidden distance d_forbid (hard zero)
      - exponential falloff exp(-gamma * d) * EPS_ALIGN**d
    """
    K = np.zeros((N_SITES, N_SITES), dtype=float)
    for i in range(N_SITES):
        for j in range(N_SITES):
            if i == j:
                K[i, j] = 1.0
            else:
                d = min(abs(i - j), N_SITES - abs(i - j))
                if d == d_forbid:
                    K[i, j] = 0.0
                else:
                    K[i, j] = math.exp(-gamma * d) * (EPS_ALIGN ** d)
    return K

def build_kernel(d_forbid: int = 2) -> np.ndarray:
    """
    Legacy uniform kernel (no sector dependence, just forbidden distance).
    Kept for debugging/reference.
    """
    K = np.zeros((N_SITES, N_SITES), dtype=float)
    for i in range(N_SITES):
        for j in range(N_SITES):
            if i == j:
                K[i, j] = 1.0
            else:
                d = min(abs(i - j), N_SITES - abs(i - j))
                if d == d_forbid:
                    K[i, j] = 0.0
                else:
                    K[i, j] = EPS_ALIGN ** d
    return K

# Sector coherence parameters (γ_x)
# Smaller gamma = more coherent; larger gamma = more decoherent
GAMMA_U   = 0.00  # up quarks: most coherent
GAMMA_D   = 0.03  # down quarks: slightly less coherent
GAMMA_E   = 0.05  # charged leptons
GAMMA_NU  = 0.08  # neutrinos: least coherent
GAMMA_MAJ = GAMMA_NU  # heavy Majorana sector aligned with neutrino window

# Precompute sector kernels
K_UP   = build_kernel_sector(GAMMA_U,  d_forbid=FORBIDDEN_D)
K_DOWN = build_kernel_sector(GAMMA_D,  d_forbid=FORBIDDEN_D)
K_E    = build_kernel_sector(GAMMA_E,  d_forbid=FORBIDDEN_D)
K_NU   = build_kernel_sector(GAMMA_NU, d_forbid=FORBIDDEN_D)
K_MAJ  = build_kernel_sector(GAMMA_MAJ, d_forbid=FORBIDDEN_D)

# A "default" kernel if needed (not used in core alignment now)
KERNEL = build_kernel(d_forbid=FORBIDDEN_D)

# ==================================
# HARMONIC PHASES (D_360, N_eff)
# ==================================
def align_9x9(M: np.ndarray) -> np.ndarray:
    """
    Apply the alignment kernel elementwise to a 9x9 matrix.
    In your setup, KERNEL is the 9x9 geometric alignment matrix.
    """
    return KERNEL * M

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
    """
    Implement a 9 → 6 → 3 structure:

      - Start from the full 9×9 aligned Majorana proto-matrix.
      - Restrict to the 6 heavy proto-sites (HEAVY_SITES): 9 → 6.
      - Project these 6 heavy modes onto 3 triadic heavy modes (B_H): 6 → 3.
      - Use the resulting 3×3 heavy Majorana matrix in a type-I seesaw with the
        3×3 effective Dirac Yukawa Ynu_eff.
    """

    # 9 → 6: take the heavy 6×6 block using the explicit heavy-site list
    M_H = M9_aligned[np.ix_(HEAVY_SITES, HEAVY_SITES)]  # shape (6, 6)

    # Label heavy states internally as indices 0..5 in this 6D heavy subspace
    h_indices = np.arange(len(HEAVY_SITES))  # = 0..5

    # 6 → 3 triadic heavy modes:
    # B_H is a 6×3 matrix whose columns are orthonormal triadic combinations
    # (discrete Fourier modes k = 1, 2, 3 on the 6 heavy sites).
    B_H = np.zeros((len(HEAVY_SITES), 3), dtype=complex)
    for col, k in enumerate([1, 2, 3]):
        B_H[:, col] = np.exp(2j * math.pi * k * h_indices / len(HEAVY_SITES)) / math.sqrt(len(HEAVY_SITES))

    # Project 6×6 → 3×3 heavy Majorana
    M_R_dimless = B_H.conj().T @ M_H @ B_H
    M_R_dimless = 0.5 * (M_R_dimless + M_R_dimless.T)  # ensure symmetric
    M_R_dimless += 1e-9 * np.eye(3)
    M_R = Lambda_Maj * M_R_dimless  # give it a physical scale

    # Dirac mass matrix (3×3) from effective Yukawa
    m_D = (v / math.sqrt(2.0)) * Ynu_eff

    # Seesaw: 3×3 light neutrino mass matrix
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
                                         phases_nu=(0,24),
                                         delta_exp_u=None,
                                         delta_exp_d=None,
                                         delta_exp_e=None,
                                         delta_exp_nu=None):
    """
    Build effective 3x3 Yukawas and light-neutrino mass matrix from:
      - generation phase shifts (delta_gen_*)
      - exponent shifts (delta_exp_*)
    Exponents are: EXP_*_BASE + delta_exp_* (with delta defaulting to 0).
    """

    # Unpack phase seeds
    n0u, delu = phases_u
    n0d, deld = phases_d
    n0e, dele = phases_e
    n0n, deln = phases_nu

    # Default exponent shifts = 0 (no deformation)
    if delta_exp_u is None:
        delta_exp_u = np.zeros(3, dtype=float)
    if delta_exp_d is None:
        delta_exp_d = np.zeros(3, dtype=float)
    if delta_exp_e is None:
        delta_exp_e = np.zeros(3, dtype=float)
    if delta_exp_nu is None:
        delta_exp_nu = np.zeros(3, dtype=float)

    # Effective exponents = base + shift
    exp_up_eff   = EXP_UP_BASE   + np.array(delta_exp_u, dtype=float)
    exp_down_eff = EXP_DOWN_BASE + np.array(delta_exp_d, dtype=float)
    exp_lep_eff  = EXP_LEP_BASE  + np.array(delta_exp_e, dtype=float)
    exp_nu_eff   = EXP_NU_BASE   + np.array(delta_exp_nu, dtype=float)

    # Sector-specific N_eff (as before)
    N_eff_u  = 360
    N_eff_d  = 180
    N_eff_e  = 120
    N_eff_nu = 90  # or 60

    # Build proto Yukawas with the *effective* exponent triples
    Yu0  = generate_proto_Y_from_gen_shifts(exp_up_eff,   n0u, delu, delta_gen_u,  N_eff_u)
    Yd0  = generate_proto_Y_from_gen_shifts(exp_down_eff, n0d, deld, delta_gen_d,  N_eff_d)
    Ye0  = generate_proto_Y_from_gen_shifts(exp_lep_eff,  n0e, dele, delta_gen_e,  N_eff_e)
    Ynu0 = generate_proto_Y_from_gen_shifts(exp_nu_eff,   n0n, deln, delta_gen_nu, N_eff_nu)

    # Geometric alignment to 9x9 and Schur reduction
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
                                lambda_geom=0.0,
                                delta_exp_u=None,
                                delta_exp_d=None,
                                delta_exp_e=None,
                                delta_exp_nu=None,
                                lambda_exp=0.1):
    """
    Evaluate cost given:
      - Δ-gen shifts (delta_gen_*)
      - exponent shifts (delta_exp_*)
    Cost = chi2 + lambda_geom * ||Δ-gen||^2 + lambda_exp * ||δ-exp||^2
    """

    Yu_eff, Yd_eff, Ye_eff, Mnu_eff = build_underlying_eff_from_gen_shifts(
        delta_gen_u, delta_gen_d, delta_gen_e, delta_gen_nu, M0, N_eff,
        delta_exp_u=delta_exp_u,
        delta_exp_d=delta_exp_d,
        delta_exp_e=delta_exp_e,
        delta_exp_nu=delta_exp_nu
    )

    # FN textures (unchanged)
    QL_u   = (0, 0, 0)
    QL_d   = (1, 0.5, 0)
    QL_e   = (0.5, 0, 0)
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

    # Δ-gen penalty
    geom_penalty = (np.sum(np.array(delta_gen_u)**2) +
                    np.sum(np.array(delta_gen_d)**2) +
                    np.sum(np.array(delta_gen_e)**2) +
                    np.sum(np.array(delta_gen_nu)**2))

    # exponent-shift penalty
    if delta_exp_u is None:
        delta_exp_u = np.zeros(3, dtype=float)
    if delta_exp_d is None:
        delta_exp_d = np.zeros(3, dtype=float)
    if delta_exp_e is None:
        delta_exp_e = np.zeros(3, dtype=float)
    if delta_exp_nu is None:
        delta_exp_nu = np.zeros(3, dtype=float)

    exp_penalty = (np.sum(np.array(delta_exp_u)**2) +
                   np.sum(np.array(delta_exp_d)**2) +
                   np.sum(np.array(delta_exp_e)**2) +
                   np.sum(np.array(delta_exp_nu)**2))

    cost = chi2 + lambda_geom * geom_penalty + lambda_exp * exp_penalty

    return cost, chi2, geom_penalty, exp_penalty, obs, pulls, Vckm, Upmns

# ==================================
# Best Δ-gen and test run
# ==================================

def run_best_with_geom(lambda_geom=0.00):
    # Best Δ-gen from some earlier optimization (can be updated later)
    delta_gen_u  = np.array([-0.02188686,  0.34068281, -2.26590014])
    delta_gen_d  = np.array([-4.03359769, -0.28489373,  0.50558494])
    delta_gen_e  = np.array([ 0.38241537,  1.68425423, -0.74604266])
    delta_gen_nu = np.array([ 2.33108129, -0.20603719, -0.02546716])

    rng_M = np.random.default_rng(9)
    M0 = generate_proto_Majorana(rng_M)

    cost, chi2, geom_pen, exp_pen, obs, pulls, Vckm, Upmns = evaluate_phase_shift_config(
        delta_gen_u, delta_gen_d, delta_gen_e, delta_gen_nu, M0,
        N_eff=180, lambda_geom=lambda_geom
    )
    return cost, chi2, geom_pen, exp_pen, obs, pulls, Vckm, Upmns



# ===============================
# Δ-gen vector utilities
# ===============================

# ===============================
# Parameter vector utilities
# ===============================

def pack_params(du, dd, de, dn,
                deu, ded, dee, den):
    """
    Flatten (4×3 Δ-gen + 4×3 δ-exp) → 24-vector.
    Order:
      [Δu, Δd, Δe, Δnu, δexp_u, δexp_d, δexp_e, δexp_nu]
    """
    return np.concatenate([du, dd, de, dn, deu, ded, dee, den])


def unpack_params(X):
    """
    24-vector → (Δu, Δd, Δe, Δnu, δexp_u, δexp_d, δexp_e, δexp_nu)
    """
    X = np.array(X, dtype=float)
    du   = X[0:3]
    dd   = X[3:6]
    de   = X[6:9]
    dn   = X[9:12]
    deu  = X[12:15]
    ded  = X[15:18]
    dee  = X[18:21]
    den  = X[21:24]
    return du, dd, de, dn, deu, ded, dee, den

# ===============================
# Cost wrapper for optimizers
# ===============================

def params_cost_vectorized(X,
                           M0,
                           lambda_geom=0.05,
                           lambda_exp=0.1,
                           N_eff=180):
    """
    X = 24-dimensional parameter vector:
        [Δ-gen (12), δ-exp (12)].
    Returns scalar cost used by CMA-ES.
    """
    du, dd, de, dn, deu, ded, dee, den = unpack_params(X)

    try:
        cost, chi2, geom_pen, exp_pen, obs, pulls, Vckm, Upmns = evaluate_phase_shift_config(
            du, dd, de, dn, M0,
            N_eff=N_eff,
            lambda_geom=lambda_geom,
            delta_exp_u=deu,
            delta_exp_d=ded,
            delta_exp_e=dee,
            delta_exp_nu=den,
            lambda_exp=lambda_exp
        )
    except Exception:
        return 1e9

    if not math.isfinite(cost):
        return 1e9
    return cost


# ===============================
# Run CMA-ES with random restarts
# ===============================

def optimize_params_CMA(num_restarts=6,
                        sigma_init=0.3,
                        lambda_geom=0.05,
                        lambda_exp=0.1,
                        N_eff=180,
                        seed=42):

    rng = np.random.default_rng(seed)
    M0 = generate_proto_Majorana(rng)  # fix proto-Majorana for the scan

    best_cost = np.inf
    best_X = None

    for r in range(num_restarts):
        print(f"\n====== CMA-ES Restart {r+1}/{num_restarts} ======")

        # Initial guess: Δ-gen ~ N(0,0.4), δ-exp ~ N(0,0.2)
        X0 = np.zeros(24, dtype=float)
        X0[0:12]  = rng.normal(scale=0.4, size=12)   # Δ-gen
        X0[12:24] = rng.normal(scale=0.2, size=12)   # δ-exp (smaller)

        es = cma.CMAEvolutionStrategy(X0, sigma_init,
                                      {'popsize': 20, 'maxiter': 200})

        while not es.stop():
            solutions = es.ask()
            costs = [params_cost_vectorized(x, M0, lambda_geom, lambda_exp, N_eff)
                     for x in solutions]
            es.tell(solutions, costs)
            es.disp()

        Xbest = es.best.x
        costbest = es.best.f

        print(f"Restart {r+1}: best cost = {costbest:.4f}")

        if costbest < best_cost:
            best_cost = costbest
            best_X = Xbest.copy()

    print("\n======= GLOBAL BEST FOUND =======")
    print(f"Cost = {best_cost:.5f}")
    print(f"Xbest = {best_X}")

    du, dd, de, dn, deu, ded, dee, den = unpack_params(best_X)
    return du, dd, de, dn, deu, ded, dee, den, best_cost, M0


# ==================================
# Main: single evaluation with fixed Δ-gen
# ==================================

fits = []

def main():
    lambda_geom = 0.00
    cost, chi2_val, geom_pen, exp_pen, obs_val, pulls_val, Vckm_val, Upmns_val = run_best_with_geom(lambda_geom)

    fit = f"{chi2_val:.3f}"
    fits.append(fit)
    print(f"Best (fixed Δ-gen) with λ_geom={lambda_geom}:")
    print(f"  χ²        ≈ {chi2_val:.3f}")
    print(f"  geom_pen  ≈ {geom_pen:.3f}")
    print(f"  exp_pen   ≈ {exp_pen:.3f}")
    print(f"  cost      ≈ {cost:.3f}\n")

    print("Observables (model vs target, pull in σ):")
    for key in exp_targets.keys():
        print(f"  {key:12s}  model = {obs_val[key]:.6g}   "
              f"target = {exp_targets[key]:.6g}   pull = {pulls_val[key]: .3f}")

    print("\n|V_CKM| ≈")
    print(np.abs(Vckm_val))

    print("\n|U_PMNS| ≈")
    print(np.abs(Upmns_val))


# if __name__ == "__main__":
#     main()
if __name__ == "__main__":
    # First, optional: baseline check
    FORBIDDEN_D = 2      # if you want the d=2 geometry you liked
    KERNEL = build_kernel()

    lambda_geom = 0.05
    lambda_exp  = 0.1

    # Run joint optimizer over Δ-gen + δ-exponents
    du, dd, de, dn, deu, ded, dee, den, best_cost, M0 = optimize_params_CMA(
        num_restarts=8,
        sigma_init=0.3,
        lambda_geom=lambda_geom,
        lambda_exp=lambda_exp,
        N_eff=180,
        seed=9
    )

    print("\nOptimized parameters:")
    print("delta_gen_u  =", du)
    print("delta_gen_d  =", dd)
    print("delta_gen_e  =", de)
    print("delta_gen_nu =", dn)
    print("delta_exp_u  =", deu)
    print("delta_exp_d  =", ded)
    print("delta_exp_e  =", dee)
    print("delta_exp_nu =", den)
    print("best_cost =", best_cost)

    # Evaluate observables at this optimum
    cost, chi2, geom_pen, exp_pen, obs, pulls, Vckm, Upmns = evaluate_phase_shift_config(
        du, dd, de, dn, M0,
        N_eff=180,
        lambda_geom=lambda_geom,
        delta_exp_u=deu,
        delta_exp_d=ded,
        delta_exp_e=dee,
        delta_exp_nu=den,
        lambda_exp=lambda_exp
    )

    print(f"\nAt optimized params: chi2 ≈ {chi2:.3f}, "
          f"geom_pen ≈ {geom_pen:.3f}, exp_pen ≈ {exp_pen:.3f}, cost ≈ {cost:.3f}\n")

    print("Observables (model vs target, pull in σ):")
    for key in exp_targets.keys():
        print(f"  {key:12s}  model = {obs[key]:.6g}   "
              f"target = {exp_targets[key]:.6g}   pull = {pulls[key]: .3f}")

    print("\n|V_CKM| ≈")
    print(np.abs(Vckm))

    print("\n|U_PMNS| ≈")
    print(np.abs(Upmns))

