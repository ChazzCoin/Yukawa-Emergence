import numpy as np
import math
import cma

# ==================================
# GEOMETRIC AXIOMS (FIXED)
# ==================================

# KAPPA = 360.0 / 89.0
KAPPA = 0.242
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

# ============================================================
#  Projection-tweak utilities (true 9→3 geometric projection)
# ============================================================
# Global leakage parameters (optimized by CMA-ES)
# These MUST exist before build_default_projection is called.
build_default_projection_eps12 = 0.0
build_default_projection_eps21 = 0.0
#
# ============================================================
# FULL TRIADIC FOURIER PROJECTION BASIS (RECOMMENDED)
# ============================================================

# Global leakage parameters (for CMA optimization)
proj_eps_03 = 0.0   # leakage from k=0 → k=3 mode
proj_eps_30 = 0.0   # leakage from k=3 → k=0
proj_eps_36 = 0.0   # mixing between k=3 and k=6

def triadic_fourier_modes():
    """
    Returns (v0, v3, v6) — the triadic DFT basis vectors on Z9.
    """
    n = 9
    j = np.arange(n)
    v0 = np.exp(2j*np.pi*0*j/n) / np.sqrt(n)
    v3 = np.exp(2j*np.pi*3*j/n) / np.sqrt(n)
    v6 = np.exp(2j*np.pi*6*j/n) / np.sqrt(n)
    return v0, v3, v6


def build_default_projection(proj_eps_03: float = 0.0,
                             proj_eps_30: float = 0.0,
                             proj_eps_36: float = 0.0) -> np.ndarray:
    """
    3×9 triadic Fourier projection with small mode mixing.
    proj_eps_03, proj_eps_30, proj_eps_36 are *local* leakage parameters
    (no globals), so CMA-ES can safely steer them.

    We mix using the original Fourier modes (v0,v3,v6) to avoid runaway
    feedback before QR.
    """
    v0, v3, v6 = triadic_fourier_modes()

    # Start from original triadic modes
    b0 = v0.copy()
    b3 = v3.copy()
    b6 = v6.copy()

    # Small mode mixing using the *original* modes
    b0 = b0 + proj_eps_03 * v3
    b3 = b3 + proj_eps_30 * v0 + proj_eps_36 * v6
    b6 = b6 + proj_eps_36 * v3  # symmetric-ish coupling

    # Orthonormalize rows via QR on the transpose
    B = np.vstack([b0, b3, b6])          # shape (3, 9)
    Q_cols, _ = np.linalg.qr(B.conj().T) # (9, 3) with orthonormal columns
    P = Q_cols.T                         # (3, 9) with orthonormal rows

    return P




def tweak_projection(P_base, g, s, a, phi):
    """
    Add a small complex deformation to row g, site s, then re-orthonormalize.
    """
    P = P_base.copy()
    P[g, s] += a * np.exp(1j * phi)

    Q = np.zeros_like(P, dtype=complex)

    Q[0] = P[0] / np.linalg.norm(P[0])

    v1 = P[1] - np.vdot(Q[0], P[1]) * Q[0]
    Q[1] = v1 / np.linalg.norm(v1)

    v2 = P[2] \
         - np.vdot(Q[0], P[2]) * Q[0] \
         - np.vdot(Q[1], P[2]) * Q[1]
    Q[2] = v2 / np.linalg.norm(v2)

    return Q


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

# Precompute sector kernels (currently not used directly in align_9x9, but kept for later refinement)
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

def align_9x9_sector(M: np.ndarray, sector: str) -> np.ndarray:
    """
    Sector-dependent geometric alignment:
      sector = 'u', 'd', 'e', 'nu', or 'maj'

    Uses the precomputed 9×9 kernels:
      K_UP, K_DOWN, K_E, K_NU, K_MAJ

    Falls back to the legacy KERNEL if an unknown sector is passed.
    """
    if sector == "u":
        K = K_UP
    elif sector == "d":
        K = K_DOWN
    elif sector == "e":
        K = K_E
    elif sector == "nu":
        K = K_NU
    elif sector == "maj":
        K = K_MAJ
    else:
        # fallback: legacy universal kernel
        K = KERNEL

    return K * M


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
# Small left-handed projection tweak (3×3) for down sector
# ==================================

def build_left_projection_tweak(a_tq: float, phi_tq: float) -> np.ndarray:
    """
    Build a small unitary acting on generations 1 and 2 in the left-handed space.
    This effectively mimics a tiny misalignment in the 9→3 projection basis:
        U_tq ~ exp(i * a_tq * generator_12)
    We parametrize it as a 1–2 rotation with a complex phase.
    """
    theta = a_tq  # we treat 'a_tq' directly as a small angle
    c = math.cos(theta)
    s = math.sin(theta)

    # e^{i phi}, e^{-i phi}
    eip = complex(math.cos(phi_tq), math.sin(phi_tq))
    eim = complex(math.cos(-phi_tq), math.sin(-phi_tq))

    U = np.eye(3, dtype=complex)
    # 1-2 block
    U[0, 0] = c
    U[0, 1] = s * eip
    U[1, 0] = -s * eim
    U[1, 1] = c
    # 3rd generation untouched
    return U

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

def build_underlying_eff_from_gen_shifts(
        delta_gen_u,
        delta_gen_d,
        delta_gen_e,
        delta_gen_nu,
        M0,
        N_eff=180,
        phases_u=(0, 6),
        phases_d=(0, 3),
        phases_e=(0, 8),
        phases_nu=(0, 24),
        delta_exp_u=None,
        delta_exp_d=None,
        delta_exp_e=None,
        delta_exp_nu=None,
        P=None,
):
    """
    Build effective 3×3 Yukawas and light-neutrino mass matrix from:
      - generation phase shifts (delta_gen_*)
      - exponent shifts (delta_exp_*)
      - optional 9→3 projection matrix P (used only in neutrino sector)

    Geometry:
      • up, down, charged leptons: Schur 9→3 (as in your best-fit runs)
      • neutrinos: 9→3 projection with P, then triadic seesaw

    Exponents:
        exp_eff = EXP_*_BASE + delta_exp_*
    (delta_exp_* default to 0 if None)
    """

    # -----------------------------
    # Phase seeds per sector
    # -----------------------------
    n0u, delu = phases_u
    n0d, deld = phases_d
    n0e, dele = phases_e
    n0n, deln = phases_nu

    # -----------------------------
    # Default exponent shifts
    # -----------------------------
    if delta_exp_u is None:
        delta_exp_u = np.zeros(3, dtype=float)
    if delta_exp_d is None:
        delta_exp_d = np.zeros(3, dtype=float)
    if delta_exp_e is None:
        delta_exp_e = np.zeros(3, dtype=float)
    if delta_exp_nu is None:
        delta_exp_nu = np.zeros(3, dtype=float)

    # -----------------------------
    # Effective exponents (base + shift)
    # -----------------------------
    exp_up_eff   = EXP_UP_BASE   + np.array(delta_exp_u, dtype=float)
    exp_down_eff = EXP_DOWN_BASE + np.array(delta_exp_d, dtype=float)
    exp_lep_eff  = EXP_LEP_BASE  + np.array(delta_exp_e, dtype=float)
    exp_nu_eff   = EXP_NU_BASE   + np.array(delta_exp_nu, dtype=float)

    # -----------------------------
    # Sector-specific N_eff
    # -----------------------------
    N_eff_u  = 360
    N_eff_d  = 180
    N_eff_e  = 120
    N_eff_nu = 90  # or 60, as you like

    # -----------------------------
    # Build proto Yukawas (9×9) with effective exponents
    # -----------------------------
    Yu0  = generate_proto_Y_from_gen_shifts(exp_up_eff,   n0u, delu, delta_gen_u,  N_eff_u)
    Yd0  = generate_proto_Y_from_gen_shifts(exp_down_eff, n0d, deld, delta_gen_d,  N_eff_d)
    Ye0  = generate_proto_Y_from_gen_shifts(exp_lep_eff,  n0e, dele, delta_gen_e,  N_eff_e)
    Ynu0 = generate_proto_Y_from_gen_shifts(exp_nu_eff,   n0n, deln, delta_gen_nu, N_eff_nu)

    # -----------------------------
    # Geometric alignment 9×9 (sector-dependent)
    # -----------------------------
    Yu9  = align_9x9_sector(Yu0,  "u")
    Yd9  = align_9x9_sector(Yd0,  "d")
    Ye9  = align_9x9_sector(Ye0,  "e")
    Ynu9 = align_9x9_sector(Ynu0, "nu")
    M9   = align_9x9_sector(M0,   "maj")


    # -----------------------------
    # 9→3 reduction:
    #   • quarks & charged leptons: Schur 9→3 (LIGHT_SITES vs HEAVY_SITES)
    #   • neutrinos: 9→3 projection with P (triadic/Fourier)
    # -----------------------------
    # Quark and charged-lepton sectors: keep the successful Schur geometry
    Yu_eff = schur_9_to_3(Yu9)   # 3×3
    Yd_eff = schur_9_to_3(Yd9)   # 3×3
    Ye_eff = schur_9_to_3(Ye9)   # 3×3

    # Neutrino sector: triadic Fourier 9→3 projection
    if P is None:
        P = build_default_projection()  # 3×9

    Ynu_eff = P @ Ynu9 @ P.conj().T     # 3×3

    # -----------------------------
    # Triadic Majorana seesaw with 3×3 Dirac neutrino Yukawa
    # -----------------------------
    Mnu_eff = triadic_Majorana_seesaw(M9, Ynu_eff)

    return Yu_eff, Yd_eff, Ye_eff, Mnu_eff




def evaluate_phase_shift_config(
        delta_gen_u,
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
        lambda_exp=0.1,
        a_tq=0.0,
        phi_tq=0.0,
        eps12=0.0,
        eps21=0.0,
        tweak_row=1,
        tweak_site=1,
        lambda_proj=1.0,
        proj_eps_03=0.0,
        proj_eps_30=0.0,
        proj_eps_36=0.0,
):
    """
    Evaluate cost given:
      - Δ-gen shifts
      - δ-exp shifts
      - projection tweaks in the neutrino sector (triadic Fourier leakage)
      - a small left-handed 1–2 twist in the *down* sector (a_tq, phi_tq)

    Geometry summary:
      • Yu, Yd, Ye: 9×9 → (sector kernels) → Schur 9→3
      • Yν:         9×9 → (kernel) → 9→3 triadic-Fourier projection P
      • M_R:        9×9 → heavy 6×6 → triadic 6→3 → seesaw
    """

    # ============================================================
    # 9→3 projection matrix for neutrinos:
    # triadic Fourier + local mode leakage, but *no* a_tq here.
    # ============================================================
    P_tq = build_default_projection(
        proj_eps_03=proj_eps_03,
        proj_eps_30=proj_eps_30,
        proj_eps_36=proj_eps_36,
    )

    # ============================================================
    # ALIGNMENT + EFFECTIVE YUKAWA SECTORS
    # ============================================================
    Yu_eff, Yd_eff, Ye_eff, Mnu_eff = build_underlying_eff_from_gen_shifts(
        delta_gen_u, delta_gen_d, delta_gen_e, delta_gen_nu,
        M0, N_eff,
        delta_exp_u=delta_exp_u,
        delta_exp_d=delta_exp_d,
        delta_exp_e=delta_exp_e,
        delta_exp_nu=delta_exp_nu,
        P=P_tq
    )

    # ============================================================
    # OPTIONAL left-handed 1–2 generation twist for the down sector
    # (this is your Cabibbo handle in the quark sector)
    # ============================================================
    if abs(a_tq) > 0.0:
        U_tq = build_left_projection_tweak(a_tq, phi_tq)
        Yd_eff = U_tq @ Yd_eff

    # ============================================================
    # FN dressing
    # ============================================================
    QL_u   = (0, 0, 0)
    QL_d   = (1, 0.5, 0)
    QL_e   = (0.5, 0, 0)
    QL_nu  = (0, 0, 0)
    eps_L  = 0.3

    Yu_FN  = apply_left_FN_3x3(Yu_eff,  QL_u, eps_L)
    Yd_FN  = apply_left_FN_3x3(Yd_eff,  QL_d, eps_L)
    Ye_FN  = apply_left_FN_3x3(Ye_eff,  QL_e, eps_L)
    Mnu_FN = apply_left_FN_Majorana_3x3(Mnu_eff, QL_nu, eps_L)

    # ============================================================
    # RG evolution
    # ============================================================
    Yu_EW  = stub_rge_run(Yu_FN)
    Yd_EW  = stub_rge_run(Yd_FN)
    Ye_EW  = stub_rge_run(Ye_FN)
    Mnu_EW = stub_rge_run(Mnu_FN)

    # ============================================================
    # Observables + χ²
    # ============================================================
    obs, Vckm, Upmns = compute_observables_from_matrices(
        Yu_EW, Yd_EW, Ye_EW, Mnu_EW
    )
    chi2, pulls = chi2_from_obs(obs)

    # ============================================================
    # Penalties
    # ============================================================
    geom_penalty = (
        np.sum(np.array(delta_gen_u)**2) +
        np.sum(np.array(delta_gen_d)**2) +
        np.sum(np.array(delta_gen_e)**2) +
        np.sum(np.array(delta_gen_nu)**2)
    )

    if delta_exp_u is None: delta_exp_u = np.zeros(3)
    if delta_exp_d is None: delta_exp_d = np.zeros(3)
    if delta_exp_e is None: delta_exp_e = np.zeros(3)
    if delta_exp_nu is None: delta_exp_nu = np.zeros(3)

    exp_penalty = (
        np.sum(delta_exp_u**2) +
        np.sum(delta_exp_d**2) +
        np.sum(delta_exp_e**2) +
        np.sum(delta_exp_nu**2)
    )

    # Projection / leakage penalties
    # (eps12/eps21 are kept as small regularized knobs; they currently do not
    #  enter the projection explicitly, only the cost.)
    proj_penalty = (
        a_tq**2 +
        eps12**2 + eps21**2 +
        proj_eps_03**2 + proj_eps_30**2 + proj_eps_36**2
    )

    # ============================================================
    # Total cost
    # ============================================================
    cost = (
        chi2
        + lambda_geom * geom_penalty
        + lambda_exp  * exp_penalty
        + lambda_proj * proj_penalty
    )

    return (
        cost, chi2,
        geom_penalty, exp_penalty, proj_penalty,
        obs, pulls, Vckm, Upmns
    )


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

    cost, chi2, geom_pen, exp_pen, proj_pen, obs, pulls, Vckm, Upmns = evaluate_phase_shift_config(
        delta_gen_u, delta_gen_d, delta_gen_e, delta_gen_nu, M0,
        N_eff=180,
        lambda_geom=lambda_geom,
        lambda_exp=0.0,   # baseline: no exp penalty
        a_tq=0.0,         # no projection tweak
        phi_tq=0.0,
        lambda_proj=0.0   # no projection regularization
    )
    return cost, chi2, geom_pen, exp_pen, proj_pen, obs, pulls, Vckm, Upmns



# ===============================
# Parameter vector utilities
# ===============================

def pack_params(du, dd, de, dn,
                deu, ded, dee, den,
                a_tq, phi_tq,
                eps12, eps21):
    return np.concatenate([
        du, dd, de, dn,
        deu, ded, dee, den,
        np.array([a_tq, phi_tq, eps12, eps21])
    ])




def unpack_params(X):
    """
    31-dimensional parameter vector X ->
      (Δu, Δd, Δe, Δnu,
       δexp_u, δexp_d, δexp_e, δexp_nu,
       a_tq, phi_tq,
       eps12, eps21,
       proj_eps_03, proj_eps_30, proj_eps_36)
    """
    X = np.array(X, dtype=float)

    du  = X[0:3]
    dd  = X[3:6]
    de  = X[6:9]
    dn  = X[9:12]

    deu = X[12:15]
    ded = X[15:18]
    dee = X[18:21]
    den = X[21:24]

    a_tq   = X[24]
    phi_tq = X[25]

    eps12 = X[26]
    eps21 = X[27]

    proj_eps_03 = X[28]
    proj_eps_30 = X[29]
    proj_eps_36 = X[30]

    return (
        du, dd, de, dn,
        deu, ded, dee, den,
        a_tq, phi_tq,
        eps12, eps21,
        proj_eps_03, proj_eps_30, proj_eps_36
    )




# ===============================
# Cost wrapper for optimizers
# ===============================

def params_cost_vectorized(X,
                           M0,
                           lambda_geom=0.05,
                           lambda_exp=0.1,
                           lambda_proj=1.0,
                           N_eff=180):
    """
    X = 31-dimensional parameter vector:
        [Δ-gen (12), δ-exp (12),
         a_tq, phi_tq,
         eps12, eps21,
         proj_eps_03, proj_eps_30, proj_eps_36]

    Returns scalar cost used by CMA-ES.
    """
    (
        du, dd, de, dn,
        deu, ded, dee, den,
        a_tq, phi_tq,
        eps12, eps21,
        proj_eps_03, proj_eps_30, proj_eps_36
    ) = unpack_params(X)

    try:
        cost, chi2, geom_pen, exp_pen, proj_pen, obs, pulls, Vckm, Upmns = \
            evaluate_phase_shift_config(
                delta_gen_u=du,
                delta_gen_d=dd,
                delta_gen_e=de,
                delta_gen_nu=dn,
                M0=M0,
                N_eff=N_eff,
                lambda_geom=lambda_geom,
                delta_exp_u=deu,
                delta_exp_d=ded,
                delta_exp_e=dee,
                delta_exp_nu=den,
                lambda_exp=lambda_exp,
                a_tq=a_tq,
                phi_tq=phi_tq,
                eps12=eps12,
                eps21=eps21,
                proj_eps_03=proj_eps_03,
                proj_eps_30=proj_eps_30,
                proj_eps_36=proj_eps_36,
                lambda_proj=lambda_proj
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
                        lambda_proj=1.0,
                        N_eff=180,
                        seed=42):

    rng = np.random.default_rng(seed)
    M0 = generate_proto_Majorana(rng)  # fix proto-Majorana for the scan

    best_cost = np.inf
    best_X = None

    for r in range(num_restarts):
        print(f"\n====== CMA-ES Restart {r+1}/{num_restarts} ======")

        # ------------------------------------------------------------------
        # PARAMETER VECTOR LAYOUT (31)
        #   0–11: Δ-gen (12)
        #  12–23: δ-exp (12)
        #     24: a_tq
        #     25: phi_tq
        #     26: eps12       (projection leakage 1→2)
        #     27: eps21       (projection leakage 2→1)
        #     28: proj_eps_03 (Fourier leakage between sites 0 and 3)
        #     29: proj_eps_30
        #     30: proj_eps_36
        # ------------------------------------------------------------------
        X0 = np.zeros(31)
        X0[0:12]  = rng.normal(scale=0.4, size=12)   # Δ-gen
        X0[12:24] = rng.normal(scale=0.2, size=12)   # δ-exp

        # Projection / leakage knobs start at 0 (no deformation)
        X0[24] = 0.0   # a_tq
        X0[25] = 0.0   # phi_tq
        X0[26] = 0.0   # eps12
        X0[27] = 0.0   # eps21
        X0[28] = 0.0   # proj_eps_03
        X0[29] = 0.0   # proj_eps_30
        X0[30] = 0.0   # proj_eps_36

        es = cma.CMAEvolutionStrategy(
            X0,
            sigma_init,
            {
                'popsize': 20,
                'maxiter': 200,
                'CMA_diagonal': False
            }
        )

        while not es.stop():
            solutions = es.ask()
            costs = [
                params_cost_vectorized(
                    x, M0,
                    lambda_geom=lambda_geom,
                    lambda_exp=lambda_exp,
                    lambda_proj=lambda_proj,
                    N_eff=N_eff
                )
                for x in solutions
            ]
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
    print("Xbest =", best_X)

    (
        du, dd, de, dn,
        deu, ded, dee, den,
        a_tq, phi_tq,
        eps12, eps21,
        proj_eps_03, proj_eps_30, proj_eps_36
    ) = unpack_params(best_X)

    print("a_tq         =", a_tq)
    print("phi_tq       =", phi_tq)
    print("eps12        =", eps12)
    print("eps21        =", eps21)
    print("proj_eps_03  =", proj_eps_03)
    print("proj_eps_30  =", proj_eps_30)
    print("proj_eps_36  =", proj_eps_36)

    return (
        du, dd, de, dn,
        deu, ded, dee, den,
        a_tq, phi_tq,
        eps12, eps21,
        proj_eps_03, proj_eps_30, proj_eps_36,
        best_cost, M0
    )



# ==================================
# Main: optional single evaluation with fixed Δ-gen
# ==================================

fits = []

def main():
    lambda_geom = 0.00
    cost, chi2_val, geom_pen, exp_pen, proj_pen, obs_val, pulls_val, Vckm_val, Upmns_val = run_best_with_geom(lambda_geom)

    fit = f"{chi2_val:.3f}"
    fits.append(fit)
    print(f"Best (fixed Δ-gen) with λ_geom={lambda_geom}:")
    print(f"  χ²        ≈ {chi2_val:.3f}")
    print(f"  geom_pen  ≈ {geom_pen:.3f}")
    print(f"  exp_pen   ≈ {exp_pen:.3f}")
    print(f"  proj_pen  ≈ {proj_pen:.3f}")
    print(f"  cost      ≈ {cost:.3f}\n")

    print("Observables (model vs target, pull in σ):")
    for key in exp_targets.keys():
        print(f"  {key:12s}  model = {obs_val[key]:.6g}   "
              f"target = {exp_targets[key]:.6g}   pull = {pulls_val[key]: .3f}")

    print("\n|V_CKM| ≈")
    print(np.abs(Vckm_val))

    print("\n|U_PMNS| ≈")
    print(np.abs(Upmns_val))


if __name__ == "__main__":

    # ------------------------------------------------------------
    # Core geometric kernel for 9-site ring
    # ------------------------------------------------------------
    FORBIDDEN_D = 2
    KERNEL = build_kernel(FORBIDDEN_D)

    # Penalty strengths
    lambda_geom = 0.05
    lambda_exp  = 0.10
    lambda_proj = 1.00

    # ------------------------------------------------------------
    # Run joint optimizer over:
    #   Δ-gen (12)
    #   δ-exp  (12)
    #   projection tweaks: a_tq, phi_tq
    #   leakage terms: eps12, eps21, proj_eps_03, proj_eps_30, proj_eps_36
    # ------------------------------------------------------------
    results = optimize_params_CMA(
        num_restarts=8,
        sigma_init=0.3,
        lambda_geom=lambda_geom,
        lambda_exp=lambda_exp,
        lambda_proj=lambda_proj,
        N_eff=180,
        seed=9
    )

    # Unpack optimizer output
    (
        du, dd, de, dn,
        deu, ded, dee, den,
        a_tq, phi_tq,
        eps12, eps21,
        proj_eps_03, proj_eps_30, proj_eps_36,
        best_cost,
        M0
    ) = results

    # ------------------------------------------------------------
    # Print optimized parameters
    # ------------------------------------------------------------
    print("\nOptimized parameters:")
    print("δ_gen_u       =", du)
    print("δ_gen_d       =", dd)
    print("δ_gen_e       =", de)
    print("δ_gen_ν       =", dn)
    print("δ_exp_u       =", deu)
    print("δ_exp_d       =", ded)
    print("δ_exp_e       =", dee)
    print("δ_exp_ν       =", den)
    print("a_tq (proj)   =", a_tq)
    print("phi_tq        =", phi_tq)
    print("eps12         =", eps12)
    print("eps21         =", eps21)
    print("proj_eps_03   =", proj_eps_03)
    print("proj_eps_30   =", proj_eps_30)
    print("proj_eps_36   =", proj_eps_36)
    print("best_cost     =", best_cost)

    # ------------------------------------------------------------
    # Evaluate at optimal parameters
    # ------------------------------------------------------------
    cost, chi2, geom_pen, exp_pen, proj_pen, obs, pulls, Vckm, Upmns = \
        evaluate_phase_shift_config(
            du, dd, de, dn,
            M0,
            N_eff=180,
            lambda_geom=lambda_geom,
            delta_exp_u=deu,
            delta_exp_d=ded,
            delta_exp_e=dee,
            delta_exp_nu=den,
            lambda_exp=lambda_exp,
            a_tq=a_tq,
            phi_tq=phi_tq,
            eps12=eps12,
            eps21=eps21,
            proj_eps_03=proj_eps_03,
            proj_eps_30=proj_eps_30,
            proj_eps_36=proj_eps_36,
            lambda_proj=lambda_proj
        )

    # ------------------------------------------------------------
    # Print evaluation summary
    # ------------------------------------------------------------
    print(f"\nAt optimized params:")
    print(f"  χ²         = {chi2:.3f}")
    print(f"  geom_pen   = {geom_pen:.3f}")
    print(f"  exp_pen    = {exp_pen:.3f}")
    print(f"  proj_pen   = {proj_pen:.3f}")
    print(f"  total cost = {cost:.3f}\n")

    print("Observables (model vs target, pull in σ):")
    for key in exp_targets.keys():
        print(f"  {key:12s}  model = {obs[key]:.6g}   "
              f"target = {exp_targets[key]:.6g}   pull = {pulls[key]: .3f}")

    print("\n|V_CKM| ≈")
    print(np.abs(Vckm))

    print("\n|U_PMNS| ≈")
    print(np.abs(Upmns))

"""
====== CMA-ES Restart 1/8 ======
(10_w,20)-aCMA-ES (mu_w=5.9,w_1=27%) in dimension 31 (seed=467213, Tue Dec  9 13:53:04 2025)
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
    1     20 2.710790409087219e+03 1.0e+00 2.80e-01  3e-01  3e-01 0:00.0
    2     40 1.078155933078272e+03 1.1e+00 2.78e-01  3e-01  3e-01 0:00.0
    3     60 5.951967958143414e+02 1.1e+00 2.79e-01  3e-01  3e-01 0:00.0
  100   2000 6.733919156747852e+01 3.4e+00 1.82e-01  9e-02  2e-01 0:01.7
  200   4000 3.840849717871453e+01 7.8e+00 7.15e-02  2e-02  8e-02 0:03.3
Restart 1: best cost = 37.0025

====== CMA-ES Restart 2/8 ======
(10_w,20)-aCMA-ES (mu_w=5.9,w_1=27%) in dimension 31 (seed=383436, Tue Dec  9 13:53:08 2025)
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
    1     20 5.175344666124461e+03 1.0e+00 2.87e-01  3e-01  3e-01 0:00.0
    2     40 1.502163513194671e+03 1.1e+00 2.86e-01  3e-01  3e-01 0:00.0
    3     60 2.358094688359610e+03 1.1e+00 2.92e-01  3e-01  3e-01 0:00.0
  100   2000 3.397734471238520e+01 2.6e+00 7.00e-02  4e-02  8e-02 0:01.9
  200   4000 2.232859597366760e+01 7.0e+00 3.12e-02  2e-02  4e-02 0:03.6
Restart 2: best cost = 22.0763

====== CMA-ES Restart 3/8 ======
(10_w,20)-aCMA-ES (mu_w=5.9,w_1=27%) in dimension 31 (seed=445805, Tue Dec  9 13:53:11 2025)
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
    1     20 5.614005547673511e+03 1.0e+00 2.84e-01  3e-01  3e-01 0:00.0
    2     40 1.312654792812980e+03 1.1e+00 2.79e-01  3e-01  3e-01 0:00.0
    3     60 4.506966049247283e+02 1.1e+00 2.76e-01  3e-01  3e-01 0:00.1
  100   2000 1.465452433987385e+02 3.7e+00 3.88e-01  2e-01  4e-01 0:01.5
  200   4000 6.217162835883452e+01 5.8e+00 2.08e-01  6e-02  3e-01 0:03.1
Restart 3: best cost = 59.9983

====== CMA-ES Restart 4/8 ======
(10_w,20)-aCMA-ES (mu_w=5.9,w_1=27%) in dimension 31 (seed=460031, Tue Dec  9 13:53:14 2025)
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
    1     20 2.615656395677563e+03 1.0e+00 2.86e-01  3e-01  3e-01 0:00.0
    2     40 1.244841948173937e+03 1.1e+00 2.80e-01  3e-01  3e-01 0:00.0
    3     60 1.072465528836756e+03 1.1e+00 2.85e-01  3e-01  3e-01 0:00.0
  100   2000 7.246413838418371e+01 3.4e+00 3.34e-01  2e-01  4e-01 0:01.9
  200   4000 4.241524061949765e+01 5.5e+00 1.18e-01  5e-02  2e-01 0:03.6
Restart 4: best cost = 40.2619

====== CMA-ES Restart 5/8 ======
(10_w,20)-aCMA-ES (mu_w=5.9,w_1=27%) in dimension 31 (seed=431315, Tue Dec  9 13:53:18 2025)
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
    1     20 7.944043884832687e+02 1.0e+00 2.79e-01  3e-01  3e-01 0:00.0
    2     40 2.256724715040257e+03 1.1e+00 2.69e-01  3e-01  3e-01 0:00.0
    3     60 1.197116017231805e+03 1.1e+00 2.61e-01  3e-01  3e-01 0:00.1
  100   2000 2.798506469441742e+01 2.4e+00 8.17e-02  5e-02  9e-02 0:02.0
  200   4000 1.995828076825566e+01 4.8e+00 2.96e-02  2e-02  4e-02 0:03.3
Restart 5: best cost = 19.9583

====== CMA-ES Restart 6/8 ======
(10_w,20)-aCMA-ES (mu_w=5.9,w_1=27%) in dimension 31 (seed=448506, Tue Dec  9 13:53:21 2025)
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
    1     20 8.316306007141295e+03 1.0e+00 2.81e-01  3e-01  3e-01 0:00.0
    2     40 3.733690437181528e+03 1.1e+00 2.71e-01  3e-01  3e-01 0:00.0
    3     60 1.428730964549957e+03 1.1e+00 2.65e-01  3e-01  3e-01 0:00.0
  100   2000 4.430225549329395e+01 2.7e+00 9.59e-02  5e-02  1e-01 0:01.3
  200   4000 2.028284425164039e+01 7.8e+00 3.92e-02  1e-02  5e-02 0:03.4
Restar

"""