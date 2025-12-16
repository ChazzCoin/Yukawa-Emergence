#!/usr/bin/env python3
import numpy as np
from numpy.linalg import svd, eigvalsh, eigh, cond, solve, pinv

# ============================================================
# Fully emergent pipeline: N = 360, forbidden distance D = 7
# ============================================================
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
    0.5 * x_exp[0], 0.5 * x_exp[1], 0.5 * x_exp[2], 0.5 * x_exp[3],
    0.5 * x_exp[4], 0.5 * x_exp[5],
    0.1 * x_exp[6], 0.1 * x_exp[7], 0.1 * x_exp[8],
    0.1 * x_exp[9], 0.1 * x_exp[10], 0.1 * x_exp[11],
    0.3 * x_exp[12], 0.3 * x_exp[13]
])
observable_names = [
    # mass ratios
    "m_c/m_t",
    "m_u/m_t",
    "m_s/m_b",
    "m_d/m_b",
    "m_mu/m_tau",
    "m_e/m_tau",
    # CKM angles
    "theta12_q (rad)",
    "theta23_q (rad)",
    "theta13_q (rad)",
    # PMNS angles
    "theta12_l (rad)",
    "theta23_l (rad)",
    "theta13_l (rad)",
    # neutrino splittings
    "Delta m2_21 (eV^2)",
    "Delta m2_31 (eV^2)",
]

class Config:
    v = 246.0        # GeV
    mu0 = 1.0e12     # GeV
    mu_EW = 91.1876  # GeV
    Lambda_Maj = 7.0e13  # GeV (overall heavy Majorana scale)

    # Alignment scale: Fibonacci / 360
    # kappa = 360/89, eps = 1/kappa
    kappa = 360.0 / 89.0
    eps = 1.0 / kappa  # decay factor in (0,1)

    seed = 12345  # overwritten per run

    t0 = np.log(mu0)
    t1 = np.log(mu_EW)
    dt = -0.01  # log-scale step size

    # Higgs quartic (approx EW value, treated constant here)
    lam = 0.13

    # Optional: interpret κ as κ(t_align) = exp(-λ_align t_align)
    lambda_align = 1.0  # arbitrary units
    t_align = -np.log(eps) / lambda_align  # so exp(-λ_align t_align) = eps

class EmergentConfig:
    """
    Minimal config for the emergent pipeline.
    No eps, no κ, no alignment decay. Only physical scales.
    """
    v = 246.0        # GeV
    mu0 = 1.0e12     # GeV (seesaw / high scale)
    mu_EW = 91.1876  # GeV
    Lambda_Maj = 7.0e13  # GeV, overall heavy Majorana scale

    # log-scale RG evolution
    t0 = np.log(mu0)
    t1 = np.log(mu_EW)
    dt = -0.01  # < 0 to run downwards

    # Higgs quartic (kept constant)
    lam = 0.13

    # seed is not used internally here except for reproducibility
    seed = 12345


# ------------------------------------------------------------
# D = 7 structure: mask and adjacency
# ------------------------------------------------------------

def build_D7_mask(N: int, D: int = 7) -> np.ndarray:
    """
    Mask implementing the *only* structural rule:
      - entries with |i-j| == D are forbidden (set to 0)
      - all other entries (including diagonal) are allowed (set to 1)

    This mask is applied elementwise to Yukawa and Majorana proto matrices.
    """
    idx = np.arange(N)
    dist = np.abs(idx[:, None] - idx[None, :])
    M = np.ones((N, N), dtype=float)
    M[dist == D] = 0.0
    return M


def build_D7_adjacency(N: int, D: int = 7) -> np.ndarray:
    """
    Canonical Hermitian adjacency using only N and D:

      A_ij = 1  if i != j and |i-j| != D
            = 0  if i == j  or |i-j| == D

    This is the universal "geometry" operator whose spectrum we use
    to define the emergent flavor subspace.
    """
    idx = np.arange(N)
    dist = np.abs(idx[:, None] - idx[None, :])

    A = np.ones((N, N), dtype=float)
    A[dist == 0] = 0.0     # no self-edges
    A[dist == D] = 0.0     # forbidden distance D
    return A


# ------------------------------------------------------------
# Emergent flavor basis from D = 7 adjacency
# ------------------------------------------------------------

def emergent_flavor_basis(A: np.ndarray, n_flavors: int = 3):
    """
    Given a Hermitian adjacency A (N×N), define the emergent flavor basis as
    the n_flavors eigenvectors with the *largest* eigenvalues of A.

    This uses only:
      - the D=7 structure encoded in A
      - the canonical spectral ordering of a Hermitian operator.

    Returns:
      B: N×n_flavors matrix with orthonormal columns (flavor basis)
      evals_sel: the selected eigenvalues (for diagnostics)
    """
    evals, evecs = eigh(A)  # evals ascending, evecs columns
    idx = np.argsort(evals)[-n_flavors:]  # indices of largest eigenvalues
    B = evecs[:, idx]                     # N×3
    return B, evals[idx]


def project_to_flavor(Y_full: np.ndarray, B: np.ndarray) -> np.ndarray:
    """
    Project a full N×N matrix Y_full onto the emergent flavor subspace
    spanned by columns of B (N×3):

      Y^(3) = B† Y_full B.

    This is the emergent 3×3 Yukawa (or Majorana) matrix.
    """
    return B.conj().T @ Y_full @ B

def rescale_yukawa_sector(Y, v, m_target_heaviest):
    """
    Rescale Y so that the heaviest mass eigenvalue (v/√2 * max singular value)
    matches m_target_heaviest. Returns (Y_rescaled, alpha).
    """
    U_L, s, U_Rh = svd(Y)
    m_current = (v / np.sqrt(2.0)) * np.max(s)
    if m_current == 0:
        return Y, 1.0
    alpha = m_target_heaviest / m_current
    return alpha * Y, alpha

def random_complex_matrix(shape, rng):
    real = rng.normal(0.0, 1.0, size=shape)
    imag = rng.normal(0.0, 1.0, size=shape)
    return (real + 1j * imag) / np.sqrt(2.0)
def normalize_by_largest_singular_value(X: np.ndarray) -> np.ndarray:
    s = svd(X, compute_uv=False)
    s_max = np.max(s)
    if s_max == 0:
        return X
    return X / s_max

def seesaw_light_neutrinos(Ynu_eff: np.ndarray,
                           M_R: np.ndarray,
                           v: float,
                           cond_tol: float = 1e12) -> np.ndarray:
    """
    Type-I seesaw:
      m_D = v/√2 Ynu_eff,
      m_ν = - m_D M_R^{-1} m_D^T (symmetric 3x3, in GeV).
    """
    m_D = (v / np.sqrt(2.0)) * Ynu_eff
    if cond(M_R) > cond_tol:
        M_R_inv = pinv(M_R)
        m_nu = -m_D @ M_R_inv @ m_D.T
    else:
        X = solve(M_R, m_D.T)
        m_nu = -m_D @ X

    m_nu = 0.5 * (m_nu + m_nu.T)
    return m_nu

def has_bad(x: np.ndarray) -> bool:
    """Check for NaN or Inf in an array."""
    return np.any(np.isnan(x)) or np.any(np.isinf(x))

def run_RGE_full(Yu0, Yd0, Ye0, Ynu0, kappa_L0,
                 g1_const, g2_const, g3_const,
                 cfg: Config):
    """
    Evolve Yukawas and Weinberg operator from μ0 down to μ_EW
    with fixed gauge couplings.
    """
    Yu, Yd, Ye, Ynu = Yu0.copy(), Yd0.copy(), Ye0.copy(), Ynu0.copy()
    kappa_L = kappa_L0.copy()
    g1, g2, g3 = g1_const, g2_const, g3_const
    lam = cfg.lam

    t = cfg.t0
    step = 0
    while (cfg.dt < 0 and t > cfg.t1) or (cfg.dt > 0 and t < cfg.t1):
        step += 1

        Yu, Yd, Ye, Ynu, kappa_L = rk4_step_full(
            Yu, Yd, Ye, Ynu, kappa_L, g1, g2, g3, lam, cfg.dt
        )
        t += cfg.dt

        # Safeguard: Break if NaN/Inf detected (prevents crash)
        if any(has_bad(M) for M in (Yu, Yd, Ye, Ynu, kappa_L)):
            print(f"Warning: NaN/Inf detected at RGE step {step}, halting evolution.")
            break

    return Yu, Yd, Ye, Ynu, kappa_L, g1, g2, g3

def beta_Yukawas(Yu, Yd, Ye, Ynu, g1, g2, g3):
    """
    1-loop SM Yukawa RGEs (in matrix form), with fixed gauge couplings.
    """
    # Safeguard: Return zero betas if inputs contain NaN or Inf
    if any(has_bad(M) for M in (Yu, Yd, Ye, Ynu)):
        Z = np.zeros_like(Yu)
        return Z, Z, Z, Z

    Yu_dagYu = Yu.conj().T @ Yu
    Yd_dagYd = Yd.conj().T @ Yd
    Ye_dagYe = Ye.conj().T @ Ye
    Ynu_dagYnu = Ynu.conj().T @ Ynu

    T = np.trace(3 * Yu_dagYu + 3 * Yd_dagYd + Ye_dagYe)

    factor_u = T - (17 / 20 * g1 ** 2 + 9 / 4 * g2 ** 2 + 8 * g3 ** 2)
    factor_d = T - (1 / 4 * g1 ** 2 + 9 / 4 * g2 ** 2 + 8 * g3 ** 2)
    factor_e = T - (9 / 4 * g1 ** 2 + 9 / 4 * g2 ** 2)
    factor_nu = T - (9 / 20 * g1 ** 2 + 9 / 4 * g2 ** 2)

    dYu = Yu * factor_u + (3 / 2) * (Yu @ Yu_dagYu - Yd @ (Yd_dagYd @ Yu))
    dYd = Yd * factor_d + (3 / 2) * (Yd @ Yd_dagYd - Yu @ (Yu_dagYu @ Yd))
    dYe = Ye * factor_e + (3 / 2) * (Ye @ Ye_dagYe)
    dYnu = Ynu * factor_nu + (3 / 2) * (Ynu @ Ynu_dagYnu - Ye @ (Ye_dagYe @ Ynu))

    dYu /= (16 * np.pi ** 2)
    dYd /= (16 * np.pi ** 2)
    dYe /= (16 * np.pi ** 2)
    dYnu /= (16 * np.pi ** 2)

    return dYu, dYd, dYe, dYnu

def beta_kappa_L(kappa_L, Yu, Ye, g2, lam):
    """
    16π² dκ_L/dt = (-3 g2² + 2λ + 6 Tr(Yu†Yu)) κ_L
                   - 3/2 (Ye†Ye κ_L + κ_L (Ye†Ye)^T).

    We treat λ as constant, and ignore g1,g3 in this operator.
    """
    if any(has_bad(M) for M in (kappa_L, Yu, Ye)):
        return np.zeros_like(kappa_L)

    Yu_dagYu = Yu.conj().T @ Yu
    Ye_dagYe = Ye.conj().T @ Ye
    T_u = np.trace(Yu_dagYu)

    pref = (-3 * g2 ** 2 + 2 * lam + 6 * T_u)
    term1 = pref * kappa_L
    term2 = -1.5 * (Ye_dagYe @ kappa_L + kappa_L @ Ye_dagYe.T.conj())

    dkappa = (term1 + term2) / (16 * np.pi ** 2)
    return dkappa

def rk4_step_full(Yu, Yd, Ye, Ynu, kappa_L, g1, g2, g3, lam, dt):
    """
    RK4 step evolving Yukawas + κ_L with fixed (g1,g2,g3,lam).
    """
    # k1
    dYu1, dYd1, dYe1, dYnu1 = beta_Yukawas(Yu, Yd, Ye, Ynu, g1, g2, g3)
    dkL1 = beta_kappa_L(kappa_L, Yu, Ye, g2, lam)

    # k2
    Yu2 = Yu + 0.5 * dt * dYu1
    Yd2 = Yd + 0.5 * dt * dYd1
    Ye2 = Ye + 0.5 * dt * dYe1
    Ynu2 = Ynu + 0.5 * dt * dYnu1
    kL2 = kappa_L + 0.5 * dt * dkL1

    dYu2, dYd2, dYe2, dYnu2 = beta_Yukawas(Yu2, Yd2, Ye2, Ynu2, g1, g2, g3)
    dkL2 = beta_kappa_L(kL2, Yu2, Ye2, g2, lam)

    # k3
    Yu3 = Yu + 0.5 * dt * dYu2
    Yd3 = Yd + 0.5 * dt * dYd2
    Ye3 = Ye + 0.5 * dt * dYe2
    Ynu3 = Ynu + 0.5 * dt * dYnu2
    kL3 = kappa_L + 0.5 * dt * dkL2

    dYu3, dYd3, dYe3, dYnu3 = beta_Yukawas(Yu3, Yd3, Ye3, Ynu3, g1, g2, g3)
    dkL3 = beta_kappa_L(kL3, Yu3, Ye3, g2, lam)

    # k4
    Yu4 = Yu + dt * dYu3
    Yd4 = Yd + dt * dYd3
    Ye4 = Ye + dt * dYe3
    Ynu4 = Ynu + dt * dYnu3
    kL4 = kappa_L + dt * dkL3

    dYu4, dYd4, dYe4, dYnu4 = beta_Yukawas(Yu4, Yd4, Ye4, Ynu4, g1, g2, g3)
    dkL4 = beta_kappa_L(kL4, Yu4, Ye4, g2, lam)

    Yu_next = Yu + (dt / 6.0) * (dYu1 + 2 * dYu2 + 2 * dYu3 + dYu4)
    Yd_next = Yd + (dt / 6.0) * (dYd1 + 2 * dYd2 + 2 * dYd3 + dYd4)
    Ye_next = Ye + (dt / 6.0) * (dYe1 + 2 * dYe2 + 2 * dYe3 + dYe4)
    Ynu_next = Ynu + (dt / 6.0) * (dYnu1 + 2 * dYnu2 + 2 * dYnu3 + dYnu4)
    kL_next = kappa_L + (dt / 6.0) * (dkL1 + 2 * dkL2 + 2 * dkL3 + dkL4)

    return Yu_next, Yd_next, Ye_next, Ynu_next, kL_next

def extract_angles_and_phase(V: np.ndarray):
    """
    Extract approximate mixing angles (in radians) and Dirac phase
    from a 3x3 unitary matrix V, assuming a PDG-like parameterization.
    """
    s13 = np.abs(V[0, 2])
    theta13 = np.arcsin(np.clip(s13, 0.0, 1.0))

    s12 = np.abs(V[0, 1])
    c12 = np.abs(V[0, 0])
    theta12 = np.arctan2(s12, c12)

    s23 = np.abs(V[1, 2])
    c23 = np.abs(V[2, 2])
    theta23 = np.arctan2(s23, c23)

    J = np.imag(V[0, 0] * V[1, 1] * np.conj(V[0, 1]) * np.conj(V[1, 0]))
    denom = (np.sin(2 * theta12) * np.sin(2 * theta23) *
             np.sin(2 * theta13) * np.cos(theta13))
    if np.abs(denom) < 1e-12:
        delta = 0.0
    else:
        x = 8 * J / denom
        x = np.clip(x, -1.0, 1.0)
        delta = np.arcsin(x)

    return theta12, theta23, theta13, delta

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
    obs.append(mu_sorted[1] / mu_sorted[2])  # m_c/m_t
    obs.append(mu_sorted[0] / mu_sorted[2])  # m_u/m_t
    obs.append(md_sorted[1] / md_sorted[2])  # m_s/m_b
    obs.append(md_sorted[0] / md_sorted[2])  # m_d/m_b
    obs.append(me_sorted[1] / me_sorted[2])  # m_mu/m_tau
    obs.append(me_sorted[0] / me_sorted[2])  # m_e/m_tau

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

def chi2_breakdown(res):
    """
    Return per-observable χ² contributions as a list of dicts:
      {"name", "theory", "exp", "sigma", "chi2_i"}
    """
    x_th = make_observables(res)
    diffs = x_th - x_exp
    chi2_i = (diffs / sigma) ** 2

    breakdown = []
    for name, th, exp, sig, c2 in zip(observable_names, x_th, x_exp, sigma, chi2_i):
        breakdown.append({
            "name": name,
            "theory": th,
            "exp": exp,
            "sigma": sig,
            "chi2_i": c2,
        })
    return breakdown


def chi2(observed, expected, sigma):
    return np.sum(((observed - expected) / sigma) ** 2)

def chi2_from_res(res):
    x_th = make_observables(res)
    return chi2(x_th, x_exp, sigma)

def neutrino_splittings(mnu_masses: np.ndarray):
    """
    Compute Δm²_21 and Δm²_31 in GeV² from the (non-negative) Takagi singular values.
    """
    m_sorted = np.sort(mnu_masses)
    m1, m2, m3 = m_sorted
    dm2_21 = m2 ** 2 - m1 ** 2
    dm2_31 = m3 ** 2 - m1 ** 2
    return dm2_21, dm2_31  # GeV^2

def diag_dirac_Y(Y: np.ndarray, v: float):
    """
    SVD for Dirac Yukawa:
      Y = U_L diag(s) U_R†,  masses = v/√2 * s.
    """
    U_L, s, U_Rh = svd(Y)
    masses = (v / np.sqrt(2.0)) * s
    return U_L, np.diag(s), U_Rh, masses

def takagi_symmetric(m: np.ndarray):
    """
    Takagi factorization via SVD for complex symmetric 3x3:
      m = U diag(s) U^T, with s ≥ 0.
    """
    U, s, Vh = svd(m)
    return U, s

def diagonalize_all(Yu, Yd, Ye, mnu, v):
    UuL, Yu_diag, UuR, mu = diag_dirac_Y(Yu, v)
    UdL, Yd_diag, UdR, md = diag_dirac_Y(Yd, v)
    UlL, Ye_diag, UlR, me = diag_dirac_Y(Ye, v)

    U_nu, mnu_vals = takagi_symmetric(mnu)
    mnu_masses = mnu_vals  # in GeV

    Vckm = UuL.conj().T @ UdL
    Vpmns = UlL.conj().T @ U_nu

    return mu, md, me, mnu_masses, Vckm, Vpmns
# ------------------------------------------------------------
# Emergent pipeline
# ------------------------------------------------------------

def run_pipeline_emergent(seed: int,
                          cfg: EmergentConfig,
                          N_sites: int = 360,
                          D_forbidden: int = 7,
                          use_RGE: bool = True):
    """
    Fully emergent pipeline:

      - internal proto space: C^N with N = N_sites
      - only structural input: forbidden distance |i-j| = D_forbidden
      - random O(1) proto Yukawas and Majorana on N sites
      - enforce D=7 mask on these matrices
      - build D=7 adjacency A and find its top-3 eigenmodes → emergent 3D flavor basis B
      - project all sectors: Y^(3) = B† Y_full B, M_R^(3) = B† M_full B
      - seesaw to get mν(μ0)
      - build Weinberg operator κ_L(μ0)
      - optionally run 1-loop RGEs down to μ_EW
      - rescale Yukawas to match m_t, m_b, m_τ
      - diagonalize and extract masses, mixings, Δm²
      - compute χ² vs your existing x_exp, sigma

    Assumes the following functions from your original code are available:
      - random_complex_matrix, normalize_by_largest_singular_value
      - seesaw_light_neutrinos
      - beta_Yukawas, beta_kappa_L, rk4_step_full, run_RGE_full
      - diag_dirac_Y, takagi_symmetric, diagonalize_all
      - neutrino_splittings, make_observables, chi2_from_res
    """
    rng = np.random.default_rng(seed)

    N = N_sites
    D = D_forbidden

    # 1. D=7 structure: mask for couplings, adjacency for emergent modes
    mask_D7 = build_D7_mask(N, D)
    A_D7 = build_D7_adjacency(N, D)

    # 2. Random proto matrices (O(1) complex, normalized)
    Yu_full = random_complex_matrix((N, N), rng)
    Yu_full = normalize_by_largest_singular_value(Yu_full)

    Yd_full = random_complex_matrix((N, N), rng)
    Yd_full = normalize_by_largest_singular_value(Yd_full)

    Ye_full = random_complex_matrix((N, N), rng)
    Ye_full = normalize_by_largest_singular_value(Ye_full)

    Ynu_full = random_complex_matrix((N, N), rng)
    Ynu_full = normalize_by_largest_singular_value(Ynu_full)

    M_full = random_complex_matrix((N, N), rng)
    M_full = normalize_by_largest_singular_value(M_full)
    M_full = 0.5 * (M_full + M_full.T)  # symmetric for Majorana

    # 3. Enforce the D=7 constraint on all sectors
    Yu_full *= mask_D7
    Yd_full *= mask_D7
    Ye_full *= mask_D7
    Ynu_full *= mask_D7
    M_full *= mask_D7

    # 4. Emergent 3D flavor basis from D=7 adjacency
    B, evals_flavor = emergent_flavor_basis(A_D7, n_flavors=3)
    # B: N×3, columns orthonormal

    # 5. Project all sectors to 3×3 in flavor space
    Yu_eff = project_to_flavor(Yu_full, B)
    Yd_eff = project_to_flavor(Yd_full, B)
    Ye_eff = project_to_flavor(Ye_full, B)
    Ynu_eff = project_to_flavor(Ynu_full, B)

    M3 = project_to_flavor(M_full, B)
    M3 = 0.5 * (M3 + M3.T)  # enforce symmetry at the 3×3 level
    M_R = cfg.Lambda_Maj * M3

    # 6. Seesaw at μ0 → mν(μ0) in GeV (3×3)
    m_nu_0 = seesaw_light_neutrinos(Ynu_eff, M_R, cfg.v)

    # 7. Weinberg operator κ_L(μ0) (dimensionful, GeV^-1)
    kappa_L_0 = (2.0 / cfg.v ** 2) * m_nu_0

    # 8. RG evolution (3×3 Yukawas + κ_L), same as your original code
    g1_0, g2_0, g3_0 = 0.46, 0.63, 0.88

    if use_RGE:
        Yu_EW, Yd_EW, Ye_EW, Ynu_EW, kappa_L_EW, g1_EW, g2_EW, g3_EW = run_RGE_full(
            Yu_eff, Yd_eff, Ye_eff, Ynu_eff, kappa_L_0, g1_0, g2_0, g3_0, cfg
        )

        # If something went wrong in RGE, bail out with a bad χ²
        if any(has_bad(M) for M in (Yu_EW, Yd_EW, Ye_EW, Ynu_EW, kappa_L_EW)):
            res = {
                "mu": np.full(3, np.nan),
                "md": np.full(3, np.nan),
                "me": np.full(3, np.nan),
                "mnu": np.full(3, np.nan),
                "th_q": (np.nan, np.nan, np.nan),
                "th_l": (np.nan, np.nan, np.nan),
                "dm2_eV2": (np.nan, np.nan),
                "chi2": np.inf,
                "evals_flavor": evals_flavor,
            }
            return res

        # reconstruct mν(μ_EW) from κ_L(μ_EW)
        m_nu_EW = 0.5 * cfg.v ** 2 * kappa_L_EW
    else:
        Yu_EW, Yd_EW, Ye_EW, Ynu_EW = Yu_eff, Yd_eff, Ye_eff, Ynu_eff
        m_nu_EW = m_nu_0
        g1_EW, g2_EW, g3_EW = g1_0, g2_0, g3_0

    # 9. Rescale sectors to fix heavy Dirac masses at EW scale
    m_t_target = 173.0   # GeV
    m_b_target = 4.18    # GeV
    m_tau_target = 1.777 # GeV

    Yu_EW, alpha_u = rescale_yukawa_sector(Yu_EW, cfg.v, m_t_target)
    Yd_EW, alpha_d = rescale_yukawa_sector(Yd_EW, cfg.v, m_b_target)
    Ye_EW, alpha_e = rescale_yukawa_sector(Ye_EW, cfg.v, m_tau_target)

    # 10. Diagonalize at μ_EW (3×3)
    mu, md, me, mnu_masses, Vckm, Vpmns = diagonalize_all(
        Yu_EW, Yd_EW, Ye_EW, m_nu_EW, cfg.v
    )

    # 11. Angles and Δm²
    th12_q, th23_q, th13_q, delta_q = extract_angles_and_phase(Vckm)
    th12_l, th23_l, th13_l, delta_l = extract_angles_and_phase(Vpmns)

    dm2_21_GeV2, dm2_31_GeV2 = neutrino_splittings(mnu_masses)
    dm2_21_eV2 = dm2_21_GeV2 * 1e18
    dm2_31_eV2 = dm2_31_GeV2 * 1e18

    res = {
        "mu": mu,
        "md": md,
        "me": me,
        "mnu": mnu_masses,  # GeV
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
        "evals_flavor": evals_flavor,  # the three selected eigenvalues of A_D7
    }

    res["chi2"] = chi2_from_res(res)
    return res


# ------------------------------------------------------------
# Simple emergent scan over seeds
# ------------------------------------------------------------

def emergent_seed_scan(N_seeds: int = 10,
                       N_sites: int = 360,
                       D_forbidden: int = 7):
    """
    Scan over random seeds using the fully emergent pipeline.
    Prints χ² and basic spectra, and reports the best seed.
    """
    cfg = EmergentConfig()

    all_results = []
    chi2_vals = []

    for seed in range(N_seeds):
        res = run_pipeline_emergent(
            seed, cfg,
            N_sites=N_sites,
            D_forbidden=D_forbidden,
            use_RGE=True,
        )
        chi2_vals.append(res["chi2"])
        all_results.append(res)
        print(f"[emergent] seed {seed}: chi2 = {res['chi2']:.3g}")

    best_idx = int(np.argmin(chi2_vals))
    best = all_results[best_idx]

    print("\n[emergent] Best seed:", best_idx)
    print("chi2 =", best["chi2"])
    print("Up masses (GeV):      ", best["mu"])
    print("Down masses (GeV):    ", best["md"])
    print("Lepton masses (GeV):  ", best["me"])
    print("Neutrino masses (GeV):", best["mnu"])
    print("Neutrino masses (eV): ", best["mnu"] * 1e9)
    print("Δm² (eV²):            ", best["dm2_eV2"])
    print("CKM angles (rad):     ", best["th_q"], "δq:", best["delta_q"])
    print("PMNS angles (rad):    ", best["th_l"], "δℓ:", best["delta_l"])
    print("Flavor eigenvalues of A_D7:", best["evals_flavor"])

    print("\n=== χ² breakdown for best emergent seed ===")
    breakdown = chi2_breakdown(best)
    for entry in breakdown:
        name = entry["name"]
        th = entry["theory"]
        exp = entry["exp"]
        sig = entry["sigma"]
        c2 = entry["chi2_i"]
        pull = (th - exp) / sig
        print(f"{name:20s}  th = {th: .4e},  exp = {exp: .4e},  "
              f"sigma = {sig: .4e},  pull = {pull: .2f},  chi2_i = {c2: .2f}")


if __name__ == "__main__":
    # You can comment out the old scans if you like.
    # simple_seed_scan()
    # structural_scan()

    # New fully emergent scan:
    emergent_seed_scan()
