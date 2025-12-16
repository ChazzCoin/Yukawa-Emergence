#!/usr/bin/env python3
import numpy as np
import cma
from scipy.integrate import solve_ivp

# =========================
# Global settings / targets
# =========================

# Use RGE in the cost function. You can turn this off for speed while debugging.
USE_RGE = True

# Scales (GeV)
MU_HIGH = 2.0e14   # "flavor / seesaw" scale
MU_LOW  = 1.0e2    # electroweak-ish scale

# Rough gauge couplings and Higgs quartic at EW scale
g1_EW    = 0.36
g2_EW    = 0.65
g3_EW    = 1.17
lambda_H = 0.13

# Higgs vev
v_higgs = 246.0  # GeV

# Base exponents for KAPPA hierarchy
EXP_U_BASE  = np.array([4.0, 2.0, 0.0])   # up
EXP_D_BASE  = np.array([3.0, 2.0, 0.0])   # down
EXP_E_BASE  = np.array([3.0, 2.0, 0.0])   # charged leptons
EXP_NU_BASE = np.array([1.0, 0.0, 0.0])   # Dirac neutrinos

KAPPA = 0.24

# Experimental targets (your previous 14 + CP phases and Jarlskogs)
exp_targets = {
    "m_c/m_t"    : 0.007,
    "m_u/m_t"    : 1.0e-5,
    "m_s/m_b"    : 0.02,
    "m_d/m_b"    : 0.001,
    "m_mu/m_tau" : 0.06,
    "m_e/m_tau"  : 3.0e-4,

    "theta12_q"  : 0.226,
    "theta23_q"  : 0.041,
    "theta13_q"  : 0.0035,

    "theta12_l"  : 0.59,
    "theta23_l"  : 0.84,
    "theta13_l"  : 0.15,

    "Delta m2_21": 7.4e-5,
    "Delta m2_31": 2.5e-3,

    # CP phases and Jarlskogs (rough targets; refine as you like)
    "delta_CKM"  : 1.20,     # ~69°
    "J_CKM"      : 3.0e-5,
    "delta_PMNS" : -1.57,    # ~ -π/2
    "J_PMNS"     : 0.03,
}

# Default 30% fractional uncertainty for dimensionful/ratio observables,
# overridden with absolute sigmas for phases and J.
sigma_targets = {}
for key, val in exp_targets.items():
    sigma_targets[key] = 0.3 * abs(val) if val != 0 else 0.3

# Override sigmas for phases and Jarlskogs
sigma_targets["delta_CKM"]  = 0.10      # ~6 deg
sigma_targets["J_CKM"]      = 0.5e-5
sigma_targets["delta_PMNS"] = 0.30
sigma_targets["J_PMNS"]     = 0.01


# =========================
# Utility: flavor geometry
# =========================

def build_phase_profile(A, B):
    """Phase per generation: phi_g = A + B*g, g=0,1,2."""
    return np.array([A + B * g for g in range(3)], dtype=float)


def build_site_phases(phi_gen, n_sites=9):
    """Assign phase per site by generation index g(i) = i mod 3."""
    phases = np.zeros(n_sites, dtype=float)
    for i in range(n_sites):
        phases[i] = phi_gen[i % 3]
    return phases


def build_phase_matrix(site_phases):
    """P_ij = exp(i (phi_i - phi_j))."""
    phi_i = site_phases.reshape(-1, 1)
    phi_j = site_phases.reshape(1, -1)
    return np.exp(1j * (phi_i - phi_j))


def site_scales(base_exponents, shifts):
    """
    Given base exponents [a0,a1,a2] and shifts [X,Y], return scale per site.
    Generation exponents:
      e0 = a0
      e1 = a1 + X
      e2 = a2 + Y
    Then s_g = KAPPA^{e_g}, s_i = s_{g(i)}.
    """
    e0 = base_exponents[0]
    e1 = base_exponents[1] + shifts[0]
    e2 = base_exponents[2] + shifts[1]
    e = np.array([e0, e1, e2], dtype=float)
    s_g = KAPPA**e
    s_sites = np.zeros(9, dtype=float)
    for i in range(9):
        s_sites[i] = s_g[i % 3]
    return s_sites


def build_kernel(n_sites=9, gamma=0.0):
    """
    Toeplitz kernel on 9-ring with forbidden distance d=2.
    K_ii = 1
    K_ij = 0 if d(i,j)==2
         = exp(-gamma * d) otherwise
    """
    K = np.zeros((n_sites, n_sites), dtype=float)
    for i in range(n_sites):
        for j in range(n_sites):
            if i == j:
                K[i, j] = 1.0
            else:
                d = abs(i - j)
                d = min(d, n_sites - d)
                if d == 2:
                    K[i, j] = 0.0
                else:
                    K[i, j] = np.exp(-gamma * d)
    return K


def build_Y_sector(A, B, base_exponents, shifts, gamma):
    """
    Build 9×9 Yukawa-like matrix for one sector:
      Y_ij ~ s_i s_j * exp(i(phi_i - phi_j)) * K_ij
    Then rescale so largest singular value = 1 (dimensionless hierarchy).
    """
    phi_gen = build_phase_profile(A, B)
    site_phi = build_site_phases(phi_gen, n_sites=9)
    P = build_phase_matrix(site_phi)
    s = site_scales(base_exponents, shifts)
    mag = np.outer(s, s)  # s_i s_j
    K = build_kernel(n_sites=9, gamma=gamma)

    Y = mag * P * K  # elementwise

    # Normalize by largest singular value
    u, sing, vh = np.linalg.svd(Y)
    max_sv = np.max(sing)
    if max_sv > 0:
        Y /= max_sv
    return Y


# =========================
# Schur 9→3 & Cabibbo
# =========================

def schur_9to3(Y9):
    """
    Reduce 9×9 Yukawa to effective 3×3 via Schur complement
    over heavy sites (last 6 indices).
    """
    A = Y9[0:3, 0:3]
    B = Y9[0:3, 3:9]
    C = Y9[3:9, 0:3]
    D = Y9[3:9, 3:9]

    # Schur complement: A_eff = A - B D^{-1} C
    # Add small ridge if D is nearly singular
    D_reg = D.copy()
    D_reg += 1e-8 * np.eye(D_reg.shape[0])

    D_inv = np.linalg.inv(D_reg)
    Y3 = A - B @ D_inv @ C
    return Y3


def cabibbo_rotation(theta_C):
    """
    Real rotation in 1–2 plane acting on left-handed down quarks.
    """
    c = np.cos(theta_C)
    s = np.sin(theta_C)
    U = np.eye(3, dtype=float)
    U[0, 0] = c
    U[0, 1] = s
    U[1, 0] = -s
    U[1, 1] = c
    return U


# =========================
# Neutrino projection & seesaw
# =========================

def build_projection_resonant(lambda_nu):
    """
    3×9 projector onto triadic modes with a resonance angle lambda_nu.
    Start from three "block" vectors supported on sites (0,3,6), (1,4,7), (2,5,8),
    then rotate first two by lambda_nu in their 3d subspace.
    """
    n = 9
    B0 = np.zeros((3, n), dtype=float)

    # Row 0: sites 0,3,6
    for k in (0, 3, 6):
        B0[0, k] = 1.0
    # Row 1: sites 1,4,7
    for k in (1, 4, 7):
        B0[1, k] = 1.0
    # Row 2: sites 2,5,8
    for k in (2, 5, 8):
        B0[2, k] = 1.0

    # Normalize rows
    B0 /= np.sqrt(3.0)

    # Rotate rows 0 and 1 by lambda_nu
    R = np.eye(3)
    c = np.cos(lambda_nu)
    s = np.sin(lambda_nu)
    R[0, 0] = c
    R[0, 1] = s
    R[1, 0] = -s
    R[1, 1] = c

    P = R @ B0  # 3×9
    # Rows remain orthonormal
    return P


def triadic_seesaw(M0, Ynu_eff, P):
    """
    Simple triadic seesaw:
      - M0: 9×9 heavy Majorana matrix
      - P : 3×9 projector
      - Ynu_eff: 3×3 Dirac Yukawa in triadic subspace
      => M_R_eff = P M0 P^T
      => M_nu = - (v^2 / 2) Ynu_eff M_R_eff^{-1} Ynu_eff^T
    """
    M_R = P @ M0 @ P.T
    # Regularize if needed
    M_R_reg = M_R + 1e-5 * np.eye(3)
    M_R_inv = np.linalg.inv(M_R_reg)

    Mnu = -0.5 * (v_higgs**2) * (Ynu_eff @ M_R_inv @ Ynu_eff.T)
    return Mnu


# =========================
# Proto Majorana
# =========================

def proto_majorana(rng, scale=7e13):
    """
    Random 9×9 complex symmetric Majorana matrix with overall scale ~ 'scale'.
    """
    M_re = rng.normal(size=(9, 9))
    M_im = rng.normal(size=(9, 9))
    M = M_re + 1j * M_im
    M = 0.5 * (M + M.T)  # symmetrize

    # Normalize spectral norm and rescale to 'scale'
    u, s, vh = np.linalg.svd(M)
    max_sv = np.max(s)
    if max_sv > 0:
        M *= (scale / max_sv)
    return M


# =========================
# Observables: angles, CP
# =========================

def mixing_angles_from_U(U):
    """
    Extract (theta12, theta23, theta13) from a 3×3 unitary U
    assuming PDG-like parameterization.
    """
    Uabs = np.abs(U)
    s13 = Uabs[0, 2]
    c13 = np.sqrt(max(0.0, 1.0 - s13**2))
    s12 = Uabs[0, 1] / c13
    s23 = Uabs[1, 2] / c13

    # Clip to avoid numerical nonsense
    s12 = np.clip(s12, -1.0, 1.0)
    s23 = np.clip(s23, -1.0, 1.0)

    theta12 = np.arcsin(s12)
    theta23 = np.arcsin(s23)
    theta13 = np.arcsin(np.clip(s13, -1.0, 1.0))

    return theta12, theta23, theta13


def jarlskog_from_U(U):
    return np.imag(U[0, 0] * U[1, 1] * np.conj(U[0, 1]) * np.conj(U[1, 0]))


def delta_from_U(U):
    """
    Extract Dirac CP phase δ from U via Jarlskog invariant.
    """
    Uabs = np.abs(U)

    s13 = Uabs[0, 2]
    c13 = np.sqrt(max(0.0, 1.0 - s13**2))

    s12 = Uabs[0, 1] / c13
    c12 = Uabs[0, 0] / c13

    s23 = Uabs[1, 2] / c13
    c23 = Uabs[2, 2] / c13

    J = jarlskog_from_U(U)

    denom = s12 * s23 * s13 * c12 * c23 * (c13**2)
    if np.abs(denom) < 1e-12:
        return 0.0, J

    sin_delta = J / denom
    sin_delta = np.clip(sin_delta, -1.0, 1.0)
    delta = np.arcsin(sin_delta)  # branch [-π/2, π/2]

    return delta, J


def compute_observables(Yu, Yd, Ye, Mnu):
    """
    Diagonalize Yukawas and Majorana matrix, compute mass ratios,
    mixing angles, mass splittings, and CP phases.
    """
    # SVD for Yukawas: Y = U diag(s) V†
    Uu, su, Vu = np.linalg.svd(Yu)
    Ud, sd, Vd = np.linalg.svd(Yd)
    Ue, se, Ve = np.linalg.svd(Ye)

    # Sort singular values ascending (lightest→heaviest)
    su_sorted = np.sort(su)
    sd_sorted = np.sort(sd)
    se_sorted = np.sort(se)

    # Mass ratios (independent of overall normalization)
    mu, mc, mt = su_sorted  # assume su_sorted[2] heaviest
    md, ms, mb = sd_sorted
    me, mmu, mtau = se_sorted

    obs = {}
    obs["m_c/m_t"]    = mc / mt
    obs["m_u/m_t"]    = mu / mt
    obs["m_s/m_b"]    = ms / mb
    obs["m_d/m_b"]    = md / mb
    obs["m_mu/m_tau"] = mmu / mtau
    obs["m_e/m_tau"]  = me / mtau

    # CKM
    Vckm = Uu.conj().T @ Ud
    th12q, th23q, th13q = mixing_angles_from_U(Vckm)
    obs["theta12_q"] = float(th12q)
    obs["theta23_q"] = float(th23q)
    obs["theta13_q"] = float(th13q)

    # Neutrino Majorana: diagonalize symmetric Mnu
    M_sym = 0.5 * (Mnu + Mnu.T)
    evals, U_nu = np.linalg.eigh(M_sym)
    # Ensure non-negative masses
    m_nu = np.abs(evals)
    # Sort ascending: m1 <= m2 <= m3
    idx = np.argsort(m_nu)
    m1, m2, m3 = m_nu[idx]
    U_nu = U_nu[:, idx]

    # PMNS = Ue† U_nu
    U_pmns = Ue.conj().T @ U_nu

    th12l, th23l, th13l = mixing_angles_from_U(U_pmns)
    obs["theta12_l"] = float(th12l)
    obs["theta23_l"] = float(th23l)
    obs["theta13_l"] = float(th13l)

    # Neutrino splittings
    obs["Delta m2_21"] = float(m2**2 - m1**2)
    obs["Delta m2_31"] = float(m3**2 - m1**2)

    # CP phases & Jarlskogs
    delta_ckm, J_ckm = delta_from_U(Vckm)
    delta_pmns, J_pmns = delta_from_U(U_pmns)

    obs["delta_CKM"]  = float(delta_ckm)
    obs["J_CKM"]      = float(J_ckm)
    obs["delta_PMNS"] = float(delta_pmns)
    obs["J_PMNS"]     = float(J_pmns)

    return obs, Vckm, U_pmns


def chi2_from_obs(obs):
    chi2 = 0.0
    pulls = {}
    for key, target in exp_targets.items():
        val = obs[key]
        sigma = sigma_targets[key]
        pull = (val - target) / sigma
        pulls[key] = pull
        chi2 += pull**2
    return chi2, pulls


# =========================
# RGE machinery (1-loop SM)
# =========================

def pack_matrices(Yu, Yd, Ye, kappa):
    """
    Flatten complex 3×3 matrices into a real vector for ODE solver.
    """
    def pack(M):
        return np.concatenate([M.real.flatten(), M.imag.flatten()])
    return np.concatenate([pack(Yu), pack(Yd), pack(Ye), pack(kappa)])


def unpack_matrices(vec):
    n = 3
    N = n * n
    def unpack_block(start):
        re = vec[start:start+N].reshape((n, n))
        im = vec[start+N:start+2*N].reshape((n, n))
        return re + 1j * im

    Yu = unpack_block(0)
    Yd = unpack_block(2*N)
    Ye = unpack_block(4*N)
    kappa = unpack_block(6*N)
    return Yu, Yd, Ye, kappa


def beta_Yukawa_system(t, y, g1, g2, g3, lam_H):
    Yu, Yd, Ye, kappa = unpack_matrices(y)

    # Trace T
    T = np.trace(3 * Yu.conj().T @ Yu +
                 3 * Yd.conj().T @ Yd +
                 Ye.conj().T @ Ye).real

    pref = 1.0 / (16.0 * np.pi**2)

    dYu = (
        Yu * (T - (17.0/20.0)*g1**2 - (9.0/4.0)*g2**2 - 8.0*g3**2)
        + 1.5 * Yu @ Yu.conj().T @ Yu
        - 1.5 * Yd @ Yd.conj().T @ Yu
    )

    dYd = (
        Yd * (T - (1.0/4.0)*g1**2 - (9.0/4.0)*g2**2 - 8.0*g3**2)
        + 1.5 * Yd @ Yd.conj().T @ Yd
        - 1.5 * Yu @ Yu.conj().T @ Yd
    )

    dYe = (
        Ye * (T - (9.0/4.0)*g1**2 - (9.0/4.0)*g2**2)
        + 1.5 * Ye @ Ye.conj().T @ Ye
    )

    YeYe = Ye @ Ye.conj().T
    dKappa = (
        (-3.0 * g2**2 + lam_H) * kappa
        + (YeYe @ kappa + kappa @ YeYe.T)
    )

    return pack_matrices(pref * dYu, pref * dYd, pref * dYe, pref * dKappa)


def run_rge(mu_high, mu_low, Yu_high, Yd_high, Ye_high, kappa_high,
            g1, g2, g3, lam_H):
    """
    Run 1-loop SM RGEs (Yukawas + Weinberg operator) from mu_high down to mu_low.
    Gauge couplings and lambda_H are held fixed here.
    """
    t_high = np.log(mu_high)
    t_low  = np.log(mu_low)

    y0 = pack_matrices(Yu_high, Yd_high, Ye_high, kappa_high)

    sol = solve_ivp(
        beta_Yukawa_system,
        t_span=(t_high, t_low),
        y0=y0,
        args=(g1, g2, g3, lam_H),
        rtol=1e-5,
        atol=1e-7,
        method='RK45'
    )

    Yu_low, Yd_low, Ye_low, kappa_low = unpack_matrices(sol.y[:, -1])
    return Yu_low, Yd_low, Ye_low, kappa_low


# =========================
# Parameters: unpack / cost
# =========================

def unpack_params_16C(X):
    """
    X (len 17):
      [A_u, B_u,
       A_d, B_d,
       A_nu, B_nu,
       shift_u1, shift_u2,
       shift_d1, shift_d2,
       shift_nu1, shift_nu2,
       shift_e1, shift_e2,
       lambda_nu,
       theta_C,
       gamma_l]
    """
    X = np.asarray(X, dtype=float)

    A_u, B_u = X[0], X[1]
    A_d, B_d = X[2], X[3]
    A_nu, B_nu = X[4], X[5]

    shifts_u  = np.array([X[6],  X[7]])
    shifts_d  = np.array([X[8],  X[9]])
    shifts_nu = np.array([X[10], X[11]])
    shifts_e  = np.array([X[12], X[13]])

    lambda_nu = X[14]
    theta_C   = X[15]
    gamma_l   = X[16]

    return (A_u, B_u,
            A_d, B_d,
            A_nu, B_nu,
            shifts_u, shifts_d, shifts_nu, shifts_e,
            lambda_nu, theta_C, gamma_l)


def resonant16C_cost(X, M0):
    """
    Full cost:
      - build flavor at MU_HIGH,
      - project & seesaw,
      - optionally run RGE to MU_LOW,
      - compute observables + chi²,
      - add regularization on shifts, lambda_nu, theta_C, gamma_l.
    """
    try:
        (A_u, B_u,
         A_d, B_d,
         A_nu, B_nu,
         shifts_u, shifts_d, shifts_nu, shifts_e,
         lambda_nu, theta_C, gamma_l) = unpack_params_16C(X)

        gamma_q = 0.0  # keep quark coherence fully on

        # 9×9 Yukawas at high scale
        Yu9  = build_Y_sector(A_u,  B_u,  EXP_U_BASE,  shifts_u,  gamma_q)
        Yd9  = build_Y_sector(A_d,  B_d,  EXP_D_BASE,  shifts_d,  gamma_q)
        Ye9  = build_Y_sector(A_d,  B_d,  EXP_E_BASE,  shifts_e,  gamma_l)
        Ynu9 = build_Y_sector(A_nu, B_nu, EXP_NU_BASE, shifts_nu, gamma_l)

        # 9→3 Schur at high scale
        Yu_high = schur_9to3(Yu9)
        Yd_high = schur_9to3(Yd9)
        Ye_high = schur_9to3(Ye9)

        # Cabibbo twist
        U_C = cabibbo_rotation(theta_C)
        Yd_high = U_C @ Yd_high

        # Neutrino projection + seesaw
        P = build_projection_resonant(lambda_nu)
        Ynu_eff_high = P @ Ynu9 @ P.conj().T
        Mnu_high = triadic_seesaw(M0, Ynu_eff_high, P)

        # Convert to Weinberg operator at high scale
        kappa_high = Mnu_high / (v_higgs**2)

        if USE_RGE:
            Yu_low, Yd_low, Ye_low, kappa_low = run_rge(
                MU_HIGH, MU_LOW,
                Yu_high, Yd_high, Ye_high, kappa_high,
                g1_EW, g2_EW, g3_EW, lambda_H
            )
            Mnu_low = kappa_low * (v_higgs**2)
        else:
            Yu_low, Yd_low, Ye_low = Yu_high, Yd_high, Ye_high
            Mnu_low = Mnu_high

        # Observables and chi²
        obs, Vckm, U_pmns = compute_observables(Yu_low, Yd_low, Ye_low, Mnu_low)
        chi2, pulls = chi2_from_obs(obs)

        # Regularization penalties
        shift_penalty = (
            np.sum(shifts_u**2) +
            np.sum(shifts_d**2) +
            np.sum(shifts_nu**2) +
            np.sum(shifts_e**2)
        )

        reg = (
            0.2  * shift_penalty +
            0.05 * (lambda_nu**2) +
            1.0  * (theta_C**2) +
            0.1  * (gamma_l**2)
        )

        return chi2 + reg

    except Exception:
        return 1e9


# =========================
# Optimizer
# =========================

def optimize_resonant16C(num_restarts=4, seed=9):
    rng = np.random.default_rng(seed)
    M0 = proto_majorana(rng)

    best_cost = np.inf
    best_X = None

    for r in range(num_restarts):
        print(f"\n=== Resonant-16C Restart {r+1}/{num_restarts} ===")

        # 17 parameters
        X0 = np.zeros(17)

        # phases
        X0[0:6] = rng.uniform(-0.5, 0.5, size=6)

        # exponent shifts
        X0[6:14] = rng.normal(scale=0.4, size=8)

        # lambda_nu, theta_C
        X0[14] = 0.15
        X0[15] = 0.10

        # gamma_l initial
        X0[16] = rng.uniform(-0.2, 0.2)

        es = cma.CMAEvolutionStrategy(
            X0,
            0.3,
            {
                'popsize': 20,
                'maxiter': 600,
                'CMA_diagonal': False,
                'seed': int(rng.integers(1, 2**31 - 1)),
            }
        )

        while not es.stop():
            xs = es.ask()
            cs = [resonant16C_cost(x, M0) for x in xs]
            es.tell(xs, cs)
            es.disp()

        if es.best.f < best_cost:
            best_cost = es.best.f
            best_X = es.best.x.copy()

    print("\nBEST Resonant-16C FIT:")
    print(best_X)
    print("cost =", best_cost)

    # Diagnostics at best fit
    (
        A_u, B_u,
        A_d, B_d,
        A_nu, B_nu,
        shifts_u, shifts_d, shifts_nu, shifts_e,
        lambda_nu, theta_C, gamma_l
    ) = unpack_params_16C(best_X)

    print("\nUnpacked parameters:")
    print("A_u, B_u    =", A_u, B_u)
    print("A_d, B_d    =", A_d, B_d)
    print("A_nu, B_nu  =", A_nu, B_nu)
    print("shifts_u    =", shifts_u)
    print("shifts_d    =", shifts_d)
    print("shifts_nu   =", shifts_nu)
    print("shifts_e    =", shifts_e)
    print("lambda_nu   =", lambda_nu)
    print("theta_C     =", theta_C)
    print("gamma_l     =", gamma_l)

    # Reconstruct best-fit low-scale matrices exactly as in cost
    gamma_q = 0.0

    Yu9  = build_Y_sector(A_u,  B_u,  EXP_U_BASE,  shifts_u,  gamma_q)
    Yd9  = build_Y_sector(A_d,  B_d,  EXP_D_BASE,  shifts_d,  gamma_q)
    Ye9  = build_Y_sector(A_d,  B_d,  EXP_E_BASE,  shifts_e,  gamma_l)
    Ynu9 = build_Y_sector(A_nu, B_nu, EXP_NU_BASE, shifts_nu, gamma_l)

    Yu_high = schur_9to3(Yu9)
    Yd_high = schur_9to3(Yd9)
    Ye_high = schur_9to3(Ye9)

    U_C = cabibbo_rotation(theta_C)
    Yd_high = U_C @ Yd_high

    P = build_projection_resonant(lambda_nu)
    Ynu_eff_high = P @ Ynu9 @ P.conj().T
    Mnu_high = triadic_seesaw(M0, Ynu_eff_high, P)
    kappa_high = Mnu_high / (v_higgs**2)

    if USE_RGE:
        Yu_low, Yd_low, Ye_low, kappa_low = run_rge(
            MU_HIGH, MU_LOW,
            Yu_high, Yd_high, Ye_high, kappa_high,
            g1_EW, g2_EW, g3_EW, lambda_H
        )
        Mnu_low = kappa_low * (v_higgs**2)
    else:
        Yu_low, Yd_low, Ye_low = Yu_high, Yd_high, Ye_high
        Mnu_low = Mnu_high

    obs, Vckm, U_pmns = compute_observables(Yu_low, Yd_low, Ye_low, Mnu_low)
    chi2, pulls = chi2_from_obs(obs)

    print(f"\nFinal evaluation at best fit:")
    print(f"  χ² = {chi2:.3f}")
    print(f"  total cost (χ² + reg) = {best_cost:.3f}\n")

    print("Observables (model vs target, pull in σ):")
    for key in exp_targets.keys():
        print(f"  {key:12s}: model={obs[key]:.6g}, "
              f"target={exp_targets[key]:.6g}, pull={pulls[key]: .3f}")

    print("\n|V_CKM| ≈")
    print(np.abs(Vckm))

    print("\n|U_PMNS| ≈")
    print(np.abs(U_pmns))

    return best_X, best_cost, M0


if __name__ == "__main__":
    best_X, best_cost, M0 = optimize_resonant16C(num_restarts=4, seed=9)