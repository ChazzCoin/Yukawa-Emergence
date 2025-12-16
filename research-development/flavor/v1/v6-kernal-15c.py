#!/usr/bin/env python3
import numpy as np
import math
import cma

# ===========================================================
#   Experimental targets & sigmas
# ===========================================================

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
    'm_c/m_t': 0.35 * exp_targets['m_c/m_t'],
    'm_u/m_t': 0.35 * exp_targets['m_u/m_t'],
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

# ===========================================================
#   Linear algebra helpers
# ===========================================================

def diagonalize_dirac(Y):
    U_L, S, U_Rh = np.linalg.svd(Y)
    return U_L, np.diag(S), U_Rh.conj().T

def diagonalize_majorana(M):
    H = 0.5 * (M + M.conj().T)
    vals, U = np.linalg.eigh(H)
    idx = np.argsort(np.abs(vals))[::-1]
    return U[:, idx], vals[idx]

def extract_angles_from_U(U):
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

def stub_rge_run(M,
                 alpha=0.1,
                 mu_high=1e14,
                 mu_EW=173.0):
    factor = math.log(mu_EW / mu_high) * alpha
    return M * math.exp(factor)

# ===========================================================
#   FN-like Yukawa builder
# ===========================================================

Q_L_UP_BASE   = np.array([4.0, 2.0, 0.0])
Q_R_UP_BASE   = np.array([4.0, 2.0, 0.0])

Q_L_DOWN_BASE = np.array([3.0, 2.0, 0.0])
Q_R_DOWN_BASE = np.array([3.0, 2.0, 0.0])

Q_L_LEP_BASE  = np.array([3.0, 2.0, 0.0])
Q_R_LEP_BASE  = np.array([3.0, 2.0, 0.0])

Q_L_NU_BASE   = np.array([1.0, 0.0, 0.0])
Q_R_NU_BASE   = np.array([1.0, 0.0, 0.0])

EPS_FN = 0.23

def build_yukawa_FNar_like(A, B, shifts, QL_base, QR_base, eps=EPS_FN):
    d1, d2 = shifts
    QL_eff = np.array(QL_base, dtype=float)
    QL_eff[0] += d1
    QL_eff[1] += d2
    power = A * QL_eff[:, None] + B * QR_base[None, :]
    return eps ** power

def build_Cabibbo_rotation(theta_C):
    c = math.cos(theta_C)
    s = math.sin(theta_C)
    U = np.eye(3, dtype=complex)
    U[0, 0] = c; U[0, 1] = s
    U[1, 0] = -s; U[1, 1] = c
    return U

# Fixed resonant parameter (from successful 16C-like geometry)
LAMBDA_NU0 = 0.06

def build_majorana_resonant(Ynu, lambda_nu=LAMBDA_NU0):
    lam = float(lambda_nu)
    inv_MR_diag = np.array([
        1.0,
        1.0 + lam,
        1.0 + 2.0 * lam
    ], dtype=float)
    inv_MR_diag = np.where(inv_MR_diag <= 1e-3, 1e-3, inv_MR_diag)
    mnu = - Ynu @ np.diag(inv_MR_diag) @ Ynu.T
    return mnu + 1e-9 * np.eye(3)

# ===========================================================
#   Observables and χ²
# ===========================================================

def compute_observables_from_matrices(Yu_EW, Yd_EW, Ye_EW, Mnu_EW):
    Uu_L, Su, _ = diagonalize_dirac(Yu_EW)
    Ud_L, Sd, _ = diagonalize_dirac(Yd_EW)
    Ue_L, Se, _ = diagonalize_dirac(Ye_EW)

    mu = np.sort(np.abs(np.diag(Su)))[::-1]
    md = np.sort(np.abs(np.diag(Sd)))[::-1]
    me = np.sort(np.abs(np.diag(Se)))[::-1]

    Vckm = Uu_L.conj().T @ Ud_L

    U_nu, mnu_vals = diagonalize_majorana(Mnu_EW)
    mnu = np.sort(np.abs(mnu_vals))[::-1]

    U_pmns = Ue_L.conj().T @ U_nu

    # Normalize absolute scales
    if mu[0] != 0:  mu *= 173.0 / mu[0]
    if md[0] != 0:  md *= 4.18 / md[0]
    if me[0] != 0:  me *= 1.77686 / me[0]
    if mnu[0] != 0: mnu *= 0.058 / mnu[0]

    th12_q, th23_q, th13_q = extract_angles_from_U(Vckm)
    th12_l, th23_l, th13_l = extract_angles_from_U(U_pmns)

    mnu_asc = np.sort(mnu)
    dm21 = mnu_asc[1]**2 - mnu_asc[0]**2
    dm31 = mnu_asc[2]**2 - mnu_asc[0]**2

    obs = {
        'm_c/m_t':     mu[1]/mu[0] if mu[0] != 0 else 0.0,
        'm_u/m_t':     mu[2]/mu[0] if mu[0] != 0 else 0.0,
        'm_s/m_b':     md[1]/md[0] if md[0] != 0 else 0.0,
        'm_d/m_b':     md[2]/md[0] if md[0] != 0 else 0.0,
        'm_mu/m_tau':  me[1]/me[0] if me[0] != 0 else 0.0,
        'm_e/m_tau':   me[2]/me[0] if me[0] != 0 else 0.0,
        'theta12_q':   th12_q,
        'theta23_q':   th23_q,
        'theta13_q':   th13_q,
        'theta12_l':   th12_l,
        'theta23_l':   th23_l,
        'theta13_l':   th13_l,
        'Delta m2_21': dm21,
        'Delta m2_31': dm31,
    }

    return obs, Vckm, U_pmns

def chi2_from_obs(obs):
    chi2 = 0.0
    pulls = {}
    for key, xexp in exp_targets.items():
        xth = obs[key]
        sig = sigma_targets[key]
        pull = (xth - xexp) / sig
        pull = max(min(pull, 20.0), -20.0)
        chi2 += pull**2
        pulls[key] = pull
    return chi2, pulls

# ===========================================================
#   Parameter unpacking: Resonant-15C′
# ===========================================================

def unpack_params_resonant15Cprime(X):
    """
    X has 15 entries:

      [A_u, B_u,
       A_d, B_d,
       A_nu, B_nu,
       su1, su2,
       sd1, sd2,
       se1, se2,
       sn1, sn2,
       theta_C]
    """
    X = np.array(X, dtype=float)

    A_u, B_u   = X[0], X[1]
    A_d, B_d   = X[2], X[3]
    A_nu, B_nu = X[4], X[5]

    su1, su2 = X[6], X[7]
    sd1, sd2 = X[8], X[9]
    se1, se2 = X[10], X[11]
    sn1, sn2 = X[12], X[13]

    theta_C = X[14]

    shifts_u  = (su1, su2)
    shifts_d  = (sd1, sd2)
    shifts_e  = (se1, se2)     # fully independent
    shifts_nu = (sn1, sn2)

    lambda_nu = LAMBDA_NU0     # fixed resonant parameter

    return (
        A_u, B_u,
        A_d, B_d,
        A_nu, B_nu,
        shifts_u,
        shifts_d,
        shifts_e,
        shifts_nu,
        lambda_nu,
        theta_C
    )

# ===========================================================
#   Cost function
# ===========================================================

def evaluate_resonant15Cprime_config(
        X,
        lambda_reg_exp=0.05,
        lambda_reg_shift=0.08
):
    (
        A_u, B_u,
        A_d, B_d,
        A_nu, B_nu,
        shifts_u,
        shifts_d,
        shifts_e,
        shifts_nu,
        lambda_nu,
        theta_C
    ) = unpack_params_resonant15Cprime(X)

    # Yukawas
    Yu = build_yukawa_FNar_like(
        A_u, B_u, shifts_u,
        Q_L_UP_BASE, Q_R_UP_BASE
    )
    Yd = build_yukawa_FNar_like(
        A_d, B_d, shifts_d,
        Q_L_DOWN_BASE, Q_R_DOWN_BASE
    )
    Ye = build_yukawa_FNar_like(
        A_d, B_d, shifts_e,
        Q_L_LEP_BASE, Q_R_LEP_BASE
    )
    Ynu = build_yukawa_FNar_like(
        A_nu, B_nu, shifts_nu,
        Q_L_NU_BASE, Q_R_NU_BASE
    )

    if abs(theta_C) > 1e-12:
        U_C = build_Cabibbo_rotation(theta_C)
        Yd = U_C @ Yd

    Mnu = build_majorana_resonant(Ynu, lambda_nu=lambda_nu)

    # RGE stub
    Yu_EW  = stub_rge_run(Yu)
    Yd_EW  = stub_rge_run(Yd)
    Ye_EW  = stub_rge_run(Ye)
    Mnu_EW = stub_rge_run(Mnu)

    obs, Vckm, Upmns = compute_observables_from_matrices(
        Yu_EW, Yd_EW, Ye_EW, Mnu_EW
    )
    chi2, pulls = chi2_from_obs(obs)

    exp_norm = (
        A_u**2 + B_u**2 +
        A_d**2 + B_d**2 +
        A_nu**2 + B_nu**2 +
        theta_C**2
        # lambda_nu is fixed; leaving it out of regularization
    )

    shift_norm = (
        shifts_u[0]**2 + shifts_u[1]**2 +
        shifts_d[0]**2 + shifts_d[1]**2 +
        shifts_e[0]**2 + shifts_e[1]**2 +
        shifts_nu[0]**2 + shifts_nu[1]**2
    )

    cost = chi2 + lambda_reg_exp * exp_norm + lambda_reg_shift * shift_norm

    return cost, chi2, exp_norm, shift_norm, obs, pulls, Vckm, Upmns

def params_cost_resonant15Cprime(X):
    try:
        cost, *_ = evaluate_resonant15Cprime_config(X)
        if not math.isfinite(cost):
            return 1e9
        return cost
    except Exception:
        return 1e9

# ===========================================================
#   CMA-ES optimizer wrapper
# ===========================================================

def optimize_resonant15Cprime(
        num_restarts=4,
        sigma_init=0.3,
        seed=42
):
    best_cost = np.inf
    best_X = None

    for r in range(num_restarts):
        restart_seed = seed + r
        print(f"\n=== Resonant-15C′ Restart {r+1}/{num_restarts} (seed={restart_seed}) ===")

        rngA = np.random.default_rng(restart_seed)
        rngS = np.random.default_rng(restart_seed + 1)

        X0 = np.zeros(15)
        X0[0:6]  = rngA.normal(scale=0.5, size=6)   # A_u..B_nu
        X0[6:14] = rngS.normal(scale=0.5, size=8)   # all shifts
        X0[14]   = 0.2                              # theta_C initial guess

        es = cma.CMAEvolutionStrategy(
            X0,
            sigma_init,
            {
                'popsize': 20,
                'maxiter': 600,
                'seed': restart_seed,
            }
        )

        while not es.stop():
            solutions = es.ask()
            costs = [params_cost_resonant15Cprime(x) for x in solutions]
            es.tell(solutions, costs)
            es.disp()

        Xbest = es.best.x
        costbest = es.best.f
        print(f"Restart {r+1}: best cost = {costbest:.6g}")

        if costbest < best_cost:
            best_cost = costbest
            best_X = Xbest.copy()

    return best_X, best_cost

# ===========================================================
#   Main
# ===========================================================

if __name__ == "__main__":
    best_X, best_cost = optimize_resonant15Cprime(
        num_restarts=4,
        sigma_init=0.3,
        seed=11
    )

    cost, chi2, exp_norm, shift_norm, obs, pulls, Vckm, Upmns = \
        evaluate_resonant15Cprime_config(best_X)

    (
        A_u, B_u,
        A_d, B_d,
        A_nu, B_nu,
        shifts_u,
        shifts_d,
        shifts_e,
        shifts_nu,
        lambda_nu,
        theta_C
    ) = unpack_params_resonant15Cprime(best_X)

    print("\n=== BEST Resonant-15C′ FIT ===")
    print("Parameters X* =", best_X)

    print("\nUnpacked:")
    print(f"A_u, B_u    = {A_u}, {B_u}")
    print(f"A_d, B_d    = {A_d}, {B_d}")
    print(f"A_nu, B_nu  = {A_nu}, {B_nu}")
    print(f"shifts_u    = {shifts_u}")
    print(f"shifts_d    = {shifts_d}")
    print(f"shifts_e    = {shifts_e}  (fully independent)")
    print(f"shifts_nu   = {shifts_nu}")
    print(f"lambda_nu   = {lambda_nu}  (fixed)")
    print(f"theta_C     = {theta_C}")

    print("\nFit quality:")
    print(f"  χ²         = {chi2:.3f}")
    print(f"  exp_norm   = {exp_norm:.3f}")
    print(f"  shift_norm = {shift_norm:.3f}")
    print(f"  total cost = {cost:.3f}\n")

    print("Observables (model vs target, pull):")
    for key in exp_targets.keys():
        print(f"  {key:12s}  model={obs[key]:.6g}  "
              f"target={exp_targets[key]:.6g}  pull={pulls[key]: .3f}")

    print("\n|V_CKM| ≈")
    print(np.abs(Vckm))

    print("\n|U_PMNS| ≈")
    print(np.abs(Upmns))