import numpy as np
import math
import cma

# ==================================
# Experimental targets & sigmas
# ==================================

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
sigma_targets['m_c/m_t'] *= 0.7
sigma_targets['m_u/m_t'] *= 0.7
...
# ==================================
# Linear algebra helpers
# ==================================

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
    """
    Standard PDG-like extraction of three mixing angles from a unitary matrix
    (ignoring CP phases).
    """
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

def stub_rge_run(M: np.ndarray,
                 alpha: float = 0.1,
                 mu_high: float = 1e14,
                 mu_EW: float = 173.0) -> np.ndarray:
    """
    Placeholder RGE: a common multiplicative factor.
    Ratios and mixings are essentially unchanged.
    """
    factor = math.log(mu_EW / mu_high) * alpha
    return M * math.exp(factor)

# ==================================
# Observables and χ²
# ==================================

def compute_observables_from_matrices(Yu_EW, Yd_EW, Ye_EW, Mnu_EW):
    # Quark Dirac diagonalization
    Uu_L, Su, _ = diagonalize_dirac(Yu_EW)
    Ud_L, Sd, _ = diagonalize_dirac(Yd_EW)

    # Leptons
    Ue_L, Se, _ = diagonalize_dirac(Ye_EW)

    mu_vals = np.diag(Su)
    md_vals = np.diag(Sd)
    me_vals = np.diag(Se)

    mu_sorted = np.sort(np.abs(mu_vals))[::-1]
    md_sorted = np.sort(np.abs(md_vals))[::-1]
    me_sorted = np.sort(np.abs(me_vals))[::-1]

    # CKM
    Vckm = Uu_L.conj().T @ Ud_L
    th12_q, th23_q, th13_q = extract_angles_from_U(Vckm)

    # Neutrinos: Majorana diagonalization
    U_nu, mnu_eig = diagonalize_majorana(Mnu_EW)
    mnu_sorted = np.sort(np.abs(mnu_eig))[::-1]

    # PMNS
    U_pmns = Ue_L.conj().T @ U_nu
    th12_l, th23_l, th13_l = extract_angles_from_U(U_pmns)

    # Fix absolute scales: normalize to physical top, bottom, tau, heaviest ν
    if mu_sorted[0] != 0:
        mu_sorted *= 173.0 / mu_sorted[0]
    if md_sorted[0] != 0:
        md_sorted *= 4.18 / md_sorted[0]
    if me_sorted[0] != 0:
        me_sorted *= 1.77686 / me_sorted[0]
    if mnu_sorted[0] != 0:
        mnu_sorted *= 0.058 / mnu_sorted[0]

    obs = {}
    # Mass ratios
    obs['m_c/m_t']      = mu_sorted[1] / mu_sorted[0] if mu_sorted[0] != 0 else 0.0
    obs['m_u/m_t']      = mu_sorted[2] / mu_sorted[0] if mu_sorted[0] != 0 else 0.0
    obs['m_s/m_b']      = md_sorted[1] / md_sorted[0] if md_sorted[0] != 0 else 0.0
    obs['m_d/m_b']      = md_sorted[2] / md_sorted[0] if md_sorted[0] != 0 else 0.0
    obs['m_mu/m_tau']   = me_sorted[1] / me_sorted[0] if me_sorted[0] != 0 else 0.0
    obs['m_e/m_tau']    = me_sorted[2] / me_sorted[0] if me_sorted[0] != 0 else 0.0

    # Mixings
    obs['theta12_q']    = th12_q
    obs['theta23_q']    = th23_q
    obs['theta13_q']    = th13_q
    obs['theta12_l']    = th12_l
    obs['theta23_l']    = th23_l
    obs['theta13_l']    = th13_l

    # Neutrino Δm²
    mnu_asc = np.sort(mnu_sorted)
    dm21 = mnu_asc[1]**2 - mnu_asc[0]**2
    dm31 = mnu_asc[2]**2 - mnu_asc[0]**2
    obs['Delta m2_21'] = dm21
    obs['Delta m2_31'] = dm31

    return obs, Vckm, U_pmns

# def chi2_from_obs(obs):
#     chi2 = 0.0
#     pulls = {}
#     for key, xexp in exp_targets.items():
#         xth = obs.get(key, np.nan)
#         if not np.isfinite(xth):
#             continue
#         sig = sigma_targets[key]
#         pull = (xth - xexp) / sig
#         chi2 += pull**2
#         pulls[key] = pull
#     return chi2, pulls

# ==================================
# Resonant-14C geometry
# ==================================

# Baseline "FN-like" charges per sector (left & right)
Q_L_UP_BASE   = np.array([4.0, 2.0, 0.0])
Q_R_UP_BASE   = np.array([4.0, 2.0, 0.0])

Q_L_DOWN_BASE = np.array([3.0, 2.0, 0.0])
Q_R_DOWN_BASE = np.array([3.0, 2.0, 0.0])

Q_L_LEP_BASE  = np.array([3.0, 2.0, 0.0])
Q_R_LEP_BASE  = np.array([3.0, 2.0, 0.0])

Q_L_NU_BASE   = np.array([1.0, 0.0, 0.0])
Q_R_NU_BASE   = np.array([1.0, 0.0, 0.0])

def build_yukawa_FNar_like(A, B, shifts, QL_base, QR_base, eps=0.23, kappa=0.1):
    d1, d2 = shifts
    QL_eff = np.array(QL_base, dtype=float)
    QL_eff[0] += d1
    QL_eff[1] += d2

    power = A * QL_eff[:, None] + B * QR_base[None, :]
    Y0 = eps ** power  # rank-1-ish

    # Minimal geometric deformation: small 2–3 “resonant” twist on the left
    D = np.diag([1.0, 1.0 + kappa, 1.0 + 2.0 * kappa])
    Y = D @ Y0
    return Y

def build_Cabibbo_rotation(theta_C: float) -> np.ndarray:
    """Left-handed 1–2 rotation for the down sector."""
    c = math.cos(theta_C)
    s = math.sin(theta_C)
    U = np.eye(3, dtype=complex)
    U[0, 0] = c
    U[0, 1] = s
    U[1, 0] = -s
    U[1, 1] = c
    return U

def build_majorana_resonant(Ynu: np.ndarray, lambda_nu: float) -> np.ndarray:
    """
    Simple resonant Majorana ansatz:

      M_R^{-1} = diag(1, 1 + λν, 1 + 2 λν)

    and m_ν ~ - Yν M_R^{-1} Yν^T (overall scale fixed later).
    """
    # Avoid singular behaviour for huge negative lambda
    lam = float(lambda_nu)
    inv_MR_diag = np.array([
        1.0,
        1.0 + lam,
        1.0 + 2.0 * lam
    ], dtype=float)

    # Make sure denominators are positive-ish
    inv_MR_diag = np.where(inv_MR_diag <= 1e-3, 1e-3, inv_MR_diag)

    inv_MR = np.diag(inv_MR_diag)
    # Up to an overall scale; we fix later via neutrino spectrum normalization
    mnu = - Ynu @ inv_MR @ Ynu.T
    # Add tiny regulator on the diagonal
    mnu = mnu + 1e-9 * np.eye(3)
    return mnu

# ==================================
# Parameter packing / unpacking
# ==================================

def unpack_params_resonant14C(X):
    """
    X has length 14:
      [A_u, B_u,
       A_d, B_d,
       A_nu, B_nu,
       su1, su2,
       sd1, sd2,      <-- ALSO used for charged leptons
       sn1, sn2,
       lambda_nu,
       theta_C]
    """
    X = np.array(X, dtype=float)

    A_u, B_u = X[0], X[1]
    A_d, B_d = X[2], X[3]
    A_nu, B_nu = X[4], X[5]

    su1, su2 = X[6], X[7]
    sd1, sd2 = X[8], X[9]
    sn1, sn2 = X[10], X[11]

    # lambda_nu = X[12]
    raw_lambda_nu = X[12]
    lambda_nu = raw_lambda_nu ** 2  # always ≥ 0
    theta_C   = X[13]

    shifts_u  = (su1, su2)
    shifts_d  = (sd1, sd2)
    shifts_e  = shifts_d              # locked: parameter reduction
    shifts_nu = (sn1, sn2)

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

# ==================================
# Cost function (Resonant-14C)
# ==================================

def evaluate_resonant14C_config(
        X,
        lambda_reg_exp=0.05,    # was 0.01
        lambda_reg_shift=0.1    # was 0.05
):
    """
    Compute cost = χ² + regularization for the 14-parameter resonant model.
    """
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
    ) = unpack_params_resonant14C(X)

    # ---------------------------------------
    # 1. Build Yukawas
    # ---------------------------------------
    Yu = build_yukawa_FNar_like(
        A_u, B_u, shifts_u,
        Q_L_UP_BASE, Q_R_UP_BASE
    )

    Yd = build_yukawa_FNar_like(
        A_d, B_d, shifts_d,
        Q_L_DOWN_BASE, Q_R_DOWN_BASE
    )

    # Charged leptons share down-sector shifts (locked)
    Ye = build_yukawa_FNar_like(
        A_d, B_d, shifts_e,
        Q_L_LEP_BASE, Q_R_LEP_BASE
    )

    Ynu = build_yukawa_FNar_like(
        A_nu, B_nu, shifts_nu,
        Q_L_NU_BASE, Q_R_NU_BASE
    )

    # Apply Cabibbo twist in left-handed down sector
    if abs(theta_C) > 0.0:
        U_C = build_Cabibbo_rotation(theta_C)
        Yd = U_C @ Yd

    # ---------------------------------------
    # 2. Majorana + seesaw-like structure
    # ---------------------------------------
    Mnu = build_majorana_resonant(Ynu, lambda_nu)

    # ---------------------------------------
    # 3. RG evolution
    # ---------------------------------------
    Yu_EW  = stub_rge_run(Yu)
    Yd_EW  = stub_rge_run(Yd)
    Ye_EW  = stub_rge_run(Ye)
    Mnu_EW = stub_rge_run(Mnu)

    # ---------------------------------------
    # 4. Observables & χ²
    # ---------------------------------------
    obs, Vckm, Upmns = compute_observables_from_matrices(
        Yu_EW, Yd_EW, Ye_EW, Mnu_EW
    )
    chi2, pulls = chi2_from_obs(obs)

    # ---------------------------------------
    # 5. Regularization (keep parameters in a natural-ish window)
    # ---------------------------------------
    exp_norm = (
        A_u**2 + B_u**2 +
        A_d**2 + B_d**2 +
        A_nu**2 + B_nu**2 +
        lambda_nu**2 + theta_C**2
    )

    shift_norm = (
        shifts_u[0]**2 + shifts_u[1]**2 +
        shifts_d[0]**2 + shifts_d[1]**2 +
        shifts_nu[0]**2 + shifts_nu[1]**2
    )

    cost = chi2 + lambda_reg_exp * exp_norm + lambda_reg_shift * shift_norm

    return cost, chi2, exp_norm, shift_norm, obs, pulls, Vckm, Upmns

def params_cost_resonant14C(X):
    try:
        cost, chi2, _, _, _, _, _, _ = evaluate_resonant14C_config(X)
    except Exception:
        return 1e9
    if not math.isfinite(cost):
        return 1e9
    return cost

# ==================================
# CMA-ES optimizer wrapper
# ==================================

def optimize_resonant14C(
        num_restarts=6,
        sigma_init=0.3,
        seed=42
):
    best_cost = np.inf
    best_X = None

    for r in range(num_restarts):
        restart_seed = seed + r
        print(f"\n=== Resonant-14C Restart {r+1}/{num_restarts} (seed={restart_seed}) ===")

        X0 = np.zeros(14)
        # Rough spreads for the initial cloud
        X0[0:6] = np.random.default_rng(restart_seed).normal(scale=0.5, size=6)  # A_* and B_*
        X0[6:12] = np.random.default_rng(restart_seed+1).normal(scale=0.5, size=6)  # shifts
        X0[12] = 0.0   # lambda_nu
        X0[13] = 0.2   # theta_C initial guess

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
            costs = [params_cost_resonant14C(x) for x in solutions]
            es.tell(solutions, costs)
            es.disp()

        Xbest = es.best.x
        costbest = es.best.f

        print(f"Restart {r+1}: best cost = {costbest:.6g}")

        if costbest < best_cost:
            best_cost = costbest
            best_X = Xbest.copy()

    return best_X, best_cost

def chi2_from_obs(obs):
    chi2 = 0.0
    pulls = {}
    for key, xexp in exp_targets.items():
        xth = obs.get(key, np.nan)
        if not np.isfinite(xth):
            continue
        sig = sigma_targets[key]
        pull = (xth - xexp) / sig
        # pull = max(min(pull, 10.0), -10.0)
        pull = max(min(pull, 20.0), -20.0)
        # IMPORTANT: no clamp here
        chi2 += pull**2
        pulls[key] = pull
    return chi2, pulls

# ==================================
# Main
# ==================================

if __name__ == "__main__":
    best_X, best_cost = optimize_resonant14C(
        num_restarts=4,
        sigma_init=0.3,
        seed=11
    )

    cost, chi2, exp_norm, shift_norm, obs, pulls, Vckm, Upmns = \
        evaluate_resonant14C_config(best_X)

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
    ) = unpack_params_resonant14C(best_X)

    print("\n=== BEST Resonant-14C FIT ===")
    print("Parameters X* =", best_X)
    print("\nUnpacked:")
    print(f"A_u, B_u    = {A_u} {B_u}")
    print(f"A_d, B_d    = {A_d} {B_d}")
    print(f"A_nu, B_nu  = {A_nu} {B_nu}")
    print(f"shifts_u    = {shifts_u}")
    print(f"shifts_d    = {shifts_d}")
    print(f"shifts_e    = {shifts_e}  (locked to shifts_d)")
    print(f"shifts_nu   = {shifts_nu}")
    print(f"lambda_nu   = {lambda_nu}")
    print(f"theta_C     = {theta_C}")

    print("\nFit quality:")
    print(f"  χ²         = {chi2:.3f}")
    print(f"  exp_norm   = {exp_norm:.3f}")
    print(f"  shift_norm = {shift_norm:.3f}")
    print(f"  total cost = {cost:.3f}\n")

    print("Observables (model vs target, pull in σ):")
    for key in exp_targets.keys():
        print(f"  {key:12s}  model = {obs[key]:.6g}   "
              f"target = {exp_targets[key]:.6g}   pull = {pulls[key]: .3f}")

    print("\n|V_CKM| ≈")
    print(np.abs(Vckm))

    print("\n|U_PMNS| ≈")
    print(np.abs(Upmns))