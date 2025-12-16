import numpy as np
import math
import cma

# ================================================================
#  RESONANT-16C FLAVOR MODEL
#  - 16 parameters
#  - Quarks and leptons maximally coherent (gamma_q = gamma_l = 0)
#  - Explicit Cabibbo 1–2 left-handed twist in the down sector
#  - Independent exponent shaping for u, d, e, nu
# ================================================================

N_SITES = 9
LIGHT_SITES = [0, 1, 2]
HEAVY_SITES = [3, 4, 5, 6, 7, 8]
PI = math.pi

# ------------------------------------------------
# Basic helpers
# ------------------------------------------------

def clamp(x, lo, hi):
    return np.minimum(np.maximum(x, lo), hi)

def generation_index(i: int) -> int:
    return i % 3

# ------------------------------------------------
# Triadic base exponents
# ------------------------------------------------

EXP_U_BASE  = np.array([4.0, 2.0, 0.0])
EXP_D_BASE  = np.array([3.0, 2.0, 0.0])
EXP_E_BASE  = np.array([3.0, 2.0, 0.0])  # same pattern, different shifts
EXP_NU_BASE = np.array([1.0, 0.0, 0.0])

# ================================================================
#  RESONANT PHASE WHEELS
# ================================================================

def build_phase_profile(A: float, B: float):
    """
    φ_g = A + B*g, g=0,1,2
    """
    return np.array([A, A + B, A + 2*B], dtype=float)

def build_site_phases(phi_gen):
    phi_site = np.zeros(N_SITES, dtype=float)
    for i in range(N_SITES):
        phi_site[i] = phi_gen[generation_index(i)]
    return phi_site

def build_phase_matrix(phi_site):
    P = np.zeros((N_SITES, N_SITES), dtype=complex)
    for i in range(N_SITES):
        for j in range(N_SITES):
            P[i, j] = np.exp(1j * (phi_site[i] - phi_site[j]))
    return P

# ================================================================
#  EXPONENT HIERARCHIES (X,Y shifts)
# ================================================================

def site_scales(base_exp, shifts):
    """
    base_exp: [a,b,c]
    shifts: (X,Y) applied to 2nd & 3rd gen exponents
      e0 = base[0]
      e1 = base[1] + X
      e2 = base[2] + Y

    s_g = KAPPA^e_g, with KAPPA ≈ 0.24
    """
    X, Y = shifts
    eff = np.array([
        base_exp[0],
        base_exp[1] + X,
        base_exp[2] + Y
    ], dtype=float)

    KAPPA = 0.24
    s_gen = np.power(KAPPA, eff)

    s = np.zeros(N_SITES, dtype=float)
    for i in range(N_SITES):
        s[i] = s_gen[generation_index(i)]
    return s

# ================================================================
#  COHERENCE KERNEL (γ) — HERE γ = 0 (COHERENT LIMIT)
# ================================================================

def build_kernel_gamma(gamma: float, forbidden_d: int = 2):
    """
    Toeplitz kernel on the ring, but we'll use gamma=0
    so K becomes almost trivial (except for forbidden distance, if kept).
    """
    K = np.zeros((N_SITES, N_SITES), dtype=float)
    for i in range(N_SITES):
        for j in range(N_SITES):
            if i == j:
                K[i, j] = 1.0
            else:
                d = min(abs(i - j), N_SITES - abs(i - j))
                if d == forbidden_d:
                    K[i, j] = 0.0
                else:
                    K[i, j] = math.exp(-gamma * d)
    return K

# ================================================================
#  NEUTRINO PROJECTION RESONANCE (λ_ν)
# ================================================================

def triadic_modes():
    n = 9
    j = np.arange(n)
    v0 = np.exp(2j*np.pi*0*j/n) / np.sqrt(n)
    v3 = np.exp(2j*np.pi*3*j/n) / np.sqrt(n)
    v6 = np.exp(2j*np.pi*6*j/n) / np.sqrt(n)
    return v0, v3, v6

def build_projection_resonant(lambda_nu: float):
    """
    λ_ν in [0, 0.3]: small mixing between triadic modes before QR.
    """
    v0, v3, v6 = triadic_modes()

    b0 = v0 + lambda_nu * v3
    b3 = v3 + lambda_nu * v6
    b6 = v6 + lambda_nu * v3

    B = np.vstack([b0, b3, b6])
    Q, _ = np.linalg.qr(B.conj().T)   # 9×3
    return Q.T                        # 3×9

# ================================================================
#  PROTO-MAJORANA (fixed per run)
# ================================================================

def proto_majorana(rng: np.random.Generator):
    M = rng.normal(size=(N_SITES, N_SITES)) + 1j * rng.normal(size=(N_SITES, N_SITES))
    M = 0.5 * (M + M.T)
    _, s, _ = np.linalg.svd(M)
    M /= s[0]
    return M

# ================================================================
#  9→3 SCHUR REDUCTION
# ================================================================

def schur_9to3(Y9: np.ndarray) -> np.ndarray:
    ls = LIGHT_SITES
    hs = HEAVY_SITES
    A = Y9[np.ix_(ls, ls)]
    B = Y9[np.ix_(ls, hs)]
    C = Y9[np.ix_(hs, ls)]
    D = Y9[np.ix_(hs, hs)]
    Dinv = np.linalg.pinv(D)
    Y_eff = A - B @ Dinv @ C
    return Y_eff + 1e-9 * np.eye(3)

# ================================================================
#  SECTOR YUKAWAS
# ================================================================

def build_Y_sector(A, B, base_exp, shifts, gamma):
    phi_gen = build_phase_profile(A, B)
    phi_site = build_site_phases(phi_gen)
    P = build_phase_matrix(phi_site)

    s = site_scales(base_exp, shifts)
    mag = np.outer(s, s)

    Y0 = mag * P
    K = build_kernel_gamma(gamma)
    Y = K * Y0

    # SVD normalization
    _, sv, _ = np.linalg.svd(Y)
    if sv[0] != 0:
        Y /= sv[0]
    return Y

# ================================================================
#  TRIADIC MAJORANA SEESAW
# ================================================================

def triadic_seesaw(M9: np.ndarray, Ynu_eff: np.ndarray) -> np.ndarray:
    M_H = M9[np.ix_(HEAVY_SITES, HEAVY_SITES)]
    h = np.arange(len(HEAVY_SITES))

    B = np.zeros((len(HEAVY_SITES), 3), dtype=complex)
    for col, k in enumerate([1, 2, 3]):
        B[:, col] = np.exp(2j * np.pi * k * h / len(HEAVY_SITES)) / math.sqrt(len(HEAVY_SITES))

    M_R = B.conj().T @ M_H @ B
    M_R = 0.5 * (M_R + M_R.T)
    M_R += 1e-9 * np.eye(3)
    M_R *= 7e13

    v = 174.0 / math.sqrt(2.0)
    mD = v * Ynu_eff
    M_Rinv = np.linalg.inv(M_R)
    return - mD @ M_Rinv @ mD.T

# ================================================================
#  OBSERVABLES
# ================================================================

def diagonalize_dirac(Y: np.ndarray):
    UL, S, URh = np.linalg.svd(Y)
    return UL, np.diag(S), URh.conj().T

def diag_majorana(M: np.ndarray):
    H = 0.5 * (M + M.conj().T)
    vals, U = np.linalg.eigh(H)
    idx = np.argsort(np.abs(vals))[::-1]
    return U[:, idx], vals[idx]

def extract_angles(U: np.ndarray):
    U = np.array(U, dtype=complex)
    s13 = abs(U[0, 2])
    s13 = min(max(s13, 0.0), 1.0)
    th13 = math.asin(s13)
    c13 = math.cos(th13)
    if c13 == 0:
        return 0.0, 0.0, th13
    s12 = abs(U[0, 1]) / c13
    s23 = abs(U[1, 2]) / c13
    s12 = min(max(s12, 0.0), 1.0)
    s23 = min(max(s23, 0.0), 1.0)
    return math.asin(s12), math.asin(s23), th13

exp_targets = {
    'm_c/m_t':      7e-3,
    'm_u/m_t':      1e-5,
    'm_s/m_b':      2e-2,
    'm_d/m_b':      1e-3,
    'm_mu/m_tau':   6e-2,
    'm_e/m_tau':    3e-4,
    'theta12_q':    0.226,
    'theta23_q':    0.041,
    'theta13_q':    0.0035,
    'theta12_l':    0.59,
    'theta23_l':    0.84,
    'theta13_l':    0.15,
    'Delta m2_21':  7.4e-5,
    'Delta m2_31':  2.5e-3,
}

sigma_targets = {k: 0.3 * v for k, v in exp_targets.items()}

def compute_observables(Yu, Yd, Ye, Mnu):
    Uu, Su, _ = diagonalize_dirac(Yu)
    Ud, Sd, _ = diagonalize_dirac(Yd)
    Ue, Se, _ = diagonalize_dirac(Ye)

    mu = np.sort(np.abs(np.diag(Su)))[::-1]
    md = np.sort(np.abs(np.diag(Sd)))[::-1]
    me = np.sort(np.abs(np.diag(Se)))[::-1]

    Vckm = Uu.conj().T @ Ud
    th12_q, th23_q, th13_q = extract_angles(Vckm)

    U_nu, mnu_vals = diag_majorana(Mnu)
    mnu = np.sort(np.abs(mnu_vals))[::-1]

    U_pmns = Ue.conj().T @ U_nu
    th12_l, th23_l, th13_l = extract_angles(U_pmns)

    # Rescale masses to physical scales
    mu *= 173.0   / mu[0]
    md *= 4.18    / md[0]
    me *= 1.77686 / me[0]
    mnu *= 0.058  / mnu[0]

    mnu_asc = np.sort(mnu)
    dm21 = mnu_asc[1]**2 - mnu_asc[0]**2
    dm31 = mnu_asc[2]**2 - mnu_asc[0]**2

    obs = {
        'm_c/m_t':      mu[1] / mu[0],
        'm_u/m_t':      mu[2] / mu[0],
        'm_s/m_b':      md[1] / md[0],
        'm_d/m_b':      md[2] / md[0],
        'm_mu/m_tau':   me[1] / me[0],
        'm_e/m_tau':    me[2] / me[0],
        'theta12_q':    th12_q,
        'theta23_q':    th23_q,
        'theta13_q':    th13_q,
        'theta12_l':    th12_l,
        'theta23_l':    th23_l,
        'theta13_l':    th13_l,
        'Delta m2_21':  dm21,
        'Delta m2_31':  dm31,
    }
    return obs, Vckm, U_pmns

def chi2_from_obs(obs):
    chi2 = 0.0
    pulls = {}
    for key, target in exp_targets.items():
        th = obs[key]
        sig = sigma_targets[key]
        pull = (th - target) / sig
        chi2 += pull**2
        pulls[key] = pull
    return chi2, pulls

# ================================================================
#  CABIBBO TWIST (DOWN-SECTOR LEFT ROTATION)
# ================================================================

def cabibbo_rotation(theta_C: float) -> np.ndarray:
    """
    Real 1–2 rotation acting on left-handed down states.
    """
    c = math.cos(theta_C)
    s = math.sin(theta_C)
    U = np.eye(3, dtype=complex)
    U[0, 0] = c
    U[0, 1] = s
    U[1, 0] = -s
    U[1, 1] = c
    return U

# ================================================================
#  PARAMETER UNPACKING (16 params)
# ================================================================

def unpack_params_16C(X):
    """
    X: array-like of length 17
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

# ================================================================
#  COST FUNCTION
# ================================================================

def resonant16C_cost(X, M0):
    """
    Cost function with:
      - free lepton coherence gamma_l (quarks fixed at gamma_q = 0),
      - softened penalty on lambda_nu so neutrino resonance can work.
    """
    try:
        (A_u, B_u,
         A_d, B_d,
         A_nu, B_nu,
         shifts_u, shifts_d, shifts_nu, shifts_e,
         lambda_nu, theta_C, gamma_l) = unpack_params_16C(X)

        # Quark kernel: keep fully coherent (except forbidden distance)
        gamma_q = 0.0

        # 9×9 Yukawas
        Yu9  = build_Y_sector(A_u,  B_u,  EXP_U_BASE,  shifts_u,  gamma_q)
        Yd9  = build_Y_sector(A_d,  B_d,  EXP_D_BASE,  shifts_d,  gamma_q)

        # Leptons feel their own coherence length gamma_l
        Ye9  = build_Y_sector(A_d,  B_d,  EXP_E_BASE,  shifts_e,   gamma_l)
        Ynu9 = build_Y_sector(A_nu, B_nu, EXP_NU_BASE, shifts_nu,  gamma_l)

        # 9→3 Schur reduction
        Yu = schur_9to3(Yu9)
        Yd = schur_9to3(Yd9)
        Ye = schur_9to3(Ye9)

        # Apply Cabibbo twist to down sector
        U_C = cabibbo_rotation(theta_C)
        Yd = U_C @ Yd

        # Neutrino projection + seesaw
        P        = build_projection_resonant(lambda_nu)
        Ynu_eff  = P @ Ynu9 @ P.conj().T
        Mnu      = triadic_seesaw(M0, Ynu_eff)

        # Observables
        obs, Vckm, U_pmns = compute_observables(Yu, Yd, Ye, Mnu)
        chi2, pulls       = chi2_from_obs(obs)

        # Soft penalties: keep shifts modest, allow lambda_nu and gamma_l to work but not blow up
        shift_penalty = (
            np.sum(np.array(shifts_u)**2) +
            np.sum(np.array(shifts_d)**2) +
            np.sum(np.array(shifts_nu)**2) +
            np.sum(np.array(shifts_e)**2)
        )

        reg = (
            0.2 * shift_penalty      # keep exponent shifts O(1)
            + 0.05 * (lambda_nu**2)  # much softer neutrino resonance penalty
            + 1.0 * (theta_C**2)     # Cabibbo stays modest
            + 0.1 * (gamma_l**2)     # lepton coherence not too extreme
        )

        return chi2 + reg

    except Exception:
        # Any numerical failure (e.g. bad Schur, singular seesaw) gets a huge penalty
        return 1e9

# ================================================================
#  OPTIMIZER DRIVER
# ================================================================

def optimize_resonant16C(num_restarts=4, seed=9):
    rng = np.random.default_rng(seed)
    M0 = proto_majorana(rng)

    best_cost = np.inf
    best_X = None

    for r in range(num_restarts):
        print(f"\n=== Resonant-16C Restart {r+1}/{num_restarts} ===")

        # 17 parameters now: 6 phases, 8 shifts, lambda_nu, theta_C, gamma_l
        X0 = np.zeros(17)

        # phases: A_u, B_u, A_d, B_d, A_nu, B_nu
        X0[0:6] = rng.uniform(-0.5, 0.5, size=6)

        # exponent shifts: (shift_u1, shift_u2, shift_d1, shift_d2,
        #                   shift_nu1, shift_nu2, shift_e1, shift_e2)
        X0[6:14] = rng.normal(scale=0.4, size=8)

        # lambda_nu, theta_C
        X0[14] = 0.15
        X0[15] = 0.1

        # initial lepton coherence gamma_l (small positive)
        X0[16] = rng.uniform(0.0, 0.2)

        es = cma.CMAEvolutionStrategy(
            X0,
            0.3,
            {
                'popsize': 20,
                'maxiter': 600,
                'CMA_diagonal': False,
                # Per-restart CMA seed, reproducible for fixed outer 'seed'
                'seed': int(rng.integers(1, 2 ** 31 - 1)),
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

    # Diagnostics at best fit: mirror the cost pipeline exactly
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

    # Quark coherence fixed, lepton coherence from fit
    gamma_q = 0.0

    # 9×9 Yukawas (same as in resonant16C_cost)
    Yu9  = build_Y_sector(A_u,  B_u,  EXP_U_BASE,  shifts_u,  gamma_q)
    Yd9  = build_Y_sector(A_d,  B_d,  EXP_D_BASE,  shifts_d,  gamma_q)
    Ye9  = build_Y_sector(A_d,  B_d,  EXP_E_BASE,  shifts_e,  gamma_l)
    Ynu9 = build_Y_sector(A_nu, B_nu, EXP_NU_BASE, shifts_nu, gamma_l)

    # 9→3 Schur
    Yu = schur_9to3(Yu9)
    Yd = schur_9to3(Yd9)
    Ye = schur_9to3(Ye9)

    # Cabibbo twist
    U_C = cabibbo_rotation(theta_C)
    Yd = U_C @ Yd

    # Neutrino projection + seesaw
    P       = build_projection_resonant(lambda_nu)
    Ynu_eff = P @ Ynu9 @ P.conj().T
    Mnu     = triadic_seesaw(M0, Ynu_eff)

    # Observables and chi² (without reg)
    obs, Vckm, U_pmns = compute_observables(Yu, Yd, Ye, Mnu)
    chi2, pulls       = chi2_from_obs(obs)

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
    optimize_resonant16C(num_restarts=4, seed=9)

"""
=== Resonant-16C Restart 1/4 ===
(10_w,20)-aCMA-ES (mu_w=5.9,w_1=27%) in dimension 17 (seed=1750994040, Mon Dec  8 22:17:22 2025)
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
    1     20 2.787083241766091e+02 1.0e+00 2.96e-01  3e-01  3e-01 0:00.0
    2     40 2.536720599946552e+02 1.2e+00 3.11e-01  3e-01  3e-01 0:00.0
    3     60 1.389187284007203e+02 1.3e+00 3.14e-01  3e-01  4e-01 0:00.0
  100   2000 6.555330932798485e+01 7.8e+00 8.34e-02  2e-02  1e-01 0:01.5
  200   4000 4.745717219840245e+01 3.6e+01 1.12e-01  2e-02  2e-01 0:02.9
  300   6000 3.558346539489390e+01 6.1e+01 4.85e-02  6e-03  8e-02 0:04.3
  400   8000 1.055058614334735e+01 8.3e+01 1.36e-01  1e-02  2e-01 0:05.7
  500  10000 9.173887312189702e+00 1.9e+02 2.95e-02  9e-04  6e-02 0:07.1
  600  12000 9.157897982590979e+00 4.5e+02 1.46e-03  2e-05  3e-03 0:08.5

=== Resonant-16C Restart 2/4 ===
(10_w,20)-aCMA-ES (mu_w=5.9,w_1=27%) in dimension 17 (seed=1096649525, Mon Dec  8 22:17:30 2025)
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
    1     20 5.873021166107926e+02 1.0e+00 2.75e-01  3e-01  3e-01 0:00.0
    2     40 3.468896054595806e+02 1.1e+00 2.67e-01  3e-01  3e-01 0:00.0
    3     60 1.049018217679392e+02 1.2e+00 2.55e-01  2e-01  3e-01 0:00.0
  100   2000 2.495067450433601e+01 9.1e+00 6.53e-02  2e-02  9e-02 0:01.4
  200   4000 8.513983490042873e+00 3.0e+01 4.18e-02  4e-03  5e-02 0:02.8
  300   6000 8.130697097402338e+00 1.0e+02 9.13e-03  4e-04  2e-02 0:04.1
  400   8000 8.098387841597075e+00 1.6e+02 9.71e-03  3e-04  2e-02 0:05.5
  500  10000 8.073391425843983e+00 2.4e+02 1.19e-02  3e-04  2e-02 0:07.1
  600  12000 8.069307318150425e+00 3.6e+02 1.60e-04  2e-06  3e-04 0:08.6

=== Resonant-16C Restart 3/4 ===
(10_w,20)-aCMA-ES (mu_w=5.9,w_1=27%) in dimension 17 (seed=1407148388, Mon Dec  8 22:17:39 2025)
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
    1     20 5.949783507537333e+02 1.0e+00 2.85e-01  3e-01  3e-01 0:00.0
    2     40 9.545029472985387e+01 1.2e+00 2.70e-01  3e-01  3e-01 0:00.0
    3     60 1.189628043814013e+02 1.2e+00 2.65e-01  3e-01  3e-01 0:00.0
  100   2000 2.162972419784850e+01 1.3e+01 1.11e-01  2e-02  2e-01 0:01.9
  200   4000 8.246742473323266e+00 2.4e+01 4.47e-02  4e-03  7e-02 0:03.8
  300   6000 7.290721723576441e+00 6.3e+01 6.53e-03  5e-04  1e-02 0:05.8
  400   8000 7.279366448513139e+00 1.4e+02 9.84e-04  4e-05  2e-03 0:07.8
  500  10000 7.279352948323840e+00 2.5e+02 7.15e-06  1e-07  1e-05 0:09.7
  600  12000 7.279352948161210e+00 4.0e+02 9.13e-07  1e-08  2e-06 0:11.8

=== Resonant-16C Restart 4/4 ===
(10_w,20)-aCMA-ES (mu_w=5.9,w_1=27%) in dimension 17 (seed=1477784359, Mon Dec  8 22:17:51 2025)
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
    1     20 2.835572268924469e+02 1.0e+00 2.76e-01  3e-01  3e-01 0:00.0
    2     40 4.539730024088169e+02 1.1e+00 2.66e-01  3e-01  3e-01 0:00.0
    3     60 8.793886231628213e+01 1.2e+00 2.59e-01  2e-01  3e-01 0:00.0
  100   2000 4.287940411408322e+01 1.0e+01 1.42e-01  4e-02  2e-01 0:01.5
  200   4000 9.817168040498007e+00 1.5e+01 4.56e-02  8e-03  6e-02 0:03.2
  300   6000 7.218234151876592e+00 5.2e+01 3.39e-03  2e-04  5e-03 0:04.6
  400   8000 7.216626050959563e+00 2.7e+02 1.38e-03  5e-05  3e-03 0:06.2
  500  10000 7.216556789378717e+00 4.2e+02 1.58e-05  3e-07  3e-05 0:08.2
  581  11620 7.216556787881613e+00 4.5e+02 3.90e-07  5e-09  5e-07 0:09.8

BEST Resonant-16C FIT:
[ 0.8547712  -1.21900827  0.42595456  0.0291953   0.20844277 -0.64194008
  0.89470344  1.4603491   0.83393488  0.84890353 -0.11832564  0.07085401
 -0.91775019  0.92356035  0.50104073 -0.12438562 -0.56027617]
cost = 7.21655678787628

Unpacked parameters:
A_u, B_u    = 0.8547711968519791 -1.2190082701856326
A_d, B_d    = 0.42595455736950455 0.02919529621174785
A_nu, B_nu  = 0.20844277469945544 -0.6419400847965141
shifts_u    = [0.89470344 1.4603491 ]
shifts_d    = [0.83393488 0.84890353]
shifts_nu   = [-0.11832564  0.07085401]
shifts_e    = [-0.91775019  0.92356035]
lambda_nu   = 0.5010407279220283
theta_C     = -0.12438562067937092
gamma_l     = -0.5602761722687236

Final evaluation at best fit:
  χ² = 5.944
  total cost (χ² + reg) = 7.217

Observables (model vs target, pull in σ):
  m_c/m_t     : model=0.0069334, target=0.007, pull=-0.032
  m_u/m_t     : model=9.9841e-06, target=1e-05, pull=-0.005
  m_s/m_b     : model=0.0238752, target=0.02, pull= 0.646
  m_d/m_b     : model=0.000385721, target=0.001, pull=-2.048
  m_mu/m_tau 
"""