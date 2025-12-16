import numpy as np
import math
import cma

# ================================================================
#  RESONANT-16+ FLAVOR MODEL
#  - 16 tunable parameters
#  - Quarks maximally coherent (gamma_q = 0)
#  - No explicit Cabibbo twist (CKM from geometry)
#  - Independent exponent shaping for charged leptons (X_e, Y_e)
# ================================================================

N_SITES = 9
LIGHT_SITES = [0, 1, 2]
HEAVY_SITES = [3, 4, 5, 6, 7, 8]
PI = math.pi


# ------------------------------------------------
# Basic utils
# ------------------------------------------------

def clamp(x, lo, hi):
    return np.minimum(np.maximum(x, lo), hi)

def generation_index(i: int) -> int:
    return i % 3


# ------------------------------------------------
# Triadic base exponents
# ------------------------------------------------

EXP_U_BASE  = np.array([4.0, 2.0, 0.0])   # up-type
EXP_D_BASE  = np.array([3.0, 2.0, 0.0])   # down-type
EXP_E_BASE  = np.array([3.0, 2.0, 0.0])   # charged leptons (same pattern, different shifts)
EXP_NU_BASE = np.array([1.0, 0.0, 0.0])   # neutrinos


# ================================================================
#  RESONANT PHASE WHEELS
# ================================================================

def build_phase_profile(A: float, B: float):
    """
    Triadic resonant phase:
      φ_g = A + B*g,  g = 0,1,2
    A,B ∈ [-π,π]
    """
    return np.array([A, A + B, A + 2*B], dtype=float)

def build_site_phases(phi_gen):
    """
    Lift generation phases to 9 sites via triadic indexing.
    """
    phi_site = np.zeros(N_SITES, dtype=float)
    for i in range(N_SITES):
        phi_site[i] = phi_gen[generation_index(i)]
    return phi_site

def build_phase_matrix(phi_site):
    """
    9×9 phase difference matrix: P_ij = exp(i (φ_i - φ_j)).
    """
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
    base_exp: length-3 (a,b,c)
    shifts: (X,Y) applied to generations 2 and 3 exponents.

      effective exponents:
        e0 = base[0]
        e1 = base[1] + X
        e2 = base[2] + Y

      magnitude per generation:
        s_g = (kappa)^e_g, with kappa ≈ 0.24
    """
    X, Y = shifts
    eff = np.array([
        base_exp[0],
        base_exp[1] + X,
        base_exp[2] + Y
    ], dtype=float)

    # alignment scale; chosen empirically from your good runs
    KAPPA = 0.24
    s_gen = np.power(KAPPA, eff)

    s = np.zeros(N_SITES, dtype=float)
    for i in range(N_SITES):
        s[i] = s_gen[generation_index(i)]
    return s


# ================================================================
#  GEOMETRIC COHERENCE KERNEL (γ)
# ================================================================

def build_kernel_gamma(gamma: float, forbidden_d: int = 2):
    """
    Toeplitz on the 9-site ring with a forbidden distance d=2.
      K_ij = exp(-γ * d_ij) for allowed distances,
      K_ij = 0 if d_ij = forbidden_d,
      K_ii = 1.
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
    """
    Triadic Fourier modes on Z9 for k=0,3,6.
    """
    n = 9
    j = np.arange(n)
    v0 = np.exp(2j*np.pi*0*j/n) / np.sqrt(n)
    v3 = np.exp(2j*np.pi*3*j/n) / np.sqrt(n)
    v6 = np.exp(2j*np.pi*6*j/n) / np.sqrt(n)
    return v0, v3, v6

def build_projection_resonant(lambda_nu: float):
    """
    λ_ν ∈ [0, 0.2]: small mixing of triadic Fourier modes before QR.

    Returns a 3×9 matrix with orthonormal rows.
    """
    v0, v3, v6 = triadic_modes()

    b0 = v0 + lambda_nu * v3
    b3 = v3 + lambda_nu * v6
    b6 = v6 + lambda_nu * v3

    B = np.vstack([b0, b3, b6])
    Q, _ = np.linalg.qr(B.conj().T)   # (9×3)
    return Q.T                        # (3×9)


# ================================================================
#  PROTO-MAJORANA (fixed per optimizer run)
# ================================================================

def proto_majorana(rng: np.random.Generator):
    M = rng.normal(size=(N_SITES, N_SITES)) + 1j * rng.normal(size=(N_SITES, N_SITES))
    M = 0.5 * (M + M.T)   # symmetric
    _, s, _ = np.linalg.svd(M)
    M /= s[0]
    return M


# ================================================================
#  9→3 REDUCTION (SCHUR)
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
#  YUKAWA SECTORS FROM RESONANT GEOMETRY
# ================================================================

def build_Y_sector(A, B, base_exp, shifts, gamma):
    """
    Build 9×9 Yukawa texture from:
      - resonant phases (A,B)
      - triadic exponents base_exp
      - exponent shifts shifts = (X,Y)
      - geometric coherence gamma
    """
    phi_gen = build_phase_profile(A, B)
    phi_site = build_site_phases(phi_gen)
    P = build_phase_matrix(phi_site)

    s = site_scales(base_exp, shifts)
    mag = np.outer(s, s)

    Y0 = mag * P
    K = build_kernel_gamma(gamma)
    Y = K * Y0

    # normalize to max singular value 1
    _, sv, _ = np.linalg.svd(Y)
    if sv[0] != 0:
        Y /= sv[0]
    return Y


# ================================================================
#  TRIADIC MAJORANA SEESAW
# ================================================================

def triadic_seesaw(M9: np.ndarray, Ynu_eff: np.ndarray) -> np.ndarray:
    """
    9D Majorana → heavy 6×6 → triadic 3×3 → seesaw with Dirac Yν_eff.
    """
    M_H = M9[np.ix_(HEAVY_SITES, HEAVY_SITES)]
    h = np.arange(len(HEAVY_SITES))

    B = np.zeros((len(HEAVY_SITES), 3), dtype=complex)
    for col, k in enumerate([1, 2, 3]):
        B[:, col] = np.exp(2j * np.pi * k * h / len(HEAVY_SITES)) / math.sqrt(len(HEAVY_SITES))

    M_R = B.conj().T @ M_H @ B
    M_R = 0.5 * (M_R + M_R.T)
    M_R += 1e-9 * np.eye(3)
    M_R *= 7e13    # heavy Majorana scale

    v = 174.0 / math.sqrt(2.0)
    mD = v * Ynu_eff
    M_Rinv = np.linalg.inv(M_R)

    return - mD @ M_Rinv @ mD.T


# ================================================================
#  OBSERVABLES + χ²
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

    # Normalize mass scales
    mu *= 173.0  / mu[0]
    md *= 4.18   / md[0]
    me *= 1.77686 / me[0]
    mnu *= 0.058 / mnu[0]

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
#  PARAMETER UNPACKING (16 PARAMETERS)
# ================================================================

def unpack_params_16plus(X):
    """
    Resonant-16+ layout:

      0  : A_u
      1  : B_u
      2  : A_d
      3  : B_d
      4  : A_nu
      5  : B_nu
      6  : X_u
      7  : Y_u
      8  : X_d
      9  : Y_d
      10 : X_nu
      11 : Y_nu
      12 : X_e
      13 : Y_e
      14 : gamma_l
      15 : lambda_nu
    """
    X = np.array(X, dtype=float)

    A_u, B_u, A_d, B_d, A_nu, B_nu = X[0:6]
    X_u, Y_u = X[6:8]
    X_d, Y_d = X[8:10]
    X_nu, Y_nu = X[10:12]
    X_e, Y_e = X[12:14]
    gamma_l, lambda_nu = X[14:16]

    # Bound phases to [-π, π]
    A_u = clamp(A_u, -PI, PI)
    B_u = clamp(B_u, -PI, PI)
    A_d = clamp(A_d, -PI, PI)
    B_d = clamp(B_d, -PI, PI)
    A_nu = clamp(A_nu, -PI, PI)
    B_nu = clamp(B_nu, -PI, PI)

    # Exponent shifts in a modest range
    for v in [X_u, Y_u, X_d, Y_d, X_nu, Y_nu, X_e, Y_e]:
        pass
    X_u, Y_u, X_d, Y_d, X_nu, Y_nu, X_e, Y_e = [
        clamp(v, -1.0, 1.0)
        for v in (X_u, Y_u, X_d, Y_d, X_nu, Y_nu, X_e, Y_e)
    ]

    gamma_l   = clamp(gamma_l, 0.0, 0.3)
    lambda_nu = clamp(lambda_nu, 0.0, 0.2)

    return (
        A_u, B_u, A_d, B_d, A_nu, B_nu,
        (X_u, Y_u), (X_d, Y_d), (X_nu, Y_nu), (X_e, Y_e),
        gamma_l, lambda_nu
    )


# ================================================================
#  COST FUNCTION (Resonant-16+)
# ================================================================

def resonant16plus_cost(X, M0):
    try:
        (A_u, B_u, A_d, B_d, A_nu, B_nu,
         shifts_u, shifts_d, shifts_nu, shifts_e,
         gamma_l, lambda_nu) = unpack_params_16plus(X)

        gamma_q = 0.0  # maximally coherent quark sectors

        # 9×9 Yukawas
        Yu9  = build_Y_sector(A_u,  B_u,  EXP_U_BASE,  shifts_u,  gamma_q)
        Yd9  = build_Y_sector(A_d,  B_d,  EXP_D_BASE,  shifts_d,  gamma_q)
        Ye9  = build_Y_sector(A_d,  B_d,  EXP_E_BASE,  shifts_e,  gamma_l)
        Ynu9 = build_Y_sector(A_nu, B_nu, EXP_NU_BASE, shifts_nu, gamma_l)

        # 9→3 Schur for Dirac sectors
        Yu = schur_9to3(Yu9)
        Yd = schur_9to3(Yd9)
        Ye = schur_9to3(Ye9)

        # Neutrino 9→3 projection
        P = build_projection_resonant(lambda_nu)
        Ynu_eff = P @ Ynu9 @ P.conj().T

        # Seesaw
        Mnu = triadic_seesaw(M0, Ynu_eff)

        # Observables
        obs, Vckm, U_pmns = compute_observables(Yu, Yd, Ye, Mnu)
        chi2, pulls = chi2_from_obs(obs)

        # Mild regularization: keep exponent shifts and gamma_l, lambda_nu moderate
        shift_penalty = (
            np.sum(np.array(shifts_u)**2) +
            np.sum(np.array(shifts_d)**2) +
            np.sum(np.array(shifts_nu)**2) +
            np.sum(np.array(shifts_e)**2)
        )
        reg = 0.1 * shift_penalty + 2.0 * (gamma_l**2) + 2.0 * (lambda_nu**2)

        return chi2 + reg

    except Exception:
        return 1e9


# ================================================================
#  OPTIMIZER DRIVER
# ================================================================

def optimize_resonant16plus(num_restarts=4, seed=7):
    rng = np.random.default_rng(seed)
    M0 = proto_majorana(rng)

    best_cost = np.inf
    best_X = None

    for r in range(num_restarts):
        print(f"\n=== Resonant-16+ Restart {r+1}/{num_restarts} ===")

        X0 = np.zeros(16)
        # randomize phases & exponent shifts a bit
        X0[0:6]  = rng.uniform(-0.5, 0.5, size=6)   # phases (smallish start)
        X0[6:14] = rng.normal(scale=0.3, size=8)    # exponent shifts
        X0[14]   = 0.1                              # gamma_l
        X0[15]   = 0.1                              # lambda_nu

        es = cma.CMAEvolutionStrategy(
            X0,
            0.3,
            {
                'popsize': 20,
                'maxiter': 600,
                'CMA_diagonal': False
            }
        )

        while not es.stop():
            xs = es.ask()
            cs = [resonant16plus_cost(x, M0) for x in xs]
            es.tell(xs, cs)
            es.disp()

        if es.best.f < best_cost:
            best_cost = es.best.f
            best_X = es.best.x.copy()

    print("\nBEST Resonant-16+ FIT:")
    print(best_X)
    print("cost =", best_cost)

    # Final diagnostic at best fit
    (
        A_u, B_u, A_d, B_d, A_nu, B_nu,
        shifts_u, shifts_d, shifts_nu, shifts_e,
        gamma_l, lambda_nu
    ) = unpack_params_16plus(best_X)

    print("\nUnpacked parameters:")
    print("A_u, B_u    =", A_u, B_u)
    print("A_d, B_d    =", A_d, B_d)
    print("A_nu, B_nu  =", A_nu, B_nu)
    print("shifts_u    =", shifts_u)
    print("shifts_d    =", shifts_d)
    print("shifts_nu   =", shifts_nu)
    print("shifts_e    =", shifts_e)
    print("gamma_l     =", gamma_l)
    print("lambda_nu   =", lambda_nu)

    # Recompute observables to show pulls
    gamma_q = 0.0
    Yu9  = build_Y_sector(A_u,  B_u,  EXP_U_BASE,  shifts_u,  gamma_q)
    Yd9  = build_Y_sector(A_d,  B_d,  EXP_D_BASE,  shifts_d,  gamma_q)
    Ye9  = build_Y_sector(A_d,  B_d,  EXP_E_BASE,  shifts_e,  gamma_l)
    Ynu9 = build_Y_sector(A_nu, B_nu, EXP_NU_BASE, shifts_nu, gamma_l)

    Yu = schur_9to3(Yu9)
    Yd = schur_9to3(Yd9)
    Ye = schur_9to3(Ye9)

    P = build_projection_resonant(lambda_nu)
    Ynu_eff = P @ Ynu9 @ P.conj().T
    Mnu = triadic_seesaw(M0, Ynu_eff)

    obs, Vckm, U_pmns = compute_observables(Yu, Yd, Ye, Mnu)
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
    optimize_resonant16plus(num_restarts=4, seed=9)

"""

=== Resonant-16+ Restart 1/4 ===
(10_w,20)-aCMA-ES (mu_w=5.9,w_1=27%) in dimension 16 (seed=438244, Tue Dec  9 13:50:42 2025)
/Users/chazzromeo/rAI/Flavor/.venv/lib/python3.11/site-packages/cma/s.py:15: UserWarning: Could not import matplotlib.pyplot, therefore ``cma.plot()`` etc. is not available
  _warnings.warn('Could not import matplotlib.pyplot, therefore'
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
    1     20 1.322553869213175e+02 1.0e+00 2.87e-01  3e-01  3e-01 0:00.0
    2     40 1.000788237730388e+02 1.2e+00 2.87e-01  3e-01  3e-01 0:00.0
    3     60 1.112157036009610e+02 1.2e+00 2.83e-01  3e-01  3e-01 0:00.0
  100   2000 2.957959861285722e+01 2.3e+01 1.04e-01  3e-02  2e-01 0:01.6
  200   4000 2.900139946630192e+01 1.1e+02 4.89e-03  4e-04  8e-03 0:03.1
  300   6000 2.900125085547902e+01 2.9e+02 2.39e-04  6e-06  5e-04 0:04.9
  400   8000 2.900125077452089e+01 8.2e+02 1.17e-05  5e-08  2e-05 0:06.6
  500  10000 2.900125077438054e+01 1.3e+03 2.16e-06  6e-09  4e-06 0:08.3
  600  12000 2.900125077437808e+01 1.8e+03 2.11e-06  5e-09  4e-06 0:10.1

=== Resonant-16+ Restart 2/4 ===
(10_w,20)-aCMA-ES (mu_w=5.9,w_1=27%) in dimension 16 (seed=361351, Tue Dec  9 13:50:53 2025)
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
    1     20 1.204101607240687e+02 1.0e+00 2.74e-01  3e-01  3e-01 0:00.0
    2     40 1.087290347182405e+02 1.1e+00 2.70e-01  3e-01  3e-01 0:00.0
    3     60 1.209213003175906e+02 1.2e+00 2.69e-01  3e-01  3e-01 0:00.1
  100   2000 3.077790981753149e+01 1.3e+01 1.01e-01  4e-02  1e-01 0:01.9
  200   4000 2.947116827458095e+01 6.2e+01 3.79e-03  5e-04  6e-03 0:03.9
  300   6000 2.947109401927158e+01 3.7e+02 3.01e-04  9e-06  6e-04 0:05.4
  400   8000 2.947109398905001e+01 1.2e+03 9.84e-06  1e-07  3e-05 0:06.8
  500  10000 2.947109398896676e+01 1.3e+03 7.12e-06  8e-08  2e-05 0:08.3
  600  12000 2.947109398892366e+01 2.3e+03 6.34e-06  7e-08  2e-05 0:09.7

=== Resonant-16+ Restart 3/4 ===
(10_w,20)-aCMA-ES (mu_w=5.9,w_1=27%) in dimension 16 (seed=458162, Tue Dec  9 13:51:02 2025)
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
    1     20 1.921359827385857e+02 1.0e+00 3.02e-01  3e-01  3e-01 0:00.0
    2     40 8.841803956936850e+01 1.2e+00 2.90e-01  3e-01  3e-01 0:00.0
    3     60 1.102926542140527e+02 1.3e+00 2.82e-01  3e-01  3e-01 0:00.0
  100   2000 3.287385698448514e+01 2.1e+01 2.53e-01  8e-02  5e-01 0:01.4
  200   4000 2.947127583762240e+01 1.1e+02 1.07e-02  1e-03  2e-02 0:02.8
  300   6000 2.947109400886301e+01 7.3e+02 3.20e-04  7e-06  1e-03 0:04.2
  400   8000 2.947109398902386e+01 2.0e+03 1.52e-05  2e-07  5e-05 0:05.9
  500  10000 2.947109398894364e+01 3.1e+03 1.15e-05  1e-07  4e-05 0:07.7
  600  12000 2.947109398884091e+01 3.4e+03 8.70e-06  7e-08  3e-05 0:09.5

=== Resonant-16+ Restart 4/4 ===
(10_w,20)-aCMA-ES (mu_w=5.9,w_1=27%) in dimension 16 (seed=415662, Tue Dec  9 13:51:12 2025)
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
    1     20 1.003759160697665e+02 1.0e+00 2.77e-01  3e-01  3e-01 0:00.0
    2     40 1.009031457442193e+02 1.1e+00 2.85e-01  3e-01  3e-01 0:00.0
    3     60 1.003767084082229e+02 1.2e+00 2.88e-01  3e-01  3e-01 0:00.0
  100   2000 3.165018351583797e+01 2.4e+01 5.06e-02  2e-02  8e-02 0:01.8
  200   4000 3.156833123192398e+01 2.0e+02 4.62e-04  3e-05  1e-03 0:03.6
  300   6000 3.156833103257412e+01 1.2e+03 1.77e-05  2e-07  6e-05 0:05.4
  400   8000 3.156833103234565e+01 2.1e+03 6.06e-06  7e-08  2e-05 0:07.1
  500  10000 3.156833103226636e+01 2.1e+03 6.81e-06  7e-08  2e-05 0:08.6
  600  12000 3.156833103249056e+01 2.2e+03 2.69e-06  2e-08  7e-06 0:10.2

BEST Resonant-16+ FIT:
[-2.06318929e-01 -5.48084118e-06 -7.72289311e-01 -8.90802082e-07
 -1.01428504e+00  3.53905063e-01  1.87376617e+00  1.85676271e+00
  5.85003324e-01  9.04143913e-01 -9.60109686e-01 -1.00003832e+00
  3.92958144e-01  1.49829531e+00 -8.37749133e-01  9.40334596e-01]
cost = 29.00125077421206

Unpacked parameters:
A_u, B_u    = -0.20631892921575462 -5.480841181215517e-06
A_d, B_d    = -0.7722893105723955 -8.908020818799611e-07
A_nu, B_nu  = -1.0142850425453311 0.3539050628891003
shifts_u    = (np.float64(1.0), np.float64(1.0))
shifts_d    = (np.float64(0.58500332375431), np.float64(0.9041439126373768))
shifts_nu   = (np.float64(-0.9601096859597213), np.float64(-1.0))
shifts_e    = (np.float64(0.39295814401315565), np.float64(1.0))
gamma_l     = 0.0
lambda_nu   = 0.2

Final evaluation at best fit:
  χ² = 28.298
  total cost (χ² + reg) = 29.001

Observables (model vs target, pull in σ):
  m_c/m_t     : model=0.00180779, target=0.007, pull=-2.472
  m_u/m_t     : model=8.3716e-06, target=1e-05, pull=-0.543
  m_s/m_b     : model=0.0260497, target=0.02, pull= 1.008
  m_d/m_b     : model=0.000356539, target=0.001, pull=-2.145
  m_mu/m_tau  : model=0.0332924, target=0.06, pull=-1.484
  m_e/m_tau   : model=0.000295307, target=0.0003, pull=-0.052
  theta12_q   : model=0.118164, target=0.226, pull=-1.591
  theta23_q   : model=0.0469116, target=0.041, pull= 0.481
  theta13_q   : model=0.00356571, target=0.0035, pull= 0.063
  theta12_l   : model=1.14375, target=0.59, pull= 3.129
  theta23_l   : model=0.947738, target=0.84, pull= 0.428
  theta13_l   : model=0.151935, target=0.15, pull= 0.043
  Delta m2_21 : model=7.54883e-05, target=7.4e-05, pull= 0.067
  Delta m2_31 : model=0.00336396, target=0.0025, pull= 1.152

|V_CKM| ≈
[[0.99302047 0.11788823 0.0035657 ]
 [0.11759324 0.99195402 0.04689409]
 [0.00906527 0.04614748 0.9988935 ]]

|U_PMNS| ≈
[[0.40941129 0.89970843 0.15135103]
 [0.58202562 0.12981091 0.80274237]
 [0.70258706 0.41674196 0.57679941]]
"""