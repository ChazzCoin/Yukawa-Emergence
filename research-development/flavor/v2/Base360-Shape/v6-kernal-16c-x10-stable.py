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
    Layout:
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
      14 : lambda_nu
      15 : theta_C   (Cabibbo twist)
    """
    X = np.array(X, dtype=float)

    A_u, B_u, A_d, B_d, A_nu, B_nu = X[0:6]
    X_u, Y_u = X[6:8]
    X_d, Y_d = X[8:10]
    X_nu, Y_nu = X[10:12]
    X_e, Y_e = X[12:14]
    lambda_nu, theta_C = X[14:16]

    # Phases bounded to [-π, π]
    A_u = clamp(A_u, -PI, PI)
    B_u = clamp(B_u, -PI, PI)
    A_d = clamp(A_d, -PI, PI)
    B_d = clamp(B_d, -PI, PI)
    A_nu = clamp(A_nu, -PI, PI)
    B_nu = clamp(B_nu, -PI, PI)

    # Exponent shifts: extended range [-1.5, 1.5]
    X_u, Y_u, X_d, Y_d, X_nu, Y_nu, X_e, Y_e = [
        clamp(v, -1.5, 1.5)
        for v in (X_u, Y_u, X_d, Y_d, X_nu, Y_nu, X_e, Y_e)
    ]

    # lambda_nu = clamp(lambda_nu, 0.0, 0.3)
    # theta_C   = clamp(theta_C,  -0.5, 0.5)  # up to ~30 degrees
    lambda_nu = clamp(lambda_nu, 0.0, 0.5)  # allow stronger resonance
    theta_C = clamp(theta_C, -0.5, 0.5)

    return (
        A_u, B_u, A_d, B_d, A_nu, B_nu,
        (X_u, Y_u), (X_d, Y_d), (X_nu, Y_nu), (X_e, Y_e),
        lambda_nu, theta_C
    )

# ================================================================
#  COST FUNCTION
# ================================================================

def resonant16C_cost(X, M0):
    try:
        (A_u, B_u, A_d, B_d, A_nu, B_nu,
         shifts_u, shifts_d, shifts_nu, shifts_e,
         lambda_nu, theta_C) = unpack_params_16C(X)

        gamma_q = 0.0
        gamma_l = 0.0

        # 9×9 Yukawas
        Yu9  = build_Y_sector(A_u,  B_u,  EXP_U_BASE,  shifts_u,  gamma_q)
        Yd9  = build_Y_sector(A_d,  B_d,  EXP_D_BASE,  shifts_d,  gamma_q)
        Ye9  = build_Y_sector(A_d,  B_d,  EXP_E_BASE,  shifts_e,  gamma_l)
        Ynu9 = build_Y_sector(A_nu, B_nu, EXP_NU_BASE, shifts_nu, gamma_l)

        # 9→3 Schur
        Yu = schur_9to3(Yu9)
        Yd = schur_9to3(Yd9)
        Ye = schur_9to3(Ye9)

        # Apply Cabibbo twist to down sector
        U_C = cabibbo_rotation(theta_C)
        Yd = U_C @ Yd

        # Neutrino projection + seesaw
        P = build_projection_resonant(lambda_nu)
        Ynu_eff = P @ Ynu9 @ P.conj().T
        Mnu = triadic_seesaw(M0, Ynu_eff)

        # Observables
        obs, Vckm, U_pmns = compute_observables(Yu, Yd, Ye, Mnu)
        chi2, pulls = chi2_from_obs(obs)

        # Soft penalties
        shift_penalty = (
            np.sum(np.array(shifts_u)**2) +
            np.sum(np.array(shifts_d)**2) +
            np.sum(np.array(shifts_nu)**2) +
            np.sum(np.array(shifts_e)**2)
        )
        # reg = 0.2 * shift_penalty + 2.0 * (lambda_nu**2) + 1.0 * (theta_C**2)
        reg = (
                0.2 * shift_penalty  # keep shifts modest
                + 0.5 * (lambda_nu ** 2)  # softer penalty: let neutrino resonance work
                + 1.0 * (theta_C ** 2)
        )
        return chi2 + reg

    except Exception:
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

        X0 = np.zeros(16)
        # phases
        X0[0:6]  = rng.uniform(-0.5, 0.5, size=6)
        # exponent shifts
        X0[6:14] = rng.normal(scale=0.4, size=8)
        # lambda_nu, theta_C
        X0[14]   = 0.15
        X0[15]   = 0.1

        es = cma.CMAEvolutionStrategy(
            X0,
            0.3,
            {
                'popsize': 20,
                'maxiter': 600,
                'CMA_diagonal': False,
                'seed': int(rng.integers(1, 2 ** 31 - 1)),  # or just a fixed integer, e.g. 12345
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
        A_u, B_u, A_d, B_d, A_nu, B_nu,
        shifts_u, shifts_d, shifts_nu, shifts_e,
        lambda_nu, theta_C
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

    gamma_q = gamma_l = 0.0

    Yu9  = build_Y_sector(A_u,  B_u,  EXP_U_BASE,  shifts_u,  gamma_q)
    Yd9  = build_Y_sector(A_d,  B_d,  EXP_D_BASE,  shifts_d,  gamma_q)
    Ye9  = build_Y_sector(A_d,  B_d,  EXP_E_BASE,  shifts_e,  gamma_l)
    Ynu9 = build_Y_sector(A_nu, B_nu, EXP_NU_BASE, shifts_nu, gamma_l)

    Yu = schur_9to3(Yu9)
    Yd = schur_9to3(Yd9)
    Ye = schur_9to3(Ye9)

    U_C = cabibbo_rotation(theta_C)
    Yd = U_C @ Yd

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
    optimize_resonant16C(num_restarts=4, seed=9)

"""


"""