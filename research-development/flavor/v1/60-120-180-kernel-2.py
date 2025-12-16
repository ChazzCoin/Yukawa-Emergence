#!/usr/bin/env python3
# =============================================================================
#   TRIADIC GEOMETRIC ALIGNMENT — D360 TRIAD (60–120–180) IMPLEMENTATION
#
#   9 continuous parameters:
#      A_u, B_u,
#      A_d, B_d,
#      A_e, B_e,
#      A_nu, B_nu,
#      kappa      ← triadic coherence
#
#   Geometry: 9-site ring (triadic indexing)
#   Phases:   120-cycle phase wheels φ(i) = A + B (i mod 3)
#   Kernel:   D360 triadic seed n=60 → (60,120,180), lifted to 9 sites
#
#   This script:
#      - builds 9×9 Yukawas from the triadic kernel + phase wheels
#      - downfolds to 3×3 (Schur complement)
#      - constructs a 3×3 seesaw neutrino mass matrix
#      - runs 1-loop SM RGEs (Yu,Yd,Ye,κ)
#      - fits 14 flavor observables with CMA-ES
# =============================================================================

import numpy as np
from scipy.integrate import solve_ivp
import cma

# --------------------------------------------------------------------
#  Physics constants
# --------------------------------------------------------------------
MU_HIGH = 2.0e14   # flavor / seesaw scale
MU_LOW  = 1.0e2    # EW scale
V_HIGGS = 246.0    # Higgs vev in GeV

# EW couplings at MU_LOW
g1_EW, g2_EW, g3_EW = 0.36, 0.65, 1.17
lam_H = 0.13

N_SITES = 9

# --------------------------------------------------------------------
#  Experimental targets (14 observables)
# --------------------------------------------------------------------
targets = {
    "m_c/m_t":   0.007,
    "m_u/m_t":   1e-5,
    "m_s/m_b":   0.02,
    "m_d/m_b":   0.001,
    "m_mu/m_tau":0.06,
    "m_e/m_tau": 3e-4,

    "theta12_q": 0.226,  # Cabibbo
    "theta23_q": 0.041,
    "theta13_q": 0.0035,

    "theta12_l": 0.59,
    "theta23_l": 0.84,
    "theta13_l": 0.15,

    "Delta_m2_21": 7.4e-5,
    "Delta_m2_31": 2.5e-3,
}
sigmas = {k: 0.3*abs(v) for k, v in targets.items()}

# ====================================================================
# 1a. (Optional) Original D120 kernel (kept for reference / debugging)
# ====================================================================
def kernel_D120_axiom(kappa, forbid_long_arc=False):
    """
    D120 axiom-projector kernel on 9 sites.

    Sites are at angles theta_i = 40° * i on a 360° circle,
    then folded into a 120° fundamental domain via mod 120.

    Distances on the 120° circle are 0, 40, or 80 degrees.
    """
    theta_360 = 40.0 * np.arange(N_SITES)     # 0,40,...,320
    theta = np.mod(theta_360, 120.0)          # fold to 120°

    K = np.zeros((N_SITES, N_SITES), dtype=float)

    for i in range(N_SITES):
        for j in range(N_SITES):
            d = abs(theta[i] - theta[j])
            d = min(d, 120.0 - d)  # geodesic distance on 120° circle

            if d < 1e-12:
                K[i, j] = 1.0
            elif abs(d - 40.0) < 1e-12:
                K[i, j] = kappa
            elif abs(d - 80.0) < 1e-12:
                K[i, j] = 0.0 if forbid_long_arc else kappa**2
            else:
                K[i, j] = 0.0
    return K

# ====================================================================
# 1b. D360 Triadic 60–120–180 Kernel lifted to 9 sites
# ====================================================================
def triad_kernel_60_120_180(kappa: float) -> np.ndarray:
    """
    3x3 triadic kernel in harmonic space for modes (60, 120, 180).

    Distances in degrees: 0, 60, 120.
    Map them to powers of kappa: 1, kappa, kappa**2.
    """
    modes = np.array([60.0, 120.0, 180.0], dtype=float)
    dist = np.abs(modes[:, None] - modes[None, :])  # 3x3 matrix

    K3 = np.zeros((3, 3), dtype=float)
    for i in range(3):
        for j in range(3):
            d = dist[i, j]
            if d < 1e-9:
                K3[i, j] = 1.0
            elif np.isclose(d, 60.0):
                K3[i, j] = kappa
            elif np.isclose(d, 120.0):
                K3[i, j] = kappa**2
            else:
                # Should not occur for (60,120,180); kept for safety
                K3[i, j] = 0.0
    return K3


def triad_to_9site_projector() -> np.ndarray:
    """
    3x9 projector from triad space (3 modes) to a 9-site ring.

    Mode 0 (60)  -> sites (0,3,6)
    Mode 1 (120) -> sites (1,4,7)
    Mode 2 (180) -> sites (2,5,8)
    """
    P = np.zeros((3, 9), dtype=complex)
    triads = [(0, 3, 6), (1, 4, 7), (2, 5, 8)]
    for a, sites in enumerate(triads):
        for s in sites:
            P[a, s] = 1.0 / np.sqrt(3.0)
    return P


def triad_kernel_9x9(kappa: float) -> np.ndarray:
    """
    Lift the 3x3 triadic kernel K3 (60,120,180) to a 9x9 kernel.

    K9 = P^† K3 P, where P is the triad_to_9site_projector.
    """
    K3 = triad_kernel_60_120_180(kappa)
    P = triad_to_9site_projector()
    return P.conj().T @ K3 @ P

# ====================================================================
# 2. Phase matrix — triadic phase wheel
# ====================================================================
def phase_matrix(A, B):
    """
    φ(i) = A + B*(i mod 3)
    Triadic periodicity → 120° geometric phase unit.
    """
    phi_gen = np.array([A + B*g for g in range(3)])
    phi = np.array([phi_gen[i % 3] for i in range(N_SITES)])
    # matrix of e^{i(φ_i - φ_j)}
    return np.exp(1j * (phi[:, None] - phi[None, :]))

# ====================================================================
# 3. Yukawa builder (triadic D360 60–120–180 kernel)
# ====================================================================
def build_Yukawa(A, B, kappa, alpha, forbid_long_arc=False):
    """
    Aligned 9x9 Yukawa using D360 triadic 60-120-180 kernel.

    Parameters
    ----------
    A, B : float
        Phase-wheel parameters for φ(i) = A + B * (i mod 3).
    kappa : float
        Triadic coherence parameter.
    alpha : float
        Overall sector normalization (fixed per sector).
    forbid_long_arc : bool, optional
        Kept for API compatibility; ignored in the triadic kernel version.
    """
    phases = phase_matrix(A, B)      # 9x9 complex phase differences
    K = triad_kernel_9x9(kappa)      # 9x9 real triadic kernel from (60,120,180)
    Y9 = phases * K

    # Normalize largest singular value to 1 and rescale to alpha
    s_max = np.linalg.svd(Y9, compute_uv=False)[0]
    if s_max == 0:
        raise RuntimeError("Degenerate Yukawa: largest singular value is zero.")
    Y9 /= s_max
    Y9 *= alpha
    return Y9

# ====================================================================
# 4. Schur 9→3
# ====================================================================
def schur_9to3(Y9: np.ndarray) -> np.ndarray:
    """
    Downfold 9x9 matrix to 3x3 via Schur complement.

    Y9 = [[A, B],
          [B†, D]]

    Yeff = A - B D^{-1} B†
    """
    A = Y9[:3, :3]
    B = Y9[:3, 3:]
    D = Y9[3:, 3:]
    Dinv = np.linalg.pinv(D + 1e-10*np.eye(6))
    return A - B @ Dinv @ B.conj().T

# ====================================================================
# 5. Proto-Majorana 9x9 matrix
# ====================================================================
def proto_majorana(rng, Lambda_Maj=7e13):
    """
    Random complex symmetric 9x9 matrix with largest singular value = Lambda_Maj.
    """
    M = rng.normal(size=(9, 9)) + 1j * rng.normal(size=(9, 9))
    M = 0.5 * (M + M.T.conj())
    M /= np.linalg.svd(M, compute_uv=False)[0]
    M *= Lambda_Maj
    return M

def summarize_results(Yu, Yd, Ye, Mnu, label=""):
    su = np.sort(np.linalg.svd(Yu, compute_uv=False))
    sd = np.sort(np.linalg.svd(Yd, compute_uv=False))
    se = np.sort(np.linalg.svd(Ye, compute_uv=False))

    print(f"\n=== FLAVOR SUMMARY {label} ===")

    print("\nMass hierarchies:")
    print(f"  up:    mu/mc/mt = {su[0]:.3e}, {su[1]:.3e}, {su[2]:.3e}")
    print(f"  down:  md/ms/mb = {sd[0]:.3e}, {sd[1]:.3e}, {sd[2]:.3e}")
    print(f"  lepton:me/mm/mt = {se[0]:.3e}, {se[1]:.3e}, {se[2]:.3e}")

    Uu = np.linalg.svd(Yu)[0]
    Ud = np.linalg.svd(Yd)[0]
    Vckm = Uu.conj().T @ Ud

    def angles(U):
        a = np.abs(U)
        s13 = a[0, 2]
        c13 = np.sqrt(max(1e-12, 1 - s13**2))
        s12 = a[0, 1]/c13
        s23 = a[1, 2]/c13
        return np.degrees(np.arcsin([s12, s23, s13]))

    print("\nCKM angles [deg]:")
    print(" ", angles(Vckm))

    evals, U_nu = np.linalg.eigh(0.5*(Mnu + Mnu.T))
    mnu = np.sort(np.abs(evals))

    Ue = np.linalg.svd(Ye)[0]
    Upmns = Ue.conj().T @ U_nu

    print("\nPMNS angles [deg]:")
    print(" ", angles(Upmns))

    print("\nNeutrino splittings:")
    print(f"  Δm21² = {mnu[1]**2 - mnu[0]**2:.3e}")
    print(f"  Δm31² = {mnu[2]**2 - mnu[0]**2:.3e}")

def robustness_scan_simple(M0, n=30):
    rng = np.random.default_rng(2025)
    ok = 0

    for _ in range(n):
        X = rng.uniform(-0.5, 0.5, 9)
        X[-1] = rng.uniform(0.2, 0.9)
        try:
            if cost(X, M0) < 2e5:
                ok += 1
        except:
            pass

    print("\n=== ROBUSTNESS / EXISTENCE CHECK ===")
    print(f"Viable points: {ok} / {n}")

def publication_proof(best_x, M0):
    A_u, B_u, A_d, B_d, A_e, B_e, A_nu, B_nu, kappa = best_x
    alpha_u, alpha_d, alpha_e, alpha_nu = 0.71, 0.095, 0.082, 0.13

    # Build full 9×9 Yukawas
    Yu9  = build_Yukawa(A_u,  B_u,  kappa, alpha_u)
    Yd9  = build_Yukawa(A_d,  B_d,  kappa, alpha_d)
    Ye9  = build_Yukawa(A_e,  B_e,  kappa, alpha_e)
    Ynu9 = build_Yukawa(A_nu, B_nu, kappa, alpha_nu)

    # Downfold
    Yu = schur_9to3(Yu9)
    Yd = schur_9to3(Yd9)
    Ye = schur_9to3(Ye9)

    # Neutrinos
    P = np.zeros((3, 9), dtype=complex)
    for c, sites in enumerate([(0,3,6),(1,4,7),(2,5,8)]):
        P[c, sites] = 1/np.sqrt(3)

    Ynu = P @ Ynu9 @ P.conj().T
    MR  = P @ M0   @ P.conj().T
    Mnu = -0.5 * V_HIGGS**2 * (
        Ynu @ np.linalg.pinv(MR + 1e-8*np.eye(3)) @ Ynu.T
    )

    # Run RG
    Yu_l, Yd_l, Ye_l, kappa_l = run_rge(Yu, Yd, Ye, Mnu / V_HIGGS**2)
    Mnu_l = kappa_l * V_HIGGS**2

    # ---- PRINT PROOF ----
    print("\n==============================")
    print("   PUBLICATION PROOF OUTPUT")
    print("==============================")

    print("\n--- Emergence: 9 → 3 modes ---")
    print("9x9 Yu singular values:", np.linalg.svd(Yu9, compute_uv=False))
    print("3x3 Yu singular values:", np.linalg.svd(Yu_l, compute_uv=False))

    summarize_results(Yu_l, Yd_l, Ye_l, Mnu_l, label="(FINAL LOW SCALE)")

    obs = get_obs(Yu_l, Yd_l, Ye_l, Mnu_l)
    print_observable_table(obs)

# ====================================================================
# 6. RGE System (1-loop SM Yukawas + κ operator)
# ====================================================================
def pack(Yu, Yd, Ye, kappa):
    def flat(M): return np.concatenate([M.real.ravel(), M.imag.ravel()])
    return np.concatenate([flat(Yu), flat(Yd), flat(Ye), flat(kappa)])

def unpack(v):
    n = 3
    N = n * n
    def block(i):
        re = v[i:i+N].reshape((n, n))
        im = v[i+N:i+2*N].reshape((n, n))
        return re + 1j * im
    Yu = block(0)
    Yd = block(2*N)
    Ye = block(4*N)
    kappa = block(6*N)
    return Yu, Yd, Ye, kappa

def beta(t, v, g1, g2, g3, lam):
    Yu, Yd, Ye, kappa = unpack(v)

    # Safety clipping to avoid extreme runaway
    for M in (Yu, Yd, Ye):
        np.clip(M, -20, 20, out=M)

    T = np.trace(3*Yu@Yu.conj().T + 3*Yd@Yd.conj().T + Ye@Ye.conj().T).real
    pref = 1.0 / (16.0 * np.pi**2)

    dYu = pref * (
        Yu * (T - (17.0/20.0)*g1**2 - (9.0/4.0)*g2**2 - 8.0*g3**2)
        + 1.5 * (Yu @ Yu.conj().T @ Yu - Yd @ Yd.conj().T @ Yu)
    )

    dYd = pref * (
        Yd * (T - (1.0/4.0)*g1**2 - (9.0/4.0)*g2**2 - 8.0*g3**2)
        + 1.5 * (Yd @ Yd.conj().T @ Yd - Yu @ Yu.conj().T @ Yd)
    )

    dYe = pref * (
        Ye * (T - (9.0/4.0)*g1**2 - (9.0/4.0)*g2**2)
        + 1.5 * (Ye @ Ye.conj().T @ Ye)
    )

    YeT = Ye @ Ye.conj().T
    dkappa = pref * ((-3.0*g2**2 + lam)*kappa + YeT @ kappa + kappa @ YeT.T)

    return pack(dYu, dYd, dYe, dkappa)

def run_rge(Yu, Yd, Ye, kappa_high):
    """
    Run RGEs from MU_HIGH down to MU_LOW.
    """
    sol = solve_ivp(
        beta,
        [np.log(MU_HIGH), np.log(MU_LOW)],
        pack(Yu, Yd, Ye, kappa_high),
        args=(g1_EW, g2_EW, g3_EW, lam_H),
        method="RK45",
        rtol=1e-5,
        atol=1e-8,
        max_step=0.4,
    )
    return unpack(sol.y[:, -1])

# ====================================================================
# 7. Observables extraction
# ====================================================================
def get_obs(Yu, Yd, Ye, Mnu):
    """
    Compute 14 observables from Yukawas + neutrino mass matrix.
    """
    def angles(U):
        a = np.abs(U)
        s13 = a[0, 2]
        c13 = np.sqrt(max(0.0, 1.0 - s13**2))
        s12 = a[0, 1]/c13 if c13 > 1e-12 else 0.0
        s23 = a[1, 2]/c13 if c13 > 1e-12 else 0.0
        # clamp to [-1,1] to avoid NaNs
        s12 = min(1.0, max(0.0, s12))
        s23 = min(1.0, max(0.0, s23))
        s13 = min(1.0, max(0.0, s13))
        return np.arcsin(s12), np.arcsin(s23), np.arcsin(s13)

    su = np.sort(np.linalg.svd(Yu, compute_uv=False))
    sd = np.sort(np.linalg.svd(Yd, compute_uv=False))
    se = np.sort(np.linalg.svd(Ye, compute_uv=False))

    obs = {
        "m_c/m_t":    su[1]/su[2],
        "m_u/m_t":    su[0]/su[2],
        "m_s/m_b":    sd[1]/sd[2],
        "m_d/m_b":    sd[0]/sd[2],
        "m_mu/m_tau": se[1]/se[2],
        "m_e/m_tau":  se[0]/se[2],
    }

    # CKM
    Uu = np.linalg.svd(Yu)[0]
    Ud = np.linalg.svd(Yd)[0]
    Vckm = Uu.conj().T @ Ud
    obs["theta12_q"], obs["theta23_q"], obs["theta13_q"] = angles(Vckm)

    # Neutrinos
    evals, U_nu = np.linalg.eigh(0.5 * (Mnu + Mnu.T))
    mnu = np.sort(np.abs(evals))

    Ue = np.linalg.svd(Ye)[0]
    Upmns = Ue.conj().T @ U_nu
    obs["theta12_l"], obs["theta23_l"], obs["theta13_l"] = angles(Upmns)

    obs["Delta_m2_21"] = mnu[1]**2 - mnu[0]**2
    obs["Delta_m2_31"] = mnu[2]**2 - mnu[0]**2

    return obs

# ====================================================================
# 8. Cost function
# ====================================================================
def cost(X, M0):
    """
    X = (A_u,B_u, A_d,B_d, A_e,B_e, A_nu,B_nu, kappa)
    """
    A_u, B_u, A_d, B_d, A_e, B_e, A_nu, B_nu, kappa = X

    # Fixed physical normalizations
    alpha_u, alpha_d, alpha_e, alpha_nu = 0.71, 0.095, 0.082, 0.13

    # Build 9×9 Yukawas (triadic D360 kernel)
    Yu9  = build_Yukawa(A_u,  B_u,  kappa, alpha_u)
    Yd9  = build_Yukawa(A_d,  B_d,  kappa, alpha_d)
    Ye9  = build_Yukawa(A_e,  B_e,  kappa, alpha_e)
    Ynu9 = build_Yukawa(A_nu, B_nu, kappa, alpha_nu)

    # Schur 9→3 for quarks and charged leptons
    Yu_h = schur_9to3(Yu9)
    Yd_h = schur_9to3(Yd9)
    Ye_h = schur_9to3(Ye9)

    if np.random.rand() < 0.001:  # occasional diagnostic
        print("9x9 Yu singular values:", np.linalg.svd(Yu9, compute_uv=False))
        print("3x3 Yu singular values:", np.linalg.svd(Yu_h, compute_uv=False))

    # Triadic projector for neutrinos (3 generations ← 9 sites)
    P = np.zeros((3, 9), dtype=complex)
    for c, sites in enumerate([(0, 3, 6), (1, 4, 7), (2, 5, 8)]):
        P[c, sites] = 1.0 / np.sqrt(3.0)

    Ynu_eff = P @ Ynu9 @ P.conj().T
    MR      = P @ M0   @ P.conj().T

    # Type-I seesaw at high scale
    Mnu_h = -0.5 * V_HIGGS**2 * (
        Ynu_eff @ np.linalg.pinv(MR + 1e-8*np.eye(3)) @ Ynu_eff.T
    )
    kappa_h = Mnu_h / V_HIGGS**2

    # RGE evolution to low scale
    Yu_l, Yd_l, Ye_l, kappa_l = run_rge(Yu_h, Yd_h, Ye_h, kappa_h)
    Mnu_l = kappa_l * V_HIGGS**2

    # Observables
    obs = get_obs(Yu_l, Yd_l, Ye_l, Mnu_l)

    chi2 = sum(((obs[k] - targets[k]) / sigmas[k])**2 for k in targets)
    reg  = 0.05 * np.sum(X**2)   # mild regularization on parameters
    if np.random.rand() < 0.001:  # occasional diagnostic
        print("9x9 Yu singular values:", np.linalg.svd(Yu9, compute_uv=False))
        print("3x3 Yu singular values:", np.linalg.svd(Yu_h, compute_uv=False))

    return chi2 + reg

# ====================================================================
# 9. RUN OPTIMIZATION
# ====================================================================
if __name__ == "__main__":
    rng = np.random.default_rng(777)
    M0 = proto_majorana(rng)

    # Initial guess: A,B in [-0.4,0.4], kappa ~ 0.27
    x0 = np.concatenate([rng.uniform(-0.4, 0.4, 8), [0.27]])

    es = cma.CMAEvolutionStrategy(
        x0, 0.35,
        {'popsize': 80, 'maxiter': 4000, 'verb_disp': 1, 'seed': 42}
    )

    print("Starting D360 triadic 60–120–180 alignment optimization...")
    while not es.stop():
        xs = es.ask()
        es.tell(xs, [cost(x, M0) for x in xs])
        es.disp()

    print("\n=== FINAL TRIADIC D360 RESULT ===")
    print("Best χ² + reg =", es.best.f)
    print("Best parameters:", es.best.x)

    # ---- PRODUCE PUBLICATION OUTPUT ----
    publication_proof(es.best.x, M0)

    # ---- ROBUSTNESS / EXISTENCE ----
    robustness_scan_simple(M0, n=40)

"""
RESULTS:
(40_w,80)-aCMA-ES (mu_w=21.8,w_1=9%) in dimension 9 (seed=42, Mon Dec 15 21:57:13 2025)
Starting D360 triadic 60–120–180 alignment optimization...
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
    1     80 1.890384674331038e+06 1.0e+00 4.32e-01  4e-01  5e-01 0:01.8
    2    160 1.150369321155479e+06 1.4e+00 5.42e-01  5e-01  6e-01 0:03.7
    3    240 1.224781811774196e+06 1.5e+00 5.95e-01  5e-01  7e-01 0:05.5
    4    320 1.177923438241167e+06 1.7e+00 6.25e-01  5e-01  7e-01 0:07.3
    5    400 9.216713354399954e+05 2.0e+00 6.15e-01  5e-01  7e-01 0:09.2
    6    480 1.319293360023705e+06 2.0e+00 5.73e-01  4e-01  7e-01 0:11.1
    7    560 4.518454070058948e+05 2.3e+00 5.94e-01  4e-01  7e-01 0:12.9
    8    640 4.101303454502988e+05 2.3e+00 5.86e-01  4e-01  7e-01 0:14.7
    9    720 4.361710906851634e+05 2.5e+00 5.06e-01  3e-01  6e-01 0:16.5
   10    800 1.012539515645331e+06 2.6e+00 5.63e-01  3e-01  6e-01 0:18.3
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
   11    880 4.648458760663532e+05 2.9e+00 5.77e-01  3e-01  7e-01 0:20.2
   12    960 8.454162229922885e+05 3.3e+00 5.40e-01  2e-01  6e-01 0:22.0
   13   1040 7.175095840682744e+05 3.9e+00 5.19e-01  2e-01  6e-01 0:23.8
   14   1120 8.949331757795989e+05 4.3e+00 5.72e-01  2e-01  6e-01 0:25.6
   15   1200 4.070832521538258e+05 4.8e+00 5.69e-01  1e-01  6e-01 0:27.5
   16   1280 5.577178147790811e+05 5.5e+00 5.79e-01  1e-01  6e-01 0:29.3
   17   1360 4.352035014696820e+05 6.4e+00 5.55e-01  1e-01  6e-01 0:31.1
   18   1440 5.990224525173972e+05 8.1e+00 5.39e-01  8e-02  6e-01 0:32.9
   19   1520 6.737230961236510e+05 9.0e+00 5.47e-01  7e-02  7e-01 0:34.8
   20   1600 6.965991729066645e+05 1.0e+01 5.71e-01  7e-02  8e-01 0:36.6
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
   21   1680 6.549982072014991e+05 1.3e+01 5.67e-01  6e-02  8e-01 0:38.4
   22   1760 6.154136210757585e+05 1.6e+01 6.17e-01  6e-02  9e-01 0:40.2
   23   1840 5.237943334630505e+05 2.0e+01 6.14e-01  6e-02  9e-01 0:42.1
   24   1920 4.508230247776448e+05 2.1e+01 5.78e-01  5e-02  8e-01 0:43.9
   25   2000 2.968928270409642e+05 2.3e+01 5.19e-01  4e-02  8e-01 0:45.7
   26   2080 2.695570345202658e+05 3.0e+01 5.49e-01  4e-02  8e-01 0:47.6
   27   2160 5.507347958882614e+05 3.2e+01 6.46e-01  4e-02  9e-01 0:49.4
   28   2240 5.710588055798565e+05 3.3e+01 5.57e-01  3e-02  8e-01 0:51.2
   29   2320 6.819989978171875e+05 3.8e+01 5.60e-01  3e-02  8e-01 0:53.0
   30   2400 6.502386763115245e+05 4.1e+01 5.78e-01  3e-02  8e-01 0:54.9
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
   31   2480 2.505747804304289e+05 4.8e+01 5.54e-01  2e-02  8e-01 0:56.7
   32   2560 5.307385505196790e+05 5.1e+01 5.68e-01  2e-02  9e-01 0:58.5
   33   2640 4.687790618079154e+05 5.5e+01 6.02e-01  2e-02  9e-01 1:00.3
   34   2720 3.904404883036885e+05 5.8e+01 6.20e-01  2e-02  9e-01 1:02.2
   35   2800 4.930440231919982e+05 6.3e+01 5.86e-01  2e-02  8e-01 1:04.0
   36   2880 3.319113900955228e+05 6.8e+01 6.06e-01  2e-02  9e-01 1:05.8
   37   2960 1.163125596983890e+05 7.7e+01 6.01e-01  2e-02  9e-01 1:07.6
   38   3040 4.678932193109416e+05 8.4e+01 6.25e-01  2e-02  9e-01 1:09.5
   39   3120 2.559995885166379e+05 9.1e+01 6.07e-01  1e-02  8e-01 1:11.3
   40   3200 5.300835372119243e+05 9.5e+01 6.18e-01  1e-02  9e-01 1:13.1
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
   41   3280 3.369937878427983e+05 1.1e+02 6.46e-01  1e-02  9e-01 1:14.9
   42   3360 6.842531361879989e+05 1.2e+02 6.13e-01  1e-02  9e-01 1:16.8
   43   3440 3.631962753164368e+05 1.4e+02 6.15e-01  1e-02  9e-01 1:18.6
   44   3520 2.610636445113499e+05 1.4e+02 6.00e-01  1e-02  9e-01 1:20.5
   45   3600 3.030424656372985e+05 1.4e+02 6.15e-01  1e-02  9e-01 1:22.3
   46   3680 6.922019517793762e+05 1.5e+02 6.03e-01  1e-02  8e-01 1:24.2
   47   3760 2.847546503623854e+05 1.5e+02 5.49e-01  9e-03  8e-01 1:26.1
   48   3840 2.960272220581251e+05 1.5e+02 5.76e-01  9e-03  9e-01 1:28.0
   49   3920 5.199388713865072e+05 1.5e+02 5.91e-01  9e-03  9e-01 1:29.9
   50   4000 4.997560964751794e+05 1.8e+02 6.77e-01  1e-02  1e+00 1:31.9
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
   51   4080 3.826259211408684e+05 2.1e+02 6.32e-01  8e-03  1e+00 1:33.8
   52   4160 5.841148883404121e+05 2.3e+02 6.34e-01  8e-03  1e+00 1:35.6
   53   4240 4.261370596583476e+05 2.3e+02 6.17e-01  7e-03  1e+00 1:37.5
   54   4320 4.490560179386796e+05 2.5e+02 5.69e-01  7e-03  9e-01 1:39.3
   55   4400 5.669445994916620e+05 2.4e+02 5.76e-01  7e-03  1e+00 1:41.1
   56   4480 5.611108209039365e+05 2.7e+02 6.46e-01  8e-03  1e+00 1:42.9
   57   4560 2.025207969557795e+05 3.0e+02 6.31e-01  7e-03  1e+00 1:44.8
   58   4640 2.826428529063871e+05 3.1e+02 6.72e-01  7e-03  1e+00 1:46.6
   59   4720 1.912578688208579e+05 3.2e+02 5.49e-01  6e-03  1e+00 1:48.4
   60   4800 5.853333260540696e+05 3.4e+02 5.98e-01  6e-03  1e+00 1:50.2
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
   61   4880 1.771689610937646e+05 3.5e+02 5.87e-01  6e-03  1e+00 1:52.0
   62   4960 2.295589135846584e+05 3.8e+02 5.39e-01  5e-03  9e-01 1:53.9
   63   5040 3.123911404618290e+05 3.5e+02 5.04e-01  5e-03  9e-01 1:55.7
   64   5120 5.075201745965198e+05 3.4e+02 4.73e-01  5e-03  8e-01 1:57.5
   65   5200 5.630213868621056e+05 3.4e+02 4.91e-01  5e-03  8e-01 1:59.3
   66   5280 4.634588785397759e+05 3.4e+02 4.84e-01  5e-03  8e-01 2:01.2
   67   5360 4.211808152041871e+05 3.2e+02 4.57e-01  5e-03  7e-01 2:03.0
   68   5440 2.946187445634573e+05 3.1e+02 4.78e-01  5e-03  8e-01 2:04.9
   69   5520 2.397363285800804e+05 3.2e+02 5.51e-01  6e-03  9e-01 2:06.7
   70   5600 3.398574928430300e+05 3.1e+02 4.83e-01  5e-03  8e-01 2:08.6
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
   71   5680 6.491643193275478e+05 3.2e+02 4.87e-01  5e-03  8e-01 2:10.4
   72   5760 4.218654804183816e+05 3.2e+02 4.89e-01  5e-03  7e-01 2:12.2
   73   5840 4.673227467908993e+05 3.1e+02 5.19e-01  5e-03  8e-01 2:14.0
   74   5920 2.989585987072120e+05 3.2e+02 5.79e-01  5e-03  9e-01 2:15.9
   75   6000 2.739230320686988e+05 3.3e+02 5.82e-01  6e-03  9e-01 2:17.7
   76   6080 2.451476504847514e+05 3.6e+02 5.60e-01  5e-03  9e-01 2:19.5
   77   6160 5.030699914655369e+05 3.6e+02 5.17e-01  5e-03  8e-01 2:21.3
   78   6240 5.335846391162868e+05 3.4e+02 4.85e-01  5e-03  7e-01 2:23.2
   79   6320 2.652146287333085e+05 3.5e+02 5.46e-01  5e-03  8e-01 2:25.0
   80   6400 3.118152117204899e+05 3.2e+02 5.16e-01  5e-03  7e-01 2:26.8
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
   81   6480 6.006623873470754e+05 3.1e+02 5.13e-01  5e-03  7e-01 2:28.6
   82   6560 1.466895664464155e+05 3.2e+02 5.15e-01  4e-03  7e-01 2:30.5
   83   6640 2.908963387528905e+05 3.1e+02 5.95e-01  5e-03  8e-01 2:32.3
   84   6720 3.984740241085584e+05 3.1e+02 6.22e-01  5e-03  9e-01 2:34.1
   85   6800 2.473177710566951e+05 3.0e+02 6.56e-01  6e-03  9e-01 2:36.0
   86   6880 2.391990205699504e+05 2.7e+02 6.53e-01  6e-03  1e+00 2:37.8
   87   6960 4.063787807794606e+05 2.8e+02 6.89e-01  7e-03  1e+00 2:39.6
   88   7040 4.131676672187257e+05 2.7e+02 6.64e-01  7e-03  9e-01 2:41.4
   89   7120 4.802385390872010e+05 2.4e+02 6.98e-01  7e-03  1e+00 2:43.3
   90   7200 2.827839549614418e+05 2.5e+02 6.69e-01  7e-03  9e-01 2:45.1
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
   91   7280 5.806860296003490e+05 2.7e+02 6.08e-01  6e-03  8e-01 2:46.9
   92   7360 3.967753417663190e+05 2.8e+02 6.13e-01  6e-03  8e-01 2:48.7
   93   7440 5.627415565297224e+05 2.8e+02 6.04e-01  6e-03  8e-01 2:50.6
   94   7520 4.122935324608524e+05 2.6e+02 5.75e-01  6e-03  7e-01 2:52.4
   95   7600 6.282267307616397e+05 2.6e+02 5.74e-01  5e-03  7e-01 2:54.2
   96   7680 2.830756428404745e+05 2.8e+02 5.35e-01  5e-03  7e-01 2:56.1
   97   7760 4.856154045476310e+05 2.8e+02 5.39e-01  5e-03  6e-01 2:57.9
   98   7840 3.798074094232084e+05 2.7e+02 5.74e-01  5e-03  7e-01 2:59.7
   99   7920 3.568090951815454e+05 2.8e+02 5.54e-01  5e-03  6e-01 3:01.6
  100   8000 2.745326123915605e+05 2.6e+02 5.55e-01  5e-03  6e-01 3:03.4
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  101   8080 2.731510668817981e+05 2.7e+02 4.72e-01  4e-03  5e-01 3:05.2
  102   8160 1.390797083935810e+05 2.5e+02 4.54e-01  4e-03  4e-01 3:07.0
  103   8240 4.045465445698603e+05 2.6e+02 4.34e-01  4e-03  4e-01 3:08.9
  104   8320 4.171375367324962e+05 2.5e+02 4.41e-01  4e-03  4e-01 3:10.7
  105   8400 3.972458983002806e+05 2.6e+02 4.84e-01  4e-03  5e-01 3:12.6
  106   8480 2.631515619453858e+05 2.4e+02 5.30e-01  4e-03  5e-01 3:14.5
  107   8560 3.771269856959054e+05 2.1e+02 4.71e-01  4e-03  4e-01 3:16.4
  108   8640 6.209500492021323e+05 2.2e+02 4.13e-01  3e-03  4e-01 3:18.2
  109   8720 4.087696453932335e+05 2.2e+02 4.36e-01  3e-03  4e-01 3:20.2
  110   8800 3.845884443572421e+05 2.3e+02 4.26e-01  3e-03  4e-01 3:22.1
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  111   8880 3.880485966074669e+05 2.5e+02 4.07e-01  3e-03  4e-01 3:24.2
  112   8960 4.099329715677807e+05 2.5e+02 3.84e-01  3e-03  4e-01 3:26.1
  113   9040 4.324002418509594e+05 2.5e+02 3.74e-01  3e-03  3e-01 3:28.0
  114   9120 4.302941908560772e+05 2.5e+02 3.83e-01  3e-03  3e-01 3:29.9
  115   9200 4.525866475739368e+05 2.3e+02 3.51e-01  2e-03  3e-01 3:31.8
  116   9280 4.255474938055625e+05 2.4e+02 3.64e-01  3e-03  3e-01 3:33.7
  117   9360 3.028901548052295e+05 2.0e+02 3.83e-01  3e-03  3e-01 3:35.6
  118   9440 2.795811013035139e+05 1.8e+02 3.72e-01  3e-03  3e-01 3:37.5
  119   9520 2.763019935779089e+05 1.7e+02 3.36e-01  2e-03  3e-01 3:39.4
  120   9600 2.149276353206184e+05 2.0e+02 3.50e-01  3e-03  3e-01 3:41.2
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  121   9680 3.369073418438304e+05 1.9e+02 3.45e-01  3e-03  3e-01 3:43.1
  122   9760 4.062292193702740e+05 1.7e+02 3.43e-01  3e-03  2e-01 3:44.9
  123   9840 2.983994028120908e+05 1.6e+02 3.36e-01  3e-03  2e-01 3:46.8
  124   9920 2.165551470776061e+05 1.7e+02 3.24e-01  3e-03  2e-01 3:48.7
  125  10000 6.141890459175524e+05 1.6e+02 3.56e-01  3e-03  2e-01 3:50.5
  126  10080 4.064132157073910e+05 1.7e+02 3.58e-01  3e-03  2e-01 3:52.4
  127  10160 4.607205857847249e+05 1.7e+02 3.69e-01  3e-03  2e-01 3:54.2
  128  10240 1.942250175223374e+05 1.8e+02 3.65e-01  3e-03  2e-01 3:56.0
  129  10320 5.137394170615328e+05 1.8e+02 3.63e-01  3e-03  2e-01 3:57.9
  130  10400 5.685290360465306e+05 1.8e+02 3.32e-01  3e-03  2e-01 3:59.7
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  131  10480 1.627916643042467e+05 1.8e+02 3.59e-01  3e-03  2e-01 4:01.6
  132  10560 2.474846040399669e+05 1.9e+02 3.39e-01  3e-03  2e-01 4:03.5
  133  10640 2.228191207625461e+05 2.0e+02 3.35e-01  3e-03  2e-01 4:05.3
  134  10720 6.763253109870419e+05 2.3e+02 3.17e-01  3e-03  2e-01 4:07.1
  135  10800 9.041366155506260e+04 2.3e+02 3.20e-01  3e-03  2e-01 4:09.0
  136  10880 1.124483836090275e+05 2.3e+02 3.49e-01  3e-03  3e-01 4:10.8
  137  10960 1.751548494985363e+05 2.5e+02 3.25e-01  3e-03  2e-01 4:12.6
  138  11040 3.551678820813823e+05 2.4e+02 2.92e-01  2e-03  2e-01 4:14.6
  139  11120 5.834200528483111e+05 2.5e+02 3.03e-01  3e-03  2e-01 4:16.5
  140  11200 2.112032719948488e+05 2.5e+02 3.36e-01  3e-03  2e-01 4:18.3
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  141  11280 2.899876344143704e+05 2.2e+02 3.25e-01  3e-03  2e-01 4:20.1
  142  11360 4.575807415953702e+05 2.1e+02 3.01e-01  3e-03  2e-01 4:22.0
  143  11440 1.630028233323093e+05 2.2e+02 2.80e-01  2e-03  2e-01 4:23.8
  144  11520 4.798366235443731e+05 2.2e+02 2.73e-01  2e-03  2e-01 4:25.6
  145  11600 3.126920448502619e+05 2.4e+02 2.76e-01  2e-03  2e-01 4:27.5
  146  11680 6.577182869426538e+05 2.6e+02 2.54e-01  2e-03  2e-01 4:29.3
  147  11760 3.341615799731213e+05 2.7e+02 2.52e-01  2e-03  2e-01 4:31.1
  148  11840 3.845306065695729e+05 2.5e+02 2.41e-01  2e-03  1e-01 4:33.0
  149  11920 4.010881998011053e+05 2.5e+02 2.71e-01  2e-03  2e-01 4:34.8
  150  12000 6.030191265287810e+05 2.5e+02 2.79e-01  2e-03  2e-01 4:36.6
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  151  12080 1.994925832844276e+05 2.8e+02 2.93e-01  2e-03  2e-01 4:38.5
  152  12160 4.340244846338289e+05 2.9e+02 2.65e-01  2e-03  2e-01 4:40.3
  153  12240 5.096376010539670e+05 3.0e+02 2.49e-01  2e-03  1e-01 4:42.1
  154  12320 1.543750430308851e+05 2.9e+02 2.66e-01  2e-03  2e-01 4:43.9
  155  12400 3.900282895185075e+05 3.1e+02 2.83e-01  2e-03  2e-01 4:45.8
  156  12480 1.956408495251638e+05 3.4e+02 3.09e-01  3e-03  2e-01 4:47.6
  157  12560 2.646914850317949e+05 3.3e+02 3.02e-01  2e-03  2e-01 4:49.5
  158  12640 4.965771377857791e+05 3.4e+02 2.96e-01  2e-03  2e-01 4:51.3
  159  12720 2.933835789763730e+05 3.6e+02 2.97e-01  2e-03  2e-01 4:53.1
  160  12800 4.417797653002021e+05 3.6e+02 2.95e-01  2e-03  2e-01 4:54.9
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  161  12880 4.522362470198605e+05 3.8e+02 3.07e-01  2e-03  2e-01 4:56.8
  162  12960 6.755896156690600e+05 3.9e+02 3.03e-01  2e-03  2e-01 4:58.6
  163  13040 2.479756737224330e+05 4.1e+02 3.06e-01  2e-03  2e-01 5:00.5
  164  13120 2.005781503515428e+05 4.1e+02 2.77e-01  2e-03  2e-01 5:02.3
  165  13200 3.555289097351115e+05 4.2e+02 2.58e-01  2e-03  2e-01 5:04.1
  166  13280 3.283735072711596e+05 4.5e+02 2.48e-01  2e-03  2e-01 5:06.0
  167  13360 3.415238676336622e+05 4.5e+02 2.45e-01  2e-03  2e-01 5:07.8
  168  13440 4.733414699425767e+05 4.6e+02 2.23e-01  1e-03  2e-01 5:09.6
  169  13520 2.914711830526897e+05 4.6e+02 2.28e-01  1e-03  2e-01 5:11.5
  170  13600 2.036134673368848e+05 4.8e+02 2.16e-01  1e-03  2e-01 5:13.3
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  171  13680 4.518265332898386e+05 4.9e+02 1.97e-01  1e-03  1e-01 5:15.1
  172  13760 1.374756067231793e+05 4.9e+02 1.96e-01  1e-03  1e-01 5:17.0
  173  13840 4.461448729051448e+05 5.7e+02 2.03e-01  1e-03  2e-01 5:18.8
  174  13920 5.052115970978581e+05 5.5e+02 2.08e-01  1e-03  2e-01 5:20.6
  175  14000 3.124631917751452e+05 5.4e+02 2.19e-01  1e-03  2e-01 5:22.4
  176  14080 2.669941755664687e+05 5.6e+02 2.33e-01  1e-03  2e-01 5:24.3
  177  14160 1.406120455813452e+05 6.0e+02 2.44e-01  1e-03  2e-01 5:26.1
  178  14240 3.370134564567520e+05 5.8e+02 2.52e-01  1e-03  2e-01 5:27.9
  179  14320 1.704383539197660e+05 6.0e+02 2.60e-01  1e-03  2e-01 5:29.8
  180  14400 4.266192066242368e+05 5.8e+02 2.59e-01  1e-03  2e-01 5:31.6
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  181  14480 3.470209737574017e+05 5.7e+02 2.52e-01  1e-03  2e-01 5:33.4
  182  14560 4.873430119513627e+05 5.7e+02 2.38e-01  1e-03  2e-01 5:35.2
  183  14640 3.037439356723921e+05 5.6e+02 2.66e-01  2e-03  2e-01 5:37.1
  184  14720 4.738142361579736e+05 5.4e+02 2.78e-01  2e-03  2e-01 5:38.9
  185  14800 5.554541896005075e+05 4.8e+02 2.59e-01  2e-03  2e-01 5:40.7
  186  14880 4.711111180138337e+05 5.0e+02 2.39e-01  2e-03  2e-01 5:42.6
  187  14960 1.418206893661525e+05 5.2e+02 2.32e-01  1e-03  2e-01 5:44.4
  188  15040 2.547594074951859e+05 5.3e+02 2.31e-01  1e-03  2e-01 5:46.2
  189  15120 3.694520699169599e+05 4.9e+02 2.27e-01  1e-03  2e-01 5:48.1
  190  15200 2.488426809817426e+05 5.2e+02 2.41e-01  1e-03  2e-01 5:49.9
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  191  15280 3.495462456832656e+05 5.2e+02 2.33e-01  1e-03  2e-01 5:51.7
  192  15360 2.914991439494549e+05 5.5e+02 2.20e-01  1e-03  2e-01 5:53.5
  193  15440 3.052305083590542e+05 5.5e+02 2.15e-01  1e-03  1e-01 5:55.4
  194  15520 2.755214550109520e+05 5.1e+02 2.18e-01  1e-03  1e-01 5:57.2
  195  15600 1.722597265276794e+05 5.0e+02 2.01e-01  1e-03  1e-01 5:59.0
  196  15680 3.092147421938565e+05 5.5e+02 1.93e-01  1e-03  1e-01 6:00.9
  197  15760 1.087282923394148e+05 5.6e+02 1.97e-01  1e-03  1e-01 6:02.7
  198  15840 1.835514061640070e+05 6.6e+02 2.21e-01  1e-03  2e-01 6:04.5
  199  15920 3.356462480922865e+05 6.9e+02 2.30e-01  1e-03  2e-01 6:06.3
  200  16000 4.223773499779609e+05 7.3e+02 2.13e-01  1e-03  2e-01 6:08.2
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  201  16080 2.421619952322700e+05 8.3e+02 2.04e-01  1e-03  2e-01 6:10.0
  202  16160 5.033476498180840e+05 8.7e+02 2.04e-01  1e-03  2e-01 6:11.8
  203  16240 4.810214148976148e+05 8.5e+02 2.01e-01  1e-03  1e-01 6:13.7
  204  16320 2.157586047649338e+05 8.7e+02 2.14e-01  1e-03  1e-01 6:15.5
  205  16400 4.857972528012389e+05 8.4e+02 2.09e-01  1e-03  1e-01 6:17.3
  206  16480 3.784769926996476e+05 8.7e+02 2.31e-01  1e-03  2e-01 6:19.1
  207  16560 4.199327638918217e+05 8.7e+02 2.36e-01  1e-03  2e-01 6:21.0
  208  16640 6.058691267308843e+05 9.4e+02 2.47e-01  1e-03  2e-01 6:22.8
  209  16720 2.606948331512090e+05 1.1e+03 2.43e-01  1e-03  2e-01 6:24.6
  210  16800 4.098288851362380e+05 1.2e+03 2.65e-01  1e-03  2e-01 6:26.5
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  211  16880 2.535688245320961e+05 1.2e+03 2.69e-01  1e-03  2e-01 6:28.3
  212  16960 4.895200204676670e+05 1.1e+03 2.46e-01  9e-04  2e-01 6:30.1
  213  17040 4.224238484018007e+05 1.1e+03 2.40e-01  9e-04  2e-01 6:31.9
  214  17120 2.150036796687361e+05 1.2e+03 2.08e-01  8e-04  1e-01 6:33.8
  215  17200 3.212306009796541e+05 1.2e+03 2.20e-01  8e-04  2e-01 6:35.6
  216  17280 4.797429620070554e+05 1.5e+03 2.06e-01  7e-04  1e-01 6:37.5
  217  17360 1.789031925435951e+05 1.4e+03 2.40e-01  7e-04  2e-01 6:39.3
  218  17440 2.698212892652859e+05 1.4e+03 2.71e-01  8e-04  2e-01 6:41.1
  219  17520 3.992370222097218e+05 1.6e+03 2.64e-01  8e-04  2e-01 6:42.9
  220  17600 2.592453928633541e+05 1.5e+03 2.52e-01  8e-04  2e-01 6:44.8
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  221  17680 2.032463547614080e+05 1.6e+03 2.32e-01  7e-04  2e-01 6:46.6
  222  17760 1.701106701284685e+05 1.7e+03 2.08e-01  6e-04  1e-01 6:48.4
  223  17840 5.079660121465717e+05 1.7e+03 2.00e-01  6e-04  1e-01 6:50.3
  224  17920 1.095770799796029e+05 1.8e+03 1.88e-01  5e-04  1e-01 6:52.1
  225  18000 3.202078406635054e+05 1.8e+03 1.94e-01  6e-04  1e-01 6:53.9
  226  18080 2.410212223477276e+05 1.9e+03 1.81e-01  5e-04  1e-01 6:55.7
  227  18160 3.422653106621394e+05 2.0e+03 1.92e-01  5e-04  1e-01 6:57.6
  228  18240 6.009855977784152e+05 2.0e+03 1.91e-01  6e-04  1e-01 6:59.4
  229  18320 2.988266786086850e+05 1.9e+03 1.76e-01  5e-04  1e-01 7:01.2
  230  18400 1.474017479097760e+05 1.9e+03 1.71e-01  5e-04  1e-01 7:03.1
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  231  18480 4.241556197213951e+05 2.0e+03 1.54e-01  5e-04  9e-02 7:04.9
  232  18560 1.372708790515507e+05 1.9e+03 1.57e-01  5e-04  9e-02 7:06.7
  233  18640 2.747383940314941e+05 2.0e+03 1.70e-01  5e-04  1e-01 7:08.6
  234  18720 1.721861894506034e+05 1.8e+03 1.50e-01  5e-04  8e-02 7:10.4
  235  18800 1.561229276614191e+05 1.7e+03 1.39e-01  4e-04  7e-02 7:12.2
  236  18880 3.995262821680196e+05 1.7e+03 1.46e-01  4e-04  7e-02 7:14.0
  237  18960 3.999831393701014e+05 1.6e+03 1.48e-01  4e-04  7e-02 7:15.9
  238  19040 3.889642169918415e+05 1.5e+03 1.53e-01  4e-04  7e-02 7:17.7
  239  19120 3.321677619840100e+05 1.5e+03 1.53e-01  4e-04  6e-02 7:19.5
  240  19200 1.817698263546127e+05 1.4e+03 1.59e-01  4e-04  7e-02 7:21.4
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  241  19280 2.408908650554652e+05 1.6e+03 1.61e-01  4e-04  7e-02 7:23.2
  242  19360 2.446984021098613e+05 1.8e+03 1.62e-01  4e-04  7e-02 7:25.1
  243  19440 4.550225776243071e+05 1.9e+03 1.68e-01  5e-04  7e-02 7:26.9
  244  19520 4.093781352996993e+05 1.9e+03 1.57e-01  4e-04  6e-02 7:28.8
  245  19600 2.219408465101299e+05 1.8e+03 1.54e-01  4e-04  6e-02 7:30.6
  246  19680 3.325110188720328e+05 1.9e+03 1.61e-01  4e-04  6e-02 7:32.4
  247  19760 4.556824473109522e+05 1.9e+03 1.60e-01  5e-04  6e-02 7:34.3
  248  19840 4.978686535290500e+05 1.9e+03 1.53e-01  4e-04  6e-02 7:36.1
  249  19920 6.081379218722436e+05 1.9e+03 1.50e-01  4e-04  6e-02 7:38.0
  250  20000 2.962071431557625e+05 2.1e+03 1.55e-01  4e-04  6e-02 7:39.8
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  251  20080 4.336263288158929e+05 2.0e+03 1.41e-01  4e-04  5e-02 7:41.7
  252  20160 1.990489659294161e+05 1.9e+03 1.40e-01  4e-04  5e-02 7:43.5
  253  20240 4.115626712083206e+05 2.0e+03 1.65e-01  5e-04  6e-02 7:45.3
  254  20320 3.299425909643659e+05 1.8e+03 1.62e-01  5e-04  6e-02 7:47.2
  255  20400 2.910304928199787e+05 1.9e+03 1.48e-01  4e-04  6e-02 7:49.0
  256  20480 1.787283222251109e+05 2.0e+03 1.36e-01  4e-04  5e-02 7:50.8
  257  20560 6.134446238575511e+05 1.8e+03 1.27e-01  3e-04  5e-02 7:52.6
  258  20640 2.407392014756129e+05 1.8e+03 1.21e-01  3e-04  4e-02 7:54.5
  259  20720 2.353992090414827e+05 1.9e+03 1.32e-01  3e-04  5e-02 7:56.3
  260  20800 1.785412136040612e+05 2.1e+03 1.38e-01  4e-04  5e-02 7:58.1
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  261  20880 1.929783964187899e+05 2.2e+03 1.43e-01  4e-04  5e-02 7:60.0
  262  20960 1.249679496796180e+05 2.3e+03 1.47e-01  4e-04  5e-02 8:01.8
  263  21040 1.533242066262002e+05 2.3e+03 1.53e-01  4e-04  6e-02 8:03.6
  264  21120 4.111891490246241e+05 2.5e+03 1.72e-01  4e-04  7e-02 8:05.5
  265  21200 1.868504336396179e+05 2.5e+03 1.64e-01  4e-04  6e-02 8:07.3
  266  21280 2.359338536647971e+05 2.6e+03 1.61e-01  4e-04  6e-02 8:09.1
  267  21360 4.101580848733981e+05 2.4e+03 1.47e-01  4e-04  5e-02 8:11.0
  268  21440 3.510092590397359e+05 2.4e+03 1.27e-01  3e-04  4e-02 8:12.8
  269  21520 3.073674136645034e+05 2.6e+03 1.15e-01  3e-04  4e-02 8:14.6
  270  21600 4.291105572335332e+05 2.6e+03 1.13e-01  3e-04  4e-02 8:16.4
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  271  21680 4.822400559427417e+05 2.6e+03 1.13e-01  3e-04  4e-02 8:18.3
  272  21760 1.580921666263223e+05 2.8e+03 1.16e-01  3e-04  4e-02 8:20.1
  273  21840 3.106345126775934e+05 2.6e+03 1.24e-01  3e-04  4e-02 8:21.9
  274  21920 4.825969245688906e+05 2.5e+03 1.21e-01  3e-04  4e-02 8:23.7
  275  22000 3.247258644832086e+05 2.6e+03 1.24e-01  3e-04  4e-02 8:25.6
  276  22080 1.966073029455554e+05 2.8e+03 1.35e-01  3e-04  4e-02 8:27.4
  277  22160 3.118776824597737e+05 2.8e+03 1.47e-01  3e-04  4e-02 8:29.2
  278  22240 3.371907648912471e+05 2.8e+03 1.56e-01  3e-04  5e-02 8:31.1
  279  22320 2.738179498596721e+05 3.0e+03 1.61e-01  3e-04  5e-02 8:32.9
  280  22400 3.389769499077631e+05 2.8e+03 1.62e-01  4e-04  5e-02 8:34.7
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  281  22480 3.122305599849589e+05 2.8e+03 1.64e-01  3e-04  5e-02 8:36.6
  282  22560 6.159104754112107e+05 2.5e+03 1.47e-01  3e-04  4e-02 8:38.4
  283  22640 3.465007568380904e+05 2.3e+03 1.49e-01  3e-04  4e-02 8:40.2
  284  22720 3.985700510154743e+05 2.4e+03 1.64e-01  4e-04  4e-02 8:42.1
  285  22800 4.227943827758851e+05 2.1e+03 1.78e-01  4e-04  5e-02 8:43.9
  286  22880 3.887049376226822e+05 2.0e+03 1.76e-01  4e-04  5e-02 8:45.7
  287  22960 1.055626190058610e+05 2.1e+03 1.74e-01  4e-04  5e-02 8:47.5
  288  23040 5.850003683571231e+05 2.2e+03 1.72e-01  3e-04  5e-02 8:49.4
  289  23120 5.125580620964894e+05 2.4e+03 1.77e-01  3e-04  5e-02 8:51.2
  290  23200 2.823239681546185e+05 2.4e+03 1.71e-01  3e-04  5e-02 8:53.0
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  291  23280 4.279242568279767e+05 2.3e+03 1.52e-01  3e-04  4e-02 8:54.9
  292  23360 3.510197440123855e+05 2.5e+03 1.49e-01  3e-04  4e-02 8:56.7
  293  23440 2.487446105113147e+05 2.4e+03 1.62e-01  3e-04  5e-02 8:58.5
  294  23520 2.789376705315134e+05 2.7e+03 1.66e-01  3e-04  5e-02 9:00.3
  295  23600 3.236504253513793e+05 2.8e+03 1.79e-01  3e-04  5e-02 9:02.2
  296  23680 2.659542107064610e+05 2.7e+03 1.65e-01  3e-04  5e-02 9:04.0
  297  23760 4.678161374903045e+05 2.7e+03 1.60e-01  3e-04  5e-02 9:05.8
  298  23840 1.835451921347342e+05 2.8e+03 1.59e-01  3e-04  5e-02 9:07.6
  299  23920 2.578924143248776e+05 3.0e+03 1.68e-01  3e-04  5e-02 9:09.5
  300  24000 7.150175817912375e+04 3.1e+03 1.61e-01  3e-04  4e-02 9:11.3
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  301  24080 2.048580211726793e+05 3.0e+03 1.63e-01  3e-04  4e-02 9:13.1
  302  24160 3.581085422995071e+05 3.1e+03 1.68e-01  3e-04  5e-02 9:15.0
  303  24240 4.443717196821170e+05 3.4e+03 1.55e-01  3e-04  4e-02 9:16.8
  304  24320 4.027290903579830e+05 3.6e+03 1.59e-01  3e-04  4e-02 9:18.6
  305  24400 2.425847221702788e+05 3.5e+03 1.54e-01  3e-04  4e-02 9:20.5
  306  24480 4.673468565301230e+05 3.8e+03 1.68e-01  3e-04  5e-02 9:22.3
  307  24560 2.779299369755411e+05 4.5e+03 1.63e-01  3e-04  5e-02 9:24.1
  308  24640 4.933854276539086e+05 4.3e+03 1.54e-01  3e-04  5e-02 9:26.0
  309  24720 2.390119477266982e+05 4.2e+03 1.47e-01  3e-04  4e-02 9:27.8
  310  24800 3.940570174029972e+05 4.5e+03 1.42e-01  2e-04  4e-02 9:29.6
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  311  24880 6.878600475358664e+05 4.9e+03 1.49e-01  2e-04  5e-02 9:31.5
  312  24960 3.873301995160466e+05 5.0e+03 1.56e-01  3e-04  5e-02 9:33.3
  313  25040 2.635858284260974e+05 5.6e+03 1.54e-01  3e-04  5e-02 9:35.1
  314  25120 4.325539183259477e+05 5.8e+03 1.55e-01  3e-04  5e-02 9:37.0
  315  25200 1.163279034028127e+05 6.1e+03 1.57e-01  3e-04  5e-02 9:38.8
  316  25280 5.656964360710818e+05 6.2e+03 1.58e-01  3e-04  5e-02 9:40.6
  317  25360 1.911352535765436e+05 6.1e+03 1.67e-01  3e-04  4e-02 9:42.5
  318  25440 2.082800337866559e+05 5.8e+03 1.52e-01  3e-04  4e-02 9:44.3
  319  25520 2.034832973070574e+05 5.2e+03 1.47e-01  3e-04  4e-02 9:46.1
  320  25600 3.592294887038206e+05 5.3e+03 1.39e-01  2e-04  4e-02 9:48.0
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  321  25680 3.354615228304304e+05 5.8e+03 1.33e-01  2e-04  3e-02 9:49.8
  322  25760 2.955686373139617e+05 5.9e+03 1.40e-01  3e-04  3e-02 9:51.6
  323  25840 1.979464510003033e+05 6.4e+03 1.23e-01  2e-04  3e-02 9:53.5
  324  25920 2.573378214250650e+05 6.3e+03 1.20e-01  2e-04  3e-02 9:55.3
  325  26000 3.697073025181435e+05 6.4e+03 1.17e-01  2e-04  3e-02 9:57.1
  326  26080 3.686209475628145e+05 6.3e+03 1.14e-01  2e-04  3e-02 9:58.9
  327  26160 1.680255092814707e+05 6.3e+03 1.31e-01  2e-04  3e-02 10:00.8
  328  26240 2.003673299304363e+05 5.7e+03 1.29e-01  2e-04  3e-02 10:02.6
  329  26320 4.039565237730753e+05 5.4e+03 1.30e-01  2e-04  3e-02 10:04.4
  330  26400 4.167737597900786e+05 6.0e+03 1.28e-01  2e-04  3e-02 10:06.3
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  331  26480 1.962859271104860e+05 6.1e+03 1.21e-01  2e-04  3e-02 10:08.1
  332  26560 4.277674308511436e+05 6.4e+03 1.05e-01  2e-04  2e-02 10:09.9
  333  26640 1.099618280341160e+05 6.3e+03 9.90e-02  2e-04  2e-02 10:11.8
  334  26720 1.765940803423911e+05 6.3e+03 1.02e-01  2e-04  2e-02 10:13.6
  335  26800 4.502998918550844e+05 5.9e+03 1.12e-01  2e-04  2e-02 10:15.4
  336  26880 3.323653211366158e+05 5.9e+03 1.15e-01  2e-04  3e-02 10:17.3
  337  26960 3.804092011037405e+05 6.0e+03 1.16e-01  2e-04  3e-02 10:19.1
  338  27040 1.511580370746691e+05 6.7e+03 1.12e-01  2e-04  3e-02 10:20.9
  339  27120 7.852533964762170e+04 7.0e+03 1.07e-01  2e-04  2e-02 10:22.7
  340  27200 2.824107314672819e+05 7.1e+03 1.16e-01  2e-04  2e-02 10:24.6
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  341  27280 2.966252808185979e+05 7.1e+03 1.21e-01  2e-04  2e-02 10:26.4
  342  27360 1.826554115207839e+05 6.4e+03 1.18e-01  2e-04  2e-02 10:28.2
  343  27440 4.600776160311333e+05 6.2e+03 1.24e-01  2e-04  2e-02 10:30.1
  344  27520 3.688278853537621e+05 5.9e+03 1.36e-01  2e-04  3e-02 10:31.9
  345  27600 4.123831129012252e+05 5.7e+03 1.38e-01  3e-04  3e-02 10:33.7
  346  27680 3.766182638582618e+05 5.5e+03 1.51e-01  3e-04  4e-02 10:35.6
  347  27760 2.590068781356228e+05 5.2e+03 1.53e-01  3e-04  4e-02 10:37.4
  348  27840 4.280668541328363e+05 5.1e+03 1.57e-01  3e-04  4e-02 10:39.2
  349  27920 3.557996095934351e+05 5.7e+03 1.53e-01  3e-04  4e-02 10:41.1
  350  28000 3.709471148501500e+05 6.3e+03 1.59e-01  3e-04  4e-02 10:43.0
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  351  28080 4.812799429261802e+05 6.7e+03 1.70e-01  4e-04  5e-02 10:44.8
  352  28160 3.519788135610974e+05 6.9e+03 1.82e-01  4e-04  5e-02 10:46.6
  353  28240 4.369750190456790e+05 7.4e+03 1.81e-01  4e-04  5e-02 10:48.5
  354  28320 3.919019689949603e+05 6.9e+03 1.91e-01  4e-04  5e-02 10:50.3
  355  28400 2.929027266997316e+05 7.4e+03 1.97e-01  5e-04  6e-02 10:52.1
  356  28480 5.008668736677389e+05 7.6e+03 2.00e-01  5e-04  6e-02 10:54.0
  357  28560 3.205129207730124e+05 7.4e+03 1.96e-01  4e-04  6e-02 10:55.8
  358  28640 1.554691015964491e+05 7.8e+03 1.93e-01  4e-04  6e-02 10:57.6
  359  28720 1.784454760355613e+05 8.3e+03 1.78e-01  4e-04  5e-02 10:59.5
  360  28800 2.864228706046575e+05 8.9e+03 1.85e-01  5e-04  6e-02 11:01.3
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  361  28880 2.452790815485439e+05 9.0e+03 1.84e-01  4e-04  6e-02 11:03.2
  362  28960 2.360813611199718e+05 9.5e+03 1.82e-01  4e-04  6e-02 11:05.0
  363  29040 2.619692731592062e+05 1.0e+04 1.79e-01  5e-04  7e-02 11:06.8
  364  29120 4.742404288051384e+05 1.2e+04 1.69e-01  5e-04  7e-02 11:08.6
  365  29200 6.469316380884362e+05 1.3e+04 1.67e-01  4e-04  6e-02 11:10.5
  366  29280 1.810850853573070e+05 1.4e+04 1.68e-01  5e-04  7e-02 11:12.3
  367  29360 2.650550472291533e+05 1.4e+04 1.85e-01  5e-04  8e-02 11:14.2
  368  29440 2.879766080361611e+05 1.6e+04 1.59e-01  4e-04  7e-02 11:16.0
  369  29520 3.039454562080991e+05 1.6e+04 1.50e-01  4e-04  6e-02 11:17.8
  370  29600 2.987763249771442e+05 1.6e+04 1.39e-01  4e-04  6e-02 11:19.6
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  371  29680 4.211116691374737e+05 1.6e+04 1.27e-01  3e-04  5e-02 11:21.5
  372  29760 4.414832711867482e+05 1.7e+04 1.15e-01  3e-04  4e-02 11:23.3
  373  29840 2.299718927574377e+05 1.8e+04 1.16e-01  3e-04  4e-02 11:25.2
  374  29920 2.270247348523468e+05 1.7e+04 1.09e-01  3e-04  4e-02 11:27.0
  375  30000 3.910726927312200e+05 1.6e+04 1.12e-01  3e-04  4e-02 11:28.8
  376  30080 1.084017565421644e+05 1.7e+04 1.22e-01  3e-04  4e-02 11:30.7
  377  30160 2.667095385742468e+05 1.8e+04 1.27e-01  3e-04  5e-02 11:32.5
  378  30240 1.235389802843539e+05 1.9e+04 1.28e-01  4e-04  5e-02 11:34.3
  379  30320 1.809057666572898e+05 2.1e+04 1.27e-01  4e-04  6e-02 11:36.2
  380  30400 4.131077656401533e+05 2.3e+04 1.28e-01  4e-04  5e-02 11:38.0
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  381  30480 2.675197498430882e+05 2.2e+04 1.33e-01  4e-04  6e-02 11:39.9
  382  30560 2.770775747529483e+05 2.4e+04 1.52e-01  5e-04  7e-02 11:41.7
  383  30640 3.916664880021139e+05 2.4e+04 1.53e-01  5e-04  7e-02 11:43.5
  384  30720 2.183310561335137e+05 2.6e+04 1.56e-01  5e-04  7e-02 11:45.4
  385  30800 2.917514776253565e+05 2.6e+04 1.64e-01  5e-04  8e-02 11:47.2
  386  30880 4.192808470573973e+05 2.9e+04 1.54e-01  5e-04  8e-02 11:49.0
  387  30960 8.691398449400981e+04 3.1e+04 1.64e-01  6e-04  8e-02 11:50.9
  388  31040 2.920294958870423e+05 3.3e+04 1.66e-01  5e-04  8e-02 11:52.7
  389  31120 2.548100300817123e+05 3.4e+04 1.85e-01  6e-04  8e-02 11:54.5
  390  31200 3.501440740323898e+05 3.3e+04 1.79e-01  5e-04  8e-02 11:56.4
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  391  31280 1.085961874679935e+05 3.5e+04 1.73e-01  5e-04  7e-02 11:58.2
  392  31360 4.401749141924077e+05 3.3e+04 1.67e-01  5e-04  7e-02 12:00.0
  393  31440 9.305991579879266e+04 3.3e+04 1.71e-01  5e-04  8e-02 12:01.9
  394  31520 5.647450601345352e+04 3.7e+04 1.85e-01  5e-04  8e-02 12:03.7
  395  31600 2.089805500278963e+05 3.7e+04 1.99e-01  5e-04  8e-02 12:05.5
  396  31680 1.996171000195488e+05 3.7e+04 1.75e-01  5e-04  8e-02 12:07.4
  397  31760 2.388873931559220e+05 3.9e+04 1.67e-01  5e-04  7e-02 12:09.2
  398  31840 1.586500622540720e+05 4.1e+04 1.93e-01  6e-04  9e-02 12:11.1
  399  31920 2.389583419364855e+05 4.3e+04 2.36e-01  7e-04  1e-01 12:12.9
  400  32000 4.329872465422895e+05 4.1e+04 2.54e-01  8e-04  1e-01 12:14.7
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  401  32080 2.068380303726364e+05 4.5e+04 2.71e-01  8e-04  1e-01 12:16.6
  402  32160 3.412004409783705e+05 4.4e+04 2.65e-01  8e-04  1e-01 12:18.4
  403  32240 2.693872064409831e+05 4.3e+04 3.05e-01  8e-04  1e-01 12:20.2
  404  32320 1.738993194514262e+05 3.9e+04 3.26e-01  9e-04  1e-01 12:22.1
  405  32400 2.491006114408485e+05 4.1e+04 2.84e-01  8e-04  1e-01 12:23.9
  406  32480 1.762895477326445e+05 4.1e+04 2.89e-01  8e-04  1e-01 12:25.8
  407  32560 3.796914272072202e+05 4.6e+04 2.84e-01  8e-04  1e-01 12:27.7
  408  32640 2.114675863344006e+05 4.6e+04 2.59e-01  7e-04  1e-01 12:29.5
  409  32720 1.749218144217771e+05 4.6e+04 2.52e-01  7e-04  1e-01 12:31.4
  410  32800 3.024413087055433e+05 4.5e+04 2.74e-01  7e-04  1e-01 12:33.2
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  411  32880 4.877842454112544e+05 4.2e+04 2.69e-01  7e-04  1e-01 12:35.1
  412  32960 4.165663256865978e+05 4.3e+04 2.92e-01  7e-04  1e-01 12:36.9
  413  33040 3.907760368485189e+05 3.9e+04 2.83e-01  7e-04  1e-01 12:38.8
  414  33120 3.327705262442791e+05 4.3e+04 2.83e-01  6e-04  9e-02 12:40.6
  415  33200 2.471959743045244e+05 4.1e+04 3.23e-01  8e-04  1e-01 12:42.4
  416  33280 2.513143601694418e+05 4.4e+04 3.49e-01  9e-04  1e-01 12:44.3
  417  33360 4.485802214649814e+05 4.5e+04 3.58e-01  8e-04  1e-01 12:46.2
  418  33440 7.807797202186013e+04 4.6e+04 4.17e-01  9e-04  1e-01 12:48.0
  419  33520 4.503616140639437e+05 4.1e+04 4.01e-01  9e-04  1e-01 12:49.9
  420  33600 1.383537976921835e+05 4.3e+04 3.74e-01  8e-04  1e-01 12:51.7
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  421  33680 2.443406592745304e+05 4.4e+04 3.71e-01  8e-04  1e-01 12:53.6
  422  33760 4.481680169042789e+05 4.6e+04 3.65e-01  8e-04  1e-01 12:55.4
  423  33840 1.344905098192809e+05 4.6e+04 3.29e-01  7e-04  1e-01 12:57.2
  424  33920 1.193464243279705e+05 4.8e+04 3.35e-01  8e-04  1e-01 12:59.1
  425  34000 4.054805270206101e+05 4.8e+04 3.41e-01  8e-04  1e-01 13:00.9
  426  34080 3.351731483460080e+05 5.1e+04 3.20e-01  8e-04  1e-01 13:02.8
  427  34160 2.643681593390380e+05 4.8e+04 3.23e-01  8e-04  1e-01 13:04.6
  428  34240 1.458380241058348e+05 5.1e+04 3.14e-01  8e-04  1e-01 13:06.5
  429  34320 3.032531310678441e+05 5.3e+04 3.31e-01  8e-04  1e-01 13:08.3
  430  34400 2.243481047994588e+05 5.0e+04 3.37e-01  7e-04  1e-01 13:10.2
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  431  34480 2.181747399818132e+05 4.8e+04 3.34e-01  8e-04  1e-01 13:12.0
  432  34560 3.991272272888055e+05 5.0e+04 3.58e-01  8e-04  1e-01 13:13.9
  433  34640 1.245629600025246e+05 5.0e+04 3.70e-01  8e-04  1e-01 13:15.7
  434  34720 2.527305221395790e+05 5.0e+04 3.43e-01  7e-04  1e-01 13:17.6
  435  34800 3.424244107566258e+05 4.8e+04 3.69e-01  8e-04  1e-01 13:19.4
  436  34880 3.483687868169111e+05 4.9e+04 3.85e-01  8e-04  1e-01 13:21.3
  437  34960 2.727059072916437e+05 4.8e+04 3.46e-01  7e-04  1e-01 13:23.1
  438  35040 3.594665873151278e+05 4.9e+04 3.40e-01  7e-04  1e-01 13:24.9
  439  35120 4.941385443197579e+05 5.3e+04 3.49e-01  7e-04  1e-01 13:26.8
  440  35200 2.682309251979704e+05 5.3e+04 3.37e-01  7e-04  1e-01 13:28.6
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  441  35280 8.245242697592237e+04 5.6e+04 3.40e-01  7e-04  1e-01 13:30.5
  442  35360 5.477020950743391e+05 6.1e+04 3.44e-01  7e-04  1e-01 13:32.3
  443  35440 2.836125215307894e+05 5.9e+04 3.56e-01  7e-04  1e-01 13:34.2
  444  35520 1.394111551172188e+05 5.5e+04 3.31e-01  6e-04  9e-02 13:36.0
  445  35600 1.693489691823403e+05 5.5e+04 2.97e-01  5e-04  8e-02 13:37.8
  446  35680 4.248023646815997e+05 5.5e+04 2.87e-01  5e-04  7e-02 13:39.7
  447  35760 3.658408639638016e+05 4.9e+04 3.00e-01  5e-04  8e-02 13:41.5
  448  35840 3.082981176156037e+05 4.9e+04 3.42e-01  6e-04  1e-01 13:43.3
  449  35920 4.588512760719704e+05 6.0e+04 3.23e-01  5e-04  9e-02 13:45.2
  450  36000 3.666363086477173e+05 6.2e+04 3.20e-01  5e-04  8e-02 13:47.0
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  451  36080 3.003786703979125e+05 6.4e+04 3.16e-01  5e-04  8e-02 13:48.8
  452  36160 1.619676623339925e+05 6.0e+04 3.41e-01  5e-04  8e-02 13:50.7
  453  36240 3.107272801827011e+05 6.0e+04 3.47e-01  5e-04  8e-02 13:52.6
  454  36320 2.288627131652387e+05 5.9e+04 3.44e-01  5e-04  8e-02 13:54.5
  455  36400 4.351071496986391e+05 6.1e+04 3.17e-01  4e-04  7e-02 13:56.3
  456  36480 2.656242834396215e+05 5.7e+04 3.14e-01  4e-04  7e-02 13:58.1
  457  36560 3.148825269936358e+05 6.1e+04 2.73e-01  4e-04  6e-02 13:60.0
  458  36640 3.057820309818357e+05 6.1e+04 2.86e-01  4e-04  7e-02 14:01.8
  459  36720 3.167500928246034e+05 6.4e+04 3.13e-01  4e-04  7e-02 14:03.7
  460  36800 2.347507928645731e+05 6.1e+04 3.38e-01  5e-04  8e-02 14:05.5
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  461  36880 1.994192518841733e+05 6.4e+04 3.25e-01  5e-04  8e-02 14:07.4
  462  36960 1.939507776210237e+05 6.9e+04 2.93e-01  4e-04  7e-02 14:09.2
  463  37040 5.802834268187883e+05 7.1e+04 2.93e-01  4e-04  6e-02 14:11.0
  464  37120 2.198870203659660e+05 6.7e+04 2.82e-01  4e-04  6e-02 14:12.9
  465  37200 2.001891565977314e+05 6.8e+04 2.64e-01  3e-04  5e-02 14:14.7
  466  37280 2.545504874689082e+05 6.8e+04 2.53e-01  3e-04  5e-02 14:16.6
  467  37360 2.836251520591929e+05 6.3e+04 2.30e-01  3e-04  4e-02 14:18.4
  468  37440 4.062366365136030e+05 5.6e+04 2.23e-01  3e-04  4e-02 14:20.2
  469  37520 2.019331666460193e+05 5.3e+04 2.12e-01  2e-04  4e-02 14:22.1
  470  37600 2.964048280815321e+05 5.4e+04 2.14e-01  2e-04  4e-02 14:23.9
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  471  37680 2.342618628737253e+05 5.5e+04 2.13e-01  2e-04  3e-02 14:25.8
  472  37760 2.619029255837608e+05 4.9e+04 2.61e-01  3e-04  4e-02 14:27.6
  473  37840 3.205538454484168e+05 4.8e+04 2.97e-01  3e-04  4e-02 14:29.5
  474  37920 5.561358190680464e+05 4.5e+04 3.15e-01  3e-04  5e-02 14:31.3
  475  38000 4.232402384040941e+05 4.7e+04 3.08e-01  3e-04  5e-02 14:33.1
  476  38080 2.508404747617227e+05 4.7e+04 3.27e-01  3e-04  5e-02 14:35.0
  477  38160 4.299975862798023e+05 4.4e+04 3.30e-01  3e-04  5e-02 14:36.8
  478  38240 4.164913008834230e+05 4.8e+04 3.24e-01  3e-04  5e-02 14:38.7
  479  38320 3.348304191643098e+05 4.9e+04 2.90e-01  3e-04  4e-02 14:40.5
  480  38400 1.830881863926781e+05 4.9e+04 2.98e-01  3e-04  4e-02 14:42.3
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  481  38480 4.034137557365547e+05 4.7e+04 3.20e-01  3e-04  5e-02 14:44.2
  482  38560 5.959321034744284e+05 5.5e+04 3.47e-01  3e-04  6e-02 14:46.0
  483  38640 3.299009016399105e+05 6.2e+04 4.00e-01  4e-04  7e-02 14:47.9
  484  38720 2.178517707955661e+05 5.7e+04 3.85e-01  3e-04  6e-02 14:49.7
  485  38800 3.158939691981182e+05 5.7e+04 3.58e-01  3e-04  6e-02 14:51.5
  486  38880 1.089513354485553e+05 5.6e+04 3.45e-01  3e-04  5e-02 14:53.4
  487  38960 3.817067835594862e+05 5.7e+04 3.39e-01  3e-04  5e-02 14:55.2
  488  39040 3.661589726727645e+05 5.9e+04 3.50e-01  3e-04  6e-02 14:57.1
  489  39120 1.502603837959366e+05 6.0e+04 3.75e-01  3e-04  6e-02 14:58.9
  490  39200 3.884340491423417e+05 6.4e+04 3.85e-01  3e-04  6e-02 15:00.8
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  491  39280 3.461808657185764e+05 6.4e+04 4.16e-01  3e-04  6e-02 15:02.8
  492  39360 6.736411303102117e+04 6.2e+04 3.85e-01  3e-04  6e-02 15:04.6
  493  39440 2.582690865515452e+05 6.5e+04 3.99e-01  3e-04  6e-02 15:06.5
  494  39520 2.985379311989335e+05 7.3e+04 4.50e-01  3e-04  7e-02 15:08.4
  495  39600 1.899118579301476e+05 7.1e+04 4.55e-01  3e-04  7e-02 15:10.2
  496  39680 3.189164186964625e+05 6.8e+04 4.61e-01  3e-04  7e-02 15:12.0
  497  39760 2.600451961900954e+05 6.9e+04 4.52e-01  3e-04  7e-02 15:13.9
  498  39840 3.710020992885957e+05 7.0e+04 4.63e-01  3e-04  7e-02 15:15.7
  499  39920 1.412358881168602e+05 7.2e+04 4.55e-01  3e-04  7e-02 15:17.6
  500  40000 1.928243840275437e+05 6.9e+04 4.84e-01  3e-04  7e-02 15:19.4
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  501  40080 3.259266219495906e+05 6.8e+04 4.73e-01  4e-04  8e-02 15:21.2
  502  40160 2.851533468990953e+05 7.6e+04 4.91e-01  4e-04  9e-02 15:23.1
  503  40240 1.672797440594927e+05 8.6e+04 5.00e-01  4e-04  9e-02 15:24.9
  504  40320 3.231884142394990e+05 9.1e+04 4.76e-01  3e-04  8e-02 15:26.7
  505  40400 2.819976172302521e+05 8.2e+04 4.76e-01  3e-04  8e-02 15:28.6
  506  40480 2.231406306415711e+05 8.1e+04 4.59e-01  3e-04  8e-02 15:30.4
  507  40560 2.953435860442707e+05 7.9e+04 4.85e-01  4e-04  8e-02 15:32.3
  508  40640 3.137543635992020e+05 8.6e+04 5.00e-01  4e-04  8e-02 15:34.1
  509  40720 3.231628968019930e+05 8.4e+04 4.60e-01  3e-04  8e-02 15:36.0
  510  40800 2.653082990624747e+05 8.6e+04 4.42e-01  3e-04  8e-02 15:37.8

=== FINAL TRIADIC D360 RESULT ===
Best χ² + reg = 56474.506013453516
Best parameters: [-1.31897057 -0.04963103  0.75232724  0.35006886  0.12823382  0.00741369
 -0.21451311  0.14247704  0.99539674]

"""