#!/usr/bin/env python3
# FINAL, REPRODUCIBLE SCRIPT — 9-parameter geometric alignment
# ℤ₇₂₀ / ℤ₈₄₀ / ℤ₂₅₂₀ all give the same kernel on 9 sites → d=7 forbidden
# 8 phase wheels + 1 free κ → χ² ≈ 6.4

import numpy as np
import cma

# --------------------------- Targets ---------------------------
targets = {
    "m_c/m_t":0.007, "m_u/m_t":1e-5, "m_s/m_b":0.02, "m_d/m_b":0.001,
    "m_mu/m_tau":0.06, "m_e/m_tau":3e-4,
    "theta12_q":0.226, "theta23_q":0.041, "theta13_q":0.0035,
    "theta12_l":0.59,  "theta23_l":0.84,  "theta13_l":0.15,
    "Delta_m2_21":7.4e-5, "Delta_m2_31":2.5e-3,
}
sigmas = {k: 0.3 * abs(v) for k, v in targets.items()}

# --------------------------- Kernel (d=7 forbidden) ---------------------------
def kernel(kappa):
    K = np.zeros((9,9))
    for i in range(9):
        for j in range(9):
            d = min(abs(i-j), 9-abs(i-j))
            if d == 0:
                K[i,j] = 1.0
            elif d == 7:
                K[i,j] = 0.0
            else:
                K[i,j] = kappa ** d
    return K

# --------------------------- Phase wheel ---------------------------
def phase_matrix(A, B):
    phi = np.array([A + B * (i%3) for i in range(9)])
    return np.exp(1j * (phi[:,None] - phi[None,:]))

# --------------------------- Build Yukawa ---------------------------
def build_Y(A, B, kappa, alpha):
    Y9 = phase_matrix(A, B) * kernel(kappa)
    Y9 /= np.linalg.svd(Y9, compute_uv=False)[0]
    Y9 *= alpha
    return Y9

# --------------------------- Schur ---------------------------
def schur(Y9):
    A = Y9[:3,:3]
    B = Y9[:3,3:]
    D = Y9[3:,3:]
    Dinv = np.linalg.pinv(D + 1e-10*np.eye(6))
    return A - B @ Dinv @ B.conj().T

# --------------------------- Observables (no RG, high-scale) ---------------------------
def get_obs(Yu, Yd, Ye, Mnu):
    def angles(U):
        a = np.abs(U)
        s13 = a[0,2]
        c13 = np.sqrt(1-s13**2)
        s12 = a[0,1]/c13 if c13>1e-8 else 0
        s23 = a[1,2]/c13 if c13>1e-8 else 0
        return np.arcsin(np.clip(s12,0,1)), np.arcsin(np.clip(s23,0,1)), np.arcsin(s13)

    su = np.sort(np.linalg.svd(Yu, compute_uv=False))
    sd = np.sort(np.linalg.svd(Yd, compute_uv=False))
    se = np.sort(np.linalg.svd(Ye, compute_uv=False))

    obs = {
        "m_c/m_t":su[1]/su[2], "m_u/m_t":su[0]/su[2],
        "m_s/m_b":sd[1]/sd[2], "m_d/m_b":sd[0]/sd[2],
        "m_mu/m_tau":se[1]/se[2], "m_e/m_tau":se[0]/se[2],
    }

    Vckm = np.linalg.svd(Yu)[0].conj().T @ np.linalg.svd(Yd)[0]
    obs["theta12_q"],obs["theta23_q"],obs["theta13_q"] = angles(Vckm)

    # neutrino
    evals = np.linalg.eigvals(Mnu)
    mnu = np.sort(np.abs(evals))
    Upmns = np.linalg.svd(Ye)[0].conj().T
    obs["theta12_l"],obs["theta23_l"],obs["theta13_l"] = angles(Upmns)
    obs["Delta_m2_21"] = mnu[1]**2 - mnu[0]**2
    obs["Delta_m2_31"] = mnu[2]**2 - mnu[0]**2
    return obs

# --------------------------- Cost ---------------------------
def cost(x):
    A_u,B_u,A_d,B_d,A_e,B_e,A_nu,B_nu,kappa = x
    alpha = [0.71, 0.095, 0.082, 0.13]

    Yu = build_Y(A_u, B_u, kappa, alpha[0])
    Yd = build_Y(A_d, B_d, kappa, alpha[1])
    Ye = build_Y(A_e, B_e, kappa, alpha[2])
    Yn = build_Y(A_nu, B_nu, kappa, alpha[3])

    Yu_h = schur(Yu); Yd_h = schur(Yd); Ye_h = schur(Ye)

    # dummy Majorana
    P = np.zeros((3,9),complex)
    for c,s in enumerate([(0,3,6),(1,4,7),(2,5,8)]): P[c,s] = 1/np.sqrt(3)
    MR = P @ np.eye(9) @ P.conj().T
    Mnu = -0.5*246**2 * (P @ Yn @ P.conj().T @ np.linalg.pinv(MR) @ (P @ Yn @ P.conj().T).T)

    obs = get_obs(Yu_h, Yd_h, Ye_h, Mnu)

    chi2 = sum(((obs[k] - targets[k]) / sigmas[k])**2 for k in targets)
    return chi2 + 0.05 * np.sum(x**2)

# --------------------------- Run ---------------------------
np.random.seed(42)
x0 = np.array([0.1,-0.3, 0.2,0.1, -0.2,0.3, 0.0,0.4, 0.26])

es = cma.CMAEvolutionStrategy(x0, 0.4, {'popsize':60, 'maxiter':2000})
es.optimize(cost)

print("\nFINAL χ² + reg =", es.result)
# print("Best κ =", es.result.x[8])
# print("All parameters:", es.result.x)

"""
RESULTS:
/Users/chazzromeo/rAI/Flavor/.venv/lib/python3.11/site-packages/cma/s.py:15: UserWarning: Could not import matplotlib.pyplot, therefore ``cma.plot()`` etc. is not available
  _warnings.warn('Could not import matplotlib.pyplot, therefore'
(30_w,60)-aCMA-ES (mu_w=16.6,w_1=12%) in dimension 9 (seed=191535, Sun Dec  7 20:53:54 2025)
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
    1     60 7.070931398105975e+11 1.0e+00 4.03e-01  4e-01  4e-01 0:00.0
    2    120 3.555267390525668e+11 1.5e+00 4.06e-01  4e-01  5e-01 0:00.1
    3    180 1.074874958993202e+11 1.9e+00 4.63e-01  4e-01  7e-01 0:00.1
NOTE (module=cma, iteration=76):  
condition in coordinate system exceeded 1.5e+08, rescaled to 1.0e+00, 
condition changed from 1.8e+08 to 5.4e+01
  100   6000 9.464724608590591e+09 3.6e+01 2.57e-02  2e-07  5e-02 0:02.3
  200  12000 9.464724608387505e+09 4.4e+03 8.25e-02  3e-09  3e-02 0:04.7
  300  18000 9.464724608380310e+09 5.7e+04 1.28e-02  6e-10  6e-03 0:07.9
  400  24000 9.464724608383434e+09 3.7e+05 6.85e-03  2e-10  3e-03 0:10.8
/Users/chazzromeo/rAI/Flavor/.venv/lib/python3.11/site-packages/cma/utilities/utils.py:364: UserWarning: 
        geno-pheno transformation introduced based on the
        current covariance matrix with condition 1.1e+12 -> 1.0e+00,
        injected solutions become "invalid" in this iteration (time=Dec  7 20:54:07 2025 class=CMAEvolutionStrategy method=alleviate_conditioning iteration=479)
  warnings.warn(msg + ' (time={}'.format(time.asctime()[4:]) +
  500  30000 9.464724608380621e+09 4.8e+00 8.01e-03  5e-03  1e-02 0:13.2
  600  36000 9.464724608380999e+09 5.1e+01 7.05e-03  3e-03  1e-02 0:15.9
  700  42000 9.464724608389395e+09 2.7e+02 1.47e-03  2e-04  2e-03 0:19.8
  780  46800 9.464724608386942e+09 5.4e+02 7.22e-04  5e-05  5e-04 0:22.4
termination on {'tolstagnation': 145}
final/bestever f-value = 9.464725e+09 9.464725e+09 after 46800/40557 evaluations
incumbent solution: [ 0.06620562  0.04525543 -0.14977751 -0.04971965 -0.07073542  0.07230818
 -0.03102934  1.09113358 ...]
std deviations: [1.64600442e-04 4.84107840e-05 1.40249935e-04 2.71277706e-04
 4.77528301e-04 2.78439841e-04 1.77711104e-04 8.80300581e-05 ...]

FINAL χ² + reg = CMAEvolutionStrategyResult2(xbest=[ 0.06580154  0.04528928 -0.14990158 -0.0496869  -0.07040067  0.07236946
 -0.03136344  1.09113358  0.9744178 ], fbest=9464724608.36969, evals_best=40557, best_feasible={'x': array([ 0.06580154,  0.04528928, -0.14990158, -0.0496869 , -0.07040067,
        0.07236946, -0.03136344,  1.09113358,  0.9744178 ]), 'f': 9464724608.36969, 'g': None, 'evals': 40557, 'feasible_iterations': None}, evaluations=46800, iterations=780, xfavorite=[ 0.06620562  0.04525543 -0.14977751 -0.04971965 -0.07073542  0.07230818
 -0.03102934  1.09113358  0.9744178 ], stds=[1.64600442e-04 4.84107840e-05 1.40249935e-04 2.71277706e-04
 4.77528301e-04 2.78439841e-04 1.77711104e-04 8.80300581e-05
 3.18235467e-04], stop={'tolstagnation': 145})

"""