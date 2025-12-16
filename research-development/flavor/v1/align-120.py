#!/usr/bin/env python3
# =============================================================================
#   TRIADIC GEOMETRIC ALIGNMENT — FINAL THEORY (FULL IMPLEMENTATION)
#
#   9 continuous parameters:
#      A_u, B_u,
#      A_d, B_d,
#      A_e, B_e,
#      A_nu, B_nu,
#      kappa      ← triadic coherence
#
#   Kernel:   triadic forbidden set {2,4}
#   Phases:   120-cycle phase wheels
#   Geometry: 9-site ring (triadic indexing)
#
#   This version is complete and ready for validation runs.
# =============================================================================

import numpy as np
import cma
from scipy.integrate import solve_ivp

# ------------------------- Physics constants -------------------------
MU_HIGH = 2.0e14   # flavor / seesaw scale
MU_LOW  = 1.0e2    # EW scale
V_HIGGS = 246.0    # Higgs vev

# EW couplings at MU_LOW
g1_EW, g2_EW, g3_EW = 0.36, 0.65, 1.17
lam_H = 0.13

N_SITES = 9
# --------------------------------------------------------------------
#  Experimental targets (14 observables)
# --------------------------------------------------------------------
targets = {
    "m_c/m_t":0.007,   "m_u/m_t":1e-5,
    "m_s/m_b":0.02,    "m_d/m_b":0.001,
    "m_mu/m_tau":0.06, "m_e/m_tau":3e-4,

    "theta12_q":0.226, "theta23_q":0.041, "theta13_q":0.0035,
    "theta12_l":0.59,  "theta23_l":0.84,  "theta13_l":0.15,

    "Delta_m2_21":7.4e-5,
    "Delta_m2_31":2.5e-3,
}
sigmas = {k:0.3*abs(v) for k,v in targets.items()}

# ====================================================================
# 1. TRIADIC KERNEL — Forbidden distances {2,4}
# ====================================================================
def kernel_D120_axiom(kappa, forbid_long_arc=False):
    """
    D120 axiom-projector kernel on 9 sites.

    Sites are at angles theta_i = 40° * i on a 360° circle,
    then folded into a 120° fundamental domain via mod 120.

    Distances on the 120° circle are 0, 40, or 80 degrees.
    """
    # 9 sites on 360°, as in the original D360 setup
    theta_360 = 40.0 * np.arange(N_SITES)
    # fold to 120°
    theta = np.mod(theta_360, 120.0)

    K = np.zeros((N_SITES, N_SITES), dtype=float)

    for i in range(N_SITES):
        for j in range(N_SITES):
            d = abs(theta[i] - theta[j])
            d = min(d, 120.0 - d)  # D120 geodesic distance

            if d < 1e-12:
                K[i, j] = 1.0
            elif abs(d - 40.0) < 1e-12:
                K[i, j] = kappa  # nearest neighbor on the 120° circle
            elif abs(d - 80.0) < 1e-12:
                if forbid_long_arc:
                    K[i, j] = 0.0    # splitter: kill the 80° distance
                else:
                    K[i, j] = kappa**2
            else:
                # In principle we shouldn't reach here with the chosen geometry,
                # but keep it safe:
                K[i, j] = 0.0
    return K

# ====================================================================
# 2. 120-PHASE WHEEL
# ====================================================================
def phase_matrix(A, B):
    """
    φ(i) = A + B*(i mod 3)
    Triadic periodicity → 120° geometric phase unit.
    """
    phi_gen = np.array([A + B*g for g in range(3)])
    phi = np.array([phi_gen[i % 3] for i in range(9)])
    return np.exp(1j * (phi[:,None] - phi[None,:]))

# ====================================================================
# 3. Yukawa builder
# ====================================================================
def build_Yukawa(A, B, kappa, alpha, forbid_long_arc=False):
    phases = phase_matrix(A, B)  # your existing phase wheels
    K = kernel_D120_axiom(kappa, forbid_long_arc=forbid_long_arc)
    Y9 = phases * K
    # normalize largest singular value to 1 and rescale to alpha
    Y9 /= np.linalg.svd(Y9, compute_uv=False)[0]
    Y9 *= alpha
    return Y9

# ====================================================================
# 4. Schur 9→3
# ====================================================================
def schur_9to3(Y9):
    A = Y9[:3,:3]
    B = Y9[:3,3:]
    D = Y9[3:,3:]
    Dinv = np.linalg.pinv(D + 1e-10*np.eye(6))
    return A - B @ Dinv @ B.conj().T

# ====================================================================
# 5. Proto-Majorana
# ====================================================================
def proto_majorana(rng, Lambda_Maj=7e13):
    M = rng.normal(size=(9,9)) + 1j*rng.normal(size=(9,9))
    M = 0.5*(M + M.T.conj())
    M /= np.linalg.svd(M, compute_uv=False)[0]
    M *= Lambda_Maj
    return M

# ====================================================================
# 6. RGE System (1-loop SM Yukawas + kappa operator)
# ====================================================================
def pack(Yu,Yd,Ye,kappa):
    def flat(M): return np.concatenate([M.real.ravel(), M.imag.ravel()])
    return np.concatenate([flat(Yu), flat(Yd), flat(Ye), flat(kappa)])

def unpack(v):
    n=3; N=n*n
    def block(i):
        re = v[i:i+N].reshape((3,3))
        im = v[i+N:i+2*N].reshape((3,3))
        return re + 1j*im
    return block(0), block(2*N), block(4*N), block(6*N)

def beta(t, v, g1,g2,g3,lam):
    Yu,Yd,Ye,kappa = unpack(v)

    # Safety clipping
    for M in (Yu,Yd,Ye):
        np.clip(M, -20, 20, out=M)

    T = np.trace(3*Yu@Yu.conj().T + 3*Yd@Yd.conj().T + Ye@Ye.conj().T).real
    pref = 1/(16*np.pi**2)

    dYu = pref*(Yu*(T-(17/20)*g1**2-(9/4)*g2**2-8*g3**2)
            +1.5*(Yu@Yu.conj().T@Yu - Yd@Yd.conj().T@Yu))

    dYd = pref*(Yd*(T-(1/4)*g1**2-(9/4)*g2**2-8*g3**2)
            +1.5*(Yd@Yd.conj().T@Yd - Yu@Yu.conj().T@Yd))

    dYe = pref*(Ye*(T-(9/4)*g1**2-(9/4)*g2**2) +1.5*(Ye@Ye.conj().T@Ye))

    YeT = Ye@Ye.conj().T
    dkappa = pref*((-3*g2**2 + lam)*kappa + YeT@kappa + kappa@YeT.T)

    return pack(dYu,dYd,dYe,dkappa)

def run_rge(Yu,Yd,Ye,kappa_high):
    sol = solve_ivp(
        beta,
        [np.log(MU_HIGH), np.log(MU_LOW)],
        pack(Yu,Yd,Ye,kappa_high),
        args=(g1_EW,g2_EW,g3_EW,lam_H),
        method="RK45",
        rtol=1e-5,
        atol=1e-8,
        max_step=0.4
    )
    return unpack(sol.y[:,-1])

# ====================================================================
# 7. Observables
# ====================================================================
def get_obs(Yu,Yd,Ye,Mnu):
    def angles(U):
        a = np.abs(U)
        s13 = a[0,2]
        c13 = np.sqrt(max(0.,1-s13**2))
        s12 = a[0,1]/c13 if c13>1e-12 else 0
        s23 = a[1,2]/c13 if c13>1e-12 else 0
        return np.arcsin(min(1,max(0,s12))), \
               np.arcsin(min(1,max(0,s23))), \
               np.arcsin(s13)

    su = np.sort(np.linalg.svd(Yu, compute_uv=False))
    sd = np.sort(np.linalg.svd(Yd, compute_uv=False))
    se = np.sort(np.linalg.svd(Ye, compute_uv=False))

    obs = {
        "m_c/m_t":su[1]/su[2],  "m_u/m_t":su[0]/su[2],
        "m_s/m_b":sd[1]/sd[2],  "m_d/m_b":sd[0]/sd[2],
        "m_mu/m_tau":se[1]/se[2],"m_e/m_tau":se[0]/se[2],
    }

    Uu = np.linalg.svd(Yu)[0]
    Ud = np.linalg.svd(Yd)[0]
    Vckm = Uu.conj().T @ Ud
    obs["theta12_q"],obs["theta23_q"],obs["theta13_q"] = angles(Vckm)

    evals,U_nu = np.linalg.eigh(0.5*(Mnu+Mnu.T))
    mnu = np.sort(np.abs(evals))
    Ue = np.linalg.svd(Ye)[0]
    Upmns = Ue.conj().T @ U_nu
    obs["theta12_l"],obs["theta23_l"],obs["theta13_l"] = angles(Upmns)

    obs["Delta_m2_21"] = mnu[1]**2 - mnu[0]**2
    obs["Delta_m2_31"] = mnu[2]**2 - mnu[0]**2

    return obs

# ====================================================================
# 8. COST FUNCTION — 9 parameters
# ====================================================================
def cost(X, M0):
    A_u,B_u,A_d,B_d,A_e,B_e,A_nu,B_nu,kappa = X

    # Fixed physical normalizations
    alpha_u,alpha_d,alpha_e,alpha_nu = 0.71,0.095,0.082,0.13

    # Build 9×9 Yukawas
    Yu9 = build_Yukawa(A_u,B_u,kappa,alpha_u)
    Yd9 = build_Yukawa(A_d,B_d,kappa,alpha_d)
    Ye9 = build_Yukawa(A_e,B_e,kappa,alpha_e)
    Ynu9= build_Yukawa(A_nu,B_nu,kappa,alpha_nu)

    # Schur 9→3
    Yu_h = schur_9to3(Yu9)
    Yd_h = schur_9to3(Yd9)
    Ye_h = schur_9to3(Ye9)

    # Triadic projector for neutrinos
    P = np.zeros((3,9), dtype=complex)
    for c,sites in enumerate([(0,3,6),(1,4,7),(2,5,8)]):
        P[c,sites] = 1/np.sqrt(3)

    Ynu_eff = P @ Ynu9 @ P.conj().T
    MR = P @ M0 @ P.conj().T

    Mnu_h = -0.5 * V_HIGGS**2 * (Ynu_eff @ np.linalg.pinv(MR + 1e-8*np.eye(3)) @ Ynu_eff.T)
    kappa_h = Mnu_h / V_HIGGS**2

    # RGE evolution
    Yu_l,Yd_l,Ye_l,kappa_l = run_rge(Yu_h,Yd_h,Ye_h,kappa_h)
    Mnu_l = kappa_l * V_HIGGS**2

    # Observables
    obs = get_obs(Yu_l,Yd_l,Ye_l,Mnu_l)

    chi2 = sum(((obs[k]-targets[k])/sigmas[k])**2 for k in targets)
    reg  = 0.05 * np.sum(X**2)   # mild regularization

    return chi2 + reg

# ====================================================================
# 9. RUN OPTIMIZATION
# ====================================================================
rng = np.random.default_rng(777)
M0 = proto_majorana(rng)

x0 = np.concatenate([rng.uniform(-0.4,0.4,8), [0.27]])

es = cma.CMAEvolutionStrategy(
    x0, 0.35,
    {'popsize':80, 'maxiter':100, 'verb_disp':1, 'seed':42}
)

print("Starting Z120 alignment optimization...")
while not es.stop():
    xs = es.ask()
    es.tell(xs, [cost(x, M0) for x in xs])
    es.disp()

print("\n=== FINAL Z120 RESULT ===")
print("Best χ² + reg =", es.best.f)
print("Best parameters:", es.best.x)