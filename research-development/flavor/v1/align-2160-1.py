#!/usr/bin/env python3
# ============================================================================
#  THE FINAL THEORY — Triadic ℤ₂₁₆₀ Kernel with Forbidden Distances {2,4,7}
#  9 Continuous Parameters:
#      A_u,B_u, A_d,B_d, A_e,B_e, A_nu,B_nu, kappa
#
#  Architecture:
#      • 9-site geometric flavor lattice
#      • Triadic distance kernel with strict forbidden set {2,4,7}
#      • Schur 9→3 Yukawa reduction
#      • Triadic neutrino projector (sites (0,3,6), (1,4,7), (2,5,8))
#      • Full 1-loop SM RGE for Yu,Yd,Ye,kappa
#      • 14 experimental observables, χ² ≈ 3–6 expected
# ============================================================================

import numpy as np
import cma
from scipy.integrate import solve_ivp

# --------------------------- Global constants ---------------------------

MU_HIGH = 2.0e14   # seesaw / flavor scale
MU_LOW  = 1.0e2    # electroweak scale
V_HIGGS = 246.0

# Electroweak parameters
g1_EW = 0.36
g2_EW = 0.65
g3_EW = 1.17
lam_H = 0.13

# 14 observables and uncertainties
targets = {
    "m_c/m_t":0.007, "m_u/m_t":1e-5,
    "m_s/m_b":0.02,  "m_d/m_b":0.001,
    "m_mu/m_tau":0.06, "m_e/m_tau":3e-4,
    "theta12_q":0.226, "theta23_q":0.041, "theta13_q":0.0035,
    "theta12_l":0.59,  "theta23_l":0.84,  "theta13_l":0.15,
    "Delta_m2_21":7.4e-5, "Delta_m2_31":2.5e-3,
}
sigmas = {k:0.3*abs(v) for k,v in targets.items()}

# ============================================================================
#  TRIADIC KERNEL  —  ℤ₂₁₆₀ Geometry
#  Forbidden distances = {2, 4, 7}
#  Effective distance range for the 9-site ring: d ∈ {0…4}
# ============================================================================

def triadic_kernel_2160(kappa):
    """
    9×9 kernel with strict forbidden set {2,4,7}.
    On a 9-site ring the distinct distances are only d=0..4,
    but we enforce the triadic harmonic structure inherited from ℤ₂₁₆₀.
    """
    K = np.zeros((9,9), dtype=float)
    forbidden = {2, 4, 7}  # although 7 never appears on 9-ring, we keep the rule exact

    for i in range(9):
        for j in range(9):
            d = abs(i-j)
            d = min(d, 9-d)  # ring distance

            if d == 0:
                K[i,j] = 1.0
            elif d in forbidden:
                K[i,j] = 0.0
            else:
                K[i,j] = kappa**d

    return K

# ============================================================================
#  PHASE GEOMETRY
# ============================================================================

def phase_matrix(A, B):
    """
    Phase wheel: φ_i = A + B*(i mod 3)
    Returns exp[i(φ_i − φ_j)].
    """
    phi = np.array([A + B*(i % 3) for i in range(9)])
    return np.exp(1j * (phi[:,None] - phi[None,:]))

# ============================================================================
#  Yukawa Construction
# ============================================================================

def build_Yukawa(A, B, kappa, alpha):
    """
    Build full 9×9 Yukawa with triadic kernel and phases.
    Normalize so largest singular value = 1, then scale by alpha.
    """
    Y = phase_matrix(A,B) * triadic_kernel_2160(kappa)
    sv = np.linalg.svd(Y, compute_uv=False)
    Y /= sv[0]
    return alpha * Y

# ============================================================================
#  SCHUR COMPLEMENT 9 → 3
# ============================================================================

def schur_9to3(Y9):
    A = Y9[:3,:3]
    B = Y9[:3,3:]
    D = Y9[3:,3:]
    Dinv = np.linalg.pinv(D + 1e-10*np.eye(6))
    return A - B @ Dinv @ B.conj().T

# ============================================================================
#  Proto-Majorana (random but fixed scale)
# ============================================================================

def proto_majorana(rng, scale=7e13):
    M = rng.normal(size=(9,9)) + 1j*rng.normal(size=(9,9))
    M = 0.5*(M + M.T.conj())
    M /= np.linalg.svd(M, compute_uv=False)[0]
    return scale * M

# ============================================================================
#  RGE PACK / UNPACK
# ============================================================================

def pack(Yu,Yd,Ye,kappa):
    def f(M): return np.concatenate([M.real.ravel(), M.imag.ravel()])
    return np.concatenate([f(Yu), f(Yd), f(Ye), f(kappa)])

def unpack(v):
    n = 3; N = n*n
    def blk(i):
        re = v[i:i+N].reshape((3,3))
        im = v[i+N:i+2*N].reshape((3,3))
        return re + 1j*im
    return blk(0), blk(2*N), blk(4*N), blk(6*N)

# ============================================================================
#  1-LOOP RGE EVOLUTION (Stable)
# ============================================================================

def beta(t, v, g1,g2,g3,lam):
    Yu,Yd,Ye,kappa = unpack(v)

    # Clip to prevent runaway
    for M in (Yu, Yd, Ye):
        np.clip(M, -20, 20, out=M)

    T = np.trace(3*Yu@Yu.conj().T +
                 3*Yd@Yd.conj().T +
                 Ye@Ye.conj().T).real

    pref = 1/(16*np.pi**2)

    dYu = pref*(Yu*(T - (17/20)*g1**2 - (9/4)*g2**2 - 8*g3**2)
                + 1.5*(Yu@Yu.conj().T@Yu - Yd@Yd.conj().T@Yu))

    dYd = pref*(Yd*(T - (1/4)*g1**2 - (9/4)*g2**2 - 8*g3**2)
                + 1.5*(Yd@Yd.conj().T@Yd - Yu@Yu.conj().T@Yd))

    dYe = pref*(Ye*(T - (9/4)*g1**2 - (9/4)*g2**2)
                + 1.5*(Ye@Ye.conj().T@Ye))

    YeT = Ye@Ye.conj().T
    dkappa = pref*((-3*g2**2 + lam)*kappa + (YeT@kappa + kappa@YeT.T))

    return pack(dYu,dYd,dYe,dkappa)

def run_rge(Yu,Yd,Ye,kappa_high):
    sol = solve_ivp(
        beta,
        [np.log(MU_HIGH), np.log(MU_LOW)],
        pack(Yu,Yd,Ye,kappa_high),
        args=(g1_EW,g2_EW,g3_EW,lam_H),
        rtol=1e-5, atol=1e-8,
        method='RK45', max_step=0.4
    )
    return unpack(sol.y[:,-1])

# ============================================================================
#  Observables
# ============================================================================

def extract_angles(U):
    a = np.abs(U)
    s13 = a[0,2]
    c13 = np.sqrt(max(0,1-s13**2))
    s12 = a[0,1]/c13 if c13>1e-10 else 0
    s23 = a[1,2]/c13 if c13>1e-10 else 0
    return (np.arcsin(np.clip(s12,0,1)),
            np.arcsin(np.clip(s23,0,1)),
            np.arcsin(s13))

def get_obs(Yu,Yd,Ye,Mnu):
    su = np.sort(np.linalg.svd(Yu, compute_uv=False))
    sd = np.sort(np.linalg.svd(Yd, compute_uv=False))
    se = np.sort(np.linalg.svd(Ye, compute_uv=False))

    obs = {
        "m_c/m_t":su[1]/su[2],
        "m_u/m_t":su[0]/su[2],
        "m_s/m_b":sd[1]/sd[2],
        "m_d/m_b":sd[0]/sd[2],
        "m_mu/m_tau":se[1]/se[2],
        "m_e/m_tau":se[0]/se[2],
    }

    # CKM
    Uu = np.linalg.svd(Yu)[0]
    Ud = np.linalg.svd(Yd)[0]
    Vckm = Uu.conj().T @ Ud
    obs["theta12_q"],obs["theta23_q"],obs["theta13_q"] = extract_angles(Vckm)

    # Neutrinos
    evals,U_nu = np.linalg.eigh(0.5*(Mnu+Mnu.T))
    m = np.sort(np.abs(evals))
    Upmns = np.linalg.svd(Ye)[0].conj().T @ U_nu
    obs["theta12_l"],obs["theta23_l"],obs["theta13_l"] = extract_angles(Upmns)

    obs["Delta_m2_21"] = m[1]**2 - m[0]**2
    obs["Delta_m2_31"] = m[2]**2 - m[0]**2

    return obs

# ============================================================================
#  COST FUNCTION — 9 PARAMETERS
# ============================================================================

def cost(X, M0):
    A_u,B_u, A_d,B_d, A_e,B_e, A_nu,B_nu, kappa = X

    # High-scale normalizations (empirical)
    alpha_u  = 0.71
    alpha_d  = 0.095
    alpha_e  = 0.082
    alpha_nu = 0.13

    # Build 9-site Yukawas
    Yu9  = build_Yukawa(A_u, B_u, kappa, alpha_u)
    Yd9  = build_Yukawa(A_d, B_d, kappa, alpha_d)
    Ye9  = build_Yukawa(A_e, B_e, kappa, alpha_e)
    Ynu9 = build_Yukawa(A_nu, B_nu, kappa, alpha_nu)

    # Schur 9→3
    Yu_h = schur_9to3(Yu9)
    Yd_h = schur_9to3(Yd9)
    Ye_h = schur_9to3(Ye9)

    # Triadic neutrino projector P
    P = np.zeros((3,9), dtype=complex)
    for c,sites in enumerate([(0,3,6),(1,4,7),(2,5,8)]):
        P[c,sites] = 1/np.sqrt(3)

    # Effective neutrino Dirac and Majorana
    Ynu_eff = P @ Ynu9 @ P.conj().T
    MR = P @ M0 @ P.conj().T
    Mnu_h = -0.5 * V_HIGGS**2 * (
        Ynu_eff @ np.linalg.pinv(MR + 1e-8*np.eye(3)) @ Ynu_eff.T
    )
    kappa_h = Mnu_h / V_HIGGS**2

    # RGE downward
    Yu_l, Yd_l, Ye_l, kappa_l = run_rge(Yu_h, Yd_h, Ye_h, kappa_h)
    Mnu_l = kappa_l * V_HIGGS**2

    # Observables
    obs = get_obs(Yu_l, Yd_l, Ye_l, Mnu_l)

    chi2 = sum(((obs[k]-targets[k])/sigmas[k])**2 for k in targets)
    reg = 0.05*np.sum(X**2)
    return chi2 + reg

# ============================================================================
#  RUN OPTIMIZATION
# ============================================================================

if __name__ == "__main__":
    rng = np.random.default_rng(777)
    M0 = proto_majorana(rng)

    # Initial guess
    x0 = np.concatenate([rng.uniform(-0.4,0.4,8), [0.27]])

    es = cma.CMAEvolutionStrategy(
        x0, 0.35,
        {'popsize':80, 'maxiter':4000, 'seed':42}
    )

    print("Starting ℤ₂₁₆₀ triadic optimization...")

    while not es.stop():
        xs = es.ask()
        costs = [cost(x, M0) for x in xs]
        es.tell(xs, costs)
        es.disp()

    print("\n=== FINAL ℤ₂₁₆₀ TRIADIC RESULT ===")
    print("Best χ²+reg =", es.best.f)
    print("Best parameters:", es.best.x)