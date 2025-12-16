import numpy as np
import math
import cma

# ============================================
# RESONANT ALIGNMENT KERNEL MODEL (RAKM-1)
# 14-PARAMETER PURE RESONANT FLAVOR ENGINE
# ============================================

N = 9
LIGHT = [0,1,2]
HEAVY = [3,4,5,6,7,8]

# ==============
# Fixed geometry
# ==============

def idx(i): return i % 3

# Fixed resonant baselines (triadic)
EXP_UP0   = np.array([4.,2.,0.])
EXP_DOWN0 = np.array([3.,2.,0.])
EXP_LEP0  = np.array([3.,2.,0.])
EXP_NU0   = np.array([1.,0.,0.])

# ============================================================
# Resonant kernels (2 parameters only)
# ============================================================

def build_kernel(gamma):
    """Sector coherence kernel: exp(-γ d) with forbidden d=2."""
    K = np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            if i==j:
                K[i,j] = 1.
            else:
                d = min(abs(i-j), N-abs(i-j))
                if d==2:
                    K[i,j] = 0.
                else:
                    K[i,j] = math.exp(-gamma*d)
    return K

# ============================================================
# Resonant phases (Δ-shifts)
# ============================================================

def phi_gen(A,B):
    """Triad phase generator: φ_g = A + B*g   (3-vector)"""
    return np.array([A, A+B, A+2*B])

def site_phases(phi):
    """Expand 3 generation phases to 9 sites."""
    return np.array([phi[idx(i)] for i in range(N)])

def phase_matrix(phi_site):
    P = np.zeros((N,N),dtype=complex)
    for i in range(N):
        for j in range(N):
            P[i,j] = np.exp(1j*(phi_site[i]-phi_site[j]))
    return P

# ============================================================
# Resonant Yukawa construction
# ============================================================

def site_scales(exp_eff):
    """Triadic magnitude from exponents."""
    a,b,c = exp_eff
    mags = np.array([ (10**-a), (10**-b), (10**-c) ])
    return np.array([ mags[idx(i)] for i in range(N) ])

def proto_Y(exp_eff, A, B):
    phi = phi_gen(A,B)
    phs = site_phases(phi)
    P = phase_matrix(phs)
    s = site_scales(exp_eff)
    Y = np.outer(s,s) * P
    # normalize
    u,sv,vh = np.linalg.svd(Y)
    return Y / np.max(sv)

# ============================================================
# Resonant alignment
# ============================================================

def align(Y, gamma):
    K = build_kernel(gamma)
    return K * Y

# ============================================================
# 9→3 Schur reduction (quarks + charged leptons)
# ============================================================

def schur_9_to_3(Y9):
    A = Y9[np.ix_(LIGHT,LIGHT)]
    B = Y9[np.ix_(LIGHT,HEAVY)]
    C = Y9[np.ix_(HEAVY,LIGHT)]
    D = Y9[np.ix_(HEAVY,HEAVY)]
    Dinv = np.linalg.pinv(D)
    Yeff = A - B @ Dinv @ C
    return Yeff + 1e-12*np.eye(3)

# ============================================================
# Neutrino 9→3 triadic projection + seesaw
# ============================================================

def triadic_proj_9_to_3(Y):
    j = np.arange(N)
    v0 = np.exp(2j*np.pi*0*j/N)/np.sqrt(N)
    v3 = np.exp(2j*np.pi*3*j/N)/np.sqrt(N)
    v6 = np.exp(2j*np.pi*6*j/N)/np.sqrt(N)
    P = np.vstack([v0,v3,v6]).conj()
    return P @ Y @ P.conj().T

def seesaw(M9, Ynu_eff):
    H = M9[np.ix_(HEAVY,HEAVY)]
    # heavy triadic (6→3) projection
    j = np.arange(6)
    B = np.zeros((6,3),dtype=complex)
    for k in range(3):
        B[:,k] = np.exp(2j*np.pi*(k+1)*j/6)/np.sqrt(6)
    MR = B.conj().T @ H @ B
    MR = 0.5*(MR+MR.T) + 1e-12*np.eye(3)
    mD = (174/math.sqrt(2))*Ynu_eff
    return - mD @ np.linalg.inv(MR) @ mD.T

# ============================================================
# Diagonalization + observables
# ============================================================

def diag_dirac(Y):
    U,S,Vh = np.linalg.svd(Y)
    return U, S, Vh.conj().T

def diag_majorana(M):
    H = 0.5*(M+M.conj().T)
    vals,U = np.linalg.eigh(H)
    idx = np.argsort(np.abs(vals))[::-1]
    return U[:,idx], vals[idx]

def mix_angles(U):
    s13 = abs(U[0,2])
    s13 = min(max(s13,0),1)
    t13 = math.asin(s13)
    c13 = math.cos(t13)
    s12 = abs(U[0,1])/c13
    s23 = abs(U[1,2])/c13
    return math.asin(s12), math.asin(s23), t13

# targets
T = {
 'm_c/m_t':7e-3,'m_u/m_t':1e-5,
 'm_s/m_b':2e-2,'m_d/m_b':1e-3,
 'theta12_q':0.226,'theta23_q':0.041,'theta13_q':0.0035,
 'theta12_l':0.59,'theta23_l':0.84,'theta13_l':0.15,
 'Delta21':7.4e-5,'Delta31':2.5e-3
}

σ = {k:0.3*T[k] for k in T}  # wide priors

# ============================================================
# COST FUNCTION
# ============================================================

def cost_func(X, M0):
    # unpack 14 parameters
    A_u,B_u,A_d,B_d,A_nu,B_nu = X[0:6]
    X_u,Y_u,X_d,Y_d,X_nu,Y_nu = X[6:12]
    gamma_q, gamma_l = X[12:14]

    # exponents
    exp_u  = EXP_UP0   + np.array([X_u,Y_u,0])
    exp_d  = EXP_DOWN0 + np.array([X_d,Y_d,0])
    exp_l  = EXP_LEP0  # fixed
    exp_nu = EXP_NU0   + np.array([X_nu,Y_nu,0])

    # proto Yukawas
    Yu0 = proto_Y(exp_u,  A_u, B_u)
    Yd0 = proto_Y(exp_d,  A_d, B_d)
    Ye0 = proto_Y(exp_l,  0,   0  )
    Ynu0= proto_Y(exp_nu, A_nu,B_nu)

    # alignment
    Yu9  = align(Yu0,  gamma_q)
    Yd9  = align(Yd0,  gamma_q)
    Ye9  = align(Ye0,  gamma_l)
    Ynu9 = align(Ynu0, gamma_l)
    M9   = align(M0,   gamma_l)

    # 9→3
    Yu = schur_9_to_3(Yu9)
    Yd = schur_9_to_3(Yd9)
    Ye = schur_9_to_3(Ye9)

    Ynu_eff = triadic_proj_9_to_3(Ynu9)
    Mnu = seesaw(M9, Ynu_eff)

    # observables
    Uu,Su,_ = diag_dirac(Yu)
    Ud,Sd,_ = diag_dirac(Yd)
    Ue,Se,_ = diag_dirac(Ye)
    Unu,mn = diag_majorana(Mnu)

    V = Uu.conj().T @ Ud
    Upm = Ue.conj().T @ Unu

    th12q,th23q,th13q = mix_angles(V)
    th12l,th23l,th13l = mix_angles(Upm)

    mc,mu = Su[1],Su[2]
    ms,md = Sd[1],Sd[2]
    mm,me = Se[1],Se[2]

    mν = np.sort(np.abs(mn))
    dm21 = mν[1]**2 - mν[0]**2
    dm31 = mν[2]**2 - mν[0]**2

    O = {
        'm_c/m_t':mc/Su[0],'m_u/m_t':mu/Su[0],
        'm_s/m_b':ms/Sd[0],'m_d/m_b':md/Sd[0],
        'theta12_q':th12q,'theta23_q':th23q,'theta13_q':th13q,
        'theta12_l':th12l,'theta23_l':th23l,'theta13_l':th13l,
        'Delta21':dm21,'Delta31':dm31
    }

    χ2 = 0.
    for k in T:
        χ2 += ((O[k]-T[k])/σ[k])**2
    return χ2

# ============================================================
# OPTIMIZER
# ============================================================

def run():
    rng = np.random.default_rng(7)
    M0 = rng.normal(size=(N,N)) + 1j*rng.normal(size=(N,N))
    M0 = 0.5*(M0+M0.T)

    X0 = np.zeros(14)
    X0[:] = rng.normal(scale=0.3,size=14)

    es = cma.CMAEvolutionStrategy(X0,0.3,{'popsize':20})
    while not es.stop():
        xs = es.ask()
        cs = [cost_func(x,M0) for x in xs]
        es.tell(xs,cs)
        es.disp()

    xbest = es.best.x
    print("\nBEST 14-PARAMETER RESONANT FIT:")
    print(xbest)
    print("cost =",es.best.f)

if __name__=="__main__":
    run()