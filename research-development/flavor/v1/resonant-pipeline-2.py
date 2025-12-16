import numpy as np
import math

# ============================================================
# FULL HARMONIC RECODE — CLEAN END-TO-END PIPELINE
# ============================================================
# Architecture:
#   1. Harmonic Seed (triadic + divisor + φ-bridge)
#   2. Harmonic Triadic Lattice (analytic, no randomness)
#   3. Divisor/GCD Kernel
#   4. Projection to 3D (geometric eigenmodes)
#   5. Sector Yukawas from lattice kernel
#   6. Majorana kernel + seesaw
#   7. PMNS / CKM extraction
# ============================================================

# ============================================================
# 0. Constants & harmonic math
# ============================================================
N_CYCLE = 360
phi = 97/60      # Fibonacci harmonic bridge

def divisors(n):
    return [k for k in range(1,n+1) if n%k==0]

D360 = divisors(N_CYCLE)

# ============================================================
# 1. Harmonic seed construction (no randomness)
# ============================================================
# Generation → triad mapping
SITE_MAP = {0:0,1:1,2:2,3:0,4:1,5:2,6:0,7:1,8:2}

def harmonic_triad_magnitudes(n0):
    s1 = phi**-1 * np.exp(1j*n0)
    s2 = phi**-2 * np.exp(1j*2*n0)
    s3 = phi**-3 * np.exp(1j*3*n0)
    return [s1,s2,s3]

def harmonic_phase_gradient(n0, N=9):
    base = 2*np.pi/N_CYCLE
    return np.array([n0*base*i for i in range(N)])

def build_harmonic_seed_Y(n0=3, N=9):
    svals = harmonic_triad_magnitudes(n0)
    phases = harmonic_phase_gradient(n0,N)
    Y = np.zeros((N,N),dtype=complex)
    for i in range(N):
        for j in range(N):
            si = svals[SITE_MAP[i]]
            sj = np.conj(svals[SITE_MAP[j]])
            Y[i,j] = si*sj*np.exp(1j*(phases[i]-phases[j]))
    return Y

# ============================================================
# 2. Harmonic Lattice Kernel
# ============================================================

def cyclic_distance(a,b):
    d = abs(a-b)
    return d if d<=N_CYCLE//2 else N_CYCLE-d

def build_harmonic_lattice(n0, positions):
    N = len(positions)
    svals = harmonic_triad_magnitudes(n0)
    phases = harmonic_phase_gradient(n0,N)
    # Proto lattice from pure harmonic triads
    L = np.zeros((N,N),dtype=float)
    for i in range(N):
        for j in range(N):
            d = cyclic_distance(positions[i],positions[j])
            theta = 2*np.pi*d/N_CYCLE
            # triadic cosines form harmonic kernel
            L[i,j] = np.cos(n0*theta)+np.cos(2*n0*theta)+np.cos(3*n0*theta)
    # normalize diagonal to 1
    L /= np.max(np.diag(L))
    return L

# ============================================================
# 3. Divisor/GCD Kernel
# ============================================================

def build_gcd_kernel(L, positions):
    N = len(positions)
    gcd_sum = {}
    gcd_count = {}
    for i in range(N):
        for j in range(i+1,N):
            d = cyclic_distance(positions[i],positions[j])
            if d==0: continue
            g = math.gcd(d, N_CYCLE)
            val = abs(L[i,j])
            gcd_sum[g]=gcd_sum.get(g,0)+val
            gcd_count[g]=gcd_count.get(g,0)+1
    gcd_mean = {g:gcd_sum[g]/gcd_count[g] for g in gcd_sum}
    LG = np.zeros_like(L)
    for i in range(N):
        LG[i,i]=1.0
        for j in range(i+1,N):
            d=cyclic_distance(positions[i],positions[j])
            g=math.gcd(d,N_CYCLE)
            val=gcd_mean.get(g,0)
            LG[i,j]=LG[j,i]=val
    LG/=np.max(np.abs(LG))
    return LG

# ============================================================
# 4. Projection to 3D
# ============================================================

def project_to_3D(K):
    evals, evecs = np.linalg.eigh(K)
    idx = np.argsort(evals)[::-1]
    P = evecs[:,idx[:3]]
    return P

def project_matrix(M,P):
    return P.conj().T @ M @ P

# ============================================================
# 5. Sector Yukawas from harmonic kernels
# ============================================================

def build_sector_Y(K, lam, phase_shift):
    N = K.shape[0]
    core = np.exp(-lam*(1-K))
    PH = np.zeros((N,N),dtype=complex)
    for i in range(N):
        for j in range(N):
            PH[i,j]=np.exp(1j*phase_shift*(i-j))
    return core*PH

# ============================================================
# 6. Majorana + Seesaw
# ============================================================

def seesaw(Y, MR, v=1.0):
    MInv=np.linalg.inv(MR)
    M = -v*v*(Y.T@MInv@Y)
    return 0.5*(M+M.conj().T)

def diag_hermitian(M):
    m,U = np.linalg.eigh(M)
    idx=np.argsort(np.abs(m))
    return m[idx],U[:,idx]

# ============================================================
# 7. PMNS angles
# ============================================================

def pmns_angles(U):
    A=np.abs(U)
    s13=A[0,2]
    t13=np.arcsin(np.clip(s13,0,1))
    c13=np.cos(t13)
    s12=A[0,1]/max(c13,1e-12)
    s23=A[1,2]/max(c13,1e-12)
    return (np.degrees(np.arcsin(np.clip(s12,0,1))),
            np.degrees(np.arcsin(np.clip(s23,0,1))),
            np.degrees(t13))

# ============================================================
# MASTER PIPELINE
# ============================================================

def run_harmonic_pipeline(n0=3, positions=None):
    if positions is None:
        positions = np.array([0,40,80,120,160,200,240,280,320])

    # Step 1: seed Yukawa
    Y9 = build_harmonic_seed_Y(n0, N=len(positions))

    # Step 2: triadic harmonic lattice
    L_tri = build_harmonic_lattice(n0, positions)

    # Step 3: gcd kernel
    L_gcd = build_gcd_kernel(L_tri, positions)

    # Step 4: choose kernel (triadic-only for now)
    K = L_tri

    # Step 5: sector Yukawas
    lam_up, lam_down, lam_e, lam_nu, lam_M = 1.2,1.0,0.9,0.4,1.1
    Yu  = build_sector_Y(K,lam_up,0.0)
    Yd  = build_sector_Y(K,lam_down,np.pi/6)
    Ye  = build_sector_Y(K,lam_e,np.pi/3)
    Ynu = build_sector_Y(K,lam_nu,np.pi/2)
    MR  = build_sector_Y(K,lam_M,0.0) + np.eye(len(positions))

    # Step 6: seesaw & 3D projection
    P = project_to_3D(K)
    Yu3  = project_matrix(Yu,P)
    Ye3  = project_matrix(Ye,P)
    Ynu3 = project_matrix(Ynu,P)
    MR3  = project_matrix(MR,P)

    mnu3 = seesaw(Ynu3,MR3)
    He3  = Ye3.conj().T@Ye3

    me2,Ue = diag_hermitian(He3)
    mnu,Uν = diag_hermitian(mnu3)

    Upmns = Ue.conj().T @ Uν
    th12,th23,th13 = pmns_angles(Upmns)

    ratios = np.abs(mnu)/np.max(np.abs(mnu))

    return {
        "Y9":Y9,
        "K":K,
        "Yu3":Yu3,
        "Ye3":Ye3,
        "Ynu3":Ynu3,
        "MR3":MR3,
        "angles":(th12,th23,th13),
        "mnu_ratios":ratios
    }

# ============================================================
if __name__ == "__main__":
    result = run_harmonic_pipeline()
    print("PMNS angles: ", np.round(result["angles"],2))
    print("Neutrino ratios:", np.round(result["mnu_ratios"],4))