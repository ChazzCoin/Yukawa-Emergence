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
    ang = 2*np.pi*n0/360
    s1 = phi**-1 * np.exp(1j*ang)
    s2 = phi**-2 * np.exp(1j*2*ang)
    s3 = phi**-3 * np.exp(1j*3*ang)
    return [s1, s2, s3]


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

def build_full_divisor_kernel(positions, weight_mode="fibonacci"):
    N = len(positions)
    K = np.zeros((N, N), dtype=float)
    base = 2*np.pi/360

    # weight assignment
    if weight_mode == "uniform":
        weights = {k:1.0 for k in D360}
    elif weight_mode == "fibonacci":
        # Fibonacci scaling: w_k = φ^(-rank(k))
        # rank = position of k in sorted D360
        weights = {}
        for r,k in enumerate(D360):
            weights[k] = phi**(-r)
    else:
        raise ValueError("unknown weight mode")

    # build divisor-harmonic kernel
    for i in range(N):
        for j in range(N):
            d = cyclic_distance(positions[i], positions[j])
            θ = d * base
            val = 0.0
            for k in D360:
                val += weights[k] * math.cos(k * θ)
            K[i,j] = val

    # normalize
    K /= np.max(np.abs(K))
    return K
