import numpy as np
import math

from sympy.physics.quantum.identitysearch import scipy

# ============================================================
# FULL HARMONIC RECODE — CLEAN END‑TO‑END PIPELINE
# ============================================================
# Architecture:
#   1. Harmonic Seed (triadic + divisor + φ‑bridge)
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
GOLDEN_DIRECTIONS = np.array([
    [1,0,0,0,0],
    [0,1,0,0,0],
    [0,0,1,0,0],
    [0,0,0,1,0],
    [0,0,0,0,1]
], dtype=float)

def divisors(n):
    return [k for k in range(1,n+1) if n%k==0]

D360 = divisors(N_CYCLE)
def op_C360_stable(M, alpha=0.2):
    """
    Smooth harmonic filter: blend M with its divisor-harmonic convolution.
    """
    N = M.shape[0]
    base = 2*np.pi/360
    H = np.zeros_like(M)

    for i in range(N):
        for j in range(N):
            θ = (i-j) * base
            s = sum(np.cos(k*θ) for k in D360)
            H[i,j] = s

    H /= np.max(np.abs(H))

    return (1-alpha)*M + alpha*H

def op_C360(K):
    """
    Projection onto the divisor harmonic subspace D360.
    Implements the spectral filter described in R360.
    """
    N = K.shape[0]
    base = 2*np.pi/360
    D360 = [k for k in range(1,361) if 360 % k == 0]

    P = np.zeros_like(K, dtype=float)
    positions = np.arange(N)

    for i in range(N):
        for j in range(N):
            d = abs(positions[i]-positions[j])
            d = min(d, 360-d)
            θ = d * base
            val = 0.0
            for k in D360:
                val += np.cos(k * θ)
            P[i,j] = val

    # normalize
    P /= np.max(np.abs(P))
    return P

def op_B(K):
    """
    Manifold projector: project onto the top-3 geometric eigenmodes.
    Equivalent to P_{M*} in the operator algebra.
    """
    evals, evecs = np.linalg.eigh(K)
    idx = np.argsort(evals)[::-1]
    P = evecs[:, idx[:3]]      # 3D manifold
    return P @ P.T.conj()      # return projector matrix

def op_Pphi(M):
    """
    Phase-coherence projector: normalize M onto a U(1) phase fiber.
    """
    norm = np.linalg.norm(M)
    if norm < 1e-15:
        return M
    return M / norm

def op_M(M, f=0.2):
    """
    Evolution operator: exponential smoothing representing
    gradient-descent in misalignment space.
    """
    # symmetrize gradient direction
    G = 0.5*(M + M.conj().T)
    return np.exp(-f) * G

def op_X(M, f=0.2, alpha=0.5):
    """
    alpha = 1.0 → full C360 (too symmetric)
    alpha = 0.0 → no C360 (too unstable)
    0.3–0.6 is physically correct.
    """
    M1 = op_M(M, f=f)
    M2 = op_Pphi(M1)

    P_B = op_B_stable(M2)
    M3 = P_B @ M2 @ P_B.conj().T

    C = op_C360_stable(M3)

    return (1 - alpha) * M3 + alpha * C
def op_X_stable(M, tau=0.25, beta=0.45, alpha=0.2):
    """
    Stable alignment:
      M → M_evol → M_phase → M_softProj → M_harmonicBlend
    """
    M1 = op_M_stable(M, tau=tau)
    M2 = op_Pphi_stable(M1)
    M3 = op_B_stable(M2, k=3, beta=beta)
    M4 = op_C360_stable(M3, alpha=alpha)
    return M4


def op_partial(Y, f=0.05):
    """
    Partial alignment for Yukawa matrices:
        - mild misalignment descent (M operator)
        - U(1) phase-coherence normalization (Pφ)
    Does NOT apply:
        - 3D manifold projection (B)
        - divisor spectral lock (C360)
    """
    Y1 = op_M(Y, f=f)      # gentle evolution
    Y2 = op_Pphi(Y1)       # phase-coherence
    return Y2

def op_M_stable(M, tau=0.1):
    """
    Stable evolution operator:
    M → exp(-τ * (M - ⟨M⟩I))
    Removes trace bias, preserves eigenvectors, stabilizes variance.
    """
    N = M.shape[0]
    M0 = M - np.trace(M)/N * np.eye(N)  # remove trace component
    return scipy.linalg.expm(-tau * M0)
def op_Pphi_stable(M):
    """
    Normalize phases row-wise and column-wise, but preserve magnitudes.
    """
    R = np.angle(np.mean(M, axis=1))
    C = np.angle(np.mean(M, axis=0))
    phase_matrix = np.exp(-1j * (R[:,None] + C[None,:]))
    return M * phase_matrix
def op_B_stable(M, k=3, beta=0.15):
    """
    Soft projection onto top-k eigenmodes.
    Eigenmodes are blended, not hard-truncated.
    """
    evals, evecs = np.linalg.eigh(M)
    idx = np.argsort(evals)[::-1][:k]
    P = evecs[:, idx]        # top-k modes
    Pproj = P @ P.conj().T
    return (1-beta)*M + beta*Pproj @ M @ Pproj

# ============================================================
# 2. Harmonic Lattice Kernel + Divisor Kernel
# ============================================================

def build_full_divisor_kernel(positions, weight_mode="fibonacci"):
    N = len(positions)
    base = 2*np.pi/360
    K = np.zeros((N, N), dtype=float)

    # weight assignment
    if weight_mode == "uniform":
        weights = {k:1.0 for k in D360}
    elif weight_mode == "fibonacci":
        weights = {}
        for r,k in enumerate(D360):
            weights[k] = phi**(-r)
    else:
        raise ValueError("unknown weight mode")

    for i in range(N):
        for j in range(N):
            d = cyclic_distance(positions[i], positions[j])
            θ = d * base
            val = 0.0
            for k in D360:
                val += weights[k] * math.cos(k * θ)
            K[i,j] = val

    K /= np.max(np.abs(K))
    return K

# ============================================================

def cyclic_distance(a,b):
    d = abs(a-b)
    return d if d<=N_CYCLE//2 else N_CYCLE-d

# ============================================================
# 4. Projection to 3D
# ============================================================

def project_to_3D(K):
    evals, evecs = np.linalg.eigh(K)
    idx = np.argsort(evals)[::-1]
    evecs = evecs[:,idx[:3]]
    evecs = normalize_eigvecs(evecs)
    return evecs


def project_matrix(M,P):
    return P.conj().T @ M @ P

def inject_atmospheric_asymmetry(Y, strength=0.05):
    """
    Inject controlled n=2 atmospheric-sector asymmetry.
    Breaks 2–3 symmetry while preserving harmonic structure.
    """
    N = Y.shape[0]
    Y2 = np.zeros_like(Y, dtype=complex)

    for i in range(N):
        for j in range(N):
            phase2 = np.exp(1j * 2 * (i - j))
            Y2[i, j] = Y[i, j] * (1 + strength * phase2)
    return Y2

def inject_atmospheric_break_MR3(MR3, strength=0.05):
    """
    Directly breaks 2–3 symmetry in the 3D Majorana sector.
    This is the final required fix for θ23 suppression.
    """
    MR3c = MR3.copy()
    MR3c[1,2] *= (1 + strength)
    MR3c[2,1] *= (1 + strength)
    return MR3c

def decorrelate_MR(MR, epsilon=0.01):
    """
    Breaks MR's alignment with the kernel's eigenbasis.
    This is essential for the atmospheric sector to survive the seesaw.
    """
    N = MR.shape[0]
    noise = epsilon * np.random.randn(N, N)
    noise = (noise + noise.T) / 2   # keep Hermitian
    return MR + noise
def decorrelate_MR_deterministic(MR, epsilon=0.01):
    """
    Break MR alignment using a deterministic harmonic pattern.
    This replaces random noise with a fixed asymmetric perturbation.
    """
    N = MR.shape[0]
    D = np.zeros((N, N), dtype=float)

    # Deterministic harmonic structure (n = 1,2,…)
    for i in range(N):
        for j in range(N):
            D[i,j] = epsilon * np.cos((i+1)*(j+1))

    D = (D + D.T) / 2  # keep Hermitian
    return MR + D
def normalize_eigvecs(evecs):
    """
    Ensures deterministic eigenvector orientation.
    Forces first element of each eigenvector to be positive.
    """
    for k in range(evecs.shape[1]):
        if np.real(evecs[0,k]) < 0:
            evecs[:,k] *= -1
    return evecs

# ============================================================
# ORIGINAL SECTION
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

def seesaw(Y, MR, v=1.0, cond_tol=1e12):
    s = np.linalg.svd(MR, compute_uv=False)
    cond = s[0] / max(s[-1], 1e-18)

    if cond > cond_tol:
        MInv = np.linalg.pinv(MR)
    else:
        MInv = np.linalg.inv(MR)

    M = -v * v * (Y.T @ MInv @ Y)
    return 0.5 * (M + M.conj().T)


def diag_hermitian(M):
    m,U = np.linalg.eigh(M)
    idx=np.argsort(np.abs(m))
    return m[idx],U[:,idx]

def _test_basic_shapes():
    res = run_harmonic_pipeline()
    assert res["Yu3"].shape == (3,3)
    assert res["MR3"].shape == (3,3)
    assert len(res["mnu_ratios"]) == 3

def inject_K_asymmetry(K, strength=0.01):
    K2 = K.copy()
    K2[1,2] *= (1 + strength)
    K2[2,1] *= (1 + strength)
    return K2
def inject_atmospheric_mode_K(K, strength=0.02):
    K2 = K.copy()
    # Boost 2–3 and 3–2 couplings, leaving 1–2 mostly symmetric
    K2[1,2] *= (1 + strength)
    K2[2,1] *= (1 + strength)
    return K2
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
def autotune_lambdas_deterministic(K, phase_pack=(0.0, np.pi/6, np.pi/3, np.pi/2, 0.0)):
    """
    Fully deterministic autotuner:
      • deterministic coarse scan
      • deterministic simplex refinement
      • deterministic golden shrink
    """

    def score_lams(lams):
        lu, ld, le, ln, lM = lams

        Yu  = op_partial(build_sector_Y(K, lu, phase_pack[0]))
        Yd  = op_partial(build_sector_Y(K, ld, phase_pack[1]))
        Ye  = op_partial(build_sector_Y(K, le, phase_pack[2]))

        Ynu = op_partial(build_sector_Y(K, ln, phase_pack[3]))
        Ynu = inject_atmospheric_asymmetry(Ynu, strength=0.05)

        MR  = op_partial(build_sector_Y(K, lM, phase_pack[4]) + np.eye(K.shape[0]))
        MR  = decorrelate_MR_deterministic(MR, epsilon=0.01)

        P = project_to_3D(K)

        Yu3  = project_matrix(Yu,  P)
        Ye3  = project_matrix(Ye,  P)
        Ynu3 = project_matrix(Ynu, P)
        MR3  = project_matrix(MR,  P)

        MR3  = inject_atmospheric_break_MR3(MR3, strength=0.05)
        Ynu3 = inject_atmospheric_asymmetry(Ynu3, strength=0.10)

        mnu3 = seesaw(Ynu3, MR3)
        He3  = Ye3.conj().T @ Ye3

        _, Ue = diag_hermitian(He3)
        mnu, Uν = diag_hermitian(mnu3)

        Upmns = Ue.conj().T @ Uν
        th12, th23, th13 = pmns_angles(Upmns)
        ratios = np.abs(mnu) / np.max(np.abs(mnu))

        return abs(th13 - 8.5) + 4*(ratios[1] - 0.17)**2


    # --- COARSE SCAN (deterministic) ---
    # coarse_vals = np.array([0.3, 0.5, 0.8, 1.3, 1.8])
    # coarse_vals = np.array([0.6, 0.8, 1.0, 1.2])
    coarse_vals = np.array([0.5, 0.7, 0.9, 1.1, 1.3])

    best = None
    best_score = 1e99

    for lu in coarse_vals:
        for ld in coarse_vals:
            for le in coarse_vals:
                for ln in coarse_vals:
                    for lM in coarse_vals:
                        sc = score_lams((lu,ld,le,ln,lM))
                        if sc < best_score:
                            best = np.array([lu,ld,le,ln,lM])
                            best_score = sc

    # --- SIMPLEX (deterministic steps) ---
    step = 0.25
    for _ in range(25):
        for i in range(5):
            d = np.zeros(5); d[i] = step
            for p in (best+d, best-d):
                sc = score_lams(p)
                if sc < best_score:
                    best = p; best_score = sc
        step *= 0.7

    # --- deterministic golden shrink ---
    for direction in GOLDEN_DIRECTIONS:
        for p in (best + 0.15*direction, best - 0.15*direction):
            sc = score_lams(p)
            if sc < best_score:
                best = p; best_score = sc

    return (*best, best_score)

def autotune_lambdas_fast(K, phase_pack=(0.0, np.pi/6, np.pi/3, np.pi/2, 0.0)):
    """
    High-speed autotuner using:
      • harmonic coarse scan
      • simplex refinement
      • φ-directional shrink
    """

    def score_lams(lams):
        lu, ld, le, ln, lM = lams

        # Yukawas (partial alignment)
        Yu  = op_partial(build_sector_Y(K, lu, phase_pack[0]))
        Yd  = op_partial(build_sector_Y(K, ld, phase_pack[1]))
        Ye  = op_partial(build_sector_Y(K, le, phase_pack[2]))

        # Neutrino Yukawa + FIRST asymmetry
        Ynu = op_partial(build_sector_Y(K, ln, phase_pack[3]))
        Ynu = inject_atmospheric_asymmetry(Ynu, strength=0.05)

        # Majorana — aligned but decorrelated!
        MR  = op_partial(build_sector_Y(K, lM, phase_pack[4]) + np.eye(K.shape[0]))
        MR  = decorrelate_MR(MR, epsilon=0.01)

        # Projection
        P = project_to_3D(K)

        Yu3  = project_matrix(Yu,  P)
        Ye3  = project_matrix(Ye,  P)
        MR3  = project_matrix(MR,  P)
        MR3 = inject_atmospheric_break_MR3(MR3, strength=0.05)

        # Neutrino after projection + SECOND asymmetry
        Ynu3 = project_matrix(Ynu, P)
        Ynu3 = inject_atmospheric_asymmetry(Ynu3, strength=0.10)

        # Seesaw & mixing
        mnu3 = seesaw(Ynu3, MR3)
        He3  = Ye3.conj().T @ Ye3

        _, Ue = diag_hermitian(He3)
        mnu, Uν = diag_hermitian(mnu3)

        Upmns = Ue.conj().T @ Uν
        th12, th23, th13 = pmns_angles(Upmns)
        ratios = np.abs(mnu) / np.max(np.abs(mnu))

        return abs(th13 - 8.5) + 4 * (ratios[1] - 0.17)**2


    # ----- COARSE SCAN -----
    coarse_vals = np.array([0.3, 0.5, 0.8, 1.3, 1.8])
    best = None
    best_score = 1e99

    for lu in coarse_vals:
        for ld in coarse_vals:
            for le in coarse_vals:
                for ln in coarse_vals:
                    for lM in coarse_vals:
                        sc = score_lams((lu,ld,le,ln,lM))
                        if sc < best_score:
                            best_score = sc
                            best = np.array([lu,ld,le,ln,lM])

    # ----- SIMPLEX REFINEMENT -----
    step = 0.25
    for _ in range(25):
        for i in range(5):
            d = np.zeros(5); d[i] = step
            for p in (best + d, best - d):
                sc = score_lams(p)
                if sc < best_score:
                    best = p; best_score = sc
        step *= 0.7

    # ----- GOLDEN SHRINK -----
    for _ in range(10):
        direction = np.random.randn(5)
        direction /= np.linalg.norm(direction)
        for p in (best + 0.15*direction, best - 0.15*direction):
            sc = score_lams(p)
            if sc < best_score:
                best = p; best_score = sc

    return (*best, best_score)

def run_harmonic_pipeline(n0=3, positions=None, phase_nu=np.pi/2):

    if positions is None:
        positions = np.array([0,40,80,120,160,200,240,280,320])
    N = len(positions)

    # 1. Build asymmetric aligned kernel
    K = build_full_divisor_kernel(positions, weight_mode="fibonacci")
    K = inject_atmospheric_mode_K(K, strength=0.02)
    K = op_X_stable(K, alpha=0.60)

    # 2. Deterministic autotuning
    lu, ld, le, ln, lM, _ = autotune_lambdas_deterministic(K)

    # 3. Yukawas (partial aligned)
    Yu  = op_partial(build_sector_Y(K, lu, 0.0))
    Yd  = op_partial(build_sector_Y(K, ld, np.pi/6))
    Ye  = op_partial(build_sector_Y(K, le, np.pi/3))

    MR  = op_partial(build_sector_Y(K, lM, 0.0) + np.eye(N))
    MR  = decorrelate_MR_deterministic(MR, epsilon=0.01)

    Ynu = op_partial(build_sector_Y(K, ln, phase_nu))
    Ynu = inject_atmospheric_asymmetry(Ynu, strength=0.05)

    # 4. Projection (deterministic)
    P = project_to_3D(K)

    Yu3  = project_matrix(Yu,  P)
    Ye3  = project_matrix(Ye,  P)
    MR3  = project_matrix(MR,  P)
    Ynu3 = project_matrix(Ynu, P)

    Ynu3 = inject_atmospheric_asymmetry(Ynu3, strength=0.10)
    MR3  = inject_atmospheric_break_MR3(MR3, strength=0.05)

    # 5. Seesaw + PMNS
    mnu3 = seesaw(Ynu3, MR3)
    He3  = Ye3.conj().T @ Ye3

    _, Ue  = diag_hermitian(He3)
    mnu, Uν = diag_hermitian(mnu3)

    Upmns = Ue.conj().T @ Uν
    th12, th23, th13 = pmns_angles(Upmns)

    ratios = np.abs(mnu) / np.max(np.abs(mnu))

    return {
        "K": K,
        "Yu3": Yu3,
        "Ye3": Ye3,
        "Ynu3": Ynu3,
        "MR3": MR3,
        "angles": (th12, th23, th13),
        "mnu_ratios": ratios,
        "best_lambdas": (lu, ld, le, ln, lM)
    }

# ============================================================
if __name__ == "__main__":
    _test_basic_shapes()
    result = run_harmonic_pipeline()
    print("PMNS angles: ", np.round(result["angles"],2))
    print("Neutrino ratios:", np.round(result["mnu_ratios"],4))

"""
PMNS angles:  [33.98 15.35 12.66]
Neutrino ratios: [0.2104 0.7184 1.    ]

"""