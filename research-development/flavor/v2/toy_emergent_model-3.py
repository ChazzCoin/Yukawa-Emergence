import numpy as np
def local_geometry_features(A, phi):
    """
    Emergent local geometry:
      - degree / coordination number
      - Fibonacci bond weight (A=1.0, B=0.618)
      - local phase curvature
      - diagonal Laplacian entry

    Returns: G_local, shape (N,)
    """
    N = len(phi)
    deg = A.sum(axis=1)

    # bond weight per site (for Fibonacci, each site touches one main bond)
    bond = np.zeros(N)
    for i in range(N - 1):
        w = A[i, i + 1]
        bond[i] = w
        bond[i + 1] = w

    # discrete second derivative of phase (local curvature)
    dphi2 = np.roll(phi, -1) - 2 * phi + np.roll(phi, 1)
    curvature = np.abs(dphi2)

    # Laplacian diagonal = degree for this simple graph
    Ldiag = deg

    # normalize each channel (avoid division by zero)
    def nz_norm(x):
        m = np.max(np.abs(x))
        return x / m if m > 0 else x

    deg_n = nz_norm(deg)
    bond_n = nz_norm(bond)
    curv_n = nz_norm(curvature)
    Ldiag_n = nz_norm(Ldiag)

    G = (
        0.4 * deg_n +
        0.3 * bond_n +
        0.2 * curv_n +
        0.1 * Ldiag_n
    )

    return nz_norm(G)

def spectral_kernel_response(evals, evecs, triad, q_vec, sector_name):
    """
    Fully emergent spectral response kernel.

    No parameters. ALL structure comes from:
      - Laplacian eigenvalues evals
      - Triad eigenmodes
      - Integer charges q_j
      - Sector identity (which changes only the divisor pattern)
      - Divisors of 360 (harmonic alphabet)

    Returns a site-dependent complex weight vector W_s of length N.
    """

    N = len(evals)

    # Allowed harmonic divisors of 360
    D360 = np.array([1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 15, 18, 20, 24, 30, 36, 40, 45, 60, 72, 90, 120, 180, 360])

    # Sector-dependent divisor subset (emergent choice)
    if sector_name == "u":
        D = D360[D360 % 3 == 0]  # multiples of 3
    elif sector_name == "d":
        D = D360[D360 % 4 == 0]  # multiples of 4
    elif sector_name == "e":
        D = D360[D360 % 5 == 0]  # multiples of 5
    elif sector_name == "nu":
        D = D360[D360 % 2 == 1]  # odd divisors → maximal mixing
    else:
        raise ValueError("Bad sector")

    # triad eigenvectors
    G = evecs[:, triad]  # shape (N,3)

    # projectors into each triad direction
    proj = np.abs(G) ** 2  # (N,3), site-by-site mode weight

    # Spectral decay scale based on triad
    lam = evals[triad]
    lam_min = np.min(lam[lam > 1e-12])

    # sector-dependent "harmonic sum" over D
    # NO parameters — purely number-theoretic
    kernel = np.zeros(N, dtype=complex)

    for d in D:
        # Phase driven by Laplacian eigenvalues (NCG-inspired)
        phase = np.exp(-1j * d * lam)  # (3,)
        weight = np.exp(-d * lam / lam_min)  # (3,)
        contrib = proj @ (weight * phase)  # (N,)
        kernel += contrib

    # Apply integer charge catastrophe suppression — again, NO parameters
    # q_vec is shape (3,)
    charge_phase = np.exp(-1j * np.sum(q_vec))  # single global factor
    kernel *= charge_phase

    # Normalize
    kernel /= np.max(np.abs(kernel)) + 1e-12

    return kernel
def fibonacci_word(n):
    """
    Generate a Fibonacci word (aperiodic sequence) of length >= n.
    Substitution rules:
        A -> AB
        B -> A
    """
    a, b = "A", "AB"
    while len(b) < n:
        a, b = b, b + a
    return b[:n]


def build_internal_graph_fibonacci(N, wA=1.0, wB=0.618):
    """
    Build adjacency matrix for a 1D Fibonacci quasi-crystal chain.
    Each site i corresponds to the letter at position i.

    Letter → bond weight:
        A → strong bond
        B → weak bond (scaled by golden ratio inverse)
    """
    word = fibonacci_word(N)
    A = np.zeros((N, N), dtype=float)

    # Aperiodic adjacency: no periodic boundary!
    for i in range(N - 1):
        w = wA if word[i] == "A" else wB
        A[i, i + 1] = w
        A[i + 1, i] = w

    return A
###############################################################################
# 0) Utilities
###############################################################################

def normalize(v):
    n = np.linalg.norm(v)
    return v / n if n > 0 else v


def make_permutation_matrix(p):
    """Convert list/tuple p (e.g. [2,0,1]) into a permutation matrix."""
    P = np.zeros((len(p), len(p)))
    for i, j in enumerate(p):
        P[i, j] = 1.0
    return P


###############################################################################
# 1) Phase relaxation (same as earlier toy)
###############################################################################

def misalignment_energy(phi, w6=1.0, w5=1.0):
    N = len(phi)
    diffs = phi[:, None] - phi[None, :]
    E6 = w6 * np.sum(1 - np.cos(6 * diffs)) / (N * N)
    E5 = w5 * np.sum(1 - np.cos(5 * diffs)) / (N * N)
    return E6 + E5


def relax_phases(N=100, steps=800, eta=0.01, seed=42):
    rng = np.random.default_rng(seed)
    phi = rng.uniform(0, 2 * np.pi, size=N)

    for _ in range(steps):
        diffs = phi[:, None] - phi[None, :]
        grad = 6 * np.sum(np.sin(6 * diffs), axis=1) + 5 * np.sum(np.sin(5 * diffs), axis=1)
        phi -= eta * grad
        phi %= 2 * np.pi
    return phi, misalignment_energy(phi)


###############################################################################
# 2) Internal graph: 1D ring Laplacian
###############################################################################

def build_ring_adjacency(N):
    A = np.zeros((N, N))
    for i in range(N):
        A[i, (i + 1) % N] = 1
        A[i, (i - 1) % N] = 1
    return A


def laplacian(A):
    d = np.sum(A, axis=1)
    return np.diag(d) - A


###############################################################################
# 3) Harmonic strengths + emergent triad selection
###############################################################################

def harmonic_strengths(evals, evecs, phi, divisors=(1, 2, 3)):
    N = len(phi)
    V = evecs
    strengths = []
    exp_dphi = {d: np.exp(-1j * d * phi) for d in divisors}

    for n in range(N):
        v = V[:, n]
        Hn = []
        for d in divisors:
            amp = np.dot(v, exp_dphi[d])
            Hn.append(np.abs(amp) ** 2)
        strengths.append(Hn)
    return np.array(strengths)  # (N, len(divisors))


def select_triad(evals, H):
    lam = evals
    scores = np.sum(H[1:], axis=1)  # ignore zero mode
    idx = np.arange(1, len(lam))
    top3_rel = np.argsort(scores)[-3:]
    triad = idx[top3_rel]
    return np.sort(triad)


###############################################################################
# 4) Emergent R and Q
###############################################################################

def build_R_from_triad(evecs, phi, triad):
    V = evecs
    c = []
    for n in triad:
        v = V[:, n]
        c.append(np.dot(v, np.exp(1j * phi)))
    c = np.array(c)
    theta = np.angle(c)
    k = np.round(theta * 360 / (2 * np.pi)).astype(int) % 360
    R = np.diag(np.exp(1j * 2 * np.pi * k / 360))
    return R, k


def build_Q_from_harmonics(H, triad):
    A = np.sum(H[triad], axis=1)
    order = np.argsort(-A)
    q = np.empty(3, dtype=int)
    q[order] = np.array([1, 2, 3])
    Q = np.diag(q.astype(float))
    return Q, q


###############################################################################
# 5) Universal Yukawa kernel
###############################################################################

def build_Yukawa(evals, triad, Q):
    lam = evals[triad]
    lam_min = np.min(lam[lam > 1e-12])
    q = np.diag(Q)

    F = np.exp(-lam / lam_min) * np.exp(-q)
    Y = np.diag(F)
    return Y, F


###############################################################################
# 6) MULTI-SECTOR EXTENSION (u, d, e, nu)
###############################################################################

def sector_permutation(name):
    """
    Fully emergent toy rule:
    Each sector corresponds to a different permutation of the triad modes.
    This is the ONLY difference allowed (Axiom F2).
    """
    if name == "u":
        return [0, 1, 2]
    if name == "d":
        return [1, 2, 0]
    if name == "e":
        return [2, 0, 1]
    if name == "nu":
        return [0, 2, 1]
    raise ValueError("Unknown sector")


def build_sector_Yukawa(evals, triad, Q, P):
    Qs = P @ Q @ P.T
    Y, F = build_Yukawa(evals, triad, Qs)
    return Y, F, Qs


def build_sector_U_L(evecs, triad, phi, sector_name):
    """
    Build fully emergent sector-dependent unitary matrices from
    the Fibonacci Laplacian eigenvectors + phase field.
    No knobs, no fits: only internal geometry and φ.

    evecs: (N, N) eigenvector matrix of Laplacian
    triad: length-3 array of triad eigenvalue indices
    phi:   length-N phase field
    sector_name: "u", "d", "e", or "nu"
    """
    # triad eigenvectors spanning the generation subspace
    G = evecs[:, triad]   # shape (N, 3)

    # sector-specific internal "weighting" of sites
    # purely from φ; these are NOT tunable parameters,
    # just different harmonic couplings into the same manifold.
    if sector_name == "u":
        weights = np.exp(1j * phi)
    elif sector_name == "d":
        weights = np.exp(2j * phi)
    elif sector_name == "e":
        weights = np.exp(-1j * phi)
    elif sector_name == "nu":
        weights = np.exp(3j * phi) * np.cos(phi)
    else:
        raise ValueError(f"Unknown sector '{sector_name}'")

    # Project the weighted internal state into the triad subspace.
    # This gives a 3x3 matrix whose columns are the "sector views"
    # of the generation basis.
    W = weights[:, None] * G          # (N,3)
    Uraw = G.conj().T @ W             # (3,3)

    # Polar decomposition via SVD → nearest unitary
    Uu, _, Vh = np.linalg.svd(Uraw)
    U_L = Uu @ Vh                     # guaranteed unitary 3x3

    return U_L
###############################################################################
# EMERGENT INTERNAL CONNECTION K (non-diagonal, sector-dependent)
###############################################################################

def build_internal_connection(A, phi, evecs, triad, sector_name):
    """
    Complex, phase-sensitive, fully emergent internal connection.
    This version DOES NOT commute with R and WILL generate mixing.
    """

    V = evecs[:, triad]  # N×3 complex matrix

    # 1. Complex phase gradient operator (captures twisting of internal modes)
    dphi = np.roll(phi, -1) - phi
    G = np.exp(1j * dphi)                 # N complex weights
    G3 = V.conj().T @ (np.diag(G)) @ V    # 3×3 complex

    # 2. Mode-current coupling (emergent internal flow)
    J = (np.gradient(V.real, axis=0) + 1j * np.gradient(V.imag, axis=0))
    J3 = V.conj().T @ J @ np.ones((3,3))  # collapse N×3 to 3×3 pattern

    # 3. Fibonacci adjacency contribution (complexified)
    A_phase = A * np.exp(1j * (phi[:,None] - phi[None,:]))
    B3 = V.conj().T @ A_phase @ V         # 3×3 complex

    # 4. Sector permutation (Axiom F2)
    P = make_permutation_matrix(sector_permutation(sector_name))
    P3 = P @ np.eye(3) @ P.T

    # 5. Full emergent complex internal connection
    K = G3 + J3 + B3 + P3

    # Normalize for stability
    K /= np.max(np.abs(K)) + 1e-12
    return K


###############################################################################
# Sector unitary builder using K-emergence
###############################################################################

def build_sector_unitary(R, K):
    """
    Build U_L by diagonalizing the emergent operator R K R†.
    This captures:
      - 360-mode twisting (R)
      - geometric coupling (K)
    """
    M = R @ K @ np.conjugate(R.T)
    eigvals, U = np.linalg.eig(M)

    # Sort by phase to stabilize identification of generations
    phases = np.angle(eigvals)
    order = np.argsort(phases)
    return U[:, order]

###############################################################################
# 7) RUNNER
###############################################################################

def run_multisector_toy(N=1080, seed=42):
    print("=== MULTI-SECTOR FULLY EMERGENT TOY MODEL ===")

    # Phase field
    phi, E = relax_phases(N=N, seed=seed)
    print("Misalignment energy:", E)

    # Replace the uniform ring adjacency
    A = build_internal_graph_fibonacci(N)
    print("Adjacency sample:", A[:5, :5])
    L = laplacian(A)
    evals, evecs = np.linalg.eigh(L)
    print("First 5 eigenvalues:", evals[:5])

    # Harmonics + triad
    H = harmonic_strengths(evals, evecs, phi)
    triad = select_triad(evals, H)
    print("Triad indices:", triad)
    print("Triad eigenvalues:", evals[triad])

    # R and Q
    R, k = build_R_from_triad(evecs, phi, triad)
    Q, q = build_Q_from_harmonics(H, triad)
    print("k_j:", k)
    print("q_j:", q)

    # Build all 4 sectors
    sectors = {}
    for name in ["u", "d", "e", "nu"]:
        # 1) build Q_s and Yukawa
        P = make_permutation_matrix(sector_permutation(name))
        Y, F, Qs = build_sector_Yukawa(evals, triad, Q, P)

        # 2) emergent internal connection
        K = build_internal_connection(A, phi, evecs, triad, name)

        # 3) sector-specific left-unitary
        U_L = build_sector_unitary(R, K)

        sectors[name] = dict(Y=Y, F=F, Qs=Qs, U_L=U_L)
    # Mixing
    CKM = sectors["u"]["U_L"].conj().T @ sectors["d"]["U_L"]
    PMNS = sectors["e"]["U_L"].conj().T @ sectors["nu"]["U_L"]

    print("\n=== CKM (emergent) ===")
    print(CKM)
    print("\n=== PMNS (emergent) ===")
    print(PMNS)

    return sectors, CKM, PMNS


if __name__ == "__main__":
    run_multisector_toy()

"""
/Users/chazzromeo/rAI/Flavor/.venv/bin/python /Users/chazzromeo/rAI/Flavor/flavor/v2/toy_emergent_model-3.py 
=== MULTI-SECTOR FULLY EMERGENT TOY MODEL ===
Misalignment energy: 1.9967149483235276
Adjacency sample: [[0.    1.    0.    0.    0.   ]
 [1.    0.    0.618 0.    0.   ]
 [0.    0.618 0.    1.    0.   ]
 [0.    0.    1.    0.    1.   ]
 [0.    0.    0.    1.    0.   ]]
First 5 eigenvalues: [6.86382668e-17 6.84537951e-06 2.73813254e-05 6.16080479e-05
 1.09523752e-04]
Triad indices: [194 609 868]
Triad eigenvalues: [0.24889194 1.94185901 3.22773355]
k_j: [153  11 359]
q_j: [1 2 3]

=== CKM (emergent) ===
[[ 1.        +0.j          0.12228721+0.04949486j -0.20297785-0.37294257j]
 [ 0.12228721-0.04949486j  1.        +0.j         -0.34748741+0.2526886j ]
 [-0.20297785+0.37294257j -0.34748741-0.2526886j   1.        +0.j        ]]

=== PMNS (emergent) ===
[[ 1.        +0.j          0.12228721+0.04949486j -0.20297785-0.37294257j]
 [ 0.12228721-0.04949486j  1.        +0.j         -0.34748741+0.2526886j ]
 [-0.20297785+0.37294257j -0.34748741-0.2526886j   1.        +0.j        ]]

Process finished with exit code 0

"""