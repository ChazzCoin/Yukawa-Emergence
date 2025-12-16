import numpy as np


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
# 7) RUNNER
###############################################################################

def run_multisector_toy(N=360, seed=42):
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
        P = make_permutation_matrix(sector_permutation(name))
        Y, F, Qs = build_sector_Yukawa(evals, triad, Q, P)

        # NEW: sector-dependent emergent U_L from (evecs, triad, phi)
        U_L = build_sector_U_L(evecs, triad, phi, name)

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
/Users/chazzromeo/rAI/Flavor/.venv/bin/python /Users/chazzromeo/rAI/Flavor/flavor/v2/toy_emergent_model-2.py 
=== MULTI-SECTOR FULLY EMERGENT TOY MODEL ===
Misalignment energy: 1.99591686211232
Adjacency sample: [[0.    1.    0.    0.    0.   ]
 [1.    0.    0.618 0.    0.   ]
 [0.    0.618 0.    1.    0.   ]
 [0.    0.    1.    0.    1.   ]
 [0.    0.    0.    1.    0.   ]]
First 5 eigenvalues: [-6.71220520e-17  6.16083772e-05  2.46419490e-04  5.54474052e-04
  9.85638822e-04]
Triad indices: [ 41 161 308]
Triad eigenvalues: [0.10211546 1.40460122 3.36836579]
k_j: [244 250 219]
q_j: [1 3 2]

=== CKM (emergent) ===
[[-0.84338122-0.50966253j  0.09604095+0.06880655j -0.10691793+0.05968735j]
 [-0.04481193+0.10408094j -0.61576273+0.75732538j -0.14103791+0.12067264j]
 [-0.11390117+0.05601609j  0.00503122+0.18251304j  0.22365806-0.94896231j]]

=== PMNS (emergent) ===
[[ 5.81370079e-01+0.76206102j -2.78762947e-01-0.05514558j
  -6.30481937e-04-0.02283895j]
 [-1.24819872e-01-0.2534427j  -8.66936050e-01-0.35356195j
   5.88932924e-02-0.20033524j]
 [-6.19800724e-03-0.0376857j  -1.50330669e-01-0.141657j
  -5.92060290e-01+0.77803597j]]

Process finished with exit code 0

"""