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


def build_sector_U_L(R, K, P):
    """
    Sectoral left-unitary:
        U_L^{(s)} = eigvectors( P R K R^\dagger P^\dagger )
    """

    # Rotate R into sector frame
    R_s = P @ R @ P.T

    # Sector-specific internal resonance operator
    M_s = R_s @ K @ R_s.conj().T

    # Extract the unitary from eigen-decomposition
    evals, U = np.linalg.eig(M_s)

    # Sort by phase for stability
    order = np.argsort(np.angle(evals))
    U = U[:, order]

    return U
###############################################################################
# EMERGENT INTERNAL CONNECTION K (non-diagonal, sector-dependent)
###############################################################################

def build_internal_connection(A, phi, evecs, triad, q):
    """
    Fully emergent internal connection K:
      - graph + curvature (B, C),
      - triad overlap,
      - charge–phase commutator [G3, Q3].

    A: adjacency (N×N)
    phi: phase field (N,)
    evecs: eigenvectors of Laplacian (N×N)
    triad: 3 indices of chosen modes
    q: length-3 integer charges associated to the triad (e.g. [1,2,3])
    """

    N = len(phi)
    V = evecs[:, triad]                     # N×3 triad eigenvectors (complex)

    # --- 1. Local graph structure (quasi-crystal adjacency) ---
    B = A / (np.max(A) + 1e-12)

    # --- 2. Phase curvature kernel ---
    dphi = np.roll(phi, -1) - phi
    C = np.outer(dphi, dphi)
    C /= (np.max(np.abs(C)) + 1e-12)

    # --- 3. Triad overlap in physical space ---
    overlap = V.T @ V                        # 3×3, non-diagonal in quasi-crystals

    # --- 4. Geometry-projected curvature ---
    G3 = V.T @ (B + C) @ V                   # 3×3 complex
    G3 /= (np.max(np.abs(G3)) + 1e-12)

    # --- 5. Charge operator projected to triad ---
    q_triad = np.asarray(q, dtype=float)     # already length-3 triad charges
    Q3 = np.diag(q_triad)
    Q3 /= (np.sum(q_triad) + 1e-12)          # normalize → no knobs

    # --- 6. Charge–phase commutator (essential asymmetry) ---
    comm = G3 @ Q3 - Q3 @ G3                 # [G3, Q3]

    # --- 7. Assemble full K ---
    Ftriad = np.fft.fft(V, axis=0)
    Sspec = (np.abs(Ftriad).T @ np.abs(Ftriad))  # simple spectral coupling

    K = (
        0.4 * overlap +
        0.3 * Sspec +
        0.3 * G3 +
        1.0j * comm
    )

    # Normalize for stability
    K /= (np.max(np.abs(K)) + 1e-12)

    return K

def evolve_RQ_along_K(R, Q, K, n_steps=40, d_beta=0.05):
    """
    Emergent β-flow on the triad manifold:

        dR/dβ = [K, R]
        dQ/dβ = [K, Q]

    R: 3x3 unitary-ish generation operator
    Q: 3x3 Hermitian-ish charge matrix
    K: 3x3 internal connection from build_internal_connection

    n_steps, d_beta are algorithmic (integration resolution), not physical knobs.
    """
    R_flow = R.copy().astype(complex)
    Q_flow = Q.copy().astype(complex)

    for _ in range(n_steps):
        dR = K @ R_flow - R_flow @ K
        dQ = K @ Q_flow - Q_flow @ K

        R_flow = R_flow + d_beta * dR
        Q_flow = Q_flow + d_beta * dQ

        # Re-unitarize R_flow via polar-like rephasing
        eigvals, eigvecs = np.linalg.eig(R_flow)
        phases = np.exp(1j * np.angle(eigvals))
        R_flow = eigvecs @ np.diag(phases) @ np.linalg.inv(eigvecs)

        # Re-Hermitize Q_flow (keep it self-adjoint)
        Q_flow = 0.5 * (Q_flow + Q_flow.conj().T)

    return R_flow, Q_flow
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

    # Fibonacci internal graph
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

    # Emergent R and Q from triad
    R, k = build_R_from_triad(evecs, phi, triad)
    Q, q = build_Q_from_harmonics(H, triad)
    print("k_j:", k)
    print("q_j:", q)

    # Build internal connection K (3x3) from graph + phases + triad + q
    K = build_internal_connection(A, phi, evecs, triad, q)

    # Emergent β-flow on (R, Q) along K
    R_flow, Q_flow = evolve_RQ_along_K(R, Q, K, n_steps=40, d_beta=0.05)

    # Build all 4 sectors
    sectors = {}
    for name in ["u", "d", "e", "nu"]:
        # Sector permutation (the only explicit sector split)
        P = make_permutation_matrix(sector_permutation(name))

        # Sectoral left-unitary from FLOWED R and K
        U_L = build_sector_U_L(R_flow, K, P)

        # Sector Yukawa from FLOWED Q
        Y, F, Qs = build_sector_Yukawa(evals, triad, Q_flow, P)

        sectors[name] = dict(Y=Y, F=F, Qs=Qs, U_L=U_L)

    # Mixing matrices
    CKM = sectors["u"]["U_L"].conj().T @ sectors["d"]["U_L"]
    PMNS = sectors["e"]["U_L"].conj().T @ sectors["nu"]["U_L"]

    print("\n=== CKM (emergent) ===")
    print(CKM)
    print("\n=== PMNS (emergent) ===")
    print(PMNS)
    print("CKM |V|:")
    print(np.abs(CKM))
    print("PMNS |U|:")
    print(np.abs(PMNS))

    return sectors, CKM, PMNS

if __name__ == "__main__":
    run_multisector_toy()

"""
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
[[ 0.42909751-0.27528765j -0.43034941-0.46962462j  0.49714132+0.29528737j]
 [-0.24209022-0.82492627j  0.33031693-0.14335271j -0.25193414+0.2603054j ]
 [ 0.03022732-0.00820507j  0.63361682+0.251242j    0.72979292+0.04276045j]]

=== PMNS (emergent) ===
[[ 0.19308435+0.36061081j  0.71012836+0.46704683j -0.10741257+0.31420659j]
 [-0.68890183+0.58640121j  0.17844112-0.37379451j  0.06646076-0.07461425j]
 [ 0.11105503-0.04350798j  0.32545058-0.01016437j -0.36565186-0.86374253j]]
CKM |V|:
[[0.5098117  0.63698344 0.57822498]
 [0.85971566 0.36008231 0.36225642]
 [0.03132114 0.68161046 0.73104457]]
PMNS |U|:
[[0.40904978 0.84995002 0.33205909]
 [0.90468343 0.41420233 0.09992156]
 [0.11927349 0.32560927 0.93795119]]

"""