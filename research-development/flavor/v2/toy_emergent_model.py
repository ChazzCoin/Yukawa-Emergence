import numpy as np

# ============================================================
# Toy Emergent Cymatic Model on an Internal Ring
# ============================================================

# This script:
#  1) Relaxes an internal phase field φ on N sites via a misalignment functional.
#  2) Builds a simple ring graph and its Laplacian.
#  3) Computes Laplacian eigenmodes (internal "cymatic patterns").
#  4) Measures divisor-harmonic strength H_n(d) for d ∈ {1,2,3}.
#  5) Selects a triad of modes purely from harmonic structure (FULL EMERGENCE).
#  6) Builds emergent R (360-cycle) and Q (integer charges) on generation space.
#  7) Builds a toy Yukawa spectrum from (λ_triad, Q) with a single universal rule.


# -----------------------------
# 1. Phase relaxation (misalignment)
# -----------------------------

def misalignment_energy(phi, w6=1.0, w5=1.0):
    """
    Misalignment functional over all pairs (all-to-all),
    encoding 6-fold and 5-fold phase preferences.

    E = w6 * <1 - cos(6Δφ)> + w5 * <1 - cos(5Δφ)>
    """
    N = len(phi)
    diffs = phi[:, None] - phi[None, :]
    cos6 = np.cos(6 * diffs)
    cos5 = np.cos(5 * diffs)
    E6 = w6 * np.sum(1.0 - cos6) / (N * N)
    E5 = w5 * np.sum(1.0 - cos5) / (N * N)
    return E6 + E5


def relax_phases(N=100, n_steps=1000, eta=0.01, w6=1.0, w5=1.0, seed=42):
    """
    Gradient descent on the misalignment functional.
    This is the engine that 'rotates the mirrors' (phases) to reduce strain.
    """
    rng = np.random.default_rng(seed)
    phi = rng.uniform(0, 2 * np.pi, size=N)

    for step in range(n_steps):
        diffs = phi[:, None] - phi[None, :]
        sin6 = np.sin(6 * diffs)
        sin5 = np.sin(5 * diffs)
        grad = 6 * w6 * np.sum(sin6, axis=1) + 5 * w5 * np.sum(sin5, axis=1)
        phi = phi - eta * grad
        phi = (phi + 2 * np.pi) % (2 * np.pi)

    E_final = misalignment_energy(phi, w6=w6, w5=w5)
    return phi, E_final


# -----------------------------
# 2. Simple internal graph: 1D ring + Laplacian
# -----------------------------

def build_ring_adjacency(N):
    """
    Build a simple ring graph: each node connected to its two neighbors.
    """
    A = np.zeros((N, N), dtype=float)
    for i in range(N):
        j_next = (i + 1) % N
        j_prev = (i - 1) % N
        A[i, j_next] = 1.0
        A[i, j_prev] = 1.0
    return A


def laplacian_from_adjacency(A):
    d = np.sum(A, axis=1)
    L = np.diag(d) - A
    return L


# -----------------------------
# 3. Divisor-harmonic strengths on the "living frequency manifold"
# -----------------------------

def harmonic_strengths(evals, evecs, phi, divisors=(1, 2, 3)):
    """
    For each eigenmode n, compute H_n(d) = |sum_j v_n(j) e^{-i d φ_j}|^2
    for d in given divisors. This measures how strongly each mode
    resonates with base-360 harmonic cycles (like cymatic corridors).
    """
    N = len(phi)
    lam = np.asarray(evals)
    V = np.asarray(evecs)  # columns V[:, n]

    # Precompute exponentials for speed
    strengths = []
    for n in range(N):
        v = V[:, n]
        H_d = []
        for d in divisors:
            phase_factor = np.exp(-1j * d * phi)
            amp = np.dot(v, phase_factor)
            H_d.append(np.abs(amp) ** 2)
        strengths.append(H_d)
    strengths = np.array(strengths)  # shape (N, len(divisors))
    return strengths  # H[n, k] corresponds to mode n, divisor divisors[k]


def select_generation_triad(evals, strengths, divisors=(1, 2, 3)):
    """
    FULLY EMERGENT triad selection:

    1. Ignore the lowest eigenvalue (zero mode).
    2. For each mode n, define a harmonic score:
         A_n = sum_d H_n(d)
    3. Pick the top 3 modes by A_n as the "generation triad".

    This is the simplest honest realization of:
      - modes are chosen by how well they align with divisor harmonics.
    """
    lam = np.asarray(evals)
    H = np.asarray(strengths)
    N, K = H.shape

    # ignore mode 0 (Laplacian zero mode)
    start = 1
    idx = np.arange(N)[start:]
    scores = np.sum(H[start:, :], axis=1)  # A_n

    # pick top 3 by score (no hand tuning)
    top3_rel = np.argsort(scores)[-3:]
    triad_indices = idx[top3_rel]

    # sort triad by eigenvalue ascending (just for readability)
    triad_indices = triad_indices[np.argsort(lam[triad_indices])]

    return triad_indices


# -----------------------------
# 4. Emergent R and Q on 3D generation space
# -----------------------------

def build_R_from_triad(evecs, phi, triad_indices):
    """
    Build a 3x3 R operator purely from emergent phase-mode structure.

    For each triad mode n_j, we compute a complex overlap:
        c_j = sum_i v_n_j(i) * e^{i φ_i}

    Then we extract an angle θ_j = arg(c_j) and map it to
    a nearest base-360 cycle index k_j.

    Finally, R = diag(exp(i 2π k_j / 360)).
    """
    V = np.asarray(evecs)
    triad_indices = np.asarray(triad_indices, dtype=int)

    c_list = []
    for n in triad_indices:
        v = V[:, n]
        c = np.dot(v, np.exp(1j * phi))
        c_list.append(c)
    c_list = np.array(c_list)
    theta = np.angle(c_list)  # angles in [-π, π]

    # Map to nearest 360-cycle
    k_list = np.round(theta * 360.0 / (2 * np.pi)).astype(int) % 360

    # Build R as a diagonal in this simple toy
    R_diag = np.exp(1j * 2 * np.pi * k_list / 360.0)
    R = np.diag(R_diag)

    return R, k_list


def build_Q_from_harmonics(strengths, triad_indices, divisors=(1, 2, 3)):
    """
    Emergent integer charge operator Q from harmonic strengths:

    - For each triad mode, define total harmonic score A_n.
    - Rank triad modes by A_n.
    - Assign charges q = (1,2,3) by rank (strongest → 1, weakest → 3).
    """
    H = np.asarray(strengths)
    triad_indices = np.asarray(triad_indices, dtype=int)

    # total score A_n for triad modes
    A_triad = np.sum(H[triad_indices, :], axis=1)
    # rank: largest score -> lowest charge
    order = np.argsort(-A_triad)  # descending by A
    q_values = np.empty(3, dtype=int)
    q_values[order] = np.array([1, 2, 3])

    Q = np.diag(q_values.astype(float))
    return Q, q_values


# -----------------------------
# 5. Toy emergent Yukawa spectrum
# -----------------------------

def build_toy_yukawa(evals, triad_indices, Q):
    """
    Build a toy Yukawa spectrum purely from:

      - the Laplacian eigenvalues at the triad (λ_triad),
      - the emergent integer charges Q.

    Universal rule (no sector knobs):
      F_j = exp(-λ_j / λ_min) * exp(- q_j)

    where λ_min is the smallest eigenvalue in the triad (excluding exact zero).
    This gives a simple hierarchical spectrum:
      - larger λ and larger q ⇒ smaller F_j.

    Return F (vector of 3 "masses") and Y = diag(F).
    """
    lam = np.asarray(evals)
    triad_indices = np.asarray(triad_indices, dtype=int)

    lam_triad = lam[triad_indices]
    lam_min = np.min(lam_triad[lam_triad > 1e-14])  # avoid zero
    lam_min = lam_min if lam_min > 0 else 1.0

    # charges from Q
    q_vals = np.diag(Q)

    # universal, emergent rule
    F = np.exp(-lam_triad / lam_min) * np.exp(-q_vals)

    Y = np.diag(F.astype(float))
    return F, Y, lam_triad


# -----------------------------
# 6. Main runner
# -----------------------------

def run_toy_model(N=100, seed=42):
    print(f"=== Toy Emergent Cymatic Model (N={N}, seed={seed}) ===")

    # 1) Relax phases
    phi, E_final = relax_phases(N=N, n_steps=1000, eta=0.01, w6=1.0, w5=1.0, seed=seed)
    print(f"Final misalignment energy: {E_final:.6e}")

    # 2) Build ring graph and Laplacian
    A = build_ring_adjacency(N)
    L = laplacian_from_adjacency(A)

    # 3) Laplacian spectrum
    evals, evecs = np.linalg.eigh(L)
    print("First 10 Laplacian eigenvalues:")
    print(evals[:10])

    # 4) Harmonic strengths H_n(d)
    divisors = (1, 2, 3)
    H = harmonic_strengths(evals, evecs, phi, divisors=divisors)

    # 5) Emergent triad selection
    triad_indices = select_generation_triad(evals, H, divisors=divisors)
    print("\nSelected generation triad indices:", triad_indices)
    print("Triad eigenvalues:", evals[triad_indices])

    print("\nHarmonic strengths H_n(d) for triad modes:")
    for idx in triad_indices:
        print(f"  mode {idx}: λ={evals[idx]:.6f}, H(d)={H[idx, :]}")

    # 6) Emergent R and Q
    R, k_list = build_R_from_triad(evecs, phi, triad_indices)
    Q, q_list = build_Q_from_harmonics(H, triad_indices, divisors=divisors)

    print("\nEmergent R (3x3):")
    print(R)
    print("Cycle indices k_j (mod 360):", k_list)

    print("\nEmergent Q (3x3):")
    print(Q)
    print("Integer charges q_j:", q_list)

    # 7) Toy Yukawa from (λ_triad, Q)
    F, Y, lam_triad = build_toy_yukawa(evals, triad_indices, Q)
    print("\nToy emergent 'mass' spectrum F_j:")
    for j, (lam_j, q_j, F_j) in enumerate(zip(lam_triad, q_list, F)):
        print(f"  gen {j}: λ={lam_j:.6f}, q={q_j}, F={F_j:.6e}")

    print("\nToy Yukawa matrix Y (diagonal in this toy):")
    print(Y)
    print("\n=== Done ===")


if __name__ == "__main__":
    run_toy_model(N=120, seed=42)