import numpy as np

# ============================================================
# Harmonic, fully aligned pipeline with a distance-7 gap
# Parent on Z_360 -> Selection Operator S^ -> Embedding search
# -> Holographic Yukawas with explicit gap at d=7
# ============================================================

N_CYCLE = 360
NUM_SITES = 9
RNG_SEED = 123


# ----------------------------
# 1. Divisors and parent modes
# ----------------------------

def divisors(n: int):
    return [k for k in range(1, n + 1) if n % k == 0]


D360 = divisors(N_CYCLE)


# ---------------------------------------------------
# 2. Parent state |Psi> with triadic closure on Z_360
# ---------------------------------------------------

def build_parent_state(gamma: float = 0.02):
    """
    Build a parent state |Psi> = sum_{n in D360} a_n |n>
    with:
      - triadic closure on chosen seeds (n,2n,3n)
      - exponential falloff in |a_n| with n (gamma)
      - simple linear phase pattern within triads
    """
    rng = np.random.default_rng(RNG_SEED)

    # Seed frequencies that divide 360 and admit 2n,3n also dividing 360
    seed_candidates = [1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 15, 18, 20, 24, 30]
    seeds = []
    for s in seed_candidates:
        if (2 * s in D360) and (3 * s in D360):
            seeds.append(s)

    active = set()
    triads = []
    for s in seeds:
        triad = [s, 2 * s, 3 * s]
        triads.append(triad)
        active.update(triad)

    freqs = sorted(active)
    amps = np.zeros(len(freqs), dtype=np.complex128)

    for triad in triads:
        base_mag = np.exp(-gamma * triad[0])
        mags = base_mag * (1.0 + 0.1 * rng.normal(size=3))
        base_phase = 2.0 * np.pi * rng.random()
        delta_phase = 2.0 * np.pi / 360.0
        phases = [
            base_phase,
            base_phase + delta_phase,
            base_phase + 2.0 * delta_phase,
        ]
        for n, mag, phi in zip(triad, mags, phases):
            idx = freqs.index(n)
            amps[idx] = mag * np.exp(1j * phi)

    norm = np.linalg.norm(amps)
    if norm == 0:
        raise RuntimeError("Parent amplitudes vanished; adjust gamma or seeds.")
    amps /= norm
    return freqs, amps


# ---------------------------------------------------
# 3. Selection Operator S^ = C^360 B^ P^phi
# ---------------------------------------------------

def apply_C360(freqs, amps):
    """
    Harmonic closure C^360:
    Here freqs are already divisors of 360; just renormalize.
    """
    amps = amps / np.linalg.norm(amps)
    return freqs, amps


def apply_P_phi(freqs, amps):
    """
    Phase-coherence projector P^phi:
    Enforce equal phase spacing within each triad (n, 2n, 3n).
    Magnitudes fixed, phases rebuilt linearly.
    """
    amps_out = amps.copy()
    freq_to_idx = {n: i for i, n in enumerate(freqs)}
    processed = set()

    for n in freqs:
        if n in processed:
            continue
        if (2 * n in freq_to_idx) and (3 * n in freq_to_idx):
            i1, i2, i3 = freq_to_idx[n], freq_to_idx[2 * n], freq_to_idx[3 * n]
            mags = np.abs([amps[i1], amps[i2], amps[i3]])

            base_phase = np.angle(amps[i1])
            delta_phase = 2.0 * np.pi / 360.0
            new_phases = [
                base_phase,
                base_phase + delta_phase,
                base_phase + 2.0 * delta_phase,
            ]
            for idx, mag, phi in zip([i1, i2, i3], mags, new_phases):
                amps_out[idx] = mag * np.exp(1j * phi)

            processed.update([n, 2 * n, 3 * n])

    amps_out /= np.linalg.norm(amps_out)
    return freqs, amps_out


def apply_B(freqs, amps, alpha=0.5):
    """
    Geometric selector B^:
    Smooth magnitudes within each triad towards their average
    (one "misalignment-minimizing" step).
    """
    amps_out = amps.copy()
    freq_to_idx = {n: i for i, n in enumerate(freqs)}
    processed = set()

    for n in freqs:
        if n in processed:
            continue
        if (2 * n in freq_to_idx) and (3 * n in freq_to_idx):
            i1, i2, i3 = freq_to_idx[n], freq_to_idx[2 * n], freq_to_idx[3 * n]
            mags = np.abs([amps[i1], amps[i2], amps[i3]])
            phases = np.angle([amps[i1], amps[i2], amps[i3]])

            avg_mag = np.mean(mags)
            new_mags = (1 - alpha) * mags + alpha * avg_mag

            for idx, mag, phi in zip([i1, i2, i3], new_mags, phases):
                amps_out[idx] = mag * np.exp(1j * phi)

            processed.update([n, 2 * n, 3 * n])

    amps_out /= np.linalg.norm(amps_out)
    return freqs, amps_out


def apply_selection_operator(freqs, amps, alpha=0.5):
    """
    Full S^ = C^360 B^ P^phi
    """
    freqs, amps = apply_C360(freqs, amps)
    freqs, amps = apply_P_phi(freqs, amps)
    freqs, amps = apply_B(freqs, amps, alpha=alpha)
    return freqs, amps


# ---------------------------------------------------
# 4. Embedding 9 boundary sites with a gap at d = 7
# ---------------------------------------------------

def cyclic_distance(a, b, N=N_CYCLE):
    d = abs(a - b)
    return d if d <= N // 2 else N - d


def embedding_score(positions):
    """
    Score an embedding by:
      - rewarding coverage of distances {1,2,3,4,5,6,8}
      - heavily penalizing appearance of distance 7
    """
    num = len(positions)
    distances = set()
    has7 = False

    for i in range(num):
        for j in range(i + 1, num):
            d = cyclic_distance(positions[i], positions[j])
            if d == 7:
                has7 = True
            if d <= 8:
                distances.add(d)

    allowed = {1, 2, 3, 4, 5, 6, 8}
    coverage = len(distances & allowed)
    score = coverage - (10 if has7 else 0)
    return score, distances, has7


def search_embedding(num_sites=NUM_SITES, max_trials=20000):
    """
    Random search for an embedding of num_sites points on Z_360
    that maximizes coverage of {1,2,3,4,5,6,8} with no d=7.
    """
    rng = np.random.default_rng(RNG_SEED)
    best_score = -1e9
    best_positions = None
    best_distances = None

    for _ in range(max_trials):
        positions = np.sort(rng.choice(N_CYCLE, size=num_sites, replace=False))
        score, dists, has7 = embedding_score(positions)
        if score > best_score:
            best_score, best_positions, best_distances = score, positions, dists
            # Perfect coverage (except 7) – early exit
            if score == len({1, 2, 3, 4, 5, 6, 8}):
                break

    return best_positions, best_distances, best_score


def boundary_distances(positions):
    num = len(positions)
    D = np.zeros((num, num), dtype=int)
    for i in range(num):
        for j in range(num):
            D[i, j] = cyclic_distance(positions[i], positions[j])
    return D


# ---------------------------------------------------
# 5. Holographic map with explicit distance gap d = 7
# ---------------------------------------------------

def holographic_kernel_with_gap(dist_matrix, lambd, gap_distance=7):
    """
    Exponential distance kernel with a strict gap at a chosen distance:
        K_ij = exp(-lambda * d_ij) for d_ij != gap_distance
        K_ij = 0 for d_ij == gap_distance
    """
    K = np.exp(-lambd * dist_matrix)
    K[dist_matrix == gap_distance] = 0.0
    return K


def build_sector_yukawa(freqs, amps, dist_matrix, lambd,
                        sector_phase_shift=0.0, gap_distance=7):
    """
    Toy Yukawa for a sector:
        Y_ij = (sum_n |a_n|^2) * K(d_ij) * phase_ij
    with K having an explicit gap at d=gap_distance.
    """
    weights = np.abs(amps) ** 2
    parent_scale = np.sum(weights)

    K = holographic_kernel_with_gap(dist_matrix, lambd, gap_distance)

    base_phase = np.angle(amps[0]) + sector_phase_shift
    num = dist_matrix.shape[0]
    phases = np.zeros((num, num), dtype=np.complex128)
    for i in range(num):
        for j in range(num):
            phi_ij = base_phase * (i - j)
            phases[i, j] = np.exp(1j * phi_ij)

    Y = parent_scale * K * phases
    return Y


def build_all_sectors(freqs, amps, dist_matrix, gap_distance=7):
    """
    Build Yukawa-like matrices for four sectors with different decay constants.
    """
    sectors = {}

    # Chosen so up-sector decays strongest (largest hierarchy), neutrinos weakest
    lambd_u = 0.25
    lambd_d = 0.18
    lambd_e = 0.21
    lambd_nu = 0.05

    sectors["up"] = build_sector_yukawa(freqs, amps, dist_matrix, lambd_u,
                                        sector_phase_shift=0.0,
                                        gap_distance=gap_distance)
    sectors["down"] = build_sector_yukawa(freqs, amps, dist_matrix, lambd_d,
                                          sector_phase_shift=np.pi / 6.0,
                                          gap_distance=gap_distance)
    sectors["charged_lepton"] = build_sector_yukawa(freqs, amps, dist_matrix, lambd_e,
                                                    sector_phase_shift=np.pi / 3.0,
                                                    gap_distance=gap_distance)
    sectors["neutrino"] = build_sector_yukawa(freqs, amps, dist_matrix, lambd_nu,
                                              sector_phase_shift=np.pi / 2.0,
                                              gap_distance=gap_distance)
    return sectors


# ---------------------------------------------------
# 6. Simple diagnostics
# ---------------------------------------------------

def summarize_matrix(name, M):
    print(f"--- {name} ---")
    print("shape:", M.shape)
    svals = np.linalg.svd(M, compute_uv=False)
    print("singular values (approx):", np.round(svals, 4))
    print("top-left 3x3 block (real):")
    print(np.round(M.real[:3, :3], 4))
    print("top-left 3x3 block (imag):")
    print(np.round(M.imag[:3, :3], 4))
    print()


def run_pipeline():
    # Parent
    print("=== Building parent state |Psi> with triadic closure on Z_360 ===")
    freqs, amps = build_parent_state(gamma=0.02)
    print("Active parent frequencies:", freqs)
    print("Number of modes:", len(freqs))
    print()

    # Selection Operator
    print("=== Applying Selection Operator S^ = C^360 B^ P^phi ===")
    freqs_sel, amps_sel = apply_selection_operator(freqs, amps, alpha=0.7)
    print("Norm after selection:", np.linalg.norm(amps_sel))
    print()

    # Embedding with distance-7 gap
    print("=== Searching for 9-site embedding with unique gap at distance 7 ===")
    positions, dists, score = search_embedding()
    print("Best embedding positions (mod 360):", positions)
    print("Realized small distances (<=8):", sorted(dists))
    print("Embedding score:", score)
    print()

    D = boundary_distances(positions)
    print("Boundary distance matrix D_ij:")
    print(D)
    print()

    # Holographic Yukawas
    print("=== Building holographic Yukawa-like matrices with distance gap d=7 ===")
    sectors = build_all_sectors(freqs_sel, amps_sel, D, gap_distance=7)
    for name, Y in sectors.items():
        summarize_matrix(f"Y_{name}", Y)


if __name__ == "__main__":
    run_pipeline()
