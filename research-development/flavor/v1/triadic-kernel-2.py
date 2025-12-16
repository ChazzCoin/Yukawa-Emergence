import numpy as np
import math

# ============================================================
# Harmonic alignment pipeline (triad-driven, consistent triad partition)
# Parent on Z_360 -> Selection S^ -> parent moments -> sector lambdas
# -> triad-based embedding -> emergent proto lattice L
# -> Yukawas & Majorana from L -> seesaw + PMNS-like mixing
#
# Triad partition is explicit and used consistently in:
# - misalignment functionals
# - P^phi and B operators
# - proto lattice construction
# ============================================================

N_CYCLE = 360
NUM_SITES = 9
RNG_SEED = 123

def set_rng_seed(seed: int):
    """
    Update the global RNG seed used in build_parent_state and search_embedding.
    """
    global RNG_SEED
    RNG_SEED = seed


# ----------------------------
# 1. Divisors and parent modes
# ----------------------------

def divisors(n: int):
    return [k for k in range(1, n + 1) if n % k == 0]


D360 = divisors(N_CYCLE)
def build_gcd_magnitude_lattice(L_proto, positions, N=N_CYCLE):
    """
    Given a proto lattice L_proto and boundary positions, build a
    gcd–projected magnitude lattice L_gcd:

      - For each gcd g = gcd(d_ij, N), compute mean |L_ij| over all pairs
        (i,j) with that gcd.
      - Then define L_gcd[i,j] = mean |L|_{gcd(d_ij,N)}, symmetric.
      - Diagonal entries set to 1 (self-alignment).

    This compresses L_proto onto the divisor lattice of N.
    """
    import math

    num = len(positions)

    # 1) Collect |L_ij| statistics per gcd
    gcd_sum = {}
    gcd_count = {}

    for i in range(num):
        for j in range(i + 1, num):
            d = cyclic_distance(positions[i], positions[j], N=N)
            if d == 0:
                continue
            g = math.gcd(d, N)
            val = abs(L_proto[i, j])

            gcd_sum[g] = gcd_sum.get(g, 0.0) + val
            gcd_count[g] = gcd_count.get(g, 0) + 1

    # 2) Build mapping g -> mean |L|
    gcd_mean = {}
    for g in gcd_sum:
        gcd_mean[g] = gcd_sum[g] / gcd_count[g]

    # 3) Construct L_gcd using only gcd-dependent magnitudes
    L_gcd = np.zeros_like(L_proto, dtype=float)

    for i in range(num):
        L_gcd[i, i] = 1.0  # self-alignment
        for j in range(i + 1, num):
            d = cyclic_distance(positions[i], positions[j], N=N)
            if d == 0:
                mag = 1.0
            else:
                g = math.gcd(d, N)
                mag = gcd_mean.get(g, 0.0)
            L_gcd[i, j] = mag
            L_gcd[j, i] = mag

    # 4) Normalize to [0,1] like L_norm
    max_val = np.max(np.abs(L_gcd))
    if max_val > 0:
        L_gcd /= max_val

    return L_gcd

def compute_distance_and_gcd_stats(L, positions, N=N_CYCLE, max_d=None, top_k=5):
    """
    Compute compact statistics for:
      - distance spectrum (mean |L_ij| vs distance),
      - gcd-based spectrum (mean |L_ij| vs gcd(d, N)).

    Returns:
      top_distances: list of (d, mean|L|) sorted by mean|L| desc, length <= top_k
      top_gcds:      list of (g, mean|L|, count) sorted by mean|L| desc, length <= top_k
    """
    # --- distance spectrum ---
    dist_spectrum = distance_alignment_spectrum(L, positions, max_d=max_d)
    top_distances = dist_spectrum[:top_k]

    # --- gcd spectrum (same logic as analyze_gcd_alignment, but no prints) ---
    import math
    num = len(positions)
    gcd_sum = {}
    gcd_count = {}

    for i in range(num):
        for j in range(i + 1, num):
            d = cyclic_distance(positions[i], positions[j], N=N)
            if d == 0:
                continue
            if max_d is not None and d > max_d:
                continue
            g = math.gcd(d, N)
            val = abs(L[i, j])
            gcd_sum[g] = gcd_sum.get(g, 0.0) + val
            gcd_count[g] = gcd_count.get(g, 0) + 1

    gcd_spectrum = []
    for g in gcd_sum:
        mean_val = gcd_sum[g] / gcd_count[g]
        gcd_spectrum.append((g, mean_val, gcd_count[g]))

    gcd_spectrum.sort(key=lambda x: x[1], reverse=True)
    top_gcds = gcd_spectrum[:top_k]

    return top_distances, top_gcds


def analyze_gcd_alignment(L, positions, N=N_CYCLE, max_d=None):
    """
    Analyze alignment as a function of gcd(d, N), where
      d = cyclic distance between boundary sites,
      N = 360 here.

    For each gcd g = gcd(d, N), we compute:
      - mean |L_ij| over all pairs (i,j) with gcd(d_ij, N) = g,
      - number of such pairs,
      - the list of distinct distances d in that gcd-class (for diagnostics).

    This reveals how the triadic kernel respects the divisor structure of N.
    """
    num = len(positions)
    # accumulate per gcd
    gcd_sum = {}
    gcd_count = {}
    gcd_dists = {}

    for i in range(num):
        for j in range(i + 1, num):
            d = cyclic_distance(positions[i], positions[j], N=N)
            if d == 0:
                continue
            if max_d is not None and d > max_d:
                continue

            g = math.gcd(d, N)
            val = abs(L[i, j])

            gcd_sum[g] = gcd_sum.get(g, 0.0) + val
            gcd_count[g] = gcd_count.get(g, 0) + 1
            if g not in gcd_dists:
                gcd_dists[g] = set()
            gcd_dists[g].add(d)

    # build spectrum: (g, mean |L|, count, sorted distances)
    spectrum_g = []
    for g in gcd_sum:
        mean_val = gcd_sum[g] / gcd_count[g]
        d_list = sorted(gcd_dists[g])
        spectrum_g.append((g, mean_val, gcd_count[g], d_list))

    # sort by mean alignment descending (most aligned gcd-classes first)
    spectrum_g.sort(key=lambda x: x[1], reverse=True)

    print("=== GCD-based alignment spectrum (grouped by gcd(d, 360)) ===")
    for g, mean_val, count, d_list in spectrum_g:
        print(
            f"  gcd = {g:3d}: mean |L| ~ {mean_val:.4f}, "
            f"count = {count:2d}, distances = {d_list}"
        )
    print()


# ---------------------------------------------------
# 2. Parent state |Psi> with triadic closure on Z_360
#     + explicit triad partition
# ---------------------------------------------------

def build_parent_state(gamma: float = 0.02):
    """
    |Psi_raw> = sum_{n in freqs} a_n |n>, with:

      - freqs: chosen subset of D_360 (same as before, 20 modes),
      - triads: non-overlapping (n,2n,3n) partition on freqs,
      - triad nodes get correlated magnitudes (set by seed n),
      - all phases initially RANDOM (no coherence),
      - non-triad nodes get their own magnitudes and random phases.

    Returns:
      freqs:       sorted list of active frequencies
      amps:        complex amplitudes on freqs (raw, misaligned phases)
      triads:      list of disjoint (n,2n,3n)
    """
    rng = np.random.default_rng(RNG_SEED)

    # Keep the same active frequencies you’ve been using:
    freqs = [
        1, 2, 3, 4, 5, 6, 8, 9, 10, 12,
        15, 18, 20, 24, 30, 36, 40, 45, 60, 90
    ]
    freqs = sorted(freqs)

    # Build a non-overlapping triad partition
    triads, triad_nodes = build_triad_partition(freqs)

    amps = np.zeros(len(freqs), dtype=np.complex128)
    idx_map = {n: i for i, n in enumerate(freqs)}

    # 1) Assign amplitudes for triad nodes
    for (n1, n2, n3) in triads:
        # magnitude profile still tied to the seed n1
        base_mag = np.exp(-gamma * n1)
        # we allow slight magnitude variation inside the triad
        mags = base_mag * (1.0 + 0.1 * rng.normal(size=3))
        # phases are fully random initially
        phases = 2.0 * np.pi * rng.random(size=3)

        for n, mag, phi in zip((n1, n2, n3), mags, phases):
            i = idx_map[n]
            amps[i] = mag * np.exp(1j * phi)

    # 2) Assign amplitudes for non-triad nodes
    for n in freqs:
        if n in triad_nodes:
            continue
        i = idx_map[n]
        base_mag = np.exp(-gamma * n)
        mag = base_mag * (1.0 + 0.1 * rng.normal())
        phi = 2.0 * np.pi * rng.random()
        amps[i] = mag * np.exp(1j * phi)

    norm = np.linalg.norm(amps)
    if norm == 0:
        raise RuntimeError("Parent amplitudes vanished; adjust gamma or freqs.")
    amps /= norm

    return freqs, amps, triads



# ---------------------------------------------------
# 3. Misalignment functionals + Selection Operator S^
# ---------------------------------------------------

def phase_misalignment(freqs, amps, triads):
    """
    Phase misalignment functional M_phi:
        sum_over_triads [ (Δφ_2 - Δφ_1)^2 ],
    where Δφ_1 = φ_2n - φ_n, Δφ_2 = φ_3n - φ_2n
    for triads (n,2n,3n) in the canonical triad list.

    Vanishes iff each triad has perfectly equal phase spacing.
    """
    phases = np.angle(amps)
    idx_map = {n: i for i, n in enumerate(freqs)}

    M_phi = 0.0
    for (n1, n2, n3) in triads:
        i1, i2, i3 = idx_map[n1], idx_map[n2], idx_map[n3]
        p1, p2, p3 = phases[i1], phases[i2], phases[i3]
        d1 = (p2 - p1 + np.pi) % (2*np.pi) - np.pi
        d2 = (p3 - p2 + np.pi) % (2*np.pi) - np.pi
        M_phi += (d2 - d1)**2
    return float(M_phi)


def magnitude_misalignment(freqs, amps, triads):
    """
    Magnitude misalignment functional M_B:
        sum_over_triads Var(|a_n|, |a_2n|, |a_3n|).
    Vanishes iff each triad has equal magnitudes.
    """
    mags = np.abs(amps)
    idx_map = {n: i for i, n in enumerate(freqs)}

    M_B = 0.0
    for (n1, n2, n3) in triads:
        i1, i2, i3 = idx_map[n1], idx_map[n2], idx_map[n3]
        m = np.array([mags[i1], mags[i2], mags[i3]])
        var = np.var(m)
        M_B += var
    return float(M_B)

def build_triad_partition(freqs):
    """
    Build a non-overlapping triad partition from a given freq list.

    For each n in sorted(freqs), if (n,2n,3n) are all in freqs and
    none of them have been used yet, we make a triad (n,2n,3n).
    Each frequency belongs to at most one triad.

    Returns:
      triads: list of disjoint (n,2n,3n)
      triad_nodes: set of all frequencies that appear in some triad
    """
    freq_set = set(freqs)
    used = set()
    triads = []

    for n in sorted(freqs):
        if n in used:
            continue
        if (2*n in freq_set) and (3*n in freq_set) and \
           (2*n not in used) and (3*n not in used):
            triads.append((n, 2*n, 3*n))
            used.update({n, 2*n, 3*n})

    triad_nodes = used
    return triads, triad_nodes


def apply_C360(freqs, amps):
    # freqs already in D360; just renormalize
    amps = amps / np.linalg.norm(amps)
    return freqs, amps


def apply_P_phi(freqs, amps, triads):
    """
    Phase-coherence projector:
    enforce equal phase spacing in each canonical triad (n,2n,3n).
    """
    amps_out = amps.copy()
    idx_map = {n: i for i, n in enumerate(freqs)}

    for (n1, n2, n3) in triads:
        i1, i2, i3 = idx_map[n1], idx_map[n2], idx_map[n3]
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

    amps_out /= np.linalg.norm(amps_out)
    return freqs, amps_out



def apply_B(freqs, amps, triads, alpha=0.5):
    """
    Geometric selector:
    smooth magnitudes in each canonical triad towards their average.

    new_mags = (1 - alpha)*mags + alpha*mean(mags)
    with 0 < alpha <= 1 ensures magnitude variance decreases
    unless already equal.
    """
    amps_out = amps.copy()
    idx_map = {n: i for i, n in enumerate(freqs)}

    for (n1, n2, n3) in triads:
        i1, i2, i3 = idx_map[n1], idx_map[n2], idx_map[n3]
        mags = np.abs([amps[i1], amps[i2], amps[i3]])
        phases = np.angle([amps[i1], amps[i2], amps[i3]])

        avg_mag = np.mean(mags)
        new_mags = (1 - alpha) * mags + alpha * avg_mag

        for idx, mag, phi in zip([i1, i2, i3], new_mags, phases):
            amps_out[idx] = mag * np.exp(1j * phi)

    amps_out /= np.linalg.norm(amps_out)
    return freqs, amps_out


# ---------------------------------------------------
# 4. Parent moments -> sector decay constants (lambdas)
# ---------------------------------------------------

def parent_moment(freqs, amps, k=1):
    """
    <n^k> with respect to |Psi_sel|^2 on D_360.
    """
    weights = np.abs(amps) ** 2
    ns = np.array(freqs, dtype=float)
    return np.sum(weights * (ns ** k))


def derive_sector_lambdas(freqs, amps_sel):
    """
    Derive decay constants (lambda's) from parent moments.
    Only fixed *ratios* are chosen; absolute scale from <n>, <n^2>.
    """
    n1 = parent_moment(freqs, amps_sel, k=1)
    n2 = parent_moment(freqs, amps_sel, k=2)
    n_max = max(freqs)

    base1 = n1 / n_max          # first moment scale
    base2 = np.sqrt(n2) / n_max # RMS scale

    # Sector weights as simple rational-ish factors
    c_up   = 6/5      # 1.2
    c_down = 1.0
    c_e    = 9/10     # 0.9
    c_nu   = 4/10     # 0.4
    c_M    = 11/10    # 1.1

    lambdas = {}
    lambdas["up"]   = c_up   * base1
    lambdas["down"] = c_down * base1
    lambdas["e"]    = c_e    * base1
    lambdas["nu"]   = c_nu   * base1
    lambdas["M"]    = c_M    * base2

    return lambdas


# ---------------------------------------------------
# 5. Triad-based embedding & proto lattice
# ---------------------------------------------------

def cyclic_distance(a, b, N=N_CYCLE):
    d = abs(a - b)
    return d if d <= N // 2 else N - d


def build_proto_lattice(freqs, amps, triads, positions):
    """
    Build proto lattice L_ij from canonical triads on Z_360:

        L_ij = Sum_{triads (n,2n,3n)} |a_n|^2
               [cos(n*theta) + cos(2n*theta) + cos(3n*theta)],

    where theta = 2π * d_ij / 360.
    """
    weights = np.abs(amps) ** 2
    idx_map = {n: i for i, n in enumerate(freqs)}

    num = len(positions)
    L = np.zeros((num, num), dtype=float)

    for i in range(num):
        for j in range(num):
            d = cyclic_distance(positions[i], positions[j])
            theta = 2.0 * np.pi * d / N_CYCLE
            s = 0.0
            for (n1, n2, n3) in triads:
                w = weights[idx_map[n1]]
                s += w * (np.cos(n1 * theta) +
                          np.cos(n2 * theta) +
                          np.cos(n3 * theta))
            L[i, j] = s

    # Normalize so that average diagonal ~ 1
    diag_mean = np.mean(np.diag(L))
    if abs(diag_mean) > 1e-12:
        L = L / diag_mean

    return L


def embedding_score(positions, freqs, amps, triads):
    """
    Score embedding using proto lattice L:
    - build L from canonical triads,
    - prefer Toeplitz-like structure (entries depend mainly on distance),
    - reward more distinct realized distances.
    """
    num = len(positions)
    L = build_proto_lattice(freqs, amps, triads, positions)

    # Collect means by distance (Toeplitz target)
    dist_sums = {}
    dist_counts = {}
    for i in range(num):
        for j in range(num):
            d = cyclic_distance(positions[i], positions[j])
            dist_sums[d] = dist_sums.get(d, 0.0) + L[i, j]
            dist_counts[d] = dist_counts.get(d, 0) + 1
    mean_by_d = {d: dist_sums[d] / dist_counts[d] for d in dist_sums}

    # Toeplitz error
    toeplitz_err = 0.0
    for i in range(num):
        for j in range(num):
            d = cyclic_distance(positions[i], positions[j])
            toeplitz_err += (L[i, j] - mean_by_d[d]) ** 2

    # Variety: more distinct nonzero distances is better
    distinct_d = len([d for d in mean_by_d if d > 0])

    score = -toeplitz_err + 0.1 * distinct_d
    return score, L


def search_embedding(freqs, amps, triads, num_sites=NUM_SITES, max_trials=20000):
    """
    Random search for an embedding of num_sites points on Z_360
    that optimizes triad-based lattice coherence.
    """
    rng = np.random.default_rng(RNG_SEED)
    best_score = -1e18
    best_positions = None
    best_L = None

    for _ in range(max_trials):
        positions = np.sort(rng.choice(N_CYCLE, size=num_sites, replace=False))
        score, L = embedding_score(positions, freqs, amps, triads)
        if score > best_score:
            best_score = score
            best_positions = positions
            best_L = L

    return best_positions, best_L, best_score


def boundary_distances(positions):
    num = len(positions)
    D = np.zeros((num, num), dtype=int)
    for i in range(num):
        for j in range(num):
            D[i, j] = cyclic_distance(positions[i], positions[j])
    return D


def rescale_distances(D, max_scale=8.0):
    """
    Compress raw distances to [0, max_scale] (for diagnostics only).
    """
    d_max = np.max(D)
    if d_max == 0:
        return D.astype(float)
    return (D / d_max) * max_scale


def distance_alignment_spectrum(L, positions, max_d=None):
    """
    Compute mean |L_ij| vs distance and return a sorted list:
    [(d1, mean|L|), ...] sorted from most aligned (largest mean|L|)
    to most entropic (smallest mean|L|).
    """
    num = len(positions)
    dist_sum_abs = {}
    dist_count = {}

    for i in range(num):
        for j in range(i + 1, num):
            d = cyclic_distance(positions[i], positions[j])
            if max_d is not None and d > max_d:
                continue
            if d == 0:
                continue
            dist_sum_abs[d] = dist_sum_abs.get(d, 0.0) + abs(L[i, j])
            dist_count[d] = dist_count.get(d, 0) + 1

    mean_abs = {d: dist_sum_abs[d] / dist_count[d] for d in dist_sum_abs}
    spectrum = sorted(mean_abs.items(), key=lambda kv: kv[1], reverse=True)
    return spectrum


def normalize_proto_lattice(L):
    """
    Normalize proto lattice to [0,1] with diag=1, for use as base kernel.
    """
    Lmin = np.min(L)
    Lmax = np.max(L)
    if Lmax > Lmin:
        L_norm = (L - Lmin) / (Lmax - Lmin)
    else:
        L_norm = np.ones_like(L)

    n = L_norm.shape[0]
    for i in range(n):
        L_norm[i, i] = 1.0

    return L_norm


# ---------------------------------------------------
# 6. Yukawas & Majorana from proto lattice
# ---------------------------------------------------

def build_sector_yukawa_from_L(L_norm, lambd_S, sector_phase_shift, amps):
    """
    Build sector Yukawa from normalized proto lattice L_norm:

        K_S = exp(-lambda_S * (1 - L_norm))

    A simple coherent phase pattern is added on top.
    """
    K = np.exp(-lambd_S * (1.0 - L_norm))

    base_phase = np.angle(amps[0]) + sector_phase_shift
    num = L_norm.shape[0]
    phases = np.zeros((num, num), dtype=np.complex128)
    for i in range(num):
        for j in range(num):
            phi_ij = base_phase * (i - j)
            phases[i, j] = np.exp(1j * phi_ij)

    Y = K * phases
    return Y


def build_all_sectors(freqs, amps, L_norm, lambdas):
    """
    Build Yukawa-like matrices for four sectors using proto lattice L_norm
    and parent-derived lambdas.
    """
    sectors = {}
    sectors["up"] = build_sector_yukawa_from_L(L_norm, lambdas["up"], 0.0, amps)
    sectors["down"] = build_sector_yukawa_from_L(L_norm, lambdas["down"], np.pi / 6.0, amps)
    sectors["charged_lepton"] = build_sector_yukawa_from_L(L_norm, lambdas["e"], np.pi / 3.0, amps)
    sectors["neutrino_D"] = build_sector_yukawa_from_L(L_norm, lambdas["nu"], np.pi / 2.0, amps)
    return sectors


def build_majorana_from_L(L_norm, lambda_M):
    """
    Heavy Majorana matrix from same proto lattice:

        M_R = exp(-lambda_M * (1 - L_norm)) + I
    """
    K_M = np.exp(-lambda_M * (1.0 - L_norm))
    M_R = K_M + np.eye(L_norm.shape[0])
    return M_R


# ---------------------------------------------------
# 7. Seesaw + mixing
# ---------------------------------------------------

def diagonalize_hermitian(M):
    """
    Diagonalize Hermitian M: M = U diag(m) U^\dagger
    Return eigenvalues sorted ascending and corresponding U.
    """
    m, U = np.linalg.eigh(M)
    idx = np.argsort(m)
    m_sorted = m[idx]
    U_sorted = U[:, idx]
    return m_sorted, U_sorted


def seesaw_light_neutrinos(Y_nu, M_R, v=1.0):
    """
    Type-I seesaw:
        m_nu = -v^2 * Y_nu^T M_R^{-1} Y_nu
    """
    M_R_inv = np.linalg.inv(M_R)
    m_nu = -v ** 2 * Y_nu.T @ M_R_inv @ Y_nu
    m_nu = 0.5 * (m_nu + m_nu.conj().T)
    return m_nu


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

def run_single_seed(seed: int):
    """
    Run the core pipeline for a single RNG seed and print a compact summary:
      - misalignment before / after S^,
      - top distance and gcd alignment classes,
      - leading singular values for Yukawas,
      - rough PMNS-like mixing matrix.
    """
    print(f"\n==================== Seed = {seed} ====================\n")
    set_rng_seed(seed)

    # --- Parent + triads ---
    freqs, amps, triads = build_parent_state(gamma=0.02)

    M_phi_0 = phase_misalignment(freqs, amps, triads)
    M_B_0   = magnitude_misalignment(freqs, amps, triads)

    # Apply S^
    freqs_C, amps_C = apply_C360(freqs, amps)
    freqs_P, amps_P = apply_P_phi(freqs_C, amps_C, triads)
    freqs_sel, amps_sel = apply_B(freqs_P, amps_P, triads, alpha=0.7)

    M_phi_P = phase_misalignment(freqs_P, amps_P, triads)
    M_B_P   = magnitude_misalignment(freqs_P, amps_P, triads)
    M_phi_S = phase_misalignment(freqs_sel, amps_sel, triads)
    M_B_S   = magnitude_misalignment(freqs_sel, amps_sel, triads)

    print("Misalignment diagnostics:")
    print(f"  M_phi (initial)    = {M_phi_0:.6f}")
    print(f"  M_phi (after P^)   = {M_phi_P:.6f}")
    print(f"  M_phi (after S^)   = {M_phi_S:.6f}")
    print(f"  M_B   (initial)    = {M_B_0:.6e}")
    print(f"  M_B   (after P^)   = {M_B_P:.6e}")
    print(f"  M_B   (after S^)   = {M_B_S:.6e}")
    print()

    # --- Lambdas from parent moments ---
    lambdas = derive_sector_lambdas(freqs_sel, amps_sel)

    # --- Embedding + proto lattice ---
    positions, L_proto, score = search_embedding(freqs_sel, amps_sel, triads)
    D_raw = boundary_distances(positions)

    print("Embedding summary:")
    print(f"  positions (mod 360): {positions}")
    print(f"  embedding score:     {score:.3f}")
    print()

    # Distance / gcd spectra (compact)
    top_distances, top_gcds = compute_distance_and_gcd_stats(L_proto, positions, N=N_CYCLE, top_k=5)

    print("Top distance alignment classes (by mean |L_ij|):")
    for d, m in top_distances:
        print(f"  d = {d:3d}: mean |L| ~ {m:.4f}")
    print()

    print("Top gcd(d,360) alignment classes (by mean |L_ij|):")
    for g, m, cnt in top_gcds:
        print(f"  gcd = {g:3d}: mean |L| ~ {m:.4f}, count = {cnt}")
    print()

    # --- Yukawas / Majorana / seesaw ---
    L_norm = normalize_proto_lattice(L_proto)
    sectors = build_all_sectors(freqs_sel, amps_sel, L_norm, lambdas)
    M_R = build_majorana_from_L(L_norm, lambdas["M"])
    Y_nu = sectors["neutrino_D"]
    m_nu = seesaw_light_neutrinos(Y_nu, M_R, v=1.0)

    # --- GCD–projected kernel ---
    L_gcd = build_gcd_magnitude_lattice(L_proto, positions, N=N_CYCLE)
    sectors_gcd = build_all_sectors(freqs_sel, amps_sel, L_gcd, lambdas)
    M_R_gcd = build_majorana_from_L(L_gcd, lambdas["M"])
    Y_nu_gcd = sectors_gcd["neutrino_D"]
    m_nu_gcd = seesaw_light_neutrinos(Y_nu_gcd, M_R_gcd, v=1.0)

    # Leading singular values per sector (just top few)
    print("Leading singular values per sector (triad kernel):")
    for name in ["up", "down", "charged_lepton", "neutrino_D"]:
        svals = np.linalg.svd(sectors[name], compute_uv=False)
        svals_sorted = np.sort(svals)[::-1]
        print(f"  {name:14s}: {np.round(svals_sorted[:4], 4)}")
    print()

    print("Leading singular values per sector (gcd-projected kernel):")
    for name in ["up", "down", "charged_lepton", "neutrino_D"]:
        svals = np.linalg.svd(sectors_gcd[name], compute_uv=False)
        svals_sorted = np.sort(svals)[::-1]
        print(f"  {name:14s}: {np.round(svals_sorted[:4], 4)}")
    print()

    # Toy 3x3 mixing: triad kernel
    Y_e = sectors["charged_lepton"][:3, :3]
    H_e = Y_e.conj().T @ Y_e
    H_nu = m_nu[:3, :3]
    m_e2, U_e = diagonalize_hermitian(H_e)
    m_nu_vals, U_nu = diagonalize_hermitian(H_nu)
    U_PMNS = U_e.conj().T @ U_nu

    # Toy 3x3 mixing: gcd kernel
    Y_e_g = sectors_gcd["charged_lepton"][:3, :3]
    H_e_g = Y_e_g.conj().T @ Y_e_g
    H_nu_g = m_nu_gcd[:3, :3]
    m_e2_g, U_e_g = diagonalize_hermitian(H_e_g)
    m_nu_vals_g, U_nu_g = diagonalize_hermitian(H_nu_g)
    U_PMNS_g = U_e_g.conj().T @ U_nu_g

    print("Toy mass spectra (3x3, triad kernel):")
    print("  m_e^2  :", np.round(m_e2, 6))
    print("  m_nu   :", np.round(m_nu_vals, 6))
    print()

    print("Toy mass spectra (3x3, gcd kernel):")
    print("  m_e^2  :", np.round(m_e2_g, 6))
    print("  m_nu   :", np.round(m_nu_vals_g, 6))
    print()

    print("PMNS-like |U|  (triad kernel):")
    print(np.round(np.abs(U_PMNS), 3))
    print()

    print("PMNS-like |U|  (gcd kernel):")
    print(np.round(np.abs(U_PMNS_g), 3))
    print()


# ---------------------------------------------------
# 8. Full pipeline with consistent triads
# ---------------------------------------------------
def run_multi_seed(seed_list):
    """
    Run a robustness scan over multiple seeds.

    For each seed:
      - rerun the core pipeline,
      - print compact diagnostics.

    Example:
      run_multi_seed([1, 2, 3, 10, 42, 123])
    """
    for seed in seed_list:
        run_single_seed(seed)

def run_pipeline():
    # Parent and explicit triads
    print("=== Parent state |Psi> with triadic closure on Z_360 ===")
    freqs, amps, triads = build_parent_state(gamma=0.02)
    print("Active parent frequencies:", freqs)
    print("Number of modes:", len(freqs))
    print("Number of triads in partition:", len(triads))
    print()

    # Misalignment of raw parent
    M_phi_0 = phase_misalignment(freqs, amps, triads)
    M_B_0 = magnitude_misalignment(freqs, amps, triads)
    print("Initial phase misalignment M_phi:", round(M_phi_0, 6))
    print("Initial magnitude misalignment M_B:", round(M_B_0, 6))
    print()

    print("=== Applying C^360 (normalization only) ===")
    freqs_C, amps_C = apply_C360(freqs, amps)
    M_phi_C = phase_misalignment(freqs_C, amps_C, triads)
    M_B_C = magnitude_misalignment(freqs_C, amps_C, triads)
    print("After C^360:  M_phi =", round(M_phi_C, 6), " M_B =", round(M_B_C, 6))
    print()

    print("=== Applying P^phi (phase-coherence projector) ===")
    freqs_P, amps_P = apply_P_phi(freqs_C, amps_C, triads)
    M_phi_P = phase_misalignment(freqs_P, amps_P, triads)
    M_B_P = magnitude_misalignment(freqs_P, amps_P, triads)
    print("After P^phi:  M_phi =", round(M_phi_P, 6), " M_B =", round(M_B_P, 6))
    print()

    print("=== Applying B (magnitude smoother) ===")
    freqs_sel, amps_sel = apply_B(freqs_P, amps_P, triads, alpha=0.7)
    M_phi_S = phase_misalignment(freqs_sel, amps_sel, triads)
    M_B_S = magnitude_misalignment(freqs_sel, amps_sel, triads)
    print("After B:      M_phi =", round(M_phi_S, 6), " M_B =", round(M_B_S, 6))
    print("Norm after full S^:", np.linalg.norm(amps_sel))
    print()

    # Parent-derived lambdas
    print("=== Deriving sector decay constants from parent moments ===")
    lambdas = derive_sector_lambdas(freqs_sel, amps_sel)
    for key, val in lambdas.items():
        print(f"lambda_{key} =", round(float(val), 4))
    print()

    # Embedding from triads
    print("=== Searching 9-site embedding via triad-based coherence ===")
    positions, L_proto, score = search_embedding(freqs_sel, amps_sel, triads)
    print("Embedding positions (mod 360):", positions)
    print("Embedding score:", score)
    print()

    D_raw = boundary_distances(positions)
    print("Boundary distance matrix D_ij (raw):")
    print(D_raw)
    print()

    D_scaled = rescale_distances(D_raw, max_scale=8.0)
    print("Scaled distance matrix D_ij (approx in [0,8]) (diagnostic):")
    print(np.round(D_scaled, 3))
    print()

    print("Proto lattice L_ij from triads (top-left 3x3, real):")
    print(np.round(L_proto[:3, :3], 4))
    print()

    # Distance alignment spectrum (diagnostic)
    print("=== Distance alignment spectrum (mean |L_ij| vs distance) ===")
    spectrum = distance_alignment_spectrum(L_proto, positions)
    for d, m in spectrum:
        print(f"  d = {d}: mean |L| ~ {m:.4f}")
    print()

    # GCD-based grouping: how alignment organizes by gcd(d, 360)
    analyze_gcd_alignment(L_proto, positions)

    # Normalized proto lattice as base kernel
    L_norm = normalize_proto_lattice(L_proto)
    print("Normalized proto lattice L_norm (top-left 3x3, real):")
    print(np.round(L_norm[:3, :3], 4))
    print()

    # Holographic Yukawas from L_norm
    print("=== Yukawa-like matrices from emergent proto lattice ===")
    sectors = build_all_sectors(freqs_sel, amps_sel, L_norm, lambdas)
    for name, Y in sectors.items():
        summarize_matrix(f"Y_{name}", Y)

    # Heavy Majorana from L_norm
    print("=== Heavy Majorana matrix M_R from same proto lattice ===")
    M_R = build_majorana_from_L(L_norm, lambdas["M"])
    summarize_matrix("M_R", M_R)

    # Seesaw: light neutrinos
    print("=== Seesaw light neutrino mass matrix m_nu ===")
    Y_nu = sectors["neutrino_D"]
    m_nu = seesaw_light_neutrinos(Y_nu, M_R, v=1.0)
    summarize_matrix("m_nu", m_nu)

    # Toy 3x3 mixing
    print("=== Toy 3x3 mixing from charged lepton and neutrino sectors ===")
    Y_e = sectors["charged_lepton"][:3, :3]
    H_e = Y_e.conj().T @ Y_e
    H_nu = m_nu[:3, :3]

    m_e2, U_e = diagonalize_hermitian(H_e)
    m_nu_light, U_nu = diagonalize_hermitian(H_nu)

    U_PMNS = U_e.conj().T @ U_nu

    print("Charged-lepton squared masses (toy units):", np.round(m_e2, 4))
    print("Light neutrino masses (toy units):", np.round(m_nu_light, 6))
    print("PMNS-like |U| matrix (absolute values):")
    print(np.round(np.abs(U_PMNS), 3))


if __name__ == "__main__":
    run_pipeline()
    # Robustness scan over multiple seeds
    seeds = [1, 2, 3, 10, 42, 123, 999]
    run_multi_seed(seeds)