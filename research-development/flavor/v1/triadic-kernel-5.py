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
# ---------------------------------------------------
# Phenomenological targets (rough neutrino sector)
# ---------------------------------------------------
THETA12_TARGET = 33.0  # degrees
THETA23_TARGET = 45.0  # degrees
THETA13_TARGET = 9.0   # degrees

# For a normal hierarchy: m2/m3 ~ sqrt(Δm21^2 / Δm31^2) ~ 0.17
R_NU21_TARGET = 0.17   # target ratio of m_nu2 / m_nu3 (by abs)
# Rough “data-ish” targets (you can tweak these)
TH12_EXP = 33.4   # deg
TH23_EXP = 49.0   # deg
TH13_EXP = 8.6    # deg
R21_EXP  = 0.17   # dimensionless (m2/m3)

def phenom_cost_weighted(theta12, theta23, theta13, r21,
                         w12=1.0, w23=1.0, w13=1.0, wr=1.0):
    """
    Tunable χ^2-like cost around approximate experimental values.

    w12, w23, w13, wr are relative weights; they don't have to be physical.
    """
    # Very rough denominators ('sigmas') to set the scale; tune as you like
    sig12 = 3.0
    sig23 = 4.0
    sig13 = 1.0
    sigr  = 0.05

    c12 = w12 * ((theta12 - TH12_EXP) / sig12) ** 2
    c23 = w23 * ((theta23 - TH23_EXP) / sig23) ** 2
    c13 = w13 * ((theta13 - TH13_EXP) / sig13) ** 2
    cr  = wr  * ((r21     - R21_EXP) / sigr ) ** 2

    return c12 + c23 + c13 + cr


def set_rng_seed(seed: int):
    """
    Update the global RNG seed used in build_parent_state and search_embedding.
    """
    global RNG_SEED
    RNG_SEED = seed

def phenom_cost(theta12, theta23, theta13, r21,
                t12=THETA12_TARGET,
                t23=THETA23_TARGET,
                t13=THETA13_TARGET,
                r21_target=R_NU21_TARGET,
                w_theta=1.0,
                w_ratio=1.0):
    """
    Simple 'alignment cost' combining angle deviations and one mass-ratio deviation.

    Lower is better. w_theta and w_ratio let you emphasize angles vs mass ratios.
    """
    d12 = theta12 - t12
    d23 = theta23 - t23
    d13 = theta13 - t13
    dr  = r21 - r21_target

    return (w_theta * (d12*d12 + d23*d23 + d13*d13) +
            w_ratio * (dr*dr))
# ============================================================
# 9. Direct 3x3 geometric benchmark with kappa and sites {1,2,5}
# ============================================================

def build_kappa_kernel_3x3(kappa: float,
                           sites=(1, 2, 5),
                           N: int = N_CYCLE) -> np.ndarray:
    """
    Direct 3x3 geometric kernel on a subset of boundary sites.

    K_ij = exp(-kappa * d_ij),  where d_ij is cyclic distance on Z_N
    between flavor sites (e.g. 1, 2, 5). This is the toy 'kappa-model'
    kernel you describe in the text.

    Returns a 3x3 real symmetric matrix.
    """
    num = len(sites)
    K = np.zeros((num, num), dtype=float)
    for i in range(num):
        for j in range(num):
            d = cyclic_distance(sites[i], sites[j], N=N)
            K[i, j] = math.exp(-kappa * d)
    return K

# ----------------------------
# 1. Divisors and parent modes
# ----------------------------

def divisors(n: int):
    return [k for k in range(1, n + 1) if n % k == 0]

def build_projection_matrix(L_ref: np.ndarray, k: int = 3) -> np.ndarray:
    """
    Build a k-dimensional projection P from a symmetric alignment matrix L_ref.

    L_ref is (N,N), real symmetric (e.g. L_norm or L_gcd).
    We take the k eigenvectors with largest eigenvalues:

        L_ref v_i = λ_i v_i,   λ_1 >= λ_2 >= ... >= λ_N

    and form P = [v_1, v_2, ..., v_k], so that P has shape (N, k) and
    P^† P = I_k (orthonormal columns).
    """
    evals, evecs = np.linalg.eigh(L_ref)
    idx = np.argsort(evals)[::-1]  # descending by eigenvalue
    top_vecs = evecs[:, idx[:k]]   # shape (N, k)
    return top_vecs


def project_matrix_to_3(M: np.ndarray, P: np.ndarray) -> np.ndarray:
    """
    Project a matrix M in boundary space (N x N) into the k-dim eigenmode
    subspace defined by P (N x k):

        M_3 = P^† M P

    For k = 3, this gives a 3x3 effective matrix (e.g. Yukawa or M_R)
    in the emergent generation basis.
    """
    return P.conj().T @ M @ P

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

def sort_spectrum_by_abs(m, U):
    """
    Given eigenvalues m and eigenvectors U (columns), reorder them
    by ascending |m| so that index 0 = lightest, 2 = heaviest.

    Returns (m_sorted, U_sorted).
    """
    idx = np.argsort(np.abs(m))
    return m[idx], U[:, idx]


def pmns_angles(U):
    """
    Extract approximate mixing angles (theta12, theta23, theta13) in degrees
    from a unitary 3x3 matrix U, using the standard PDG-style parameterization
    without CP phase (we ignore phases and just use absolute values).

    |U_e3| = s13
    |U_e2| = s12 c13
    |U_mu3| = s23 c13
    """
    U_abs = np.abs(U)

    s13 = np.clip(U_abs[0, 2], 0.0, 1.0)
    theta13 = np.arcsin(s13)
    c13 = np.cos(theta13) if theta13 < np.pi / 2 else 1e-12  # avoid divide-by-zero

    s12 = np.clip(U_abs[0, 1] / max(c13, 1e-12), 0.0, 1.0)
    s23 = np.clip(U_abs[1, 2] / max(c13, 1e-12), 0.0, 1.0)

    theta12 = np.arcsin(s12)
    theta23 = np.arcsin(s23)

    return np.degrees(theta12), np.degrees(theta23), np.degrees(theta13)


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

def effective_kernel_3x3_from_seed(summary, kernel="triad"):
    """
    Given a summary from run_single_seed and a choice of kernel ("triad" or "gcd"),
    construct the 3x3 *geometric* kernel in the eigenmode basis:

      L_eff3 = P^T L P

    where P is the 3x3 projection matrix built from that kernel.
    """
    L = summary["L_norm"] if kernel == "triad" else summary["L_gcd"]
    P = build_projection_matrix(L, k=3)
    L3 = P.conj().T @ L @ P
    # normalize diagonal to 1
    d = np.diag(L3)
    # avoid zero
    d[d == 0] = 1.0
    L3_norm = L3 / d.max()
    return L3_norm

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

def build_effective_kernel(L_triad, L_gcd, beta):
    return np.cos(beta) * L_triad + np.sin(beta) * L_gcd

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

def scan_beta_for_seed(seed: int,
                       beta_grid=None,
                       gamma: float = 0.02,
                       alpha: float = 0.7):
    """
    For a fixed seed:
      - build the parent state and proto lattice,
      - construct triad (L_norm) and gcd (L_gcd) kernels,
      - scan beta in K_eff = cos(beta) L_norm + sin(beta) L_gcd,
      - compute phenomenology cost for each beta,
      - report the best beta and corresponding angles / r21.

    This reuses the same pipeline as run_single_seed but with an
    interpolated geometric kernel.
    """
    print(f"\n******** beta-scan for seed = {seed} ********\n")
    set_rng_seed(seed)

    if beta_grid is None:
        beta_grid = np.linspace(0.0, 0.5*np.pi, 16)

    # --- Parent + triads ---
    freqs, amps, triads = build_parent_state(gamma=gamma)

    freqs_C, amps_C = apply_C360(freqs, amps)
    freqs_P, amps_P = apply_P_phi(freqs_C, amps_C, triads)
    freqs_sel, amps_sel = apply_B(freqs_P, amps_P, triads, alpha=alpha)

    lambdas = derive_sector_lambdas(freqs_sel, amps_sel)
    positions, L_proto, score = search_embedding(freqs_sel, amps_sel, triads)

    L_triad = normalize_proto_lattice(L_proto)
    L_gcd = build_gcd_magnitude_lattice(L_proto, positions, N=N_CYCLE)

    sector_names = ["up", "down", "charged_lepton", "neutrino_D"]

    best = {
        "beta": None,
        "cost": float("inf"),
        "angles": None,
        "r21": None,
        "mnu_ratios": None,
    }

    for beta in beta_grid:
        # Effective kernel
        L_eff = np.cos(beta) * L_triad + np.sin(beta) * L_gcd

        sectors_eff = build_all_sectors(freqs_sel, amps_sel, L_eff, lambdas)
        M_R_eff = build_majorana_from_L(L_eff, lambdas["M"])

        # 9D seesaw just for completeness (not strictly used)
        Y_nu_eff = sectors_eff["neutrino_D"]
        m_nu_eff = seesaw_light_neutrinos(Y_nu_eff, M_R_eff, v=1.0)

        # Projection to 3 generations from geometric kernel
        P_eff = build_projection_matrix(L_eff, k=3)

        Y3_eff = {name: project_matrix_to_3(sectors_eff[name], P_eff)
                  for name in sector_names}
        M_R3_eff = project_matrix_to_3(M_R_eff, P_eff)

        # 3D seesaw
        m_nu3_eff = seesaw_light_neutrinos(Y3_eff["neutrino_D"], M_R3_eff, v=1.0)
        H_e3_eff = Y3_eff["charged_lepton"].conj().T @ Y3_eff["charged_lepton"]

        m_e2_3_eff, U_e3_eff_raw = diagonalize_hermitian(H_e3_eff)
        m_nu_3_eff, U_nu3_eff_raw = diagonalize_hermitian(m_nu3_eff)

        m_e2_3_eff, U_e3_eff = sort_spectrum_by_abs(m_e2_3_eff, U_e3_eff_raw)
        m_nu_3_eff, U_nu3_eff = sort_spectrum_by_abs(m_nu_3_eff, U_nu3_eff_raw)

        U_PMNS3_eff = U_e3_eff.conj().T @ U_nu3_eff
        th12, th23, th13 = pmns_angles(U_PMNS3_eff)

        m_nu_ratios_eff = np.abs(m_nu_3_eff) / np.max(np.abs(m_nu_3_eff))
        r21_eff = float(m_nu_ratios_eff[1])

        cost_eff = phenom_cost(th12, th23, th13, r21_eff)

        print(f"beta = {beta:.3f}: cost = {cost_eff:.3f}, "
              f"angles = ({th12:.2f}, {th23:.2f}, {th13:.2f}), r21 = {r21_eff:.4f}")

        if cost_eff < best["cost"]:
            best["beta"] = float(beta)
            best["cost"] = float(cost_eff)
            best["angles"] = (float(th12), float(th23), float(th13))
            best["r21"] = float(r21_eff)
            best["mnu_ratios"] = m_nu_ratios_eff.copy()

    print("\n=== Best beta for this seed (triad–gcd interpolation) ===")
    print(f"Seed = {seed}")
    print(f"beta_best = {best['beta']:.3f}")
    th12_b, th23_b, th13_b = best["angles"]
    print(f"Angles (deg): theta12 = {th12_b:.2f}, "
          f"theta23 = {th23_b:.2f}, theta13 = {th13_b:.2f}")
    print(f"r21 = {best['r21']:.4f}")
    print("m_nu ratios:", np.round(best["mnu_ratios"], 6))
    print(f"Phenomenology cost at best beta: {best['cost']:.3f}")
    print()

    return best
def run_beta_scan_for_seed(seed: int,
                           betas=None,
                           verbose=True):
    """
    For a given seed:
      * build the triad and gcd kernels,
      * form an interpolated kernel L_eff(beta),
      * project to 3x3, do seesaw,
      * extract angles, r21,
      * evaluate the main phenomenology cost vs beta.
    """
    if betas is None:
        # match your existing grid: 0 to pi/2 in ~16 steps
        betas = np.linspace(0.0, 0.5 * np.pi, 16)

    if verbose:
        print(f"\n******** beta-scan for seed = {seed} ********\n")

    # -------------------------------
    # 1. Rebuild everything for seed
    # -------------------------------
    set_rng_seed(seed)

    # Parent + triads
    freqs, amps, triads = build_parent_state(gamma=0.02)

    # Apply S^ (C360, P_phi, B)
    freqs_C, amps_C     = apply_C360(freqs, amps)
    freqs_P, amps_P     = apply_P_phi(freqs_C, amps_C, triads)
    freqs_sel, amps_sel = apply_B(freqs_P, amps_P, triads, alpha=0.7)

    # Embedding + proto lattice
    positions, L_proto, score = search_embedding(freqs_sel, amps_sel, triads)

    # Triad and gcd kernels
    L_triad = normalize_proto_lattice(L_proto)
    L_gcd   = build_gcd_magnitude_lattice(L_proto, positions, N=N_CYCLE)

    # Sector Yukawas & Majorana for *both* kernels
    lambdas     = derive_sector_lambdas(freqs_sel, amps_sel)
    sectors_tri = build_all_sectors(freqs_sel, amps_sel, L_triad, lambdas)
    sectors_gcd = build_all_sectors(freqs_sel, amps_sel, L_gcd,   lambdas)

    M_R_triad = build_majorana_from_L(L_triad, lambdas["M"])
    M_R_gcd   = build_majorana_from_L(L_gcd,   lambdas["M"])

    # ----------------------------------------
    # 2. β interpolation & 3x3 flavor analysis
    # ----------------------------------------
    best_cost    = np.inf
    best_beta    = None
    best_summary = None

    for beta in betas:
        # Interpolation between kernels
        # If you used a different interpolation before, adapt this line.
        L_eff = np.cos(beta) * L_triad + np.sin(beta) * L_gcd

        # Build sectors & Majorana for this L_eff
        sectors_eff = build_all_sectors(freqs_sel, amps_sel, L_eff, lambdas)
        M_R_eff     = build_majorana_from_L(L_eff, lambdas["M"])

        # 3x3 projection using eigenmodes of L_eff
        P_eff  = build_projection_matrix(L_eff, k=3)
        Y3_e   = project_matrix_to_3(sectors_eff["charged_lepton"], P_eff)
        Y3_nuD = project_matrix_to_3(sectors_eff["neutrino_D"],     P_eff)
        M_R3   = project_matrix_to_3(M_R_eff,                       P_eff)

        # Seesaw in 3D
        m_nu3 = seesaw_light_neutrinos(Y3_nuD, M_R3, v=1.0)

        # Charged-lepton Hermitian mass-squared
        H_e3 = Y3_e.conj().T @ Y3_e

        # Diagonalize & sort
        m_e2, U_e_raw  = diagonalize_hermitian(H_e3)
        m_nu, U_nu_raw = diagonalize_hermitian(m_nu3)

        m_e2, U_e  = sort_spectrum_by_abs(m_e2,  U_e_raw)
        m_nu, U_nu = sort_spectrum_by_abs(m_nu,  U_nu_raw)

        # PMNS and angles
        U_PMNS = U_e.conj().T @ U_nu
        th12, th23, th13 = pmns_angles(U_PMNS)

        # Neutrino hierarchy ratios
        m_nu_ratios = np.abs(m_nu) / np.max(np.abs(m_nu))
        r21 = float(m_nu_ratios[1])

        # Your original cost
        cost = phenom_cost(th12, th23, th13, r21)

        if verbose:
            print(f"beta = {beta:5.3f}: "
                  f"cost = {cost:7.3f}, "
                  f"angles = ({th12:5.2f}, {th23:5.2f}, {th13:5.2f}), "
                  f"r21 = {r21:5.4f}")

        # Track best β
        if cost < best_cost:
            best_cost = cost
            best_beta = beta
            best_summary = (th12, th23, th13, r21, m_nu_ratios.copy())

    if verbose:
        print(f"\n=== Best beta for this seed (phenom_cost) ===")
        th12, th23, th13, r21, mrat = best_summary
        print(f"Seed = {seed}")
        print(f"beta_best = {best_beta:.3f}")
        print(f"Angles (deg): theta12 = {th12:.2f}, "
              f"theta23 = {th23:.2f}, theta13 = {th13:.2f}")
        print(f"r21 = {r21:.4f}")
        print(f"m_nu ratios: {np.round(mrat, 4)}")
        print(f"Phenomenology cost at best beta: {best_cost:.3f}")

    return {
        "seed": seed,
        "best_beta": best_beta,
        "best_cost": best_cost,
        "best_summary": best_summary,
    }

def build_sector_yukawa_from_K3(K_norm, lam_S, sector_phase_shift):
    """
    3x3 analogue of build_sector_yukawa_from_L:
        Y_S = exp(-lam_S * (1 - K_norm)) * phase_pattern
    """
    K = np.exp(-lam_S * (1.0 - K_norm))

    base_phase = sector_phase_shift
    num = K.shape[0]
    phases = np.zeros_like(K, dtype=np.complex128)
    for i in range(num):
        for j in range(num):
            phi_ij = base_phase * (i - j)
            phases[i, j] = np.exp(1j * phi_ij)

    return K * phases

def analyze_kappa_benchmark(kappa: float,
                            sites=(1, 2, 5),
                            lambdas_scale=(1.2, 1.0, 0.9, 0.4, 1.1),
                            v=1.0):
    """
    Analyze a direct 3x3 geometric benchmark with parameter kappa and
    flavor sites `sites`.

    We:
      - build K = K(kappa, sites),
      - define sector Yukawas as exp(-lambda_S * (1 - K_norm)),
      - build M_R similarly,
      - perform 3x3 seesaw and extract mixing angles, mass ratios, etc.

    lambdas_scale: (c_up, c_down, c_e, c_nu, c_M) are relative factors
    analogous to derive_sector_lambdas, but now we just treat them as
    dimensionless multipliers times a common base ~ kappa.
    """
    print("\n========== κ-benchmark analysis ==========")
    print(f"kappa = {kappa}, sites = {sites}")

    # 1) Basic geometric kernel (3x3), then normalize to [0,1] with diag=1
    K = build_kappa_kernel_3x3(kappa, sites=sites, N=N_CYCLE)

    K_norm = K / np.max(K)
    np.fill_diagonal(K_norm, 1.0)

    # 2) Sector "lambdas" from kappa (very simple ansatz: lambda_S = c_S * kappa)
    c_up, c_down, c_e, c_nu, c_M = lambdas_scale
    lam_up   = c_up   * kappa
    lam_down = c_down * kappa
    lam_e    = c_e    * kappa
    lam_nu   = c_nu   * kappa
    lam_M    = c_M    * kappa

    def yuk_from_K(lam):
        return np.exp(-lam * (1.0 - K_norm))

    Y_u = build_sector_yukawa_from_K3(K_norm, lam_up, 0.0)
    Y_d = build_sector_yukawa_from_K3(K_norm, lam_down, np.pi / 6)
    Y_e = build_sector_yukawa_from_K3(K_norm, lam_e, np.pi / 3)
    Y_nuD = build_sector_yukawa_from_K3(K_norm, lam_nu, np.pi / 2)
    M_R = np.exp(-lam_M * (1.0 - K_norm)) + np.eye(3)  # can also phase-twist if you like

    # 3) Seesaw + mixing
    m_nu = seesaw_light_neutrinos(Y_nuD, M_R, v=v)
    H_e  = Y_e.conj().T @ Y_e

    m_e2, U_e_raw  = diagonalize_hermitian(H_e)
    m_nu_eig, U_nu_raw = diagonalize_hermitian(m_nu)

    m_e2,  U_e  = sort_spectrum_by_abs(m_e2,  U_e_raw)
    m_nu_eig, U_nu = sort_spectrum_by_abs(m_nu_eig, U_nu_raw)

    U_PMNS = U_e.conj().T @ U_nu
    th12, th23, th13 = pmns_angles(U_PMNS)

    m_nu_ratios = np.abs(m_nu_eig) / np.max(np.abs(m_nu_eig))
    r21 = float(m_nu_ratios[1])

    # Charged lepton ratios
    m_e2_ratios = m_e2 / np.max(np.abs(m_e2))

    cost_simple = phenom_cost(th12, th23, th13, r21)
    cost_weight = phenom_cost_weighted(th12, th23, th13, r21,
                                       w12=1.0, w23=1.0, w13=1.0, wr=3.0)

    print("Kappa kernel K_norm:")
    print(np.round(K_norm, 4))
    print("\nCharged-lepton eigenvalues (m_e^2):", np.round(m_e2, 6))
    print("m_e^2 ratios (to max):", np.round(m_e2_ratios, 6))
    print("Neutrino eigenvalues (|m_nu|):", np.round(np.abs(m_nu_eig), 6))
    print("m_nu ratios (to max):", np.round(m_nu_ratios, 6))
    print("\nPMNS-like |U| (kappa benchmark):")
    print(np.round(np.abs(U_PMNS), 3))
    print(f"Angles (deg): theta12 = {th12:.2f}, "
          f"theta23 = {th23:.2f}, theta13 = {th13:.2f}")
    print(f"r21 = {r21:.4f}")
    print(f"phenom_cost       = {cost_simple:.3f}")
    print(f"phenom_cost_weight= {cost_weight:.3f}")
    print("===========================================\n")

    return {
        "kappa": kappa,
        "sites": sites,
        "K_norm": K_norm,
        "Y_u": Y_u,
        "Y_d": Y_d,
        "Y_e": Y_e,
        "Y_nuD": Y_nuD,
        "M_R": M_R,
        "m_e2": m_e2,
        "m_e2_ratios": m_e2_ratios,
        "m_nu": m_nu_eig,
        "m_nu_ratios": m_nu_ratios,
        "angles": (th12, th23, th13),
        "r21": r21,
        "cost_simple": cost_simple,
        "cost_weight": cost_weight,
    }

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

    Returns
    -------
    summary : dict
        Dictionary collecting key diagnostics for this seed:
        - seed, positions, L_proto, L_norm, L_gcd
        - angles_triad / angles_gcd
        - r21_triad / r21_gcd
        - m_nu_ratios_triad / m_nu_ratios_gcd
        - m_e2_ratios_triad / m_e2_ratios_gcd
        - cost_triad / cost_gcd
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

    # --- Yukawas / Majorana / seesaw (triad kernel) ---
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

    # -------------------------------------------------------
    # Eigenmode-based 9 -> 3 projection (triad vs gcd kernels)
    # -------------------------------------------------------

    # Build projection matrices from the *geometric* kernels
    P_triad = build_projection_matrix(L_norm, k=3)
    P_gcd = build_projection_matrix(L_gcd, k=3)

    # Project all sector Yukawas and M_R into 3D generation space
    sector_names = ["up", "down", "charged_lepton", "neutrino_D"]

    Y3_triad = {}
    Y3_gcd = {}

    for name in sector_names:
        Y3_triad[name] = project_matrix_to_3(sectors[name], P_triad)
        Y3_gcd[name] = project_matrix_to_3(sectors_gcd[name], P_gcd)

    M_R3_triad = project_matrix_to_3(M_R, P_triad)
    M_R3_gcd = project_matrix_to_3(M_R_gcd, P_gcd)

    # Seesaw in 3D generation basis
    m_nu3_triad = seesaw_light_neutrinos(Y3_triad["neutrino_D"], M_R3_triad, v=1.0)
    m_nu3_gcd = seesaw_light_neutrinos(Y3_gcd["neutrino_D"], M_R3_gcd, v=1.0)

    # Hermitian mass-squared matrices for charged leptons in 3D
    H_e3_triad = Y3_triad["charged_lepton"].conj().T @ Y3_triad["charged_lepton"]
    H_e3_gcd = Y3_gcd["charged_lepton"].conj().T @ Y3_gcd["charged_lepton"]

    # Diagonalize to get spectra + mixing (unsorted -> sorted by |m|)
    m_e2_3_triad, U_e3_triad_raw = diagonalize_hermitian(H_e3_triad)
    m_nu_3_triad, U_nu3_triad_raw = diagonalize_hermitian(m_nu3_triad)

    m_e2_3_gcd, U_e3_gcd_raw = diagonalize_hermitian(H_e3_gcd)
    m_nu_3_gcd, U_nu3_gcd_raw = diagonalize_hermitian(m_nu3_gcd)

    # Sort spectra by ascending |m| to define generation ordering
    m_e2_3_triad, U_e3_triad = sort_spectrum_by_abs(m_e2_3_triad, U_e3_triad_raw)
    m_nu_3_triad, U_nu3_triad = sort_spectrum_by_abs(m_nu_3_triad, U_nu3_triad_raw)

    m_e2_3_gcd, U_e3_gcd = sort_spectrum_by_abs(m_e2_3_gcd, U_e3_gcd_raw)
    m_nu_3_gcd, U_nu3_gcd = sort_spectrum_by_abs(m_nu_3_gcd, U_nu3_gcd_raw)

    # PMNS in sorted generation basis
    U_PMNS3_triad = U_e3_triad.conj().T @ U_nu3_triad
    U_PMNS3_gcd = U_e3_gcd.conj().T @ U_nu3_gcd

    # Extract mixing angles (degrees)
    th12_t, th23_t, th13_t = pmns_angles(U_PMNS3_triad)
    th12_g, th23_g, th13_g = pmns_angles(U_PMNS3_gcd)

    # Neutrino ratios (normalized to heaviest)
    m_nu_ratios_triad = np.abs(m_nu_3_triad) / np.max(np.abs(m_nu_3_triad))
    m_nu_ratios_gcd = np.abs(m_nu_3_gcd) / np.max(np.abs(m_nu_3_gcd))
    r21_triad = float(m_nu_ratios_triad[1])
    r21_gcd = float(m_nu_ratios_gcd[1])

    # Charged-lepton ratios (normalized to heaviest)
    m_e2_ratios_triad = m_e2_3_triad / np.max(np.abs(m_e2_3_triad))
    m_e2_ratios_gcd = m_e2_3_gcd / np.max(np.abs(m_e2_3_gcd))

    # Phenomenology cost
    cost_triad = phenom_cost(th12_t, th23_t, th13_t, r21_triad)
    cost_gcd = phenom_cost(th12_g, th23_g, th13_g, r21_gcd)
    alt_cost_triad = phenom_cost_weighted(th12_t, th23_t, th13_t, r21_triad,
                                          w12=1.0, w23=1.0, w13=1.0, wr=3.0)
    alt_cost_gcd = phenom_cost_weighted(th12_g, th23_g, th13_g, r21_gcd,
                                        w12=1.0, w23=1.0, w13=1.0, wr=3.0)
    print("alt_cost_triad")
    print(alt_cost_triad)
    print("alt_cost_gcd")
    print(alt_cost_gcd)
    # -------------------------
    # Print eigenmode 3x3 story
    # -------------------------

    print("=== Eigenmode-based 3x3 flavor (triad kernel) ===")
    print("m_e^2 (3-gen, triad):", np.round(m_e2_3_triad, 6))
    print("m_e^2 ratios (to max):", np.round(m_e2_ratios_triad, 6))
    print("m_nu  (3-gen, triad, abs):", np.round(np.abs(m_nu_3_triad), 6))
    print("m_nu ratios (to max):", np.round(m_nu_ratios_triad, 6))
    print("PMNS-like |U| (triad kernel, eigenmode 3x3):")
    print(np.round(np.abs(U_PMNS3_triad), 3))
    print(f"Angles (deg, triad): theta12 = {th12_t:.2f}, "
          f"theta23 = {th23_t:.2f}, theta13 = {th13_t:.2f}")
    print(f"r21 (triad) = {r21_triad:.4f},  phenom cost (triad) = {cost_triad:.3f}")
    print()

    print("=== Eigenmode-based 3x3 flavor (gcd kernel) ===")
    print("m_e^2 (3-gen, gcd):", np.round(m_e2_3_gcd, 6))
    print("m_e^2 ratios (to max):", np.round(m_e2_ratios_gcd, 6))
    print("m_nu  (3-gen, gcd, abs):", np.round(np.abs(m_nu_3_gcd), 6))
    print("m_nu ratios (to max):", np.round(m_nu_ratios_gcd, 6))
    print("PMNS-like |U| (gcd kernel, eigenmode 3x3):")
    print(np.round(np.abs(U_PMNS3_gcd), 3))
    print(f"Angles (deg, gcd):   theta12 = {th12_g:.2f}, "
          f"theta23 = {th23_g:.2f}, theta13 = {th13_g:.2f}")
    print(f"r21 (gcd)   = {r21_gcd:.4f},  phenom cost (gcd)   = {cost_gcd:.3f}")
    print()

    # -------------------------
    # Return a compact summary
    # -------------------------
    summary = {
        "seed": seed,
        "positions": positions,
        "L_proto": L_proto,
        "L_norm": L_norm,
        "L_gcd": L_gcd,
        "angles_triad": (th12_t, th23_t, th13_t),
        "angles_gcd": (th12_g, th23_g, th13_g),
        "r21_triad": r21_triad,
        "r21_gcd": r21_gcd,
        "m_nu_ratios_triad": m_nu_ratios_triad,
        "m_nu_ratios_gcd": m_nu_ratios_gcd,
        "m_e2_ratios_triad": m_e2_ratios_triad,
        "m_e2_ratios_gcd": m_e2_ratios_gcd,
        "cost_triad": float(cost_triad),
        "cost_gcd": float(cost_gcd),
    }
    return summary

def run_multi_seed(seed_list):
    """
    Run a robustness scan over multiple seeds.

    For each seed:
      - rerun the core pipeline,
      - print compact diagnostics,
      - collect a summary dict.

    At the end:
      - report the best triad-kernel configuration (minimal phenom cost),
      - report the best gcd-kernel configuration (minimal phenom cost).

    Example:
      run_multi_seed([1, 2, 3, 10, 42, 123])
    """
    best_triad = None
    best_triad_cost = float("inf")

    best_gcd = None
    best_gcd_cost = float("inf")

    summaries = []

    for seed in seed_list:
        summary = run_single_seed(seed)
        summaries.append(summary)

        # Update best triad
        if summary["cost_triad"] < best_triad_cost:
            best_triad_cost = summary["cost_triad"]
            best_triad = summary

        # Update best gcd
        if summary["cost_gcd"] < best_gcd_cost:
            best_gcd_cost = summary["cost_gcd"]
            best_gcd = summary

    # -------------------------
    # Global bests over seeds
    # -------------------------
    print("\n===================================================")
    print("=== Best triad & gcd configurations over seeds ===")
    print("===================================================\n")

    if best_triad is not None:
        th12_t, th23_t, th13_t = best_triad["angles_triad"]
        print(">>> Best triad kernel configuration:")
        print(f"  Seed: {best_triad['seed']}")
        print(f"  Positions (mod 360): {best_triad['positions']}")
        print(f"  Phenomenology cost (triad): {best_triad_cost:.3f}")
        print(f"  Angles (deg): theta12 = {th12_t:.2f}, "
              f"theta23 = {th23_t:.2f}, theta13 = {th13_t:.2f}")
        print(f"  r21 (triad): {best_triad['r21_triad']:.4f}")
        print("  m_nu ratios (triad):",
              np.round(best_triad["m_nu_ratios_triad"], 6))
        print("  m_e^2 ratios (triad):",
              np.round(best_triad["m_e2_ratios_triad"], 6))
        print()

    if best_gcd is not None:
        th12_g, th23_g, th13_g = best_gcd["angles_gcd"]
        print(">>> Best gcd kernel configuration:")
        print(f"  Seed: {best_gcd['seed']}")
        print(f"  Positions (mod 360): {best_gcd['positions']}")
        print(f"  Phenomenology cost (gcd): {best_gcd_cost:.3f}")
        print(f"  Angles (deg): theta12 = {th12_g:.2f}, "
              f"theta23 = {th23_g:.2f}, theta13 = {th13_g:.2f}")
        print(f"  r21 (gcd): {best_gcd['r21_gcd']:.4f}")
        print("  m_nu ratios (gcd):",
              np.round(best_gcd["m_nu_ratios_gcd"], 6))
        print("  m_e^2 ratios (gcd):",
              np.round(best_gcd["m_e2_ratios_gcd"], 6))
        print()
    return summaries, best_triad, best_gcd

if __name__ == "__main__":
    summaries, best_triad, best_gcd = run_multi_seed([1, 2, 3, 10, 42, 123, 999])

    # Compare κ-benchmark to your 9-site model
    kappa_bench = 0.24
    kappa_result = analyze_kappa_benchmark(kappa_bench, sites=(1, 2, 5))

    # Existing β-scans (optional)
    _res2   = run_beta_scan_for_seed(2)
    _res3   = run_beta_scan_for_seed(3)
    _res123 = run_beta_scan_for_seed(123)
