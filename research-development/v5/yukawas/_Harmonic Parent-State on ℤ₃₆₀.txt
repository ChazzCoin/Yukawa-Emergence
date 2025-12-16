import numpy as np

# ============================================================
# Harmonic alignment pipeline (triad-driven, no special distance)
# Parent on Z_360 -> Selection S^ -> parent moments -> sector lambdas
# -> triad-based embedding -> emergent proto lattice L
# -> Yukawas & Majorana from L -> seesaw + PMNS-like mixing
#
# Distances are not hard-coded: alignment and entropy patterns
# emerge solely from triadic structure on D_360 and the embedding.
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
    |Psi> = sum_{n in D360} a_n |n>
    - triadic closure on seeds (n,2n,3n)
    - exponential falloff |a_n| ~ exp(-gamma * n)
    - linear phase pattern within triads (step 2π/360)
    """
    rng = np.random.default_rng(RNG_SEED)

    seed_candidates = [1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 15, 18, 20, 24, 30, 36, 40]
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
    # freqs already in D360; just renormalize
    amps = amps / np.linalg.norm(amps)
    return freqs, amps


def apply_P_phi(freqs, amps):
    """
    Phase-coherence projector:
    enforce equal phase spacing in each triad (n,2n,3n).
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
    Geometric selector:
    smooth magnitudes in each triad towards their average
    (one gradient-flow step towards triad magnitude alignment).
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
    freqs, amps = apply_C360(freqs, amps)
    freqs, amps = apply_P_phi(freqs, amps)
    freqs, amps = apply_B(freqs, amps, alpha=alpha)
    return freqs, amps


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


def build_triads_from_freqs(freqs):
    """
    Return list of triads (n,2n,3n) present in freqs.
    """
    triads = []
    freq_set = set(freqs)
    for n in freqs:
        if (2 * n in freq_set) and (3 * n in freq_set):
            triads.append((n, 2 * n, 3 * n))
    return triads


def build_proto_lattice(freqs, amps, positions):
    """
    Build proto lattice L_ij from triads on Z_360:

        L_ij = Sum_{triads (n,2n,3n)} |a_n|^2
               [cos(n*theta) + cos(2n*theta) + cos(3n*theta)],

    where theta = 2π * d_ij / 360.

    High |L_ij| = aligned distance (constructive interference).
    Low |L_ij| = entropic distance (destructive interference).
    """
    triads = build_triads_from_freqs(freqs)
    weights = np.abs(amps) ** 2
    idx_map = {n: i for i, n in enumerate(freqs)}

    num = len(positions)
    L = np.zeros((num, num), dtype=float)

    for i in range(num):
        for j in range(num):
            d = cyclic_distance(positions[i], positions[j])
            theta = 2.0 * np.pi * d / N_CYCLE
            s = 0.0
            for (n, n2, n3) in triads:
                w = weights[idx_map[n]]
                s += w * (np.cos(n * theta) +
                          np.cos(n2 * theta) +
                          np.cos(n3 * theta))
            L[i, j] = s

    # Normalize so that average diagonal ~ 1
    diag_mean = np.mean(np.diag(L))
    if abs(diag_mean) > 1e-12:
        L = L / diag_mean

    return L


def embedding_score(positions, freqs, amps):
    """
    Score embedding using proto lattice L:
    - build L from triads,
    - prefer Toeplitz-like structure (entries depend mainly on distance),
    - reward more distinct realized distances.

    No distance is singled out here; the spectrum is emergent.
    """
    num = len(positions)
    L = build_proto_lattice(freqs, amps, positions)

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


def search_embedding(freqs, amps, num_sites=NUM_SITES, max_trials=20000):
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
        score, L = embedding_score(positions, freqs, amps)
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

    This is purely diagnostic; no special distance is picked.
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

    so that:
      - K_S(ii) = 1,
      - off-diagonals are suppressed according to both:
            * triad-induced alignment L_norm(d),
            * sector alignment strength lambda_S.

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

    It shares the same harmonic skeleton but with its own alignment strength.
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


# ---------------------------------------------------
# 8. Full pipeline
# ---------------------------------------------------

def run_pipeline():
    # Parent and selection
    print("=== Parent state |Psi> with triadic closure on Z_360 ===")
    freqs, amps = build_parent_state(gamma=0.02)
    print("Active parent frequencies:", freqs)
    print("Number of modes:", len(freqs))
    print()

    print("=== Selection Operator S^ = C^360 B^ P^phi ===")
    freqs_sel, amps_sel = apply_selection_operator(freqs, amps, alpha=0.7)
    print("Norm after selection:", np.linalg.norm(amps_sel))
    print()

    # Parent-derived lambdas
    print("=== Deriving sector decay constants from parent moments ===")
    lambdas = derive_sector_lambdas(freqs_sel, amps_sel)
    for key, val in lambdas.items():
        print(f"lambda_{key} =", round(float(val), 4))
    print()

    # Embedding from triads
    print("=== Searching 9-site embedding via triad-based coherence ===")
    positions, L_proto, score = search_embedding(freqs_sel, amps_sel)
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

    # Distance alignment spectrum (diagnostic only)
    print("=== Distance alignment spectrum (mean |L_ij| vs distance) ===")
    spectrum = distance_alignment_spectrum(L_proto, positions)
    for d, m in spectrum:
        print(f"  d = {d}: mean |L| ~ {m:.4f}")
    print()

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

"""
RESULTS:
=== Parent state |Psi> with triadic closure on Z_360 ===
Active parent frequencies: [1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 15, 18, 20, 24, 30, 36, 40, 45, 60, 90]
Number of modes: 20

=== Selection Operator S^ = C^360 B^ P^phi ===
Norm after selection: 1.0

=== Deriving sector decay constants from parent moments ===
lambda_up = 0.2164
lambda_down = 0.1803
lambda_e = 0.1623
lambda_nu = 0.0721
lambda_M = 0.2972

=== Searching 9-site embedding via triad-based coherence ===
Embedding positions (mod 360): [ 75  81  83  88 184 260 283 284 295]
Embedding score: 3.6

Boundary distance matrix D_ij (raw):
[[  0   6   8  13 109 175 152 151 140]
 [  6   0   2   7 103 179 158 157 146]
 [  8   2   0   5 101 177 160 159 148]
 [ 13   7   5   0  96 172 165 164 153]
 [109 103 101  96   0  76  99 100 111]
 [175 179 177 172  76   0  23  24  35]
 [152 158 160 165  99  23   0   1  12]
 [151 157 159 164 100  24   1   0  11]
 [140 146 148 153 111  35  12  11   0]]

Scaled distance matrix D_ij (approx in [0,8]) (diagnostic):
[[0.    0.268 0.358 0.581 4.872 7.821 6.793 6.749 6.257]
 [0.268 0.    0.089 0.313 4.603 8.    7.061 7.017 6.525]
 [0.358 0.089 0.    0.223 4.514 7.911 7.151 7.106 6.615]
 [0.581 0.313 0.223 0.    4.291 7.687 7.374 7.33  6.838]
 [4.872 4.603 4.514 4.291 0.    3.397 4.425 4.469 4.961]
 [7.821 8.    7.911 7.687 3.397 0.    1.028 1.073 1.564]
 [6.793 7.061 7.151 7.374 4.425 1.028 0.    0.045 0.536]
 [6.749 7.017 7.106 7.33  4.469 1.073 0.045 0.    0.492]
 [6.257 6.525 6.615 6.838 4.961 1.564 0.536 0.492 0.   ]]

Proto lattice L_ij from triads (top-left 3x3, real):
[[1.     0.3132 0.2385]
 [0.3132 1.     0.7638]
 [0.2385 0.7638 1.    ]]

=== Distance alignment spectrum (mean |L_ij| vs distance) ===
  d = 1: mean |L| ~ 0.9269
  d = 2: mean |L| ~ 0.7638
  d = 179: mean |L| ~ 0.3942
  d = 5: mean |L| ~ 0.3816
  d = 6: mean |L| ~ 0.3132
  d = 7: mean |L| ~ 0.2761
  d = 8: mean |L| ~ 0.2385
  d = 177: mean |L| ~ 0.2093
  d = 158: mean |L| ~ 0.1814
  d = 101: mean |L| ~ 0.1729
  d = 103: mean |L| ~ 0.1640
  d = 157: mean |L| ~ 0.1607
  d = 100: mean |L| ~ 0.1575
  d = 159: mean |L| ~ 0.1569
  d = 11: mean |L| ~ 0.1418
  d = 153: mean |L| ~ 0.1376
  d = 164: mean |L| ~ 0.1353
  d = 99: mean |L| ~ 0.1344
  d = 165: mean |L| ~ 0.1249
  d = 12: mean |L| ~ 0.1187
  d = 160: mean |L| ~ 0.1127
  d = 76: mean |L| ~ 0.0941
  d = 146: mean |L| ~ 0.0905
  d = 152: mean |L| ~ 0.0872
  d = 24: mean |L| ~ 0.0766
  d = 140: mean |L| ~ 0.0719
  d = 175: mean |L| ~ 0.0690
  d = 151: mean |L| ~ 0.0613
  d = 96: mean |L| ~ 0.0551
  d = 13: mean |L| ~ 0.0479
  d = 23: mean |L| ~ 0.0384
  d = 109: mean |L| ~ 0.0369
  d = 111: mean |L| ~ 0.0363
  d = 172: mean |L| ~ 0.0327
  d = 148: mean |L| ~ 0.0244
  d = 35: mean |L| ~ 0.0162

Normalized proto lattice L_norm (top-left 3x3, real):
[[1.     0.4186 0.3554]
 [0.4186 1.     0.8001]
 [0.3554 0.8001 1.    ]]

=== Yukawa-like matrices from emergent proto lattice ===
--- Y_up ---
shape: (9, 9)
singular values (approx): [7.7426 0.39   0.2325 0.1812 0.163  0.148  0.0938 0.0365 0.0124]
top-left 3x3 block (real):
[[ 1.      0.3534 -0.5904]
 [ 0.3534  1.      0.3838]
 [-0.5904  0.3838  1.    ]]
top-left 3x3 block (imag):
[[ 0.     -0.8079 -0.6388]
 [ 0.8079  0.     -0.8774]
 [ 0.6388  0.8774  0.    ]]

--- Y_down ---
shape: (9, 9)
singular values (approx): [7.9364 0.3313 0.1974 0.1533 0.1378 0.125  0.0782 0.0304 0.0103]
top-left 3x3 block (real):
[[ 1.     -0.1    -0.8683]
 [-0.1     1.     -0.1071]
 [-0.8683 -0.1071  1.    ]]
top-left 3x3 block (imag):
[[ 0.     -0.8949  0.1964]
 [ 0.8949  0.     -0.9586]
 [-0.1964  0.9586  0.    ]]

--- Y_charged_lepton ---
shape: (9, 9)
singular values (approx): [8.0355 0.301  0.1793 0.139  0.1249 0.1132 0.0705 0.0274 0.0092]
top-left 3x3 block (real):
[[ 1.     -0.5397 -0.2671]
 [-0.5397  1.     -0.5741]
 [-0.2671 -0.5741  1.    ]]
top-left 3x3 block (imag):
[[ 0.     -0.7327  0.8602]
 [ 0.7327  0.     -0.7795]
 [-0.8602  0.7795  0.    ]]

--- Y_neutrino_D ---
shape: (9, 9)
singular values (approx): [8.5547e+00 1.4050e-01 8.3500e-02 6.4200e-02 5.7600e-02 5.2000e-02
 3.1400e-02 1.2100e-02 4.1000e-03]
top-left 3x3 block (real):
[[ 1.     -0.8786  0.6479]
 [-0.8786  1.     -0.9031]
 [ 0.6479 -0.9031  1.    ]]
top-left 3x3 block (imag):
[[ 0.     -0.3843  0.701 ]
 [ 0.3843  0.     -0.395 ]
 [-0.701   0.395   0.    ]]

=== Heavy Majorana matrix M_R from same proto lattice ===
--- M_R ---
shape: (9, 9)
singular values (approx): [8.3295 1.5133 1.3065 1.2406 1.2169 1.1976 1.1282 1.0502 1.0171]
top-left 3x3 block (real):
[[2.     0.8413 0.8257]
 [0.8413 2.     0.9423]
 [0.8257 0.9423 2.    ]]
top-left 3x3 block (imag):
[[0. 0. 0.]
 [0. 0. 0.]
 [0. 0. 0.]]

=== Seesaw light neutrino mass matrix m_nu ===
--- m_nu ---
shape: (9, 9)
singular values (approx): [6.9991e+00 4.1135e+00 1.6100e-02 3.9000e-03 2.6000e-03 1.8000e-03
 1.7000e-03 1.2000e-03 1.0000e-04]
top-left 3x3 block (real):
[[-1.2761  1.2558 -1.0364]
 [ 1.2558 -1.0758  0.7082]
 [-1.0364  0.7082 -0.2379]]
top-left 3x3 block (imag):
[[ 0.  0.  0.]
 [-0.  0. -0.]
 [-0.  0.  0.]]

=== Toy 3x3 mixing from charged lepton and neutrino sectors ===
Charged-lepton squared masses (toy units): [1.0000e-03 1.3400e-02 8.1386e+00]
Light neutrino masses (toy units): [-2.996893e+00 -6.490000e-04  4.076680e-01]
PMNS-like |U| matrix (absolute values):
[[0.318 0.787 0.528]
 [0.602 0.593 0.536]
 [0.733 0.171 0.659]]
"""
