import numpy as np

# ============================================================
# Fully derived harmonic alignment pipeline (upgraded)
# Parent on Z_360 -> Selection S^ -> parent moments -> sector lambdas
# -> triad-based embedding -> emergent proto lattice L
# -> emergent entropic gap -> Yukawas & Majorana from L
# -> seesaw + PMNS-like mixing
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
    - reward more distinct nonzero distances.

    No explicit mention of distance 7 here.
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

    # Toeplitz error: how far L_ij deviates from mean_by_d(distance)
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
    that optimizes triad-based lattice coherence (no explicit d=7 logic).
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


def find_entropic_gap(L, positions):
    """
    From the emergent proto lattice L, compute average |L_ij| vs distance
    and identify the distance d_gap with minimal average amplitude
    (candidate "entropic gap").
    """
    num = len(positions)
    dist_sum_abs = {}
    dist_count = {}

    for i in range(num):
        for j in range(i + 1, num):
            d = cyclic_distance(positions[i], positions[j])
            if d == 0:
                continue
            dist_sum_abs[d] = dist_sum_abs.get(d, 0.0) + abs(L[i, j])
            dist_count[d] = dist_count.get(d, 0) + 1

    mean_abs = {d: dist_sum_abs[d] / dist_count[d] for d in dist_sum_abs}
    d_gap = min(mean_abs, key=lambda d: mean_abs[d])

    return d_gap, mean_abs


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

    # Force exact diag = 1
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

    so that K_S has diag = 1, off-diagonals suppressed by both alignment and L_norm.
    Then multiply by a simple coherent phase pattern.
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

    so it shares the same harmonic skeleton but with its own alignment strength.
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

    # Emergent entropic gap
    d_gap, mean_abs = find_entropic_gap(L_proto, positions)
    print("Average |L_ij| vs distance (emergent):")
    for d in sorted(mean_abs):
        print(f"  d = {d}: mean |L| ~ {mean_abs[d]:.4f}")
    print()
    print("Emergent entropic gap candidate distance d_gap =", d_gap)
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
"""