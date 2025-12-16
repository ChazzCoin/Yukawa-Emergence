import numpy as np

# ============================================================
# Fully derived harmonic alignment pipeline
# Parent on Z_360 -> Selection S^ -> parent moments -> sector lambdas
# -> embedding (gap d=7) -> holographic Yukawas & Majorana
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
    - linear phase pattern within triads
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
    equal phase spacing in each triad (n,2n,3n).
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
    smooth magnitudes in each triad towards their average.
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
    Derive all decay constants (lambda's) from parent moments.
    Only fixed *ratios* are chosen; absolute scale comes from <n>, <n^2>.
    """
    n1 = parent_moment(freqs, amps_sel, k=1)
    n2 = parent_moment(freqs, amps_sel, k=2)
    n_max = max(freqs)

    # Base scales from moments:
    base1 = n1 / n_max          # first moment scale
    base2 = np.sqrt(n2) / n_max # "RMS" scale

    # Sector weights as simple rational-ish factors
    # Up > down > lepton > neutrino, Majorana comparable to lepton/down
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
# 5. Embedding 9 boundary sites with a gap at d = 7
# ---------------------------------------------------

def cyclic_distance(a, b, N=N_CYCLE):
    d = abs(a - b)
    return d if d <= N // 2 else N - d


def embedding_score(positions):
    """
    Reward coverage of {1,2,3,4,5,6,8}, penalize distance 7.
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
    rng = np.random.default_rng(RNG_SEED)
    best_score = -1e9
    best_positions = None
    best_distances = None

    for _ in range(max_trials):
        positions = np.sort(rng.choice(N_CYCLE, size=num_sites, replace=False))
        score, dists, has7 = embedding_score(positions)
        if score > best_score:
            best_score, best_positions, best_distances = score, positions, dists
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


def rescale_distances(D, max_scale=8.0):
    """
    Compress raw distances to [0, max_scale] for richer exponential structure.
    """
    d_max = np.max(D)
    if d_max == 0:
        return D.astype(float)
    return (D / d_max) * max_scale


# ---------------------------------------------------
# 6. Holographic map with explicit d = 7 gap
# ---------------------------------------------------

def holographic_kernel_with_gap(dist_matrix, lambd, gap_distance=7.0, tol=1e-6):
    """
    Exponential kernel with explicit gap at d ~ gap_distance (scaled).
    """
    K = np.exp(-lambd * dist_matrix)
    mask_gap = np.abs(dist_matrix - gap_distance) < tol
    K[mask_gap] = 0.0
    return K


def build_sector_yukawa(freqs, amps, dist_matrix_scaled, lambd,
                        sector_phase_shift=0.0,
                        gap_distance=7.0):
    """
    Yukawa-like matrix:
        Y_ij = (sum |a_n|^2) * K(d_ij) * phase_ij
    with K having a gap at d = gap_distance in the scaled metric.
    """
    weights = np.abs(amps) ** 2
    parent_scale = np.sum(weights)

    K = holographic_kernel_with_gap(dist_matrix_scaled, lambd, gap_distance)

    base_phase = np.angle(amps[0]) + sector_phase_shift
    num = dist_matrix_scaled.shape[0]
    phases = np.zeros((num, num), dtype=np.complex128)
    for i in range(num):
        for j in range(num):
            phi_ij = base_phase * (i - j)
            phases[i, j] = np.exp(1j * phi_ij)

    Y = parent_scale * K * phases
    return Y


def build_all_sectors(freqs, amps, dist_matrix_scaled, lambdas, gap_distance=7.0):
    """
    Build Yukawa-like matrices for four sectors using parent-derived lambdas.
    """
    sectors = {}
    sectors["up"] = build_sector_yukawa(freqs, amps, dist_matrix_scaled,
                                        lambdas["up"], 0.0, gap_distance)
    sectors["down"] = build_sector_yukawa(freqs, amps, dist_matrix_scaled,
                                          lambdas["down"], np.pi / 6.0, gap_distance)
    sectors["charged_lepton"] = build_sector_yukawa(freqs, amps, dist_matrix_scaled,
                                                    lambdas["e"], np.pi / 3.0, gap_distance)
    sectors["neutrino_D"] = build_sector_yukawa(freqs, amps, dist_matrix_scaled,
                                                lambdas["nu"], np.pi / 2.0, gap_distance)
    return sectors


def build_majorana_matrix(dist_matrix_scaled, lambda_M, gap_distance=7.0):
    """
    Heavy Majorana matrix M_R with same geometry + gap, using parent-derived lambda_M.
    """
    K = holographic_kernel_with_gap(dist_matrix_scaled, lambda_M, gap_distance)
    # Ensure positive definiteness by adding diagonal mass scale 1
    M_R = K + np.eye(K.shape[0])
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

    # Embedding
    print("=== 9-site embedding with gap at distance 7 (raw distances) ===")
    positions, dists, score = search_embedding()
    print("Embedding positions (mod 360):", positions)
    print("Realized small distances (<=8):", sorted(dists))
    print("Embedding score:", score)
    print()

    D_raw = boundary_distances(positions)
    print("Boundary distance matrix D_ij (raw):")
    print(D_raw)
    print()

    # Rescale distances
    D_scaled = rescale_distances(D_raw, max_scale=8.0)
    print("Scaled distance matrix D_ij (approx in [0,8]):")
    print(np.round(D_scaled, 3))
    print()

    # Holographic Yukawas
    print("=== Holographic Yukawa-like matrices with d=7 gap (scaled) ===")
    sectors = build_all_sectors(freqs_sel, amps_sel, D_scaled, lambdas, gap_distance=7.0)
    for name, Y in sectors.items():
        summarize_matrix(f"Y_{name}", Y)

    # Heavy Majorana matrix from parent-derived lambda_M
    print("=== Heavy Majorana matrix M_R from same geometry ===")
    M_R = build_majorana_matrix(D_scaled, lambdas["M"], gap_distance=7.0)
    summarize_matrix("M_R", M_R)

    # Seesaw: light neutrinos
    print("=== Seesaw light neutrino mass matrix m_nu ===")
    Y_nu = sectors["neutrino_D"]
    m_nu = seesaw_light_neutrinos(Y_nu, M_R, v=1.0)
    summarize_matrix("m_nu", m_nu)

    # Toy 3x3 mixing
    print("=== Toy 3x3 mixing from charged lepton and neutrino sectors ===")
    Y_e = sectors["charged_lepton"][:3, :3]
    Y_nu_light = Y_nu[:3, :3]  # just to see structure

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
