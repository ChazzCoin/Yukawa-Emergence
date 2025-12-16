import numpy as np

# ============================================================
# A360 / Alignment-Math toy pipeline:
# Parent -> Selection Operator S^ -> Holographic Yukawas
# ============================================================

# ----------------------------
# 1. Global harmonic settings
# ----------------------------

N_CYCLE = 360         # parent cycle size
NUM_SITES = 9         # boundary / flavor sites
RNG_SEED = 123        # reproducibility


def divisors(n: int):
    """Return sorted list of positive divisors of n."""
    ds = []
    for k in range(1, n + 1):
        if n % k == 0:
            ds.append(k)
    return ds


D360 = divisors(N_CYCLE)


# ---------------------------------------------------
# 2. Parent state |Psi> and triadic closure building
# ---------------------------------------------------

def build_parent_state(gamma: float = 0.02):
    """
    Build a parent state |Psi> = sum_{n in D360} a_n |n>
    with:
      - triadic closure on a chosen seed set
      - exponential falloff in |a_n| with n (controlled by gamma)
      - simple linear phase pattern along triads
    Returns:
        freqs: list of n in D360 with nonzero amplitude
        amps:  np.array of complex amplitudes a_n
    """
    rng = np.random.default_rng(RNG_SEED)

    # Choose some "seed" frequencies whose triads we will activate.
    # These should divide 360 and give meaningful 2n, 3n also dividing 360.
    seed_candidates = [1, 2, 3, 4, 5, 6, 8, 9, 10]
    seeds = []
    for s in seed_candidates:
        if (2 * s) in D360 and (3 * s) in D360:
            seeds.append(s)

    active = set()
    triads = []
    for s in seeds:
        triad = [s, 2 * s, 3 * s]
        triads.append(triad)
        active.update(triad)

    freqs = sorted(active)
    amps = np.zeros(len(freqs), dtype=np.complex128)

    # Assign amplitudes triad by triad
    for triad in triads:
        # base magnitude and phase for this triad
        base_mag = np.exp(-gamma * triad[0])
        # small random spread in magnitudes to keep things non-degenerate
        mags = base_mag * (1.0 + 0.1 * rng.normal(size=3))

        base_phase = 2.0 * np.pi * rng.random()
        delta_phase = 2.0 * np.pi / 360.0  # one "phase tick" in base-360
        phases = [
            base_phase,
            base_phase + delta_phase,
            base_phase + 2.0 * delta_phase,
        ]

        for n, mag, phi in zip(triad, mags, phases):
            idx = freqs.index(n)
            amps[idx] = mag * np.exp(1j * phi)

    # Normalize to unit norm
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
    Project onto divisor-closed subset of D360.
    In this toy implementation, we simply keep everything in D360
    and renormalize (since freqs are already divisors of 360).
    """
    # All freqs should already be in D360; we just normalize.
    amps = amps / np.linalg.norm(amps)
    return freqs, amps


def apply_P_phi(freqs, amps):
    """
    Phase-coherence projector P^phi:
    Enforce equal phase spacing within each triad (n, 2n, 3n).
    We rebuild phases to be exactly linear, keeping magnitudes fixed.
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

            # new coherent phases
            base_phase = np.angle(amps[i1])
            delta_phase = 2.0 * np.pi / 360.0
            new_phases = [
                base_phase,
                base_phase + delta_phase,
                base_phase + 2.0 * delta_phase,
            ]
            amps_out[i1] = mags[0] * np.exp(1j * new_phases[0])
            amps_out[i2] = mags[1] * np.exp(1j * new_phases[1])
            amps_out[i3] = mags[2] * np.exp(1j * new_phases[2])

            processed.update([n, 2 * n, 3 * n])

    # Normalize
    amps_out /= np.linalg.norm(amps_out)
    return freqs, amps_out


def apply_B(freqs, amps, alpha=0.5):
    """
    Geometric / misalignment selector B^:
    Minimize variance of magnitudes within each triad (n, 2n, 3n).
    We do a single "smoothing" step towards equal magnitudes.
    alpha in [0,1] controls the strength of the smoothing.
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

            amps_out[i1] = new_mags[0] * np.exp(1j * phases[0])
            amps_out[i2] = new_mags[1] * np.exp(1j * phases[1])
            amps_out[i3] = new_mags[2] * np.exp(1j * phases[2])

            processed.update([n, 2 * n, 3 * n])

    # Normalize
    amps_out /= np.linalg.norm(amps_out)
    return freqs, amps_out


def apply_selection_operator(freqs, amps, alpha=0.5):
    """
    Full selection operator S^ = C^360 B^ P^phi.
    """
    freqs, amps = apply_C360(freqs, amps)
    freqs, amps = apply_P_phi(freqs, amps)
    freqs, amps = apply_B(freqs, amps, alpha=alpha)
    return freqs, amps


# ---------------------------------------------------
# 4. Boundary / flavor sites on the 360-cycle
# ---------------------------------------------------

def build_embedding(num_sites=NUM_SITES, step=40):
    """
    Simple harmonic embedding of boundary sites into Z_360:
    site i -> position i * step mod 360.
    For step=40 and num_sites=9, we get 9 evenly spaced points.
    """
    positions = [(i * step) % N_CYCLE for i in range(num_sites)]
    return np.array(positions, dtype=int)


def cyclic_distance(a, b, N=N_CYCLE):
    """
    Minimal distance on a cycle of size N.
    """
    delta = abs(a - b)
    return min(delta, N - delta)


def boundary_distances(positions):
    """
    Build a matrix D_ij of cyclic distances between boundary sites.
    """
    num = len(positions)
    D = np.zeros((num, num), dtype=int)
    for i in range(num):
        for j in range(num):
            D[i, j] = cyclic_distance(positions[i], positions[j], N_CYCLE)
    return D


# ---------------------------------------------------
# 5. Holographic map: parent -> Yukawa matrices
# ---------------------------------------------------

def holographic_kernel(dist_matrix, lambd):
    """
    Build a simple exponential distance kernel:
        K_ij = exp(-lambda * d_ij)
    """
    return np.exp(-lambd * dist_matrix)


def build_sector_yukawa(freqs, amps, dist_matrix, lambd, sector_phase_shift=0.0):
    """
    Construct a toy Yukawa matrix for one sector as:
        Y_ij = (sum_n |a_n|^2) * exp(-lambda * d_ij) * phase_ij
    where phase_ij encodes a simple sector-dependent pattern derived
    from the parent phases.
    This keeps things simple but aligned with the parent.
    """
    weights = np.abs(amps) ** 2
    # effective overall scale from parent spectrum
    parent_scale = np.sum(weights)

    K = holographic_kernel(dist_matrix, lambd)

    # Build a simple coherent phase pattern on the boundary
    # using one representative phase from the parent.
    # This is just to make Y complex and sector-distinct.
    base_phase = np.angle(amps[0]) + sector_phase_shift
    num = dist_matrix.shape[0]
    phases = np.zeros((num, num), dtype=np.complex128)
    for i in range(num):
        for j in range(num):
            # Use a linear function of (i-j) as a crude phase pattern.
            phi_ij = base_phase * (i - j)
            phases[i, j] = np.exp(1j * phi_ij)

    Y = parent_scale * K * phases
    return Y


def build_all_sectors(freqs, amps, dist_matrix):
    """
    Build Yukawa-like matrices for four sectors:
    up, down, charged lepton, neutrino.
    Each gets its own decay constant and phase shift.
    """
    sectors = {}

    # Decay constants: stronger alignment for heavier hierarchy.
    lambd_u = 0.03
    lambd_d = 0.02
    lambd_e = 0.025
    lambd_nu = 0.005

    sectors["up"] = build_sector_yukawa(freqs, amps, dist_matrix, lambd_u, sector_phase_shift=0.0)
    sectors["down"] = build_sector_yukawa(freqs, amps, dist_matrix, lambd_d, sector_phase_shift=np.pi / 6.0)
    sectors["charged_lepton"] = build_sector_yukawa(freqs, amps, dist_matrix, lambd_e, sector_phase_shift=np.pi / 3.0)
    sectors["neutrino"] = build_sector_yukawa(freqs, amps, dist_matrix, lambd_nu, sector_phase_shift=np.pi / 2.0)

    return sectors


# ---------------------------------------------------
# 6. Simple diagnostics / "test run"
# ---------------------------------------------------

def summarize_matrix(name, M):
    """
    Print a compact summary: shape, singular values, and a few entries.
    """
    print(f"--- {name} ---")
    print("shape:", M.shape)
    # singular values
    svals = np.linalg.svd(M, compute_uv=False)
    print("singular values (approx):", np.round(svals, 4))
    # first 3x3 block
    print("top-left 3x3 block (real part):")
    print(np.round(M.real[:3, :3], 4))
    print("top-left 3x3 block (imag part):")
    print(np.round(M.imag[:3, :3], 4))
    print()


def run_pipeline():
    """
    Run the full aligned pipeline:
      1. Build parent state |Psi>
      2. Apply selection operator S^
      3. Embed boundary sites and build distance matrix
      4. Build sector Yukawa matrices holographically
      5. Print quick diagnostics
    """
    print("=== Building parent state |Psi> on Z_360 with triadic closure ===")
    freqs, amps = build_parent_state(gamma=0.02)
    print("Active parent frequencies (divisors of 360):", freqs)
    print("Number of modes:", len(freqs))
    print()

    print("=== Applying Selection Operator S^ = C^360 B^ P^phi ===")
    freqs_sel, amps_sel = apply_selection_operator(freqs, amps, alpha=0.7)
    print("Norm after selection:", np.linalg.norm(amps_sel))
    print()

    print("=== Embedding boundary / flavor sites into Z_360 ===")
    positions = build_embedding(num_sites=NUM_SITES, step=40)
    print("Boundary positions (mod 360):", positions)
    D = boundary_distances(positions)
    print("Boundary distance matrix D_ij (in units on Z_360):")
    print(D)
    print()

    print("=== Building holographic Yukawa-like matrices for each sector ===")
    sectors = build_all_sectors(freqs_sel, amps_sel, D)

    for name, Y in sectors.items():
        summarize_matrix(f"Y_{name}", Y)


if __name__ == "__main__":
    run_pipeline()
