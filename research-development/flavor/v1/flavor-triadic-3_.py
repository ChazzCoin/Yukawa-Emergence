#!/usr/bin/env python3
import numpy as np
from numpy.linalg import svd, solve, cond, pinv
import pandas as pd


# =========================
# Global config
# =========================

class Config:
    v = 246.0        # GeV
    mu0 = 1.0e12     # GeV
    mu_EW = 91.1876  # GeV
    Lambda_Maj = 1.0e14  # GeV (overall heavy Majorana scale)

    # Alignment scale: Fibonacci / 360
    # kappa = 360/89, eps = 1/kappa = 89/360
    kappa = 360.0 / 89.0
    eps = 1.0 / kappa

    seed = 12345  # overwritten per run

    t0 = np.log(mu0)
    t1 = np.log(mu_EW)
    dt = -0.01  # log-scale step size

    # Higgs quartic (approx EW value, treated constant here)
    lam = 0.13

    # Contextual alignment toggle
    use_contextual_kernel = True


# Global indices for light / heavy sites (9 = 3 + 6)
LIGHT = slice(0, 3)
HEAVY = slice(3, 9)

# Path to Base-60 harmonic table (parent structure)
# NOTE: update this path on your machine if needed.
BASE60_TABLE_PATH = "/Users/chazzromeo/Downloads/Base-60_Harmonic_Table_with_Human-Readable_Breakdown.xlsx"
_base60_df = None  # cached table

# Caches for sector triads / exponents
_sector_triads_cache = None
_sector_exponents_cache = None
_allowed_digits_cache = None

# Sector ordering and theta dimension
SECTORS = ["up", "down", "lepton", "neutrino"]
N_SECTORS = len(SECTORS)
N_GEN = 3


THETA_SECTORS = ["up", "down", "lepton", "neutrino"]
THETA_DIM = len(THETA_SECTORS) * N_GEN   # now 12

# Pi-vortex growth constants in D360 encoding
PI_D360 = 377.0 / 120.0      # ≈ 3.1416
EPS_PI  = 120.0 / 377.0      # ≈ 0.318 (vortex damping factor)


# Names for each observable in the same order as make_observables()
observable_names = [
    # mass ratios
    "m_c/m_t",
    "m_u/m_t",
    "m_s/m_b",
    "m_d/m_b",
    "m_mu/m_tau",
    "m_e/m_tau",
    # CKM angles
    "theta12_q (rad)",
    "theta23_q (rad)",
    "theta13_q (rad)",
    # PMNS angles
    "theta12_l (rad)",
    "theta23_l (rad)",
    "theta13_l (rad)",
    # neutrino splittings
    "Delta m2_21 (eV^2)",
    "Delta m2_31 (eV^2)",
]


# experimental targets (rough)
x_exp = np.array([
    # mass ratios
    0.007,    # m_c/m_t
    1e-5,     # m_u/m_t
    0.02,     # m_s/m_b
    0.001,    # m_d/m_b
    0.06,     # m_mu/m_tau
    0.0003,   # m_e/m_tau
    # CKM angles (rad)
    0.226, 0.041, 0.0035,
    # PMNS angles (rad)
    0.59, 0.84, 0.15,
    # Delta m^2 (eV^2)
    7.4e-5, 2.5e-3
])

sigma = np.array([
    0.5 * x_exp[0], 0.5 * x_exp[1], 0.5 * x_exp[2], 0.5 * x_exp[3],
    0.5 * x_exp[4], 0.5 * x_exp[5],
    0.1 * x_exp[6], 0.1 * x_exp[7], 0.1 * x_exp[8],
    0.1 * x_exp[9], 0.1 * x_exp[10], 0.1 * x_exp[11],
    0.3 * x_exp[12], 0.3 * x_exp[13]
])


# =========================
# Utility
# =========================

def has_bad(x: np.ndarray) -> bool:
    """Check for NaN or Inf in an array."""
    return np.any(np.isnan(x)) or np.any(np.isinf(x))


def chi2(observed, expected, sigma_arr):
    return np.sum(((observed - expected) / sigma_arr) ** 2)


# =========================
# Base-60 parent structure → sector triads → generation patterns
# =========================

K_DISTANCES = [1, 2, 3, 4, 5, 6, 8]

def kernel_params_to_weights(phi):
    """
    Map unconstrained real params phi_d to [0,1] via logistic sigmoids,
    approximating 0/1 selection for each distance.
    """
    phi = np.asarray(phi, float)
    assert phi.size == len(K_DISTANCES)
    return 1.0 / (1.0 + np.exp(-phi))  # σ(phi)


def get_base60_table() -> pd.DataFrame:
    """Load and cache the Base-60 harmonic table from disk."""
    global _base60_df
    if _base60_df is None:
        _base60_df = pd.read_excel(BASE60_TABLE_PATH, sheet_name=0)
    return _base60_df

def C360_exponents(e_raw):
    """
    A360 harmonic closure on exponent triads:
      C360 ∘ B ∘ Pϕ  acting on exponent-vector e_raw.

    Here:
      - Pϕ: project onto triadic line (E, E+Δ, E+2Δ)
      - B: minimal-misalignment “rounding” to integer step
      - C360: identification of integer step with D360 lattice
    """
    e_tri = project_to_triad(e_raw)   # Pϕ / B in exponent space
    e_h = project_to_D360(e_tri)      # D360 step locking
    return e_h

def divisors(N):
    return [d for d in range(1, N+1) if N % d == 0]

def C360_distances(N=9):
    """
    D360-based allowed distances for N-site chain.
    For N=9, distances run 1..8. We keep those that
    survive both D360 and D9 logic, with the 7-gap.
    """
    # Distances available on the N-chain
    d_chain = list(range(1, N))
    # D360 harmonics we want to sample
    D360 = divisors(360)
    allowed = [d for d in d_chain if d in D360]

    # Enforce the known 7-gap for N=9
    if N == 9 and 7 in allowed:
        allowed.remove(7)

    return tuple(sorted(allowed))

def gap_penalty(theta: np.ndarray,
                cfg: Config,
                lambda_gap: float = 1.0) -> float:
    """
    Penalize distortion of inter-generation exponent gaps
    away from the Base-60 triad gaps:

        M_gap[θ] = λ_gap * Σ_sectors [ (Δ12 - Δ12^0)^2 + (Δ23 - Δ23^0)^2 ]

    where Δij = e_i - e_j, e = e_base + δ(θ).
    """
    theta = np.asarray(theta, dtype=float)
    _, base_exps, _ = get_sector_harmonic_data()
    deltas = theta_to_deltas(theta)

    total = 0.0
    for s in SECTORS:
        e0 = np.array(base_exps[s], dtype=float)
        d  = np.array(deltas[s],    dtype=float)
        e  = e0 + d

        # current gaps
        d12 = e[0] - e[1]
        d23 = e[1] - e[2]

        # base gaps
        d12_0 = e0[0] - e0[1]
        d23_0 = e0[1] - e0[2]

        total += (d12 - d12_0)**2 + (d23 - d23_0)**2

    return lambda_gap * float(total)

def find_allowed_triads(df: pd.DataFrame):
    """
    Find all triads (n,2n,3n) with:
      - 'Allowed harmonic' status for each digit,
      - 3n < 60.
    """
    allowed = df[df["Harmonic Status"] == "Allowed harmonic"]["Base-60 Digit"].tolist()
    triads = []
    for n in allowed:
        tri = [n, 2 * n, 3 * n]
        if tri[-1] < 60 and all(t in allowed for t in tri):
            triads.append(tri)
    return triads, sorted(allowed)


def derive_sector_triads() -> dict:
    """
    Choose sector-dependent triads using Base-60 symbolic meanings.

      up-type quarks:    triad [1,2,3]
      down-type quarks:  triad [2,4,6]
      charged leptons:   triad [3,6,9]
      neutrinos:         triad [12,24,36]
    """
    df = get_base60_table()
    triads, allowed_digits = find_allowed_triads(df)

    triad_set = {tuple(t) for t in triads}

    sector_triads = {
        "up": [1, 2, 3],
        "down": [2, 4, 6],
        "lepton": [3, 6, 9],
        "neutrino": [12, 24, 36],
    }

    # Sanity check: ensure these are valid allowed triads
    for name, tri in sector_triads.items():
        if tuple(tri) not in triad_set:
            raise ValueError(f"Sector {name} triad {tri} is not an allowed (n,2n,3n) triad.")

    return sector_triads, allowed_digits


def exponents_from_triads(tri: list, allowed_digits: list) -> list:
    """
    Map a Base-60 triad to eps-exponents using a harmonic-depth rule.

    Define a "harmonic depth" index hd(d) for each allowed digit using
    its position in the sorted allowed list. For a triad (d1,d2,d3) with
    d3 = max digit, set exponents

        e_g = hd(d_max) - hd(d_g)

    so that the heaviest generation (largest digit) has exponent 0 and
    lighter generations have positive integer exponents.
    """
    hd = {d: i for i, d in enumerate(sorted(allowed_digits))}
    dmax = max(tri)
    exponents = [hd[dmax] - hd[d] for d in tri]
    return exponents  # [e1,e2,e3] with e3 = 0 by construction


def _init_sector_harmonics():
    """Initialize cached sector triads and base exponents from Base-60 table."""
    global _sector_triads_cache, _sector_exponents_cache, _allowed_digits_cache
    if _sector_triads_cache is None:
        sector_triads, allowed_digits = derive_sector_triads()
        base_exps = {
            sector: exponents_from_triads(tri, allowed_digits)
            for sector, tri in sector_triads.items()
        }
        _sector_triads_cache = sector_triads
        _sector_exponents_cache = base_exps
        _allowed_digits_cache = allowed_digits


def get_sector_harmonic_data():
    """Return (sector_triads, base_exponents, allowed_digits) from cache."""
    _init_sector_harmonics()
    return _sector_triads_cache, _sector_exponents_cache, _allowed_digits_cache

def theta_scale_penalty(theta: np.ndarray,
                        lambda_scale: float = 1.0) -> float:
    theta = np.asarray(theta, dtype=float)
    return lambda_scale * float(np.dot(theta, theta))


def theta_to_deltas(theta):
    theta = np.asarray(theta, dtype=float)
    if theta.size != THETA_DIM:
        raise ValueError(f"theta must have length {THETA_DIM}, got {theta.size}")
    deltas = {s: [0.0, 0.0, 0.0] for s in SECTORS}
    idx = 0
    for s in THETA_SECTORS:
        deltas[s] = theta[idx:idx + N_GEN].tolist()
        idx += N_GEN
    return deltas

def clamp_exponent(e, e_min=-10, e_max=+10):
    """
    Prevent eps**e from overflowing.
    e_min = -10 corresponds to eps^{-10} ≈ (1/0.247)^10 ≈ 1e6 (safe)
    e_max =  +10 corresponds to eps^{+10} ≈ 1e-6 (safe)

    Adjust the limits if needed, but [-10, +10] is a good starting band.
    """
    return max(e_min, min(e_max, e))
# ============================================================
# Triadic Projection Operator  P_tri
# ============================================================

def project_to_triad(exponents):
    """
    Project a 3-vector [e1,e2,e3] onto a triadic manifold:

        (e1, e2, e3) = (E, E+Δ, E+2Δ)     with Δ ∈ R

    This enforces 3n triadic closure, guarantees harmonic coherence,
    and prevents exponent drift into a collapsed or unphysical region.
    """
    e1, e2, e3 = exponents

    # Fit Δ by least squares under weights (1,1,1)
    # Solve min || (e1,e2,e3) - (E, E+Δ, E+2Δ) ||
    # Gives Δ = (e2 - e1 + e3 - e2)/2 = (e3 - e1)/2.
    Delta = 0.5 * (e3 - e1)

    # Then E = e1
    E = e1

    # Construct triadic vector
    e_tri = np.array([E, E + Delta, E + 2 * Delta], dtype=float)
    return e_tri


def project_to_D360(exponents):
    """
    Enforce D360 divisor alignment:
      - Exponents correspond to harmonic positions modulo repeating 360-lattice.
      - Practically: round to nearest allowed fractional harmonic step.

    Here we use the fundamental step 1/kappa = eps.

    So exponent space is restricted to integer multiples of 1.

    This keeps the structure aligned to the Base-60 / 360 harmonic cycle.
    """
    e = np.asarray(exponents, float)
    # Round each exponent to nearest integer (D360 harmonic step)
    return np.round(e)


def project_exponents_aligned(e_raw):
    """
    Soft alignment for evolution:
      e_raw → P_tri(e_raw)

    D360 locking will be penalized, not hard-enforced.
    """
    e_tri = project_to_triad(e_raw)   # only triadic, no rounding
    return e_tri



def sector_generation_patterns(cfg: Config, theta=None):
    eps = cfg.eps
    _, base_exps, _ = get_sector_harmonic_data()

    if theta is None:
        deltas = {s: [0.0, 0.0, 0.0] for s in SECTORS}
    else:
        deltas = theta_to_deltas(theta)

    patterns = {}
    for sector in SECTORS:
        e_base = np.array(base_exps[sector], dtype=float)
        d = np.array(deltas[sector], dtype=float)
        e_eff_raw = e_base + d
        e_eff = project_exponents_aligned(e_eff_raw)

        # TEMP DEBUG:
        print(f"[sector_generation_patterns] sector={sector}, "
              f"e_base={e_base}, d={d}, e_eff={e_eff}")

        patterns[sector] = [eps ** ee for ee in e_eff]

    return patterns



def harmonic_penalty(theta: np.ndarray,
                     cfg: Config,
                     w_tri: float = 1.0,
                     w_D360: float = 1.0) -> float:
    """
    A360 harmonic misalignment:

      M_harm[theta] = w_tri * Σ_sectors || e_tri(sector; θ) - e_tri_base(sector) ||^2
                    + w_D360 * Σ_sectors || e_tri(sector; θ) - round(e_tri(sector; θ)) ||^2

    where:
      - e_tri_base(sector) = P_tri(e_base) from the Base-60 table,
      - e_tri(sector; θ)   = P_tri(e_base + δ(θ)).

    This penalizes:
      • drifting away from the original Base-60 triads (triadic deformation),
      • drifting away from the D360 integer exponent lattice.
    """
    theta = np.asarray(theta, dtype=float)

    # Base exponents per sector from the Base-60 table
    _, base_exps, _ = get_sector_harmonic_data()
    deltas = theta_to_deltas(theta)

    penalty_tri = 0.0
    penalty_D360 = 0.0

    for sector in SECTORS:
        e0 = np.array(base_exps[sector], dtype=float)   # base exponents
        d  = np.array(deltas[sector],    dtype=float)   # θ-shifts

        # Effective exponents with θ deformation
        e_raw = e0 + d

        # Triadic projection for current and base state
        e_tri  = project_to_triad(e_raw)
        e_tri0 = project_to_triad(e0)

        diff_tri = e_tri - e_tri0
        penalty_tri += float(np.dot(diff_tri, diff_tri))

        # Distance from nearest D360 integer lattice
        e_round = np.round(e_tri)
        diff_D  = e_tri - e_round
        penalty_D360 += float(np.dot(diff_D, diff_D))

    return w_tri * penalty_tri + w_D360 * penalty_D360





# =========================
# Alignment kernel K (contextual)
# =========================
def build_alignment_kernel_parametric(eps: float, phi, N: int = 9) -> np.ndarray:
    """
    Parametric alignment kernel:
      K_ij = w_d * eps^{|i-j|} for d = |i-j| in K_DISTANCES,
      K_ij = 0 otherwise, K_ii = 1.

    Here w_d ∈ (0,1) is emergent via gradient descent.
    """
    weights = kernel_params_to_weights(phi)  # len=7
    w_map = {d: w for d, w in zip(K_DISTANCES, weights)}

    K = np.zeros((N, N), dtype=float)
    for i in range(N):
        for j in range(N):
            d = abs(i - j)
            if d == 0:
                K[i, j] = 1.0
            elif d in w_map:
                K[i, j] = w_map[d] * (eps ** d)
            else:
                K[i, j] = 0.0
    return K

def build_alignment_kernel(eps: float, N: int = 9,
                           allowed_distances=None) -> np.ndarray:

    """
    Build an NxN alignment kernel:
      K_ij = eps^{|i-j|} for |i-j| in allowed_distances,
      K_ij = 0 for other off-diagonals,
      K_ii = 1.

    Allowed distances are fixed by the D_360 divisor law, with some
    distances (e.g. 7 in the 9-site case) forbidden.
    """
    if allowed_distances is None:
        allowed_distances = C360_distances(N)

    allowed_distances = set(allowed_distances)
    K = np.zeros((N, N), dtype=float)
    for i in range(N):
        for j in range(N):
            d = abs(i - j)
            if d == 0:
                K[i, j] = 1.0
            elif d in allowed_distances:
                K[i, j] = eps ** d
            else:
                K[i, j] = 0.0
    return K


def eps_at_scale(t: float, cfg: Config) -> float:
    """
    Simple example of a scale-dependent epsilon.

    Maps t in [t1, t0] (EW → high scale) to a mild change in eps.
    """
    if not getattr(cfg, "use_contextual_kernel", False):
        return cfg.eps

    # Normalize t to x in [0,1]
    x = (t - cfg.t1) / (cfg.t0 - cfg.t1)
    x = max(0.0, min(1.0, x))
    # Interpolate between eps_low and eps_high
    eps_high = cfg.eps          # at high scale
    eps_low = cfg.eps * 0.8     # slightly smaller at EW
    return eps_low + (eps_high - eps_low) * x


def allowed_distances_at_scale(t: float, cfg: Config, N: int = 9):
    if not getattr(cfg, "use_contextual_kernel", False):
        return C360_distances(N)

    if t > np.log(1e10):
        # full C360 pattern
        return C360_distances(N)
    elif t > np.log(1e4):
        # shave off the outermost harmonics (example)
        d_full = list(C360_distances(N))
        return tuple(d for d in d_full if d <= 6)
    else:
        # ultra-local triad core
        return (1, 2, 3)


def build_alignment_kernel_contextual(cfg: Config, t: float, N: int = 9) -> np.ndarray:
    """
    Contextual alignment kernel K(t) in N-dimensional site space.
    """
    eps_t = eps_at_scale(t, cfg)
    allowed = allowed_distances_at_scale(t, cfg)
    return build_alignment_kernel(eps_t, N=N, allowed_distances=allowed)


def build_alignment_kernel_3(cfg: Config, t: float) -> np.ndarray:
    """
    3x3 contextual kernel for the effective 3-generation Yukawa / Weinberg sector.
    """
    return build_alignment_kernel_contextual(cfg, t, N=3)


# =========================
# Proto-matrices
# =========================

def random_complex_matrix(shape, rng):
    real = rng.normal(0.0, 1.0, size=shape)
    imag = rng.normal(0.0, 1.0, size=shape)
    return (real + 1j * imag) / np.sqrt(2.0)


def normalize_by_largest_singular_value(X: np.ndarray) -> np.ndarray:
    s = svd(X, compute_uv=False)
    s_max = np.max(s)
    if s_max == 0:
        return X
    return X / s_max


def random_weighted_proto(shape, rng, site_scales):
    """
    Gaussian proto-matrix with site-dependent variances:
      Var[X_ij] ~ site_scales[i] * site_scales[j].
    Then normalized so largest singular value = 1.
    """
    site_scales = np.asarray(site_scales, dtype=float)
    S = np.outer(site_scales, site_scales)
    real = rng.normal(0.0, S)
    imag = rng.normal(0.0, S)
    X = (real + 1j * imag) / np.sqrt(2.0)
    return normalize_by_largest_singular_value(X)


def build_site_scales_from_generations(gen_pattern):
    """
    gen_pattern: length-3 array [s1, s2, s3] for 'generation' (1,2,3).

    We assign site_scales[i] = gen_pattern[i % 3] for a 9-site chain.
    This enforces triadic repetition:
      sites (0,3,6)->s1, (1,4,7)->s2, (2,5,8)->s3.
    """
    gen_pattern = np.array(gen_pattern, dtype=float)
    scales = np.zeros(9, dtype=float)
    for i in range(9):
        g = i % 3
        scales[i] = gen_pattern[g]
    return scales


def generate_proto_matrices(cfg: Config, theta=None):
    """
    Generate proto Yukawa and Majorana matrices on the 9-site proto-flavor space,
    with sector-dependent structure baked in from the parent harmonic table.

    - For each sector (u,d,e,nu), we:
        * pick a Base-60 triad,
        * map it to exponents via a universal harmonic rule,
        * optionally shift exponents using theta,
        * build a 9-site scale profile by triadic repetition.

    - All randomness is in the complex Gaussian proto entries; hierarchy comes
      purely from these harmonic scale profiles and the alignment kernel.
    """
    rng = np.random.default_rng(cfg.seed)
    patterns = sector_generation_patterns(cfg, theta)

    site_scales_u = build_site_scales_from_generations(patterns["up"])
    site_scales_d = build_site_scales_from_generations(patterns["down"])
    site_scales_e = build_site_scales_from_generations(patterns["lepton"])
    site_scales_nu = build_site_scales_from_generations(patterns["neutrino"])

    # Draw weighted proto-matrices
    Yu0 = random_weighted_proto((9, 9), rng, site_scales_u)
    Yd0 = random_weighted_proto((9, 9), rng, site_scales_d)
    Ye0 = random_weighted_proto((9, 9), rng, site_scales_e)
    Ynu0 = random_weighted_proto((9, 9), rng, site_scales_nu)

    # Majorana proto: O(1) and symmetric, no extra site hierarchy yet
    M0 = random_complex_matrix((9, 9), rng)
    M0 = normalize_by_largest_singular_value(M0)
    M0 = 0.5 * (M0 + M0.T)

    return Yu0, Yd0, Ye0, Ynu0, M0


# =========================
# Alignment Φ: K ⊙ X
# =========================

def apply_alignment(K: np.ndarray, X: np.ndarray) -> np.ndarray:
    """Hadamard (elementwise) alignment: Φ(X) = K ⊙ X."""
    return K * X


def align_all(K, Yu0, Yd0, Ye0, Ynu0, M0):
    Yu9 = apply_alignment(K, Yu0)
    Yd9 = apply_alignment(K, Yd0)
    Ye9 = apply_alignment(K, Ye0)
    Ynu9 = apply_alignment(K, Ynu0)
    M9 = apply_alignment(K, M0)
    return Yu9, Yd9, Ye9, Ynu9, M9


def selection_operator(cfg: Config, proto_state, t: float = None):
    """
    Selection operator S^: apply contextual kernel K(t) to proto state.
    """
    if t is None:
        t = cfg.t0
    K9 = build_alignment_kernel_contextual(cfg, t, N=9)
    Yu0, Yd0, Ye0, Ynu0, M0 = proto_state
    Yu9, Yd9, Ye9, Ynu9, M9 = align_all(K9, Yu0, Yd0, Ye0, Ynu0, M0)
    return (Yu9, Yd9, Ye9, Ynu9, M9), K9


# =========================
# Schur complement 9→3
# =========================

def schur_9_to_3(Y9: np.ndarray, cond_tol: float = 1e12) -> np.ndarray:
    """
    Y9 is 9x9. Light sites: 0,1,2; heavy: 3..8.
    Effective 3x3 Yukawa via Schur complement:
      Y_eff = A - B D^{-1} B†.

    If D is ill-conditioned, uses pseudo-inverse.
    """
    A = Y9[LIGHT, LIGHT]
    B = Y9[LIGHT, HEAVY]
    D = Y9[HEAVY, HEAVY]

    if cond(D) > cond_tol:
        # fall back to pseudo-inverse
        D_inv = pinv(D)
        Y_eff = A - B @ D_inv @ B.conj().T
    else:
        X = solve(D, B.conj().T)  # D X = B†
        Y_eff = A - B @ X
    return Y_eff


def reduce_all_dirac(Yu9, Yd9, Ye9, Ynu9):
    Yu_eff = schur_9_to_3(Yu9)
    Yd_eff = schur_9_to_3(Yd9)
    Ye_eff = schur_9_to_3(Ye9)
    Ynu_eff = schur_9_to_3(Ynu9)
    return Yu_eff, Yd_eff, Ye_eff, Ynu_eff


# =========================
# Majorana sector: triadic projection 6→3
# =========================

def heavy_block(M9: np.ndarray) -> np.ndarray:
    """Extract 6x6 heavy block (sites 3..8, 0-based)."""
    return M9[HEAVY, HEAVY]


def triad_heavy_basis(Nh=6):
    """
    Build a 6x3 triadic basis in heavy space using DFT modes k = 0,1,2.
    Columns are normalized.

    k=0 is the uniform (democratic) mode, k=1,2 are the first two harmonics.
    """
    ks = np.array([0, 1, 2])
    i = np.arange(Nh)
    basis = []
    for k in ks:
        vec = np.exp(2j * np.pi * k * i / Nh)
        vec /= np.linalg.norm(vec)
        basis.append(vec)
    return np.stack(basis, axis=1)


def build_M_R_triadic(M9_aligned: np.ndarray,
                      Lambda_Maj: float) -> np.ndarray:
    """
    9x9 aligned Majorana → 6x6 heavy block → triadic 3x3 projection.

    M_R = Λ_Maj * B_H† M_H B_H, symmetrized.
    """
    M_H = heavy_block(M9_aligned)    # 6x6
    B_H = triad_heavy_basis(6)       # 6x3 (harmonic ks from Nh)
    M3 = B_H.conj().T @ M_H @ B_H    # 3x3
    M3 = 0.5 * (M3 + M3.T)           # enforce symmetry
    M_R = Lambda_Maj * M3
    return M_R


def seesaw_light_neutrinos(Ynu_eff: np.ndarray,
                           M_R: np.ndarray,
                           v: float,
                           cond_tol: float = 1e12) -> np.ndarray:
    """
    Type-I seesaw:
      m_D = v/√2 Ynu_eff,
      m_ν = - m_D M_R^{-1} m_D^T (symmetric 3x3, in GeV).
    """
    m_D = (v / np.sqrt(2.0)) * Ynu_eff
    if cond(M_R) > cond_tol:
        M_R_inv = pinv(M_R)
        m_nu = -m_D @ M_R_inv @ m_D.T
    else:
        X = solve(M_R, m_D.T)
        m_nu = -m_D @ X

    m_nu = 0.5 * (m_nu + m_nu.T)
    return m_nu


# =========================
# 1-loop Yukawa RGEs (g frozen)
# + Weinberg operator RGE
# =========================

def beta_Yukawas(Yu, Yd, Ye, Ynu, g1, g2, g3):
    """
    1-loop SM Yukawa RGEs (in matrix form), with fixed gauge couplings.
    """
    if any(has_bad(M) for M in (Yu, Yd, Ye, Ynu)):
        Z = np.zeros_like(Yu)
        return Z, Z, Z, Z

    Yu_dagYu = Yu.conj().T @ Yu
    Yd_dagYd = Yd.conj().T @ Yd
    Ye_dagYe = Ye.conj().T @ Ye
    Ynu_dagYnu = Ynu.conj().T @ Ynu

    T = np.trace(3 * Yu_dagYu + 3 * Yd_dagYd + Ye_dagYe)

    factor_u = T - (17 / 20 * g1 ** 2 + 9 / 4 * g2 ** 2 + 8 * g3 ** 2)
    factor_d = T - (1 / 4 * g1 ** 2 + 9 / 4 * g2 ** 2 + 8 * g3 ** 2)
    factor_e = T - (9 / 4 * g1 ** 2 + 9 / 4 * g2 ** 2)
    factor_nu = T - (9 / 20 * g1 ** 2 + 9 / 4 * g2 ** 2)

    dYu = Yu * factor_u + (3 / 2) * (Yu @ Yu_dagYu - Yd @ (Yd_dagYd @ Yu))
    dYd = Yd * factor_d + (3 / 2) * (Yd @ Yd_dagYd - Yu @ (Yu_dagYu @ Yd))
    dYe = Ye * factor_e + (3 / 2) * (Ye @ Ye_dagYe)
    dYnu = Ynu * factor_nu + (3 / 2) * (Ynu @ Ynu_dagYnu - Ye @ (Ye_dagYe @ Ynu))

    dYu /= (16 * np.pi ** 2)
    dYd /= (16 * np.pi ** 2)
    dYe /= (16 * np.pi ** 2)
    dYnu /= (16 * np.pi ** 2)

    return dYu, dYd, dYe, dYnu


def beta_kappa_L(kappa_L, Yu, Ye, g2, lam):
    """
    16π² dκ_L/dt = (-3 g2² + 2λ + 6 Tr(Yu†Yu)) κ_L
                   - 3/2 (Ye†Ye κ_L + κ_L (Ye†Ye)^T).

    We treat λ as constant, and ignore g1,g3 in this operator.
    """
    if any(has_bad(M) for M in (kappa_L, Yu, Ye)):
        return np.zeros_like(kappa_L)

    Yu_dagYu = Yu.conj().T @ Yu
    Ye_dagYe = Ye.conj().T @ Ye
    T_u = np.trace(Yu_dagYu)

    pref = (-3 * g2 ** 2 + 2 * lam + 6 * T_u)
    term1 = pref * kappa_L
    term2 = -1.5 * (Ye_dagYe @ kappa_L + kappa_L @ Ye_dagYe.T.conj())

    dkappa = (term1 + term2) / (16 * np.pi ** 2)
    return dkappa


def rk4_step_full(Yu, Yd, Ye, Ynu, kappa_L, g1, g2, g3, lam, dt):
    """
    RK4 step evolving Yukawas + κ_L with fixed (g1,g2,g3,lam).
    """
    # k1
    dYu1, dYd1, dYe1, dYnu1 = beta_Yukawas(Yu, Yd, Ye, Ynu, g1, g2, g3)
    dkL1 = beta_kappa_L(kappa_L, Yu, Ye, g2, lam)

    # k2
    Yu2 = Yu + 0.5 * dt * dYu1
    Yd2 = Yd + 0.5 * dt * dYd1
    Ye2 = Ye + 0.5 * dt * dYe1
    Ynu2 = Ynu + 0.5 * dt * dYnu1
    kL2 = kappa_L + 0.5 * dt * dkL1

    dYu2, dYd2, dYe2, dYnu2 = beta_Yukawas(Yu2, Yd2, Ye2, Ynu2, g1, g2, g3)
    dkL2 = beta_kappa_L(kL2, Yu2, Ye2, g2, lam)

    # k3
    Yu3 = Yu + 0.5 * dt * dYu2
    Yd3 = Yd + 0.5 * dt * dYd2
    Ye3 = Ye + 0.5 * dt * dYe2
    Ynu3 = Ynu + 0.5 * dt * dYnu2
    kL3 = kappa_L + 0.5 * dt * dkL2

    dYu3, dYd3, dYe3, dYnu3 = beta_Yukawas(Yu3, Yd3, Ye3, Ynu3, g1, g2, g3)
    dkL3 = beta_kappa_L(kL3, Yu3, Ye3, g2, lam)

    # k4
    Yu4 = Yu + dt * dYu3
    Yd4 = Yd + dt * dYd3
    Ye4 = Ye + dt * dYe3
    Ynu4 = Ynu + dt * dYnu3
    kL4 = kappa_L + dt * dkL3

    dYu4, dYd4, dYe4, dYnu4 = beta_Yukawas(Yu4, Yd4, Ye4, Ynu4, g1, g2, g3)
    dkL4 = beta_kappa_L(kL4, Yu4, Ye4, g2, lam)

    Yu_next = Yu + (dt / 6.0) * (dYu1 + 2 * dYu2 + 2 * dYu3 + dYu4)
    Yd_next = Yd + (dt / 6.0) * (dYd1 + 2 * dYd2 + 2 * dYd3 + dYd4)
    Ye_next = Ye + (dt / 6.0) * (dYe1 + 2 * dYe2 + 2 * dYe3 + dYe4)
    Ynu_next = Ynu + (dt / 6.0) * (dYnu1 + 2 * dYnu2 + 2 * dYnu3 + dYnu4)
    kL_next = kappa_L + (dt / 6.0) * (dkL1 + 2 * dkL2 + 2 * dkL3 + dkL4)

    return Yu_next, Yd_next, Ye_next, Ynu_next, kL_next


def run_RGE_full(Yu0, Yd0, Ye0, Ynu0, kappa_L0,
                 g1_const, g2_const, g3_const,
                 cfg: Config):
    """
    Evolve Yukawas and Weinberg operator from μ0 down to μ_EW
    with fixed gauge couplings.

    Contextual selection is applied at the proto 9x9 level only;
    the effective 3x3 theory is allowed to develop mixing freely.
    """
    Yu, Yd, Ye, Ynu = Yu0.copy(), Yd0.copy(), Ye0.copy(), Ynu0.copy()
    kappa_L = kappa_L0.copy()
    g1, g2, g3 = g1_const, g2_const, g3_const
    lam = cfg.lam

    t = cfg.t0
    step = 0
    while (cfg.dt < 0 and t > cfg.t1) or (cfg.dt > 0 and t < cfg.t1):
        step += 1

        Yu, Yd, Ye, Ynu, kappa_L = rk4_step_full(
            Yu, Yd, Ye, Ynu, kappa_L, g1, g2, g3, lam, cfg.dt
        )
        t += cfg.dt

        # NO 3x3 K3-projection here

        if any(has_bad(M) for M in (Yu, Yd, Ye, Ynu, kappa_L)):
            print(f"Warning: NaN/Inf detected at RGE step {step}, halting evolution.")
            break

    return Yu, Yd, Ye, Ynu, kappa_L, g1, g2, g3

# =========================
# Diagonalization and angles
# =========================

def diag_dirac_Y(Y: np.ndarray, v: float):
    """
    SVD for Dirac Yukawa:
      Y = U_L diag(s) U_R†,  masses = v/√2 * s.
    """
    U_L, s, U_Rh = svd(Y)
    masses = (v / np.sqrt(2.0)) * s
    return U_L, np.diag(s), U_Rh, masses


def takagi_symmetric(m: np.ndarray):
    """
    Takagi factorization via SVD for complex symmetric 3x3:
      m = U diag(s) U^T, with s ≥ 0.
    """
    U, s, Vh = svd(m)
    return U, s


def diagonalize_all(Yu, Yd, Ye, mnu, v):
    UuL, Yu_diag, UuR, mu = diag_dirac_Y(Yu, v)
    UdL, Yd_diag, UdR, md = diag_dirac_Y(Yd, v)
    UlL, Ye_diag, UlR, me = diag_dirac_Y(Ye, v)

    U_nu, mnu_vals = takagi_symmetric(mnu)
    mnu_masses = mnu_vals  # in GeV

    Vckm = UuL.conj().T @ UdL
    Vpmns = UlL.conj().T @ U_nu

    return mu, md, me, mnu_masses, Vckm, Vpmns


def extract_angles_and_phase(V: np.ndarray):
    """
    Extract approximate mixing angles (in radians) and Dirac phase
    from a 3x3 unitary matrix V, assuming a PDG-like parameterization.
    """
    s13 = np.abs(V[0, 2])
    theta13 = np.arcsin(np.clip(s13, 0.0, 1.0))

    s12 = np.abs(V[0, 1])
    c12 = np.abs(V[0, 0])
    theta12 = np.arctan2(s12, c12)

    s23 = np.abs(V[1, 2])
    c23 = np.abs(V[2, 2])
    theta23 = np.arctan2(s23, c23)

    J = np.imag(V[0, 0] * V[1, 1] * np.conj(V[0, 1]) * np.conj(V[1, 0]))
    denom = (np.sin(2 * theta12) * np.sin(2 * theta23) *
             np.sin(2 * theta13) * np.cos(theta13))
    if np.abs(denom) < 1e-12:
        delta = 0.0
    else:
        x = 8 * J / denom
        x = np.clip(x, -1.0, 1.0)
        delta = np.arcsin(x)

    return theta12, theta23, theta13, delta


def neutrino_splittings(mnu_masses: np.ndarray):
    """
    Compute Δm²_21 and Δm²_31 in GeV² from the (non-negative) Takagi singular values.
    """
    m_sorted = np.sort(mnu_masses)
    m1, m2, m3 = m_sorted
    dm2_21 = m2 ** 2 - m1 ** 2
    dm2_31 = m3 ** 2 - m1 ** 2
    return dm2_21, dm2_31  # GeV^2


# =========================
# χ² and observables
# =========================

def make_observables(res):
    mu, md, me, mnu = res["mu"], res["md"], res["me"], res["mnu"]
    th12_q, th23_q, th13_q = res["th_q"]
    th12_l, th23_l, th13_l = res["th_l"]
    dm2_21, dm2_31 = res["dm2_eV2"]

    # sort ascending so index 2 is heaviest
    mu_sorted = np.sort(mu)
    md_sorted = np.sort(md)
    me_sorted = np.sort(me)

    obs = []

    # mass ratios
    obs.append(mu_sorted[1] / mu_sorted[2])  # m_c/m_t
    obs.append(mu_sorted[0] / mu_sorted[2])  # m_u/m_t
    obs.append(md_sorted[1] / md_sorted[2])  # m_s/m_b
    obs.append(md_sorted[0] / md_sorted[2])  # m_d/m_b
    obs.append(me_sorted[1] / me_sorted[2])  # m_mu/m_tau
    obs.append(me_sorted[0] / me_sorted[2])  # m_e/m_tau

    # CKM
    obs.append(th12_q)
    obs.append(th23_q)
    obs.append(th13_q)

    # PMNS
    obs.append(th12_l)
    obs.append(th23_l)
    obs.append(th13_l)

    # neutrino splittings (eV²)
    obs.append(dm2_21)
    obs.append(dm2_31)

    return np.array(obs)


def chi2_from_res(res, w_mix: float = 3.0, w_dm2: float = 2.0):
    x_th = make_observables(res)
    diffs = x_th - x_exp
    # copies of sigma
    sigma_eff = sigma.copy()

    # indices: [0..5 mass ratios, 6..8 CKM, 9..11 PMNS, 12..13 Δm²]
    mix_idx = np.arange(6, 12)
    dm2_idx = np.arange(12, 14)

    # effectively reduce sigma → increases weight
    sigma_eff[mix_idx] /= np.sqrt(w_mix)
    sigma_eff[dm2_idx] /= np.sqrt(w_dm2)

    return np.sum((diffs / sigma_eff) ** 2)



def chi2_breakdown(res):
    """
    Return per-observable χ² contributions as a list of dicts:
      {"name", "theory", "exp", "sigma", "chi2_i"}
    """
    x_th = make_observables(res)
    diffs = x_th - x_exp
    chi2_i = (diffs / sigma) ** 2

    breakdown = []
    for name, th, exp, sig, c2 in zip(observable_names, x_th, x_exp, sigma, chi2_i):
        breakdown.append({
            "name": name,
            "theory": th,
            "exp": exp,
            "sigma": sig,
            "chi2_i": c2,
        })
    return breakdown


def rescale_yukawa_sector(Y, v, m_target_heaviest):
    """
    Rescale Y so that the heaviest mass eigenvalue (v/√2 * max singular value)
    matches m_target_heaviest. Returns (Y_rescaled, alpha).
    """
    U_L, s, U_Rh = svd(Y)
    m_current = (v / np.sqrt(2.0)) * np.max(s)
    if m_current == 0:
        return Y, 1.0
    alpha = m_target_heaviest / m_current
    return alpha * Y, alpha


# =========================
# Manifestation operator (extract + score)
# =========================

def manifestation_operator(res):
    """
    Apply 'measurement': extract observables and compute χ².
    """
    obs = make_observables(res)
    return obs, chi2_from_res(res)


# =========================
# run_pipeline with theta (evolution + selection + manifestation)
# =========================

def run_pipeline(seed: int,
                 cfg: Config,
                 use_RGE: bool = True,
                 theta=None):
    """
    Full pipeline with harmonically derived, sector-dependent proto structure.
    """
    cfg.seed = seed

    # Proto state from Base-60-driven patterns (+theta)
    proto_state = generate_proto_matrices(cfg, theta=theta)

    # Selection at high scale
    (Yu9, Yd9, Ye9, Ynu9, M9), K9 = selection_operator(cfg, proto_state, t=cfg.t0)

    # Schur 9→3 for Dirac Yukawas
    Yu_eff, Yd_eff, Ye_eff, Ynu_eff = reduce_all_dirac(Yu9, Yd9, Ye9, Ynu9)

    # M_R (triadic heavy projection)
    M_R = build_M_R_triadic(M9, cfg.Lambda_Maj)

    # Seesaw at μ0 → mν(μ0) in GeV
    m_nu_0 = seesaw_light_neutrinos(Ynu_eff, M_R, cfg.v)

    # Weinberg operator κ_L(μ0) (dimensionful, GeV^-1)
    kappa_L_0 = (2.0 / cfg.v ** 2) * m_nu_0

    # RG evolution
    g1_0, g2_0, g3_0 = 0.46, 0.63, 0.88
    if use_RGE:
        Yu_EW, Yd_EW, Ye_EW, Ynu_EW, kappa_L_EW, g1_EW, g2_EW, g3_EW = run_RGE_full(
            Yu_eff, Yd_eff, Ye_eff, Ynu_eff, kappa_L_0, g1_0, g2_0, g3_0, cfg
        )
        # reconstruct mν(μ_EW) from κ_L(μ_EW)
        m_nu_EW = 0.5 * cfg.v ** 2 * kappa_L_EW
    else:
        Yu_EW, Yd_EW, Ye_EW, Ynu_EW = Yu_eff, Yd_eff, Ye_eff, Ynu_eff
        m_nu_EW = m_nu_0
        g1_EW, g2_EW, g3_EW = g1_0, g2_0, g3_0

    # Rescale sectors to fix heavy masses
    m_t_target = 173.0
    m_b_target = 4.18
    m_tau_target = 1.777

    Yu_EW, alpha_u = rescale_yukawa_sector(Yu_EW, cfg.v, m_t_target)
    Yd_EW, alpha_d = rescale_yukawa_sector(Yd_EW, cfg.v, m_b_target)
    Ye_EW, alpha_e = rescale_yukawa_sector(Ye_EW, cfg.v, m_tau_target)

    # Diagonalize at μ_EW
    mu, md, me, mnu_masses, Vckm, Vpmns = diagonalize_all(
        Yu_EW, Yd_EW, Ye_EW, m_nu_EW, cfg.v
    )

    # Angles, Δm²
    th12_q, th23_q, th13_q, delta_q = extract_angles_and_phase(Vckm)
    th12_l, th23_l, th13_l, delta_l = extract_angles_and_phase(Vpmns)
    dm2_21_GeV2, dm2_31_GeV2 = neutrino_splittings(mnu_masses)
    dm2_21_eV2 = dm2_21_GeV2 * 1e18
    dm2_31_eV2 = dm2_31_GeV2 * 1e18

    res = {
        "mu": mu,
        "md": md,
        "me": me,
        "mnu": mnu_masses,  # GeV
        "Vckm": Vckm,
        "Vpmns": Vpmns,
        "th_q": (th12_q, th23_q, th13_q),
        "delta_q": delta_q,
        "th_l": (th12_l, th23_l, th13_l),
        "delta_l": delta_l,
        "dm2_GeV2": (dm2_21_GeV2, dm2_31_GeV2),
        "dm2_eV2": (dm2_21_eV2, dm2_31_eV2),
        "alphas": (alpha_u, alpha_d, alpha_e),
        "g_EW": (g1_EW, g2_EW, g3_EW),
        "K9": K9,
        "theta": theta,
    }

    res["chi2"] = chi2_from_res(res)
    return res


# =========================
# Misalignment functional & gradient flow
# =========================
def misalignment(theta: np.ndarray,
                 cfg: Config,
                 seed: int = 0,
                 use_RGE: bool = True,
                 alpha_phys: float = 1.0,
                 alpha_harm: float = 1.0,
                 alpha_gap: float = 1.0,
                 w_tri: float = 1.0,
                 w_D360: float = 1.0,
                 lambda_scale: float = 1.0,
                 lambda_gap: float = 1.0,
                 w_mix: float = 3.0,
                 w_dm2: float = 2.0) -> float:
    """
    Full A360 misalignment functional:

    M[θ] = α_phys * χ²_phys(θ; w_mix, w_dm2)
         + α_harm * M_harm(θ)
         + α_gap  * M_gap(θ)
         + M_scale(θ)
    """
    theta = np.asarray(theta, dtype=float)

    # Run pipeline with current θ
    res = run_pipeline(seed, cfg, use_RGE=use_RGE, theta=theta)

    # Physical χ² with enhanced weight on mixing & Δm²
    M_phys = chi2_from_res(res, w_mix=w_mix, w_dm2=w_dm2)

    # Harmonic penalties (triad alignment + D360 lattice)
    M_harm = harmonic_penalty(theta, cfg, w_tri=w_tri, w_D360=w_D360)

    # Gap penalty to prevent hyper-hierarchical exponents
    M_gap = gap_penalty(theta, cfg, lambda_gap=lambda_gap)

    # θ-scale regularization
    M_scale = theta_scale_penalty(theta, lambda_scale=lambda_scale)

    return alpha_phys * M_phys + alpha_harm * M_harm + alpha_gap * M_gap + M_scale




def grad_misalignment(theta, cfg: Config, seed: int = 0,
                      eps_theta: float = 1e-2,
                      use_RGE: bool = True,
                      alpha_phys: float = 1.0,
                      alpha_harm: float = 1.0,
                      w_tri: float = 1.0,
                      w_D360: float = 1.0,
                      lambda_scale: float = 0.1):
    """
    Finite-difference gradient of the full misalignment functional M[theta].
    """
    theta = np.asarray(theta, dtype=float)
    base = misalignment(theta, cfg, seed=seed, use_RGE=use_RGE,
                        alpha_phys=alpha_phys, alpha_harm=alpha_harm,
                        w_tri=w_tri, w_D360=w_D360,
                        lambda_scale=lambda_scale)
    grad = np.zeros_like(theta)
    for i in range(theta.size):
        th2 = theta.copy()
        th2[i] += eps_theta
        grad[i] = (
            misalignment(th2, cfg, seed=seed, use_RGE=use_RGE,
                         alpha_phys=alpha_phys, alpha_harm=alpha_harm,
                         w_tri=w_tri, w_D360=w_D360,
                         lambda_scale=lambda_scale)
            - base
        ) / eps_theta
    return grad



def optimize_alignment(cfg: Config,
                       seed: int = 0,
                       n_steps: int = 10,
                       eta0: float = 0.1,
                       eps_theta: float = 1e-2,
                       use_RGE: bool = True,
                       alpha_phys: float = 1.0,
                       alpha_harm: float = 1.0,
                       lambda_scale: float = 0.1):
    """
    A360 pi-vortex evolution in theta-space:
      theta_{n+1} = theta_n - eta_n ∂M/∂theta
      eta_n = eta0 * EPS_PI**n
    """
    theta = np.zeros(THETA_DIM)
    history = []

    for step in range(n_steps):
        eta_n = eta0 * (EPS_PI ** step)  # π-vortex growth schedule

        chi2_val = misalignment(theta, cfg, seed=seed,
                                use_RGE=use_RGE,
                                alpha_phys=alpha_phys,
                                alpha_harm=alpha_harm,
                                lambda_scale=lambda_scale)
        history.append((step, eta_n, chi2_val, theta.copy()))
        print(f"[opt step {step}] eta = {eta_n:.3g}, M = {chi2_val:.3g}")

        grad = grad_misalignment(theta, cfg, seed=seed,
                                 eps_theta=eps_theta, use_RGE=use_RGE)

        theta = theta - eta_n * grad
        theta = np.clip(theta, -3.0, +3.0)

    res_final = run_pipeline(seed, cfg, use_RGE=use_RGE, theta=theta)
    return res_final, theta, history



# =========================
# Scan driver
# =========================
if __name__ == "__main__":
    cfg = Config()
    N_seeds = 10

    print("=== Seed scan with Base-60-driven patterns (theta=0) ===")
    all_results = []
    chi2_vals = []

    for seed in range(N_seeds):
        r = run_pipeline(seed, cfg, use_RGE=True, theta=None)
        all_results.append(r)
        chi2_vals.append(r["chi2"])
        print(f"seed {seed}: chi2 = {r['chi2']:.3g}")

    best_idx = int(np.argmin(chi2_vals))
    best = all_results[best_idx]

    print("\nBest seed (theta=0 base state):", best_idx)
    print("chi2 =", best["chi2"])
    print("Up masses (GeV):      ", best["mu"])
    print("Down masses (GeV):    ", best["md"])
    print("Lepton masses (GeV):  ", best["me"])
    print("Neutrino masses (GeV):", best["mnu"])
    print("Neutrino masses (eV): ", best["mnu"] * 1e9)
    print("Δm² (eV²):            ", best["dm2_eV2"])
    print("CKM angles (rad):     ", best["th_q"], "δq:", best["delta_q"])
    print("PMNS angles (rad):    ", best["th_l"], "δℓ:", best["delta_l"])

    print("\n=== χ² breakdown for best seed (theta=0) ===")
    breakdown = chi2_breakdown(best)
    for entry in breakdown:
        name = entry["name"]
        th = entry["theory"]
        exp = entry["exp"]
        sig = entry["sigma"]
        c2 = entry["chi2_i"]
        pull = (th - exp) / sig
        print(f"{name:20s}  th = {th: .4e},  exp = {exp: .4e},  "
              f"sigma = {sig: .4e},  pull = {pull: .2f},  chi2_i = {c2: .2f}")

    # Now turn on the derivative harmonic layer: optimize theta
    print("\n=== Gradient-flow alignment in theta-space ===")
    # You can start with a small n_steps because each step is expensive.
    # If runtime is okay, increase n_steps (e.g. 5–10).
    res_opt, theta_opt, history = optimize_alignment(
        cfg,
        seed=best_idx,
        n_steps=3,       # try 3 first; increase if it's fast enough
        eps_theta=1e-2,  # finite-difference step
        use_RGE=True,
    )

    print("\nOptimized theta (exponent shifts):")
    print(theta_opt)
    print("Optimized chi2:", res_opt["chi2"])

    print("\n=== χ² breakdown after theta-alignment ===")
    breakdown_opt = chi2_breakdown(res_opt)
    for entry in breakdown_opt:
        name = entry["name"]
        th = entry["theory"]
        exp = entry["exp"]
        sig = entry["sigma"]
        c2 = entry["chi2_i"]
        pull = (th - exp) / sig
        print(f"{name:20s}  th = {th: .4e},  exp = {exp: .4e},  "
              f"sigma = {sig: .4e},  pull = {pull: .2f},  chi2_i = {c2: .2f}")

"""


"""