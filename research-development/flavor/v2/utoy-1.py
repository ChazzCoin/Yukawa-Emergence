import itertools
import math

import numpy as np

# N_SITES = 1080

def dynamic_light_heavy_sites(Y9):
    """
    Determine LIGHT_SITES and HEAVY_SITES dynamically
    from the magnitude profile of the proto-Yukawa.

    Strategy:
      - Compute per-site norm: sum_j |Y9[i,j]|
      - Take 3 smallest as LIGHT
      - Remaining 6 as HEAVY
    """
    mags = np.sum(np.abs(Y9), axis=1)   # length-9

    order = np.argsort(mags)
    light = order[:3].tolist()
    heavy = order[3:].tolist()

    return light, heavy
def build_proto_basis(evecs, triad):
    """
    Build a 360×9 proto-flavor basis matrix B
    using the triad + next 6 eigenvectors.
    """
    # triad eigenvectors
    v1, v2, v3 = [evecs[:, i] for i in triad]

    # next 6 modes sorted by eigenvalue
    # (exclude triad and zero mode)
    all_modes = [i for i in range(len(evecs[0])) if i not in triad and i != 0]
    next6 = sorted(all_modes)[:6]  # simple and effective

    extra = [evecs[:, i] for i in next6]

    # stack 9 basis vectors
    B = np.stack([v1, v2, v3] + extra, axis=1)  # shape (360,9)

    # Orthonormalize
    Q, _ = np.linalg.qr(B)
    return Q[:, :9]
def project_to_9x9(Y360, B):
    """
    Project full 360×360 Yukawa onto 9×9 proto-flavor basis.
    """
    return B.conj().T @ Y360 @ B
def sector_exponents_from_Q(q_vals, perm):
    """
    q_vals = [1,2,3] (sorted harmonic strengths)
    perm = tuple of 3 indices giving sector ordering
    """
    return q_vals[list(perm)]
def apply_alignment(Y9, kappa, forbidden_d):
    N = Y9.shape[0]
    K = np.zeros((N,N), dtype=complex)

    for i in range(N):
        for j in range(N):
            d = min(abs(i-j), N - abs(i-j))
            if d == forbidden_d:
                K[i,j] = 0
            else:
                K[i,j] = kappa**d
    return K * Y9
def dynamic_schur(Y9):
    LIGHT, HEAVY = dynamic_light_heavy_sites(Y9)

    A = Y9[np.ix_(LIGHT, LIGHT)]
    B = Y9[np.ix_(LIGHT, HEAVY)]
    C = Y9[np.ix_(HEAVY, LIGHT)]
    D = Y9[np.ix_(HEAVY, HEAVY)]

    Dinv = np.linalg.pinv(D)
    return A - B @ Dinv @ C

# ----------------------------------------------------------------------
# 1. Misalignment functional M[theta] and its gradient
# ----------------------------------------------------------------------

def misalignment_energy(theta, w6=1.0, w5=1.0, J=None):
    """
    Compute total misalignment energy M[theta].

    theta: array shape (N,)
    J: optional coupling matrix shape (N,N); if None, J_ij = 1/N.
    """
    theta = theta.reshape(-1)
    N = theta.shape[0]
    if J is None:
        # uniform coupling J_ij = 1/N (not critical in this toy)
        J = np.ones((N, N), dtype=float) / N

    # pairwise differences Δ_ij
    # we can vectorize using broadcasting
    dtheta = theta[:, None] - theta[None, :]  # shape (N, N)

    E6 = 1.0 - np.cos(6.0 * dtheta)
    E5 = 1.0 - np.cos(5.0 * dtheta)

    M = 0.5 * np.sum(J * (w6 * E6 + w5 * E5))  # 1/2 to avoid double-count
    return M

def misalignment_grad(theta, w6=1.0, w5=1.0, J=None):
    """
    Gradient dM/dtheta_i.

    d/dθ_i (1 - cos(kΔ_ij)) = k sin(kΔ_ij), where Δ_ij = θ_i - θ_j.

    So

        ∂M/∂θ_i = sum_j J_ij [ w6*6 sin(6Δ_ij) + w5*5 sin(5Δ_ij) ].
    """
    theta = theta.reshape(-1)
    N = theta.shape[0]
    if J is None:
        J = np.ones((N, N), dtype=float) / N

    dtheta = theta[:, None] - theta[None, :]  # shape (N, N)
    term6 = 6.0 * np.sin(6.0 * dtheta)
    term5 = 5.0 * np.sin(5.0 * dtheta)

    # sum over j: (J_ij * (w6*term6 + w5*term5))
    grad = np.sum(J * (w6 * term6 + w5 * term5), axis=1)
    return grad

def relax_phases(N=200, n_steps=500, eta=0.01, w6=1.0, w5=1.0, random_seed=42):
    """
    Perform gradient descent on M[theta] starting from random phases.
    Returns final theta and a history of energies.
    """
    rng = np.random.default_rng(random_seed)
    theta = rng.uniform(0.0, 2.0 * math.pi, size=N)
    J = None  # uniform couplings

    energy_hist = []

    for step in range(n_steps):
        M = misalignment_energy(theta, w6=w6, w5=w5, J=J)
        energy_hist.append(M)
        grad = misalignment_grad(theta, w6=w6, w5=w5, J=J)
        theta -= eta * grad
        # wrap back into [0, 2π)
        theta = np.mod(theta, 2.0 * math.pi)

    return theta, np.array(energy_hist)

def build_emergent_adjacency(theta, w6=1.0, w5=1.0, keep_fraction=0.1):
    """
    From a relaxed configuration theta, build an emergent adjacency A_int.

    For each pair (i,j):

        Δ_ij = θ_i - θ_j
        S_ij = w6*cos(6Δ_ij) + w5*cos(5Δ_ij)
        W_ij = max(0, S_ij)

    Then keep only the top 'keep_fraction' of W_ij (i<j) as edges.
    """
    theta = theta.reshape(-1)
    N = theta.shape[0]

    dtheta = theta[:, None] - theta[None, :]  # (N,N)
    S = w6 * np.cos(6.0 * dtheta) + w5 * np.cos(5.0 * dtheta)
    W = np.maximum(0.0, S)

    # Zero out diagonal
    np.fill_diagonal(W, 0.0)

    # Threshold
    # flatten upper triangle (i<j), pick top fraction
    iu, ju = np.triu_indices(N, k=1)
    weights = W[iu, ju]
    if keep_fraction <= 0.0:
        keep_fraction = 0.1
    n_edges = max(1, int(keep_fraction * weights.size))

    # indices of top weights
    top_idx = np.argpartition(weights, -n_edges)[-n_edges:]
    mask = np.zeros_like(weights, dtype=bool)
    mask[top_idx] = True

    # build adjacency
    A = np.zeros((N, N), dtype=float)
    A[iu[mask], ju[mask]] = 1.0
    A[ju[mask], iu[mask]] = 1.0

    return A

def laplacian(A):
    d = np.sum(A, axis=1)
    return np.diag(d) - A


###############################################################################
# 3. Harmonics, triad selection, and R,Q on H_gen  (B1,B2)
###############################################################################


def harmonic_strengths(evals, evecs, phi, divisors=(1, 2, 3)):
    """Overlap of each Laplacian mode with e^{-i d φ} for a few d."""
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
    return np.array(strengths)


def select_triad(evals, H):
    """Select 3 modes (≠ zero) with largest harmonic strength as generation triad."""
    scores = np.sum(H[1:], axis=1)
    idx = np.arange(1, len(evals))
    top3_rel = np.argsort(scores)[-3:]
    triad = np.sort(idx[top3_rel])
    return triad

def build_R_from_triad(evecs, phi, triad):
    """Base-360 internal symmetry R on 3D triad subspace (B1)."""
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
    """Integer charges q_j from triad harmonic strengths (B2, F2)."""
    A = np.sum(H[triad], axis=1)
    order = np.argsort(-A)
    q = np.empty(3, dtype=int)
    q[order] = np.array([1, 2, 3])
    Q = np.diag(q.astype(float))
    return Q, q

# ================================================================
#  COHERENCE KERNEL (γ) — HERE γ = 0 (COHERENT LIMIT)
# ================================================================

def build_kernel_gamma(gamma: float, n_sites, forbidden_d: int = 2):
    """
    Toeplitz kernel on the ring, but we'll use gamma=0
    so K becomes almost trivial (except for forbidden distance, if kept).
    """
    K = np.zeros((n_sites, n_sites), dtype=float)
    for i in range(n_sites):
        for j in range(n_sites):
            if i == j:
                K[i, j] = 1.0
            else:
                d = min(abs(i - j), n_sites - abs(i - j))
                if d == forbidden_d:
                    K[i, j] = 0.0
                else:
                    K[i, j] = math.exp(-gamma * d)
    return K


# ================================================================
#  9→3 SCHUR REDUCTION
# ================================================================

# def schur_9to3(Y9: np.ndarray) -> np.ndarray:
#     ls = LIGHT_SITES
#     hs = HEAVY_SITES
#     A = Y9[np.ix_(ls, ls)]
#     B = Y9[np.ix_(ls, hs)]
#     C = Y9[np.ix_(hs, ls)]
#     D = Y9[np.ix_(hs, hs)]
#     Dinv = np.linalg.pinv(D)
#     Y_eff = A - B @ Dinv @ C
#     return Y_eff + 1e-9 * np.eye(3)

# ================================================================
#  RESONANT PHASE WHEELS
# ================================================================

def build_phase_matrix(phi_site, n_sites):
    P = np.zeros((n_sites, n_sites), dtype=complex)
    for i in range(n_sites):
        for j in range(n_sites):
            P[i, j] = np.exp(1j * (phi_site[i] - phi_site[j]))
    return P

# ================================================================
#  EXPONENT HIERARCHIES (X,Y shifts)
# ================================================================
def generation_index(i: int) -> int:
    return i % 3
def site_scales(base_exp, shifts, n_sites):
    """
    base_exp: [a,b,c]
    shifts: (X,Y) applied to 2nd & 3rd gen exponents
      e0 = base[0]
      e1 = base[1] + X
      e2 = base[2] + Y

    s_g = KAPPA^e_g, with KAPPA ≈ 0.24
    """
    X, Y = shifts
    eff = np.array([
        base_exp[0],
        base_exp[1] + X,
        base_exp[2] + Y
    ], dtype=float)

    KAPPA = 0.24
    s_gen = np.power(KAPPA, eff)

    s = np.zeros(n_sites, dtype=float)
    for i in range(n_sites):
        s[i] = s_gen[generation_index(i)]
    return s

# ================================================================
#  SECTOR YUKAWAS
# ================================================================
def build_phase_profile(A: float, B: float):
    """
    φ_g = A + B*g, g=0,1,2
    """
    return np.array([A, A + B, A + 2*B], dtype=float)

def build_site_phases(phi_gen, n_sites):
    phi_site = np.zeros(n_sites, dtype=float)
    for i in range(n_sites):
        phi_site[i] = phi_gen[generation_index(i)]
    return phi_site
def build_Y_sector(A, B, base_exp, shifts, gamma, n_sites):
    phi_gen = build_phase_profile(A, B)
    phi_site = build_site_phases(phi_gen, n_sites)
    P = build_phase_matrix(phi_site, n_sites)

    s = site_scales(base_exp, shifts, n_sites)
    mag = np.outer(s, s)

    Y0 = mag * P
    K = build_kernel_gamma(gamma, n_sites)
    Y = K * Y0

    # SVD normalization
    _, sv, _ = np.linalg.svd(Y)
    if sv[0] != 0:
        Y /= sv[0]
    return Y

def build_geometric_unitary(gen_vecs: np.ndarray, region_list):
    """
    Given:
      gen_vecs: shape (N_sites, 3) = eigenvectors for the generation triad
      region_list: list of 3 index arrays (sites in each region for this sector)

    Construct 3 vectors in generation space by summing gen_vecs over each region,
    then orthonormalize them to get a 3x3 unitary-ish matrix U_geom.
    """
    cols = []
    for inds in region_list:
        # sum over sites in this region
        v = np.sum(gen_vecs[inds, :], axis=0)
        cols.append(v)
    M = np.stack(cols, axis=1)  # shape (3,3) with each col a vector in generation space

    # QR decomposition to orthonormalize columns
    Q, R = np.linalg.qr(M)
    # Optional: enforce det(Q) ~ +1 by flipping a column sign if needed
    if np.linalg.det(Q) < 0:
        Q[:, 0] = -Q[:, 0]
    return Q  # unitary (up to numerical noise)

TARGETS = {
    "mu/mt":     2.2e-5,
    "mc/mt":     7.5e-3,
    "md/mb":     1.1e-3,
    "ms/mb":     2.2e-2,
    "me/mtau":   2.9e-4,
    "mmu/mtau":  0.059,
    "theta12_q": 0.227,
    "theta23_q": 0.041,
    "theta13_q": 0.0036,
    "theta12_l": 0.584,
    "theta23_l": 0.785,
    "theta13_l": 0.150,
}

def search_best_lepton_regions(
    gen_vecs,
    regions,
    U_geom_u, U_geom_d,
    F_u, F_d, F_e, F_n,
    P_phi_12, P_phi_23, C_12,
    N_SOLAR=36, N_REACTOR=45
):
    """
    Brute-force search over all permutations of the 3 regions for
    charged leptons and neutrinos, keeping:
      - up/down geometry fixed (U_geom_u, U_geom_d),
      - masses F_s fixed,
      - golden/Cabibbo/neutrino-dressing operators fixed.

    Returns (best_assign_e, best_assign_nu, best_chi2, best_results),
    where best_results includes the mixing matrices and angles.
    """
    R0, R1, R2 = regions
    region_list = [R0, R1, R2]
    perms = list(itertools.permutations(range(3)))

    best_chi2 = None
    best_assign_e = None
    best_assign_nu = None
    best_dat = None

    for pe in perms:
        for pn in perms:
            assign_e  = [region_list[i] for i in pe]
            assign_nu = [region_list[i] for i in pn]

            U_geom_e  = build_geometric_unitary(gen_vecs, assign_e)
            U_geom_nu = build_geometric_unitary(gen_vecs, assign_nu)

            U_geom = {
                "u":  U_geom_u,
                "d":  U_geom_d,
                "e":  U_geom_e,
                "nu": U_geom_nu,
            }

            sector_bases = build_sector_bases(
                P_phi_12, P_phi_23, C_12,
                U_geom,
                use_neutrino_dressing=True,
                N_SOLAR=N_SOLAR,
                N_REACTOR=N_REACTOR
            )

            U_L_u,  U_R_u  = sector_bases["u"]
            U_L_d,  U_R_d  = sector_bases["d"]
            U_L_e,  U_R_e  = sector_bases["e"]
            U_L_nu, U_R_nu = sector_bases["nu"]

            # Mixing matrices
            V_ckm  = mixing_matrix(U_L_u, U_L_d)
            U_pmns = mixing_matrix(U_L_e, U_L_nu)

            theta12_q, theta23_q, theta13_q = mixing_angles_from_U(V_ckm)
            theta12_l, theta23_l, theta13_l = mixing_angles_from_U(U_pmns)

            # Mass ratios from F_s (unchanged per iteration)
            mu_mt, mc_mt   = mass_ratios(F_u)
            md_mb, ms_mb   = mass_ratios(F_d)
            me_mt, mmu_mt  = mass_ratios(F_e)

            obs = compute_observables(
                mu_mt, mc_mt, md_mb, ms_mb, me_mt, mmu_mt,
                theta12_q, theta23_q, theta13_q,
                theta12_l, theta23_l, theta13_l
            )
            chi2_val, chi2_details = chi2(obs, TARGETS)

            if (best_chi2 is None) or (chi2_val < best_chi2):
                best_chi2 = chi2_val
                best_assign_e = pe
                best_assign_nu = pn
                best_dat = {
                    "chi2": chi2_val,
                    "chi2_details": chi2_details,
                    "V_ckm": V_ckm,
                    "U_pmns": U_pmns,
                    "angles_q": (theta12_q, theta23_q, theta13_q),
                    "angles_l": (theta12_l, theta23_l, theta13_l),
                }

    return best_assign_e, best_assign_nu, best_chi2, best_dat

def rot12(theta: float) -> np.ndarray:
    c = math.cos(theta)
    s = math.sin(theta)
    return np.array([
        [ c,  s, 0.0],
        [-s,  c, 0.0],
        [0.0, 0.0, 1.0]
    ], dtype=complex)

def rot23(theta: float) -> np.ndarray:
    c = math.cos(theta)
    s = math.sin(theta)
    return np.array([
        [1.0, 0.0, 0.0],
        [0.0,  c,  s],
        [0.0, -s,  c]
    ], dtype=complex)

def rot13(theta: float) -> np.ndarray:
    c = math.cos(theta)
    s = math.sin(theta)
    return np.array([
        [ c, 0.0,  s],
        [0.0, 1.0, 0.0],
        [-s, 0.0,  c]
    ], dtype=complex)

def build_sector_bases(P_phi_12, P_phi_23, C_12,
                       U_geom,
                       use_neutrino_dressing=True,
                       N_SOLAR=36, N_REACTOR=45):
    """
    Build the sector L/R bases for quarks + leptons.
    This updated version gracefully handles None inputs for
    P_phi_12, P_phi_23, and C_12 by replacing them with identity.
    """

    I3 = np.eye(3, dtype=complex)

    # ------------------------------------------------------------
    # FIX: Replace None with identity matrices
    # ------------------------------------------------------------
    if P_phi_12 is None:
        P_phi_12 = I3
    if P_phi_23 is None:
        P_phi_23 = I3
    if C_12 is None:
        C_12 = I3

    # ------------------------------------------------------------
    # Quark sectors
    # ------------------------------------------------------------
    U_L_u = U_geom["u"] @ P_phi_12
    U_L_d = U_geom["d"] @ P_phi_12 @ C_12

    # ------------------------------------------------------------
    # Charged leptons
    # ------------------------------------------------------------
    U_L_e = U_geom["e"]

    # ------------------------------------------------------------
    # Neutrino dressing
    # ------------------------------------------------------------
    if use_neutrino_dressing:
        theta12_nu = 2.0 * math.pi / float(N_SOLAR)
        theta13_nu = 2.0 * math.pi / float(N_REACTOR)

        U_L_nu = (
            U_geom["nu"]
            @ rot12(theta12_nu)
            @ P_phi_23
            @ rot13(theta13_nu)
        )
    else:
        U_L_nu = U_geom["nu"] @ P_phi_23

    # R-bases are trivial in this framework
    U_R_u  = I3
    U_R_d  = I3
    U_R_e  = I3
    U_R_nu = I3

    return {
        "u":  (U_L_u,  U_R_u),
        "d":  (U_L_d,  U_R_d),
        "e":  (U_L_e,  U_R_e),
        "nu": (U_L_nu, U_R_nu),
    }

def mixing_matrix(U_L_up: np.ndarray, U_L_down: np.ndarray) -> np.ndarray:
    return U_L_up.conj().T @ U_L_down

def mixing_angles_from_U(U: np.ndarray):
    U_abs = np.abs(U)
    s13 = U_abs[0, 2]
    c13 = math.sqrt(max(0.0, 1.0 - s13**2))
    if c13 < 1e-12:
        return 0.0, 0.0, math.pi / 2.0
    s12 = U_abs[0, 1] / c13
    s23 = U_abs[1, 2] / c13
    s12 = max(-1.0, min(1.0, s12))
    s23 = max(-1.0, min(1.0, s23))
    theta12 = math.asin(s12)
    theta23 = math.asin(s23)
    theta13 = math.asin(s13)
    return theta12, theta23, theta13

def mass_ratios(F_s: np.ndarray):
    s_sorted = np.sort(F_s)
    m1, m2, m3 = s_sorted
    return m1 / m3, m2 / m3

def compute_observables(mu_mt, mc_mt, md_mb, ms_mb, me_mt, mmu_mt,
                        theta12_q, theta23_q, theta13_q,
                        theta12_l, theta23_l, theta13_l):
    return {
        "mu/mt":      mu_mt,
        "mc/mt":      mc_mt,
        "md/mb":      md_mb,
        "ms/mb":      ms_mb,
        "me/mtau":    me_mt,
        "mmu/mtau":   mmu_mt,
        "theta12_q":  theta12_q,
        "theta23_q":  theta23_q,
        "theta13_q":  theta13_q,
        "theta12_l":  theta12_l,
        "theta23_l":  theta23_l,
        "theta13_l":  theta13_l,
    }

def chi2(observables, targets):
    """
    Robust, warning-free χ² evaluator.
    Converts all observables into clean real scalars.
    Handles tiny imaginary parts and numpy dtypes gracefully.
    """

    def to_scalar(x):
        """
        Convert any input (scalar, numpy scalar, array, complex)
        to a *real float*, safely stripping tiny imaginary parts.
        """
        arr = np.asarray(x).reshape(-1)
        val = arr[0]

        # Strip tiny imag parts (floating-point noise)
        if np.iscomplex(val):
            val = np.real(val)

        # If still nan or inf, convert to zero to avoid χ² explosions
        if not np.isfinite(val):
            return 0.0

        return float(val)

    chi2_total = 0.0
    details = []

    ratio_keys = [
        "mu/mt", "mc/mt",
        "md/mb", "ms/mb",
        "me/mtau", "mmu/mtau"
    ]

    angle_keys = [
        "theta12_q", "theta23_q", "theta13_q",
        "theta12_l", "theta23_l", "theta13_l"
    ]

    # ============================================================
    # MASS RATIOS — evaluated in log10-space
    # ============================================================
    for k in ratio_keys:
        m = to_scalar(observables[k])
        t = to_scalar(targets[k])

        # Skip impossible or zero values
        if m <= 0 or t <= 0:
            contrib = 1000.0   # penalize pathological cases
            chi2_total += contrib
            details.append((k, m, t, contrib))
            continue

        logm = math.log10(m)
        logt = math.log10(t)
        sigma_log = 0.30

        contrib = ((logm - logt) / sigma_log) ** 2
        chi2_total += contrib
        details.append((k, m, t, contrib))

    # ============================================================
    # MIXING ANGLES — evaluated in linear space
    # ============================================================
    for k in angle_keys:
        m = to_scalar(observables[k])
        t = to_scalar(targets[k])

        sigma = 0.20
        contrib = ((m - t) / sigma) ** 2

        chi2_total += contrib
        details.append((k, m, t, contrib))

    return chi2_total, details


# ===============================================================
# Helper: Build phase-based geometric regions
# ===============================================================

def build_regions(phi, n_regions=3):
    """
    Partition phase field φ into 3 disjoint masks (R0,R1,R2).
    These are used by build_geometric_unitary().
    """
    phi = np.asarray(phi).reshape(-1)
    N = len(phi)

    # Normalize phases into [0, 2π)
    phi_norm = (phi - np.min(phi)) % (2*np.pi)

    # Sort indices by φ
    idx_sorted = np.argsort(phi_norm)

    block_size = N // n_regions
    regions = []

    for r in range(n_regions):
        start = r * block_size
        end   = (r + 1) * block_size if r < n_regions - 1 else N

        mask = np.zeros(N, dtype=bool)
        mask[idx_sorted[start:end]] = True
        regions.append(mask)

    return regions


# ===============================================================
# Unified Yukawa Builder Using R, Q, K (Your True Alignment Theory)
# ===============================================================

def build_Y_sector_unified(phi, R, Q_vals, K, kappa=0.24):
    """
    Unified Yukawa:
        magnitude from integer charges Q
        phase from φ and the diagonal phases of R
        alignment via Schur kernel K
    """

    N = len(phi)

    # ------------------------------------------------------------
    # FIX: R is a 3×3 diag in your code; convert to a length-3 vector
    # ------------------------------------------------------------
    R_vec = np.diag(R)  # extract diagonal entries e^{iθ_g}

    # magnitude from Q charges
    mag_gen = kappa ** Q_vals
    mag_site = np.array([mag_gen[i % 3] for i in range(N)])

    # geometric phase field
    phi_site = phi
    Y0 = np.zeros((N, N), dtype=complex)
    P = np.zeros((N, N), dtype=complex)

    for i in range(N):
        gi = i % 3
        for j in range(N):
            gj = j % 3

            geom_phase = np.exp(1j * (phi_site[i] - phi_site[j]))

            # 360-harmonic twist using diagonal R phases
            rtwist = R_vec[gi] * np.conj(R_vec[gj])

            P[i, j] = geom_phase * rtwist
            Y0[i, j] = mag_site[i] * mag_site[j]

    # final alignment: Φ(Y) = K ∘ (Y0 ∘ P)
    Y = K * (Y0 * P)
    return Y


# ===============================================================
# MAIN FULL PIPELINE — test_run()
# ===============================================================

def test_run(
    N=360,
    kappa=10,
    forbidden_d=2,
    seed=42
):
    print("\n=== UTM: Unified Toy Model Test Run ===")
    print("N = %d" % N)
    print("forbidden_d = %d" % forbidden_d)
    print("seed = %d" % seed)
    print("Kappa = %f" % kappa)

    # ------------------------------------------------------------
    # MODULE 1: Emergent Internal Geometry
    # ------------------------------------------------------------
    print("→ Relaxing misalignment phases...")
    phi, E = relax_phases(N=N, n_steps=800, eta=0.01, random_seed=seed)

    print("→ Building adjacency and Laplacian...")
    A_int = build_emergent_adjacency(phi)
    L_int = laplacian(A_int)

    evals, evecs = np.linalg.eigh(L_int)

    # ------------------------------------------------------------
    # MODULE 2: Generational Triad
    # ------------------------------------------------------------
    print("→ Extracting spectral triad...")
    H = harmonic_strengths(evals, evecs, phi)
    triad = select_triad(evals, H)

    # ------------------------------------------------------------
    # MODULE 3: Build R, Q
    # ------------------------------------------------------------
    print("→ Building harmonic operator R and charge operator Q...")
    R, k_vals = build_R_from_triad(evecs, phi, triad)
    Q, q_vals = build_Q_from_harmonics(H, triad)

    # sector-specific Q permutations
    # (example assignments — can be refined)
    Q_u  = q_vals
    Q_d  = q_vals[[1,0,2]]
    Q_e  = q_vals[[2,1,0]]
    Q_nu = q_vals

    # ------------------------------------------------------------
    # MODULE 4: Build Alignment Kernel Φ
    # ------------------------------------------------------------
    print("→ Building alignment kernel K...")
    K = build_kernel_gamma(gamma=0.0, n_sites=N, forbidden_d=forbidden_d)

    # ------------------------------------------------------------
    # MODULE 5: Build Sector Yukawas (Unified)
    # ------------------------------------------------------------
    print("→ Building unified Yukawa matrices...")

    Y_u9  = build_Y_sector_unified(phi, R, Q_u,  K, kappa)
    Y_d9  = build_Y_sector_unified(phi, R, Q_d,  K, kappa)
    Y_e9  = build_Y_sector_unified(phi, R, Q_e,  K, kappa)
    Y_nu9 = build_Y_sector_unified(phi, R, Q_nu, K, kappa)

    # 9→3 reduction
    Y_u  = dynamic_schur(Y_u9)
    Y_d  = dynamic_schur(Y_d9)
    Y_e  = dynamic_schur(Y_e9)
    Y_nu = dynamic_schur(Y_nu9)

    # ------------------------------------------------------------
    # MODULE 6: Regions → CKM & PMNS
    # ------------------------------------------------------------
    print("→ Building geometric mixing unitaries...")
    regions = build_regions(phi)

    # Up/down geometry (can be enhanced later)
    U_geom_u = build_geometric_unitary(evecs[:, triad], [regions[0], regions[1], regions[2]])
    U_geom_d = U_geom_u

    print("→ Computing CKM & PMNS via region permutations...")
    assign_e, assign_n, best_chi2_mix, data = search_best_lepton_regions(
        gen_vecs = evecs[:, triad],
        regions  = regions,
        U_geom_u = U_geom_u,
        U_geom_d = U_geom_d,
        F_u      = Y_u,
        F_d      = Y_d,
        F_e      = Y_e,
        F_n      = Y_nu,
        P_phi_12 = None,
        P_phi_23 = None,
        C_12     = None
    )

    V_ckm  = data["V_ckm"]
    U_pmns = data["U_pmns"]

    # ------------------------------------------------------------
    # MODULE 7: FULL χ² EVALUATION
    # ------------------------------------------------------------
    print("→ Computing SM χ² score...")

    mu_mt, mc_mt   = mass_ratios(Y_u)
    md_mb, ms_mb   = mass_ratios(Y_d)
    me_mt, mmu_mt  = mass_ratios(Y_e)

    obs = compute_observables(
        mu_mt, mc_mt,
        md_mb, ms_mb,
        me_mt, mmu_mt,
        *data["angles_q"],
        *data["angles_l"]
    )

    chi2_total, details = chi2(obs, TARGETS)

    # ------------------------------------------------------------
    # OUTPUT
    # ------------------------------------------------------------

    print("\n=== RESULTS ===")
    print("CKM:\n", V_ckm)
    print("PMNS:\n", U_pmns)
    print("Total χ² =", chi2_total)
    print("Mixing χ² =", best_chi2_mix)
    print("χ² Details:", details)

    return dict(
        chi2=chi2_total,
        V_ckm=V_ckm,
        U_pmns=U_pmns,
        Yukawas=dict(Y_u=Y_u, Y_d=Y_d, Y_e=Y_e, Y_nu=Y_nu),
        details=details
    )


# ===============================================================
# RUN IF EXECUTED DIRECTLY
# ===============================================================

if __name__ == "__main__":
    out = test_run()
    print("\nCompleted UTM Test Run.")