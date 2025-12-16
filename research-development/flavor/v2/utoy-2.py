###############################################################
# Unified Toy Model v2 (UTM_v2)
# Emergent Geometry → Proto-Flavor Projection → Alignment → SM
# Using Harmonic Neighborhood Projection (Resonant Operator)
###############################################################

import numpy as np
import math
import itertools
###############################################################
# 0. Utility
###############################################################

def normalize_vec(v):
    return v / np.linalg.norm(v)

###############################################################
# 1. Emergent Geometry (you already have these)
###############################################################




# N_SITES = 1080
#
# def dynamic_light_heavy_sites(Y9):
#     """
#     Determine LIGHT_SITES and HEAVY_SITES dynamically
#     from the magnitude profile of the proto-Yukawa.
#
#     Strategy:
#       - Compute per-site norm: sum_j |Y9[i,j]|
#       - Take 3 smallest as LIGHT
#       - Remaining 6 as HEAVY
#     """
#     mags = np.sum(np.abs(Y9), axis=1)   # length-9
#
#     order = np.argsort(mags)
#     light = order[:3].tolist()
#     heavy = order[3:].tolist()
#
#     return light, heavy
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

def sector_exponents_from_Q(q_vals, perm):
    """
    q_vals = [1,2,3] (sorted harmonic strengths)
    perm = tuple of 3 indices giving sector ordering
    """
    return q_vals[list(perm)]
# def apply_alignment(Y9, kappa, forbidden_d):
#     N = Y9.shape[0]
#     K = np.zeros((N,N), dtype=complex)
#
#     for i in range(N):
#         for j in range(N):
#             d = min(abs(i-j), N - abs(i-j))
#             if d == forbidden_d:
#                 K[i,j] = 0
#             else:
#                 K[i,j] = kappa**d
#     return K * Y9
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


def harmonic_strengths(evals, evecs, phi):
    """
    Compute harmonic strength for ALL eigenmodes.
    H[k] = |sum_i evecs[i,k] * exp(i phi_i)|^2
    """
    N = len(phi)
    H = np.zeros(N)

    phase = np.exp(1j * phi)

    for k in range(N):
        H[k] = np.abs(np.sum(evecs[:, k] * phase))**2

    return H


def select_triad(evals, H):
    """
    Select 3 generational modes using harmonic strength H.

    H is a 1D array of length N with H[k] = strength of mode k.
    We select the top 3 nonzero harmonic modes.
    """

    N = len(H)

    # Exclude trivial zero-mode k=0
    indices = np.argsort(-H)   # sorted by strength
    indices = [i for i in indices if i != 0]

    # Simply pick the top three modes
    triad = indices[:3]

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
    """
    Construct Q-values (hierarchy exponents) from harmonic strengths.

    Input:
        H: 1D array of length N (harmonic strengths for all modes)
        triad: list of 3 mode indices

    Output:
        Q_matrix: 3×3 diagonal matrix diag(q_vals)
        q_vals: array([q0, q1, q2]) where q_i reflects relative strength
    """

    # Extract strengths for the three generational modes
    triad_strengths = np.array([H[t] for t in triad], dtype=float)

    # Rank-order them (small → large)
    order = np.argsort(triad_strengths)

    # Assign integer-like hierarchy exponents: 1, 2, 3
    q_vals = np.zeros(3, dtype=float)
    q_vals[order[0]] = 1.0   # light gen
    q_vals[order[1]] = 2.0   # middle gen
    q_vals[order[2]] = 3.0   # heavy gen

    Q_matrix = np.diag(q_vals)
    return Q_matrix, q_vals

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
      region_list: list of 3 index arrays (boolean masks) defining regions

    Construct 3 vectors in generation space by summing gen_vecs over each region,
    then orthonormalize them to get a 3x3 unitary-ish matrix U_geom.
    """
    cols = []

    for inds in region_list:
        # inds is a boolean mask over sites
        if inds.dtype != bool:
            inds = np.asarray(inds, dtype=bool)

        # All sites in this region
        idxs = np.where(inds)[0]
        if idxs.size == 0:
            raise ValueError("Encountered an empty region in build_geometric_unitary.")

        # Sum over sites in this region
        v = np.sum(gen_vecs[idxs, :], axis=0)

        # Guard against destructive cancellation
        norm = np.linalg.norm(v)
        if norm < 1e-12:
            # Fallback: just take the first site's vector in this region
            v = gen_vecs[idxs[0], :]
            norm = np.linalg.norm(v)

        v = v / norm
        cols.append(v)

    M = np.stack(cols, axis=1)  # shape (3,3), columns are the 3 region vectors

    # QR decomposition to orthonormalize columns
    Q, _ = np.linalg.qr(M)

    # Enforce det(Q) ≈ +1 for a proper rotation (flip one column if needed)
    if np.linalg.det(Q) < 0:
        Q[:, 0] = -Q[:, 0]

    return Q  # 3x3 unitary (up to numerical noise)
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
            mu_mt, mc_mt = mass_ratios(F_u)  # now uses SVD of Yu3
            md_mb, ms_mb = mass_ratios(F_d)  # SVD of Yd3
            me_mt, mmu_mt = mass_ratios(F_e)  # SVD of Ye3

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
    """
    Compute (m1/m3, m2/m3) from either:
      - a 3×3 Yukawa matrix (use its singular values), or
      - a length-3 vector of masses/eigenvalues.

    This makes all existing call sites (Yu3, Yd3, Ye3, or vectors)
    behave as intended.
    """
    arr = np.asarray(F_s)

    if arr.ndim == 1:
        # Already an eigenvalue/mass vector
        s_sorted = np.sort(np.real_if_close(arr))
    elif arr.ndim == 2 and arr.shape == (3, 3):
        # Treat as Yukawa matrix: use singular values as "masses"
        _, s, _ = np.linalg.svd(arr)
        s_sorted = np.sort(np.real_if_close(s))
    else:
        raise ValueError(
            f"mass_ratios expects a length-3 vector or a 3x3 matrix, got shape {arr.shape}"
        )

    m1, m2, m3 = s_sorted
    if m3 == 0:
        # Degenerate / pathological: avoid division by zero
        return 0.0, 0.0

    return float(m1 / m3), float(m2 / m3)

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

# def build_regions(phi, n_regions=3):
#     """
#     Partition phase field φ into 3 disjoint masks (R0,R1,R2).
#     These are used by build_geometric_unitary().
#     """
#     phi = np.asarray(phi).reshape(-1)
#     N = len(phi)
#
#     # Normalize phases into [0, 2π)
#     phi_norm = (phi - np.min(phi)) % (2*np.pi)
#
#     # Sort indices by φ
#     idx_sorted = np.argsort(phi_norm)
#
#     block_size = N // n_regions
#     regions = []
#
#     for r in range(n_regions):
#         start = r * block_size
#         end   = (r + 1) * block_size if r < n_regions - 1 else N
#
#         mask = np.zeros(N, dtype=bool)
#         mask[idx_sorted[start:end]] = True
#         regions.append(mask)
#
#     return regions


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





###############################################################
# 2. Harmonic-Neighborhood Proto-Basis Construction (9 modes)
###############################################################

def build_proto_basis_harmonic(evecs, triad, H, num_extra=6):
    """
    Build a 360×9 proto-flavor basis using:
        - 3 triad eigenvectors
        - next-6 modes with highest harmonic strength
    This matches Resonant Model operator theory.
    """

    # Extract triad eigenvectors
    triad_vecs = [normalize_vec(evecs[:, i]) for i in triad]

    # Sort all modes by harmonic strength (descending)
    idx_sorted = np.argsort(-H)   # high → low

    # Remove triad indices
    idx_sorted = [i for i in idx_sorted if i not in triad]

    # Pick next 6 "harmonic neighbors"
    extra_modes = idx_sorted[:num_extra]

    extra_vecs = [normalize_vec(evecs[:, i]) for i in extra_modes]

    # Stack into matrix: shape (360, 9)
    B = np.stack(triad_vecs + extra_vecs, axis=1)

    # Orthonormalize basis
    Q, _ = np.linalg.qr(B)

    # Return first 9 orthonormal columns
    return Q[:, :9]


###############################################################
# 3. Project 360×360 → 9×9 Proto-Yukawa
###############################################################

def project_to_9x9(Y360, B):
    """
    Y360: full geometric Yukawa (360×360)
    B:   proto-basis (360×9)
    """
    return B.conj().T @ Y360 @ B


###############################################################
# 4. Build Full 360×360 Yukawa Seed From Q and R
###############################################################

# def build_Y360(phi, q_vals, R_vec, kappa):
#     """
#     Build the initial 360×360 Yukawa matrix BEFORE Φ-alignment and projection.
#     Structure:
#         Y_ij = kappa^(q[g(i)] + q[g(j)]) * exp(i(φ_i - φ_j)) * R_g(i) R_g(j)^*
#     """
#
#     N = len(phi)
#     Y = np.zeros((N, N), dtype=complex)
#
#     # generation index mapping
#     def g(i):
#         return i % 3
#
#     mag_gen = kappa ** q_vals
#     mag = np.array([mag_gen[g(i)] for i in range(N)])
#
#     for i in range(N):
#         for j in range(N):
#             geom_phase = np.exp(1j * (phi[i] - phi[j]))
#             rtwist     = R_vec[g(i)] * np.conj(R_vec[g(j)])
#             Y[i, j]    = mag[i] * mag[j] * geom_phase * rtwist
#
#     return Y


###############################################################
# 5. Alignment Operator Φ applied to 9×9
###############################################################

def build_alignment_kernel_9(kappa, forbidden_d):
    N = 9
    K = np.zeros((N, N), dtype=complex)

    for i in range(N):
        for j in range(N):
            d = min(abs(i-j), N - abs(i-j))
            if d == forbidden_d:
                K[i, j] = 0
            else:
                K[i, j] = kappa**d
    return K


def apply_alignment(Y9, kappa, forbidden_d):
    K = build_alignment_kernel_9(kappa, forbidden_d)
    return K * Y9


###############################################################
# 6. Dynamic LIGHT / HEAVY Selection + Schur Reduction
###############################################################

def dynamic_light_heavy_sites(Y9):
    mags = np.sum(np.abs(Y9), axis=1)
    order = np.argsort(mags)
    return order[:3].tolist(), order[3:].tolist()


def schur_reduce_dynamic(Y9):
    LIGHT, HEAVY = dynamic_light_heavy_sites(Y9)

    A = Y9[np.ix_(LIGHT, LIGHT)]
    B = Y9[np.ix_(LIGHT, HEAVY)]
    C = Y9[np.ix_(HEAVY, LIGHT)]
    D = Y9[np.ix_(HEAVY, HEAVY)]

    Dinv = np.linalg.pinv(D)
    Y3 = A - B @ Dinv @ C
    return Y3


def dynamic_generation_map(B, triad_indices=[0,1,2]):
    """
    Determine generation index per site using projection onto the
    triad basis inside the 9D proto-flavor space.

    For each of the N internal sites, we compute:
        g(i) = argmax(|B[i, triad_indices]|)
    where B is the 360×9 proto-basis matrix.

    Returns:
        g: integer array of shape (N,), values ∈ {0,1,2}
    """
    N = B.shape[0]
    g = np.zeros(N, dtype=int)

    for i in range(N):
        # Look only at components corresponding to the triad directions
        components = np.abs(B[i, triad_indices])
        g[i] = int(np.argmax(components))

    return g
def build_Y360(phi, q_vals, R_vec, kappa, g_map):
    """
    Full geometric Yukawa on the 360-site space.

    Structure:
        Y_ij = kappa^(q[g(i)] + q[g(j)])
               * exp(i(φ_i − φ_j))
               * R[g(i)] * R[g(j)]^*
    """
    N = len(phi)
    Y = np.zeros((N, N), dtype=complex)

    for i in range(N):
        gi = g_map[i]
        mag_i = kappa ** q_vals[gi]

        for j in range(N):
            gj = g_map[j]
            mag_j = kappa ** q_vals[gj]

            geom_phase = np.exp(1j * (phi[i] - phi[j]))
            rtwist     = R_vec[gi] * np.conj(R_vec[gj])

            Y[i, j] = mag_i * mag_j * geom_phase * rtwist

    return Y

###############################################################
# 7. Build Sector Yukawas Using Q Permutations
###############################################################

# Sector orderings (tunable)
PERM_U  = (0, 1, 2)
PERM_D  = (1, 2, 0)
PERM_E  = (2, 1, 0)
PERM_NU = (0, 2, 1)


def exponents_from_perm(q_vals, perm):
    return np.array([q_vals[i] for i in perm])


###############################################################
# 8. CKM / PMNS (use your existing functions)
###############################################################



def build_regions(phi):
    phi = np.asarray(phi)
    N = len(phi)
    idx = np.argsort(phi % (2*np.pi))
    thirds = N // 3
    regions = []
    for k in range(3):
        s = k * thirds
        e = (k+1)*thirds if k < 2 else N
        mask = np.zeros(N, dtype=bool)
        mask[idx[s:e]] = True
        regions.append(mask)
    return regions
def build_ckm(Yu3, Yd3, U_geom_u, U_geom_d):
    """
    Construct CKM matrix using both:
    - geometric quark unitaries (U_geom_u, U_geom_d)
    - Yukawa diagonalization unitaries

    CKM = (U_geom_u · U_L_u)† · (U_geom_d · U_L_d)
    """

    # Left singular vectors of Yukawas (analog of left-handed diagonalization)
    Uu, _, _ = np.linalg.svd(Yu3)
    Ud, _, _ = np.linalg.svd(Yd3)

    # Compose geometric and Yukawa rotations
    U_L_u = U_geom_u @ Uu
    U_L_d = U_geom_d @ Ud

    # CKM = misalignment between up and down left-handed bases
    CKM = U_L_u.conj().T @ U_L_d
    return CKM
def build_pmns(Ye3, Yn3, U_geom_e, U_geom_nu):
    Ue, _, _ = np.linalg.svd(Ye3)
    Uv, _, _ = np.linalg.svd(Yn3)
    U_L_e  = U_geom_e  @ Ue
    U_L_nu = U_geom_nu @ Uv
    return U_L_e.conj().T @ U_L_nu
###############################################################
# 9. Main Pipeline — UTM v2
###############################################################

def test_run(
    N=60,
    kappa=8,
    forbidden_d=2,
    seed=42,
    debug=False
):
    print("\n=== UTM v2: Unified Toy Model (Harmonic-Neighborhood) ===")
    print(f"N = {N}")
    print(f"kappa = {kappa}")
    print(f"forbidden_d = {forbidden_d}")
    print(f"seed = {seed}")

    # ------------------------------------------------------------
    # Step 1: Emergent geometry
    # ------------------------------------------------------------
    print("→ Relaxing phases...")
    phi, E = relax_phases(N=N, n_steps=1000, eta=0.01, random_seed=seed)

    print("→ Building adjacency + Laplacian...")
    A = build_emergent_adjacency(phi)
    L = laplacian(A)
    evals, evecs = np.linalg.eigh(L)

    # ------------------------------------------------------------
    # Step 2: Triad + harmonic strengths
    # ------------------------------------------------------------
    print("→ Extracting harmonic triad...")
    H = harmonic_strengths(evals, evecs, phi)
    triad = select_triad(evals, H)

    # ------------------------------------------------------------
    # Step 3: Build R and Q
    # ------------------------------------------------------------
    print("→ Building R and Q operators...")
    R, _ = build_R_from_triad(evecs, phi, triad)
    q_matrix, q_vals = build_Q_from_harmonics(H, triad)
    R_vec = np.diag(R)

    # Sector Q-permutations
    Q_u  = exponents_from_perm(q_vals, PERM_U)
    Q_d  = exponents_from_perm(q_vals, PERM_D)
    Q_e  = exponents_from_perm(q_vals, PERM_E)
    Q_nu = exponents_from_perm(q_vals, PERM_NU)
    # ------------------------------------------------------------
    # Step 5: Build proto basis (harmonic neighborhood)
    # ------------------------------------------------------------
    print("→ Building harmonic-neighborhood proto basis...")
    B = build_proto_basis_harmonic(evecs, triad, H)
    # ------------------------------------------------------------
    # Step 4: Build 360×360 Yukawa seed
    # ------------------------------------------------------------
    print("→ Building 360×360 Yukawa seeds...")
    # Create dynamic generation map from proto-basis
    g_map = dynamic_generation_map(B)

    Yu360 = build_Y360(phi, Q_u, R_vec, kappa, g_map)
    Yd360 = build_Y360(phi, Q_d, R_vec, kappa, g_map)
    Ye360 = build_Y360(phi, Q_e, R_vec, kappa, g_map)
    Yn360 = build_Y360(phi, Q_nu, R_vec, kappa, g_map)



    # ------------------------------------------------------------
    # Step 6: Project 360→9
    # ------------------------------------------------------------
    Yu9 = project_to_9x9(Yu360, B)
    Yd9 = project_to_9x9(Yd360, B)
    Ye9 = project_to_9x9(Ye360, B)
    Yn9 = project_to_9x9(Yn360, B)

    # ------------------------------------------------------------
    # Step 7: Apply Φ to each 9×9
    # ------------------------------------------------------------
    print("→ Applying alignment operator Φ...")
    Yu9a = apply_alignment(Yu9, kappa, forbidden_d)
    Yd9a = apply_alignment(Yd9, kappa, forbidden_d)
    Ye9a = apply_alignment(Ye9, kappa, forbidden_d)
    Yn9a = apply_alignment(Yn9, kappa, forbidden_d)

    # ------------------------------------------------------------
    # Step 8: Schur reduction to 3×3 Yukawas
    # ------------------------------------------------------------
    print("→ Schur reducing to 3×3 Yukawas...")
    Yu3 = schur_reduce_dynamic(Yu9a)
    Yd3 = schur_reduce_dynamic(Yd9a)
    Ye3 = schur_reduce_dynamic(Ye9a)
    Yn3 = schur_reduce_dynamic(Yn9a)

    # ------------------------------------------------------------
    # Step 9: CKM & PMNS
    # ------------------------------------------------------------
    print("→ Computing CKM & PMNS...")
    # ------------------------------------------------------------
    # Build generational eigenvector matrix (360 × 3)
    # ------------------------------------------------------------
    gen_vecs = evecs[:, triad]  # triad eigenvectors define generation subspace

    # ------------------------------------------------------------
    # Build region masks for leptons and quarks
    # Leptons use pure phi-regions
    # Quarks use shifted versions to induce CKM misalignment
    # ------------------------------------------------------------
    regions_l = build_regions(phi)  # leptons

    phi_shifted = (phi + 0.0137) % (2 * np.pi)  # tiny geometric phase offset
    regions_u = build_regions(phi)  # up quarks use unshifted
    regions_d = build_regions(phi_shifted)  # down quarks use shifted → CKM ≠ I

    # ------------------------------------------------------------
    # Build geometric unitaries from generational vectors + regions
    # ------------------------------------------------------------
    U_geom_u = build_geometric_unitary(gen_vecs, regions_u)
    U_geom_d = build_geometric_unitary(gen_vecs, regions_d)
    gen_vecs = evecs[:, triad]  # this is the correct generational basis: (360×3)

    assign_e, assign_nu, mixchi2, mix_data = search_best_lepton_regions(
        gen_vecs,
        regions_l,  # leptons use unshifted φ-regions
        U_geom_u,
        U_geom_d,
        Yu3, Yd3, Ye3, Yn3,
        None, None, None
    )

    # Use the CKM/PMNS matrices that were actually used in the χ² search
    V_ckm = mix_data["V_ckm"]
    U_pmns = mix_data["U_pmns"]

    # ------------------------------------------------------------
    # Step 10: Full chi²
    # ------------------------------------------------------------
    # Compute mass ratios directly from Yukawas

    ru = mass_ratios(Yu3)
    rd = mass_ratios(Yd3)
    re = mass_ratios(Ye3)

    # Extract angles
    theta12_q, theta23_q, theta13_q = mix_data["angles_q"]
    theta12_l, theta23_l, theta13_l = mix_data["angles_l"]

    # Build the observable dictionary
    observables = {
        "mu/mt": ru[0],
        "mc/mt": ru[1],
        "md/mb": rd[0],
        "ms/mb": rd[1],
        "me/mtau": re[0],
        "mmu/mtau": re[1],
        "theta12_q": theta12_q,
        "theta23_q": theta23_q,
        "theta13_q": theta13_q,
        "theta12_l": theta12_l,
        "theta23_l": theta23_l,
        "theta13_l": theta13_l,
    }

    chi2_total, details = chi2(observables, TARGETS)

    print("\n=== RESULTS ===")
    print("CKM =\n", V_ckm)
    print("PMNS =\n", U_pmns)
    print("χ² =", chi2_total)
    print("Mixing contribution =", mixchi2)

    return dict(
        Y3=dict(u=Yu3, d=Yd3, e=Ye3, nu=Yn3),
        Vckm=V_ckm,
        Upmns=U_pmns,
        chi2=chi2_total,
        details=details
    )


if __name__ == "__main__":
    test_run()

"""
=== UTM v2: Unified Toy Model (Harmonic-Neighborhood) ===
N = 60
kappa = 8
forbidden_d = 2
seed = 42
→ Relaxing phases...
→ Building adjacency + Laplacian...
→ Extracting harmonic triad...
→ Building R and Q operators...
→ Building harmonic-neighborhood proto basis...
→ Building 360×360 Yukawa seeds...
→ Applying alignment operator Φ...
→ Schur reducing to 3×3 Yukawas...
→ Computing CKM & PMNS...

=== RESULTS ===
CKM =
 [[1.00000000e+00+0.j 1.98321845e-16+0.j 1.00201270e-16+0.j]
 [1.98321845e-16+0.j 1.00000000e+00+0.j 1.55563096e-16+0.j]
 [1.11022302e-16+0.j 1.80411242e-16+0.j 1.00000000e+00+0.j]]
PMNS =
 [[ 0.97522367+0.j  0.17364818+0.j  0.13705875+0.j]
 [-0.12644838+0.j  0.94592532+0.j -0.2987241 +0.j]
 [-0.18152024+0.j  0.27399196+0.j  0.9444463 +0.j]]
χ² = 93.437308012988
Mixing contribution = 93.437308012988
"""