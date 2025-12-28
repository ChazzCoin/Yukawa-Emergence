import numpy as np
import math

"""
Emergent aether toy:
====================

1) Start with N sites, each with a phase theta_i ∈ [0, 2π).
2) Define a misalignment functional

       M[theta] = sum_{i<j} J_ij [ (1 - cos(6Δ_ij)) + (1 - cos(5Δ_ij)) ]

   where Δ_ij = theta_i - theta_j.

   - The cos(6Δ) term encodes a 6-fold (60°) alignment preference.
   - The cos(5Δ) term encodes a 5-fold (72°) golden alignment preference.
   - Competing 5- and 6-fold preferences lead to frustrated, quasi-crystal-like order.

3) Perform gradient descent on {theta_i} to minimize M.

4) From the relaxed configuration, build an emergent adjacency matrix

       S_ij = cos(6Δ_ij) + cos(5Δ_ij)
       W_ij = max(0, S_ij)

   and keep only the strongest edges to define an unweighted graph A_int.

5) Build the Laplacian L_int from A_int.

6) Plug this L_int into the operator-first flavor machinery:

   - extract a 3-mode "generation triad" from its spectrum,
   - build F_base(λ), integer-charge hierarchies F_s,
   - apply golden P_phi and Cabibbo C_12 to get CKM & PMNS,
   - compute rough chi^2 vs SM-inspired targets.

This is still a toy, but now the internal graph is *emergent from operator-like rules*,
not chosen a priori as fib2d or 24-cell.
"""
import itertools

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


# ----------------------------------------------------------------------
# 2. Build emergent adjacency and Laplacian from relaxed phases
# ----------------------------------------------------------------------
def build_geometric_regions(theta: np.ndarray, n_regions: int = 3):
    """
    Partition sites into n_regions contiguous blocks in phase-order.
    This uses only the emergent phase field: no coordinates assumed.
    """
    theta = theta.reshape(-1)
    N = theta.shape[0]
    order = np.argsort(theta)  # sites sorted by phase
    # Split the sorted list into n_regions nearly equal chunks
    base = N // n_regions
    extra = N % n_regions
    regions = []
    start = 0
    for r in range(n_regions):
        size = base + (1 if r < extra else 0)
        idx = order[start:start+size]
        regions.append(idx)
        start += size
    return regions  # list of arrays of site indices

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

def laplacian_from_adjacency(A: np.ndarray) -> np.ndarray:
    D = np.diag(A.sum(axis=1))
    return D - A


# ----------------------------------------------------------------------
# 3. Operator-first flavor machinery (reused structure)
# ----------------------------------------------------------------------

def spectral_triad(L_int: np.ndarray):
    """
    Extract a 3-mode generation triad from L_int:
    the three lowest nonzero eigenvalues.
    """
    eigvals, eigvecs = np.linalg.eigh(L_int)
    eps = 1e-10
    nonzero_indices = np.where(eigvals > eps)[0]
    gen_indices = nonzero_indices[:3]
    lam_gen = eigvals[gen_indices]
    return lam_gen, gen_indices, eigvals

def base_kernel(lams: np.ndarray, alpha: float = 3.0, form: str = "lambda_sq") -> np.ndarray:
    """
    Spectral kernel F(lambda). We keep the same choices as before:
      "lambda_sq":  F = exp(-alpha * lambda^2)
    """
    if form == "lambda_sq":
        return np.exp(-alpha * (lams ** 2))
    elif form == "lambda":
        return np.exp(-alpha * lams)
    else:
        raise ValueError(f"Unknown kernel form: {form}")

def build_sector_charges():
    """
    Integer charges q_{s,g} for sector+generation hierarchies.

    Indices g = 0,1,2 correspond to the three internal modes in lam_gen
    (here ~[0.98, 1.82, 2.0]). Physical generations (1st,2nd,3rd) are
    determined by sorting the resulting F_s.

    These q_{s,g} are small integers chosen so that, given the fixed
    emergent F_base(lambda_gen), the sorted mass ratios (m1/m3, m2/m3)
    in each sector approximate the observed SM hierarchies:

      - Up:   mu/mt ~ 2.2e-5, mc/mt ~ 7.5e-3
      - Down: md/mb ~ 1.1e-3, ms/mb ~ 2.2e-2
      - E:    me/mtau ~ 2.9e-4, mmu/mtau ~ 5.9e-2

    No continuous Yukawa parameters are introduced; only discrete
    integer exponents acting on the emergent 3-mode triad.
    """
    sector_charges_gen = {
        # Up-type quarks
        "u":  np.array([4.0, 8.0, 0.0]),
        # Down-type quarks
        "d":  np.array([5.0, 5.0, 0.0]),
        # Charged leptons
        "e":  np.array([4.0, 0.0, 3.0]),
        # Neutrinos (kept as a simple, more-suppressed pattern for now)
        "nu": np.array([6.0, 5.0, 4.0]),
    }
    return sector_charges_gen

def sector_weights(F_base: np.ndarray, q_vec: np.ndarray, beta: float = 1.0) -> np.ndarray:
    return F_base * np.exp(-beta * q_vec)

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

def build_generation_operators(phi_order=5, cab_denom=28):
    theta_phi = 2.0 * math.pi / float(phi_order)
    theta_C   = 2.0 * math.pi / float(cab_denom)
    P_phi_12 = rot12(theta_phi)
    P_phi_23 = rot23(theta_phi)
    C_12     = rot12(theta_C)
    return P_phi_12, P_phi_23, C_12, theta_phi, theta_C

def build_sector_bases(P_phi_12, P_phi_23, C_12,
                       U_geom,
                       use_neutrino_dressing=True,
                       N_SOLAR=36, N_REACTOR=45):
    I3 = np.eye(3, dtype=complex)

    # Quarks
    U_L_u  = U_geom["u"]  @ P_phi_12
    U_L_d  = U_geom["d"]  @ P_phi_12 @ C_12

    # Charged leptons
    U_L_e  = U_geom["e"]

    # Neutrinos
    if use_neutrino_dressing:
        theta12_nu = 2.0 * math.pi / float(N_SOLAR)
        theta13_nu = 2.0 * math.pi / float(N_REACTOR)
        U_L_nu = U_geom["nu"] @ rot12(theta12_nu) @ P_phi_23 @ rot13(theta13_nu)
    else:
        U_L_nu = U_geom["nu"] @ P_phi_23

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

def yukawa_from_F_and_UL(F_s: np.ndarray, U_L: np.ndarray, U_R: np.ndarray) -> np.ndarray:
    F_diag = np.diag(F_s.astype(complex))
    return U_L.conj().T @ F_diag @ U_R

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

def chi2(observables, targets):
    chi2_total = 0.0
    details = []
    ratio_keys = ["mu/mt", "mc/mt", "md/mb", "ms/mb", "me/mtau", "mmu/mtau"]
    angle_keys = ["theta12_q", "theta23_q", "theta13_q",
                  "theta12_l", "theta23_l", "theta13_l"]

    for k in ratio_keys:
        m = observables[k]
        t = targets[k]
        if m <= 0 or t <= 0:
            continue
        logm = math.log10(m)
        logt = math.log10(t)
        sigma_log = 0.3
        contrib = ((logm - logt) / sigma_log)**2
        chi2_total += contrib
        details.append((k, m, t, contrib))

    for k in angle_keys:
        m = observables[k]
        t = targets[k]
        sigma = 0.2
        contrib = ((m - t) / sigma)**2
        chi2_total += contrib
        details.append((k, m, t, contrib))

    return chi2_total, details

def largest_connected_component(A):
    N = A.shape[0]
    visited = np.zeros(N, dtype=bool)
    components = []

    for i in range(N):
        if not visited[i]:
            stack = [i]
            comp = []
            visited[i] = True
            while stack:
                u = stack.pop()
                comp.append(u)
                neighbors = np.where(A[u] > 0)[0]
                for v in neighbors:
                    if not visited[v]:
                        visited[v] = True
                        stack.append(v)
            components.append(comp)

    # pick largest
    comp_sizes = [len(c) for c in components]
    largest_idx = np.argmax(comp_sizes)
    nodes = np.array(components[largest_idx], dtype=int)

    # induced subgraph
    A_sub = A[np.ix_(nodes, nodes)]
    return A_sub, nodes
# ----------------------------------------------------------------------
# 4. Main: emergent graph → L_int → flavor operators
# ----------------------------------------------------------------------

def main():
    # Step 1: relax phases under misalignment functional
    N = 200
    theta_final, energy_hist = relax_phases(
        N=N,
        n_steps=600,
        eta=0.01,
        w6=1.0,
        w5=1.0,
        random_seed=42
    )
    print("Relaxation complete.")
    print(f"Final misalignment energy: {energy_hist[-1]:.6f}")
    print()

    # Step 2: build emergent adjacency and Laplacian
    A_int_full = build_emergent_adjacency(
        theta_final,
        w6=1.0,
        w5=1.0,
        keep_fraction=0.05
    )
    A_int, nodes = largest_connected_component(A_int_full)
    L_int = laplacian_from_adjacency(A_int)

    # Spectrum and generation triad
    lam_gen, gen_indices, eigvals = spectral_triad(L_int)
    F_base = base_kernel(lam_gen, alpha=3.0, form="lambda_sq")

    print("=== Emergent internal graph ===")
    print(f"Number of sites: {A_int.shape[0]}")
    print("First 10 eigenvalues of L_int:")
    print(eigvals[:10])
    print()
    print("Generation eigenvalue indices:", gen_indices)
    print("Generation triad lam_gen:", lam_gen)
    print("Base kernel F_base(lam_gen):", F_base)
    print()

    # Generation eigenvectors restricted to the triad
    eigvals_full, eigvecs_full = np.linalg.eigh(L_int)
    gen_vecs = eigvecs_full[:, gen_indices]  # shape (N_sub, 3)

    # Build 3 geometric regions from the emergent phase field,
    # restricted to the nodes in the largest connected component.
    theta_sub = theta_final[nodes]
    regions = build_geometric_regions(theta_sub, n_regions=3)
    R0, R1, R2 = regions

    # Quarks: share the same geometric basis so CKM stays Cabibbo-like
    assign_u = [R0, R1, R2]
    assign_d = [R0, R1, R2]
    U_geom_u = build_geometric_unitary(gen_vecs, assign_u)
    U_geom_d = build_geometric_unitary(gen_vecs, assign_d)

    # Sector charges & F_s (fixed integer Q pattern)
    sector_charges_gen = build_sector_charges()
    F_u = sector_weights(F_base, sector_charges_gen["u"],  beta=1.0)
    F_d = sector_weights(F_base, sector_charges_gen["d"],  beta=1.0)
    F_e = sector_weights(F_base, sector_charges_gen["e"],  beta=1.0)
    F_n = sector_weights(F_base, sector_charges_gen["nu"], beta=1.0)

    print("=== Yukawa-like mass scales F_s ===")
    print("Up-type (F_u):        ", F_u)
    print("Down-type (F_d):      ", F_d)
    print("Charged leptons (F_e):", F_e)
    print("Neutrino (F_n):       ", F_n)
    print()

    # Generation-space operators (golden + Cabibbo)
    P_phi_12, P_phi_23, C_12, theta_phi, theta_C = build_generation_operators(
        phi_order=5,
        cab_denom=28
    )

    # Step 3: search over geometric assignments for e and nu
    best_pe, best_pn, best_chi2, best_dat = search_best_lepton_regions(
        gen_vecs,
        regions,
        U_geom_u, U_geom_d,
        F_u, F_d, F_e, F_n,
        P_phi_12, P_phi_23, C_12,
        N_SOLAR=36, N_REACTOR=45
    )

    print("Best lepton region permutations:")
    print("  pe (e sectors)  =", best_pe)
    print("  pn (nu sectors) =", best_pn)
    print(f"Best total chi^2  ≈ {best_chi2:.2f}")
    print()

    # Reconstruct the best U_geom using that assignment
    region_list = [R0, R1, R2]
    assign_e  = [region_list[i] for i in best_pe]
    assign_nu = [region_list[i] for i in best_pn]

    U_geom = {
        "u":  U_geom_u,
        "d":  U_geom_d,
        "e":  build_geometric_unitary(gen_vecs, assign_e),
        "nu": build_geometric_unitary(gen_vecs, assign_nu),
    }

    # Build sector bases using both geometry and flavor operators
    sector_bases = build_sector_bases(
        P_phi_12, P_phi_23, C_12,
        U_geom,
        use_neutrino_dressing=True,
        N_SOLAR=36,
        N_REACTOR=45
    )

    U_L_u,  U_R_u  = sector_bases["u"]
    U_L_d,  U_R_d  = sector_bases["d"]
    U_L_e,  U_R_e  = sector_bases["e"]
    U_L_nu, U_R_nu = sector_bases["nu"]

    # Yukawa-like operators (not strictly needed for ratios, but kept for completeness)
    Y_u  = yukawa_from_F_and_UL(F_u, U_L_u,  U_R_u)
    Y_d  = yukawa_from_F_and_UL(F_d, U_L_d,  U_R_d)
    Y_e  = yukawa_from_F_and_UL(F_e, U_L_e,  U_R_e)
    Y_nu = yukawa_from_F_and_UL(F_n, U_L_nu, U_R_nu)

    # Mass ratios from F_s (eigenvalues of Y†Y will be very close to these)
    mu_mt, mc_mt   = mass_ratios(F_u)
    md_mb, ms_mb   = mass_ratios(F_d)
    me_mt, mmu_mt  = mass_ratios(F_e)

    print("Mass ratios (m1/m3, m2/m3) from F_s:")
    print(f"mu/mt:     {mu_mt:.3e}, mc/mt:     {mc_mt:.3e}")
    print(f"md/mb:     {md_mb:.3e}, ms/mb:     {ms_mb:.3e}")
    print(f"me/mtau:   {me_mt:.3e}, mmu/mtau:  {mmu_mt:.3e}")
    print()

    # Step 5: mixing matrices
    V_ckm  = mixing_matrix(U_L_u, U_L_d)
    U_pmns = mixing_matrix(U_L_e, U_L_nu)

    theta12_q, theta23_q, theta13_q = mixing_angles_from_U(V_ckm)
    theta12_l, theta23_l, theta13_l = mixing_angles_from_U(U_pmns)

    print("=== CKM-like mixing matrix (geometry + operator) ===")
    print(V_ckm)
    print(f"theta12_q ≈ {theta12_q:.3f} rad, theta23_q ≈ {theta23_q:.3f}, theta13_q ≈ {theta13_q:.3e}")
    print(f"(Cabibbo operator angle = 2π/28 ≈ {theta_C:.3f} rad)")
    print()

    print("=== PMNS-like mixing matrix (geometry + operator) ===")
    print(U_pmns)
    print(f"theta12_l ≈ {theta12_l:.3f} rad, theta23_l ≈ {theta23_l:.3f}, theta13_l ≈ {theta13_l:.3e}")
    print(f"(Golden operator angle = 2π/5 ≈ {theta_phi:.3f} rad)")
    print()

    # Step 6: chi^2 vs rough targets (recompute for the final configuration)
    obs = compute_observables(
        mu_mt, mc_mt, md_mb, ms_mb, me_mt, mmu_mt,
        theta12_q, theta23_q, theta13_q,
        theta12_l, theta23_l, theta13_l
    )
    chi2_value, chi2_details = chi2(obs, TARGETS)

    print("=== Observables vs rough targets ===")
    for k, m, t, contrib in chi2_details:
        print(f"{k:12s}: model={m:.3e}, target={t:.3e}, chi2_contrib={contrib:.2f}")
    print()
    print(f"Total chi^2 ≈ {chi2_value:.2f}")
    print()
    print("NOTES:")
    print("- The internal graph is emergent from the misalignment functional M[theta],")
    print("  which encodes 6-fold (C_360) and 5-fold (golden) alignment preferences.")
    print("- We then restrict to the largest connected component to define a single,")
    print("  coherent aether vacuum, and build its Laplacian L_int.")
    print("- The generation triad and F_base(lambda) come from the spectrum of L_int,")
    print("  sector hierarchies from discrete integer charges Q_{s,g}, and mixing")
    print("  from a combination of geometry-derived U_geom[s] and fixed operators")
    print("  P_phi (golden) and C_12 (Cabibbo).")
    print("- No random Yukawas or continuous per-sector fits are used; everything")
    print("  comes from the emergent graph, a universal kernel, integer exponents,")
    print("  and discrete 2π/n phase rotations.")

if __name__ == "__main__":
    main()

"""
RESULTS:

Relaxation complete.
Final misalignment energy: 99.972531

=== Emergent internal graph ===
Number of sites: 39
First 10 eigenvalues of L_int:
[7.04495820e-17 9.80951287e-01 1.81564639e+00 2.00000000e+00
 2.00000000e+00 2.00000000e+00 2.00000000e+00 2.00000000e+00
 4.38302718e+00 5.00000000e+00]

Generation eigenvalue indices: [1 2 3]
Generation triad lam_gen: [0.98095129 1.81564639 2.        ]
Base kernel F_base(lam_gen): [5.57545487e-02 5.06933717e-05 6.14421235e-06]

=== Yukawa-like mass scales F_s ===
Up-type (F_u):         [1.02118018e-03 1.70057317e-08 6.14421235e-06]
Down-type (F_d):       [3.75671194e-04 3.41569251e-07 6.14421235e-06]
Charged leptons (F_e): [1.02118018e-03 5.06933717e-05 3.05902321e-07]
Neutrino (F_n):        [1.38201709e-04 3.41569251e-07 1.12535175e-07]

Best lepton region permutations:
  pe (e sectors)  = (1, 2, 0)
  pn (nu sectors) = (2, 0, 1)
Best total chi^2  ≈ 11.90

Mass ratios (m1/m3, m2/m3) from F_s:
mu/mt:     1.665e-05, mc/mt:     6.017e-03
md/mb:     9.092e-04, ms/mb:     1.636e-02
me/mtau:   2.996e-04, mmu/mtau:  4.964e-02

=== CKM-like mixing matrix (geometry + operator) ===
[[ 9.74927912e-01+0.j  2.22520934e-01+0.j -3.35683117e-17+0.j]
 [-2.22520934e-01+0.j  9.74927912e-01+0.j -5.38467983e-17+0.j]
 [-5.55111512e-17+0.j -5.55111512e-17+0.j  1.00000000e+00+0.j]]
theta12_q ≈ 0.224 rad, theta23_q ≈ 0.000, theta13_q ≈ 3.357e-17
(Cabibbo operator angle = 2π/28 ≈ 0.224 rad)

=== PMNS-like mixing matrix (geometry + operator) ===
[[-0.82852367+0.j  0.20511282+0.j  0.52103481+0.j]
 [-0.55830005+0.j -0.23112818+0.j -0.79679409+0.j]
 [-0.04300685+0.j -0.95105652+0.j  0.30600966+0.j]]
theta12_l ≈ 0.243 rad, theta23_l ≈ 1.204, theta13_l ≈ 5.481e-01
(Golden operator angle = 2π/5 ≈ 1.257 rad)

=== Observables vs rough targets ===
mu/mt       : model=1.665e-05, target=2.200e-05, chi2_contrib=0.16
mc/mt       : model=6.017e-03, target=7.500e-03, chi2_contrib=0.10
md/mb       : model=9.092e-04, target=1.100e-03, chi2_contrib=0.08
ms/mb       : model=1.636e-02, target=2.200e-02, chi2_contrib=0.18
me/mtau     : model=2.996e-04, target=2.900e-04, chi2_contrib=0.00
mmu/mtau    : model=4.964e-02, target=5.900e-02, chi2_contrib=0.06
theta12_q   : model=2.244e-01, target=2.270e-01, chi2_contrib=0.00
theta23_q   : model=5.385e-17, target=4.100e-02, chi2_contrib=0.04
theta13_q   : model=3.357e-17, target=3.600e-03, chi2_contrib=0.00
theta12_l   : model=2.427e-01, target=5.840e-01, chi2_contrib=2.91
theta23_l   : model=1.204e+00, target=7.850e-01, chi2_contrib=4.39
theta13_l   : model=5.481e-01, target=1.500e-01, chi2_contrib=3.96

Total chi^2 ≈ 11.90
"""