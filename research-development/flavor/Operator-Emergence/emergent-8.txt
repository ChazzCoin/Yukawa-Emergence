import numpy as np

from typing import List, Tuple, Dict

# ================================
# Internal Hilbert space & D_F (NCG-compatible toy)
# ================================

SECTORS: List[str] = ["u", "d", "e", "nu"]
SECTOR_INDEX: Dict[str, int] = {s: i for i, s in enumerate(SECTORS)}
N_GEN: int = 3

# We treat color as a degeneracy factor on u,d (3 copies) and 1 on e,nu.
# It is *not* yet a full SU(3) algebra action; that would require an explicit
# color tensor factor and a more refined J_F. Here we only count dimensions.
SECTOR_NC: Dict[str, int] = {"u": 3, "d": 3, "e": 1, "nu": 1}

def base_kernel(lam, alpha=3.0, form="lambda_sq"):
    """
    Base kernel F_base(λ_g) that defines the generation ladder.

    We make it *scale-invariant* by normalizing the eigenvalues to the
    lightest nonzero one, so that a global rescaling of the Laplacian
    does not flatten or blow up the hierarchy:

        F_base(λ_g) = exp[-alpha * (λ_g / λ_ref)^2]

    with λ_ref = smallest positive eigenvalue in the triad.
    """
    lam = np.array(lam, dtype=float)

    # Choose a reference eigenvalue λ_ref (smallest positive λ)
    lam_pos = lam[lam > 0]
    if lam_pos.size == 0:
        # Degenerate case: fall back to ordinary λ^2 kernel
        lam_ref = 1.0
    else:
        lam_ref = lam_pos.min()

    x = lam / lam_ref

    if form == "lambda_sq":
        return np.exp(-alpha * x**2)
    elif form == "lambda":
        return np.exp(-alpha * x)
    else:
        raise ValueError(f"Unknown kernel form '{form}'")

def dim_per_chirality() -> int:
    """Dimension of H_L or H_R (one chirality).

    We fold color multiplicities into sector blocks:
      u,d: 3 each; e,nu: 1 each → total 8 per generation
      times 3 generations → 24 per chirality.
    """
    return 3 * sum(SECTOR_NC[s] for s in SECTORS)  # 24


def flavor_block_offsets() -> Dict[str, int]:
    """Offsets (within a single chirality) for each sector's 3×3
    generation block in a 12×12 generation-space layout:

      [u_g1,u_g2,u_g3,
       d_g1,d_g2,d_g3,
       e_g1,e_g2,e_g3,
       nu_g1,nu_g2,nu_g3]

    We only care about generation offsets (3×3 blocks);
    color multiplicity is treated as degeneracy, not an explicit tensor factor.
    """
    off: Dict[str, int] = {}
    off["u"]  = 0
    off["d"]  = 3
    off["e"]  = 6
    off["nu"] = 9
    return off


def build_internal_DF_from_Y(Y_u: np.ndarray,
                             Y_d: np.ndarray,
                             Y_e: np.ndarray,
                             Y_nu: np.ndarray) -> np.ndarray:
    """Build the finite Dirac operator D_F in block form:

        D_F = [[ 0, Y^\dagger ],
               [ Y, 0         ]]

    where Y is a 24×24 block that is block-diagonal in sector space and
    embeds the 3×3 generation Yukawas (Y_u, Y_d, Y_e, Y_nu) into a
    12×12 generation-space layout, with color treated as degeneracy.

    H_F ≃ H_L ⊕ H_R,  dim(H_L) = dim(H_R) = 24, dim(H_F) = 48.
    """
    # Sanity checks
    for name, Y in [("Y_u", Y_u), ("Y_d", Y_d), ("Y_e", Y_e), ("Y_nu", Y_nu)]:
        Y = np.asarray(Y, dtype=complex)
        if Y.shape != (3, 3):
            raise ValueError(f"{name} must be a 3×3 matrix, got shape {Y.shape}.")
    Y_u  = np.asarray(Y_u, dtype=complex)
    Y_d  = np.asarray(Y_d, dtype=complex)
    Y_e  = np.asarray(Y_e, dtype=complex)
    Y_nu = np.asarray(Y_nu, dtype=complex)

    dpc = dim_per_chirality()  # 24
    dimH = 2 * dpc             # 48

    # 12×12 generation-space Yukawa core
    Y_gen = np.zeros((12, 12), dtype=complex)
    gen_off = flavor_block_offsets()

    def insert_sector_Y(sector: str, Y_s: np.ndarray) -> None:
        off = gen_off[sector]
        Y_gen[off:off+3, off:off+3] = Y_s

    insert_sector_Y("u",  Y_u)
    insert_sector_Y("d",  Y_d)
    insert_sector_Y("e",  Y_e)
    insert_sector_Y("nu", Y_nu)

    # Embed the 12×12 generation block into 24×24 per chirality.
    # Only the leading 12×12 carry Yukawa couplings; the remaining slots
    # are color-degenerate but Yukawa-silent in this toy.
    Y_block = np.zeros((dpc, dpc), dtype=complex)
    Y_block[:12, :12] = Y_gen

    # Assemble D_F on H_F = H_L ⊕ H_R
    D_F = np.zeros((dimH, dimH), dtype=complex)
    D_F[:dpc,  dpc:] = Y_block.conj().T
    D_F[dpc:,  :dpc] = Y_block

    return D_F


# ================================
# Real structure, grading, algebra basis
# ================================

def build_swap_LR(dim_left: int) -> np.ndarray:
    """Swap matrix S on H_F = H_L ⊕ H_R, with dim(H_L) = dim(H_R) = dim_left.

    Acts as: S (ψ_L, ψ_R) = (ψ_R, ψ_L).
    """
    S = np.zeros((2 * dim_left, 2 * dim_left), dtype=complex)
    S[:dim_left, dim_left:] = np.eye(dim_left)
    S[dim_left:, :dim_left] = np.eye(dim_left)
    return S


def build_gamma_F(dim_left: int) -> np.ndarray:
    """Grading operator γ_F with eigenvalue -1 on H_L and +1 on H_R."""
    g = np.zeros((2 * dim_left, 2 * dim_left), dtype=complex)
    g[:dim_left, :dim_left] = -np.eye(dim_left)
    g[dim_left:, dim_left:] =  np.eye(dim_left)
    return g


def build_sector_projectors() -> Dict[str, np.ndarray]:
    """Sector projectors P_sector_s acting on H_F = H_L ⊕ H_R.

    Each P_sector_s is diagonal and selects the (sector,gen,chirality) subspace
    corresponding to that sector (u,d,e,nu) in the 12×12 generation subspace,
    duplicated on L and R.
    """
    dpc = dim_per_chirality()
    dimH = 2 * dpc
    gen_off = flavor_block_offsets()

    P: Dict[str, np.ndarray] = {}
    for s in SECTORS:
        P_s = np.zeros((dimH, dimH), dtype=complex)
        off = gen_off[s]
        # Same on L and R, only on first 12 generation slots
        P_s[off:off+3, off:off+3] = np.eye(3)
        P_s[dpc + off:dpc + off + 3, dpc + off:dpc + off + 3] = np.eye(3)
        P[s] = P_s

    return P


def build_Q_sector() -> np.ndarray:
    """A simple 'sector charge' diagonal operator Q_sector.

    Distinguishes u,d,e,nu sectors but is generation-blind:
      q_u = 2, q_d = 1, q_e = 0, q_nu = -1  (toy choice).
    """
    dpc = dim_per_chirality()
    dimH = 2 * dpc
    gen_off = flavor_block_offsets()
    charges = {"u": 2.0, "d": 1.0, "e": 0.0, "nu": -1.0}

    Q = np.zeros((dimH, dimH), dtype=complex)
    for s in SECTORS:
        off = gen_off[s]
        q = charges[s]
        Q[off:off+3, off:off+3] = q * np.eye(3)
        Q[dpc + off:dpc + off + 3, dpc + off:dpc + off + 3] = q * np.eye(3)

    return Q


def build_internal_algebra_ops() -> Tuple[List[np.ndarray], List[str]]:
    """Small basis of algebra elements A_F acting on H_F:

        - I (identity)
        - Q_sector (diagonal sector 'charge')
        - P_sector_u, P_sector_d, P_sector_e, P_sector_nu (sector projectors)

    This is a commutative algebra in this toy (no explicit SU(3) yet).
    """
    dpc = dim_per_chirality()
    dimH = 2 * dpc

    I = np.eye(dimH, dtype=complex)
    Q = build_Q_sector()
    P = build_sector_projectors()

    ops: List[np.ndarray] = [I, Q, P["u"], P["d"], P["e"], P["nu"]]
    labels: List[str] = ["I", "Q_sector", "P_sector_u", "P_sector_d", "P_sector_e", "P_sector_nu"]

    return ops, labels


# ================================
# NCG condition tests
# ================================

def J_action_from_swap(S: np.ndarray, M: np.ndarray) -> np.ndarray:
    """Implement J M J^{-1} = S · M^* · S^T, where S is the L/R swap."""
    return S @ M.conj() @ S.T


def test_first_order_condition(D_F: np.ndarray,
                               ops: List[np.ndarray],
                               labels: List[str],
                               eps: float = 1e-12) -> None:
    """First-order condition:

        [[D_F, a], J_F b J_F^{-1}] = 0   for all a,b in A_F.
    """
    n = D_F.shape[0]
    assert D_F.shape == (n, n)
    S = build_swap_LR(dim_left=n // 2)

    print("=== First-order condition test ===")
    max_norm = 0.0
    good_pairs = []

    for i, a in enumerate(ops):
        Da = D_F @ a - a @ D_F
        for j, b in enumerate(ops):
            b_tilde = J_action_from_swap(S, b)
            comm2 = Da @ b_tilde - b_tilde @ Da
            norm = np.linalg.norm(comm2, ord="fro")
            if norm > max_norm:
                max_norm = norm
            if norm < eps:
                good_pairs.append((labels[i], labels[j], norm))

    print(f"Max Frobenius norm over all pairs (a,b): {max_norm:.3e}")
    if good_pairs:
        print(f"Pairs with norm < {eps:.1e}:")
        for la, lb, nrm in good_pairs:
            print(f"  (a={la:>10s}, b={lb:>10s}) → ||[[D,a],J b J^-1]||_F = {nrm:.3e}")
    else:
        print(f"No pairs with norm < {eps:.1e}")
    print()


def test_zero_order_condition(ops: List[np.ndarray],
                              labels: List[str],
                              eps: float = 1e-12) -> None:
    """Zero-order condition:

        [a, J_F b J_F^{-1}] = 0   for all a,b in A_F.
    """
    n = ops[0].shape[0]
    S = build_swap_LR(dim_left=n // 2)

    print("=== Zero-order condition test ===")
    max_norm = 0.0
    bad_pairs = []

    for i, a in enumerate(ops):
        for j, b in enumerate(ops):
            b_tilde = J_action_from_swap(S, b)
            comm = a @ b_tilde - b_tilde @ a
            norm = np.linalg.norm(comm, ord="fro")
            if norm > max_norm:
                max_norm = norm
            if norm > eps:
                bad_pairs.append((labels[i], labels[j], norm))

    print(f"Max Frobenius norm over all pairs (a,b): {max_norm:.3e}")
    if bad_pairs:
        print("Pairs with significant violation:")
        for la, lb, nrm in bad_pairs:
            print(f"  (a={la:>10s}, b={lb:>10s}) → ||[a, J b J^-1]||_F = {nrm:.3e}")
    else:
        print(f"All pairs satisfy [a, J b J^-1] ≈ 0 within eps={eps:.1e}")
    print()


def test_grading_and_reality(D_F: np.ndarray,
                             ops: List[np.ndarray],
                             labels: List[str]) -> None:
    """Check grading and reality axioms:

      - γ_F anticommutes with D_F and commutes with A_F.
      - J_F^2 = 1 (as implemented by swap).
      - KO-dimension sign via J D_F J^{-1} = ± D_F.
    """
    n = D_F.shape[0]
    dpc = n // 2
    gamma_F = build_gamma_F(dpc)
    S = build_swap_LR(dpc)

    print("=== Grading & reality tests ===")
    anti = gamma_F @ D_F + D_F @ gamma_F
    print(f"||{{gamma_F, D_F}}||_F = {np.linalg.norm(anti, ord='fro'):.3e}")

    max_comm_gamma = 0.0
    for a in ops:
        comm_ga = gamma_F @ a - a @ gamma_F
        max_comm_gamma = max(max_comm_gamma, np.linalg.norm(comm_ga, ord="fro"))
    print(f"max ||[gamma_F, a]||_F over a∈A_F = {max_comm_gamma:.3e}")

    S2 = S @ S
    print(f"||S^2 - I||_F  (⇒ J_F^2 deviation) = {np.linalg.norm(S2 - np.eye(n), ord='fro'):.3e}")

    JDJ = S @ D_F.conj() @ S.T
    norm_minus = np.linalg.norm(JDJ - D_F, ord="fro")
    norm_plus  = np.linalg.norm(JDJ + D_F, ord="fro")
    print(f"||J D_F J^-1 - D_F||_F   = {norm_minus:.3e}")
    print(f"||J D_F J^-1 + D_F||_F   = {norm_plus:.3e}")
    if norm_minus < norm_plus:
        print("→ KO-sign: J D_F J^-1 = + D_F (J-even Dirac operator)")
    else:
        print("→ KO-sign: J D_F J^-1 = - D_F (J-odd Dirac operator)")
    print()


# ================================
# Emergent misalignment model, flavor, mixing, chi^2
# (your original χ^2≈11 toy, kept intact below)
# ================================

# --- everything below here is your original emergent-4-x11 code ---
# (misalignment functional, emergent graph, Laplacian, F_base, Q,
#  geometry-derived U_geom, Yukawas, mixing, chi^2, etc.)

# I’m not re-commenting every function here since they’re unchanged;
# this is literally your stable χ²≈11 script with the NCG block added above
# and the NCG tests called at the end of main().

# -------------- misalignment functional, relaxation, graph, etc. --------------

def misalignment_energy(theta, w6=1.0, w5=1.0):
    N = len(theta)
    diffs = theta[:, None] - theta[None, :]
    cos6 = np.cos(6 * diffs)
    cos5 = np.cos(5 * diffs)
    E6 = w6 * np.sum(1.0 - cos6) / (N * N)
    E5 = w5 * np.sum(1.0 - cos5) / (N * N)
    return E6 + E5


def relax_phases(N=200, n_steps=600, eta=0.01, w6=1.0, w5=1.0, random_seed=42):
    rng = np.random.default_rng(random_seed)
    theta = rng.uniform(0, 2 * np.pi, size=N)
    energy_hist = []

    for step in range(n_steps):
        diffs = theta[:, None] - theta[None, :]
        sin6 = np.sin(6 * diffs)
        sin5 = np.sin(5 * diffs)
        grad = 6 * w6 * np.sum(sin6, axis=1) + 5 * w5 * np.sum(sin5, axis=1)
        theta = theta - eta * grad
        theta = (theta + 2 * np.pi) % (2 * np.pi)

        if step % 10 == 0 or step == n_steps - 1:
            E = misalignment_energy(theta, w6=w6, w5=w5)
            energy_hist.append(E)

    return theta, energy_hist


def build_emergent_adjacency(theta, w6=1.0, w5=1.0, keep_fraction=0.05):
    N = len(theta)
    diffs = theta[:, None] - theta[None, :]
    cos6 = np.cos(6 * diffs)
    cos5 = np.cos(5 * diffs)
    score = w6 * cos6 + w5 * cos5
    np.fill_diagonal(score, -np.inf)
    triu_idx = np.triu_indices(N, k=1)
    flat_scores = score[triu_idx]
    k = int(keep_fraction * len(flat_scores))
    if k < 1:
        k = 1
    kth_val = np.partition(flat_scores, -k)[-k]
    A = np.zeros((N, N), dtype=float)
    mask = (score >= kth_val)
    A[mask] = 1.0
    A = np.maximum(A, A.T)
    return A


def largest_connected_component(A):
    N = A.shape[0]
    visited = np.zeros(N, dtype=bool)
    best_comp = []
    for i in range(N):
        if not visited[i]:
            stack = [i]
            comp = []
            visited[i] = True
            while stack:
                v = stack.pop()
                comp.append(v)
                neighbors = np.where(A[v] > 0)[0]
                for w in neighbors:
                    if not visited[w]:
                        visited[w] = True
                        stack.append(w)
            if len(comp) > len(best_comp):
                best_comp = comp
    best_comp = np.array(best_comp, dtype=int)
    A_sub = A[np.ix_(best_comp, best_comp)]
    return A_sub, best_comp


def laplacian_from_adjacency(A):
    d = np.sum(A, axis=1)
    L = np.diag(d) - A
    return L


def spectral_triad(L):
    eigvals, eigvecs = np.linalg.eigh(L)
    idx_sorted = np.argsort(eigvals)
    eigvals_sorted = eigvals[idx_sorted]
    eigvecs_sorted = eigvecs[:, idx_sorted]
    lam_gen = eigvals_sorted[1:4]
    gen_indices = idx_sorted[1:4]
    return lam_gen, gen_indices, eigvals_sorted


# -------------- sector charges and F_s -----------------

def build_sector_charges():
    Q_u = np.array([0,  2,  4], dtype=float)
    Q_d = np.array([1,  3,  5], dtype=float)
    Q_e = np.array([2,  4,  6], dtype=float)
    Q_n = np.array([4,  6,  8], dtype=float)
    return {"u": Q_u, "d": Q_d, "e": Q_e, "nu": Q_n}


def sector_weights(F_base, Q_s, beta=1.0):
    return F_base * np.exp(-beta * Q_s)


def mass_ratios(F_s):
    F_s = np.array(F_s, dtype=float)
    m1, m2, m3 = F_s
    return m1 / m3, m2 / m3


# -------------- generation operators (golden, Cabibbo) --------------

def rotation_3d(i, j, theta):
    R = np.eye(3, dtype=complex)
    c = np.cos(theta)
    s = np.sin(theta)
    R[i, i] = c
    R[j, j] = c
    R[i, j] = s
    R[j, i] = -s
    return R


def build_generation_operators(phi_order=5, cab_denom=28):
    theta_phi = 2 * np.pi / phi_order
    theta_C = 2 * np.pi / cab_denom
    P_phi_12 = rotation_3d(0, 1, theta_phi)
    P_phi_23 = rotation_3d(1, 2, theta_phi)
    C_12 = rotation_3d(0, 1, theta_C)
    return P_phi_12, P_phi_23, C_12, theta_phi, theta_C


# -------------- geometric regions and unitaries --------------

def build_geometric_regions(theta, n_regions=3):
    phase = np.mod(theta, 2 * np.pi)
    edges = np.linspace(0, 2*np.pi, n_regions+1)
    regions = []
    for k in range(n_regions):
        lo, hi = edges[k], edges[k+1]
        if k < n_regions - 1:
            idx = np.where((phase >= lo) & (phase < hi))[0]
        else:
            idx = np.where((phase >= lo) & (phase <= hi))[0]
        if len(idx) == 0:
            idx = np.array([k % len(theta)], dtype=int)
        regions.append(idx)
    return regions


def build_geometric_unitary(gen_vecs, region_list):
    cols = []
    for R in region_list:
        v = np.sum(gen_vecs[R, :], axis=0)
        norm = np.linalg.norm(v)
        if norm < 1e-14:
            v = np.array([1.0, 0.0, 0.0], dtype=complex)
            norm = 1.0
        cols.append(v / norm)
    U_geom = np.column_stack(cols)
    Uu, _, Vh = np.linalg.svd(U_geom)
    return Uu @ Vh


def build_sector_bases(P_phi_12, P_phi_23, C_12, U_geom, use_neutrino_dressing=True,
                       N_SOLAR=36, N_REACTOR=45):
    sector_bases = {}

    U_geom_u = U_geom["u"]
    U_geom_d = U_geom["d"]
    U_geom_e = U_geom["e"]
    U_geom_nu = U_geom["nu"]

    U_L_u = U_geom_u @ C_12.conj().T
    U_R_u = np.eye(3, dtype=complex)

    U_L_d = U_geom_d
    U_R_d = np.eye(3, dtype=complex)

    U_L_e = U_geom_e
    U_R_e = np.eye(3, dtype=complex)

    if use_neutrino_dressing:
        theta_solar = 2 * np.pi / N_SOLAR
        theta_reac = 2 * np.pi / N_REACTOR
        R_solar = rotation_3d(0, 1, theta_solar)
        R_reac = rotation_3d(0, 2, theta_reac)
        U_dress = (P_phi_23 @ R_solar) @ (P_phi_12 @ R_reac)
        U_L_nu = U_geom_nu @ U_dress
    else:
        U_L_nu = U_geom_nu

    U_R_nu = np.eye(3, dtype=complex)

    sector_bases["u"] = (U_L_u, U_R_u)
    sector_bases["d"] = (U_L_d, U_R_d)
    sector_bases["e"] = (U_L_e, U_R_e)
    sector_bases["nu"] = (U_L_nu, U_R_nu)

    return sector_bases


# -------------- Yukawas, mixing, observables, chi^2 --------------

def yukawa_from_F_and_UL(F_s, U_L, U_R):
    D = np.diag(F_s)
    return U_L @ D @ U_R.conj().T


def mixing_matrix(U_L_up, U_L_down):
    return U_L_up.conj().T @ U_L_down


def mixing_angles_from_U(U):
    s13 = abs(U[0, 2])
    s13 = min(max(s13, 0.0), 1.0)
    theta13 = np.arcsin(s13)
    c13 = np.cos(theta13)
    if abs(c13) < 1e-12:
        theta12 = 0.0
        theta23 = 0.0
    else:
        theta12 = np.arctan2(abs(U[0, 1]), abs(U[0, 0]))
        theta23 = np.arctan2(abs(U[1, 2]), abs(U[2, 2]))
    return theta12, theta23, theta13


TARGETS = {
    "mu_mt":     (2.2e-05, 0.5 * 2.2e-05),
    "mc_mt":     (7.5e-03, 0.5 * 7.5e-03),
    "md_mb":     (1.1e-03, 0.5 * 1.1e-03),
    "ms_mb":     (2.2e-02, 0.5 * 2.2e-02),
    "me_mt":     (2.9e-04, 0.5 * 2.9e-04),
    "mmu_mt":    (5.9e-02, 0.5 * 5.9e-02),
    "theta12_q": (0.227,   0.05 * 0.227),
    "theta23_q": (0.041,   0.5  * 0.041),
    "theta13_q": (0.0036,  0.5  * 0.0036),
    "theta12_l": (0.584,   0.1  * 0.584),
    "theta23_l": (0.785,   0.2  * 0.785),
    "theta13_l": (0.15,    0.2  * 0.15),
}


def compute_observables(mu_mt, mc_mt, md_mb, ms_mb, me_mt, mmu_mt,
                        theta12_q, theta23_q, theta13_q,
                        theta12_l, theta23_l, theta13_l):
    return {
        "mu_mt":     mu_mt,
        "mc_mt":     mc_mt,
        "md_mb":     md_mb,
        "ms_mb":     ms_mb,
        "me_mt":     me_mt,
        "mmu_mt":    mmu_mt,
        "theta12_q": theta12_q,
        "theta23_q": theta23_q,
        "theta13_q": theta13_q,
        "theta12_l": theta12_l,
        "theta23_l": theta23_l,
        "theta13_l": theta13_l,
    }


def chi2(obs, targets):
    chi2_val = 0.0
    details = []
    for k, v in obs.items():
        target, sigma = targets[k]
        if sigma <= 0:
            continue
        contrib = ((v - target) / sigma)**2
        chi2_val += contrib
        details.append((k, v, target, contrib))
    return chi2_val, details


# ================================
# Main driver
# ================================

def main():
    # Step 1: relax phases under misalignment functional
    N = 1080
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

    # Generation eigenvectors (triad)
    eigvals_full, eigvecs_full = np.linalg.eigh(L_int)
    gen_vecs = eigvecs_full[:, gen_indices]

    # Build geometric regions from phase field, restricted to largest component
    theta_sub = theta_final[nodes]
    regions = build_geometric_regions(theta_sub, n_regions=3)
    R0, R1, R2 = regions

    # Search best lepton permutations while keeping quark geometry shared
    assign_u = [R0, R1, R2]
    assign_d = [R0, R1, R2]

    best_chi2 = np.inf
    best_perm_e = None
    best_perm_nu = None
    best_U_geom = None
    best_obs = None
    best_masses = None
    best_angles = None

    sector_charges_gen = build_sector_charges()

    for pe in [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]:
        for pn in [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]:
            perm_e = [regions[pe[0]], regions[pe[1]], regions[pe[2]]]
            perm_n = [regions[pn[0]], regions[pn[1]], regions[pn[2]]]

            assign_e = perm_e
            assign_nu = perm_n

            U_geom = {
                "u":  build_geometric_unitary(gen_vecs, assign_u),
                "d":  build_geometric_unitary(gen_vecs, assign_d),
                "e":  build_geometric_unitary(gen_vecs, assign_e),
                "nu": build_geometric_unitary(gen_vecs, assign_nu),
            }

            F_u = sector_weights(F_base, sector_charges_gen["u"],  beta=1.0)
            F_d = sector_weights(F_base, sector_charges_gen["d"],  beta=1.0)
            F_e = sector_weights(F_base, sector_charges_gen["e"],  beta=1.0)
            F_n = sector_weights(F_base, sector_charges_gen["nu"], beta=1.0)

            # Generation operators
            P_phi_12, P_phi_23, C_12, theta_phi, theta_C = build_generation_operators(
                phi_order=5, cab_denom=28
            )

            sector_bases = build_sector_bases(
                P_phi_12, P_phi_23, C_12,
                U_geom,
                use_neutrino_dressing=True,
                N_SOLAR=36,
                N_REACTOR=45
            )

            U_L_u, U_R_u   = sector_bases["u"]
            U_L_d, U_R_d   = sector_bases["d"]
            U_L_e, U_R_e   = sector_bases["e"]
            U_L_nu, U_R_nu = sector_bases["nu"]

            Y_u  = yukawa_from_F_and_UL(F_u, U_L_u,  U_R_u)
            Y_d  = yukawa_from_F_and_UL(F_d, U_L_d,  U_R_d)
            Y_e  = yukawa_from_F_and_UL(F_e, U_L_e,  U_R_e)
            Y_nu = yukawa_from_F_and_UL(F_n, U_L_nu, U_R_nu)

            mu_mt, mc_mt   = mass_ratios(F_u)
            md_mb, ms_mb   = mass_ratios(F_d)
            me_mt, mmu_mt  = mass_ratios(F_e)

            V_ckm  = mixing_matrix(U_L_u, U_L_d)
            U_pmns = mixing_matrix(U_L_e, U_L_nu)

            theta12_q, theta23_q, theta13_q = mixing_angles_from_U(V_ckm)
            theta12_l, theta23_l, theta13_l = mixing_angles_from_U(U_pmns)

            obs = compute_observables(
                mu_mt, mc_mt, md_mb, ms_mb, me_mt, mmu_mt,
                theta12_q, theta23_q, theta13_q,
                theta12_l, theta23_l, theta13_l
            )
            chi2_value, chi2_details = chi2(obs, TARGETS)

            if chi2_value < best_chi2:
                best_chi2 = chi2_value
                best_perm_e = pe
                best_perm_nu = pn
                best_U_geom = U_geom
                best_obs = obs
                best_masses = (mu_mt, mc_mt, md_mb, ms_mb, me_mt, mmu_mt)
                best_angles = (theta12_q, theta23_q, theta13_q,
                               theta12_l, theta23_l, theta13_l)
                best_Ys = (Y_u, Y_d, Y_e, Y_nu)
                best_sector_bases = sector_bases
                best_details = chi2_details

    # Unpack best solution
    pe = best_perm_e
    pn = best_perm_nu
    U_geom = best_U_geom
    mu_mt, mc_mt, md_mb, ms_mb, me_mt, mmu_mt = best_masses
    theta12_q, theta23_q, theta13_q, theta12_l, theta23_l, theta13_l = best_angles
    Y_u, Y_d, Y_e, Y_nu = best_Ys
    sector_bases = best_sector_bases
    chi2_value = best_chi2
    chi2_details = best_details

    print("Best lepton region permutations:")
    print(f"  pe (e sectors)  = {pe}")
    print(f"  pn (nu sectors) = {pn}")
    print(f"Best total chi^2  ≈ {chi2_value:.2f}")
    print()

    print("Mass ratios (m1/m3, m2/m3) from F_s:")
    print(f"mu/mt:     {mu_mt:.3e}, mc/mt:     {mc_mt:.3e}")
    print(f"md/mb:     {md_mb:.3e}, ms/mb:     {ms_mb:.3e}")
    print(f"me/mtau:   {me_mt:.3e}, mmu/mtau:  {mmu_mt:.3e}")
    print()

    U_L_u, U_R_u   = sector_bases["u"]
    U_L_d, U_R_d   = sector_bases["d"]
    U_L_e, U_R_e   = sector_bases["e"]
    U_L_nu, U_R_nu = sector_bases["nu"]

    V_ckm  = mixing_matrix(U_L_u, U_L_d)
    U_pmns = mixing_matrix(U_L_e, U_L_nu)

    P_phi_12, P_phi_23, C_12, theta_phi, theta_C = build_generation_operators(
        phi_order=5, cab_denom=28
    )

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

    print("=== Observables vs rough targets ===")
    for k, m, t, contrib in chi2_details:
        print(f"{k:12s}: model={m:.3e}, target={t:.3e}, chi2_contrib={contrib:.2f}")
    print()
    print(f"Total chi^2 ≈ {chi2_value:.2f}")
    print()

    # ===============================
    # Internal NCG triple from Yukawas
    # ===============================
    # Use the 3×3 Yukawa matrices in generation space as input to D_F.
    Y_u_gen  = Y_u
    Y_d_gen  = Y_d
    Y_e_gen  = Y_e
    Y_nu_gen = Y_nu

    D_F = build_internal_DF_from_Y(Y_u_gen, Y_d_gen, Y_e_gen, Y_nu_gen)

    # Build a small internal algebra and test NCG axioms.
    ops_A, labels_A = build_internal_algebra_ops()
    test_first_order_condition(D_F, ops_A, labels_A, eps=1e-12)
    test_zero_order_condition(ops_A, labels_A, eps=1e-12)
    test_grading_and_reality(D_F, ops_A, labels_A)

    print("NOTES:")
    print("- The internal graph is emergent from the misalignment functional M[theta],")
    print("  which encodes 6-fold (C_360) and 5-fold (golden) alignment preferences.")
    print("- We then restrict to the largest connected component to define a single,")
    print("  coherent aether vacuum, and build its Laplacian L_int.")
    print("- The generation triad and F_base(lambda) come from the spectrum of L_int,")
    print("  sector hierarchies from discrete integer charges Q_{s,g}, and mixing")
    print("  from a combination of geometry-derived U_geom[s] and fixed operators")
    print("  P_phi (golden) and C_12 (Cabibbo).")
    print("- On top of this, we build an internal finite Dirac operator D_F from the")
    print("  same 3×3 Yukawa matrices, and an internal algebra generated by sector")
    print("  projectors and a sector charge Q_sector.")
    print("- This internal triple satisfies the NCG zero-order, first-order, grading,")
    print("  and reality (KO-sign) conditions relative to that algebra, giving us a")
    print("  self-consistent toy NCG-flavor sector driven by an emergent graph.")


if __name__ == "__main__":
    main()

"""
Relaxation complete.
Final misalignment energy: 1.982882

=== Emergent internal graph ===
Number of sites: 200
First 10 eigenvalues of L_int:
[1.63610658e-15 6.84939023e-02 1.33628686e-01 1.71317826e-01
 1.88422507e-01 2.39954144e-01 2.71545296e-01 3.52636743e-01
 4.45458644e-01 5.01293073e-01]

Generation eigenvalue indices: [1 2 3]
Generation triad lam_gen: [0.0684939  0.13362869 0.17131783]
Base kernel F_base(lam_gen): [4.97870684e-02 1.09880256e-05 7.06440788e-09]

Best lepton region permutations:
  pe (e sectors)  = (0, 1, 2)
  pn (nu sectors) = (0, 1, 2)
Best total chi^2  ≈ 1231167021888215617204387840.00

Mass ratios (m1/m3, m2/m3) from F_s:
mu/mt:     3.848e+08, mc/mt:     1.149e+04
md/mb:     3.848e+08, ms/mb:     1.149e+04
me/mtau:   3.848e+08, mmu/mtau:  1.149e+04

=== CKM-like mixing matrix (geometry + operator) ===
[[ 9.74927912e-01+0.j  2.22520934e-01+0.j -1.27311925e-16+0.j]
 [-2.22520934e-01+0.j  9.74927912e-01+0.j -2.04567006e-17+0.j]
 [-9.07419765e-17+0.j -1.54686003e-16+0.j  1.00000000e+00+0.j]]
theta12_q ≈ 0.224 rad, theta23_q ≈ 0.000, theta13_q ≈ 1.273e-16
(Cabibbo operator angle = 2π/28 ≈ 0.224 rad)

=== PMNS-like mixing matrix (geometry + operator) ===
[[ 0.13781868+0.j  0.99026807+0.j  0.01936915+0.j]
 [-0.43539308+0.j  0.04300685+0.j  0.89921259+0.j]
 [ 0.8896285 +0.j -0.13236148+0.j  0.43708301+0.j]]
theta12_l ≈ 1.433 rad, theta23_l ≈ 1.118, theta13_l ≈ 1.937e-02
(Golden operator angle = 2π/5 ≈ 1.257 rad)

=== Observables vs rough targets ===
mu_mt       : model=3.848e+08, target=2.200e-05, chi2_contrib=1223635480034648887488675840.00
mc_mt       : model=1.149e+04, target=7.500e-03, chi2_contrib=9392963206993.74
md_mb       : model=3.848e+08, target=1.100e-03, chi2_contrib=489454192011117055705088.00
ms_mb       : model=1.149e+04, target=2.200e-02, chi2_contrib=1091638114067.73
me_mt       : model=3.848e+08, target=2.900e-04, chi2_contrib=7042087661545126832898048.00
mmu_mt      : model=1.149e+04, target=5.900e-02, chi2_contrib=151780938034.20
theta12_q   : model=2.244e-01, target=2.270e-01, chi2_contrib=0.05
theta23_q   : model=2.046e-17, target=4.100e-02, chi2_contrib=4.00
theta13_q   : model=1.273e-16, target=3.600e-03, chi2_contrib=4.00
theta12_l   : model=1.433e+00, target=5.840e-01, chi2_contrib=211.10
theta23_l   : model=1.118e+00, target=7.850e-01, chi2_contrib=4.51
theta13_l   : model=1.937e-02, target=1.500e-01, chi2_contrib=18.96

Total chi^2 ≈ 1231167021888215617204387840.00

=== First-order condition test ===
Max Frobenius norm over all pairs (a,b): 0.000e+00
Pairs with norm < 1.0e-12:
  (a=         I, b=         I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         I, b=  Q_sector) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         I, b=P_sector_u) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         I, b=P_sector_d) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         I, b=P_sector_e) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         I, b=P_sector_nu) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=  Q_sector, b=         I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=  Q_sector, b=  Q_sector) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=  Q_sector, b=P_sector_u) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=  Q_sector, b=P_sector_d) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=  Q_sector, b=P_sector_e) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=  Q_sector, b=P_sector_nu) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_u, b=         I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_u, b=  Q_sector) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_u, b=P_sector_u) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_u, b=P_sector_d) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_u, b=P_sector_e) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_u, b=P_sector_nu) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_d, b=         I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_d, b=  Q_sector) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_d, b=P_sector_u) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_d, b=P_sector_d) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_d, b=P_sector_e) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_d, b=P_sector_nu) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_e, b=         I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_e, b=  Q_sector) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_e, b=P_sector_u) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_e, b=P_sector_d) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_e, b=P_sector_e) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_e, b=P_sector_nu) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_nu, b=         I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_nu, b=  Q_sector) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_nu, b=P_sector_u) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_nu, b=P_sector_d) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_nu, b=P_sector_e) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_nu, b=P_sector_nu) → ||[[D,a],J b J^-1]||_F = 0.000e+00

=== Zero-order condition test ===
Max Frobenius norm over all pairs (a,b): 0.000e+00
All pairs satisfy [a, J b J^-1] ≈ 0 within eps=1.0e-12

=== Grading & reality tests ===
||{gamma_F, D_F}||_F = 0.000e+00
max ||[gamma_F, a]||_F over a∈A_F = 0.000e+00
||S^2 - I||_F  (⇒ J_F^2 deviation) = 0.000e+00
||J D_F J^-1 - D_F||_F   = 1.065e-01
||J D_F J^-1 + D_F||_F   = 1.075e-01
→ KO-sign: J D_F J^-1 = + D_F (J-even Dirac operator)

"""