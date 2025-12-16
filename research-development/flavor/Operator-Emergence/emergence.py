import numpy as np
from typing import List, Tuple, Dict

# =========================
# Internal Laplacian scaling
# =========================

from itertools import combinations


class InternalSpectrumConfig:
    """
    Config for rescaling the internal Laplacian eigenvalues
    before feeding them into the universal spectral kernel F_base.
    """
    # Global, sector-independent rescale factor for \hat L:
    # lam_scaled = L_rescale_factor * lam_raw.
    #
    # Start with something like 0.3–0.5 and tune later.
    L_rescale_factor: float = 0.3

    # How many low-lying nonzero eigenvalues we allow triads to pick from
    max_triad_index: int = 20  # use lam[1..20] as the "triad pool"

def rescale_laplacian_evals(lam_raw: np.ndarray,
                            cfg: InternalSpectrumConfig) -> np.ndarray:
    """
    Apply a global, sector-independent rescale factor to the internal Laplacian
    eigenvalues. This is equivalent to choosing the 'unit of internal distance'.

    lam_raw : 1D array of eigenvalues of L_int (sorted ascending).
    returns : lam_scaled (same shape).
    """
    return cfg.L_rescale_factor * lam_raw

def choose_quark_and_lepton_triads(lam: np.ndarray,
                                   max_triad_index: int = 20):
    """
    Given Laplacian eigenvalues lam (sorted ascending), choose two triads
    of indices (i<j<k) among the first max_triad_index nonzero eigenvalues:

      - quark triad  : 'nicely spaced' low-lying modes (small mixing)
      - lepton triad : two nearly-degenerate low modes + one higher (large mixing)

    Returns:
      triad_quark  : np.array([i_q1, i_q2, i_q3])
      triad_lepton : np.array([i_l1, i_l2, i_l3])
    """
    # ignore exact zero mode at index 0
    start = 1
    stop = min(max_triad_index, len(lam))
    nonzero_indices = np.arange(start, stop, dtype=int)

    if len(nonzero_indices) < 3:
        raise ValueError("Not enough nonzero eigenvalues to form triads.")

    triads = list(combinations(nonzero_indices, 3))

    def quark_score(triad):
        i, j, k = triad
        span = lam[k] - lam[i]
        middle_gap = lam[j] - lam[i]
        upper_gap = lam[k] - lam[j]
        skew = abs(middle_gap - upper_gap)
        # prefer large span, relatively even spacing
        return span - 0.5 * skew

    def lepton_score(triad):
        i, j, k = triad
        deg_pair_gap = abs(lam[j] - lam[i])
        upper_gap = lam[k] - lam[j]
        # prefer a nearly-degenerate pair (i,j) and a separated k
        return -deg_pair_gap + 0.5 * upper_gap

    best_q, best_q_score = None, -np.inf
    best_l, best_l_score = None, -np.inf

    for triad in triads:
        s_q = quark_score(triad)
        if s_q > best_q_score:
            best_q_score = s_q
            best_q = triad

        s_l = lepton_score(triad)
        if s_l > best_l_score:
            best_l_score = s_l
            best_l = triad

    triad_quark = np.array(best_q, dtype=int)
    triad_lepton = np.array(best_l, dtype=int)
    return triad_quark, triad_lepton

class FlavorNCGOperators:
    """
    Master class collecting:
      - Emergent misalignment + internal graph machinery
      - Flavor hierarchy and mixing operators
      - Internal NCG finite Dirac operator + algebra + tests
    """

    # ---------------------------
    # Global conventions
    # ---------------------------
    SECTORS: List[str] = ["u", "d", "e", "nu"]
    N_GEN: int = 3
    # Color multiplicities (degeneracies) per sector (not yet full SU(3) action)
    SECTOR_NC: Dict[str, int] = {"u": 3, "d": 3, "e": 1, "nu": 1}

    # Rough SM targets
    TARGETS: Dict[str, Tuple[float, float]] = {
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

    # ===========================
    # 1. Internal Hilbert space & D_F (NCG side)
    # ===========================

    def dim_per_chirality(self) -> int:
        """
        For now we only keep the 4×3 generation degrees of freedom explicitly.
        Color is *not* yet an explicit tensor factor here.
        """
        return len(self.SECTORS) * self.N_GEN  # 4 * 3 = 12

    def flavor_block_offsets(self) -> Dict[str, int]:
        """
        Offsets (within a single chirality) for each sector's 3×3
        generation block in a 12×12 generation-space layout:

          [u_g1,u_g2,u_g3,
           d_g1,d_g2,d_g3,
           e_g1,e_g2,e_g3,
           nu_g1,nu_g2,nu_g3]
        """
        off: Dict[str, int] = {}
        off["u"]  = 0
        off["d"]  = 3
        off["e"]  = 6
        off["nu"] = 9
        return off

    def build_internal_DF_from_Y(self, Y_u, Y_d, Y_e, Y_nu):
        """
        Build the finite Dirac operator D_F in block form:

            D_F = [[ 0, Y^\dagger ],
                   [ Y, 0         ]]

        where Y is a 12×12 block that is block-diagonal in sector space and
        embeds the 3×3 generation Yukawas (Y_u, Y_d, Y_e, Y_nu) into a
        4×3 generation-space layout:

            [u_g1,u_g2,u_g3,
             d_g1,d_g2,d_g3,
             e_g1,e_g2,e_g3,
             nu_g1,nu_g2,nu_g3]

        H_F ≃ H_L ⊕ H_R,  dim(H_L) = dim(H_R) = 12, dim(H_F) = 24.
        """
        # --- sanity checks on shapes ---
        for name, Y in [("Y_u", Y_u), ("Y_d", Y_d), ("Y_e", Y_e), ("Y_nu", Y_nu)]:
            arr = np.asarray(Y, dtype=complex)
            if arr.shape != (3, 3):
                raise ValueError(f"{name} must be a 3×3 matrix, got shape {arr.shape}.")

        # Promote to complex NumPy arrays
        Y_u = np.asarray(Y_u, dtype=complex)
        Y_d = np.asarray(Y_d, dtype=complex)
        Y_e = np.asarray(Y_e, dtype=complex)
        Y_nu = np.asarray(Y_nu, dtype=complex)

        # One chirality dimension: 4 sectors × 3 generations = 12
        dpc = self.dim_per_chirality()  # expected 12
        dimH = 2 * dpc  # 24

        # 12×12 generation-space Yukawa core
        Y_gen = np.zeros((dpc, dpc), dtype=complex)
        gen_off = self.flavor_block_offsets()

        def insert_sector_Y(sector: str, Y_s: np.ndarray) -> None:
            off = gen_off[sector]
            Y_gen[off:off + 3, off:off + 3] = Y_s

        insert_sector_Y("u", Y_u)
        insert_sector_Y("d", Y_d)
        insert_sector_Y("e", Y_e)
        insert_sector_Y("nu", Y_nu)

        # Here Y_block is just the 12×12 generation Yukawa
        Y_block = Y_gen

        # Assemble D_F on H_F = H_L ⊕ H_R
        D_F = np.zeros((dimH, dimH), dtype=complex)
        D_F[:dpc, dpc:] = Y_block.conj().T
        D_F[dpc:, :dpc] = Y_block

        return D_F

    # --- Real structure & grading ---

    def build_swap_LR(self, dim_left: int) -> np.ndarray:
        """Swap matrix S on H_F = H_L ⊕ H_R, dim(H_L) = dim(H_R) = dim_left."""
        S = np.zeros((2 * dim_left, 2 * dim_left), dtype=complex)
        S[:dim_left, dim_left:] = np.eye(dim_left)
        S[dim_left:, :dim_left] = np.eye(dim_left)
        return S

    def build_gamma_F(self, dim_left: int) -> np.ndarray:
        """Grading operator γ_F with eigenvalue -1 on H_L and +1 on H_R."""
        g = np.zeros((2 * dim_left, 2 * dim_left), dtype=complex)
        g[:dim_left, :dim_left] = -np.eye(dim_left)
        g[dim_left:, dim_left:] =  np.eye(dim_left)
        return g

    def build_sector_projectors(self):
        dpc = self.dim_per_chirality()  # 12
        dimH = 2 * dpc  # 24
        gen_off = self.flavor_block_offsets()

        P: Dict[str, np.ndarray] = {}
        for s in self.SECTORS:
            P_s = np.zeros((dimH, dimH), dtype=complex)
            off = gen_off[s]
            # Same on L and R, only on first 12 gen slots
            P_s[off:off+3, off:off+3] = np.eye(3)
            P_s[dpc + off:dpc + off + 3, dpc + off:dpc + off + 3] = np.eye(3)
            P[s] = P_s

        return P

    def build_Q_sector(self) -> np.ndarray:
        """
        A simple 'sector charge' diagonal operator Q_sector.

        Distinguishes u,d,e,nu sectors but is generation-blind:
          q_u = 2, q_d = 1, q_e = 0, q_nu = -1 (toy choice).
        """
        dpc = self.dim_per_chirality()
        dimH = 2 * dpc
        gen_off = self.flavor_block_offsets()
        charges = {"u": 2.0, "d": 1.0, "e": 0.0, "nu": -1.0}

        Q = np.zeros((dimH, dimH), dtype=complex)
        for s in self.SECTORS:
            off = gen_off[s]
            q = charges[s]
            Q[off:off+3, off:off+3] = q * np.eye(3)
            Q[dpc + off:dpc + off + 3, dpc + off:dpc + off + 3] = q * np.eye(3)

        return Q

    def build_internal_algebra_ops(self) -> Tuple[List[np.ndarray], List[str]]:
        """
        Small basis of algebra elements A_F acting on H_F:

          - I (identity)
          - Q_sector (diagonal sector 'charge')
          - P_sector_u, P_sector_d, P_sector_e, P_sector_nu (sector projectors)

        This is a commutative algebra in this toy (no explicit SU(3) yet).
        """
        dpc = self.dim_per_chirality()
        dimH = 2 * dpc

        I = np.eye(dimH, dtype=complex)
        Q = self.build_Q_sector()
        P = self.build_sector_projectors()

        ops: List[np.ndarray] = [I, Q, P["u"], P["d"], P["e"], P["nu"]]
        labels: List[str] = ["I", "Q_sector", "P_sector_u", "P_sector_d", "P_sector_e", "P_sector_nu"]

        return ops, labels

    # --- NCG tests ---

    def J_action_from_swap(self, S: np.ndarray, M: np.ndarray) -> np.ndarray:
        """Implement J M J^{-1} = S · M^* · S^T, where S is the L/R swap."""
        return S @ M.conj() @ S.T

    def test_first_order_condition(
        self, D_F: np.ndarray, ops: List[np.ndarray], labels: List[str], eps: float = 1e-12
    ) -> None:
        """First-order condition: [[D_F, a], J_F b J_F^{-1}] = 0 for all a,b."""
        n = D_F.shape[0]
        assert D_F.shape == (n, n)
        S = self.build_swap_LR(dim_left=n // 2)

        print("=== First-order condition test ===")
        max_norm = 0.0
        good_pairs = []

        for i, a in enumerate(ops):
            Da = D_F @ a - a @ D_F
            for j, b in enumerate(ops):
                b_tilde = self.J_action_from_swap(S, b)
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

    def test_zero_order_condition(
        self, ops: List[np.ndarray], labels: List[str], eps: float = 1e-12
    ) -> None:
        """Zero-order condition: [a, J_F b J_F^{-1}] = 0 for all a,b."""
        n = ops[0].shape[0]
        S = self.build_swap_LR(dim_left=n // 2)

        print("=== Zero-order condition test ===")
        max_norm = 0.0
        bad_pairs = []

        for i, a in enumerate(ops):
            for j, b in enumerate(ops):
                b_tilde = self.J_action_from_swap(S, b)
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

    def test_grading_and_reality(
        self, D_F: np.ndarray, ops: List[np.ndarray], labels: List[str]
    ) -> None:
        """
        Check grading and reality axioms:

          - γ_F anticommutes with D_F and commutes with A_F.
          - J_F^2 = 1 (swap^2 = I).
          - KO-dimension sign via J D_F J^{-1} = ± D_F.
        """
        n = D_F.shape[0]
        dpc = n // 2
        gamma_F = self.build_gamma_F(dpc)
        S = self.build_swap_LR(dpc)

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

    # ===========================
    # 2. Emergent misalignment model, graph, spectrum
    # ===========================

    def misalignment_energy(self, theta, w6=1.0, w5=1.0):
        N = len(theta)
        diffs = theta[:, None] - theta[None, :]
        cos6 = np.cos(6 * diffs)
        cos5 = np.cos(5 * diffs)
        E6 = w6 * np.sum(1.0 - cos6) / (N * N)
        E5 = w5 * np.sum(1.0 - cos5) / (N * N)
        return E6 + E5

    def relax_phases(self, N=200, n_steps=600, eta=0.01, w6=1.0, w5=1.0, random_seed=42):
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
                E = self.misalignment_energy(theta, w6=w6, w5=w5)
                energy_hist.append(E)

        return theta, energy_hist

    def build_emergent_adjacency(self, theta, w6=1.0, w5=1.0, keep_fraction=0.05):
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

    def largest_connected_component(self, A):
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

    def laplacian_from_adjacency(self, A):
        d = np.sum(A, axis=1)
        L = np.diag(d) - A
        return L

    def base_kernel(self, lam, alpha=3.0, form="lambda_sq"):
        """
        Scale-invariant base kernel F_base(λ_g):

            F_base(λ_g) = exp[-alpha * (λ_g / λ_ref)^2]

        with λ_ref = smallest positive eigenvalue in the triad.
        """
        lam = np.array(lam, dtype=float)

        lam_pos = lam[lam > 0]
        if lam_pos.size == 0:
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

    def spectral_triad(self, L):
        """
        Extract the 'generation triad' of internal Laplacian eigenvalues/eigenvectors.

        We:
          - diagonalize L,
          - take the first three NON-zero eigenvalues,
          - THEN reorder this triad so that generation index g=2 (3rd index)
            corresponds to the *smallest* λ, i.e. the heaviest generation
            after passing through a decreasing base kernel.
        """
        eigvals, eigvecs = np.linalg.eigh(L)
        idx_sorted = np.argsort(eigvals)
        eigvals_sorted = eigvals[idx_sorted]
        eigvecs_sorted = eigvecs[:, idx_sorted]

        # Triad of first three non-zero eigenvalues (ascending initially)
        triad_idx  = idx_sorted[1:4]
        triad_vals = eigvals_sorted[1:4]

        # Reorder so that lam_gen[2] is the smallest λ (→ largest Yukawa)
        order = np.argsort(triad_vals)[::-1]  # DESC by λ
        lam_gen = triad_vals[order]
        gen_indices = triad_idx[order]

        return lam_gen, gen_indices, eigvals_sorted

    # ===========================
    # 3. Sector charges, Yukawas, mixing
    # ===========================

    def build_sector_charges(self):
        # Example pattern that breaks degeneracy across sectors
        # (Just a *toy* suggestion; you can tune these later.)
        Q_u = np.array([0, 2, 4], dtype=float)  # keep as is for up
        Q_d = np.array([2, 4, 5], dtype=float)  # smaller gap to 3rd
        Q_e = np.array([4, 5, 6], dtype=float)  # even smaller gaps
        Q_nu = np.array([6, 7, 8], dtype=float)  # very compressed

        return {"u": Q_u, "d": Q_d, "e": Q_e, "nu": Q_nu}

    def sector_weights(self, F_base, Q_s, beta=1.0):
        return F_base * np.exp(-beta * Q_s)

    def mass_ratios(self, F_s):
        """
        Given a 3-component 'spectrum' F_s ~ (m1, m2, m3) up to an overall factor,
        return the ratios (m1/m3, m2/m3) with the convention that m3 is the
        largest eigenvalue. We:

          - take absolute values (masses are positive),
          - regularize very small entries to avoid zeros,
          - sort the three values so m1 <= m2 <= m3.

        This prevents enormous '1e16' ratios that are purely an ordering artifact.
        """
        F_s = np.array(F_s, dtype=float)
        F_s = np.abs(F_s)

        max_val = np.max(F_s)
        if max_val <= 0.0 or not np.isfinite(max_val):
            # Completely degenerate/vanishing or pathological spectrum:
            # treat as (1,1,1) so the fit will assign a large chi^2 and reject it.
            return 1.0, 1.0

        # Floor very small entries so none are exactly zero
        eps = 1e-16 * max_val
        F_s[F_s < eps] = eps

        # Sort so that m3 is the largest
        m1, m2, m3 = np.sort(F_s)

        return m1 / m3, m2 / m3
    # --- generation operators ---

    def rotation_3d(self, i, j, theta):
        R = np.eye(3, dtype=complex)
        c = np.cos(theta)
        s = np.sin(theta)
        R[i, i] = c
        R[j, j] = c
        R[i, j] = s
        R[j, i] = -s
        return R

    def build_generation_operators(self, phi_order=5, cab_denom=28):
        theta_phi = 2 * np.pi / phi_order
        theta_C = 2 * np.pi / cab_denom
        P_phi_12 = self.rotation_3d(0, 1, theta_phi)
        P_phi_23 = self.rotation_3d(1, 2, theta_phi)
        C_12 = self.rotation_3d(0, 1, theta_C)
        return P_phi_12, P_phi_23, C_12, theta_phi, theta_C

    # --- geometric regions and unitaries ---

    def build_geometric_regions(self, theta, n_regions=3):
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

    def build_geometric_unitary(self, gen_vecs, region_list):
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

    def build_sector_bases(self, P_phi_12, P_phi_23, C_12, U_geom,
                           use_neutrino_dressing: bool = True,
                           N_SOLAR: int = 36,
                           N_REACTOR: int = 45,
                           N_ATM: int = 24):
        """
        Build left/right generation bases for each sector.

        N_SOLAR, N_REACTOR, N_ATM set discrete angles for solar, reactor,
        and atmospheric-like rotations in the neutrino sector:
          θ_solar  = 2π / N_SOLAR
          θ_reac   = 2π / N_REACTOR
          θ_atm    = 2π / N_ATM
        """
        sector_bases = {}

        U_geom_u = U_geom["u"]
        U_geom_d = U_geom["d"]
        U_geom_e = U_geom["e"]
        U_geom_nu = U_geom["nu"]

        # Quarks: Cabibbo on up-type only
        U_L_u = U_geom_u @ C_12.conj().T
        U_R_u = np.eye(3, dtype=complex)

        U_L_d = U_geom_d
        U_R_d = np.eye(3, dtype=complex)

        # Charged leptons: pure geometry
        U_L_e = U_geom_e
        U_R_e = np.eye(3, dtype=complex)

        # Neutrinos: geometry + golden + 3 discrete rotations
        if use_neutrino_dressing:
            theta_solar = 2 * np.pi / N_SOLAR
            theta_reac = 2 * np.pi / N_REACTOR
            theta_atm = 2 * np.pi / N_ATM

            R_solar = self.rotation_3d(0, 1, theta_solar)
            R_reac = self.rotation_3d(0, 2, theta_reac)
            R_atm = self.rotation_3d(1, 2, theta_atm)

            U_dress = R_atm @ P_phi_23 @ R_solar @ P_phi_12 @ R_reac
            U_L_nu = U_geom_nu @ U_dress
        else:
            U_L_nu = U_geom_nu

        U_R_nu = np.eye(3, dtype=complex)

        sector_bases["u"] = (U_L_u, U_R_u)
        sector_bases["d"] = (U_L_d, U_R_d)
        sector_bases["e"] = (U_L_e, U_R_e)
        sector_bases["nu"] = (U_L_nu, U_R_nu)

        return sector_bases

    # --- Yukawas, mixing, chi² ---

    def yukawa_from_F_and_UL(self, F_s, U_L, U_R):
        D = np.diag(F_s)
        return U_L @ D @ U_R.conj().T

    def mixing_matrix(self, U_L_up, U_L_down):
        return U_L_up.conj().T @ U_L_down

    def mixing_angles_from_U(self, U):
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

    def compute_observables(
        self,
        mu_mt, mc_mt, md_mb, ms_mb, me_mt, mmu_mt,
        theta12_q, theta23_q, theta13_q,
        theta12_l, theta23_l, theta13_l
    ):
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

    def chi2(self, obs, targets=None):
        if targets is None:
            targets = self.TARGETS
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


class EmergentFlavorNCGModel(FlavorNCGOperators):
    """
    Subclass that:
      - configures the emergent graph / misalignment parameters
      - runs the full pipeline
      - constructs the internal NCG triple from Yukawas
      - runs NCG tests
    """

    def __init__(
        self,
        N_sites: int = 200,
        n_steps: int = 600,
        eta: float = 0.01,
        w6: float = 1.0,
        w5: float = 1.0,
        keep_fraction: float = 0.05,
    ):
        super().__init__()
        self.N_sites = N_sites
        self.n_steps = n_steps
        self.eta = eta
        self.w6 = w6
        self.w5 = w5
        self.keep_fraction = keep_fraction

    def run(self):
        # Step 1: relax phases under misalignment functional
        theta_final, energy_hist = self.relax_phases(
            N=self.N_sites,
            n_steps=self.n_steps,
            eta=self.eta,
            w6=self.w6,
            w5=self.w5,
            random_seed=42,
        )
        print("Relaxation complete.")
        print(f"Final misalignment energy: {energy_hist[-1]:.6f}")
        print()

        # Step 2: build emergent adjacency and Laplacian
        A_int_full = self.build_emergent_adjacency(
            theta_final,
            w6=self.w6,
            w5=self.w5,
            keep_fraction=self.keep_fraction,
        )
        A_int, nodes = self.largest_connected_component(A_int_full)
        L_int = self.laplacian_from_adjacency(A_int)

        # ---------------------------
        # Spectrum, rescaling, triads
        # ---------------------------

        # Full spectrum of the internal Laplacian
        eigvals_full_raw, eigvecs_full = np.linalg.eigh(L_int)

        # Global rescale factor for eigenvalues (unit of internal distance).
        # You can tune this; 1.0 reproduces the raw spectrum.
        L_rescale_factor = 1.0
        lam = L_rescale_factor * eigvals_full_raw

        # Choose distinct generation triads for quarks vs leptons
        triad_quark, triad_lepton = choose_quark_and_lepton_triads(
            lam, max_triad_index=90
        )

        lam_gen_quark = lam[triad_quark]
        lam_gen_lepton = lam[triad_lepton]

        # Universal base kernel evaluated on each triad
        # F_base_quark = self.base_kernel(lam_gen_quark, alpha=3.0, form="lambda_sq")
        # F_base_lepton = self.base_kernel(lam_gen_lepton, alpha=3.0, form="lambda_sq")
        alpha_internal = 300.0
        F_base_quark = self.base_kernel(lam_gen_quark, alpha=alpha_internal, form="lambda_sq")
        F_base_lepton = self.base_kernel(lam_gen_lepton, alpha=alpha_internal, form="lambda_sq")

        # Regularize extremely small kernel entries to avoid exact zeros
        def regularize_F_base(F):
            F = np.array(F, dtype=float)
            max_val = np.max(F)
            if max_val <= 0.0 or not np.isfinite(max_val):
                # pathological; return something tiny but nonzero
                return np.full_like(F, 1e-16)
            eps = 1e-16 * max_val
            F[F < eps] = eps
            return F

        F_base_quark = regularize_F_base(F_base_quark)
        F_base_lepton = regularize_F_base(F_base_lepton)

        print("=== Emergent internal graph ===")
        print(f"Number of sites: {A_int.shape[0]}")
        print("First 10 eigenvalues of L_int (raw, unscaled):")
        print(eigvals_full_raw[:10])
        print()
        print("Laplacian rescale factor L_rescale_factor =", L_rescale_factor)
        print("Quark generation eigenvalue indices:", triad_quark)
        print("Lepton generation eigenvalue indices:", triad_lepton)
        print("Quark lam_gen:", lam_gen_quark)
        print("Lepton lam_gen:", lam_gen_lepton)
        print("Base kernel F_base_quark(lam_gen):", F_base_quark)
        print("Base kernel F_base_lepton(lam_gen):", F_base_lepton)
        print()

        # Generation eigenvectors for each sector family
        gen_vecs_quark = eigvecs_full[:, triad_quark]
        gen_vecs_lepton = eigvecs_full[:, triad_lepton]

        # Step 3: build geometric regions from phase field, restricted to largest component
        theta_sub = theta_final[nodes]
        regions = self.build_geometric_regions(theta_sub, n_regions=3)
        R0, R1, R2 = regions

        # Quark assignments keep shared region geometry
        assign_u = [R0, R1, R2]
        assign_d = [R0, R1, R2]

        # Best-fit bookkeeping
        best_chi2 = np.inf
        best_perm_e = None
        best_perm_nu = None
        best_U_geom = None
        best_obs = None
        best_masses = None
        best_angles = None
        best_Ys = None
        best_sector_bases = None
        best_details = None

        sector_charges_gen = self.build_sector_charges()

        perms = [
            (0, 1, 2),
            (0, 2, 1),
            (1, 0, 2),
            (1, 2, 0),
            (2, 0, 1),
            (2, 1, 0),
        ]

        # Scan over lepton region permutations, with quark geometry fixed
        for pe in perms:
            for pn in perms:
                perm_e = [regions[pe[0]], regions[pe[1]], regions[pe[2]]]
                perm_n = [regions[pn[0]], regions[pn[1]], regions[pn[2]]]

                assign_e = perm_e
                assign_nu = perm_n

                # ---------------------------
                # Geometric unitaries per sector
                #   - quarks ⟶ quark triad
                #   - leptons ⟶ lepton triad
                # ---------------------------
                U_geom = {
                    "u": self.build_geometric_unitary(gen_vecs_quark, assign_u),
                    "d": self.build_geometric_unitary(gen_vecs_quark, assign_d),
                    "e": self.build_geometric_unitary(gen_vecs_lepton, assign_e),
                    "nu": self.build_geometric_unitary(gen_vecs_lepton, assign_nu),
                }

                # Sector weights from their respective base kernels
                F_u = self.sector_weights(F_base_quark, sector_charges_gen["u"], beta=1.0)
                F_d = self.sector_weights(F_base_quark, sector_charges_gen["d"], beta=1.0)
                F_e = self.sector_weights(F_base_lepton, sector_charges_gen["e"], beta=1.0)
                F_n = self.sector_weights(F_base_lepton, sector_charges_gen["nu"], beta=1.0)

                # generation operators (Cabibbo + golden)
                P_phi_12, P_phi_23, C_12, theta_phi, theta_C = self.build_generation_operators(
                    phi_order=5, cab_denom=28
                )

                # Sector bases: left/right unitaries including operator dressing
                sector_bases = self.build_sector_bases(
                    P_phi_12, P_phi_23, C_12,
                    U_geom,
                    use_neutrino_dressing=True,
                    N_SOLAR=36,
                    N_REACTOR=45,
                    N_ATM=24,
                )

                U_L_u, U_R_u = sector_bases["u"]
                U_L_d, U_R_d = sector_bases["d"]
                U_L_e, U_R_e = sector_bases["e"]
                U_L_nu, U_R_nu = sector_bases["nu"]

                # Yukawas from F_s and sector bases
                Y_u = self.yukawa_from_F_and_UL(F_u, U_L_u, U_R_u)
                Y_d = self.yukawa_from_F_and_UL(F_d, U_L_d, U_R_d)
                Y_e = self.yukawa_from_F_and_UL(F_e, U_L_e, U_R_e)
                Y_nu = self.yukawa_from_F_and_UL(F_n, U_L_nu, U_R_nu)

                # Mass ratios from F_s
                mu_mt, mc_mt = self.mass_ratios(F_u)
                md_mb, ms_mb = self.mass_ratios(F_d)
                me_mt, mmu_mt = self.mass_ratios(F_e)

                # Mixing matrices
                V_ckm = self.mixing_matrix(U_L_u, U_L_d)
                U_pmns = self.mixing_matrix(U_L_e, U_L_nu)

                theta12_q, theta23_q, theta13_q = self.mixing_angles_from_U(V_ckm)
                theta12_l, theta23_l, theta13_l = self.mixing_angles_from_U(U_pmns)

                obs = self.compute_observables(
                    mu_mt, mc_mt, md_mb, ms_mb, me_mt, mmu_mt,
                    theta12_q, theta23_q, theta13_q,
                    theta12_l, theta23_l, theta13_l,
                )
                chi2_value, chi2_details = self.chi2(obs)

                # Skip pathological points (NaN / inf chi^2)
                if not np.isfinite(chi2_value):
                    continue

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

        # Sanity check: did we actually find at least one finite-chi^2 point?
        if best_masses is None:
            raise RuntimeError(
                "No finite-chi^2 configuration found. "
                "Check base_kernel steepness and Laplacian rescale factor."
            )

        # ---------------------------
        # Unpack best solution
        # ---------------------------
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

        U_L_u, U_R_u = sector_bases["u"]
        U_L_d, U_R_d = sector_bases["d"]
        U_L_e, U_R_e = sector_bases["e"]
        U_L_nu, U_R_nu = sector_bases["nu"]

        V_ckm = self.mixing_matrix(U_L_u, U_L_d)
        U_pmns = self.mixing_matrix(U_L_e, U_L_nu)

        P_phi_12, P_phi_23, C_12, theta_phi, theta_C = self.build_generation_operators(
            phi_order=5, cab_denom=28
        )

        theta12_q, theta23_q, theta13_q = self.mixing_angles_from_U(V_ckm)
        theta12_l, theta23_l, theta13_l = self.mixing_angles_from_U(U_pmns)

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
        D_F = self.build_internal_DF_from_Y(Y_u, Y_d, Y_e, Y_nu)

        # Build a small internal algebra and test NCG axioms.
        ops_A, labels_A = self.build_internal_algebra_ops()
        self.test_first_order_condition(D_F, ops_A, labels_A, eps=1e-12)
        self.test_zero_order_condition(ops_A, labels_A, eps=1e-12)
        self.test_grading_and_reality(D_F, ops_A, labels_A)

        print("NOTES:")
        print("- The internal graph is emergent from the misalignment functional M[theta],")
        print("  which encodes 6-fold (C_360) and 5-fold (golden) alignment preferences.")
        print("- We then restrict to the largest connected component to define a single,")
        print("  coherent aether vacuum, and build its Laplacian L_int.")
        print("- The Laplacian spectrum is globally rescaled by L_rescale_factor before")
        print("  feeding its low modes into the universal base kernel F_base.")
        print("- Quark and lepton sectors select different generation triads from the same")
        print("  internal spectrum: quarks on one triad (small mixing), leptons on another")
        print("  (with potentially quasi-degenerate modes for large mixing).")
        print("- Sector hierarchies come from discrete integer charges Q_{s,g} acting on")
        print("  these F_base(lambda) values, and mixing from geometry-derived U_geom[s]")
        print("  combined with fixed operators P_phi (golden) and C_12 (Cabibbo).")
        print("- On top of this, we build an internal finite Dirac operator D_F from the")
        print("  same 3×3 Yukawa matrices, and an internal algebra generated by sector")
        print("  projectors and a sector charge Q_sector.")
        print("- This internal triple satisfies the NCG zero-order, first-order, grading,")
        print("  and reality (KO-sign) conditions relative to that algebra, giving us a")
        print("  self-consistent toy NCG-flavor sector driven by an emergent graph.")


if __name__ == "__main__":
    model = EmergentFlavorNCGModel()
    model.run()