import numpy as np
from typing import List, Tuple, Dict
from itertools import combinations

# ============================================================
#  Helpers: D_360, angle quantization, harmonic triad scoring
# ============================================================

def divisors_360() -> np.ndarray:
    """Divisors of 360 used as the allowed harmonic alphabet."""
    return np.array([1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 15, 18,
                     20, 24, 30, 36, 40, 45, 60, 72, 90, 120, 180, 360])


def nearest_divisor_angle(theta: float, divisors: np.ndarray = None) -> Tuple[float, int]:
    """
    Project a mixing angle theta onto the nearest divisor angle 2π/N
    with N ∈ D_360. Returns (theta_proj, N).
    """
    if divisors is None:
        divisors = divisors_360()
    theta = float(theta)
    candidates = 2.0 * np.pi / divisors
    idx = int(np.argmin(np.abs(candidates - theta)))
    return candidates[idx], int(divisors[idx])


def triad_harmonic_score(lam: np.ndarray, triad: Tuple[int, int, int],
                         mode: str = "quark") -> float:
    """
    Harmonic triad score based purely on spectral ratios.

    mode="quark": favor simple rational spacing (1:2, 2:3, etc).
    mode="lepton": favor near-degenerate pair + separated third, with φ-like ratio.
    """
    i, j, k = triad
    lam_i, lam_j, lam_k = lam[i], lam[j], lam[k]
    if lam_k <= 0:
        return -np.inf

    # Sort by value to get ordered gaps
    vals = np.array([lam_i, lam_j, lam_k])
    order = np.argsort(vals)
    li, lj, lk = vals[order]
    span = lk - li
    if span <= 0:
        return -np.inf

    g1 = lj - li
    g2 = lk - lj
    r1 = g1 / span
    r2 = g2 / span

    # Simple rational targets for "even" spacing
    simple_ratios = np.array([1/3, 1/2, 2/3])
    # φ-like target for leptons
    phi = (1 + np.sqrt(5)) / 2.0
    phi_ratio = 1.0 / phi  # ≈ 0.618

    if mode == "quark":
        # Smallest distance of r1,r2 to simple rational fractions
        d1 = np.min(np.abs(simple_ratios - r1))
        d2 = np.min(np.abs(simple_ratios - r2))
        return - (d1 + d2)  # smaller distance = higher score

    elif mode == "lepton":
        # Favor near-degenerate pair (smallest gap) and φ-like split of span
        dgaps = np.array([g1, g2])
        small_gap = np.min(dgaps)
        large_gap = np.max(dgaps)
        if span <= 0:
            return -np.inf
        # degeneracy measure: prefer small small_gap/span
        deg_measure = small_gap / span
        # φ-like measure on normalized large gap
        large_ratio = large_gap / span
        dphi = np.abs(large_ratio - phi_ratio)
        # score: want deg_measure small and dphi small
        return - (deg_measure + dphi)

    else:
        return -np.inf


# ============================================================
# Internal Laplacian scaling config (kept for compatibility)
# ============================================================

class InternalSpectrumConfig:
    """
    Config for rescaling the internal Laplacian eigenvalues
    before feeding them into the universal spectral kernel F_base.
    In the fully emergent version, the effective rescale is derived
    directly from the spectrum and this class is not used for tuning.
    """
    L_rescale_factor: float = 0.3
    max_triad_index: int = 20


def rescale_laplacian_evals(lam_raw: np.ndarray,
                            cfg: InternalSpectrumConfig) -> np.ndarray:
    """Legacy helper; not used in the emergent run(), kept for compatibility."""
    return cfg.L_rescale_factor * lam_raw


# ============================================================
# Emergent triad chooser
# ============================================================

def choose_quark_and_lepton_triads(lam: np.ndarray,
                                   max_triad_index: int = 20):
    """
    Choose quark- and lepton-like triads purely by harmonic spectral criteria.

    - Quark triad: near-rational spacing (simple gap ratios).
    - Lepton triad: near-degenerate pair + separated third with φ-like ratio.
    """
    start = 1  # ignore exact zero mode at 0
    stop = min(max_triad_index, len(lam))
    nonzero_indices = np.arange(start, stop, dtype=int)

    if len(nonzero_indices) < 3:
        raise ValueError("Not enough nonzero eigenvalues to form triads.")

    triads = list(combinations(nonzero_indices, 3))

    best_q, best_q_score = None, -np.inf
    best_l, best_l_score = None, -np.inf

    for triad in triads:
        s_q = triad_harmonic_score(lam, triad, mode="quark")
        if s_q > best_q_score:
            best_q_score = s_q
            best_q = triad

        s_l = triad_harmonic_score(lam, triad, mode="lepton")
        if s_l > best_l_score:
            best_l_score = s_l
            best_l = triad

    triad_quark = np.array(best_q, dtype=int)
    triad_lepton = np.array(best_l, dtype=int)
    return triad_quark, triad_lepton


# ============================================================
#  FlavorNCGOperators: NCG side (largely unchanged)
# ============================================================

class FlavorNCGOperators:
    """
    Master class collecting:
      - Emergent misalignment + internal graph machinery
      - Flavor hierarchy and mixing operators
      - Internal NCG finite Dirac operator + algebra + tests
    """

    SECTORS: List[str] = ["u", "d", "e", "nu"]
    N_GEN: int = 3
    SECTOR_NC: Dict[str, int] = {"u": 3, "d": 3, "e": 1, "nu": 1}

    # Rough SM targets kept ONLY as an external diagnostic (not used for selection)
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
        return len(self.SECTORS) * self.N_GEN  # 4 * 3 = 12

    def flavor_block_offsets(self) -> Dict[str, int]:
        off: Dict[str, int] = {}
        off["u"]  = 0
        off["d"]  = 3
        off["e"]  = 6
        off["nu"] = 9
        return off

    def build_internal_DF_from_Y(self, Y_u, Y_d, Y_e, Y_nu):
        for name, Y in [("Y_u", Y_u), ("Y_d", Y_d), ("Y_e", Y_e), ("Y_nu", Y_nu)]:
            arr = np.asarray(Y, dtype=complex)
            if arr.shape != (3, 3):
                raise ValueError(f"{name} must be a 3×3 matrix, got shape {arr.shape}.")

        Y_u = np.asarray(Y_u, dtype=complex)
        Y_d = np.asarray(Y_d, dtype=complex)
        Y_e = np.asarray(Y_e, dtype=complex)
        Y_nu = np.asarray(Y_nu, dtype=complex)

        dpc = self.dim_per_chirality()
        dimH = 2 * dpc

        Y_gen = np.zeros((dpc, dpc), dtype=complex)
        gen_off = self.flavor_block_offsets()

        def insert_sector_Y(sector: str, Y_s: np.ndarray) -> None:
            off = gen_off[sector]
            Y_gen[off:off + 3, off:off + 3] = Y_s

        insert_sector_Y("u", Y_u)
        insert_sector_Y("d", Y_d)
        insert_sector_Y("e", Y_e)
        insert_sector_Y("nu", Y_nu)

        Y_block = Y_gen

        D_F = np.zeros((dimH, dimH), dtype=complex)
        D_F[:dpc, dpc:] = Y_block.conj().T
        D_F[dpc:, :dpc] = Y_block
        return D_F

    # --- Real structure & grading ---

    def build_swap_LR(self, dim_left: int) -> np.ndarray:
        S = np.zeros((2 * dim_left, 2 * dim_left), dtype=complex)
        S[:dim_left, dim_left:] = np.eye(dim_left)
        S[dim_left:, :dim_left] = np.eye(dim_left)
        return S

    def build_gamma_F(self, dim_left: int) -> np.ndarray:
        g = np.zeros((2 * dim_left, 2 * dim_left), dtype=complex)
        g[:dim_left, :dim_left] = -np.eye(dim_left)
        g[dim_left:, dim_left:] =  np.eye(dim_left)
        return g

    def build_sector_projectors(self):
        dpc = self.dim_per_chirality()
        dimH = 2 * dpc
        gen_off = self.flavor_block_offsets()

        P: Dict[str, np.ndarray] = {}
        for s in self.SECTORS:
            P_s = np.zeros((dimH, dimH), dtype=complex)
            off = gen_off[s]
            P_s[off:off+3, off:off+3] = np.eye(3)
            P_s[dpc + off:dpc + off + 3, dpc + off:dpc + off + 3] = np.eye(3)
            P[s] = P_s
        return P

    def build_Q_sector(self) -> np.ndarray:
        """
        Simple sector charge operator distinguishing u,d,e,nu.
        In the emergent scheme, these sector charges are fixed only
        at this operator level; generation-wise structure is emergent.
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
        dpc = self.dim_per_chirality()
        dimH = 2 * dpc

        I = np.eye(dimH, dtype=complex)
        Q = self.build_Q_sector()
        P = self.build_sector_projectors()

        ops: List[np.ndarray] = [I, Q, P["u"], P["d"], P["e"], P["nu"]]
        labels: List[str] = ["I", "Q_sector", "P_sector_u", "P_sector_d", "P_sector_e", "P_sector_nu"]
        return ops, labels

    # --- NCG tests and alignment score ---

    def J_action_from_swap(self, S: np.ndarray, M: np.ndarray) -> np.ndarray:
        return S @ M.conj() @ S.T

    def test_first_order_condition(
        self, D_F: np.ndarray, ops: List[np.ndarray], labels: List[str], eps: float = 1e-12
    ) -> None:
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

    def ncg_alignment_score(self, D_F: np.ndarray, ops: List[np.ndarray]) -> float:
        """
        Scalar NCG coherence measure (smaller = more aligned).
        Combines grading, zero-order, and first-order deviations.
        """
        n = D_F.shape[0]
        dpc = n // 2
        gamma_F = self.build_gamma_F(dpc)
        S = self.build_swap_LR(dpc)

        # Grading: {γ, D_F} ≈ 0
        anti = gamma_F @ D_F + D_F @ gamma_F
        norm_anti = np.linalg.norm(anti, ord='fro')

        # gamma commutators with algebra
        max_comm_gamma = 0.0
        for a in ops:
            comm_ga = gamma_F @ a - a @ gamma_F
            max_comm_gamma = max(max_comm_gamma, np.linalg.norm(comm_ga, ord="fro"))

        # zero- and first-order
        max_zero = 0.0
        max_first = 0.0
        for i, a in enumerate(ops):
            Da = D_F @ a - a @ D_F
            for j, b in enumerate(ops):
                b_tilde = self.J_action_from_swap(S, b)
                # zero-order
                comm0 = a @ b_tilde - b_tilde @ a
                max_zero = max(max_zero, np.linalg.norm(comm0, ord="fro"))
                # first-order
                comm2 = Da @ b_tilde - b_tilde @ Da
                max_first = max(max_first, np.linalg.norm(comm2, ord="fro"))

        # Simple linear combo; no tunable weights beyond unity
        return norm_anti + max_comm_gamma + max_zero + max_first

    # ===========================
    # 2. Emergent misalignment model, graph, spectrum
    # ===========================

    def allowed_harmonics(self) -> np.ndarray:
        """Allowed global harmonic set D_360."""
        return divisors_360()

    def contextual_harmonics(self, step: int, total_steps: int) -> np.ndarray:
        """
        Contextual selection of subset of D_360 as relaxation proceeds.
        Early time: small subset; late time: full set.
        """
        D = self.allowed_harmonics()
        frac = step / max(total_steps, 1)
        k = int(1 + frac * (len(D) - 1))
        return D[:k]

    def misalignment_energy(self, theta, ns: np.ndarray = None):
        if ns is None:
            ns = self.allowed_harmonics()
        N = len(theta)
        diffs = theta[:, None] - theta[None, :]
        E = 0.0
        for n in ns:
            w_n = 1.0 / n
            E += w_n * np.sum(1.0 - np.cos(n * diffs)) / (N * N)
        return E

    def relax_phases(self, N=200, n_steps=600, eta=0.01, random_seed=42):
        rng = np.random.default_rng(random_seed)
        theta = rng.uniform(0, 2 * np.pi, size=N)
        energy_hist = []

        for step in range(n_steps):
            ns = self.contextual_harmonics(step, n_steps)
            diffs = theta[:, None] - theta[None, :]
            grad = np.zeros(N, dtype=float)

            for n in ns:
                w_n = 1.0 / n
                sin_n = np.sin(n * diffs)
                grad += w_n * n * np.sum(sin_n, axis=1)

            theta = theta - eta * grad
            theta = (theta + 2 * np.pi) % (2 * np.pi)

            if step % 10 == 0 or step == n_steps - 1:
                E = self.misalignment_energy(theta, ns=ns)
                energy_hist.append(E)

        return theta, energy_hist

    def build_emergent_adjacency(self, theta, ns: np.ndarray = None, keep_fraction: float = 0.05):
        """
        Adjacency from the same harmonic set ns used at late-time misalignment.
        Score_ij = Σ_n (1/n) cos(n(θ_i - θ_j)).
        """
        if ns is None:
            ns = self.allowed_harmonics()

        N = len(theta)
        diffs = theta[:, None] - theta[None, :]
        score = np.zeros((N, N), dtype=float)

        for n in ns:
            w_n = 1.0 / n
            score += w_n * np.cos(n * diffs)

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
        Universal base kernel F_base(λ_g):

            F_base(λ_g) = exp[-alpha * (λ_g / λ_ref)^p]

        with λ_ref = smallest positive eigenvalue in the triad.
        alpha will be emergently set from triad spread.
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

    def emergent_alpha_for_triad(self, lam_triad: np.ndarray) -> float:
        """
        Derive kernel steepness from the triad itself.
        Use spread in log(λ) to set alpha ~ 1 / Var(log λ).
        """
        lam = np.array(lam_triad, dtype=float)
        lam_pos = lam[lam > 0]
        if lam_pos.size <= 1:
            return 1.0
        logs = np.log(lam_pos)
        var = np.var(logs)
        eps = 1e-6
        alpha = 1.0 / (var + eps)
        return alpha

    def spectral_triad(self, L):
        eigvals, eigvecs = np.linalg.eigh(L)
        idx_sorted = np.argsort(eigvals)
        eigvals_sorted = eigvals[idx_sorted]
        eigvecs_sorted = eigvecs[:, idx_sorted]

        triad_idx = idx_sorted[1:4]
        triad_vals = eigvals_sorted[1:4]

        order = np.argsort(triad_vals)[::-1]  # DESC by λ
        lam_gen = triad_vals[order]
        gen_indices = triad_idx[order]
        return lam_gen, gen_indices, eigvals_sorted

    # ===========================
    # 3. Sector charges, Yukawas, mixing
    # ===========================

    def build_sector_charges_from_spectrum(self, lam: np.ndarray,
                                           triad_quark: np.ndarray,
                                           triad_lepton: np.ndarray) -> Dict[str, np.ndarray]:
        """
        Emergent sector/generation charges from local spectral density
        around each triad eigenvalue.
        """
        lam = np.array(lam, dtype=float)

        def local_density(idx: int) -> float:
            v = lam[idx]
            if v <= 0:
                return 1.0
            window = (lam >= 0.9 * v) & (lam <= 1.1 * v)
            return float(np.sum(window))

        def triad_charges(triad: np.ndarray) -> np.ndarray:
            qs = np.array([local_density(int(i)) for i in triad], dtype=float)
            # log compress
            return np.log1p(qs)

        Q_quark = triad_charges(triad_quark)
        Q_lepton = triad_charges(triad_lepton)

        charges = {
            "u":  Q_quark,
            "d":  Q_quark,
            "e":  Q_lepton,
            "nu": Q_lepton,
        }
        return charges

    def sector_weights(self, F_base: np.ndarray, Q_s: np.ndarray):
        """
        No free β: masses ~ F_base * exp(-Q_s).
        """
        return F_base * np.exp(-Q_s)

    def mass_ratios(self, F_s):
        F_s = np.array(F_s, dtype=float)
        F_s = np.abs(F_s)
        max_val = np.max(F_s)
        if max_val <= 0.0 or not np.isfinite(max_val):
            return 1.0, 1.0
        eps = 1e-16 * max_val
        F_s[F_s < eps] = eps
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
        """
        In the emergent scheme, phi_order and cab_denom are NOT free:
        they are derived from geometric mixing and projected to the
        nearest divisor-based angles before calling this.
        """
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

    def emergent_neutrino_denominators(self, lam_gen_lepton: np.ndarray) -> Tuple[int, int, int]:
        """
        Set N_SOLAR, N_REACTOR, N_ATM from lepton triad degeneracies.
        Smaller gap -> larger N (finer angle).
        """
        lam = np.array(lam_gen_lepton, dtype=float)
        if lam.size != 3:
            return 36, 45, 24

        gaps = np.abs(np.diff(np.sort(lam)))
        # Protect against zero
        gaps = gaps + 1e-8
        inv_gaps = 1.0 / gaps
        inv_gaps /= np.max(inv_gaps)

        # Map to a subset of divisors
        D = divisors_360()
        candidates = D[D <= 90]  # keep it modest

        def map_val(v):
            # v in [0,1] -> candidate index
            idx = int(np.clip(round(v * (len(candidates)-1)), 0, len(candidates)-1))
            return int(candidates[idx])

        N_SOLAR = map_val(inv_gaps[0])   # g12
        N_ATM   = map_val(inv_gaps[-1])  # g23
        N_REACTOR = map_val(0.5 * (inv_gaps[0] + inv_gaps[-1]))
        return N_SOLAR, N_REACTOR, N_ATM

    # --- Yukawas, mixing, diagnostics ---

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


# ============================================================
# EmergentFlavorNCGModel: FULL EMERGENCE RUN PIPELINE
# ============================================================

class EmergentFlavorNCGModel(FlavorNCGOperators):
    def __init__(
        self,
        N_sites: int = 200,
        n_steps: int = 600,
        eta: float = 0.01,
        keep_fraction: float = 0.05,
    ):
        super().__init__()
        self.N_sites = N_sites
        self.n_steps = n_steps
        self.eta = eta
        self.keep_fraction = keep_fraction

    def run(self):
        # Step 1: relax phases under D_360-driven misalignment
        theta_final, energy_hist = self.relax_phases(
            N=self.N_sites,
            n_steps=self.n_steps,
            eta=self.eta,
            random_seed=42,
        )
        print("Relaxation complete.")
        print(f"Final misalignment energy: {energy_hist[-1]:.6f}")
        print()

        # Harmonics active at final time
        ns_final = self.contextual_harmonics(self.n_steps - 1, self.n_steps)

        # Step 2: emergent adjacency & Laplacian
        A_int_full = self.build_emergent_adjacency(
            theta_final,
            ns=ns_final,
            keep_fraction=self.keep_fraction,
        )
        A_int, nodes = self.largest_connected_component(A_int_full)
        L_int = self.laplacian_from_adjacency(A_int)

        # Spectrum and emergent rescaling
        eigvals_full_raw, eigvecs_full = np.linalg.eigh(L_int)
        pos = eigvals_full_raw[eigvals_full_raw > 1e-12]
        if pos.size > 0:
            L_rescale_factor = 1.0 / pos[0]
        else:
            L_rescale_factor = 1.0
        lam = L_rescale_factor * eigvals_full_raw

        # Emergent triads from harmonic scoring
        triad_quark, triad_lepton = choose_quark_and_lepton_triads(
            lam, max_triad_index=min(90, len(lam))
        )
        lam_gen_quark = lam[triad_quark]
        lam_gen_lepton = lam[triad_lepton]

        # Emergent alpha from triad spread
        alpha_quark = self.emergent_alpha_for_triad(lam_gen_quark)
        alpha_lepton = self.emergent_alpha_for_triad(lam_gen_lepton)

        F_base_quark = self.base_kernel(lam_gen_quark, alpha=alpha_quark, form="lambda_sq")
        F_base_lepton = self.base_kernel(lam_gen_lepton, alpha=alpha_lepton, form="lambda_sq")

        def regularize_F_base(F):
            F = np.array(F, dtype=float)
            max_val = np.max(F)
            if max_val <= 0.0 or not np.isfinite(max_val):
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
        print("Quark triad indices:", triad_quark, "lam_gen_quark:", lam_gen_quark)
        print("Lepton triad indices:", triad_lepton, "lam_gen_lepton:", lam_gen_lepton)
        print("Alpha_quark (emergent):", alpha_quark)
        print("Alpha_lepton (emergent):", alpha_lepton)
        print("Base kernel F_base_quark:", F_base_quark)
        print("Base kernel F_base_lepton:", F_base_lepton)
        print()

        # Generation eigenvectors
        gen_vecs_quark = eigvecs_full[:, triad_quark]
        gen_vecs_lepton = eigvecs_full[:, triad_lepton]

        # Step 3: geometric regions from phase field (restricted to largest component)
        theta_sub = theta_final[nodes]
        regions = self.build_geometric_regions(theta_sub, n_regions=3)
        R0, R1, R2 = regions

        # Quark assignments share region geometry
        assign_u = [R0, R1, R2]
        assign_d = [R0, R1, R2]

        # Sector charges from spectrum
        sector_charges_gen = self.build_sector_charges_from_spectrum(
            lam,
            triad_quark=triad_quark,
            triad_lepton=triad_lepton,
        )

        # Emergent neutrino denominators from lepton triad degeneracies
        N_SOLAR, N_REACTOR, N_ATM = self.emergent_neutrino_denominators(lam_gen_lepton)
        print("Emergent neutrino denominators (SOLAR, REACTOR, ATM):", N_SOLAR, N_REACTOR, N_ATM)
        print()

        # Permutations for leptons (internal alignment selection only)
        perms = [
            (0, 1, 2),
            (0, 2, 1),
            (1, 0, 2),
            (1, 2, 0),
            (2, 0, 1),
            (2, 1, 0),
        ]

        best_align_score = np.inf
        best_perm_e = None
        best_perm_nu = None
        best_U_geom = None
        best_masses = None
        best_angles = None
        best_Ys = None
        best_sector_bases = None
        best_chi2 = None
        best_chi2_details = None

        # Build algebra once for NCG scoring (size known: 24x24)
        ops_A, labels_A = self.build_internal_algebra_ops()

        for pe in perms:
            for pn in perms:
                perm_e = [regions[pe[0]], regions[pe[1]], regions[pe[2]]]
                perm_n = [regions[pn[0]], regions[pn[1]], regions[pn[2]]]

                assign_e = perm_e
                assign_nu = perm_n

                # Geometric unitaries
                U_geom = {
                    "u": self.build_geometric_unitary(gen_vecs_quark, assign_u),
                    "d": self.build_geometric_unitary(gen_vecs_quark, assign_d),
                    "e": self.build_geometric_unitary(gen_vecs_lepton, assign_e),
                    "nu": self.build_geometric_unitary(gen_vecs_lepton, assign_nu),
                }

                # Pure geometric mixing
                V_ckm_geom = self.mixing_matrix(U_geom["u"], U_geom["d"])
                U_pmns_geom = self.mixing_matrix(U_geom["e"], U_geom["nu"])
                theta12_q_geom, theta23_q_geom, theta13_q_geom = self.mixing_angles_from_U(V_ckm_geom)
                theta12_l_geom, theta23_l_geom, theta13_l_geom = self.mixing_angles_from_U(U_pmns_geom)

                # Emergent Cabibbo and golden angles via divisor projection
                theta_C_proj, cab_denom = nearest_divisor_angle(theta12_q_geom)
                theta_phi_proj, phi_order = nearest_divisor_angle(theta12_l_geom)

                P_phi_12, P_phi_23, C_12, theta_phi, theta_C = self.build_generation_operators(
                    phi_order=phi_order, cab_denom=cab_denom
                )

                # Sector weights from spectrum and charges
                F_u = self.sector_weights(F_base_quark, sector_charges_gen["u"])
                F_d = self.sector_weights(F_base_quark, sector_charges_gen["d"])
                F_e = self.sector_weights(F_base_lepton, sector_charges_gen["e"])
                F_n = self.sector_weights(F_base_lepton, sector_charges_gen["nu"])

                # Sector bases: geometry + emergent operators
                sector_bases = self.build_sector_bases(
                    P_phi_12, P_phi_23, C_12,
                    U_geom,
                    use_neutrino_dressing=True,
                    N_SOLAR=N_SOLAR,
                    N_REACTOR=N_REACTOR,
                    N_ATM=N_ATM,
                )

                U_L_u, U_R_u = sector_bases["u"]
                U_L_d, U_R_d = sector_bases["d"]
                U_L_e, U_R_e = sector_bases["e"]
                U_L_nu, U_R_nu = sector_bases["nu"]

                # Yukawas from emergent F_s
                Y_u = self.yukawa_from_F_and_UL(F_u, U_L_u, U_R_u)
                Y_d = self.yukawa_from_F_and_UL(F_d, U_L_d, U_R_d)
                Y_e = self.yukawa_from_F_and_UL(F_e, U_L_e, U_R_e)
                Y_nu = self.yukawa_from_F_and_UL(F_n, U_L_nu, U_R_nu)

                # Mass ratios from F_s
                mu_mt, mc_mt = self.mass_ratios(F_u)
                md_mb, ms_mb = self.mass_ratios(F_d)
                me_mt, mmu_mt = self.mass_ratios(F_e)

                # Mixing matrices with dressed U_L
                V_ckm = self.mixing_matrix(U_L_u, U_L_d)
                U_pmns = self.mixing_matrix(U_L_e, U_L_nu)

                theta12_q, theta23_q, theta13_q = self.mixing_angles_from_U(V_ckm)
                theta12_l, theta23_l, theta13_l = self.mixing_angles_from_U(U_pmns)

                # Emergent alignment: angles close to divisor angles + NCG coherence
                # Angle errors to nearest divisor angles
                def angle_error(theta):
                    _, _N = nearest_divisor_angle(theta)
                    theta_proj, _ = nearest_divisor_angle(theta)
                    return abs(theta - theta_proj)

                angle_errors = (
                    angle_error(theta12_q) +
                    angle_error(theta23_q) +
                    angle_error(theta13_q) +
                    angle_error(theta12_l) +
                    angle_error(theta23_l) +
                    angle_error(theta13_l)
                )

                # NCG alignment score
                D_F = self.build_internal_DF_from_Y(Y_u, Y_d, Y_e, Y_nu)
                ncg_score = self.ncg_alignment_score(D_F, ops_A)

                align_score = angle_errors + ncg_score  # no external data used

                if align_score < best_align_score:
                    best_align_score = align_score
                    best_perm_e = pe
                    best_perm_nu = pn
                    best_U_geom = U_geom
                    best_masses = (mu_mt, mc_mt, md_mb, ms_mb, me_mt, mmu_mt)
                    best_angles = (theta12_q, theta23_q, theta13_q,
                                   theta12_l, theta23_l, theta13_l)
                    best_Ys = (Y_u, Y_d, Y_e, Y_nu)
                    best_sector_bases = sector_bases

                    # External diagnostic: SM χ² (NOT used to select)
                    obs = self.compute_observables(
                        mu_mt, mc_mt, md_mb, ms_mb, me_mt, mmu_mt,
                        theta12_q, theta23_q, theta13_q,
                        theta12_l, theta23_l, theta13_l,
                    )
                    chi2_value, chi2_details = self.chi2(obs)
                    best_chi2 = chi2_value
                    best_chi2_details = chi2_details

        if best_masses is None:
            raise RuntimeError("No emergent alignment configuration found.")

        # ---------------------------
        # Unpack best emergent solution
        # ---------------------------
        pe = best_perm_e
        pn = best_perm_nu
        U_geom = best_U_geom
        mu_mt, mc_mt, md_mb, ms_mb, me_mt, mmu_mt = best_masses
        theta12_q, theta23_q, theta13_q, theta12_l, theta23_l, theta13_l = best_angles
        Y_u, Y_d, Y_e, Y_nu = best_Ys
        sector_bases = best_sector_bases
        chi2_value = best_chi2
        chi2_details = best_chi2

        print("=== Emergent lepton region permutations (internal alignment only) ===")
        print(f"  pe (e sectors)  = {pe}")
        print(f"  pn (nu sectors) = {pn}")
        print(f"Best internal alignment score  ≈ {best_align_score:.3e}")
        print()

        print("Mass ratios (m1/m3, m2/m3) from emergent F_s:")
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

        # Reconstruct emergent Cabibbo / golden parameters for reporting
        V_ckm_geom = self.mixing_matrix(U_geom["u"], U_geom["d"])
        U_pmns_geom = self.mixing_matrix(U_geom["e"], U_geom["nu"])
        theta12_q_geom, theta23_q_geom, theta13_q_geom = self.mixing_angles_from_U(V_ckm_geom)
        theta12_l_geom, theta23_l_geom, theta13_l_geom = self.mixing_angles_from_U(U_pmns_geom)
        theta_C_proj, cab_denom = nearest_divisor_angle(theta12_q_geom)
        theta_phi_proj, phi_order = nearest_divisor_angle(theta12_l_geom)

        print("=== CKM-like mixing matrix (emergent geometry + operators) ===")
        print(V_ckm)
        print(f"theta12_q ≈ {theta12_q:.3f} rad, theta23_q ≈ {theta23_q:.3f}, theta13_q ≈ {theta13_q:.3e}")
        print(f"(Emergent Cabibbo: 2π/{cab_denom} ≈ {theta_C_proj:.3f} rad)")
        print()

        print("=== PMNS-like mixing matrix (emergent geometry + operators) ===")
        print(U_pmns)
        print(f"theta12_l ≈ {theta12_l:.3f} rad, theta23_l ≈ {theta23_l:.3f}, theta13_l ≈ {theta13_l:.3e}")
        print(f"(Emergent golden-like: 2π/{phi_order} ≈ {theta_phi_proj:.3f} rad)")
        print()

        # External diagnostic only
        obs = self.compute_observables(
            mu_mt, mc_mt, md_mb, ms_mb, me_mt, mmu_mt,
            theta12_q, theta23_q, theta13_q,
            theta12_l, theta23_l, theta13_l,
        )
        chi2_value, chi2_details = self.chi2(obs)

        print("=== Observables vs rough SM targets (diagnostic ONLY) ===")
        for k, m, t, contrib in chi2_details:
            print(f"{k:12s}: model={m:.3e}, target={t:.3e}, chi2_contrib={contrib:.2f}")
        print()
        print(f"Total diagnostic chi^2 ≈ {chi2_value:.2f}")
        print()

        # ===============================
        # Internal NCG triple from emergent Yukawas
        # ===============================
        D_F = self.build_internal_DF_from_Y(Y_u, Y_d, Y_e, Y_nu)

        # Internal algebra and NCG axiom checks (now emergent-consistent)
        ops_A, labels_A = self.build_internal_algebra_ops()
        self.test_first_order_condition(D_F, ops_A, labels_A, eps=1e-12)
        self.test_zero_order_condition(ops_A, labels_A, eps=1e-12)
        self.test_grading_and_reality(D_F, ops_A, labels_A)

        print("NOTES:")
        print("- Misalignment uses a context-dependent subset of D_360 harmonics only.")
        print("- The internal graph, Laplacian, and rescaling are entirely emergent from that harmonic engine.")
        print("- Quark and lepton triads are chosen by harmonic spectral criteria (rational vs φ-like spacing).")
        print("- Sector/generation charges Q_{s,g} come from local spectral density near each triad eigenvalue.")
        print("- Base-kernel steepness alpha is derived from the triad's log-spectrum variance.")
        print("- Cabibbo, golden, and neutrino rotation denominators are read off from geometric mixing")
        print("  and projected onto nearest divisor-based 2π/N angles.")
        print("- Region assignments for leptons are selected by internal alignment score (divisor-angle match")
        print("  + NCG coherence), not by fitting external SM data.")
        print("- SM targets are retained only as an external diagnostic chi^2 and do not feed back into")
        print("  the emergent vacuum selection.")
        print("- The internal NCG triple is built from the same emergent Yukawas and tested against the")
        print("  zero-order, first-order, grading, and reality axioms, providing a fully emergent,")
        print("  self-consistent toy NCG-flavor sector.")


if __name__ == "__main__":
    model = EmergentFlavorNCGModel()
    model.run()

"""
RESULTS:

Relaxation complete.
Final misalignment energy: 2.653294

=== Emergent internal graph ===
Number of sites: 177
First 10 eigenvalues of L_int (raw, unscaled):
[-1.51899334e-15  4.47196842e-03  1.39832968e-02  4.30280833e-02
  8.14410445e-02  1.05968431e-01  1.43320623e-01  2.52640915e-01
  4.08466369e-01  5.01017328e-01]

Laplacian rescale factor L_rescale_factor = 223.61517478622167
Quark triad indices: [23 50 79] lam_gen_quark: [ 894.46069914 1565.3062235  2236.15174786]
Lepton triad indices: [16 50 53] lam_gen_lepton: [ 469.45491797 1565.3062235  1569.78274107]
Alpha_quark (emergent): 7.03133486459314
Alpha_lepton (emergent): 3.09553932586928
Base kernel F_base_quark: [8.83751305e-04 4.44770355e-10 8.83751305e-20]
Base kernel F_base_lepton: [4.52506011e-02 1.13181752e-15 9.29323040e-16]

Emergent neutrino denominators (SOLAR, REACTOR, ATM): 1 15 90

=== Emergent lepton region permutations (internal alignment only) ===
  pe (e sectors)  = (2, 0, 1)
  pn (nu sectors) = (2, 0, 1)
Best internal alignment score  ≈ 4.440e-02

Mass ratios (m1/m3, m2/m3) from emergent F_s:
mu/mt:     1.000e-16, mc/mt:     2.745e-07
md/mb:     1.000e-16, ms/mb:     2.745e-07
me/mtau:   3.734e-15, mmu/mtau:  4.548e-15

=== CKM-like mixing matrix (emergent geometry + operators) ===
[[ 9.99847695e-01+0.j  1.74524064e-02+0.j -6.29704910e-18+0.j]
 [-1.74524064e-02+0.j  9.99847695e-01+0.j -2.08809152e-16+0.j]
 [ 0.00000000e+00+0.j -2.77555756e-16+0.j  1.00000000e+00+0.j]]
theta12_q ≈ 0.017 rad, theta23_q ≈ 0.000, theta13_q ≈ 6.297e-18
(Emergent Cabibbo: 2π/360 ≈ 0.017 rad)

=== PMNS-like mixing matrix (emergent geometry + operators) ===
[[ 0.91340632+0.j  0.01745241+0.j  0.4066747 +0.j]
 [-0.05133233+0.j  0.99604297+0.j  0.07254921+0.j]
 [-0.40379931+0.j -0.08714247+0.j  0.91068782+0.j]]
theta12_l ≈ 0.019 rad, theta23_l ≈ 0.079, theta13_l ≈ 4.188e-01
(Emergent golden-like: 2π/360 ≈ 0.017 rad)

=== Observables vs rough SM targets (diagnostic ONLY) ===
mu_mt       : model=1.000e-16, target=2.200e-05, chi2_contrib=4.00
mc_mt       : model=2.745e-07, target=7.500e-03, chi2_contrib=4.00
md_mb       : model=1.000e-16, target=1.100e-03, chi2_contrib=4.00
ms_mb       : model=2.745e-07, target=2.200e-02, chi2_contrib=4.00
me_mt       : model=3.734e-15, target=2.900e-04, chi2_contrib=4.00
mmu_mt      : model=4.548e-15, target=5.900e-02, chi2_contrib=4.00
theta12_q   : model=1.745e-02, target=2.270e-01, chi2_contrib=340.86
theta23_q   : model=2.088e-16, target=4.100e-02, chi2_contrib=4.00
theta13_q   : model=6.297e-18, target=3.600e-03, chi2_contrib=4.00
theta12_l   : model=1.910e-02, target=5.840e-01, chi2_contrib=93.56
theta23_l   : model=7.950e-02, target=7.850e-01, chi2_contrib=20.19
theta13_l   : model=4.188e-01, target=1.500e-01, chi2_contrib=80.29

Total diagnostic chi^2 ≈ 566.90

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
||J D_F J^-1 - D_F||_F   = 5.488e-02
||J D_F J^-1 + D_F||_F   = 7.196e-02
→ KO-sign: J D_F J^-1 = + D_F (J-even Dirac operator)

"""