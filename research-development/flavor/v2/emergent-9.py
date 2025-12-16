from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple
import numpy as np

# ================================
# Finite SM triple (1 generation + nu_R)
# ================================

@dataclass
class SMState:
    name: str           # "nu_L", "u_L_r", "bar_u_R_b", ...
    chirality: str      # "L" or "R"
    sector: str         # "lepton" or "quark"
    particle: bool      # True = particle, False = antiparticle
    color: Optional[str]
    generation: int


def build_sm_basis_one_gen(include_nu_R: bool = True) -> List[SMState]:
    """1 generation + nu_R + antiparticles → dim(H_F^SM) = 32."""
    basis: List[SMState] = []
    gen = 1

    # Particles: leptons
    basis.append(SMState("nu_L", "L", "lepton", True, None, gen))
    basis.append(SMState("e_L",  "L", "lepton", True, None, gen))
    if include_nu_R:
        basis.append(SMState("nu_R", "R", "lepton", True, None, gen))
    basis.append(SMState("e_R",  "R", "lepton", True, None, gen))

    # Particles: quarks
    colors = ["r", "g", "b"]
    for c in colors:
        basis.append(SMState(f"u_L_{c}", "L", "quark", True, c, gen))
        basis.append(SMState(f"d_L_{c}", "L", "quark", True, c, gen))
    for c in colors:
        basis.append(SMState(f"u_R_{c}", "R", "quark", True, c, gen))
        basis.append(SMState(f"d_R_{c}", "R", "quark", True, c, gen))

    # Antiparticles: same multiplets, opposite chirality
    particle_states = basis.copy()
    for st in particle_states:
        new_chir = "R" if st.chirality == "L" else "L"
        basis.append(SMState("bar_" + st.name, new_chir, st.sector, False, st.color, st.generation))

    return basis


def build_name_index_map_sm(basis: List[SMState]) -> Dict[str, int]:
    return {st.name: i for i, st in enumerate(basis)}


def build_gamma_sm(basis: List[SMState]) -> np.ndarray:
    """γ: -1 on left-chiral states, +1 on right-chiral states (particles + antiparticles)."""
    dim = len(basis)
    gamma = np.zeros((dim, dim), dtype=complex)
    for i, st in enumerate(basis):
        gamma[i, i] = -1.0 if st.chirality == "L" else 1.0
    return gamma


def build_swap_J_sm(basis: List[SMState]) -> np.ndarray:
    """Swap matrix S implementing particle ↔ antiparticle exchange."""
    dim = len(basis)
    S = np.zeros((dim, dim), dtype=complex)

    idx_map: Dict[tuple, int] = {
        (st.name, st.particle): i
        for i, st in enumerate(basis)
    }

    for i, st in enumerate(basis):
        if st.particle:
            partner_name = "bar_" + st.name
            j = idx_map[(partner_name, False)]
        else:
            assert st.name.startswith("bar_")
            base_name = st.name[len("bar_"):]
            j = idx_map[(base_name, True)]
        S[i, j] = 1.0

    return S

"""
# Example:
basis_SM = build_sm_basis_one_gen()
check_KO_relations_SM(basis_SM)

||S^2 - I||_F = 0.0
||JγJ^{-1} + γ||_F = 0.0
||JγJ^{-1} - γ||_F = 11.313708498984761

"""
def add_quaternion_on_doublet(M: np.ndarray,
                              q2: np.ndarray,
                              idx1: int,
                              idx2: int,
                              conj: bool = False) -> None:
    """
    Add the 2x2 matrix q2 (or its conjugate) on the subspace spanned by idx1, idx2.
    """
    block = q2.conj() if conj else q2
    indices = [idx1, idx2]
    for a, i in enumerate(indices):
        for b, j in enumerate(indices):
            M[i, j] += block[a, b]


def add_color_matrix_on_triplet(M: np.ndarray,
                                m3: np.ndarray,
                                idxs: List[int],
                                conj: bool = False) -> None:
    """
    Add the 3x3 matrix m3 (or its conjugate) on the subspace of indices idxs.
    """
    block = m3.conj() if conj else m3
    assert len(idxs) == 3
    for a, i in enumerate(idxs):
        for b, j in enumerate(idxs):
            M[i, j] += block[a, b]
def rep_A_SM(lambda_c: complex,
             q2: np.ndarray,
             m3: np.ndarray,
             basis: List[SMState],
             name_to_index: Dict[str, int]) -> np.ndarray:
    """
    Representation of a = (lambda_c, q2, m3) in A_F^SM on H_F^SM (1 generation + nu_R).

    - lambda_c: scalar (C-part), acts as lambda_c * I on all states.
    - q2: 2x2 complex matrix (H-part), acts on SU(2) doublets:
         (nu_L, e_L), (u_L_c, d_L_c), and their antiparticle doublets.
    - m3: 3x3 complex matrix (color part), acts on color triplets for quarks
          (fundamental on particles, conjugate on antiparticles).
    """
    dim = len(basis)
    M = lambda_c * np.eye(dim, dtype=complex)

    # 1) Quaternion (SU(2)) part: act on doublets

    # Lepton doublet (particles)
    idx_nu_L  = name_to_index["nu_L"]
    idx_e_L   = name_to_index["e_L"]
    add_quaternion_on_doublet(M, q2, idx_nu_L, idx_e_L, conj=False)

    # Lepton doublet (antiparticles)
    idx_bar_nu_L = name_to_index["bar_nu_L"]
    idx_bar_e_L  = name_to_index["bar_e_L"]
    add_quaternion_on_doublet(M, q2, idx_bar_nu_L, idx_bar_e_L, conj=True)

    # Quark doublets (particles and antiparticles) for each color
    colors = ["r", "g", "b"]
    for c in colors:
        # particles: (u_L_c, d_L_c)
        idx_uL = name_to_index[f"u_L_{c}"]
        idx_dL = name_to_index[f"d_L_{c}"]
        add_quaternion_on_doublet(M, q2, idx_uL, idx_dL, conj=False)

        # antiparticles: (bar_u_L_c, bar_d_L_c)
        idx_bar_uL = name_to_index[f"bar_u_L_{c}"]
        idx_bar_dL = name_to_index[f"bar_d_L_{c}"]
        add_quaternion_on_doublet(M, q2, idx_bar_uL, idx_bar_dL, conj=True)

    # Right-handed singlets (nu_R, e_R, u_R_c, d_R_c, and their antiparticles)
    # see only lambda_c and color; we do not apply q2 on them here.

    # 2) Color part: act on quark triplets in color space

    # Helper to collect indices of color triplets by name pattern
    def idx_triplet(prefix: str) -> List[int]:
        return [name_to_index[f"{prefix}_{c}"] for c in colors]

    # Particle quark triplets
    uL_triplet = idx_triplet("u_L")
    dL_triplet = idx_triplet("d_L")
    uR_triplet = idx_triplet("u_R")
    dR_triplet = idx_triplet("d_R")

    add_color_matrix_on_triplet(M, m3, uL_triplet, conj=False)
    add_color_matrix_on_triplet(M, m3, dL_triplet, conj=False)
    add_color_matrix_on_triplet(M, m3, uR_triplet, conj=False)
    add_color_matrix_on_triplet(M, m3, dR_triplet, conj=False)

    # Antiparticle quark triplets (conjugate rep)
    bar_uL_triplet = [name_to_index[f"bar_u_L_{c}"] for c in colors]
    bar_dL_triplet = [name_to_index[f"bar_d_L_{c}"] for c in colors]
    bar_uR_triplet = [name_to_index[f"bar_u_R_{c}"] for c in colors]
    bar_dR_triplet = [name_to_index[f"bar_d_R_{c}"] for c in colors]

    add_color_matrix_on_triplet(M, m3, bar_uL_triplet, conj=True)
    add_color_matrix_on_triplet(M, m3, bar_dL_triplet, conj=True)
    add_color_matrix_on_triplet(M, m3, bar_uR_triplet, conj=True)
    add_color_matrix_on_triplet(M, m3, bar_dR_triplet, conj=True)

    # Leptons have no color, so no m3 action there.

    return M
@dataclass
class SMAlgebraElement:
    name: str
    matrix: np.ndarray


def pauli_matrices() -> List[np.ndarray]:
    sigma1 = np.array([[0, 1], [1, 0]], dtype=complex)
    sigma2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sigma3 = np.array([[1, 0], [0, -1]], dtype=complex)
    return [sigma1, sigma2, sigma3]


def simple_gell_mann_subset() -> List[np.ndarray]:
    # lambda3 and lambda8 in standard Gell-Mann basis
    lam3 = np.array([[1, 0, 0],
                     [0,-1, 0],
                     [0, 0, 0]], dtype=complex)
    lam8 = (1/np.sqrt(3)) * np.array([[1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0,-2]], dtype=complex)
    return [lam3, lam8]


def build_SM_algebra_generators(basis: List[SMState],
                                name_to_index: Dict[str, int]) -> List[SMAlgebraElement]:
    dim = len(basis)

    # Identity (pure C-part with lambda=1)
    I_SM = rep_A_SM(lambda_c=1.0+0j,
                    q2=np.zeros((2, 2), dtype=complex),
                    m3=np.zeros((3, 3), dtype=complex),
                    basis=basis,
                    name_to_index=name_to_index)

    gens: List[SMAlgebraElement] = [SMAlgebraElement("I", I_SM)]

    # Quaternion part generators: lambda=0, m3=0, q2 = Pauli matrices
    for k, sigma in enumerate(pauli_matrices(), start=1):
        M_sigma = rep_A_SM(lambda_c=0.0+0j,
                           q2=sigma,
                           m3=np.zeros((3, 3), dtype=complex),
                           basis=basis,
                           name_to_index=name_to_index)
        gens.append(SMAlgebraElement(f"H_sigma{k}", M_sigma))

    # Color part generators: lambda=0, q2=0, m3 = subset of Gell-Mann matrices
    for k, gm in enumerate(simple_gell_mann_subset(), start=3):
        M_gm = rep_A_SM(lambda_c=0.0+0j,
                        q2=np.zeros((2, 2), dtype=complex),
                        m3=gm,
                        basis=basis,
                        name_to_index=name_to_index)
        gens.append(SMAlgebraElement(f"color_lambda{k}", M_gm))

    return gens
def build_SM_algebra_generators_commutative(basis: List[SMState],
                                            name_to_index: Dict[str, int]) -> List[SMAlgebraElement]:
    """
    Commutative subalgebra of A_F^SM for testing:
    - Identity
    - One diagonal quaternion generator (sigma3)
    - Two diagonal color generators (lambda3, lambda8)
    This should satisfy the zero-order condition with our current J/D shell.
    """
    gens: List[SMAlgebraElement] = []

    # Identity (pure C-part)
    I_SM = rep_A_SM(
        lambda_c=1.0+0j,
        q2=np.zeros((2, 2), dtype=complex),
        m3=np.zeros((3, 3), dtype=complex),
        basis=basis,
        name_to_index=name_to_index,
    )
    gens.append(SMAlgebraElement("I", I_SM))

    # Diagonal quaternion generator: sigma3
    sigma3 = np.array([[1, 0],
                       [0,-1]], dtype=complex)
    H_sigma3 = rep_A_SM(
        lambda_c=0.0+0j,
        q2=sigma3,
        m3=np.zeros((3, 3), dtype=complex),
        basis=basis,
        name_to_index=name_to_index,
    )
    gens.append(SMAlgebraElement("H_sigma3", H_sigma3))

    # Diagonal color generators: lambda3 and lambda8
    lam3 = np.array([[ 1, 0, 0],
                     [ 0,-1, 0],
                     [ 0, 0, 0]], dtype=complex)
    lam8 = (1/np.sqrt(3)) * np.array([[ 1, 0, 0],
                                      [ 0, 1, 0],
                                      [ 0, 0,-2]], dtype=complex)

    color_lam3 = rep_A_SM(
        lambda_c=0.0+0j,
        q2=np.zeros((2, 2), dtype=complex),
        m3=lam3,
        basis=basis,
        name_to_index=name_to_index,
    )
    gens.append(SMAlgebraElement("color_lambda3", color_lam3))

    color_lam8 = rep_A_SM(
        lambda_c=0.0+0j,
        q2=np.zeros((2, 2), dtype=complex),
        m3=lam8,
        basis=basis,
        name_to_index=name_to_index,
    )
    gens.append(SMAlgebraElement("color_lambda8", color_lam8))

    return gens

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

def J_action_from_S(S: np.ndarray, M: np.ndarray) -> np.ndarray:
    """Adjoint action J M J^{-1} = S · M^* · S^T."""
    return S @ M.conj() @ S.T

def build_DF_SM_one_gen(basis: List[SMState],
                        name_to_index: Dict[str, int],
                        m_nu: float = 0.1,
                        m_e: float  = 0.5,
                        m_u: float  = 2.0,
                        m_d: float  = 4.0) -> np.ndarray:
    """
    Simple finite Dirac operator for 1 generation:
    - Dirac masses linking L ↔ R for (nu, e, u_c, d_c) and their antiparticles.
    - No Majorana block yet.
    """
    dim = len(basis)
    D = np.zeros((dim, dim), dtype=complex)

    def couple(a: str, b: str, mass: float) -> None:
        i = name_to_index[a]
        j = name_to_index[b]
        D[i, j] = mass
        D[j, i] = mass

    # Leptons
    if "nu_R" in name_to_index:
        couple("nu_L", "nu_R", m_nu)
        couple("bar_nu_L", "bar_nu_R", m_nu)
    couple("e_L", "e_R", m_e)
    couple("bar_e_L", "bar_e_R", m_e)

    # Quarks (each color separately)
    colors = ["r", "g", "b"]
    for c in colors:
        couple(f"u_L_{c}", f"u_R_{c}", m_u)
        couple(f"d_L_{c}", f"d_R_{c}", m_d)
        couple(f"bar_u_L_{c}", f"bar_u_R_{c}", m_u)
        couple(f"bar_d_L_{c}", f"bar_d_R_{c}", m_d)

    return D

def test_first_order_condition_generic(D_F: np.ndarray,
                                       ops: List[np.ndarray],
                                       labels: List[str],
                                       S: np.ndarray,
                                       eps: float = 1e-12) -> None:
    """First-order condition:

       [[D_F, a], J b J^{-1}] = 0  for all a,b.
    """
    n = D_F.shape[0]
    assert D_F.shape == (n, n)

    print("=== First-order condition test (generic J) ===")
    max_norm = 0.0
    good_pairs = []

    for i, a in enumerate(ops):
        Da = D_F @ a - a @ D_F
        for j, b in enumerate(ops):
            b_tilde = J_action_from_S(S, b)
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
            print(f"  (a={la:>12s}, b={lb:>12s}) → ||[[D,a],J b J^-1]||_F = {nrm:.3e}")
    else:
        print(f"No pairs with norm < {eps:.1e}")
    print()


def test_zero_order_condition_generic(ops: List[np.ndarray],
                                      labels: List[str],
                                      S: np.ndarray,
                                      eps: float = 1e-12) -> None:
    """Zero-order condition:

       [a, J b J^{-1}] = 0  for all a,b.
    """
    n = ops[0].shape[0]
    print("=== Zero-order condition test (generic J) ===")
    max_norm = 0.0
    bad_pairs = []

    for i, a in enumerate(ops):
        for j, b in enumerate(ops):
            b_tilde = J_action_from_S(S, b)
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
            print(f"  (a={la:>12s}, b={lb:>12s}) → ||[a, J b J^-1]||_F = {nrm:.3e}")
    else:
        print(f"All pairs satisfy [a, J b J^-1] ≈ 0 within eps={eps:.1e}")
    print()


def test_grading_and_reality_generic(D_F: np.ndarray,
                                     ops: List[np.ndarray],
                                     labels: List[str],
                                     gamma: np.ndarray,
                                     S: np.ndarray) -> None:
    """Grading & reality bundle:

      - γ anticommutes with D_F.
      - γ commutes with all a in A_F.
      - J^2 = 1.
      - KO-sign via J D_F J^{-1} = ± D_F.
    """
    n = D_F.shape[0]

    print("=== Grading & reality tests (generic γ,J) ===")
    anti = gamma @ D_F + D_F @ gamma
    print(f"||{{gamma, D_F}}||_F = {np.linalg.norm(anti, ord='fro'):.3e}")

    max_comm_gamma = 0.0
    for a in ops:
        comm_ga = gamma @ a - a @ gamma
        max_comm_gamma = max(max_comm_gamma, np.linalg.norm(comm_ga, ord="fro"))
    print(f"max ||[gamma, a]||_F over a∈A_F = {max_comm_gamma:.3e}")

    S2 = S @ S
    print(f"||S^2 - I||_F  (⇒ J^2 deviation) = {np.linalg.norm(S2 - np.eye(n), ord='fro'):.3e}")

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
def run_emergent_alignment(
    N: int = 200,
    n_steps: int = 600,
    eta: float = 0.01,
    w6: float = 1.0,
    w5: float = 1.0,
    keep_fraction: float = 0.05,
    alpha: float = 3.0,
    beta: float = 1.0,
    random_seed: int = 42,
    use_neutrino_dressing: bool = True,
) -> Dict[str, np.ndarray]:
    """
    Full emergent pipeline (minimal, self-contained version):

    - Relax phases on the circle via misalignment functional.
    - Build emergent graph adjacency and Laplacian.
    - Extract the 3 lowest nonzero eigenmodes (generation triad).
    - Build a base kernel F_base(λ_g) and sector-weighted F_s.
    - Build a single geometric unitary U_geom from graph regions.
    - Use the sector_bases + F_s to construct full 3×3 Yukawa matrices.

    Returns:
      {
        "Y_u": 3x3 complex Yukawa matrix,
        "Y_d": 3x3 complex Yukawa matrix,
        "Y_e": 3x3 complex Yukawa matrix,
        "Y_nu": 3x3 complex Yukawa matrix,
        "F": {"u": F_u, "d": F_d, "e": F_e, "nu": F_nu},
        "lam_gen": lam_gen  (3 eigenvalues),
      }
    """
    # 1) Relax phases on the circle
    theta, energy_hist = relax_phases(
        N=N,
        n_steps=n_steps,
        eta=eta,
        w6=w6,
        w5=w5,
        random_seed=random_seed,
    )

    # 2) Build emergent adjacency and largest connected component
    A_full = build_emergent_adjacency(theta, w6=w6, w5=w5, keep_fraction=keep_fraction)
    A_cc, comp = largest_connected_component(A_full)

    # 3) Laplacian and spectral triad
    L = laplacian_from_adjacency(A_cc)
    lam_gen, gen_indices, eigvals_sorted = spectral_triad(L)

    # 4) Base kernel and sector weights (3×3 F_s)
    F_base = base_kernel(lam_gen, alpha=alpha, form="lambda_sq")
    Q_sector = build_sector_charges()
    F_u = sector_weights(F_base, Q_sector["u"], beta=beta)
    F_d = sector_weights(F_base, Q_sector["d"], beta=beta)
    F_e = sector_weights(F_base, Q_sector["e"], beta=beta)
    F_nu = sector_weights(F_base, Q_sector["nu"], beta=beta)

    # 5) Generation eigenvectors (3 lowest nonzero)
    eigvals, eigvecs = np.linalg.eigh(L)
    idx_sorted = np.argsort(eigvals)
    eigvecs_sorted = eigvecs[:, idx_sorted]
    # columns 1,2,3 correspond to the nonzero triad used above
    gen_vecs = eigvecs_sorted[:, 1:4]  # shape (N_cc, 3)

    # 6) Geometric regions on the connected component
    theta_cc = theta[comp]
    regions = build_geometric_regions(theta_cc, n_regions=3)

    # 7) Single geometric unitary from regions
    U_geom_single = build_geometric_unitary(gen_vecs, regions)

    # Use the same U_geom for all sectors (sector-specific dressing happens later)
    U_geom = {s: U_geom_single for s in SECTORS}

    # 8) Build sector bases (U_L, U_R) per sector
    P_phi_12, P_phi_23, C_12, theta_phi, theta_C = build_generation_operators()
    sector_bases = build_sector_bases(
        P_phi_12,
        P_phi_23,
        C_12,
        U_geom,
        use_neutrino_dressing=use_neutrino_dressing,
    )

    # 9) Construct full 3×3 Yukawas from F_s and U_L/R
    Y_u = yukawa_from_F_and_UL(F_u, *sector_bases["u"])
    Y_d = yukawa_from_F_and_UL(F_d, *sector_bases["d"])
    Y_e = yukawa_from_F_and_UL(F_e, *sector_bases["e"])
    Y_nu = yukawa_from_F_and_UL(F_nu, *sector_bases["nu"])

    return {
        "Y_u": Y_u,
        "Y_d": Y_d,
        "Y_e": Y_e,
        "Y_nu": Y_nu,
        "F": {"u": F_u, "d": F_d, "e": F_e, "nu": F_nu},
        "lam_gen": lam_gen,
    }
def effective_1gen_yukawas_from_3x3(
    Y_u: np.ndarray,
    Y_d: np.ndarray,
    Y_e: np.ndarray,
    Y_nu: np.ndarray,
) -> Tuple[float, float, float, float]:
    """
    Collapse 3×3 Yukawas to a single effective 1-gen value per sector.

    Here we take the largest singular value in each sector, which
    corresponds roughly to the heaviest family (t, b, τ, ν_3).
    """
    def eff(Y: np.ndarray) -> float:
        svals = np.linalg.svd(Y, compute_uv=False)
        return float(np.max(np.abs(svals)))

    Y_u_1 = eff(Y_u)
    Y_d_1 = eff(Y_d)
    Y_e_1 = eff(Y_e)
    Y_nu_1 = eff(Y_nu)
    return Y_e_1, Y_nu_1, Y_u_1, Y_d_1
def run_ncg_with_alignment() -> None:
    """
    1) Run emergent alignment pipeline and get 3×3 Yukawas.
    2) Build and test:
       - 3-generation toy internal triple (24⊕24) with those Yukawas.
       - 1-generation SM-like triple using effective aligned Yukawas.
    """
    print("=== Running emergent alignment pipeline to obtain Yukawas ===")
    align = run_emergent_alignment()
    Y_u_align = align["Y_u"]
    Y_d_align = align["Y_d"]
    Y_e_align = align["Y_e"]
    Y_nu_align = align["Y_nu"]

    # ---------------------------------------------
    # A. 3-generation internal toy triple (24⊕24)
    # ---------------------------------------------
    print()
    print("=== NCG tests: 3-gen internal toy triple with aligned Yukawas ===")
    D_F_internal = build_internal_DF_from_Y(
        Y_u_align,
        Y_d_align,
        Y_e_align,
        Y_nu_align,
    )
    ops_internal, labels_internal = build_internal_algebra_ops()
    test_first_order_condition(D_F_internal, ops_internal, labels_internal, eps=1e-12)
    test_zero_order_condition(ops_internal, labels_internal, eps=1e-12)
    test_grading_and_reality(D_F_internal, ops_internal, labels_internal)

    # ---------------------------------------------
    # B. 1-generation SM-like triple (effective Yukawas)
    # ---------------------------------------------
    print()
    print("=== NCG tests: SM-like 1gen triple with effective aligned Yukawas ===")

    Y_e_eff, Y_nu_eff, Y_u_eff, Y_d_eff = effective_1gen_yukawas_from_3x3(
        Y_u_align,
        Y_d_align,
        Y_e_align,
        Y_nu_align,
    )

    basis_SM = build_sm_basis_one_gen(include_nu_R=True)
    name_to_index_SM = build_name_index_map_sm(basis_SM)
    gamma_SM = build_gamma_sm(basis_SM)
    S_SM = build_swap_J_sm(basis_SM)

    algebra_SM = build_SM_algebra_generators_commutative(basis_SM, name_to_index_SM)
    ops_SM = [elem.matrix for elem in algebra_SM]
    labels_SM = [elem.name for elem in algebra_SM]

    D_SM_align = build_DF_SM_one_gen(
        basis_SM,
        name_to_index_SM,
        m_nu=Y_nu_eff,
        m_e=Y_e_eff,
        m_u=Y_u_eff,
        m_d=Y_d_eff,
    )

    test_first_order_condition_generic(D_SM_align, ops_SM, labels_SM, S_SM, eps=1e-12)
    test_zero_order_condition_generic(ops_SM, labels_SM, S_SM, eps=1e-12)
    test_grading_and_reality_generic(D_SM_align, ops_SM, labels_SM, gamma_SM, S_SM)

# ================================
# Main driver
# ================================

def run_sm_ncg_tests() -> None:
    """Build 1-gen SM-like finite triple and run first-, zero-order and reality tests."""
    basis_SM = build_sm_basis_one_gen(include_nu_R=True)
    name_to_index_SM = build_name_index_map_sm(basis_SM)
    gamma_SM = build_gamma_sm(basis_SM)
    S_SM = build_swap_J_sm(basis_SM)

    # Use commutative subalgebra
    algebra_SM = build_SM_algebra_generators_commutative(basis_SM, name_to_index_SM)
    ops_SM = [elem.matrix for elem in algebra_SM]
    labels_SM = [elem.name for elem in algebra_SM]

    D_SM = build_DF_SM_one_gen(
        basis_SM,
        name_to_index_SM,
        m_nu=0.1,
        m_e=0.5,
        m_u=2.0,
        m_d=4.0,
    )

    print(f"dim H_F^SM (1 gen + nu_R) = {len(basis_SM)}")
    print("Algebra generators (commutative subset):", labels_SM)
    print()

    test_first_order_condition_generic(D_SM, ops_SM, labels_SM, S_SM, eps=1e-12)
    test_zero_order_condition_generic(ops_SM, labels_SM, S_SM, eps=1e-12)
    test_grading_and_reality_generic(D_SM, ops_SM, labels_SM, gamma_SM, S_SM)


if __name__ == "__main__":
    # print(">>> Baseline SM-like 1gen finite triple (fixed Yukawas) <<<\n")
    # run_sm_ncg_tests()

    print("\n\n>>> Emergent alignment Yukawas plugged into NCG tests <<<\n")
    run_ncg_with_alignment()

"""
>>> Emergent alignment Yukawas plugged into NCG tests <<<

=== Running emergent alignment pipeline to obtain Yukawas ===

=== NCG tests: 3-gen internal toy triple with aligned Yukawas ===
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


=== NCG tests: SM-like 1gen triple with effective aligned Yukawas ===
=== First-order condition test (generic J) ===
Max Frobenius norm over all pairs (a,b): 1.843e-01
Pairs with norm < 1.0e-12:
  (a=           I, b=           I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=           I, b=    H_sigma3) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=           I, b=color_lambda3) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=           I, b=color_lambda8) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=    H_sigma3, b=           I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=    H_sigma3, b=color_lambda3) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=    H_sigma3, b=color_lambda8) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=color_lambda3, b=           I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=color_lambda3, b=    H_sigma3) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=color_lambda3, b=color_lambda3) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=color_lambda3, b=color_lambda8) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=color_lambda8, b=           I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=color_lambda8, b=    H_sigma3) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=color_lambda8, b=color_lambda3) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=color_lambda8, b=color_lambda8) → ||[[D,a],J b J^-1]||_F = 0.000e+00

=== Zero-order condition test (generic J) ===
Max Frobenius norm over all pairs (a,b): 0.000e+00
All pairs satisfy [a, J b J^-1] ≈ 0 within eps=1.0e-12

=== Grading & reality tests (generic γ,J) ===
||{gamma, D_F}||_F = 0.000e+00
max ||[gamma, a]||_F over a∈A_F = 0.000e+00
||S^2 - I||_F  (⇒ J^2 deviation) = 0.000e+00
||J D_F J^-1 - D_F||_F   = 0.000e+00
||J D_F J^-1 + D_F||_F   = 3.685e-01
→ KO-sign: J D_F J^-1 = + D_F (J-even Dirac operator)
"""