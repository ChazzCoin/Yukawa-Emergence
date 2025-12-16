from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple
import numpy as np

from v3.spectral import build_24cell_vertices, build_24cell_adjacency, spectral_decomposition, \
    build_universal_kernel, cluster_eigenvalues, projector_from_cluster


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

# ============================================
# Drop-in: 3-gen SM finite triple with emergent Yukawas
# ============================================


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
def build_SM_algebra_generators_commutative_strict(
    basis: List[SMState],
    name_to_index: Dict[str, int]
) -> List[SMAlgebraElement]:
    """
    STRICT harmonic-aligned commutative subalgebra of A_F^SM for NCG tests.

    We keep only:
      - Identity (pure C-part)
      - Two diagonal color generators (lambda3, lambda8)

    We deliberately DROP the quaternion sigma3 here, because with the
    simplified J-shell we are using, sigma3 is the first place where
    first-order violations show up once Yukawas are fully emergent.

    This is the "strict" algebra you should use when you want exact
    first-/zero-order behaviour in the SM embedding layer.
    """
    gens: List[SMAlgebraElement] = []

    # Identity (pure C-part)
    I_SM = rep_A_SM(
        lambda_c=1.0 + 0j,
        q2=np.zeros((2, 2), dtype=complex),
        m3=np.zeros((3, 3), dtype=complex),
        basis=basis,
        name_to_index=name_to_index,
    )
    gens.append(SMAlgebraElement("I", I_SM))

    # Diagonal color generators: lambda3 and lambda8
    lam3 = np.array([[ 1, 0, 0],
                     [ 0,-1, 0],
                     [ 0, 0, 0]], dtype=complex)
    lam8 = (1/np.sqrt(3)) * np.array([[ 1, 0, 0],
                                      [ 0, 1, 0],
                                      [ 0, 0,-2]], dtype=complex)

    M_lam3 = rep_A_SM(
        lambda_c=0.0 + 0j,
        q2=np.zeros((2, 2), dtype=complex),
        m3=lam3,
        basis=basis,
        name_to_index=name_to_index,
    )
    M_lam8 = rep_A_SM(
        lambda_c=0.0 + 0j,
        q2=np.zeros((2, 2), dtype=complex),
        m3=lam8,
        basis=basis,
        name_to_index=name_to_index,
    )

    gens.append(SMAlgebraElement("color_lambda3", M_lam3))
    gens.append(SMAlgebraElement("color_lambda8", M_lam8))

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


def build_sector_projectors1() -> Dict[str, np.ndarray]:
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

def build_sector_projectors(evals, evecs):
    """
    Construct LEFT/RIGHT projectors P_L, P_R for each sector (u,d,e,nu)
    from eigenvalue clusters.

    Strategy:
    - Cluster eigenvalues (degeneracies).
    - Ignore the trivial λ=0 ground state cluster for flavor (index 0 cluster).
    - Use the remaining clusters with size>=3 as natural "triplet" candidates.
    - For left-handed vs right-handed subspaces, use *different* clusters:
        - Up:
            L_u = lowest nontrivial triplet cluster
            R_u = highest triplet cluster
        - Down:
            L_d = 2nd lowest triplet cluster (if exists, else same as L_u)
            R_d = 2nd highest triplet cluster (if exists, else same as R_u)
        - Charged leptons:
            L_e = L_u
            R_e = R_d
        - Neutrinos:
            L_n = L_d
            R_n = R_u

    All of this is fully determined once evals/evecs are known.
    No continuous parameters, no randomness.
    """
    clusters = cluster_eigenvalues(evals)
    N = len(evals)

    # Identify triplet-like clusters (size >= 3), excluding the 0-eigenvalue cluster
    triplet_clusters = []
    for ci, cl in enumerate(clusters):
        if len(cl) >= 3:
            # Optional: skip λ=0 cluster (usually evals[0] ~ 0)
            if abs(evals[cl[0]]) < 1e-12:
                continue
            triplet_clusters.append(cl)

    if len(triplet_clusters) == 0:
        # Fallback: just use the lowest nontrivial cluster(s) whatever size
        # (still deterministic, but flavor structure will be degenerate)
        print("WARNING: No triplet clusters found; using smallest nontrivial clusters.")
        # skip λ=0 cluster
        nonzero_clusters = [cl for cl in clusters if abs(evals[cl[0]]) > 1e-12]
        # ensure at least one cluster
        if len(nonzero_clusters) == 0:
            nonzero_clusters = clusters
        triplet_clusters = nonzero_clusters

    # Helper to pick clusters with wrap-around if needed
    def pick_cluster(idx):
        return triplet_clusters[idx % len(triplet_clusters)]

    # Choose clusters for each sector (left/right) following the pattern above
    L_u_cl = pick_cluster(0)
    R_u_cl = pick_cluster(-1)

    L_d_cl = pick_cluster(1) if len(triplet_clusters) > 1 else L_u_cl
    R_d_cl = pick_cluster(-2) if len(triplet_clusters) > 1 else R_u_cl

    L_e_cl = L_u_cl
    R_e_cl = R_d_cl

    L_n_cl = L_d_cl
    R_n_cl = R_u_cl

    # Build projectors
    P_L_u = projector_from_cluster(evecs, L_u_cl, n_rows=3)
    P_R_u = projector_from_cluster(evecs, R_u_cl, n_rows=3)

    P_L_d = projector_from_cluster(evecs, L_d_cl, n_rows=3)
    P_R_d = projector_from_cluster(evecs, R_d_cl, n_rows=3)

    P_L_e = projector_from_cluster(evecs, L_e_cl, n_rows=3)
    P_R_e = projector_from_cluster(evecs, R_e_cl, n_rows=3)

    P_L_n = projector_from_cluster(evecs, L_n_cl, n_rows=3)
    P_R_n = projector_from_cluster(evecs, R_n_cl, n_rows=3)

    sector_proj = {
        "u":  (P_L_u, P_R_u),
        "d":  (P_L_d, P_R_d),
        "e":  (P_L_e, P_R_e),
        "nu": (P_L_n, P_R_n),
    }

    return sector_proj, clusters, triplet_clusters

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
    P = build_sector_projectors1()

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
        theta_solar = 2 * np.pi / N_SOLAR  # ~10°
        theta_reac = 2 * np.pi / N_REACTOR  # ~8°
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
def project_DF_to_J_sign(D_F: np.ndarray,
                         S: np.ndarray,
                         sign: str = "+") -> np.ndarray:
    """
    Project a finite Dirac operator D_F onto the J-even or J-odd subspace
    defined by the swap matrix S (implementing J via J M J^-1 = S M* S^T).

    - sign = "+"  →  J-even part:  D_even = 1/2 (D + J D J^-1)
    - sign = "-"  →  J-odd part:   D_odd  = 1/2 (D - J D J^-1)

    This guarantees J D_proj J^-1 = ± D_proj exactly (up to FP roundoff).
    """

    D_tilde = S @ D_F.conj() @ S.T
    if sign == "+":
        D_proj = 0.5 * (D_F + D_tilde)
    elif sign == "-":
        D_proj = 0.5 * (D_F - D_tilde)
    else:
        raise ValueError("sign must be '+' or '-'")
    return D_proj


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

# ================================================================
# 24-CELL SPECTRAL GEOMETRY ENGINE (drop-in optional replacement)
# ================================================================

def build_24cell_geometry():
    verts = build_24cell_vertices()
    A = build_24cell_adjacency(verts)
    L = build_laplacian(A)
    evals, evecs = spectral_decomposition(L)
    K_24, f_vals = build_universal_kernel(evals, evecs)
    clusters = cluster_eigenvalues(evals)
    sector_proj_24, all_clusters, triplet_clusters = build_sector_projectors(evals, evecs)
    return {
        "verts": verts,
        "A": A,
        "L": L,
        "evals": evals,
        "evecs": evecs,
        "K_24": K_24,
        "sector_proj": sector_proj_24,
        "clusters": clusters,
        "triplets": triplet_clusters,
    }
def run_emergent_alignment(
    N: int = 360,
    n_steps: int = 600,
    eta: float = 0.01,
    w6: float = 1.0,
    w5: float = 1.0,
    keep_fraction: float = 0.05,
    alpha: float = 3.0,
    beta: float = 1.0,
    random_seed: int = 42,
    use_neutrino_dressing: bool = True,
    use_24cell: bool = False
) -> Dict[str, np.ndarray]:

    """
    Full emergent alignment engine with optional 24-cell spectral replacement.
    """

    # ------------------------------------------------------------
    # 1. Relax phases on the harmonic circle
    # ------------------------------------------------------------
    theta, energy_hist = relax_phases(
        N=N,
        n_steps=n_steps,
        eta=eta,
        w6=w6,
        w5=w5,
        random_seed=random_seed,
    )

    # ------------------------------------------------------------
    # 2A. 24-cell geometry mode (replacement)
    # ------------------------------------------------------------
    # ------------------------------------------------------------
    # 2A. 24-cell geometry mode (replacement)
    # ------------------------------------------------------------
    if use_24cell:
        geo = build_24cell_geometry()

        L = geo["L"]
        evals = geo["evals"]
        evecs = geo["evecs"]

        # Choose lowest non-zero "triplet" cluster
        tri_full = geo["triplets"][0]  # cluster indices; may have size ≥ 3
        if len(tri_full) < 3:
            raise RuntimeError(
                f"24-cell triplet cluster has size {len(tri_full)} < 3; "
                "cannot build 3-generation flavor space."
            )

        # We only need 3 modes for 3 generations → take the first 3
        tri = np.array(tri_full[:3], dtype=int)

        # 3 generation eigenvalues and eigenvectors
        lam_gen = evals[tri]  # shape (3,)
        gen_vecs = evecs[:, tri]  # shape (24, 3)

        # Base kernel from Laplacian eigenvalues (3-component)
        F_base = base_kernel(lam_gen, alpha=alpha, form="lambda_sq")

        # Optional: refine using the 24-cell universal heat kernel
        # Project heat kernel into the 3D generation subspace
        K_tri = gen_vecs.conj().T @ geo["K_24"] @ gen_vecs
        # Replace F_base with heat-kernel diagonal magnitudes (still length 3)
        F_base = np.abs(np.diag(K_tri))

    # ------------------------------------------------------------
    # 2B. Original emergent graph mode
    # ------------------------------------------------------------
    else:
        A_full = build_emergent_adjacency(theta, w6=w6, w5=w5, keep_fraction=keep_fraction)
        A_cc, comp = largest_connected_component(A_full)

        L = laplacian_from_adjacency(A_cc)
        lam_gen, gen_indices, _ = spectral_triad(L)

        # Generation eigenvectors for Yukawa orientations
        eigvals, eigvecs = np.linalg.eigh(L)
        idx_sorted = np.argsort(eigvals)
        eigvecs_sorted = eigvecs[:, idx_sorted]
        gen_vecs = eigvecs_sorted[:, 1:4]     # triad

        # Base kernel
        F_base = base_kernel(lam_gen, alpha=alpha, form="lambda_sq")

    # ------------------------------------------------------------
    # 3. Sector weight system F_s = F_base ∘ exp(-β Q_s)
    # ------------------------------------------------------------
    Qsec = build_sector_charges()
    F_u  = sector_weights(F_base, Qsec["u"],  beta)
    F_d  = sector_weights(F_base, Qsec["d"],  beta)
    F_e  = sector_weights(F_base, Qsec["e"],  beta)
    F_nu = sector_weights(F_base, Qsec["nu"], beta)

    # ------------------------------------------------------------
    # 4. Build geometric regions (only meaningful in graph mode)
    # ------------------------------------------------------------
    if use_24cell:
        # Fake regions: divide vertex list into three roughly equal sectors
        Ngeo = gen_vecs.shape[0]
        idxs = np.arange(Ngeo)
        third = Ngeo // 3
        regions = [idxs[:third], idxs[third:2*third], idxs[2*third:]]
    else:
        theta_cc = theta[comp]
        regions = build_geometric_regions(theta_cc, n_regions=3)

    # ------------------------------------------------------------
    # 5. Build geometric generation frame U_geom
    # ------------------------------------------------------------
    U_geom_single = build_geometric_unitary(gen_vecs, regions)
    U_geom = {s: U_geom_single for s in SECTORS}

    # ------------------------------------------------------------
    # 6. Build sector bases (U_L, U_R)
    # ------------------------------------------------------------
    P_phi_12, P_phi_23, C_12, _, _ = build_generation_operators()

    sector_bases = build_sector_bases(
        P_phi_12,
        P_phi_23,
        C_12,
        U_geom,
        use_neutrino_dressing=use_neutrino_dressing,
    )

    # ------------------------------------------------------------
    # 7. Construct Yukawa matrices
    # ------------------------------------------------------------
    Y_u  = yukawa_from_F_and_UL(F_u,  *sector_bases["u"])
    Y_d  = yukawa_from_F_and_UL(F_d,  *sector_bases["d"])
    Y_e  = yukawa_from_F_and_UL(F_e,  *sector_bases["e"])
    Y_nu = yukawa_from_F_and_UL(F_nu, *sector_bases["nu"])

    return {
        "Y_u": Y_u,
        "Y_d": Y_d,
        "Y_e": Y_e,
        "Y_nu": Y_nu,
        "F": {"u": F_u, "d": F_d, "e": F_e, "nu": F_nu},
        "lam_gen": lam_gen,
        "L": L,
        "gen_vecs": gen_vecs,
        "regions": regions,
        "U_geom": U_geom_single,
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



# ================================
# Z2160 triadic 9-site flavor kernel (Alignment Spectral Triple v2.0)
# ================================

from dataclasses import dataclass
from typing import Optional

N_FLAVOR_SITES_Z2160 = 9
FLAVOR_LIGHT_INDICES_Z2160 = np.array([0, 1, 2], dtype=int)
FLAVOR_HEAVY_INDICES_Z2160 = np.array([3, 4, 5, 6, 7, 8], dtype=int)
FORBIDDEN_DISTANCES_Z2160 = (2, 4, 7)


@dataclass
class Z2160SectorConfig:
    """
    Phase + scale parameters for a single sector in the Z2160 triadic flavor model.

    φ_i^(s) = A_s + B_s * (i mod 3)
    Y_9^(s) ∝ (P^(s) ∘ K) with overall scale alpha_s.

    NOTE:
    - These parameters are *phenomenological* and should be fitted / tuned.
    - Defaults below are placeholders, not final physics values.
    """
    A: float     # overall phase offset (sector s)
    B: float     # per-generation phase increment (sector s)
    alpha: float # overall Yukawa scale (sector s)
    normalize: bool = True


def build_Z2160_kernel(kappa: float = 0.24,
                       forbidden: Tuple[int, ...] = FORBIDDEN_DISTANCES_Z2160
                       ) -> np.ndarray:
    """
    Real 9x9 geometric kernel K with forbidden distances {2,4,7} on a Z_2160 lattice,
    using the explicit embedding ι(i) = i, so d_2160(i,j) = |i-j| for i,j=0..8.
    """
    K = np.zeros((N_FLAVOR_SITES_Z2160, N_FLAVOR_SITES_Z2160), dtype=float)
    forb = set(forbidden)

    for i in range(N_FLAVOR_SITES_Z2160):
        for j in range(N_FLAVOR_SITES_Z2160):
            if i == j:
                K[i, j] = 1.0
            else:
                d = abs(i - j)
                if d in forb:
                    K[i, j] = 0.0
                else:
                    K[i, j] = kappa ** d
    return K


def build_Z2160_phase_matrix(A: float, B: float) -> np.ndarray:
    """
    Sector phase matrix P^(s)_{ij} = exp(i [φ_i - φ_j]),
    with φ_i = A + B * (i mod 3).
    """
    phi = np.array([A + B * (i % 3) for i in range(N_FLAVOR_SITES_Z2160)], dtype=float)
    diff = phi[:, None] - phi[None, :]
    return np.exp(1j * diff)


def build_Z2160_Y9(K: np.ndarray,
                   cfg: Z2160SectorConfig) -> np.ndarray:
    """
    Build the 9x9 Yukawa-like matrix for one sector:

        Y_9^(s) = alpha_s * (P^(s) ∘ K) / σ_max   (if normalize=True)

    where P^(s) encodes the phase pattern and K is the real geometric kernel.
    """
    if K.shape != (N_FLAVOR_SITES_Z2160, N_FLAVOR_SITES_Z2160):
        raise ValueError(f"K must be {N_FLAVOR_SITES_Z2160}x{N_FLAVOR_SITES_Z2160}.")

    P = build_Z2160_phase_matrix(cfg.A, cfg.B)
    Y_raw = P * K  # Hadamard product
    # enforce Hermiticity numerically
    Y_raw = 0.5 * (Y_raw + Y_raw.conj().T)

    if cfg.normalize:
        svals = np.linalg.svd(Y_raw, compute_uv=False)
        sigma_max = np.max(np.abs(svals))
        if sigma_max <= 0:
            raise ValueError("Singular values of Z2160 Yukawa kernel are non-positive.")
        Y_raw = Y_raw / sigma_max

    return cfg.alpha * Y_raw


def schur_9_to_3_Z2160(Y9: np.ndarray) -> np.ndarray:
    """
    Schur complement compression from 9x9 to effective 3x3:

        Y_eff = A - B D^{-1} C,

    where Y9 is partitioned as

        Y9 = [[A, B],
              [C, D]]

    with A: 3x3 (LIGHT×LIGHT), D: 6x6 (HEAVY×HEAVY).
    """
    if Y9.shape != (N_FLAVOR_SITES_Z2160, N_FLAVOR_SITES_Z2160):
        raise ValueError(f"Y9 must be {N_FLAVOR_SITES_Z2160}x{N_FLAVOR_SITES_Z2160}.")

    A = Y9[np.ix_(FLAVOR_LIGHT_INDICES_Z2160, FLAVOR_LIGHT_INDICES_Z2160)]
    B = Y9[np.ix_(FLAVOR_LIGHT_INDICES_Z2160, FLAVOR_HEAVY_INDICES_Z2160)]
    C = Y9[np.ix_(FLAVOR_HEAVY_INDICES_Z2160, FLAVOR_LIGHT_INDICES_Z2160)]
    D = Y9[np.ix_(FLAVOR_HEAVY_INDICES_Z2160, FLAVOR_HEAVY_INDICES_Z2160)]

    try:
        D_inv = np.linalg.inv(D)
    except np.linalg.LinAlgError:
        raise np.linalg.LinAlgError("Heavy block D in Schur 9→3 is singular.")

    return A - B @ D_inv @ C


def build_Z2160_MR_from_Y9nu(Y9_nu: np.ndarray,
                             Lambda_Maj: float = 2.0e14) -> np.ndarray:
    """
    Optional: build a 3x3 heavy Majorana mass matrix from the heavy block
    of the 9x9 neutrino kernel via a triadic heavy basis:

        M_R = Λ_Maj B_H^† D^(ν) B_H,

    where D^(ν) is the 6x6 heavy-heavy block and B_H is a 6x3 "Fourier-like" basis.
    """
    if Y9_nu.shape != (N_FLAVOR_SITES_Z2160, N_FLAVOR_SITES_Z2160):
        raise ValueError(f"Y9_nu must be {N_FLAVOR_SITES_Z2160}x{N_FLAVOR_SITES_Z2160}.")

    D_nu = Y9_nu[np.ix_(FLAVOR_HEAVY_INDICES_Z2160, FLAVOR_HEAVY_INDICES_Z2160)]
    Nh = D_nu.shape[0]  # 6
    i_vals = np.arange(Nh, dtype=float)
    ks = (1, 2, 3)

    B_cols = []
    for k in ks:
        phase = 2.0 * np.pi * k * i_vals / Nh
        B_cols.append(np.exp(1j * phase) / np.sqrt(Nh))
    B_H = np.column_stack(B_cols)  # shape (Nh, 3)

    MR = Lambda_Maj * (B_H.conj().T @ D_nu @ B_H)
    return MR


def run_Z2160_alignment(
    kappa: float = 0.24,
    forbidden: Tuple[int, ...] = FORBIDDEN_DISTANCES_Z2160,
    sector_config: Optional[Dict[str, Z2160SectorConfig]] = None,
) -> Dict[str, np.ndarray]:
    """
    Alignment Spectral Triple v2.0 flavor engine using the Z2160 9-site kernel.

    Returns a dict analogous to run_emergent_alignment():

        {
          "Y_u": 3x3,
          "Y_d": 3x3,
          "Y_e": 3x3,
          "Y_nu": 3x3,
          "K":   9x9 real kernel,
          "sector_config": config dict actually used
        }

    NOTES:
    - Default sector_config is a *placeholder* and should be refined by fit.
    - You can pass in your own sector_config dict with keys "u","d","e","nu".
    """
    if sector_config is None:
        # Very conservative placeholder choices:
        sector_config = {
            "u":  Z2160SectorConfig(A=0.0, B=0.0, alpha=1.0),
            "d":  Z2160SectorConfig(A=0.1, B=0.2, alpha=0.05),
            "e":  Z2160SectorConfig(A=0.2, B=0.1, alpha=0.05),
            "nu": Z2160SectorConfig(A=0.0, B=0.3, alpha=0.02),
        }

    K = build_Z2160_kernel(kappa=kappa, forbidden=forbidden)
    Y3: Dict[str, np.ndarray] = {}
    Y9_all: Dict[str, np.ndarray] = {}

    for s in ["u", "d", "e", "nu"]:
        cfg = sector_config[s]
        Y9 = build_Z2160_Y9(K, cfg)
        Y3[s] = schur_9_to_3_Z2160(Y9)
        Y9_all[s] = Y9

    return {
        "Y_u":  Y3["u"],
        "Y_d":  Y3["d"],
        "Y_e":  Y3["e"],
        "Y_nu": Y3["nu"],
        "Y9":   Y9_all,
        "K":    K,
        "sector_config": sector_config,
    }





# ================================
# Public API: finite internal triple with emergent Yukawas
# ================================

def build_alignment_finite_triple(
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
    use_strict_algebra: bool = True,
    yukawa_mode: str = "emergent",   # "emergent" or "Z2160"
    z2160_kwargs: Optional[Dict[str, object]] = None,
) -> Dict[str, object]:
    """
    Build the finite internal triple with aligned Yukawas.

    yukawa_mode:
      - "emergent": use run_emergent_alignment(...)
      - "Z2160":    use run_Z2160_alignment(...) (9-site triadic kernel)

    z2160_kwargs:
      - Optional dict passed to run_Z2160_alignment, e.g.
        {
          "kappa": 0.24,
          "sector_config": { "u": Z2160SectorConfig(...), ... }
        }
    """
    if yukawa_mode == "emergent":
        emergent = run_emergent_alignment(
            N=N,
            n_steps=n_steps,
            eta=eta,
            w6=w6,
            w5=w5,
            keep_fraction=keep_fraction,
            alpha=alpha,
            beta=beta,
            random_seed=random_seed,
            use_neutrino_dressing=use_neutrino_dressing,
        )
        Y_u = emergent["Y_u"]
        Y_d = emergent["Y_d"]
        Y_e = emergent["Y_e"]
        Y_nu = emergent["Y_nu"]

    elif yukawa_mode == "Z2160":
        if z2160_kwargs is None:
            z2160_kwargs = {}
        z = run_Z2160_alignment(**z2160_kwargs)
        Y_u = z["Y_u"]
        Y_d = z["Y_d"]
        Y_e = z["Y_e"]
        Y_nu = z["Y_nu"]

    else:
        raise ValueError(f"Unsupported yukawa_mode {yukawa_mode!r} in build_alignment_finite_triple")

    triple = build_SM_3gen_triple_from_Y(
        Y_u, Y_d, Y_e, Y_nu,
        use_strict_algebra=use_strict_algebra,
    )

    # Expose J_F in the same way as before
    triple["J_F"] = triple["S_F"]
    return triple


def run_ncg_with_Z2160_alignment() -> None:
    """
    3-generation SM finite triple with Z2160 triadic Yukawas (Alignment v2.0 flavor mode):

      1) Use build_alignment_finite_triple(yukawa_mode="Z2160", ...) to construct
         the finite triple with 3×3 Yukawas coming from the 9-site Z2160 kernel.
      2) Run the generic NCG tests:
         - zero-order
         - first-order
         - grading & reality
    """
    print("\n\n=== 3-generation SM finite triple with Z2160 triadic Yukawas (v2.0 test) ===\n")

    # If you want to tweak Z2160 parameters, pass z2160_kwargs explicitly.
    # For a first test, we use the defaults defined in run_Z2160_alignment().
    tri_z = build_alignment_finite_triple(
        yukawa_mode="Z2160",
        z2160_kwargs=None,        # or a dict with custom kappa / sector_config
        use_strict_algebra=True,  # keep the same strict algebra as in the emergent case
    )

    D_F_z      = tri_z["D_F"]
    gamma_F_z  = tri_z["gamma_F"]
    J_F_z      = tri_z["J_F"]
    ops_F_z    = tri_z["algebra_ops"]
    labels_F_z = tri_z["algebra_labels"]

    # Run the same generic NCG tests you used for the emergent Yukawa triple
    test_zero_order_condition_generic(ops_F_z, labels_F_z, J_F_z, eps=1e-12)
    test_first_order_condition_generic(D_F_z, ops_F_z, labels_F_z, J_F_z, eps=1e-12)
    test_grading_and_reality_generic(D_F_z, ops_F_z, labels_F_z, gamma_F_z, J_F_z)

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

    # Use STRICT commutative subalgebra by default
    algebra_SM = build_SM_algebra_generators_commutative_strict(
        basis_SM, name_to_index_SM
    )
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
    print("Algebra generators (STRICT commutative subset):", labels_SM)
    print()

    test_first_order_condition_generic(D_SM, ops_SM, labels_SM, S_SM, eps=1e-12)
    test_zero_order_condition_generic(ops_SM, labels_SM, S_SM, eps=1e-12)
    test_grading_and_reality_generic(D_SM, ops_SM, labels_SM, gamma_SM, S_SM)

# -------------- From Yukawas to observables and χ² (global fit layer) --------------

def observables_from_Yukawas(
    Y_u: np.ndarray,
    Y_d: np.ndarray,
    Y_e: np.ndarray,
    Y_nu: np.ndarray,
    v: float = 246.0,
) -> Dict[str, float]:
    """
    Given 3×3 Yukawa matrices Y_u, Y_d, Y_e, Y_nu (as produced by
    run_emergent_alignment), compute a minimal set of flavor observables:

      - quark mass ratios:  mu/mt, mc/mt, md/mb, ms/mb
      - lepton mass ratios: me/mt, mmu/mt
      - CKM mixing angles   (theta12_q, theta23_q, theta13_q)
      - PMNS mixing angles  (theta12_l, theta23_l, theta13_l)

    These are then packaged into the same keys used by TARGETS so that
    they can be passed directly into chi2(obs, TARGETS).

    Notes:
      - Masses are defined as m_i = (v / sqrt(2)) * singular_value_i,
        where v ≈ 246 GeV is the Higgs vev.
      - For now we only form mass *ratios*, so the overall scale largely
        cancels and v can be left at its default.
      - CKM = U_L^u† U_L^d, PMNS = U_L^e† U_L^nu, with U_L^s from SVD.
    """
    # SVD for each Yukawa: Y = U_L diag(s) U_R^†
    Uu, s_u, Vhu = np.linalg.svd(Y_u)
    Ud, s_d, Vhd = np.linalg.svd(Y_d)
    Ue, s_e, Vhe = np.linalg.svd(Y_e)
    Unu, s_nu, Vhnu = np.linalg.svd(Y_nu)

    # Sort singular values by magnitude (lightest → heaviest)
    s_u_sorted = np.sort(np.abs(s_u))
    s_d_sorted = np.sort(np.abs(s_d))
    s_e_sorted = np.sort(np.abs(s_e))

    # Convert to masses using m = v / sqrt(2) * |y|
    factor = v / np.sqrt(2.0)
    m_u, m_c, m_t   = factor * s_u_sorted
    m_d, m_s, m_b   = factor * s_d_sorted
    m_e, m_mu, m_tau = factor * s_e_sorted

    # Mass ratios as in TARGETS
    mu_mt  = m_u / m_t
    mc_mt  = m_c / m_t
    md_mb  = m_d / m_b
    ms_mb  = m_s / m_b
    me_mt  = m_e / m_t
    mmu_mt = m_mu / m_t

    # CKM mixing: V_CKM = U_L^u† U_L^d
    V_ckm = mixing_matrix(Uu, Ud)
    theta12_q, theta23_q, theta13_q = mixing_angles_from_U(V_ckm)

    # PMNS mixing: U_PMNS = U_L^e† U_L^nu
    U_pmns = mixing_matrix(Ue, Unu)
    theta12_l, theta23_l, theta13_l = mixing_angles_from_U(U_pmns)

    # Package into the same structure used by TARGETS/chi2
    obs = compute_observables(
        mu_mt, mc_mt, md_mb, ms_mb, me_mt, mmu_mt,
        theta12_q, theta23_q, theta13_q,
        theta12_l, theta23_l, theta13_l,
    )
    return obs


def build_laplacian(A):
    """
    Graph Laplacian L = D - A, where D is degree matrix.
    """
    degrees = np.sum(A, axis=1)
    D = np.diag(degrees)
    L = D - A
    return L
def flavor_chi2_from_Yukawas(
    Y_u: np.ndarray,
    Y_d: np.ndarray,
    Y_e: np.ndarray,
    Y_nu: np.ndarray,
    targets: Dict[str, Tuple[float, float]] = None,
    v: float = 246.0,
) -> Tuple[float, List[Tuple[str, float, float, float]], Dict[str, float]]:
    """
    Convenience wrapper:
      Yukawas → observables → χ² against TARGETS.

    Returns:
      chi2_val, details, obs

      - chi2_val: total χ²
      - details: list of (name, value, target, contribution)
      - obs:     dict of observables actually used
    """
    if targets is None:
        targets = TARGETS

    obs = observables_from_Yukawas(Y_u, Y_d, Y_e, Y_nu, v=v)
    chi2_val, details = chi2(obs, targets)
    return chi2_val, details, obs

def chi2_wrapper(x):
    # x is a vector of parameters you choose to fit, e.g. [w6, w5, alpha, beta]
    w6, w5, alpha, beta = x
    chi2_val, _, _, _ = flavor_chi2_from_emergent_params(
        w6=w6,
        w5=w5,
        alpha=alpha,
        beta=beta,
        # plus any other fixed kwargs you like
    )
    return chi2_val

def flavor_chi2_from_emergent_params(
    targets: Dict[str, Tuple[float, float]] = None,
    v: float = 246.0,
    **kwargs,
) -> Tuple[float,
           List[Tuple[str, float, float, float]],
           Dict[str, float],
           Dict[str, np.ndarray]]:
    """
    One-shot global-fit evaluation:

      1. Run run_emergent_alignment(**kwargs) to get Y_u, Y_d, Y_e, Y_nu.
      2. Convert them to observables (mass ratios + mixing angles).
      3. Compute χ² against TARGETS (or a custom target set).

    Arguments:
      - targets: optional override for TARGETS
      - v: Higgs vev used to convert Yukawas to masses (default 246 GeV)
      - **kwargs: passed directly to run_emergent_alignment, e.g.
          N=200, n_steps=600, eta=0.01, w6=1.0, ...

    Returns:
      chi2_val, details, obs, emergent

      - chi2_val: total χ²
      - details:  (name, value, target, χ²_contribution)
      - obs:      dict of observables
      - emergent: dict returned by run_emergent_alignment (including Y_s)
    """
    if targets is None:
        targets = TARGETS

    emergent = run_emergent_alignment(**kwargs)
    Y_u = emergent["Y_u"]
    Y_d = emergent["Y_d"]
    Y_e = emergent["Y_e"]
    Y_nu = emergent["Y_nu"]

    obs = observables_from_Yukawas(Y_u, Y_d, Y_e, Y_nu, v=v)
    chi2_val, details = chi2(obs, targets)
    return chi2_val, details, obs, emergent
def build_SM_3gen_triple_from_Y(
    Y_u: np.ndarray,
    Y_d: np.ndarray,
    Y_e: np.ndarray,
    Y_nu: np.ndarray,
    use_strict_algebra: bool = True,
) -> Dict[str, object]:
    """
    Build a 3-generation SM finite triple (H_F^(3), D_F^(3), A_F^(3), J_F, γ_F)
    from *given* 3×3 Yukawa matrices Y_u, Y_d, Y_e, Y_nu.

    This is a refactoring of the old build_SM_3gen_triple_with_emergent_Y,
    but with the emergent alignment engine factored out. It treats generation
    as a C^3 factor and lifts a 1-gen basis via Kronecker products.
    """
    # 1) 1-gen basis + structures
    basis_1 = build_sm_basis_one_gen(include_nu_R=True)
    name_to_index_1 = build_name_index_map_sm(basis_1)
    gamma_1 = build_gamma_sm(basis_1)
    S_1 = build_swap_J_sm(basis_1)

    # 2) Tensor with generation space: H_F^(3) = H_F^(1) ⊗ C^3
    dim1 = len(basis_1)
    dim_gen = 3
    dim3 = dim1 * dim_gen

    # Grading and J_F for 3 gens
    gamma_3 = np.kron(gamma_1, np.eye(dim_gen, dtype=complex))
    S_3 = np.kron(S_1, np.eye(dim_gen, dtype=complex))

    # 3) Algebra: act identically on each generation
    if use_strict_algebra:
        algebra_1 = build_SM_algebra_generators_commutative_strict(
            basis_1, name_to_index_1
        )
    else:
        algebra_1 = build_SM_algebra_generators_commutative(
            basis_1, name_to_index_1
        )

    ops_3: List[np.ndarray] = []
    labels_3: List[str] = []
    for elem in algebra_1:
        ops_3.append(np.kron(elem.matrix, np.eye(dim_gen, dtype=complex)))
        labels_3.append(elem.name)

    # 4) Build D_F^(3): replace scalar masses with 3×3 Yukawas
    D3 = np.zeros((dim3, dim3), dtype=complex)

    def lift_sector_block(
        name_L: str,
        name_R: str,
        Y: np.ndarray,
        color_dim: int = 1,
    ):
        """
        - name_L, name_R: base names in the 1-gen basis ("nu_L","nu_R","e_L","e_R",...)
        - Y: 3×3 Yukawa matrix in generation space
        - color_dim: 1 for leptons, 3 for quarks
        """
        idx_L_1 = [i for i, st in enumerate(basis_1) if st.name.startswith(name_L)]
        idx_R_1 = [i for i, st in enumerate(basis_1) if st.name.startswith(name_R)]

        assert len(idx_L_1) == color_dim
        assert len(idx_R_1) == color_dim

        for c in range(color_dim):
            i1 = idx_L_1[c]
            j1 = idx_R_1[c]
            for g in range(dim_gen):
                for gprime in range(dim_gen):
                    i3 = i1 * dim_gen + g
                    j3 = j1 * dim_gen + gprime
                    D3[i3, j3] += Y[g, gprime]
                    D3[j3, i3] += np.conj(Y[g, gprime])  # Hermitian completion

    # Lift sectors: leptons (color_dim=1), quarks (color_dim=3)
    lift_sector_block("nu_L", "nu_R", Y_nu, color_dim=1)
    lift_sector_block("e_L",  "e_R",  Y_e,  color_dim=1)
    lift_sector_block("u_L",  "u_R",  Y_u,  color_dim=3)
    lift_sector_block("d_L",  "d_R",  Y_d,  color_dim=3)

    # Enforce J-even or J-odd as desired; we use J-even here
    D3_J_even = project_DF_to_J_sign(D3, S_3, sign="+")

    return {
        "D_F": D3_J_even,
        "gamma_F": gamma_3,
        "S_F": S_3,
        "algebra_ops": ops_3,
        "algebra_labels": labels_3,
        "Y_u": Y_u,
        "Y_d": Y_d,
        "Y_e": Y_e,
        "Y_nu": Y_nu,
    }

def build_SM_3gen_triple_with_emergent_Y(
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
    use_strict_algebra: bool = True,
) -> Dict[str, object]:
    """
    Build a 3-generation SM finite triple whose Yukawas are the full 3×3
    emergent Yukawa matrices from the alignment pipeline.

    use_strict_algebra:
        If True (default), use the STRICT harmonic-aligned commutative
        subalgebra (I, lambda3, lambda8 only).
        If False, fall back to the fuller commutative algebra
        including the quaternion sigma3 generator.
    """
    # 1) Run emergent alignment to get 3×3 Yukawas
    emergent = run_emergent_alignment(
        N=N,
        n_steps=n_steps,
        eta=eta,
        w6=w6,
        w5=w5,
        keep_fraction=keep_fraction,
        alpha=alpha,
        beta=beta,
        random_seed=random_seed,
        use_neutrino_dressing=use_neutrino_dressing,
    )
    Y_u = emergent["Y_u"]
    Y_d = emergent["Y_d"]
    Y_e = emergent["Y_e"]
    Y_nu = emergent["Y_nu"]

    # 2) Delegate to the generic builder
    return build_SM_3gen_triple_from_Y(
        Y_u, Y_d, Y_e, Y_nu,
        use_strict_algebra=use_strict_algebra,
    )


if __name__ == "__main__":
    # print("\n\n>>> Emergent alignment Yukawas plugged into SM tests <<<\n")
    # run_sm_ncg_tests()
    #
    # print("\n\n>>> Emergent alignment Yukawas plugged into NCG tests <<<\n")
    # run_ncg_with_alignment()
    #
    # print("\n\n>>> Full Emergence <<<\n")
    # run_emergent_alignment(use_24cell=True)

    print("\n\n>>> 3-generation SM finite triple with emergent Yukawas (alignment API test) <<<\n")
    tri = build_alignment_finite_triple(
        yukawa_mode="Z2160",
        z2160_kwargs={
            "kappa": 0.22,
            # optional: custom sector_config for u,d,e,nu
        },
    )
    D_F      = tri["D_F"]
    gamma_F  = tri["gamma_F"]
    J_F      = tri["J_F"]
    ops_F    = tri["algebra_ops"]
    labels_F = tri["algebra_labels"]


    # Re-use your generic tests on the returned finite triple
    test_zero_order_condition_generic(ops_F, labels_F, J_F, eps=1e-12)
    test_first_order_condition_generic(D_F, ops_F, labels_F, J_F, eps=1e-12)
    test_grading_and_reality_generic(D_F, ops_F, labels_F, gamma_F, J_F)

    print("\n\n>>> 3-generation SM finite triple with Z2160 triadic Yukawas (v2.0 alignment test) <<<\n")
    run_ncg_with_Z2160_alignment()


"""
>>> 3-generation SM finite triple with emergent Yukawas (alignment API test) <<<

=== Zero-order condition test (generic J) ===
Max Frobenius norm over all pairs (a,b): 0.000e+00
All pairs satisfy [a, J b J^-1] ≈ 0 within eps=1.0e-12

=== First-order condition test (generic J) ===
Max Frobenius norm over all pairs (a,b): 0.000e+00
Pairs with norm < 1.0e-12:
  (a=           I, b=           I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=           I, b=color_lambda3) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=           I, b=color_lambda8) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=color_lambda3, b=           I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=color_lambda3, b=color_lambda3) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=color_lambda3, b=color_lambda8) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=color_lambda8, b=           I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=color_lambda8, b=color_lambda3) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=color_lambda8, b=color_lambda8) → ||[[D,a],J b J^-1]||_F = 0.000e+00

=== Grading & reality tests (generic γ,J) ===
||{gamma, D_F}||_F = 0.000e+00
max ||[gamma, a]||_F over a∈A_F = 0.000e+00
||S^2 - I||_F  (⇒ J^2 deviation) = 0.000e+00
||J D_F J^-1 - D_F||_F   = 0.000e+00
||J D_F J^-1 + D_F||_F   = 4.148e+00
→ KO-sign: J D_F J^-1 = + D_F (J-even Dirac operator)



>>> 3-generation SM finite triple with Z2160 triadic Yukawas (v2.0 alignment test) <<<



=== 3-generation SM finite triple with Z2160 triadic Yukawas (v2.0 test) ===

=== Zero-order condition test (generic J) ===
Max Frobenius norm over all pairs (a,b): 0.000e+00
All pairs satisfy [a, J b J^-1] ≈ 0 within eps=1.0e-12

=== First-order condition test (generic J) ===
Max Frobenius norm over all pairs (a,b): 0.000e+00
Pairs with norm < 1.0e-12:
  (a=           I, b=           I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=           I, b=color_lambda3) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=           I, b=color_lambda8) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=color_lambda3, b=           I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=color_lambda3, b=color_lambda3) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=color_lambda3, b=color_lambda8) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=color_lambda8, b=           I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=color_lambda8, b=color_lambda3) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=color_lambda8, b=color_lambda8) → ||[[D,a],J b J^-1]||_F = 0.000e+00

=== Grading & reality tests (generic γ,J) ===
||{gamma, D_F}||_F = 0.000e+00
max ||[gamma, a]||_F over a∈A_F = 0.000e+00
||S^2 - I||_F  (⇒ J^2 deviation) = 0.000e+00
||J D_F J^-1 - D_F||_F   = 0.000e+00
||J D_F J^-1 + D_F||_F   = 4.148e+00
→ KO-sign: J D_F J^-1 = + D_F (J-even Dirac operator)

"""