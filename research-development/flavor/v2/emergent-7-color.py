# ================================
# Internal Hilbert space & D_F
# ================================
import numpy as np
from typing import List, Tuple, Dict

SECTORS = ["u", "d", "e", "nu"]
SECTOR_INDEX: Dict[str, int] = {s: i for i, s in enumerate(SECTORS)}
N_GEN = 3
# Color multiplicities (degeneracies). In this toy, we don't explicitly
# tensor out color, we just keep track of the full 24-dim per chirality.
SECTOR_NC = {"u": 3, "d": 3, "e": 1, "nu": 1}


def dim_per_chirality() -> int:
    """
    Dimension of H_L or H_R (one chirality).

    We treat color multiplicities as degeneracy, not as an explicit tensor
    factor for now, but the total number of internal states per chirality
    still comes out as:

        dim(H_L) = N_GEN * sum_s N_c(s) = 3 * (3+3+1+1) = 24

    This matches your emergent-5 setup.
    """
    return N_GEN * sum(SECTOR_NC[s] for s in SECTORS)  # 24


def flavor_block_offsets() -> Dict[str, int]:
    """
    Return offsets (within the *generation subspace*) for each sector's 3×3
    generation block in a 12×12 layout:

      [u_g1,u_g2,u_g3,
       d_g1,d_g2,d_g3,
       e_g1,e_g2,e_g3,
       nu_g1,nu_g2,nu_g3]

    We only care about the leading 12 entries as "generation space".
    The remaining 12 (per chirality) are currently unused / reserved
    for future color-explicit extensions.
    """
    off: Dict[str, int] = {}
    off["u"]  = 0
    off["d"]  = 3
    off["e"]  = 6
    off["nu"] = 9
    return off


# -------------------------------------------------------------------
# F-based D_F builder (diagonal Yukawas) – optional fallback
# -------------------------------------------------------------------
def build_internal_DF(F_u: np.ndarray,
                      F_d: np.ndarray,
                      F_e: np.ndarray,
                      F_n: np.ndarray) -> np.ndarray:
    """
    Build the finite Dirac operator D_F in block form:

      D_F = [[ 0, Y^\dagger ],
             [ Y, 0       ]]

    where Y is a 24×24 block that is block-diagonal in sector space and
    uses 3×3 diagonal generation Yukawas diag(F_s) in a 12×12
    generation-space layout (color folded in as degeneracy).

    H_F ≃ H_L ⊕ H_R,  dim(H_L) = dim(H_R) = 24, dim(H_F) = 48.
    """
    # Sanity checks
    for name, F in [("F_u", F_u), ("F_d", F_d), ("F_e", F_e), ("F_n", F_n)]:
        F_arr = np.asarray(F, dtype=float)
        if F_arr.shape != (3,):
            raise ValueError(f"{name} must be a length-3 array, got shape {F_arr.shape}.")

    dpc = dim_per_chirality()  # 24
    dimH = 2 * dpc             # 48

    # 12×12 generation-space Yukawa core; then embedded into 24×24 per chirality
    Y_gen = np.zeros((12, 12), dtype=complex)

    Y_u  = np.diag(F_u)
    Y_d  = np.diag(F_d)
    Y_e  = np.diag(F_e)
    Y_nu = np.diag(F_n)

    gen_off = flavor_block_offsets()

    def insert_sector_Y(sector: str, Y_s: np.ndarray) -> None:
        off = gen_off[sector]
        Y_gen[off:off+3, off:off+3] = Y_s

    insert_sector_Y("u",  Y_u)
    insert_sector_Y("d",  Y_d)
    insert_sector_Y("e",  Y_e)
    insert_sector_Y("nu", Y_nu)

    # Embed this 12×12 generation block into 24×24 per chirality.
    Y_block = np.zeros((dpc, dpc), dtype=complex)
    Y_block[:12, :12] = Y_gen

    # Assemble D_F
    D_F = np.zeros((dimH, dimH), dtype=complex)
    D_F[:dpc,  dpc:] = Y_block.conj().T
    D_F[dpc:,  :dpc] = Y_block

    return D_F


# -------------------------------------------------------------------
# Y-based D_F builder (full Yukawa matrices, with mixing)
# -------------------------------------------------------------------
def build_internal_DF_from_Y(Y_u: np.ndarray,
                             Y_d: np.ndarray,
                             Y_e: np.ndarray,
                             Y_nu: np.ndarray) -> np.ndarray:
    """
    Build the finite Dirac operator D_F in block form:

      D_F = [[ 0, Y^\dagger ],
             [ Y, 0        ]]

    where Y is a 24×24 block, block-diagonal in sector space, with
    3×3 Yukawa matrices per sector (not necessarily diagonal):

      Y_gen = diag( Y_u, Y_d, Y_e, Y_nu ) in the leading 12×12 generation space.

    The remaining 12 entries per chirality are unused (reserved for future
    explicit color structure). For now they are set to zero.

    H_F ≃ H_L ⊕ H_R, dim(H_L)=dim(H_R)=24, dim(H_F)=48.
    """
    # Sanity checks on shapes
    for name, Y in [("Y_u", Y_u), ("Y_d", Y_d), ("Y_e", Y_e), ("Y_nu", Y_nu)]:
        Y_arr = np.asarray(Y, dtype=complex)
        if Y_arr.shape != (3, 3):
            raise ValueError(f"{name} must be a 3×3 matrix, got shape {Y_arr.shape}.")

    dpc = dim_per_chirality()  # 24
    dimH = 2 * dpc             # 48

    # Build 12×12 generation block first
    Y_gen = np.zeros((12, 12), dtype=complex)
    gen_off = flavor_block_offsets()

    def insert_sector_Y(sector: str, Y_s: np.ndarray) -> None:
        off = gen_off[sector]
        Y_gen[off:off+3, off:off+3] = Y_s

    insert_sector_Y("u",  Y_u)
    insert_sector_Y("d",  Y_d)
    insert_sector_Y("e",  Y_e)
    insert_sector_Y("nu", Y_nu)

    # Embed generation block into 24×24 per chirality
    Y_block = np.zeros((dpc, dpc), dtype=complex)
    Y_block[:12, :12] = Y_gen

    # Assemble D_F
    D_F = np.zeros((dimH, dimH), dtype=complex)
    D_F[:dpc,  dpc:] = Y_block.conj().T
    D_F[dpc:,  :dpc] = Y_block

    return D_F


# ================================
# Real structure, grading, algebra basis
# ================================

def build_swap_LR(dim_left: int) -> np.ndarray:
    """
    Build the swap matrix S on H_F = H_L ⊕ H_R, where dim(H_L) = dim(H_R) = dim_left.
    Acts as:
      S ( ψ_L, ψ_R ) = ( ψ_R, ψ_L )
    """
    S = np.zeros((2*dim_left, 2*dim_left), dtype=complex)
    S[:dim_left, dim_left:] = np.eye(dim_left)
    S[dim_left:, :dim_left] = np.eye(dim_left)
    return S


def build_gamma_F(dim_left: int) -> np.ndarray:
    """
    Grading operator γ_F with eigenvalue -1 on H_L and +1 on H_R.
    """
    g = np.zeros((2*dim_left, 2*dim_left), dtype=complex)
    g[:dim_left, :dim_left] = -np.eye(dim_left)
    g[dim_left:, dim_left:] =  np.eye(dim_left)
    return g


def build_sector_projectors() -> Dict[str, np.ndarray]:
    """
    Build sector projectors P_sector_s acting on H_F = H_L ⊕ H_R.
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
        # Act the same on L and R (block-diagonal), and only on the first 12 gen slots
        P_s[off:off+3, off:off+3] = np.eye(3)
        P_s[dpc + off:dpc + off + 3, dpc + off:dpc + off + 3] = np.eye(3)
        P[s] = P_s

    return P  # dict with keys "u","d","e","nu"


def build_Q_sector() -> np.ndarray:
    """
    Build a simple 'sector charge' diagonal operator Q_sector which distinguishes
    u,d,e,nu sectors but is generation-blind.
    Example charges:
      q_u = 2, q_d = 1, q_e = 0, q_nu = -1
    """
    dpc = dim_per_chirality()
    dimH = 2 * dpc

    gen_off = flavor_block_offsets()
    charges = {"u": 2.0, "d": 1.0, "e": 0.0, "nu": -1.0}

    Q = np.zeros((dimH, dimH), dtype=complex)
    for s in SECTORS:
        off = gen_off[s]
        q   = charges[s]
        Q[off:off+3, off:off+3] = q * np.eye(3)
        Q[dpc + off:dpc + off + 3, dpc + off:dpc + off + 3] = q * np.eye(3)

    return Q


def build_internal_algebra_ops() -> Tuple[List[np.ndarray], List[str]]:
    """
    Build a small basis of algebra elements A_F acting on H_F:
      - I (identity)
      - Q_sector (diagonal sector 'charge')
      - P_sector_u, P_sector_d, P_sector_e, P_sector_nu (sector projectors)

    This is a commutative algebra in this toy (we are not yet including full SU(3)).
    """
    dpc = dim_per_chirality()
    dimH = 2 * dpc

    I = np.eye(dimH, dtype=complex)
    Q = build_Q_sector()
    P = build_sector_projectors()

    ops: List[np.ndarray] = [I, Q, P["u"], P["d"], P["e"], P["nu"]]
    labels: List[str] = ["I", "Q_sector", "P_sector_u", "P_sector_d",
                         "P_sector_e", "P_sector_nu"]

    return ops, labels


# ================================
# NCG condition tests
# ================================

def J_action_from_swap(S: np.ndarray, M: np.ndarray) -> np.ndarray:
    """
    Implement J M J^{-1} = S * M^* * S^T, where S is the L/R swap.
    """
    return S @ M.conj() @ S.T


def test_first_order_condition(D_F: np.ndarray,
                               ops: List[np.ndarray],
                               labels: List[str],
                               eps: float = 1e-12) -> None:
    """
    First-order condition:
      [[D_F, a], J_F b J_F^{-1}] = 0
    for all a,b in algebra.
    """
    n = D_F.shape[0]
    assert D_F.shape == (n, n)
    S = build_swap_LR(dim_left=n//2)

    print("=== First-order condition test ===")
    max_norm = 0.0
    good_pairs = []

    for i, a in enumerate(ops):
        Da = D_F @ a - a @ D_F
        for j, b in enumerate(ops):
            b_tilde = J_action_from_swap(S, b)
            comm2 = Da @ b_tilde - b_tilde @ Da
            norm = np.linalg.norm(comm2, ord='fro')

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
    """
    Zero-order condition:
      [a, J_F b J_F^{-1}] = 0
    for all a,b in algebra.
    """
    n = ops[0].shape[0]
    S = build_swap_LR(dim_left=n//2)

    print("=== Zero-order condition test ===")
    max_norm = 0.0
    bad_pairs = []

    for i, a in enumerate(ops):
        for j, b in enumerate(ops):
            b_tilde = J_action_from_swap(S, b)
            comm = a @ b_tilde - b_tilde @ a
            norm = np.linalg.norm(comm, ord='fro')
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
    """
    - Check γ_F anticommutes with D_F and commutes with A_F.
    - Check J_F^2 = 1 (as implemented by swap).
    - Detect the KO-dimension sign via:
        J D_F J^{-1} = ± D_F
    """
    n = D_F.shape[0]
    dpc = n // 2
    gamma_F = build_gamma_F(dpc)
    S = build_swap_LR(dpc)

    print("=== Grading & reality tests ===")
    # γ_F anti-commutes with D_F
    anti = gamma_F @ D_F + D_F @ gamma_F
    print(f"||{{gamma_F, D_F}}||_F = {np.linalg.norm(anti, ord='fro'):.3e}")

    # γ_F commutes with algebra
    max_comm_gamma = 0.0
    for a in ops:
        comm_ga = gamma_F @ a - a @ gamma_F
        max_comm_gamma = max(max_comm_gamma, np.linalg.norm(comm_ga, ord='fro'))
    print(f"max ||[gamma_F, a]||_F over a∈A_F = {max_comm_gamma:.3e}")

    # J_F^2 = 1 (swap^2 = I)
    S2 = S @ S
    print(f"||S^2 - I||_F  (⇒ J_F^2 deviation) = {np.linalg.norm(S2 - np.eye(n), ord='fro'):.3e}")

    # J D_F J^{-1} vs ± D_F
    JDJ = S @ D_F.conj() @ S.T
    norm_plus  = np.linalg.norm(JDJ - D_F, ord='fro')
    norm_minus = np.linalg.norm(JDJ + D_F, ord='fro')

    print(f"||J D_F J^-1 - D_F||_F   = {norm_plus:.3e}")
    print(f"||J D_F J^-1 + D_F||_F   = {norm_minus:.3e}")

    if norm_plus < 1e-12 and norm_minus > norm_plus:
        print("→ KO-sign: J D_F J^-1 = + D_F (J-even Dirac operator)")
    elif norm_minus < 1e-12 and norm_plus > norm_minus:
        print("→ KO-sign: J D_F J^-1 = - D_F (J-odd Dirac operator)")
    else:
        print("→ KO-sign ambiguous or not clean at numerical precision.")
    print()


# ================================
# Minimal example / entry point
# ================================

def example_Fs() -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Provide a simple example set of F_s triads so this file can be run
    standalone. In your full emergent model, replace these with the
    F_u, F_d, F_e, F_n you compute from the internal graph + Q.
    """
    # Example: slightly hierarchical triad (not meant to match SM)
    F_base = np.array([0.05, 0.005, 0.0005], dtype=float)

    # Simple integer exponents per sector (toy, not SM):
    q_u  = np.array([0,  2, 4], dtype=float)
    q_d  = np.array([1,  3, 5], dtype=float)
    q_e  = np.array([2,  4, 6], dtype=float)
    q_nu = np.array([4,  6, 8], dtype=float)

    beta = 1.0

    def sector_weights(Fb: np.ndarray, q: np.ndarray, beta_val: float) -> np.ndarray:
        return Fb * np.exp(-beta_val * q)

    F_u  = sector_weights(F_base, q_u,  beta)
    F_d  = sector_weights(F_base, q_d,  beta)
    F_e  = sector_weights(F_base, q_e,  beta)
    F_n  = sector_weights(F_base, q_nu, beta)

    return F_u, F_d, F_e, F_n


def main() -> None:
    # Example Yukawa triads (replace with your emergent F_s in the full model)
    F_u, F_d, F_e, F_n = example_Fs()

    # For this minimal example, we build Y_s as simple diagonal matrices.
    # In your full emergent pipeline, REPLACE these with your actual emergent
    # Yukawa matrices, e.g. Y_u = U_L_u @ np.diag(F_u) @ U_R_u.conj().T, etc.
    Y_u  = np.diag(F_u)
    Y_d  = np.diag(F_d)
    Y_e  = np.diag(F_e)
    Y_nu = np.diag(F_n)

    # --- Build internal Dirac from full Yukawas ---
    D_F = build_internal_DF_from_Y(Y_u, Y_d, Y_e, Y_nu)

    # --- Build algebra basis (same as before) ---
    ops_A, labels_A = build_internal_algebra_ops()

    # --- Run NCG tests ---
    test_first_order_condition(D_F, ops_A, labels_A, eps=1e-12)
    test_zero_order_condition(ops_A, labels_A, eps=1e-12)
    test_grading_and_reality(D_F, ops_A, labels_A)


if __name__ == "__main__":
    main()

"""
RESULTS:
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
||J D_F J^-1 - D_F||_F   = 0.000e+00
||J D_F J^-1 + D_F||_F   = 1.519e-01
→ KO-sign: J D_F J^-1 = + D_F (J-even Dirac operator)
"""