# =================================================================
# triple.py
# Full Emergence Architecture — Module 6
# =================================================================
# Builds the fully emergent internal spectral triple:
#
#    (A_int, H_int, D_int)
#
# All components arise from:
#   - the internal quasi-crystal graph
#   - L (Laplacian)
#   - R (360-rotation)
#   - Q (integer charge operator)
#   - Yukawas Y_s
#   - Emergent symmetry blocks from manifestation operator X
#
# Absolutely no external parameters.
# =================================================================

import numpy as np
import scipy


# ---------------------------------------------------------------
# 1. Internal Hilbert space  H_int  (Axiom A1)
# ---------------------------------------------------------------
class InternalHilbertSpace:
    """
    Finite-dimensional internal Hilbert space representing:
        - quasi-crystal modes
        - generation subspace
        - flavor charge structure

    FULL EMERGENCE:
    Dim(H_int) = dimension of the largest connected component
    of the quasi-crystal graph.
    """

    def __init__(self, N):
        self.N = N
        self.H = np.eye(N, dtype=complex)

    def dim(self):
        return self.N

    def basis(self):
        return np.eye(self.N, dtype=complex)


# ---------------------------------------------------------------
# 2. Internal algebra A_int  (Axiom B + D + F)
# ---------------------------------------------------------------
def build_internal_algebra(R, Q, C360, B, Pphi):
    """
    Axiom G1 minimality:
    The operator set {R, Q, C360, B, Pphi} generates A_int.

    We construct the smallest unital *-algebra containing them:
         A_int = alg{ I, R, Q, C360, B, Pphi }.
    """

    gens = [np.eye(R.shape[0]), R, Q, C360, B, Pphi]

    # Closure under Hermitian conjugation and composition
    # (finite system → brute-force stabilization)
    changed = True
    A = gens.copy()

    while changed:
        changed = False

        # star closure
        for M in list(A):
            Mdag = M.conj().T
            if not any(np.allclose(Mdag, X) for X in A):
                A.append(Mdag)
                changed = True

        # multiplication closure
        for M in list(A):
            for N in list(A):
                MN = M @ N
                if not any(np.allclose(MN, X) for X in A):
                    A.append(MN)
                    changed = True

    return A


# ---------------------------------------------------------------
# 3. Finite Dirac operator  D_int  (Axiom F + NCG structure)
# ---------------------------------------------------------------
def build_internal_Dirac(Y, R, Q):
    """
    Axiom F1-F3: Yukawa structure supplies the Dirac operator.

    We follow the NCG prescription:
        D_int = direct_sum_s ( Y_s + Y_s^† )

    PLUS emergent geometric terms built from commutators with R and Q:
        D_R = [R, Y_s]
        D_Q = [Q, Y_s]

    All components are summed with unit weight (FULL EMERGENCE).
    """

    # Identify sector list
    sectors = list(Y.keys())

    # Start with block diagonal Yukawa contribution
    blocks = []
    for s in sectors:
        Y_s = Y[s]
        Ds = Y_s + Y_s.conj().T
        blocks.append(Ds)

    D0 = scipy.linalg.block_diag(*blocks)

    # Add emergent commutator contributions
    D_R = np.zeros_like(D0)
    D_Q = np.zeros_like(D0)

    offset = 0
    for s in sectors:
        Y_s = Y[s]
        dim = Y_s.shape[0]

        R_s = R[:dim, :dim]
        Q_s = Q[:dim, :dim]

        D_R[offset:offset+dim, offset:offset+dim] = R_s @ Y_s - Y_s @ R_s
        D_Q[offset:offset+dim, offset:offset+dim] = Q_s @ Y_s - Y_s @ Q_s

        offset += dim

    D_int = D0 + D_R + D_Q
    return D_int


# ---------------------------------------------------------------
# 4. Full internal triple
# ---------------------------------------------------------------
class InternalSpectralTriple:
    """
    The full triple (A_int, H_int, D_int) entirely emergent.
    """

    def __init__(self, H_int, A_int, D_int):
        self.H = H_int
        self.A = A_int
        self.D = D_int

    # Zero-order condition [A, J B J^{-1}] = 0 equivalent
    def test_zero_order(self, J):
        print("Testing zero order...")
        max_norm = 0
        for a in self.A:
            for b in self.A:
                comm = a @ (J @ b @ J.T) - (J @ b @ J.T) @ a
                n = np.linalg.norm(comm)
                max_norm = max(max_norm, n)
        print("Zero-order max norm:", max_norm)
        return max_norm

    # First-order [[D,a], J b J^{-1}] = 0
    def test_first_order(self, J):
        print("Testing first order...")
        max_norm = 0
        for a in self.A:
            Da = self.D @ a - a @ self.D
            for b in self.A:
                Jb = J @ b @ J.T
                comm = Da @ Jb - Jb @ Da
                n = np.linalg.norm(comm)
                max_norm = max(max_norm, n)
        print("First-order max norm:", max_norm)
        return max_norm

    # Grading & reality tests
    def test_grading_reality(self, J, gamma):
        print("Testing grading & reality...")
        anti = gamma @ self.D + self.D @ gamma
        print("||{γ, D}|| =", np.linalg.norm(anti))

        S2 = J @ J
        print("||J^2 - I|| =", np.linalg.norm(S2 - np.eye(S2.shape[0])))

        JDJ = J @ self.D @ J.T
        print("||J D J^-1 - D|| =", np.linalg.norm(JDJ - self.D))
        print("||J D J^-1 + D|| =", np.linalg.norm(JDJ + self.D))