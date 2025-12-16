import numpy as np
import cmath

class GeometricTriple:
    """
    Geometric spectral triple (A_geom, H_geom, D_geom)
    with divisor projector C360 and phase projector Pphi.
    """

    def __init__(self, N_max=360):
        self.N_max = N_max
        self.divisors = self._compute_divisors(360)

    @staticmethod
    def _compute_divisors(n):
        return [d for d in range(1, n+1) if n % d == 0]

    def C360(self, n):
        return 1 if n in self.divisors else 0

    def B(self, n):
        return 1

    def Pphi(self, n, phi):
        return cmath.exp(1j * phi * n)

    def A_projector(self, n, phi):
        return self.C360(n) * self.B(n) * self.Pphi(n, phi)

    def D_geom(self, modes):
        """
        Build geometric Dirac operator diag(n).
        """
        return np.diag(modes)

    def gamma_geom(self, modes):
        """
        Gamma = sign(n)
        """
        return np.diag([1 if m >= 0 else -1 for m in modes])

    def J_geom(self, n):
        """
        Complex conjugation: e^{inx} -> e^{-inx}
        """
        return -n

class FiniteTriple:
    """
    Finite triadic sector: (A_F, H_F, D_F)
    D_F is the 9×9 triadic Z2160 kernel K with phases.
    """

    def __init__(self, kernel, phases):
        self.kernel = kernel          # real K_{ij}
        self.phases = np.array(phases)
        self.phase_matrix = self._build_phase_matrix()
        self.D_F = self.phase_matrix * self.kernel

    def _build_phase_matrix(self):
        phi = self.phases
        return np.exp(1j * (phi[:,None] - phi[None,:]))

    def J_F(self, vec):
        return np.conjugate(vec)

    def gamma_F(self, signs):
        """
        signs: list of +1/-1 of length 9
        """
        return np.diag(signs)

class TriadicKernel:
    """
    Build K_{ij} from:
    - forbidden distances {2,4,7}
    - kappa^{distance}
    - Z_2160 distance metric
    """

    def __init__(self, kappa=0.24, forbidden=(2,4,7), N=2160, sites=None):
        self.kappa = kappa
        self.forbidden = set(forbidden)
        self.N = N
        if sites is None:
            raise ValueError("Must specify 9 embedded sites in Z_2160.")
        self.sites = sites           # list of 9 integers in [0,N-1]
        self.K = self._build_kernel()

    def dist(self, a, b):
        """
        distance on the cycle Z_N
        """
        d = abs(a-b)
        return min(d, self.N - d)

    def _build_kernel(self):
        K = np.zeros((9,9), dtype=float)
        for i in range(9):
            for j in range(9):
                if i == j:
                    K[i,j] = 1.0
                else:
                    d = self.dist(self.sites[i], self.sites[j])
                    if d in self.forbidden:
                        K[i,j] = 0.0
                    else:
                        K[i,j] = self.kappa ** d
        return K

class TriadicCompression:
    """
    Triadic 9→3 projection:
    - triads T1, T2, T3 (each size 3)
    - build S (3×9)
    - compress Y^(s) = S K^(s) S†
    """

    def __init__(self, triads):
        self.triads = triads
        self.S = self._build_S()

    def _build_S(self):
        S = np.zeros((3,9), dtype=float)
        for a, triad in enumerate(self.triads):
            for idx in triad:
                S[a,idx] = 1/np.sqrt(3)
        return S

    def compress(self, K_sector):
        S = self.S
        return S @ K_sector @ S.conjugate().T

    def sector_kernel(self, K, R):
        """
        K: 9×9 complex
        R: 9-element weight vector
        """
        Rm = np.diag(R)
        return Rm @ K @ Rm

class AlignmentSpectralTriple:
    def __init__(self, geom, finite):
        self.geom = geom
        self.finite = finite

    def build_product_dirac(self, modes):
        """
        D = D_geom ⊗ 1 + gamma_geom ⊗ D_F
        """
        Dg = self.geom.D_geom(modes)
        gamma = self.geom.gamma_geom(modes)
        DF = self.finite.D_F

        term1 = np.kron(Dg, np.eye(9))
        term2 = np.kron(gamma, DF)

        return term1 + term2

    def inner_fluctuation(self, A):
        """
        Placeholder for A + JAJ^{-1}.
        A can be geometric or internal fluctuation.
        """
        return A

class SpectralAction:
    """
    Computes:
    S = Tr f(D_A^2/Λ^2) + <Ψ, D_A Ψ>
    Uses truncated heat-kernel expansion.
    """

    def __init__(self, cutoff=1e14):
        self.Lambda = cutoff

    def D_squared(self, D):
        return D @ D

    def heat_kernel_coeffs(self, D2):
        a0 = np.trace(np.eye(D2.shape[0]))
        a2 = np.trace(D2)
        a4 = np.trace(D2 @ D2)
        return a0, a2, a4

    def bosonic_action(self, D):
        D2 = self.D_squared(D)
        a0, a2, a4 = self.heat_kernel_coeffs(D2)
        f4, f2, f0 = 1, 1, 1  # can later be user inputs
        L = self.Lambda
        S = f4*L**4*a0 + f2*L**2*a2 + f0*a4
        return S.real


class GaugeSector:
    """
    Placeholder SU(3)_A, SU(2)_A, U(1)_A fields
    from inner fluctuations of the algebra.
    """

    def __init__(self):
        pass

    def gauge_field_U1(self, size):
        return np.zeros((size,size), dtype=complex)

    def gauge_field_SU2(self, size):
        return np.zeros((size,size), dtype=complex)

    def gauge_field_SU3(self, size):
        return np.zeros((size,size), dtype=complex)

def normalize_matrix(M):
    s = np.linalg.svd(M, compute_uv=False)
    if s[0] != 0:
        return M / s[0]
    return M
