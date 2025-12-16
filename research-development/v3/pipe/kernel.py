import numpy as np
from typing import Iterable, Sequence, List, Set


class TriadicKernel:
    """
    Triadic Z_2160 Finite Kernel for the Alignment Spectral Triple v2.0.

    This class constructs the 9×9 real kernel K_{ij} encoding:
        • geometric distances on the lattice Z_N (default N=2160),
        • a forbidden-distance set 𝓕 = {2,4,7},
        • geometric decay K_{ij} = κ^d for allowed d,
        • normalized diagonal terms K_{ii} = 1.

    The kernel is used as the base finite Dirac data before adding phases.
    It is guaranteed to be real, symmetric, and positive along the diagonal.
    """

    # ------------------------------------------------------------
    # Initialization
    # ------------------------------------------------------------
    def __init__(
        self,
        kappa: float = 0.24,
        forbidden: Sequence[int] = (2, 4, 7),
        N: int = 2160,
        sites: Sequence[int] = None
    ) -> None:
        """
        Parameters
        ----------
        kappa : float
            Decay parameter controlling geometric suppression (0 < κ < 1).
        forbidden : sequence of int
            Forbidden distances; any pair separated by these distances has K_ij = 0.
        N : int
            Size of the cyclic lattice Z_N.
        sites : sequence of int
            The 9 embedded site positions in Z_N. Must be length 9, each in [0, N-1].

        Notes
        -----
        This class constructs only the *real* kernel K_{ij}.
        Phases are applied in the FiniteTriple class to construct D_F.
        """
        self._validate_params(kappa, forbidden, N, sites)

        self.kappa: float = float(kappa)
        self.forbidden: Set[int] = set(int(d) for d in forbidden)
        self.N: int = int(N)
        self.sites: List[int] = [int(s) for s in sites]

        # Build real kernel
        self.K: np.ndarray = self._build_kernel()

    # ------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------
    @staticmethod
    def _validate_params(
        kappa: float,
        forbidden: Sequence[int],
        N: int,
        sites: Sequence[int]
    ) -> None:
        if not (0 < kappa < 1):
            raise ValueError("kappa must satisfy 0 < kappa < 1.")

        if any(d <= 0 for d in forbidden):
            raise ValueError("Forbidden distances must be positive integers.")

        if N <= 0:
            raise ValueError("N must be a positive integer.")

        if sites is None:
            raise ValueError("Must specify 9 embedded sites in Z_N.")
        if len(sites) != 9:
            raise ValueError("sites must contain exactly 9 integers.")
        if any((s < 0 or s >= N) for s in sites):
            raise ValueError("sites must lie in the interval [0, N-1].")

    # ------------------------------------------------------------
    # Z_N distance
    # ------------------------------------------------------------
    def dist(self, a, b):
        d = abs(a - b)
        d = min(d, self.N - d)
        return max(0, min(d, 15))  # clip to 0–15 for intra-triad stability

    # ------------------------------------------------------------
    # Kernel construction
    # ------------------------------------------------------------
    def _build_kernel(self) -> np.ndarray:
        """
        Construct the real kernel K_{ij}.

        Returns
        -------
        np.ndarray
            9×9 real symmetric matrix K.
        """
        K = np.zeros((9, 9), dtype=float)

        for i in range(9):
            for j in range(9):
                if i == j:
                    # Normalized diagonal
                    K[i, j] = 1.0
                else:
                    d = self.dist(self.sites[i], self.sites[j])
                    if d in self.forbidden:
                        K[i, j] = 0.0
                    else:
                        K[i, j] = self.kappa ** d

        # Symmetry sanity check (should hold exactly)
        if not np.allclose(K, K.T, atol=1e-15):
            raise ValueError("Kernel is not symmetric; internal error detected.")

        return K

    # ------------------------------------------------------------
    # Convenience
    # ------------------------------------------------------------
    def as_array(self) -> np.ndarray:
        """
        Return the kernel as a numpy array (copy).

        Returns
        -------
        np.ndarray
            Real symmetric 9×9 kernel.
        """
        return np.array(self.K, copy=True)

    def print_summary(self) -> None:
        """
        Print a human-readable summary of the kernel details.
        """
        print("TriadicKernel Summary")
        print("---------------------")
        print(f"  κ       : {self.kappa}")
        print(f"  N       : {self.N}")
        print(f"  Forbidden distances : {sorted(self.forbidden)}")
        print(f"  Sites   : {self.sites}")
        print("  Kernel shape:", self.K.shape)
        print()

"""
Concrete upgrades I’d make to this class (so it’s “v5-ready”)
1) Parameterize d_max + clarify forbidden semantics
def __init__(..., d_max: int | None = 15, forbid_on: str = "raw"):
    ...


d_max=None means no clipping

forbid_on="raw" checks forbidden on raw distance, then optionally clips for decay

forbid_on="clipped" keeps your current behavior

2) Vectorize (faster + simpler)

Right now it’s nested loops. Fine for 9×9, but if you later generalize to bigger N-site models you’ll want vectorization.

3) Export an adjacency / Laplacian for anchoring

Your newer anchor logic wants a clock/Laplacian 
𝐿
𝑁
L
N
	​

. This kernel is already “weights on edges”, so you can define:

𝑊
𝑖
𝑗
=
𝐾
𝑖
𝑗
W
ij
	​

=K
ij
	​

 for 
𝑖
≠
𝑗
i

=j (maybe threshold small values)

Graph Laplacian:

𝐿
=
d
i
a
g
(
𝑊
1
)
−
𝑊
L=diag(W1)−W

Then your 
Π
360
Π
360
	​

 logic can be built from the spectrum of 
𝐿
L (or of a clock operator derived from it).
"""