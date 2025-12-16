import numpy as np
from typing import Sequence, List


class TriadicCompression:
    """
    Triadic 9→3 projection operator for the Alignment Spectral Triple v2.0.

    This class implements:

        • A triadic partition of the 9 internal sites:
              I = T1 ⊔ T2 ⊔ T3,  |T_a| = 3,
        • The projection matrix S : ℂ^9 → ℂ^3 defined by
              S_{a,i} = 1/sqrt(3) if i ∈ T_a else 0,
        • Sector-dependent kernels:
              K^(s) = R^(s) ⋅ K ⋅ R^(s),
        • The effective 3×3 Yukawa matrix:
              Y^(s) = S ⋅ K^(s) ⋅ S†.

    All dimensionality and partitioning is strictly validated.
    """

    # ------------------------------------------------------------
    # Initialization
    # ------------------------------------------------------------
    def __init__(self, triads: Sequence[Sequence[int]]) -> None:
        """
        Parameters
        ----------
        triads : sequence of three sequences
            Each element must be a sequence of exactly 3 site indices.
            Together they must form a partition of {0,…,8}.

        Example
        -------
        triads = [
            [0,1,2],
            [3,4,5],
            [6,7,8]
        ]
        """
        self._validate_triads(triads)
        self.triads: List[List[int]] = [list(t) for t in triads]

        # Build S ∈ ℝ^{3×9}
        self.S: np.ndarray = self._build_S()

    # ------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------
    @staticmethod
    def _validate_triads(triads: Sequence[Sequence[int]]) -> None:
        if len(triads) != 3:
            raise ValueError("triads must contain exactly 3 triads (3 groups).")

        flat = []
        for t in triads:
            if len(t) != 3:
                raise ValueError("Each triad must have exactly 3 indices.")
            for x in t:
                flat.append(x)

        # Must partition the 9 sites
        if sorted(flat) != list(range(9)):
            raise ValueError(
                "Triads must form a partition of {0,…,8} with no duplicates."
            )

    # ------------------------------------------------------------
    # Build S matrix
    # ------------------------------------------------------------
    def _build_S(self) -> np.ndarray:
        """
        Construct the projection matrix S_{ai}, mapping ℂ^9 → ℂ^3:

            |t_a⟩ = (1/√3) ∑_{i∈T_a} |i⟩

        Returns
        -------
        np.ndarray
            Real 3×9 projection matrix.
        """
        S = np.zeros((3, 9), dtype=float)
        for a, triad in enumerate(self.triads):
            for idx in triad:
                S[a, idx] = 1.0 / np.sqrt(3.0)
        return S

    # ------------------------------------------------------------
    # Sector Kernel (R K R)
    # ------------------------------------------------------------
    @staticmethod
    def sector_kernel(K: np.ndarray, R: Sequence[float]) -> np.ndarray:
        """
        Build the sector-specific 9×9 finite kernel:

            K^(s) = R ⋅ K ⋅ R

        Parameters
        ----------
        K : np.ndarray (9×9)
            Base Hermitian kernel (after phase multiplication).
        R : sequence of float, length 9
            Sector weight vector (positive weights).

        Returns
        -------
        np.ndarray
            The 9×9 sector-dependent kernel.

        Raises
        ------
        ValueError if inputs are malformed.
        """
        if not isinstance(K, np.ndarray):
            raise TypeError("K must be a numpy array.")
        if K.shape != (9, 9):
            raise ValueError("K must be a 9×9 matrix.")
        if len(R) != 9:
            raise ValueError("R must be a 9-element weight vector.")

        Rm = np.diag(np.array(R, dtype=float))
        K_sector = Rm @ K @ Rm

        return K_sector

    # ------------------------------------------------------------
    # Compression Y^(s) = S K^(s) S†
    # ------------------------------------------------------------
    def compress(self, K_sector: np.ndarray) -> np.ndarray:
        """
        Compress a 9×9 sector kernel to an effective 3×3 Yukawa matrix:

            Y^(s) = S ⋅ K^(s) ⋅ S†

        Parameters
        ----------
        K_sector : np.ndarray (9×9)
            Sector-tempered kernel for one sector (u, d, e, ν).

        Returns
        -------
        np.ndarray
            3×3 complex Yukawa matrix.

        Raises
        ------
        ValueError for incorrect matrix shape.
        """
        if not isinstance(K_sector, np.ndarray):
            raise TypeError("K_sector must be a numpy array.")
        if K_sector.shape != (9, 9):
            raise ValueError("K_sector must be a 9×9 matrix.")

        S = self.S
        Y = S @ K_sector @ S.conjugate().T

        # Hermiticity should hold up to machine precision.
        if not np.allclose(Y, Y.conjugate().T, atol=1e-12):
            raise ValueError(
                "Compressed Yukawa matrix is not Hermitian — check K and R."
            )

        return Y

    # ------------------------------------------------------------
    # Convenience
    # ------------------------------------------------------------
    def print_summary(self) -> None:
        """Print a summary of the triadic configuration."""
        print("TriadicCompression Summary")
        print("--------------------------")
        print(f"Triads: {self.triads}")
        print("Projection matrix S (3×9):")
        print(self.S)
        print()
