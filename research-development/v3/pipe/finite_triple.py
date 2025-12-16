import numpy as np
from typing import Sequence, List


class FiniteTriple:
    """
    Finite Triadic Sector (A_F, H_F, D_F) of the Alignment Spectral Triple v2.0.

    This class implements:
        • A_F = ℂ^9  acting diagonally on |i⟩ for i=0,…,8,
        • H_F = ℂ^9  with canonical orthonormal basis,
        • D_F = 𝒦 = P ∘ K  where:
              K is the real Z_2160 triadic kernel with forbidden distances,
              P is the phase matrix P_ij = exp(i(φ_i - φ_j)).

    The resulting operator 𝒦 is Hermitian and serves as the finite Dirac operator.
    """

    # ------------------------------------------------------------
    # Initialization
    # ------------------------------------------------------------
    def __init__(self, kernel: np.ndarray, phases: Sequence[float]) -> None:
        """
        Parameters
        ----------
        kernel : np.ndarray
            The real 9×9 triadic kernel K_{ij} (with forbidden distances applied).
        phases : sequence of float
            Phase profile {φ_0, …, φ_8} defining the matrix P_ij = e^{i(φ_i - φ_j)}.

        Notes
        -----
        Validates that:
            • kernel is 9×9 real matrix,
            • phases is length 9,
            • D_F = P * K (entrywise) is constructed.
        """
        self._validate_kernel(kernel)
        self._validate_phases(phases)

        self.kernel: np.ndarray = kernel.astype(float)
        self.phases: np.ndarray = np.array(phases, dtype=float)

        # Build P_ij = exp(i(φ_i - φ_j))
        self.phase_matrix: np.ndarray = self._build_phase_matrix()

        # Finite Dirac operator D_F = P ∘ K
        self.D_F: np.ndarray = self.phase_matrix * self.kernel

        # Hermiticity check for safety (not required every time)
        self._validate_hermiticity(self.D_F)

    # ------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------
    @staticmethod
    def _validate_kernel(kernel: np.ndarray) -> None:
        if not isinstance(kernel, np.ndarray):
            raise TypeError("kernel must be a numpy ndarray.")
        if kernel.shape != (9, 9):
            raise ValueError("kernel must be a 9×9 matrix.")
        if not np.isrealobj(kernel):
            raise ValueError("kernel must be real-valued (forbidden distances have zeros).")

    @staticmethod
    def _validate_phases(phases: Sequence[float]) -> None:
        if len(phases) != 9:
            raise ValueError("phases must be a sequence of 9 real values.")

    @staticmethod
    def _validate_hermiticity(M: np.ndarray) -> None:
        """Ensure finite Dirac operator is Hermitian: M = M† to machine precision."""
        if not np.allclose(M, M.conjugate().T, atol=1e-12):
            raise ValueError("Finite Dirac operator D_F must be Hermitian.")

    # ------------------------------------------------------------
    # Phase Matrix
    # ------------------------------------------------------------
    def _build_phase_matrix(self) -> np.ndarray:
        """
        Compute P_ij = exp(i(φ_i - φ_j)).

        Returns
        -------
        np.ndarray
            A unitary 9×9 matrix with |P_ij| = 1.
        """
        phi = self.phases
        return np.exp(1j * (phi[:, None] - phi[None, :]))

    # ------------------------------------------------------------
    # Real Structure
    # ------------------------------------------------------------
    @staticmethod
    def J_F(vec: np.ndarray) -> np.ndarray:
        """
        Finite real structure J_F acting by complex conjugation.

        Parameters
        ----------
        vec : np.ndarray
            A complex vector in H_F.

        Returns
        -------
        np.ndarray
            Componentwise complex conjugate.
        """
        if not isinstance(vec, np.ndarray):
            raise TypeError("J_F expects a numpy array.")
        return np.conjugate(vec)

    # ------------------------------------------------------------
    # Finite Grading
    # ------------------------------------------------------------
    @staticmethod
    def gamma_F(signs: Sequence[int]) -> np.ndarray:
        """
        Finite grading γ_F: diagonal matrix with entries ±1.

        Parameters
        ----------
        signs : sequence of int
            Must be a list of 9 integers, each +1 or −1.

        Returns
        -------
        np.ndarray
            9×9 diagonal matrix diag(signs).

        Raises
        ------
        ValueError if signs does not have length 9 or contains invalid values.
        """
        if len(signs) != 9:
            raise ValueError("gamma_F requires exactly 9 signs.")
        if any(s not in (-1, 1) for s in signs):
            raise ValueError("gamma_F signs must be ±1.")

        return np.diag(np.array(signs, dtype=float))

"""
What I’d do next, concretely (to “fully upgrade things”)
Step A: Add a triad op (using your TriadicCompression)

loads 9×9 D_F

applies S D_F S†

stores 3×3 D into state

Step B: Decide what diag means

If treating as “Dirac-like Hermitian”: diagonalize D directly (what you’re doing now)

If treating as “Yukawa”: diagonalize D D† (left) and/or D† D (right)

This choice will matter when you reintroduce CKM/PMNS extraction.
"""