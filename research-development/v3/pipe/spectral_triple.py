import numpy as np
from typing import Sequence


class AlignmentSpectralTriple:
    """
    Full Product Spectral Triple:
        (A, H, D, J, Γ)

    combining:
        • (A_geom, H_geom, D_geom, J_geom, γ_geom)
        • (A_F,   H_F,   D_F,   J_F,   γ_F)

    with product rules:
        D = D_geom ⊗ 1_9 + γ_geom ⊗ D_F
        J = J_geom ⊗ J_F
        Γ = γ_geom ⊗ γ_F
    """

    # ------------------------------------------------------------
    # Initialization
    # ------------------------------------------------------------
    def __init__(self, geom, finite, gamma_F=None) -> None:
        """
        Parameters
        ----------
        geom : GeometricTriple
            Instance containing geometric Dirac, grading, real structure.
        finite : FiniteTriple
            Instance containing D_F, J_F, γ_F.
        gamma_F : np.ndarray, optional
            9×9 finite grading. If None, must be provided later.

        Notes
        -----
        gamma_F is not required to build the Dirac operator, but is required for:
            • full spectral-triple structure,
            • constructing the full Γ = γ_geom ⊗ γ_F.
        """
        self.geom = geom
        self.finite = finite

        # Finite grading (if provided)
        self.gamma_F = gamma_F

    # ------------------------------------------------------------
    # Product Dirac
    # ------------------------------------------------------------
    def build_product_dirac(self, modes: Sequence[int]) -> np.ndarray:
        """
        Construct the product Dirac operator:

            D = D_geom ⊗ I_9  +  γ_geom ⊗ D_F.

        Parameters
        ----------
        modes : sequence of int
            Harmonic modes defining the geometric basis.

        Returns
        -------
        np.ndarray
            Full Dirac operator of dimension (len(modes)*9) × (len(modes)*9).

        Raises
        ------
        ValueError for invalid shapes.
        """
        # Geometric pieces
        Dg = self.geom.D_geom(modes)            # (N × N)
        gamma_g = self.geom.gamma_geom(modes)   # (N × N)

        # Finite piece
        DF = self.finite.D_F                    # (9 × 9)

        if DF.shape != (9, 9):
            raise ValueError("Finite Dirac operator D_F must be 9×9.")

        # Kronecker construction
        term1 = np.kron(Dg, np.eye(9))
        term2 = np.kron(gamma_g, DF)

        D = term1 + term2

        # Hermiticity test
        if not np.allclose(D, D.conjugate().T, atol=1e-12):
            raise ValueError("Product Dirac operator is not Hermitian.")

        return D

    # ------------------------------------------------------------
    # Product Grading Γ = γ_geom ⊗ γ_F
    # ------------------------------------------------------------
    def build_product_grading(self, modes: Sequence[int]) -> np.ndarray:
        """
        Construct Γ = γ_geom ⊗ γ_F.

        Parameters
        ----------
        modes : sequence of int

        Returns
        -------
        np.ndarray
            Grading operator.

        Raises
        ------
        RuntimeError if finite grading is not provided.
        """
        if self.gamma_F is None:
            raise RuntimeError(
                "Finite grading γ_F not provided. "
                "Pass gamma_F to constructor or call set_finite_grading()."
            )

        gamma_g = self.geom.gamma_geom(modes)
        return np.kron(gamma_g, self.gamma_F)

    # ------------------------------------------------------------
    # Product Real Structure J = J_geom ⊗ J_F
    # ------------------------------------------------------------
    def J(self, vec: np.ndarray, modes: Sequence[int]) -> np.ndarray:
        """
        Apply the real structure:
            J = J_geom ⊗ J_F

        Parameters
        ----------
        vec : np.ndarray
            Vector in H_geom ⊗ H_F.
        modes : sequence of int
            Harmonic modes.

        Returns
        -------
        np.ndarray
            Componentwise conjugated vector.

        Notes
        -----
        J_geom acts by complex conjugation on the geometric coefficients.
        J_F acts by complex conjugation on finite coefficients.

        Thus the tensor product action reduces to:
            J(vec) = conjugate(vec)
        """
        if not isinstance(vec, np.ndarray):
            raise TypeError("vec must be a numpy array.")

        return np.conjugate(vec)

    # ------------------------------------------------------------
    # Inner Fluctuations: D_A = D + A + J A J^{-1}
    # ------------------------------------------------------------
    def fluctuated_dirac(
        self,
        D: np.ndarray,
        A: np.ndarray
    ) -> np.ndarray:
        """
        Construct inner-fluctuated Dirac operator:

            D_A = D + A + J A J^{-1}.

        Parameters
        ----------
        D : np.ndarray
            Base Dirac operator.
        A : np.ndarray
            Inner-fluctuation field.

        Returns
        -------
        np.ndarray
            The fluctuated Dirac operator.

        Raises
        ------
        ValueError for shape mismatch.
        """
        if D.shape != A.shape:
            raise ValueError("D and A must have the same shape for fluctuation.")

        # J acts as complex conjugation for this representation
        JAJ = A.conjugate().T  # if you're treating J as conjugation in this basis
        DA = D + A + JAJ

        # Hermiticity check: D_A must remain Hermitian
        if not np.allclose(DA, DA.conjugate().T, atol=1e-12):
            raise ValueError("Fluctuated Dirac operator is not Hermitian.")

        return DA

    # ------------------------------------------------------------
    # Helper: Build D² for spectral action
    # ------------------------------------------------------------
    @staticmethod
    def D_squared(D: np.ndarray) -> np.ndarray:
        """
        Compute D², required for spectral action evaluation.

        Parameters
        ----------
        D : np.ndarray

        Returns
        -------
        np.ndarray
            D @ D
        """
        return D @ D

    # ------------------------------------------------------------
    # Optional setter for γ_F
    # ------------------------------------------------------------
    def set_finite_grading(self, gamma_F: np.ndarray) -> None:
        """
        Set the finite grading γ_F after initialization.

        Parameters
        ----------
        gamma_F : np.ndarray
            9×9 diagonal matrix of ±1.

        Raises
        ------
        ValueError if the grading is malformed.
        """
        if gamma_F.shape != (9, 9):
            raise ValueError("gamma_F must be 9×9.")
        if not np.allclose(gamma_F, np.diag(np.diag(gamma_F))):
            raise ValueError("gamma_F must be a diagonal matrix.")
        if not all(el in (-1, 1) for el in np.diag(gamma_F)):
            raise ValueError("gamma_F diagonal entries must be ±1.")

        self.gamma_F = gamma_F
