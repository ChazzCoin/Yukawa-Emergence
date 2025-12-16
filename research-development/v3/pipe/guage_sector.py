import numpy as np
from typing import Optional, Sequence


class GaugeSector:
    """
    Gauge Sector for the Alignment Spectral Triple v2.0

    Based on inner fluctuations of the internal algebra:
        A_F = C ⊕ H ⊕ M_3(C)

    which generates the gauge group:
        U(1)_A × SU(2)_A × SU(3)_A.

    This class provides:
        • Lie algebra generators for U(1), SU(2), SU(3)
        • Constructed gauge fields A_μ in matrix form
        • Validation of Hermiticity (required in spectral triple)
        • Validation of tracelessness (for SU(n))
    """

    # ------------------------------------------------------------
    # Initialization
    # ------------------------------------------------------------
    def __init__(self) -> None:
        """
        Prepare generator sets for SU(2) and SU(3), with identity for U(1).
        """
        self.U1_gen = 1.0  # scalar generator; represented as identity of the embedding space

        # Pauli matrices (SU(2))
        self.SU2_gens = self._load_su2_generators()

        # Gell-Mann matrices (SU(3))
        self.SU3_gens = self._load_su3_generators()

    # ------------------------------------------------------------
    # SU(2) generators
    # ------------------------------------------------------------
    @staticmethod
    def _load_su2_generators() -> Sequence[np.ndarray]:
        """
        Return Pauli matrices σ₁, σ₂, σ₃ scaled to be Hermitian generators of su(2).
        """
        σ1 = np.array([[0, 1], [1, 0]], dtype=complex)
        σ2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
        σ3 = np.array([[1, 0], [0, -1]], dtype=complex)
        return [σ1, σ2, σ3]

    # ------------------------------------------------------------
    # SU(3) generators (Gell-Mann)
    # ------------------------------------------------------------
    @staticmethod
    def _load_su3_generators() -> Sequence[np.ndarray]:
        """
        Return the eight Hermitian, traceless Gell-Mann matrices λ1,...,λ8.
        """
        λ1 = np.array([[0,1,0],[1,0,0],[0,0,0]], dtype=complex)
        λ2 = np.array([[0,-1j,0],[1j,0,0],[0,0,0]], dtype=complex)
        λ3 = np.array([[1,0,0],[0,-1,0],[0,0,0]], dtype=complex)
        λ4 = np.array([[0,0,1],[0,0,0],[1,0,0]], dtype=complex)
        λ5 = np.array([[0,0,-1j],[0,0,0],[1j,0,0]], dtype=complex)
        λ6 = np.array([[0,0,0],[0,0,1],[0,1,0]], dtype=complex)
        λ7 = np.array([[0,0,0],[0,0,-1j],[0,1j,0]], dtype=complex)
        λ8 = (1/np.sqrt(3)) * np.array([[1,0,0],[0,1,0],[0,0,-2]], dtype=complex)

        return [λ1, λ2, λ3, λ4, λ5, λ6, λ7, λ8]

    # ------------------------------------------------------------
    # Gauge field constructors
    # ------------------------------------------------------------
    def gauge_field_U1(self, size: int, coeff: float = 0.0) -> np.ndarray:
        """
        Construct a U(1)_A gauge field A = coeff * I_size.

        Parameters
        ----------
        size : int
            Size of the representation matrix (depends on field type).
        coeff : float
            U(1) coupling multiplier.

        Returns
        -------
        np.ndarray
            A Hermitian U(1) gauge field matrix.
        """
        if size <= 0:
            raise ValueError("Matrix size must be positive.")
        return coeff * np.eye(size, dtype=complex)

    def gauge_field_SU2(
        self,
        size: int,
        params: Optional[Sequence[float]] = None
    ) -> np.ndarray:
        """
        Construct an SU(2)_A gauge field:

            A = Σ_k params[k] * σ_k

        Parameters
        ----------
        size : int
            Must be 2 for SU(2) doublet representation.
        params : sequence of length 3
            Coefficients for the σ₁,σ₂,σ₃ generators.

        Returns
        -------
        np.ndarray
            2×2 Hermitian matrix.
        """
        if size != 2:
            raise ValueError("SU(2) gauge fields require size = 2.")

        if params is None:
            params = [0.0, 0.0, 0.0]
        if len(params) != 3:
            raise ValueError("SU(2) gauge field requires 3 parameters.")

        A = sum(params[k] * self.SU2_gens[k] for k in range(3))

        if not np.allclose(A, A.conjugate().T, atol=1e-12):
            raise ValueError("Constructed SU(2) gauge field is not Hermitian.")

        # SU(2) fields need not be traceless individually; generators are traceless.
        return A

    def gauge_field_SU3(
        self,
        size: int,
        params: Optional[Sequence[float]] = None
    ) -> np.ndarray:
        """
        Construct an SU(3)_A gauge field:

            A = Σ_k params[k] * λ_k

        Parameters
        ----------
        size : int
            Must be 3 for SU(3) triplet representation.
        params : sequence of length 8
            Coefficients for the λ₁,…,λ₈ generators.

        Returns
        -------
        np.ndarray
            3×3 Hermitian, traceless gauge field.
        """
        if size != 3:
            raise ValueError("SU(3) gauge fields require size = 3.")

        if params is None:
            params = [0.0] * 8
        if len(params) != 8:
            raise ValueError("SU(3) gauge field requires 8 parameters.")

        A = sum(params[k] * self.SU3_gens[k] for k in range(8))

        # Hermiticity check
        if not np.allclose(A, A.conjugate().T, atol=1e-12):
            raise ValueError("Constructed SU(3) gauge field is not Hermitian.")

        # Tracelessness check (required for su(3))
        if abs(np.trace(A)) > 1e-12:
            raise ValueError("Constructed SU(3) gauge field is not traceless.")

        return A

    # ------------------------------------------------------------
    # Utility
    # ------------------------------------------------------------
    @staticmethod
    def embed_block(A: np.ndarray, block_size: int, total_size: int, position: int) -> np.ndarray:
        """
        Embed a gauge block (like SU(2) or SU(3) field) into a larger matrix
        for fields in combined representations.

        This is essential for building the fluctuation matrix A
        on the full Hilbert space.

        Parameters
        ----------
        A : np.ndarray
            Block matrix to embed.
        block_size : int
            Dimension of the block.
        total_size : int
            Dimension of full matrix.
        position : int
            Starting row/column index where block is placed.

        Returns
        -------
        np.ndarray
            total_size × total_size matrix with A embedded inside.
        """
        M = np.zeros((total_size, total_size), dtype=complex)
        M[position:position+block_size, position:position+block_size] = A
        return M

"""
How to upgrade this cleanly (best practice)
Upgrade 1: Make factor-aware embedding the default

Right now embed_block() is positional indexing (fragile and easy to misuse). In product triples you really want Kronecker-factor embedding:

Gauge acts on the internal/gauge factor

Identity on geometry and flavor factors (unless explicitly intended)

Add utilities like:

def kron_embed(A_block, dims, on):
    dims: list like [dim_geom, dim_internal, dim_flavor]
    on: index of factor to apply A_block to
    Returns: full operator on tensor product
    


This prevents “oops I gauged flavor” mistakes.

Upgrade 2: Split “toy gauge fields” from “inner fluctuations”

Keep your current class but rename intent:

ToyGaugeSector (manual A matrices)

InnerFluctuationBuilder (builds A from algebra elements and commutators)

Then your pipeline can choose:

demo mode: toy

structural mode: inner-fluctuation

Upgrade 3: enforce commutant-boundary safety

If you’re in your v5-style “textures in commutant” mode, you can add a gate:

if something is acting nontrivially on the flavor factor as a “gauge field”, reject it unless explicitly allowed

This is literally one check once you represent operators in factorized form.
"""