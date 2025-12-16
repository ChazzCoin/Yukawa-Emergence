
import numpy as np
from dataclasses import dataclass
from typing import Literal, Optional, Dict, Any

Mode = Literal["gaussian", "moments"]

@dataclass
class SpectralActionV3:
    """
    Numerically-stable spectral action / score for development (v5-friendly).

    - gaussian:  Tr exp(-(D^2)/Lambda^2) = sum_i exp(-(lambda_i^2)/Lambda^2)
      (bounded in [0, dim], smooth, great for tuning)

    - moments:   uses eigenvalue moments as a proxy:
          a0 = n
          a2 = sum lambda^2
          a4 = sum lambda^4
      with normalization and optional baseline subtraction so it’s informative.
    """
    Lambda: float = 1.0
    mode: Mode = "gaussian"

    # Moments weights (only used in mode="moments")
    f4: float = 1.0
    f2: float = 1.0
    f0: float = 1.0

    # If True, remove constant term f4*Lambda^4*a0 from Srel (recommended)
    subtract_a0: bool = True

    # Normalize moments by dimension to reduce scaling with matrix size
    normalize_by_dim: bool = True

    # Eigen computation
    hermitian_tol: float = 1e-10
    eig_preview: int = 8  # store a few eigenvalues for debugging

    def _eigvals(self, D: np.ndarray) -> np.ndarray:
        if D.ndim != 2 or D.shape[0] != D.shape[1]:
            raise ValueError("D must be a square matrix.")
        if not np.allclose(D, D.conjugate().T, atol=self.hermitian_tol):
            raise ValueError("D must be Hermitian (within tolerance) for stable eigvalsh().")
        return np.linalg.eigvalsh(D)

    def evaluate(self, D: np.ndarray, psi: Optional[np.ndarray] = None) -> Dict[str, Any]:
        """
        Returns a dictionary with:
          - a0,a2,a4 (moments)
          - S0 (raw bosonic score)
          - Srel (baseline-subtracted / normalized score used for tuning)
          - Sferm (if psi provided)
          - Stotal (if psi provided)
          - eigvals_preview
        """
        lam = self._eigvals(D)
        n = lam.size
        L = float(self.Lambda)
        if L <= 0:
            raise ValueError("Lambda must be positive.")

        # moments from eigenvalues (stable)
        a0 = float(n)
        a2 = float(np.sum(lam * lam))
        a4 = float(np.sum((lam * lam) ** 2))

        # optional normalization
        norm = float(n) if (self.normalize_by_dim and n > 0) else 1.0
        a2n = a2 / norm
        a4n = a4 / norm

        out: Dict[str, Any] = {
            "mode": self.mode,
            "Lambda": L,
            "dim": n,
            "a0": a0,
            "a2": a2,
            "a4": a4,
            "a2_norm": a2n,
            "a4_norm": a4n,
            "eigvals_preview": lam[: min(self.eig_preview, n)].tolist(),
        }

        # Bosonic score
        if self.mode == "gaussian":
            # S0 is in [0, n]
            x = (lam * lam) / (L * L)
            S0 = float(np.sum(np.exp(-x)))
            out["S0"] = S0
            out["Srel"] = S0  # gaussian already baseline-free and comparable

        elif self.mode == "moments":
            S0 = float(self.f4 * (L ** 4) * a0 + self.f2 * (L ** 2) * a2n + self.f0 * a4n)
            out["S0"] = S0
            Srel = S0
            if self.subtract_a0:
                Srel -= float(self.f4 * (L ** 4) * a0)
            out["Srel"] = float(Srel)

        else:
            raise ValueError(f"Unknown mode: {self.mode}")

        # Fermionic score (optional)
        if psi is not None:
            psi = np.asarray(psi)
            if psi.ndim != 1 or psi.shape[0] != n:
                raise ValueError("psi must be a vector with the same dimension as D.")
            Sferm = np.vdot(psi, D @ psi)
            out["Sferm"] = float(np.real(Sferm))
            out["Stotal"] = float(out["Srel"] + out["Sferm"])

        return out

class SpectralAction:
    """
    Connes–Chamseddine Spectral Action for the Alignment Spectral Triple v2.0.

    Implements the truncated heat-kernel expansion:

        S_bos(D_A) = f4 * Λ^4 * a0
                    + f2 * Λ^2 * a2
                    + f0 * a4

    where {a0, a2, a4} are the first Seeley–DeWitt coefficients of D_A^2.

    The fermionic part of the action is:

        S_ferm(Ψ, D_A) = <Ψ, D_A Ψ>.
    """

    # ------------------------------------------------------------
    # Initialization
    # ------------------------------------------------------------
    def __init__(
        self,
        cutoff: float = 1e14,
        f4: float = 1.0,
        f2: float = 1.0,
        f0: float = 1.0,
    ) -> None:
        """
        Parameters
        ----------
        cutoff : float
            The energy scale Λ.
        f4, f2, f0 : float
            Cutoff-function moments f_k = ∫₀^∞ f(u) u^{k/2 - 1} du.
            Defaults are placeholders for user-controlled physics inputs.
        """
        if cutoff <= 0:
            raise ValueError("cutoff Λ must be positive.")

        self.Lambda: float = cutoff
        self.f4: float = f4
        self.f2: float = f2
        self.f0: float = f0

    # ------------------------------------------------------------
    # D²
    # ------------------------------------------------------------
    @staticmethod
    def D_squared(D: np.ndarray) -> np.ndarray:
        """
        Compute D² = D @ D.

        Parameters
        ----------
        D : np.ndarray
            Hermitian matrix.

        Returns
        -------
        np.ndarray
            The squared operator.

        Raises
        ------
        ValueError if D is not square.
        """
        if D.ndim != 2 or D.shape[0] != D.shape[1]:
            raise ValueError("D must be a square matrix.")
        return D @ D

    # ------------------------------------------------------------
    # Heat-kernel coefficients
    # ------------------------------------------------------------
    @staticmethod
    def heat_kernel_coeffs(D2: np.ndarray) -> tuple[float, float, float]:
        """
        Compute approximate heat-kernel coefficients:

            a0 = Tr(1)
            a2 = Tr(D²)
            a4 = Tr(D⁴)

        Parameters
        ----------
        D2 : np.ndarray
            D², must be square.

        Returns
        -------
        (a0, a2, a4)
        """
        if D2.ndim != 2 or D2.shape[0] != D2.shape[1]:
            raise ValueError("D2 must be a square matrix.")

        n = D2.shape[0]
        I = np.eye(n)

        a0 = np.trace(I)
        a2 = np.trace(D2)
        a4 = np.trace(D2 @ D2)

        return float(a0), float(a2), float(a4)

    # ------------------------------------------------------------
    # Bosonic spectral action
    # ------------------------------------------------------------
    def bosonic_action(self, D: np.ndarray) -> float:
        """
        Evaluate the bosonic spectral action:

            S_bos = f4 * Λ^4 * a0
                  + f2 * Λ^2 * a2
                  + f0 * a4

        Parameters
        ----------
        D : np.ndarray
            Dirac operator (Hermitian).

        Returns
        -------
        float
            Value of the truncated spectral action.
        """
        if not np.allclose(D, D.conjugate().T, atol=1e-12):
            raise ValueError("Dirac operator D must be Hermitian.")

        D2 = self.D_squared(D)
        a0, a2, a4 = self.heat_kernel_coeffs(D2)

        L = self.Lambda
        S = (
            self.f4 * L**4 * a0 +
            self.f2 * L**2 * a2 +
            self.f0 * a4
        )
        return float(np.real(S))

    # ------------------------------------------------------------
    # Fermionic action  <Ψ, D_A Ψ>
    # ------------------------------------------------------------
    @staticmethod
    def fermionic_action(psi: np.ndarray, D: np.ndarray) -> float:
        """
        Compute the fermionic spectral action:

            S_ferm = <Ψ, D Ψ> = Ψ† D Ψ

        Parameters
        ----------
        psi : np.ndarray
            Fermion Hilbert-space vector.
        D : np.ndarray
            Dirac operator.

        Returns
        -------
        float
            Real part of Ψ† D Ψ.
        """
        if psi.ndim != 1:
            raise ValueError("psi must be a 1D vector.")
        if D.ndim != 2:
            raise ValueError("D must be a matrix.")
        if D.shape[0] != D.shape[1]:
            raise ValueError("D must be square.")
        if psi.shape[0] != D.shape[0]:
            raise ValueError("Dimension mismatch between psi and D.")

        val = np.vdot(psi, D @ psi)
        return float(np.real(val))

    # ------------------------------------------------------------
    # Full spectral action
    # ------------------------------------------------------------
    def total_action(
        self,
        D: np.ndarray,
        psi: Optional[np.ndarray] = None
    ) -> float:
        """
        Compute the total spectral action:

            S = S_bos + S_ferm

        Parameters
        ----------
        D : np.ndarray
            Dirac operator.
        psi : np.ndarray or None
            Fermion field. If None, fermionic term omitted.

        Returns
        -------
        float
            Total spectral action value.
        """
        S_bos = self.bosonic_action(D)

        if psi is None:
            return S_bos

        S_ferm = self.fermionic_action(psi, D)
        return S_bos + S_ferm
