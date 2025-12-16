from dataclasses import dataclass
import numpy as np
from .hilbert import HilbertSpace

@dataclass
class DiracOperator:
    H: HilbertSpace
    D: np.ndarray  # finite-model matrix

    def __post_init__(self):
        if self.D.shape != (self.H.dim, self.H.dim):
            raise ValueError(f"D shape {self.D.shape} != ({self.H.dim},{self.H.dim})")

    def is_self_adjoint(self, tol: float = 1e-10) -> bool:
        return np.allclose(self.D, self.D.conj().T, atol=tol)

    def spectrum(self) -> np.ndarray:
        # self-adjoint => real eigvals; use eigh for stability
        w = np.linalg.eigvalsh(self.D)
        return w
