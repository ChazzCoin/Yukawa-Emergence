from dataclasses import dataclass
from typing import Dict, Any
import numpy as np
from ..core.hilbert import HilbertSpace
from ..core.algebra import Algebra, Representation
from ..core.dirac import DiracOperator
from ..core.triple import SpectralTriple

@dataclass
class V3Backend:
    name: str = "v3"

    def build(self, cfg: Dict[str, Any]) -> SpectralTriple:
        dim = int(cfg.get("dim", 9))
        H = HilbertSpace(dim=dim, label="H_flav")

        # v3-ish: commutative site algebra placeholder
        A = Algebra(name="C^G (stub)", elements={"G": cfg.get("G", f"Z_{dim}")})

        # rep: diagonal action placeholder
        def pi(a: Any) -> np.ndarray:
            return np.eye(dim, dtype=np.complex128)

        rep = Representation(name="diag (stub)", pi=pi)

        # Dirac/kernel: start with a Hermitian random seed or configured matrix
        seed = int(cfg.get("seed", 0))
        rng = np.random.default_rng(seed)
        X = rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim))
        Dm = (X + X.conj().T) / 2.0

        D = DiracOperator(H=H, D=Dm)
        return SpectralTriple(A=A, H=H, rep=rep, D=D)
