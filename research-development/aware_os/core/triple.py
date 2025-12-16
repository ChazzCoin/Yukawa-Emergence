from dataclasses import dataclass
from .algebra import Algebra, Representation
from .hilbert import HilbertSpace
from .dirac import DiracOperator

@dataclass
class SpectralTriple:
    A: Algebra
    H: HilbertSpace
    rep: Representation
    D: DiracOperator

    def check(self) -> dict:
        return {
            "D_self_adjoint": self.D.is_self_adjoint(),
            "dim": self.H.dim,
            "algebra": self.A.name,
            "rep": self.rep.name,
        }
