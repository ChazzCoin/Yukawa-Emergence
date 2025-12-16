from dataclasses import dataclass
from typing import Dict, Any
import numpy as np
from ..core.hilbert import HilbertSpace
from ..core.algebra import Algebra, Representation
from ..core.dirac import DiracOperator
from ..core.triple import SpectralTriple

@dataclass
class V5Backend:
    name: str = "v5"

    def build(self, cfg: Dict[str, Any]) -> SpectralTriple:
        dim = int(cfg.get("dim", 9))
        H = HilbertSpace(dim=dim, label="H_flav")

        # v5-ish: "pure multiplicity" flavor factor — algebra acts trivially on flavor
        A = Algebra(name="A_SM (acts trivially on flav) (stub)", elements={"note": "commutant insertion model"})

        def pi(a: Any) -> np.ndarray:
            # acts as identity on flavor multiplicity
            return np.eye(dim, dtype=np.complex128)

        rep = Representation(name="trivial-on-flavor (stub)", pi=pi)

        # base D seed (later you’ll inject textures via ops.anchor / commutant modules)
        seed = int(cfg.get("seed", 1))
        rng = np.random.default_rng(seed)
        X = rng.normal(size=(dim, dim))
        Dm = (X + X.T) / 2.0

        D = DiracOperator(H=H, D=Dm.astype(np.complex128))
        return SpectralTriple(A=A, H=H, rep=rep, D=D)
