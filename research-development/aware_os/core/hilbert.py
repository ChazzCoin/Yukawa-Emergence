from dataclasses import dataclass

@dataclass(frozen=True)
class HilbertSpace:
    dim: int
    label: str = "H"
