from dataclasses import dataclass
from typing import Callable, Any, Dict
import numpy as np

@dataclass(frozen=True)
class Algebra:
    name: str
    elements: Dict[str, Any]  # symbolic payload

@dataclass(frozen=True)
class Representation:
    name: str
    pi: Callable[[Any], np.ndarray]  # maps algebra element -> operator on H
