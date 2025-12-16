from __future__ import annotations
from dataclasses import dataclass
from typing import Protocol, Dict, Any
from ..core.triple import SpectralTriple

class Backend(Protocol):
    name: str
    def build(self, cfg: Dict[str, Any]) -> SpectralTriple: ...

def get_backend(name: str) -> Backend:
    name = name.lower().strip()
    if name == "v3":
        from .v3 import V3Backend
        return V3Backend()
    if name == "v5":
        from .v5 import V5Backend
        return V5Backend()
    raise ValueError(f"Unknown backend: {name!r}")
