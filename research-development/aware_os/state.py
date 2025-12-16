from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict
import numpy as np

@dataclass
class StateStore:
    workdir: Path

    @property
    def path(self) -> Path:
        return self.workdir / "state.npz"

    def save(self, **payload: Any) -> None:
        self.workdir.mkdir(parents=True, exist_ok=True)
        np.savez(self.path, **payload)

    def load(self) -> Dict[str, Any]:
        if not self.path.exists():
            raise FileNotFoundError(f"No state found at {self.path}")
        data = np.load(self.path, allow_pickle=True)
        return {k: data[k] for k in data.files}
