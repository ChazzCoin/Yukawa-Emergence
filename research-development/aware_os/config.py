from pathlib import Path
from typing import Dict, Any
import yaml

def load_config(path: Path="/Users/chazzromeo/rAI/Flavor/aware_os/config.yml") -> Dict[str, Any]:
    if isinstance(path, str):
        path = Path(path)
    if not path.exists():
        raise FileNotFoundError(str(path))
    return yaml.safe_load(path.read_text()) or {}
