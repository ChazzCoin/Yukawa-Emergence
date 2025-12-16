# main.py
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

# Import your Typer app (this is your existing CLI module)
from cli import build as cli_build
from cli import anchor as cli_anchor
from cli import schur as cli_schur
from cli import diag as cli_diag
from cli import score as cli_score

from F import OS

current_path = OS.get_path(__file__=__file__)
config_path = f"{current_path}/aware_os/config.yml"
# ---------------------------
# Edit these and hit Run
# ---------------------------
@dataclass(frozen=True)
class RunConfig:
    # workspace
    workdir: str = Path("/aware_os")
    config_path: str = Path("/aware_os/config.yml")

    # backend/build
    backend: Literal["v3", "v5"] = "v5"

    # pipeline toggles (turn on/off steps)
    do_build: bool = True
    do_anchor: bool = True
    do_schur: bool = True
    do_diag: bool = True
    do_score: bool = True

    # anchor params
    N: int = 2160
    alpha: float = 0.08
    light_dim: int = 3

    # diag params
    sector: str = "up"

    # score params
    action: str = "spectral"
    Lambda: float = 1.0


def run(cfg: RunConfig) -> None:
    """
    Runs the pipeline using your existing Typer command functions directly.
    This is ideal for PyCharm: edit cfg defaults, click Run.
    """
    # 1) Build
    if cfg.do_build:
        cli_build(
            backend=cfg.backend,
            config=cfg.config_path,
            workdir=cfg.workdir,
        )

    # 2) Anchor
    if cfg.do_anchor:
        cli_anchor(
            alpha=cfg.alpha,
            N=cfg.N,
            light_dim=cfg.light_dim,
            workdir=cfg.workdir,
        )

    # 3) Schur
    if cfg.do_schur:
        cli_schur(
            light_dim=cfg.light_dim,
            workdir=cfg.workdir,
        )

    # 4) Diag
    if cfg.do_diag:
        cli_diag(
            sector=cfg.sector,
            workdir=cfg.workdir,
        )

    # 5) Score
    if cfg.do_score:
        cli_score(
            action=cfg.action,
            Lambda=cfg.Lambda,
            workdir=cfg.workdir,
        )


if __name__ == "__main__":
    cfg = RunConfig(
        # <- optionally override defaults here, or just edit the dataclass defaults above
        # alpha=0.1,
        # backend="v3",
    )
    run(cfg)
