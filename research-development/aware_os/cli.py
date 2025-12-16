from pathlib import Path
import json
import numpy as np
import typer
from rich import print

from .config import load_config
from .state import StateStore
from .backends.base import get_backend
from .ops.anchor import anchor as op_anchor
from .ops.schur import schur_complement
from .ops.diag import diagonalize_hermitian
from .ops.score import spectral_action_gaussian

app = typer.Typer(add_completion=False)

@app.command()
def build(
    backend: str = typer.Option("v5", help="Backend: v3 or v5"),
    config: Path = typer.Option(Path("config.yml"), help="YAML config path"),
    workdir: Path = typer.Option(Path(".alignncg"), help="Workspace dir"),
):
    cfg = load_config(config)
    be = get_backend(backend)
    triple = be.build(cfg)
    checks = triple.check()

    store = StateStore(workdir)
    store.save(
        backend=np.array([be.name]),
        cfg=np.array([json.dumps(cfg)]),
        D=triple.D.D,
        dim=np.array([triple.H.dim]),
        checks=np.array([json.dumps(checks)]),
    )
    print({"status": "ok", "built": be.name, "checks": checks, "state": str(store.path)})

@app.command()
def anchor(
    alpha: float = typer.Option(..., help="Filter strength"),
    N: int = typer.Option(2160, help="Anchor label (kept for your 360 story; not used yet)"),
    light_dim: int = typer.Option(3, help="Light subspace dim"),
    workdir: Path = typer.Option(Path(".alignncg"), help="Workspace dir"),
):
    store = StateStore(workdir)
    st = store.load()

    D = st["D"]
    out = op_anchor(D, alpha=alpha, light_dim=light_dim)

    payload = dict(st)  # copy
    payload.update(
        D=out["D_eff"],
        Pi=out["Pi"],
        L=out["L"],
        K=out["K"],
        last=np.array([json.dumps({"op": "anchor", "alpha": alpha, "N": N, "light_dim": light_dim})]),
    )
    store.save(**payload)

    print({"status": "ok", "op": "anchor", "alpha": alpha, "N": N, "light_dim": light_dim})

@app.command()
def schur(
    light_dim: int = typer.Option(3, help="Light dim for Schur complement"),
    workdir: Path = typer.Option(Path(".alignncg"), help="Workspace dir"),
):
    store = StateStore(workdir)
    st = store.load()
    D = st["D"]
    D_eff = schur_complement(D, light_dim=light_dim)

    payload = dict(st)
    payload.update(D=D_eff, last=np.array([json.dumps({"op": "schur", "light_dim": light_dim})]))
    store.save(**payload)

    print({"status":"ok","op":"schur","light_dim":light_dim,"new_dim":int(D_eff.shape[0])})

@app.command()
def diag(
    sector: str = typer.Option("kernel", help="Label only (up/down/leptons/etc later)"),
    workdir: Path = typer.Option(Path(".alignncg"), help="Workspace dir"),
):
    store = StateStore(workdir)
    st = store.load()
    D = st["D"]
    out = diagonalize_hermitian(D)
    w = out["eigenvalues"]

    payload = dict(st)
    payload.update(eigvals=w, last=np.array([json.dumps({"op": "diag", "sector": sector})]))
    store.save(**payload)

    print({"status":"ok","op":"diag","sector":sector,"eigs_preview":w[:min(10,len(w))].tolist()})

@app.command()
def score(
    action: str = typer.Option("spectral", help="spectral (gaussian proxy)"),
    Lambda: float = typer.Option(1.0, help="Scale parameter"),
    workdir: Path = typer.Option(Path(".alignncg"), help="Workspace dir"),
):
    store = StateStore(workdir)
    st = store.load()
    D = st["D"]
    if action != "spectral":
        raise typer.BadParameter("Only action=spectral is implemented in this skeleton.")
    val = spectral_action_gaussian(D, Lambda=Lambda)
    payload = dict(st)
    payload.update(score=np.array([val]),
                   last=np.array([json.dumps({"op": "score", "action": action, "Lambda": Lambda})]))
    store.save(**payload)

    print({"status":"ok","op":"score","action":action,"Lambda":Lambda,"value":val})
