#!/usr/bin/env python3
"""
cycle_clock_emergence_yukawa_test.py

Copy/paste into JetBrains and click Run.

What it does
------------
1) Builds a cycle-based parent clock on Z_N (a step-matched Laplacian using shift^s).
2) Embeds a child cycle Z_d (d | N) by selecting every s = N/d site.
3) Verifies:
   - leakage/invariance of the child subspace inside the parent
   - spectrum match: induced child clock vs direct child clock
   - heat-kernel match: exp(-alpha L) operator match (probe-vector based, fast)
4) Optionally tests a Yukawa-style texture map preservation:
     Y(alpha) = Pi * K(alpha) * Y0 * K(alpha) * Pi
   comparing induced-vs-direct on probe vectors (fast) or via full matrices (slow).
5) Optionally shows failure when using the WRONG parent generator (parent step=1).

No CLI. Edit CONFIG below and run.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import scipy.sparse as sp
import scipy.linalg as la

# =============================================================================
# CONFIG (edit these and run)
# =============================================================================

@dataclass(frozen=True)
class Config:
    # Parent/child sizes
    parent_N: int = 2160
    child_d: int = 360  # must divide parent_N

    # Heat-kernel scales to test
    alphas: Tuple[float, ...] = (0.1, 0.5, 1.0, 2.0, 4.0)

    # Run mismatch demo (compressing parent step=1 clock)
    mismatch_demo: bool = True

    # Kernel comparison method (fast): number of random probe vectors
    kernel_probes: int = 8

    # Yukawa stage
    run_yukawa: bool = True
    yukawa_compare_mode: str = "probe"   # "probe" (fast) or "matrix" (slow)
    yukawa_probes: int = 8               # used if yukawa_compare_mode="probe"

    # Y0 seed choices: "distance", "random_sym", "triadic"
    y0_kind: str = "distance"
    kappa: float = 0.24                  # used for distance/triadic
    seed: int = 12345                    # RNG seed for random_sym AND probe vectors
    scale: float = 1.0                   # used for random_sym
    triad_phase: bool = False            # triadic complex phase wheel

    # Pi projector choices: "none", "lowk", "band"
    pi_mode: str = "none"
    pi_lowk: int = 24
    pi_band: Tuple[float, float] = (0.0, 0.5)

    # Tolerances
    eig_abs_tol: float = 1e-10
    eig_rel_tol: float = 1e-10
    leakage_rel_tol: float = 1e-12
    kernel_rel_tol: float = 1e-10        # probe-relative tolerance for kernel operator
    yukawa_rel_tol: float = 1e-10        # probe-relative or Frobenius-relative (if matrix mode)

    # Output
    write_detailed_json: bool = False
    detailed_json_path: str = "cycle_clock_results.json"

    # Logging
    log_level: int = logging.INFO


CONFIG = Config()

# =============================================================================
# Logging
# =============================================================================

LOG = logging.getLogger("cycle_clock_yukawa_test")


def setup_logging(level: int) -> None:
    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


# =============================================================================
# Cycle operators on Z_N
# =============================================================================

def make_shift(N: int, step: int = 1, dtype=np.float64) -> sp.csr_matrix:
    """
    Sparse permutation matrix U such that:
        U |x> = |x + step (mod N)>
    In vector form, (U v)[row] = v[col] where col = row - step (mod N).
    """
    step = step % N
    rows = np.arange(N, dtype=np.int64)
    cols = (rows - step) % N
    data = np.ones(N, dtype=dtype)
    return sp.csr_matrix((data, (rows, cols)), shape=(N, N))


def laplacian_from_shift(U: sp.csr_matrix) -> sp.csr_matrix:
    """
    Step-Laplacian clock:
        L = 2I - U - U*
    For real permutation U, U* = U^T.
    """
    N = U.shape[0]
    I = sp.identity(N, format="csr", dtype=U.dtype)
    return (2 * I) - U - U.transpose()


def selection_matrix_every_s(N: int, s: int, dtype=np.float64) -> sp.csr_matrix:
    """
    Selection matrix S (d x N), d = N/s, that picks indices:
        0, s, 2s, ..., (d-1)s
    So (S v)[j] = v[j*s].
    """
    if N % s != 0:
        raise ValueError(f"Expected child_d | parent_N. Got N={N}, s={s} (N%s != 0).")
    d = N // s
    rows = np.arange(d, dtype=np.int64)
    cols = (rows * s) % N
    data = np.ones(d, dtype=dtype)
    return sp.csr_matrix((data, (rows, cols)), shape=(d, N))


# =============================================================================
# Utilities: norms, leakage, eigendecomp, projectors
# =============================================================================

def fro_norm(A: np.ndarray) -> float:
    return float(np.linalg.norm(A, ord="fro"))


def op_norm_sparse(A: sp.spmatrix) -> float:
    """
    Operator norm estimate. Exact 2-norm for moderate sizes; conservative otherwise.
    """
    n = A.shape[0]
    if n <= 600:
        return float(np.linalg.norm(A.toarray(), ord=2))
    return float(np.linalg.norm(A.data))


def leakage_ratio(L: sp.csr_matrix, S: sp.csr_matrix) -> float:
    """
    leakage = ||(I-P) L P|| / ||L|| where P = S^T S.
    Zero means the selected child subspace is invariant under L.
    """
    N = L.shape[0]
    I = sp.identity(N, format="csr", dtype=L.dtype)
    P = (S.transpose() @ S).tocsr()
    leak_op = (I - P) @ (L @ P)
    denom = op_norm_sparse(L)
    num = op_norm_sparse(leak_op)
    return 0.0 if denom == 0 else (num / denom)


def eig_decomp_sym(L: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    For symmetric real L, returns (w, V) with L = V diag(w) V^T, w ascending.
    """
    w, V = la.eigh(L)
    return w, V


def compare_eigs(w1: np.ndarray, w2: np.ndarray, abs_tol: float, rel_tol: float) -> Dict[str, Any]:
    if w1.shape != w2.shape:
        raise ValueError(f"Eigen arrays mismatch: {w1.shape} vs {w2.shape}")
    abs_err = float(np.max(np.abs(w1 - w2)))
    rel_err = float(np.max(np.abs(w1 - w2) / np.maximum(1.0, np.abs(w2))))
    return {
        "eig_abs_max": abs_err,
        "eig_rel_max": rel_err,
        "eig_abs_ok": abs_err <= abs_tol,
        "eig_rel_ok": rel_err <= rel_tol,
    }


def projector_mask_from_eigs(
    w: np.ndarray,
    mode: str,
    lowk: int,
    band: Tuple[float, float],
) -> np.ndarray:
    """
    Returns p in {0,1}^n that defines Pi = V diag(p) V^T (basis-independent).
    """
    n = w.shape[0]
    if mode == "none":
        return np.ones(n, dtype=np.float64)

    if mode == "lowk":
        if not (1 <= lowk <= n):
            raise ValueError(f"pi_lowk must be in [1,{n}] (got {lowk}).")
        idx = np.arange(lowk, dtype=np.int64)
    elif mode == "band":
        lo, hi = float(band[0]), float(band[1])
        if hi < lo:
            raise ValueError(f"pi_band must satisfy hi >= lo (got {band}).")
        idx = np.where((w >= lo) & (w <= hi))[0]
        if idx.size == 0:
            raise ValueError(f"pi_band {band} selects no eigenmodes.")
    else:
        raise ValueError(f"Unknown pi_mode: {mode}")

    p = np.zeros(n, dtype=np.float64)
    p[idx] = 1.0
    return p


# =============================================================================
# Fast heat-kernel application (no dense expm)
# =============================================================================

def apply_heat_kernel(w: np.ndarray, V: np.ndarray, alpha: float, x: np.ndarray) -> np.ndarray:
    """
    y = exp(-alpha L) x, with L=V diag(w) V^T.
    """
    e = np.exp(-alpha * w)
    return V @ (e * (V.T @ x))


# =============================================================================
# Y0 constructors (seed coupling kernels)
# =============================================================================

def cyclic_geodesic_dist(i: int, j: int, N: int) -> int:
    d = abs(i - j)
    return min(d, N - d)


def y0_distance_decay(d: int, kappa: float) -> np.ndarray:
    Y0 = np.empty((d, d), dtype=np.float64)
    for i in range(d):
        for j in range(d):
            Y0[i, j] = np.exp(-kappa * cyclic_geodesic_dist(i, j, d))
    return Y0


def y0_random_symmetric(d: int, seed: int, scale: float = 1.0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    A = rng.normal(loc=0.0, scale=scale, size=(d, d))
    return 0.5 * (A + A.T)


def y0_triadic(d: int, kappa: float, phase: bool = False) -> np.ndarray:
    """
    Triadic coupling based on residue classes mod 3.
    Magnitudes by dc = (ci - cj) mod 3:
      dc=0 -> 1
      dc=1 -> kappa
      dc=2 -> kappa^2
    If phase=True, multiply by a mod-3 phase wheel.
    """
    if d % 3 != 0:
        LOG.warning("triadic Y0 is cleanest when d divisible by 3 (got d=%d).", d)

    mags = np.array([1.0, kappa, kappa**2], dtype=np.float64)

    if phase:
        ph = np.array([0.0, 2*np.pi/3, 4*np.pi/3], dtype=np.float64)
        Y0 = np.empty((d, d), dtype=np.complex128)
        for i in range(d):
            ci = i % 3
            for j in range(d):
                cj = j % 3
                dc = (ci - cj) % 3
                Y0[i, j] = mags[dc] * np.exp(1j * (ph[ci] - ph[cj]))
        return Y0

    Y0 = np.empty((d, d), dtype=np.float64)
    for i in range(d):
        ci = i % 3
        for j in range(d):
            cj = j % 3
            dc = (ci - cj) % 3
            Y0[i, j] = mags[dc]
    return Y0


def build_y0(cfg: Config, d: int) -> np.ndarray:
    if cfg.y0_kind == "distance":
        return y0_distance_decay(d, cfg.kappa)
    if cfg.y0_kind == "random_sym":
        return y0_random_symmetric(d, cfg.seed, cfg.scale)
    if cfg.y0_kind == "triadic":
        return y0_triadic(d, cfg.kappa, cfg.triad_phase)
    raise ValueError(f"Unknown y0_kind: {cfg.y0_kind}")


# =============================================================================
# Yukawa operator: apply to vector without forming Y
# =============================================================================

def apply_yukawa_operator(
    w: np.ndarray,
    V: np.ndarray,
    p_mask: np.ndarray,
    Y0: np.ndarray,
    alpha: float,
    x: np.ndarray,
) -> np.ndarray:
    """
    Apply Y(alpha) = Pi K Y0 K Pi to x, where:
      K = exp(-alpha L)
      Pi = spectral projector of L with mask p_mask

    Since Pi and K are spectral functions of L, they commute:
      G(alpha) := Pi K = V diag(g) V^T, g = p_mask * exp(-alpha w)

    Then:
      y = G Y0 G x
    """
    g = p_mask * np.exp(-alpha * w)
    x_hat = V.T @ x
    u = V @ (g * x_hat)
    v = Y0 @ u
    v_hat = V.T @ v
    y = V @ (g * v_hat)
    return y


def build_full_yukawa_matrix(
    w: np.ndarray,
    V: np.ndarray,
    p_mask: np.ndarray,
    Y0: np.ndarray,
    alpha: float,
) -> np.ndarray:
    """
    Builds dense Y(alpha) using eigenbasis (O(d^3)).
    Only use for small d or explicit Frobenius comparisons.
    """
    g = p_mask * np.exp(-alpha * w)
    Y0_hat = V.T @ (Y0 @ V)
    Y_hat = (g[:, None] * Y0_hat) * g[None, :]
    return V @ (Y_hat @ V.T)


# =============================================================================
# Core tests
# =============================================================================

def emergence_tests(cfg: Config) -> Dict[str, Any]:
    N = cfg.parent_N
    d = cfg.child_d
    if N % d != 0:
        raise ValueError(f"child_d must divide parent_N (got N={N}, d={d}).")
    s = N // d

    # Step-matched parent clock: uses shift^s
    U_parent_s = make_shift(N, step=s)
    L_parent_s = laplacian_from_shift(U_parent_s).tocsr()

    # Sublattice selection
    S = selection_matrix_every_s(N, s)

    # Leakage test (invariance)
    leak = leakage_ratio(L_parent_s, S)

    # Induced child clock
    L_child_induced = (S @ (L_parent_s @ S.transpose())).toarray()

    # Direct child clock (one-step on Z_d)
    U_child_1 = make_shift(d, step=1)
    L_child_direct = laplacian_from_shift(U_child_1).toarray()

    # Spectra compare
    w_ind_sorted = np.sort(np.linalg.eigvalsh(L_child_induced))
    w_dir_sorted = np.sort(np.linalg.eigvalsh(L_child_direct))
    eig_cmp = compare_eigs(w_ind_sorted, w_dir_sorted, cfg.eig_abs_tol, cfg.eig_rel_tol)

    # Eigendecompositions for fast operator application
    w_dir, V_dir = eig_decomp_sym(L_child_direct)
    w_ind, V_ind = eig_decomp_sym(L_child_induced)

    # Kernel (heat) operator compare via probe vectors
    rng = np.random.default_rng(20250101)
    probes = max(int(cfg.kernel_probes), 1)

    def sample_unit() -> np.ndarray:
        x = rng.normal(size=d)
        nrm = float(np.linalg.norm(x))
        return x / (nrm if nrm > 0 else 1.0)

    kernel_compares: List[Dict[str, Any]] = []
    for a in cfg.alphas:
        worst = 0.0
        for _ in range(probes):
            x = sample_unit()
            y_ind = apply_heat_kernel(w_ind, V_ind, a, x)
            y_dir = apply_heat_kernel(w_dir, V_dir, a, x)
            rel = float(np.linalg.norm(y_ind - y_dir) / max(np.linalg.norm(y_dir), 1e-300))
            worst = max(worst, rel)
        kernel_compares.append({"alpha": float(a), "probe_rel_max": float(worst), "ok": worst <= cfg.kernel_rel_tol})

    mismatch = None
    mismatch_child_L = None
    if cfg.mismatch_demo:
        # WRONG parent generator: step=1 (not step-matched to sublattice)
        U_parent_1 = make_shift(N, step=1)
        L_parent_1 = laplacian_from_shift(U_parent_1).tocsr()
        leak_1 = leakage_ratio(L_parent_1, S)
        mismatch_child_L = (S @ (L_parent_1 @ S.transpose())).toarray()

        w_bad_sorted = np.sort(np.linalg.eigvalsh(mismatch_child_L))
        bad_eig_cmp = compare_eigs(w_bad_sorted, w_dir_sorted, cfg.eig_abs_tol, cfg.eig_rel_tol)

        mismatch = {
            "leakage_one_step_parent": float(leak_1),
            "eig_compare_one_step_parent": bad_eig_cmp,
        }

    passed = True
    if leak > cfg.leakage_rel_tol:
        passed = False
    if not (eig_cmp["eig_abs_ok"] and eig_cmp["eig_rel_ok"]):
        passed = False
    if not all(k["ok"] for k in kernel_compares):
        passed = False

    return {
        "passed": bool(passed),
        "parent_N": N,
        "child_d": d,
        "step_s": s,
        "leakage_step_matched": float(leak),
        "eigs_compare": eig_cmp,
        "kernel_compares": kernel_compares,
        "mismatch_demo": mismatch,
        # carry forward
        "eig_child_direct": (w_dir, V_dir),
        "eig_child_induced": (w_ind, V_ind),
        "mismatch_child_L": mismatch_child_L,
    }


def yukawa_tests(cfg: Config, emerg: Dict[str, Any]) -> Dict[str, Any]:
    w_dir, V_dir = emerg["eig_child_direct"]
    w_ind, V_ind = emerg["eig_child_induced"]
    d = int(emerg["child_d"])

    # Build Y0 on child
    Y0 = build_y0(cfg, d)

    # Spectral projector masks (basis-independent selection by eigenvalues)
    p_dir = projector_mask_from_eigs(w_dir, cfg.pi_mode, cfg.pi_lowk, cfg.pi_band)
    p_ind = projector_mask_from_eigs(w_ind, cfg.pi_mode, cfg.pi_lowk, cfg.pi_band)

    results: Dict[str, Any] = {
        "passed": True,
        "mode": cfg.yukawa_compare_mode,
        "y0_kind": cfg.y0_kind,
        "kappa": cfg.kappa,
        "seed": cfg.seed,
        "scale": cfg.scale,
        "triad_phase": cfg.triad_phase,
        "pi_mode": cfg.pi_mode,
        "pi_lowk": cfg.pi_lowk,
        "pi_band": [float(cfg.pi_band[0]), float(cfg.pi_band[1])],
        "comparisons": [],
        "mismatch_demo": None,
    }

    if cfg.yukawa_compare_mode not in ("probe", "matrix"):
        raise ValueError("yukawa_compare_mode must be 'probe' or 'matrix'.")

    if cfg.yukawa_compare_mode == "matrix":
        # Full dense Yukawa comparisons (slow for large d)
        for a in cfg.alphas:
            Y_dir = build_full_yukawa_matrix(w_dir, V_dir, p_dir, Y0, a)
            Y_ind = build_full_yukawa_matrix(w_ind, V_ind, p_ind, Y0, a)
            num = fro_norm(Y_ind - Y_dir)
            den = max(fro_norm(Y_dir), 1e-300)
            rel = float(num / den)
            ok = rel <= cfg.yukawa_rel_tol
            results["comparisons"].append({"alpha": float(a), "fro_rel": rel, "ok": ok})
            if not ok:
                results["passed"] = False

        if cfg.mismatch_demo and emerg["mismatch_child_L"] is not None:
            w_bad, V_bad = eig_decomp_sym(emerg["mismatch_child_L"])
            p_bad = projector_mask_from_eigs(w_bad, cfg.pi_mode, cfg.pi_lowk, cfg.pi_band)
            mm = []
            for a in cfg.alphas:
                Y_dir = build_full_yukawa_matrix(w_dir, V_dir, p_dir, Y0, a)
                Y_bad = build_full_yukawa_matrix(w_bad, V_bad, p_bad, Y0, a)
                rel = float(fro_norm(Y_bad - Y_dir) / max(fro_norm(Y_dir), 1e-300))
                mm.append({"alpha": float(a), "fro_rel": rel})
            results["mismatch_demo"] = mm

        return results

    # Probe-vector operator comparison (fast)
    probes = max(int(cfg.yukawa_probes), 1)
    rng = np.random.default_rng(cfg.seed + 999)

    def sample_unit() -> np.ndarray:
        x = rng.normal(size=d)
        nrm = float(np.linalg.norm(x))
        return x / (nrm if nrm > 0 else 1.0)

    for a in cfg.alphas:
        worst = 0.0
        for _ in range(probes):
            x = sample_unit()
            y_dir = apply_yukawa_operator(w_dir, V_dir, p_dir, Y0, a, x)
            y_ind = apply_yukawa_operator(w_ind, V_ind, p_ind, Y0, a, x)
            rel = float(np.linalg.norm(y_ind - y_dir) / max(np.linalg.norm(y_dir), 1e-300))
            worst = max(worst, rel)
        ok = worst <= cfg.yukawa_rel_tol
        results["comparisons"].append({"alpha": float(a), "probe_rel_max": float(worst), "ok": ok, "probes": probes})
        if not ok:
            results["passed"] = False

    if cfg.mismatch_demo and emerg["mismatch_child_L"] is not None:
        w_bad, V_bad = eig_decomp_sym(emerg["mismatch_child_L"])
        p_bad = projector_mask_from_eigs(w_bad, cfg.pi_mode, cfg.pi_lowk, cfg.pi_band)
        mm = []
        for a in cfg.alphas:
            worst = 0.0
            for _ in range(probes):
                x = sample_unit()
                y_dir = apply_yukawa_operator(w_dir, V_dir, p_dir, Y0, a, x)
                y_bad = apply_yukawa_operator(w_bad, V_bad, p_bad, Y0, a, x)
                rel = float(np.linalg.norm(y_bad - y_dir) / max(np.linalg.norm(y_dir), 1e-300))
                worst = max(worst, rel)
            mm.append({"alpha": float(a), "probe_rel_max": float(worst), "probes": probes})
        results["mismatch_demo"] = mm

    return results


# =============================================================================
# Main
# =============================================================================

def main() -> None:
    setup_logging(CONFIG.log_level)

    emerg = emergence_tests(CONFIG)
    yk = yukawa_tests(CONFIG, emerg) if CONFIG.run_yukawa else None

    summary: Dict[str, Any] = {
        "passed_emergence": emerg["passed"],
        "parent_N": emerg["parent_N"],
        "child_d": emerg["child_d"],
        "step_s": emerg["step_s"],
        "leakage_step_matched": emerg["leakage_step_matched"],
        "eig_abs_max": emerg["eigs_compare"]["eig_abs_max"],
        "eig_rel_max": emerg["eigs_compare"]["eig_rel_max"],
        "kernel_probe_rel_max": max(k["probe_rel_max"] for k in emerg["kernel_compares"]),
    }

    if yk is not None:
        summary["passed_yukawa"] = yk["passed"]
        summary["yukawa_mode"] = yk["mode"]
        summary["y0_kind"] = yk["y0_kind"]
        summary["pi_mode"] = yk["pi_mode"]
        if yk["mode"] == "probe":
            summary["yukawa_probe_rel_max"] = max(c["probe_rel_max"] for c in yk["comparisons"])
        else:
            summary["yukawa_fro_rel_max"] = max(c["fro_rel"] for c in yk["comparisons"])

    print(json.dumps(summary, indent=2))

    if CONFIG.write_detailed_json:
        detailed = {
            "config": CONFIG.__dict__,
            "emergence": {
                k: v for k, v in emerg.items()
                if k not in ("eig_child_direct", "eig_child_induced", "mismatch_child_L")
            },
            "yukawa": yk,
        }
        out_path = Path(CONFIG.detailed_json_path).expanduser().resolve()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(detailed, indent=2))
        LOG.info("Wrote detailed results to %s", str(out_path))


if __name__ == "__main__":
    main()
