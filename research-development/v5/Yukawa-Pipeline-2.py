#!/usr/bin/env python3
"""
cycle_clock_emergence_yukawa_test.py

Production-ready harness (NO CLI): edit CONFIG and click Run.

Tests:
  (1) Emergence of a child cycle clock (size d) from a parent cycle clock (size N),
      via sublattice selection (every s = N/d site), using the step-matched parent clock.
  (2) Preservation of Yukawa-style texture maps on the child sector.
  (3) Optional mismatch demo: compressing the WRONG parent clock (step=1) to show failure.

Design:
  - Efficient: uses eigendecompositions and probe-vector operator tests by default.
  - Robust: complex-safe algebra, strong validation, structured outputs.
  - Stable metrics: uses symmetric relative error to avoid denominator blow-ups.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import scipy.sparse as sp
import scipy.linalg as la

# =============================================================================
# CONFIG (edit and run)
# =============================================================================

@dataclass(frozen=True)
class Config:
    # Sizes
    parent_N: int = 2160
    child_d: int = 360  # must divide parent_N

    # Heat-kernel scales to test
    alphas: Tuple[float, ...] = (0.1, 0.5, 1.0, 2.0, 4.0)

    # Mismatch demo: compress parent step=1 clock (should fail)
    mismatch_demo: bool = True

    # Kernel operator comparison: number of random probe vectors
    kernel_probes: int = 8

    # Yukawa stage
    run_yukawa: bool = True
    yukawa_compare_mode: str = "probe"  # "probe" (fast) or "matrix" (slow)
    yukawa_probes: int = 8             # used if yukawa_compare_mode="probe"

    # Y0 seed types: "distance", "random_sym", "triadic"
    y0_kind: str = "distance"
    kappa: float = 0.24
    seed: int = 12345
    scale: float = 1.0
    triad_phase: bool = False

    # Pi projector: "none", "lowk", "band"
    pi_mode: str = "none"
    pi_lowk: int = 24
    pi_band: Tuple[float, float] = (0.0, 0.5)

    # Tolerances
    eig_abs_tol: float = 1e-10
    eig_rel_tol: float = 1e-10
    leakage_rel_tol: float = 1e-12
    kernel_rel_tol: float = 1e-10
    yukawa_rel_tol: float = 1e-10

    # Output
    print_summary: bool = True
    write_detailed_json: bool = False
    detailed_json_path: str = "cycle_clock_detailed_results.json"

    # Behavior
    raise_on_fail: bool = False

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
# Validation
# =============================================================================

def validate_config(cfg: Config) -> None:
    if cfg.parent_N <= 1 or cfg.child_d <= 1:
        raise ValueError("parent_N and child_d must be > 1.")
    if cfg.parent_N % cfg.child_d != 0:
        raise ValueError(f"child_d must divide parent_N (got {cfg.child_d} | {cfg.parent_N} is false).")

    if len(cfg.alphas) == 0 or any(a <= 0 for a in cfg.alphas):
        raise ValueError("alphas must be non-empty and all positive.")

    if cfg.kernel_probes < 1 or cfg.yukawa_probes < 1:
        raise ValueError("kernel_probes and yukawa_probes must be >= 1.")

    if cfg.yukawa_compare_mode not in ("probe", "matrix"):
        raise ValueError("yukawa_compare_mode must be 'probe' or 'matrix'.")

    if cfg.y0_kind not in ("distance", "random_sym", "triadic"):
        raise ValueError("y0_kind must be 'distance', 'random_sym', or 'triadic'.")

    if cfg.pi_mode not in ("none", "lowk", "band"):
        raise ValueError("pi_mode must be 'none', 'lowk', or 'band'.")

    if cfg.pi_mode == "lowk" and cfg.pi_lowk < 1:
        raise ValueError("pi_lowk must be >= 1 when pi_mode='lowk'.")
    if cfg.pi_mode == "lowk" and cfg.pi_lowk > cfg.child_d:
        raise ValueError("pi_lowk must be <= child_d.")

    if cfg.pi_mode == "band":
        lo, hi = cfg.pi_band
        if hi < lo:
            raise ValueError("pi_band must satisfy hi >= lo.")

    for name, val in [
        ("eig_abs_tol", cfg.eig_abs_tol),
        ("eig_rel_tol", cfg.eig_rel_tol),
        ("leakage_rel_tol", cfg.leakage_rel_tol),
        ("kernel_rel_tol", cfg.kernel_rel_tol),
        ("yukawa_rel_tol", cfg.yukawa_rel_tol),
    ]:
        if val < 0:
            raise ValueError(f"{name} must be non-negative.")


# =============================================================================
# Stable error metric
# =============================================================================

def rel_err_symmetric(y_a: np.ndarray, y_b: np.ndarray, eps: float = 1e-12) -> float:
    """
    Symmetric relative error:
      ||y_a - y_b|| / max(||y_a||, ||y_b||, eps)

    This avoids artificial blow-ups when one side is near zero.
    """
    na = float(np.linalg.norm(y_a))
    nb = float(np.linalg.norm(y_b))
    denom = max(na, nb, eps)
    return float(np.linalg.norm(y_a - y_b) / denom)


# =============================================================================
# Cycle operators on Z_N
# =============================================================================

def make_shift(N: int, step: int = 1, dtype=np.float64) -> sp.csr_matrix:
    """Sparse permutation matrix U such that U|x> = |x + step mod N>."""
    step = step % N
    rows = np.arange(N, dtype=np.int64)
    cols = (rows - step) % N  # (U v)[row] = v[col]
    data = np.ones(N, dtype=dtype)
    return sp.csr_matrix((data, (rows, cols)), shape=(N, N))


def laplacian_from_shift(U: sp.csr_matrix) -> sp.csr_matrix:
    """L = 2I - U - U^* (for real permutation U, U^* = U^T)."""
    N = U.shape[0]
    I = sp.identity(N, format="csr", dtype=U.dtype)
    return (2 * I) - U - U.transpose()


def selection_matrix_every_s(N: int, s: int, dtype=np.float64) -> sp.csr_matrix:
    """Selection matrix S (d x N), d=N/s, picks indices [0, s, 2s, ..., (d-1)s]."""
    if N % s != 0:
        raise ValueError(f"N must be divisible by s (got N={N}, s={s}).")
    d = N // s
    rows = np.arange(d, dtype=np.int64)
    cols = (rows * s) % N
    data = np.ones(d, dtype=dtype)
    return sp.csr_matrix((data, (rows, cols)), shape=(d, N))


# =============================================================================
# Utilities
# =============================================================================

def op_norm_sparse(A: sp.spmatrix) -> float:
    """Operator norm estimate; exact 2-norm for moderate sizes, conservative otherwise."""
    n = A.shape[0]
    if n <= 600:
        return float(np.linalg.norm(A.toarray(), ord=2))
    return float(np.linalg.norm(A.data))


def leakage_ratio(L: sp.csr_matrix, S: sp.csr_matrix) -> float:
    """leakage = ||(I-P) L P|| / ||L||, where P = S^T S."""
    N = L.shape[0]
    I = sp.identity(N, format="csr", dtype=L.dtype)
    P = (S.transpose() @ S).tocsr()
    leak = (I - P) @ (L @ P)
    denom = op_norm_sparse(L)
    num = op_norm_sparse(leak)
    return 0.0 if denom == 0 else (num / denom)


def eig_decomp_sym(L: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Symmetric eigen-decomposition: L = V diag(w) V^H, w ascending."""
    w, V = la.eigh(L)
    return w, V


def compare_eigs(w1: np.ndarray, w2: np.ndarray, abs_tol: float, rel_tol: float) -> Dict[str, Any]:
    abs_err = float(np.max(np.abs(w1 - w2)))
    rel_err = float(np.max(np.abs(w1 - w2) / np.maximum(1.0, np.abs(w2))))
    return {
        "eig_abs_max": abs_err,
        "eig_rel_max": rel_err,
        "eig_abs_ok": abs_err <= abs_tol,
        "eig_rel_ok": rel_err <= rel_tol,
    }


def projector_mask_from_eigs(w: np.ndarray, mode: str, lowk: int, band: Tuple[float, float]) -> np.ndarray:
    """
    Returns p in {0,1}^n such that Pi = V diag(p) V^H (basis-independent selection).
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
        idx = np.where((w >= lo) & (w <= hi))[0]
        if idx.size == 0:
            raise ValueError(f"pi_band {band} selects no eigenmodes.")
    else:
        raise ValueError(f"Unknown pi_mode: {mode}")

    p = np.zeros(n, dtype=np.float64)
    p[idx] = 1.0
    return p


def sample_unit_vector(rng: np.random.Generator, d: int) -> np.ndarray:
    x = rng.normal(size=d)
    nrm = float(np.linalg.norm(x))
    return x / (nrm if nrm > 0 else 1.0)


# =============================================================================
# Heat-kernel application: K(alpha) = exp(-alpha L)
# =============================================================================

def apply_heat_kernel(w: np.ndarray, V: np.ndarray, alpha: float, x: np.ndarray) -> np.ndarray:
    """
    Apply K(alpha) to x using L = V diag(w) V^H:
      Kx = V (exp(-alpha*w) * (V^H x))
    """
    e = np.exp(-alpha * w)
    return V @ (e * (V.conj().T @ x))


# =============================================================================
# Y0 constructors
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


def y0_random_symmetric(d: int, seed: int, scale: float) -> np.ndarray:
    rng = np.random.default_rng(seed)
    A = rng.normal(loc=0.0, scale=scale, size=(d, d))
    return 0.5 * (A + A.T)


def y0_triadic(d: int, kappa: float, phase: bool) -> np.ndarray:
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
# Yukawa operator: Y(alpha) = Pi K Y0 K Pi
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
    Apply Y(alpha) to x without building Y:
      g = p_mask * exp(-alpha*w)
      y = V diag(g) V^H  *  Y0  *  V diag(g) V^H x
    """
    g = p_mask * np.exp(-alpha * w)
    x_hat = V.conj().T @ x
    u = V @ (g * x_hat)
    v = Y0 @ u
    v_hat = V.conj().T @ v
    y = V @ (g * v_hat)
    return y


def build_full_yukawa_matrix(
    w: np.ndarray,
    V: np.ndarray,
    p_mask: np.ndarray,
    Y0: np.ndarray,
    alpha: float,
) -> np.ndarray:
    g = p_mask * np.exp(-alpha * w)
    Y0_hat = V.conj().T @ (Y0 @ V)
    Y_hat = (g[:, None] * Y0_hat) * g[None, :]
    return V @ (Y_hat @ V.conj().T)


# =============================================================================
# Spectral diagnostics (cycle k/phase)
# =============================================================================

def cycle_k_from_lambda(d: int, lam: float) -> float:
    """
    Approximate k in [0, d/2] from cycle Laplacian eigenvalue:
      lam = 2 - 2 cos(2πk/d)
    => cos = 1 - lam/2
    => k = d/(2π) * arccos(1 - lam/2)
    """
    c = 1.0 - 0.5 * float(lam)
    c = max(-1.0, min(1.0, c))
    theta = float(np.arccos(c))
    return float(d * theta / (2.0 * np.pi))


def band_mode_diagnostics(w: np.ndarray, lo: float, hi: float, d: int) -> Dict[str, Any]:
    """
    Given child eigenvalues w and band [lo, hi], report:
      - number of selected eigenmodes
      - eigen-index range in sorted w
      - approximate k-range implied by endpoints
      - phase range in degrees
    """
    w = np.asarray(w, dtype=float)
    idx = np.where((w >= lo) & (w <= hi))[0]
    k_lo = cycle_k_from_lambda(d, lo)
    k_hi = cycle_k_from_lambda(d, hi)
    deg_lo = 360.0 * k_lo / d
    deg_hi = 360.0 * k_hi / d

    return {
        "modes_selected": int(idx.size),
        "eig_index_min": int(idx.min()) if idx.size else None,
        "eig_index_max": int(idx.max()) if idx.size else None,
        "k_range_est": [float(k_lo), float(k_hi)],
        "phase_deg_est": [float(deg_lo), float(deg_hi)],
    }


# =============================================================================
# Tests
# =============================================================================

def run_emergence(cfg: Config) -> Dict[str, Any]:
    N = cfg.parent_N
    d = cfg.child_d
    s = N // d

    # Step-matched parent clock L_N^(s)
    U_parent_s = make_shift(N, step=s)
    L_parent_s = laplacian_from_shift(U_parent_s).tocsr()

    # Selection onto every s-th site
    S = selection_matrix_every_s(N, s)

    # Leakage test (invariance)
    leak = leakage_ratio(L_parent_s, S)

    # Induced child clock (d x d): S L S^T
    L_child_induced = (S @ (L_parent_s @ S.transpose())).toarray()

    # Direct child clock (d x d): step-1 on Z_d
    U_child_1 = make_shift(d, step=1)
    L_child_direct = laplacian_from_shift(U_child_1).toarray()

    # Eigen spectrum compare
    w_ind_sorted = np.sort(np.linalg.eigvalsh(L_child_induced))
    w_dir_sorted = np.sort(np.linalg.eigvalsh(L_child_direct))
    eig_cmp = compare_eigs(w_ind_sorted, w_dir_sorted, cfg.eig_abs_tol, cfg.eig_rel_tol)

    # Eigendecomps for operator tests
    w_dir, V_dir = eig_decomp_sym(L_child_direct)
    w_ind, V_ind = eig_decomp_sym(L_child_induced)

    # Kernel operator compare via probes (symmetric error)
    rng_k = np.random.default_rng(20250101)  # fixed stream for kernel stage
    kernel_compares: List[Dict[str, Any]] = []
    for a in cfg.alphas:
        worst = 0.0
        for _ in range(cfg.kernel_probes):
            x = sample_unit_vector(rng_k, d)
            y_ind = apply_heat_kernel(w_ind, V_ind, a, x)
            y_dir = apply_heat_kernel(w_dir, V_dir, a, x)
            rel = rel_err_symmetric(y_ind, y_dir, eps=1e-12)
            worst = max(worst, rel)
        kernel_compares.append(
            {"alpha": float(a), "probe_rel_max": float(worst), "ok": worst <= cfg.kernel_rel_tol}
        )

    # Mismatch demo
    mismatch = None
    mismatch_child_L = None
    if cfg.mismatch_demo:
        U_parent_1 = make_shift(N, step=1)
        L_parent_1 = laplacian_from_shift(U_parent_1).tocsr()
        leak_1 = leakage_ratio(L_parent_1, S)
        mismatch_child_L = (S @ (L_parent_1 @ S.transpose())).toarray()
        w_bad_sorted = np.sort(np.linalg.eigvalsh(mismatch_child_L))
        bad_eig_cmp = compare_eigs(w_bad_sorted, w_dir_sorted, cfg.eig_abs_tol, cfg.eig_rel_tol)
        mismatch = {"leakage_one_step_parent": float(leak_1), "eig_compare_one_step_parent": bad_eig_cmp}

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
        # carry for yukawa stage
        "eig_child_direct": (w_dir, V_dir),
        "eig_child_induced": (w_ind, V_ind),
        "mismatch_child_L": mismatch_child_L,
    }


def run_yukawa(cfg: Config, emerg: Dict[str, Any]) -> Dict[str, Any]:
    w_dir, V_dir = emerg["eig_child_direct"]
    w_ind, V_ind = emerg["eig_child_induced"]
    d = int(emerg["child_d"])

    Y0 = build_y0(cfg, d)

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

    if cfg.yukawa_compare_mode == "matrix":
        # Full matrix comparisons (slow for large d)
        for a in cfg.alphas:
            Y_dir = build_full_yukawa_matrix(w_dir, V_dir, p_dir, Y0, a)
            Y_ind = build_full_yukawa_matrix(w_ind, V_ind, p_ind, Y0, a)
            num = float(np.linalg.norm(Y_ind - Y_dir, ord="fro"))
            den = float(max(np.linalg.norm(Y_dir, ord="fro"), 1e-300))
            rel = num / den
            ok = rel <= cfg.yukawa_rel_tol
            results["comparisons"].append({"alpha": float(a), "fro_rel": float(rel), "ok": bool(ok)})
            if not ok:
                results["passed"] = False

        if cfg.mismatch_demo and emerg["mismatch_child_L"] is not None:
            w_bad, V_bad = eig_decomp_sym(emerg["mismatch_child_L"])
            p_bad = projector_mask_from_eigs(w_bad, cfg.pi_mode, cfg.pi_lowk, cfg.pi_band)
            mm = []
            for a in cfg.alphas:
                Y_dir = build_full_yukawa_matrix(w_dir, V_dir, p_dir, Y0, a)
                Y_bad = build_full_yukawa_matrix(w_bad, V_bad, p_bad, Y0, a)
                rel = float(np.linalg.norm(Y_bad - Y_dir, ord="fro") / max(np.linalg.norm(Y_dir, ord="fro"), 1e-300))
                mm.append({"alpha": float(a), "fro_rel": rel})
            results["mismatch_demo"] = mm

        return results

    # Probe operator comparisons (fast) — use symmetric error metric
    rng_y = np.random.default_rng(cfg.seed + 999)  # deterministic stream for yukawa stage
    for a in cfg.alphas:
        worst = 0.0
        for _ in range(cfg.yukawa_probes):
            x = sample_unit_vector(rng_y, d)
            y_dir = apply_yukawa_operator(w_dir, V_dir, p_dir, Y0, a, x)
            y_ind = apply_yukawa_operator(w_ind, V_ind, p_ind, Y0, a, x)
            rel = rel_err_symmetric(y_ind, y_dir, eps=1e-12)
            worst = max(worst, rel)

        ok = worst <= cfg.yukawa_rel_tol
        results["comparisons"].append(
            {"alpha": float(a), "probe_rel_max": float(worst), "ok": bool(ok), "probes": int(cfg.yukawa_probes)}
        )
        if not ok:
            results["passed"] = False

    if cfg.mismatch_demo and emerg["mismatch_child_L"] is not None:
        w_bad, V_bad = eig_decomp_sym(emerg["mismatch_child_L"])
        p_bad = projector_mask_from_eigs(w_bad, cfg.pi_mode, cfg.pi_lowk, cfg.pi_band)
        mm = []
        for a in cfg.alphas:
            worst = 0.0
            for _ in range(cfg.yukawa_probes):
                x = sample_unit_vector(rng_y, d)
                y_dir = apply_yukawa_operator(w_dir, V_dir, p_dir, Y0, a, x)
                y_bad = apply_yukawa_operator(w_bad, V_bad, p_bad, Y0, a, x)
                rel = rel_err_symmetric(y_bad, y_dir, eps=1e-12)
                worst = max(worst, rel)
            mm.append({"alpha": float(a), "probe_rel_max": float(worst), "probes": int(cfg.yukawa_probes)})
        results["mismatch_demo"] = mm

    return results


# =============================================================================
# Child candidates (divisors)
# =============================================================================

def divisors(n: int) -> List[int]:
    ds = set()
    for k in range(1, int(np.sqrt(n)) + 1):
        if n % k == 0:
            ds.add(k)
            ds.add(n // k)
    return sorted(ds)


def enumerate_child_candidates(parent_N: int) -> List[int]:
    """Candidate child sizes d are all divisors of parent_N (excluding 1 and parent itself)."""
    return [d for d in divisors(parent_N) if 1 < d < parent_N]


# =============================================================================
# Alignment scans
# =============================================================================

def scan_pi_lowk(cfg: Config, lowk_values: List[int]) -> List[Dict[str, Any]]:
    """
    Scan pi_mode='lowk' across a list of k values and report stability/mismatch.
    Returns list of records sorted by score (best first).
    """
    records: List[Dict[str, Any]] = []
    emerg = run_emergence(cfg)
    if not emerg["passed"]:
        raise RuntimeError("Emergence failed; scan_pi_lowk requires a valid emergent child sector.")

    for k in lowk_values:
        trial = Config(**{**asdict(cfg), "pi_mode": "lowk", "pi_lowk": int(k)})
        yk = run_yukawa(trial, emerg)

        score = float(max(c["probe_rel_max"] for c in yk["comparisons"])) if yk["mode"] == "probe" \
            else float(max(c["fro_rel"] for c in yk["comparisons"]))

        rec = {"pi_mode": "lowk", "pi_lowk": int(k), "score": score, "passed_yukawa": bool(yk["passed"])}

        if trial.mismatch_demo and yk.get("mismatch_demo") is not None:
            if yk["mode"] == "probe":
                rec["mismatch_score"] = float(max(m["probe_rel_max"] for m in yk["mismatch_demo"]))
            else:
                rec["mismatch_score"] = float(max(m["fro_rel"] for m in yk["mismatch_demo"]))

        records.append(rec)

    records.sort(key=lambda r: r["score"])
    return records


def make_band_grid_from_child_spectrum(
    w: np.ndarray,
    num_lo: int = 12,
    num_hi: int = 12,
    min_width: float = 1e-9,
) -> List[Tuple[float, float]]:
    """
    Build a grid of candidate [lo, hi] bands from child eigenvalues w (ascending),
    using quantiles of the spectrum.
    """
    w = np.asarray(w, dtype=float)
    w_sorted = np.sort(w)
    q_lo = np.linspace(0.0, 0.9, num_lo)
    q_hi = np.linspace(0.1, 1.0, num_hi)

    bands: List[Tuple[float, float]] = []
    for a in q_lo:
        lo = float(np.quantile(w_sorted, a))
        for b in q_hi:
            hi = float(np.quantile(w_sorted, b))
            if hi >= lo + min_width:
                bands.append((lo, hi))

    uniq: Dict[Tuple[float, float], Tuple[float, float]] = {}
    for lo, hi in bands:
        key = (round(lo, 8), round(hi, 8))
        uniq[key] = (lo, hi)

    return list(uniq.values())


def scan_pi_band_auto(cfg: Config, num_lo: int = 12, num_hi: int = 12) -> List[Dict[str, Any]]:
    """
    Automatically scan many band projectors based on the DIRECT child spectrum.
    Returns records sorted by (score asc, mismatch_score desc).
    """
    emerg = run_emergence(cfg)
    if not emerg["passed"]:
        raise RuntimeError("Emergence failed; band auto-scan requires a valid emergent child sector.")

    w_dir, _V_dir = emerg["eig_child_direct"]
    bands = make_band_grid_from_child_spectrum(w_dir, num_lo=num_lo, num_hi=num_hi)

    records: List[Dict[str, Any]] = []
    for lo, hi in bands:
        trial = Config(**{**asdict(cfg), "pi_mode": "band", "pi_band": (float(lo), float(hi))})
        try:
            yk = run_yukawa(trial, emerg)
        except ValueError:
            continue

        score = float(max(c["probe_rel_max"] for c in yk["comparisons"])) if yk["mode"] == "probe" \
            else float(max(c["fro_rel"] for c in yk["comparisons"]))

        rec: Dict[str, Any] = {
            "pi_mode": "band",
            "pi_band": [float(lo), float(hi)],
            "score": score,
            "passed_yukawa": bool(yk["passed"]),
            "band_diag": band_mode_diagnostics(w_dir, lo, hi, d=cfg.child_d),
        }

        if trial.mismatch_demo and yk.get("mismatch_demo") is not None:
            if yk["mode"] == "probe":
                rec["mismatch_score"] = float(max(m["probe_rel_max"] for m in yk["mismatch_demo"]))
            else:
                rec["mismatch_score"] = float(max(m["fro_rel"] for m in yk["mismatch_demo"]))

        records.append(rec)

    records.sort(key=lambda r: (r["score"], -r.get("mismatch_score", 0.0)))
    return records


def refine_band_around_lambda(
    cfg: Config,
    lam0: float = 2.0,
    widths: Optional[List[float]] = None,
    min_modes: int = 2,
) -> List[Dict[str, Any]]:
    """
    Refine scan of Pi = band[lam0 - w, lam0 + w] over widths.
    Sort by (score asc, mismatch_score desc, modes_selected asc).
    """
    if widths is None:
        widths = [0.50, 0.30, 0.20, 0.15, 0.10, 0.08, 0.06, 0.05, 0.04,
                  0.03, 0.025, 0.02, 0.015, 0.01, 0.0075, 0.005]

    emerg = run_emergence(cfg)
    if not emerg["passed"]:
        raise RuntimeError("Emergence failed; refine requires a valid emergent child sector.")

    w_dir, _V_dir = emerg["eig_child_direct"]
    d = cfg.child_d

    records: List[Dict[str, Any]] = []
    for hw in widths:
        lo = float(lam0 - hw)
        hi = float(lam0 + hw)
        trial = Config(**{**asdict(cfg), "pi_mode": "band", "pi_band": (lo, hi)})

        try:
            yk = run_yukawa(trial, emerg)
        except ValueError:
            continue

        score = float(max(c["probe_rel_max"] for c in yk["comparisons"])) if yk["mode"] == "probe" \
            else float(max(c["fro_rel"] for c in yk["comparisons"]))

        mismatch_score = None
        if trial.mismatch_demo and yk.get("mismatch_demo") is not None:
            if yk["mode"] == "probe":
                mismatch_score = float(max(m["probe_rel_max"] for m in yk["mismatch_demo"]))
            else:
                mismatch_score = float(max(m["fro_rel"] for m in yk["mismatch_demo"]))

        diag = band_mode_diagnostics(w_dir, lo, hi, d=d)
        if diag["modes_selected"] < min_modes:
            continue

        rec: Dict[str, Any] = {
            "pi_mode": "band",
            "pi_band": [lo, hi],
            "half_width": float(hw),
            "score": float(score),
            "passed_yukawa": bool(yk["passed"]),
            "mismatch_score": mismatch_score,
            "band_diag": diag,
        }
        if mismatch_score is not None:
            rec["mismatch_per_mode"] = float(mismatch_score / max(diag["modes_selected"], 1))

        records.append(rec)

    records.sort(
        key=lambda r: (
            r["score"],
            -(r["mismatch_score"] if r["mismatch_score"] is not None else 0.0),
            r["band_diag"]["modes_selected"],
        )
    )
    return records


def pick_best_discriminant(records: List[Dict[str, Any]], score_floor: float = 1e-12) -> Optional[Dict[str, Any]]:
    """
    Stable discriminator when score can be exactly 0:
      discriminant = mismatch_score / max(score, score_floor)

    Tie-break: higher mismatch_score.
    """
    best = None
    best_disc = -1.0
    best_mismatch = -1.0

    for r in records:
        if "mismatch_score" not in r:
            continue
        score = float(r.get("score", 0.0))
        mismatch = float(r.get("mismatch_score", 0.0))
        disc = mismatch / max(score, score_floor)

        if (disc > best_disc) or (disc == best_disc and mismatch > best_mismatch):
            best_disc = disc
            best_mismatch = mismatch
            best = {**r, "discriminant": float(disc), "score_floor": float(score_floor)}

    return best


# =============================================================================
# Runner
# =============================================================================

def run(cfg: Config) -> Dict[str, Any]:
    validate_config(cfg)

    emerg = run_emergence(cfg)
    yk = run_yukawa(cfg, emerg) if cfg.run_yukawa else None

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

    if cfg.mismatch_demo and emerg["mismatch_demo"] is not None:
        summary["mismatch_leakage_one_step_parent"] = emerg["mismatch_demo"]["leakage_one_step_parent"]
        summary["mismatch_eig_abs_max"] = emerg["mismatch_demo"]["eig_compare_one_step_parent"]["eig_abs_max"]
        summary["mismatch_eig_rel_max"] = emerg["mismatch_demo"]["eig_compare_one_step_parent"]["eig_rel_max"]

    if yk is not None:
        summary["passed_yukawa"] = yk["passed"]
        summary["yukawa_mode"] = yk["mode"]
        summary["y0_kind"] = yk["y0_kind"]
        summary["pi_mode"] = yk["pi_mode"]

        if yk["mode"] == "probe":
            summary["yukawa_probe_rel_max"] = max(c["probe_rel_max"] for c in yk["comparisons"])
        else:
            summary["yukawa_fro_rel_max"] = max(c["fro_rel"] for c in yk["comparisons"])

        if cfg.mismatch_demo and yk.get("mismatch_demo") is not None:
            if yk["mode"] == "probe":
                summary["mismatch_yukawa_probe_rel_max"] = max(m["probe_rel_max"] for m in yk["mismatch_demo"])
            else:
                summary["mismatch_yukawa_fro_rel_max"] = max(m["fro_rel"] for m in yk["mismatch_demo"])

    if cfg.print_summary:
        print(json.dumps(summary, indent=2))

    if cfg.write_detailed_json:
        detailed = {
            "config": asdict(cfg),
            "emergence": {
                k: v for k, v in emerg.items()
                if k not in ("eig_child_direct", "eig_child_induced", "mismatch_child_L")
            },
            "yukawa": yk,
        }
        out_path = Path(cfg.detailed_json_path).expanduser().resolve()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(detailed, indent=2))
        LOG.info("Wrote detailed results to %s", str(out_path))

    passed_all = emerg["passed"] and (True if yk is None else yk["passed"])
    if cfg.raise_on_fail and not passed_all:
        raise RuntimeError("One or more tests failed. See printed summary/detailed JSON.")

    return summary


def main() -> None:
    setup_logging(CONFIG.log_level)

    # 1) Run core test
    summary = run(CONFIG)

    # 2) Alignment scans
    lowk_scan = scan_pi_lowk(CONFIG, lowk_values=[6, 12, 18, 24, 36, 48, 60, 72])
    best_lowk = lowk_scan[0] if lowk_scan else None

    band_scan = scan_pi_band_auto(CONFIG, num_lo=14, num_hi=14)
    best_band = band_scan[0] if band_scan else None

    report = {
        "core_summary": summary,
        "best_pi_lowk": best_lowk,
        "best_pi_band": best_band,
        "top5_pi_lowk": lowk_scan[:5],
        "top5_pi_band": band_scan[:5],
    }

    print("\nALIGNMENT REPORT")
    print(json.dumps(report, indent=2))

    best_lowk_disc = pick_best_discriminant(lowk_scan, score_floor=1e-12)
    best_band_disc = pick_best_discriminant(band_scan, score_floor=1e-12)

    print("\nBEST DISCRIMINANTS")
    print(json.dumps({
        "best_lowk_discriminant": best_lowk_disc,
        "best_band_discriminant": best_band_disc,
    }, indent=2))

    refined = refine_band_around_lambda(CONFIG, lam0=2.0, min_modes=2)
    best_refined = refined[0] if refined else None

    print("\nREFINED QUARTER-TURN BAND")
    print(json.dumps({
        "best_refined_band": best_refined,
        "top5_refined_bands": refined[:5],
    }, indent=2))


if __name__ == "__main__":
    main()
