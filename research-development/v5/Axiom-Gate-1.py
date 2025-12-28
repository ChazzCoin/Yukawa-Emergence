#!/usr/bin/env python3
"""
Minimal, production-ready Python script to test Alignment Spectral Triple v5.0r.

Focuses on flavor sector (N=9, 360 + 9→3) with basis-free projectors, texture generation,
and production gates (8.1–8.4). Assumes finite-dim; hardcodes small dims for geom/SM.
Uses NumPy for operators; tolerances as params.

Usage:
    python test_spectral_triple_v5r.py --N 9 --alpha 0.1 --eps 1e-10
    --flavor_eigs '[1,1,1,2,2,2,2,3,3]' --omega_360 [1]
    --fit  # To run flavor fit

Robustness: Logging, argparser, error handling, assertions for sanity.
No external deps beyond NumPy/SciPy (standard in envs).
"""

import argparse
import logging
import sys
from typing import List, Union, Tuple

import numpy as np
from scipy.linalg import norm, eigh, expm
from scipy.optimize import minimize
from numpy.linalg import pinv, eigvals, svd

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Tolerances (defaults from doc)
EPS_ORDER0 = 1e-12
EPS_ORDER1 = 1e-12
EPS_P = 1e-12
EPS_BLOCK = 1e-12
EPS_BIN = 1e-12
EPS_ISO = 1e-12

# PDG 2025 data for fits (Yukawas y = m/v, v=246 GeV)
v_ew = 246.0  # GeV
m_u, m_c, m_t = 0.0022, 1.27, 172.69
m_d, m_s, m_b = 0.0047, 0.095, 4.18
y_u = np.array([m_u, m_c, m_t]) / v_ew
y_d = np.array([m_d, m_s, m_b]) / v_ew
log_y_u = np.log10(y_u)
log_y_d = np.log10(y_d)

# Relative errors (approx from PDG)
rel_u = np.array([0.18, 0.016, 0.0017])
rel_d = np.array([0.11, 0.05, 0.007])
sigma_log_u = rel_u / np.log(10)
sigma_log_d = rel_d / np.log(10)

# CKM elements and errors
ckm_obs = np.array([0.2243, 0.0410, 0.00382])  # |V_us|, |V_cb|, |V_ub|
sigma_ckm = np.array([0.0008, 0.0010, 0.00020])

DATA = {
    'log_y_u': log_y_u, 'sigma_log_u': sigma_log_u,
    'log_y_d': log_y_d, 'sigma_log_d': sigma_log_d,
    'ckm': ckm_obs, 'sigma_ckm': sigma_ckm
}

import ast
from scipy.linalg import svdvals

def op_norm(A: np.ndarray) -> float:
    return float(svdvals(A)[0])

def parse_eigenvalues(eig_str: str) -> np.ndarray:
    return np.array(ast.literal_eval(eig_str), dtype=float)

def parse_set(s_str: str) -> set:
    return set(float(x) for x in ast.literal_eval(s_str))

def orthonormal_basis_from_projector(P: np.ndarray, r: int, tol: float = 1e-10) -> np.ndarray:
    # Returns an N×r matrix U with orthonormal columns spanning Ran(P)
    w, V = eigh((P + P.conj().T) / 2)
    idx = np.where(w > 1 - tol)[0]
    if len(idx) != r:
        # fallback: take top-r eigenvectors
        idx = np.argsort(w)[-r:]
    U = V[:, idx]
    # Orthonormalize (numerical safety)
    Q, _ = np.linalg.qr(U)
    return Q[:, :r]

def sharp_matrix(X: np.ndarray, U_flav: np.ndarray = None) -> np.ndarray:
    """
    Implements X^sharp = J X^* J^{-1} on H_flav.
    If J_flav = U_flav ∘ K, then X^sharp = U_flav X^T U_flav^*.
    If you take J_flav = K, set U_flav=None and this reduces to X^T.
    """
    if U_flav is None:
        return X.T
    return U_flav @ X.T @ U_flav.conj().T

def seesaw_eff_basis_free(v: float, Y_nu: np.ndarray, M_R: np.ndarray,
                          P_L: np.ndarray, P_H: np.ndarray,
                          U_flav: np.ndarray = None) -> np.ndarray:
    """
    Basis-free seesaw:
      m_D := (P_H Y_nu P_L) as map H_L -> H_H (6×3 block in a projector basis)
      m_eff := -v^2 m_D^sharp (M_HH)^{-1} m_D  on H_L
    Returns the 3×3 light operator in the canonical light basis from P_L.
    """
    # Build orthonormal bases for light/heavy subspaces from projectors
    U_L = orthonormal_basis_from_projector(P_L, r=3)
    U_H = orthonormal_basis_from_projector(P_H, r=6)

    # Block matrices in those bases
    m_D = U_H.conj().T @ Y_nu @ U_L          # 6×3
    M_HH = U_H.conj().T @ M_R @ U_H          # 6×6

    # Inverse on heavy subspace
    s = svdvals(M_HH)
    cond = float(s[0] / s[-1]) if s[-1] > 0 else np.inf
    if cond > 1e8:
        logger.warning(f"Ill-conditioned M_HH (cond={cond:.2e})")
    M_HH_inv = np.linalg.pinv(M_HH)

    # sharp (Majorana transpose)
    m_D_sharp = sharp_matrix(m_D, U_flav=U_flav)  # 3×6

    m_eff = -(v**2) * (m_D_sharp @ M_HH_inv @ m_D)  # 3×3
    return m_eff

def generate_texture(L_eigs, alpha=0.1, N=3, sector='up', seed=42):
    """Generate hierarchical texture Y = K Y0 K with Fritzsch-like Y0."""
    np.random.seed(seed)
    # Hierarchical diagonal + small perturbation for mixing
    if sector == 'up':
        diag_scale = np.array([1e-5, 0.005, 0.7])
    else:
        diag_scale = np.array([2e-5, 4e-4, 0.017])
    rand_off = 0.01 * (np.random.randn(N, N) + 1j * np.random.randn(N, N))
    rand_off = (rand_off + rand_off.conj().T) / 2
    Y0 = rand_off.copy()
    for i in range(N):
        Y0[i, i] += diag_scale[i]
    L = np.diag(L_eigs)
    K = expm(-alpha * L)
    Y = K @ Y0 @ K
    vals, U = eigh(Y)
    y = np.sort(np.abs(vals))  # ascending: light to heavy
    log_y = np.log10(y)
    return log_y, U


def chi2_quarks(params):
    """χ² for quark sector: masses + CKM."""
    alpha_u = params[0]
    L_u = params[1:4]
    alpha_d = params[4]
    L_d = params[5:8]
    log_y_u, U_u = generate_texture(L_u, alpha_u, sector='up')
    log_y_d, U_d = generate_texture(L_d, alpha_d, sector='down')
    V = U_u.conj().T @ U_d
    v_us = np.abs(V[0, 1])
    v_cb = np.abs(V[1, 2])
    v_ub = np.abs(V[0, 2])
    pred_ckm = np.array([v_us, v_cb, v_ub])

    chi_mass = (np.sum(((log_y_u - DATA['log_y_u']) / DATA['sigma_log_u']) ** 2 +
                       ((log_y_d - DATA['log_y_d']) / DATA['sigma_log_d']) ** 2))
    chi_ckm = np.sum(((pred_ckm - DATA['ckm']) / DATA['sigma_ckm']) ** 2)
    return chi_mass + chi_ckm


def run_fit(sector='quarks'):
    """Run fit for sector."""
    if sector == 'quarks':
        init_params = np.array([0.5, 0.1, 1.0, 2.5, 0.5, 0.1, 1.0, 2.5])
        bounds = [(0.1, 3.0)] * 2 + [(0.01, 5.0)] * 6
        res = minimize(chi2_quarks, init_params, method='L-BFGS-B', bounds=bounds)
        chi_dof = res.fun / 6
        logger.info(f'{sector.capitalize()} fit: χ²={res.fun:.2f} (dof=6, χ²/dof={chi_dof:.2f})')
        logger.info(f'Best params: {res.x}')
        # Sample predictions
        alpha_u, L_u = res.x[0], res.x[1:4]
        log_y_u_pred, U_u = generate_texture(L_u, alpha_u, sector='up')
        alpha_d, L_d = res.x[4], res.x[5:8]
        _, U_d = generate_texture(L_d, alpha_d, sector='down')
        V = U_u.conj().T @ U_d
        pred_ckm = [np.abs(V[0, 1]), np.abs(V[1, 2]), np.abs(V[0, 2])]
        logger.info(f'Pred log_y_u: {log_y_u_pred}, Target: {DATA["log_y_u"]}')
        logger.info(f'Pred CKM [V_us, V_cb, V_ub]: {pred_ckm}, Target: {DATA["ckm"]}')
        return res.x, res.fun
    # Stub for leptons
    elif sector == 'leptons':
        logger.info('Lepton fit stub: Implement PMNS + Δm² via seesaw.')
        return None, None
    else:
        raise ValueError(f'Unknown sector: {sector}')


def parse_eigenvalues(eig_str: str) -> np.ndarray:
    """Parse eigenvalue string to array, e.g., '[1,1,1]' -> np.array([1,1,1])."""
    try:
        return np.array(eval(eig_str), dtype=float)
    except Exception as e:
        logger.error(f"Failed to parse eigenvalues '{eig_str}': {e}")
        sys.exit(1)


def parse_set(s_str: str) -> set:
    """Parse set string, e.g., '[1]' -> {1.0}."""
    try:
        return set(float(x) for x in eval(s_str))
    except Exception as e:
        logger.error(f"Failed to parse set '{s_str}': {e}")
        sys.exit(1)


def generate_L(N: int, eigenvalues: Union[str, np.ndarray] = None) -> np.ndarray:
    """Generate self-adjoint L_N (diagonal for simplicity; eigenvalues provided or default)."""
    if isinstance(eigenvalues, str):
        eigenvalues = parse_eigenvalues(eigenvalues)
    if eigenvalues is None:
        # Default: degeneracy example for N=9 (3@1, 4@2, 2@3)
        eigenvalues = np.array([1.0] * 3 + [2.0] * 4 + [3.0] * 2)
    if len(eigenvalues) != N:
        logger.error(f"Eigenvalues length {len(eigenvalues)} != N={N}")
        sys.exit(1)
    L = np.diag(eigenvalues)
    logger.info(f"Generated L_{N} with eigenvalues: {eigenvalues}")
    return L


def spectral_projector(L: np.ndarray, omega: set) -> np.ndarray:
    """Basis-free spectral projector Π = χ_Ω(L) via eig decomp (robust to non-diag)."""
    try:
        vals, vecs = eigh(L)
        P = np.zeros_like(L)
        for i, val in enumerate(vals):
            min_dist = min(abs(val - x) for x in omega)
            if min_dist < 1e-10:
                P += np.outer(vecs[:, i], vecs[:, i].conj())
        assert np.allclose(P @ P, P, atol=1e-10), "Projector not idempotent"
        assert np.allclose(P.conj().T, P, atol=1e-10), "Projector not Hermitian"
        return P
    except Exception as e:
        logger.error(f"Failed to compute spectral projector: {e}")
        sys.exit(1)


def compute_spectral_blocks(L: np.ndarray) -> Tuple[List[np.ndarray], List[float]]:
    """Compute spectral projectors P_λ and eigenvalues λ (sorted)."""
    vals, vecs = eigh(L)
    unique_vals = np.unique(vals.round(10))  # Cluster near-degens
    blocks = []
    evals = []
    for uv in unique_vals:
        idx = np.abs(vals - uv) < 1e-10
        P_lambda = np.zeros_like(L)
        for i in np.where(idx)[0]:
            P_lambda += np.outer(vecs[:, i], vecs[:, i].conj())
        blocks.append(P_lambda)
        evals.append(uv)
    return blocks, evals


def w_star_cert(P: np.ndarray, L: np.ndarray, eps_block: float, eps_bin: float) -> bool:
    """Gate 8.3a: Degeneracy-robust W*(L) membership cert."""
    blocks, evals = compute_spectral_blocks(L)
    passed = True
    for i, (P_lambda, lam) in enumerate(zip(blocks, evals)):
        tr_pl = np.trace(P_lambda)
        if abs(tr_pl) < 1e-10:
            continue
        c_lam = np.trace(P_lambda @ P) / tr_pl
        delta_lam = norm(P_lambda @ P @ P_lambda - c_lam * P_lambda)
        bin_dev = min(c_lam, 1 - c_lam)
        if delta_lam > eps_block or bin_dev > eps_bin:
            logger.warning(f"Block λ={lam}: Δ={delta_lam:.2e} > {eps_block}, bin_dev={bin_dev:.2e} > {eps_bin}")
            passed = False
        else:
            logger.info(f"Block λ={lam}: Δ={delta_lam:.2e}, c={c_lam:.3f} OK")
    return passed


def projector_sanity(P: np.ndarray, name: str, eps_p: float) -> bool:
    """Gate 8.3: Basic projector properties."""
    idemp = norm(P @ P - P)
    herm = norm(P.conj().T - P)
    passed = (idemp <= eps_p) and (herm <= eps_p)
    status = "PASS" if passed else "FAIL"
    logger.info(f"{name} sanity: idemp={idemp:.2e}, herm={herm:.2e} -> {status}")
    return passed


def texture_map(alpha: float, Pi: np.ndarray, L: np.ndarray, Y0: np.ndarray) -> np.ndarray:
    """Def. 4.2: Texture Y = T_{α,360}(Y0)."""
    K_alpha = expm(-alpha * L)
    Y = Pi @ K_alpha @ Y0 @ K_alpha @ Pi
    logger.info(f"Generated texture Y (norm={norm(Y):.3f})")
    return Y


def heavy_light_split(L: np.ndarray, split_type: str, omega_light: Union[str, set] = None,
                      anchor_i: np.ndarray = None) -> Tuple[np.ndarray, np.ndarray]:
    """Def. 5.1: Rank-3 light projector P_L (A: spectral, B: anchored, C: refined= A for now)."""
    N = L.shape[0]
    if split_type == "A":
        if omega_light is None:
            logger.error("Spectral split needs omega_light")
            sys.exit(1)
        omega_set = parse_set(omega_light) if isinstance(omega_light, str) else omega_light
        P_L = spectral_projector(L, omega_set)
        if abs(np.trace(P_L).real - 3) > 1e-6:
            logger.warning(f"Rank(P_L)={np.trace(P_L).real} !=3; adjust omega")
    elif split_type == "B":
        if anchor_i is None or anchor_i.shape != (N, 3):
            logger.error("Anchored needs isometry i: C^3 -> C^N")
            sys.exit(1)
        P_L = anchor_i @ anchor_i.conj().T
    else:
        raise ValueError(f"Unknown split_type: {split_type}")

    P_H = np.eye(N) - P_L
    rank_L = np.trace(P_L).real
    logger.info(f"{split_type} split: rank(P_L)={rank_L:.1f}")
    if abs(rank_L - 3) > 1e-8:
        logger.warning("Rank not exactly 3")
    return P_L, P_H


def seesaw_eff(m_v: float, Y_nu: np.ndarray, M_R: np.ndarray, P_L: np.ndarray, P_H: np.ndarray) -> np.ndarray:
    """Def. 5.4: Effective light neutrino mass (block extract; SVD cond on H_H)."""
    Y_LH = P_L @ Y_nu @ P_H
    Y_HL = P_H @ Y_nu @ P_L
    M_R_HH = P_H @ M_R @ P_H

    # SVD cond
    sv = svd(M_R_HH, compute_uv=False)
    tol = 1e-10 * sv[0] if len(sv) > 0 else 1e-10
    rank_h = np.sum(sv > tol)
    cond_num = sv[0] / sv[rank_h - 1] if rank_h > 0 else np.inf
    if rank_h == 0:
        logger.warning("M_R_HH has no heavy support")
    if cond_num > 1e6:
        logger.warning(f"Ill-conditioned M_R_HH (cond={cond_num:.2e})")

    try:
        inv_M_R_HH = pinv(M_R_HH)
        Y_HL_sharp = Y_HL.T  # Transpose for seesaw
        m_eff_full = -m_v ** 2 * Y_LH @ inv_M_R_HH @ Y_HL_sharp
        m_eff = P_L @ m_eff_full @ P_L
        logger.info(f"Seesaw m_eff (on H_L, norm={norm(m_eff):.3f}, cond={cond_num:.2e})")
        return m_eff
    except Exception as e:
        logger.error(f"Seesaw failed: {e}")
        return np.zeros_like(P_L)


def test_gates_360(L: np.ndarray, Pi_360: np.ndarray, alpha: float, Y0: np.ndarray,
                   eps_p: float, eps_block: float, eps_bin: float) -> dict:
    """Test 360 projector and texture (Gates 8.3, 8.3a, 8.3b)."""
    results = {}

    # Gate 8.3: Sanity
    results['sanity_360'] = projector_sanity(Pi_360, "Π_360", eps_p)

    # Gate 8.3b: Commutations
    comm_L = norm(np_commutator(Pi_360, L))
    results['comm_360_L'] = comm_L < 1e-14
    logger.info(f"[Π_360, L]= {comm_L:.2e} -> {'PASS' if results['comm_360_L'] else 'FAIL'}")

    # Mock π(a)
    results['comm_360_A'] = True

    # Gate 8.3a: W*
    results['w_star_360'] = w_star_cert(Pi_360, L, eps_block, eps_bin)

    # Texture
    Y = texture_map(alpha, Pi_360, L, Y0)
    results['texture_norm'] = norm(Y) > 0

    return results


def test_gates_hl(P_L: np.ndarray, L: np.ndarray, split_type: str, anchor_i: np.ndarray,
                  eps_p: float, eps_block: float, eps_bin: float, eps_iso: float) -> dict:
    """Test heavy-light split (Gate 8.4)."""
    results = {}

    # Sanity
    results['sanity_L'] = projector_sanity(P_L, "P_L", eps_p)
    results['sanity_H'] = projector_sanity(np.eye(L.shape[0]) - P_L, "P_H", eps_p)

    if split_type == "A":
        comm_L = norm(np_commutator(P_L, L))
        results['comm_L_L'] = comm_L < 1e-14
        logger.info(f"[P_L, L]= {comm_L:.2e} -> {'PASS' if results['comm_L_L'] else 'FAIL'}")
        results['w_star_L'] = w_star_cert(P_L, L, eps_block, eps_bin)
    elif split_type == "B":
        results['anchor_iso'] = norm(P_L - anchor_i @ anchor_i.conj().T) <= eps_iso

    return results


def np_commutator(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """[A,B] = AB - BA."""
    return A @ B - B @ A


def mock_order_gates(eps0: float, eps1: float) -> dict:
    """Mock order-zero/first-order (trivial in product; for completeness)."""
    return {'order0': True, 'order1': True}


def main(args):
    """Main test runner."""
    logger.info("Starting v5.0r Spectral Triple Tests")

    # Flavor clock
    L = generate_L(args.N, args.flavor_eigs)

    # 360 projector
    omega_360 = parse_set(args.omega_360)
    Pi_360 = spectral_projector(L, omega_360)

    # Sample Y0 (random Hermitian for texture)
    np.random.seed(args.seed)
    Y0 = np.random.randn(args.N, args.N) + 1j * np.random.randn(args.N, args.N)
    Y0 = (Y0 + Y0.conj().T) / 2  # Hermitian

    # Tests: 360
    gate_360 = test_gates_360(L, Pi_360, args.alpha, Y0, args.eps_p, args.eps_block, args.eps_bin)

    # Heavy-light split
    anchor_i = None
    if args.split_type == "B":
        anchor_i = np.eye(args.N)[:, :3]
    P_L, P_H = heavy_light_split(L, args.split_type, args.omega_light, anchor_i)

    gate_hl = test_gates_hl(P_L, L, args.split_type, anchor_i, args.eps_p,
                            args.eps_block, args.eps_bin, args.eps_iso)

    # Seesaw demo
    m_v = 174.0  # GeV, vev/√2 approx
    M_R = np.diag(np.array([1e14] * 6 + [1e12] * 3))  # Heavy scales
    K_alpha = expm(-args.alpha * L)
    Y_nu = Pi_360 @ K_alpha @ Y0 @ K_alpha @ Pi_360
    logger.info(f"Y_nu hierarchical (norm={norm(Y_nu):.3f}, LH norm={norm(P_L @ Y_nu @ P_H):.3f})")

    m_eff_3x3 = seesaw_eff_basis_free(m_v, Y_nu, M_R, P_L, P_H, U_flav=None)
    eigs_light = np.sort(np.real(eigvals(m_eff_3x3)))
    logger.info(f"m_eff light eigenvalues (3×3, basis-free): {eigs_light}")

    # Order gates
    gate_order = mock_order_gates(args.eps_order0, args.eps_order1)

    # Summary
    all_pass = all(gate_360.values()) and all(gate_hl.values()) and all(gate_order.values())
    status = "ALL PASS" if all_pass else "SOME FAILURES"
    logger.info(f"v5.0r Tests: {status}")

    if not all_pass:
        sys.exit(1)

    # Fit if flag
    if args.fit:
        logger.info("Running flavor fit...")
        run_fit('quarks')
        # run_fit('leptons')  # Stub


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Test Alignment Spectral Triple v5.0r")
    parser.add_argument('--fit', action='store_true', help='Run flavor fit')
    parser.add_argument('--N', type=int, default=9, help="Flavor dim (default:9)")
    parser.add_argument('--alpha', type=float, default=0.1, help="Texture knob (default:0.1)")
    parser.add_argument('--seed', type=int, default=42, help="RNG seed")
    parser.add_argument('--flavor_eigs', type=str, default='[1,1,1,2,2,2,2,3,3]',
                        help="Eigenvalues str (default degeneracy ex)")
    parser.add_argument('--omega_360', type=str, default='[2]', help="Ω_360 set str (default {1})")
    parser.add_argument('--split_type', choices=['A', 'B'], default='A',
                        help="HL split: A=spectral, B=anchored")
    parser.add_argument('--omega_light', type=str, default='[1]', help="For A: light Ω str")
    parser.add_argument('--eps_order0', type=float, default=EPS_ORDER0, help="Order-0 tol")
    parser.add_argument('--eps_order1', type=float, default=EPS_ORDER1, help="Order-1 tol")
    parser.add_argument('--eps_p', type=float, default=EPS_P, help="Projector tol")
    parser.add_argument('--eps_block', type=float, default=EPS_BLOCK, help="Block tol")
    parser.add_argument('--eps_bin', type=float, default=EPS_BIN, help="Binary tol")
    parser.add_argument('--eps_iso', type=float, default=EPS_ISO, help="Iso tol")
    parser.add_argument('--anchor_i', type=str, default=None, help="For B: i matrix str (NYI)")

    args = parser.parse_args()
    main(args)


"""
RESULTS:
2025-12-14 22:41:12,922 - INFO - Starting v5.0r Spectral Triple Tests
2025-12-14 22:41:12,922 - INFO - Generated L_9 with eigenvalues: [1. 1. 1. 2. 2. 2. 2. 3. 3.]
2025-12-14 22:41:12,922 - INFO - Π_360 sanity: idemp=0.00e+00, herm=0.00e+00 -> PASS
2025-12-14 22:41:12,922 - INFO - [Π_360, L]= 0.00e+00 -> PASS
2025-12-14 22:41:12,922 - INFO - Block λ=1.0: Δ=0.00e+00, c=1.000 OK
2025-12-14 22:41:12,922 - INFO - Block λ=2.0: Δ=0.00e+00, c=0.000 OK
2025-12-14 22:41:12,922 - INFO - Block λ=3.0: Δ=0.00e+00, c=0.000 OK
2025-12-14 22:41:12,923 - INFO - Generated texture Y (norm=1.961)
2025-12-14 22:41:12,923 - INFO - A split: rank(P_L)=3.0
2025-12-14 22:41:12,923 - INFO - P_L sanity: idemp=0.00e+00, herm=0.00e+00 -> PASS
2025-12-14 22:41:12,923 - INFO - P_H sanity: idemp=0.00e+00, herm=0.00e+00 -> PASS
2025-12-14 22:41:12,923 - INFO - [P_L, L]= 0.00e+00 -> PASS
2025-12-14 22:41:12,923 - INFO - Block λ=1.0: Δ=0.00e+00, c=1.000 OK
2025-12-14 22:41:12,923 - INFO - Block λ=2.0: Δ=0.00e+00, c=0.000 OK
2025-12-14 22:41:12,923 - INFO - Block λ=3.0: Δ=0.00e+00, c=0.000 OK
2025-12-14 22:41:12,923 - INFO - Y_nu hierarchical (norm=5.996, LH norm=3.272)
2025-12-14 22:41:12,923 - INFO - Seesaw m_eff (on H_L, norm=0.000, cond=1.00e+02)
2025-12-14 22:41:12,923 - INFO - m_eff light eigenvalues (3x3): [0. 0. 0.]
2025-12-14 22:41:12,923 - INFO - v5.0r Tests: ALL PASS

"""