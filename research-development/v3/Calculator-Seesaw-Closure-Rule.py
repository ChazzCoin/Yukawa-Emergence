#!/usr/bin/env python3
# ------------------------------------------------------------
#  Alignment Spectral Triple v3.3 — Geometry-Selected Flavor
#  Implements:
#   • non-convolution kernel (B operator)
#   • phase + magnitude separation
#   • broken μ–τ symmetry
#   • orthonormal but asymmetric triad compression
# ------------------------------------------------------------

import numpy as np
from numpy.linalg import eigh
from scipy.linalg import expm

def apply_misalignment_flow(K, beta):
    """
    Implements the Alignment evolution operator:

        M = exp(-beta * K)

    where K is Hermitian.
    """
    return expm(-beta * K)
def emergent_C360_projector(K, beta, rel_cut=0.15):
    """
    Build the emergent C360 projector from stability under misalignment flow.

    Parameters
    ----------
    K : Hermitian matrix
        Alignment kernel.
    beta : float
        Misalignment flow strength.
    rel_cut : float
        Relative cutoff: keep modes with
        exp(-beta * lambda_i) >= rel_cut * max(exp(-beta * lambda))

    Returns
    -------
    P : projector matrix (same size as K)
    kept_indices : list of kept eigenmode indices
    """

    # Diagonalize K
    evals, evecs = np.linalg.eigh(K)

    # Apply one-step flow to eigenvalues
    flowed = np.exp(-beta * evals)

    # Stability criterion (relative)
    max_val = flowed.max()
    keep = flowed >= rel_cut * max_val

    kept_indices = np.where(keep)[0]

    # Build projector
    P = np.zeros_like(K, dtype=complex)
    for i in kept_indices:
        v = evecs[:, i:i+1]
        P += v @ v.conj().T

    return P, kept_indices

# ============================================================
# 1. CONFIGURATION
# ============================================================

class AlignmentV33Config:

    # -------- group --------
    group_elements = [(i, j) for i in range(3) for j in range(3)]
    subgroup_H = [(0, 0), (1, 1), (2, 2)]
    triad_shifts = [(0, 0), (1, 0), (0, 1)]

    # -------- spectral kernel characters (phase only) --------
    kernel_characters = [
        (1, 0,  1.0),
        (0, 1,  0.6),
        (1, 1,  0.35),
    ]

    # -------- geometry weights α(g)  (implements B) --------
    # small, smooth, nonuniform → breaks translation invariance
    geometry_weights = {
        (0,0): 1.00, (0,1): 0.92, (0,2): 0.85,
        (1,0): 0.95, (1,1): 1.10, (1,2): 0.88,
        (2,0): 0.80, (2,1): 0.90, (2,2): 1.05,
    }

    # -------- geometric damping W(g,h) --------
    # misalignment functional (distance on torus)
    damping_strength = 0.35

    # -------- compression --------
    compression_characters = [
        (0, 0),
        (1, 0),
        (1, 1),
    ]

    higgs_vev = 174.0


# ============================================================
# 2. GROUP OPS
# ============================================================

def add_g(a, b):
    return ((a[0]+b[0]) % 3, (a[1]+b[1]) % 3)

def sub_g(a, b):
    return ((a[0]-b[0]) % 3, (a[1]-b[1]) % 3)


# ============================================================
# 3. CHARACTERS
# ============================================================

def chi(g, p, q):
    i, j = g
    return np.exp(2j*np.pi*(p*i+q*j)/3.0)


# ============================================================
# 4. KERNEL WITH GEOMETRY SELECTION
# ============================================================

def build_kernel(cfg):

    G = cfg.group_elements
    n = len(G)
    K = np.zeros((n,n), dtype=complex)

    for a,g in enumerate(G):
        for b,h in enumerate(G):

            # spectral part
            F = sum(
                w * chi(sub_g(g,h), p, q)
                for (p,q,w) in cfg.kernel_characters
            )

            # geometry weights (B operator)
            alpha_g = cfg.geometry_weights[g]
            alpha_h = cfg.geometry_weights[h]

            # misalignment damping
            dist = min(
                abs(g[0]-h[0]), 3-abs(g[0]-h[0])
            ) + min(
                abs(g[1]-h[1]), 3-abs(g[1]-h[1])
            )
            W = np.exp(-cfg.damping_strength * dist)

            K[a,b] = alpha_g * F * np.conj(alpha_h) * W

    # enforce Hermitian
    return 0.5*(K + K.conj().T)


# ============================================================
# 5. TRIADS AND COMPRESSION
# ============================================================

def build_triads(cfg):
    index = {g:i for i,g in enumerate(cfg.group_elements)}
    triads = []
    for s in cfg.triad_shifts:
        triads.append([index[add_g(h,s)] for h in cfg.subgroup_H])
    return triads


def build_S(cfg, triads):

    G = cfg.group_elements
    S = np.zeros((3,9), dtype=complex)

    for i, triad in enumerate(triads):
        p,q = cfg.compression_characters[i]
        for idx in triad:
            g = G[idx]
            S[i,idx] = chi(g,p,q) / np.sqrt(3)

    return S


# ============================================================
# 6. YUKAWA + DIAGONALIZATION
# ============================================================

def effective_yukawa(K,S):
    return S @ K @ S.conj().T


def diagonalize(Y, vev):
    vals, U = eigh(Y)
    return vals, U, np.abs(vals)*vev

def harmonic_blocks(K, beta, tol=1e-3):
    """
    Identify C360-stable harmonic blocks by eigenvalue degeneracy.
    """
    evals, evecs = np.linalg.eigh(K)
    flowed = np.exp(-beta * evals)

    blocks = []
    used = set()

    for i, val in enumerate(flowed):
        if i in used:
            continue
        block = [i]
        for j in range(i+1, len(flowed)):
            if abs(flowed[j] - val) < tol:
                block.append(j)
        for j in block:
            used.add(j)
        blocks.append(block)

    return blocks
def harmonic_blocks_by_degeneracy(A, tol_rel=0.02):
    """
    Group eigenmodes into blocks if eigenvalues are close *relative* to scale.
    Works best on the effective 3x3 Y (or any small Hermitian operator).
    """
    evals, _ = np.linalg.eigh(A)
    evals = np.sort(np.real(evals))

    blocks = []
    block = [0]
    for i in range(1, len(evals)):
        scale = max(1.0, abs(evals[i-1]), abs(evals[i]))
        if abs(evals[i] - evals[i-1]) <= tol_rel * scale:
            block.append(i)
        else:
            blocks.append(block)
            block = [i]
    blocks.append(block)
    return blocks, evals
def apply_seesaw_to_block(Y, block, M):
    """
    Apply a minimal seesaw to a given harmonic block of Y.

    Parameters
    ----------
    Y : Hermitian matrix (effective Yukawa, e.g. 3x3)
    block : list of indices (e.g. [1,2])
    M : heavy alignment scale

    Returns
    -------
    light_eigs : light eigenvalues after seesaw
    """

    # Extract the block
    Yb = Y[np.ix_(block, block)]

    n = len(block)
    zero = np.zeros_like(Yb)
    MR = M * np.eye(n)

    # Seesaw extension
    big = np.block([
        [zero, Yb],
        [Yb.conj().T, MR]
    ])

    eigvals, _ = np.linalg.eigh(big)

    # Sort by absolute value
    eigvals = np.sort(np.abs(eigvals))

    # Return the light modes only
    return eigvals[:n]

# ============================================================
# 7. MAIN
# ============================================================

def run():

    cfg = AlignmentV33Config()

    triads = build_triads(cfg)
    K = build_kernel(cfg)
    S = build_S(cfg, triads)
    beta = 1.5
    print(harmonic_blocks(K, beta))
    # 1. Misalignment flow
    K_flow = apply_misalignment_flow(K, beta)

    # 2. Emergent C360 projector
    P_C360, kept = emergent_C360_projector(K, beta)

    # 3. Project kernel onto harmonic subspace
    K_proj = P_C360 @ K_flow @ P_C360

    # 4. Effective Yukawa
    Y = effective_yukawa(K_proj, S)

    vals, U, masses = diagonalize(Y, cfg.higgs_vev)

    print("\nEigenvalues:", vals)
    print("\nMasses [GeV]:", masses)
    print("\n|U|:\n", np.abs(U))

    # Identify harmonic blocks (already done)
    blocks, _ = harmonic_blocks_by_degeneracy(Y, tol_rel=0.03)

    print("\nHarmonic blocks:", blocks)

    # Apply seesaw only to blocks with dim >= 2
    for block in blocks:
        if len(block) >= 2:
            light = apply_seesaw_to_block(Y, block, M=1e6)
            print(f"\nSeesaw applied to block {block}")
            print("Light eigenvalues after seesaw:", light)


if __name__ == "__main__":
    run()

"""
RESULTS:
[[0], [1], [2], [3], [4], [5], [6], [7], [8]]

Eigenvalues: [0.20443467 0.30696381 0.31326296]

Masses [GeV]: [35.57163219 53.41170349 54.50775485]

|U|:
 [[0.99858093 0.04133058 0.03358422]
 [0.02739089 0.89481925 0.44558752]
 [0.04567117 0.44451129 0.89460822]]

Harmonic blocks: [[0], [1, 2]]

Seesaw applied to block [1, 2]
Light eigenvalues after seesaw: [9.41276264e-08 9.80979446e-08]
"""