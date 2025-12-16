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
import matplotlib.pyplot as plt

def plot_U_heatmaps_vs_beta(
    beta_vals=(1.1, 1.3, 1.5, 1.7),
    save_path="fig_U_heatmaps_vs_beta.pdf",
    show=True
):
    """
    FIGURE 3:
    |U| heatmaps vs β for all four sectors.
    """

    sectors = make_sector_configs()
    n_beta = len(beta_vals)
    n_sec = len(sectors)

    fig, axes = plt.subplots(
        n_sec, n_beta,
        figsize=(3.0 * n_beta, 2.6 * n_sec),
        constrained_layout=True
    )

    for i, (sec, cfg) in enumerate(sectors.items()):
        for j, beta in enumerate(beta_vals):

            absU = run_sector(cfg, beta)
            ax = axes[i, j]

            im = ax.imshow(absU, cmap="viridis", vmin=0, vmax=1)

            if i == 0:
                ax.set_title(rf"$\beta = {beta}$")
            if j == 0:
                ax.set_ylabel(sec, rotation=0, labelpad=25, fontsize=12)

            ax.set_xticks([])
            ax.set_yticks([])

    cbar = fig.colorbar(im, ax=axes, fraction=0.02, pad=0.02)
    cbar.set_label(r"$|U_{ij}|$")

    fig.suptitle(
        r"Emergent mixing matrices $|U|$ vs misalignment strength $\beta$",
        fontsize=14
    )

    plt.savefig(save_path)
    if show:
        plt.show()
    plt.close()

def plot_yukawa_spectrum_and_mixing(
    cfg,
    beta=1.5,
    save_path="fig_yukawa_spectrum_and_mixing.pdf",
    show=True
):
    """
    FIGURE 2:
    (a) Normalized Yukawa eigenvalue spectrum
    (b) |U| heatmap of the emergent mixing matrix
    """

    # --- build objects ---
    triads = build_triads(cfg)
    K = build_kernel(cfg)
    S = build_S(cfg, triads)

    # alignment flow + projector
    K_flow = apply_misalignment_flow(K, beta)
    P, _ = emergent_C360_projector(K, beta)
    K_proj = P @ K_flow @ P

    # effective Yukawa
    Y = effective_yukawa(K_proj, S)

    # diagonalize
    vals, U, _ = diagonalize(Y, cfg.higgs_vev)
    vals = np.abs(vals)
    vals /= vals.max()  # normalize

    absU = np.abs(U)

    # --- plotting ---
    fig, axes = plt.subplots(
        1, 2, figsize=(9.5, 4.0),
        gridspec_kw={"width_ratios": [1.0, 1.2]}
    )

    # ---- Panel (a): spectrum ----
    ax = axes[0]
    ax.plot(
        range(1, 4), vals,
        "o-", lw=2, ms=7
    )
    ax.set_xticks([1, 2, 3])
    ax.set_xlabel("Eigenmode index")
    ax.set_ylabel(r"Normalized $|\lambda_i|$")
    ax.set_title("Emergent Yukawa spectrum")
    ax.grid(alpha=0.3)

    # ---- Panel (b): |U| heatmap ----
    ax = axes[1]
    im = ax.imshow(
        absU,
        cmap="viridis",
        vmin=0.0,
        vmax=1.0
    )

    ax.set_xticks([0, 1, 2])
    ax.set_yticks([0, 1, 2])
    ax.set_xticklabels([r"$1$", r"$2$", r"$3$"])
    ax.set_yticklabels([r"$1$", r"$2$", r"$3$"])
    ax.set_xlabel("Right-handed index")
    ax.set_ylabel("Left-handed index")
    ax.set_title(r"Mixing matrix $|U|$")

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(r"$|U_{ij}|$")

    plt.tight_layout()
    plt.savefig(save_path)
    if show:
        plt.show()
    plt.close()

def plot_flowed_eigenvalues_vs_beta(
    cfg,
    beta_vals=np.linspace(0.8, 2.2, 120),
    save_path="fig_flowed_eigenvalues_vs_beta.pdf",
    show=True
):
    """
    FIGURE 1:
    Flowed eigenvalues exp(-β λ_i(K)) vs β.

    Reveals stability-selected harmonic degeneracies.
    """

    # Build kernel once
    K = build_kernel(cfg)

    # Eigenvalues of the *bare* kernel
    lam, _ = eigh(K)
    lam = np.sort(np.real(lam))

    # Compute flowed eigenvalues for each β
    flowed = np.zeros((len(beta_vals), len(lam)))
    for i, beta in enumerate(beta_vals):
        flowed[i, :] = np.exp(-beta * lam)

    # ---- plotting ----
    plt.figure(figsize=(6.5, 4.5))

    for j in range(len(lam)):
        plt.plot(
            beta_vals,
            flowed[:, j],
            lw=2,
            label=fr"$e^{{-\beta \lambda_{{{j+1}}}}}$"
        )

    plt.xlabel(r"Misalignment strength $\beta$", fontsize=12)
    plt.ylabel(r"Flowed eigenvalues $e^{-\beta \lambda}$", fontsize=12)
    plt.title("Stability-selected harmonic spectrum", fontsize=13)

    plt.legend(frameon=False)
    plt.grid(alpha=0.25)

    plt.tight_layout()
    plt.savefig(save_path)
    if show:
        plt.show()
    plt.close()

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
"""
Axiomatic inputs:
• finite internal space Z3 × Z3
• spectral kernel characters
• misalignment flow exp(-βK)
• stability-selected projector
• triadic compression

Conventional choices (not tuned):
• β = O(1) misalignment strength
• rel_cut = 0.15 stability threshold
• smooth O(1) geometry weights
"""
# ============================================================

class AlignmentV33Config:

    # -------- group --------
    group_elements = [(i, j) for i in range(3) for j in range(3)]
    subgroup_H = [(0, 0), (1, 1), (2, 2)]
    triad_shifts = [(0, 0), (1, 0), (0, 1)]

    # -------- spectral kernel characters (phase only) --------
    kernel_characters = [
        (1, 0, 1.0),
        (0, 1, 0.7),
        (1, 1, 0.4),
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
    damping_strength = 0.45

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

def add_g(a, b): return ((a[0]+b[0]) % 3, (a[1]+b[1]) % 3)
def sub_g(a, b): return ((a[0]-b[0]) % 3, (a[1]-b[1]) % 3)

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

def effective_yukawa(K,S): return S @ K @ S.conj().T

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

def run_sector(cfg, beta):
    """
    Run the alignment pipeline for a single sector configuration.
    Returns |U|.
    """
    triads = build_triads(cfg)
    K = build_kernel(cfg)
    S = build_S(cfg, triads)

    K_flow = apply_misalignment_flow(K, beta)
    P, _ = emergent_C360_projector(K, beta)
    K_proj = P @ K_flow @ P

    Y = effective_yukawa(K_proj, S)
    _, U, _ = diagonalize(Y, cfg.higgs_vev)

    return np.abs(U)

def make_sector_configs():
    """
    Returns dict of sector -> AlignmentV33Config
    """

    base = AlignmentV33Config()

    sectors = {}

    # up-type quarks: very rigid, almost diagonal
    cfg_u = AlignmentV33Config()
    cfg_u.geometry_weights = {
        k: v * (1.05 if k[0] == k[1] else 0.95)
        for k, v in base.geometry_weights.items()
    }
    sectors["u"] = cfg_u

    # down-type quarks: slightly softer
    cfg_d = AlignmentV33Config()
    cfg_d.geometry_weights = {
        k: v * (1.02 if k[0] == k[1] else 0.98)
        for k, v in base.geometry_weights.items()
    }
    sectors["d"] = cfg_d

    # charged leptons: similar to down, marginally more mixing
    cfg_e = AlignmentV33Config()
    cfg_e.geometry_weights = {
        k: v * (1.00 if k[0] == k[1] else 1.00)
        for k, v in base.geometry_weights.items()
    }
    sectors["e"] = cfg_e

    # neutrinos: slightly softened diagonal rigidity
    cfg_nu = AlignmentV33Config()
    cfg_nu.geometry_weights = {
        k: v * (0.95 if k[0] == k[1] else 1.05)
        for k, v in base.geometry_weights.items()
    }
    sectors["nu"] = cfg_nu

    return sectors

# ============================================================
# 7. MAIN
# ============================================================

def run():

    cfg = AlignmentV33Config()

    triads = build_triads(cfg)
    mean_geom = np.mean(list(cfg.geometry_weights.values()))
    cfg.geometry_weights = {k: v / mean_geom for k, v in cfg.geometry_weights.items()}

    K = build_kernel(cfg)
    S = build_S(cfg, triads)
    beta = 1.5
    for beta in [1.2, 1.5, 1.8]:
        K_flow = apply_misalignment_flow(K, beta)
        P, _ = emergent_C360_projector(K, beta)
        Y = effective_yukawa(P @ K_flow @ P, S)
        blocks, _ = harmonic_blocks_by_degeneracy(Y)
        print(f"beta={beta}: blocks={blocks}")

    print("β acts as a control parameter for harmonic stability,")
    print("not a fitted quantity.")
    print("\nInterpretation:")
    print("• Fully split spectrum → Dirac-like sector")
    print("• Degenerate harmonic block → requires Majorana completion")
    print("• Intermediate β yields neutrino-like structure")

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
        if len(block) == 1:
            print(f"Block {block}: Dirac-like (no seesaw)")
        elif len(block) >= 2:
            light = apply_seesaw_to_block(Y, block, M=1e6)
            print(f"\nSeesaw applied to block {block}")
            print("Light eigenvalues after seesaw:", light)
            print("Seesaw suppression factor ~", light.max() / masses.max())
        else:
            print(f"Block {block}: non-closed → Majorana seesaw")

    print("\n=== Emergent Yukawa Summary ===")
    print("Eigenvalue ratios:", vals / vals.max())
    print("Mixing pattern |U|:")
    print(np.round(np.abs(U), 3))
    print("Harmonic blocks:", blocks)
    # ---- FIGURE 1: spectral flow ----
    plot_flowed_eigenvalues_vs_beta(cfg)
    # ---- FIGURE 2: Yukawa spectrum + mixing ----
    plot_yukawa_spectrum_and_mixing(cfg, beta=1.5)
    # ---- FIGURE 3: |U| vs beta across sectors ----
    plot_U_heatmaps_vs_beta()

if __name__ == "__main__":
    run()

"""
RESULTS:
beta=1.2: blocks=[[0], [1, 2]]
beta=1.5: blocks=[[0], [1], [2]]
beta=1.8: blocks=[[0], [1], [2]]
β acts as a control parameter for harmonic stability,
not a fitted quantity.

Interpretation:
• Fully split spectrum → Dirac-like sector
• Degenerate harmonic block → requires Majorana completion
• Intermediate β yields neutrino-like structure
[[0], [1], [2], [3], [4], [5], [6], [7, 8]]

Eigenvalues: [0.01868867 0.16377328 0.24118664]

Masses [GeV]: [ 3.25182833 28.49655043 41.96647612]

|U|:
 [[0.99846476 0.04280896 0.03514993]
 [0.04148383 0.99836285 0.03937913]
 [0.03670444 0.03793438 0.99860591]]

Harmonic blocks: [[0], [1], [2]]
Block [0]: Dirac-like (no seesaw)
Block [1]: Dirac-like (no seesaw)
Block [2]: Dirac-like (no seesaw)

=== Emergent Yukawa Summary ===
Eigenvalue ratios: [0.07748633 0.67903129 1.        ]
Mixing pattern |U|:
[[0.998 0.043 0.035]
 [0.041 0.998 0.039]
 [0.037 0.038 0.999]]
Harmonic blocks: [[0], [1], [2]]
"""