# FE/selection_apply.py

import numpy as np


def projector_B_from_L(L, n_light=3):
    """
    Sublattice selector B in site space.

    FULL EMERGENCE:
      - Take the n_light lowest Laplacian modes.
      - Project onto their span.
    """
    evals, evecs = np.linalg.eigh(L)
    idx = np.argsort(evals)[:n_light]
    V = evecs[:, idx]          # N × n_light
    P_B = V @ V.T              # N × N projector
    return P_B


def projector_Pphi_from_phi(phi, k_list):
    """
    Phase-coherence projector P_phi in site space.

    FULL EMERGENCE:
      - For each divisor harmonic k in k_list, build a phase vector
            v_k(i) = exp(i k φ_i).
      - Sum them and project onto the resulting coherent direction.
    """
    phi = np.asarray(phi).reshape(-1)
    N = phi.shape[0]

    v_total = np.zeros(N, dtype=complex)
    for k in k_list:
        v_total += np.exp(1j * k * phi)

    # normalize
    norm = np.linalg.norm(v_total)
    if norm < 1e-15:
        # fall back to trivial identity if phases cancel too much
        return np.eye(N, dtype=complex)

    v = v_total / norm
    P_phi = np.outer(v, v.conj())   # rank-1 projector
    return P_phi


def build_selection_operator(phi, L, k_list):
    """
    FULL-EMERGENT SELECTION OPERATOR IN SITE SPACE:

        S_site = B * P_phi * B

    where:
      - B      : sublattice (light) projector from Laplacian.
      - P_phi  : phase-coherence projector built from φ and k_list.

    All are N×N; no mixing with 3×3 generation space here.
    """
    # N×N
    P_B   = projector_B_from_L(L)
    P_phi = projector_Pphi_from_phi(phi, k_list)

    S = P_B @ P_phi @ P_B

    # If it's effectively real, strip tiny imaginary parts
    if np.allclose(S.imag, 0.0, atol=1e-12):
        S = S.real

    return S