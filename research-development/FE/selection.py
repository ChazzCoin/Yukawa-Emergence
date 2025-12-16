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
    Phase-coherence projector P_phi in site space (single, non-sector).

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

    norm = np.linalg.norm(v_total)
    if norm < 1e-15:
        return np.eye(N, dtype=complex)

    v = v_total / norm
    P_phi = np.outer(v, v.conj())   # rank-1 projector
    return P_phi


def projector_Pphi_sector(phi, k_list, multiplier):
    """
    Sector-specific phase-coherence projector.

    FULL EMERGENCE:
      - Use harmonic multiples of the same triad k_list:
            v_k^sector(i) = exp(i * multiplier * k * φ_i).
      - multiplier is an integer from D_360 (1,2,3,4,...)
        and we choose different ones per sector (u,d,e,nu).

    This gives distinct coherent directions in site space
    from the same underlying φ and k_list.
    """
    phi = np.asarray(phi).reshape(-1)
    N = phi.shape[0]

    v_total = np.zeros(N, dtype=complex)
    for k in k_list:
        k_eff = (multiplier * k) % 360
        v_total += np.exp(1j * k_eff * phi)

    norm = np.linalg.norm(v_total)
    if norm < 1e-15:
        return np.eye(N, dtype=complex)

    v = v_total / norm
    P_phi = np.outer(v, v.conj())
    return P_phi


def build_selection_operator(phi, L, k_list):
    """
    Single selection operator in site space (legacy).

        S_site = B * P_phi * B
    """
    P_B   = projector_B_from_L(L)
    P_phi = projector_Pphi_from_phi(phi, k_list)

    S = P_B @ P_phi @ P_B
    if np.allclose(S.imag, 0.0, atol=1e-12):
        S = S.real
    return S


def build_sector_selection_operators(phi, L, k_list):
    """
    FULL-EMERGENT SECTOR-RESOLVED SELECTION:

        For each sector s ∈ {u,d,e,nu}, build:

            S_s = B * P_phi^(s) * B

        where P_phi^(s) differs only by an integer harmonic multiplier
        drawn from D_360:

            u  → multiplier 1
            d  → multiplier 2
            e  → multiplier 3
            nu → multiplier 4

        No continuous parameters, no tuning. All sector
        differences come from discrete harmonic structure.
    """
    P_B = projector_B_from_L(L)

    multipliers = {
        "u": 1,
        "d": 2,
        "e": 3,
        "nu": 4,
    }

    S_sectors = {}
    for s, m in multipliers.items():
        P_phi_s = projector_Pphi_sector(phi, k_list, multiplier=m)
        S_s = P_B @ P_phi_s @ P_B
        if np.allclose(S_s.imag, 0.0, atol=1e-12):
            S_s = S_s.real
        S_sectors[s] = S_s

    return S_sectors