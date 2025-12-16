import math

import numpy as np


class Geometry:
    N_SITES = 9
    LIGHT_SITES = [0, 1, 2]
    HEAVY_SITES = [3, 4, 5, 6, 7, 8]

    @staticmethod
    def generation_index(i: int) -> int:
        return i % 3

def chi2_from_obs(obs):
    chi2 = 0.0
    pulls = {}
    for key, xexp in exp_targets.items():
        xth = obs.get(key, np.nan)
        if not np.isfinite(xth):
            continue
        sig = sigma_targets[key]
        pull = (xth - xexp) / sig
        chi2 += pull**2
        pulls[key] = pull
    return chi2, pulls

def triadic_fourier_modes():
    n = 9
    j = np.arange(n)
    v0 = np.exp(2j*np.pi*0*j/n) / np.sqrt(n)
    v3 = np.exp(2j*np.pi*3*j/n) / np.sqrt(n)
    v6 = np.exp(2j*np.pi*6*j/n) / np.sqrt(n)
    return v0, v3, v6

def build_default_projection(proj_eps_03: float = 0.0,
                             proj_eps_30: float = 0.0,
                             proj_eps_36: float = 0.0) -> np.ndarray:
    v0, v3, v6 = triadic_fourier_modes()

    b0 = v0.copy()
    b3 = v3.copy()
    b6 = v6.copy()

    b0 = b0 + proj_eps_03 * v3
    b3 = b3 + proj_eps_30 * v0 + proj_eps_36 * v6
    b6 = b6 + proj_eps_36 * v3

    B = np.vstack([b0, b3, b6])          # (3, 9)
    Q_cols, _ = np.linalg.qr(B.conj().T) # (9, 3)
    P = Q_cols.T                         # (3, 9)
    return P

def tweak_projection(P_base, g, s, a, phi):
    P = P_base.copy()
    P[g, s] += a * np.exp(1j * phi)

    Q = np.zeros_like(P, dtype=complex)

    Q[0] = P[0] / np.linalg.norm(P[0])

    v1 = P[1] - np.vdot(Q[0], P[1]) * Q[0]
    Q[1] = v1 / np.linalg.norm(v1)

    v2 = P[2] \
         - np.vdot(Q[0], P[2]) * Q[0] \
         - np.vdot(Q[1], P[2]) * Q[1]
    Q[2] = v2 / np.linalg.norm(v2)

    return Q

""""""

def build_phase_profile_gen(n0_tilde: int, delta_tilde: int, N_eff: int = 360) -> np.ndarray:
    q = 360 // N_eff
    n0_eff = q * n0_tilde
    delta_eff = q * delta_tilde
    phi_gen = np.zeros(3, dtype=float)
    for g in range(3):
        angle_deg = n0_eff + g * delta_eff
        phi_gen[g] = 2.0 * math.pi * angle_deg / 360.0
    return phi_gen

def build_site_phases(phi_gen: np.ndarray) -> np.ndarray:
    phi_site = np.zeros(Geometry.N_SITES, dtype=float)
    for i in range(Geometry.N_SITES):
        phi_site[i] = phi_gen[Geometry.generation_index(i)]
    return phi_site

def build_phase_matrix(phi_site: np.ndarray) -> np.ndarray:
    N = len(phi_site)
    P = np.zeros((N, N), dtype=complex)
    for i in range(N):
        for j in range(N):
            P[i, j] = np.exp(1j * (phi_site[i] - phi_site[j]))
    return P