# FE/yukawa.py

import numpy as np


def normalize(v):
    v = np.asarray(v)
    n = np.linalg.norm(v)
    if n < 1e-15:
        return v
    return v / n


# def build_generation_basis_from_psi_sector(psi_s, R, triad_vecs):
#     """
#     Sector-resolved generation basis.
#
#     psi_s      : (N,) manifested internal state for sector s (site space).
#     triad_vecs : (N,3) Laplacian eigenvectors for triad modes.
#     R          : (3,3) generation-space harmonic operator.
#
#     Steps:
#       1) Project psi_s into the 3D triad subspace:
#            psi_gen = V^T psi_s  (3,)
#       2) Use psi_gen as seed in 3D generation space.
#       3) Apply R repeatedly to build a triad basis.
#     """
#     psi_s = np.asarray(psi_s).reshape(-1)
#     triad_vecs = np.asarray(triad_vecs)  # (N,3)
#
#     psi_gen = triad_vecs.T @ psi_s       # (3,)
#
#     g0 = normalize(psi_gen)
#     g1 = normalize(R @ g0)
#     g2 = normalize(R @ g1)
#
#     G_s = np.column_stack([g0, g1, g2])  # (3,3)
#     return G_s


def universal_kernel_F(R):
    """
    Universal spectral kernel F(R).

    For diagonal R with entries e^{iθ_j}, we use:
        F_j = exp(-0.5 * (1 - Re(e^{iθ_j})))
    """
    lam = np.diag(R)
    Re = np.real(lam)
    Fvals = np.exp(-0.5 * (1.0 - Re))
    return Fvals


def build_yukawa_sector(Fvals, Q, U):
    """
    Yukawa in generation space:

        Y = U^\dagger  diag(Fvals)  exp(-Q_diag)  U

    Fvals : (3,) kernel values from R
    Q     : (3,3) diagonal charge operator
    U     : (3,3) sector-specific generation basis
    """
    Fvals = np.asarray(Fvals).reshape(-1)
    q_diag = np.diag(Q)
    E = np.diag(np.exp(-q_diag))

    D_F = np.diag(Fvals)
    Y = U.conj().T @ D_F @ E @ U
    return Y


def build_yukawas_from_manifested_state(psi_sectors, R, Q, evecs, mode_indices):
    """
    FULL EMERGENCE, SECTOR-RESOLVED:

    Inputs:
      psi_sectors : dict with keys "u","d","e","nu", each (N,) vector
      R           : (3,3) generation harmonic operator
      Q           : (3,3) integer charge operator
      evecs       : (N,N) Laplacian eigenvectors
      mode_indices: length-3 triad indices from build_R_three

    Steps:
      1) Build sector-specific generation bases G_s from psi_s.
      2) Compute universal kernel F(R).
      3) Construct sector Yukawas from (F, Q, G_s) with discrete patterns.
      4) Extract spectra and mixing matrices.
    """
    triad_vecs = evecs[:, mode_indices]  # (N,3)

    # 1) sector-specific generation bases
    G_u  = build_generation_basis_from_psi_sector(psi_sectors["u"],  R, triad_vecs)
    G_d  = build_generation_basis_from_psi_sector(psi_sectors["d"],  R, triad_vecs)
    G_e  = build_generation_basis_from_psi_sector(psi_sectors["e"],  R, triad_vecs)
    G_nu = build_generation_basis_from_psi_sector(psi_sectors["nu"], R, triad_vecs)

    # 2) universal kernel from R
    F_base = universal_kernel_F(R)

    # 3) sector Yukawas (discrete patterns, no tuning)
    Y_u  = build_yukawa_sector(F_base,             Q, G_u)
    Y_d  = build_yukawa_sector(F_base[::-1],       Q, G_d)
    Y_e  = build_yukawa_sector(F_base**2,          Q, G_e)
    Y_nu = build_yukawa_sector(np.abs(F_base),     Q, G_nu)

    # 4) singular values (mass spectra)
    def sv(Y):
        return np.sort(np.linalg.svd(Y, compute_uv=False))[::-1]

    spectra = {
        "u":  sv(Y_u),
        "d":  sv(Y_d),
        "e":  sv(Y_e),
        "nu": sv(Y_nu),
    }

    # 5) emergent mixing matrices:
    #    left-handed mixing only depends on left bases G_s

    V_ckm  = G_u.conj().T  @ G_d
    U_pmns = G_e.conj().T  @ G_nu

    mixings = {
        "CKM":  V_ckm,
        "PMNS": U_pmns,
    }

    return Y_u, Y_d, Y_e, Y_nu, spectra, mixings

def build_generation_basis_from_psi_sector(psi_s, R, triad_vecs):
    """
    Sector-resolved generation basis with orthonormalization.

    psi_s      : (N,) manifested internal state for sector s (site space).
    triad_vecs : (N,3) Laplacian eigenvectors for triad modes.
    R          : (3,3) generation-space harmonic operator.

    Steps:
      1) Project psi_s into 3D triad subspace: psi_gen = V^T psi_s.
      2) Build a 3×3 seed matrix from {psi_gen, R psi_gen, R^2 psi_gen}.
      3) Orthonormalize via QR to get a sharp unitary-like G_s.
    """
    psi_s = np.asarray(psi_s).reshape(-1)
    triad_vecs = np.asarray(triad_vecs)  # (N,3)

    # 1) Project into generation subspace
    psi_gen = triad_vecs.T @ psi_s       # (3,)

    # If psi_gen is essentially zero, fall back to a trivial basis
    if np.linalg.norm(psi_gen) < 1e-15:
        return np.eye(3, dtype=complex)

    # 2) Build seed matrix in generation space
    v0 = psi_gen
    v1 = R @ v0
    v2 = R @ v1

    W = np.column_stack([v0, v1, v2])    # 3×3

    # 3) Orthonormalize (QR decomposition)
    #    This enforces an orthonormal generation basis without adding knobs.
    Q, _ = np.linalg.qr(W)

    return Q  # (3×3) sector-specific, orthonormal basis