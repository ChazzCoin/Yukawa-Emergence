"""
flavor_common.py -- shared helpers for the theta23-paper audit scripts.

What it provides: PDG parametrization U(t12,t13,t23,delta), extraction of the PDG
angles back out of a unitary matrix, Hermitian-from-spectrum builders, the 2x2 Bloch
form used by theta23_obstruction_paper.tex, and a Fourier-content helper.

Inputs: none (library). Angles are in radians unless a name ends in _deg.
Source of the PDG convention: Review of Particle Physics, "Neutrino Mixing" (standard
parametrization U = R23 . U13(delta) . R12).
Run: imported by the sibling scripts; `python3 flavor_common.py` runs a self-check.
Made by: crew flavor-obstruction-paper, item #1 (R1), data-analyst.
"""
import numpy as np

SX = np.array([[0, 1], [1, 0]], complex)
SZ = np.array([[1, 0], [0, -1]], complex)


def rot2(t):
    return np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]], complex)


def embed_E(R2):
    """1 (+) R2 : acts on the (2,3) plane, identity on index 1."""
    W = np.eye(3, dtype=complex)
    W[1:, 1:] = R2
    return W


def pdg_U(t12, t13, t23, delta):
    c12, s12, c13, s13, c23, s23 = (np.cos(t12), np.sin(t12), np.cos(t13),
                                    np.sin(t13), np.cos(t23), np.sin(t23))
    e = np.exp(1j * delta)
    return np.array([
        [c12 * c13, s12 * c13, s13 / e],
        [-s12 * c23 - c12 * s23 * s13 * e, c12 * c23 - s12 * s23 * s13 * e, s23 * c13],
        [s12 * s23 - c12 * c23 * s13 * e, -c12 * s23 - s12 * c23 * s13 * e, c23 * c13]])


def pdg_angles(U):
    """(t12, t13, t23, |delta|) in radians from |U_ij| (rephasing invariant).
    |delta| in [0, pi] from |U_tau1|^2; sign of delta from the Jarlskog sign."""
    a = np.abs(U)
    t13 = np.arcsin(min(a[0, 2], 1.0))
    t12 = np.arctan2(a[0, 1], a[0, 0])
    t23 = np.arctan2(a[1, 2], a[2, 2])
    s12, c12, s23, c23, s13 = np.sin(t12), np.cos(t12), np.sin(t23), np.cos(t23), np.sin(t13)
    den = 2 * s12 * c12 * s23 * c23 * s13
    if den < 1e-12:
        return t12, t13, t23, float('nan')
    cd = (s12**2 * s23**2 + c12**2 * c23**2 * s13**2 - a[2, 0]**2) / den
    d = np.arccos(np.clip(cd, -1, 1))
    J = np.imag(U[0, 0] * U[1, 1] * np.conj(U[0, 1]) * np.conj(U[1, 0]))
    return t12, t13, t23, (d if J >= 0 else -d)


def herm_from(eigs, V):
    return V @ np.diag(eigs) @ V.conj().T


def bloch2(gap, xi, c):
    """2x2 Hermitian, eigenvalues c +- gap/2, larger-eigenvalue eigenvector (cos xi, sin xi).
    (Independent of the paper's construction: built from the projector, not from sigma.)"""
    v = np.array([np.cos(xi), np.sin(xi)], complex)
    P = np.outer(v, v.conj())
    return (c - gap / 2) * np.eye(2) + gap * P


def fourier_amps(vals, kmax=8):
    """vals sampled uniformly over [0, 2pi). Returns amplitudes of cos/sin pairs k=0..kmax."""
    f = np.fft.rfft(vals) / len(vals)
    amps = np.abs(f[:kmax + 1]) * 2.0
    amps[0] = abs(f[0])
    return amps


if __name__ == "__main__":
    # self-check against known cases
    U = pdg_U(0.58, 0.15, 0.85, 1.2)
    assert np.allclose(U @ U.conj().T, np.eye(3), atol=1e-14)
    r = pdg_angles(U)
    assert np.allclose(r, (0.58, 0.15, 0.85, 1.2), atol=1e-12), r
    print("flavor_common self-check OK: unitary, PDG round-trip to 1e-12", r)
