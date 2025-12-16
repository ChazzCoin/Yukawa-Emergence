# FE/evolution.py (or evolution_operator.py)

import numpy as np

def build_misalignment_operator(L, Pphi):
    """
    FULL EMERGENCE, site-space only:

    Misalignment operator acts on the internal crystal (site basis).
    It should not mix in the 3x3 generation operator R directly.

    We build:
        M_op = aL * L + aP * (I - Pphi)

    where:
      - L     is the graph Laplacian (N x N)
      - Pphi  is a phase-coherence projector in site space (N x N)
      - aL,aP are emergent weights from spectral variances, no tuning.
    """
    n = L.shape[0]

    # Emergent weights (no knobs)
    var_L = np.var(np.diag(L))
    var_P = np.var(np.diag(Pphi))

    denom = var_L + var_P + 1e-12
    aL = var_L / denom if denom > 0 else 1.0
    aP = var_P / denom if denom > 0 else 0.0

    M_op = aL * L + aP * (np.eye(n) - Pphi)
    return M_op


def evolve_to_fixed_point(M_op, psi0, tol=1e-12, max_iter=20000):
    """
    Continuous semigroup evolution in site space:
        psi(t+dt) = psi(t) - dt * M_op psi(t)

    dt is chosen emergently from the spectral radius of M_op.
    """
    psi = psi0.astype(float).copy()

    # Emergent stable time step
    eigvals = np.linalg.eigvals(M_op)
    lam_max = np.max(np.abs(eigvals))
    dt = 1.0 / (lam_max + 1e-12)

    for _ in range(max_iter):
        psi_new = psi - dt * (M_op @ psi)
        if np.linalg.norm(psi_new - psi) < tol:
            return psi_new
        psi = psi_new

    return psi