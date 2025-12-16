import numpy as np

def schur_complement(M: np.ndarray, light_dim: int = 3, ridge: float = 1e-12) -> np.ndarray:
    """
    Schur complement reduction to a light_dim×light_dim effective operator.

    Robust behavior:
      - If M is already light_dim×light_dim: return M (Hermitian-symmetrized) as a no-op.
      - If light_dim >= n: error (unless equal, handled above).
    """
    n = M.shape[0]

    # Already reduced -> no-op
    if n == light_dim:
        H = (M + M.conj().T) / 2.0
        return H

    if not (1 <= light_dim < n):
        raise ValueError(f"light_dim must be in [1, n-1], got light_dim={light_dim}, n={n}")

    L = slice(0, light_dim)
    Hh = slice(light_dim, n)

    M_LL = M[L, L]
    M_LH = M[L, Hh]
    M_HL = M[Hh, L]
    M_HH = M[Hh, Hh]

    # Stabilize inversion
    M_HH_reg = M_HH + ridge * np.eye(n - light_dim, dtype=M.dtype)

    # Solve M_HH_reg * X = M_HL  (more stable than explicit inverse)
    X = np.linalg.solve(M_HH_reg, M_HL)

    M_eff = M_LL - M_LH @ X
    M_eff = (M_eff + M_eff.conj().T) / 2.0  # keep Hermitian
    return M_eff
