import numpy as np
from scipy.linalg import expm

def cycle_laplacian(n: int) -> np.ndarray:
    L = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        L[i, i] = 2.0
        L[i, (i - 1) % n] = -1.0
        L[i, (i + 1) % n] = -1.0
    return L

def topk_subspace_from_operator(X: np.ndarray, k: int, which: str = "lowest") -> tuple[np.ndarray, np.ndarray]:
    """
    Returns (Vk, w_sel) where Vk is n×k with orthonormal columns spanning the chosen subspace.
    Assumes X is Hermitian (or close; caller should symmetrize).
    """
    w, v = np.linalg.eigh(X)
    idx = np.argsort(w)
    if which == "highest":
        idx = idx[::-1]
    sel = idx[:k]
    Vk = v[:, sel]         # n×k
    w_sel = w[sel]         # k
    return Vk, w_sel

def anchor(D: np.ndarray, alpha: float, light_dim: int = 3, which: str = "lowest") -> dict:
    """
    Heat-kernel filter + subspace selection + compression.

    Returns:
      - D_eff: k×k effective operator (what you usually want to carry forward)
      - Pi: n×n projector onto chosen subspace
      - V: n×k isometry (basis for subspace)
      - X: filtered operator (n×n)
      - L, K: supporting operators
    """
    n = D.shape[0]
    L = cycle_laplacian(n).astype(np.complex128)
    K = expm((-alpha) * L)  # heat-kernel style filter

    # Filter
    X = K @ D @ K
    Xh = (X + X.conj().T) / 2.0  # enforce Hermitian for stable spectral selection

    # Choose subspace
    V, w_sel = topk_subspace_from_operator(Xh, k=light_dim, which=which)
    Pi = V @ V.conj().T  # n×n projector

    # Compress to k×k (no crushing)
    D_eff = V.conj().T @ Xh @ V
    D_eff = (D_eff + D_eff.conj().T) / 2.0

    return {
        "D_eff": D_eff,
        "Pi": Pi,
        "V": V,
        "w_sel": w_sel,
        "X": Xh,
        "L": L,
        "K": K,
    }
