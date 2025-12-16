import numpy as np

def diagonalize_hermitian(M: np.ndarray) -> dict:
    H = (M + M.conj().T) / 2.0
    w, v = np.linalg.eigh(H)
    return {"eigenvalues": w, "eigenvectors": v}
