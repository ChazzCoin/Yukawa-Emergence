import numpy as np

def spectral_action_gaussian(D: np.ndarray, Lambda: float = 1.0) -> float:
    # simple, stable proxy: Tr exp(-(D^2)/Lambda^2)
    H = (D + D.conj().T) / 2.0
    w = np.linalg.eigvalsh(H)
    return float(np.sum(np.exp(-(w*w)/(Lambda*Lambda))))
