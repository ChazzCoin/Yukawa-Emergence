# FE/harmonics.py
import numpy as np

# ------------------------------------------------
# Divisors of 360 (D_360)
# ------------------------------------------------
def divisors_360():
    n = 360
    divs = []
    for k in range(1, n + 1):
        if n % k == 0:
            divs.append(k)
    return np.array(divs, dtype=int)

D360 = divisors_360()


# ------------------------------------------------
# Helper: choose triad periods from eigenvalue triad
# ------------------------------------------------
def _best_divisor_triad_for_evals(triad_evals):
    """
    Given three positive eigenvalues λ1 < λ2 < λ3,
    choose a triad (T1,T2,T3) from D360 that best matches
    the ratios λ2/λ1 and λ3/λ1.

    This enforces Axiom B1: triad of distinct divisors of 360.
    """
    lam = np.array(triad_evals, dtype=float)
    lam_sorted = np.sort(lam)
    lam1, lam2, lam3 = lam_sorted

    r2 = lam2 / lam1
    r3 = lam3 / lam1

    divs = D360
    best_err = np.inf
    best_triad = None

    # scan ordered triples T1 < T2 < T3
    for i in range(len(divs)):
        for j in range(i + 1, len(divs)):
            for k in range(j + 1, len(divs)):
                T1, T2, T3 = float(divs[i]), float(divs[j]), float(divs[k])
                R2 = T2 / T1
                R3 = T3 / T1
                # compare in log-space to respect ratios
                err = (np.log(r2 + 1e-12) - np.log(R2 + 1e-12))**2 \
                    + (np.log(r3 + 1e-12) - np.log(R3 + 1e-12))**2
                if err < best_err:
                    best_err = err
                    best_triad = (int(T1), int(T2), int(T3))

    return best_triad  # (T1, T2, T3)


# ------------------------------------------------
# 1. Build emergent 3×3 R operator from spectrum
# ------------------------------------------------
def build_R_three(evals: np.ndarray,
                  evecs: np.ndarray,
                  k_gen: int = 3):
    """
    Construct the emergent base-360 operator R on a 3D generation
    subspace, using only the Laplacian eigenvalues (no hand labels).

    Steps:
      1. Sort eigenvalues, skip zero mode.
      2. Take first three positive eigenvalues as the 'generation triad'.
      3. Find triad of periods (T1,T2,T3) ∈ D360^3 that best match
         the eigenvalue ratios.
      4. Set k_j = 360 / T_j and define
         R = diag(exp(2π i k_j / 360)).

    Returns:
      R        : 3x3 unitary
      k_list   : length-3 int array of cycle indices
      triad_ix : length-3 int array of eigenvalue indices used
    """

    # 1. sort and find first three positive eigenvalues
    idx_sorted = np.argsort(evals)
    evals_sorted = evals[idx_sorted]

    # ignore the zero (or smallest) mode
    # pick next three as the generation triad
    triad_ix_local = np.arange(1, 4)  # 1,2,3 in sorted space
    triad_evals = evals_sorted[triad_ix_local]

    # 2. map eigenvalue triad to divisor triad
    T1, T2, T3 = _best_divisor_triad_for_evals(triad_evals)
    T_list = np.array([T1, T2, T3], dtype=int)

    # 3. compute k_j from T_j
    # Axiom B1: T_j = 360 / gcd(360,k_j)  ⇒ choose k_j = 360 / T_j
    k_list = (360 // T_list).astype(int)

    # 4. build R
    phases = np.exp(2j * np.pi * k_list / 360.0)
    R = np.diag(phases)

    # map triad_ix_local back to original eigenvalue indices
    triad_ix = idx_sorted[triad_ix_local]

    return R, k_list, triad_ix


# ------------------------------------------------
# 2. Build emergent integer charge operator Q
# ------------------------------------------------
def build_charge_operator(evals: np.ndarray,
                          evecs: np.ndarray,
                          mode_indices: np.ndarray = None):
    """
    Build a 3x3 diagonal integer charge operator Q in the same
    generation subspace as R, using only spectral data and D_360.

    FULL EMERGENCE:

      1. Use the same triad of eigenvalues as build_R_three.
      2. Get (T1,T2,T3) for that triad via the same mapping.
      3. First attempt a ratio-based mapping:
            q_j = floor(T_max / T_j)
         If that yields a trivial, fully-degenerate vector (all equal),
         fall back to a pure ordering-based assignment:
            - smallest T → q = 1
            - middle  T → q = 2
            - largest T → q = 3

      This ensures:
        - integer charges
        - triadic hierarchy
        - no hand-crafted table or continuous fit
    """

    # If caller didn’t pass mode_indices, reuse build_R_three’s eigenvalue triad
    if mode_indices is None:
        idx_sorted = np.argsort(evals)
        evals_sorted = evals[idx_sorted]
        triad_ix_local = np.arange(1, 4)  # first three positive modes
        triad_evals = evals_sorted[triad_ix_local]
    else:
        triad_evals = evals[mode_indices]

    # Get best divisor triad (T1,T2,T3)
    T1, T2, T3 = _best_divisor_triad_for_evals(triad_evals)
    T_list = np.array([T1, T2, T3], dtype=float)

    # --- 1) ratio-based mapping (primary emergent rule) ---
    T_max = np.max(T_list)
    q_ratio = np.floor(T_max / T_list).astype(int)
    q_ratio[q_ratio < 1] = 1  # enforce positivity

    # Check if this is nontrivial (at least two distinct values)
    if len(np.unique(q_ratio)) >= 2:
        q_list = q_ratio
    else:
        # --- 2) fallback: pure ordering-based hierarchy ---
        # smaller T → faster cycle → treat as "heavier" (less suppressed)
        # assign:
        #   smallest T → q=1
        #   middle   T → q=2
        #   largest  T → q=3
        order = np.argsort(T_list)  # indices of T in ascending order
        q_temp = np.arange(1, len(T_list) + 1)  # [1,2,3]
        q_list = np.zeros_like(q_temp)
        # place 1,2,3 according to sorted T
        q_list[order] = q_temp

    Q = np.diag(q_list.astype(float))
    return Q, q_list