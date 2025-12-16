import numpy as np

# ===========================
# 1. Geometry and alignment kernel
# ===========================

def circle_distance(a, b, N=360):
    """
    Minimal distance on a circle of size N between lattice points a and b.
    """
    d = abs(a - b) % N
    return min(d, N - d)


def build_alignment_kernel(positions, kappa, forbidden_distance=None, epsilon=0.0, N=360):
    """
    Build the exponential Toeplitz-like alignment kernel K_ij = kappa^{d(i,j)}
    on a subset of Z_N, with the option to suppress a single forbidden distance.

    positions : list of ints
        Embedding of the internal sites into Z_N.
    kappa : float
        Geometric decay parameter (e.g. ~0.24).
    forbidden_distance : int or None
        If not None, entries with distance == forbidden_distance are set to epsilon.
    epsilon : float
        Value to use at the forbidden distance (typically 0 or very small).
    N : int
        Size of the ambient cycle (360 in the A360 construction).
    """
    n = len(positions)
    K = np.zeros((n, n), dtype=float)
    for i, pi in enumerate(positions):
        for j, pj in enumerate(positions):
            d = circle_distance(pi, pj, N=N)
            if forbidden_distance is not None and d == forbidden_distance:
                K[i, j] = epsilon
            else:
                K[i, j] = kappa**d
    return K


# ===========================
# 2. Alignment map on Yukawa / Majorana blocks
# ===========================

def align_matrix(M0, K):
    # Φ(X) = K ⊙ X  (elementwise)
    return K * M0


# ===========================
# 3. Diagonalization and mixing extraction
# ===========================

def hermitian_diagonalization(M):
    """
    Diagonalize a Yukawa (or mass) matrix M via its hermitian square:

        H = M^\dagger M

    Returns eigenvalues (masses >= 0) and the unitary matrix U such that

        U^\dagger H U = diag(m_i^2),

    and in the flavor basis the left-rotation is U_L = U.
    """
    H = M.conj().T @ M
    vals, vecs = np.linalg.eigh(H)
    # sort by eigenvalue
    idx = np.argsort(vals)
    vals = vals[idx]
    vecs = vecs[:, idx]
    masses = np.sqrt(np.clip(vals, 0, None))
    U = vecs  # columns are eigenvectors
    return masses, U


def ckm_or_pmns(U_up_like, U_down_like):
    """
    Given left-unitary matrices U_f (from diagonalizing M_f),
    build the relative mixing matrix

        V = U_up^\dagger U_down

    CKM:  V = U_u^\dagger U_d
    PMNS: U = U_e^\dagger U_nu
    """
    return U_up_like.conj().T @ U_down_like


def standard_mixing_angles(U):
    """
    Extract PDG-like mixing angles (theta12, theta23, theta13) and
    an estimate of delta_CP from a 3x3 unitary matrix U.

    Angles are returned in degrees.
    """
    # s13 = |U_13|
    s13 = abs(U[0, 2])
    theta13 = np.arcsin(np.clip(s13, 0.0, 1.0))

    c13 = np.cos(theta13)
    # s12 = |U_12| / c13
    s12 = abs(U[0, 1]) / max(c13, 1e-12)
    s12 = np.clip(s12, 0.0, 1.0)
    theta12 = np.arcsin(s12)

    # s23 = |U_23| / c13
    s23 = abs(U[1, 2]) / max(c13, 1e-12)
    s23 = np.clip(s23, 0.0, 1.0)
    theta23 = np.arcsin(s23)

    # Jarlskog invariant: J = Im(U11 U22 U12* U21*)
    J = np.imag(U[0, 0]*U[1, 1]*np.conj(U[0, 1])*np.conj(U[1, 0]))

    c12 = np.cos(theta12)
    c23 = np.cos(theta23)
    denom = (s12*c12*s23*c23*s13*c13*c13)
    if abs(denom) < 1e-12:
        delta = 0.0
    else:
        sdelta = np.clip(J/denom, -1.0, 1.0)
        delta = np.arcsin(sdelta)

    deg = 180.0/np.pi
    return theta12*deg, theta23*deg, theta13*deg, delta*deg


# ===========================
# 4. Seesaw reduction for 9x9 -> 3x3
# ===========================

def seesaw_effective_block(M_full, n_light=3):
    """
    Given a full (n_light + n_heavy) x (n_light + n_heavy) matrix M_full
    in block form

        [ M_LL  M_LH ]
        [ M_HL  M_HH ]

    return the effective light 3x3 block after integrating out the heavy modes:

        M_eff = M_LL - M_LH M_HH^{-1} M_HL

    This is the standard Schur-complement / seesaw formula.
    """
    n = M_full.shape[0]
    assert M_full.shape[1] == n
    n_heavy = n - n_light
    M_LL = M_full[:n_light, :n_light]
    M_LH = M_full[:n_light, n_light:]
    M_HL = M_full[n_light:, :n_light]
    M_HH = M_full[n_light:, n_light:]

    M_HH_inv = np.linalg.inv(M_HH)
    M_eff = M_LL - M_LH @ M_HH_inv @ M_HL
    return M_eff


# ===========================
# 5. Random proto-data and a simple driver
# ===========================

def random_proto_matrix(n, complex_entries=False, scale=1.0, rng=None):
    """
    Generate a random n x n proto-matrix with O(1) entries.

    complex_entries : if True, generates complex entries with random phases.
    """
    if rng is None:
        rng = np.random.default_rng()
    if complex_entries:
        mag = rng.normal(scale=scale, size=(n, n))
        phase = rng.uniform(0.0, 2*np.pi, size=(n, n))
        return mag * np.exp(1j*phase)
    else:
        return rng.normal(scale=scale, size=(n, n))


def run_single_sector(positions_light, kappa, forbidden_distance=None,
                      N_cycle=360, complex_entries=False, seed=None):
    """
    Minimal example: take the 3 light sites, build their 3x3 alignment kernel,
    align a random 3x3 proto-Yukawa, and extract mass ratios and mixing angles.

    This is the 'light-only' version; in the full 9-site model you would:

      - build K_9x9 on all sites
      - align a 9x9 proto-matrix
      - integrate out heavy states via seesaw
      - then restrict to the 3x3 light block

    but the interface is essentially the same.
    """
    if seed is not None:
        rng = np.random.default_rng(seed)
    else:
        rng = np.random.default_rng()

    # Build 3x3 K for the light sites alone
    K_light = build_alignment_kernel(
        positions_light,
        kappa=kappa,
        forbidden_distance=forbidden_distance,
        N=N_cycle
    )

    # Proto-Yukawa in the light 3x3 sector
    Y0 = random_proto_matrix(3, complex_entries=complex_entries, rng=rng)

    # Alignment map
    Y = align_matrix(Y0, K_light)

    # Masses and mixing (for a single sector we just get the left rotation)
    masses, U = hermitian_diagonalization(Y)
    ratios = masses / masses[-1]  # normalize to the heaviest

    return {
        "K_light": K_light,
        "Y_proto": Y0,
        "Y_aligned": Y,
        "masses": masses,
        "mass_ratios": ratios,
        "U_left": U
    }


def example_flavor_pipeline():
    """
    Example of a *full* CKM / PMNS pipeline:

      1. Choose internal positions for the light sites (1,2,5 here).
      2. Choose kappa and forbidden distance.
      3. Build K_light.
      4. Align independent proto-Yukawas for up, down, charged lepton, neutrino.
      5. Diagonalize each; compute CKM and PMNS.
      6. Extract mixing angles and CP phase.
    """
    positions_light = [1, 2, 5]
    kappa = 0.24
    forbidden_distance = 7  # can be None if you only use the 3x3 light subgraph

    # For reproducibility
    rng = np.random.default_rng(42)

    # Build K_light once
    K_light = build_alignment_kernel(
        positions_light,
        kappa=kappa,
        forbidden_distance=forbidden_distance,
        N=360
    )

    # Up, down, charged lepton, neutrino proto-Yukawas
    Y0_u = random_proto_matrix(3, complex_entries=True, rng=rng)
    Y0_d = random_proto_matrix(3, complex_entries=True, rng=rng)
    Y0_e = random_proto_matrix(3, complex_entries=True, rng=rng)
    Y0_n = random_proto_matrix(3, complex_entries=True, rng=rng)

    # Alignment
    Y_u = align_matrix(Y0_u, K_light)
    Y_d = align_matrix(Y0_d, K_light)
    Y_e = align_matrix(Y0_e, K_light)
    Y_n = align_matrix(Y0_n, K_light)

    # Diagonalization
    m_u, U_u = hermitian_diagonalization(Y_u)
    m_d, U_d = hermitian_diagonalization(Y_d)
    m_e, U_e = hermitian_diagonalization(Y_e)
    m_n, U_n = hermitian_diagonalization(Y_n)

    # Normalize mass spectra
    ru = m_u / m_u[-1]
    rd = m_d / m_d[-1]
    re = m_e / m_e[-1]
    rn = m_n / m_n[-1]

    # CKM and PMNS
    V_ckm = ckm_or_pmns(U_u, U_d)
    U_pmns = ckm_or_pmns(U_e, U_n)

    # Extract angles
    th12_q, th23_q, th13_q, dcp_q = standard_mixing_angles(V_ckm)
    th12_l, th23_l, th13_l, dcp_l = standard_mixing_angles(U_pmns)

    results = {
        "K_light": K_light,
        "masses": {
            "up": m_u,
            "down": m_d,
            "charged_lepton": m_e,
            "neutrino": m_n,
        },
        "mass_ratios": {
            "up": ru,
            "down": rd,
            "charged_lepton": re,
            "neutrino": rn,
        },
        "CKM": V_ckm,
        "PMNS": U_pmns,
        "angles_quarks_deg": {
            "theta12": th12_q,
            "theta23": th23_q,
            "theta13": th13_q,
            "delta_CP": dcp_q,
        },
        "angles_leptons_deg": {
            "theta12": th12_l,
            "theta23": th23_l,
            "theta13": th13_l,
            "delta_CP": dcp_l,
        },
    }

    return results


if __name__ == "__main__":
    # Run the example pipeline and print a compact summary
    res = example_flavor_pipeline()

    print("=== Alignment kernel K_light (3x3) ===")
    print(np.array_str(res["K_light"], precision=6, suppress_small=True))

    print("\n=== Mass ratios (normalized to heaviest) ===")
    for sector, ratios in res["mass_ratios"].items():
        print(f"{sector:>15s} : {np.array_str(ratios, precision=4, suppress_small=True)}")

    print("\n=== Quark mixing angles (deg) ===")
    for k, v in res["angles_quarks_deg"].items():
        print(f"{k:>10s} = {v:7.3f}")

    print("\n=== Lepton mixing angles (deg) ===")
    for k, v in res["angles_leptons_deg"].items():
        print(f"{k:>10s} = {v:7.3f}")
