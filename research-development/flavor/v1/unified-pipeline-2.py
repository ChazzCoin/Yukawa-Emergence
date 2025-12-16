import numpy as np

# =========================
# Basic config
# =========================

v_HIGGS = 246.0              # GeV
Lambda_Maj = 7.0e13          # GeV (overall heavy RH neutrino scale)
kappa = 360.0 / 89.0         # geometric scale in paper
eps = 1.0 / kappa            # decay factor in (0,1), ≈ 0.247
LIGHT = slice(0, 3)
HEAVY = slice(3, 9)

# =========================
# Utilities
# =========================

def random_complex_matrix(shape, rng):
    X = rng.normal(size=shape)
    Y = rng.normal(size=shape)
    return X + 1j * Y

def normalize_by_largest_singular_value(M):
    s_max = np.linalg.svd(M, compute_uv=False)[0]
    if s_max == 0:
        return M
    return M / s_max

def generation_pattern(eps, exponents):
    a, b, c = exponents
    return np.array([eps**a, eps**b, eps**c], dtype=float)

def build_site_scales_from_generations(gen3):
    """
    gen3: array-like of length 3 (generation scales).
    Sites (0,3,6) -> gen3[0]
    Sites (1,4,7) -> gen3[1]
    Sites (2,5,8) -> gen3[2]
    """
    s = np.zeros(9, dtype=float)
    s[[0, 3, 6]] = gen3[0]
    s[[1, 4, 7]] = gen3[1]
    s[[2, 5, 8]] = gen3[2]
    return s

def random_weighted_proto(shape, rng, site_scales):
    """
    9x9 complex Gaussian with variance ~ s_i * s_j, then normalized.
    """
    N = shape[0]
    assert shape == (9, 9)
    var = np.zeros((N, N), dtype=float)
    for i in range(N):
        for j in range(N):
            var[i, j] = site_scales[i] * site_scales[j]
    # draw real and imag with variance var
    X = rng.normal(scale=np.sqrt(var), size=shape)
    Y = rng.normal(scale=np.sqrt(var), size=shape)
    M = X + 1j * Y
    return normalize_by_largest_singular_value(M)

# =========================
# Alignment kernel K (9x9)
# =========================

def build_alignment_kernel(eps, N=9, d_star=7):
    """
    K_ij = eps^{|i-j|} for |i-j| != d_star, |i-j|>0
         = 1 for |i-j| = 0
         = 0 for |i-j| = d_star.
    """
    K = np.zeros((N, N), dtype=float)
    for i in range(N):
        for j in range(N):
            d = abs(i - j)
            if d == 0:
                K[i, j] = 1.0
            elif d == d_star:
                K[i, j] = 0.0
            else:
                K[i, j] = eps**d
    return K

def apply_alignment(K, X):
    """Hadamard / Schur alignment: Φ(X) = K ⊙ X."""
    return K * X

# =========================
# Proto matrices (9x9)
# =========================

def generate_proto_matrices(seed,
                            use_site_hierarchy=True,
                            exponents_u=(4, 2, 0),
                            exponents_d=(3, 1, 0),
                            exponents_e=(4, 2, 0),
                            exponents_nu=(1, 0, 0)):
    rng = np.random.default_rng(seed)

    if use_site_hierarchy:
        gen_u = generation_pattern(eps, exponents_u)
        gen_d = generation_pattern(eps, exponents_d)
        gen_e = generation_pattern(eps, exponents_e)
        gen_nu = generation_pattern(eps, exponents_nu)
    else:
        gen_u = gen_d = gen_e = gen_nu = np.array([1.0, 1.0, 1.0])

    site_scales_u = build_site_scales_from_generations(gen_u)
    site_scales_d = build_site_scales_from_generations(gen_d)
    site_scales_e = build_site_scales_from_generations(gen_e)
    site_scales_nu = build_site_scales_from_generations(gen_nu)

    Yu0 = random_weighted_proto((9, 9), rng, site_scales_u)
    Yd0 = random_weighted_proto((9, 9), rng, site_scales_d)
    Ye0 = random_weighted_proto((9, 9), rng, site_scales_e)
    Ynu0 = random_weighted_proto((9, 9), rng, site_scales_nu)

    # Majorana proto: symmetric, O(1), no extra site hierarchy
    M0 = random_complex_matrix((9, 9), rng)
    M0 = normalize_by_largest_singular_value(M0)
    M0 = 0.5 * (M0 + M0.T)

    return Yu0, Yd0, Ye0, Ynu0, M0

def align_all(K, Yu0, Yd0, Ye0, Ynu0, M0):
    Yu9 = apply_alignment(K, Yu0)
    Yd9 = apply_alignment(K, Yd0)
    Ye9 = apply_alignment(K, Ye0)
    Ynu9 = apply_alignment(K, Ynu0)
    M9 = apply_alignment(K, M0)
    return Yu9, Yd9, Ye9, Ynu9, M9

# =========================
# Schur complement 9→3 (Dirac)
# =========================

def schur_9_to_3(Y9, cond_tol=1e12):
    """
    Light sites 0,1,2; heavy 3..8.
    Y_eff = A - B D^{-1} B† (or pseudo-inverse if ill-conditioned).
    """
    A = Y9[LIGHT, LIGHT]
    B = Y9[LIGHT, HEAVY]
    D = Y9[HEAVY, HEAVY]

    # condition number in 2-norm
    s = np.linalg.svd(D, compute_uv=False)
    cond = s[0] / max(s[-1], 1e-18)

    if cond > cond_tol:
        D_inv = np.linalg.pinv(D)
        Y_eff = A - B @ D_inv @ B.conj().T
    else:
        # solve D X = B†
        X = np.linalg.solve(D, B.conj().T)
        Y_eff = A - B @ X
    return Y_eff

def reduce_all_dirac(Yu9, Yd9, Ye9, Ynu9):
    Yu_eff = schur_9_to_3(Yu9)
    Yd_eff = schur_9_to_3(Yd9)
    Ye_eff = schur_9_to_3(Ye9)
    Ynu_eff = schur_9_to_3(Ynu9)
    return Yu_eff, Yd_eff, Ye_eff, Ynu_eff

# =========================
# Majorana: triadic heavy projection and seesaw
# =========================

def heavy_block(M9):
    return M9[HEAVY, HEAVY]   # 6x6

def triad_heavy_basis(Nh=6, ks=(1, 2, 3)):
    """
    Nh x len(ks) matrix with orthonormal columns built from DFT modes.
    """
    i = np.arange(Nh)
    cols = []
    for k in ks:
        vec = np.exp(2j * np.pi * k * i / Nh)
        vec = vec / np.linalg.norm(vec)
        cols.append(vec)
    return np.stack(cols, axis=1)   # 6 x 3

def build_M_R_triadic(M9_aligned, Lambda_Maj, ks=(1, 2, 3)):
    M_H = heavy_block(M9_aligned)      # 6x6
    B_H = triad_heavy_basis(6, ks)     # 6x3
    M3 = B_H.conj().T @ M_H @ B_H      # 3x3
    M3 = 0.5 * (M3 + M3.T)             # enforce symmetry
    return Lambda_Maj * M3

def seesaw_light_neutrinos(Ynu_eff, M_R, v=v_HIGGS, cond_tol=1e12):
    """
    m_D = v / sqrt(2) * Ynu_eff
    m_ν = - m_D M_R^{-1} m_D^T
    """
    m_D = (v / np.sqrt(2.0)) * Ynu_eff
    s = np.linalg.svd(M_R, compute_uv=False)
    cond = s[0] / max(s[-1], 1e-18)
    if cond > cond_tol:
        M_R_inv = np.linalg.pinv(M_R)
    else:
        M_R_inv = np.linalg.inv(M_R)
    m_nu = - m_D @ M_R_inv @ m_D.T
    # enforce symmetry numerically
    m_nu = 0.5 * (m_nu + m_nu.T)
    return m_nu

# =========================
# Diagonalization and mixing
# =========================

def diag_dirac_Y(Y, v=v_HIGGS):
    U_L, s, U_Rh = np.linalg.svd(Y)
    masses = (v / np.sqrt(2.0)) * s
    return U_L, np.diag(s), U_Rh, masses

def takagi_symmetric(m):
    U, s, Vh = np.linalg.svd(m)
    return U, s

def diagonalize_all(Yu, Yd, Ye, mnu, v=v_HIGGS):
    UuL, Yu_diag, UuR, mu = diag_dirac_Y(Yu, v)
    UdL, Yd_diag, UdR, md = diag_dirac_Y(Yd, v)
    UlL, Ye_diag, UlR, me = diag_dirac_Y(Ye, v)
    U_nu, mnu_vals = takagi_symmetric(mnu)

    Vckm = UuL.conj().T @ UdL
    Vpmns = UlL.conj().T @ U_nu

    return mu, md, me, mnu_vals, Vckm, Vpmns

def extract_angles_and_phase(V):
    """
    CKM/PMNS-like parameterization (approx).
    Returns (theta12, theta23, theta13, delta) in radians.
    """
    s13 = np.abs(V[0, 2])
    theta13 = np.arcsin(np.clip(s13, 0.0, 1.0))

    s12 = np.abs(V[0, 1])
    c12 = np.abs(V[0, 0])
    theta12 = np.arctan2(s12, c12)

    s23 = np.abs(V[1, 2])
    c23 = np.abs(V[2, 2])
    theta23 = np.arctan2(s23, c23)

    # Jarlskog invariant
    J = np.imag(V[0, 0] * V[1, 1] * np.conj(V[0, 1]) * np.conj(V[1, 0]))
    denom = (np.sin(2 * theta12) * np.sin(2 * theta23) *
             np.sin(2 * theta13) * np.cos(theta13))

    if np.abs(denom) < 1e-12:
        delta = 0.0
    else:
        x = 8 * J / denom
        x = np.clip(x, -1.0, 1.0)
        delta = np.arcsin(x)

    return theta12, theta23, theta13, delta

# =========================
# Top-level pipeline (no RGE)
# =========================

def run_pipeline_minimal(seed=0,
                         exponents_u=(4, 2, 0),
                         exponents_d=(3, 1, 0),
                         exponents_e=(4, 2, 0),
                         exponents_nu=(1, 0, 0),
                         triad_ks=(1, 2, 3),
                         use_site_hierarchy=True):
    # 1. Alignment kernel
    K = build_alignment_kernel(eps, N=9, d_star=7)

    # 2. Proto matrices
    Yu0, Yd0, Ye0, Ynu0, M0 = generate_proto_matrices(
        seed,
        use_site_hierarchy=use_site_hierarchy,
        exponents_u=exponents_u,
        exponents_d=exponents_d,
        exponents_e=exponents_e,
        exponents_nu=exponents_nu,
    )

    # 3. Alignment (Schur)
    Yu9, Yd9, Ye9, Ynu9, M9 = align_all(K, Yu0, Yd0, Ye0, Ynu0, M0)

    # 4. Schur complement 9→3
    Yu_eff, Yd_eff, Ye_eff, Ynu_eff = reduce_all_dirac(Yu9, Yd9, Ye9, Ynu9)

    # 5. Majorana sector: triadic heavy projection
    M_R = build_M_R_triadic(M9, Lambda_Maj, ks=triad_ks)

    # 6. Seesaw → light neutrinos
    m_nu = seesaw_light_neutrinos(Ynu_eff, M_R, v=v_HIGGS)

    # 7. Diagonalize
    mu, md, me, mnu_vals, Vckm, Vpmns = diagonalize_all(
        Yu_eff, Yd_eff, Ye_eff, m_nu, v=v_HIGGS
    )

    # sort masses (for ratios)
    mu_sorted = np.sort(mu)
    md_sorted = np.sort(md)
    me_sorted = np.sort(me)
    mnu_sorted = np.sort(mnu_vals)

    # mixing angles
    th12_q, th23_q, th13_q, delta_q = extract_angles_and_phase(Vckm)
    th12_l, th23_l, th13_l, delta_l = extract_angles_and_phase(Vpmns)

    results = {
        "mu": mu_sorted,
        "md": md_sorted,
        "me": me_sorted,
        "mnu": mnu_sorted,
        "Vckm": Vckm,
        "Vpmns": Vpmns,
        "th_q": (th12_q, th23_q, th13_q, delta_q),
        "th_l": (th12_l, th23_l, th13_l, delta_l),
    }
    return results

if __name__ == "__main__":
    res = run_pipeline_minimal(seed=0)
    mu, md, me, mnu = res["mu"], res["md"], res["me"], res["mnu"]
    th12_q, th23_q, th13_q, delta_q = res["th_q"]
    th12_l, th23_l, th13_l, delta_l = res["th_l"]

    print("=== Mass ratios (normalized to heaviest) ===")
    print("up:",   mu / mu[-1])
    print("down:", md / md[-1])
    print("lep:",  me / me[-1])
    print("nu:",   mnu / mnu[-1])

    print("\n=== Quark mixing angles (deg) ===")
    print("theta12 =", np.degrees(th12_q))
    print("theta23 =", np.degrees(th23_q))
    print("theta13 =", np.degrees(th13_q))
    print("delta_CP (q) =", np.degrees(delta_q))

    print("\n=== Lepton mixing angles (deg) ===")
    print("theta12 =", np.degrees(th12_l))
    print("theta23 =", np.degrees(th23_l))
    print("theta13 =", np.degrees(th13_l))
    print("delta_CP (ℓ) =", np.degrees(delta_l))
