import numpy as np
import math

"""
Operator-first quasi-crystal flavor toy (geometry-driven mixing)
================================================================

Upgrades vs previous version:
  - Internal space: same 2D Fibonacci quasi-crystal patch (Cartesian product of
    two 1D Fibonacci chains).
  - We now use the *eigenvectors* of the internal Laplacian to define:
      * A 3D "generation" subspace from the three smallest nonzero eigenmodes.
      * Sector-dependent right-handed flavor bases U_R^s built from how simple,
        purely geometric region masks project into that subspace.
  - No F3, no hand-chosen 30° rotations: mixing arises from geometry + charges.

Flavor construction (axiom-driven):
  - Generation subspace H_gen spanned by v1, v2, v3 (lowest nonzero eigenmodes).
  - Spectral kernel F(lambda) = exp(-alpha * lambda) on these modes.
  - Discrete charge operator Q_s: integer charges q_{s,g} per sector & generation.
  - For each sector s:
      * define 3 geometrical region masks on the internal graph,
      * project them into H_gen and orthonormalize to get U_R^s (3x3 unitary),
      * build Yukawa Y_s = F_s_diag * U_R^s, with F_s_diag = diag(F_base * e^-beta q_s).
  - Diagonalize Y_s via SVD to get singular values (mass hierarchies) and
    left-handed unitaries (mixing matrices).
  - Compute a simple chi^2 vs rough SM-like targets.

Still a toy, but now *all* flavor structure comes from:
  - internal quasi-crystal geometry,
  - integer charges,
  - universal spectral kernel.
"""

# ----------------------------------------------------------------------
# 1. 1D Fibonacci chains and Laplacians
# ----------------------------------------------------------------------

def fibonacci_word(n_iter: int = 6, seed: str = "A") -> str:
    """
    Generate a Fibonacci substitution word:
      A -> AB
      B -> A
    After n_iter substitutions starting from 'seed'.
    """
    s = seed
    for _ in range(n_iter):
        s = "".join("AB" if c == "A" else "A" for c in s)
    return s

def fibonacci_chain_adjacency(word: str) -> np.ndarray:
    """
    Build adjacency matrix for a 1D Fibonacci chain with nearest-neighbor edges.
    Edge weights depend slightly on local A/B pattern to encode a weak
    quasi-crystal modulation along the chain.

    Nodes are indexed 0..N-1 corresponding to characters in 'word'.
    """
    N = len(word)
    A = np.zeros((N, N), dtype=float)

    # base connectivity
    for i in range(N - 1):
        A[i, i+1] = A[i+1, i] = 1.0

    # modulate weights by local pattern (AA, AB, BA)
    for i in range(N - 1):
        pair = word[i:i+2]
        if pair == "AA":
            w = 1.0
        elif pair == "AB":
            w = 1.1
        elif pair == "BA":
            w = 0.9
        else:  # "BB" (rare)
            w = 1.0
        A[i, i+1] = A[i+1, i] = w

    return A

def laplacian_from_adjacency(A: np.ndarray) -> np.ndarray:
    D = np.diag(A.sum(axis=1))
    return D - A

# Build 1D Fibonacci chains along x and y
word_x = fibonacci_word(6)   # length ~ 21
word_y = fibonacci_word(5)   # length ~ 13

A_x = fibonacci_chain_adjacency(word_x)
A_y = fibonacci_chain_adjacency(word_y)

L_x = laplacian_from_adjacency(A_x)
L_y = laplacian_from_adjacency(A_y)

N_x = A_x.shape[0]
N_y = A_y.shape[0]

# ----------------------------------------------------------------------
# 2. 2D quasi-crystal patch: Cartesian product of chains
# ----------------------------------------------------------------------

def cartesian_product_adjacency(A1: np.ndarray, A2: np.ndarray) -> np.ndarray:
    """
    Adjacency of the Cartesian product graph G = G1 □ G2.
    Nodes are pairs (i, j). Two nodes (i1, j1) and (i2, j2) are adjacent if:
      - i1 == i2 and j1, j2 are adjacent in G2, OR
      - j1 == j2 and i1, i2 are adjacent in G1.

    Edge weights inherit from 1D chains.
    """
    N1 = A1.shape[0]
    N2 = A2.shape[0]
    N  = N1 * N2
    A  = np.zeros((N, N), dtype=float)

    def idx(i1: int, i2: int) -> int:
        return i1 * N2 + i2

    # edges along x (vary i1, fixed i2)
    for i1 in range(N1):
        for i1p in range(N1):
            if A1[i1, i1p] != 0.0:
                for j in range(N2):
                    u = idx(i1,  j)
                    v = idx(i1p, j)
                    w = A1[i1, i1p]
                    A[u, v] = A[v, u] = max(A[u, v], w)

    # edges along y (vary i2, fixed i1)
    for j1 in range(N2):
        for j2 in range(N2):
            if A2[j1, j2] != 0.0:
                for i in range(N1):
                    u = idx(i, j1)
                    v = idx(i, j2)
                    w = A2[j1, j2]
                    A[u, v] = A[v, u] = max(A[u, v], w)

    return A

A_int = cartesian_product_adjacency(A_x, A_y)
L_int = laplacian_from_adjacency(A_int)
N_sites = A_int.shape[0]

# Diagonalize internal Laplacian
eigvals, eigvecs = np.linalg.eigh(L_int)

# ----------------------------------------------------------------------
# 3. Generation subspace from eigenvectors
# ----------------------------------------------------------------------

# Skip near-zero mode, take three smallest nonzero eigenvalues as "generations"
eps = 1e-10
nonzero_indices = np.where(eigvals > eps)[0]
gen_indices = nonzero_indices[:3]         # indices of λ1, λ2, λ3
lam_gen = eigvals[gen_indices]           # shape (3,)
V_gen   = eigvecs[:, gen_indices]        # shape (N_sites, 3); columns are eigenvectors

# Orthonormalize V_gen just in case (should already be orthonormal from eigh)
# but we can QR it to be safe
Q_gen, _ = np.linalg.qr(V_gen)
V_gen = Q_gen  # use orthonormalized generation basis

# ----------------------------------------------------------------------
# 4. Spectral kernel F(lambda) from quasi-crystal
# ----------------------------------------------------------------------

def base_kernel(lams: np.ndarray, alpha: float = 3.0) -> np.ndarray:
    """
    Spectral kernel from the internal quasi-crystal:
        F(lambda) = exp(-alpha * lambda).
    """
    return np.exp(-alpha * lams)

alpha = 3.0  # fixed, not fitted
F_base = base_kernel(lam_gen, alpha=alpha)  # shape (3,)

# ----------------------------------------------------------------------
# 5. Integer charges Q_s: sector + generation hierarchy
# ----------------------------------------------------------------------

beta = 1.0  # fixed, not fitted

sector_charges_gen = {
    "u":  np.array([2.0, 1.0, 0.0], dtype=float),  # u, c, t
    "d":  np.array([3.0, 2.0, 1.0], dtype=float),  # d, s, b
    "e":  np.array([4.0, 3.0, 2.0], dtype=float),  # e, mu, tau
    "nu": np.array([6.0, 5.0, 4.0], dtype=float),  # nu1, nu2, nu3
}

def sector_weights(F_base: np.ndarray, q_vec: np.ndarray, beta: float = 1.0) -> np.ndarray:
    """
    F_s(g) = F_base(g) * exp(-beta * q_vec[g]).
    """
    return F_base * np.exp(-beta * q_vec)

# ----------------------------------------------------------------------
# 6. Geometry-based region masks and U_R^s from eigenvectors
# ----------------------------------------------------------------------

def site_coords(n: int, Ny: int) -> tuple[int, int]:
    """Return (ix, iy) for flattened index n in a grid of width Ny."""
    ix = n // Ny
    iy = n % Ny
    return ix, iy

def build_region_masks_LR(Nx: int, Ny: int, N_sites: int):
    """
    Build simple geometric masks for LEFT and RIGHT flavor spaces.

    LEFT masks (V^L): emphasize different coarse regions of the patch:
      - L0: left half
      - L1: right half
      - L2: central band

    RIGHT masks (V^R): emphasize different, shifted regions:
      - R0: bottom half
      - R1: top half
      - R2: diagonal-ish strip
    """
    def site_coords(n: int, Ny: int):
        ix = n // Ny
        iy = n % Ny
        return ix, iy

    def make_mask_left_half():
        m = np.zeros(N_sites)
        for n in range(N_sites):
            ix, iy = site_coords(n, Ny)
            if ix < Nx // 2:
                m[n] = 1.0
        return m

    def make_mask_right_half():
        m = np.zeros(N_sites)
        for n in range(N_sites):
            ix, iy = site_coords(n, Ny)
            if ix >= Nx // 2:
                m[n] = 1.0
        return m

    def make_mask_central_band():
        m = np.zeros(N_sites)
        for n in range(N_sites):
            ix, iy = site_coords(n, Ny)
            if Nx // 3 <= ix < 2 * Nx // 3:
                m[n] = 1.0
        return m

    def make_mask_bottom_half():
        m = np.zeros(N_sites)
        for n in range(N_sites):
            ix, iy = site_coords(n, Ny)
            if iy < Ny // 2:
                m[n] = 1.0
        return m

    def make_mask_top_half():
        m = np.zeros(N_sites)
        for n in range(N_sites):
            ix, iy = site_coords(n, Ny)
            if iy >= Ny // 2:
                m[n] = 1.0
        return m

    def make_mask_diagonal_strip():
        m = np.zeros(N_sites)
        for n in range(N_sites):
            ix, iy = site_coords(n, Ny)
            # crude diagonal condition: |ix - scaled iy| small
            if abs(ix - (Nx * iy / max(1, Ny))) < Nx / 4:
                m[n] = 1.0
        return m

    # Base masks
    L_masks_base = [
        make_mask_left_half(),
        make_mask_right_half(),
        make_mask_central_band(),
    ]
    R_masks_base = [
        make_mask_bottom_half(),
        make_mask_top_half(),
        make_mask_diagonal_strip(),
    ]

    # Assign different permutations per sector
    masks_L = {}
    masks_R = {}

    # Up: (L0,L1,L2), (R0,R1,R2)
    masks_L["u"]  = [L_masks_base[0], L_masks_base[1], L_masks_base[2]]
    masks_R["u"]  = [R_masks_base[0], R_masks_base[1], R_masks_base[2]]

    # Down: permute right masks
    masks_L["d"]  = [L_masks_base[1], L_masks_base[2], L_masks_base[0]]
    masks_R["d"]  = [R_masks_base[1], R_masks_base[2], R_masks_base[0]]

    # Charged leptons: permute left masks
    masks_L["e"]  = [L_masks_base[2], L_masks_base[0], L_masks_base[1]]
    masks_R["e"]  = [R_masks_base[0], R_masks_base[2], R_masks_base[1]]

    # Neutrinos: swap roles more strongly
    masks_L["nu"] = [L_masks_base[1], L_masks_base[0], L_masks_base[2]]
    masks_R["nu"] = [R_masks_base[2], R_masks_base[1], R_masks_base[0]]

    return masks_L, masks_R

mask_L, mask_R = build_region_masks_LR(N_x, N_y, N_sites)

def normalize_column(v: np.ndarray) -> np.ndarray:
    nrm = np.linalg.norm(v)
    if nrm < 1e-12:
        return v
    return v / nrm

def geometry_unitaries_LR_from_eigenvectors(V_gen: np.ndarray, masks_L, masks_R):
    """
    For each sector s:
      - LEFT: project its 3 L-masks into generation subspace to get a 3x3 matrix,
              QR-orthonormalize → U_L^s.
      - RIGHT: same with R-masks → U_R^s.
    """
    sectors = {}
    for s in masks_L.keys():
        # LEFT
        cols_L = [normalize_column(m.astype(float)) for m in masks_L[s]]
        W_L = np.stack(cols_L, axis=1)            # (N_sites, 3)
        C_L = V_gen.conj().T @ W_L                # (3,3)
        Q_L, _ = np.linalg.qr(C_L)
        U_L_s = Q_L.astype(complex)

        # RIGHT
        cols_R = [normalize_column(m.astype(float)) for m in masks_R[s]]
        W_R = np.stack(cols_R, axis=1)            # (N_sites, 3)
        C_R = V_gen.conj().T @ W_R                # (3,3)
        Q_R, _ = np.linalg.qr(C_R)
        U_R_s = Q_R.astype(complex)

        sectors[s] = (U_L_s, U_R_s)
    return sectors

sector_unitaries_geom = geometry_unitaries_LR_from_eigenvectors(
    V_gen, mask_L, mask_R
)

# ----------------------------------------------------------------------
# 7. Build Yukawa matrices from operators (geometry-driven U_R)
# ----------------------------------------------------------------------

def build_yukawas(F_base: np.ndarray, sector_charges_gen, beta: float = 1.0):
    """
    Build Yukawa matrices for all sectors:
        F_s_diag = diag(F_base(g) * exp(-beta * q_{s,g}))
        Y_s      = U_L^s†  F_s_diag  U_R^s
    Here both U_L^s and U_R^s are geometry-derived.
    """
    Y = {}
    for name, (U_L, U_R) in sector_unitaries_geom.items():
        q_vec = sector_charges_gen[name]
        weights = sector_weights(F_base, q_vec, beta=beta)
        F_s_diag = np.diag(weights.astype(complex))
        Y[name] = U_L.conj().T @ F_s_diag @ U_R
    return Y

Yukawas = build_yukawas(F_base, sector_charges_gen, beta=beta)

# ----------------------------------------------------------------------
# 8. Diagonalization and mixing
# ----------------------------------------------------------------------

def diagonalize_dirac(Y: np.ndarray):
    """
    Dirac-like Yukawa diagonalization via SVD:
        Y = U_L diag(s) U_R†
    Returns U_L, singular values s, U_R.
    """
    U_L, s_vals, U_Rh = np.linalg.svd(Y)
    U_R = U_Rh.conj().T
    return U_L, s_vals, U_R

Uu_L, su, Uu_R = diagonalize_dirac(Yukawas["u"])
Ud_L, sd, Ud_R = diagonalize_dirac(Yukawas["d"])
Ue_L, se, Ue_R = diagonalize_dirac(Yukawas["e"])
Un_L, sn, Un_R = diagonalize_dirac(Yukawas["nu"])

def mixing_matrix(U_L_up: np.ndarray, U_L_down: np.ndarray) -> np.ndarray:
    """CKM/PMNS-like mixing matrix: V = U_L_up† U_L_down."""
    return U_L_up.conj().T @ U_L_down

def mixing_angles_from_U(U: np.ndarray):
    """
    Extract approximate (theta12, theta23, theta13) from a unitary 3x3 matrix U
    using PDG-like conventions on |U|, ignoring CP phases.
    """
    U_abs = np.abs(U)
    s13 = U_abs[0, 2]
    c13 = math.sqrt(max(0.0, 1.0 - s13**2))
    if c13 < 1e-12:
        return 0.0, 0.0, math.pi / 2.0
    s12 = U_abs[0, 1] / c13
    s23 = U_abs[1, 2] / c13
    s12 = max(-1.0, min(1.0, s12))
    s23 = max(-1.0, min(1.0, s23))
    theta12 = math.asin(s12)
    theta23 = math.asin(s23)
    theta13 = math.asin(s13)
    return theta12, theta23, theta13

V_ckm  = mixing_matrix(Uu_L, Ud_L)
U_pmns = mixing_matrix(Ue_L, Un_L)

theta12_q, theta23_q, theta13_q = mixing_angles_from_U(V_ckm)
theta12_l, theta23_l, theta13_l = mixing_angles_from_U(U_pmns)

# ----------------------------------------------------------------------
# 9. Compare to rough SM targets via chi^2
# ----------------------------------------------------------------------

def compute_observables(su, sd, se, sn,
                        theta12_q, theta23_q, theta13_q,
                        theta12_l, theta23_l, theta13_l):
    """
    Build a dictionary of dimensionless observables:
      - mass ratios m1/m3, m2/m3 per sector (using singular values)
      - mixing angles (radians) for CKM- and PMNS-like matrices
    """
    def ratios(s):
        s_sorted = np.sort(s)  # ascending: [light, mid, heavy]
        m1, m2, m3 = s_sorted
        return m1 / m3, m2 / m3

    mu_mt, mc_mt   = ratios(su)
    md_mb, ms_mb   = ratios(sd)
    me_mt, mmu_mt  = ratios(se)

    return {
        "mu/mt":      mu_mt,
        "mc/mt":      mc_mt,
        "md/mb":      md_mb,
        "ms/mb":      ms_mb,
        "me/mtau":    me_mt,
        "mmu/mtau":   mmu_mt,
        "theta12_q":  theta12_q,
        "theta23_q":  theta23_q,
        "theta13_q":  theta13_q,
        "theta12_l":  theta12_l,
        "theta23_l":  theta23_l,
        "theta13_l":  theta13_l,
    }

targets = {
    "mu/mt":     2.2e-5,
    "mc/mt":     7.5e-3,
    "md/mb":     1.1e-3,
    "ms/mb":     2.2e-2,
    "me/mtau":   2.9e-4,
    "mmu/mtau":  0.059,
    "theta12_q": 0.227,
    "theta23_q": 0.041,
    "theta13_q": 0.0036,
    "theta12_l": 0.584,
    "theta23_l": 0.785,
    "theta13_l": 0.150,
}

def chi2(observables, targets):
    """
    Simple chi^2:
      - for ratios, use log10 error with sigma_log = 0.3 dex (~ factor 2)
      - for angles, use sigma = 0.2 rad
    """
    chi2_total = 0.0
    details = []

    ratio_keys = ["mu/mt", "mc/mt", "md/mb", "ms/mb", "me/mtau", "mmu/mtau"]
    angle_keys = ["theta12_q", "theta23_q", "theta13_q",
                  "theta12_l", "theta23_l", "theta13_l"]

    for k in ratio_keys:
        m = observables[k]
        t = targets[k]
        if m <= 0 or t <= 0:
            continue
        logm = math.log10(m)
        logt = math.log10(t)
        sigma_log = 0.3
        contrib = ((logm - logt) / sigma_log)**2
        chi2_total += contrib
        details.append((k, m, t, contrib))

    for k in angle_keys:
        m = observables[k]
        t = targets[k]
        sigma = 0.2
        contrib = ((m - t) / sigma)**2
        chi2_total += contrib
        details.append((k, m, t, contrib))

    return chi2_total, details

obs = compute_observables(su, sd, se, sn,
                          theta12_q, theta23_q, theta13_q,
                          theta12_l, theta23_l, theta13_l)
chi2_value, chi2_details = chi2(obs, targets)

# ----------------------------------------------------------------------
# 10. Print summary
# ----------------------------------------------------------------------

def main():
    print("=== 2D Fibonacci quasi-crystal internal graph ===")
    print(f"Chain lengths: Nx = {N_x}, Ny = {N_y}, total sites = {N_sites}")
    print("First 10 Laplacian eigenvalues (internal L_int):")
    print(eigvals[:10])
    print()
    print("Generation eigenvalue indices:", gen_indices)
    print("Generation triad lam_gen:", lam_gen)
    print("Base kernel F_base(lam_gen):", F_base)
    print()

    print("=== Yukawa singular values (up to overall scale) ===")
    print("Up-type (su):        ", su)
    print("Down-type (sd):      ", sd)
    print("Charged leptons (se):", se)
    print("Neutrino Dirac (sn): ", sn)
    print()

    print("=== CKM-like mixing matrix (geometry-driven) ===")
    print(V_ckm)
    print("Mixing angles (radians):")
    print(f"theta12_q ≈ {theta12_q:.3f}, theta23_q ≈ {theta23_q:.3f}, theta13_q ≈ {theta13_q:.3e}")
    print()

    print("=== PMNS-like mixing matrix (geometry-driven) ===")
    print(U_pmns)
    print("Mixing angles (radians):")
    print(f"theta12_l ≈ {theta12_l:.3f}, theta23_l ≈ {theta23_l:.3f}, theta13_l ≈ {theta13_l:.3e}")
    print()

    print("=== Observables vs rough targets ===")
    for k, m, t, contrib in chi2_details:
        print(f"{k:12s}: model={m:.3e}, target={t:.3e}, chi2_contrib={contrib:.2f}")
    print()
    print(f"Total chi^2 ≈ {chi2_value:.2f}")
    print()

    print("NOTES:")
    print("- Generation subspace is spanned by the three lowest nonzero eigenmodes")
    print("  of the internal Laplacian; this is the 'R+L' spectral data in axiom form.")
    print("- Sector right-handed flavor bases U_R^s are built purely from geometry:")
    print("  we project simple region masks on the quasi-crystal onto the generation")
    print("  subspace and orthonormalize; no F3, no arbitrary rotations.")
    print("- Yukawa structure Y_s = diag(F_base * exp(-beta * q_{s,g})) * U_R^s.")
    print("- Mass hierarchies come from F_base(lambda) and integer charges q_{s,g}.")
    print("- Mixing comes entirely from the geometry-derived U_R^s and the subsequent")
    print("  SVD diagonalization; there are no continuous flavor-tuning parameters.")
    print("- The chi^2 will still be large, but now whatever mixing you see is a")
    print("  genuinely emergent property of the quasi-crystal + operator setup.")

if __name__ == "__main__":
    main()