#!/usr/bin/env python3
"""
Resonant 24-cell spectral flavor toy model.

Core principles (no cheating):
- ONE parent shape: the regular 24-cell in 4D (24 vertices).
- ONE geometric operator: the graph Laplacian Δ on the 24-cell's vertex graph.
- ONE universal kernel: K = exp(-Δ), same for all sectors.
- All flavor structure (up, down, charged lepton, neutrino) emerges from:
    * the spectrum (eigenvalues + eigenvectors) of Δ, and
    * discrete choices of *which eigenvalue clusters* define left/right subspaces.

NO:
- Random matrices,
- Sector-specific scaling parameters,
- Hand-tuned exponent tables,
- Continuous fit parameters.

Everything is determined by the geometry + pure linear algebra.
"""

import numpy as np
import math

# ---------------------------------------------------------------------------
# 1. 24-cell geometry: vertices and Laplacian
# ---------------------------------------------------------------------------

def build_24cell_vertices():
    """
    Construct the 24 vertices of the regular 24-cell in R^4.

    A standard coordinate realization:
    - 8 vertices of type A: (±1, 0, 0, 0) and permutations over coordinates.
    - 16 vertices of type B: (±1/2, ±1/2, ±1/2, ±1/2) for all 16 sign choices.

    This set is invariant under the symmetry group of the 24-cell and
    sits on a sphere in R^4.
    """
    verts = []

    # Type A: permutations of (±1, 0, 0, 0)
    for axis in range(4):
        for sign in (+1.0, -1.0):
            v = np.zeros(4)
            v[axis] = sign
            verts.append(v)

    # Type B: all sign combinations of (±1/2, ±1/2, ±1/2, ±1/2)
    for s0 in (+0.5, -0.5):
        for s1 in (+0.5, -0.5):
            for s2 in (+0.5, -0.5):
                for s3 in (+0.5, -0.5):
                    v = np.array([s0, s1, s2, s3])
                    verts.append(v)

    verts = np.array(verts)   # shape (24, 4)
    assert verts.shape == (24, 4)
    return verts


def build_24cell_adjacency(vertices, tol=1e-8):
    """
    Build the adjacency matrix A (24x24) of the 24-cell graph.

    Two vertices are connected by an edge if their Euclidean distance
    equals the minimal nonzero distance between any pair.

    This is purely geometric, no arbitrary thresholds beyond numerical tol.
    """
    N = vertices.shape[0]
    A = np.zeros((N, N), dtype=int)

    # Compute squared distances between all pairs
    d2 = np.zeros((N, N))
    for i in range(N):
        diff = vertices[i] - vertices
        d2[i, :] = np.sum(diff * diff, axis=1)

    # Find the smallest non-zero squared distance
    d2_flat = d2.flatten()
    nonzero = d2_flat[d2_flat > tol]
    d2_min = np.min(nonzero)

    # Connect vertices at minimal distance
    for i in range(N):
        for j in range(i+1, N):
            if abs(d2[i, j] - d2_min) < tol:
                A[i, j] = 1
                A[j, i] = 1

    return A


def build_laplacian(A):
    """
    Graph Laplacian L = D - A, where D is degree matrix.
    """
    degrees = np.sum(A, axis=1)
    D = np.diag(degrees)
    L = D - A
    return L


# ---------------------------------------------------------------------------
# 2. Spectral decomposition and universal kernel
# ---------------------------------------------------------------------------

def spectral_decomposition(L):
    """
    Diagonalize symmetric Laplacian:

        L v_i = λ_i v_i

    Returns:
        evals : eigenvalues sorted ascending (shape (N,))
        evecs : eigenvectors in columns (shape (N, N))
    """
    evals, evecs = np.linalg.eigh(L)
    # eigh already returns sorted evals for symmetric matrices
    return evals, evecs


def build_universal_kernel(evals, evecs):
    """
    Build universal kernel K = exp(-L) in the vertex basis.

    In spectral form:
        K = V diag(exp(-λ_i)) V^T

    with V columns = eigenvectors.

    This is the discrete heat kernel with "time" t=1, no free parameter.
    """
    f_vals = np.exp(-evals)
    # evecs: shape (N, N), columns are eigenvectors
    # K = V * diag(f_vals) * V^T
    K = (evecs * f_vals) @ evecs.T.conj()
    return K, f_vals


# ---------------------------------------------------------------------------
# 3. Eigenvalue clustering (degeneracies)
# ---------------------------------------------------------------------------

def cluster_eigenvalues(evals, tol=1e-8):
    """
    Group eigenvalues into clusters of (approximately) equal values.

    Returns:
        clusters: list of lists of indices, e.g.
                  [[0], [1,2,3], [4,5,6,7], ...]

    This is purely spectral: it reads off degeneracies from geometry.
    """
    clusters = []
    current_cluster = [0]

    for i in range(1, len(evals)):
        if abs(evals[i] - evals[i-1]) < tol:
            current_cluster.append(i)
        else:
            clusters.append(current_cluster)
            current_cluster = [i]
    clusters.append(current_cluster)

    return clusters


# ---------------------------------------------------------------------------
# 4. Build projectors from spectral clusters
# ---------------------------------------------------------------------------

def projector_from_cluster(evecs, cluster_indices, n_rows=3):
    """
    Build an n_rows × N projector from a given eigenvalue cluster:

        P[r, :] = eigenvector^T

    for the first n_rows eigenvectors in that cluster.

    If the cluster has fewer than n_rows eigenvectors, we pad by taking
    additional eigenvectors from the global spectrum (next indices).
    This is still purely deterministic and spectral.
    """
    N = evecs.shape[0]
    P = np.zeros((n_rows, N), dtype=complex)

    # Flatten cluster indices into an ordered list
    idx_list = list(cluster_indices)

    # If fewer than n_rows, pad with additional indices
    if len(idx_list) < n_rows:
        # Find all indices 0..N-1 not in cluster
        all_idx = list(range(N))
        remaining = [i for i in all_idx if i not in idx_list]
        # Append as many as needed
        idx_list = idx_list + remaining[:(n_rows - len(idx_list))]

    # Take first n_rows eigenvectors
    for r in range(n_rows):
        idx = idx_list[r]
        v = evecs[:, idx]          # shape (N,)
        P[r, :] = v.conj().T       # row = eigenvector^T
    return P


def build_sector_projectors(evals, evecs):
    """
    Construct LEFT/RIGHT projectors P_L, P_R for each sector (u,d,e,nu)
    from eigenvalue clusters.

    Strategy:
    - Cluster eigenvalues (degeneracies).
    - Ignore the trivial λ=0 ground state cluster for flavor (index 0 cluster).
    - Use the remaining clusters with size>=3 as natural "triplet" candidates.
    - For left-handed vs right-handed subspaces, use *different* clusters:
        - Up:
            L_u = lowest nontrivial triplet cluster
            R_u = highest triplet cluster
        - Down:
            L_d = 2nd lowest triplet cluster (if exists, else same as L_u)
            R_d = 2nd highest triplet cluster (if exists, else same as R_u)
        - Charged leptons:
            L_e = L_u
            R_e = R_d
        - Neutrinos:
            L_n = L_d
            R_n = R_u

    All of this is fully determined once evals/evecs are known.
    No continuous parameters, no randomness.
    """
    clusters = cluster_eigenvalues(evals)
    N = len(evals)

    # Identify triplet-like clusters (size >= 3), excluding the 0-eigenvalue cluster
    triplet_clusters = []
    for ci, cl in enumerate(clusters):
        if len(cl) >= 3:
            # Optional: skip λ=0 cluster (usually evals[0] ~ 0)
            if abs(evals[cl[0]]) < 1e-12:
                continue
            triplet_clusters.append(cl)

    if len(triplet_clusters) == 0:
        # Fallback: just use the lowest nontrivial cluster(s) whatever size
        # (still deterministic, but flavor structure will be degenerate)
        print("WARNING: No triplet clusters found; using smallest nontrivial clusters.")
        # skip λ=0 cluster
        nonzero_clusters = [cl for cl in clusters if abs(evals[cl[0]]) > 1e-12]
        # ensure at least one cluster
        if len(nonzero_clusters) == 0:
            nonzero_clusters = clusters
        triplet_clusters = nonzero_clusters

    # Helper to pick clusters with wrap-around if needed
    def pick_cluster(idx):
        return triplet_clusters[idx % len(triplet_clusters)]

    # Choose clusters for each sector (left/right) following the pattern above
    L_u_cl = pick_cluster(0)
    R_u_cl = pick_cluster(-1)

    L_d_cl = pick_cluster(1) if len(triplet_clusters) > 1 else L_u_cl
    R_d_cl = pick_cluster(-2) if len(triplet_clusters) > 1 else R_u_cl

    L_e_cl = L_u_cl
    R_e_cl = R_d_cl

    L_n_cl = L_d_cl
    R_n_cl = R_u_cl

    # Build projectors
    P_L_u = projector_from_cluster(evecs, L_u_cl, n_rows=3)
    P_R_u = projector_from_cluster(evecs, R_u_cl, n_rows=3)

    P_L_d = projector_from_cluster(evecs, L_d_cl, n_rows=3)
    P_R_d = projector_from_cluster(evecs, R_d_cl, n_rows=3)

    P_L_e = projector_from_cluster(evecs, L_e_cl, n_rows=3)
    P_R_e = projector_from_cluster(evecs, R_e_cl, n_rows=3)

    P_L_n = projector_from_cluster(evecs, L_n_cl, n_rows=3)
    P_R_n = projector_from_cluster(evecs, R_n_cl, n_rows=3)

    sector_proj = {
        "u":  (P_L_u, P_R_u),
        "d":  (P_L_d, P_R_d),
        "e":  (P_L_e, P_R_e),
        "nu": (P_L_n, P_R_n),
    }

    return sector_proj, clusters, triplet_clusters


# ---------------------------------------------------------------------------
# 5. Yukawas, diagonalization, mixing
# ---------------------------------------------------------------------------

def build_yukawa(P_L, P_R, K):
    """
    Yukawa matrix from geometry:

        Y = P_L @ K @ P_R^†

    P_L: 3×N, P_R: 3×N, K: N×N
    => Y: 3×3
    """
    return P_L @ K @ P_R.conj().T


def diagonalize_dirac(Y):
    """
    SVD for Dirac-like Yukawa:

        Y = U_L diag(s) U_R^†

    Returns:
        U_L, s_vals, U_R
    """
    U_L, s_vals, U_Rh = np.linalg.svd(Y)
    U_R = U_Rh.conj().T
    return U_L, s_vals, U_R


def mixing_matrix(U_L_up, U_L_down):
    """
    CKM/PMNS-like mixing matrix:

        V = U_L_up^† U_L_down
    """
    return U_L_up.conj().T @ U_L_down


def mixing_angles_from_U(U):
    """
    Extract approximate mixing angles (θ12, θ23, θ13) from a unitary 3×3 matrix U,
    ignoring CP phase, using standard PDG-like formulae on |U|:

        s13 = |U_13|
        c13 = sqrt(1 - s13^2)
        s12 = |U_12| / c13
        s23 = |U_23| / c13
    """
    U_abs = np.abs(U)
    s13 = U_abs[0, 2]
    c13 = math.sqrt(max(0.0, 1.0 - s13**2))
    if c13 < 1e-12:
        # pathological corner
        return 0.0, 0.0, math.pi / 2.0

    s12 = U_abs[0, 1] / c13
    s23 = U_abs[1, 2] / c13

    # Clamp to [-1,1] for safety
    s12 = max(-1.0, min(1.0, s12))
    s23 = max(-1.0, min(1.0, s23))

    theta12 = math.asin(s12)
    theta23 = math.asin(s23)
    theta13 = math.asin(s13)

    return theta12, theta23, theta13


# ---------------------------------------------------------------------------
# 6. Main: put it all together
# ---------------------------------------------------------------------------

def main():
    # 1) Geometry: 24-cell
    verts = build_24cell_vertices()
    A = build_24cell_adjacency(verts)
    L = build_laplacian(A)

    # 2) Spectrum & kernel
    evals, evecs = spectral_decomposition(L)
    K, f_vals = build_universal_kernel(evals, evecs)

    # 3) Spectral clusters
    clusters = cluster_eigenvalues(evals)
    sector_proj, all_clusters, triplet_clusters = build_sector_projectors(evals, evecs)

    print("=== 24-cell spectral data ===")
    print("Number of vertices:", verts.shape[0])
    print("Eigenvalues of Laplacian (sorted):")
    print(evals)
    print()

    print("Eigenvalue clusters (indices):")
    for i, cl in enumerate(all_clusters):
        lam = evals[cl[0]]
        print(f"  Cluster {i}: size={len(cl)}, λ≈{lam:.6f}, indices={cl}")
    print()

    print("Triplet-like clusters (size >= 3, excluding λ≈0):")
    for i, cl in enumerate(triplet_clusters):
        lam = evals[cl[0]]
        print(f"  Triplet cluster {i}: size={len(cl)}, λ≈{lam:.6f}, indices={cl}")
    print()

    # 4) Build Yukawas for each sector
    P_L_u, P_R_u = sector_proj["u"]
    P_L_d, P_R_d = sector_proj["d"]
    P_L_e, P_R_e = sector_proj["e"]
    P_L_n, P_R_n = sector_proj["nu"]

    Yu  = build_yukawa(P_L_u, P_R_u, K)
    Yd  = build_yukawa(P_L_d, P_R_d, K)
    Ye  = build_yukawa(P_L_e, P_R_e, K)
    Ynu = build_yukawa(P_L_n, P_R_n, K)

    # 5) Diagonalize Yukawas
    Uu_L, su, Uu_R = diagonalize_dirac(Yu)
    Ud_L, sd, Ud_R = diagonalize_dirac(Yd)
    Ue_L, se, Ue_R = diagonalize_dirac(Ye)
    Un_L, sn, Un_R = diagonalize_dirac(Ynu)

    print("=== Yukawa singular values (up to overall scale) ===")
    print("Up-type (su):        ", su)
    print("Down-type (sd):      ", sd)
    print("Charged leptons (se):", se)
    print("Neutrino Dirac (sn): ", sn)
    print()

    # 6) Mixing matrices
    V_ckm  = mixing_matrix(Uu_L, Ud_L)
    U_pmns = mixing_matrix(Ue_L, Un_L)

    theta12_q, theta23_q, theta13_q = mixing_angles_from_U(V_ckm)
    theta12_l, theta23_l, theta13_l = mixing_angles_from_U(U_pmns)

    print("=== CKM-like mixing matrix (quarks) ===")
    print(V_ckm)
    print("Approx CKM mixing angles (radians):")
    print(f"theta12_q ≈ {theta12_q:.3f}, theta23_q ≈ {theta23_q:.3f}, theta13_q ≈ {theta13_q:.3f}")
    print()

    print("=== PMNS-like mixing matrix (leptons) ===")
    print(U_pmns)
    print("Approx PMNS mixing angles (radians):")
    print(f"theta12_l ≈ {theta12_l:.3f}, theta23_l ≈ {theta23_l:.3f}, theta13_l ≈ {theta13_l:.3f}")
    print()

    print("NOTES:")
    print("- No randomness, no sector-specific scales, no hand-tuned exponents.")
    print("- Parent object: the 24-cell (24 vertices in R^4).")
    print("- Geometry → Laplacian spectrum λ_i, eigenvectors v_i.")
    print("- Universal kernel: K = exp(-Δ) built purely from {λ_i, v_i}.")
    print("- Flavor sectors: (u, d, e, ν) defined via spectral clusters (degeneracies).")
    print("- Left/right projectors pick 3D subspaces from different eigenvalue clusters.")
    print("- Yukawas: Y_s = P_L^(s) K P_R^(s)†.")
    print("- Mixing arises solely from misalignment of left subspaces under the same K.")
    print("- This is a testbed: not designed to match SM data, but to explore how")
    print("  a single spectral object (Δ on the 24-cell) can generate hierarchical")
    print("  patterns and mixing without any arbitrary continuous parameters.")


if __name__ == "__main__":
    main()

"""
RESULTS:
=== 24-cell spectral data ===
Number of vertices: 24
Eigenvalues of Laplacian (sorted):
[5.76557525e-16 4.00000000e+00 4.00000000e+00 4.00000000e+00
 4.00000000e+00 8.00000000e+00 8.00000000e+00 8.00000000e+00
 8.00000000e+00 8.00000000e+00 8.00000000e+00 8.00000000e+00
 8.00000000e+00 8.00000000e+00 1.00000000e+01 1.00000000e+01
 1.00000000e+01 1.00000000e+01 1.00000000e+01 1.00000000e+01
 1.00000000e+01 1.00000000e+01 1.20000000e+01 1.20000000e+01]

Eigenvalue clusters (indices):
  Cluster 0: size=1, λ≈0.000000, indices=[0]
  Cluster 1: size=4, λ≈4.000000, indices=[1, 2, 3, 4]
  Cluster 2: size=9, λ≈8.000000, indices=[5, 6, 7, 8, 9, 10, 11, 12, 13]
  Cluster 3: size=8, λ≈10.000000, indices=[14, 15, 16, 17, 18, 19, 20, 21]
  Cluster 4: size=2, λ≈12.000000, indices=[22, 23]

Triplet-like clusters (size >= 3, excluding λ≈0):
  Triplet cluster 0: size=4, λ≈4.000000, indices=[1, 2, 3, 4]
  Triplet cluster 1: size=9, λ≈8.000000, indices=[5, 6, 7, 8, 9, 10, 11, 12, 13]
  Triplet cluster 2: size=8, λ≈10.000000, indices=[14, 15, 16, 17, 18, 19, 20, 21]

=== Yukawa singular values (up to overall scale) ===
Up-type (su):         [2.48234330e-17 2.04201431e-17 5.05511849e-18]
Down-type (sd):       [0.00033546 0.00033546 0.00033546]
Charged leptons (se): [1.44578138e-17 7.50327625e-18 6.53644303e-18]
Neutrino Dirac (sn):  [1.28640670e-17 8.84288118e-18 2.75481512e-18]

=== CKM-like mixing matrix (quarks) ===
[[ 0.6761714 +0.j  0.61757366+0.j -0.40173997+0.j]
 [ 0.61905558+0.j -0.18060801+0.j  0.76429768+0.j]
 [-0.39945266+0.j  0.7654956 +0.j  0.50443439+0.j]]
Approx CKM mixing angles (radians):
theta12_q ≈ 0.740, theta23_q ≈ 0.987, theta13_q ≈ 0.413

=== PMNS-like mixing matrix (leptons) ===
[[ 0.46470253+0.j -0.24086518+0.j -0.85207718+0.j]
 [-0.7985389 +0.j -0.52980016+0.j -0.28574012+0.j]
 [-0.38260578+0.j  0.81320093+0.j -0.43853969+0.j]]
Approx PMNS mixing angles (radians):
theta12_l ≈ 0.478, theta23_l ≈ 0.577, theta13_l ≈ 1.020
"""