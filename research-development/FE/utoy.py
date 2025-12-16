###############################################################
# Unified Toy Model v2 (UTM_v2)
# Emergent Geometry → Proto-Flavor Projection → Alignment → SM
# Using Harmonic Neighborhood Projection (Resonant Operator)
###############################################################

import numpy as np
import math

###############################################################
# 0. Utility
###############################################################

def normalize_vec(v):
    return v / np.linalg.norm(v)

###############################################################
# 1. Emergent Geometry (you already have these)
###############################################################


###############################################################
# 2. Harmonic-Neighborhood Proto-Basis Construction (9 modes)
###############################################################

def build_proto_basis_harmonic(evecs, triad, H, num_extra=6):
    """
    Build a 360×9 proto-flavor basis using:
        - 3 triad eigenvectors
        - next-6 modes with highest harmonic strength
    This matches Resonant Model operator theory.
    """

    # Extract triad eigenvectors
    triad_vecs = [normalize_vec(evecs[:, i]) for i in triad]

    # Sort all modes by harmonic strength (descending)
    idx_sorted = np.argsort(-H)   # high → low

    # Remove triad indices
    idx_sorted = [i for i in idx_sorted if i not in triad]

    # Pick next 6 "harmonic neighbors"
    extra_modes = idx_sorted[:num_extra]

    extra_vecs = [normalize_vec(evecs[:, i]) for i in extra_modes]

    # Stack into matrix: shape (360, 9)
    B = np.stack(triad_vecs + extra_vecs, axis=1)

    # Orthonormalize basis
    Q, _ = np.linalg.qr(B)

    # Return first 9 orthonormal columns
    return Q[:, :9]


###############################################################
# 3. Project 360×360 → 9×9 Proto-Yukawa
###############################################################

def project_to_9x9(Y360, B):
    """
    Y360: full geometric Yukawa (360×360)
    B:   proto-basis (360×9)
    """
    return B.conj().T @ Y360 @ B


###############################################################
# 4. Build Full 360×360 Yukawa Seed From Q and R
###############################################################

def build_Y360(phi, q_vals, R_vec, kappa):
    """
    Build the initial 360×360 Yukawa matrix BEFORE Φ-alignment and projection.
    Structure:
        Y_ij = kappa^(q[g(i)] + q[g(j)]) * exp(i(φ_i - φ_j)) * R_g(i) R_g(j)^*
    """

    N = len(phi)
    Y = np.zeros((N, N), dtype=complex)

    # generation index mapping
    def g(i):
        return i % 3

    mag_gen = kappa ** q_vals
    mag = np.array([mag_gen[g(i)] for i in range(N)])

    for i in range(N):
        for j in range(N):
            geom_phase = np.exp(1j * (phi[i] - phi[j]))
            rtwist     = R_vec[g(i)] * np.conj(R_vec[g(j)])
            Y[i, j]    = mag[i] * mag[j] * geom_phase * rtwist

    return Y


###############################################################
# 5. Alignment Operator Φ applied to 9×9
###############################################################

def build_alignment_kernel_9(kappa, forbidden_d):
    N = 9
    K = np.zeros((N, N), dtype=complex)

    for i in range(N):
        for j in range(N):
            d = min(abs(i-j), N - abs(i-j))
            if d == forbidden_d:
                K[i, j] = 0
            else:
                K[i, j] = kappa**d
    return K


def apply_alignment(Y9, kappa, forbidden_d):
    K = build_alignment_kernel_9(kappa, forbidden_d)
    return K * Y9


###############################################################
# 6. Dynamic LIGHT / HEAVY Selection + Schur Reduction
###############################################################

def dynamic_light_heavy_sites(Y9):
    mags = np.sum(np.abs(Y9), axis=1)
    order = np.argsort(mags)
    return order[:3].tolist(), order[3:].tolist()


def schur_reduce_dynamic(Y9):
    LIGHT, HEAVY = dynamic_light_heavy_sites(Y9)

    A = Y9[np.ix_(LIGHT, LIGHT)]
    B = Y9[np.ix_(LIGHT, HEAVY)]
    C = Y9[np.ix_(HEAVY, LIGHT)]
    D = Y9[np.ix_(HEAVY, HEAVY)]

    Dinv = np.linalg.pinv(D)
    Y3 = A - B @ Dinv @ C
    return Y3


###############################################################
# 7. Build Sector Yukawas Using Q Permutations
###############################################################

# Sector orderings (tunable)
PERM_U  = (0, 1, 2)
PERM_D  = (1, 2, 0)
PERM_E  = (2, 1, 0)
PERM_NU = (0, 2, 1)


def exponents_from_perm(q_vals, perm):
    return np.array([q_vals[i] for i in perm])


###############################################################
# 8. CKM / PMNS (use your existing functions)
###############################################################



def build_regions(phi):
    phi = np.asarray(phi)
    N = len(phi)
    idx = np.argsort(phi % (2*np.pi))
    thirds = N // 3
    regions = []
    for k in range(3):
        s = k * thirds
        e = (k+1)*thirds if k < 2 else N
        mask = np.zeros(N, dtype=bool)
        mask[idx[s:e]] = True
        regions.append(mask)
    return regions


###############################################################
# 9. Main Pipeline — UTM v2
###############################################################

def test_run(
    N=360,
    kappa=0.24,
    forbidden_d=5,
    seed=42,
    debug=False
):
    print("\n=== UTM v2: Unified Toy Model (Harmonic-Neighborhood) ===")
    print(f"N = {N}")
    print(f"kappa = {kappa}")
    print(f"forbidden_d = {forbidden_d}")
    print(f"seed = {seed}")

    # ------------------------------------------------------------
    # Step 1: Emergent geometry
    # ------------------------------------------------------------
    print("→ Relaxing phases...")
    phi, E = relax_phases(N=N, steps=1000, eta=0.01, seed=seed)

    print("→ Building adjacency + Laplacian...")
    A = build_adjacency_matrix(phi)
    L = laplacian(A)
    evals, evecs = np.linalg.eigh(L)

    # ------------------------------------------------------------
    # Step 2: Triad + harmonic strengths
    # ------------------------------------------------------------
    print("→ Extracting harmonic triad...")
    H = harmonic_strengths(evals, evecs, phi)
    triad = select_triad(evals, H)

    # ------------------------------------------------------------
    # Step 3: Build R and Q
    # ------------------------------------------------------------
    print("→ Building R and Q operators...")
    R, _ = build_R_from_triad(evecs, phi, triad)
    q_matrix, q_vals = build_Q_from_harmonics(H, triad)
    R_vec = np.diag(R)

    # Sector Q-permutations
    Q_u  = exponents_from_perm(q_vals, PERM_U)
    Q_d  = exponents_from_perm(q_vals, PERM_D)
    Q_e  = exponents_from_perm(q_vals, PERM_E)
    Q_nu = exponents_from_perm(q_vals, PERM_NU)

    # ------------------------------------------------------------
    # Step 4: Build 360×360 Yukawa seed
    # ------------------------------------------------------------
    print("→ Building 360×360 Yukawa seeds...")
    Yu360 = build_Y360(phi, Q_u,  R_vec, kappa)
    Yd360 = build_Y360(phi, Q_d,  R_vec, kappa)
    Ye360 = build_Y360(phi, Q_e,  R_vec, kappa)
    Yn360 = build_Y360(phi, Q_nu, R_vec, kappa)

    # ------------------------------------------------------------
    # Step 5: Build proto basis (harmonic neighborhood)
    # ------------------------------------------------------------
    print("→ Building harmonic-neighborhood proto basis...")
    B = build_proto_basis_harmonic(evecs, triad, H)

    # ------------------------------------------------------------
    # Step 6: Project 360→9
    # ------------------------------------------------------------
    Yu9 = project_to_9x9(Yu360, B)
    Yd9 = project_to_9x9(Yd360, B)
    Ye9 = project_to_9x9(Ye360, B)
    Yn9 = project_to_9x9(Yn360, B)

    # ------------------------------------------------------------
    # Step 7: Apply Φ to each 9×9
    # ------------------------------------------------------------
    print("→ Applying alignment operator Φ...")
    Yu9a = apply_alignment(Yu9, kappa, forbidden_d)
    Yd9a = apply_alignment(Yd9, kappa, forbidden_d)
    Ye9a = apply_alignment(Ye9, kappa, forbidden_d)
    Yn9a = apply_alignment(Yn9, kappa, forbidden_d)

    # ------------------------------------------------------------
    # Step 8: Schur reduction to 3×3 Yukawas
    # ------------------------------------------------------------
    print("→ Schur reducing to 3×3 Yukawas...")
    Yu3 = schur_reduce_dynamic(Yu9a)
    Yd3 = schur_reduce_dynamic(Yd9a)
    Ye3 = schur_reduce_dynamic(Ye9a)
    Yn3 = schur_reduce_dynamic(Yn9a)

    # ------------------------------------------------------------
    # Step 9: CKM & PMNS
    # ------------------------------------------------------------
    print("→ Computing CKM & PMNS...")
    regions = build_regions(phi)

    # geometric unitaries (temporary identity)
    U_geom_u = np.eye(3)
    U_geom_d = np.eye(3)

    assign_e, assign_nu, mixchi2, mix_data = search_best_lepton_regions(
        triad,
        regions,
        U_geom_u,
        U_geom_d,
        Yu3, Yd3, Ye3, Yn3,
        None, None, None
    )

    V_ckm  = mix_data["V_ckm"]
    U_pmns = mix_data["U_pmns"]

    # ------------------------------------------------------------
    # Step 10: Full chi²
    # ------------------------------------------------------------
    print("→ Computing χ²...")
    obs = compute_observables(
        *mix_data["mass_ratios_u"],
        *mix_data["mass_ratios_d"],
        *mix_data["mass_ratios_e"],
        *mix_data["angles_q"],
        *mix_data["angles_l"],
    )

    chi2_total, details = chi2(obs, TARGETS)

    print("\n=== RESULTS ===")
    print("CKM =\n", V_ckm)
    print("PMNS =\n", U_pmns)
    print("χ² =", chi2_total)
    print("Mixing contribution =", mixchi2)

    return dict(
        Y3=dict(u=Yu3, d=Yd3, e=Ye3, nu=Yn3),
        Vckm=V_ckm,
        Upmns=U_pmns,
        chi2=chi2_total,
        details=details
    )


if __name__ == "__main__":
    test_run()
