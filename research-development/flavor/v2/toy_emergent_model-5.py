import numpy as np

"""
AF-2S Unified Fully Emergent Flavor Pipeline

This file implements a self-contained toy realization of the AF-2S
(Alignment Framework, strong internal–external spectral interlock)
within the constraints of fully emergent, real-physics-inspired axioms:

  - A: Internal quasi-crystal graph + Laplacian L on H_int
  - B: Core internal operators R (base-360), Q (integer charges)
  - C: Misalignment functional M[phi] and its gradient flow
  - D: Selection operator S = C_360 B P_phi on H_int
  - E: Manifested internal state Psi_phys as fixed point of S∘M
  - F: Flavor / Yukawa structure from (R,Q) plus geometric sector embeddings
  - G: Minimality + SM chi^2 as EXTERNAL diagnostic only

All sector structure (u,d,e,nu), hierarchies and mixings are emergent
from (L, phi) and discrete operators derived from them. No sector-tuned
continuous parameters are introduced; SM data only enters through a
χ^2 diagnostic printed at the end.
"""

###############################################################################
# 0. Generic utilities
###############################################################################


def normalize(v):
    n = np.linalg.norm(v)
    return v / n if n > 0 else v


def make_projector(cols):
    """Given a list/array of column vectors, return projector onto their span."""
    if isinstance(cols, list):
        M = np.column_stack(cols)
    else:
        M = np.asarray(cols)
    gram = M.conj().T @ M
    try:
        gram_inv = np.linalg.inv(gram)
    except np.linalg.LinAlgError:
        gram_inv = np.linalg.pinv(gram)
    return M @ gram_inv @ M.conj().T


###############################################################################
# 1. Internal Fibonacci quasi-crystal graph  (A2)
###############################################################################


def fibonacci_word(n):
    """Generate a Fibonacci word (aperiodic sequence) of length >= n."""
    a, b = "A", "AB"  # substitution: A→AB, B→A
    while len(b) < n:
        a, b = b, b + a
    return b[:n]


def build_internal_graph_fibonacci(N, wA=1.0, wB=0.618):
    """Finite 1D Fibonacci quasi-crystal segment as weighted graph."""
    word = fibonacci_word(N)
    A = np.zeros((N, N), dtype=float)
    for i in range(N - 1):
        w = wA if word[i] == "A" else wB
        A[i, i + 1] = w
        A[i + 1, i] = w
    return A


def laplacian(A):
    d = np.sum(A, axis=1)
    return np.diag(d) - A


###############################################################################
# 2. Misalignment functional on phases (C1)
###############################################################################


def misalignment_energy(phi, w6=1.0, w5=1.0):
    """5- and 6-fold phase misalignment energy on the complete graph."""
    N = len(phi)
    diffs = phi[:, None] - phi[None, :]
    E6 = w6 * np.sum(1.0 - np.cos(6 * diffs)) / (N * N)
    E5 = w5 * np.sum(1.0 - np.cos(5 * diffs)) / (N * N)
    return E6 + E5


def relax_phases(N=360, steps=800, eta=0.01, seed=42):
    """Gradient-flow for M[phi] → relaxed phase configuration (C2)."""
    rng = np.random.default_rng(seed)
    phi = rng.uniform(0, 2 * np.pi, size=N)
    for _ in range(steps):
        diffs = phi[:, None] - phi[None, :]
        grad = 6 * np.sum(np.sin(6 * diffs), axis=1) + 5 * np.sum(
            np.sin(5 * diffs), axis=1
        )
        phi -= eta * grad
        phi %= 2 * np.pi
    return phi, misalignment_energy(phi)


###############################################################################
# 3. Harmonics, triad selection, and R,Q on H_gen  (B1,B2)
###############################################################################


def harmonic_strengths(evals, evecs, phi, divisors=(1, 2, 3)):
    """Overlap of each Laplacian mode with e^{-i d φ} for a few d."""
    N = len(phi)
    V = evecs
    strengths = []
    exp_dphi = {d: np.exp(-1j * d * phi) for d in divisors}
    for n in range(N):
        v = V[:, n]
        Hn = []
        for d in divisors:
            amp = np.dot(v, exp_dphi[d])
            Hn.append(np.abs(amp) ** 2)
        strengths.append(Hn)
    return np.array(strengths)


def select_triad(evals, H):
    """Select 3 modes (≠ zero) with largest harmonic strength as generation triad."""
    scores = np.sum(H[1:], axis=1)
    idx = np.arange(1, len(evals))
    top3_rel = np.argsort(scores)[-3:]
    triad = np.sort(idx[top3_rel])
    return triad


def build_R_from_triad(evecs, phi, triad):
    """Base-360 internal symmetry R on 3D triad subspace (B1)."""
    V = evecs
    c = []
    for n in triad:
        v = V[:, n]
        c.append(np.dot(v, np.exp(1j * phi)))
    c = np.array(c)
    theta = np.angle(c)
    k = np.round(theta * 360 / (2 * np.pi)).astype(int) % 360
    R = np.diag(np.exp(1j * 2 * np.pi * k / 360))
    return R, k


def build_Q_from_harmonics(H, triad):
    """Integer charges q_j from triad harmonic strengths (B2, F2)."""
    A = np.sum(H[triad], axis=1)
    order = np.argsort(-A)
    q = np.empty(3, dtype=int)
    q[order] = np.array([1, 2, 3])
    Q = np.diag(q.astype(float))
    return Q, q


###############################################################################
# 4. Selection operator pieces C360, B, P_phi on H_int  (D1–D4)
###############################################################################


def projector_C360_from_triad(evecs, triad):
    """C_360: projector onto triad subspace inside H_int (D1)."""
    cols = [evecs[:, n] for n in triad]
    P = make_projector(cols)
    return 0.5 * (P + P.conj().T)


def local_phase_curvature(A, phi):
    """Discrete phase curvature kappa_i = sum_j A_ij (phi_i - phi_j)."""
    N = len(phi)
    kappa = np.zeros(N, dtype=float)
    for i in range(N):
        nbrs = np.where(A[i] != 0.0)[0]
        if len(nbrs) == 0:
            kappa[i] = 0.0
        else:
            kappa[i] = np.sum(A[i, nbrs] * (phi[i] - phi[nbrs]))
    return kappa


def projector_B_from_curvature(A, phi, frac_light=0.5):
    """B: projector onto 'light' low-curvature subspace (D2)."""
    N = len(phi)
    kappa = local_phase_curvature(A, phi)
    order = np.argsort(np.abs(kappa))
    keep = order[: int(frac_light * N)]
    mask = np.zeros(N, dtype=float)
    mask[keep] = 1.0
    return np.diag(mask)


def phase_regions(phi, n_regions=3):
    """Partition sites into n_regions by phase angle."""
    phase = np.mod(phi, 2 * np.pi)
    edges = np.linspace(0, 2 * np.pi, n_regions + 1)
    labels = np.zeros(len(phi), dtype=int)
    for k in range(n_regions):
        lo, hi = edges[k], edges[k + 1]
        if k < n_regions - 1:
            idx = np.where((phase >= lo) & (phase < hi))[0]
        else:
            idx = np.where((phase >= lo) & (phase <= hi))[0]
        labels[idx] = k
    return labels


def projector_Pphi_from_phase(phi):
    """P_phi: phase-coherence projector from a finite Z_3 symmetry (D3)."""
    N = len(phi)
    labels = phase_regions(phi, n_regions=3)
    U = []
    for g in range(3):
        phase_g = np.exp(2j * np.pi * g * labels / 3.0)
        U_g = np.diag(phase_g)
        U.append(U_g)
    P = sum(U) / 3.0
    return 0.5 * (P + P.conj().T)


def build_selection_operator(A, phi, evecs, triad):
    """Full selection operator S = C_360 B P_phi on H_int (D4)."""
    C = projector_C360_from_triad(evecs, triad)
    B = projector_B_from_curvature(A, phi, frac_light=0.5)
    P_phi = projector_Pphi_from_phase(phi)

    S = C @ B @ P_phi
    P_eff = S @ S.conj().T
    P_eff = 0.5 * (P_eff + P_eff.conj().T)
    return P_eff


###############################################################################
# 5. Manifested internal state Ψ_phys (E1,E2)
###############################################################################


def manifested_state(phi, S):
    """Ψ_phys = S Ψ_0 / ||S Ψ_0||, with Ψ_0 uniform on H_int (E2)."""
    N = len(phi)
    psi0 = np.ones(N, dtype=complex) / np.sqrt(N)
    psi = S @ psi0
    norm = np.linalg.norm(psi)
    if norm < 1e-14:
        return psi0
    return psi / norm


def projected_generation_state(psi_phys, evecs, triad):
    """Coordinates of Ψ_phys in the generation triad basis (3-vector)."""
    V_gen = evecs[:, triad]
    c = V_gen.conj().T @ psi_phys
    return normalize(c)


###############################################################################
# 6. Yukawa kernel and sector embeddings (F1–F3)
###############################################################################


def yukawa_kernel_from_RQ(R, Q):
    """Universal Yukawa kernel F(R,Q) on H_gen (F1)."""
    alpha = 1.0
    beta = np.log(5.0)  # λ ~ 0.2 hierarchy scale (global, no sector tuning)

    ReR = (R + R.conj().T) / 2.0
    F_R = np.exp(-alpha * (np.eye(3) - ReR.real))
    Q_diag = np.diag(Q)
    F_Q = np.diag(np.exp(-beta * Q_diag))
    return F_R @ F_Q


def sector_masks_from_curvature_and_phase(A, phi):
    """Four geometric masks on H_int, one per sector (F2)."""
    N = len(phi)
    kappa = local_phase_curvature(A, phi)
    labels_phi = phase_regions(phi, n_regions=2)
    kappa_n = kappa / (np.max(np.abs(kappa)) + 1e-12)

    masks = {}
    masks["u"] = ((np.abs(kappa_n) < 0.5) & (labels_phi == 0)).astype(float)
    masks["d"] = ((np.abs(kappa_n) < 0.5) & (labels_phi == 1)).astype(float)
    masks["e"] = ((np.abs(kappa_n) >= 0.5) & (labels_phi == 0)).astype(float)
    masks["nu"] = ((np.abs(kappa_n) >= 0.5) & (labels_phi == 1)).astype(float)
    return masks


def sector_U_L_from_masks(evecs, triad, masks):
    """Emergent U_L^{(s)} from geometric Gram metrics G_s (F2,F3)."""
    V_gen = evecs[:, triad]
    sectors_U = {}
    for name, m in masks.items():
        D_s = np.diag(m)
        G_s = V_gen.conj().T @ D_s @ V_gen
        G_s = 0.5 * (G_s + G_s.conj().T)
        w, U = np.linalg.eigh(G_s)
        order = np.argsort(-w.real)
        U_L_s = U[:, order]
        sectors_U[name] = U_L_s
    return sectors_U


def build_sector_Yukawas(F_Y, sectors_U):
    """Sector Yukawas Y_s = U_L^{(s)} F_Y; masses = singular values (F3)."""
    sectors = {}
    for name, U_L in sectors_U.items():
        Y_s = U_L @ F_Y
        sv = np.linalg.svd(Y_s, compute_uv=False)
        sectors[name] = dict(Y=Y_s, masses=sv, U_L=U_L)
    return sectors


###############################################################################
# 7. SM-like diagnostic: mixing angles + χ² (external only)
###############################################################################

# Rough SM targets (at some reference scale); used ONLY for diagnostics.
SM_TARGETS = {
    # mass ratios (m1/m3, m2/m3) per sector
    "u_m2/m3": (2.2e-05, 0.5 * 2.2e-05),
    "u_m1/m3": (7.5e-03, 0.5 * 7.5e-03),
    "d_m2/m3": (1.1e-03, 0.5 * 1.1e-03),
    "d_m1/m3": (2.2e-02, 0.5 * 2.2e-02),
    "e_m2/m3": (2.9e-04, 0.5 * 2.9e-04),
    "e_m1/m3": (5.9e-02, 0.5 * 5.9e-02),
    # CKM
    "CKM_theta12": (0.227, 0.05 * 0.227),
    "CKM_theta23": (0.041, 0.5 * 0.041),
    "CKM_theta13": (0.0036, 0.5 * 0.0036),
    # PMNS
    "PMNS_theta12": (0.584, 0.1 * 0.584),
    "PMNS_theta23": (0.785, 0.2 * 0.785),
    "PMNS_theta13": (0.15, 0.2 * 0.15),
}


def mass_ratio_observables(sectors):
    """Compute (m1/m3, m2/m3) per sector from singular values."""
    obs = {}
    for name in ["u", "d", "e", "nu"]:
        sv = np.sort(np.abs(sectors[name]["masses"]))
        if sv[-1] <= 0:
            r1, r2 = 1.0, 1.0
        else:
            r1, r2 = sv[0] / sv[-1], sv[1] / sv[-1]
        if name != "nu":  # only diagnose charged sectors vs SM
            obs[f"{name}_m1/m3"] = r1
            obs[f"{name}_m2/m3"] = r2
    return obs


def mixing_angles_from_U(U):
    """Standard parameterization angles from a 3×3 unitary U."""
    s13 = abs(U[0, 2])
    s13 = min(max(s13, 0.0), 1.0)
    theta13 = np.arcsin(s13)
    c13 = np.cos(theta13)
    if abs(c13) < 1e-12:
        theta12 = 0.0
        theta23 = 0.0
    else:
        theta12 = np.arctan2(abs(U[0, 1]), abs(U[0, 0]))
        theta23 = np.arctan2(abs(U[1, 2]), abs(U[2, 2]))
    return theta12, theta23, theta13


def mixing_observables(CKM, PMNS):
    """Extract CKM/PMNS angles for χ² diagnostics."""
    t12_q, t23_q, t13_q = mixing_angles_from_U(CKM)
    t12_l, t23_l, t13_l = mixing_angles_from_U(PMNS)
    return {
        "CKM_theta12": t12_q,
        "CKM_theta23": t23_q,
        "CKM_theta13": t13_q,
        "PMNS_theta12": t12_l,
        "PMNS_theta23": t23_l,
        "PMNS_theta13": t13_l,
    }


def chi2_from_observables(obs, targets=SM_TARGETS):
    """External χ² diagnostic: how SM-like is this emergent vacuum?"""
    chi2 = 0.0
    contribs = {}
    for k, v in obs.items():
        if k not in targets:
            continue
        mean, sigma = targets[k]
        if sigma <= 0:
            continue
        c = ((v - mean) / sigma) ** 2
        chi2 += c
        contribs[k] = (v, mean, c)
    return chi2, contribs


###############################################################################
# 8. Full AF-2S internal pipeline + diagnostics
###############################################################################


def run_emergent_pipeline(N=1080, seed=42):
    print("=== MULTI-SECTOR FULLY EMERGENT (CURVATURE-SPLIT) MODEL ===")

    # 1) Misalignment flow on phases (C1,C2)
    phi, E = relax_phases(N=N, seed=seed)
    print(f"Misalignment energy: {E}")

    # 2) Internal quasi-crystal graph + Laplacian (A2)
    A = build_internal_graph_fibonacci(N)
    L = laplacian(A)
    evals, evecs = np.linalg.eigh(L)
    print("First 5 eigenvalues:", evals[:5])

    # 3) Harmonic analysis + triad (B1,B2,B3)
    H = harmonic_strengths(evals, evecs, phi)
    triad = select_triad(evals, H)
    print("Triad selection →")
    print("Triad:", triad)
    print("Triad λ:", evals[triad])

    R, k = build_R_from_triad(evecs, phi, triad)
    Q, q = build_Q_from_harmonics(H, triad)
    print("k_j:", k)
    print("q_j:", q)

    # 4) Selection operator on H_int (D1–D4)
    S = build_selection_operator(A, phi, evecs, triad)

    # 5) Manifested internal state and its projection to H_gen (E1,E2)
    psi_phys = manifested_state(phi, S)
    c_gen = projected_generation_state(psi_phys, evecs, triad)
    print("||Ψ_phys||:", np.linalg.norm(psi_phys))
    print("Projected generation state (components in triad basis):", c_gen)

    # 6) Universal Yukawa kernel (F1)
    F_Y = yukawa_kernel_from_RQ(R, Q)

    # 7) Sector masks + emergent U_L^{(s)} (F2,F3)
    masks = sector_masks_from_curvature_and_phase(A, phi)
    sectors_U = sector_U_L_from_masks(evecs, triad, masks)

    # 8) Build sector Yukawas and mixing matrices (F3)
    sectors = build_sector_Yukawas(F_Y, sectors_U)

    U_L_u = sectors["u"]["U_L"]
    U_L_d = sectors["d"]["U_L"]
    U_L_e = sectors["e"]["U_L"]
    U_L_nu = sectors["nu"]["U_L"]

    CKM = U_L_u.conj().T @ U_L_d
    PMNS = U_L_e.conj().T @ U_L_nu

    print("\n=== Emergent mass spectra (singular values, normalized) ===")
    for name in ["u", "d", "e", "nu"]:
        sv = sectors[name]["masses"]
        sv = sv / sv.max() if sv.max() > 0 else sv
        print(f"{name}: {sv}")

    print("\n=== CKM (emergent) ===")
    print(CKM)
    print("CKM |V|:")
    print(np.abs(CKM))

    print("\n=== PMNS (emergent) ===")
    print(PMNS)
    print("PMNS |U|:")
    print(np.abs(PMNS))

    # 9) External SM-like χ² diagnostic (G2)
    obs = {}
    obs.update(mass_ratio_observables(sectors))
    obs.update(mixing_observables(CKM, PMNS))
    chi2, contribs = chi2_from_observables(obs)

    print("\n=== Global χ² Diagnostic (external, non-emergent) ===")
    print("Total χ²:", chi2)
    print("Per-term contributions:")
    for k, (val, tgt, c) in contribs.items():
        print(f"  {k:12s}: model={val:.3e}, target={tgt:.3e}, χ²={c:.2f}")

    return dict(
        A=A,
        L=L,
        evals=evals,
        evecs=evecs,
        phi=phi,
        triad=triad,
        R=R,
        Q=Q,
        S=S,
        psi_phys=psi_phys,
        sectors=sectors,
        CKM=CKM,
        PMNS=PMNS,
        chi2=chi2,
        chi2_contribs=contribs,
    )


if __name__ == "__main__":
    run_emergent_pipeline()

"""
/Users/chazzromeo/rAI/Flavor/.venv/bin/python /Users/chazzromeo/rAI/Flavor/toy_emergent_model-5.py 
=== MULTI-SECTOR FULLY EMERGENT (CURVATURE-SPLIT) MODEL ===
Misalignment energy: 1.9967149483235276
First 5 eigenvalues: [6.86382668e-17 6.84537951e-06 2.73813254e-05 6.16080479e-05
 1.09523752e-04]
Triad selection →
Triad: [194 609 868]
Triad λ: [0.24889194 1.94185901 3.22773355]
k_j: [153  11 359]
q_j: [1 2 3]
||Ψ_phys||: 0.9999999999999969
Projected generation state (components in triad basis): [ 0.87847349+0.j  0.4610587 +0.j -0.12533635+0.j]

=== Emergent mass spectra (singular values, normalized) ===
u: [1.00000000e+00 1.15841670e-01 3.49296761e-04]
d: [1.00000000e+00 1.15841670e-01 3.49296761e-04]
e: [1.00000000e+00 1.15841670e-01 3.49296761e-04]
nu: [1.00000000e+00 1.15841670e-01 3.49296761e-04]

=== CKM (emergent) ===
[[ 0.1083333   0.19494879 -0.97481222]
 [ 0.32881092 -0.93241871 -0.1499291 ]
 [ 0.93816165  0.30428658  0.16511329]]
CKM |V|:
[[0.1083333  0.19494879 0.97481222]
 [0.32881092 0.93241871 0.1499291 ]
 [0.93816165 0.30428658 0.16511329]]

=== PMNS (emergent) ===
[[ 0.75446421 -0.24947939 -0.60707807]
 [-0.62537881  0.00748843 -0.78028537]
 [ 0.19921118  0.96835115 -0.15036937]]
PMNS |U|:
[[0.75446421 0.24947939 0.60707807]
 [0.62537881 0.00748843 0.78028537]
 [0.19921118 0.96835115 0.15036937]]

=== Global χ² Diagnostic (external, non-emergent) ===
Total χ²: 112102704.74822819
Per-term contributions:
  u_m1/m3     : model=3.493e-04, target=7.500e-03, χ²=3.64
  u_m2/m3     : model=1.158e-01, target=2.200e-05, χ²=110861123.29
  d_m1/m3     : model=3.493e-04, target=2.200e-02, χ²=3.87
  d_m2/m3     : model=1.158e-01, target=1.100e-03, χ²=43522.81
  e_m1/m3     : model=3.493e-04, target=5.900e-02, χ²=3.95
  e_m2/m3     : model=1.158e-01, target=2.900e-04, χ²=635062.47
  CKM_theta12 : model=1.064e+00, target=2.270e-01, χ²=5432.88
  CKM_theta23 : model=7.372e-01, target=4.100e-02, χ²=1153.47
  CKM_theta13 : model=1.346e+00, target=3.600e-03, χ²=556083.02
  PMNS_theta12: model=3.194e-01, target=5.840e-01, χ²=20.54
  PMNS_theta23: model=1.380e+00, target=7.850e-01, χ²=14.38
  PMNS_theta13: model=6.524e-01, target=1.500e-01, χ²=280.43

"""