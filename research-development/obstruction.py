import numpy as np
import math

# ========= Core model implementation =========

def build_alignment_kernel(kappa1=1.0, kappa4=0.3, phi=0.2, N=9):
    """
    Hermitian cyclic kernel: nearest neighbors ±1 with phase, and ±4 without phase.
    """
    K = np.zeros((N, N), dtype=complex)
    for i in range(N):
        # ±1 with phase
        K[i, (i+1) % N] += kappa1 * np.exp(1j * phi)
        K[i, (i-1) % N] += kappa1 * np.exp(-1j * phi)
        # ±4 real
        K[i, (i+4) % N] += kappa4
        K[i, (i-4) % N] += kappa4
    assert np.allclose(K, K.conj().T)
    return K


def build_orientation_operator(N=9, phases=None):
    """
    O = Phi @ X where X is cyclic shift, Phi is diagonal phase.
    """
    X = np.zeros((N, N), dtype=complex)
    for i in range(N):
        X[i, (i+1) % N] = 1.0

    if phases is None:
        phases = np.linspace(0.0, 0.6, N)

    Phi = np.diag(np.exp(1j * phases))
    O = Phi @ X
    return O


def build_hermitian_orientation(O, mu):
    """
    Hermitian orientation combination.
    NOTE: as run, this uses mu twice (mu^2 factor on Im part):
        O_H = 0.5*(O + O†) + mu * (O - O†)/(2j) * mu
    """
    return (0.5 * (O + O.conj().T)
            + mu * (O - O.conj().T) / (2j)
            + 0.06 * np.diag([0, 1, -1, 0, 1, -1, 0, 1, -1])
            + 0.015 * np.diag([0, 1, 1, 0, 1, 1, 0, 1, 1]))


def sector_selector(triad, epsilon, N=9):
    """
    Sector mask: triad entries weight 1, others epsilon.
    """
    w = np.ones(N) * epsilon
    for i in triad:
        w[i] = 1.0
    return np.diag(w)


def generation_projector_tilted(K, O_full, alpha, mu, mode='lowest'):
    """
    Generation projector as lowest (or highest) eigenmodes of K + alpha * O_H(mu).
    """
    O_H = build_hermitian_orientation(O_full, mu)
    Keff = K + alpha * O_H
    evals, evecs = np.linalg.eigh(Keff)
    if mode == 'lowest':
        idx = np.argsort(evals)[:3]
    else:
        idx = np.argsort(evals)[-3:]
    return evecs[:, idx]


def expm_from_hermitian(K, beta):
    """
    e^{-beta K} for Hermitian K via spectral decomposition.
    """
    evals, evecs = np.linalg.eigh(K)
    return evecs @ np.diag(np.exp(-beta * evals)) @ evecs.conj().T


def generate_yukawa(K, O_full, Wgen, beta, triad, epsilon, c_f):
    """
    Y_f = Wgen† (R_f e^{-beta K} R_f + c_f R_f O R_f) Wgen, normalized.
    """
    N = K.shape[0]
    R = sector_selector(triad, epsilon, N)
    Kflow = expm_from_hermitian(K, beta)
    M_hier = R @ Kflow @ R
    M_orient = c_f * R @ (O_full + 0.12 * np.diag([0, 1, -1, 0, 1, -1, 0, 1, -1])) @ R

    M = M_hier + M_orient
    Y = Wgen.conj().T @ M @ Wgen
    return Y / np.linalg.norm(Y)


def diagonalize_yukawa(Y):
    """
    SVD of Yukawa: Y = U diag(s) V†. We use U (left) and s (singular values).
    """
    U, s, Vh = np.linalg.svd(Y)
    return U, s


def mixing_angles(U):
    """
    PDG-like extraction:
    s13 = |U_{e3}|
    s12 = |U_{e2}|/c13
    s23 = |U_{μ3}|/c13
    """
    s13 = abs(U[0, 2])
    c13 = np.sqrt(max(1.0 - s13**2, 0.0))
    s12 = abs(U[0, 1]) / c13 if c13 > 1e-12 else 0.0
    s23 = abs(U[1, 2]) / c13 if c13 > 1e-12 else 0.0
    s12 = np.clip(s12, 0.0, 1.0)
    s23 = np.clip(s23, 0.0, 1.0)
    th12 = np.degrees(np.arcsin(s12))
    th23 = np.degrees(np.arcsin(s23))
    th13 = np.degrees(np.arcsin(s13))
    return np.array([th12, th23, th13])


def jarlskog(U):
    """
    Simple Jarlskog invariant from first 2×2 block.
    """
    return np.imag(U[0, 0] * U[1, 1] * np.conj(U[0, 1]) * np.conj(U[1, 0]))

def relative_spectral_entropy(Ue, Unu):
    A = Ue.conj().T @ Unu
    p = np.abs(A)**2
    p = p / np.sum(p)
    return -np.sum(p * np.log(p + 1e-12))

# ========= PDG-like neutrino labeling and seesaw =========

def pdg_label_neutrinos(masses):
    """
    Given 3 real masses (non-negative), return permutation indices (i1,i2,i3)
    corresponding to (1,2,3) = (solar1, solar2, remaining) for normal ordering.
    """
    m = np.abs(masses)
    best_pair = None
    best_val = None
    for i in range(3):
        for j in range(i + 1, 3):
            dm2 = m[j]**2 - m[i]**2
            v = abs(dm2)
            if v < 1e-12:
                continue
            if best_val is None or v < best_val:
                best_val = v
                best_pair = (i, j)
    if best_pair is None:
        return np.arange(3)
    i1, i2 = best_pair
    if m[i2] < m[i1]:
        i1, i2 = i2, i1
    remaining = [k for k in range(3) if k not in (i1, i2)][0]
    return np.array([i1, i2, remaining])


def seesaw_neutrino_mass_interpolating(Ynu, Ogen, eta, LambdaR=1.0, rcond=1e-12):
    """
    m_nu = Ynu^T M_R^{-1} Ynu, with
    M_R(eta) = LambdaR * [ (1-eta) I + eta (Ogen† Ogen) ].
    """
    dim = Ynu.shape[0]
    I = np.eye(dim, dtype=complex)
    H = Ogen.conj().T @ Ogen + 0.08 * np.diag([0, 1, -1])

    MR = LambdaR * ((1 - eta) * I + eta * H + 0.04 * np.diag([0, 1, 1]))

    U, s, Vh = np.linalg.svd(MR)
    s_inv = np.array([1 / x if x > rcond else 0.0 for x in s])
    MR_inv = (Vh.conj().T * s_inv) @ U.conj().T
    mnu = Ynu.T @ MR_inv @ Ynu
    mnu = 0.5 * (mnu + mnu.T)
    return mnu


def run_pipeline(alpha_e, mu_e, alpha_nu, mu_nu, eta):
    """
    Compute full pipeline for given lepton embedding + seesaw params.

    Returns dict containing:
      - pe: |U_ei|^2 row after PDG labeling
      - r:  Δm²_21 / |Δm²_31|
      - th12, th23, th13
      - J:  Jarlskog of PMNS
      - Vckm: CKM matrix (for sanity)
    """
    N = 9
    # core objects
    K = build_alignment_kernel()
    O = build_orientation_operator(N=N)

    # projectors
    W_Q = generation_projector_tilted(K, O, 0.0, 0.0, mode='lowest')
    W_e = generation_projector_tilted(K, O, alpha_e, mu_e, mode='lowest')
    W_nu = generation_projector_tilted(K, O, alpha_nu, mu_nu, mode='lowest')

    # triads
    T0 = [0, 3, 6]
    T1 = [1, 4, 7]
    T2 = [2, 5, 8]

    # sector definitions: (beta, triad, epsilon, c_f)
    sectors = {
        'u': (1.00, T0, 0.15, 0.12 * np.exp(1j * 0.30)),
        'd': (1.15, T1, 0.20, 0.10 * np.exp(1j * 1.10)),
        'e': (1.10, T1, 0.30, 0.08 * np.exp(1j * 0.70)),
        'nu': (0.65, T2, 0.70, 0.25 * np.exp(1j * 2.00))
    }

    Y = {}
    U = {}
    masses = {}

    # quarks share W_Q, leptons use their own
    for f, (beta, triad, eps, c_f) in sectors.items():
        if f in ['u', 'd']:
            W = W_Q
        elif f == 'e':
            W = W_e
        else:
            W = W_nu
        Y[f] = generate_yukawa(K, O, W, beta, triad, eps, c_f)
        U[f], masses[f] = diagonalize_yukawa(Y[f])

    # CKM
    Vckm = U['u'].conj().T @ U['d']

    # PMNS from Yukawas (pre-seesaw)
    Upmns = U['e'].conj().T @ U['nu']

    # neutrino seesaw refinement
    Ogen_nu = W_nu.conj().T @ O @ W_nu
    mnu = seesaw_neutrino_mass_interpolating(Y['nu'], Ogen_nu, eta)
    m_evals, U_nu = np.linalg.eigh(mnu)

    # PDG relabel
    perm = pdg_label_neutrinos(m_evals)
    m_sorted = m_evals[perm]
    U_nu_sorted = U_nu[:, perm]

    Upmns = U['e'].conj().T @ U_nu_sorted

    # electron row
    pe = np.abs(Upmns[0, :])**2

    # mass-squared ratio r
    m = np.abs(m_sorted)
    dm2_21 = m[1]**2 - m[0]**2
    dm2_31 = m[2]**2 - m[0]**2
    r = abs(dm2_21) / max(abs(dm2_31), 1e-15)

    th12, th23, th13 = mixing_angles(Upmns)
    J = jarlskog(Upmns)

    return dict(
        pe=pe,
        r=r,
        th12=th12,
        th23=th23,
        th13=th13,
        J=J,
        Vckm=Vckm
    )


# Quick sanity test of pipeline
test_res = run_pipeline(alpha_e=0.05, mu_e=1.0, alpha_nu=0.4, mu_nu=2.0, eta=0.1)


# ========= Loss functions for the scan =========

t = np.array([0.67, 0.30, 0.02])

def H(p):
    p = np.clip(p, 1e-12, 1.0)
    return -np.sum(p * np.log(p))

H_t = H(t)


def electron_row_loss(pe):
    # pe = (|U_e1|^2, |U_e2|^2, |U_e3|^2)
    Le = np.sum((pe - t)**2)
    Lpeak = (np.max(pe) - 0.67)**2
    LH = (H(pe) - H_t)**2
    return Le + 0.5 * Lpeak + 0.5 * LH


def r_window_loss(r):
    # encourage r in [0.015, 0.06]
    if 0.015 <= r <= 0.06:
        return 0.0
    return (r - 0.03)**2


def total_loss(pe, r):
    return electron_row_loss(pe) + 5.0 * r_window_loss(r)


# ========= Random search over parameter space =========

rng = np.random.default_rng(0)
results = []

n_samples = 400

for _ in range(n_samples):
    alpha_e = rng.uniform(0.0, 0.4)
    mu_e = rng.choice([0.5, 1.0, 2.0])
    alpha_nu = rng.uniform(0.1, 0.7)
    mu_nu = rng.choice([1.0, 2.0, 3.0])

    # log-uniform eta, with some exact-zero points
    if rng.random() < 0.1:
        eta = 0.0
    else:
        eta = 10**rng.uniform(-3, math.log10(0.8))

    try:
        res = run_pipeline(alpha_e, mu_e, alpha_nu, mu_nu, eta)
    except Exception:
        continue

    pe = res['pe']
    r = res['r']
    L = total_loss(pe, r)

    results.append((
        L, alpha_e, mu_e, alpha_nu, mu_nu, eta,
        pe, r, res['th12'], res['th23'], res['th13'], res['J']
    ))

# sort by loss and print top 10
results.sort(key=lambda x: x[0])

print(f"Total valid points: {len(results)}\n")

top = results[:10]
for i, (L, a_e, mu_e, a_nu, mu_nu, eta, pe, r, th12, th23, th13, J) in enumerate(top):
    print(f"#{i}: L={L:.4g}, a_e={a_e:.3f}, mu_e={mu_e}, "
          f"a_nu={a_nu:.3f}, mu_nu={mu_nu}, eta={eta:.4f}")
    print(f"    pe={pe}, r={r:.4f}")
    print(f"    th12={th12:.2f}, th23={th23:.2f}, th13={th13:.2f}, J={J:.3e}")
    print()

"""
RESULTS:

Total valid points: 400

#0: L=0.01523, a_e=0.363, mu_e=1.0, a_nu=0.689, mu_nu=1.0, eta=0.0047
    pe=[0.6012254  0.3893162  0.00945839], r=0.0290
    th12=38.82, th23=15.97, th13=5.58, J=8.840e-03

#1: L=0.03013, a_e=0.237, mu_e=1.0, a_nu=0.282, mu_nu=2.0, eta=0.0599
    pe=[0.57114671 0.42347312 0.00538017], r=0.0318
    th12=40.73, th23=80.90, th13=4.21, J=6.830e-04

#2: L=0.04437, a_e=0.129, mu_e=2.0, a_nu=0.695, mu_nu=1.0, eta=0.0069
    pe=[0.54764493 0.44713275 0.00522232], r=0.0298
    th12=42.10, th23=68.85, th13=4.14, J=2.047e-03

#3: L=0.04872, a_e=0.299, mu_e=2.0, a_nu=0.296, mu_nu=2.0, eta=0.2841
    pe=[0.53850723 0.44732131 0.01417146], r=0.0565
    th12=42.35, th23=4.17, th13=6.84, J=3.442e-03

#4: L=0.06981, a_e=0.098, mu_e=0.5, a_nu=0.386, mu_nu=3.0, eta=0.0018
    pe=[6.13211049e-01 3.86394085e-01 3.94866419e-04], r=0.1362
    th12=38.44, th23=80.24, th13=1.14, J=1.219e-03

#5: L=0.07087, a_e=0.127, mu_e=0.5, a_nu=0.698, mu_nu=1.0, eta=0.0933
    pe=[0.51348098 0.48391067 0.00260835], r=0.0300
    th12=44.15, th23=86.11, th13=2.93, J=-6.013e-04

#6: L=0.08021, a_e=0.061, mu_e=0.5, a_nu=0.690, mu_nu=1.0, eta=0.0000
    pe=[0.503342   0.49539247 0.00126554], r=0.0279
    th12=44.77, th23=83.47, th13=2.04, J=1.557e-03

#7: L=0.08754, a_e=0.390, mu_e=2.0, a_nu=0.389, mu_nu=3.0, eta=0.0061
    pe=[0.73526616 0.26056667 0.00416717], r=0.1512
    th12=30.77, th23=3.46, th13=3.70, J=-1.590e-03

#8: L=0.09744, a_e=0.250, mu_e=2.0, a_nu=0.333, mu_nu=2.0, eta=0.5250
    pe=[5.96896545e-01 4.02573246e-01 5.30209302e-04], r=0.1550
    th12=39.39, th23=1.57, th13=1.32, J=2.775e-04

#9: L=0.1097, a_e=0.302, mu_e=1.0, a_nu=0.690, mu_nu=2.0, eta=0.2694
    pe=[0.63061022 0.35993445 0.00945532], r=0.1740
    th12=37.07, th23=76.96, th13=5.58, J=3.926e-03

"""