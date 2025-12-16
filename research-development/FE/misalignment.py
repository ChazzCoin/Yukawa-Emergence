# =================================================================
# misalignment.py
# Full Emergence Architecture — Module 2
# =================================================================
# Implements:
#   - Full divisor-360 misalignment functional M[Ψ]
#   - Local strain, bend, phase, defect energies (Axiom C1)
#   - Gradient-flow evolution operator 𝕄(t)
#   - Emergent dissipative semigroup (Axiom C2)
# =================================================================

import numpy as np
from numpy.linalg import eigh
from FE.harmonics import D360

# ---------------------------------------------------------------
# Utility: normalized phases for a state Ψ
# ---------------------------------------------------------------
def extract_phases(Psi):
    """
    Each internal site has a complex amplitude Psi[i].
    Extract normalized phase φ_i = arg(Psi[i]) mod 2π.
    """
    return np.angle(Psi + 1e-12)


# ---------------------------------------------------------------
# Local misalignment components (Axiom C1)
# ---------------------------------------------------------------
def strain_energy(phi, A):
    """
    M_strain = Σ_edges (φ_i - φ_j)^2 weighted by adjacency.
    Local elasticity of the quasi-crystal.
    """
    N = len(phi)
    E = 0.0
    for i in range(N):
        for j in range(N):
            if A[i, j] > 0:
                d = phi[i] - phi[j]
                E += d * d
    return 0.5 * E


def bend_energy(phi, A):
    """
    M_bend = Σ i (Δφ_i)^2
    where Δ is graph Laplacian acting on phases.
    """
    N = len(phi)
    d = np.sum(A, axis=1)
    Lphi = d * phi - A.dot(phi)
    return np.sum(Lphi * Lphi)


def phase_energy(phi):
    """
    M_phase = Σ_{n,m∈D360} Σ_{i,j} (1 - cos(n φ_i - m φ_j))

    The full harmonic misalignment functional demanded by Mode A.
    """
    N = len(phi)
    E = 0.0

    for i in range(N):
        for j in range(N):
            for n in D360:
                for m in D360:
                    E += 1.0 - np.cos(n * phi[i] - m * phi[j])

    return E


def defect_energy(phi, A):
    """
    M_defect = Σ_i |coord(i) - median_coord| * local phase incoherence

    Integer-coordination mismatches produce defect energy,
    matching the Axiom C1 requirement for defect contributions.
    """
    N = len(phi)
    degrees = np.sum(A, axis=1)
    med = np.median(degrees)

    E = 0.0
    for i in range(N):
        mismatch = abs(degrees[i] - med)
        incoh = 1.0 - np.abs(np.mean(np.exp(1j * (phi[i] - phi))))
        E += mismatch * incoh

    return E


# ---------------------------------------------------------------
# Full misalignment functional (Axiom C1)
# ---------------------------------------------------------------
def misalignment(Psi, A):
    """
    M[Ψ] = M_strain + M_bend + M_phase + M_defect

    No coefficients. No α, β, λ parameters.
    Full emergence requires all contributions be equally fundamental.
    """
    phi = extract_phases(Psi)

    M_s = strain_energy(phi, A)
    M_b = bend_energy(phi, A)
    M_p = phase_energy(phi)
    M_d = defect_energy(phi, A)

    return M_s + M_b + M_p + M_d


# ---------------------------------------------------------------
# Gradient of misalignment wrt phases
# ---------------------------------------------------------------
def misalignment_grad(Psi, A):
    """
    Compute δM / δφ_i using numerical differentiation.
    Because the misalignment is highly nonlinear (full D360 tensor),
    explicit analytic gradients would be enormous.

    Numerical gradient is acceptable under FULL EMERGENCE because
    it is not a model parameter — just a computational representation
    of the exact functional derivative.
    """
    eps = 1e-6
    phi = extract_phases(Psi)
    grad = np.zeros_like(phi, dtype=float)

    for i in range(len(phi)):
        phi_perturb = phi.copy()

        # Forward variation
        phi_perturb[i] = phi[i] + eps
        Psi_f = np.exp(1j * phi_perturb)
        M_f = misalignment(Psi_f, A)

        # Backward variation
        phi_perturb[i] = phi[i] - eps
        Psi_b = np.exp(1j * phi_perturb)
        M_b = misalignment(Psi_b, A)

        grad[i] = (M_f - M_b) / (2 * eps)

    return grad


# ---------------------------------------------------------------
# Evolution operator 𝕄(t) (Axiom C2)
# ---------------------------------------------------------------
def evolve(Psi0, A, n_steps=200, eta=0.01):
    """
    Evolution operator:
        Ψ(t+dt) = Ψ(t) - η * δM/δΨ*
    but since Ψ lives in U(1)^N, we evolve only phases:

        φ_i ← φ_i - η * δM/δφ_i

    The semigroup property 𝕄(t+s) = 𝕄(t)𝕄(s) is satisfied
    by explicit Euler integration of a gradient flow.

    Contractivity emerges automatically from minimization of M.

    Full emergence: no tuning. eta is not a physics parameter,
    just discretization step size of gradient descent.
    """

    Psi = Psi0.copy()

    for _ in range(n_steps):
        grad_phi = misalignment_grad(Psi, A)
        phi = extract_phases(Psi)
        phi -= eta * grad_phi
        Psi = np.exp(1j * phi)

    return Psi


# ---------------------------------------------------------------
# Master function: compute fully relaxed internal state
# ---------------------------------------------------------------
def relax_internal_state(A, N_sites):
    """
    FULL EMERGENCE MASTER FUNCTION:
    Relax an internal state Ψ until misalignment is minimized.

    Returns:
        Ψ_relaxed  — fixed-point of internal evolution
        M_final    — final misalignment value
    """

    # Initial random internal phases
    rng = np.random.default_rng(1234)
    Psi0 = np.exp(1j * rng.uniform(0, 2 * np.pi, size=N_sites))

    Psi_rel = evolve(Psi0, A, n_steps=150, eta=0.01)
    M_final = misalignment(Psi_rel, A)

    return Psi_rel, M_final

# ------------------------------------------------
# 1. Relax internal phases (FULL EMERGENCE)
# ------------------------------------------------
def relax_internal_phases(
    N_sites: int = 200,
    n_steps: int = 300,
    eta: float = 0.01,
    seed: int = 42,
) -> np.ndarray:
    """
    Fully emergent phase relaxation using divisor-360 harmonics.

    We evolve phases φ_i via gradient descent on

        M[φ] = Σ_{i,j} Σ_{n∈D360} (1 - cos(n (φ_i - φ_j))).

    No phenomenological parameters; eta is just a numerical step.
    """
    rng = np.random.default_rng(seed)
    phi = rng.uniform(0.0, 2.0 * np.pi, size=N_sites)

    for _ in range(n_steps):
        # pairwise differences
        dphi = phi[:, None] - phi[None, :]

        # misalignment gradient (simplified but fully divisor-based)
        grad = np.zeros_like(phi)
        for n in D360:
            grad += n * np.sum(np.sin(n * dphi), axis=1)

        phi = (phi - eta * grad) % (2.0 * np.pi)

    return phi


# ------------------------------------------------
# 2. Build internal quasi-crystal graph from phases
# ------------------------------------------------
def build_internal_graph(phi: np.ndarray) -> np.ndarray:
    """
    Build emergent adjacency matrix A from relaxed phases φ using
    divisor-harmonic-induced sparsity (FULL EMERGENCE).

    A_ij = 1 if C_ij >= mean(C) + std(C)
    """

    N = len(phi)
    C = np.zeros((N, N), dtype=float)

    # Compute harmonic correlation using D360
    for i in range(N):
        for j in range(N):
            d = phi[i] - phi[j]
            C_ij = 0.0
            for n in D360:
                C_ij += np.cos(n * d)
            C[i, j] = C_ij

    # Compute emergent threshold
    mask_off = ~np.eye(N, dtype=bool)
    flat = C[mask_off]
    mean_C = np.mean(flat)
    std_C = np.std(flat)

    threshold = mean_C + std_C

    # Build adjacency
    A = (C >= threshold).astype(float)
    np.fill_diagonal(A, 0.0)
    A = np.maximum(A, A.T)

    return A

# ------------------------------------------------
# 3. Laplacian and spectrum
# ------------------------------------------------
def laplacian_from_adjacency(A: np.ndarray) -> np.ndarray:
    d = np.sum(A, axis=1)
    return np.diag(d) - A


def laplacian_spectrum(A: np.ndarray):
    """
    Return eigenvalues and eigenvectors of the internal Laplacian.
    """
    L = laplacian_from_adjacency(A)
    evals, evecs = eigh(L)
    return evals, evecs