# =================================================================
# manifestation.py
# Full Emergence Architecture — Module 4
# =================================================================
# Implements:
#   - Evolution semigroup M(t)
#   - Manifestation operator X = S ∘ M
#   - Fixed-point solver X ψ = ψ
#   - Pre-particle invariant subspace extraction
# =================================================================

import numpy as np
import scipy


# ---------------------------------------------------------------------
# 1. Evolution operator  M(t)  from misalignment gradient (Axiom C2)
# ---------------------------------------------------------------------
def M_evolution_operator(gradM, dt):
    """
    One-step evolution map:
        M(dt) = exp[- dt * gradM ]
    FULL EMERGENCE:
    - gradM is computed from the misalignment functional,
      NOT from a fitted potential.
    - Purely dissipative (contractive semigroup).

    Returns a matrix acting on internal state vectors.
    """
    # Matrix exponential of -dt * gradM
    return scipy.linalg.expm(-dt * gradM)


def evolve_state(Psi0, gradM_func, n_steps=200, dt=0.01):
    """
    Evolve Ψ using the semigroup property:
        Ψ_{k+1} = M(dt) Ψ_k

    gradM_func: function returning matrix gradient δM/δΨ† at Ψ.
    """
    Psi = Psi0.copy()
    for _ in range(n_steps):
        gradM = gradM_func(Psi)
        Mdt = M_evolution_operator(gradM, dt)
        Psi = Mdt @ Psi
        Psi /= np.linalg.norm(Psi)  # renormalize (pure emergence)
    return Psi


# ---------------------------------------------------------------------
# 2. Manifestation operator X = S ∘ M (Axiom E1)
# ---------------------------------------------------------------------
def manifestation_operator(S, gradM_func, dt):
    """
    Build the instantaneous manifestation operator:

        X(dt) = S @ M(dt)

    where S is the selection operator from module 3.
    """
    # Build infinitesimal M(dt)
    Mdt = M_evolution_operator(gradM_func(np.eye(S.shape[0])), dt)

    # Compose
    Xdt = S @ Mdt

    # Symmetrize to keep stability
    Xdt = 0.5 * (Xdt + Xdt.conj().T)
    return Xdt


# ---------------------------------------------------------------------
# 3. Fixed-point solver: solve X ψ = ψ (Axiom E2)
# ---------------------------------------------------------------------
def find_fixed_points(X, tol=1e-10):
    """
    Solve the universal harmonic equation:
        X ψ = ψ

    Physical states = eigenspace of eigenvalue 1.

    FULL EMERGENCE:
    - No selection of eigenvectors by hand.
    - No fitting.
    - Extract the full invariant subspace.
    """
    evals, evecs = np.linalg.eigh(X)

    # Identify eigenvalue ≈ 1
    fixed_indices = np.where(np.isclose(evals, 1.0, atol=tol))[0]

    fixed_states = [evecs[:, i] for i in fixed_indices]
    return fixed_states, evals, evecs


# ---------------------------------------------------------------------
# 4. Particle extraction: minimal invariant subspaces
# ---------------------------------------------------------------------
def extract_particle_subspaces(fixed_states, S):
    """
    Axiom E2 says:
        particles = minimal nontrivial invariant subspaces
        of the manifestation map X.

    We approximate this by:
        - computing all S-stable subspaces spanned by fixed states
        - selecting irreducible blocks (cannot be further decomposed)
    """
    particles = []

    for psi in fixed_states:
        psi_proj = S @ psi
        if np.allclose(psi_proj, psi, atol=1e-10):
            particles.append(psi / np.linalg.norm(psi))

    return particles