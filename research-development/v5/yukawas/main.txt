"""
Main execution pipeline for the Alignment Spectral Triple v2.0.

This script:
    1. Builds the geometric triple
    2. Builds the triadic Z_2160 kernel
    3. Builds the finite spectral triple
    4. Applies phases to get D_F
    5. Builds 9→3 triadic compression operator
    6. Builds the full product Dirac operator D
    7. Applies a sample inner fluctuation
    8. Computes the bosonic + fermionic spectral action

Everything here is modular; users can insert custom site choices,
custom weight vectors, and custom fluctuations A.
"""

import numpy as np

# ------------------------------------------------------------
# Import your modules
# ------------------------------------------------------------
from geometry import GeometricTriple
from kernel import TriadicKernel
from finite_triple import FiniteTriple
from compression import TriadicCompression
from spectral_triple import AlignmentSpectralTriple
from spectral_action import SpectralAction, SpectralActionV3
from guage_sector import GaugeSector


# ------------------------------------------------------------
# Utility: Create sample data for testing
# ------------------------------------------------------------
def sample_sites():
    """
    Example choice of 9 distinct sites embedded in Z_2160.
    These should be tuned based on phenomenology.
    """
    return [0, 10, 25, 100, 200, 350, 777, 1024, 1500]


def sample_phases():
    """
    Example phase vector φ_i for finite Dirac operator.
    In real fits one adjusts these phases.
    """
    return np.linspace(0, np.pi, 9)


def sample_weights():
    """
    Example R^(s) sector weights.
    Real model tuning will choose different weights per sector.
    """
    return np.ones(9)


def sample_modes():
    """
    Example geometric harmonic modes.
    These must be divisors of 360 to preserve alignment.
    """
    return [1, 2, 3, 4, 5, 6]


# ------------------------------------------------------------
# Full pipeline
# ------------------------------------------------------------
def main():
    print("\n=== Alignment Spectral Triple v2.0 Pipeline ===\n")

    # ------------------------------------------------------------
    # 1. Geometric Triple
    # ------------------------------------------------------------
    geom = GeometricTriple(N_max=360)
    modes = sample_modes()

    # ------------------------------------------------------------
    # 2. Triadic Kernel (real)
    # ------------------------------------------------------------
    sites = sample_sites()
    triadic_kernel = TriadicKernel(
        kappa=0.24,
        forbidden=(2, 4, 7),
        N=2160,
        sites=sites
    )
    K_real = triadic_kernel.K

    # ------------------------------------------------------------
    # 3. Finite triple (add phases)
    # ------------------------------------------------------------
    phases = sample_phases()
    finite = FiniteTriple(K_real, phases)

    # Optional: Provide finite grading
    gamma_F = np.diag([1, 1, 1, -1, -1, -1, 1, -1, 1])
    finite_gamma = gamma_F

    # ------------------------------------------------------------
    # 4. Triadic compression (9 → 3)
    # ------------------------------------------------------------
    triads = [
        [0, 1, 2],
        [3, 4, 5],
        [6, 7, 8],
    ]
    compressor = TriadicCompression(triads)

    # Example sector compression (for Yukawa calculation)
    R = sample_weights()
    K_sector = compressor.sector_kernel(finite.D_F, R)
    Y_sector = compressor.compress(K_sector)

    print("Effective 3×3 Yukawa (example):")
    print(Y_sector, "\n")

    # ------------------------------------------------------------
    # 5. Product triple
    # ------------------------------------------------------------
    product = AlignmentSpectralTriple(geom=geom, finite=finite, gamma_F=finite_gamma)
    D = product.build_product_dirac(modes)

    print("Product Dirac operator D shape:", D.shape)

    # ------------------------------------------------------------
    # 6. Example inner fluctuation
    # ------------------------------------------------------------
    gauge = GaugeSector()
    size = D.shape[0]

    # A U(1) block (placeholder)
    A = gauge.gauge_field_U1(size=size, coeff=0.01)

    D_fluct = product.fluctuated_dirac(D, A)
    print("Fluctuated Dirac constructed.\n")

    # ------------------------------------------------------------
    # 7. Spectral Action (v5 dev: gaussian / heat-kernel score)
    # ------------------------------------------------------------
    spectral = SpectralActionV3(Lambda=1.0, mode="gaussian")

    # Bosonic-only metrics
    metrics_bos = spectral.evaluate(D_fluct)
    print("\n--- SpectralActionV3 (bosonic) ---")
    print("dim:", metrics_bos["dim"])
    print("Srel (gaussian):", metrics_bos["Srel"])
    print("eigvals_preview:", metrics_bos["eigvals_preview"])
    print("a2_norm:", metrics_bos["a2_norm"])
    print("a4_norm:", metrics_bos["a4_norm"])

    # Fermionic vector for demonstration (optional)
    psi = np.random.randn(size) + 1j * np.random.randn(size)

    metrics_tot = spectral.evaluate(D_fluct, psi=psi)
    print("\n--- SpectralActionV3 (bos + ferm) ---")
    print("Srel (bosonic):", metrics_tot["Srel"])
    print("Sferm:", metrics_tot["Sferm"])
    print("Stotal:", metrics_tot["Stotal"])

    print("\n=== Pipeline Complete ===\n")


# ------------------------------------------------------------
# Entry point
# ------------------------------------------------------------
if __name__ == "__main__":
    main()

"""
=== Alignment Spectral Triple v2.0 Pipeline ===

Effective 3×3 Yukawa (example):
[[ 1.00000039e+00+1.29032335e-23j  5.22268121e-10-1.26086678e-09j
  -9.65025655e-10-9.65025655e-10j]
 [ 5.22268121e-10+1.26086678e-09j  1.00000000e+00+4.38579393e-27j
   5.22268121e-10-1.26086678e-09j]
 [-9.65025655e-10+9.65025655e-10j  5.22268121e-10+1.26086678e-09j
   1.00000000e+00+2.97145844e-27j]] 

Product Dirac operator D shape: (54, 54)
Fluctuated Dirac constructed.

Bosonic spectral action: 5.4000000000000004e+57
Fermionic action: 497.0713623386307
Total action: 5.4000000000000004e+57

=== Pipeline Complete ===
"""