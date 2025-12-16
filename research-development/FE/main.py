# main.py — sector-resolved unitaries

import time
import numpy as np

from FE.chi2 import compute_global_chi2
from FE.misalignment import (
    relax_internal_phases,
    build_internal_graph,
    laplacian_spectrum,
)

from FE.harmonics import (
    build_R_three,
    build_charge_operator,
)

from FE.evolution import (
    build_misalignment_operator,
    evolve_to_fixed_point,
)

from FE.selection import (
    build_sector_selection_operators,
)

from FE.yukawa import (
    build_yukawas_from_manifested_state,
)


def run_once(N_sites=360, seed=42):
    print(f"=== Running full-emergent internal pipeline (N_sites={N_sites}, seed={seed}) ===")
    t0 = time.perf_counter()

    # 1) Relax phases
    phi = relax_internal_phases(N_sites=N_sites, seed=seed)

    # 2) Internal graph + Laplacian
    A = build_internal_graph(phi)
    evals, evecs = laplacian_spectrum(A)
    L = np.diag(A.sum(axis=1)) - A

    # 3) Emergent R, Q, triad modes
    R, k_list, mode_indices = build_R_three(evals, evecs)
    Q, q_list = build_charge_operator(evals, evecs, mode_indices=mode_indices)

    # 4) Misalignment operator & evolution (site space)
    Pphi = np.eye(L.shape[0])
    M_op = build_misalignment_operator(L, Pphi)
    psi0 = np.random.randn(L.shape[0])
    psi_inf = evolve_to_fixed_point(M_op, psi0)

    # 5) Sector-resolved selection operators
    S_sectors = build_sector_selection_operators(phi, L, k_list)

    # 6) Manifested internal states per sector
    psi_sectors = {}
    for s in ["u", "d", "e", "nu"]:
        psi_s = S_sectors[s] @ psi_inf
        n = np.linalg.norm(psi_s)
        if n > 1e-15:
            psi_s = psi_s / n
        psi_sectors[s] = psi_s

    print("=== Manifested Internal State Norms (per sector) ===")
    for s in psi_sectors:
        print(f"{s}: {np.linalg.norm(psi_sectors[s])}")

    # 7) Fully emergent Yukawas from manifested sector states
    Y_u, Y_d, Y_e, Y_nu, spectra, mixings = \
        build_yukawas_from_manifested_state(psi_sectors, R, Q, evecs, mode_indices)

    print("\n=== Fully Emergent Mass Spectra ===")
    for k, v in spectra.items():
        print(f"{k}: {v}")

    print("\n=== Fully Emergent Mixing Matrices ===")
    for name, M in mixings.items():
        print(f"{name}:")
        print(M)
        print()

    chi2_total, chi2_terms = compute_global_chi2(spectra, mixings)

    print("\n=== Global χ² Diagnostic ===")
    print("Total χ²:", chi2_total)
    print("Per-term contributions:")
    for k, v in chi2_terms.items():
        print(f"  {k:12s}: {v:.4f}")

    t1 = time.perf_counter()
    print(f"Total run time: {t1 - t0:.3f} s")
    print("=== Done ===")

if __name__ == "__main__":
    run_once()

"""
RESULTS:
=== Running full-emergent internal pipeline (N_sites=360, seed=42) ===
=== Manifested Internal State Norms (per sector) ===
u: 1.0
d: 1.0000000000000002
e: 1.0
nu: 1.0

=== Fully Emergent Mass Spectra ===
u: [0.36509557 0.13450474 0.04954539]
d: [0.3660937  0.13450474 0.04941031]
e: [0.36233277 0.1336793  0.04930489]
nu: [0.36509557 0.13450474 0.04954539]

=== Fully Emergent Mixing Matrices ===
CKM:
[[-4.61143074e-02-5.01065375e-01j -8.21601839e-01-2.67913008e-01j
  -1.60624128e-15+3.37750128e-15j]
 [ 8.46893847e-01-1.71981898e-01j -1.03238664e-01+4.92478241e-01j
   2.75386549e-15+6.27953025e-16j]
 [ 2.88509593e-16+5.34749291e-16j -4.04288935e-15+2.29161018e-15j
  -4.08579684e-01-9.12722653e-01j]]

PMNS:
[[ 8.17795846e-01+3.00403524e-01j  4.73738656e-01-1.28605453e-01j
  -9.64587809e-17+5.83696855e-16j]
 [-4.24154736e-01-2.47104099e-01j  8.67729737e-01-7.79578507e-02j
   4.46176924e-16+6.90244868e-16j]
 [ 1.63348309e-16+3.84081468e-16j -4.52978049e-16+8.03820741e-16j
   9.60867910e-01+2.77006965e-01j]]


=== Global χ² Diagnostic ===
Total χ²: 20378402373.778545
Per-term contributions:
  u_m2/m3     : 132066.5936
  u_m1/m3     : 287579168.3803
  d_m2/m3     : 1007.4608
  d_m1/m3     : 48375.2902
  e_m2/m3     : 91476.3442
  e_m1/m3     : 20090545867.6485
  nu_m2/m3    : 34.0540
  nu_m1/m3    : 3950.4540
  CKM_θ12     : 327.1283
  CKM_θ23     : 25.0000
  CKM_θ13     : 25.0000
  PMNS_θ12    : 0.4247
  PMNS_θ23    : 25.0000
  PMNS_θ13    : 25.0000
Total run time: 9.291 s
=== Done ===

"""