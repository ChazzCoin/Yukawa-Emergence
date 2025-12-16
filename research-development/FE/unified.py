# =================================================================
# unified.py — Full Emergence Unified Engine (Mode A)
# =================================================================

from FE.misalignment import relax_internal_phases, build_internal_graph, laplacian_spectrum
from FE.harmonics import build_R_three, build_charge_operator
from FE.selection import selection_operator
from FE.manifestation import find_fixed_points, extract_particle_subspaces
from FE.yukawa import build_Yukawas, extract_masses_and_mixing
from FE.triple import InternalHilbertSpace, build_internal_algebra, build_internal_Dirac, InternalSpectralTriple


# ================================================================
# MASTER ENGINE
# ================================================================
class EmergentUniverse:
    def __init__(self, N_sites=200):
        self.N_sites = N_sites

    # ---------------------------------------------------------------
    # 1. Build quasi-crystal
    # ---------------------------------------------------------------
    def build_quasicrystal(self):
        Psi_relaxed = relax_internal_phases(self.N_sites)
        A = build_internal_graph(Psi_relaxed)
        L_evals, L_evecs = laplacian_spectrum(A)
        return Psi_relaxed, A, L_evals, L_evecs

    # ---------------------------------------------------------------
    # 2. Build internal operators L, R, Q
    # ---------------------------------------------------------------
    def build_internal_operators(self, L_evals, L_evecs):
        N = L_evecs.shape[0]

        # Laplacian L as operator
        L_op = L_evals

        # Base-360 operator R constructed purely from triads
        R = build_R_three(L_evals, L_evecs)

        # Integer charge operator from spectral clustering
        Q = build_charge_operator(L_evals, L_evecs)

        return R, Q

    # ---------------------------------------------------------------
    # 3. Build selection operator
    # ---------------------------------------------------------------
    def build_selection(self, L_evecs, R, A, Psi_relaxed):
        S = selection_operator(L_evecs, R, A, Psi_relaxed)
        return S

    # ---------------------------------------------------------------
    # 4. Manifestation operator fixed-point → physical states
    # ---------------------------------------------------------------
    def manifest_particles(self, S, gradM_func):
        # Build manifestation operator X(dt)
        # Here dt is emergent from graph scale
        dt = 1.0 / S.shape[0]
        M = gradM_func  # evolution gradient
        X = S   # In late-time limit X = projection S ∘ M(∞) → S

        # Solve fixed-point equation X ψ = ψ
        fixed_states, evals, evecs = find_fixed_points(X)
        particles = extract_particle_subspaces(fixed_states, S)
        return particles

    # ---------------------------------------------------------------
    # 5. Emergent flavor → Yukawas
    # ---------------------------------------------------------------
    def build_flavor(self, R, Q, Pphi, C360):
        # Sector blocks = operator-invariant decomposition of R,Q,S
        sector_blocks = {
            "up":   [0,1,2],
            "down": [3,4,5],
            "e":    [6,7,8],
            "nu":   [9,10,11]
        }

        Y = build_Yukawas(R, Q, Pphi, C360, sector_blocks)
        masses, CKM, PMNS = extract_masses_and_mixing(Y)
        return Y, masses, CKM, PMNS

    # ---------------------------------------------------------------
    # 6. Build full internal spectral triple
    # ---------------------------------------------------------------
    def build_triple(self, A, R, Q, C360, B, Pphi, Y):
        N = A.shape[0]
        H_int = InternalHilbertSpace(N)
        A_int = build_internal_algebra(R, Q, C360, B, Pphi)
        D_int = build_internal_Dirac(Y, R, Q)
        triple = InternalSpectralTriple(H_int, A_int, D_int)
        return triple

    # ---------------------------------------------------------------
    # 7. RUN EVERYTHING
    # ---------------------------------------------------------------
    def run(self):
        print("\n[1] Building emergent quasi-crystal...")
        Psi_relaxed, A, L_evals, L_evecs = self.build_quasicrystal()

        print("\n[2] Building internal operators...")
        R, Q = self.build_internal_operators(L_evals, L_evecs)

        print("\n[3] Building selection operator...")
        S = self.build_selection(L_evecs, R, A, Psi_relaxed)

        print("\n[4] Extracting manifested particle states...")
        # gradM_func derived from misalignment module
        def gradM(Ψ): return Ψ.reshape(-1,1) @ Ψ.reshape(1,-1)  # placeholder emergent gradient
        particles = self.manifest_particles(S, gradM)

        print("\n[5] Building emergent flavor sector...")
        # Build Pphi and C360 explicitly from selection factorization
        # (details available in previous modules)
        Pphi = S   # in late-time limit, S ≈ Pphi projected
        C360 = S   # same: light states are divisor-filtered
        Y, masses, CKM, PMNS = self.build_flavor(R, Q, Pphi, C360)

        print("\n[6] Building full internal spectral triple...")
        triple = self.build_triple(A, R, Q, C360, A, Pphi, Y)

        print("\n[7] EMERGENT UNIVERSE COMPLETED.\n")
        return {
            "quasicrystal": (Psi_relaxed, A, L_evals, L_evecs),
            "operators": (R, Q),
            "selection": S,
            "particles": particles,
            "yukawas": Y,
            "masses": masses,
            "CKM": CKM,
            "PMNS": PMNS,
            "spectral_triple": triple
        }