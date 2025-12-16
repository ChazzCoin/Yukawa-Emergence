#!/usr/bin/env python3
import numpy as np
from numpy.linalg import eigh
from scipy.linalg import expm

def scan_seesaw_scale(pipe, M_values):
    """
    Scan heavy scale M and record light neutrino masses.
    """
    results = []

    for M in M_values:
        pipe.seesaw_M = M
        out = pipe.run()

        # extract Majorana-like light eigenvalues
        majorana = out["majorana_blocks"]
        assert len(majorana) == 1, "Expected exactly one Majorana block"

        block, light_eigs = majorana[0]
        results.append((M, light_eigs.copy()))

    return results
def set_geometry_weights(cfg, g_diag=1.0, g_off=1.0):
    """
    Controlled two-parameter geometry deformation.

    g_diag > 1  : stronger pinning of diagonal triad (Dirac rigidity)
    g_off  < 1  : softer off-diagonal geometry (PMNS enhancement)
    """
    base = {
        (0,0): 1.00, (0,1): 0.92, (0,2): 0.85,
        (1,0): 0.95, (1,1): 1.10, (1,2): 0.88,
        (2,0): 0.80, (2,1): 0.90, (2,2): 1.05,
    }

    diag = {(0,0), (1,1), (2,2)}

    for g in base:
        if g in diag:
            base[g] *= g_diag
        else:
            base[g] *= g_off

    cfg.geometry_weights = base


def angles_from_absU(absU):
    s13 = absU[0,2]
    c13 = np.sqrt(max(0.0, 1.0 - s13**2))

    s12 = absU[0,1] / c13 if c13 > 0 else np.nan
    s23 = absU[1,2] / c13 if c13 > 0 else np.nan

    th12 = np.degrees(np.arcsin(np.clip(s12, -1, 1)))
    th13 = np.degrees(np.arcsin(np.clip(s13, -1, 1)))
    th23 = np.degrees(np.arcsin(np.clip(s23, -1, 1)))

    return th12, th13, th23
def scan_geometry(pipe, g_diag_vals, g_off_vals):
    """
    Perform a controlled geometry scan and extract mixing angles.
    """
    results = []

    for g_diag in g_diag_vals:
        for g_off in g_off_vals:
            set_geometry_weights(pipe.cfg, g_diag=g_diag, g_off=g_off)

            out = pipe.run()
            th12, th13, th23 = angles_from_absU(out["absU"])

            results.append({
                "g_diag": g_diag,
                "g_off": g_off,
                "theta_12": th12,
                "theta_13": th13,
                "theta_23": th23,
                "blocks": out["blocks"],
            })

    return results
def apply_full_seesaw_3x3(self, Y):
    """
    Apply a full 3x3 seesaw to a single 3D harmonic block.
    Assumes Y is already the effective 3x3 lepton Yukawa.
    """
    n = Y.shape[0]
    assert n == 3, "Full seesaw requires 3x3 Y."

    zero = np.zeros_like(Y)
    MR = self.seesaw_M * np.eye(3)

    big = np.block([
        [zero, Y],
        [Y.conj().T, MR]
    ])

    eigvals, eigvecs = eigh(big)
    eigvals = np.real(eigvals)

    # sort by absolute value
    order = np.argsort(np.abs(eigvals))
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]

    # light neutrinos = first 3 eigenvalues
    light_masses = np.abs(eigvals[:3])
    light_vectors = eigvecs[:3, :3]

    return light_masses, light_vectors

# ============================================================
#  CONFIG (same as yours, unchanged)
# ============================================================

class AlignmentV33Config:
    group_elements = [(i, j) for i in range(3) for j in range(3)]
    subgroup_H = [(0, 0), (1, 1), (2, 2)]
    triad_shifts = [(0, 0), (1, 0), (0, 1)]

    kernel_characters = [
        (1, 0,  1.0),
        (0, 1,  0.6),
        (1, 1,  0.35),
    ]

    geometry_weights = {
        (0,0): 1.00, (0,1): 0.92, (0,2): 0.85,
        (1,0): 0.95, (1,1): 1.10, (1,2): 0.88,
        (2,0): 0.80, (2,1): 0.90, (2,2): 1.05,
    }

    damping_strength = 0.35

    compression_characters = [
        (0, 0),
        (1, 0),
        (1, 1),
    ]

    higgs_vev = 174.0


# ============================================================
#  PIPELINE
# ============================================================

class AlignmentPipeline:
    """
    AlignmentPipeline encapsulates the operator sequence:

      K --(misalignment flow)--> K_flow
        --(emergent C360 projector)--> K_proj
        --(compression)--> Y
        --(block diagnosis)--> blocks on Y
        --(conditional seesaw)--> neutrino-like light eigenvalues

    This makes quark/lepton emergence automatic:
      - 1D harmonic blocks => Dirac-like (no seesaw)
      - >=2D harmonic blocks => closure failure => seesaw extension
    """

    def __init__(
        self,
        cfg: AlignmentV33Config,
        beta: float = 1.5,
        rel_cut: float = 0.15,
        tol_rel_blocks: float = 0.03,
        seesaw_M: float = 1e6,
    ):
        self.cfg = cfg
        self.beta = float(beta)
        self.rel_cut = float(rel_cut)
        self.tol_rel_blocks = float(tol_rel_blocks)
        self.seesaw_M = float(seesaw_M)

        # caches (filled by run())
        self.triads = None
        self.K = None
        self.S = None
        self.K_flow = None
        self.P_C360 = None
        self.kept_indices = None
        self.K_proj = None
        self.Y = None

    # ---------- group + characters ----------

    @staticmethod
    def add_g(a, b):
        return ((a[0]+b[0]) % 3, (a[1]+b[1]) % 3)

    @staticmethod
    def sub_g(a, b):
        return ((a[0]-b[0]) % 3, (a[1]-b[1]) % 3)

    @staticmethod
    def chi(g, p, q):
        i, j = g
        return np.exp(2j*np.pi*(p*i+q*j)/3.0)

    # ---------- build K (geometry-selected, non-convolution) ----------

    def build_kernel(self):
        cfg = self.cfg
        G = cfg.group_elements
        n = len(G)
        K = np.zeros((n, n), dtype=complex)

        for a, g in enumerate(G):
            for b, h in enumerate(G):
                # spectral part
                F = sum(
                    w * self.chi(self.sub_g(g, h), p, q)
                    for (p, q, w) in cfg.kernel_characters
                )

                # geometry weights (B)
                alpha_g = cfg.geometry_weights[g]
                alpha_h = cfg.geometry_weights[h]

                # misalignment damping W(g,h)
                dist = min(abs(g[0]-h[0]), 3-abs(g[0]-h[0])) + \
                       min(abs(g[1]-h[1]), 3-abs(g[1]-h[1]))
                W = np.exp(-cfg.damping_strength * dist)

                K[a, b] = alpha_g * F * np.conj(alpha_h) * W

        # enforce Hermitian
        return 0.5 * (K + K.conj().T)

    # ---------- triads + compression S ----------

    def build_triads(self):
        cfg = self.cfg
        index = {g: i for i, g in enumerate(cfg.group_elements)}
        triads = []
        for s in cfg.triad_shifts:
            triads.append([index[self.add_g(h, s)] for h in cfg.subgroup_H])
        return triads

    def build_S(self, triads):
        cfg = self.cfg
        G = cfg.group_elements
        S = np.zeros((3, 9), dtype=complex)

        for i, triad in enumerate(triads):
            p, q = cfg.compression_characters[i]
            for idx in triad:
                g = G[idx]
                S[i, idx] = self.chi(g, p, q) / np.sqrt(3)
        return S

    # ---------- operators: flow, C360 projector, compression ----------

    def misalignment_flow(self, K):
        return expm(-self.beta * K)

    def emergent_C360_projector(self, K):
        evals, evecs = eigh(K)
        flowed = np.exp(-self.beta * evals)

        max_val = flowed.max()
        keep = flowed >= self.rel_cut * max_val
        kept = np.where(keep)[0]

        P = np.zeros_like(K, dtype=complex)
        for i in kept:
            v = evecs[:, i:i+1]
            P += v @ v.conj().T

        return P, kept.tolist()

    @staticmethod
    def effective_yukawa(K_like, S):
        return S @ K_like @ S.conj().T

    # ---------- block diagnosis on Y ----------

    def harmonic_blocks_on_Y(self, Y):
        evals, _ = eigh(Y)
        evals = np.sort(np.real(evals))

        blocks = []
        block = [0]
        for i in range(1, len(evals)):
            scale = max(1.0, abs(evals[i-1]), abs(evals[i]))
            if abs(evals[i] - evals[i-1]) <= self.tol_rel_blocks * scale:
                block.append(i)
            else:
                blocks.append(block)
                block = [i]
        blocks.append(block)
        return blocks, evals

    # ---------- conditional seesaw on >=2D blocks ----------

    @staticmethod
    def _seesaw_light_eigs_for_block(Y, block, M):
        Yb = Y[np.ix_(block, block)]
        n = len(block)

        zero = np.zeros_like(Yb)
        MR = M * np.eye(n)

        big = np.block([
            [zero, Yb],
            [Yb.conj().T, MR]
        ])

        eigvals, _ = eigh(big)
        eigvals = np.sort(np.abs(eigvals))
        return eigvals[:n]

    def apply_conditional_seesaw(self, Y, blocks):
        """
        Returns:
          - dirac_blocks: list of (block, eigenvalues) for 1D blocks
          - majorana_blocks: list of (block, light_eigs_after_seesaw) for >=2D blocks
        """
        dirac = []
        majorana = []

        # we use eigenvalues of the sub-block as "Dirac-like" for 1D
        # (you can later map these to charged lepton vs quark readouts)
        for block in blocks:
            if len(block) == 1:
                i = block[0]
                dirac.append((block, np.array([np.real(eigh(Y)[0][i])])))
            else:
                light = self._seesaw_light_eigs_for_block(Y, block, self.seesaw_M)
                majorana.append((block, light))

        return dirac, majorana

    def apply_full_seesaw_3x3(self, Y):
        """
        Apply a full 3x3 seesaw to a single 3D harmonic block.
        Assumes Y is already the effective 3x3 lepton Yukawa.
        """
        n = Y.shape[0]
        assert n == 3, "Full seesaw requires 3x3 Y."

        zero = np.zeros_like(Y)
        MR = self.seesaw_M * np.eye(3)

        big = np.block([
            [zero, Y],
            [Y.conj().T, MR]
        ])

        eigvals, eigvecs = eigh(big)
        eigvals = np.real(eigvals)

        # sort by absolute value
        order = np.argsort(np.abs(eigvals))
        eigvals = eigvals[order]
        eigvecs = eigvecs[:, order]

        # light neutrinos = first 3 eigenvalues
        light_masses = np.abs(eigvals[:3])
        light_vectors = eigvecs[:3, :3]

        return light_masses, light_vectors

    # ---------- run the whole pipeline ----------

    def run(self):
        cfg = self.cfg

        self.triads = self.build_triads()
        self.K = self.build_kernel()
        self.S = self.build_S(self.triads)

        self.K_flow = self.misalignment_flow(self.K)
        self.P_C360, self.kept_indices = self.emergent_C360_projector(self.K)
        self.K_proj = self.P_C360 @ self.K_flow @ self.P_C360

        self.Y = self.effective_yukawa(self.K_proj, self.S)

        # Y spectrum
        y_vals, U = eigh(self.Y)
        masses = np.abs(y_vals) * cfg.higgs_vev

        # block structure (this is the “emergence diagnosis”)
        blocks, y_sorted = self.harmonic_blocks_on_Y(self.Y)

        # seesaw only where needed
        dirac_blocks, majorana_blocks = self.apply_conditional_seesaw(self.Y, blocks)

        return {
            "Y": self.Y,
            "Y_eigvals": y_vals,
            "masses_GeV": masses,
            "U": U,
            "absU": np.abs(U),
            "blocks": blocks,
            "Y_sorted_eigvals": y_sorted,
            "dirac_blocks": dirac_blocks,
            "majorana_blocks": majorana_blocks,
            "kept_indices_C360": self.kept_indices,
        }



# ============================================================
#  Example usage
# ============================================================

def main_seesaw_scan():
    # ---------------------------------------
    # Scan seesaw scale
    # ---------------------------------------

    cfg = AlignmentV33Config()
    pipe = AlignmentPipeline(cfg, beta=1.5, rel_cut=0.15, tol_rel_blocks=0.03)

    out = pipe.run()

    print("\nY eigenvalues:", out["Y_eigvals"])
    print("Masses [GeV]:", out["masses_GeV"])
    print("\n|U|:\n", out["absU"])

    print("\nHarmonic blocks (on Y):", out["blocks"])
    print("C360 kept indices (on K):", out["kept_indices_C360"])

    print("\nDirac-like blocks (1D):", out["dirac_blocks"])
    print("Majorana-like blocks (>=2D) light eigs:", out["majorana_blocks"])

    M_values = [1e3, 1e4, 1e5, 1e6, 1e7]
    scan = scan_seesaw_scale(pipe, M_values)

    print("\nSeesaw scale scan:")
    for M, eigs in scan:
        print(f"M = {M:.1e}  -> light eigenvalues = {eigs}")
        print(f"   M * m_light ≈ {M * eigs}")

    g_diag_vals = [0.9, 1.0, 1.1, 1.2]
    g_off_vals = [0.7, 0.85, 1.0, 1.15]
    results = scan_geometry(pipe, g_diag_vals, g_off_vals)
    for r in results:
        print(
            f"g_diag={r['g_diag']:.2f}, g_off={r['g_off']:.2f} | "
            f"θ12={r['theta_12']:.1f}°, "
            f"θ13={r['theta_13']:.1f}°, "
            f"θ23={r['theta_23']:.1f}° | "
            f"blocks={r['blocks']}"
        )
    print("\n--Main Seesaw Scan Complete--\n")

def main():
    cfg = AlignmentV33Config()
    pipe = AlignmentPipeline(cfg, beta=1.5, rel_cut=0.15, tol_rel_blocks=0.03, seesaw_M=1e6)

    out = pipe.run()

    print("\nY eigenvalues:", out["Y_eigvals"])
    print("Masses [GeV]:", out["masses_GeV"])
    print("\n|U|:\n", out["absU"])

    print("\nHarmonic blocks (on Y):", out["blocks"])
    print("C360 kept indices (on K):", out["kept_indices_C360"])

    print("\nDirac-like blocks (1D):", out["dirac_blocks"])
    print("Majorana-like blocks (>=2D) light eigs:", out["majorana_blocks"])

def run_lepton_sector_full_3nu():
    cfg = AlignmentV33Config()

    # Lock geometry into 3D harmonic regime
    set_geometry_weights(cfg, g_diag=0.90, g_off=1.00)

    pipe = AlignmentPipeline(
        cfg,
        beta=1.5,
        rel_cut=0.15,
        tol_rel_blocks=0.03,
        seesaw_M=1e6
    )

    out = pipe.run()

    print("\nLepton-sector harmonic blocks:", out["blocks"])

    # sanity check
    assert out["blocks"] == [[0,1,2]], "Not in 3D harmonic regime."

    # apply full seesaw
    light_masses, light_vectors = pipe.apply_full_seesaw_3x3(out["Y"])

    print("\nLight neutrino masses (dimensionless):")
    print(light_masses)

    print("\nNeutrino mixing matrix |U_PMNS|:")
    print(np.abs(light_vectors))

    return light_masses, light_vectors

if __name__ == "__main__":
    main_seesaw_scan()
    run_lepton_sector_full_3nu()

"""

Y eigenvalues: [0.20443467 0.30696381 0.31326296]
Masses [GeV]: [35.57163219 53.41170349 54.50775485]

|U|:
 [[0.99858093 0.04133058 0.03358422]
 [0.02739089 0.89481925 0.44558752]
 [0.04567117 0.44451129 0.89460822]]

Harmonic blocks (on Y): [[0], [1, 2]]
C360 kept indices (on K): [0, 1, 2, 3, 4, 5]

Dirac-like blocks (1D): [([0], array([0.20443467]))]
Majorana-like blocks (>=2D) light eigs: [([1, 2], array([9.41276264e-08, 9.80979446e-08]))]

Seesaw scale scan:
M = 1.0e+03  -> light eigenvalues = [9.41172042e-05 9.80589011e-05]
   M * m_light ≈ [0.0941172 0.0980589]
M = 1.0e+04  -> light eigenvalues = [9.41172093e-06 9.80588972e-06]
   M * m_light ≈ [0.09411721 0.0980589 ]
M = 1.0e+05  -> light eigenvalues = [9.41176667e-07 9.80605911e-07]
   M * m_light ≈ [0.09411767 0.09806059]
M = 1.0e+06  -> light eigenvalues = [9.41276264e-08 9.80979446e-08]
   M * m_light ≈ [0.09412763 0.09809794]
M = 1.0e+07  -> light eigenvalues = [9.48507732e-09 1.19789497e-08]
   M * m_light ≈ [0.09485077 0.1197895 ]
g_diag=0.90, g_off=0.70 | θ12=1.4°, θ13=0.1°, θ23=87.2° | blocks=[[0], [1], [2]]
g_diag=0.90, g_off=0.85 | θ12=2.8°, θ13=0.7°, θ23=85.6° | blocks=[[0], [1], [2]]
g_diag=0.90, g_off=1.00 | θ12=11.1°, θ13=3.4°, θ23=60.8° | blocks=[[0, 1, 2]]
g_diag=0.90, g_off=1.15 | θ12=76.1°, θ13=78.4°, θ23=8.9° | blocks=[[0, 1, 2]]
g_diag=1.00, g_off=0.70 | θ12=1.0°, θ13=0.3°, θ23=87.3° | blocks=[[0], [1], [2]]
g_diag=1.00, g_off=0.85 | θ12=1.6°, θ13=0.7°, θ23=86.0° | blocks=[[0], [1], [2]]
g_diag=1.00, g_off=1.00 | θ12=2.4°, θ13=1.9°, θ23=26.5° | blocks=[[0], [1, 2]]
g_diag=1.00, g_off=1.15 | θ12=6.5°, θ13=4.4°, θ23=3.1° | blocks=[[0, 1], [2]]
g_diag=1.10, g_off=0.70 | θ12=0.9°, θ13=0.4°, θ23=87.4° | blocks=[[0], [1], [2]]
g_diag=1.10, g_off=0.85 | θ12=1.3°, θ13=0.7°, θ23=86.3° | blocks=[[0], [1], [2]]
g_diag=1.10, g_off=1.00 | θ12=1.4°, θ13=1.6°, θ23=10.7° | blocks=[[0], [1, 2]]
g_diag=1.10, g_off=1.15 | θ12=2.3°, θ13=2.5°, θ23=2.5° | blocks=[[0], [1], [2]]
g_diag=1.20, g_off=0.70 | θ12=0.9°, θ13=0.8°, θ23=87.5° | blocks=[[0], [1], [2]]
g_diag=1.20, g_off=0.85 | θ12=1.2°, θ13=1.2°, θ23=86.6° | blocks=[[0], [1], [2]]
g_diag=1.20, g_off=1.00 | θ12=1.8°, θ13=1.4°, θ23=5.8° | blocks=[[0], [1, 2]]
g_diag=1.20, g_off=1.15 | θ12=2.3°, θ13=1.8°, θ23=2.3° | blocks=[[0], [1], [2]]

--Main Seesaw Scan Complete--


Lepton-sector harmonic blocks: [[0, 1, 2]]

Light neutrino masses (dimensionless):
[7.55617404e-08 8.72520849e-08 9.13414306e-08]

Neutrino mixing matrix |U_PMNS|:
[[0.97898577 0.19406889 0.0626428 ]
 [0.0904764  0.48103731 0.872019  ]
 [0.18275907 0.85495051 0.48544696]]
"""