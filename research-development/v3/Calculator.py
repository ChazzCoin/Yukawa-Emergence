#!/usr/bin/env python3
# ------------------------------------------------------------
#  Alignment Spectral Triple v3.2 — Finite Flavor Geometry
#  Z3 x Z3 heap-based internal group, character-sum kernel K(g,h),
#  triadic cosets, universal D_F, compressed Yukawa Y, mixing.
# ------------------------------------------------------------

import numpy as np
from numpy.linalg import eigh  # for Hermitian eigendecomposition

# ============================================================
# 1. CONFIGURATION: ALL MODEL CHOICES LIVE HERE
# ============================================================

class AlignmentV32Config:
    """
    Configuration for the v3.2 heap-based finite flavor sector.

    - Internal group: G = Z_3 x Z_3
    - Triads: cosets of a subgroup H of order 3
    - Kernel: K(g,h) = sum_a w_a chi_{p_a,q_a}(g - h)
    - Compression: flat or character-weighted triad modes
    """

    # ----- group structure -----
    # Z3 x Z3 = {(i,j) | i,j in {0,1,2}}
    group_elements = [(i, j) for i in range(3) for j in range(3)]

    # subgroup H of order 3
    # we choose H = {(0,0), (1,1), (2,2)} ⊂ Z3 x Z3
    subgroup_H = [(0, 0), (1, 1), (2, 2)]

    # coset shifts to generate triads: H, H+(1,0), H+(0,1)
    triad_shifts = [(0, 0), (1, 0), (0, 1)]

    # --------------------------------------------------------
    # CHARACTER SET FOR THE KERNEL
    # --------------------------------------------------------
    # General character:
    #   chi_{p,q}(i,j) = exp(2π i / 3 * (p*i + q*j)),  (p,q in {0,1,2})
    #
    # We now allow a SUM of characters:
    #   K(g,h) = sum_a w_a * chi_{p_a,q_a}(g - h)
    #
    # Listed as: kernel_characters = [(p, q, weight), ...]
    #
    # Example below: a superposition of two nontrivial characters (1,0) and (0,1)
    # This breaks the extreme symmetry of a single-character kernel.
    kernel_characters = [
        (1, 0, 1.0),
        (0, 1, 1.0),
    ]

    # --------------------------------------------------------
    # COMPRESSION MODE FOR TRIADS
    # --------------------------------------------------------
    # Options:
    #   "flat"               : (S psi)_i = (1/sqrt(3)) sum_{g ∈ T_i} psi(g)
    #   "character_weighted" : (S psi)_i = (1/√N_i) ∑_{g ∈ T_i} chi_i(g) psi(g)
    #
    # character_weighted uses one character per triad, specified below.
    compression_mode = "character_weighted"

    # Characters used for each triad row in "character_weighted" mode.
    # compression_characters[i] = (p_i, q_i) for triad T_i.
    #
    # Example choice:
    #   T1 uses trivial character (0,0) → flat sum but normalized;
    #   T2 uses (1,0), T3 uses (0,1) → internal phases in those triads.
    compression_characters = [
        (0, 0),  # triad 1
        (1, 0),  # triad 2
        (0, 1),  # triad 3
    ]

    # ----- physical scale (optional) -----
    higgs_vev = 174.0  # GeV, used if you want to convert eigenvalues to masses


# ============================================================
# 2. GROUP ARITHMETIC ON Z3 x Z3
# ============================================================

def add_g(a, b):
    """Add two elements of Z3 x Z3: a, b are (i,j) pairs."""
    return ((a[0] + b[0]) % 3, (a[1] + b[1]) % 3)

def neg_g(a):
    """Additive inverse in Z3 x Z3."""
    return ((-a[0]) % 3, (-a[1]) % 3)

def sub_g(a, b):
    """Difference a - b in Z3 x Z3."""
    return add_g(a, neg_g(b))


# ============================================================
# 3. CHARACTERS AND KERNEL
# ============================================================

def character_chi_pq(g, p, q):
    """
    Character of G = Z3 x Z3 with parameters (p,q):

        chi_{p,q}(i,j) = exp(2π i / 3 * (p*i + q*j)).
    """
    i, j = g
    phase = (p * i + q * j) % 3
    return np.exp(2j * np.pi * phase / 3.0)


def build_flavor_kernel(cfg):
    """
    Build the universal 9x9 flavor kernel K with entries:

        K_{g,h} = sum_a w_a * chi_{p_a,q_a}(g - h),

    where g,h ∈ G = Z3 x Z3, indexed in the order given by cfg.group_elements.
    """
    G = cfg.group_elements
    n = len(G)
    K = np.zeros((n, n), dtype=complex)

    chars = cfg.kernel_characters  # list of (p,q,weight)

    for i, g in enumerate(G):
        for j, h in enumerate(G):
            diff = sub_g(g, h)
            val = 0.0 + 0.0j
            for p, q, w in chars:
                val += w * character_chi_pq(diff, p, q)
            K[i, j] = val

    # Optional: enforce Hermitian symmetry numerically
    K = 0.5 * (K + K.conj().T)
    return K


# ============================================================
# 4. TRIADS AS COSETS
# ============================================================

def build_triad_indices(cfg):
    """
    Build triads T1, T2, T3 as cosets of H in G.

    We use cfg.subgroup_H and cfg.triad_shifts. Each triad is a list
    of indices into cfg.group_elements.
    """
    G = cfg.group_elements
    H = cfg.subgroup_H
    shifts = cfg.triad_shifts

    # maps element -> index in G list
    index_of = {g: idx for idx, g in enumerate(G)}

    triads = []
    for shift in shifts:
        triad_elems = [add_g(h, shift) for h in H]
        triad_indices = [index_of[g] for g in triad_elems]
        triads.append(triad_indices)

    return triads  # list of 3 lists, each of length 3


def build_compression_matrix(cfg, triads):
    """
    Build S : C^9 -> C^3, which maps site basis to generation basis.

    Mode = "flat":
        (S psi)_i = (1/sqrt(3)) sum_{g ∈ T_i} psi(g).

    Mode = "character_weighted":
        (S psi)_i = (1/√N_i) ∑_{g ∈ T_i} chi_i(g) psi(g),
        where chi_i is the character with parameters (p_i,q_i) given
        in cfg.compression_characters[i], and N_i is the number of sites
        in T_i (here 3).
    """
    G = cfg.group_elements
    n_sites = len(G)
    n_triads = len(triads)
    S = np.zeros((n_triads, n_sites), dtype=complex)

    mode = cfg.compression_mode

    if mode == "flat":
        # Original uniform coset averages (leads to Y≈0 for pure nontrivial chi)
        for i, triad in enumerate(triads):
            for idx in triad:
                S[i, idx] = 1.0 / np.sqrt(len(triad))

    elif mode == "character_weighted":
        # Each row i uses its own character (p_i, q_i)
        comp_chars = cfg.compression_characters
        assert len(comp_chars) == n_triads, \
            "Need one compression character (p,q) per triad."

        for i, triad in enumerate(triads):
            p_i, q_i = comp_chars[i]
            # normalization: sqrt(sum_{g∈T_i} |chi_i(g)|^2) = sqrt(3)
            norm = np.sqrt(len(triad))
            for idx in triad:
                g = G[idx]
                chi_val = character_chi_pq(g, p_i, q_i)
                S[i, idx] = chi_val / norm
    else:
        raise ValueError(f"Unknown compression_mode: {mode}")

    return S


# ============================================================
# 5. EFFECTIVE YUKAWA AND SPECTRUM
# ============================================================

def build_effective_yukawa(K, S):
    """
    Effective 3x3 Yukawa operator:

        Y = S K S^\dagger.

    K : 9x9 kernel, S : 3x9 compression matrix.
    """
    return S @ K @ S.conj().T


def diagonalize_yukawa(Y, vev=None):
    """
    Diagonalize Y (3x3 Hermitian):

        Y v_i = λ_i v_i,

    where λ_i are eigenvalues and v_i columns of U.

    If vev is not None, we also compute "masses" = |λ_i| * vev.
    """
    eigvals, U = eigh(Y)  # Hermitian eigendecomposition
    masses = np.abs(eigvals) * vev if vev is not None else None
    return eigvals, U, masses


# ============================================================
# 6. MAIN ROUTINE
# ============================================================

def run_alignment_v32_demo():
    cfg = AlignmentV32Config()

    print("\n==============================================")
    print(" Alignment Spectral Triple v3.2 — Flavor Demo ")
    print(" Internal group: Z3 x Z3 with heap-based kernel")
    print("==============================================\n")

    # 1) group and triads
    print("Group elements G = Z3 x Z3 (in order):")
    for idx, g in enumerate(cfg.group_elements):
        print(f"  {idx}: {g}")
    print()

    triads = build_triad_indices(cfg)
    print("Triads (cosets of H) by indices:")
    for i, triad in enumerate(triads, start=1):
        elems = [cfg.group_elements[idx] for idx in triad]
        print(f"  T_{i}: indices {triad}, elements {elems}")
    print()

    # 2) universal kernel
    K = build_flavor_kernel(cfg)
    print("Universal flavor kernel K (9x9):")
    print(K)
    print()

    # 3) compression matrix
    S = build_compression_matrix(cfg, triads)
    print("Compression matrix S (3x9), mode =", cfg.compression_mode)
    print(S)
    print("Check S S^† (should be identity on C^3):")
    print(S @ S.conj().T)
    print()

    # 4) effective Yukawa
    Y = build_effective_yukawa(K, S)
    print("Effective Yukawa operator Y = S K S^† (3x3):")
    print(Y)
    print()

    # 5) spectrum and mixing
    eigvals, U, masses = diagonalize_yukawa(Y, vev=cfg.higgs_vev)

    print("Eigenvalues λ_i of Y (dimensionless):")
    print(eigvals)
    print()

    print(f"Effective 'masses' m_i = |λ_i| * v with v = {cfg.higgs_vev} GeV:")
    print(masses)
    print()

    print("Unitary matrix U (columns are eigenvectors of Y):")
    print(U)
    print()

    print("Absolute values |U| (mixing pattern in generation space):")
    print(np.abs(U))
    print()

    print("Demo complete.\n")


# ============================================================
# ENTRY POINT
# ============================================================

if __name__ == "__main__":
    run_alignment_v32_demo()

"""
==============================================
 Alignment Spectral Triple v3.2 — Flavor Demo 
 Internal group: Z3 x Z3 with heap-based kernel
==============================================

Group elements G = Z3 x Z3 (in order):
  0: (0, 0)
  1: (0, 1)
  2: (0, 2)
  3: (1, 0)
  4: (1, 1)
  5: (1, 2)
  6: (2, 0)
  7: (2, 1)
  8: (2, 2)

Triads (cosets of H) by indices:
  T_1: indices [0, 4, 8], elements [(0, 0), (1, 1), (2, 2)]
  T_2: indices [3, 7, 2], elements [(1, 0), (2, 1), (0, 2)]
  T_3: indices [1, 5, 6], elements [(0, 1), (1, 2), (2, 0)]

Universal flavor kernel K (9x9):
[[ 2. +0.j          0.5-0.8660254j   0.5+0.8660254j   0.5-0.8660254j
  -1. -1.73205081j -1. +0.j          0.5+0.8660254j  -1. +0.j
  -1. +1.73205081j]
 [ 0.5+0.8660254j   2. +0.j          0.5-0.8660254j  -1. +0.j
   0.5-0.8660254j  -1. -1.73205081j -1. +1.73205081j  0.5+0.8660254j
  -1. +0.j        ]
 [ 0.5-0.8660254j   0.5+0.8660254j   2. +0.j         -1. -1.73205081j
  -1. +0.j          0.5-0.8660254j  -1. +0.j         -1. +1.73205081j
   0.5+0.8660254j ]
 [ 0.5+0.8660254j  -1. +0.j         -1. +1.73205081j  2. +0.j
   0.5-0.8660254j   0.5+0.8660254j   0.5-0.8660254j  -1. -1.73205081j
  -1. +0.j        ]
 [-1. +1.73205081j  0.5+0.8660254j  -1. +0.j          0.5+0.8660254j
   2. +0.j          0.5-0.8660254j  -1. +0.j          0.5-0.8660254j
  -1. -1.73205081j]
 [-1. +0.j         -1. +1.73205081j  0.5+0.8660254j   0.5-0.8660254j
   0.5+0.8660254j   2. +0.j         -1. -1.73205081j -1. +0.j
   0.5-0.8660254j ]
 [ 0.5-0.8660254j  -1. -1.73205081j -1. +0.j          0.5+0.8660254j
  -1. +0.j         -1. +1.73205081j  2. +0.j          0.5-0.8660254j
   0.5+0.8660254j ]
 [-1. +0.j          0.5-0.8660254j  -1. -1.73205081j -1. +1.73205081j
   0.5+0.8660254j  -1. +0.j          0.5+0.8660254j   2. +0.j
   0.5-0.8660254j ]
 [-1. -1.73205081j -1. +0.j          0.5-0.8660254j  -1. +0.j
  -1. +1.73205081j  0.5+0.8660254j   0.5-0.8660254j   0.5+0.8660254j
   2. +0.j        ]]

Compression matrix S (3x9), mode = character_weighted
[[ 0.57735027+0.j   0.        +0.j   0.        +0.j   0.        +0.j
   0.57735027+0.j   0.        +0.j   0.        +0.j   0.        +0.j
   0.57735027+0.j ]
 [ 0.        +0.j   0.        +0.j   0.57735027+0.j  -0.28867513+0.5j
   0.        +0.j   0.        +0.j   0.        +0.j  -0.28867513-0.5j
   0.        +0.j ]
 [ 0.        +0.j  -0.28867513+0.5j  0.        +0.j   0.        +0.j
   0.        +0.j  -0.28867513-0.5j  0.57735027+0.j   0.        +0.j
   0.        +0.j ]]
Check S S^† (should be identity on C^3):
[[1.+0.j 0.+0.j 0.+0.j]
 [0.+0.j 1.+0.j 0.+0.j]
 [0.+0.j 0.+0.j 1.+0.j]]

Effective Yukawa operator Y = S K S^† (3x3):
[[-4.04424761e-16-4.41217410e-17j  2.47902802e-18-3.98279385e-17j
   2.47902802e-18-3.98279385e-17j]
 [ 6.40987562e-17-5.58611185e-18j  5.76888806e-16+6.26978536e-17j
  -2.24345647e-16-3.13489268e-17j]
 [ 6.40987562e-17-5.58611185e-18j -2.24345647e-16-3.13489268e-17j
   5.76888806e-16+6.26978536e-17j]]

Eigenvalues λ_i of Y (dimensionless):
[-4.15220315e-16  3.61125548e-16  8.03447618e-16]

Effective 'masses' m_i = |λ_i| * v with v = 174.0 GeV:
[7.22483349e-14 6.28358453e-14 1.39799885e-13]

Unitary matrix U (columns are eigenvectors of Y):
[[ 0.99303093+0.j          0.11773453-0.j          0.00530508+0.j        ]
 [-0.08280674+0.00937152j  0.69343814-0.11051167j  0.11084486+0.69835474j]
 [-0.08317953+0.00509394j  0.70210107-0.01110751j -0.01162983-0.70700118j]]

Absolute values |U| (mixing pattern in generation space):
[[0.99303093 0.11773453 0.00530508]
 [0.08333536 0.70218893 0.70709683]
 [0.08333536 0.70218893 0.70709683]]

Demo complete.
"""