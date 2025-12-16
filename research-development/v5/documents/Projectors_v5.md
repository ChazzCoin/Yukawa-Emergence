## Canonical (H_{\rm SM}) basis ordering (16 particle + 16 conjugate)

Fix **one-generation** internal particle space (H_{\rm SM}^{\rm part}\cong \mathbb C^{16}) and conjugate space (H_{\rm SM}^{\rm conj}\cong \mathbb C^{16}). Define
[
H_{\rm SM}=H_{\rm SM}^{\rm part}\oplus H_{\rm SM}^{\rm conj},
\qquad
H_{\rm int}=H_{\rm SM}\otimes H_{\rm flav},
\qquad \dim(H_{\rm int})=32N.
]

### Particle ordering (indices (0)–(15))

I will use the standard “8 left + 8 right” ordering inside the particle block:

**Left (8):**

1. (u_L^r)
2. (u_L^g)
3. (u_L^b)
4. (d_L^r)
5. (d_L^g)
6. (d_L^b)
7. (\nu_L)
8. (e_L)

**Right (8):**
9.  (u_R^r)
10. (u_R^g)
11. (u_R^b)
12. (d_R^r)
13. (d_R^g)
14. (d_R^b)
15. (e_R)
16. (\nu_R)

So: (Q_L) is (0)–(5), (L_L) is (6)–(7), (u_R) is (8)–(10), (d_R) is (11)–(13), (e_R) is (14), (\nu_R) is (15).

### Conjugate ordering (indices (16)–(31))

Take the conjugates in the **same order**:

* (16)–(31) correspond to ((u_L^{r})^c,\dots,(\nu_R)^c) in the same sequence as above.

Equivalently, for each particle basis vector (e_i) (with (i=0,\dots,15)), the conjugate partner is (e_{i+16}).

---

## Exact projector matrices on (H_{\rm SM})

Let (E_{ii}) denote the (32\times 32) matrix with a (1) at ((i,i)) and zeros elsewhere. For any index set (S\subset{0,\dots,31}), define
[
P_S := \sum_{i\in S} E_{ii}.
]

### Particle-sector projectors (acting on indices (0)–(15))

[
P_{Q_L}=P_{{0,1,2,3,4,5}},\qquad
P_{L_L}=P_{{6,7}},
]
[
P_{u_R}=P_{{8,9,10}},\qquad
P_{d_R}=P_{{11,12,13}},
]
[
P_{e_R}=P_{{14}},\qquad
P_{\nu_R}=P_{{15}}.
]

### Conjugate-sector projectors (indices (16)–(31))

[
P_{Q_L^c}=P_{{16,17,18,19,20,21}},\qquad
P_{L_L^c}=P_{{22,23}},
]
[
P_{u_R^c}=P_{{24,25,26}},\qquad
P_{d_R^c}=P_{{27,28,29}},
]
[
P_{e_R^c}=P_{{30}},\qquad
P_{\nu_R^c}=P_{{31}}.
]

### Chiral projectors on (H_{\rm SM})

Define the left/right projectors on the **particle** subspace:
[
P_L^{\rm part}=P_{{0,1,2,3,4,5,6,7}},\qquad
P_R^{\rm part}=P_{{8,9,10,11,12,13,14,15}}.
]
And on conjugates:
[
P_L^{\rm conj}=P_{{16,17,18,19,20,21,22,23}},\qquad
P_R^{\rm conj}=P_{{24,25,26,27,28,29,30,31}}.
]

---

## Exact (\gamma_{\rm SM}) in this basis

To match the v4.0 template (H_{\rm SM}=H_L\oplus H_R\oplus H_R^c\oplus H_L^c) (so that the “conjugate-right” sits with the (+) grading), set
[
\gamma_{\rm SM}
===============

\mathrm{diag}(\underbrace{+1,\dots,+1}*{8\ \text{(particle L)}},
\underbrace{-1,\dots,-1}*{8\ \text{(particle R)}},
\underbrace{+1,\dots,+1}*{8\ \text{(conj R)}},
\underbrace{-1,\dots,-1}*{8\ \text{(conj L)}}).
]
In indices:

* (+1) on ({0,\dots,7}\cup{24,\dots,31}),
* (-1) on ({8,\dots,15}\cup{16,\dots,23}).

Then an internal operator (D_{\rm int}) is **odd** iff it only connects (+\leftrightarrow-) index blocks.

---

## Exact (J_{\rm SM}) as a permutation (antiunitary (=U\circ K))

In this concrete basis, take
[
J_{\rm SM} = U_{\rm SM}\circ K,
\qquad
U_{\rm SM}e_i = e_{i+16},\ \ U_{\rm SM}e_{i+16}=e_i.
]
So (U_{\rm SM}) is the (32\times 32) permutation matrix swapping particle (\leftrightarrow) conjugate with the same internal label. (If you later need SM-specific phases/signs, you multiply (U_{\rm SM}) by a diagonal unitary; the swap structure stays.)

---

## Lift to (H_{\rm SM}\otimes H_{\rm flav}) (dimension (32N))

Every SM projector becomes
[
\widetilde P = P\otimes 1_N,
]
and every flavor texture becomes
[
\widetilde Y = 1_{32}\otimes Y,\qquad Y\in\mathcal B(H_{\rm flav}),
]
which makes the commutant condition (Y\in\pi(A)') manifest in code.

---

## Literal sparse blueprint for (D_{\rm int})

With the ordering above, define the LR Yukawa map (acting on (H_{\rm SM}\otimes H_{\rm flav})) as
[
M_D
===

\big(P_{u_R}P_{u_L}\big)\otimes Y_u
+
\big(P_{d_R}P_{d_L}\big)\otimes Y_d
+
\big(P_{\nu_R}P_{\nu_L}\big)\otimes Y_\nu
+
\big(P_{e_R}P_{e_L}\big)\otimes Y_e,
]
where the “component projectors” are the diagonals
[
P_{u_L}=P_{{0,1,2}},\quad P_{d_L}=P_{{3,4,5}},\quad
P_{\nu_L}=P_{{6}},\quad P_{e_L}=P_{{7}}.
]

Then a v4.0-compatible internal Dirac template is the (4\times 4) block operator on
[
H_{\rm SM}\otimes H_{\rm flav}
==============================

(H_L\otimes H_{\rm flav})\oplus(H_R\otimes H_{\rm flav})\oplus(H_R^c\otimes H_{\rm flav})\oplus(H_L^c\otimes H_{\rm flav}),
]
[
D_{\rm int}
===========

\begin{pmatrix}
0 & M_D^\dagger & 0 & 0 \
M_D & 0 & M_R & 0 \
0 & M_R^\dagger & 0 & \overline{M_D} \
0 & 0 & \overline{M_D}^\dagger & 0
\end{pmatrix},
]
with (M_R=(P_{\nu_R}\otimes 1)\otimes M_R^{\rm flav}) (Takagi-symmetric on flavor) if you include seesaw, else (M_R=0).

This template is now **fully unambiguous**: every nonzero entry corresponds to a named projector channel and a specific flavor texture.

---

## Minimal code to generate these projectors (NumPy)

```python
import numpy as np

def proj_from_indices(dim, idx):
    P = np.zeros((dim, dim), dtype=float)
    P[idx, idx] = 1.0
    return P

# H_SM dimension
DIM_SM = 32

# Particle indices
uL = [0,1,2]; dL = [3,4,5]; nuL = [6]; eL = [7]
uR = [8,9,10]; dR = [11,12,13]; eR = [14]; nuR = [15]

# Conjugates: +16
uL_c = [i+16 for i in uL]; dL_c = [i+16 for i in dL]; nuL_c = [i+16 for i in nuL]; eL_c = [i+16 for i in eL]
uR_c = [i+16 for i in uR]; dR_c = [i+16 for i in dR]; eR_c = [i+16 for i in eR]; nuR_c = [i+16 for i in nuR]

P_uL  = proj_from_indices(DIM_SM, uL)
P_dL  = proj_from_indices(DIM_SM, dL)
P_nuL = proj_from_indices(DIM_SM, nuL)
P_eL  = proj_from_indices(DIM_SM, eL)

P_uR  = proj_from_indices(DIM_SM, uR)
P_dR  = proj_from_indices(DIM_SM, dR)
P_eR  = proj_from_indices(DIM_SM, eR)
P_nuR = proj_from_indices(DIM_SM, nuR)

# gamma_SM: + on particle L and conjugate R; - on particle R and conjugate L
plus  = list(range(0,8)) + list(range(24,32))
minus = list(range(8,16)) + list(range(16,24))
gamma_SM = proj_from_indices(DIM_SM, plus) - proj_from_indices(DIM_SM, minus)

# U_SM for J_SM = U_SM ∘ K : swap i <-> i+16
U_SM = np.zeros((DIM_SM, DIM_SM), dtype=complex)
for i in range(16):
    U_SM[i, i+16] = 1.0
    U_SM[i+16, i] = 1.0
```

Once you tensor these with (1_N) on flavor, you have a literal sparse blueprint for assembling (D_{\rm int}) and for computing the spectral-action traces.

Not really — your **diagonal projectors (P_S=\sum_{i\in S}E_{ii})** are fine **as long as you keep exactly the basis ordering you stated**. What *does* need updating is **how you use them to build the off-diagonal maps** (your (M_D) blueprint) and, depending on taste, the **4×4 block layout** (because your conjugate ordering is (L^c) then (R^c), not (R^c) then (L^c)).

### 1) Your projectors are consistent with your stated ordering

* Particle indices (0!-!15) and conjugates (16!-!31) (same order) ⇒ all the listed (P_{Q_L},P_{u_R},\dots,P_{\nu_R^c}) are correct.
* Your (\gamma_{\rm SM}) choice (“(+)” on particle (L) and conjugate (R), “(-)” on particle (R) and conjugate (L)) is consistent with that ordering.

### 2) What *must* be updated: your (M_D) formula

This line (and its siblings)
[
(P_{u_R}P_{u_L})\otimes Y_u
]
is **always zero**, because (P_{u_R}P_{u_L}=0) (orthogonal diagonal projectors onto disjoint subspaces).

What you want is a **partial isometry / matrix-unit sum** that maps the (u_L) basis vectors to the corresponding (u_R) basis vectors (and similarly for (d,e,\nu)):

[
S_u:=\sum_{c=r,g,b},|u_R^c\rangle\langle u_L^c|
=E_{8,0}+E_{9,1}+E_{10,2},
]
[
S_d:=E_{11,3}+E_{12,4}+E_{13,5},\quad
S_e:=E_{14,7},\quad
S_\nu:=E_{15,6}.
]

Then the LR Dirac/Yukawa map on (H_{\rm SM}\otimes H_{\rm flav}) is
[
\boxed{
M_D
===

S_u\otimes Y_u
+
S_d\otimes Y_d
+
S_e\otimes Y_e
+
S_\nu\otimes Y_\nu,
}
]
which is automatically **odd** w.r.t. your (\gamma_{\rm SM}) (it connects (+\to-)).

### 3) Majorana block: place it between (\nu_R) and (\nu_R^c)

With your indexing, (\nu_R) is (15) and (\nu_R^c) is (31), so the SM-side matrix unit is (E_{31,15}) (and adjoint (E_{15,31})). Thus
[
\boxed{
M_R = E_{31,15}\otimes M_R^{\rm flav}
}
]
(with (M_R^{\rm flav}) “Takagi-symmetric” in your (\sharp)-sense).

### 4) Your 4×4 block template needs a permutation **or** a block reorder

Your current *basis order* is
[
(H_L)\ \oplus\ (H_R)\ \oplus\ (H_L^c)\ \oplus\ (H_R^c),
]
but the “v4 template” you wrote assumes
[
(H_L)\ \oplus\ (H_R)\ \oplus\ (H_R^c)\ \oplus\ (H_L^c).
]

So either:

* **(Preferred for code)** keep your basis as-is, and write the block matrix in the ((L,R,L^c,R^c)) order; **or**
* introduce a **permutation** that swaps the conjugate blocks (16!-!23 \leftrightarrow 24!-!31) and conjugate (D_{\rm int}) by it to recover the v4 ordering.

Your projectors themselves don’t need to change for the first option.

---

## Minimal NumPy fix: build the intertwiners (S_u,S_d,S_e,S_\nu) (not projector products)

```python
import numpy as np

DIM_SM = 32

def partial_isometry(dim, tgt, src):
    """Sum_k |tgt_k><src_k| as a dim×dim matrix."""
    M = np.zeros((dim, dim), dtype=complex)
    for t, s in zip(tgt, src):
        M[t, s] = 1.0
    return M

# indices (your convention)
uL = [0,1,2]; dL = [3,4,5]; nuL = [6]; eL = [7]
uR = [8,9,10]; dR = [11,12,13]; eR = [14]; nuR = [15]
nuR_c = [31]

S_u  = partial_isometry(DIM_SM, uR,  uL)
S_d  = partial_isometry(DIM_SM, dR,  dL)
S_e  = partial_isometry(DIM_SM, eR,  eL)
S_nu = partial_isometry(DIM_SM, nuR, nuL)

# SM-side Majorana unit: |nu_R^c><nu_R|
S_MR = partial_isometry(DIM_SM, nuR_c, nuR)

# Then on H_SM ⊗ H_flav:
# M_D = kron(S_u, Y_u) + kron(S_d, Y_d) + kron(S_e, Y_e) + kron(S_nu, Y_nu)
# M_R = kron(S_MR, M_R_flav)
```

---

### Bottom line

* **No**, your *diagonal index projectors* don’t need updating.
* **Yes**, update the **construction of (M_D)** (replace (P_{u_R}P_{u_L}) by matrix-unit/partial-isometry maps), and ensure your **4×4 block layout** matches your conjugate ordering (or permute the basis).
