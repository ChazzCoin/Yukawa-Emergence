## Explicit block template for (D_{\rm int}(Y[\mathcal K])) in v5

Goal: a **self-adjoint**, **odd** internal operator on (H_{\rm SM}\otimes H_{\rm flav}) whose LR blocks carry the sector Yukawas ({Y_u,Y_d,Y_e,Y_\nu}) (acting on (H_{\rm flav})), and whose optional Majorana piece implements the seesaw **without violating the v5 evenness seal** (i.e. without any separate (\gamma_{\rm geom}\otimes 1\otimes D_{\rm flav}) term in the *total* Dirac).

---

### 1) Internal Hilbert space and grading (v5-compatible)

Take one-generation SM internal particle content (generations live only in (H_{\rm flav})):

[
Q_L \cong \mathbb C^2\otimes \mathbb C^3,\qquad
L_L \cong \mathbb C^2,
]
[
u_R\cong \mathbb C^3,\quad
d_R\cong \mathbb C^3,\quad
e_R\cong \mathbb C,\quad
\nu_R\cong \mathbb C\ (\text{optional}).
]

Define
[
H_L := Q_L\oplus L_L,\qquad
H_R := u_R\oplus d_R\oplus e_R\oplus \nu_R.
]

Use the **Majorana-compatible ordering**
[
H_{\rm SM}=H_L\oplus H_R\oplus H_R^{c}\oplus H_L^{c},
]
and define the grading on (H_{\rm SM}) by
[
\gamma_{\rm SM}
===============

\mathrm{diag}(+1,,-1,,+1,,-1)
\quad\text{on}\quad
(H_L,\ H_R,\ H_R^{c},\ H_L^{c}).
]
Then (D_{\rm int}) is **odd** iff it connects only the (+) and (-) blocks:
[
\boxed{\ {\gamma_{\rm SM}\otimes 1_{\rm flav},,D_{\rm int}}=0\ }.
]

Finally,
[
H_{\rm int}:=H_{\rm SM}\otimes H_{\rm flav}.
]

---

### 2) Sector projectors on (H_{\rm SM}) (fixed once the basis is fixed)

Let
[
P_{Q_L},,P_{L_L},,P_{u_R},,P_{d_R},,P_{e_R},,P_{\nu_R}
]
be the orthogonal projectors onto the indicated subspaces of (H_{\rm SM}), and similarly
(P_{Q_L^c},\dots,P_{\nu_R^c}) onto their conjugates.

Inside each (SU(2)) doublet (\mathbb C^2) with basis ((\uparrow,\downarrow)), define
[
p_\uparrow=\begin{pmatrix}1&0\0&0\end{pmatrix},\qquad
p_\downarrow=\begin{pmatrix}0&0\0&1\end{pmatrix}.
]
On (Q_L=\mathbb C^2\otimes\mathbb C^3) use (p_{\uparrow,\downarrow}\otimes 1_3), and on (L_L=\mathbb C^2) use (p_{\uparrow,\downarrow}).

---

### 3) Yukawa textures as **commutant-inserted** flavor operators (v5 invariant)

Let
[
Y_u,\ Y_d,\ Y_e,\ Y_\nu\in\mathcal B(H_{\rm flav})
]
be the alignment-produced textures. In v5, the production invariant is **literal commutant insertion**:

[
\boxed{
\text{All flavor/Yukawa/Majorana dependence appearing inside }D_{\rm int}
\text{ lies in }
\pi(A_{\rm SM})'\otimes\mathcal B(H_{\rm flav})
\subset \pi(A)'.
}
]

A sufficient (and typical) implementation is the pure-flavor lift
[
\widetilde Y_s:=1_{H_{\rm SM}}\otimes Y_s\in \pi(A)'.
]

---

### 4) The Dirac LR map (M_D) (use channel maps, not projector products)

Define (M_D:H_L\otimes H_{\rm flav}\to H_R\otimes H_{\rm flav}) as a sum of **channel partial isometries** tensored with flavor textures:
[
\boxed{
M_D
===

V_u\otimes Y_u
+
V_d\otimes Y_d
+
V_\nu\otimes Y_\nu
+
V_e\otimes Y_e,
}
]
where the SM channel maps (V_s\in\mathcal B(H_{\rm SM})) are the canonical injections from the appropriate left submultiplet into the appropriate right singlet:

[
V_u := P_{u_R},(p_\uparrow\otimes 1_3),P_{Q_L},
\qquad
V_d := P_{d_R},(p_\downarrow\otimes 1_3),P_{Q_L},
]
[
V_\nu := P_{\nu_R},p_\uparrow,P_{L_L},
\qquad
V_e := P_{e_R},p_\downarrow,P_{L_L}.
]

> Note: this is the v5-correct way to write it. Expressions like (P_{u_R}P_{u_L}) vanish identically; the *channel map* (V_u) is the right object.

---

### 5) Optional Majorana map (M_R) (seesaw) in a basis-free way

Let the heavy Majorana texture live on flavor space:
[
M_R^{\rm flav}\in\mathcal B(H_{\rm flav}).
]

For basis-free reality, it is cleaner to use the (J_{\rm flav})-transpose (\sharp) (Majorana adjoint). The symmetry condition is
[
\boxed{\ (M_R^{\rm flav})^\sharp = M_R^{\rm flav}\ }
]
(which reduces to (M_R^{T}=M_R) in a basis with (J_{\rm flav}=K)).

Embed it only on the (\nu_R) channel as a map (H_R\otimes H_{\rm flav}\to H_R^c\otimes H_{\rm flav}):
[
\boxed{
M_R := W_R\otimes M_R^{\rm flav},
\qquad
W_R:=P_{\nu_R^c},J_{\rm SM},P_{\nu_R}.
}
]
(Equivalently, in a fixed basis one may take (W_R=|\nu_R^c\rangle\langle \nu_R|). If seesaw is omitted, set (M_R=0).)

This placement keeps (D_{\rm int}) odd with respect to (\gamma_{\rm SM}) and is compatible with the v5 evenness seal because it lives *inside* (D_{\rm int}).

---

### 6) Full block template for (D_{\rm int}) (v5 form)

With the ordering
[
H_{\rm SM}=H_L\oplus H_R\oplus H_R^c\oplus H_L^c,
]
define
[
\boxed{
D_{\rm int}
===========

\begin{pmatrix}
0 & M_D^\dagger & 0 & 0 \
M_D & 0 & M_R & 0 \
0 & M_R^\dagger & 0 & J_{\rm int} M_D J_{\rm int}^{-1} \
0 & 0 & (J_{\rm int} M_D J_{\rm int}^{-1})^\dagger & 0
\end{pmatrix},
}
]
where
[
J_{\rm int}:=J_{\rm SM}\otimes J_{\rm flav}.
]

In a basis where (J_{\rm SM}) is “swap-with-conjugate” and (J_{\rm flav}=K), the conjugate Dirac block is entrywise complex conjugation:
[
J_{\rm int} M_D J_{\rm int}^{-1}=\overline{M_D}.
]

This template satisfies by construction:

* **self-adjointness:** (D_{\rm int}^\dagger=D_{\rm int}),
* **oddness:** ({\gamma_{\rm SM}\otimes 1_{\rm flav},D_{\rm int}}=0),
* **sector clarity:** each (Y_s) appears only in its named projector-defined channel (V_s),
* **seesaw clarity:** the only Majorana insertion is (M_R) on the (\nu_R\leftrightarrow \nu_R^c) channel,
* **v5 evenness seal:** there is **no** separate (\gamma_{\rm geom}\otimes 1\otimes D_{\rm flav}) term in the *total* Dirac.

---

### 7) Why the spectral-action Yukawa dependence is unambiguous (v5 wording)

All Yukawa/Majorana dependence enters through explicit traces over flavor space (times fixed SM multiplicities, e.g. color).

Define the alignment invariants
[
\mathcal Y_2 := \mathrm{Tr}*{\rm flav}(Y_u^\dagger Y_u)+\mathrm{Tr}*{\rm flav}(Y_d^\dagger Y_d)+\mathrm{Tr}*{\rm flav}(Y_e^\dagger Y_e)+\mathrm{Tr}*{\rm flav}(Y_\nu^\dagger Y_\nu),
]
[
\mathcal Y_4 := \mathrm{Tr}*{\rm flav}!\big((Y_u^\dagger Y_u)^2\big)+\mathrm{Tr}*{\rm flav}!\big((Y_d^\dagger Y_d)^2\big)+\mathrm{Tr}*{\rm flav}!\big((Y_e^\dagger Y_e)^2\big)+\mathrm{Tr}*{\rm flav}!\big((Y_\nu^\dagger Y_\nu)^2\big),
]
and, if seesaw is present,
[
\mathcal M_2 := \mathrm{Tr}*{\rm flav}(M_R^{\rm flav,\dagger}M_R^{\rm flav}),\qquad
\mathcal M_4 := \mathrm{Tr}*{\rm flav}!\big((M_R^{\rm flav,\dagger}M_R^{\rm flav})^2\big).
]
Then every Yukawa-dependent coefficient in the heat-kernel/spectral-action expansion is a polynomial in (\mathcal Y_2,\mathcal Y_4) (and (\mathcal M_2,\mathcal M_4) if present) with fixed SM multiplicities determined by the projector-defined channels (V_s). No ambiguity remains because the SM projectors and channel maps specify exactly which internal degrees of freedom contribute.

## v5 sparse template in the 32-index (H_{\rm SM}) basis

### Basis (as you fixed it)

Particle indices (0)–(15):

* (u_L^{r,g,b}): (0,1,2)
* (d_L^{r,g,b}): (3,4,5)
* (\nu_L): (6)
* (e_L): (7)
* (u_R^{r,g,b}): (8,9,10)
* (d_R^{r,g,b}): (11,12,13)
* (e_R): (14)
* (\nu_R): (15)

Conjugates (16)–(31) in the same order: (i\mapsto i+16).

Matrix units: (E_{ij}:=|e_i\rangle\langle e_j|\in M_{32}(\mathbb C)).

You also fixed
[
\gamma_{\rm SM}=\mathrm{diag}(+ \text{ on } {0,\dots,7}\cup{24,\dots,31},\ \ - \text{ on } {8,\dots,15}\cup{16,\dots,23}),
]
so “odd” means: **only connect (+\leftrightarrow-)**.

---

## 1) Correct partial isometries (S_u,S_d,S_e,S_\nu) (particle LR)

These implement the SM channel maps (L\to R) **without any bogus (P_{u_R}P_{u_L}) products**:

[
\boxed{
\begin{aligned}
S_u &:= E_{8,0}+E_{9,1}+E_{10,2},\
S_d &:= E_{11,3}+E_{12,4}+E_{13,5},\
S_e &:= E_{14,7},\
S_\nu &:= E_{15,6}.
\end{aligned}}
]

Each (S_s) maps particle-left ((+)) to particle-right ((-)), hence is (\gamma_{\rm SM})-odd.

---

## 2) Conjugate-side LR partial isometries (S_u^c,S_d^c,S_e^c,S_\nu^c)

These implement (L^c\to R^c) (which in your grading is (-\to+)):

[
\boxed{
\begin{aligned}
S_u^{c} &:= E_{24,16}+E_{25,17}+E_{26,18},\
S_d^{c} &:= E_{27,19}+E_{28,20}+E_{29,21},\
S_e^{c} &:= E_{30,23},\
S_\nu^{c} &:= E_{31,22}.
\end{aligned}}
]

---

## 3) Flavor insertions (v5 commutant channel)

Take sector textures
[
Y_u,Y_d,Y_e,Y_\nu\in\mathcal B(H_{\rm flav})\cong M_N(\mathbb C).
]

On the conjugate side, the v5-correct thing is **antiunitary conjugation by (J_{\rm flav})** (not (\sharp)):
[
Y_s^{c}:=J_{\rm flav},Y_s,J_{\rm flav}^{-1}.
]
In the common coding basis (J_{\rm flav}=K), this is just
[
Y_s^{c}=\overline{Y_s}.
]

(And all these lifted operators live in the commutant channel, e.g. (1_{32}\otimes Y_s\in\pi(A)').)

---

## 4) Optional Majorana channel (seesaw) in this basis

Let (M_R^{\rm flav}\in\mathcal B(H_{\rm flav})) satisfy the **Majorana symmetry**
[
\boxed{(M_R^{\rm flav})^\sharp=M_R^{\rm flav}}
\qquad\text{(so if }J_{\rm flav}=K,\ \ M_R^{T}=M_R\text{)}.
]

Embed only on (\nu_R\leftrightarrow \nu_R^c):
[
\boxed{
S_M:=E_{31,15}
\quad(\nu_R\to \nu_R^c),\qquad
S_M^\dagger=E_{15,31}.
}
]
This is also ( - \leftrightarrow +), hence odd.

---

## 5) The fully sparse v5 (D_{\rm int}) on (H_{\rm SM}\otimes H_{\rm flav})

Define the particle Dirac map and its conjugate copy:
[
M_D := S_u\otimes Y_u+S_d\otimes Y_d+S_e\otimes Y_e+S_\nu\otimes Y_\nu,
]
[
M_D^{c} := S_u^{c}\otimes Y_u^{c}+S_d^{c}\otimes Y_d^{c}+S_e^{c}\otimes Y_e^{c}+S_\nu^{c}\otimes Y_\nu^{c}.
]

Then the v5 internal operator is simply
[
\boxed{
D_{\rm int}
===========

\Big(M_D+M_D^\dagger\Big)
+
\Big(M_D^{c}+(M_D^{c})^\dagger\Big)
+
\Big(S_M\otimes M_R^{\rm flav}+S_M^\dagger\otimes (M_R^{\rm flav})^\dagger\Big).
}
]

### Why this is v5-correct

* **Self-adjoint:** built as sum of (X+X^\dagger).
* **Odd:** every term connects (+) to (-) under your (\gamma_{\rm SM}).
* **Seesaw without breaking evenness:** Majorana lives *inside* (D_{\rm int}); you still take total
  [
  D=D_{\rm geom}\otimes 1\otimes 1+\gamma_{\rm geom}\otimes D_{\rm int},
  ]
  and **you do not add** (\gamma_{\rm geom}\otimes 1\otimes D_{\rm flav}) (v5 evenness seal).

---

# Minimal NumPy “code-is-the-proof” builder

```python
import numpy as np

def E(i, j, dim=32, dtype=complex):
    M = np.zeros((dim, dim), dtype=dtype)
    M[i, j] = 1.0
    return M

def kron(A, B):
    return np.kron(A, B)

def build_S_maps():
    # Particle LR maps
    S_u  = E(8,0)  + E(9,1)  + E(10,2)
    S_d  = E(11,3) + E(12,4) + E(13,5)
    S_e  = E(14,7)
    S_nu = E(15,6)

    # Conjugate LR maps (L^c -> R^c)
    S_u_c  = E(24,16) + E(25,17) + E(26,18)
    S_d_c  = E(27,19) + E(28,20) + E(29,21)
    S_e_c  = E(30,23)
    S_nu_c = E(31,22)

    # Majorana channel: nu_R -> nu_R^c
    S_M = E(31,15)

    return (S_u, S_d, S_e, S_nu, S_u_c, S_d_c, S_e_c, S_nu_c, S_M)

def J_flav_conj(Y, U_flav=None):
    """
    Implements Y^c := J_flav Y J_flav^{-1}.
    If J_flav=K, then Y^c = conjugate(Y).
    If J_flav = U_flav ∘ K, then Y^c = U_flav * conjugate(Y) * U_flav^{-1}.
    """
    if U_flav is None:
        return np.conjugate(Y)
    return U_flav @ np.conjugate(Y) @ np.linalg.inv(U_flav)

def build_D_int(Yu, Yd, Ye, Ynu, MR=None, U_flav=None):
    """
    Returns D_int on H_SM ⊗ H_flav, i.e. a (32N x 32N) matrix.
    MR is M_R^{flav}. If MR is None, seesaw is omitted.
    """
    N = Yu.shape[0]
    I_SM = np.eye(32, dtype=complex)

    (S_u, S_d, S_e, S_nu,
     S_u_c, S_d_c, S_e_c, S_nu_c,
     S_M) = build_S_maps()

    # Conjugate-side flavor insertions
    Yu_c  = J_flav_conj(Yu,  U_flav=U_flav)
    Yd_c  = J_flav_conj(Yd,  U_flav=U_flav)
    Ye_c  = J_flav_conj(Ye,  U_flav=U_flav)
    Ynu_c = J_flav_conj(Ynu, U_flav=U_flav)

    # Build M_D and M_D^c
    M_D   = kron(S_u,  Yu)  + kron(S_d,  Yd)  + kron(S_e,  Ye)  + kron(S_nu,  Ynu)
    M_D_c = kron(S_u_c,Yu_c)+ kron(S_d_c,Yd_c)+ kron(S_e_c,Ye_c)+ kron(S_nu_c,Ynu_c)

    D_int = (M_D + M_D.conj().T) + (M_D_c + M_D_c.conj().T)

    if MR is not None:
        # Majorana insertion nu_R <-> nu_R^c
        D_int += kron(S_M, MR) + kron(S_M.conj().T, MR.conj().T)

    return D_int
```

