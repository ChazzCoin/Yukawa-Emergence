## Explicit block template for (D_{\rm int}(Y[\mathcal K])) in v4.0

Goal: a **self-adjoint**, **odd** internal operator on (H_{\rm SM}\otimes H_{\rm flav}) whose LR blocks carry the sector Yukawas ({Y_u,Y_d,Y_e,Y_\nu}) (each acting on (H_{\rm flav})), and whose optional Majorana piece implements the seesaw **without breaking evenness**.

### 1) Internal Hilbert space and grading

Take one-generation SM particle space (no generations here; generations live in (H_{\rm flav})):

* Quark doublet: (Q_L \cong \mathbb C^2\otimes \mathbb C^3)
* Lepton doublet: (L_L \cong \mathbb C^2)
* Singlets: (u_R\cong \mathbb C^3,; d_R\cong \mathbb C^3,; e_R\cong \mathbb C,; \nu_R\cong \mathbb C) (include (\nu_R) if you want seesaw)

Define
[
H_L := Q_L\oplus L_L,\qquad
H_R := u_R\oplus d_R\oplus e_R\oplus \nu_R.
]

To place Majorana blocks while keeping oddness, use the **Majorana-compatible charge-conjugate ordering**
[
H_{\rm SM}
==========

H_L;\oplus;H_R;\oplus;H_R^{c};\oplus;H_L^{c},
]
and set the SM grading (on this internal space) to be
[
\gamma_{\rm SM}
===============

\mathrm{diag}(+1,,-1,,+1,,-1)
\quad\text{on}\quad
(H_L,;H_R,;H_R^{c},;H_L^{c}).
]
Then any operator (D_{\rm int}) that only connects (+\leftrightarrow-) blocks automatically satisfies
[
{\gamma_{\rm SM},D_{\rm int}}=0,
]
which is exactly what v4.0 requires (since ({\gamma_{\rm SM}\otimes 1, D_{\rm int}}=0)).

Finally,
[
H_{\rm int}:=H_{\rm SM}\otimes H_{\rm flav}.
]

### 2) Sector projectors on (H_{\rm SM})

Let
[
P_{Q_L},,P_{L_L},,P_{u_R},,P_{d_R},,P_{e_R},,P_{\nu_R}
]
be the orthogonal projectors onto the indicated subspaces of (H_{\rm SM}) (and similarly (P_{Q_L^c},\dots) onto conjugates). These are fixed once you fix the basis ordering.

Inside each SU(2) doublet (\mathbb C^2) with basis ((\uparrow,\downarrow)), define the component projectors
[
p_\uparrow=\begin{pmatrix}1&0\0&0\end{pmatrix},\qquad
p_\downarrow=\begin{pmatrix}0&0\0&1\end{pmatrix}.
]
On (Q_L=\mathbb C^2\otimes \mathbb C^3), use (p_{\uparrow,\downarrow}\otimes 1_3). On (L_L=\mathbb C^2), use (p_{\uparrow,\downarrow}).

### 3) Yukawa textures as commutant operators on (H_{\rm flav})

Let
[
Y_u,;Y_d,;Y_e,;Y_\nu \in \mathcal B(H_{\rm flav})
]
be your Alignment-generated textures (from the kernel pipeline), and impose the v4.0 commutant protection:
[
Y_s \in \pi(A)'\qquad (s\in{u,d,e,\nu}).
]

### 4) The Dirac (LR) map (M_D)

Define (M_D: H_L\otimes H_{\rm flav}\to H_R\otimes H_{\rm flav}) by
[
\boxed{
M_D
===

\big(P_{u_R}(p_\uparrow\otimes 1_3)P_{Q_L}\big)\otimes Y_u
+
\big(P_{d_R}(p_\downarrow\otimes 1_3)P_{Q_L}\big)\otimes Y_d
+
\big(P_{\nu_R}p_\uparrow P_{L_L}\big)\otimes Y_\nu
+
\big(P_{e_R}p_\downarrow P_{L_L}\big)\otimes Y_e.
}
]
This formula is the “unambiguous matching” between:

* SM sector projectors ((P_{(\cdot)})),
* SU(2) component selection ((p_\uparrow,p_\downarrow)),
* flavor textures ((Y_s)).

### 5) Optional Majorana map (M_R) (seesaw)

Let the heavy Majorana scale live in a symmetric operator on flavor space:
[
M_R^{\rm flav}\in \mathcal B(H_{\rm flav}),
\qquad
(M_R^{\rm flav})^T = M_R^{\rm flav}
\quad\text{(Takagi symmetry in a fixed basis)}.
]
Embed it only on the (\nu_R) subspace:
[
\boxed{
M_R
===

\big(P_{\nu_R}\big)\otimes M_R^{\rm flav}
:; H_R\otimes H_{\rm flav}\to H_R^{c}\otimes H_{\rm flav}.
}
]
(If you omit seesaw, set (M_R=0).)

### 6) The full block template for (D_{\rm int}(Y[\mathcal K]))

With the ordering
[
H_{\rm SM}=H_L\oplus H_R\oplus H_R^c\oplus H_L^c,
]
define
[
\boxed{
D_{\rm int}(Y[\mathcal K])
==========================

\begin{pmatrix}
0 & M_D^\dagger & 0 & 0 \
M_D & 0 & M_R & 0 \
0 & M_R^\dagger & 0 & \overline{M_D} \
0 & 0 & \overline{M_D}^\dagger & 0
\end{pmatrix}.
}
]
Here (\overline{M_D}) is the charge-conjugate Dirac map on the conjugate subspaces (in a fixed basis it is the entrywise complex conjugate of (M_D); more invariantly it is (J_{\rm SM} M_D J_{\rm SM}^{-1}) restricted to the conjugate blocks).

This template satisfies, by construction:

* self-adjointness: (D_{\rm int}^\dagger=D_{\rm int}),
* oddness: ({\gamma_{\rm SM}\otimes 1,;D_{\rm int}}=0) (since only (+\leftrightarrow-) blocks are populated),
* sector clarity: every appearance of (Y_u,Y_d,Y_e,Y_\nu) is pinned to a specific SM projector channel,
* seesaw clarity: the only Majorana insertion is (M_R) on the (\nu_R) channel.

### 7) Why this makes the spectral action expansion unambiguous

All bosonic coefficients that depend on Yukawas now reduce to **explicit traces over flavor space** (times known SM multiplicities like color):

Define the Alignment invariants
[
\mathcal Y_2 := \mathrm{Tr}(Y_u^\dagger Y_u)+\mathrm{Tr}(Y_d^\dagger Y_d)+\mathrm{Tr}(Y_e^\dagger Y_e)+\mathrm{Tr}(Y_\nu^\dagger Y_\nu),
]
[
\mathcal Y_4 := \mathrm{Tr}\big((Y_u^\dagger Y_u)^2\big)+\mathrm{Tr}\big((Y_d^\dagger Y_d)^2\big)+\mathrm{Tr}\big((Y_e^\dagger Y_e)^2\big)+\mathrm{Tr}\big((Y_\nu^\dagger Y_\nu)^2\big),
]
and (if present)
[
\mathcal M_2 := \mathrm{Tr}(M_R^{\dagger}M_R),\qquad
\mathcal M_4 := \mathrm{Tr}\big((M_R^{\dagger}M_R)^2\big).
]
Then every Yukawa-dependent term in the heat-kernel coefficients is a polynomial in (\mathcal Y_2,\mathcal Y_4) (and (\mathcal M_2,\mathcal M_4) if seesaw is included), with fixed SM factors (color 3 for quarks, SU(2) doublet structure from the (p_{\uparrow,\downarrow}) insertions). No ambiguity remains because the projectors specify exactly which channels contribute.

