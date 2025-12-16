
# Alignment Spectral Triple v5.0r

## Robust, production-ready manuscript

## 0. Standing hypotheses and notation

Fix:

**(H0) Geometric triple.** A real, even spectral triple
[
(A_{\rm geom},H_{\rm geom},D_{\rm geom},J_{\rm geom},\gamma_{\rm geom})
]
with (\gamma_{\rm geom}^2=1), (\gamma_{\rm geom}^\ast=\gamma_{\rm geom}), ([\gamma_{\rm geom},\pi_{\rm geom}(A_{\rm geom})]=0), and ({\gamma_{\rm geom},D_{\rm geom}}=0).

**(H1) One-generation SM internal triple.** A real, even finite triple
[
(A_{\rm SM},H_{\rm SM},D_{\rm SM}(\cdot),J_{\rm SM},\gamma_{\rm SM})
]
where Yukawas are inserted in the standard left–right off-diagonal blocks of the internal Dirac operator.

**(H2) Flavor space.** A finite-dimensional Hilbert space (H_{\rm flav}\cong \mathbb C^N) with antiunitary real structure
[
J_{\rm flav}=U_{\rm flav}\circ K,
]
where (K) is complex conjugation in some (any) basis and (U_{\rm flav}) is unitary.

**(H3) Flavor clock.** A fixed self-adjoint operator (the “clock/Laplacian”)
[
L_N=L_N^\ast\in\mathcal B(H_{\rm flav}),
]
specified once per run (e.g. a graph Laplacian on (C_N) or (P_N)). Spectral projections are taken in the Borel functional calculus.

**Commutators and (anti)commutators.**
[
[X,Y]:=XY-YX,\qquad {X,Y}:=XY+YX.
]

**Norm conventions.** We write (|X|) for the operator norm of (X\in\mathcal B(H)), and (|z|) for the absolute value of a scalar (z\in\mathbb C).

**Majorana adjoint (basis-free transpose).** Define the (J_{\rm flav})-transpose (Majorana adjoint)
[
X^{\sharp}:=J_{\rm flav},X^\ast,J_{\rm flav}^{-1}.
]
In any basis where (J_{\rm flav}=K), one has (X^\sharp=X^T) (ordinary matrix transpose). When (X) acts only on (H_{\rm flav}), (X^\sharp) coincides with the restriction of the full conjugation (X\mapsto JX^\ast J^{-1}) (with (J=J_{\rm geom}\otimes J_{\rm SM}\otimes J_{\rm flav})) to the flavor factor.

**Trace convention (finite-dimensional).** We write (\mathrm{tr}*{\rm flav}) for the usual matrix trace on (\mathcal B(H*{\rm flav})).

---

## 1. The v5.0r spectral datum

### Definition 1.1 (Aligned product algebra, Hilbert space, and representation)

Define
[
A:=A_{\rm geom}\otimes A_{\rm SM},\qquad
H:=H_{\rm geom}\otimes H_{\rm SM}\otimes H_{\rm flav},
]
and represent (A) on (H) by
[
\pi(a_{\rm geom}\otimes a_{\rm SM})
:=\pi_{\rm geom}(a_{\rm geom})\otimes \pi_{\rm SM}(a_{\rm SM})\otimes 1_{\rm flav}.
]

### Definition 1.2 (Real structure and grading)

Define
[
J:=J_{\rm geom}\otimes J_{\rm SM}\otimes J_{\rm flav},\qquad
\Gamma:=\gamma_{\rm geom}\otimes\gamma_{\rm SM}\otimes 1_{\rm flav}.
]

### Definition 1.3 (Flavor multiplicity principle)

Flavor is pure multiplicity: the represented algebra acts trivially on (H_{\rm flav}), i.e.
[
\pi(A)\subseteq \mathcal B(H_{\rm geom}\otimes H_{\rm SM})\otimes 1_{\rm flav}.
]

### Definition 1.4 (v5.0r Dirac operator; no pure flavor term)

Let (D_{\rm int}(Y[\mathcal K],M_R,\dots)) be an internal Dirac operator on (H_{\rm SM}\otimes H_{\rm flav}) (built from LR Yukawa blocks and optionally Majorana blocks) satisfying the **oddness requirement**
[
\boxed{\ {\gamma_{\rm SM}\otimes 1_{\rm flav},D_{\rm int}}=0\ }.
]
Define the total Dirac operator by the **even-product form**
[
\boxed{\ D:=D_{\rm geom}\otimes 1\otimes 1;+;\gamma_{\rm geom}\otimes D_{\rm int}.\ }
]
**Production constraint (evenness seal).** There is no additional summand of the form (\gamma_{\rm geom}\otimes 1\otimes D_{\rm flav}).

### Definition 1.5 (Alignment Spectral Triple v5.0r)

The tuple
[
(A,H,D,J,\Gamma)
]
with Definitions 1.1–1.4 is the **Alignment Spectral Triple v5.0r**.

---

## 2. Commutant protection and exclusion of flavor gauge fields

### Definition 2.1′ (Production commutant insertion invariant)

All operators that carry Yukawa/flavor/Majorana dependence inside (D_{\rm int}) lie in
[
\pi(A_{\rm SM})'\otimes \mathcal B(H_{\rm flav})\ \subset\ \pi(A)'.
]
In particular, the sufficient pattern (\widetilde Y=1_{H_{\rm SM}}\otimes Y) is allowed, but the invariant is literal membership in (\pi(A)').

### Lemma 2.2 (Multiplicity implies commutant membership for pure flavor operators)

Under Def. 1.3, every operator of the form
[
1_{H_{\rm geom}}\otimes 1_{H_{\rm SM}}\otimes T
]
lies in (\pi(A)'). In particular, every (\widetilde Y:=1_{H_{\rm SM}}\otimes Y) lies in (\pi(A)').

### Definition 2.3 (One-forms and inner fluctuations)

Define the space of one-forms
[
\boxed{
\Omega_D^1(A):=\left{\sum_i \pi(a_i)[D,\pi(b_i)] : a_i,b_i\in A\right}\subset\mathcal B(H),
}
]
and for (A\in\Omega_D^1(A)) define the fluctuated operator
[
D_A:=D+A+JAJ^{-1}.
]

### Proposition 2.4 (Inner fluctuations do not generate flavor gauge fields)

Assume Def. 1.3 and the production invariant Def. 2.1′. Then inner fluctuations cannot produce gauge bosons with nontrivial action on (H_{\rm flav}).

**Reason (referee-tight).** Since (\pi(\cdot)) acts as (1_{\rm flav}), every commutator ([D,\pi(b)]) is of the form ((\cdots)\otimes 1_{\rm flav}), hence (\Omega_D^1(A)\subset\mathcal B(H_{\rm geom}\otimes H_{\rm SM})\otimes 1_{\rm flav}). Therefore any fluctuation (A) acts trivially on flavor, and so does (JAJ^{-1}).

---

## 3. Canonical spectral anchor for “360”

### Definition 3.1 (Spectral 360-subspace and projector)

Fix a Borel set (\boxed{\ \Omega^{(N)}*{360}\subset\sigma(L_N)\ }). Define the **360 projector**
[
\boxed{\ \Pi*{360}^{(N)}:=\chi_{\Omega^{(N)}_{360}}(L_N).\ }
]
This is basis-independent by functional calculus.

### Lemma 3.2 (Anchor commutation identities)

[
[\Pi_{360}^{(N)},L_N]=0.
]
Moreover, under Def. 1.3 (multiplicity),
[
[\Pi_{360}^{(N)},\pi(a)]=0\qquad \forall a\in A,
]
since (\Pi_{360}^{(N)}) acts only on (H_{\rm flav}).

---

## 4. Geometric texture generation via functional calculus

### Definition 4.1 (Kernel from the flavor clock)

For (\alpha>0), define the heat-kernel operator
[
K_\alpha:=e^{-\alpha L_N}\in\mathcal B(H_{\rm flav}).
]
(Any fixed functional calculus choice (K=f(L_N)) is admissible; (\alpha) is the hierarchy-control knob.)

### Definition 4.2 (v5.0r texture map)

Given a seed (Y_0\in\mathcal B(H_{\rm flav})), define the produced texture
[
\boxed{
Y:=\mathcal T_{\alpha,360}(Y_0):=\Pi_{360}^{(N)},K_\alpha,Y_0,K_\alpha,\Pi_{360}^{(N)}.
}
]
Insert (Y) into the internal Dirac operator through the commutant channel required by Def. 2.1′ (e.g. the sufficient lift (\widetilde Y:=1_{H_{\rm SM}}\otimes Y)).

### Lemma 4.3 (Covariant basis-independence)

Under a unitary change of flavor basis (U) with covariant transport
[
L_N\mapsto U L_N U^\ast,\qquad Y_0\mapsto U Y_0 U^\ast,
]
one has
[
\boxed{
\mathcal T_{\alpha,360}(U Y_0 U^\ast)=U,\mathcal T_{\alpha,360}(Y_0),U^\ast.
}
]

---

## 5. Canonical heavy–light split (9→3) without basis leakage

From here set (N=9).

### Definition 5.1 (Light projector: spectral / anchored / refined)

There are three admissible **basis-free** ways to define a rank-(3) light projector (P_L) on (H_{\rm flav}\cong\mathbb C^9):

**(A) Purely spectral split (when feasible).** Choose a Borel set (\Lambda\subset\sigma(L_9)) and set
[
P_L:=\chi_{\Lambda}(L_9),\qquad P_H:=1-P_L,
]
with (\mathrm{rank}(P_L)=3).

**(B) Anchored split (always feasible).** Fix an isometry (i:\mathbb C^3\hookrightarrow\mathbb C^9) as model data and set
[
P_L:=ii^\ast,\qquad P_H:=1-ii^\ast.
]

**(C) Refined spectral split (degeneracy-safe).** Choose an additional self-adjoint anchor (Q=Q^\ast\in\mathcal B(H_{\rm flav})) such that
[
[L_9,Q]=0,
]
and define (P_L) as a **joint spectral projector** in (W^*(L_9,Q)), i.e.
[
P_L\in W^*(L_9,Q),\qquad \mathrm{rank}(P_L)=3.
]
Since ([L_9,Q]=0) in finite dimension, (L_9) and (Q) are simultaneously diagonalizable and projections in (W^*(L_9,Q)) are sums of joint spectral projectors (P_{\lambda,\mu}=P^{(L_9)}*\lambda P^{(Q)}*\mu).

### Lemma 5.1a (Feasibility of rank-3 purely spectral splitting)

Let (L=L^\ast) with spectral decomposition (L=\sum_{\lambda}\lambda P_\lambda) and eigenspaces (E_\lambda=\mathrm{Ran}(P_\lambda)). Then
[
\mathrm{rank}\big(\chi_{\Lambda}(L)\big)=\sum_{\lambda\in\Lambda}\dim(E_\lambda).
]
Hence Option (A) is feasible iff there exists (\Lambda\subset\sigma(L_9)) such that
[
\sum_{\lambda\in\Lambda}\dim(E_\lambda)=3.
]
If this fails due to degeneracy pattern, a rank-(3) split requires either Option (B) (anchored (i)) or Option (C) (commuting refinement).

### Lemma 5.2 (Projector algebra)

For any of the options (A)–(C),
[
P_L^2=P_L=P_L^\ast,\quad P_H^2=P_H=P_H^\ast,\quad P_LP_H=0,\quad P_L+P_H=1.
]
Moreover, in the spectral options (A) and (C),
[
[P_L,L_9]=0,\qquad [P_H,L_9]=0,
]
and in option (C) additionally ([P_L,Q]=[P_H,Q]=0).

### Definition 5.3 (Block extraction by projectors)

For (X\in\mathcal B(H_{\rm flav})), define
[
X_{LL}=P_LXP_L,\quad X_{LH}=P_LXP_H,\quad X_{HL}=P_HXP_L,\quad X_{HH}=P_HXP_H.
]

### Definition 5.4 (Seesaw effective light operator; basis-free)

Let (Y_\nu\in\mathcal B(H_{\rm flav})) and (M_R\in\mathcal B(H_{\rm flav})) such that (M_{R,HH}) is invertible on (H_H). Define
[
m_{\nu,\mathrm{eff}}
:=
-v^2,Y_{\nu,LH},(M_{R,HH})^{-1},Y_{\nu,HL}^{\sharp}
\quad\text{on }H_L.
]
This is a projector-defined functional of the same operators inserted into (D_{\rm int}).

---

## 6. Evenness and the “no pure flavor term” rule

### Theorem 6.1 (Evenness of v5.0r)

Assume ({\gamma_{\rm SM}\otimes 1_{\rm flav},D_{\rm int}}=0). Then ((A,H,D,J,\Gamma)) is even:
[
\Gamma^\ast=\Gamma,\quad \Gamma^2=1,\quad [\Gamma,\pi(A)]=0,\quad {\Gamma,D}=0.
]

### Proposition 6.2 (Failure mode: adding a pure flavor term)

If one adds a separate summand (\gamma_{\rm geom}\otimes 1\otimes D_{\rm flav}) with generic (D_{\rm flav}), then ({\Gamma,D}=0) generally fails. Hence the v5.0r constraint in Def. 1.4 is structurally necessary.

---

## 7. Order-zero and first-order conditions (stability under the 360 move)

### Theorem 7.1 (Order-zero unchanged)

Order-zero depends only on (\pi) and (J):
[
[\pi(a),J\pi(b)J^{-1}]=0\qquad \forall a,b\in A.
]
Since v5.0r modifies neither (\pi) nor (J), and (\Pi_{360}^{(N)}) enters only through texture generation (not through (\pi)), the order-zero condition is unchanged.

### Theorem 7.2 (First-order stable under commutant insertion)

First-order is
[
[[D,\pi(a)],J\pi(b)J^{-1}]=0\qquad \forall a,b\in A.
]
In v5.0r, dependence on the 360 mechanism enters only through (D_{\rm int}) via operators constrained by Def. 2.1′, while flavor remains multiplicity (Def. 1.3). Therefore no new first-order obstruction is introduced by the 360 anchoring mechanism itself.

---

## 8. Production gates (operator norm), including degeneracy-robust (W^*(L)) certification

Fix a finite generator set (\mathcal G_A\subset A) for gate-testing. All gates are expressed in operator norm (|\cdot|).

### Gate 8.1 (Order-zero)

[
\max_{a,b\in\mathcal G_A}\ |[\pi(a),J\pi(b)J^{-1}]|\le \varepsilon_0.
]

### Gate 8.2 (First-order)

[
\max_{a,b\in\mathcal G_A}\ \big|\big[[D,\pi(a)],J\pi(b)J^{-1}\big]\big|\le \varepsilon_1.
]

### Gate 8.3 (Projector sanity)

For each projector (P\in{\Pi_{360}^{(N)},P_L,P_H}),
[
|P^2-P|\le \varepsilon_P,\qquad |P^\ast-P|\le \varepsilon_P.
]

### Lemma 8.4 (Finite-dimensional characterization of ({L}') and (W^*(L)))

Let (L=L^\ast\in\mathcal B(\mathbb C^N)) with spectral decomposition
[
L=\sum_{\lambda\in\sigma(L)} \lambda P_\lambda,
\qquad E_\lambda:=\mathrm{Ran}(P_\lambda).
]

1. **Commutant.**
   [
   \boxed{
   {L}'
   =
   {X\in\mathcal B(\mathbb C^N):[X,L]=0}
   =
   \bigoplus_{\lambda\in\sigma(L)} \mathcal B(E_\lambda).
   }
   ]

2. **Von Neumann algebra generated by (L).**
   [
   \boxed{
   W^*(L)={L}''=({L}')'
   =
   \left{\sum_{\lambda\in\sigma(L)} x_\lambda P_\lambda:\ x_\lambda\in\mathbb C\right}.
   }
   ]
   Equivalently: (X\in W^*(L)) iff (X|*{E*\lambda}=x_\lambda,\mathrm{id}*{E*\lambda}) is scalar on each eigenspace. In finite dimension,
   [
   W^*(L)={f(L): f \text{ Borel}}.
   ]

### Lemma 8.5 (Projections in (W^*(L)) are sums of spectral projectors)

With (L) as above, a projection (\Pi=\Pi^\ast=\Pi^2) satisfies (\Pi\in W^*(L)) iff
[
\boxed{
\Pi=\sum_{\lambda\in S} P_\lambda
\quad\text{for some subset }S\subset\sigma(L).
}
]
Equivalently: (\Pi) selects whole eigenspaces (whole spectral blocks), never a proper subspace of a degenerate block.

### Gate 8.3a (Degeneracy-robust (W^*(L)) membership gate)

Let (P_\lambda) be the spectral projectors of (L_N). For a candidate projector (\Pi) (e.g. (\Pi_{360}^{(N)}), and (P_L) when defined spectrally), define
[
c_\lambda := \frac{\mathrm{tr}*{\rm flav}(P*\lambda \Pi)}{\mathrm{tr}*{\rm flav}(P*\lambda)}\in[0,1],
\qquad
\Delta_\lambda := \big|P_\lambda \Pi P_\lambda - c_\lambda P_\lambda\big|.
]
Require
[
\boxed{
\max_{\lambda}\Delta_\lambda \le \varepsilon_{\rm block}
\quad\text{and}\quad
\max_{\lambda}\min{c_\lambda,1-c_\lambda}\le \varepsilon_{\rm bin}.
}
]
This is an approximate operational certificate that (\Pi) is scalar (hence in (W^*(L_N)), not merely in ({L_N}')) and binary on each spectral block, so (\Pi\approx \sum_{\lambda\in S}P_\lambda) by Lemma 8.5.

**Numerical implementation note (clusters).** If eigenvalues are resolved only up to clusters (\mathcal C), replace (P_\lambda) by cluster projectors (P_{\mathcal C}=\chi_{I_{\mathcal C}}(L_N)) and apply the same test over (\mathcal C).

### Gate 8.3b (Anchor commutation telemetry)

These are exact when projectors are built by functional calculus (useful as sanity checks):
[
|[\Pi_{360}^{(N)},L_N]|=0,\qquad
\max_{a\in\mathcal G_A}|[\Pi_{360}^{(N)},\pi(a)]|=0.
]
(They do not, by themselves, rule out hidden subprojections in degenerate blocks; Gate 8.3a does.)

### Gate 8.4 (Heavy–light split gates)

If (P_L) is defined spectrally (Option A or C), apply Gate 8.3a to (P_L) and require
[
|[P_L,L_9]|=0,
]
and in Option (C) also
[
|[P_L,Q]|=0.
]

If (P_L) is anchored (Option B), require the anchored identity (up to tolerance)
[
|P_L-ii^\ast|\le \varepsilon_{\rm iso},
]
along with the projector sanity gate.

If seesaw is used, require numerical stability of the heavy inverse:
[
|(M_{R,HH})^{-1}|\ \text{finite and stable (monitor condition number as telemetry).}
]

---

## 9. Closure statement (v5.0r meaning of “360 + 9→3”)

**v5.0r does not add a “360 term” to (D).** It fixes a canonical internal spectral object (L_N), defines (\Pi_{360}^{(N)}) (and (P_L) for (N=9)) as spectral/anchored projectors, and uses them only to generate commutant-protected Yukawa/Majorana textures inserted into the odd LR (and optional Majorana) channels of (D_{\rm int}). Thus “360 + 9→3” is basis-free end-to-end and auditable by gates.

A compact operator-algebra slogan is:
[
\boxed{
\Pi_{360}^{(N)}\in W^*(L_N),
\qquad
P_L\in W^*(L_9)\ \text{(Options A/C)}\ \ \text{or}\ \ P_L=ii^\ast\ \text{(Option B).}
}
]
