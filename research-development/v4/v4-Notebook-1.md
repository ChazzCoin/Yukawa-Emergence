## Essential upgrades before we call v4.0 “production-ready”

### 1) Evenness currently fails if (D) contains a pure flavor term

With your present choice
[
\Gamma=\gamma_{\rm geom}\otimes\gamma_{\rm SM}\otimes 1_{\rm flav},
\qquad
D_{\rm flav}^{\rm tot}=\gamma_{\rm geom}\otimes 1\otimes D_{\rm flav},
]
one has ([,\Gamma, D_{\rm flav}^{\rm tot},]=0), hence ({\Gamma,D}\neq 0).
So if we insist the **full triple is even**, (D_{\rm flav}) must not appear as an “even commuting add-on” term.

Production fix (canonical):

* **Flavor lives inside the LR Yukawa blocks** (odd part) as operators on (H_{\rm flav}) that lie in the commutant of (\pi(A)).
* No separate (\gamma_{\rm geom}\otimes 1\otimes D_{\rm flav}) term.

This keeps both: (i) full evenness, (ii) nontrivial mixing.

### 2) Make explicit the commutant condition that protects fluctuations

To keep inner fluctuations from generating “generation-gauge” bosons, state (once) that all flavor/Yukawa matrices satisfy
[
Y \in \pi(A)'\quad\text{(commutant)},
]
i.e. they commute with the represented algebra (they only weight the usual Higgs fluctuations).

### 3) Minor but important notation repairs

* Replace all (|\cdot|) by (|\cdot|) (operator norm) in gates.
* Fix typos: (J_{\rm flav}D_{\rm flav}=\varepsilon'*{\rm flav}D*{\rm flav}J_{\rm flav}) (no stray (*)).
* Remove accidental “==” and doubled commas.
* Use ({X,Y}=XY+YX) and ([X,Y]=XY-YX) consistently.

What follows is the **clean v4.0 text** with those upgrades, still in your Section 1–10 theorem-chain style.

## Section 1 — Definitions

### Definition 1.1 (Alignment v4.0 spectral datum)

Assume:

* A real, even spectral triple ((A_{\rm geom},H_{\rm geom},D_{\rm geom},J_{\rm geom},\gamma_{\rm geom})).
* A real, even **finite** SM triple ((A_{\rm SM},H_{\rm SM},D_{\rm SM}(Y[\mathcal K]),J_{\rm SM},\gamma_{\rm SM})) (one-generation internal geometry) where Yukawa textures are **inserted** into the LR blocks of (D_{\rm SM}).
* A finite-dimensional flavor multiplicity space (H_{\rm flav}\cong\mathbb C^N) with antiunitary (J_{\rm flav}=U_{\rm flav}\circ K).
* Yukawa/flavor operators (the matrices that carry mixing) act on (H_{\rm flav}) and lie in the commutant of the represented algebra:
  [
  Y[\mathcal K]\in \pi(A)'\quad\text{(as operators on }H_{\rm SM}\otimes H_{\rm flav}\text{)}.
  ]

Define
[
A:=A_{\rm geom}\otimes A_{\rm SM},\qquad
H:=H_{\rm geom}\otimes H_{\rm SM}\otimes H_{\rm flav},
]
[
\pi(a_{\rm geom}\otimes a_{\rm SM})
:=\pi_{\rm geom}(a_{\rm geom})\otimes \pi_{\rm SM}(a_{\rm SM})\otimes 1_{\rm flav},
]
and take the Dirac operator in the **even product form**
[
D:=D_{\rm geom}\otimes 1\otimes 1
+\gamma_{\rm geom}\otimes D_{\rm int}(Y[\mathcal K]),
]
where (D_{\rm int}(Y[\mathcal K])) acts on (H_{\rm SM}\otimes H_{\rm flav}) and is **odd** for (\gamma_{\rm SM}\otimes 1):
[
{\gamma_{\rm SM}\otimes 1,,D_{\rm int}(Y[\mathcal K])}=0.
]
Finally,
[
J:=J_{\rm geom}\otimes J_{\rm SM}\otimes J_{\rm flav},\qquad
\Gamma:=\gamma_{\rm geom}\otimes\gamma_{\rm SM}\otimes 1_{\rm flav}.
]
The tuple ((A,H,D,J,\Gamma)) is the **Alignment Spectral Triple v4.0**.

### Definition 1.2 (Flavor multiplicity principle)

Flavor is a pure multiplicity factor: the represented algebra acts trivially on (H_{\rm flav}):
[
\pi(A)\subseteq\mathcal B(H_{\rm geom}\otimes H_{\rm SM})\otimes 1_{\rm flav}.
]
Equivalently, there is **no** additional commutative “site algebra” acting diagonally on (H_{\rm flav}).

### Definition 1.3 (Reality-compatibility for flavor Yukawas)

A flavor operator (Y\in\mathcal B(H_{\rm flav})) is **(J_{\rm flav})-real** with sign (\varepsilon'*{\rm flav}\in{\pm1}) if
[
J*{\rm flav},Y = \varepsilon'*{\rm flav},Y,J*{\rm flav}.
]
In the common basis choice (J_{\rm flav}=K), this is (\overline{Y}=\varepsilon'_{\rm flav}Y).

## Section 2 — Lemma: finiteness and completeness

### Lemma 2.1 (Hilbert space structure and bounded representation)

If (H_{\rm SM}) and (H_{\rm flav}) are finite-dimensional, then (H) is a Hilbert space with the tensor-product inner product, and
(\pi:A\to\mathcal B(H)) is a bounded (*)-representation provided (\pi_{\rm geom}) and (\pi_{\rm SM}) are bounded.

## Section 3 — Theorem: self-adjointness and compact resolvent

Let
[
D_0:=D_{\rm geom}\otimes 1\otimes 1,\qquad
B:=\gamma_{\rm geom}\otimes D_{\rm int}(Y[\mathcal K]).
]
If (D_{\rm int}) is finite (or, more generally, bounded relative to (D_0) with relative bound (<1)), then (B) is (D_0)-bounded.

### Theorem 3.1 (Essential self-adjointness)

If (D_{\rm geom}) is essentially self-adjoint on (\mathrm{Dom}(D_{\rm geom})), then (D=D_0+B) is essentially self-adjoint on
[
\mathrm{Dom}(D):=\mathrm{Dom}(D_{\rm geom})\otimes H_{\rm SM}\otimes H_{\rm flav}.
]

### Theorem 3.2 (Compact resolvent)

If (D_{\rm geom}) has compact resolvent, then (D) has compact resolvent.

## Section 4 — Theorem: bounded commutators

### Theorem 4.1 (Bounded commutators)

For all simple tensors (a=a_{\rm geom}\otimes a_{\rm SM}\in A),
[
[D,\pi(a)]\in\mathcal B(H),
]
provided ([D_{\rm geom},\pi_{\rm geom}(a_{\rm geom})]) is bounded and ([D_{\rm int},\pi_{\rm SM}(a_{\rm SM})\otimes 1]) is bounded.

Explicitly,
[
[D,\pi(a)]
==========

[D_{\rm geom},\pi_{\rm geom}(a_{\rm geom})]\otimes \pi_{\rm SM}(a_{\rm SM})\otimes 1
+\gamma_{\rm geom}\pi_{\rm geom}(a_{\rm geom})\otimes [D_{\rm int},\pi_{\rm SM}(a_{\rm SM})\otimes 1].
]

## Section 5 — Theorem: grading (evenness)

### Theorem 5.1 (Evenness relations)

Assume the geometric and SM factors are even and ([\gamma_{\rm geom},\pi_{\rm geom}(A_{\rm geom})]=[\gamma_{\rm SM},\pi_{\rm SM}(A_{\rm SM})]=0).
If also ({\gamma_{\rm SM}\otimes 1, D_{\rm int}}=0), then with
(\Gamma=\gamma_{\rm geom}\otimes\gamma_{\rm SM}\otimes 1),
[
\Gamma^\ast=\Gamma,\qquad \Gamma^2=1,\qquad {\Gamma,D}=0,\qquad [\Gamma,\pi(A)]=0.
]

## Section 6 — Theorem: reality and KO-signs

### Theorem 6.1 (Reality / opposite action)

(J) is antiunitary on (H) and defines the opposite representation
[
\pi^\circ(a):=J\pi(a)J^{-1}.
]
Implementation rule (matrix form in a fixed basis):
[
JXJ^{-1}=U_J,\overline{X},U_J^\dagger,
\qquad J=U_J\circ K.
]

### Theorem 6.2 (KO-sign relations inherited from factors)

Assume each factor satisfies
[
J_i^2=\varepsilon_i,\qquad J_iD_i=\varepsilon_i' D_iJ_i,\qquad J_i\gamma_i=\varepsilon_i''\gamma_iJ_i,
]
(with (\gamma_{\rm flav}:=1)). Then the product triple satisfies
[
J^2=\varepsilon,\qquad JD=\varepsilon' DJ,\qquad J\Gamma=\varepsilon''\Gamma J,
]
with ((\varepsilon,\varepsilon',\varepsilon'')) determined by the graded tensor product convention, equivalently by total KO-dimension modulo (8).

## Section 7 — Theorem: order-zero condition

### Theorem 7.1 (Order-zero stability under flavor multiplicity)

If order-zero holds on the geometric and SM factors, then for all (a,b\in A),
[
[\pi(a),J\pi(b)J^{-1}]=0.
]
Reason: (\pi(\cdot)) acts as (1_{\rm flav}) on flavor, so no new order-zero obstruction can arise from the flavor factor.

## Section 8 — Theorem: first-order condition

### Theorem 8.1 (First-order stability)

If first-order holds on the geometric and SM factors (with the same graded-product conventions), then for all (a,b\in A),
[
[[D,\pi(a)],J\pi(b)J^{-1}]=0.
]
The flavor multiplicity does not introduce new obstructions because (\pi(A)) is trivial on (H_{\rm flav}) and the flavor/Yukawa operators are chosen in (\pi(A)').

## Section 9 — Inner fluctuations and the unified operator

### Definition 9.1 (One-forms and fluctuated Dirac)

[
\Omega^1_D(A):=\left{\sum_i \pi(a_i)[D,\pi(b_i)] : a_i,b_i\in A\right}\subset\mathcal B(H).
]
Given (A\in\Omega^1_D(A)), define its Hermitian part (A\leftarrow \tfrac12(A+A^\dagger)) and
[
D_A:=D + A + JAJ^{-1}.
]

### Proposition 9.2 (No flavor gauge fields)

Under Definitions 1.2 and the commutant condition (Y[\mathcal K]\in\pi(A)'), the inner fluctuations do not create new gauge bosons acting on flavor: the only dynamical fields generated are those already associated to the (A_{\rm geom}) and (A_{\rm SM}) factors (including Higgs in the SM sense).

## Section 10 — Production commutator gates

Fix a finite generating/spanning set (\mathcal G_A\subset A). On finite truncations, use the operator norm (|\cdot|).

### Gate 10.1 (Order-zero gate)

[
\max_{a,b\in\mathcal G_A}\ \big|[\pi(a),J\pi(b)J^{-1}]\big|\le \varepsilon_0.
]

### Gate 10.2 (First-order gate)

[
\max_{a,b\in\mathcal G_A}\ \big|[[D,\pi(a)],J\pi(b)J^{-1}]\big|\le \varepsilon_1.
]

### Gate 10.3 (Fluctuation stability gates)

For (D_A=D+A+JAJ^{-1}):
[
|D_A-D_A^\dagger|\le \varepsilon_{\rm sa},
\qquad
\max_{a,b\in\mathcal G_A}\ \big|[[D_A,\pi(a)],J\pi(b)J^{-1}]\big|\le \varepsilon_A.
]

### Gate 10.4 (Projector-defined chirality)

With chiral projectors (P_L,P_R),
[
|P_L^2-P_L|,\ |P_R^2-P_R|,\ |P_LP_R| \ \text{small},\qquad
|P_L+P_R-1|\ \text{small on the chiral subspace}.
]

## Theorem 7.1 (Order-zero) — two-line proof

**Claim.** If order-zero holds on the geometric and SM factors, then for all (a,b\in A),
[
[\pi(a),,J\pi(b)J^{-1}]=0 .
]

**Proof.** It suffices to check on simple tensors (a=a_{\rm geom}\otimes a_{\rm SM}), (b=b_{\rm geom}\otimes b_{\rm SM}) (finite sums are algebraically dense). Using
[
\pi(a)=\pi_{\rm geom}(a_{\rm geom})\otimes\pi_{\rm SM}(a_{\rm SM})\otimes 1,\quad
J\pi(b)J^{-1}=\big(J_{\rm geom}\pi_{\rm geom}(b_{\rm geom})J_{\rm geom}^{-1}\big)\otimes\big(J_{\rm SM}\pi_{\rm SM}(b_{\rm SM})J_{\rm SM}^{-1}\big)\otimes 1,
]
the commutator factorizes as
[
[\pi(a),J\pi(b)J^{-1}]
======================

[\pi_{\rm geom}(a_{\rm geom}),J_{\rm geom}\pi_{\rm geom}(b_{\rm geom})J_{\rm geom}^{-1}]
\otimes \pi_{\rm SM}(a_{\rm SM}),J_{\rm SM}\pi_{\rm SM}(b_{\rm SM})J_{\rm SM}^{-1}\otimes 1
]
[
+;
\pi_{\rm geom}(a_{\rm geom}),J_{\rm geom}\pi_{\rm geom}(b_{\rm geom})J_{\rm geom}^{-1}
\otimes[\pi_{\rm SM}(a_{\rm SM}),J_{\rm SM}\pi_{\rm SM}(b_{\rm SM})J_{\rm SM}^{-1}]\otimes 1
=0
]
by the factor order-zero hypotheses. Extend from simple tensors to finite sums by bilinearity, and to the chosen (C^\ast)-completion by norm-density and continuity of ((x,y)\mapsto [x,y]) on (\mathcal B(H)). ∎

## Theorem 8.1 (First-order) — two-line proof

**Claim.** If first-order holds on the geometric factor and on the internal factor ((A_{\rm SM},H_{\rm SM}\otimes H_{\rm flav},D_{\rm int},J_{\rm SM}\otimes J_{\rm flav})), then for all (a,b\in A),
[
[[D,\pi(a)],,J\pi(b)J^{-1}]=0 .
]

**Proof.** Check on simple tensors (a=a_{\rm geom}\otimes a_{\rm SM}), (b=b_{\rm geom}\otimes b_{\rm SM}). With
[
D=D_{\rm geom}\otimes 1\otimes 1+\gamma_{\rm geom}\otimes D_{\rm int},\qquad
\pi(a)=\pi_{\rm geom}(a_{\rm geom})\otimes\pi_{\rm SM}(a_{\rm SM})\otimes 1,
]
we have
[
[D,\pi(a)]
==========

[D_{\rm geom},\pi_{\rm geom}(a_{\rm geom})]\otimes \pi_{\rm SM}(a_{\rm SM})\otimes 1
+\gamma_{\rm geom}\pi_{\rm geom}(a_{\rm geom})\otimes [D_{\rm int},,\pi_{\rm SM}(a_{\rm SM})\otimes 1].
]
Also (J\pi(b)J^{-1}=(J_{\rm geom}\pi_{\rm geom}(b_{\rm geom})J_{\rm geom}^{-1})\otimes(J_{\rm SM}\pi_{\rm SM}(b_{\rm SM})J_{\rm SM}^{-1})\otimes 1). Hence the double commutator splits as the sum of two tensor-factorized terms:
[
\big[[D_{\rm geom},\pi_{\rm geom}(a_{\rm geom})],,J_{\rm geom}\pi_{\rm geom}(b_{\rm geom})J_{\rm geom}^{-1}\big]\otimes(\cdots)\otimes 1
;+;
(\cdots)\otimes \big[[D_{\rm int},\pi_{\rm SM}(a_{\rm SM})\otimes 1],, (J_{\rm SM}\pi_{\rm SM}(b_{\rm SM})J_{\rm SM}^{-1})\otimes 1\big],
]
and each vanishes by the respective first-order hypotheses; the flavor multiplicity is inert because (\pi(A)) acts as (1_{\rm flav}). Extend to finite sums by bilinearity and to the completion by norm-density/continuity exactly as in Theorem 7.1. ∎

## Consolidated Factorization Lemma (v4.0)

*(One lemma that yields Sections 4, 7, 8 by inspection.)*

### Lemma F.1 (Tensor–commutator factorization with trivial flavor action)

Let
[
A=A_{\rm geom}\otimes A_{\rm SM},\quad
H=H_{\rm geom}\otimes H_{\rm SM}\otimes H_{\rm flav},
]
[
\pi(a_{\rm geom}!\otimes! a_{\rm SM})
=\pi_{\rm geom}(a_{\rm geom})\otimes \pi_{\rm SM}(a_{\rm SM})\otimes 1_{\rm flav},
]
and
[
D=D_{\rm geom}\otimes 1\otimes 1
+\gamma_{\rm geom}\otimes D_{\rm SM}\otimes 1
+\gamma_{\rm geom}\otimes 1\otimes D_{\rm flav}.
]
Let (J=J_{\rm geom}\otimes J_{\rm SM}\otimes J_{\rm flav}) and write, for (b=b_{\rm geom}\otimes b_{\rm SM}),
[
b^\circ := J\pi(b)J^{-1}
=\big(J_{\rm geom}\pi_{\rm geom}(b_{\rm geom})J_{\rm geom}^{-1}\big)\otimes
\big(J_{\rm SM}\pi_{\rm SM}(b_{\rm SM})J_{\rm SM}^{-1}\big)\otimes 1_{\rm flav}.
]

Then for all simple tensors (a=a_{\rm geom}\otimes a_{\rm SM}\in A), the following identities hold (and extend to all (a\in A) by linearity and continuity):

#### (F.1) **Single commutator** decomposes, flavor term vanishes identically

[
[D,\pi(a)]
==========

\underbrace{,[D_{\rm geom},\pi_{\rm geom}(a_{\rm geom})]\otimes \pi_{\rm SM}(a_{\rm SM})\otimes 1,}*{\text{geom piece}}
+
\underbrace{,\gamma*{\rm geom}\pi_{\rm geom}(a_{\rm geom})\otimes [D_{\rm SM},\pi_{\rm SM}(a_{\rm SM})]\otimes 1,}*{\text{SM piece}},
]
and
[
[\gamma*{\rm geom}\otimes 1\otimes D_{\rm flav},,\pi(a)]=0
\qquad(\text{Definition 1.2 / flavor multiplicity}).
]

#### (F.2) **Order-zero commutator** factorizes into factor order-zero commutators

[
[\pi(a),b^\circ]
================

[\pi_{\rm geom}(a_{\rm geom}),,b_{\rm geom}^\circ]\otimes
\pi_{\rm SM}(a_{\rm SM}),b_{\rm SM}^\circ\otimes 1
+
\pi_{\rm geom}(a_{\rm geom}),b_{\rm geom}^\circ\otimes
[\pi_{\rm SM}(a_{\rm SM}),,b_{\rm SM}^\circ]\otimes 1.
]

#### (F.3) **First-order double commutator** splits into a sum of *factor* first-order/zero checks

Let
[
X_1:=[D_{\rm geom},\pi_{\rm geom}(a_{\rm geom})],\quad Y_1:=\pi_{\rm SM}(a_{\rm SM}),
]
[
X_2:=\gamma_{\rm geom}\pi_{\rm geom}(a_{\rm geom}),\quad Y_2:=[D_{\rm SM},\pi_{\rm SM}(a_{\rm SM})],
]
so that ([D,\pi(a)]=(X_1\otimes Y_1+X_2\otimes Y_2)\otimes 1). Then
[
\big[[D,\pi(a)],b^\circ\big]
============================

\big([X_1,b_{\rm geom}^\circ]\otimes Y_1 b_{\rm SM}^\circ

* b_{\rm geom}^\circ X_1\otimes [Y_1,b_{\rm SM}^\circ]\big)\otimes 1
  ]
  [
  \qquad\qquad\quad
  +\big([X_2,b_{\rm geom}^\circ]\otimes Y_2 b_{\rm SM}^\circ
* b_{\rm geom}^\circ X_2\otimes [Y_2,b_{\rm SM}^\circ]\big)\otimes 1.
  ]

Moreover, if ([\gamma_{\rm geom},\pi_{\rm geom}(A_{\rm geom})]=0) and (J_{\rm geom}\gamma_{\rm geom}=\varepsilon''\gamma_{\rm geom}J_{\rm geom}) (even real triple), then ([\gamma_{\rm geom},b_{\rm geom}^\circ]=0), hence
[
[X_2,b_{\rm geom}^\circ]=[\gamma_{\rm geom}\pi_{\rm geom}(a_{\rm geom}),b_{\rm geom}^\circ]
=\gamma_{\rm geom}[\pi_{\rm geom}(a_{\rm geom}),b_{\rm geom}^\circ].
]

---

### Proof (mechanical; this is what becomes code)

Use the elementary tensor-commutator identity
[
[A\otimes B,\ C\otimes D]=[A,C]\otimes BD + CA\otimes[B,D],
]
(and its obvious extension with an extra (\otimes 1_{\rm flav})).

* **(F.1)** Expand ([D,\pi(a)]) term-by-term; the flavor term commutes because (\pi(a)) acts as (1_{\rm flav}).
* **(F.2)** Apply the identity to (\pi_{\rm geom}(a_{\rm geom})\otimes \pi_{\rm SM}(a_{\rm SM})) versus (b_{\rm geom}^\circ\otimes b_{\rm SM}^\circ).
* **(F.3)** Write ([D,\pi(a)]) as a sum of two simple tensors (geom piece + SM piece) and apply the same identity twice; the additional observation ([\gamma_{\rm geom},b_{\rm geom}^\circ]=0) follows from (J_{\rm geom}\gamma_{\rm geom}=\varepsilon''\gamma_{\rm geom}J_{\rm geom}) and ([\gamma_{\rm geom},\pi_{\rm geom}(A_{\rm geom})]=0).

∎

---

## Immediate corollaries = Sections 4, 7, 8 “for free”

### Corollary F.2 (Section 4: bounded commutators)

If ([D_{\rm geom},\pi_{\rm geom}(a_{\rm geom})]) and ([D_{\rm SM},\pi_{\rm SM}(a_{\rm SM})]) are bounded for generators, then ([D,\pi(a)]) is bounded for all (a\in A).
*(From (F.1); flavor contributes nothing.)*

### Corollary F.3 (Section 7: order-zero)

If
[
[\pi_{\rm geom}(a_{\rm geom}),b_{\rm geom}^\circ]=0,\qquad
[\pi_{\rm SM}(a_{\rm SM}),b_{\rm SM}^\circ]=0
]
for all (a_{\rm geom},b_{\rm geom}) and (a_{\rm SM},b_{\rm SM}), then ([\pi(a),b^\circ]=0) for all (a,b\in A).
*(From (F.2).)*

### Corollary F.4 (Section 8: first-order)

If first-order holds on each factor,
[
[[D_{\rm geom},\pi_{\rm geom}(a_{\rm geom})],b_{\rm geom}^\circ]=0,\qquad
[[D_{\rm SM},\pi_{\rm SM}(a_{\rm SM})],b_{\rm SM}^\circ]=0,
]
and order-zero holds on each factor (to kill the ([Y_1,b_{\rm SM}^\circ]) and ([\pi_{\rm geom}(a_{\rm geom}),b_{\rm geom}^\circ]) terms), then
[
[[D,\pi(a)],b^\circ]=0\quad\forall,a,b\in A.
]
*(From (F.3), with ([X_2,b_{\rm geom}^\circ]=\gamma_{\rm geom}[\pi_{\rm geom}(a_{\rm geom}),b_{\rm geom}^\circ]=0) under order-zero.)*

## Lemma F.1′ (Factorization Operator Identity)

Fix (a=a_{\rm geom}\otimes a_{\rm SM}), (b=b_{\rm geom}\otimes b_{\rm SM}\in A) and set
[
\pi(a)=\pi_{\rm geom}(a_{\rm geom})\otimes\pi_{\rm SM}(a_{\rm SM})\otimes 1,\qquad
b^\circ := J\pi(b)J^{-1}=b^\circ_{\rm geom}\otimes b^\circ_{\rm SM}\otimes 1,
]
with (b^\circ_{\rm geom}:=J_{\rm geom}\pi_{\rm geom}(b_{\rm geom})J_{\rm geom}^{-1}) and (b^\circ_{\rm SM}:=J_{\rm SM}\pi_{\rm SM}(b_{\rm SM})J_{\rm SM}^{-1}).
Define the “commutator building blocks”
[
C_{\rm geom}(a_{\rm geom}):=[D_{\rm geom},\pi_{\rm geom}(a_{\rm geom})],\qquad
C_{\rm SM}(a_{\rm SM}):=[D_{\rm SM},\pi_{\rm SM}(a_{\rm SM})].
]

Then, in the v4.0 alignment triple (with (\pi(A)\subset \mathcal B(H_{\rm geom}\otimes H_{\rm SM})\otimes 1_{\rm flav})), the following **single factorization identity** holds:
[
\boxed{
\begin{aligned}
[D,\pi(a)]
&=\Big(C_{\rm geom}(a_{\rm geom})\otimes \pi_{\rm SM}(a_{\rm SM})
;+;\gamma_{\rm geom}\pi_{\rm geom}(a_{\rm geom})\otimes C_{\rm SM}(a_{\rm SM})\Big)\otimes 1,
[3pt]
[\pi(a),b^\circ]
&=\Big([\pi_{\rm geom}(a_{\rm geom}),b^\circ_{\rm geom}]\otimes \pi_{\rm SM}(a_{\rm SM}),b^\circ_{\rm SM}
;+;\pi_{\rm geom}(a_{\rm geom}),b^\circ_{\rm geom}\otimes [\pi_{\rm SM}(a_{\rm SM}),b^\circ_{\rm SM}]\Big)\otimes 1,
[3pt]
\big[[D,\pi(a)],b^\circ\big]
&=\Big(
[C_{\rm geom}(a_{\rm geom}),b^\circ_{\rm geom}]\otimes \pi_{\rm SM}(a_{\rm SM}),b^\circ_{\rm SM}
;+;b^\circ_{\rm geom}C_{\rm geom}(a_{\rm geom})\otimes [\pi_{\rm SM}(a_{\rm SM}),b^\circ_{\rm SM}]
\
&\qquad\quad
+;[\gamma_{\rm geom}\pi_{\rm geom}(a_{\rm geom}),b^\circ_{\rm geom}]\otimes C_{\rm SM}(a_{\rm SM}),b^\circ_{\rm SM}
;+;b^\circ_{\rm geom}\gamma_{\rm geom}\pi_{\rm geom}(a_{\rm geom})\otimes [C_{\rm SM}(a_{\rm SM}),b^\circ_{\rm SM}]
\Big)\otimes 1,
\end{aligned}}
]
where (1) denotes (1_{\rm flav}). In particular,
[
[\gamma_{\rm geom}\otimes 1\otimes D_{\rm flav},,\pi(a)]=0
]
is built into the first line via the final (\otimes 1) (Definition 1.2).

*Proof.* Repeated use of ([X\otimes Y,,X'\otimes Y']=[X,X']\otimes YY'+X'X\otimes [Y,Y']), together with (\pi(\cdot)\otimes 1_{\rm flav}). ∎

---

## Corollaries (hence Sections 4/7/8)

### Corollary F.2 (Hence Section 4: bounded commutators)

If (C_{\rm geom}(a_{\rm geom})) and (C_{\rm SM}(a_{\rm SM})) are bounded on their factors for generators (hence for all (a) by linearity/continuity), then ([D,\pi(a)]\in\mathcal B(H)) for all (a\in A).

### Corollary F.3 (Hence Section 7: order-zero)

If ([\pi_{\rm geom}(a_{\rm geom}),b^\circ_{\rm geom}]=0) and ([\pi_{\rm SM}(a_{\rm SM}),b^\circ_{\rm SM}]=0) for all (a_{\rm geom},b_{\rm geom},a_{\rm SM},b_{\rm SM}), then ([\pi(a),b^\circ]=0) for all (a,b\in A).

### Corollary F.4 (Hence Section 8: first-order)

Assume the factor first-order conditions
[
[C_{\rm geom}(a_{\rm geom}),b^\circ_{\rm geom}]=0,\qquad [C_{\rm SM}(a_{\rm SM}),b^\circ_{\rm SM}]=0,
]
and the factor order-zero conditions as in Corollary F.3. If moreover ([\gamma_{\rm geom},b^\circ_{\rm geom}]=0) (true in an even real triple), then the boxed third line collapses to
[
[[D,\pi(a)],b^\circ]=0\quad\forall,a,b\in A.
]

# The Unified Alignment Dirac Operator v4.0

## 0. Design constraints that v4.0 must satisfy

### Evenness

The full triple must admit a grading (\Gamma) with
[
\Gamma^\ast=\Gamma,\qquad \Gamma^2=1,\qquad {\Gamma,D}=0,\qquad [\Gamma,\pi(A)]=0.
]
Therefore **no standalone “pure flavor” term** of the form (\gamma_{\rm geom}\otimes 1\otimes D_{\rm flav}) is allowed, because it commutes with (\Gamma) and would spoil ({\Gamma,D}=0).

### Flavor multiplicity

Flavor is a **multiplicity factor**:
[
\pi(A)\subseteq \mathcal B(H_{\rm geom}\otimes H_{\rm SM})\otimes 1_{\rm flav},
]
so there is **no** extra commutative diagonal “site algebra” acting on (H_{\rm flav}). This is the v4.0 structural cure that protects order-zero and first-order.

### Commutant protection

All flavor/Yukawa textures are chosen in the commutant:
[
Y[\mathcal K]\in \pi(A)',,
]
so inner fluctuations do not generate “generation-gauge” bosons; they only produce the usual gauge/Higgs degrees of freedom.

---

## 1. The spectral datum and the even product Dirac operator

### 1.1 Spaces and algebra

Assume:

* A real, even geometric triple ((A_{\rm geom},H_{\rm geom},D_{\rm geom},J_{\rm geom},\gamma_{\rm geom})).
* A real, even finite SM triple ((A_{\rm SM},H_{\rm SM},D_{\rm SM},J_{\rm SM},\gamma_{\rm SM})) (one-generation internal geometry).
* A finite-dimensional multiplicity space (H_{\rm flav}\cong\mathbb C^N) with antiunitary (J_{\rm flav}=U_{\rm flav}\circ K).

Define
[
A:=A_{\rm geom}\otimes A_{\rm SM},\qquad
H:=H_{\rm geom}\otimes H_{\rm SM}\otimes H_{\rm flav},
]
[
\pi(a_{\rm geom}\otimes a_{\rm SM})
===================================

\pi_{\rm geom}(a_{\rm geom})\otimes \pi_{\rm SM}(a_{\rm SM})\otimes 1_{\rm flav}.
]

### 1.2 The even product Dirac operator

The v4.0 unfluctuated Dirac operator is **by definition**
[
D
=

D_{\rm geom}\otimes 1\otimes 1
;+;
\gamma_{\rm geom}\otimes D_{\rm int}(Y[\mathcal K]),
]
where (D_{\rm int}(Y[\mathcal K])) acts on (H_{\rm SM}\otimes H_{\rm flav}) and is **odd** with respect to (\gamma_{\rm SM}\otimes 1):
[
{\gamma_{\rm SM}\otimes 1,;D_{\rm int}(Y[\mathcal K])}=0.
]

The total grading and real structure are
[
\Gamma:=\gamma_{\rm geom}\otimes \gamma_{\rm SM}\otimes 1_{\rm flav},\qquad
J:=J_{\rm geom}\otimes J_{\rm SM}\otimes J_{\rm flav}.
]

This is the v4.0 “even product” lock: all mixing must appear **inside** the odd internal operator (D_{\rm int}), never as an even add-on.

---

## 2. Where the hierarchy and mixing live (without breaking axioms)

### 2.1 Flavor textures as commutant operators

The alignment hierarchy is encoded in operators on (H_{\rm flav}), for example using your (\kappa)-distance kernel on sites ({1,2,5}):
[
(Y_{\rm flav})*{ij}=\kappa^{|n_i-n_j|},
\qquad {n_i}={1,2,5}.
]
But in v4.0 this object is **not** itself “the finite algebra”; it is a **commutant texture** acting on the multiplicity factor and inserted into the LR blocks of (D*{\rm int}).

Formally, require
[
Y_{\rm flav}\in \pi(A)'\quad\text{(commutant condition)}.
]

### 2.2 The internal Dirac operator as an odd LR-coupler

Write the internal operator in chiral block form (on (H_{\rm SM}\otimes H_{\rm flav})):
[
D_{\rm int}(Y[\mathcal K])
==========================

\begin{pmatrix}
0 & \mathcal Y^\dagger \
\mathcal Y & 0
\end{pmatrix}
;\oplus;
\text{(possible Majorana blocks)}.
]
Here (\mathcal Y) is the LR coupling operator, built sector-by-sector using SM projectors and your flavor textures. Schematically:
[
\mathcal Y
==========

P_u\otimes Y_u[\mathcal K]
+
P_d\otimes Y_d[\mathcal K]
+
P_e\otimes Y_e[\mathcal K]
+
P_\nu\otimes Y_\nu[\mathcal K],
]
with (P_u,P_d,P_e,P_\nu) acting on (H_{\rm SM}) and (Y_x[\mathcal K]) acting on (H_{\rm flav}). Each (Y_x[\mathcal K]) is generated by your kernel pipeline and lies in (\pi(A)').

### 2.3 Neutrino seesaw (if included)

If you include Majorana masses, they belong in the internal operator as additional blocks (Takagi/real structure compatible), not as a “fluctuation trick”. The reality structure (J) constrains the allowed form; it does not “create” the scale by itself.

---

## 3. Inner fluctuations and the unified operator

### 3.1 One-forms and fluctuations

Define the represented one-forms
[
\Omega_D^1(A)
:=
\left{
\sum_i \pi(a_i)[D,\pi(b_i)] : a_i,b_i\in A
\right}
\subset \mathcal B(H).
]
Given (A_1\in\Omega_D^1(A)), take its Hermitian part
[
A:=\tfrac12(A_1+A_1^\dagger),
]
and define the fluctuated operator
[
D_A:=D + A + JAJ^{-1}.
]

### 3.2 What fluctuations generate (and what they do not)

* **Generated by fluctuations:** gauge bosons and Higgs fields associated to the represented algebra (A=A_{\rm geom}\otimes A_{\rm SM}).
* **Not generated by fluctuations:** Yukawa textures. Those live in (D_{\rm int}(Y[\mathcal K])) by construction; fluctuations Higgs-dress them.

### 3.3 The antiunitary implementation rule

In a fixed matrix basis, if (J=U_J\circ K), then
[
JXJ^{-1}=U_J,\overline{X},U_J^\dagger.
]
This must be used for the (JAJ^{-1}) term.

---

## 4. The final unified operator (v4.0 canonical form)

[
\boxed{
D_A
===

\underbrace{D_{\rm geom}\otimes 1\otimes 1}*{\text{geometric phase-gradient}}
+
\underbrace{\gamma*{\rm geom}\otimes D_{\rm int}(Y[\mathcal K])}*{\text{odd internal LR + textures}}
+
\underbrace{A}*{\text{gauge + Higgs from inner fluctuations}}
+
\underbrace{JAJ^{-1}}_{\text{reality completion}}
}
]

This is the single operator compatible with:

* evenness,
* flavor multiplicity,
* commutant protection,
* order-zero and first-order stability.

---

## 5. Spectral action (what becomes computable once (D_A) is fixed)

With
[
S=\operatorname{Tr},f!\left(\frac{D_A^2}{\Lambda^2}\right),
]
the standard heat-kernel expansion yields:

* geometric terms (from (D_{\rm geom})),
* gauge kinetic terms (from fluctuations),
* Higgs kinetic and potential (from internal fluctuations),
* Yukawa normalization and fermion masses (from (D_{\rm int}(Y[\mathcal K])) Higgs-dressed),
* neutrino seesaw structure (if Majorana blocks are included in (D_{\rm int})).

The unification is real: everything is encoded in one spectral datum—without violating the axioms.

---

## 6. Production commutator gates (norm-corrected)

Fix a finite spanning/generating set (\mathcal G_A\subset A). On truncations use the operator norm (|\cdot|).

[
\max_{a,b\in\mathcal G_A}\ \big|[\pi(a),J\pi(b)J^{-1}]\big|\le \varepsilon_0,
]
[
\max_{a,b\in\mathcal G_A}\ \big|[[D,\pi(a)],J\pi(b)J^{-1}]\big|\le \varepsilon_1,
]
[
|D_A-D_A^\dagger|\le \varepsilon_{\rm sa},
\qquad
\max_{a,b\in\mathcal G_A}\ \big|[[D_A,\pi(a)],J\pi(b)J^{-1}]\big|\le \varepsilon_A.
]

---

## 7. Updated implementation skeleton (v4.0 structural correctness)

This skeleton only encodes the **structural** requirements: even product form, Hermitian one-forms, antiunitary (J), and projector-defined Yukawa extraction.

```python
import numpy as np

def antiunitary_conjugation(UJ: np.ndarray, X: np.ndarray) -> np.ndarray:
    """
    J X J^{-1} for J = UJ ∘ K.
    """
    return UJ @ X.conj() @ UJ.conj().T


class UnifiedAlignmentDiracV4:
    """
    Builds the v4.0 even-product unified operator:
      D   = D_geom ⊗ 1 ⊗ 1  +  gamma_geom ⊗ D_int(Y[K])
      D_A = D + A + J A J^{-1}
    """

    def __init__(self,
                 D_geom_full: np.ndarray,
                 gamma_geom_full: np.ndarray,
                 D_int_full: np.ndarray,
                 rep_pairs: list,
                 UJ_full: np.ndarray):
        """
        All inputs are already represented on the full H = H_geom ⊗ H_SM ⊗ H_flav.

        rep_pairs: list of (a, b) represented operators used in A = Σ a [D, b]
        UJ_full: unitary implementing the antiunitary J = UJ ∘ K
        """
        self.D_geom = D_geom_full
        self.gamma_geom = gamma_geom_full
        self.D_int = D_int_full
        self.rep_pairs = rep_pairs
        self.UJ = UJ_full

        self.D = None
        self.A = None
        self.DA = None

    def build_D(self):
        # Even product form
        self.D = self.D_geom + (self.gamma_geom @ self.D_int)
        return self.D

    def build_A(self):
        if self.D is None:
            self.build_D()

        A1 = np.zeros_like(self.D, dtype=complex)
        for (a, b) in self.rep_pairs:
            comm = self.D @ b - b @ self.D
            A1 += a @ comm

        # Hermitian one-form
        self.A = 0.5 * (A1 + A1.conj().T)
        return self.A

    def build_DA(self):
        if self.D is None:
            self.build_D()
        if self.A is None:
            self.build_A()

        JAJinv = antiunitary_conjugation(self.UJ, self.A)
        self.DA = self.D + self.A + JAJinv
        return self.DA


def extract_lr(D_int: np.ndarray, P_L: np.ndarray, P_R: np.ndarray) -> np.ndarray:
    """
    Projector-defined LR extraction (no row slicing).
    Convention: L → R block is P_R D_int P_L.
    """
    return P_R @ D_int @ P_L
```
# The Alignment Spectral Action v4.0

## 1. The spectral action principle (unchanged, now v4.0-correct)

The spectral action is still the same master functional:
[
\boxed{
S_{\rm spec}[D_A,\Psi]
======================

\operatorname{Tr}, f!\left(\frac{D_A^2}{\Lambda^2}\right)
;+;
\langle \Psi,; D_A \Psi\rangle
}
]
where (f) is positive, smooth, rapidly decaying, and (D_A) is the **inner-fluctuated** Dirac operator.

The v4.0 update is not the principle—it is the **definition of (D)** (and therefore of (D_A)).

---

## 2. v4.0 unified operator input to the spectral action

### 2.1 Spectral datum (recall)

[
A=A_{\rm geom}\otimes A_{\rm SM},\qquad
H=H_{\rm geom}\otimes H_{\rm SM}\otimes H_{\rm flav},
]
[
\pi(a_{\rm geom}\otimes a_{\rm SM})
===================================

\pi_{\rm geom}(a_{\rm geom})\otimes \pi_{\rm SM}(a_{\rm SM})\otimes 1_{\rm flav}.
]

Flavor is a **multiplicity factor**, and all flavor/Yukawa textures satisfy the **commutant condition**
[
Y[\mathcal K]\in \pi(A)',.
]

### 2.2 The even product Dirac operator (the key replacement)

The v2.1 form
[
D_{\rm geom}\otimes 1\otimes 1
;+;
\gamma_{\rm geom}\otimes D_{\rm SM}\otimes 1
;+;
\gamma_{\rm geom}\otimes 1\otimes D_F
]
is replaced in v4.0 by the **even product form**
[
\boxed{
D
=

D_{\rm geom}\otimes 1\otimes 1
;+;
\gamma_{\rm geom}\otimes D_{\rm int}(Y[\mathcal K])
}
]
where (D_{\rm int}(Y[\mathcal K])) acts on (H_{\rm SM}\otimes H_{\rm flav}) and is **odd** with respect to (\gamma_{\rm SM}\otimes 1):
[
\boxed{
{\gamma_{\rm SM}\otimes 1,;D_{\rm int}(Y[\mathcal K])}=0.
}
]

This single change enforces full evenness with
[
\Gamma=\gamma_{\rm geom}\otimes\gamma_{\rm SM}\otimes 1_{\rm flav},
\qquad
{\Gamma,D}=0,
]
and removes the forbidden “pure flavor even add-on” obstruction.

### 2.3 Inner fluctuations and the unified fluctuated operator

One-forms:
[
A_1=\sum_i \pi(a_i),[D,\pi(b_i)],\qquad a_i,b_i\in A,
\qquad
A:=\tfrac12(A_1+A_1^\dagger),
]
and the fluctuated Dirac operator:
[
\boxed{
D_A = D + A + JAJ^{-1}.
}
]
Implementation rule (antiunitary (J=U_J\circ K)):
[
JXJ^{-1}=U_J,\overline{X},U_J^\dagger.
]

### 2.4 What fluctuations generate (v4.0 clarification)

* **Generated by** (\operatorname{Tr} f(D_A^2/\Lambda^2)): gauge bosons + Higgs sector (from the represented algebra (A_{\rm SM}) and the geometric factor as appropriate).
* **Not generated by** fluctuations: Yukawa textures. Those are encoded in (D_{\rm int}(Y[\mathcal K])) (LR blocks) and are Higgs-dressed by fluctuations.

This is the correct v4.0 separation: **textures in (D_{\rm int})**, **fields in (A)**.

---

## 3. Bosonic spectral action in v4.0 form

The bosonic action is
[
S_{\rm bos}=\operatorname{Tr}, f!\left(\frac{D_A^2}{\Lambda^2}\right).
]

In four dimensions, the heat-kernel expansion takes the standard form
[
\operatorname{Tr}, f!\left(\frac{D_A^2}{\Lambda^2}\right)
\sim
\sum_{k\ge 0} F_{4-k},\Lambda^{4-k}, a_k(D_A^2),
]
where (a_k(D_A^2)) are the Seeley–DeWitt coefficients and the (F_\ell) are moments of (f) (one convenient normalization is)
[
F_4=\int_0^\infty f(u),u,du,\qquad
F_2=\int_0^\infty f(u),du,\qquad
F_0=f(0),
]
with higher terms suppressed by inverse powers of (\Lambda) under suitable regularity assumptions.

Structurally, the expansion yields:
[
\mathcal L_{\rm bos}
====================

\mathcal L_{\rm grav}
+\mathcal L_{\rm gauge}
+\mathcal L_{H}
+\text{(higher-curvature / higher-dimension terms)}.
]

In v4.0 the **only** change is what sits inside (D_A): the internal flavor/mixing content is now located in (D_{\rm int}(Y[\mathcal K])) (odd block), rather than in a separate (D_F) term.

---

## 4. Fermionic spectral action in v4.0 form

The fermionic action is
[
S_{\rm ferm}=\langle \Psi,; D_A\Psi\rangle.
]
This produces:

* fermion kinetic terms,
* gauge-fermion couplings (through (A)),
* Higgs-fermion couplings (through the Higgs component of (A)),
* Yukawa interactions because the LR blocks of (D_{\rm int}(Y[\mathcal K])) carry the textures (Y[\mathcal K]),
* Majorana/seesaw structure if (D_{\rm int}) includes the appropriate Majorana blocks compatible with (J).

In short: **Yukawas appear because they are in (D_{\rm int})**, and become physical couplings because the Higgs appears in the fluctuated operator.

---

## 5. Alignment-harmonic cutoff function (f_{\rm align})

Updated to guarantee NCG positivity while encoding harmonic structure

Your original idea (divisor enhancement + triadic emphasis + phase coherence + forbidden-mode suppression) is compatible with NCG **as long as (f) remains positive and rapidly decaying**.

A safe “Alignment-coded” class is:
[
\boxed{
f_{\rm align}(x)
================

e^{-x};
\exp!\Bigg(
\sum_{d\in D_{360}}\alpha_d\cos!\frac{x}{d}
;+;
\sum_{m=1}^M \beta_m\cos(3m,x)
;-;
\gamma,(1-\cos(\varphi x))
;-;
\sum_{r\in\mathcal F}\delta_r,e^{-(x-\lambda_r)^2}
\Bigg),
}
]
with parameters chosen so the exponent is bounded above (automatic for finite sums) and (\gamma,\delta_r\ge 0). This function is:

* smooth,
* strictly positive,
* rapidly decaying (due to (e^{-x})),
* “Alignment-coded” through divisor/triadic/phase/forbidden-mode factors.

This preserves the spectral action axioms while allowing harmonic weighting in a controlled way.

---

## 6. Drop-in LaTeX replacement (your section, rewritten for v4.0)

```latex
\section{The Alignment Spectral Action (v4.0)}

In the Alignment Spectral Triple v4.0,
\[
(A,H,D,J,\Gamma)
=
\left(
A_{\mathrm{geom}} \otimes A_{\mathrm{SM}},\;
H_{\mathrm{geom}} \otimes H_{\mathrm{SM}} \otimes H_{\mathrm{flav}},\;
D,\;
J,\;
\Gamma
\right),
\]
the full physical Lagrangian is generated by the spectral action principle:
\begin{equation}
S_{\mathrm{spec}}
=
\mathrm{Tr}\, f(D_A^2/\Lambda^2)
\;+\;
\langle \Psi, D_A \Psi\rangle.
\end{equation}
Here $\Lambda$ is a cutoff scale, $f$ is smooth, positive, and rapidly decaying,
and $D_A$ is the inner-fluctuated Dirac operator.

\subsection{Even product Dirac operator (v4.0)}
The unfluctuated Dirac operator is taken in the even product form
\begin{equation}
D
=
D_{\mathrm{geom}}\otimes 1\otimes 1
+
\gamma_{\mathrm{geom}} \otimes D_{\mathrm{int}}(Y[\mathcal K]),
\end{equation}
where $D_{\mathrm{int}}(Y[\mathcal K])$ acts on $H_{\mathrm{SM}}\otimes H_{\mathrm{flav}}$
and is odd with respect to $\gamma_{\mathrm{SM}}\otimes 1$:
\begin{equation}
\{\gamma_{\mathrm{SM}}\otimes 1,\; D_{\mathrm{int}}(Y[\mathcal K])\}=0.
\end{equation}
Flavor is a pure multiplicity factor: the represented algebra acts trivially on
$H_{\mathrm{flav}}$, and the Yukawa/flavor textures are chosen in the commutant
$Y[\mathcal K]\in\pi(A)'$.

\subsection{Inner fluctuations}
Inner fluctuations arise from elements $a_i,b_i\in A$:
\begin{equation}
A_1=\sum_i \pi(a_i)[D,\pi(b_i)],\qquad A=\tfrac12(A_1+A_1^\dagger),
\qquad
D_A = D + A + JAJ^{-1}.
\end{equation}
The real structure $J$ is antiunitary and in a fixed basis satisfies
$JXJ^{-1}=U_J\,\overline{X}\,U_J^\dagger$.

\subsection{Bosonic spectral action}
The bosonic part is
\begin{equation}
S_{\mathrm{bos}} = \mathrm{Tr}\, f(D_A^2/\Lambda^2).
\end{equation}
Using the heat-kernel expansion,
\begin{equation}
\mathrm{Tr}\,f(D_A^2/\Lambda^2)
\sim
\sum_{k\ge 0}
F_{4-k}\,\Lambda^{4-k}\, a_k(D_A^2),
\end{equation}
one obtains gravitational terms, gauge kinetic terms for the
$U(1)\times SU(2)\times SU(3)$ sector encoded by $A_{\mathrm{SM}}$, and the Higgs
kinetic and potential terms arising from the internal fluctuations.

\subsection{Fermionic spectral action}
The fermionic part is
\begin{equation}
S_{\mathrm{ferm}} = \langle \Psi, D_A \Psi \rangle,
\end{equation}
producing fermion kinetic terms, gauge couplings, Higgs couplings, and Yukawa
interactions through the LR blocks of $D_{\mathrm{int}}(Y[\mathcal K])$ (with
textures supplied by the Alignment kernel pipeline).

\subsection{Alignment-harmonic cutoff}
A representative Alignment-coded cutoff is
\begin{equation}
f_{\mathrm{align}}(x)
=
e^{-x}\exp\!\Bigg(
\sum_{d\in D_{360}}\alpha_d\cos\!\frac{x}{d}
+\sum_{m=1}^M\beta_m\cos(3m\,x)
-\gamma(1-\cos(\varphi x))
-\sum_{r\in\mathcal F}\delta_r e^{-(x-\lambda_r)^2}
\Bigg),
\end{equation}
which remains smooth, positive, and rapidly decaying while encoding divisor
harmonicity, triadic closure emphasis, phase coherence, and forbidden-mode
suppression.
```
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

## Lifted operators and an explicit (32N\times 32N) sparse blueprint for (D_{\rm int})

We keep the **exact ordering** already fixed:

* Particle indices (0\ldots15): ((u_L^r,u_L^g,u_L^b,d_L^r,d_L^g,d_L^b,\nu_L,e_L,u_R^r,u_R^g,u_R^b,d_R^r,d_R^g,d_R^b,e_R,\nu_R)).
* Conjugates (16\ldots31): same order, shifted by (+16).

Thus
[
H_{\rm int}=H_{\rm SM}\otimes H_{\rm flav},\qquad \dim(H_{\rm int})=32N.
]

### 1) The lifted operators (P\otimes 1_N) and (1\otimes Y)

Let (I_N) be the (N\times N) identity.

* For any (32\times 32) SM projector (P),
  [
  \widetilde P := P\otimes I_N \in M_{32N}(\mathbb C).
  ]

* For any flavor texture (Y\in M_N(\mathbb C)),
  [
  \widetilde Y := I_{32}\otimes Y \in M_{32N}(\mathbb C).
  ]

More generally, for a matrix unit (E_{ij}) on (H_{\rm SM}) (one at ((i,j))),
[
\widetilde E_{ij}(Y) := E_{ij}\otimes Y,
]
which is exactly “place the (N\times N) block (Y)” in the ((i,j)) SM slot.

This (\widetilde E_{ij}(Y)) is the cleanest literal sparse-building primitive.

---

## 2) Grading-oddness target (so v4.0 even product holds)

With the earlier choice
[
\gamma_{\rm SM}=\mathrm{diag}(+1\ \text{on}\ {0,\dots,7}\cup{24,\dots,31},;;-1\ \text{on}\ {8,\dots,15}\cup{16,\dots,23}),
]
oddness is equivalent to: **(D_{\rm int}) only connects (+\leftrightarrow-)**.
Our template below enforces this automatically.

---

## 3) Explicit block-sparse (D_{\rm int}) (Dirac + Majorana) in this ordering

Let (Y_u,Y_d,Y_e,Y_\nu\in M_N(\mathbb C)) be your Alignment textures (in the commutant (\pi(A)')). Let (M_R\in M_N(\mathbb C)) be the heavy Majorana matrix (usually Takagi-symmetric: (M_R^T=M_R) in a (J_{\rm flav}=K) basis).

### 3.1 Dirac (LR) couplings: particles

For each color channel, couple (u_L^c\leftrightarrow u_R^c) and (d_L^c\leftrightarrow d_R^c). Couple leptons (\nu_L\leftrightarrow\nu_R), (e_L\leftrightarrow e_R).

Using the SM indices (particle):

* (u_L:{0,1,2}), (u_R:{8,9,10})
* (d_L:{3,4,5}), (d_R:{11,12,13})
* (\nu_L:{6}), (\nu_R:{15})
* (e_L:{7}), (e_R:{14})

Add the terms:
[
\sum_{c=0}^2 \Big(\widetilde E_{,8+c,;0+c}(Y_u)+\widetilde E_{,0+c,;8+c}(Y_u^\dagger)\Big)
;+;
\sum_{c=0}^2 \Big(\widetilde E_{,11+c,;3+c}(Y_d)+\widetilde E_{,3+c,;11+c}(Y_d^\dagger)\Big)
]
[
+;\widetilde E_{,14,;7}(Y_e)+\widetilde E_{,7,;14}(Y_e^\dagger)
;+;
\widetilde E_{,15,;6}(Y_\nu)+\widetilde E_{,6,;15}(Y_\nu^\dagger).
]

### 3.2 Dirac couplings: conjugates

Conjugate indices are shifted by (+16). The conjugate Dirac blocks use complex-conjugated Yukawas to match the charge-conjugate sector (and keep self-adjointness):
[
\sum_{c=0}^2 \Big(\widetilde E_{,16+(0+c),;16+(8+c)}(\overline{Y_u})
+\widetilde E_{,16+(8+c),;16+(0+c)}(Y_u^{T})\Big)
]
[
+\sum_{c=0}^2 \Big(\widetilde E_{,16+(3+c),;16+(11+c)}(\overline{Y_d})
+\widetilde E_{,16+(11+c),;16+(3+c)}(Y_d^{T})\Big)
]
[
+;\widetilde E_{,23,;30}(\overline{Y_e})
+\widetilde E_{,30,;23}(Y_e^{T})
;+;
\widetilde E_{,22,;31}(\overline{Y_\nu})
+\widetilde E_{,31,;22}(Y_\nu^{T}).
]
(Here (22=16+6), (23=16+7), (30=16+14), (31=16+15).)

### 3.3 Majorana seesaw block (right-handed neutrino only)

Couple (\nu_R) (index (15)) to (\nu_R^c) (index (31)):
[
\widetilde E_{,31,;15}(M_R);+;\widetilde E_{,15,;31}(M_R^\dagger).
]
If you want the Takagi condition explicitly: impose (M_R^T=M_R) (in the same basis where (J_{\rm flav}=K)).

### 3.4 Final definition

[
\boxed{
D_{\rm int}
===========

D_{\rm Dirac}^{\rm part}
+
D_{\rm Dirac}^{\rm conj}
+
D_{\rm Maj}.
}
]
By construction:

* (D_{\rm int}^\dagger=D_{\rm int}),
* ({\gamma_{\rm SM}\otimes 1,;D_{\rm int}}=0) (only (+\leftrightarrow-) couplings are used).

---

## 4) Literal assembly code (block-sparse by SM slots)

This constructs an explicit (32N\times 32N) matrix by writing only (N\times N) blocks into the required SM slots.

```python
import numpy as np

# ----------------------------
# Index convention (fixed)
# ----------------------------
uL = [0, 1, 2]
dL = [3, 4, 5]
nuL = [6]
eL  = [7]

uR = [8, 9, 10]
dR = [11, 12, 13]
eR = [14]
nuR = [15]

def put_block(D, i, j, B):
    """
    Place an NxN block B into the (i,j) SM slot of a 32N x 32N matrix D.
    """
    N = B.shape[0]
    r = slice(i*N, (i+1)*N)
    c = slice(j*N, (j+1)*N)
    D[r, c] += B

def assemble_D_int(N, Yu, Yd, Ye, Ynu, MR=None):
    """
    Returns D_int as an explicit (32N x 32N) matrix.

    Yu, Yd, Ye, Ynu: NxN Yukawa textures on H_flav.
    MR: NxN Majorana matrix (optional). If None, Majorana block is omitted.
    """
    D = np.zeros((32*N, 32*N), dtype=complex)

    # --- Particle Dirac blocks: L <-> R
    for k in range(3):  # colors r,g,b
        put_block(D, uR[k], uL[k], Yu)
        put_block(D, uL[k], uR[k], Yu.conj().T)

        put_block(D, dR[k], dL[k], Yd)
        put_block(D, dL[k], dR[k], Yd.conj().T)

    put_block(D, eR[0], eL[0], Ye)
    put_block(D, eL[0], eR[0], Ye.conj().T)

    put_block(D, nuR[0], nuL[0], Ynu)
    put_block(D, nuL[0], nuR[0], Ynu.conj().T)

    # --- Conjugate Dirac blocks (indices +16):
    # convention: (block from R^c -> L^c) uses conj(Y), and adjoint uses Y^T
    for k in range(3):
        put_block(D, (uL[k]+16), (uR[k]+16), Yu.conj())
        put_block(D, (uR[k]+16), (uL[k]+16), Yu.T)

        put_block(D, (dL[k]+16), (dR[k]+16), Yd.conj())
        put_block(D, (dR[k]+16), (dL[k]+16), Yd.T)

    put_block(D, (eL[0]+16), (eR[0]+16), Ye.conj())
    put_block(D, (eR[0]+16), (eL[0]+16), Ye.T)

    put_block(D, (nuL[0]+16), (nuR[0]+16), Ynu.conj())
    put_block(D, (nuR[0]+16), (nuL[0]+16), Ynu.T)

    # --- Majorana: nu_R <-> nu_R^c
    if MR is not None:
        i = nuR[0]       # 15
        ic = nuR[0] + 16 # 31
        put_block(D, ic, i, MR)
        put_block(D, i, ic, MR.conj().T)

    return D
```

### Optional sanity checks (oddness + self-adjointness)

If you also build (\gamma_{\rm SM}) as a (32\times 32) diagonal, lift it to (32N) by (\gamma_{\rm SM}\otimes I_N), then check:
[
|D_{\rm int}-D_{\rm int}^\dagger|\approx 0,\qquad
|{\gamma_{\rm SM}\otimes I_N,;D_{\rm int}}|\approx 0.
]

---

## 5) How this plugs into v4.0 immediately

Once (D_{\rm int}) is assembled as above, the v4.0 even product operator is literally:
[
D
=

D_{\rm geom}\otimes 1_{32N}
;+;
\gamma_{\rm geom}\otimes D_{\rm int},
\qquad
D_A = D + A + JAJ^{-1},
]
with (J) implemented antiunitarily.

## Lifted projectors and lifted textures

Fix (H_{\rm int}=H_{\rm SM}\otimes H_{\rm flav}) with (\dim H_{\rm SM}=32), (\dim H_{\rm flav}=N).

### Lifted SM projectors

For any projector (P\in M_{32}(\mathbb C)),
[
\widetilde P ;:=; P\otimes I_N ;\in; M_{32N}(\mathbb C).
]
In particular, with the index sets (S\subset{0,\dots,31}) from the fixed ordering,
[
P_S=\sum_{i\in S}E_{ii},\qquad \widetilde P_S=P_S\otimes I_N.
]

### Lifted flavor textures

For any flavor matrix (Y\in M_N(\mathbb C)),
[
\widetilde Y ;:=; I_{32}\otimes Y;\in;M_{32N}(\mathbb C).
]
Equivalently, on a chosen SM slot ((i,j)) we use the block-primitive
[
E_{ij}\otimes Y \in M_{32N}(\mathbb C),
]
which is exactly “place the (N\times N) block (Y) in SM position ((i,j)).”

---

## Partial isometries (matrix units) and “projector language”

Let (E_{ij}\in M_{32}(\mathbb C)) be matrix units on (H_{\rm SM}):
[
E_{ij}e_j=e_i,\qquad E_{ij}^\dagger=E_{ji},\qquad E_{ij}E_{kl}=\delta_{jk}E_{il}.
]
Interpret (E_{ij}) as a **partial isometry** from the (j)-slot to the (i)-slot with
[
E_{ij}^\dagger E_{ij}=E_{jj},\qquad E_{ij}E_{ij}^\dagger=E_{ii}.
]
Lift them by
[
\widetilde E_{ij}:=E_{ij}\otimes I_N.
]

Then a Yukawa coupling “from (j) to (i)” with texture (Y) is literally
[
(E_{ij}\otimes Y);=;(\widetilde E_{ij})(\widetilde Y).
]
This is the symbol-for-symbol match to the code primitive “put an (N\times N) block into SM slot ((i,j)).”

---

## Explicit (D_{\rm int}) as sums of partial isometries times (\widetilde Y)

Use the fixed particle indices:
[
u_L^c\in{0,1,2},\quad d_L^c\in{3,4,5},\quad \nu_L=6,\quad e_L=7,
]
[
u_R^c\in{8,9,10},\quad d_R^c\in{11,12,13},\quad e_R=14,\quad \nu_R=15,
]
and conjugates shifted by (+16).

Let (Y_u,Y_d,Y_e,Y_\nu\in M_N(\mathbb C)) and optionally (M_R\in M_N(\mathbb C)) (Majorana).

### 1) Particle Dirac sector (D_{\rm Dirac}^{\rm part})

Define the particle partial-isometry sums (color channels (c=0,1,2)):
[
V_u:=\sum_{c=0}^2 E_{,8+c,;0+c},\qquad
V_d:=\sum_{c=0}^2 E_{,11+c,;3+c},
]
[
V_e:=E_{14,7},\qquad V_\nu:=E_{15,6}.
]

Lift them: (\widetilde V_\bullet = V_\bullet\otimes I_N). Then
[
\boxed{
D_{\rm Dirac}^{\rm part}
========================

(\widetilde V_u)(I_{32}\otimes Y_u)+(\widetilde V_u)^\dagger(I_{32}\otimes Y_u^\dagger)
+
(\widetilde V_d)(I_{32}\otimes Y_d)+(\widetilde V_d)^\dagger(I_{32}\otimes Y_d^\dagger)
}
]
[
\boxed{
\qquad\quad
+
(\widetilde V_e)(I_{32}\otimes Y_e)+(\widetilde V_e)^\dagger(I_{32}\otimes Y_e^\dagger)
+
(\widetilde V_\nu)(I_{32}\otimes Y_\nu)+(\widetilde V_\nu)^\dagger(I_{32}\otimes Y_\nu^\dagger).
}
]

This is exactly the particle part of the earlier “put_block” code: each (E_{ij}\otimes Y) term is one block insertion.

### 2) Conjugate Dirac sector (D_{\rm Dirac}^{\rm conj})

Define the conjugate partial isometries:
[
V_u^c:=\sum_{c=0}^2 E_{,16+(0+c),;16+(8+c)},\qquad
V_d^c:=\sum_{c=0}^2 E_{,16+(3+c),;16+(11+c)},
]
[
V_e^c:=E_{,16+7,;16+14}=E_{23,30},\qquad
V_\nu^c:=E_{,16+6,;16+15}=E_{22,31}.
]

Then the conjugate blocks are
[
\boxed{
D_{\rm Dirac}^{\rm conj}
========================

(\widetilde V_u^c)(I_{32}\otimes \overline{Y_u})+(\widetilde V_u^c)^\dagger(I_{32}\otimes Y_u^{T})
+
(\widetilde V_d^c)(I_{32}\otimes \overline{Y_d})+(\widetilde V_d^c)^\dagger(I_{32}\otimes Y_d^{T})
}
]
[
\boxed{
\qquad\quad
+
(\widetilde V_e^c)(I_{32}\otimes \overline{Y_e})+(\widetilde V_e^c)^\dagger(I_{32}\otimes Y_e^{T})
+
(\widetilde V_\nu^c)(I_{32}\otimes \overline{Y_\nu})+(\widetilde V_\nu^c)^\dagger(I_{32}\otimes Y_\nu^{T}).
}
]

This matches the earlier code convention:

* forward conjugate blocks use (\overline{Y}),
* backward conjugate blocks use (Y^T),
  so that self-adjointness is automatic.

### 3) Majorana block (D_{\rm Maj}) (optional seesaw)

Couple (\nu_R) (index (15)) to (\nu_R^c) (index (31)):
[
W_R:=E_{31,15},\qquad \widetilde W_R=W_R\otimes I_N.
]
Then
[
\boxed{
D_{\rm Maj}
===========

# (\widetilde W_R)(I_{32}\otimes M_R)+(\widetilde W_R)^\dagger(I_{32}\otimes M_R^\dagger)

(E_{31,15}\otimes M_R) + (E_{15,31}\otimes M_R^\dagger).
}
]
(If you want Takagi symmetry in a (J_{\rm flav}=K) basis: impose (M_R^T=M_R).)

### Final internal operator

[
\boxed{
D_{\rm int}
===========

D_{\rm Dirac}^{\rm part}
+
D_{\rm Dirac}^{\rm conj}
+
D_{\rm Maj}.
}
]

---

## How this matches the code symbol-for-symbol

Your code primitive was:

* “place an (N\times N) block (B) into SM slot ((i,j))”.

That is exactly the matrix
[
E_{ij}\otimes B.
]

So each line like

```python
put_block(D, i, j, Yu)
```

is precisely the term (E_{ij}\otimes Y_u).

And each “Hermitian completion”

```python
put_block(D, j, i, Yu.conj().T)
```

is precisely the adjoint term (E_{ji}\otimes Y_u^\dagger).

## Invariant assembly of (D_{\rm int}) (no indices)

Fix (H_{\rm int}=H_{\rm SM}\otimes H_{\rm flav}) with (\dim H_{\rm SM}=32), (\dim H_{\rm flav}=N). Let the SM particle space decompose (one generation)
[
H_{\rm SM}^{\rm part}=Q_L\oplus L_L\oplus u_R\oplus d_R\oplus e_R\oplus \nu_R,
]
and take the full internal space with conjugates
[
H_{\rm SM}
==========

H_{\rm SM}^{\rm part}\oplus (H_{\rm SM}^{\rm part})^c.
]
Let (P_X) be the orthogonal projector onto the subspace (X\subset H_{\rm SM}^{\rm part}) (extended by (0) on the orthogonal complement), and similarly (P_{X^c}) onto the conjugate subspace.

### Lifted projectors and textures

[
\widetilde P_X := P_X\otimes I_N,\qquad \widetilde P_{X^c}:=P_{X^c}\otimes I_N,\qquad
\widetilde Y := I_{32}\otimes Y.
]

---

## Partial isometries as canonical channel maps

Choose fixed isometries (they are uniquely determined once a basis in each multiplet is chosen, but the *form* is basis-free):

* (V_u: Q_L \to u_R) identifies the **up-component** of the quark doublet with (u_R) (color-preserving).
* (V_d: Q_L \to d_R) identifies the **down-component** of the quark doublet with (d_R) (color-preserving).
* (V_\nu: L_L \to \nu_R) identifies the **neutrino component** of the lepton doublet with (\nu_R).
* (V_e: L_L \to e_R) identifies the **electron component** of the lepton doublet with (e_R).

These satisfy the partial-isometry relations
[
V_u^\dagger V_u = P_{Q_L^{\uparrow}},\qquad V_u V_u^\dagger = P_{u_R},
]
[
V_d^\dagger V_d = P_{Q_L^{\downarrow}},\qquad V_d V_d^\dagger = P_{d_R},
]
[
V_\nu^\dagger V_\nu = P_{L_L^{\uparrow}},\qquad V_\nu V_\nu^\dagger = P_{\nu_R},
]
[
V_e^\dagger V_e = P_{L_L^{\downarrow}},\qquad V_e V_e^\dagger = P_{e_R},
]
where (Q_L^{\uparrow,\downarrow}) and (L_L^{\uparrow,\downarrow}) denote the SU(2) components inside the doublets.

Lift them to (H_{\rm int}) by (\widetilde V_\bullet := V_\bullet\otimes I_N).

---

## Yukawa textures (Alignment output) as commutant operators

Let
[
Y_u,;Y_d,;Y_e,;Y_\nu \in M_N(\mathbb C),
\qquad
\widetilde Y_s := I_{32}\otimes Y_s,
]
and impose the v4.0 commutant condition (Y_s\in\pi(A)').

---

## Definition: Dirac (LR) part in projector/partial-isometry notation

Define the particle Dirac contribution
[
\boxed{
D_{\rm Dirac}^{\rm part}
========================

\widetilde V_u,\widetilde Y_u
+\widetilde V_u^\dagger,\widetilde Y_u^\dagger
+\widetilde V_d,\widetilde Y_d
+\widetilde V_d^\dagger,\widetilde Y_d^\dagger
+\widetilde V_e,\widetilde Y_e
+\widetilde V_e^\dagger,\widetilde Y_e^\dagger
+\widetilde V_\nu,\widetilde Y_\nu
+\widetilde V_\nu^\dagger,\widetilde Y_\nu^\dagger.
}
]
This is self-adjoint by construction.

---

## Conjugate Dirac part

Let (J_{\rm SM}=U_{\rm SM}\circ K) be the antiunitary real structure on (H_{\rm SM}), and define the conjugate channel maps invariantly by
[
V_s^c := J_{\rm SM},V_s,J_{\rm SM}^{-1},\qquad s\in{u,d,e,\nu}.
]
Lift: (\widetilde V_s^c := V_s^c\otimes I_N).

Then define
[
\boxed{
D_{\rm Dirac}^{\rm conj}
========================

\widetilde V_u^c,\widetilde{\overline{Y_u}}
+(\widetilde V_u^c)^\dagger,\widetilde{Y_u^{T}}
+
\widetilde V_d^c,\widetilde{\overline{Y_d}}
+(\widetilde V_d^c)^\dagger,\widetilde{Y_d^{T}}
+
\widetilde V_e^c,\widetilde{\overline{Y_e}}
+(\widetilde V_e^c)^\dagger,\widetilde{Y_e^{T}}
+
\widetilde V_\nu^c,\widetilde{\overline{Y_\nu}}
+(\widetilde V_\nu^c)^\dagger,\widetilde{Y_\nu^{T}}.
}
]
Again, this is self-adjoint term-by-term.

(If you prefer a basis-free statement: “the conjugate blocks are fixed by the real structure”; the (\overline{\cdot}) and ((\cdot)^T) are the matrix realizations in the chosen (J_{\rm flav}=K) basis.)

---

## Majorana (seesaw) part (optional)

Let (P_{\nu_R}) be the projector onto (\nu_R\subset H_{\rm SM}^{\rm part}) and (P_{\nu_R^c}) onto its conjugate. Choose a partial isometry (W_R:\nu_R\to \nu_R^c) with
[
W_R^\dagger W_R=P_{\nu_R},\qquad W_R W_R^\dagger=P_{\nu_R^c}.
]
Lift (\widetilde W_R=W_R\otimes I_N). Let (M_R\in M_N(\mathbb C)) (Takagi-symmetric if desired: (M_R^T=M_R)) and (\widetilde M_R=I_{32}\otimes M_R). Define
[
\boxed{
D_{\rm Maj}
===========

\widetilde W_R,\widetilde M_R
+\widetilde W_R^\dagger,\widetilde M_R^\dagger.
}
]

---

## Final internal operator

[
\boxed{
D_{\rm int}(Y[\mathcal K])
==========================

D_{\rm Dirac}^{\rm part}
+
D_{\rm Dirac}^{\rm conj}
+
D_{\rm Maj}.
}
]

This is the invariant “paper form”: sums of (lifted) partial isometries times (lifted) flavor textures.

---

# Lemma: index expansion equals block insertion

**Lemma (Index expansion).** Fix the canonical basis ordering of (H_{\rm SM}) (16 particle + 16 conjugate) and choose bases in each multiplet so that each partial isometry (V) is a sum of matrix units:
[
V=\sum_{(i,j)\in\mathcal I(V)} E_{ij}.
]
Then for any flavor matrix (Y\in M_N(\mathbb C)),
[
\widetilde V,\widetilde Y
=========================

# (V\otimes I_N)(I_{32}\otimes Y)

\sum_{(i,j)\in\mathcal I(V)} (E_{ij}\otimes Y),
]
i.e. (\widetilde V,\widetilde Y) is exactly the (32N\times 32N) matrix whose ((i,j)) SM-slot block equals (Y) and all other blocks are zero.

**Corollary (Code correspondence).** The instruction “put (Y) into SM slot ((i,j))” is precisely the term (E_{ij}\otimes Y). Therefore the invariant formula for (D_{\rm int}) expands to the earlier block-sparse assembly code verbatim, with self-adjoint completion coming from adding the adjoint terms.

*Proof.* Immediate from Kronecker product identities and (E_{ij}E_{kl}=\delta_{jk}E_{il}). ∎

---

## Minimal “paper-to-code dictionary” (one paragraph you can paste)

* **Projector lift:** (P\mapsto \widetilde P=P\otimes I_N).
* **Texture lift:** (Y\mapsto \widetilde Y=I_{32}\otimes Y).
* **Channel map:** (V=\sum E_{ij}) (partial isometry).
* **Block insertion:** (E_{ij}\otimes Y) equals “place (Y) in SM slot ((i,j))”.
* **Internal Dirac:** (D_{\rm int}=\sum \widetilde V_s\widetilde Y_s+\widetilde V_s^\dagger \widetilde Y_s^\dagger+\cdots).

## Channel maps (V_u,V_d,V_e,V_\nu,W_R) from projectors and representation data

We work on the **one-generation** internal SM Hilbert space
[
H_{\rm SM}^{\rm part}=Q_L\oplus L_L\oplus u_R\oplus d_R\oplus e_R\oplus \nu_R,
\qquad
H_{\rm SM}=H_{\rm SM}^{\rm part}\oplus (H_{\rm SM}^{\rm part})^c,
]
and then lift to (H_{\rm int}=H_{\rm SM}\otimes H_{\rm flav}) by (\widetilde{(\cdot)}=(\cdot)\otimes I_N).

### 1) Projectors and SU(2) component selection

Let
[
P_{Q_L},,P_{L_L},,P_{u_R},,P_{d_R},,P_{e_R},,P_{\nu_R}
]
be the orthogonal projectors onto the indicated summands of (H_{\rm SM}^{\rm part}) (extended by (0) on the other summands, and by (0) on the conjugate half unless explicitly stated).

Inside the doublets (Q_L\simeq \mathbb C^2\otimes \mathbb C^3) and (L_L\simeq \mathbb C^2), fix the **component projectors**
[
p_\uparrow,;p_\downarrow\in \mathrm{End}(\mathbb C^2),\qquad
p_\uparrow+p_\downarrow=1,\quad p_\uparrow p_\downarrow=0.
]
This is exactly the allowed “SU(2) component selection” datum.

Define the induced subspace projectors
[
P_{Q_L^\uparrow}:=(p_\uparrow\otimes 1_3),P_{Q_L},\qquad
P_{Q_L^\downarrow}:=(p_\downarrow\otimes 1_3),P_{Q_L},
]
[
P_{L_L^\uparrow}:=p_\uparrow P_{L_L},\qquad
P_{L_L^\downarrow}:=p_\downarrow P_{L_L}.
]

### 2) Canonical color-preserving intertwiners

Because
[
Q_L^\uparrow \cong \mathbb C^3 \cong u_R,\qquad
Q_L^\downarrow \cong \mathbb C^3 \cong d_R,
]
there exist **SU(3)-equivariant unitary identifications** (unique up to an overall phase):
[
U_u: Q_L^\uparrow \to u_R,\qquad
U_d: Q_L^\downarrow \to d_R.
]
Likewise, since (L_L^\uparrow\cong \nu_R\cong \mathbb C) and (L_L^\downarrow\cong e_R\cong \mathbb C), fix unitary identifications (phases):
[
U_\nu: L_L^\uparrow \to \nu_R,\qquad
U_e: L_L^\downarrow \to e_R.
]

These (U)’s are the minimal representation-theoretic choices: *color-preserving* on quarks and *component-preserving* on doublets.

### 3) Definition of the partial isometries (V_u,V_d,V_\nu,V_e)

Define operators on (H_{\rm SM}) by extending by zero outside the stated domains/codomains:
[
\boxed{
V_u := P_{u_R},U_u,P_{Q_L^\uparrow},\qquad
V_d := P_{d_R},U_d,P_{Q_L^\downarrow},
}
]
[
\boxed{
V_\nu := P_{\nu_R},U_\nu,P_{L_L^\uparrow},\qquad
V_e := P_{e_R},U_e,P_{L_L^\downarrow}.
}
]

They are **partial isometries** with the desired initial/final projections:
[
V_u^\dagger V_u=P_{Q_L^\uparrow},\qquad V_uV_u^\dagger=P_{u_R},
]
[
V_d^\dagger V_d=P_{Q_L^\downarrow},\qquad V_dV_d^\dagger=P_{d_R},
]
[
V_\nu^\dagger V_\nu=P_{L_L^\uparrow},\qquad V_\nu V_\nu^\dagger=P_{\nu_R},
]
[
V_e^\dagger V_e=P_{L_L^\downarrow},\qquad V_e V_e^\dagger=P_{e_R}.
]

### 4) Conjugate channel maps (V_s^c)

Let the real structure be (J_{\rm SM}=U_{\rm SM}\circ K). Define conjugate maps **representation-theoretically** by
[
\boxed{
V_s^c := J_{\rm SM},V_s,J_{\rm SM}^{-1}\quad (s\in{u,d,e,\nu}),
}
]
so their initial/final projections are automatically the conjugate projectors:
[
(V_s^c)^\dagger V_s^c=P_{(\text{domain of }V_s)^c},\qquad
V_s^c(V_s^c)^\dagger=P_{(\text{codomain of }V_s)^c}.
]

### 5) Majorana channel map (W_R) (seesaw)

The Majorana coupling is a map (\nu_R\to \nu_R^c). Define it using the **unitary part** of the real structure:
[
\boxed{
W_R := P_{\nu_R^c},U_{\rm SM},P_{\nu_R}.
}
]
Then (W_R) is a partial isometry with
[
W_R^\dagger W_R=P_{\nu_R},\qquad W_R W_R^\dagger=P_{\nu_R^c}.
]

This is the clean representation-theoretic content of “(\nu_R) pairs with its charge conjugate under (J)”.

### 6) Lift to (H_{\rm SM}\otimes H_{\rm flav})

[
\boxed{
\widetilde V_s := V_s\otimes I_N,\qquad
\widetilde V_s^c := V_s^c\otimes I_N,\qquad
\widetilde W_R := W_R\otimes I_N.
}
]
And textures lift as (\widetilde Y:=I_{32}\otimes Y).

With these definitions, the invariant internal operator assembly from the previous step,
[
D_{\rm int}= \sum_{s\in{u,d,e,\nu}}\big(\widetilde V_s\widetilde Y_s+\widetilde V_s^\dagger\widetilde Y_s^\dagger\big)
+\sum_{s\in{u,d,e,\nu}}\big(\widetilde V_s^c\widetilde{\overline{Y_s}}+(\widetilde V_s^c)^\dagger\widetilde{Y_s^T}\big)
+\big(\widetilde W_R\widetilde M_R+\widetilde W_R^\dagger\widetilde M_R^\dagger\big),
]
is now entirely **representation-theoretic**: only projectors, SU(2) component selection ((p_\uparrow,p_\downarrow)), color-preserving intertwiners, and the real structure.

## Lemma: index expansion of the channel maps equals block-slot insertion

### Setup

Fix a concrete orthonormal basis of the one-generation SM particle space
[
H_{\rm SM}^{\rm part}=Q_L\oplus L_L\oplus u_R\oplus d_R\oplus e_R\oplus \nu_R
]
compatible with the decomposition, the SU(2) component splitting ((\uparrow,\downarrow)) on (Q_L,L_L), and the color basis on (\mathbb C^3). Extend it to (H_{\rm SM}=H_{\rm SM}^{\rm part}\oplus (H_{\rm SM}^{\rm part})^c) by taking the conjugate basis.

Let ({e_i}*{i=0}^{31}) be the resulting ordered basis of (H*{\rm SM}) and let (E_{ij}\in M_{32}(\mathbb C)) be the matrix units:
[
E_{ij} e_j = e_i,\qquad E_{ij}^\dagger=E_{ji}.
]

Lift to (H_{\rm int}=H_{\rm SM}\otimes H_{\rm flav}) via the Kronecker product:
[
\widetilde E_{ij}:=E_{ij}\otimes I_N,\qquad \widetilde Y:=I_{32}\otimes Y.
]

### Lemma 1 (Channel maps are sums of matrix units)

With the basis above, each partial isometry channel map (V\in{V_u,V_d,V_e,V_\nu,W_R}) is a finite sum of matrix units:
[
\boxed{
V=\sum_{(i,j)\in \mathcal I(V)} E_{ij},
}
]
where (\mathcal I(V)) is the set of basis-pairs ((i,j)) for which (V e_j=e_i). In particular:

* (V_u) is the color-preserving identification (Q_L^\uparrow\to u_R), hence (\mathcal I(V_u)) consists of the three color-matched pairs.
* (V_d) is the color-preserving identification (Q_L^\downarrow\to d_R), hence (\mathcal I(V_d)) consists of the three color-matched pairs.
* (V_\nu) and (V_e) each contribute a single pair.
* (W_R) contributes the single pair (\nu_R\to \nu_R^c).

*Proof.* By construction, each (V) is a partial isometry mapping an orthonormal basis of its initial subspace bijectively onto an orthonormal basis of its final subspace (color-preserving and SU(2)-component preserving). Therefore (V) acts as a permutation (possibly with phases; fix phases to (1) for the canonical choice) on the chosen basis vectors in its domain, and vanishes on the orthogonal complement. Any such operator is a finite sum of matrix units (E_{ij}). ∎

### Lemma 2 (Index expansion equals ((i,j))-slot block insertion)

Let (V=\sum_{(i,j)\in\mathcal I(V)}E_{ij}). Then for any flavor matrix (Y\in M_N(\mathbb C)),
[
\boxed{
(\widetilde V)(\widetilde Y)
=(V\otimes I_N)(I_{32}\otimes Y)
=\sum_{(i,j)\in\mathcal I(V)} (E_{ij}\otimes Y).
}
]
Thus (\widetilde V,\widetilde Y) is exactly the (32N\times 32N) matrix whose ((i,j)) SM-slot block equals (Y) for each ((i,j)\in\mathcal I(V)), and is zero elsewhere.

*Proof.* Use bilinearity and the Kronecker identity
[
(E_{ij}\otimes I_N)(I_{32}\otimes Y)=E_{ij}\otimes Y,
]
then sum over (\mathcal I(V)). ∎

### Corollary (Code-level guarantee)

The earlier instruction “insert the (N\times N) block (Y) into SM slot ((i,j))” is precisely the term (E_{ij}\otimes Y). Therefore, the invariant assembly
[
D_{\rm int}=\sum_s\big(\widetilde V_s\widetilde Y_s+\widetilde V_s^\dagger\widetilde Y_s^\dagger\big)
+\sum_s\big(\widetilde V_s^c\widetilde{\overline{Y_s}}+(\widetilde V_s^c)^\dagger\widetilde{Y_s^T}\big)
+\big(\widetilde W_R\widetilde M_R+\widetilde W_R^\dagger\widetilde M_R^\dagger\big)
]
expands exactly into the block-sparse ((i,j))-slot insertion implementation, term-for-term, once the concrete basis ordering is fixed. ∎

## One-line “seal” sentence for the paper

“Upon fixing a basis compatible with the SM multiplet decomposition, the channel maps (V_u,V_d,V_e,V_\nu,W_R) become sums of matrix units (E_{ij}), and the lifted operators (E_{ij}\otimes Y) coincide with the ((i,j))-block insertion rule used in the numerical implementation; hence the code is the literal index expansion of the invariant projector/partial-isometry construction.”
## Yes — Origin Mixing should be updated for v4.0

Your current “Origin of Mixing” text is **v3.x-faithful** (mixing explained via emergent block-dimension + closure rule), but in v4.0 the *production object* is the **even internal operator**
[
D_{\rm int}(Y[\mathcal K])\in M_{32N}(\mathbb C),
\qquad
D = D_{\rm geom}\otimes I_{32N}+\gamma_{\rm geom}\otimes D_{\rm int},
\qquad
{\gamma_{\rm SM}\otimes I_N,\ D_{\rm int}}=0,
]
with textures living in the **commutant**
[
Y_s\in \pi(A)'\quad(\text{acting on }H_{\rm flav}),
]
so “block diagnosis” becomes **diagnostic telemetry**, not a structural axiom.

## 3. Alignment Origin of Mixing (v4.0)

We work in the Alignment Spectral Triple v4.0, where flavor is a pure multiplicity space (H_{\rm flav}\cong\mathbb C^N) and the represented algebra acts trivially on flavor. The sector textures (Y_u,Y_d,Y_e,Y_\nu\in M_N(\mathbb C)) are generated by the alignment kernel pipeline (misalignment flow, (C_{360}) projection, triadic compression), and then lifted to the commutant:
[
\widetilde Y_s := I_{32}\otimes Y_s \in \pi(A)' \subset \mathcal B(H_{\rm SM}\otimes H_{\rm flav}).
]

The finite internal operator is assembled as a projector/partial-isometry sum (Dirac channels plus optional Majorana channel):
[
D_{\rm int}(Y[\mathcal K])
==========================

\sum_{s\in{u,d,e,\nu}}
\big(\widetilde V_s,\widetilde Y_s+\widetilde V_s^\dagger,\widetilde Y_s^\dagger\big)
;+;
\big(\widetilde W_R,\widetilde M_R+\widetilde W_R^\dagger,\widetilde M_R^\dagger\big),
]
with (\widetilde V_s:=V_s\otimes I_N) and (\widetilde W_R:=W_R\otimes I_N). This construction is **odd** for (\gamma_{\rm SM}\otimes I_N), hence compatible with full evenness of the product triple.

### Proposition 1 (Dirac vs Majorana is a channel choice, not a block heuristic)

A sector is Dirac-type if its internal channel contains only LR Dirac couplings (\widetilde V_s\widetilde Y_s) (and conjugates). A sector becomes Majorana-type precisely when a right-handed Majorana channel (\widetilde W_R\widetilde M_R) is present, yielding an effective light mass operator after integrating out heavy modes:
[
m_\nu ;\sim; -,v^2, Y_\nu, M_R^{-1}, Y_\nu^{T}
\quad (\text{Takagi-diagonalized}).
]
Thus “closure failure (\Rightarrow) seesaw” is retained only as *diagnostic intuition* about when the kernel tends to generate structures that demand a Majorana channel, not as the structural definition.

### Proposition 2 (Mixing origin as misalignment of sector diagonalizations)

Mixing angles arise from mismatches between the unitary (or Takagi) diagonalizations of **distinct sector textures** extracted from (D_{\rm int}) by projector-defined channels:

* Quarks:
  [
  V_{\rm CKM} = U_u^\dagger U_d,
  ]
  where (U_u,U_d) diagonalize the Hermitian combinations (Y_u Y_u^\dagger) and (Y_d Y_d^\dagger).

* Leptons:
  [
  U_{\rm PMNS} = U_e^\dagger U_\nu,
  ]
  where (U_e) diagonalizes (Y_eY_e^\dagger) and (U_\nu) Takagi-diagonalizes (m_\nu).

In alignment terms, the kernel pipeline shapes the relative *phase geometry* of the sector textures. Small mixing corresponds to near-simultaneous diagonalizability; large mixing corresponds to near-degenerate planes in which rotations cost little spectral weight.

### Diagnostic note (retain the block language, but demote it)

The empirical “split vs (1\oplus 2) vs (3D)” regimes observed by harmonic block diagnosis of a compressed operator (Y_{\rm raw}) remain valuable telemetry:

* split (\Rightarrow) rigid eigenframe (\Rightarrow) small mixing,
* (1\oplus 2) (\Rightarrow) one stable axis + one soft plane (\Rightarrow) CKM-small / PMNS-large pattern,
* (3D) (\Rightarrow) fully soft frame (\Rightarrow) large PMNS-like mixing.

But in v4.0 these regimes are interpreted as **approximate degeneracy structure of the sector textures** (and therefore of the induced diagonalization frames), not as an axiom-level “closure rule” that triggers new terms in (D). The seesaw/Majorana structure is carried explicitly by the (\widetilde W_R\widetilde M_R) channel inside (D_{\rm int}).

## Alignment Origin of Mixing v4.0

### Section 1–10 theorem chain (referee-friendly)

## Section 1 — Definitions

### Definition 1.1 (v4.0 sector textures and commutant lift)

Let (H_{\rm SM}) be the fixed 32-dimensional one-generation SM finite space (16 particle (+) 16 conjugate) with grading (\gamma_{\rm SM}), and let (H_{\rm flav}\cong\mathbb C^N) be the flavor multiplicity space.
For each sector (s\in{u,d,e,\nu}) fix a texture
[
Y_s \in M_N(\mathbb C),
]
and (optionally) a Majorana texture
[
M_R\in M_N(\mathbb C)\quad\text{(complex symmetric preferred)}.
]
Assume the **commutant condition**
[
\widetilde Y_s := I_{32}\otimes Y_s \in \pi(A)'\subset \mathcal B(H_{\rm SM}\otimes H_{\rm flav}),
\qquad
\widetilde M_R := I_{32}\otimes M_R \in \pi(A)'.
]

### Definition 1.2 (Channel maps and internal operator)

Fix partial isometries (channel maps) (V_u,V_d,V_e,V_\nu,W_R\in\mathcal B(H_{\rm SM})) (color-preserving, SU(2)-component selecting) and lift them:
[
\widetilde V_s:=V_s\otimes I_N,\qquad \widetilde W_R:=W_R\otimes I_N.
]
Define the internal finite operator
[
D_{\rm int}(Y[\mathcal K])
==========================

\sum_{s\in{u,d,e,\nu}}
\big(\widetilde V_s\widetilde Y_s+\widetilde V_s^\dagger\widetilde Y_s^\dagger\big)
;+;
\big(\widetilde W_R\widetilde M_R+\widetilde W_R^\dagger\widetilde M_R^\dagger\big),
]
with the v4.0 **oddness** requirement
[
{\gamma_{\rm SM}\otimes I_N,\ D_{\rm int}}=0.
]

## Section 2 — Lemma: sector extraction is projector-defined

### Lemma 2.1 (Projector extraction equals channel insertion)

Let (P_{sL},P_{sR}) be the SM projectors onto the left/right multiplets used by the channel (V_s), and lift (\widetilde P:=P\otimes I_N). Then the LR block of (D_{\rm int}) on sector (s) is
[
\widetilde P_{sR},D_{\rm int},\widetilde P_{sL}
===============================================

\widetilde V_s,\widetilde Y_s,
\qquad
\widetilde P_{sL},D_{\rm int},\widetilde P_{sR}
===============================================

\widetilde V_s^\dagger,\widetilde Y_s^\dagger.
]
**Proof.** Expand (D_{\rm int}) as a sum over channels; all non-(s) terms are annihilated by the sector projectors, while the (s)-term survives and factorizes by construction. (\square)

## Section 3 — Definitions: diagonalization data

### Definition 3.1 (Dirac-sector Hermitian forms)

For (s\in{u,d,e}) define
[
H_s := Y_s Y_s^\dagger \in M_N(\mathbb C),\qquad H_s=H_s^\dagger\ge 0,
]
and let (U_s\in U(N)) diagonalize it:
[
U_s^\dagger H_s U_s = \mathrm{diag}(m_{s,1}^2,\dots,m_{s,N}^2).
]

### Definition 3.2 (Majorana neutrino form)

Given (Y_\nu) and (M_R) (invertible on its support), define the effective light Majorana mass matrix
[
m_\nu := -v^2,Y_\nu,M_R^{-1},Y_\nu^{T},
\qquad m_\nu = m_\nu^{T},
]
and Takagi-diagonalize:
[
U_\nu^{T} m_\nu U_\nu = \mathrm{diag}(m_1,m_2,m_3),\qquad m_i\ge 0.
]

## Section 4 — Theorem: mixing matrices are mismatch of diagonalizers

### Theorem 4.1 (CKM and PMNS as relative eigenframes)

Define
[
V_{\rm CKM}:=U_u^\dagger U_d,\qquad U_{\rm PMNS}:=U_e^\dagger U_\nu.
]
Then all quark/lepton mixing angles arise from the **relative orientation** of the diagonalizing eigenframes produced by the sector textures extracted from (D_{\rm int}).
**Proof.** Standard: left-handed rotations that diagonalize the mass (or mass-squared) operators differ between sectors; the charged-current coupling compares these bases, yielding the stated products. (\square)

## Section 5 — Definition: block regimes and eigenvalue gaps

### Definition 5.1 (Spectral clusters and gaps)

Let (H) be Hermitian with eigenvalues (\lambda_1\le \cdots\le \lambda_N). A **cluster** (I\subset{1,\dots,N}) is a set of indices whose eigenvalues lie in an interval ([a,b]). The **external gap** of (I) is
[
\Delta_I := \min{\lambda_{j}-b: j\notin I,\ \lambda_j>b}\ \wedge\ \min{a-\lambda_{j}: j\notin I,\ \lambda_j<a}.
]
A “split / (1\oplus 2) / (3D)” regime is shorthand for the cluster pattern of ({\lambda_i}) (one-by-one vs one isolated + near-degenerate plane vs near-degenerate triple).

## Section 6 — Lemma: exact degeneracy gives free rotations

### Lemma 6.1 (Rotation freedom inside a degenerate cluster)

If (H) has an exactly degenerate cluster (I) (all (\lambda_i=\lambda) for (i\in I)), then for any unitary (W) acting only on (\mathrm{span}{e_i:i\in I}),
[
W^\dagger H W = H.
]
So the eigenbasis in that subspace is not geometrically fixed, and large mixing can occur without changing masses.
**Proof.** On the degenerate subspace (H=\lambda I), hence it commutes with all unitaries supported there. (\square)

## Section 7 — Quantitative lemma: “block-regime ⇒ rotation freedom” via gaps

### Lemma 7.1 (Gap controls subspace stability: Davis–Kahan bound)

Let (H) be Hermitian and let (I) be a cluster with spectral projector (P_I). Let (H' = H + E) with projector (P'_I) onto the corresponding cluster of (H') (same index set size). If the cluster remains separated and (\Delta_I>0), then
[
|P_I - P'_I|\ \le\ \frac{|E|}{\Delta_I}.
]
In particular, when (\Delta_I) is small (near-degeneracy), the invariant subspace can rotate by an (O(1)) angle under a perturbation (E) of size comparable to (\Delta_I).
**Proof (two lines).** The Davis–Kahan sin(\Theta) theorem gives (|\sin\Theta|\le |E|/\Delta_I); (|P_I-P'_I|=|\sin\Theta|). (\square)

## Section 8 — Theorem: v4.0 origin of CKM–PMNS disparity

### Theorem 8.1 (Rigidity vs softness produces small vs large mixing)

Assume quark-sector Hermitian forms (H_u,H_d) have **large spectral gaps** separating their eigenvalues (split regime), while at least one lepton-sector operator (typically the neutrino effective form) has a **small gap** supporting a near-degenerate plane or triple. Then generically:

* (U_u) and (U_d) are **rigid** (eigenframes stable), hence (V_{\rm CKM}=U_u^\dagger U_d) has small angles.
* (U_\nu) (and/or (U_e)) is **soft** on the near-degenerate cluster, hence (U_{\rm PMNS}=U_e^\dagger U_\nu) can realize large angles without destabilizing the spectrum.

**Proof.** By Lemma 7.1, a cluster with large (\Delta) fixes its invariant subspace up to (O(|E|/\Delta)), suppressing rotation freedom; a cluster with small (\Delta) allows (O(1)) subspace rotations. Mixing is exactly relative eigenframe orientation (Theorem 4.1). (\square)

## Section 9 — Corollary: how the kernel pipeline enters in v4.0

### Corollary 9.1 (Kernel pipeline shapes gaps, gaps shape mixing)

The alignment kernel pipeline (flow, (C_{360}) projection, compression) determines the sector textures (Y_s) and hence the spectra and gaps of (H_s=Y_sY_s^\dagger) (and (m_\nu)). Therefore, the observed “split / (1\oplus 2) / (3D)” regimes remain meaningful as **diagnostic telemetry**: they predict whether the sector eigenframes are rigid or soft, and thus whether small or large mixing is natural.

## Section 10 — Production diagnostics (gap-gates for scans)

### Gate 10.1 (Gap telemetry for rigidity/softness)

For each sector (s), compute ordered eigenvalues of (H_s) (or (m_\nu m_\nu^\dagger)) and define gap ratios, e.g.
[
\rho_s := \frac{\min(\lambda_{i+1}-\lambda_i)}{\max(\lambda_i)}.
]
Small (\rho_s) flags rotation freedom (soft plane), large (\rho_s) flags rigidity.

### Gate 10.2 (Frame sensitivity test)

Apply a controlled perturbation (E) (e.g. from scan-to-scan parameter drift) and record
[
|P_I - P'_I|\quad \text{and compare to}\quad |E|/\Delta_I
]
as a numerical validation of Lemma 7.1 in your engine.

#!/usr/bin/env python3
"""
ncg_tests.py

Reusable Noncommutative Geometry test harness for finite spectral triples.

Provides:
  - AlgebraElement dataclass
  - First-order condition checker
  - Zero-order condition checker
  - Grading & reality checker
  - A convenience function to run the full test suite and pretty-print results
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional
import numpy as np
from scipy.linalg import kron

from v3.spectral import SMState, build_internal_algebra_ops

# =========================
# v4.0 EVEN PRODUCT TOGGLES
# =========================
USE_V4_EVEN_PRODUCT = True   # set False to keep your current "odd product" behavior
EPS_GAMMA = 1e-12            # tolerance for grading anti/commutation sanity
TOPK_VIOLATIONS = 10         # show worst offenders for debug


def hermitian_basis(n: int):
    """
    Full Hermitian basis {E_k} for n x n Hermitian matrices.
    """
    basis = []
    # diagonal elements
    for i in range(n):
        M = np.zeros((n, n), dtype=complex)
        M[i, i] = 1.0
        basis.append(M)
    # off-diagonal: real-symmetric and imaginary-antisymmetric parts
    for i in range(n):
        for j in range(i + 1, n):
            M_re = np.zeros((n, n), dtype=complex)
            M_im = np.zeros((n, n), dtype=complex)
            M_re[i, j] = 1.0
            M_re[j, i] = 1.0
            M_im[i, j] = 1.0j
            M_im[j, i] = -1.0j
            basis.append(M_re)
            basis.append(M_im)
    return basis


def filter_basis_gamma_J(basis_D, gamma: np.ndarray, S: np.ndarray,
                         tol: float = 1e-12):
    """
    Filter Hermitian basis elements to those that are:
      - gamma-odd:  {gamma, E} = 0
      - J-even:     S E* S^T - E = 0
    Returns filtered_basis.
    """
    filtered = []
    for E in basis_D:
        Cg = gamma @ E + E @ gamma
        CJ = S @ E.conj() @ S.T - E
        if (np.linalg.norm(Cg, ord="fro") < tol and
            np.linalg.norm(CJ, ord="fro") < tol):
            filtered.append(E)
    return filtered

# ============================================================
# 1. Basis for 1-generation internal Hilbert space
# ============================================================

@dataclass
class BasisState:
    name: str          # e.g. "nu_L", "e_L", "u_L_r", "d_R_b", ...
    chirality: str     # "L" or "R"
    particle: bool     # True = particle, False = antiparticle
    is_quark: bool     # True for quarks, False for leptons
    # before: color: Optional[str]
    color: Optional[Tuple[str, str]]  # (left_color, right_color) or None
    generation: int    # 1 (for now)


def build_sm_basis_1gen(include_nu_R: bool = True) -> Tuple[List[BasisState], Dict[str, int]]:
    """
    One generation of SM fermions + antiparticles.

    Particle sector:
      Leptons:
        L_L: (nu_L, e_L)
        R :  (nu_R (optional), e_R)

      Quarks (color = r,g,b):
        Q_L: (u_L^r, d_L^r, u_L^g, d_L^g, u_L^b, d_L^b)
        R :  (u_R^r, u_R^g, u_R^b, d_R^r, d_R^g, d_R^b)

    Antiparticle sector: charge conjugates of all above, with same chirality label
    (J will swap particle ↔ antiparticle sectors).
    """
    basis: List[BasisState] = []

    def add(name, chirality, particle, is_quark, color=None):
        basis.append(BasisState(
            name=name, chirality=chirality, particle=particle,
            is_quark=is_quark, color=color, generation=1
        ))

    # --- Particle Leptons ---
    add("nu_L", "L", True, False)
    add("e_L",  "L", True, False)
    if include_nu_R:
        add("nu_R", "R", True, False)
    add("e_R",  "R", True, False)

    # --- Particle Quarks as color bimodules: 3x3 = 9 states each species ---
    colors = ["r", "g", "b"]

    def add_quark(species: str, chir: str, particle: bool):
        for cl in colors:
            for cr in colors:
                add(f"{species}_{chir}_{cl}{cr}", chir, particle, True, color=(cl, cr))

    # Q_L doublets split by species; we keep species labels explicit
    add_quark("u", "L", True)
    add_quark("d", "L", True)
    add_quark("u", "R", True)
    add_quark("d", "R", True)


    # At this point, count particle states:
    # leptons: 2L + (1 or 2)R = 3 or 4
    # quarks   6L + 3R + 3R = 12
    # total particle states = 15 or 16
    n_particle = len(basis)

    # --- Antiparticles: one conjugate state per particle state ---
    for bs in list(basis):
        add(bs.name + "_c", bs.chirality, False, bs.is_quark, bs.color)

    idx: Dict[str, int] = {bs.name: i for i, bs in enumerate(basis)}
    return basis, idx


# ============================================================
# 2. Gamma_F and J_F (swap matrix)
# ============================================================

def build_gamma_F_SM(basis: List[BasisState]) -> np.ndarray:
    """
    gamma_F = -1 on left-handed states, +1 on right-handed states,
    for both particles and antiparticles.
    """
    dimH = len(basis)
    gamma = np.zeros((dimH, dimH), dtype=complex)
    for i, bs in enumerate(basis):
        sgn = -1.0 if bs.chirality == "L" else +1.0
        gamma[i, i] = sgn
    return gamma


def build_swap_particle_antiparticle(basis: List[BasisState], idx: Dict[str, int]) -> np.ndarray:
    """
    Build S such that:
      S |particle> = |particle_c>
      S |particle_c> = |particle>
    i.e. S^2 = I.
    """
    dimH = len(basis)
    S = np.zeros((dimH, dimH), dtype=complex)

    for bs in basis:
        if bs.particle:
            i = idx[bs.name]
            j = idx[bs.name + "_c"]
            S[i, j] = 1.0
            S[j, i] = 1.0

    return S


def J_action_SM(M: np.ndarray, S: np.ndarray, phase: complex = 1.0) -> np.ndarray:
    """
    Real structure on matrices:
        J M J^{-1} = phase * S M^* S^T
    """
    return phase * (S @ M.conj() @ S.T)


# ============================================================
# 3. Representation of A_F = C ⊕ H ⊕ M_3(C)
# ============================================================

@dataclass
class SMAlgebraElement:
    label: str
    op: np.ndarray


def quaternion_to_2x2(a: complex, b: complex) -> np.ndarray:
    """
    Represent q = a + b j as 2x2 complex matrix:
      [ a   b ]
      [-b* a* ]
    For our purposes, we only need a toy faithful representation.
    """
    a = complex(a)
    b = complex(b)
    return np.array([[a, b], [-b.conjugate(), a.conjugate()]], dtype=complex)

def gell_mann():
    """Return a small SU(3) generator set (not normalized)."""
    zero = 0.0 + 0.0j
    one  = 1.0 + 0.0j
    i    = 1.0j

    lam1 = np.array([[0,1,0],[1,0,0],[0,0,0]], dtype=complex)
    lam2 = np.array([[0,-i,0],[i,0,0],[0,0,0]], dtype=complex)
    lam3 = np.array([[1,0,0],[0,-1,0],[0,0,0]], dtype=complex)
    lam4 = np.array([[0,0,1],[0,0,0],[1,0,0]], dtype=complex)
    lam5 = np.array([[0,0,-i],[0,0,0],[i,0,0]], dtype=complex)
    lam6 = np.array([[0,0,0],[0,0,1],[0,1,0]], dtype=complex)
    lam7 = np.array([[0,0,0],[0,0,-i],[0,i,0]], dtype=complex)
    lam8 = np.array([[1,0,0],[0,1,0],[0,0,-2]], dtype=complex)

    return {"lam1": lam1, "lam2": lam2, "lam3": lam3, "lam6": lam6, "lam8": lam8}


def rep_A_SM(
    lam: complex,
    q: np.ndarray,
    m3: np.ndarray,
    basis: List[BasisState],
    idx: Dict[str, int],
) -> np.ndarray:
    """
    Represent (lam ∈ C, q ∈ H≈2x2, m3 ∈ M3(C)) on H_F.

    Simplified rules:
      - lambda (C-part) acts as scalar on all states.
      - q (H-part) acts non-trivially on SU(2)_L doublets:
          (nu_L, e_L) and (u_L^c, d_L^c) for each color,
        and acts trivially on SU(2) singlets (all R states).
      - m3 acts on color indices of quarks (3-dim rep), trivially on leptons.

    Antiparticle sector: use same representation (J takes care of conjugation).
    """
    dimH = len(basis)
    A = np.zeros((dimH, dimH), dtype=complex)

    # C-part: global scalar
    A += lam * np.eye(dimH, dtype=complex)

    # H-part: SU(2)_L doublets
    # (nu_L, e_L)
    if "nu_L" in idx and "e_L" in idx:
        i_nu = idx["nu_L"]
        i_e  = idx["e_L"]
        # insert q on the (nu_L, e_L) subspace
        A[np.ix_([i_nu, i_e], [i_nu, i_e])] += q

    # M3-part: faithful bimodule color action on quarks
    # Particle quarks:  L(m3) = m3 ⊗ I3  on (cl, cr)
    # Antiparticle quarks: R(m3) = I3 ⊗ m3^T on (cl, cr)
    I3 = np.eye(3, dtype=complex)
    cpos = {"r": 0, "g": 1, "b": 2}

    def block_indices(species: str, chir: str, conj: bool) -> Optional[List[int]]:
        suffix = "_c" if conj else ""
        names = []
        for cl in ["r", "g", "b"]:
            for cr in ["r", "g", "b"]:
                names.append(f"{species}_{chir}_{cl}{cr}{suffix}")
        if all(n in idx for n in names):
            return [idx[n] for n in names]
        return None

    def add_color_block(species: str, chir: str, conj: bool):
        ii = block_indices(species, chir, conj)
        if ii is None:
            return
        if not conj:
            B = np.kron(m3, I3)  # left action on particles
        else:
            B = np.kron(I3, m3.T)  # right action on antiparticles
        A[np.ix_(ii, ii)] += B

    for conj in (False, True):
        for chir in ("L", "R"):
            add_color_block("u", chir, conj)
            add_color_block("d", chir, conj)

    return A


def build_SM_algebra_generators(
    basis: List[BasisState],
    idx: Dict[str, int],
) -> List[SMAlgebraElement]:
    """
    Upgraded generating set for A_F = C ⊕ H ⊕ M_3(C), represented on H_F.

    What’s improved vs your draft:
      - Uses *true* direct-sum generators (separates C, H, M3 identities)
      - Adds a full Hermitian SU(2) generator set (Pauli σ1,σ2,σ3) on doublets
      - Adds several non-diagonal SU(3) (Gell-Mann) generators to force color mixing
      - Ensures every returned operator is Hermitian numerically
    """
    # Zero elements for the summands
    q0  = np.zeros((2, 2), dtype=complex)
    m30 = np.zeros((3, 3), dtype=complex)

    # SU(2) (Hermitian) generators acting on L-doublets
    qI = np.eye(2, dtype=complex)
    sigma1 = np.array([[0, 1], [1, 0]], dtype=complex)
    sigma2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sigma3 = np.array([[1, 0], [0, -1]], dtype=complex)

    # SU(3) color generators (Hermitian) acting on quark color triplets
    gm = gell_mann()
    m3I = np.eye(3, dtype=complex)

    # Pick a small but genuinely noncommuting / non-diagonal subset
    color_gens: list[tuple[str, np.ndarray]] = []
    for key in ("lam1", "lam2", "lam3", "lam6", "lam8"):
        if key not in gm:
            raise KeyError(f"gell_mann() did not provide '{key}'. Available: {list(gm.keys())}")
        color_gens.append((f"color_{key}", gm[key]))

    # Build generators as (label, lambda, q, m3)
    # Note: these are *direct-sum* generators; do NOT collapse them into one “identity triple”.
    gen_specs: list[tuple[str, complex, np.ndarray, np.ndarray]] = [
        # C-summand (scalar on all states)
        ("C_1", 1.0 + 0j, q0,  m30),

        # H-summand (acts on SU(2)_L doublets)
        ("H_1",        0.0 + 0j, qI,     m30),
        ("H_sigma1",   0.0 + 0j, sigma1, m30),
        ("H_sigma2",   0.0 + 0j, sigma2, m30),
        ("H_sigma3",   0.0 + 0j, sigma3, m30),

        # M3-summand (acts on color triplets)
        ("M3_1", 0.0 + 0j, q0,  m3I),
    ]

    # Add color generators
    for lab, m3 in color_gens:
        gen_specs.append((lab, 0.0 + 0j, q0, m3))

    ops: List[SMAlgebraElement] = []

    for label, lam, q, m3 in gen_specs:
        A = rep_A_SM(lam, q, m3, basis, idx)

        # Enforce Hermiticity numerically (your tests assume *-structure sanity)
        A = 0.5 * (A + A.conj().T)

        ops.append(SMAlgebraElement(label, A))

    # Convenience: include the literal identity on H_F (useful for debugging)
    # (This is redundant with C_1 in your simplified rep, but harmless and explicit.)
    dimH = len(basis)
    ops.append(SMAlgebraElement("I_HF", np.eye(dimH, dtype=complex)))

    return ops



# ============================================================
# 4. Dirac operator D_F for 1 generation
# ============================================================

def build_DF_SM_1gen(
    Y_e: complex,
    Y_nu: complex,
    Y_u: complex,
    Y_d: complex,
    basis: List[BasisState],
    idx: Dict[str, int],
    include_nu_R: bool = True,
) -> np.ndarray:
    """
    Very minimal 1-generation Dirac operator:

    - Acts only between particle L and R in each sector using Yukawa couplings:
        D_F |e_L>  ~ Y_e |e_R>
        D_F |nu_L> ~ Y_nu|nu_R> (if present)
        D_F |u_L_c> ~ Y_u |u_R_c>
        D_F |d_L_c> ~ Y_d |d_R_c>
    - Antiparticle block mirrors the same structure.

    This is just enough structure to let you:
      - plug in canonical SM Yukawas,
      - plug in aligned Yukawas from your pipeline,
      - and run order tests against the SM-like algebra above.
    """
    dimH = len(basis)
    D = np.zeros((dimH, dimH), dtype=complex)

    def couple(L_name: str, R_name: str, Y: complex):
        if L_name in idx and R_name in idx:
            iL = idx[L_name]
            iR = idx[R_name]
            D[iL, iR] = Y.conjugate()
            D[iR, iL] = Y

    # --- Particle sector couplings ---
    # Leptons
    if include_nu_R and "nu_L" in idx and "nu_R" in idx:
        couple("nu_L", "nu_R", Y_nu)
    couple("e_L", "e_R", Y_e)

    # Quarks (3x3 color bimodule indices)
    for cl in ["r", "g", "b"]:
        for cr in ["r", "g", "b"]:
            couple(f"u_L_{cl}{cr}", f"u_R_{cl}{cr}", Y_u)
            couple(f"d_L_{cl}{cr}", f"d_R_{cl}{cr}", Y_d)

    # Antiparticles mirror
    for cl in ["r", "g", "b"]:
        for cr in ["r", "g", "b"]:
            couple(f"u_L_{cl}{cr}_c", f"u_R_{cl}{cr}_c", Y_u)
            couple(f"d_L_{cl}{cr}_c", f"d_R_{cl}{cr}_c", Y_d)

    return D





# ===========================================
# Basic helpers
# ===========================================

def assert_square(M: np.ndarray, name: str = "matrix") -> None:
    if M.shape[0] != M.shape[1]:
        raise ValueError(f"{name} must be square, got {M.shape}.")


def assert_even_dim(M: np.ndarray, name: str = "matrix") -> None:
    assert_square(M, name)
    n = M.shape[0]
    if n % 2 != 0:
        raise ValueError(f"{name} must have even dimension, got {n}.")


def frob_norm(M: np.ndarray) -> float:
    return float(np.linalg.norm(M, ord="fro"))


# ===========================================
# Data structures
# ===========================================

@dataclass
class AlgebraElement:
    label: str
    op: np.ndarray


@dataclass
class FirstOrderResult:
    max_norm: float
    worst_pair: Optional[Tuple[str, str]]
    good_pairs: List[Tuple[str, str, float]]


@dataclass
class ZeroOrderResult:
    max_norm: float
    worst_pair: Optional[Tuple[str, str]]
    bad_pairs: List[Tuple[str, str, float]]


@dataclass
class GradingRealityResult:
    anticom_norm: float          # || {gamma_F, D_F} ||_F
    max_comm_gamma: float        # max ||[gamma_F, a]||_F over a in A
    J2_deviation: float          # ||J^2 - I||_F
    norm_plus: float             # ||J D_F J^-1 - D_F||_F
    norm_minus: float            # ||J D_F J^-1 + D_F||_F
    df_norm: float               # ||D_F||_F
    rel_plus: float              # norm_plus / ||D_F||_F
    rel_minus: float             # norm_minus / ||D_F||_F
    ko_sign: Optional[int]       # +1, -1, or None

# ===========================================
# Real structure & grading builders
# ===========================================

def build_swap_LR(dim_left: int) -> np.ndarray:
    """
    Swap matrix S on H = H_L ⊕ H_R, dim(H_L) = dim(H_R) = dim_left.
    """
    S = np.zeros((2 * dim_left, 2 * dim_left), dtype=complex)
    S[:dim_left, dim_left:] = np.eye(dim_left)
    S[dim_left:, :dim_left] = np.eye(dim_left)
    return S


def build_gamma_F(dim_left: int) -> np.ndarray:
    """
    Grading operator gamma_F with eigenvalue -1 on H_L and +1 on H_R.
    """
    g = np.zeros((2 * dim_left, 2 * dim_left), dtype=complex)
    g[:dim_left, :dim_left] = -np.eye(dim_left)
    g[dim_left:, dim_left:] = +np.eye(dim_left)
    return g


def J_action(M: np.ndarray, S: np.ndarray, phase: complex = 1.0) -> np.ndarray:
    """
    Implement J M J^{-1} = phase * S M^* S^T, where S encodes the L/R swap
    (and potentially more structure in the future).
    """
    return phase * (S @ M.conj() @ S.T)

def evenize_geom_dirac_and_grading(D_geom: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Standard evenization of an odd geometric triple:
      H -> H ⊕ H
      D_even = [[0, D],[D, 0]]
      gamma  = [[I, 0],[0,-I]]
      lift(A)= [[A, 0],[0, A]]

    Returns: (D_even, gamma_geom, lift_matrix)
      lift_matrix is a 2x2 block operator that can lift any A via:
        A_even = lift(A) = L @ kron(A, I) ??? (we implement explicit lift below)
    """
    dim = D_geom.shape[0]
    I = np.eye(dim, dtype=complex)
    Z = np.zeros((dim, dim), dtype=complex)

    D_even = np.block([[Z, D_geom.astype(complex)],
                       [D_geom.astype(complex), Z]])
    gamma = np.block([[ I, Z],
                      [ Z,-I]])

    # We return gamma; lifting is done by lift_even(A) below.
    return D_even, gamma, I  # I returned only for convenience


def lift_even(A: np.ndarray) -> np.ndarray:
    """Lift an operator A on H to diag(A, A) on H ⊕ H."""
    dim = A.shape[0]
    Z = np.zeros((dim, dim), dtype=complex)
    return np.block([[A.astype(complex), Z],
                     [Z, A.astype(complex)]])
def build_product_dirac_v4(D_geom_even: np.ndarray, gamma_geom: np.ndarray, D_F: np.ndarray) -> np.ndarray:
    """
    v4.0 even product:
      D = D_geom_even ⊗ I_F + gamma_geom ⊗ D_F
    """
    dim_geom = D_geom_even.shape[0]
    dimF = D_F.shape[0]
    I_F = np.eye(dimF, dtype=complex)
    return kron(D_geom_even, I_F) + kron(gamma_geom.astype(complex), D_F.astype(complex))

def build_geom_algebra_generators(N: int) -> dict:
    dim = 2 * N + 1
    n_vals = np.arange(-N, N + 1, dtype=int)
    I_geom = np.eye(dim, dtype=complex)

    gens = {"I_geom": I_geom}

    def proj_div(d: int) -> np.ndarray:
        mask = (n_vals % d == 0)
        return np.diag(mask.astype(float))

    for d in [2, 3, 5]:
        gens[f"P_div_{d}"] = proj_div(d)

    return gens

def build_product_algebra(N: int, *, even_geom: bool) -> tuple[list[np.ndarray], list[str], list[str]]:
    """
    Returns:
      ops_prod, labels_prod, tags_prod
    where tags_prod ∈ {"geom", "finite"} for category diagnostics.
    """
    geom_gens = build_geom_algebra_generators(N)  # acts on H_geom (undoubled)
    I_geom = geom_gens["I_geom"]

    ops_F, labels_F = build_internal_algebra_ops()
    dimF = ops_F[0].shape[0]
    I_F = np.eye(dimF, dtype=complex)

    ops_prod: list[np.ndarray] = []
    labels_prod: list[str] = []
    tags_prod: list[str] = []

    # Lift geom generators if we are in v4 (evenized geom)
    for name, A_geom in geom_gens.items():
        A = lift_even(A_geom) if even_geom else A_geom
        ops_prod.append(kron(A, I_F))
        labels_prod.append(f"{name}⊗I_F")
        tags_prod.append("geom")

    # Finite part always acts as identity on geom factor (even or odd)
    I_geom_used = lift_even(I_geom) if even_geom else I_geom

    for A_F, lab in zip(ops_F, labels_F):
        ops_prod.append(kron(I_geom_used, A_F))
        labels_prod.append(f"I_geom⊗{lab}")
        tags_prod.append("finite")

    return ops_prod, labels_prod, tags_prod

def check_v4_grading_axioms(
    D_geom_even: np.ndarray,
    gamma_geom: np.ndarray,
    ops_prod: list[np.ndarray],
    tags_prod: list[str],
    dimF: int,
    eps: float = 1e-12,
) -> None:
    """
    v4.0 core sanity:
      {gamma_geom, D_geom_even} = 0
      [gamma_geom, pi(a_geom)] = 0   (geom algebra commutes with gamma)
      [gamma_geom, I_geom] = 0       (trivial)
    We check these after tensoring to the full Hilbert space where needed.
    """
    print("=== v4.0 grading sanity checks ===")

    # 1) geometric anti-commutation
    anti = gamma_geom @ D_geom_even + D_geom_even @ gamma_geom
    n1 = np.linalg.norm(anti, ord="fro")
    print(f"||{{γ_geom, D_geom_even}}||_F = {n1:.3e}  (target ~0)")

    # 2) gamma commutes with geom algebra (lifted), checked on full product reps that are tagged "geom"
    # Build γ_full = γ_geom ⊗ I_F
    I_F = np.eye(dimF, dtype=complex)
    gamma_full = kron(gamma_geom.astype(complex), I_F)

    max_comm = 0.0
    for A, tag in zip(ops_prod, tags_prod):
        if tag != "geom":
            continue
        comm = gamma_full @ A - A @ gamma_full
        max_comm = max(max_comm, np.linalg.norm(comm, ord="fro"))

    print(f"max ||[γ_full, π(a_geom)]||_F = {max_comm:.3e}  (target ~0)")
    if n1 > eps or max_comm > eps:
        print(f"WARNING: grading sanity exceeds eps={eps:.1e}")
    print()

# ===========================================
# Order-zero and first-order checkers (finite triple)
# ===========================================

def check_first_order_condition(
    D_F: np.ndarray,
    algebra: List[AlgebraElement],
    eps: float = 1e-12,
    J_phase: complex = 1.0,
    S: np.ndarray | None = None,
    topk: int = 10,
) -> FirstOrderResult:
    """
    First-order: [[D, a], J b J^{-1}] = 0 for all a,b in algebra.

    If S is None, default to the simple H_L ⊕ H_R swap.
    """
    assert_square(D_F, "D_F")
    n = D_F.shape[0]

    if S is None:
        assert_even_dim(D_F, "D_F")
        S = build_swap_LR(n // 2)

    max_norm = 0.0
    worst_pair: Optional[Tuple[str, str]] = None

    # We store violating pairs (label_a, label_b, norm) up to topk largest.
    viols: List[Tuple[str, str, float]] = []

    # Precompute J b J^{-1}
    Jb_list = [J_action(b.op, S, phase=J_phase) for b in algebra]

    for a in algebra:
        Da = D_F @ a.op - a.op @ D_F
        for b, Jb in zip(algebra, Jb_list):
            comm2 = Da @ Jb - Jb @ Da
            nrm = frob_norm(comm2)

            if nrm > max_norm:
                max_norm = nrm
                worst_pair = (a.label, b.label)

            if nrm > eps:
                viols.append((a.label, b.label, nrm))

    # keep only topk by norm
    viols.sort(key=lambda x: x[2], reverse=True)
    viols = viols[:topk]

    return FirstOrderResult(
        max_norm=max_norm,
        worst_pair=worst_pair,
        good_pairs=viols,   # NOTE: field name kept for backward compatibility; these are VIOLATIONS.
    )


def check_zero_order_condition(
    algebra: List[AlgebraElement],
    eps: float = 1e-12,
    J_phase: complex = 1.0,
    S: np.ndarray | None = None,
    topk: int = 10,
) -> ZeroOrderResult:
    """
    Zero-order: [a, J b J^{-1}] = 0 for all a,b in algebra.

    If S is None, default to the simple H_L ⊕ H_R swap.
    """
    # infer n from first algebra element
    if not algebra:
        raise ValueError("algebra must be non-empty")

    n = algebra[0].op.shape[0]
    for a in algebra:
        assert_square(a.op, f"algebra element {a.label}")

    if S is None:
        if n % 2 != 0:
            raise ValueError("Default LR swap requires even-dimensional representation. Provide S explicitly.")
        S = build_swap_LR(n // 2)

    max_norm = 0.0
    worst_pair: Optional[Tuple[str, str]] = None

    bad_pairs: List[Tuple[str, str, float]] = []

    # Precompute J b J^{-1}
    Jb_list = [J_action(b.op, S, phase=J_phase) for b in algebra]

    for a in algebra:
        for b, Jb in zip(algebra, Jb_list):
            comm = a.op @ Jb - Jb @ a.op
            nrm = frob_norm(comm)

            if nrm > max_norm:
                max_norm = nrm
                worst_pair = (a.label, b.label)

            if nrm > eps:
                bad_pairs.append((a.label, b.label, nrm))

    bad_pairs.sort(key=lambda x: x[2], reverse=True)
    bad_pairs = bad_pairs[:topk]

    return ZeroOrderResult(
        max_norm=max_norm,
        worst_pair=worst_pair,
        bad_pairs=bad_pairs,
    )


def print_first_order_result(res: FirstOrderResult, eps: float = 1e-12) -> None:
    print("=== First-order condition ===")
    print(f"max ||[[D,a], JbJ^-1]||_F = {res.max_norm:.3e}")
    if res.worst_pair is not None:
        print(f"worst pair: a={res.worst_pair[0]}, b={res.worst_pair[1]}")
    if res.good_pairs:  # these are actually violations (see note in checker)
        print(f"Top violations (> eps={eps:.1e}):")
        for la, lb, nrm in res.good_pairs:
            print(f"  a={la:>12s}, b={lb:>12s}  -> {nrm:.3e}")
    else:
        print(f"All pairs satisfy first-order within eps={eps:.1e}")
    print()


def print_zero_order_result(res: ZeroOrderResult, eps: float = 1e-12) -> None:
    print("=== Zero-order condition ===")
    print(f"max ||[a, JbJ^-1]||_F = {res.max_norm:.3e}")
    if res.worst_pair is not None:
        print(f"worst pair: a={res.worst_pair[0]}, b={res.worst_pair[1]}")
    if res.bad_pairs:
        print(f"Top violations (> eps={eps:.1e}):")
        for la, lb, nrm in res.bad_pairs:
            print(f"  a={la:>12s}, b={lb:>12s}  -> {nrm:.3e}")
    else:
        print(f"All pairs satisfy zero-order within eps={eps:.1e}")
    print()

# ===========================================
# First-order condition
# ===========================================

def _topk_insert(lst, item, k):
    # lst is list[(norm, info)] sorted descending by norm
    lst.append(item)
    lst.sort(key=lambda x: x[0], reverse=True)
    del lst[k:]


def test_first_order_condition_product(
    D: np.ndarray,
    ops: list[np.ndarray],
    labels: list[str],
    S_prod: np.ndarray,
    eps: float = 1e-12,
    tags: list[str] | None = None,
    topk: int = 10,
) -> None:
    print("=== First-order condition test (v4 gate) ===")
    max_norm = 0.0
    worst: list[tuple[float, str]] = []

    # Optional category maxima
    cat_max = {
        ("geom","geom"): 0.0,
        ("geom","finite"): 0.0,
        ("finite","geom"): 0.0,
        ("finite","finite"): 0.0,
    } if tags is not None else None

    for i, a in enumerate(ops):
        Da = D @ a - a @ D
        for j, b in enumerate(ops):
            b_tilde = J_action(S_prod, b)
            comm2 = Da @ b_tilde - b_tilde @ Da
            norm = np.linalg.norm(comm2, ord="fro")

            if norm > max_norm:
                max_norm = norm

            if norm > eps:
                info = f"(a={labels[i]}, b={labels[j]})  ||[[D,a], JbJ^-1]||_F={norm:.3e}"
                _topk_insert(worst, (norm, info), topk)

            if cat_max is not None:
                key = (tags[i], tags[j])
                cat_max[key] = max(cat_max[key], norm)

    print(f"Max Frobenius norm over all pairs (a,b): {max_norm:.3e}")
    if cat_max is not None:
        print("Category maxima (a-tag, b-tag):")
        for k, v in cat_max.items():
            print(f"  {k}: {v:.3e}")

    if worst:
        print(f"Top violations > eps={eps:.1e}:")
        for _, info in worst:
            print(" ", info)
    else:
        print(f"All pairs satisfy first-order within eps={eps:.1e}")
    print()


def test_zero_order_condition_product(
    ops: list[np.ndarray],
    labels: list[str],
    S_prod: np.ndarray,
    eps: float = 1e-12,
    tags: list[str] | None = None,
    topk: int = 10,
) -> None:
    print("=== Zero-order condition test (v4 gate) ===")
    max_norm = 0.0
    worst: list[tuple[float, str]] = []

    cat_max = {
        ("geom","geom"): 0.0,
        ("geom","finite"): 0.0,
        ("finite","geom"): 0.0,
        ("finite","finite"): 0.0,
    } if tags is not None else None

    for i, a in enumerate(ops):
        for j, b in enumerate(ops):
            b_tilde = J_action(S_prod, b)
            comm = a @ b_tilde - b_tilde @ a
            norm = np.linalg.norm(comm, ord="fro")

            max_norm = max(max_norm, norm)

            if norm > eps:
                info = f"(a={labels[i]}, b={labels[j]})  ||[a, JbJ^-1]||_F={norm:.3e}"
                _topk_insert(worst, (norm, info), topk)

            if cat_max is not None:
                key = (tags[i], tags[j])
                cat_max[key] = max(cat_max[key], norm)

    print(f"Max Frobenius norm over all pairs (a,b): {max_norm:.3e}")
    if cat_max is not None:
        print("Category maxima (a-tag, b-tag):")
        for k, v in cat_max.items():
            print(f"  {k}: {v:.3e}")

    if worst:
        print(f"Top violations > eps={eps:.1e}:")
        for _, info in worst:
            print(" ", info)
    else:
        print(f"All pairs satisfy zero-order within eps={eps:.1e}")
    print()


# ===========================================
# Grading & reality
# ===========================================

# in ncg_tests.py

def check_grading_and_reality(
    D_F: np.ndarray,
    algebra: List[AlgebraElement],
    J_phase: complex = 1.0,
    gamma_F: np.ndarray | None = None,
    S: np.ndarray | None = None,
) -> GradingRealityResult:
    """
    - Check gamma_F anticommutes with D_F and commutes with A_F.
    - Check J^2 = 1 (via S^2).
    - Estimate KO-dimension sign via J D_F J^{-1} = ± D_F.

    If gamma_F or S are None, fall back to the simple H_L ⊕ H_R split.
    """
    assert_even_dim(D_F, "D_F")
    n = D_F.shape[0]
    dim_left = n // 2

    if gamma_F is None:
        gamma_F = build_gamma_F(dim_left)      # old behavior
    if S is None:
        S = build_swap_LR(dim_left)            # old behavior

    # {gamma_F, D_F}
    anti = gamma_F @ D_F + D_F @ gamma_F
    anticom_norm = frob_norm(anti)

    # [gamma_F, a]
    max_comm_gamma = 0.0
    for a in algebra:
        comm = gamma_F @ a.op - a.op @ gamma_F
        max_comm_gamma = max(max_comm_gamma, frob_norm(comm))

    # J^2 - I
    S2 = S @ S
    J2_deviation = frob_norm(S2 - np.eye(n))

    # KO-sign
    JDJ = J_action(D_F, S, phase=J_phase)
    norm_plus  = frob_norm(JDJ - D_F)
    norm_minus = frob_norm(JDJ + D_F)

    df_norm = frob_norm(D_F)
    scale = df_norm + 1e-16

    rel_plus = norm_plus / scale
    rel_minus = norm_minus / scale


    if rel_plus < 1e-12 and rel_minus > rel_plus:
        ko_sign: int | None = +1
    elif rel_minus < 1e-12 and rel_plus > rel_minus:
        ko_sign = -1
    else:
        ko_sign = None

    return GradingRealityResult(
        anticom_norm=anticom_norm,
        max_comm_gamma=max_comm_gamma,
        J2_deviation=J2_deviation,
        norm_plus=norm_plus,
        norm_minus=norm_minus,
        df_norm=df_norm,
        rel_plus=rel_plus,
        rel_minus=rel_minus,
        ko_sign=ko_sign,
    )


def print_grading_and_reality_result(
    result: GradingRealityResult,
    title: str = "Grading & reality tests",
) -> None:
    print(f"=== {title} ===")
    print(f"||{{gamma_F, D_F}}||_F         = {result.anticom_norm:.3e}")
    print(f"max ||[gamma_F, a]||_F         = {result.max_comm_gamma:.3e}")
    print(f"||J^2 - I||_F                  = {result.J2_deviation:.3e}")
    print(f"||D_F||_F                      = {result.df_norm:.3e}")
    print(f"||J D_F J^-1 - D_F||_F         = {result.norm_plus:.3e}   (rel {result.rel_plus:.3e})")
    print(f"||J D_F J^-1 + D_F||_F         = {result.norm_minus:.3e}  (rel {result.rel_minus:.3e})")

    if result.ko_sign == +1:
        print("→ KO-sign: J D_F J^-1 = + D_F (J-even). Expect rel_minus ≈ 2.")
    elif result.ko_sign == -1:
        print("→ KO-sign: J D_F J^-1 = - D_F (J-odd). Expect rel_plus ≈ 2.")
    else:
        print("→ KO-sign: ambiguous at numerical precision.")
    print()


# ===========================================
# Master convenience function
# ===========================================

def run_ncg_test_suite(
    D_F: np.ndarray,
    algebra: List[AlgebraElement],
    eps_first: float = 1e-12,
    eps_zero: float = 1e-12,
    J_phase: complex = 1.0,
    gamma_F: np.ndarray | None = None,
    S: np.ndarray | None = None,
    name: str = "",
) -> Tuple[FirstOrderResult, ZeroOrderResult, GradingRealityResult]:
    """
    Run first-order, zero-order, and grading/reality tests and pretty-print a summary.

    If gamma_F and/or S are provided, they are used for grading/reality;
    and S is used for order tests (so your SM particle↔antiparticle J is honored).
    """
    if name:
        print(f"=== NCG test suite for {name} ===")

    fo_res = check_first_order_condition(D_F, algebra, eps=eps_first, J_phase=J_phase, S=S, topk=TOPK_VIOLATIONS)
    print_first_order_result(fo_res, eps=eps_first)

    zo_res = check_zero_order_condition(algebra, eps=eps_zero, J_phase=J_phase, S=S, topk=TOPK_VIOLATIONS)
    print_zero_order_result(zo_res, eps=eps_zero)

    gr_res = check_grading_and_reality(D_F, algebra, J_phase=J_phase, gamma_F=gamma_F, S=S)
    print_grading_and_reality_result(gr_res)

    return fo_res, zo_res, gr_res



def hermitian_parametrization(n: int):
    """
    Return:
      - basis: list of n x n Hermitian matrices E_k
      so that any Hermitian D can be written as D = sum_k x_k E_k
      with real coefficients x_k.
    """
    basis = []
    # diagonal basis
    for i in range(n):
        M = np.zeros((n, n), dtype=complex)
        M[i, i] = 1.0
        basis.append(M)
    # off-diagonal (i<j): real-symmetric and imaginary-antisymmetric parts
    for i in range(n):
        for j in range(i+1, n):
            M_re = np.zeros((n, n), dtype=complex)
            M_im = np.zeros((n, n), dtype=complex)
            M_re[i, j] = 1.0
            M_re[j, i] = 1.0
            M_im[i, j] = 1.0j
            M_im[j, i] = -1.0j
            basis.append(M_re)
            basis.append(M_im)
    return basis  # length = n + 2 * n*(n-1)/2 = n^2
def build_first_order_constraint_matrix(basis_D_reduced,
                                        ops: list[np.ndarray],
                                        S: np.ndarray) -> np.ndarray:
    """
    Build constraint matrix V for first-order condition only:

       [[D, a_i], J a_j J^{-1}] = 0   for all i,j,

    where D = sum_k x_k E_k, with E_k in basis_D_reduced (already gamma-odd, J-even).

    Returns:
      V: complex matrix of shape (num_constraints, num_unknowns)
         such that V @ x = 0 encodes all constraints.
    """
    n = ops[0].shape[0]
    num_unknowns = len(basis_D_reduced)
    num_ops = len(ops)
    block_size = n * n
    num_blocks = num_ops * num_ops
    vec_len = num_blocks * block_size

    V = np.zeros((vec_len, num_unknowns), dtype=complex)

    # Precompute J a_j J^{-1}
    Jops = [S @ a.conj() @ S.T for a in ops]

    for k, E in enumerate(basis_D_reduced):
        offset = 0
        for a in ops:
            Da = E @ a - a @ E
            for btilde in Jops:
                C_first = Da @ btilde - btilde @ Da
                V[offset:offset + block_size, k] = C_first.reshape(-1)
                offset += block_size

    return V
def find_DF_solution_basis(basis_D_reduced, V: np.ndarray,
                           tol: float = 1e-12) -> list[np.ndarray]:
    """
    Given reduced basis {E_k} and constraint matrix V (M x K),
    find a basis of Hermitian matrices D_alpha = sum_k x_k^{(alpha)} E_k
    spanning the nullspace of V @ x = 0.

    Returns:
      DF_basis: list of n x n Hermitian matrices D_alpha.
    """
    # SVD: V = U diag(s) Vh, right-singular vectors in rows of Vh
    U, s, Vh = np.linalg.svd(V, full_matrices=False)

    null_mask = (s < tol)
    null_vectors = Vh[null_mask, :]   # shape: (num_null, K)

    DF_basis = []
    for vec in null_vectors:
        D = np.zeros_like(basis_D_reduced[0])
        for coeff, E in zip(vec, basis_D_reduced):
            D += coeff * E
        # ensure Hermitian numerically
        D = 0.5 * (D + D.conj().T)
        DF_basis.append(D)

    return DF_basis
def print_coupling_pattern(D: np.ndarray,
                           basis: List[SMState],
                           name_to_index: Dict[str, int],
                           thresh: float = 1e-6):
    """
    Print which L↔R pairs (for each species) are significantly coupled by D.
    """
    def show_pair(a, b):
        i = name_to_index[a]
        j = name_to_index[b]
        z = D[i, j]
        if abs(z) > thresh:
            print(f"  {a:>10s} <-> {b:<10s}  |D_ij| = {abs(z):.3e}")

    print("Lepton couplings:")
    if "nu_R" in name_to_index:
        show_pair("nu_L", "nu_R")
    show_pair("e_L", "e_R")

    print("Quark couplings (per color):")
    colors = ["r", "g", "b"]
    for c in colors:
        show_pair(f"u_L_{c}", f"u_R_{c}")
    for c in colors:
        show_pair(f"d_L_{c}", f"d_R_{c}")

    print("Antiparticle couplings (sanity check):")
    if "bar_nu_R" in name_to_index:
        show_pair("bar_nu_L", "bar_nu_R")
    show_pair("bar_e_L", "bar_e_R")
    for c in colors:
        show_pair(f"bar_u_L_{c}", f"bar_u_R_{c}")
    for c in colors:
        show_pair(f"bar_d_L_{c}", f"bar_d_R_{c}")

def main():
    basis, idx = build_sm_basis_1gen(include_nu_R=True)

    # SM-like Yukawas (1 generation, rough magnitudes)
    Y_e  = 2.94e-6     # me / v
    Y_nu = 1.0e-12     # tiny
    Y_u  = 1.3e-5      # up-type
    Y_d  = 2.8e-5      # down-type

    D_F_SM = build_DF_SM_1gen(Y_e, Y_nu, Y_u, Y_d, basis, idx)
    sm_generators = build_SM_algebra_generators(basis, idx)
    algebra_SM = [AlgebraElement(op.label, op.op) for op in sm_generators]

    # SM-specific gamma and J swap
    gamma_SM = build_gamma_F_SM(basis)
    S_SM = build_swap_particle_antiparticle(basis, idx)

    # Full suite using SM gamma + SM J
    run_ncg_test_suite(
        D_F_SM,
        algebra_SM,
        eps_first=1e-12,
        eps_zero=1e-12,
        gamma_F=gamma_SM,
        S=S_SM,
        name="SM-like 1gen finite triple (SM gamma,J)",
    )



if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Alignment Spectral Triple v4.0 — Engine Outline (Python)

Goal:
  Build a production-ready v4.0 engine that:
    (1) defines the internal SM basis (32 states per generation),
    (2) defines lifted projectors, gamma_SM, and real structure J,
    (3) builds channel maps V_u,V_d,V_e,V_nu,W_R,
    (4) assembles D_int(Y[K]) as a 32N×32N sparse block matrix,
    (5) assembles the even product Dirac D = D_geom⊗1 + γ_geom⊗D_int,
    (6) builds one-forms and inner fluctuations D_A = D + A + JAJ^{-1},
    (7) runs commutator gates (order-zero, first-order, fluctuation stability),
    (8) provides sector projectors + extraction utilities for Yukawas/PMNS,
    (9) provides serialization + reproducibility and scan harness hooks.

This is an OUTLINE: signatures + responsibilities + minimal stubs.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple, Callable

import numpy as np
from numpy.linalg import eigh, norm
from scipy.linalg import expm


# ============================================================
#  Core numerics helpers
# ============================================================

def op_norm(A: np.ndarray) -> float:
    """Operator norm (spectral norm): ||A||_2."""
    # robust: largest singular value
    s = np.linalg.svd(A, compute_uv=False)
    return float(s[0]) if s.size else 0.0

def hermitian_part(A: np.ndarray) -> np.ndarray:
    return 0.5 * (A + A.conj().T)

def comm(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    return A @ B - B @ A

def anti_comm(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    return A @ B + B @ A

def is_hermitian(A: np.ndarray, atol: float = 1e-12) -> bool:
    return bool(np.allclose(A, A.conj().T, atol=atol))

def is_projector(P: np.ndarray, atol: float = 1e-12) -> bool:
    return bool(np.allclose(P @ P, P, atol=atol) and is_hermitian(P, atol=atol))

def put_block(D: np.ndarray, i: int, j: int, B: np.ndarray) -> None:
    """Place an NxN block into SM slot (i,j) of a (32N)x(32N) matrix."""
    N = B.shape[0]
    r = slice(i * N, (i + 1) * N)
    c = slice(j * N, (j + 1) * N)
    D[r, c] += B


# ============================================================
#  v4.0 Configuration and parameter surfaces
# ============================================================

@dataclass(frozen=True)
class V40Params:
    """
    Parameters that vary scan-to-scan (textures, scales, cutoffs).
    Treat these as immutable per run for reproducibility.
    """
    N_flav: int

    # Higgs vev (for physical scaling, optional)
    higgs_vev_GeV: float = 174.0

    # Majorana scale (can be a matrix via builder hook)
    seesaw_M_GeV: float = 2.0e14

    # Gate tolerances
    eps_order0: float = 1e-10
    eps_first: float = 1e-10
    eps_sa: float = 1e-10
    eps_fluct: float = 1e-10

    # Optional scan knobs
    beta: float = 1.5
    rel_cut: float = 0.15
    tol_rel_blocks: float = 0.03


@dataclass(frozen=True)
class TexturePack:
    """
    Flavor textures (operators on H_flav), one per fermion sector.
    Each Y_* is NxN complex.
    """
    Yu: np.ndarray
    Yd: np.ndarray
    Ye: np.ndarray
    Ynu: np.ndarray
    MR: Optional[np.ndarray] = None  # NxN complex symmetric preferred


# ============================================================
#  SM basis (32 states) + projectors + grading
# ============================================================

class SMBasis32:
    """
    Owns the canonical 32-state ordering (16 particle + 16 conjugate),
    and provides projector matrices P_X in M_32(C) and grading gamma_SM.
    """

    def __init__(self) -> None:
        self.dim = 32
        self._index_map = self._build_index_map()  # names -> indices
        self.state_labels = self._build_state_labels()  # 32 labels in canonical order
        self._label_to_index = {lbl: i for i, lbl in enumerate(self.state_labels)}

        # cached projectors (32x32)
        self.P = self._build_projectors()

        # grading gamma_SM (32x32)
        self.gamma_SM = self._build_gamma_SM()

    # -------- basis + indices --------

    def _build_index_map(self) -> Dict[str, List[int]]:
        """
        Return dict mapping semantic labels -> index lists in the 32 ordering.
        Must match the fixed ordering used in the v4.0 paper.
        """
        # Particle: 0..15
        # (uL r,g,b, dL r,g,b, nuL, eL, uR r,g,b, dR r,g,b, eR, nuR)
        # Conjugates: +16
        uL = [0, 1, 2]
        dL = [3, 4, 5]
        nuL = [6]
        eL = [7]
        uR = [8, 9, 10]
        dR = [11, 12, 13]
        eR = [14]
        nuR = [15]

        # conjugates
        uL_c = [i + 16 for i in uL]
        dL_c = [i + 16 for i in dL]
        nuL_c = [i + 16 for i in nuL]
        eL_c = [i + 16 for i in eL]
        uR_c = [i + 16 for i in uR]
        dR_c = [i + 16 for i in dR]
        eR_c = [i + 16 for i in eR]
        nuR_c = [i + 16 for i in nuR]

        return {
            "uL": uL, "dL": dL, "nuL": nuL, "eL": eL,
            "uR": uR, "dR": dR, "eR": eR, "nuR": nuR,
            "uL_c": uL_c, "dL_c": dL_c, "nuL_c": nuL_c, "eL_c": eL_c,
            "uR_c": uR_c, "dR_c": dR_c, "eR_c": eR_c, "nuR_c": nuR_c,
            "part": list(range(0, 16)),
            "conj": list(range(16, 32)),
            "L_part": list(range(0, 8)),
            "R_part": list(range(8, 16)),
            "L_conj": list(range(24, 32)),
            "R_conj": list(range(16, 24)),
        }

    def indices(self, key: str) -> List[int]:
        return list(self._index_map[key])


    def label(self, i: int) -> str:
        """Human-readable label for basis index i (0..31)."""
        return self.state_labels[i]

    def index_of(self, label: str) -> int:
        """Inverse lookup: label -> index."""
        return int(self._label_to_index[label])

    def _build_state_labels(self) -> List[str]:
        """
        Canonical 32-state labels matching the fixed ordering:
          0..15  : particles
          16..31 : charge-conjugates (suffix _c)
        """
        # particles
        labels: List[str] = [
            "uL_r","uL_g","uL_b",
            "dL_r","dL_g","dL_b",
            "nuL","eL",
            "uR_r","uR_g","uR_b",
            "dR_r","dR_g","dR_b",
            "eR","nuR",
        ]
        # conjugates
        labels += [f"{x}_c" for x in labels]
        return labels

    # -------- projectors --------

    def _proj_from_indices(self, idx: Sequence[int]) -> np.ndarray:
        P = np.zeros((self.dim, self.dim), dtype=complex)
        for i in idx:
            P[i, i] = 1.0
        return P

    def _build_projectors(self) -> Dict[str, np.ndarray]:
        """
        Build the basic sector projectors needed for v4.0 extraction and assembly.
        """
        P: Dict[str, np.ndarray] = {}
        for k, idx in self._index_map.items():
            P[k] = self._proj_from_indices(idx)
        return P

    # -------- grading gamma_SM --------

    def _build_gamma_SM(self) -> np.ndarray:
        """
        v4.0 grading on the SM finite space (32x32).

        Convention used earlier:
          +1 on: L particles (0..7) and L conjugates (24..31)
          -1 on: R particles (8..15) and R conjugates (16..23)
        """
        g = np.zeros((self.dim,), dtype=float)
        for i in self._index_map["L_part"] + self._index_map["L_conj"]:
            g[i] = +1.0
        for i in self._index_map["R_part"] + self._index_map["R_conj"]:
            g[i] = -1.0
        return np.diag(g).astype(complex)

    # -------- lifted operators --------

    def lift_projector(self, P32: np.ndarray, N: int) -> np.ndarray:
        """Lift P (32x32) to (32N)x(32N): P ⊗ I_N."""
        return np.kron(P32, np.eye(N, dtype=complex))

    def lift_gamma(self, N: int) -> np.ndarray:
        return self.lift_projector(self.gamma_SM, N)


# ============================================================
#  Real structures (antiunitary) on finite spaces
# ============================================================

@dataclass(frozen=True)
class AntiUnitary:
    """
    Antiunitary J implemented as J = U K where K is entrywise conjugation.
    We implement conjugation action on operators: J X J^{-1} = U \bar{X} U^\dagger.
    """
    U: np.ndarray  # unitary matrix

    def conj_op(self, X: np.ndarray) -> np.ndarray:
        """Return J X J^{-1}."""
        return self.U @ np.conj(X) @ self.U.conj().T

    def square_sign(self, atol: float = 1e-12) -> Optional[int]:
        """
        Return sign of J^2 if inferable in this basis (optional).
        For J=U K: J^2 = U \bar{U}.
        """
        J2 = self.U @ np.conj(self.U)
        # if J2 ≈ ±I, return that sign
        I = np.eye(J2.shape[0], dtype=complex)
        if np.allclose(J2, I, atol=atol):
            return +1
        if np.allclose(J2, -I, atol=atol):
            return -1
        return None


class RealStructureFactory:
    """
    Factory for building J_SM, J_flav, and lifted product J_int = J_SM ⊗ J_flav.
    """

    @staticmethod
    def build_J_flav(N: int, U_flav: Optional[np.ndarray] = None) -> AntiUnitary:
        """
        Default: J_flav = K (so U_flav = I).
        """
        if U_flav is None:
            U_flav = np.eye(N, dtype=complex)
        return AntiUnitary(U=U_flav)

    @staticmethod
    def build_J_SM(basis: SMBasis32) -> AntiUnitary:
        """
        Concrete finite real structure for the 32-state basis.

        We choose the standard "swap" charge-conjugation that maps each particle
        basis vector to its conjugate partner and vice-versa:
            J |i> = |i+16>   for i=0..15
            J |i> = |i-16>   for i=16..31

        Implemented as J = U K with:
            U = [[0, I_16],
                 [I_16, 0]]

        This choice is unitary, involutive (J^2 = +1 in this basis), and matches
        the fixed ordering used by SMBasis32.
        """
        if basis.dim != 32:
            raise ValueError(f"Expected SMBasis32 dim=32, got {basis.dim}")

        I16 = np.eye(16, dtype=complex)
        Z16 = np.zeros((16, 16), dtype=complex)
        U_SM = np.block([[Z16, I16],
                         [I16, Z16]])
        return AntiUnitary(U=U_SM)
    @staticmethod
    def kron(J1: AntiUnitary, J2: AntiUnitary) -> AntiUnitary:
        """Product antiunitary: (U1⊗U2) K."""
        return AntiUnitary(U=np.kron(J1.U, J2.U))


# ============================================================
#  Channel maps V_u,V_d,V_e,V_nu and W_R
# ============================================================

class ChannelMaps:
    """
    Builds partial isometries on H_SM (32x32):
      V_u,V_d,V_e,V_nu  and W_R (nu_R -> nu_R^c).

    The minimal concrete implementation uses the fixed index ordering
    and builds each map as a sum of matrix units E_ij (color-preserving).
    """

    def __init__(self, basis: SMBasis32) -> None:
        self.basis = basis
        self.Vu = self._build_Vu()
        self.Vd = self._build_Vd()
        self.Ve = self._build_Ve()
        self.Vnu = self._build_Vnu()
        self.WR = self._build_WR()

    def _E(self, i: int, j: int) -> np.ndarray:
        E = np.zeros((self.basis.dim, self.basis.dim), dtype=complex)
        E[i, j] = 1.0
        return E

    def _sum_E(self, pairs: Sequence[Tuple[int, int]]) -> np.ndarray:
        V = np.zeros((self.basis.dim, self.basis.dim), dtype=complex)
        for (i, j) in pairs:
            V += self._E(i, j)
        return V

    def _build_Vu(self) -> np.ndarray:
        uL = self.basis.indices("uL")
        uR = self.basis.indices("uR")
        pairs = list(zip(uR, uL))  # uR <- uL (color-preserving)
        return self._sum_E(pairs)

    def _build_Vd(self) -> np.ndarray:
        dL = self.basis.indices("dL")
        dR = self.basis.indices("dR")
        pairs = list(zip(dR, dL))
        return self._sum_E(pairs)

    def _build_Ve(self) -> np.ndarray:
        eL = self.basis.indices("eL")[0]
        eR = self.basis.indices("eR")[0]
        return self._sum_E([(eR, eL)])

    def _build_Vnu(self) -> np.ndarray:
        nuL = self.basis.indices("nuL")[0]
        nuR = self.basis.indices("nuR")[0]
        return self._sum_E([(nuR, nuL)])

    def _build_WR(self) -> np.ndarray:
        nuR = self.basis.indices("nuR")[0]
        nuR_c = self.basis.indices("nuR_c")[0]
        return self._sum_E([(nuR_c, nuR)])

    # ---- conjugate maps via J_SM ----

    def conjugate_map(self, V: np.ndarray, J_SM: AntiUnitary) -> np.ndarray:
        """V^c = J_SM V J_SM^{-1}."""
        return J_SM.conj_op(V)



    def validate(self, atol: float = 1e-12) -> Dict[str, float]:
        """
        Sanity checks for partial isometries:
          Vu†Vu = P_uL, VuVu† = P_uR, etc.
        Returns operator-norm deviations for each identity.
        """
        b = self.basis
        out: Dict[str, float] = {}

        def dev(A: np.ndarray, B: np.ndarray) -> float:
            return float(op_norm(A - B))

        Vu, Vd, Ve, Vnu, WR = self.Vu, self.Vd, self.Ve, self.Vnu, self.WR

        out["Vu_dagVu"] = dev(Vu.conj().T @ Vu, b.P["uL"])
        out["VuVu_dag"] = dev(Vu @ Vu.conj().T, b.P["uR"])

        out["Vd_dagVd"] = dev(Vd.conj().T @ Vd, b.P["dL"])
        out["VdVd_dag"] = dev(Vd @ Vd.conj().T, b.P["dR"])

        out["Ve_dagVe"] = dev(Ve.conj().T @ Ve, b.P["eL"])
        out["VeVe_dag"] = dev(Ve @ Ve.conj().T, b.P["eR"])

        out["Vnu_dagVnu"] = dev(Vnu.conj().T @ Vnu, b.P["nuL"])
        out["VnuVnu_dag"] = dev(Vnu @ Vnu.conj().T, b.P["nuR"])

        # WR is a partial isometry nuR -> nuR_c
        out["WR_dagWR"] = dev(WR.conj().T @ WR, b.P["nuR"])
        out["WRWR_dag"] = dev(WR @ WR.conj().T, b.P["nuR_c"])

        # also assert they're close to partial isometries (V†V and VV† are projectors)
        # (return numeric diagnostics; raising is left to the caller)
        return out

# ============================================================
#  D_int builder (32N x 32N sparse block matrix)
# ============================================================

class DIntBuilder:
    """
    Assemble D_int(Y[K]) on H_SM ⊗ H_flav:

        D_int = Dirac(particles) + Dirac(conjugates) + Majorana

    Conventions are tied to SMBasis32.gamma_SM:
      - "left"  := gamma_SM = +1
      - "right" := gamma_SM = -1

    In the conjugate sector chirality is flipped (as encoded in SMBasis32),
    so the conjugate Dirac blocks use (V^c)† ⊗ \bar{Y} and V^c ⊗ Y^T.
    """

    def __init__(self, basis: SMBasis32, maps: ChannelMaps) -> None:
        self.basis = basis
        self.maps = maps

    # ---------------- internal validators ----------------

    @staticmethod
    def _require_square_NxN(name: str, A: np.ndarray, N: int) -> None:
        if A is None:
            raise ValueError(f"{name} is None (expected {N}x{N} complex array).")
        A = np.asarray(A)
        if A.shape != (N, N):
            raise ValueError(f"{name} must be shape ({N},{N}), got {A.shape}.")

    def _validate_textures(self, textures: TexturePack, N: int) -> None:
        self._require_square_NxN("Yu", textures.Yu, N)
        self._require_square_NxN("Yd", textures.Yd, N)
        self._require_square_NxN("Ye", textures.Ye, N)
        self._require_square_NxN("Ynu", textures.Ynu, N)
        if textures.MR is not None:
            self._require_square_NxN("MR", textures.MR, N)

    # ---------------- build ----------------

    def build_D_int(self, textures: TexturePack, J_SM: AntiUnitary, N: int) -> np.ndarray:
        """
        Return the internal Dirac operator as a dense (32N)x(32N) matrix.
        """
        self._validate_textures(textures, N)

        Yu, Yd, Ye, Ynu = (np.asarray(textures.Yu, dtype=complex),
                          np.asarray(textures.Yd, dtype=complex),
                          np.asarray(textures.Ye, dtype=complex),
                          np.asarray(textures.Ynu, dtype=complex))
        MR = None if textures.MR is None else np.asarray(textures.MR, dtype=complex)

        Vu, Vd, Ve, Vnu, WR = self.maps.Vu, self.maps.Vd, self.maps.Ve, self.maps.Vnu, self.maps.WR

        # Conjugate-channel maps on H_SM: V^c = J_SM V J_SM^{-1}.
        Vu_c = self.maps.conjugate_map(Vu, J_SM)
        Vd_c = self.maps.conjugate_map(Vd, J_SM)
        Ve_c = self.maps.conjugate_map(Ve, J_SM)
        Vnu_c = self.maps.conjugate_map(Vnu, J_SM)

        # --- particle Dirac: (R <- L)⊗Y  +  (L <- R)⊗Y†
        D_part = (
            np.kron(Vu, Yu) + np.kron(Vu.conj().T, Yu.conj().T) +
            np.kron(Vd, Yd) + np.kron(Vd.conj().T, Yd.conj().T) +
            np.kron(Ve, Ye) + np.kron(Ve.conj().T, Ye.conj().T) +
            np.kron(Vnu, Ynu) + np.kron(Vnu.conj().T, Ynu.conj().T)
        )

        # --- conjugate Dirac (lemma-consistent with your explicit +16 version)
        D_conj = (
            np.kron(Vu_c.conj().T, np.conj(Yu)) + np.kron(Vu_c, Yu.T) +
            np.kron(Vd_c.conj().T, np.conj(Yd)) + np.kron(Vd_c, Yd.T) +
            np.kron(Ve_c.conj().T, np.conj(Ye)) + np.kron(Ve_c, Ye.T) +
            np.kron(Vnu_c.conj().T, np.conj(Ynu)) + np.kron(Vnu_c, Ynu.T)
        )

        # --- Majorana: nu_R <-> nu_R^c (WR: nuR_c <- nuR)
        D_M = np.zeros_like(D_part)
        if MR is not None:
            D_M = np.kron(WR, MR) + np.kron(WR.conj().T, MR.conj().T)

        D = np.asarray(D_part + D_conj + D_M, dtype=complex)
        if D.shape != (32 * N, 32 * N):
            raise RuntimeError(f"D_int has wrong shape {D.shape}, expected {(32*N, 32*N)}.")
        return D

    # ---------------- axioms ----------------

    def check_evenness(self, D_int: np.ndarray, gamma_SM_lifted: np.ndarray, atol: float = 1e-10) -> Dict[str, float]:
        """
        v4.0 hard checks:
          (i) self-adjointness
          (ii) oddness: {gamma_SM⊗I, D_int} = 0
        """
        D_int = np.asarray(D_int, dtype=complex)
        gamma_SM_lifted = np.asarray(gamma_SM_lifted, dtype=complex)

        sa = op_norm(D_int - D_int.conj().T)
        odd = op_norm(anti_comm(gamma_SM_lifted, D_int))

        return {"sa": float(sa), "odd": float(odd), "sa_ok": sa <= atol, "odd_ok": odd <= atol}


# ============================================================
#  Geometry factor hooks (D_geom, gamma_geom, J_geom)
# ============================================================

class GeometryFactor:
    """
    Minimal hooks for the geometric triple (H_geom, D_geom, gamma_geom, J_geom).

    Defaults (backward-compatible with the template):
      - D_geom: diagonal truncation diag(modes)
      - gamma_geom: diagonal sign(modes) (zeros treated as +1)
      - J_geom: complex conjugation (AntiUnitary with U = I)

    If you want an even truncation satisfying {D_geom, gamma_geom} = 0, use
    `GeometryFactor.even_truncation(lambdas)`, which builds a chiral doubling.
    """

    def __init__(
        self,
        modes: np.ndarray,
        *,
        gamma: np.ndarray | None = None,
        J: AntiUnitary | None = None,
        zero_sign: float = +1.0,
        atol: float = 1e-12,
    ) -> None:
        self.modes = np.asarray(modes, dtype=float).reshape(-1)

        # Finite truncation diagonal Φ|n> = n|n>
        self.D_geom = np.diag(self.modes).astype(complex)

        # Chirality: default sign(modes) with a stable convention at 0.
        if gamma is None:
            s = np.sign(self.modes)
            if zero_sign not in (+1.0, -1.0):
                raise ValueError("zero_sign must be +1.0 or -1.0")
            s = np.where(s == 0.0, float(zero_sign), s)
            self.gamma_geom = np.diag(s).astype(complex)
        else:
            g = np.asarray(gamma, dtype=complex)
            if g.shape != self.D_geom.shape:
                raise ValueError(f"gamma must have shape {self.D_geom.shape}, got {g.shape}")
            self.gamma_geom = g

        # Real structure (antiunitary): default is complex conjugation in this basis.
        self.J_geom = J if J is not None else AntiUnitary(U=np.eye(self.D_geom.shape[0], dtype=complex))

        self._validate_basics(atol=atol)

    @classmethod
    def even_truncation(
        cls,
        lambdas: np.ndarray,
        *,
        J: AntiUnitary | None = None,
        atol: float = 1e-12,
    ) -> "GeometryFactor":
        """
        Build a chiral-doubled even finite triple:

            H_geom = C^n ⊕ C^n
            D_geom = [[0, Λ],
                      [Λ, 0]]   where Λ = diag(lambdas)
            gamma  = [[ I, 0],
                      [ 0,-I]]

        This satisfies {D_geom, gamma} = 0 (up to atol).
        """
        lam = np.asarray(lambdas, dtype=float).reshape(-1)
        n = lam.shape[0]
        Z = np.zeros((n, n), dtype=complex)
        Lam = np.diag(lam).astype(complex)

        obj = cls.__new__(cls)
        obj.modes = lam.copy()
        obj.D_geom = np.block([[Z, Lam], [Lam, Z]]).astype(complex)
        obj.gamma_geom = np.block([[np.eye(n, dtype=complex), Z], [Z, -np.eye(n, dtype=complex)]]).astype(complex)
        obj.J_geom = J if J is not None else AntiUnitary(U=np.eye(2 * n, dtype=complex))

        obj._validate_basics(atol=atol)

        anti = anti_comm(obj.gamma_geom, obj.D_geom)
        if op_norm(anti) > atol:
            raise ValueError(f"even_truncation: {{gamma_geom, D_geom}} not ~0 (norm={op_norm(anti)})")
        return obj

    def dim(self) -> int:
        return int(self.D_geom.shape[0])

    def lift_D(self, fin_dim: int) -> np.ndarray:
        """Return D_geom ⊗ I_fin."""
        return np.kron(self.D_geom, np.eye(fin_dim, dtype=complex))

    def lift_gamma(self, fin_dim: int) -> np.ndarray:
        """Return gamma_geom ⊗ I_fin."""
        return np.kron(self.gamma_geom, np.eye(fin_dim, dtype=complex))

    def lift_J(self, J_fin: AntiUnitary) -> AntiUnitary:
        """
        Return J_geom ⊗ J_fin as AntiUnitary with U_total = U_geom ⊗ U_fin.
        """
        return AntiUnitary(U=np.kron(self.J_geom.U, J_fin.U))

    def _validate_basics(self, *, atol: float = 1e-12) -> None:
        # D_geom should be self-adjoint.
        if op_norm(self.D_geom - self.D_geom.conj().T) > atol:
            raise ValueError("D_geom must be self-adjoint (within atol).")

        # gamma_geom should be self-adjoint and an involution.
        if op_norm(self.gamma_geom - self.gamma_geom.conj().T) > atol:
            raise ValueError("gamma_geom must be self-adjoint (within atol).")
        I = np.eye(self.gamma_geom.shape[0], dtype=complex)
        if op_norm(self.gamma_geom @ self.gamma_geom - I) > 1e-9:
            raise ValueError("gamma_geom must satisfy gamma_geom^2 = I (within tolerance).")

        # J_geom's unitary part should be unitary.
        U = self.J_geom.U
        if U.shape != self.D_geom.shape:
            raise ValueError(f"J_geom.U must have shape {self.D_geom.shape}, got {U.shape}")
        if op_norm(U.conj().T @ U - I) > 1e-9:
            raise ValueError("J_geom.U must be unitary (within tolerance).")

# ============================================================
#  Even product Dirac operator D = D_geom⊗1 + γ_geom⊗D_int
# ============================================================

class EvenProductDirac:
    """
    Build the v4.0 even product Dirac on:
        H = H_geom ⊗ (H_SM ⊗ H_flav)

    Core constructions:
        D_total  = (D_geom ⊗ I_int) + (gamma_geom ⊗ D_int)
        Gamma    = gamma_geom ⊗ (gamma_SM ⊗ I_N)

    Notes:
      - Oddness of D_total w.r.t Gamma holds if both
            {D_geom, gamma_geom} = 0   and   {D_int, gamma_SM⊗I_N} = 0.
      - This class provides hard shape checks and basic axioms checks.
    """

    def __init__(self, geom: GeometryFactor, basis: SMBasis32, params: V40Params) -> None:
        self.geom = geom
        self.basis = basis
        self.params = params

    # ---------------- builders ----------------

    def build_D(
        self,
        D_int: np.ndarray,
        *,
        scale_geom: float = 1.0,
        scale_int: float = 1.0,
    ) -> np.ndarray:
        """
        Return D_total = scale_geom*(D_geom⊗I) + scale_int*(gamma_geom⊗D_int).
        """
        D_int = np.asarray(D_int, dtype=complex)
        if D_int.ndim != 2 or D_int.shape[0] != D_int.shape[1]:
            raise ValueError(f"D_int must be square, got {D_int.shape}.")

        dG = np.asarray(self.geom.D_geom, dtype=complex)
        gG = np.asarray(self.geom.gamma_geom, dtype=complex)

        I_int = np.eye(D_int.shape[0], dtype=complex)
        D_total = (float(scale_geom) * np.kron(dG, I_int)) + (float(scale_int) * np.kron(gG, D_int))

        expected = (self.geom.dim() * D_int.shape[0], self.geom.dim() * D_int.shape[0])
        if D_total.shape != expected:
            raise RuntimeError(f"D_total has shape {D_total.shape}, expected {expected}.")
        return D_total

    def build_Gamma_total(self, N: int) -> np.ndarray:
        """
        Total grading Γ = γ_geom ⊗ (γ_SM ⊗ 1_N).
        """
        gamma_SM_lift = self.basis.lift_gamma(N)  # (32N)x(32N)
        Gamma = np.kron(np.asarray(self.geom.gamma_geom, dtype=complex), gamma_SM_lift)
        return Gamma

    def build_J_internal(self, J_SM: AntiUnitary, N: int, J_flav: AntiUnitary | None = None) -> AntiUnitary:
        """
        Internal real structure:
            J_int = J_SM ⊗ J_flav
        Default J_flav is complex conjugation on C^N (U = I_N).
        """
        if J_flav is None:
            J_flav = AntiUnitary(U=np.eye(N, dtype=complex))
        return AntiUnitary(U=np.kron(J_SM.U, J_flav.U))

    def build_J_total(self, J_int: AntiUnitary) -> AntiUnitary:
        """
        Total real structure:
            J_total = J_geom ⊗ J_int.
        """
        return self.geom.lift_J(J_int)

    # ---------------- checks ----------------

    def check_axioms(
        self,
        D_total: np.ndarray,
        D_int: np.ndarray,
        N: int,
        *,
        atol: float | None = None,
    ) -> Dict[str, float]:
        """
        Basic even-product hard checks:
          (i)  self-adjointness of D_total
          (ii) oddness: {Gamma_total, D_total} = 0
          (iii) diagnostics:
               {gamma_geom, D_geom} and {gamma_int, D_int}

        Returns numeric norms and boolean flags.
        """
        eps = float(self.params.eps_sa if atol is None else atol)

        D_total = np.asarray(D_total, dtype=complex)
        D_int = np.asarray(D_int, dtype=complex)

        Gamma = self.build_Gamma_total(N)

        sa = op_norm(D_total - D_total.conj().T)
        odd_total = op_norm(anti_comm(Gamma, D_total))

        # diagnostics: the two sufficient conditions for odd_total≈0
        dG = np.asarray(self.geom.D_geom, dtype=complex)
        gG = np.asarray(self.geom.gamma_geom, dtype=complex)
        odd_geom = op_norm(anti_comm(gG, dG))

        gamma_int = self.basis.lift_gamma(N)
        odd_int = op_norm(anti_comm(gamma_int, D_int))

        return {
            "sa_total": float(sa),
            "odd_total": float(odd_total),
            "odd_geom": float(odd_geom),
            "odd_int": float(odd_int),
            "sa_ok": sa <= eps,
            "odd_ok": odd_total <= eps,
        }


# ============================================================
#  Representations π(a) and generators for A_SM (minimal scaffold)
# ============================================================

class SMRepresentation:
    """
       Minimal representation layer needed to run commutator gates.
       In production, π is the actual NCG SM representation on H_SM.
       Here we provide an interface and placeholders.
       """

    def __init__(self, basis: SMBasis32) -> None:
        self.basis = basis

    def _build_index_map(self) -> Dict[str, List[int]]:
        uL = [0, 1, 2]
        dL = [3, 4, 5]
        nuL = [6]
        eL = [7]
        uR = [8, 9, 10]
        dR = [11, 12, 13]
        eR = [14]
        nuR = [15]

        uL_c = [i + 16 for i in uL]
        dL_c = [i + 16 for i in dL]
        nuL_c = [i + 16 for i in nuL]
        eL_c = [i + 16 for i in eL]
        uR_c = [i + 16 for i in uR]
        dR_c = [i + 16 for i in dR]
        eR_c = [i + 16 for i in eR]
        nuR_c = [i + 16 for i in nuR]

        # multiplet aliases (particle)
        Q_L = uL + dL
        L_L = nuL + eL
        Q_R = uR + dR
        L_R = eR + nuR  # (singlets)

        # multiplet aliases (conjugate)
        Q_L_c = uL_c + dL_c
        L_L_c = nuL_c + eL_c
        Q_R_c = uR_c + dR_c
        L_R_c = eR_c + nuR_c

        return {
            # fine-grained
            "uL": uL, "dL": dL, "nuL": nuL, "eL": eL,
            "uR": uR, "dR": dR, "eR": eR, "nuR": nuR,
            "uL_c": uL_c, "dL_c": dL_c, "nuL_c": nuL_c, "eL_c": eL_c,
            "uR_c": uR_c, "dR_c": dR_c, "eR_c": eR_c, "nuR_c": nuR_c,

            # aliases used by SMRepresentation / paper language
            "Q_L": Q_L,
            "L_L": L_L,
            "Q_R": Q_R,
            "L_R": L_R,
            "Q_L_c": Q_L_c,
            "L_L_c": L_L_c,
            "Q_R_c": Q_R_c,
            "L_R_c": L_R_c,

            # bookkeeping
            "part": list(range(0, 16)),
            "conj": list(range(16, 32)),
            "L_part": list(range(0, 8)),
            "R_part": list(range(8, 16)),
            "L_conj": list(range(24, 32)),
            "R_conj": list(range(16, 24)),
        }

    def generators(self) -> List[np.ndarray]:
        """
        Return a finite generating/spanning set G_A ⊂ π(A_SM).
        MUST be replaced by the true SM algebra reps for real gates.
        """
        # placeholder: projectors onto major summands
        return [
            self.basis.P["uL"] + self.basis.P["dL"],  # Q_L
            self.basis.P["nuL"] + self.basis.P["eL"],  # L_L
        ]

    # Example placeholders; in v4.0 you implement actual SM algebra action
    @property
    def P_Q_L(self) -> np.ndarray:
        # You would add these keys if you expose those subspaces explicitly.
        raise NotImplementedError

    # ----------------- full algebra representation π(λ,q,m) -----------------

    @staticmethod
    def _as2(name: str, A: np.ndarray) -> np.ndarray:
        A = np.asarray(A, dtype=complex)
        if A.shape != (2, 2):
            raise ValueError(f"{name} must be 2x2, got {A.shape}")
        return A

    @staticmethod
    def _as3(name: str, A: np.ndarray) -> np.ndarray:
        A = np.asarray(A, dtype=complex)
        if A.shape != (3, 3):
            raise ValueError(f"{name} must be 3x3, got {A.shape}")
        return A

    def pi(self, lam: complex, q: np.ndarray, m: np.ndarray, *, enforce_quaternion: bool = False) -> np.ndarray:
        """
        Algebra representation for a=(lam,q,m) ∈ C ⊕ H ⊕ M3(C).

        Minimal NCG-consistent block action:
          - On Q_L (uL,dL):         q ⊗ m
          - On L_L (nuL,eL):        q
          - On uR,dR:              lam * m
          - On nuR,eR:             lam
          - On conjugates:         complex-conjugate parameters (lam̄, q̄, m̄)

        Returns 32x32 complex matrix in the SMBasis32 ordering.
        """
        lam = complex(lam)
        q = self._as2("q", q)
        m = self._as3("m", m)

        if enforce_quaternion:
            # Optional light check: q should be in the standard complex 2x2 image of H:
            # q = [[a, b], [-b̄, ā]]
            a = q[0, 0]
            b = q[0, 1]
            target = np.array([[a, b], [-np.conj(b), np.conj(a)]], dtype=complex)
            if op_norm(q - target) > 1e-10:
                raise ValueError("q does not look like a quaternion-embedded 2x2 matrix.")

        X = np.zeros((32, 32), dtype=complex)

        # --- particles (0..15) ---
        # Q_L block: (uL_r,uL_g,uL_b,dL_r,dL_g,dL_b)
        idx_QL = self.uL + self.dL
        X[np.ix_(idx_QL, idx_QL)] = np.kron(q, m)

        # L_L block: (nuL,eL)
        idx_LL = [self.nuL[0], self.eL[0]]
        X[np.ix_(idx_LL, idx_LL)] = q

        # uR, dR: lam*m on each color triplet
        X[np.ix_(self.uR, self.uR)] = lam * m
        X[np.ix_(self.dR, self.dR)] = lam * m

        # nuR, eR: lam (singlets)
        X[self.nuR[0], self.nuR[0]] = lam
        X[self.eR[0], self.eR[0]] = lam

        # --- conjugates (16..31): use conjugate parameters ---
        lamc = np.conj(lam)
        qc = np.conj(q)
        mc = np.conj(m)

        idx_QLc = self.uL_c + self.dL_c
        X[np.ix_(idx_QLc, idx_QLc)] = np.kron(qc, mc)

        idx_LLc = [self.nuL_c[0], self.eL_c[0]]
        X[np.ix_(idx_LLc, idx_LLc)] = qc

        X[np.ix_(self.uR_c, self.uR_c)] = lamc * mc
        X[np.ix_(self.dR_c, self.dR_c)] = lamc * mc

        X[self.nuR_c[0], self.nuR_c[0]] = lamc
        X[self.eR_c[0], self.eR_c[0]] = lamc

        return X

    def pi_lifted(self, lam: complex, q: np.ndarray, m: np.ndarray, N: int, *, enforce_quaternion: bool = False) -> np.ndarray:
        """
        Lift π(a) to H_SM ⊗ H_flav: π(a) ⊗ I_N.
        """
        Pi32 = self.pi(lam, q, m, enforce_quaternion=enforce_quaternion)
        return np.kron(Pi32, np.eye(N, dtype=complex))

# ============================================================
#  One-forms Ω^1_D(A) and inner fluctuations
# ============================================================

class InnerFluctuations:
    """
    Build one-forms A = Σ π(a_i)[D, π(b_i)] and D_A = D + A + JAJ^{-1}.

    Hermitization options:
      - hermitize="A":        replace A by (A+A†)/2 before adding JAJ^{-1}
      - hermitize="total":    hermitize the full fluctuation term (A + JAJ^{-1})  [recommended]
      - hermitize=False:      do nothing (debug only)
    """

    def __init__(self, J_total: AntiUnitary) -> None:
        self.J_total = J_total

    # ---------------- core primitives ----------------

    @staticmethod
    def _require_same_shape(name: str, X: np.ndarray, Y: np.ndarray) -> None:
        if X.shape != Y.shape:
            raise ValueError(f"{name}: shape mismatch {X.shape} vs {Y.shape}")

    @staticmethod
    def _require_square(name: str, X: np.ndarray) -> None:
        if X.ndim != 2 or X.shape[0] != X.shape[1]:
            raise ValueError(f"{name} must be square, got {X.shape}")

    def one_form(
        self,
        D: np.ndarray,
        pi_a: np.ndarray,
        pi_b: np.ndarray,
        *,
        remove_trace: bool = False,
    ) -> np.ndarray:
        """
        A(a,b) = π(a)[D,π(b)].

        remove_trace=True subtracts (tr(A)/n) I to keep a traceless contribution (often useful
        when you want to drop U(1) pieces for diagnostics).
        """
        D = np.asarray(D, dtype=complex)
        pa = np.asarray(pi_a, dtype=complex)
        pb = np.asarray(pi_b, dtype=complex)

        self._require_square("D", D)
        self._require_square("pi_a", pa)
        self._require_square("pi_b", pb)
        self._require_same_shape("one_form inputs", D, pa)
        self._require_same_shape("one_form inputs", D, pb)

        A = pa @ comm(D, pb)

        if remove_trace:
            n = A.shape[0]
            tr = np.trace(A) / complex(n)
            A = A - tr * np.eye(n, dtype=complex)

        return A

    def build_A(
        self,
        D: np.ndarray,
        pairs: Sequence[Tuple[np.ndarray, np.ndarray]],
        *,
        weights: Sequence[complex] | None = None,
        remove_trace: bool = False,
    ) -> np.ndarray:
        """
        Sum A = Σ w_i π(a_i)[D,π(b_i)].
        """
        D = np.asarray(D, dtype=complex)
        self._require_square("D", D)

        if weights is not None and len(weights) != len(pairs):
            raise ValueError(f"weights must have same length as pairs: {len(weights)} vs {len(pairs)}")

        A = np.zeros_like(D, dtype=complex)
        for idx, (pa, pb) in enumerate(pairs):
            w = 1.0 if weights is None else complex(weights[idx])
            A += w * self.one_form(D, pa, pb, remove_trace=remove_trace)
        return A

    def fluctuate(
        self,
        D: np.ndarray,
        A: np.ndarray,
        *,
        hermitize: str | bool = "total",
        stabilize_output: bool = True,
    ) -> np.ndarray:
        """
        D_A = D + A + JAJ^{-1}, with optional hermitization policy.

        stabilize_output=True forces final D_A -> (D_A + D_A†)/2 to suppress tiny numerical skew.
        """
        D = np.asarray(D, dtype=complex)
        A = np.asarray(A, dtype=complex)

        self._require_square("D", D)
        self._require_square("A", A)
        self._require_same_shape("fluctuate inputs", D, A)

        if hermitize is True:
            hermitize = "A"
        if hermitize not in (False, "A", "total"):
            raise ValueError("hermitize must be False, 'A', or 'total'.")

        if hermitize == "A":
            A_use = hermitian_part(A)
            DA = D + A_use + self.J_total.conj_op(A_use)
        else:
            JA = self.J_total.conj_op(A)
            total_fluct = A + JA
            if hermitize == "total":
                total_fluct = hermitian_part(total_fluct)
            DA = D + total_fluct

        return hermitian_part(DA) if stabilize_output else DA

    # ---------------- convenience ----------------

    def fluctuate_from_pairs(
        self,
        D: np.ndarray,
        pairs: Sequence[Tuple[np.ndarray, np.ndarray]],
        *,
        weights: Sequence[complex] | None = None,
        remove_trace: bool = False,
        hermitize: str | bool = "total",
        stabilize_output: bool = True,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Convenience: build A from pairs then return (D_A, A).
        """
        A = self.build_A(D, pairs, weights=weights, remove_trace=remove_trace)
        D_A = self.fluctuate(D, A, hermitize=hermitize, stabilize_output=stabilize_output)
        return D_A, A

    # ---------------- diagnostics ----------------

    def check_self_adjoint(self, X: np.ndarray) -> float:
        X = np.asarray(X, dtype=complex)
        self._require_square("X", X)
        return float(op_norm(X - X.conj().T))

    def check_oddness(self, X: np.ndarray, Gamma_total: np.ndarray) -> float:
        X = np.asarray(X, dtype=complex)
        G = np.asarray(Gamma_total, dtype=complex)
        self._require_square("X", X)
        self._require_square("Gamma_total", G)
        self._require_same_shape("oddness inputs", X, G)
        return float(op_norm(anti_comm(G, X)))

    def check_commutator(self, D: np.ndarray, pi: np.ndarray) -> float:
        D = np.asarray(D, dtype=complex)
        pi = np.asarray(pi, dtype=complex)
        self._require_square("D", D)
        self._require_square("pi", pi)
        self._require_same_shape("commutator inputs", D, pi)
        return float(op_norm(comm(D, pi)))

# ============================================================
#  v4.0 production gates (order-zero, first-order, stability)
# ============================================================

@dataclass(frozen=True)
class GateReport:
    order0: float
    first: float
    self_adj: float
    odd_int: float
    # optional extras
    order0_pair: Optional[tuple[int, int]] = None
    first_pair: Optional[tuple[int, int]] = None




class V4Gates:
    """
    Implements Section 10 production commutator gates.

    Gate definitions:
      - order-0:  max_{a,b in Π} || [π(a), J π(b) J^{-1}] ||
      - first-order: max_{a,b in Π} || [[D, π(a)], J π(b) J^{-1}] ||
      - self-adjoint: ||X - X†|| (for whatever X you're checking)
      - oddness: ||{Γ, D}|| or ||{γ_int, D_int}|| depending on what you pass in

    Production details:
      * Precomputes Jπ(b)J^{-1} and [D,π(a)] once.
      * Optional early-exit if a gate exceeds its epsilon.
      * Optional sampling over Π for quick runs.
    """

    def __init__(self, eps0: float, eps1: float, eps_sa: float, epsA: float) -> None:
        self.eps0 = float(eps0)
        self.eps1 = float(eps1)
        self.eps_sa = float(eps_sa)
        self.epsA = float(epsA)

    # ---------- basic checks ----------

    @staticmethod
    def _require_square(name: str, X: np.ndarray) -> None:
        if X.ndim != 2 or X.shape[0] != X.shape[1]:
            raise ValueError(f"{name} must be square, got {X.shape}")

    @staticmethod
    def _require_same_dim(name: str, A: np.ndarray, B: np.ndarray) -> None:
        if A.shape != B.shape:
            raise ValueError(f"{name}: shape mismatch {A.shape} vs {B.shape}")

    def gate_selfadj(self, X: np.ndarray) -> float:
        X = np.asarray(X, dtype=complex)
        self._require_square("X", X)
        return float(op_norm(X - X.conj().T))

    def _sanitize_sample(self, n: int, sample: Optional[Sequence[int]]) -> List[int]:
        """
        Normalize and validate sample indices for a list of length n.
        - Supports negative indices like Python (-1 = last).
        - Drops out-of-range entries.
        - Preserves order and removes duplicates.
        """
        if n <= 0:
            raise ValueError("pis must be non-empty")

        if sample is None:
            return list(range(n))

        out: List[int] = []
        seen = set()
        for s in sample:
            i = int(s)
            if -n <= i < 0:
                i = i + n
            if 0 <= i < n and i not in seen:
                out.append(i)
                seen.add(i)

        if not out:
            raise ValueError(f"sample has no valid indices for pis length {n}. "
                             f"Valid range is 0..{n - 1} (or -{n}..-1). Got: {list(sample)}")
        return out

    # ---------- Section 10 gates ----------

    def gate_order0(
        self,
        pis: List[np.ndarray],
        J_total: AntiUnitary,
        *,
        early_exit: bool = False,
        eps: Optional[float] = None,
        sample: Optional[Sequence[int]] = None,
        return_witness: bool = False,
    ) -> float | tuple[float, tuple[int, int]]:
        """
        max_{a,b} || [a, J(b)] || where J(b)=J_total.conj_op(b)

        sample: optional list of indices into pis to restrict the sweep.
        """
        if len(pis) == 0:
            raise ValueError("pis must be non-empty")

        n = len(pis)
        idxs = self._sanitize_sample(n, sample)

        mats = [np.asarray(pis[i], dtype=complex) for i in idxs]
        for k, A in enumerate(mats):
            self._require_square(f"pis[{idxs[k]}]", A)
            if k > 0:
                self._require_same_dim("pis", mats[0], A)

        # Precompute J(b)
        Jpis = [J_total.conj_op(B) for B in mats]

        thresh = self.eps0 if eps is None else float(eps)
        mx = 0.0
        w = (idxs[0], idxs[0])

        for i, A in enumerate(mats):
            for j, JB in enumerate(Jpis):
                v = float(op_norm(comm(A, JB)))
                if v > mx:
                    mx = v
                    w = (idxs[i], idxs[j])
                    if early_exit and mx > thresh:
                        return (mx, w) if return_witness else mx

        return (mx, w) if return_witness else mx

    def gate_first(
        self,
        D: np.ndarray,
        pis: List[np.ndarray],
        J_total: AntiUnitary,
        *,
        early_exit: bool = False,
        eps: Optional[float] = None,
        sample: Optional[Sequence[int]] = None,
        return_witness: bool = False,
    ) -> float | tuple[float, tuple[int, int]]:
        """
        max_{a,b} || [[D,a], J(b)] ||.
        """
        D = np.asarray(D, dtype=complex)
        self._require_square("D", D)

        if len(pis) == 0:
            raise ValueError("pis must be non-empty")

        n = len(pis)
        idxs = self._sanitize_sample(n, sample)

        mats = [np.asarray(pis[i], dtype=complex) for i in idxs]
        for k, A in enumerate(mats):
            self._require_square(f"pis[{idxs[k]}]", A)
            self._require_same_dim("D vs pi", D, A)

        # Precompute [D, a] and J(b)
        Das = [comm(D, A) for A in mats]
        Jpis = [J_total.conj_op(B) for B in mats]

        thresh = self.eps1 if eps is None else float(eps)
        mx = 0.0
        w = (idxs[0], idxs[0])

        for i, Da in enumerate(Das):
            for j, JB in enumerate(Jpis):
                v = float(op_norm(comm(Da, JB)))
                if v > mx:
                    mx = v
                    w = (idxs[i], idxs[j])
                    if early_exit and mx > thresh:
                        return (mx, w) if return_witness else mx

        return (mx, w) if return_witness else mx

    # ---------- report helpers ----------

    def build_report(
        self,
        *,
        order0: float,
        first: float,
        self_adj: float,
        odd_int: float,
        order0_pair: Optional[tuple[int, int]] = None,
        first_pair: Optional[tuple[int, int]] = None,
    ) -> GateReport:
        return GateReport(
            order0=float(order0),
            first=float(first),
            self_adj=float(self_adj),
            odd_int=float(odd_int),
            order0_pair=order0_pair,
            first_pair=first_pair,
        )

    def pass_fail(self, report: GateReport) -> bool:
        # oddness is a hard axiom check: use epsA as the “hard” tolerance unless you prefer eps_sa
        return (
            report.order0 <= self.eps0 and
            report.first <= self.eps1 and
            report.self_adj <= self.eps_sa and
            report.odd_int <= self.epsA
        )
# ============================================================
#  Texture generation hooks (your current v3.3 kernel pipeline)
# ============================================================

class TextureGenerator:
    """
    Wrap your existing K→flow→P_C360→S compression machinery as a texture generator.

    Supported cfg styles (pick one):

    1) cfg.builder: callable(params) -> TexturePack | dict with keys Yu,Yd,Ye,Ynu,MR
    2) cfg.npz_path: path to .npz containing arrays named Yu,Yd,Ye,Ynu,MR
    3) cfg.mode: one of:
         - "identity_leptons" (default): Yu=Yd=0, Ye=I, Ynu=I, MR=seesaw*I
         - "identity_all": Yu=Yd=Ye=Ynu=I, MR=seesaw*I
         - "zeros_all": Yu=Yd=Ye=Ynu=0, MR=None
         - "random": random complex textures (seedable), MR symmetric

    Common cfg knobs:
      - disable_quarks: bool (forces Yu=Yd=0)
      - scales: dict with optional keys {"Yu","Yd","Ye","Ynu","MR"}
      - enforce_MR_symmetric: bool (default True)
      - seed: int (for mode="random")
      - normalize: bool (rescale each Yukawa so ||Y||_2 ~= 1 before scales)
    """

    def __init__(self, cfg: Any) -> None:
        self.cfg = cfg

    # ----------------- helpers -----------------

    @staticmethod
    def _as_mat(name: str, X: Any, N: int, allow_none: bool = False) -> np.ndarray | None:
        if X is None:
            if allow_none:
                return None
            raise ValueError(f"{name} is None, expected {N}x{N} matrix.")
        A = np.asarray(X, dtype=complex)
        if A.shape != (N, N):
            raise ValueError(f"{name} must be shape ({N},{N}), got {A.shape}.")
        return A

    @staticmethod
    def _symmetrize_majorana(M: np.ndarray) -> np.ndarray:
        # Majorana mass matrix is typically complex symmetric: M = M^T.
        return 0.5 * (M + M.T)

    @staticmethod
    def _opnorm_rescale(Y: np.ndarray, target: float = 1.0, eps: float = 1e-14) -> np.ndarray:
        n = op_norm(Y)
        if n < eps:
            return Y
        return (target / n) * Y

    def _get_cfg(self, key: str, default: Any = None) -> Any:
        # Supports dict-like cfg or object-with-attrs cfg
        if isinstance(self.cfg, dict):
            return self.cfg.get(key, default)
        return getattr(self.cfg, key, default)

    # ----------------- main entry -----------------

    def build_textures(self, params: V40Params) -> TexturePack:
        N = int(params.N_flav)

        # Highest priority: a callable builder (your v3.3 pipeline hook)
        builder = self._get_cfg("builder", None)
        if builder is not None:
            out = builder(params)
            if isinstance(out, TexturePack):
                tp = out
            elif isinstance(out, dict):
                tp = TexturePack(
                    Yu=out.get("Yu"), Yd=out.get("Yd"), Ye=out.get("Ye"),
                    Ynu=out.get("Ynu"), MR=out.get("MR")
                )
            else:
                raise ValueError("cfg.builder must return TexturePack or dict.")
        else:
            # Next: load from npz if provided
            npz_path = self._get_cfg("npz_path", None)
            if npz_path is not None:
                data = np.load(npz_path, allow_pickle=False)
                tp = TexturePack(
                    Yu=data.get("Yu"), Yd=data.get("Yd"), Ye=data.get("Ye"),
                    Ynu=data.get("Ynu"), MR=data.get("MR", None)
                )
            else:
                # Fallback: built-in modes
                mode = self._get_cfg("mode", "identity_leptons")

                if mode == "identity_leptons":
                    Yu = np.zeros((N, N), dtype=complex)
                    Yd = np.zeros((N, N), dtype=complex)
                    Ye = np.eye(N, dtype=complex)
                    Ynu = np.eye(N, dtype=complex)
                    MR = complex(params.seesaw_M_GeV) * np.eye(N, dtype=complex)
                elif mode == "identity_all":
                    Yu = np.eye(N, dtype=complex)
                    Yd = np.eye(N, dtype=complex)
                    Ye = np.eye(N, dtype=complex)
                    Ynu = np.eye(N, dtype=complex)
                    MR = complex(params.seesaw_M_GeV) * np.eye(N, dtype=complex)
                elif mode == "zeros_all":
                    Yu = np.zeros((N, N), dtype=complex)
                    Yd = np.zeros((N, N), dtype=complex)
                    Ye = np.zeros((N, N), dtype=complex)
                    Ynu = np.zeros((N, N), dtype=complex)
                    MR = None
                elif mode == "random":
                    seed = self._get_cfg("seed", 0)
                    rng = np.random.default_rng(seed)
                    def rc():
                        return rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N))
                    Yu, Yd, Ye, Ynu = rc(), rc(), rc(), rc()
                    MR = rc()
                    MR = self._symmetrize_majorana(MR)
                    MR *= complex(params.seesaw_M_GeV) / max(op_norm(MR), 1e-14)
                else:
                    raise ValueError(f"Unknown cfg.mode={mode!r}")

                tp = TexturePack(Yu=Yu, Yd=Yd, Ye=Ye, Ynu=Ynu, MR=MR)

        # Validate shapes + dtype
        Yu = self._as_mat("Yu", tp.Yu, N)
        Yd = self._as_mat("Yd", tp.Yd, N)
        Ye = self._as_mat("Ye", tp.Ye, N)
        Ynu = self._as_mat("Ynu", tp.Ynu, N)
        MR = self._as_mat("MR", tp.MR, N, allow_none=True)

        # Optional toggles / post-processing
        if bool(self._get_cfg("disable_quarks", False)):
            Yu = np.zeros_like(Yu)
            Yd = np.zeros_like(Yd)

        if bool(self._get_cfg("normalize", False)):
            Yu = self._opnorm_rescale(Yu)
            Yd = self._opnorm_rescale(Yd)
            Ye = self._opnorm_rescale(Ye)
            Ynu = self._opnorm_rescale(Ynu)

        scales = self._get_cfg("scales", {}) or {}
        Yu *= complex(scales.get("Yu", 1.0))
        Yd *= complex(scales.get("Yd", 1.0))
        Ye *= complex(scales.get("Ye", 1.0))
        Ynu *= complex(scales.get("Ynu", 1.0))
        if MR is not None:
            MR *= complex(scales.get("MR", 1.0))

        if MR is not None and bool(self._get_cfg("enforce_MR_symmetric", True)):
            MR = self._symmetrize_majorana(MR)

        return TexturePack(Yu=Yu, Yd=Yd, Ye=Ye, Ynu=Ynu, MR=MR)


# ============================================================
#  v4.0 Engine: orchestrates everything
# ============================================================

class AlignmentV40Engine:
    """
    One-stop v4.0 engine:
      - builds textures
      - builds D_int
      - builds even-product D_total
      - optionally builds fluctuations D_A
      - runs axioms + commutator gates
      - exposes sector extraction utilities
    """

    def __init__(
        self,
        geom: GeometryFactor,
        params: V40Params,
        texture_gen: TextureGenerator,
    ) -> None:
        self.params = params
        self.geom = geom
        self.texture_gen = texture_gen

        # fixed structures
        self.basis = SMBasis32()
        self.maps = ChannelMaps(self.basis)

        # internal real structures
        self.J_SM = RealStructureFactory.build_J_SM(self.basis)
        self.J_flav = RealStructureFactory.build_J_flav(params.N_flav)
        self.J_int = RealStructureFactory.kron(self.J_SM, self.J_flav)

        # total J on H_geom ⊗ H_int
        self.J_total = RealStructureFactory.kron(self.geom.J_geom, self.J_int)

        # builders
        self.dint_builder = DIntBuilder(self.basis, self.maps)
        self.product_dirac = EvenProductDirac(self.geom, self.basis, params)

        # representation layer
        self.pi_SM = SMRepresentation(self.basis)

        # gates
        self.gates = V4Gates(
            eps0=params.eps_order0,
            eps1=params.eps_first,
            eps_sa=params.eps_sa,
            epsA=params.eps_fluct,   # use as hard oddness tolerance too
        )

    # -------- build core objects --------

    def build_textures(self) -> TexturePack:
        return self.texture_gen.build_textures(self.params)

    def build_D_int(self, textures: TexturePack) -> np.ndarray:
        return self.dint_builder.build_D_int(
            textures,
            J_SM=self.J_SM,
            N=self.params.N_flav,
        )

    def build_D(self, D_int: np.ndarray) -> np.ndarray:
        return self.product_dirac.build_D(D_int)

    def build_Gamma_total(self) -> np.ndarray:
        return self.product_dirac.build_Gamma_total(self.params.N_flav)

    # -------- fluctuations --------

    def build_fluctuated(
        self,
        D: np.ndarray,
        one_form_pairs: Sequence[Tuple[np.ndarray, np.ndarray]],
        *,
        weights: Sequence[complex] | None = None,
        hermitize: str | bool = "total",
        remove_trace: bool = False,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Returns (D_A, A).
        """
        fl = InnerFluctuations(self.J_total)
        A = fl.build_A(D, one_form_pairs, weights=weights, remove_trace=remove_trace)
        D_A = fl.fluctuate(D, A, hermitize=hermitize, stabilize_output=True)
        return D_A, A

    # -------- diagnostics / axioms --------

    def run_axiom_checks_internal(self, D_int: np.ndarray) -> Dict[str, float]:
        gamma_int = self.basis.lift_gamma(self.params.N_flav)
        sa = op_norm(D_int - D_int.conj().T)
        odd_int = op_norm(anti_comm(gamma_int, D_int))
        return {"sa": float(sa), "odd": float(odd_int), "sa_ok": sa <= self.params.eps_sa, "odd_ok": odd_int <= self.params.eps_fluct}

    def run_axiom_checks_total(self, D_total: np.ndarray) -> Dict[str, float]:
        Gamma = self.build_Gamma_total()
        sa = op_norm(D_total - D_total.conj().T)
        odd_total = op_norm(anti_comm(Gamma, D_total))
        return {"sa_total": float(sa), "odd_total": float(odd_total), "sa_ok": sa <= self.params.eps_sa, "odd_ok": odd_total <= self.params.eps_fluct}

    # -------- commutator gates --------

    def run_commutator_gates(
        self,
        D_total: np.ndarray,
        *,
        return_witness: bool = True,
        early_exit: bool = False,
        sample: Sequence[int] | None = None,
        odd_int_value: float | None = None,
    ) -> GateReport:
        pis = self._lift_pi_generators_to_full()

        if return_witness:
            order0, w0 = self.gates.gate_order0(pis, self.J_total, return_witness=True, early_exit=early_exit, sample=sample)
            first, w1 = self.gates.gate_first(D_total, pis, self.J_total, return_witness=True, early_exit=early_exit, sample=sample)
        else:
            order0 = self.gates.gate_order0(pis, self.J_total, early_exit=early_exit, sample=sample)
            first = self.gates.gate_first(D_total, pis, self.J_total, early_exit=early_exit, sample=sample)
            w0 = None
            w1 = None

        sa = self.gates.gate_selfadj(D_total)

        # oddness is fundamentally internal; accept caller value or compute from D_int elsewhere.
        odd_int = float(odd_int_value) if odd_int_value is not None else 0.0

        return GateReport(
            order0=float(order0),
            first=float(first),
            self_adj=float(sa),
            odd_int=float(odd_int),
            order0_pair=w0 if return_witness else None,
            first_pair=w1 if return_witness else None,
        )

    def _lift_pi_generators_to_full(self) -> List[np.ndarray]:
        """
        Lift π(A) generators to full H_geom ⊗ H_SM ⊗ H_flav
        (π acts trivially on flavor and geometry).
        """
        N = self.params.N_flav
        gens_SM = self.pi_SM.generators()  # 32x32
        gens_int = [np.kron(g, np.eye(N, dtype=complex)) for g in gens_SM]  # 32N x 32N

        I_geom = np.eye(self.geom.dim(), dtype=complex)
        return [np.kron(I_geom, g_int) for g_int in gens_int]

    # -------- sector extraction utilities --------

    def sector_projector_int(self, key: str) -> np.ndarray:
        """(32N x 32N) projector onto an internal sector key."""
        P32 = self.basis.P[key]
        return self.basis.lift_projector(P32, self.params.N_flav)

    def sector_projector_full(self, key: str) -> np.ndarray:
        """(dim_geom*32N x dim_geom*32N) projector onto a sector in the full product space."""
        I_geom = np.eye(self.geom.dim(), dtype=complex)
        return np.kron(I_geom, self.sector_projector_int(key))

    def extract_lr_block_int(self, D_int: np.ndarray, left_key: str, right_key: str) -> np.ndarray:
        """P_L D_int P_R on H_SM ⊗ H_flav."""
        PL = self.sector_projector_int(left_key)
        PR = self.sector_projector_int(right_key)
        return PL @ D_int @ PR

    def extract_lr_block_full(self, D_total: np.ndarray, left_key: str, right_key: str) -> np.ndarray:
        """(I⊗P_L) D_total (I⊗P_R) on H_geom ⊗ H_int."""
        PL = self.sector_projector_full(left_key)
        PR = self.sector_projector_full(right_key)
        return PL @ D_total @ PR

    # -------- end-to-end run --------

    def run(
        self,
        *,
        do_fluctuate: bool = False,
        one_form_pairs: Sequence[Tuple[np.ndarray, np.ndarray]] | None = None,
        fluct_weights: Sequence[complex] | None = None,
        hermitize: str | bool = "total",
        gate_sample: Sequence[int] | None = None,
        early_exit_gates: bool = False,
    ) -> Dict[str, Any]:
        textures = self.build_textures()
        D_int = self.build_D_int(textures)
        D_total = self.build_D(D_int)

        ax_int = self.run_axiom_checks_internal(D_int)
        ax_tot = self.run_axiom_checks_total(D_total)

        D_used = D_total
        A = None
        if do_fluctuate:
            if one_form_pairs is None:
                raise ValueError("do_fluctuate=True requires one_form_pairs.")
            D_used, A = self.build_fluctuated(
                D_total,
                one_form_pairs,
                weights=fluct_weights,
                hermitize=hermitize,
            )

        gates = self.run_commutator_gates(
            D_used,
            return_witness=True,
            early_exit=early_exit_gates,
            sample=gate_sample,
            odd_int_value=ax_int["odd"],
        )

        return {
            "textures": textures,
            "D_int": D_int,
            "D_total": D_total,
            "A": A,
            "D_used": D_used,
            "axioms_internal": ax_int,
            "axioms_total": ax_tot,
            "gates": gates,
            "pass": self.gates.pass_fail(gates),
        }


# ============================================================
#  Optional: Spectral action evaluators (conceptual / numerical)
# ============================================================

class SpectralAction:
    """
    Bosonic spectral action approximations on finite truncations.

    Supported evaluation modes:
      1) exact diagonalization (dense):     Tr f(D^2/Λ^2)
      2) partial eigenvalues (k largest):   good if f decays fast and you only need top modes
      3) stochastic trace (Hutchinson):    Tr f(D^2/Λ^2) ≈ (1/M) Σ z^T f(D^2/Λ^2) z
         using Chebyshev polynomial approximation for f on [0, R].

    Notes:
      - For stability we always hermitize D first.
      - We prefer working with D^2 as a PSD operator when possible.
    """

    def __init__(
        self,
        cutoff_Lambda: float,
        f: Callable[[np.ndarray], np.ndarray],
        *,
        atol_sa: float = 1e-10,
    ) -> None:
        self.Lambda = float(cutoff_Lambda)
        self.f = f
        self.atol_sa = float(atol_sa)

    # ---------------- basics ----------------

    @staticmethod
    def _require_square(name: str, X: np.ndarray) -> None:
        if X.ndim != 2 or X.shape[0] != X.shape[1]:
            raise ValueError(f"{name} must be square, got {X.shape}")

    def _herm(self, D: np.ndarray) -> np.ndarray:
        D = np.asarray(D, dtype=complex)
        self._require_square("D", D)
        return hermitian_part(D)

    # ---------------- exact trace ----------------

    def trace_f(self, D: np.ndarray) -> float:
        """
        Exact: Tr f(D^2/Λ^2) by full diagonalization of hermitian(D).
        """
        Dh = self._herm(D)
        evals = np.linalg.eigvalsh(Dh)  # real
        x = (evals.astype(float) ** 2) / (self.Lambda ** 2)
        fx = np.asarray(self.f(x))
        return float(np.sum(fx))

    def trace_f_from_eigs(self, eigs: np.ndarray) -> float:
        """
        If you already have eigenvalues of hermitian(D), reuse them.
        """
        eigs = np.asarray(eigs, dtype=float).reshape(-1)
        x = (eigs ** 2) / (self.Lambda ** 2)
        return float(np.sum(self.f(x)))

    # ---------------- partial spectrum (optional) ----------------

    def trace_f_topk(self, D: np.ndarray, k: int) -> float:
        """
        Approximate Tr f(D^2/Λ^2) using only the k largest-magnitude eigenvalues of D.
        Useful only if f(x) decays quickly and truncation is justified.
        Dense fallback implementation (still O(n^3) worst-case), but isolates the API.

        If you want true sparse partial eigensolve, plug scipy.sparse.linalg.eigsh in here.
        """
        Dh = self._herm(D)
        n = Dh.shape[0]
        k = int(k)
        if k <= 0 or k > n:
            raise ValueError(f"k must be in 1..{n}, got {k}")

        # Dense "partial" via full eigvals then select; replace with sparse eigsh in production.
        evals = np.linalg.eigvalsh(Dh)
        # select k largest |λ|
        idx = np.argsort(np.abs(evals))[::-1][:k]
        sel = evals[idx]
        x = (sel.astype(float) ** 2) / (self.Lambda ** 2)
        return float(np.sum(self.f(x)))

    # ---------------- fermionic term ----------------

    def fermionic_term(self, psi: np.ndarray, D: np.ndarray, *, normalize: bool = False) -> complex:
        """
        <psi, D psi>. If normalize=True, uses psi/||psi||.
        """
        Dh = self._herm(D)
        psi = np.asarray(psi, dtype=complex).reshape(-1)
        if psi.shape[0] != Dh.shape[0]:
            raise ValueError(f"psi length {psi.shape[0]} does not match D dim {Dh.shape[0]}")
        if normalize:
            nrm = np.linalg.norm(psi)
            if nrm > 0:
                psi = psi / nrm
        return complex(np.vdot(psi, Dh @ psi))

    # ---------------- stochastic trace (Hutchinson) ----------------

    def trace_f_hutchinson(
        self,
        D: np.ndarray,
        *,
        num_samples: int = 64,
        rng_seed: int = 0,
    ) -> float:
        """
        Hutchinson estimator for Tr f(D^2/Λ^2) using *exact* f(D^2/Λ^2) application
        via diagonalization of D (so this is still expensive but shows the interface).

        In production, you’d replace the body with:
          - build a LinearOperator for A = D^2/Λ^2
          - approximate f(A)z via Chebyshev / Lanczos
        """
        Dh = self._herm(D)
        evals, U = np.linalg.eigh(Dh)  # Dh = U diag(evals) U†
        x = (evals.astype(float) ** 2) / (self.Lambda ** 2)
        fx = np.asarray(self.f(x), dtype=float)

        rng = np.random.default_rng(rng_seed)
        n = Dh.shape[0]
        M = int(num_samples)
        if M <= 0:
            raise ValueError("num_samples must be > 0")

        # z are Rademacher vectors for Hutchinson
        acc = 0.0
        for _ in range(M):
            z = rng.choice([-1.0, 1.0], size=(n,))
            # compute z^T f(A) z where f(A)=U diag(fx) U†
            Uz = U.conj().T @ z
            acc += float(np.sum(fx * (np.abs(Uz) ** 2)))
        return float(acc / M)

# ============================================================
#  Minimal self-tests (fast numerical sanity)
# ============================================================

def _self_test_sm_structures() -> None:
    b = SMBasis32()

    # gamma_SM sanity
    assert is_hermitian(b.gamma_SM)
    assert np.allclose(b.gamma_SM @ b.gamma_SM, np.eye(32), atol=1e-12)

    # J_SM swap sanity
    J = RealStructureFactory.build_J_SM(b)
    U = J.U
    assert np.allclose(U.conj().T @ U, np.eye(32), atol=1e-12)
    # swap property on basis vectors (up to exact equality in this basis)
    e0 = np.zeros((32,), dtype=complex); e0[0] = 1.0
    e16 = np.zeros((32,), dtype=complex); e16[16] = 1.0
    assert np.allclose(U @ e0, e16)
    assert np.allclose(U @ e16, e0)

    # ChannelMaps partial-isometry diagnostics
    maps = ChannelMaps(b)
    diags = maps.validate()
    # they should be exactly 0 in floating arithmetic for these integer matrices
    for k, v in diags.items():
        assert v < 1e-12, f"{k} deviation too large: {v}"

    # Representation generator sanity
    rep = SMRepresentation(b)
    gens = rep.generators()
    assert all(g.shape == (32, 32) for g in gens)

def _self_test_dirac_builders() -> None:
    # tiny sizes for speed
    N = 2

    # SM core
    b = SMBasis32()
    maps = ChannelMaps(b)
    J_SM = RealStructureFactory.build_J_SM(b)

    # textures: simple deterministic (nontrivial but stable)
    Yu = np.zeros((N, N), dtype=complex)
    Yd = np.zeros((N, N), dtype=complex)
    Ye = np.eye(N, dtype=complex)
    Ynu = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)  # swap
    MR = 1e9 * np.eye(N, dtype=complex)
    textures = TexturePack(Yu=Yu, Yd=Yd, Ye=Ye, Ynu=Ynu, MR=MR)

    # Build D_int
    dint = DIntBuilder(b, maps).build_D_int(textures, J_SM=J_SM, N=N)
    assert dint.shape == (32 * N, 32 * N)
    assert is_hermitian(dint)

    # Internal oddness: {gamma_SM⊗I, D_int} = 0
    gamma_int = b.lift_gamma(N)
    odd_int = op_norm(anti_comm(gamma_int, dint))
    assert odd_int < 1e-10, f"odd_int too large: {odd_int}"

    # Geometry: use even truncation so {gamma_geom, D_geom}=0 by construction
    geom = GeometryFactor.even_truncation(lambdas=np.array([1.0, 2.0]))

    # Product Dirac
    class _P:  # minimal params shim
        eps_sa = 1e-10
        eps_fluct = 1e-10
    prod = EvenProductDirac(geom, b, _P())

    Dtot = prod.build_D(dint)
    assert Dtot.shape == (geom.dim() * 32 * N, geom.dim() * 32 * N)
    assert is_hermitian(Dtot)

    # Total oddness: {Gamma_total, D_total} = 0
    Gamma = prod.build_Gamma_total(N)
    odd_tot = op_norm(anti_comm(Gamma, Dtot))
    assert odd_tot < 1e-10, f"odd_total too large: {odd_tot}"

def _self_test_engine_smoke() -> None:
    # small even geometry
    geom = GeometryFactor.even_truncation(lambdas=np.array([1.0, 3.0]))

    # minimal params object (match the fields your engine uses)
    class _Params:
        N_flav = 2
        seesaw_M_GeV = 1e9

        eps_order0 = 1e-10
        eps_first = 1e-10
        eps_sa = 1e-10
        eps_fluct = 1e-10

    params = _Params()

    # textures: deterministic identity leptons (fast)
    texgen = TextureGenerator(cfg={"mode": "identity_leptons"})

    eng = AlignmentV40Engine(geom=geom, params=params, texture_gen=texgen)

    # build core
    textures = eng.build_textures()
    D_int = eng.build_D_int(textures)
    D_tot = eng.build_D(D_int)

    # internal axioms
    ax_int = eng.run_axiom_checks_internal(D_int)
    assert ax_int["sa_ok"] and ax_int["odd_ok"], f"internal axioms failed: {ax_int}"

    # total axioms
    ax_tot = eng.run_axiom_checks_total(D_tot)
    assert ax_tot["sa_ok"] and ax_tot["odd_ok"], f"total axioms failed: {ax_tot}"

    # commutator gates (sample a few generators for speed)
    rep = eng.run_commutator_gates(
        D_tot,
        return_witness=True,
        early_exit=True,
        sample=[0, 1, 2, 3],          # I, Y, first SU2 gens...
        odd_int_value=ax_int["odd"],
    )
    assert np.isfinite(rep.order0) and np.isfinite(rep.first) and np.isfinite(rep.self_adj)

    # optional: do one tiny fluctuation using a couple of lifted generators
    pis = eng._lift_pi_generators_to_full()
    pairs = [(pis[1], pis[2]), (pis[2], pis[1])]  # simple symmetric-ish pair set
    D_A, A = eng.build_fluctuated(D_tot, pairs, hermitize="total")
    assert D_A.shape == D_tot.shape
    assert is_hermitian(D_A)

def run_self_tests() -> None:
    _self_test_sm_structures()
    _self_test_dirac_builders()
    _self_test_engine_smoke()

# ============================================================
#  Example wiring (scaffold)
# ============================================================

def example_build_engine() -> AlignmentV40Engine:
    # Params: assumes your V40Params has these defaults; override eps if desired.
    params = V40Params(N_flav=3)

    # Use an even truncation so {D_geom, gamma_geom}=0 holds.
    # lambdas are the geometric mode eigenvalues in the doubled construction.
    geom = GeometryFactor.even_truncation(lambdas=np.array([1, 2, 3, 4, 5], dtype=float))

    # Default stub textures (identity leptons) unless cfg says otherwise.
    texgen = TextureGenerator(cfg={"mode": "identity_leptons"})

    return AlignmentV40Engine(geom=geom, params=params, texture_gen=texgen)


def main() -> None:
    eng = example_build_engine()

    out = eng.run(
        do_fluctuate=False,
        gate_sample=[0, 1, 2, 3],     # quick gate smoke run (optional)
        early_exit_gates=True,
    )

    print("Axioms (internal):", out["axioms_internal"])
    print("Axioms (total):   ", out["axioms_total"])
    print("Gate report:", out["gates"])
    print("PASS:", out["pass"])


if __name__ == "__main__":
    main()

"""
/Users/chazzromeo/rAI/Flavor/.venv/bin/python /Users/chazzromeo/rAI/Flavor/alignment_v40_engine_step1.py 
Axioms (internal): {'sa': 0.0, 'odd': 0.0, 'sa_ok': True, 'odd_ok': True}
Axioms (total):    {'sa_total': 0.0, 'odd_total': 0.0, 'sa_ok': True, 'odd_ok': True}
Gate report: GateReport(order0=0.0, first=0.0, self_adj=0.0, odd_int=0.0, order0_pair=(0, 0), first_pair=(0, 0))
PASS: True
"""