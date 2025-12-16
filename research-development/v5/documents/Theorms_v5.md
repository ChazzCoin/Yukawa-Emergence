

## Alignment Spectral Triple v4 Theorem Chain, v5-Compatible (Standalone Manuscript Section)

### 0. Standing hypotheses and notation (v5-compatible)

Fix:

**(H0) Geometric triple.** A real, even spectral triple
[
(A_{\rm geom},H_{\rm geom},D_{\rm geom},J_{\rm geom},\gamma_{\rm geom})
]
with (\gamma_{\rm geom}^2=1), (\gamma_{\rm geom}^\ast=\gamma_{\rm geom}), ([\gamma_{\rm geom},\pi_{\rm geom}(A_{\rm geom})]=0), and ({\gamma_{\rm geom},D_{\rm geom}}=0).

**(H1) SM internal triple.** A real, even finite triple
[
(A_{\rm SM},H_{\rm SM},D_{\rm SM}(\cdot),J_{\rm SM},\gamma_{\rm SM})
]
(one generation) where Yukawa textures are inserted in the standard left–right off-diagonal blocks of the internal Dirac operator.

**(H2) Flavor multiplicity.** A finite-dimensional Hilbert space (H_{\rm flav}\cong\mathbb C^N) with antiunitary real structure
[
J_{\rm flav}=U_{\rm flav}\circ K
]
((K) complex conjugation in some basis, (U_{\rm flav}) unitary).

**Commutators and anticommutators.**
[
[X,Y]:=XY-YX,\qquad {X,Y}:=XY+YX.
]

**Norm.**
[
|X|\ \text{denotes the operator norm for }X\in\mathcal B(H),\qquad |z|\ \text{absolute value for }z\in\mathbb C.
]

**Set braces.** Use ({\cdot}) for sets: (\left{\cdots\right}).

---

## Section 1 — Definitions

### Definition 1.1 (Alignment spectral datum, v5-compatible)

Assume:

* A real, even geometric triple ((A_{\rm geom},H_{\rm geom},D_{\rm geom},J_{\rm geom},\gamma_{\rm geom})).
* A real, even finite SM triple ((A_{\rm SM},H_{\rm SM},D_{\rm SM}(\cdot),J_{\rm SM},\gamma_{\rm SM})).
* A flavor multiplicity space ((H_{\rm flav},J_{\rm flav})).

Define
[
A:=A_{\rm geom}\otimes A_{\rm SM},\qquad
H:=H_{\rm geom}\otimes H_{\rm SM}\otimes H_{\rm flav},
]
[
\pi(a_{\rm geom}\otimes a_{\rm SM})
:=\pi_{\rm geom}(a_{\rm geom})\otimes \pi_{\rm SM}(a_{\rm SM})\otimes 1_{\rm flav}.
]

Let (D_{\rm int}) act on (H_{\rm SM}\otimes H_{\rm flav}) and satisfy the **oddness requirement**
[
\boxed{\ {\gamma_{\rm SM}\otimes 1_{\rm flav},,D_{\rm int}}=0\ }.
]

Define the total Dirac operator in the **even product form**
[
\boxed{\ D:=D_{\rm geom}\otimes 1\otimes 1;+;\gamma_{\rm geom}\otimes D_{\rm int}\ }.
]

Define
[
J:=J_{\rm geom}\otimes J_{\rm SM}\otimes J_{\rm flav},\qquad
\Gamma:=\gamma_{\rm geom}\otimes\gamma_{\rm SM}\otimes 1_{\rm flav}.
]

**Production commutant insertion invariant (v5).** All operators carrying Yukawa/flavor/Majorana dependence inside (D_{\rm int}) lie in
[
\pi(A_{\rm SM})'\otimes \mathcal B(H_{\rm flav})\ \subset\ \pi(A)'\subset\mathcal B(H).
]
(In particular, the sufficient pattern (\widetilde Y=1_{H_{\rm SM}}\otimes Y) is allowed, but literal membership in (\pi(A)') is the invariant.)

**Production constraint (evenness seal, v5).** No additional summand of the form (\gamma_{\rm geom}\otimes 1\otimes D_{\rm flav}) is allowed.

The tuple ((A,H,D,J,\Gamma)) is the **Alignment Spectral Triple (v5-compatible v4 theorem-chain form)**.

---

### Definition 1.2 (Flavor multiplicity principle)

Flavor is a pure multiplicity factor: the represented algebra acts trivially on (H_{\rm flav}):
[
\pi(A)\subseteq\mathcal B(H_{\rm geom}\otimes H_{\rm SM})\otimes 1_{\rm flav}.
]
Equivalently, there is no additional commutative “site algebra” acting diagonally on (H_{\rm flav}).

---

### Definition 1.3 (Reality-compatibility for flavor Yukawas)

A flavor operator (Y\in\mathcal B(H_{\rm flav})) is ((J_{\rm flav}))-real with sign (\varepsilon'*{\rm flav}\in{\pm1}) if
[
\boxed{\ J*{\rm flav}Y=\varepsilon'*{\rm flav},Y,J*{\rm flav}\ }.
]
In a basis where (J_{\rm flav}=K), this is (\overline{Y}=\varepsilon'_{\rm flav}Y).

---

## Section 2 — Lemma: finiteness and completeness

### Lemma 2.1 (Hilbert space structure and bounded representation)

If (H_{\rm SM}) and (H_{\rm flav}) are finite-dimensional, then (H) is a Hilbert space with the tensor-product inner product, and (\pi:A\to\mathcal B(H)) is a bounded (^\ast)-representation provided (\pi_{\rm geom}) and (\pi_{\rm SM}) are bounded.

---

## Section 3 — Theorem: self-adjointness and compact resolvent

Let
[
D_0:=D_{\rm geom}\otimes 1\otimes 1,\qquad
B:=\gamma_{\rm geom}\otimes D_{\rm int}.
]
If (D_{\rm int}) is finite (or bounded relative to (D_0) with relative bound (<1)), then (B) is (D_0)-bounded.

### Theorem 3.1 (Essential self-adjointness)

If (D_{\rm geom}) is essentially self-adjoint on (\mathrm{Dom}(D_{\rm geom})), then (D=D_0+B) is essentially self-adjoint on
[
\mathrm{Dom}(D):=\mathrm{Dom}(D_{\rm geom})\otimes H_{\rm SM}\otimes H_{\rm flav}.
]

### Theorem 3.2 (Compact resolvent)

If (D_{\rm geom}) has compact resolvent, then (D) has compact resolvent.

---

## Section 4 — Theorem: bounded commutators

### Theorem 4.1 (Bounded commutators)

For all simple tensors (a=a_{\rm geom}\otimes a_{\rm SM}\in A),
[
[D,\pi(a)]\in\mathcal B(H),
]
provided ([D_{\rm geom},\pi_{\rm geom}(a_{\rm geom})]) is bounded and ([D_{\rm int},\pi_{\rm SM}(a_{\rm SM})\otimes 1_{\rm flav}]) is bounded.

Explicitly,
[
[D,\pi(a)]
==========

[D_{\rm geom},\pi_{\rm geom}(a_{\rm geom})]\otimes \pi_{\rm SM}(a_{\rm SM})\otimes 1_{\rm flav}
+\gamma_{\rm geom}\pi_{\rm geom}(a_{\rm geom})\otimes [D_{\rm int},\pi_{\rm SM}(a_{\rm SM})\otimes 1_{\rm flav}].
]

---

## Section 5 — Theorem: grading (evenness)

### Theorem 5.1 (Evenness relations)

Assume the geometric and SM factors are even and
[
[\gamma_{\rm geom},\pi_{\rm geom}(A_{\rm geom})]=0,\qquad [\gamma_{\rm SM},\pi_{\rm SM}(A_{\rm SM})]=0.
]
If also ({\gamma_{\rm SM}\otimes 1_{\rm flav},D_{\rm int}}=0), then with
(\Gamma=\gamma_{\rm geom}\otimes\gamma_{\rm SM}\otimes 1_{\rm flav}),
[
\Gamma^\ast=\Gamma,\qquad \Gamma^2=1,\qquad {\Gamma,D}=0,\qquad [\Gamma,\pi(A)]=0.
]

### Corollary 5.2 (Forbidden add-on term)

If one adds a separate summand (\gamma_{\rm geom}\otimes 1\otimes D_{\rm flav}) with generic (D_{\rm flav}), then typically ({\Gamma,D}\neq 0). Hence the production constraint in Definition 1.1 is structurally necessary.

---

## Section 6 — Theorem: reality and KO-signs

### Theorem 6.1 (Reality / opposite action)

(J) is antiunitary on (H) and defines the opposite representation
[
\pi^\circ(a):=J\pi(a)J^{-1}.
]
In a fixed basis with (J=U_J\circ K),
[
JXJ^{-1}=U_J,\overline{X},U_J^\dagger.
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
with ((\varepsilon,\varepsilon',\varepsilon'')) determined by the graded tensor product convention (equivalently by total KO-dimension mod (8)).

---

## Section 7 — Theorem: order-zero condition

### Theorem 7.1 (Order-zero stability under flavor multiplicity)

If order-zero holds on the geometric and SM factors, then for all (a,b\in A),
[
[\pi(a),J\pi(b)J^{-1}]=0.
]
Reason: (\pi(\cdot)) acts as (1_{\rm flav}) on flavor, so no new order-zero obstruction arises from the flavor factor.

---

## Section 8 — Theorem: first-order condition

### Theorem 8.1 (First-order stability)

If first-order holds on the geometric and SM factors (with the same graded-product conventions), then for all (a,b\in A),
[
[[D,\pi(a)],J\pi(b)J^{-1}]=0.
]
The flavor multiplicity does not introduce new obstructions because (\pi(A)) is trivial on (H_{\rm flav}) and all flavor/Yukawa dependence in (D_{\rm int}) is constrained by the commutant insertion invariant (Definition 1.1).

---

## Section 9 — Inner fluctuations and the unified operator

### Definition 9.1 (One-forms and fluctuated Dirac)

[
\Omega^1_D(A):=\left{\sum_i \pi(a_i)[D,\pi(b_i)] : a_i,b_i\in A\right}\subset\mathcal B(H).
]
Given (A\in\Omega^1_D(A)), replace (A\leftarrow \tfrac12(A+A^\dagger)) and define
[
D_A:=D + A + JAJ^{-1}.
]

### Proposition 9.2 (No flavor gauge fields, v5-compatible)

Under Definition 1.2 (multiplicity) and the commutant insertion invariant in Definition 1.1, inner fluctuations cannot generate gauge bosons with nontrivial action on (H_{\rm flav}).

**Reason (one line).** Since (\pi(\cdot)) acts as (1_{\rm flav}), every commutator ([D,\pi(b)]) is in (\mathcal B(H_{\rm geom}\otimes H_{\rm SM})\otimes 1_{\rm flav}), hence (\Omega^1_D(A)\subset(\cdots)\otimes 1_{\rm flav}), and likewise (JAJ^{-1}).

---

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
|P_L^2-P_L|,\ |P_R^2-P_R|,\ |P_LP_R|\ \text{small},\qquad
|P_L+P_R-1|\ \text{small on the chiral subspace}.
]

---

## Appendix: two-line proofs (cleaned)

### Theorem 7.1 (Order-zero) — two-line proof

**Claim.** If order-zero holds on the geometric and SM factors, then for all (a,b\in A),
[
[\pi(a),J\pi(b)J^{-1}]=0 .
]

**Proof.** Check on simple tensors (a=a_{\rm geom}\otimes a_{\rm SM}), (b=b_{\rm geom}\otimes b_{\rm SM}). Using
[
\pi(a)=\pi_{\rm geom}(a_{\rm geom})\otimes\pi_{\rm SM}(a_{\rm SM})\otimes 1_{\rm flav},
]
[
J\pi(b)J^{-1}=\big(J_{\rm geom}\pi_{\rm geom}(b_{\rm geom})J_{\rm geom}^{-1}\big)\otimes\big(J_{\rm SM}\pi_{\rm SM}(b_{\rm SM})J_{\rm SM}^{-1}\big)\otimes 1_{\rm flav},
]
the commutator factorizes into a sum of tensor-factorized terms that vanish by the factor order-zero hypotheses. Extend by bilinearity and norm continuity. (\square)

### Theorem 8.1 (First-order) — two-line proof

**Claim.** If first-order holds on the geometric factor and on the internal factor ((A_{\rm SM},H_{\rm SM}\otimes H_{\rm flav},D_{\rm int},J_{\rm SM}\otimes J_{\rm flav})), then for all (a,b\in A),
[
[[D,\pi(a)],J\pi(b)J^{-1}]=0 .
]

**Proof.** On simple tensors,
[
[D,\pi(a)]
==========

[D_{\rm geom},\pi_{\rm geom}(a_{\rm geom})]\otimes \pi_{\rm SM}(a_{\rm SM})\otimes 1_{\rm flav}
+\gamma_{\rm geom}\pi_{\rm geom}(a_{\rm geom})\otimes [D_{\rm int},\pi_{\rm SM}(a_{\rm SM})\otimes 1_{\rm flav}],
]
and (J\pi(b)J^{-1}) similarly factorizes with (\otimes 1_{\rm flav}). The double commutator splits into the sum of a geometric first-order term and an internal first-order term, each zero by hypothesis. Extend by bilinearity and norm continuity. (\square)

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

Below is a **drop-in update** of the v4.0 theorem-chain so it **fully embodies v5.0r** (evenness seal + commutant-insertion invariant + the v5 “no basis leakage” hygiene). I’m only rewriting the parts that must change; everything else can stay verbatim.

---

## Section 0 — Standing notation (v5 upgrades)

Add/replace these global lines:

**Commutators and anticommutators.**
[
[X,Y]:=XY-YX,\qquad {X,Y}:=XY+YX.
]

**Norm.**
[
|X|\ \text{denotes the operator norm for }X\in\mathcal B(H),\qquad |z|\ \text{absolute value for }z\in\mathbb C.
]

**Set braces.** Use ({\cdot}) for sets: (\left{\cdots\right}), not (\left{\cdots\right}).

---

## Section 1 — Definitions (v5.0r datum)

### Definition 1.1 (Alignment spectral datum v5.0r, replacing v4.0 Def. 1.1)

Assume:

* A real, even spectral triple ((A_{\rm geom},H_{\rm geom},D_{\rm geom},J_{\rm geom},\gamma_{\rm geom})).
* A real, even finite SM triple ((A_{\rm SM},H_{\rm SM},D_{\rm SM}(\cdot),J_{\rm SM},\gamma_{\rm SM})) where Yukawas are inserted in LR off-diagonal blocks of the internal Dirac operator.
* A finite-dimensional flavor space (H_{\rm flav}\cong\mathbb C^N) with antiunitary (J_{\rm flav}=U_{\rm flav}\circ K).
* **Production commutant insertion invariant (v5):** every operator carrying Yukawa/flavor/Majorana dependence inside (D_{\rm int}) lies in
  [
  \pi(A_{\rm SM})'\otimes \mathcal B(H_{\rm flav})\ \subset\ \pi(A)'\subset\mathcal B(H),
  ]
  where (H:=H_{\rm geom}\otimes H_{\rm SM}\otimes H_{\rm flav}) and (\pi) is defined below. (The sufficient pattern (\widetilde Y=1_{H_{\rm SM}}\otimes Y) is allowed, but literal membership in (\pi(A)') is the invariant.)

Define
[
A:=A_{\rm geom}\otimes A_{\rm SM},\qquad
H:=H_{\rm geom}\otimes H_{\rm SM}\otimes H_{\rm flav},
]
[
\pi(a_{\rm geom}\otimes a_{\rm SM})
:=\pi_{\rm geom}(a_{\rm geom})\otimes \pi_{\rm SM}(a_{\rm SM})\otimes 1_{\rm flav}.
]

Let (D_{\rm int}) be an internal operator on (H_{\rm SM}\otimes H_{\rm flav}) satisfying the **oddness requirement**
[
\boxed{\ {\gamma_{\rm SM}\otimes 1_{\rm flav},,D_{\rm int}}=0\ }.
]

Define the total Dirac operator by the **even product form**
[
\boxed{\ D:=D_{\rm geom}\otimes 1\otimes 1;+;\gamma_{\rm geom}\otimes D_{\rm int}\ }.
]

**Production constraint (evenness seal, v5):**
[
\boxed{\ \text{No additional summand of the form }\gamma_{\rm geom}\otimes 1\otimes D_{\rm flav}\text{ is allowed.}\ }
]

Finally,
[
J:=J_{\rm geom}\otimes J_{\rm SM}\otimes J_{\rm flav},\qquad
\Gamma:=\gamma_{\rm geom}\otimes\gamma_{\rm SM}\otimes 1_{\rm flav}.
]
The tuple ((A,H,D,J,\Gamma)) is the **Alignment Spectral Triple v5-compatible (v4 theorem-chain form)**.

---

### Definition 1.2 (Flavor multiplicity principle)  *unchanged*

Keep as is.

---

### Definition 1.3 (Reality-compatibility for flavor Yukawas)  *fix notation*

Replace
[
J_{\rm flav},Y = \varepsilon'*{\rm flav},Y,J*{\rm flav}
]
by the standard relation
[
\boxed{\ J_{\rm flav}Y=\varepsilon'*{\rm flav},Y,J*{\rm flav}\ }.
]
(And remove stray asterisks anywhere: (J_{\rm flav}), (\varepsilon'_{\rm flav}).)

If you want the (K)-basis statement:
[
J_{\rm flav}=K\ \Rightarrow\ \overline{Y}=\varepsilon'_{\rm flav}Y.
]

---

## Section 4 — Bounded commutators (typo fixes only)

In the displayed expansion, remove the “====” separator and keep the identity as an equation. Content is fine.

---

## Section 5 — Evenness (add the v5 failure mode)

### Theorem 5.1 (Evenness relations)  *keep, but now genuinely enforced*

Your theorem is correct once Def. 1.1 includes the **evenness seal** and the oddness condition is written correctly:
[
{\gamma_{\rm SM}\otimes 1_{\rm flav},D_{\rm int}}=0.
]

Add a one-sentence corollary right after Theorem 5.1:

**Corollary 5.2 (Forbidden add-on term).**
If a separate term (\gamma_{\rm geom}\otimes 1\otimes D_{\rm flav}) is added with generic (D_{\rm flav}), then typically ({\Gamma,D}\neq 0). Hence the production constraint in Def. 1.1 is structurally necessary.

---

## Section 9 — Inner fluctuations (v5 commutant protection + set braces)

### Definition 9.1 (One-forms and fluctuated Dirac)  *replace braces*

Replace
[
\Omega^1_D(A):=\left{\sum_i \pi(a_i)[D,\pi(b_i)] : a_i,b_i\in A\right}\subset\mathcal B(H)
]
by
[
\boxed{
\Omega^1_D(A):=\left{\sum_i \pi(a_i)[D,\pi(b_i)] : a_i,b_i\in A\right}\subset\mathcal B(H).
}
]

### Proposition 9.2 (No flavor gauge fields)  *make the v5 logic explicit*

Replace the statement with:

**Proposition 9.2 (No flavor gauge fields, v5-compatible).**
Assume the flavor multiplicity principle (Def. 1.2) and the commutant insertion invariant (Def. 1.1). Then every one-form (A\in\Omega^1_D(A)) acts trivially on (H_{\rm flav}), hence the inner fluctuation (D_A=D+A+JAJ^{-1}) cannot generate gauge bosons with nontrivial flavor action.

**Reason (one line).** Since (\pi(\cdot)) acts as (1_{\rm flav}), we have ([D,\pi(b)]\in \mathcal B(H_{\rm geom}\otimes H_{\rm SM})\otimes 1_{\rm flav}), so (\Omega^1_D(A)\subset(\cdots)\otimes 1_{\rm flav}), and likewise (JAJ^{-1}).

---

## Section 10 — Gates (norm clarity only)

Your gate formulas are fine once Section 0 declares (|\cdot|) as the operator norm. No other change needed.

---

## Proof sketches: remove doubled commas (two occurrences)

* In Theorem 7.1 claim line:
  [
  [\pi(a),J\pi(b)J^{-1}]=0.
  ]
* In Theorem 8.1 claim line:
  [
  [[D,\pi(a)],J\pi(b)J^{-1}]=0.
  ]
* In Theorem 8.1 proof, replace
  ([D_{\rm int},,\pi_{\rm SM}(\cdot)\otimes 1]) by ([D_{\rm int},\pi_{\rm SM}(\cdot)\otimes 1]).
