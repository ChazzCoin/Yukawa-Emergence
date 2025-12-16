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

