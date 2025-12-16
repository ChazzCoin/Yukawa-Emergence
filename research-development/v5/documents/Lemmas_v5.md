## Consolidated Factorization Lemma (v5-compatible)

*(One lemma that yields Sections 4, 7, 8 “by inspection,” now consistent with the v5 evenness seal and commutant-insertion invariant.)*

### Lemma F.1 (Tensor–commutator factorization with multiplicity flavor and v5 Dirac form)

Let
[
A=A_{\rm geom}\otimes A_{\rm SM},\quad
H=H_{\rm geom}\otimes H_{\rm SM}\otimes H_{\rm flav},
]
[
\pi(a_{\rm geom}\otimes a_{\rm SM})
=\pi_{\rm geom}(a_{\rm geom})\otimes \pi_{\rm SM}(a_{\rm SM})\otimes 1_{\rm flav},
]
and let the total Dirac operator be in the **v5 even product form**
[
\boxed{
D
=

D_{\rm geom}\otimes 1\otimes 1
;+;
\gamma_{\rm geom}\otimes D_{\rm int},
}
]
where (D_{\rm int}) acts on (H_{\rm SM}\otimes H_{\rm flav}) and satisfies the **oddness requirement**
[
\boxed{\ {\gamma_{\rm SM}\otimes 1_{\rm flav},,D_{\rm int}}=0\ }.
]

**(v5 evenness seal.)** No separate summand of the form (\gamma_{\rm geom}\otimes 1\otimes D_{\rm flav}) is present.

Let (J=J_{\rm geom}\otimes J_{\rm SM}\otimes J_{\rm flav}) and for (b=b_{\rm geom}\otimes b_{\rm SM}) write
[
b^\circ := J\pi(b)J^{-1}
========================

\big(J_{\rm geom}\pi_{\rm geom}(b_{\rm geom})J_{\rm geom}^{-1}\big)\otimes
\big(J_{\rm SM}\pi_{\rm SM}(b_{\rm SM})J_{\rm SM}^{-1}\big)\otimes 1_{\rm flav}.
]
(Here the final factor is (1_{\rm flav}) because (\pi(\cdot)) acts trivially on flavor.)

Then for all simple tensors (a=a_{\rm geom}\otimes a_{\rm SM}\in A), the following identities hold (and extend to all (a\in A) by linearity and continuity):

---

#### (F.1) Single commutator decomposes into geometric and internal pieces

[
\boxed{
[D,\pi(a)]
==========

[D_{\rm geom},\pi_{\rm geom}(a_{\rm geom})]\otimes \pi_{\rm SM}(a_{\rm SM})\otimes 1_{\rm flav}
;+;
\gamma_{\rm geom}\pi_{\rm geom}(a_{\rm geom})\otimes [D_{\rm int},\pi_{\rm SM}(a_{\rm SM})\otimes 1_{\rm flav}].
}
]
In particular, any **flavor/Yukawa/Majorana dependence** inside (D_{\rm int}) that lies in (\pi(A_{\rm SM})'\otimes \mathcal B(H_{\rm flav})) commutes with (\pi_{\rm SM}(a_{\rm SM})\otimes 1_{\rm flav}) and therefore does not spoil this factorization.

---

#### (F.2) Order-zero commutator factorizes into factor commutators

[
\boxed{
[\pi(a),b^\circ]
================

[\pi_{\rm geom}(a_{\rm geom}),b_{\rm geom}^\circ]\otimes \pi_{\rm SM}(a_{\rm SM}),b_{\rm SM}^\circ\otimes 1_{\rm flav}
;+;
\pi_{\rm geom}(a_{\rm geom}),b_{\rm geom}^\circ\otimes [\pi_{\rm SM}(a_{\rm SM}),b_{\rm SM}^\circ]\otimes 1_{\rm flav}.
}
]

---

#### (F.3) First-order double commutator splits into a sum of “factor checks”

Define
[
X_1:=[D_{\rm geom},\pi_{\rm geom}(a_{\rm geom})],\quad Y_1:=\pi_{\rm SM}(a_{\rm SM}),
]
[
X_2:=\gamma_{\rm geom}\pi_{\rm geom}(a_{\rm geom}),\quad Y_2:=[D_{\rm int},\pi_{\rm SM}(a_{\rm SM})\otimes 1_{\rm flav}],
]
so that ([D,\pi(a)]=(X_1\otimes Y_1 + X_2\otimes Y_2)\otimes 1_{\rm flav}). Then
[
\boxed{
\begin{aligned}
\big[[D,\pi(a)],b^\circ\big]
&=
\Big([X_1,b_{\rm geom}^\circ]\otimes Y_1 b_{\rm SM}^\circ
;+;
b_{\rm geom}^\circ X_1\otimes [Y_1,b_{\rm SM}^\circ]\Big)\otimes 1_{\rm flav}
\
&\quad+
\Big([X_2,b_{\rm geom}^\circ]\otimes Y_2 b_{\rm SM}^\circ
;+;
b_{\rm geom}^\circ X_2\otimes [Y_2,b_{\rm SM}^\circ]\Big)\otimes 1_{\rm flav}.
\end{aligned}
}
]
Moreover, if ([\gamma_{\rm geom},\pi_{\rm geom}(A_{\rm geom})]=0) and (J_{\rm geom}\gamma_{\rm geom}=\varepsilon''\gamma_{\rm geom}J_{\rm geom}) (even real triple), then ([\gamma_{\rm geom},b_{\rm geom}^\circ]=0), hence
[
[X_2,b_{\rm geom}^\circ]
========================

# [\gamma_{\rm geom}\pi_{\rm geom}(a_{\rm geom}),b_{\rm geom}^\circ]

\gamma_{\rm geom}[\pi_{\rm geom}(a_{\rm geom}),b_{\rm geom}^\circ].
]

---

### Proof (mechanical; this is what becomes code)

Use the tensor commutator identity
[
[A\otimes B,\ C\otimes D]=[A,C]\otimes BD + CA\otimes[B,D],
]
and its extension with an extra (\otimes 1_{\rm flav}). Apply it term-by-term to the v5 Dirac form
(D=D_{\rm geom}\otimes 1\otimes 1+\gamma_{\rm geom}\otimes D_{\rm int}),
and use that (b^\circ=b^\circ_{\rm geom}\otimes b^\circ_{\rm SM}\otimes 1_{\rm flav}). (\square)

---

## Immediate corollaries = Sections 4, 7, 8 “for free” (v5-compatible)

### Corollary F.2 (Section 4: bounded commutators)

If ([D_{\rm geom},\pi_{\rm geom}(a_{\rm geom})]) is bounded and ([D_{\rm int},\pi_{\rm SM}(a_{\rm SM})\otimes 1_{\rm flav}]) is bounded for generators, then ([D,\pi(a)]) is bounded for all (a\in A).

### Corollary F.3 (Section 7: order-zero)

If
[
[\pi_{\rm geom}(a_{\rm geom}),b_{\rm geom}^\circ]=0,\qquad
[\pi_{\rm SM}(a_{\rm SM}),b_{\rm SM}^\circ]=0
]
for all (a_{\rm geom},b_{\rm geom},a_{\rm SM},b_{\rm SM}), then ([\pi(a),b^\circ]=0) for all (a,b\in A).

### Corollary F.4 (Section 8: first-order)

Assume factor first-order on the geometric triple and on the **internal triple**
[
(A_{\rm SM},\ H_{\rm SM}\otimes H_{\rm flav},\ D_{\rm int},\ J_{\rm SM}\otimes J_{\rm flav}),
]
namely
[
[[D_{\rm geom},\pi_{\rm geom}(a_{\rm geom})],b_{\rm geom}^\circ]=0,\qquad
\big[[D_{\rm int},\pi_{\rm SM}(a_{\rm SM})\otimes 1_{\rm flav}],,b_{\rm SM}^\circ\otimes 1_{\rm flav}\big]=0,
]
and assume order-zero on the geometric and SM factors as in Corollary F.3. Then ([[D,\pi(a)],b^\circ]=0) for all (a,b\in A).

---

## Lemma F.1′ (Single factorization operator identity, v5 form)

Fix (a=a_{\rm geom}\otimes a_{\rm SM}), (b=b_{\rm geom}\otimes b_{\rm SM}\in A) and set
[
\pi(a)=\pi_{\rm geom}(a_{\rm geom})\otimes\pi_{\rm SM}(a_{\rm SM})\otimes 1_{\rm flav},\qquad
b^\circ=b^\circ_{\rm geom}\otimes b^\circ_{\rm SM}\otimes 1_{\rm flav},
]
with (b^\circ_{\rm geom}:=J_{\rm geom}\pi_{\rm geom}(b_{\rm geom})J_{\rm geom}^{-1}) and (b^\circ_{\rm SM}:=J_{\rm SM}\pi_{\rm SM}(b_{\rm SM})J_{\rm SM}^{-1}). Define
[
C_{\rm geom}(a_{\rm geom}):=[D_{\rm geom},\pi_{\rm geom}(a_{\rm geom})],
\qquad
C_{\rm int}(a_{\rm SM}):=[D_{\rm int},\pi_{\rm SM}(a_{\rm SM})\otimes 1_{\rm flav}].
]
Then:
[
\boxed{
\begin{aligned}
[D,\pi(a)]
&=\Big(C_{\rm geom}(a_{\rm geom})\otimes \pi_{\rm SM}(a_{\rm SM})
;+;
\gamma_{\rm geom}\pi_{\rm geom}(a_{\rm geom})\otimes C_{\rm int}(a_{\rm SM})\Big)\otimes 1_{\rm flav},
[3pt]
[\pi(a),b^\circ]
&=\Big([\pi_{\rm geom}(a_{\rm geom}),b^\circ_{\rm geom}]\otimes \pi_{\rm SM}(a_{\rm SM}),b^\circ_{\rm SM}
;+;
\pi_{\rm geom}(a_{\rm geom}),b^\circ_{\rm geom}\otimes [\pi_{\rm SM}(a_{\rm SM}),b^\circ_{\rm SM}]\Big)\otimes 1_{\rm flav},
[3pt]
\big[[D,\pi(a)],b^\circ\big]
&=\Big(
[C_{\rm geom}(a_{\rm geom}),b^\circ_{\rm geom}]\otimes \pi_{\rm SM}(a_{\rm SM}),b^\circ_{\rm SM}
;+;
b^\circ_{\rm geom}C_{\rm geom}(a_{\rm geom})\otimes [\pi_{\rm SM}(a_{\rm SM}),b^\circ_{\rm SM}]
\
&\qquad\quad
+;[\gamma_{\rm geom}\pi_{\rm geom}(a_{\rm geom}),b^\circ_{\rm geom}]\otimes C_{\rm int}(a_{\rm SM}),b^\circ_{\rm SM}
;+;
b^\circ_{\rm geom}\gamma_{\rm geom}\pi_{\rm geom}(a_{\rm geom})\otimes [C_{\rm int}(a_{\rm SM}),b^\circ_{\rm SM}\otimes 1_{\rm flav}]
\Big)\otimes 1_{\rm flav}.
\end{aligned}}
]
*Proof.* Repeated use of ([X\otimes Y, X'\otimes Y']=[X,X']\otimes YY' + X'X\otimes [Y,Y']), together with (\pi(\cdot)\otimes 1_{\rm flav}). (\square)

