3.5 Axioms for Alignment Operators on the Proto–Flavor Lattice
In this section we isolate the structural assumptions that define the class of Schur–type operators allowed to act on the flavor blocks of the finite Dirac operator. These axioms are additional finite–geometry postulates: they are not implied by the standard noncommutative geometry axioms, but encode our alignment principle for fermion masses and mixings.
They will be the only input used to single out the Alignment Operator Φ\PhiΦ used in the phenomenological analysis. In particular, in later subsections we show that, under these axioms, Φ\PhiΦ is essentially unique up to a single continuous parameter.
Proto–flavor space and Schur action
We work with a 9–dimensional proto–flavor space
F9  =  Span{∣i⟩:i=1,…,9},\mathcal{F}_9 \;=\; \mathrm{Span}\{|i\rangle : i=1,\dots,9\},F9​=Span{∣i⟩:i=1,…,9},
equipped with the standard orthonormal basis ∣i⟩|i\rangle∣i⟩. This space represents a finite internal flavor lattice on which the three light generations and their heavy partners are embedded.
Proto–Yukawa and proto–Majorana couplings are encoded as 9×99\times 99×9 matrices
Y, M∈M9(C).Y,\,M \in M_9(\mathbb{C}).Y,M∈M9​(C).
An alignment operator is a linear map
Φ:M9(C)⟶M9(C)\Phi : M_9(\mathbb{C}) \longrightarrow M_9(\mathbb{C})Φ:M9​(C)⟶M9​(C)
acting by Schur (entrywise) multiplication with a fixed kernel K∈M9(C)K\in M_9(\mathbb{C})K∈M9​(C):
 \begin{equation}
 \label{eq:Schur-action}
 \bigl(\Phi(Y)\bigr){ij} ;:=; (K\circ Y){ij} ;:=; K_{ij},Y_{ij},
 \qquad i,j=1,\dots,9.
 \end{equation}
 Our goal is to constrain the kernel KKK by a small set of axioms and then classify the resulting operators.
Axioms A1–A6
We now formulate the structural axioms for admissible alignment kernels. These axioms package the geometric, symmetry, divisor, and energy–minimization requirements that underlie the Alignment Operator used in this work.
Axiom A1 (Finite cyclic geometry)
The proto–flavor space F9\mathcal{F}_9F9​ carries the structure of a cyclic graph with 9 sites. Concretely, we fix a graph distance
d:{1,…,9}2⟶{0,1,2,3,4}d : \{1,\dots,9\}^2 \longrightarrow \{0,1,2,3,4\}d:{1,…,9}2⟶{0,1,2,3,4}
given by
d(i,j)  :=  min⁡(∣i−j∣, 9−∣i−j∣),i,j=1,…,9.d(i,j) \;:=\; \min\bigl(|i-j|,\,9-|i-j|\bigr), \qquad i,j=1,\dots,9.d(i,j):=min(∣i−j∣,9−∣i−j∣),i,j=1,…,9.
The alignment kernel KijK_{ij}Kij​ depends only on this graph distance, so that there exists a function
k:{0,1,2,3,4}⟶Ck : \{0,1,2,3,4\} \longrightarrow \mathbb{C}k:{0,1,2,3,4}⟶C
with
 \begin{equation}
 \label{eq:Toeplitz}
 K_{ij} ;=; k\bigl(d(i,j)\bigr),
 \qquad i,j=1,\dots,9.
 \end{equation}
 Thus Φ\PhiΦ is a Toeplitz (translation–invariant) Schur operator on the 9–site cycle.
Axiom A2 (Hermiticity and complete positivity)
The operator Φ\PhiΦ preserves Hermiticity and does not introduce anti–alignment.
For any Hermitian matrix Y=Y†Y=Y^\daggerY=Y†, the aligned matrix Φ(Y)\Phi(Y)Φ(Y) is also Hermitian.


The kernel KKK is real, symmetric, and positive semidefinite, with unit diagonal:


Kij=Kji∈R,K⪰0,Kii=1for all i=1,…,9.K_{ij} = K_{ji} \in \mathbb{R},\qquad K \succeq 0,\qquad K_{ii}=1\quad\text{for all }i=1,\dots,9.Kij​=Kji​∈R,K⪰0,Kii​=1for all i=1,…,9.
Equivalently, Φ\PhiΦ is a completely positive, unital Schur multiplier on M9(C)M_9(\mathbb{C})M9​(C).
Axiom A3 (Generation permutation covariance)
The three light generations are a priori indistinguishable. Let G⊂S9G\subset S_9G⊂S9​ be the subgroup of permutations implementing the flavor symmetries of the three light generations and their associated heavy partners on F9\mathcal{F}_9F9​, and let PPP denote the corresponding permutation matrices acting on M9(C)M_9(\mathbb{C})M9​(C) by conjugation.
The action of Φ\PhiΦ commutes with these flavor permutations:
 \begin{equation}
 \label{eq:permutation-covariance}
 \Phi(PYP^{-1}) ;=; P,\Phi(Y),P^{-1},
 \qquad\text{for all }Y\in M_9(\mathbb{C})\text{ and }P\in G.
 \end{equation}
 In particular, the alignment kernel does not distinguish between individual generations beyond their relative positions on the proto–flavor graph.
Axiom A4 (Ambient divisor–closed harmonic support with a single forbidden distance)
There exist an integer N∈NN \in \mathbb{N}N∈N, an injective map
ι:{1,…,9}↪ZN,\iota : \{1,\dots,9\} \hookrightarrow \mathbb{Z}_N,ι:{1,…,9}↪ZN​,
and a function k~:{0,1,…,N−1}→C\tilde k : \{0,1,\dots,N-1\} \to \mathbb{C}k~:{0,1,…,N−1}→C such that:
Ambient cyclic geometry and induced distance.
 Let dN:ZN×ZN→{0,1,…,⌊N/2⌋}d_N : \mathbb{Z}_N \times \mathbb{Z}_N \to \{0,1,\dots,\lfloor N/2 \rfloor\}dN​:ZN​×ZN​→{0,1,…,⌊N/2⌋} be the cyclic graph distance
 dN(x,y)  :=  min⁡(∣x−y∣, N−∣x−y∣),x,y∈ZN.d_N(x,y) \;:=\; \min\bigl(|x-y|,\,N - |x-y|\bigr),\qquad x,y\in\mathbb{Z}_N.dN​(x,y):=min(∣x−y∣,N−∣x−y∣),x,y∈ZN​.
 The proto–flavor graph F9\mathcal F_9F9​ is realized as the induced subgraph of the cyclic graph on ZN\mathbb{Z}_NZN​ with vertex set ι({1,…,9})\iota(\{1,\dots,9\})ι({1,…,9}), and this induced subgraph is isomorphic to the 9–cycle. The distance ddd used in Axiom A1 is the restriction
 d(i,j)  :=  dN(ι(i),ι(j)),i,j∈{1,…,9}.d(i,j) \;:=\; d_N\bigl(\iota(i),\iota(j)\bigr),\qquad i,j\in\{1,\dots,9\}.d(i,j):=dN​(ι(i),ι(j)),i,j∈{1,…,9}.
Ambient Toeplitz kernel and restriction to F9\mathcal F_9F9​.
 Define an N×NN\times NN×N Toeplitz kernel K~\tilde KK~ on CN\mathbb{C}^NCN by
 K~xy  :=  k~(dN(x,y)),x,y∈ZN.\tilde K_{xy} \;:=\; \tilde k\bigl(d_N(x,y)\bigr), \qquad x,y\in\mathbb{Z}_N.K~xy​:=k~(dN​(x,y)),x,y∈ZN​.
 The 9×99\times 99×9 alignment kernel KKK on F9\mathcal F_9F9​ is the restriction of K~\tilde KK~ to the embedded sites:
 Kij  :=  K~ι(i),ι(j)  =  k~(dN(ι(i),ι(j)))  =  k~(d(i,j)),i,j=1,…,9.K_{ij} \;:=\; \tilde K_{\iota(i),\iota(j)} \;=\; \tilde k\bigl(d_N(\iota(i),\iota(j))\bigr) \;=\; \tilde k\bigl(d(i,j)\bigr), \qquad i,j=1,\dots,9.Kij​:=K~ι(i),ι(j)​=k~(dN​(ι(i),ι(j)))=k~(d(i,j)),i,j=1,…,9.
Divisor–closed ambient harmonic support.
 Let k~^\widehat{\tilde k}k~ denote the discrete Fourier transform of k~\tilde kk~ on ZN\mathbb{Z}_NZN​, so that
 k~(d)  =  1N∑n=0N−1k~^(n) e2πind/N,d=0,1,…,N−1.\tilde k(d) \;=\; \frac{1}{N} \sum_{n=0}^{N-1} \widehat{\tilde k}(n)\,e^{2\pi i n d/N}, \qquad d=0,1,\dots,N-1.k~(d)=N1​n=0∑N−1​k~(n)e2πind/N,d=0,1,…,N−1.
 Define the spectral support
 D  :=  { n∈{0,…,N−1}:k~^(n)≠0 }.D \;:=\; \{\, n \in \{0,\dots,N-1\} : \widehat{\tilde k}(n) \neq 0 \,\}.D:={n∈{0,…,N−1}:k~(n)=0}.
 Then DDD is divisor–closed, in the sense that
 n∈D,  m∣n  ⟹  m∈D.n \in D,\; m\mid n \;\Longrightarrow\; m \in D.n∈D,m∣n⟹m∈D.
Single forbidden ambient distance.
 There exists a unique integer
 d⋆∈{1,…,N−1}d_\star \in \{1,\dots,N-1\}d⋆​∈{1,…,N−1}
 such that
 k~(d⋆)=0,k~(d)≠0    for all    d∈{1,…,N−1}∖{d⋆}.\tilde k(d_\star) = 0, \qquad \tilde k(d) \neq 0 \;\;\text{for all}\;\; d \in \{1,\dots,N-1\}\setminus\{d_\star\}.k~(d⋆​)=0,k~(d)=0for alld∈{1,…,N−1}∖{d⋆​}.
 Thus k~\tilde kk~ has a single forbidden separation scale in the ambient cyclic geometry. The 9–site profile kkk appearing in Axiom A1 is the restriction
 k(d)  :=  k~(d),d∈{0,1,2,3,4},k(d) \;:=\; \tilde k(d),\qquad d \in \{0,1,2,3,4\},k(d):=k~(d),d∈{0,1,2,3,4},
 obtained by evaluating k~\tilde kk~ only on the distances actually realized among the embedded sites. No vanishing of k(d)k(d)k(d) on {1,2,3,4}\{1,2,3,4\}{1,2,3,4} is assumed a priori.


Minimality of the ambient geometry.
 The integer NNN is minimal with this property: there is no smaller N′<NN'<NN′<N and injective map ι′:{1,…,9}↪ZN′\iota' : \{1,\dots,9\}\hookrightarrow \mathbb{Z}_{N'}ι′:{1,…,9}↪ZN′​ together with a function k~′:{0,…,N′−1}→C\tilde k' : \{0,\dots,N'-1\}\to\mathbb{C}k~′:{0,…,N′−1}→C satisfying items (1)–(4) above.
Axiom A5 (NCG compatibility)
The action of Φ\PhiΦ on flavor blocks of the finite Dirac operator DFD_FDF​ preserves the axioms of the finite spectral triple (AF,HF,DF,JF,ΓF)(A_F,H_F,D_F,J_F,\Gamma_F)(AF​,HF​,DF​,JF​,ΓF​).
Concretely:
Φ\PhiΦ acts only on flavor multiplicity indices and commutes with the representation of AFA_FAF​ on HFH_FHF​; in particular, if πF:AF→B(HF)\pi_F : A_F \to \mathcal{B}(H_F)πF​:AF​→B(HF​) denotes the representation, then
 [ πF(a), Φ(DF) ]  =  [ πF(a), DF ]for all a∈AF.[\,\pi_F(a),\,\Phi(D_F)\,] \;=\; [\,\pi_F(a),\,D_F\,] \qquad\text{for all }a\in A_F.[πF​(a),Φ(DF​)]=[πF​(a),DF​]for all a∈AF​.
The reality structure JFJ_FJF​, grading ΓF\Gamma_FΓF​, and KO–dimension constraints are unchanged:
 [JF,Φ(DF)]=[JF,DF],[ΓF,Φ(DF)]=[ΓF,DF].[J_F,\Phi(D_F)] = [J_F,D_F],\qquad [\Gamma_F,\Phi(D_F)] = [\Gamma_F,D_F].[JF​,Φ(DF​)]=[JF​,DF​],[ΓF​,Φ(DF​)]=[ΓF​,DF​].
The first–order condition remains satisfied when the Yukawa and Majorana blocks of DFD_FDF​ are replaced by their aligned versions, i.e.
 [[DF,πF(a)], πF(b)∘]=0⟹[[Φ(DF),πF(a)], πF(b)∘]=0,[[D_F,\pi_F(a)],\,\pi_F(b)^{\circ}] = 0 \quad\Longrightarrow\quad [[\Phi(D_F),\pi_F(a)],\,\pi_F(b)^{\circ}] = 0,[[DF​,πF​(a)],πF​(b)∘]=0⟹[[Φ(DF​),πF​(a)],πF​(b)∘]=0,
 for all a,b∈AFa,b\in A_Fa,b∈AF​, where πF(b)∘:=JF πF(b)∗JF−1\pi_F(b)^{\circ} := J_F\,\pi_F(b)^\ast J_F^{-1}πF​(b)∘:=JF​πF​(b)∗JF−1​.


In other words, Φ\PhiΦ refines the internal geometry without changing the algebraic content of the NCG Standard Model.
Axiom A6 (Alignment energy minimization)
Among all kernels satisfying A1–A5, the physical alignment kernel minimizes a discrete “misalignment energy” functional built from a graph Laplacian on the 9–site cycle.
More precisely, there exists a positive quadratic functional
 \begin{equation}
 \label{eq:alignment-energy}
 \mathcal{E}[k] ;=; \sum_{d=0}^{4} \Bigl[(\nabla k)(d)^2 + m^2,k(d)^2\Bigr],
 \end{equation}
 where ∇\nabla∇ is a fixed discrete derivative on the distance set {0,1,2,3,4}\{0,1,2,3,4\}{0,1,2,3,4} induced by the 9–cycle and m>0m>0m>0 sets an internal curvature scale, such that the realized profile kkk is the unique minimizer of E\mathcal{E}E subject to the constraints of A2–A4 (in particular k(0)=1k(0)=1k(0)=1, symmetry on the 9–cycle, and the single forbidden ambient distance).
Operationally, A6 expresses the requirement that alignment decays with distance in the energetically least costly way allowed by the discrete geometry. In Section 3.z we will show that this principle forces an exponential profile
 \begin{equation}
 \label{eq:exponential-profile}
 k(d) ;\propto; \kappa^{,d},\qquad d=0,1,2,3,4,
 \end{equation}
 for some κ∈(0,1)\kappa\in(0,1)κ∈(0,1).

Definition of the admissible class C0\mathcal{C}_0C0​
We can now summarize:
\boxed{ \begin{minipage}{0.9\textwidth} \textbf{Definition 3.x (Admissible alignment operators).} Let \(\mathcal{C}_0\) denote the set of linear maps \(\Phi : M_9(\mathbb{C})\to M_9(\mathbb{C})\) of the Schur form \eqref{eq:Schur-action} whose kernel \(K\) satisfies axioms A1–A6. An \emph{alignment operator} in the sense of this paper is any \(\Phi\in\mathcal{C}_0\). \end{minipage} }
In the remainder of this section we classify the elements of C0\mathcal{C}_0C0​. We will show that, under A1–A6, C0\mathcal{C}_0C0​ contains a unique element up to a single continuous parameter κ\kappaκ, and that this element coincides with the Alignment Operator Φ\PhiΦ used in the flavor sector of the NCG Standard Model.
(Immediately after A4/A6 you can insert your Proposition stating N=360,d⋆=7N=360, d_\star=7N=360,d⋆​=7 under A1–A4 + A6.)



	•	A. Kinematics & internal space
	•	B. Core operators (R, Q, L)
	•	C. Misalignment & Evolution operator \hat{\mathbb{M}}
	•	D. Selection operator \hat{\mathbb{S}}
	•	E. Manifestation & physical states \hat{\mathbb{X}}
	•	F. Flavor / Yukawa structure
	•	G. Minimality / no cheating

A. Kinematics & internal space

Axiom A1 (Total Hilbert space).
There exists a separable Hilbert space
\mathcal{H} = \mathcal{H}_{\text{ext}} \otimes \mathcal{H}_{\text{int}},
where:
	•	\mathcal{H}_{\text{ext}} carries the usual spacetime degrees of freedom (e.g. spinors over a 4D Lorentzian manifold, or a standard QFT Fock space).
	•	\mathcal{H}_{\text{int}} is a finite- or countable-dimensional internal Hilbert space encoding flavor, charges, and the “aether crystal” structure.

Axiom A2 (Internal quasi-crystal structure).
On \mathcal{H}_{\text{int}} there exists a self-adjoint “internal Laplacian” \hat{L} such that:
	1.	\hat{L} is the graph Laplacian of a bounded-degree, locally finite graph G = (V,E) (internal crystal graph).
	2.	G has finite local complexity (each finite patch appears with bounded frequency), but is not globally periodic (quasi-crystal / aperiodic allowed).
	3.	The spectrum \sigma(\hat{L}) is purely discrete or has at least a dense discrete component sufficient to define a spectral kernel F(\hat{L}).

This axiom is your “aether is a (quasi)crystalline medium” in precise form.

⸻

B. Core internal operators: R, Q, L

Axiom B1 (Base-360 internal symmetry).
There exists a unitary operator \hat{R} on a generation subspace
\mathcal{H}_{\text{gen}} \subset \mathcal{H}_{\text{int}}, \quad \dim \mathcal{H}_{\text{gen}} = 3,
such that:
	1.	\hat{R}^{360} = \mathbb{I}_{\text{gen}}.
	2.	The spectrum of \hat{R} on \mathcal{H}_{\text{gen}} consists of three distinct characters
\chi_j = e^{2\pi i k_j / 360}, \quad j=1,2,3,
with integers k_j chosen so that:
	•	the associated periods T_j = 360 / \gcd(360,k_j) form a triad of distinct divisors of 360 (e.g. 60, 120, 180),
	•	no nontrivial linear relation a_1 k_1 + a_2 k_2 + a_3 k_3 \equiv 0\ \text{mod}\ 360 holds with small integers a_i beyond the obvious ones.

This encodes your base-360 triadic generation structure as a spectral property of \hat{R}.

Axiom B2 (Charge operator Q).
There exists a self-adjoint operator \hat{Q} on \mathcal{H}_{\text{int}} such that:
	1.	\hat{Q} has pure point spectrum contained in \mathbb{Z} (integer charges).
	2.	The restriction of \hat{Q} to \mathcal{H}_{\text{gen}} has a finite set of eigenvalues \{q_{s,g}\}, where:
	•	s labels sectors (e.g. up, down, charged lepton, neutrino),
	•	g = 1,2,3 labels generation.
	3.	[\hat{Q}, \hat{R}] = 0 on \mathcal{H}_{\text{gen}} (charges and 360-phases can be simultaneously diagonalized on generation space).

This is the rigorous version of your discrete exponent / “FN-like” charge structure.

Axiom B3 (Compatibility of L, R, Q).
On \mathcal{H}_{\text{int}}:
	1.	\hat{L} commutes with \hat{R} on \mathcal{H}_{\text{gen}} up to a phase or a bounded operator that preserves the triadic subspace.
	2.	The joint spectrum of (\hat{L}, \hat{R}, \hat{Q}) is nondegenerate enough to distinguish at least three generation modes and four sector patterns, i.e. no accidental degeneracy identifies distinct generations or sectors.

This is a “no accidental degeneracy” axiom: your operators are rich enough to label what you call flavors.

⸻

C. Misalignment functional & Evolution operator \hat{\mathbb{M}}

Axiom C1 (Misalignment functional).
There exists a functional
M[\Psi]: \mathcal{D} \subset \mathcal{H} \to \mathbb{R}_{\ge 0}
defined on a dense domain \mathcal{D} of \mathcal{H}, such that:
	1.	Positivity: M[\Psi] \ge 0, and M[\Psi] = 0 iff \Psi is perfectly aligned (ground state configuration) with respect to \hat{L}, \hat{R}, \hat{Q} and external constraints (e.g. gauge, gravity).
	2.	Locality on the internal graph:
M decomposes as
M[\Psi] = M_{\text{strain}}[\Psi] + M_{\text{bend}}[\Psi] + M_{\text{phase}}[\Psi] + M_{\text{defect}}[\Psi],
where each term is a sum/integral over local contributions on the internal quasi-crystal graph G (e.g. involving nearest-neighbor differences, discrete curvature, phase differences, defect densities).
	3.	Gauge & symmetry invariance: M[\Psi] is invariant under the action of the external gauge group and under global phase rotations of Ψ; it is also invariant (or covariant in a fixed way) under the internal symmetry generated by \hat{R} and the automorphism group of G.

This is your “viscoelastic crystal energy” written as a rigorous functional.

Axiom C2 (Evolution operator).
There exists a one-parameter family of linear maps \hat{\mathbb{M}}(t) on \mathcal{H} such that:
	1.	\hat{\mathbb{M}}(0) = \mathbb{I}.
	2.	\hat{\mathbb{M}}(t+s) = \hat{\mathbb{M}}(t)\hat{\mathbb{M}}(s) for t,s \ge 0 (semigroup property).
	3.	For any \Psi \in \mathcal{D},
\left.\frac{d}{dt}\right|_{t=0} \Psi(t) = - f(0)\,\frac{\delta M}{\delta\Psi^\dagger}(\Psi),
where f(t) is a nonnegative scalar function and \Psi(t) = \hat{\mathbb{M}}(t)\Psi.
	4.	\hat{\mathbb{M}}(t) is completely positive and contractive (or unitary in an appropriate limit), ensuring physical evolution.

This encodes your idea
\hat{\mathbb{M}}(t) \sim \exp\left[-f(t)\frac{\delta M}{\delta\Psi^\dagger}\right]
as a proper dissipative/flow evolution operator.

⸻

D. Selection operator \hat{\mathbb{S}}

You wrote:
\hat{\mathbb{S}} = \hat{C}_{360}\,\hat{B}\,\hat{P}_\phi.

We now define each piece.

Axiom D1 (\hat{C}_{360}: harmonic / divisor filter).
There exists a self-adjoint operator \hat{C}_{360} on \mathcal{H}_{\text{int}} such that:
	1.	\hat{C}_{360} is a spectral projector built from \hat{R} and/or \hat{L}:
\hat{C}_{360} = \sum_{(k,\lambda) \in \mathcal{K}} |k,\lambda\rangle\langle k,\lambda|,
where |k,\lambda\rangle are joint eigenstates of \hat{R}, \hat{L} and \mathcal{K} is a subset of modes whose 360-periods and Laplacian eigenvalues satisfy a divisor/triad condition (e.g. only modes corresponding to periods {60,120,180} or similar).
	2.	\hat{C}_{360} is an orthogonal projector: \hat{C}_{360}^2 = \hat{C}_{360} = \hat{C}_{360}^\dagger.

This is the precise “only certain base-360 harmonics are allowed” rule.

Axiom D2 (\hat{B}: geometric / sublattice selector).
There exists a self-adjoint operator \hat{B} on \mathcal{H}_{\text{int}} such that:
	1.	\hat{B} is a projector onto a “light” internal subspace:
\hat{B}^2 = \hat{B} = \hat{B}^\dagger.
	2.	In a basis adapted to the internal graph G, \hat{B} acts as multiplication by 0 or 1 at each site (or block of sites), selecting sublattices / defect cores / regions that support low-energy, observable excitations.
	3.	The complement (\mathbb{I}-\hat{B}) corresponds to “heavy” mediator states that can be integrated out (e.g. via Schur complement or effective field theory).

This is the mathematically clean version of your “light vs heavy sites / sublattices” idea.

Axiom D3 (\hat{P}_\phi: phase-coherence projector).
There exists a finite group \mathcal{G}_\phi of unitary operators acting on \mathcal{H}_{\text{int}} (phase symmetries) and a projector:

\hat{P}_\phi = \frac{1}{|\mathcal{G}_\phi|} \sum_{g \in \mathcal{G}_\phi} U_g.

Such that:
	1.	\hat{P}_\phi^2 = \hat{P}_\phi = \hat{P}_\phi^\dagger.
	2.	\hat{P}_\phi projects onto the subspace of phase-locked, generation-structured states, implementing discrete relations among phases (e.g. triadic phase patterns, Cabibbo-like rotations, discrete flavor group A_4/S_4 structure).

This is the rigorous form of your “phase wheels” / “phase lock” selection.

Axiom D4 (Selection operator).
Define the Selection operator by:

\hat{\mathbb{S}} = \hat{C}_{360}\,\hat{B}\,\hat{P}_\phi.

Assume:
	1.	\hat{C}_{360}, \hat{B}, \hat{P}_\phi commute up to bounded corrections that do not change their joint eigenspaces on the physical subspace.
	2.	\hat{\mathbb{S}} is an orthogonal projector on \mathcal{H}.

⸻

E. Manifestation & physical states \hat{\mathbb{X}}

Axiom E1 (Manifestation operator).
Define:

\hat{\mathbb{X}}(t) = \hat{\mathbb{S}} \circ \hat{\mathbb{M}}(t)
on \mathcal{H}.

We will usually suppress t and consider the late-time / fixed-point limit:
\hat{\mathbb{X}} := \lim_{t\to\infty} \hat{\mathbb{S}} \circ \hat{\mathbb{M}}(t),
assuming the limit exists on the relevant domain.

Axiom E2 (Universal harmonic state equation).
A state \Psi \in \mathcal{H} is physically manifested if and only if it satisfies:

\hat{\mathbb{X}} \Psi = \Psi,
or in expanded form,
\left( \hat{C}_{360}\,\hat{B}\,\hat{P}_\phi \, e^{-f(\infty)\frac{\delta M}{\delta \Psi^\dagger}} - \mathbb{I} \right)\Psi = 0.

Interpretation:
	•	Only states that are fixed points of “evolve + select” are actually realized in the universe.
	•	Particle species correspond to minimal nontrivial invariant subspaces of \hat{\mathbb{X}}.
	•	Masses correspond to eigenvalues of the appropriate Dirac/Yukawa operators on these subspaces.

This is your central “Selection × Evolution = Reality” condition.

⸻

F. Flavor & Yukawa structure

Now we pin down what “flavor” means in this axiom system.

Axiom F1 (Yukawa operators from R and Q).
For each sector s (up, down, charged lepton, neutrino), there exists a pair of unitaries (U_L^{(s)}, U_R^{(s)}) on \mathcal{H}_{\text{gen}} such that the internal Yukawa operator for that sector is:

\hat{Y}_s = U_L^{(s)\dagger} \, F(\hat{R}) \, e^{-\beta \hat{Q}_s} \, U_R^{(s)},

where:
	1.	F(\hat{R}) is a universal spectral function of \hat{R} (or of \hat{L}, \hat{R} jointly), e.g.
F(\hat{R}) = \exp\big( - \lambda ( \mathbb{I} - \Re \hat{R} )\big) \quad\text{or}\quad F(\hat{L}) = e^{-\alpha \hat{L}},
with fixed, universal \lambda,\alpha > 0.
	2.	\hat{Q}_s is \hat{Q} restricted to the sector s acting on \mathcal{H}_{\text{gen}}, with integer eigenvalues q_{s,g}.
	3.	\beta is a fixed universal constant (no sector-dependent continuous tuning).
	4.	The only sector dependence enters via:
	•	the discrete spectrum of \hat{Q}_s (the charges q_{s,g}),
	•	and the discrete choice of unitaries (U_L^{(s)}, U_R^{(s)}), constrained by the internal symmetry group.

This is the “Yukawas = function of (R,Q) in a sector-dependent basis” axiom.

Axiom F2 (No continuous flavor knobs).
The following are not allowed as primitive input:
	•	Sector-dependent continuous parameters (no \alpha_u, \alpha_d, \dots).
	•	Hand-drawn exponent tables outside the integer charge spectrum of \hat{Q}.
	•	Arbitrary random matrices in the flavor sector.

All flavor differences must arise from:
	1.	Integer charge spectra of \hat{Q}.
	2.	Representation-theoretic data (discrete groups, irreps, and their Clebsch–Gordan coefficients) used to pick (U_L^{(s)}, U_R^{(s)}).
	3.	The spectral data of the internal operators \hat{L}, \hat{R} on the quasi-crystal.

This is the “no cheating” axiom: flavor must be emergent from the operator algebra.

Axiom F3 (Emergent mass hierarchies and mixing).

When restricted to \mathcal{H}_{\text{gen}}:
	1.	The singular values of \hat{Y}_s (up to an overall scale) are hierarchical, i.e.
m_{s,1} \ll m_{s,2} \ll m_{s,3},
and these hierarchies must be expressible as:
m_{s,g} \propto \lambda^{q_{s,g}} \times \text{(O(1) factors from }F(\hat{R})\text{)},
for some small universal \lambda = e^{-\beta} \approx 0.2\text{–}0.3.
	2.	The left-handed diagonalization matrices U_L^{(s)} for quark and lepton sectors yield:
	•	CKM and PMNS mixing matrices
V_{\text{CKM}} = U_L^{(u)\dagger} U_L^{(d)},\quad
U_{\text{PMNS}} = U_L^{(e)\dagger} U_L^{(\nu)},
with entries determined entirely by the internal symmetry group, \hat{R}, \hat{Q}, and \hat{P}_\phi, not by extra continuous fits.

This is the axiom that ties your operator algebra directly to observable flavor.

⸻

G. Minimality / “top-tier physics” sanity conditions

Axiom G1 (Spectral minimality).
The operator set \{\hat{L}, \hat{R}, \hat{Q}, \hat{C}_{360}, \hat{B}, \hat{P}_\phi\} is minimal in the sense that removing any one of them (or its role) either:
	•	collapses the hierarchy structure (e.g. generations become degenerate), or
	•	leads to vanishing Yukawas or trivial mixing.

This is the “none of these are decorative; each plays an essential role” axiom.

Axiom G2 (Consistency with external physics).
When tensored with a suitable external spectral triple (e.g. ordinary 4D geometry + SM gauge algebra), the combined structure must:
	1.	Obey standard locality, causality, and unitarity (or a well-defined generalization).
	2.	Recover known gauge charges and Lorentz representations for observed particles.
	3.	Admit a renormalizable or at least controllably effective field theory description at accessible energies.

In other words: the internal aether crystal is allowed to be exotic, but it must sit consistently under the usual QFT / GR umbrella.

⸻

Master equation (your “universal harmonic equation”, cleaned up)

All of this can be summarized in the single “master condition”:

Physical states are fixed points of Selection ∘ Evolution:
> \hat{\mathbb{X}} \Psi = \Psi,
>
with
> \hat{\mathbb{X}}
> = \hat{\mathbb{S}} \circ \hat{\mathbb{M}}
> = \hat{C}_{360}\,\hat{B}\,\hat{P}_\phi \,
> \exp\!\left[-f(\infty) \frac{\delta M}{\delta \Psi^\dagger} \right].
>

Expanded as an operator equation:
> \left(
> \hat{C}_{360}\,\hat{B}\,\hat{P}_\phi\,
> e^{-f(\infty)\frac{\delta M}{\delta \Psi^\dagger}}
> - \mathbb{I}
> \right)\Psi = 0.
>

Everything else in the axioms is just:
	•	specifying what \hat{C}_{360}, \hat{B}, \hat{P}_\phi, \hat{L}, \hat{R}, \hat{Q}, M are,
	•	and demanding that they come only from discrete spectra + quasi-crystal geometry + internal symmetry, not from arbitrary knobs.

import numpy as np
import math

"""
Emergent aether toy:
====================

1) Start with N sites, each with a phase theta_i ∈ [0, 2π).
2) Define a misalignment functional

       M[theta] = sum_{i<j} J_ij [ (1 - cos(6Δ_ij)) + (1 - cos(5Δ_ij)) ]

   where Δ_ij = theta_i - theta_j.

   - The cos(6Δ) term encodes a 6-fold (60°) alignment preference.
   - The cos(5Δ) term encodes a 5-fold (72°) golden alignment preference.
   - Competing 5- and 6-fold preferences lead to frustrated, quasi-crystal-like order.

3) Perform gradient descent on {theta_i} to minimize M.

4) From the relaxed configuration, build an emergent adjacency matrix

       S_ij = cos(6Δ_ij) + cos(5Δ_ij)
       W_ij = max(0, S_ij)

   and keep only the strongest edges to define an unweighted graph A_int.

5) Build the Laplacian L_int from A_int.

6) Plug this L_int into the operator-first flavor machinery:

   - extract a 3-mode "generation triad" from its spectrum,
   - build F_base(λ), integer-charge hierarchies F_s,
   - apply golden P_phi and Cabibbo C_12 to get CKM & PMNS,
   - compute rough chi^2 vs SM-inspired targets.

This is still a toy, but now the internal graph is *emergent from operator-like rules*,
not chosen a priori as fib2d or 24-cell.
"""
import itertools

def search_best_lepton_regions(
    gen_vecs,
    regions,
    U_geom_u, U_geom_d,
    F_u, F_d, F_e, F_n,
    P_phi_12, P_phi_23, C_12,
    N_SOLAR=36, N_REACTOR=45
):
    """
    Brute-force search over all permutations of the 3 regions for
    charged leptons and neutrinos, keeping:
      - up/down geometry fixed (U_geom_u, U_geom_d),
      - masses F_s fixed,
      - golden/Cabibbo/neutrino-dressing operators fixed.

    Returns (best_assign_e, best_assign_nu, best_chi2, best_results),
    where best_results includes the mixing matrices and angles.
    """
    R0, R1, R2 = regions
    region_list = [R0, R1, R2]
    perms = list(itertools.permutations(range(3)))

    best_chi2 = None
    best_assign_e = None
    best_assign_nu = None
    best_dat = None

    for pe in perms:
        for pn in perms:
            assign_e  = [region_list[i] for i in pe]
            assign_nu = [region_list[i] for i in pn]

            U_geom_e  = build_geometric_unitary(gen_vecs, assign_e)
            U_geom_nu = build_geometric_unitary(gen_vecs, assign_nu)

            U_geom = {
                "u":  U_geom_u,
                "d":  U_geom_d,
                "e":  U_geom_e,
                "nu": U_geom_nu,
            }

            sector_bases = build_sector_bases(
                P_phi_12, P_phi_23, C_12,
                U_geom,
                use_neutrino_dressing=True,
                N_SOLAR=N_SOLAR,
                N_REACTOR=N_REACTOR
            )

            U_L_u,  U_R_u  = sector_bases["u"]
            U_L_d,  U_R_d  = sector_bases["d"]
            U_L_e,  U_R_e  = sector_bases["e"]
            U_L_nu, U_R_nu = sector_bases["nu"]

            # Mixing matrices
            V_ckm  = mixing_matrix(U_L_u, U_L_d)
            U_pmns = mixing_matrix(U_L_e, U_L_nu)

            theta12_q, theta23_q, theta13_q = mixing_angles_from_U(V_ckm)
            theta12_l, theta23_l, theta13_l = mixing_angles_from_U(U_pmns)

            # Mass ratios from F_s (unchanged per iteration)
            mu_mt, mc_mt   = mass_ratios(F_u)
            md_mb, ms_mb   = mass_ratios(F_d)
            me_mt, mmu_mt  = mass_ratios(F_e)

            obs = compute_observables(
                mu_mt, mc_mt, md_mb, ms_mb, me_mt, mmu_mt,
                theta12_q, theta23_q, theta13_q,
                theta12_l, theta23_l, theta13_l
            )
            chi2_val, chi2_details = chi2(obs, TARGETS)

            if (best_chi2 is None) or (chi2_val < best_chi2):
                best_chi2 = chi2_val
                best_assign_e = pe
                best_assign_nu = pn
                best_dat = {
                    "chi2": chi2_val,
                    "chi2_details": chi2_details,
                    "V_ckm": V_ckm,
                    "U_pmns": U_pmns,
                    "angles_q": (theta12_q, theta23_q, theta13_q),
                    "angles_l": (theta12_l, theta23_l, theta13_l),
                }

    return best_assign_e, best_assign_nu, best_chi2, best_dat
# ----------------------------------------------------------------------
# 1. Misalignment functional M[theta] and its gradient
# ----------------------------------------------------------------------

def misalignment_energy(theta, w6=1.0, w5=1.0, J=None):
    """
    Compute total misalignment energy M[theta].

    theta: array shape (N,)
    J: optional coupling matrix shape (N,N); if None, J_ij = 1/N.
    """
    theta = theta.reshape(-1)
    N = theta.shape[0]
    if J is None:
        # uniform coupling J_ij = 1/N (not critical in this toy)
        J = np.ones((N, N), dtype=float) / N

    # pairwise differences Δ_ij
    # we can vectorize using broadcasting
    dtheta = theta[:, None] - theta[None, :]  # shape (N, N)

    E6 = 1.0 - np.cos(6.0 * dtheta)
    E5 = 1.0 - np.cos(5.0 * dtheta)

    M = 0.5 * np.sum(J * (w6 * E6 + w5 * E5))  # 1/2 to avoid double-count
    return M

def misalignment_grad(theta, w6=1.0, w5=1.0, J=None):
    """
    Gradient dM/dtheta_i.

    d/dθ_i (1 - cos(kΔ_ij)) = k sin(kΔ_ij), where Δ_ij = θ_i - θ_j.

    So

        ∂M/∂θ_i = sum_j J_ij [ w6*6 sin(6Δ_ij) + w5*5 sin(5Δ_ij) ].
    """
    theta = theta.reshape(-1)
    N = theta.shape[0]
    if J is None:
        J = np.ones((N, N), dtype=float) / N

    dtheta = theta[:, None] - theta[None, :]  # shape (N, N)
    term6 = 6.0 * np.sin(6.0 * dtheta)
    term5 = 5.0 * np.sin(5.0 * dtheta)

    # sum over j: (J_ij * (w6*term6 + w5*term5))
    grad = np.sum(J * (w6 * term6 + w5 * term5), axis=1)
    return grad

def relax_phases(N=200, n_steps=500, eta=0.01, w6=1.0, w5=1.0, random_seed=42):
    """
    Perform gradient descent on M[theta] starting from random phases.
    Returns final theta and a history of energies.
    """
    rng = np.random.default_rng(random_seed)
    theta = rng.uniform(0.0, 2.0 * math.pi, size=N)
    J = None  # uniform couplings

    energy_hist = []

    for step in range(n_steps):
        M = misalignment_energy(theta, w6=w6, w5=w5, J=J)
        energy_hist.append(M)
        grad = misalignment_grad(theta, w6=w6, w5=w5, J=J)
        theta -= eta * grad
        # wrap back into [0, 2π)
        theta = np.mod(theta, 2.0 * math.pi)

    return theta, np.array(energy_hist)


# ----------------------------------------------------------------------
# 2. Build emergent adjacency and Laplacian from relaxed phases
# ----------------------------------------------------------------------
def build_geometric_regions(theta: np.ndarray, n_regions: int = 3):
    """
    Partition sites into n_regions contiguous blocks in phase-order.
    This uses only the emergent phase field: no coordinates assumed.
    """
    theta = theta.reshape(-1)
    N = theta.shape[0]
    order = np.argsort(theta)  # sites sorted by phase
    # Split the sorted list into n_regions nearly equal chunks
    base = N // n_regions
    extra = N % n_regions
    regions = []
    start = 0
    for r in range(n_regions):
        size = base + (1 if r < extra else 0)
        idx = order[start:start+size]
        regions.append(idx)
        start += size
    return regions  # list of arrays of site indices

def build_geometric_unitary(gen_vecs: np.ndarray, region_list):
    """
    Given:
      gen_vecs: shape (N_sites, 3) = eigenvectors for the generation triad
      region_list: list of 3 index arrays (sites in each region for this sector)

    Construct 3 vectors in generation space by summing gen_vecs over each region,
    then orthonormalize them to get a 3x3 unitary-ish matrix U_geom.
    """
    cols = []
    for inds in region_list:
        # sum over sites in this region
        v = np.sum(gen_vecs[inds, :], axis=0)
        cols.append(v)
    M = np.stack(cols, axis=1)  # shape (3,3) with each col a vector in generation space

    # QR decomposition to orthonormalize columns
    Q, R = np.linalg.qr(M)
    # Optional: enforce det(Q) ~ +1 by flipping a column sign if needed
    if np.linalg.det(Q) < 0:
        Q[:, 0] = -Q[:, 0]
    return Q  # unitary (up to numerical noise)

def build_emergent_adjacency(theta, w6=1.0, w5=1.0, keep_fraction=0.1):
    """
    From a relaxed configuration theta, build an emergent adjacency A_int.

    For each pair (i,j):

        Δ_ij = θ_i - θ_j
        S_ij = w6*cos(6Δ_ij) + w5*cos(5Δ_ij)
        W_ij = max(0, S_ij)

    Then keep only the top 'keep_fraction' of W_ij (i<j) as edges.
    """
    theta = theta.reshape(-1)
    N = theta.shape[0]

    dtheta = theta[:, None] - theta[None, :]  # (N,N)
    S = w6 * np.cos(6.0 * dtheta) + w5 * np.cos(5.0 * dtheta)
    W = np.maximum(0.0, S)

    # Zero out diagonal
    np.fill_diagonal(W, 0.0)

    # Threshold
    # flatten upper triangle (i<j), pick top fraction
    iu, ju = np.triu_indices(N, k=1)
    weights = W[iu, ju]
    if keep_fraction <= 0.0:
        keep_fraction = 0.1
    n_edges = max(1, int(keep_fraction * weights.size))

    # indices of top weights
    top_idx = np.argpartition(weights, -n_edges)[-n_edges:]
    mask = np.zeros_like(weights, dtype=bool)
    mask[top_idx] = True

    # build adjacency
    A = np.zeros((N, N), dtype=float)
    A[iu[mask], ju[mask]] = 1.0
    A[ju[mask], iu[mask]] = 1.0

    return A

def laplacian_from_adjacency(A: np.ndarray) -> np.ndarray:
    D = np.diag(A.sum(axis=1))
    return D - A


# ----------------------------------------------------------------------
# 3. Operator-first flavor machinery (reused structure)
# ----------------------------------------------------------------------

def spectral_triad(L_int: np.ndarray):
    """
    Extract a 3-mode generation triad from L_int:
    the three lowest nonzero eigenvalues.
    """
    eigvals, eigvecs = np.linalg.eigh(L_int)
    eps = 1e-10
    nonzero_indices = np.where(eigvals > eps)[0]
    gen_indices = nonzero_indices[:3]
    lam_gen = eigvals[gen_indices]
    return lam_gen, gen_indices, eigvals

def base_kernel(lams: np.ndarray, alpha: float = 3.0, form: str = "lambda_sq") -> np.ndarray:
    """
    Spectral kernel F(lambda). We keep the same choices as before:
      "lambda_sq":  F = exp(-alpha * lambda^2)
    """
    if form == "lambda_sq":
        return np.exp(-alpha * (lams ** 2))
    elif form == "lambda":
        return np.exp(-alpha * lams)
    else:
        raise ValueError(f"Unknown kernel form: {form}")

def build_sector_charges():
    """
    Integer charges q_{s,g} for sector+generation hierarchies.

    Indices g = 0,1,2 correspond to the three internal modes in lam_gen
    (here ~[0.98, 1.82, 2.0]). Physical generations (1st,2nd,3rd) are
    determined by sorting the resulting F_s.

    These q_{s,g} are small integers chosen so that, given the fixed
    emergent F_base(lambda_gen), the sorted mass ratios (m1/m3, m2/m3)
    in each sector approximate the observed SM hierarchies:

      - Up:   mu/mt ~ 2.2e-5, mc/mt ~ 7.5e-3
      - Down: md/mb ~ 1.1e-3, ms/mb ~ 2.2e-2
      - E:    me/mtau ~ 2.9e-4, mmu/mtau ~ 5.9e-2

    No continuous Yukawa parameters are introduced; only discrete
    integer exponents acting on the emergent 3-mode triad.
    """
    sector_charges_gen = {
        # Up-type quarks
        "u":  np.array([4.0, 8.0, 0.0]),
        # Down-type quarks
        "d":  np.array([5.0, 5.0, 0.0]),
        # Charged leptons
        "e":  np.array([4.0, 0.0, 3.0]),
        # Neutrinos (kept as a simple, more-suppressed pattern for now)
        "nu": np.array([6.0, 5.0, 4.0]),
    }
    return sector_charges_gen

def sector_weights(F_base: np.ndarray, q_vec: np.ndarray, beta: float = 1.0) -> np.ndarray:
    return F_base * np.exp(-beta * q_vec)

def rot12(theta: float) -> np.ndarray:
    c = math.cos(theta)
    s = math.sin(theta)
    return np.array([
        [ c,  s, 0.0],
        [-s,  c, 0.0],
        [0.0, 0.0, 1.0]
    ], dtype=complex)

def rot23(theta: float) -> np.ndarray:
    c = math.cos(theta)
    s = math.sin(theta)
    return np.array([
        [1.0, 0.0, 0.0],
        [0.0,  c,  s],
        [0.0, -s,  c]
    ], dtype=complex)

def rot13(theta: float) -> np.ndarray:
    c = math.cos(theta)
    s = math.sin(theta)
    return np.array([
        [ c, 0.0,  s],
        [0.0, 1.0, 0.0],
        [-s, 0.0,  c]
    ], dtype=complex)

def build_generation_operators(phi_order=5, cab_denom=28):
    theta_phi = 2.0 * math.pi / float(phi_order)
    theta_C   = 2.0 * math.pi / float(cab_denom)
    P_phi_12 = rot12(theta_phi)
    P_phi_23 = rot23(theta_phi)
    C_12     = rot12(theta_C)
    return P_phi_12, P_phi_23, C_12, theta_phi, theta_C

def build_sector_bases(P_phi_12, P_phi_23, C_12,
                       U_geom,
                       use_neutrino_dressing=True,
                       N_SOLAR=36, N_REACTOR=45):
    I3 = np.eye(3, dtype=complex)

    # Quarks
    U_L_u  = U_geom["u"]  @ P_phi_12
    U_L_d  = U_geom["d"]  @ P_phi_12 @ C_12

    # Charged leptons
    U_L_e  = U_geom["e"]

    # Neutrinos
    if use_neutrino_dressing:
        theta12_nu = 2.0 * math.pi / float(N_SOLAR)
        theta13_nu = 2.0 * math.pi / float(N_REACTOR)
        U_L_nu = U_geom["nu"] @ rot12(theta12_nu) @ P_phi_23 @ rot13(theta13_nu)
    else:
        U_L_nu = U_geom["nu"] @ P_phi_23

    U_R_u  = I3
    U_R_d  = I3
    U_R_e  = I3
    U_R_nu = I3

    return {
        "u":  (U_L_u,  U_R_u),
        "d":  (U_L_d,  U_R_d),
        "e":  (U_L_e,  U_R_e),
        "nu": (U_L_nu, U_R_nu),
    }

def yukawa_from_F_and_UL(F_s: np.ndarray, U_L: np.ndarray, U_R: np.ndarray) -> np.ndarray:
    F_diag = np.diag(F_s.astype(complex))
    return U_L.conj().T @ F_diag @ U_R

def mixing_matrix(U_L_up: np.ndarray, U_L_down: np.ndarray) -> np.ndarray:
    return U_L_up.conj().T @ U_L_down

def mixing_angles_from_U(U: np.ndarray):
    U_abs = np.abs(U)
    s13 = U_abs[0, 2]
    c13 = math.sqrt(max(0.0, 1.0 - s13**2))
    if c13 < 1e-12:
        return 0.0, 0.0, math.pi / 2.0
    s12 = U_abs[0, 1] / c13
    s23 = U_abs[1, 2] / c13
    s12 = max(-1.0, min(1.0, s12))
    s23 = max(-1.0, min(1.0, s23))
    theta12 = math.asin(s12)
    theta23 = math.asin(s23)
    theta13 = math.asin(s13)
    return theta12, theta23, theta13

def mass_ratios(F_s: np.ndarray):
    s_sorted = np.sort(F_s)
    m1, m2, m3 = s_sorted
    return m1 / m3, m2 / m3

def compute_observables(mu_mt, mc_mt, md_mb, ms_mb, me_mt, mmu_mt,
                        theta12_q, theta23_q, theta13_q,
                        theta12_l, theta23_l, theta13_l):
    return {
        "mu/mt":      mu_mt,
        "mc/mt":      mc_mt,
        "md/mb":      md_mb,
        "ms/mb":      ms_mb,
        "me/mtau":    me_mt,
        "mmu/mtau":   mmu_mt,
        "theta12_q":  theta12_q,
        "theta23_q":  theta23_q,
        "theta13_q":  theta13_q,
        "theta12_l":  theta12_l,
        "theta23_l":  theta23_l,
        "theta13_l":  theta13_l,
    }

TARGETS = {
    "mu/mt":     2.2e-5,
    "mc/mt":     7.5e-3,
    "md/mb":     1.1e-3,
    "ms/mb":     2.2e-2,
    "me/mtau":   2.9e-4,
    "mmu/mtau":  0.059,
    "theta12_q": 0.227,
    "theta23_q": 0.041,
    "theta13_q": 0.0036,
    "theta12_l": 0.584,
    "theta23_l": 0.785,
    "theta13_l": 0.150,
}

def chi2(observables, targets):
    chi2_total = 0.0
    details = []
    ratio_keys = ["mu/mt", "mc/mt", "md/mb", "ms/mb", "me/mtau", "mmu/mtau"]
    angle_keys = ["theta12_q", "theta23_q", "theta13_q",
                  "theta12_l", "theta23_l", "theta13_l"]

    for k in ratio_keys:
        m = observables[k]
        t = targets[k]
        if m <= 0 or t <= 0:
            continue
        logm = math.log10(m)
        logt = math.log10(t)
        sigma_log = 0.3
        contrib = ((logm - logt) / sigma_log)**2
        chi2_total += contrib
        details.append((k, m, t, contrib))

    for k in angle_keys:
        m = observables[k]
        t = targets[k]
        sigma = 0.2
        contrib = ((m - t) / sigma)**2
        chi2_total += contrib
        details.append((k, m, t, contrib))

    return chi2_total, details

def largest_connected_component(A):
    N = A.shape[0]
    visited = np.zeros(N, dtype=bool)
    components = []

    for i in range(N):
        if not visited[i]:
            stack = [i]
            comp = []
            visited[i] = True
            while stack:
                u = stack.pop()
                comp.append(u)
                neighbors = np.where(A[u] > 0)[0]
                for v in neighbors:
                    if not visited[v]:
                        visited[v] = True
                        stack.append(v)
            components.append(comp)

    # pick largest
    comp_sizes = [len(c) for c in components]
    largest_idx = np.argmax(comp_sizes)
    nodes = np.array(components[largest_idx], dtype=int)

    # induced subgraph
    A_sub = A[np.ix_(nodes, nodes)]
    return A_sub, nodes
# ----------------------------------------------------------------------
# 4. Main: emergent graph → L_int → flavor operators
# ----------------------------------------------------------------------

def main():
    # Step 1: relax phases under misalignment functional
    N = 200
    theta_final, energy_hist = relax_phases(
        N=N,
        n_steps=600,
        eta=0.01,
        w6=1.0,
        w5=1.0,
        random_seed=42
    )
    print("Relaxation complete.")
    print(f"Final misalignment energy: {energy_hist[-1]:.6f}")
    print()

    # Step 2: build emergent adjacency and Laplacian
    A_int_full = build_emergent_adjacency(
        theta_final,
        w6=1.0,
        w5=1.0,
        keep_fraction=0.05
    )
    A_int, nodes = largest_connected_component(A_int_full)
    L_int = laplacian_from_adjacency(A_int)

    # Spectrum and generation triad
    lam_gen, gen_indices, eigvals = spectral_triad(L_int)
    F_base = base_kernel(lam_gen, alpha=3.0, form="lambda_sq")

    print("=== Emergent internal graph ===")
    print(f"Number of sites: {A_int.shape[0]}")
    print("First 10 eigenvalues of L_int:")
    print(eigvals[:10])
    print()
    print("Generation eigenvalue indices:", gen_indices)
    print("Generation triad lam_gen:", lam_gen)
    print("Base kernel F_base(lam_gen):", F_base)
    print()

    # Generation eigenvectors restricted to the triad
    eigvals_full, eigvecs_full = np.linalg.eigh(L_int)
    gen_vecs = eigvecs_full[:, gen_indices]  # shape (N_sub, 3)

    # Build 3 geometric regions from the emergent phase field,
    # restricted to the nodes in the largest connected component.
    theta_sub = theta_final[nodes]
    regions = build_geometric_regions(theta_sub, n_regions=3)
    R0, R1, R2 = regions

    # Quarks: share the same geometric basis so CKM stays Cabibbo-like
    assign_u = [R0, R1, R2]
    assign_d = [R0, R1, R2]
    U_geom_u = build_geometric_unitary(gen_vecs, assign_u)
    U_geom_d = build_geometric_unitary(gen_vecs, assign_d)

    # Sector charges & F_s (fixed integer Q pattern)
    sector_charges_gen = build_sector_charges()
    F_u = sector_weights(F_base, sector_charges_gen["u"],  beta=1.0)
    F_d = sector_weights(F_base, sector_charges_gen["d"],  beta=1.0)
    F_e = sector_weights(F_base, sector_charges_gen["e"],  beta=1.0)
    F_n = sector_weights(F_base, sector_charges_gen["nu"], beta=1.0)

    print("=== Yukawa-like mass scales F_s ===")
    print("Up-type (F_u):        ", F_u)
    print("Down-type (F_d):      ", F_d)
    print("Charged leptons (F_e):", F_e)
    print("Neutrino (F_n):       ", F_n)
    print()

    # Generation-space operators (golden + Cabibbo)
    P_phi_12, P_phi_23, C_12, theta_phi, theta_C = build_generation_operators(
        phi_order=5,
        cab_denom=28
    )

    # Step 3: search over geometric assignments for e and nu
    best_pe, best_pn, best_chi2, best_dat = search_best_lepton_regions(
        gen_vecs,
        regions,
        U_geom_u, U_geom_d,
        F_u, F_d, F_e, F_n,
        P_phi_12, P_phi_23, C_12,
        N_SOLAR=36, N_REACTOR=45
    )

    print("Best lepton region permutations:")
    print("  pe (e sectors)  =", best_pe)
    print("  pn (nu sectors) =", best_pn)
    print(f"Best total chi^2  ≈ {best_chi2:.2f}")
    print()

    # Reconstruct the best U_geom using that assignment
    region_list = [R0, R1, R2]
    assign_e  = [region_list[i] for i in best_pe]
    assign_nu = [region_list[i] for i in best_pn]

    U_geom = {
        "u":  U_geom_u,
        "d":  U_geom_d,
        "e":  build_geometric_unitary(gen_vecs, assign_e),
        "nu": build_geometric_unitary(gen_vecs, assign_nu),
    }

    # Build sector bases using both geometry and flavor operators
    sector_bases = build_sector_bases(
        P_phi_12, P_phi_23, C_12,
        U_geom,
        use_neutrino_dressing=True,
        N_SOLAR=36,
        N_REACTOR=45
    )

    U_L_u,  U_R_u  = sector_bases["u"]
    U_L_d,  U_R_d  = sector_bases["d"]
    U_L_e,  U_R_e  = sector_bases["e"]
    U_L_nu, U_R_nu = sector_bases["nu"]

    # Yukawa-like operators (not strictly needed for ratios, but kept for completeness)
    Y_u  = yukawa_from_F_and_UL(F_u, U_L_u,  U_R_u)
    Y_d  = yukawa_from_F_and_UL(F_d, U_L_d,  U_R_d)
    Y_e  = yukawa_from_F_and_UL(F_e, U_L_e,  U_R_e)
    Y_nu = yukawa_from_F_and_UL(F_n, U_L_nu, U_R_nu)

    # Mass ratios from F_s (eigenvalues of Y†Y will be very close to these)
    mu_mt, mc_mt   = mass_ratios(F_u)
    md_mb, ms_mb   = mass_ratios(F_d)
    me_mt, mmu_mt  = mass_ratios(F_e)

    print("Mass ratios (m1/m3, m2/m3) from F_s:")
    print(f"mu/mt:     {mu_mt:.3e}, mc/mt:     {mc_mt:.3e}")
    print(f"md/mb:     {md_mb:.3e}, ms/mb:     {ms_mb:.3e}")
    print(f"me/mtau:   {me_mt:.3e}, mmu/mtau:  {mmu_mt:.3e}")
    print()

    # Step 5: mixing matrices
    V_ckm  = mixing_matrix(U_L_u, U_L_d)
    U_pmns = mixing_matrix(U_L_e, U_L_nu)

    theta12_q, theta23_q, theta13_q = mixing_angles_from_U(V_ckm)
    theta12_l, theta23_l, theta13_l = mixing_angles_from_U(U_pmns)

    print("=== CKM-like mixing matrix (geometry + operator) ===")
    print(V_ckm)
    print(f"theta12_q ≈ {theta12_q:.3f} rad, theta23_q ≈ {theta23_q:.3f}, theta13_q ≈ {theta13_q:.3e}")
    print(f"(Cabibbo operator angle = 2π/28 ≈ {theta_C:.3f} rad)")
    print()

    print("=== PMNS-like mixing matrix (geometry + operator) ===")
    print(U_pmns)
    print(f"theta12_l ≈ {theta12_l:.3f} rad, theta23_l ≈ {theta23_l:.3f}, theta13_l ≈ {theta13_l:.3e}")
    print(f"(Golden operator angle = 2π/5 ≈ {theta_phi:.3f} rad)")
    print()

    # Step 6: chi^2 vs rough targets (recompute for the final configuration)
    obs = compute_observables(
        mu_mt, mc_mt, md_mb, ms_mb, me_mt, mmu_mt,
        theta12_q, theta23_q, theta13_q,
        theta12_l, theta23_l, theta13_l
    )
    chi2_value, chi2_details = chi2(obs, TARGETS)

    print("=== Observables vs rough targets ===")
    for k, m, t, contrib in chi2_details:
        print(f"{k:12s}: model={m:.3e}, target={t:.3e}, chi2_contrib={contrib:.2f}")
    print()
    print(f"Total chi^2 ≈ {chi2_value:.2f}")
    print()
    print("NOTES:")
    print("- The internal graph is emergent from the misalignment functional M[theta],")
    print("  which encodes 6-fold (C_360) and 5-fold (golden) alignment preferences.")
    print("- We then restrict to the largest connected component to define a single,")
    print("  coherent aether vacuum, and build its Laplacian L_int.")
    print("- The generation triad and F_base(lambda) come from the spectrum of L_int,")
    print("  sector hierarchies from discrete integer charges Q_{s,g}, and mixing")
    print("  from a combination of geometry-derived U_geom[s] and fixed operators")
    print("  P_phi (golden) and C_12 (Cabibbo).")
    print("- No random Yukawas or continuous per-sector fits are used; everything")
    print("  comes from the emergent graph, a universal kernel, integer exponents,")
    print("  and discrete 2π/n phase rotations.")

if __name__ == "__main__":
    main()

"""
RESULTS:

Relaxation complete.
Final misalignment energy: 99.972531

=== Emergent internal graph ===
Number of sites: 39
First 10 eigenvalues of L_int:
[7.04495820e-17 9.80951287e-01 1.81564639e+00 2.00000000e+00
 2.00000000e+00 2.00000000e+00 2.00000000e+00 2.00000000e+00
 4.38302718e+00 5.00000000e+00]

Generation eigenvalue indices: [1 2 3]
Generation triad lam_gen: [0.98095129 1.81564639 2.        ]
Base kernel F_base(lam_gen): [5.57545487e-02 5.06933717e-05 6.14421235e-06]

=== Yukawa-like mass scales F_s ===
Up-type (F_u):         [1.02118018e-03 1.70057317e-08 6.14421235e-06]
Down-type (F_d):       [3.75671194e-04 3.41569251e-07 6.14421235e-06]
Charged leptons (F_e): [1.02118018e-03 5.06933717e-05 3.05902321e-07]
Neutrino (F_n):        [1.38201709e-04 3.41569251e-07 1.12535175e-07]

Best lepton region permutations:
  pe (e sectors)  = (1, 2, 0)
  pn (nu sectors) = (2, 0, 1)
Best total chi^2  ≈ 11.90

Mass ratios (m1/m3, m2/m3) from F_s:
mu/mt:     1.665e-05, mc/mt:     6.017e-03
md/mb:     9.092e-04, ms/mb:     1.636e-02
me/mtau:   2.996e-04, mmu/mtau:  4.964e-02

=== CKM-like mixing matrix (geometry + operator) ===
[[ 9.74927912e-01+0.j  2.22520934e-01+0.j -3.35683117e-17+0.j]
 [-2.22520934e-01+0.j  9.74927912e-01+0.j -5.38467983e-17+0.j]
 [-5.55111512e-17+0.j -5.55111512e-17+0.j  1.00000000e+00+0.j]]
theta12_q ≈ 0.224 rad, theta23_q ≈ 0.000, theta13_q ≈ 3.357e-17
(Cabibbo operator angle = 2π/28 ≈ 0.224 rad)

=== PMNS-like mixing matrix (geometry + operator) ===
[[-0.82852367+0.j  0.20511282+0.j  0.52103481+0.j]
 [-0.55830005+0.j -0.23112818+0.j -0.79679409+0.j]
 [-0.04300685+0.j -0.95105652+0.j  0.30600966+0.j]]
theta12_l ≈ 0.243 rad, theta23_l ≈ 1.204, theta13_l ≈ 5.481e-01
(Golden operator angle = 2π/5 ≈ 1.257 rad)

=== Observables vs rough targets ===
mu/mt       : model=1.665e-05, target=2.200e-05, chi2_contrib=0.16
mc/mt       : model=6.017e-03, target=7.500e-03, chi2_contrib=0.10
md/mb       : model=9.092e-04, target=1.100e-03, chi2_contrib=0.08
ms/mb       : model=1.636e-02, target=2.200e-02, chi2_contrib=0.18
me/mtau     : model=2.996e-04, target=2.900e-04, chi2_contrib=0.00
mmu/mtau    : model=4.964e-02, target=5.900e-02, chi2_contrib=0.06
theta12_q   : model=2.244e-01, target=2.270e-01, chi2_contrib=0.00
theta23_q   : model=5.385e-17, target=4.100e-02, chi2_contrib=0.04
theta13_q   : model=3.357e-17, target=3.600e-03, chi2_contrib=0.00
theta12_l   : model=2.427e-01, target=5.840e-01, chi2_contrib=2.91
theta23_l   : model=1.204e+00, target=7.850e-01, chi2_contrib=4.39
theta13_l   : model=5.481e-01, target=1.500e-01, chi2_contrib=3.96

Total chi^2 ≈ 11.90
"""

# ================================
# Internal Hilbert space & D_F
# ================================
import numpy as np
from typing import List, Tuple, Dict

SECTORS = ["u", "d", "e", "nu"]
SECTOR_INDEX: Dict[str, int] = {s: i for i, s in enumerate(SECTORS)}
N_GEN = 3
# Color multiplicities (degeneracies). In this toy, we don't explicitly
# tensor out color, we just keep track of the full 24-dim per chirality.
SECTOR_NC = {"u": 3, "d": 3, "e": 1, "nu": 1}


def dim_per_chirality() -> int:
    """
    Dimension of H_L or H_R (one chirality).

    We treat color multiplicities as degeneracy, not as an explicit tensor
    factor for now, but the total number of internal states per chirality
    still comes out as:

        dim(H_L) = N_GEN * sum_s N_c(s) = 3 * (3+3+1+1) = 24

    This matches your emergent-5 setup.
    """
    return N_GEN * sum(SECTOR_NC[s] for s in SECTORS)  # 24


def flavor_block_offsets() -> Dict[str, int]:
    """
    Return offsets (within the *generation subspace*) for each sector's 3×3
    generation block in a 12×12 layout:

      [u_g1,u_g2,u_g3,
       d_g1,d_g2,d_g3,
       e_g1,e_g2,e_g3,
       nu_g1,nu_g2,nu_g3]

    We only care about the leading 12 entries as "generation space".
    The remaining 12 (per chirality) are currently unused / reserved
    for future color-explicit extensions.
    """
    off: Dict[str, int] = {}
    off["u"]  = 0
    off["d"]  = 3
    off["e"]  = 6
    off["nu"] = 9
    return off


# -------------------------------------------------------------------
# F-based D_F builder (diagonal Yukawas) – optional fallback
# -------------------------------------------------------------------
def build_internal_DF(F_u: np.ndarray,
                      F_d: np.ndarray,
                      F_e: np.ndarray,
                      F_n: np.ndarray) -> np.ndarray:
    """
    Build the finite Dirac operator D_F in block form:

      D_F = [[ 0, Y^\dagger ],
             [ Y, 0       ]]

    where Y is a 24×24 block that is block-diagonal in sector space and
    uses 3×3 diagonal generation Yukawas diag(F_s) in a 12×12
    generation-space layout (color folded in as degeneracy).

    H_F ≃ H_L ⊕ H_R,  dim(H_L) = dim(H_R) = 24, dim(H_F) = 48.
    """
    # Sanity checks
    for name, F in [("F_u", F_u), ("F_d", F_d), ("F_e", F_e), ("F_n", F_n)]:
        F_arr = np.asarray(F, dtype=float)
        if F_arr.shape != (3,):
            raise ValueError(f"{name} must be a length-3 array, got shape {F_arr.shape}.")

    dpc = dim_per_chirality()  # 24
    dimH = 2 * dpc             # 48

    # 12×12 generation-space Yukawa core; then embedded into 24×24 per chirality
    Y_gen = np.zeros((12, 12), dtype=complex)

    Y_u  = np.diag(F_u)
    Y_d  = np.diag(F_d)
    Y_e  = np.diag(F_e)
    Y_nu = np.diag(F_n)

    gen_off = flavor_block_offsets()

    def insert_sector_Y(sector: str, Y_s: np.ndarray) -> None:
        off = gen_off[sector]
        Y_gen[off:off+3, off:off+3] = Y_s

    insert_sector_Y("u",  Y_u)
    insert_sector_Y("d",  Y_d)
    insert_sector_Y("e",  Y_e)
    insert_sector_Y("nu", Y_nu)

    # Embed this 12×12 generation block into 24×24 per chirality.
    Y_block = np.zeros((dpc, dpc), dtype=complex)
    Y_block[:12, :12] = Y_gen

    # Assemble D_F
    D_F = np.zeros((dimH, dimH), dtype=complex)
    D_F[:dpc,  dpc:] = Y_block.conj().T
    D_F[dpc:,  :dpc] = Y_block

    return D_F


# -------------------------------------------------------------------
# Y-based D_F builder (full Yukawa matrices, with mixing)
# -------------------------------------------------------------------
def build_internal_DF_from_Y(Y_u: np.ndarray,
                             Y_d: np.ndarray,
                             Y_e: np.ndarray,
                             Y_nu: np.ndarray) -> np.ndarray:
    """
    Build the finite Dirac operator D_F in block form:

      D_F = [[ 0, Y^\dagger ],
             [ Y, 0        ]]

    where Y is a 24×24 block, block-diagonal in sector space, with
    3×3 Yukawa matrices per sector (not necessarily diagonal):

      Y_gen = diag( Y_u, Y_d, Y_e, Y_nu ) in the leading 12×12 generation space.

    The remaining 12 entries per chirality are unused (reserved for future
    explicit color structure). For now they are set to zero.

    H_F ≃ H_L ⊕ H_R, dim(H_L)=dim(H_R)=24, dim(H_F)=48.
    """
    # Sanity checks on shapes
    for name, Y in [("Y_u", Y_u), ("Y_d", Y_d), ("Y_e", Y_e), ("Y_nu", Y_nu)]:
        Y_arr = np.asarray(Y, dtype=complex)
        if Y_arr.shape != (3, 3):
            raise ValueError(f"{name} must be a 3×3 matrix, got shape {Y_arr.shape}.")

    dpc = dim_per_chirality()  # 24
    dimH = 2 * dpc             # 48

    # Build 12×12 generation block first
    Y_gen = np.zeros((12, 12), dtype=complex)
    gen_off = flavor_block_offsets()

    def insert_sector_Y(sector: str, Y_s: np.ndarray) -> None:
        off = gen_off[sector]
        Y_gen[off:off+3, off:off+3] = Y_s

    insert_sector_Y("u",  Y_u)
    insert_sector_Y("d",  Y_d)
    insert_sector_Y("e",  Y_e)
    insert_sector_Y("nu", Y_nu)

    # Embed generation block into 24×24 per chirality
    Y_block = np.zeros((dpc, dpc), dtype=complex)
    Y_block[:12, :12] = Y_gen

    # Assemble D_F
    D_F = np.zeros((dimH, dimH), dtype=complex)
    D_F[:dpc,  dpc:] = Y_block.conj().T
    D_F[dpc:,  :dpc] = Y_block

    return D_F


# ================================
# Real structure, grading, algebra basis
# ================================

def build_swap_LR(dim_left: int) -> np.ndarray:
    """
    Build the swap matrix S on H_F = H_L ⊕ H_R, where dim(H_L) = dim(H_R) = dim_left.
    Acts as:
      S ( ψ_L, ψ_R ) = ( ψ_R, ψ_L )
    """
    S = np.zeros((2*dim_left, 2*dim_left), dtype=complex)
    S[:dim_left, dim_left:] = np.eye(dim_left)
    S[dim_left:, :dim_left] = np.eye(dim_left)
    return S


def build_gamma_F(dim_left: int) -> np.ndarray:
    """
    Grading operator γ_F with eigenvalue -1 on H_L and +1 on H_R.
    """
    g = np.zeros((2*dim_left, 2*dim_left), dtype=complex)
    g[:dim_left, :dim_left] = -np.eye(dim_left)
    g[dim_left:, dim_left:] =  np.eye(dim_left)
    return g


def build_sector_projectors() -> Dict[str, np.ndarray]:
    """
    Build sector projectors P_sector_s acting on H_F = H_L ⊕ H_R.
    Each P_sector_s is diagonal and selects the (sector,gen,chirality) subspace
    corresponding to that sector (u,d,e,nu) in the 12×12 generation subspace,
    duplicated on L and R.
    """
    dpc = dim_per_chirality()
    dimH = 2 * dpc

    gen_off = flavor_block_offsets()

    P: Dict[str, np.ndarray] = {}
    for s in SECTORS:
        P_s = np.zeros((dimH, dimH), dtype=complex)
        off = gen_off[s]
        # Act the same on L and R (block-diagonal), and only on the first 12 gen slots
        P_s[off:off+3, off:off+3] = np.eye(3)
        P_s[dpc + off:dpc + off + 3, dpc + off:dpc + off + 3] = np.eye(3)
        P[s] = P_s

    return P  # dict with keys "u","d","e","nu"


def build_Q_sector() -> np.ndarray:
    """
    Build a simple 'sector charge' diagonal operator Q_sector which distinguishes
    u,d,e,nu sectors but is generation-blind.
    Example charges:
      q_u = 2, q_d = 1, q_e = 0, q_nu = -1
    """
    dpc = dim_per_chirality()
    dimH = 2 * dpc

    gen_off = flavor_block_offsets()
    charges = {"u": 2.0, "d": 1.0, "e": 0.0, "nu": -1.0}

    Q = np.zeros((dimH, dimH), dtype=complex)
    for s in SECTORS:
        off = gen_off[s]
        q   = charges[s]
        Q[off:off+3, off:off+3] = q * np.eye(3)
        Q[dpc + off:dpc + off + 3, dpc + off:dpc + off + 3] = q * np.eye(3)

    return Q


def build_internal_algebra_ops() -> Tuple[List[np.ndarray], List[str]]:
    """
    Build a small basis of algebra elements A_F acting on H_F:
      - I (identity)
      - Q_sector (diagonal sector 'charge')
      - P_sector_u, P_sector_d, P_sector_e, P_sector_nu (sector projectors)

    This is a commutative algebra in this toy (we are not yet including full SU(3)).
    """
    dpc = dim_per_chirality()
    dimH = 2 * dpc

    I = np.eye(dimH, dtype=complex)
    Q = build_Q_sector()
    P = build_sector_projectors()

    ops: List[np.ndarray] = [I, Q, P["u"], P["d"], P["e"], P["nu"]]
    labels: List[str] = ["I", "Q_sector", "P_sector_u", "P_sector_d",
                         "P_sector_e", "P_sector_nu"]

    return ops, labels


# ================================
# NCG condition tests
# ================================

def J_action_from_swap(S: np.ndarray, M: np.ndarray) -> np.ndarray:
    """
    Implement J M J^{-1} = S * M^* * S^T, where S is the L/R swap.
    """
    return S @ M.conj() @ S.T


def test_first_order_condition(D_F: np.ndarray,
                               ops: List[np.ndarray],
                               labels: List[str],
                               eps: float = 1e-12) -> None:
    """
    First-order condition:
      [[D_F, a], J_F b J_F^{-1}] = 0
    for all a,b in algebra.
    """
    n = D_F.shape[0]
    assert D_F.shape == (n, n)
    S = build_swap_LR(dim_left=n//2)

    print("=== First-order condition test ===")
    max_norm = 0.0
    good_pairs = []

    for i, a in enumerate(ops):
        Da = D_F @ a - a @ D_F
        for j, b in enumerate(ops):
            b_tilde = J_action_from_swap(S, b)
            comm2 = Da @ b_tilde - b_tilde @ Da
            norm = np.linalg.norm(comm2, ord='fro')

            if norm > max_norm:
                max_norm = norm
            if norm < eps:
                good_pairs.append((labels[i], labels[j], norm))

    print(f"Max Frobenius norm over all pairs (a,b): {max_norm:.3e}")
    if good_pairs:
        print(f"Pairs with norm < {eps:.1e}:")
        for la, lb, nrm in good_pairs:
            print(f"  (a={la:>10s}, b={lb:>10s}) → ||[[D,a],J b J^-1]||_F = {nrm:.3e}")
    else:
        print(f"No pairs with norm < {eps:.1e}")
    print()


def test_zero_order_condition(ops: List[np.ndarray],
                              labels: List[str],
                              eps: float = 1e-12) -> None:
    """
    Zero-order condition:
      [a, J_F b J_F^{-1}] = 0
    for all a,b in algebra.
    """
    n = ops[0].shape[0]
    S = build_swap_LR(dim_left=n//2)

    print("=== Zero-order condition test ===")
    max_norm = 0.0
    bad_pairs = []

    for i, a in enumerate(ops):
        for j, b in enumerate(ops):
            b_tilde = J_action_from_swap(S, b)
            comm = a @ b_tilde - b_tilde @ a
            norm = np.linalg.norm(comm, ord='fro')
            if norm > max_norm:
                max_norm = norm
            if norm > eps:
                bad_pairs.append((labels[i], labels[j], norm))

    print(f"Max Frobenius norm over all pairs (a,b): {max_norm:.3e}")
    if bad_pairs:
        print("Pairs with significant violation:")
        for la, lb, nrm in bad_pairs:
            print(f"  (a={la:>10s}, b={lb:>10s}) → ||[a, J b J^-1]||_F = {nrm:.3e}")
    else:
        print(f"All pairs satisfy [a, J b J^-1] ≈ 0 within eps={eps:.1e}")
    print()


def test_grading_and_reality(D_F: np.ndarray,
                             ops: List[np.ndarray],
                             labels: List[str]) -> None:
    """
    - Check γ_F anticommutes with D_F and commutes with A_F.
    - Check J_F^2 = 1 (as implemented by swap).
    - Detect the KO-dimension sign via:
        J D_F J^{-1} = ± D_F
    """
    n = D_F.shape[0]
    dpc = n // 2
    gamma_F = build_gamma_F(dpc)
    S = build_swap_LR(dpc)

    print("=== Grading & reality tests ===")
    # γ_F anti-commutes with D_F
    anti = gamma_F @ D_F + D_F @ gamma_F
    print(f"||{{gamma_F, D_F}}||_F = {np.linalg.norm(anti, ord='fro'):.3e}")

    # γ_F commutes with algebra
    max_comm_gamma = 0.0
    for a in ops:
        comm_ga = gamma_F @ a - a @ gamma_F
        max_comm_gamma = max(max_comm_gamma, np.linalg.norm(comm_ga, ord='fro'))
    print(f"max ||[gamma_F, a]||_F over a∈A_F = {max_comm_gamma:.3e}")

    # J_F^2 = 1 (swap^2 = I)
    S2 = S @ S
    print(f"||S^2 - I||_F  (⇒ J_F^2 deviation) = {np.linalg.norm(S2 - np.eye(n), ord='fro'):.3e}")

    # J D_F J^{-1} vs ± D_F
    JDJ = S @ D_F.conj() @ S.T
    norm_plus  = np.linalg.norm(JDJ - D_F, ord='fro')
    norm_minus = np.linalg.norm(JDJ + D_F, ord='fro')

    print(f"||J D_F J^-1 - D_F||_F   = {norm_plus:.3e}")
    print(f"||J D_F J^-1 + D_F||_F   = {norm_minus:.3e}")

    if norm_plus < 1e-12 and norm_minus > norm_plus:
        print("→ KO-sign: J D_F J^-1 = + D_F (J-even Dirac operator)")
    elif norm_minus < 1e-12 and norm_plus > norm_minus:
        print("→ KO-sign: J D_F J^-1 = - D_F (J-odd Dirac operator)")
    else:
        print("→ KO-sign ambiguous or not clean at numerical precision.")
    print()


# ================================
# Minimal example / entry point
# ================================

def example_Fs() -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Provide a simple example set of F_s triads so this file can be run
    standalone. In your full emergent model, replace these with the
    F_u, F_d, F_e, F_n you compute from the internal graph + Q.
    """
    # Example: slightly hierarchical triad (not meant to match SM)
    F_base = np.array([0.05, 0.005, 0.0005], dtype=float)

    # Simple integer exponents per sector (toy, not SM):
    q_u  = np.array([0,  2, 4], dtype=float)
    q_d  = np.array([1,  3, 5], dtype=float)
    q_e  = np.array([2,  4, 6], dtype=float)
    q_nu = np.array([4,  6, 8], dtype=float)

    beta = 1.0

    def sector_weights(Fb: np.ndarray, q: np.ndarray, beta_val: float) -> np.ndarray:
        return Fb * np.exp(-beta_val * q)

    F_u  = sector_weights(F_base, q_u,  beta)
    F_d  = sector_weights(F_base, q_d,  beta)
    F_e  = sector_weights(F_base, q_e,  beta)
    F_n  = sector_weights(F_base, q_nu, beta)

    return F_u, F_d, F_e, F_n


def main() -> None:
    # Example Yukawa triads (replace with your emergent F_s in the full model)
    F_u, F_d, F_e, F_n = example_Fs()

    # For this minimal example, we build Y_s as simple diagonal matrices.
    # In your full emergent pipeline, REPLACE these with your actual emergent
    # Yukawa matrices, e.g. Y_u = U_L_u @ np.diag(F_u) @ U_R_u.conj().T, etc.
    Y_u  = np.diag(F_u)
    Y_d  = np.diag(F_d)
    Y_e  = np.diag(F_e)
    Y_nu = np.diag(F_n)

    # --- Build internal Dirac from full Yukawas ---
    D_F = build_internal_DF_from_Y(Y_u, Y_d, Y_e, Y_nu)

    # --- Build algebra basis (same as before) ---
    ops_A, labels_A = build_internal_algebra_ops()

    # --- Run NCG tests ---
    test_first_order_condition(D_F, ops_A, labels_A, eps=1e-12)
    test_zero_order_condition(ops_A, labels_A, eps=1e-12)
    test_grading_and_reality(D_F, ops_A, labels_A)


if __name__ == "__main__":
    main()

"""
RESULTS:
=== First-order condition test ===
Max Frobenius norm over all pairs (a,b): 0.000e+00
Pairs with norm < 1.0e-12:
  (a=         I, b=         I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         I, b=  Q_sector) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         I, b=P_sector_u) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         I, b=P_sector_d) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         I, b=P_sector_e) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         I, b=P_sector_nu) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=  Q_sector, b=         I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=  Q_sector, b=  Q_sector) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=  Q_sector, b=P_sector_u) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=  Q_sector, b=P_sector_d) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=  Q_sector, b=P_sector_e) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=  Q_sector, b=P_sector_nu) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_u, b=         I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_u, b=  Q_sector) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_u, b=P_sector_u) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_u, b=P_sector_d) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_u, b=P_sector_e) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_u, b=P_sector_nu) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_d, b=         I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_d, b=  Q_sector) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_d, b=P_sector_u) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_d, b=P_sector_d) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_d, b=P_sector_e) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_d, b=P_sector_nu) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_e, b=         I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_e, b=  Q_sector) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_e, b=P_sector_u) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_e, b=P_sector_d) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_e, b=P_sector_e) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_e, b=P_sector_nu) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_nu, b=         I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_nu, b=  Q_sector) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_nu, b=P_sector_u) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_nu, b=P_sector_d) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_nu, b=P_sector_e) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_nu, b=P_sector_nu) → ||[[D,a],J b J^-1]||_F = 0.000e+00

=== Zero-order condition test ===
Max Frobenius norm over all pairs (a,b): 0.000e+00
All pairs satisfy [a, J b J^-1] ≈ 0 within eps=1.0e-12

=== Grading & reality tests ===
||{gamma_F, D_F}||_F = 0.000e+00
max ||[gamma_F, a]||_F over a∈A_F = 0.000e+00
||S^2 - I||_F  (⇒ J_F^2 deviation) = 0.000e+00
||J D_F J^-1 - D_F||_F   = 0.000e+00
||J D_F J^-1 + D_F||_F   = 1.519e-01
→ KO-sign: J D_F J^-1 = + D_F (J-even Dirac operator)
"""

from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple
import numpy as np

# ================================
# Finite SM triple (1 generation + nu_R)
# ================================

@dataclass
class SMState:
    name: str           # "nu_L", "u_L_r", "bar_u_R_b", ...
    chirality: str      # "L" or "R"
    sector: str         # "lepton" or "quark"
    particle: bool      # True = particle, False = antiparticle
    color: Optional[str]
    generation: int


def build_sm_basis_one_gen(include_nu_R: bool = True) -> List[SMState]:
    """1 generation + nu_R + antiparticles → dim(H_F^SM) = 32."""
    basis: List[SMState] = []
    gen = 1

    # Particles: leptons
    basis.append(SMState("nu_L", "L", "lepton", True, None, gen))
    basis.append(SMState("e_L",  "L", "lepton", True, None, gen))
    if include_nu_R:
        basis.append(SMState("nu_R", "R", "lepton", True, None, gen))
    basis.append(SMState("e_R",  "R", "lepton", True, None, gen))

    # Particles: quarks
    colors = ["r", "g", "b"]
    for c in colors:
        basis.append(SMState(f"u_L_{c}", "L", "quark", True, c, gen))
        basis.append(SMState(f"d_L_{c}", "L", "quark", True, c, gen))
    for c in colors:
        basis.append(SMState(f"u_R_{c}", "R", "quark", True, c, gen))
        basis.append(SMState(f"d_R_{c}", "R", "quark", True, c, gen))

    # Antiparticles: same multiplets, opposite chirality
    particle_states = basis.copy()
    for st in particle_states:
        new_chir = "R" if st.chirality == "L" else "L"
        basis.append(SMState("bar_" + st.name, new_chir, st.sector, False, st.color, st.generation))

    return basis


def build_name_index_map_sm(basis: List[SMState]) -> Dict[str, int]:
    return {st.name: i for i, st in enumerate(basis)}


def build_gamma_sm(basis: List[SMState]) -> np.ndarray:
    """γ: -1 on left-chiral states, +1 on right-chiral states (particles + antiparticles)."""
    dim = len(basis)
    gamma = np.zeros((dim, dim), dtype=complex)
    for i, st in enumerate(basis):
        gamma[i, i] = -1.0 if st.chirality == "L" else 1.0
    return gamma


def build_swap_J_sm(basis: List[SMState]) -> np.ndarray:
    """Swap matrix S implementing particle ↔ antiparticle exchange."""
    dim = len(basis)
    S = np.zeros((dim, dim), dtype=complex)

    idx_map: Dict[tuple, int] = {
        (st.name, st.particle): i
        for i, st in enumerate(basis)
    }

    for i, st in enumerate(basis):
        if st.particle:
            partner_name = "bar_" + st.name
            j = idx_map[(partner_name, False)]
        else:
            assert st.name.startswith("bar_")
            base_name = st.name[len("bar_"):]
            j = idx_map[(base_name, True)]
        S[i, j] = 1.0

    return S

"""
# Example:
basis_SM = build_sm_basis_one_gen()
check_KO_relations_SM(basis_SM)

||S^2 - I||_F = 0.0
||JγJ^{-1} + γ||_F = 0.0
||JγJ^{-1} - γ||_F = 11.313708498984761

"""
def add_quaternion_on_doublet(M: np.ndarray,
                              q2: np.ndarray,
                              idx1: int,
                              idx2: int,
                              conj: bool = False) -> None:
    """
    Add the 2x2 matrix q2 (or its conjugate) on the subspace spanned by idx1, idx2.
    """
    block = q2.conj() if conj else q2
    indices = [idx1, idx2]
    for a, i in enumerate(indices):
        for b, j in enumerate(indices):
            M[i, j] += block[a, b]


def add_color_matrix_on_triplet(M: np.ndarray,
                                m3: np.ndarray,
                                idxs: List[int],
                                conj: bool = False) -> None:
    """
    Add the 3x3 matrix m3 (or its conjugate) on the subspace of indices idxs.
    """
    block = m3.conj() if conj else m3
    assert len(idxs) == 3
    for a, i in enumerate(idxs):
        for b, j in enumerate(idxs):
            M[i, j] += block[a, b]
def rep_A_SM(lambda_c: complex,
             q2: np.ndarray,
             m3: np.ndarray,
             basis: List[SMState],
             name_to_index: Dict[str, int]) -> np.ndarray:
    """
    Representation of a = (lambda_c, q2, m3) in A_F^SM on H_F^SM (1 generation + nu_R).

    - lambda_c: scalar (C-part), acts as lambda_c * I on all states.
    - q2: 2x2 complex matrix (H-part), acts on SU(2) doublets:
         (nu_L, e_L), (u_L_c, d_L_c), and their antiparticle doublets.
    - m3: 3x3 complex matrix (color part), acts on color triplets for quarks
          (fundamental on particles, conjugate on antiparticles).
    """
    dim = len(basis)
    M = lambda_c * np.eye(dim, dtype=complex)

    # 1) Quaternion (SU(2)) part: act on doublets

    # Lepton doublet (particles)
    idx_nu_L  = name_to_index["nu_L"]
    idx_e_L   = name_to_index["e_L"]
    add_quaternion_on_doublet(M, q2, idx_nu_L, idx_e_L, conj=False)

    # Lepton doublet (antiparticles)
    idx_bar_nu_L = name_to_index["bar_nu_L"]
    idx_bar_e_L  = name_to_index["bar_e_L"]
    add_quaternion_on_doublet(M, q2, idx_bar_nu_L, idx_bar_e_L, conj=True)

    # Quark doublets (particles and antiparticles) for each color
    colors = ["r", "g", "b"]
    for c in colors:
        # particles: (u_L_c, d_L_c)
        idx_uL = name_to_index[f"u_L_{c}"]
        idx_dL = name_to_index[f"d_L_{c}"]
        add_quaternion_on_doublet(M, q2, idx_uL, idx_dL, conj=False)

        # antiparticles: (bar_u_L_c, bar_d_L_c)
        idx_bar_uL = name_to_index[f"bar_u_L_{c}"]
        idx_bar_dL = name_to_index[f"bar_d_L_{c}"]
        add_quaternion_on_doublet(M, q2, idx_bar_uL, idx_bar_dL, conj=True)

    # Right-handed singlets (nu_R, e_R, u_R_c, d_R_c, and their antiparticles)
    # see only lambda_c and color; we do not apply q2 on them here.

    # 2) Color part: act on quark triplets in color space

    # Helper to collect indices of color triplets by name pattern
    def idx_triplet(prefix: str) -> List[int]:
        return [name_to_index[f"{prefix}_{c}"] for c in colors]

    # Particle quark triplets
    uL_triplet = idx_triplet("u_L")
    dL_triplet = idx_triplet("d_L")
    uR_triplet = idx_triplet("u_R")
    dR_triplet = idx_triplet("d_R")

    add_color_matrix_on_triplet(M, m3, uL_triplet, conj=False)
    add_color_matrix_on_triplet(M, m3, dL_triplet, conj=False)
    add_color_matrix_on_triplet(M, m3, uR_triplet, conj=False)
    add_color_matrix_on_triplet(M, m3, dR_triplet, conj=False)

    # Antiparticle quark triplets (conjugate rep)
    bar_uL_triplet = [name_to_index[f"bar_u_L_{c}"] for c in colors]
    bar_dL_triplet = [name_to_index[f"bar_d_L_{c}"] for c in colors]
    bar_uR_triplet = [name_to_index[f"bar_u_R_{c}"] for c in colors]
    bar_dR_triplet = [name_to_index[f"bar_d_R_{c}"] for c in colors]

    add_color_matrix_on_triplet(M, m3, bar_uL_triplet, conj=True)
    add_color_matrix_on_triplet(M, m3, bar_dL_triplet, conj=True)
    add_color_matrix_on_triplet(M, m3, bar_uR_triplet, conj=True)
    add_color_matrix_on_triplet(M, m3, bar_dR_triplet, conj=True)

    # Leptons have no color, so no m3 action there.

    return M
@dataclass
class SMAlgebraElement:
    name: str
    matrix: np.ndarray


def pauli_matrices() -> List[np.ndarray]:
    sigma1 = np.array([[0, 1], [1, 0]], dtype=complex)
    sigma2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sigma3 = np.array([[1, 0], [0, -1]], dtype=complex)
    return [sigma1, sigma2, sigma3]


def simple_gell_mann_subset() -> List[np.ndarray]:
    # lambda3 and lambda8 in standard Gell-Mann basis
    lam3 = np.array([[1, 0, 0],
                     [0,-1, 0],
                     [0, 0, 0]], dtype=complex)
    lam8 = (1/np.sqrt(3)) * np.array([[1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0,-2]], dtype=complex)
    return [lam3, lam8]


def build_SM_algebra_generators(basis: List[SMState],
                                name_to_index: Dict[str, int]) -> List[SMAlgebraElement]:
    dim = len(basis)

    # Identity (pure C-part with lambda=1)
    I_SM = rep_A_SM(lambda_c=1.0+0j,
                    q2=np.zeros((2, 2), dtype=complex),
                    m3=np.zeros((3, 3), dtype=complex),
                    basis=basis,
                    name_to_index=name_to_index)

    gens: List[SMAlgebraElement] = [SMAlgebraElement("I", I_SM)]

    # Quaternion part generators: lambda=0, m3=0, q2 = Pauli matrices
    for k, sigma in enumerate(pauli_matrices(), start=1):
        M_sigma = rep_A_SM(lambda_c=0.0+0j,
                           q2=sigma,
                           m3=np.zeros((3, 3), dtype=complex),
                           basis=basis,
                           name_to_index=name_to_index)
        gens.append(SMAlgebraElement(f"H_sigma{k}", M_sigma))

    # Color part generators: lambda=0, q2=0, m3 = subset of Gell-Mann matrices
    for k, gm in enumerate(simple_gell_mann_subset(), start=3):
        M_gm = rep_A_SM(lambda_c=0.0+0j,
                        q2=np.zeros((2, 2), dtype=complex),
                        m3=gm,
                        basis=basis,
                        name_to_index=name_to_index)
        gens.append(SMAlgebraElement(f"color_lambda{k}", M_gm))

    return gens
def build_SM_algebra_generators_commutative(basis: List[SMState],
                                            name_to_index: Dict[str, int]) -> List[SMAlgebraElement]:
    """
    Commutative subalgebra of A_F^SM for testing:
    - Identity
    - One diagonal quaternion generator (sigma3)
    - Two diagonal color generators (lambda3, lambda8)
    This should satisfy the zero-order condition with our current J/D shell.
    """
    gens: List[SMAlgebraElement] = []

    # Identity (pure C-part)
    I_SM = rep_A_SM(
        lambda_c=1.0+0j,
        q2=np.zeros((2, 2), dtype=complex),
        m3=np.zeros((3, 3), dtype=complex),
        basis=basis,
        name_to_index=name_to_index,
    )
    gens.append(SMAlgebraElement("I", I_SM))

    # Diagonal quaternion generator: sigma3
    sigma3 = np.array([[1, 0],
                       [0,-1]], dtype=complex)
    H_sigma3 = rep_A_SM(
        lambda_c=0.0+0j,
        q2=sigma3,
        m3=np.zeros((3, 3), dtype=complex),
        basis=basis,
        name_to_index=name_to_index,
    )
    gens.append(SMAlgebraElement("H_sigma3", H_sigma3))

    # Diagonal color generators: lambda3 and lambda8
    lam3 = np.array([[ 1, 0, 0],
                     [ 0,-1, 0],
                     [ 0, 0, 0]], dtype=complex)
    lam8 = (1/np.sqrt(3)) * np.array([[ 1, 0, 0],
                                      [ 0, 1, 0],
                                      [ 0, 0,-2]], dtype=complex)

    color_lam3 = rep_A_SM(
        lambda_c=0.0+0j,
        q2=np.zeros((2, 2), dtype=complex),
        m3=lam3,
        basis=basis,
        name_to_index=name_to_index,
    )
    gens.append(SMAlgebraElement("color_lambda3", color_lam3))

    color_lam8 = rep_A_SM(
        lambda_c=0.0+0j,
        q2=np.zeros((2, 2), dtype=complex),
        m3=lam8,
        basis=basis,
        name_to_index=name_to_index,
    )
    gens.append(SMAlgebraElement("color_lambda8", color_lam8))

    return gens

# ================================
# Internal Hilbert space & D_F (NCG-compatible toy)
# ================================

SECTORS: List[str] = ["u", "d", "e", "nu"]
SECTOR_INDEX: Dict[str, int] = {s: i for i, s in enumerate(SECTORS)}
N_GEN: int = 3

# We treat color as a degeneracy factor on u,d (3 copies) and 1 on e,nu.
# It is *not* yet a full SU(3) algebra action; that would require an explicit
# color tensor factor and a more refined J_F. Here we only count dimensions.
SECTOR_NC: Dict[str, int] = {"u": 3, "d": 3, "e": 1, "nu": 1}

def base_kernel(lam, alpha=3.0, form="lambda_sq"):
    """
    Base kernel F_base(λ_g) that defines the generation ladder.

    We make it *scale-invariant* by normalizing the eigenvalues to the
    lightest nonzero one, so that a global rescaling of the Laplacian
    does not flatten or blow up the hierarchy:

        F_base(λ_g) = exp[-alpha * (λ_g / λ_ref)^2]

    with λ_ref = smallest positive eigenvalue in the triad.
    """
    lam = np.array(lam, dtype=float)

    # Choose a reference eigenvalue λ_ref (smallest positive λ)
    lam_pos = lam[lam > 0]
    if lam_pos.size == 0:
        # Degenerate case: fall back to ordinary λ^2 kernel
        lam_ref = 1.0
    else:
        lam_ref = lam_pos.min()

    x = lam / lam_ref

    if form == "lambda_sq":
        return np.exp(-alpha * x**2)
    elif form == "lambda":
        return np.exp(-alpha * x)
    else:
        raise ValueError(f"Unknown kernel form '{form}'")

def dim_per_chirality() -> int:
    """Dimension of H_L or H_R (one chirality).

    We fold color multiplicities into sector blocks:
      u,d: 3 each; e,nu: 1 each → total 8 per generation
      times 3 generations → 24 per chirality.
    """
    return 3 * sum(SECTOR_NC[s] for s in SECTORS)  # 24


def flavor_block_offsets() -> Dict[str, int]:
    """Offsets (within a single chirality) for each sector's 3×3
    generation block in a 12×12 generation-space layout:

      [u_g1,u_g2,u_g3,
       d_g1,d_g2,d_g3,
       e_g1,e_g2,e_g3,
       nu_g1,nu_g2,nu_g3]

    We only care about generation offsets (3×3 blocks);
    color multiplicity is treated as degeneracy, not an explicit tensor factor.
    """
    off: Dict[str, int] = {}
    off["u"]  = 0
    off["d"]  = 3
    off["e"]  = 6
    off["nu"] = 9
    return off


def build_internal_DF_from_Y(Y_u: np.ndarray,
                             Y_d: np.ndarray,
                             Y_e: np.ndarray,
                             Y_nu: np.ndarray) -> np.ndarray:
    """Build the finite Dirac operator D_F in block form:

        D_F = [[ 0, Y^\dagger ],
               [ Y, 0         ]]

    where Y is a 24×24 block that is block-diagonal in sector space and
    embeds the 3×3 generation Yukawas (Y_u, Y_d, Y_e, Y_nu) into a
    12×12 generation-space layout, with color treated as degeneracy.

    H_F ≃ H_L ⊕ H_R,  dim(H_L) = dim(H_R) = 24, dim(H_F) = 48.
    """
    # Sanity checks
    for name, Y in [("Y_u", Y_u), ("Y_d", Y_d), ("Y_e", Y_e), ("Y_nu", Y_nu)]:
        Y = np.asarray(Y, dtype=complex)
        if Y.shape != (3, 3):
            raise ValueError(f"{name} must be a 3×3 matrix, got shape {Y.shape}.")
    Y_u  = np.asarray(Y_u, dtype=complex)
    Y_d  = np.asarray(Y_d, dtype=complex)
    Y_e  = np.asarray(Y_e, dtype=complex)
    Y_nu = np.asarray(Y_nu, dtype=complex)

    dpc = dim_per_chirality()  # 24
    dimH = 2 * dpc             # 48

    # 12×12 generation-space Yukawa core
    Y_gen = np.zeros((12, 12), dtype=complex)
    gen_off = flavor_block_offsets()

    def insert_sector_Y(sector: str, Y_s: np.ndarray) -> None:
        off = gen_off[sector]
        Y_gen[off:off+3, off:off+3] = Y_s

    insert_sector_Y("u",  Y_u)
    insert_sector_Y("d",  Y_d)
    insert_sector_Y("e",  Y_e)
    insert_sector_Y("nu", Y_nu)

    # Embed the 12×12 generation block into 24×24 per chirality.
    # Only the leading 12×12 carry Yukawa couplings; the remaining slots
    # are color-degenerate but Yukawa-silent in this toy.
    Y_block = np.zeros((dpc, dpc), dtype=complex)
    Y_block[:12, :12] = Y_gen

    # Assemble D_F on H_F = H_L ⊕ H_R
    D_F = np.zeros((dimH, dimH), dtype=complex)
    D_F[:dpc,  dpc:] = Y_block.conj().T
    D_F[dpc:,  :dpc] = Y_block

    return D_F


# ================================
# Real structure, grading, algebra basis
# ================================

def build_swap_LR(dim_left: int) -> np.ndarray:
    """Swap matrix S on H_F = H_L ⊕ H_R, with dim(H_L) = dim(H_R) = dim_left.

    Acts as: S (ψ_L, ψ_R) = (ψ_R, ψ_L).
    """
    S = np.zeros((2 * dim_left, 2 * dim_left), dtype=complex)
    S[:dim_left, dim_left:] = np.eye(dim_left)
    S[dim_left:, :dim_left] = np.eye(dim_left)
    return S


def build_gamma_F(dim_left: int) -> np.ndarray:
    """Grading operator γ_F with eigenvalue -1 on H_L and +1 on H_R."""
    g = np.zeros((2 * dim_left, 2 * dim_left), dtype=complex)
    g[:dim_left, :dim_left] = -np.eye(dim_left)
    g[dim_left:, dim_left:] =  np.eye(dim_left)
    return g


def build_sector_projectors() -> Dict[str, np.ndarray]:
    """Sector projectors P_sector_s acting on H_F = H_L ⊕ H_R.

    Each P_sector_s is diagonal and selects the (sector,gen,chirality) subspace
    corresponding to that sector (u,d,e,nu) in the 12×12 generation subspace,
    duplicated on L and R.
    """
    dpc = dim_per_chirality()
    dimH = 2 * dpc
    gen_off = flavor_block_offsets()

    P: Dict[str, np.ndarray] = {}
    for s in SECTORS:
        P_s = np.zeros((dimH, dimH), dtype=complex)
        off = gen_off[s]
        # Same on L and R, only on first 12 generation slots
        P_s[off:off+3, off:off+3] = np.eye(3)
        P_s[dpc + off:dpc + off + 3, dpc + off:dpc + off + 3] = np.eye(3)
        P[s] = P_s

    return P


def build_Q_sector() -> np.ndarray:
    """A simple 'sector charge' diagonal operator Q_sector.

    Distinguishes u,d,e,nu sectors but is generation-blind:
      q_u = 2, q_d = 1, q_e = 0, q_nu = -1  (toy choice).
    """
    dpc = dim_per_chirality()
    dimH = 2 * dpc
    gen_off = flavor_block_offsets()
    charges = {"u": 2.0, "d": 1.0, "e": 0.0, "nu": -1.0}

    Q = np.zeros((dimH, dimH), dtype=complex)
    for s in SECTORS:
        off = gen_off[s]
        q = charges[s]
        Q[off:off+3, off:off+3] = q * np.eye(3)
        Q[dpc + off:dpc + off + 3, dpc + off:dpc + off + 3] = q * np.eye(3)

    return Q


def build_internal_algebra_ops() -> Tuple[List[np.ndarray], List[str]]:
    """Small basis of algebra elements A_F acting on H_F:

        - I (identity)
        - Q_sector (diagonal sector 'charge')
        - P_sector_u, P_sector_d, P_sector_e, P_sector_nu (sector projectors)

    This is a commutative algebra in this toy (no explicit SU(3) yet).
    """
    dpc = dim_per_chirality()
    dimH = 2 * dpc

    I = np.eye(dimH, dtype=complex)
    Q = build_Q_sector()
    P = build_sector_projectors()

    ops: List[np.ndarray] = [I, Q, P["u"], P["d"], P["e"], P["nu"]]
    labels: List[str] = ["I", "Q_sector", "P_sector_u", "P_sector_d", "P_sector_e", "P_sector_nu"]

    return ops, labels


# ================================
# NCG condition tests
# ================================

def J_action_from_swap(S: np.ndarray, M: np.ndarray) -> np.ndarray:
    """Implement J M J^{-1} = S · M^* · S^T, where S is the L/R swap."""
    return S @ M.conj() @ S.T


def test_first_order_condition(D_F: np.ndarray,
                               ops: List[np.ndarray],
                               labels: List[str],
                               eps: float = 1e-12) -> None:
    """First-order condition:

        [[D_F, a], J_F b J_F^{-1}] = 0   for all a,b in A_F.
    """
    n = D_F.shape[0]
    assert D_F.shape == (n, n)
    S = build_swap_LR(dim_left=n // 2)

    print("=== First-order condition test ===")
    max_norm = 0.0
    good_pairs = []

    for i, a in enumerate(ops):
        Da = D_F @ a - a @ D_F
        for j, b in enumerate(ops):
            b_tilde = J_action_from_swap(S, b)
            comm2 = Da @ b_tilde - b_tilde @ Da
            norm = np.linalg.norm(comm2, ord="fro")
            if norm > max_norm:
                max_norm = norm
            if norm < eps:
                good_pairs.append((labels[i], labels[j], norm))

    print(f"Max Frobenius norm over all pairs (a,b): {max_norm:.3e}")
    if good_pairs:
        print(f"Pairs with norm < {eps:.1e}:")
        for la, lb, nrm in good_pairs:
            print(f"  (a={la:>10s}, b={lb:>10s}) → ||[[D,a],J b J^-1]||_F = {nrm:.3e}")
    else:
        print(f"No pairs with norm < {eps:.1e}")
    print()


def test_zero_order_condition(ops: List[np.ndarray],
                              labels: List[str],
                              eps: float = 1e-12) -> None:
    """Zero-order condition:

        [a, J_F b J_F^{-1}] = 0   for all a,b in A_F.
    """
    n = ops[0].shape[0]
    S = build_swap_LR(dim_left=n // 2)

    print("=== Zero-order condition test ===")
    max_norm = 0.0
    bad_pairs = []

    for i, a in enumerate(ops):
        for j, b in enumerate(ops):
            b_tilde = J_action_from_swap(S, b)
            comm = a @ b_tilde - b_tilde @ a
            norm = np.linalg.norm(comm, ord="fro")
            if norm > max_norm:
                max_norm = norm
            if norm > eps:
                bad_pairs.append((labels[i], labels[j], norm))

    print(f"Max Frobenius norm over all pairs (a,b): {max_norm:.3e}")
    if bad_pairs:
        print("Pairs with significant violation:")
        for la, lb, nrm in bad_pairs:
            print(f"  (a={la:>10s}, b={lb:>10s}) → ||[a, J b J^-1]||_F = {nrm:.3e}")
    else:
        print(f"All pairs satisfy [a, J b J^-1] ≈ 0 within eps={eps:.1e}")
    print()


def test_grading_and_reality(D_F: np.ndarray,
                             ops: List[np.ndarray],
                             labels: List[str]) -> None:
    """Check grading and reality axioms:

      - γ_F anticommutes with D_F and commutes with A_F.
      - J_F^2 = 1 (as implemented by swap).
      - KO-dimension sign via J D_F J^{-1} = ± D_F.
    """
    n = D_F.shape[0]
    dpc = n // 2
    gamma_F = build_gamma_F(dpc)
    S = build_swap_LR(dpc)

    print("=== Grading & reality tests ===")
    anti = gamma_F @ D_F + D_F @ gamma_F
    print(f"||{{gamma_F, D_F}}||_F = {np.linalg.norm(anti, ord='fro'):.3e}")

    max_comm_gamma = 0.0
    for a in ops:
        comm_ga = gamma_F @ a - a @ gamma_F
        max_comm_gamma = max(max_comm_gamma, np.linalg.norm(comm_ga, ord="fro"))
    print(f"max ||[gamma_F, a]||_F over a∈A_F = {max_comm_gamma:.3e}")

    S2 = S @ S
    print(f"||S^2 - I||_F  (⇒ J_F^2 deviation) = {np.linalg.norm(S2 - np.eye(n), ord='fro'):.3e}")

    JDJ = S @ D_F.conj() @ S.T
    norm_minus = np.linalg.norm(JDJ - D_F, ord="fro")
    norm_plus  = np.linalg.norm(JDJ + D_F, ord="fro")
    print(f"||J D_F J^-1 - D_F||_F   = {norm_minus:.3e}")
    print(f"||J D_F J^-1 + D_F||_F   = {norm_plus:.3e}")
    if norm_minus < norm_plus:
        print("→ KO-sign: J D_F J^-1 = + D_F (J-even Dirac operator)")
    else:
        print("→ KO-sign: J D_F J^-1 = - D_F (J-odd Dirac operator)")
    print()


# ================================
# Emergent misalignment model, flavor, mixing, chi^2
# (your original χ^2≈11 toy, kept intact below)
# ================================

# --- everything below here is your original emergent-4-x11 code ---
# (misalignment functional, emergent graph, Laplacian, F_base, Q,
#  geometry-derived U_geom, Yukawas, mixing, chi^2, etc.)

# I’m not re-commenting every function here since they’re unchanged;
# this is literally your stable χ²≈11 script with the NCG block added above
# and the NCG tests called at the end of main().

# -------------- misalignment functional, relaxation, graph, etc. --------------

def misalignment_energy(theta, w6=1.0, w5=1.0):
    N = len(theta)
    diffs = theta[:, None] - theta[None, :]
    cos6 = np.cos(6 * diffs)
    cos5 = np.cos(5 * diffs)
    E6 = w6 * np.sum(1.0 - cos6) / (N * N)
    E5 = w5 * np.sum(1.0 - cos5) / (N * N)
    return E6 + E5


def relax_phases(N=200, n_steps=600, eta=0.01, w6=1.0, w5=1.0, random_seed=42):
    rng = np.random.default_rng(random_seed)
    theta = rng.uniform(0, 2 * np.pi, size=N)
    energy_hist = []

    for step in range(n_steps):
        diffs = theta[:, None] - theta[None, :]
        sin6 = np.sin(6 * diffs)
        sin5 = np.sin(5 * diffs)
        grad = 6 * w6 * np.sum(sin6, axis=1) + 5 * w5 * np.sum(sin5, axis=1)
        theta = theta - eta * grad
        theta = (theta + 2 * np.pi) % (2 * np.pi)

        if step % 10 == 0 or step == n_steps - 1:
            E = misalignment_energy(theta, w6=w6, w5=w5)
            energy_hist.append(E)

    return theta, energy_hist


def build_emergent_adjacency(theta, w6=1.0, w5=1.0, keep_fraction=0.05):
    N = len(theta)
    diffs = theta[:, None] - theta[None, :]
    cos6 = np.cos(6 * diffs)
    cos5 = np.cos(5 * diffs)
    score = w6 * cos6 + w5 * cos5
    np.fill_diagonal(score, -np.inf)
    triu_idx = np.triu_indices(N, k=1)
    flat_scores = score[triu_idx]
    k = int(keep_fraction * len(flat_scores))
    if k < 1:
        k = 1
    kth_val = np.partition(flat_scores, -k)[-k]
    A = np.zeros((N, N), dtype=float)
    mask = (score >= kth_val)
    A[mask] = 1.0
    A = np.maximum(A, A.T)
    return A


def largest_connected_component(A):
    N = A.shape[0]
    visited = np.zeros(N, dtype=bool)
    best_comp = []
    for i in range(N):
        if not visited[i]:
            stack = [i]
            comp = []
            visited[i] = True
            while stack:
                v = stack.pop()
                comp.append(v)
                neighbors = np.where(A[v] > 0)[0]
                for w in neighbors:
                    if not visited[w]:
                        visited[w] = True
                        stack.append(w)
            if len(comp) > len(best_comp):
                best_comp = comp
    best_comp = np.array(best_comp, dtype=int)
    A_sub = A[np.ix_(best_comp, best_comp)]
    return A_sub, best_comp


def laplacian_from_adjacency(A):
    d = np.sum(A, axis=1)
    L = np.diag(d) - A
    return L


def spectral_triad(L):
    eigvals, eigvecs = np.linalg.eigh(L)
    idx_sorted = np.argsort(eigvals)
    eigvals_sorted = eigvals[idx_sorted]
    eigvecs_sorted = eigvecs[:, idx_sorted]
    lam_gen = eigvals_sorted[1:4]
    gen_indices = idx_sorted[1:4]
    return lam_gen, gen_indices, eigvals_sorted


# -------------- sector charges and F_s -----------------

def build_sector_charges():
    Q_u = np.array([0,  2,  4], dtype=float)
    Q_d = np.array([1,  3,  5], dtype=float)
    Q_e = np.array([2,  4,  6], dtype=float)
    Q_n = np.array([4,  6,  8], dtype=float)
    return {"u": Q_u, "d": Q_d, "e": Q_e, "nu": Q_n}


def sector_weights(F_base, Q_s, beta=1.0):
    return F_base * np.exp(-beta * Q_s)


def mass_ratios(F_s):
    F_s = np.array(F_s, dtype=float)
    m1, m2, m3 = F_s
    return m1 / m3, m2 / m3


# -------------- generation operators (golden, Cabibbo) --------------

def rotation_3d(i, j, theta):
    R = np.eye(3, dtype=complex)
    c = np.cos(theta)
    s = np.sin(theta)
    R[i, i] = c
    R[j, j] = c
    R[i, j] = s
    R[j, i] = -s
    return R


def build_generation_operators(phi_order=5, cab_denom=28):
    theta_phi = 2 * np.pi / phi_order
    theta_C = 2 * np.pi / cab_denom
    P_phi_12 = rotation_3d(0, 1, theta_phi)
    P_phi_23 = rotation_3d(1, 2, theta_phi)
    C_12 = rotation_3d(0, 1, theta_C)
    return P_phi_12, P_phi_23, C_12, theta_phi, theta_C


# -------------- geometric regions and unitaries --------------

def build_geometric_regions(theta, n_regions=3):
    phase = np.mod(theta, 2 * np.pi)
    edges = np.linspace(0, 2*np.pi, n_regions+1)
    regions = []
    for k in range(n_regions):
        lo, hi = edges[k], edges[k+1]
        if k < n_regions - 1:
            idx = np.where((phase >= lo) & (phase < hi))[0]
        else:
            idx = np.where((phase >= lo) & (phase <= hi))[0]
        if len(idx) == 0:
            idx = np.array([k % len(theta)], dtype=int)
        regions.append(idx)
    return regions


def build_geometric_unitary(gen_vecs, region_list):
    cols = []
    for R in region_list:
        v = np.sum(gen_vecs[R, :], axis=0)
        norm = np.linalg.norm(v)
        if norm < 1e-14:
            v = np.array([1.0, 0.0, 0.0], dtype=complex)
            norm = 1.0
        cols.append(v / norm)
    U_geom = np.column_stack(cols)
    Uu, _, Vh = np.linalg.svd(U_geom)
    return Uu @ Vh


def build_sector_bases(P_phi_12, P_phi_23, C_12, U_geom, use_neutrino_dressing=True,
                       N_SOLAR=36, N_REACTOR=45):
    sector_bases = {}

    U_geom_u = U_geom["u"]
    U_geom_d = U_geom["d"]
    U_geom_e = U_geom["e"]
    U_geom_nu = U_geom["nu"]

    U_L_u = U_geom_u @ C_12.conj().T
    U_R_u = np.eye(3, dtype=complex)

    U_L_d = U_geom_d
    U_R_d = np.eye(3, dtype=complex)

    U_L_e = U_geom_e
    U_R_e = np.eye(3, dtype=complex)

    if use_neutrino_dressing:
        theta_solar = 2 * np.pi / N_SOLAR
        theta_reac = 2 * np.pi / N_REACTOR
        R_solar = rotation_3d(0, 1, theta_solar)
        R_reac = rotation_3d(0, 2, theta_reac)
        U_dress = (P_phi_23 @ R_solar) @ (P_phi_12 @ R_reac)
        U_L_nu = U_geom_nu @ U_dress
    else:
        U_L_nu = U_geom_nu

    U_R_nu = np.eye(3, dtype=complex)

    sector_bases["u"] = (U_L_u, U_R_u)
    sector_bases["d"] = (U_L_d, U_R_d)
    sector_bases["e"] = (U_L_e, U_R_e)
    sector_bases["nu"] = (U_L_nu, U_R_nu)

    return sector_bases


# -------------- Yukawas, mixing, observables, chi^2 --------------

def yukawa_from_F_and_UL(F_s, U_L, U_R):
    D = np.diag(F_s)
    return U_L @ D @ U_R.conj().T


def mixing_matrix(U_L_up, U_L_down):
    return U_L_up.conj().T @ U_L_down


def mixing_angles_from_U(U):
    s13 = abs(U[0, 2])
    s13 = min(max(s13, 0.0), 1.0)
    theta13 = np.arcsin(s13)
    c13 = np.cos(theta13)
    if abs(c13) < 1e-12:
        theta12 = 0.0
        theta23 = 0.0
    else:
        theta12 = np.arctan2(abs(U[0, 1]), abs(U[0, 0]))
        theta23 = np.arctan2(abs(U[1, 2]), abs(U[2, 2]))
    return theta12, theta23, theta13


TARGETS = {
    "mu_mt":     (2.2e-05, 0.5 * 2.2e-05),
    "mc_mt":     (7.5e-03, 0.5 * 7.5e-03),
    "md_mb":     (1.1e-03, 0.5 * 1.1e-03),
    "ms_mb":     (2.2e-02, 0.5 * 2.2e-02),
    "me_mt":     (2.9e-04, 0.5 * 2.9e-04),
    "mmu_mt":    (5.9e-02, 0.5 * 5.9e-02),
    "theta12_q": (0.227,   0.05 * 0.227),
    "theta23_q": (0.041,   0.5  * 0.041),
    "theta13_q": (0.0036,  0.5  * 0.0036),
    "theta12_l": (0.584,   0.1  * 0.584),
    "theta23_l": (0.785,   0.2  * 0.785),
    "theta13_l": (0.15,    0.2  * 0.15),
}


def compute_observables(mu_mt, mc_mt, md_mb, ms_mb, me_mt, mmu_mt,
                        theta12_q, theta23_q, theta13_q,
                        theta12_l, theta23_l, theta13_l):
    return {
        "mu_mt":     mu_mt,
        "mc_mt":     mc_mt,
        "md_mb":     md_mb,
        "ms_mb":     ms_mb,
        "me_mt":     me_mt,
        "mmu_mt":    mmu_mt,
        "theta12_q": theta12_q,
        "theta23_q": theta23_q,
        "theta13_q": theta13_q,
        "theta12_l": theta12_l,
        "theta23_l": theta23_l,
        "theta13_l": theta13_l,
    }


def chi2(obs, targets):
    chi2_val = 0.0
    details = []
    for k, v in obs.items():
        target, sigma = targets[k]
        if sigma <= 0:
            continue
        contrib = ((v - target) / sigma)**2
        chi2_val += contrib
        details.append((k, v, target, contrib))
    return chi2_val, details

def J_action_from_S(S: np.ndarray, M: np.ndarray) -> np.ndarray:
    """Adjoint action J M J^{-1} = S · M^* · S^T."""
    return S @ M.conj() @ S.T

def build_DF_SM_one_gen(basis: List[SMState],
                        name_to_index: Dict[str, int],
                        m_nu: float = 0.1,
                        m_e: float  = 0.5,
                        m_u: float  = 2.0,
                        m_d: float  = 4.0) -> np.ndarray:
    """
    Simple finite Dirac operator for 1 generation:
    - Dirac masses linking L ↔ R for (nu, e, u_c, d_c) and their antiparticles.
    - No Majorana block yet.
    """
    dim = len(basis)
    D = np.zeros((dim, dim), dtype=complex)

    def couple(a: str, b: str, mass: float) -> None:
        i = name_to_index[a]
        j = name_to_index[b]
        D[i, j] = mass
        D[j, i] = mass

    # Leptons
    if "nu_R" in name_to_index:
        couple("nu_L", "nu_R", m_nu)
        couple("bar_nu_L", "bar_nu_R", m_nu)
    couple("e_L", "e_R", m_e)
    couple("bar_e_L", "bar_e_R", m_e)

    # Quarks (each color separately)
    colors = ["r", "g", "b"]
    for c in colors:
        couple(f"u_L_{c}", f"u_R_{c}", m_u)
        couple(f"d_L_{c}", f"d_R_{c}", m_d)
        couple(f"bar_u_L_{c}", f"bar_u_R_{c}", m_u)
        couple(f"bar_d_L_{c}", f"bar_d_R_{c}", m_d)

    return D

def test_first_order_condition_generic(D_F: np.ndarray,
                                       ops: List[np.ndarray],
                                       labels: List[str],
                                       S: np.ndarray,
                                       eps: float = 1e-12) -> None:
    """First-order condition:

       [[D_F, a], J b J^{-1}] = 0  for all a,b.
    """
    n = D_F.shape[0]
    assert D_F.shape == (n, n)

    print("=== First-order condition test (generic J) ===")
    max_norm = 0.0
    good_pairs = []

    for i, a in enumerate(ops):
        Da = D_F @ a - a @ D_F
        for j, b in enumerate(ops):
            b_tilde = J_action_from_S(S, b)
            comm2 = Da @ b_tilde - b_tilde @ Da
            norm = np.linalg.norm(comm2, ord="fro")
            if norm > max_norm:
                max_norm = norm
            if norm < eps:
                good_pairs.append((labels[i], labels[j], norm))

    print(f"Max Frobenius norm over all pairs (a,b): {max_norm:.3e}")
    if good_pairs:
        print(f"Pairs with norm < {eps:.1e}:")
        for la, lb, nrm in good_pairs:
            print(f"  (a={la:>12s}, b={lb:>12s}) → ||[[D,a],J b J^-1]||_F = {nrm:.3e}")
    else:
        print(f"No pairs with norm < {eps:.1e}")
    print()


def test_zero_order_condition_generic(ops: List[np.ndarray],
                                      labels: List[str],
                                      S: np.ndarray,
                                      eps: float = 1e-12) -> None:
    """Zero-order condition:

       [a, J b J^{-1}] = 0  for all a,b.
    """
    n = ops[0].shape[0]
    print("=== Zero-order condition test (generic J) ===")
    max_norm = 0.0
    bad_pairs = []

    for i, a in enumerate(ops):
        for j, b in enumerate(ops):
            b_tilde = J_action_from_S(S, b)
            comm = a @ b_tilde - b_tilde @ a
            norm = np.linalg.norm(comm, ord="fro")
            if norm > max_norm:
                max_norm = norm
            if norm > eps:
                bad_pairs.append((labels[i], labels[j], norm))

    print(f"Max Frobenius norm over all pairs (a,b): {max_norm:.3e}")
    if bad_pairs:
        print("Pairs with significant violation:")
        for la, lb, nrm in bad_pairs:
            print(f"  (a={la:>12s}, b={lb:>12s}) → ||[a, J b J^-1]||_F = {nrm:.3e}")
    else:
        print(f"All pairs satisfy [a, J b J^-1] ≈ 0 within eps={eps:.1e}")
    print()


def test_grading_and_reality_generic(D_F: np.ndarray,
                                     ops: List[np.ndarray],
                                     labels: List[str],
                                     gamma: np.ndarray,
                                     S: np.ndarray) -> None:
    """Grading & reality bundle:

      - γ anticommutes with D_F.
      - γ commutes with all a in A_F.
      - J^2 = 1.
      - KO-sign via J D_F J^{-1} = ± D_F.
    """
    n = D_F.shape[0]

    print("=== Grading & reality tests (generic γ,J) ===")
    anti = gamma @ D_F + D_F @ gamma
    print(f"||{{gamma, D_F}}||_F = {np.linalg.norm(anti, ord='fro'):.3e}")

    max_comm_gamma = 0.0
    for a in ops:
        comm_ga = gamma @ a - a @ gamma
        max_comm_gamma = max(max_comm_gamma, np.linalg.norm(comm_ga, ord="fro"))
    print(f"max ||[gamma, a]||_F over a∈A_F = {max_comm_gamma:.3e}")

    S2 = S @ S
    print(f"||S^2 - I||_F  (⇒ J^2 deviation) = {np.linalg.norm(S2 - np.eye(n), ord='fro'):.3e}")

    JDJ = S @ D_F.conj() @ S.T
    norm_minus = np.linalg.norm(JDJ - D_F, ord="fro")
    norm_plus  = np.linalg.norm(JDJ + D_F, ord="fro")
    print(f"||J D_F J^-1 - D_F||_F   = {norm_minus:.3e}")
    print(f"||J D_F J^-1 + D_F||_F   = {norm_plus:.3e}")
    if norm_minus < norm_plus:
        print("→ KO-sign: J D_F J^-1 = + D_F (J-even Dirac operator)")
    else:
        print("→ KO-sign: J D_F J^-1 = - D_F (J-odd Dirac operator)")
    print()
def run_emergent_alignment(
    N: int = 200,
    n_steps: int = 600,
    eta: float = 0.01,
    w6: float = 1.0,
    w5: float = 1.0,
    keep_fraction: float = 0.05,
    alpha: float = 3.0,
    beta: float = 1.0,
    random_seed: int = 42,
    use_neutrino_dressing: bool = True,
) -> Dict[str, np.ndarray]:
    """
    Full emergent pipeline (minimal, self-contained version):

    - Relax phases on the circle via misalignment functional.
    - Build emergent graph adjacency and Laplacian.
    - Extract the 3 lowest nonzero eigenmodes (generation triad).
    - Build a base kernel F_base(λ_g) and sector-weighted F_s.
    - Build a single geometric unitary U_geom from graph regions.
    - Use the sector_bases + F_s to construct full 3×3 Yukawa matrices.

    Returns:
      {
        "Y_u": 3x3 complex Yukawa matrix,
        "Y_d": 3x3 complex Yukawa matrix,
        "Y_e": 3x3 complex Yukawa matrix,
        "Y_nu": 3x3 complex Yukawa matrix,
        "F": {"u": F_u, "d": F_d, "e": F_e, "nu": F_nu},
        "lam_gen": lam_gen  (3 eigenvalues),
      }
    """
    # 1) Relax phases on the circle
    theta, energy_hist = relax_phases(
        N=N,
        n_steps=n_steps,
        eta=eta,
        w6=w6,
        w5=w5,
        random_seed=random_seed,
    )

    # 2) Build emergent adjacency and largest connected component
    A_full = build_emergent_adjacency(theta, w6=w6, w5=w5, keep_fraction=keep_fraction)
    A_cc, comp = largest_connected_component(A_full)

    # 3) Laplacian and spectral triad
    L = laplacian_from_adjacency(A_cc)
    lam_gen, gen_indices, eigvals_sorted = spectral_triad(L)

    # 4) Base kernel and sector weights (3×3 F_s)
    F_base = base_kernel(lam_gen, alpha=alpha, form="lambda_sq")
    Q_sector = build_sector_charges()
    F_u = sector_weights(F_base, Q_sector["u"], beta=beta)
    F_d = sector_weights(F_base, Q_sector["d"], beta=beta)
    F_e = sector_weights(F_base, Q_sector["e"], beta=beta)
    F_nu = sector_weights(F_base, Q_sector["nu"], beta=beta)

    # 5) Generation eigenvectors (3 lowest nonzero)
    eigvals, eigvecs = np.linalg.eigh(L)
    idx_sorted = np.argsort(eigvals)
    eigvecs_sorted = eigvecs[:, idx_sorted]
    # columns 1,2,3 correspond to the nonzero triad used above
    gen_vecs = eigvecs_sorted[:, 1:4]  # shape (N_cc, 3)

    # 6) Geometric regions on the connected component
    theta_cc = theta[comp]
    regions = build_geometric_regions(theta_cc, n_regions=3)

    # 7) Single geometric unitary from regions
    U_geom_single = build_geometric_unitary(gen_vecs, regions)

    # Use the same U_geom for all sectors (sector-specific dressing happens later)
    U_geom = {s: U_geom_single for s in SECTORS}

    # 8) Build sector bases (U_L, U_R) per sector
    P_phi_12, P_phi_23, C_12, theta_phi, theta_C = build_generation_operators()
    sector_bases = build_sector_bases(
        P_phi_12,
        P_phi_23,
        C_12,
        U_geom,
        use_neutrino_dressing=use_neutrino_dressing,
    )

    # 9) Construct full 3×3 Yukawas from F_s and U_L/R
    Y_u = yukawa_from_F_and_UL(F_u, *sector_bases["u"])
    Y_d = yukawa_from_F_and_UL(F_d, *sector_bases["d"])
    Y_e = yukawa_from_F_and_UL(F_e, *sector_bases["e"])
    Y_nu = yukawa_from_F_and_UL(F_nu, *sector_bases["nu"])

    return {
        "Y_u": Y_u,
        "Y_d": Y_d,
        "Y_e": Y_e,
        "Y_nu": Y_nu,
        "F": {"u": F_u, "d": F_d, "e": F_e, "nu": F_nu},
        "lam_gen": lam_gen,
    }
def effective_1gen_yukawas_from_3x3(
    Y_u: np.ndarray,
    Y_d: np.ndarray,
    Y_e: np.ndarray,
    Y_nu: np.ndarray,
) -> Tuple[float, float, float, float]:
    """
    Collapse 3×3 Yukawas to a single effective 1-gen value per sector.

    Here we take the largest singular value in each sector, which
    corresponds roughly to the heaviest family (t, b, τ, ν_3).
    """
    def eff(Y: np.ndarray) -> float:
        svals = np.linalg.svd(Y, compute_uv=False)
        return float(np.max(np.abs(svals)))

    Y_u_1 = eff(Y_u)
    Y_d_1 = eff(Y_d)
    Y_e_1 = eff(Y_e)
    Y_nu_1 = eff(Y_nu)
    return Y_e_1, Y_nu_1, Y_u_1, Y_d_1
def run_ncg_with_alignment() -> None:
    """
    1) Run emergent alignment pipeline and get 3×3 Yukawas.
    2) Build and test:
       - 3-generation toy internal triple (24⊕24) with those Yukawas.
       - 1-generation SM-like triple using effective aligned Yukawas.
    """
    print("=== Running emergent alignment pipeline to obtain Yukawas ===")
    align = run_emergent_alignment()
    Y_u_align = align["Y_u"]
    Y_d_align = align["Y_d"]
    Y_e_align = align["Y_e"]
    Y_nu_align = align["Y_nu"]

    # ---------------------------------------------
    # A. 3-generation internal toy triple (24⊕24)
    # ---------------------------------------------
    print()
    print("=== NCG tests: 3-gen internal toy triple with aligned Yukawas ===")
    D_F_internal = build_internal_DF_from_Y(
        Y_u_align,
        Y_d_align,
        Y_e_align,
        Y_nu_align,
    )
    ops_internal, labels_internal = build_internal_algebra_ops()
    test_first_order_condition(D_F_internal, ops_internal, labels_internal, eps=1e-12)
    test_zero_order_condition(ops_internal, labels_internal, eps=1e-12)
    test_grading_and_reality(D_F_internal, ops_internal, labels_internal)

    # ---------------------------------------------
    # B. 1-generation SM-like triple (effective Yukawas)
    # ---------------------------------------------
    print()
    print("=== NCG tests: SM-like 1gen triple with effective aligned Yukawas ===")

    Y_e_eff, Y_nu_eff, Y_u_eff, Y_d_eff = effective_1gen_yukawas_from_3x3(
        Y_u_align,
        Y_d_align,
        Y_e_align,
        Y_nu_align,
    )

    basis_SM = build_sm_basis_one_gen(include_nu_R=True)
    name_to_index_SM = build_name_index_map_sm(basis_SM)
    gamma_SM = build_gamma_sm(basis_SM)
    S_SM = build_swap_J_sm(basis_SM)

    algebra_SM = build_SM_algebra_generators_commutative(basis_SM, name_to_index_SM)
    ops_SM = [elem.matrix for elem in algebra_SM]
    labels_SM = [elem.name for elem in algebra_SM]

    D_SM_align = build_DF_SM_one_gen(
        basis_SM,
        name_to_index_SM,
        m_nu=Y_nu_eff,
        m_e=Y_e_eff,
        m_u=Y_u_eff,
        m_d=Y_d_eff,
    )

    test_first_order_condition_generic(D_SM_align, ops_SM, labels_SM, S_SM, eps=1e-12)
    test_zero_order_condition_generic(ops_SM, labels_SM, S_SM, eps=1e-12)
    test_grading_and_reality_generic(D_SM_align, ops_SM, labels_SM, gamma_SM, S_SM)

# ================================
# Main driver
# ================================

def run_sm_ncg_tests() -> None:
    """Build 1-gen SM-like finite triple and run first-, zero-order and reality tests."""
    basis_SM = build_sm_basis_one_gen(include_nu_R=True)
    name_to_index_SM = build_name_index_map_sm(basis_SM)
    gamma_SM = build_gamma_sm(basis_SM)
    S_SM = build_swap_J_sm(basis_SM)

    # Use commutative subalgebra
    algebra_SM = build_SM_algebra_generators_commutative(basis_SM, name_to_index_SM)
    ops_SM = [elem.matrix for elem in algebra_SM]
    labels_SM = [elem.name for elem in algebra_SM]

    D_SM = build_DF_SM_one_gen(
        basis_SM,
        name_to_index_SM,
        m_nu=0.1,
        m_e=0.5,
        m_u=2.0,
        m_d=4.0,
    )

    print(f"dim H_F^SM (1 gen + nu_R) = {len(basis_SM)}")
    print("Algebra generators (commutative subset):", labels_SM)
    print()

    test_first_order_condition_generic(D_SM, ops_SM, labels_SM, S_SM, eps=1e-12)
    test_zero_order_condition_generic(ops_SM, labels_SM, S_SM, eps=1e-12)
    test_grading_and_reality_generic(D_SM, ops_SM, labels_SM, gamma_SM, S_SM)


if __name__ == "__main__":
    # print(">>> Baseline SM-like 1gen finite triple (fixed Yukawas) <<<\n")
    # run_sm_ncg_tests()

    print("\n\n>>> Emergent alignment Yukawas plugged into NCG tests <<<\n")
    run_ncg_with_alignment()

"""
>>> Emergent alignment Yukawas plugged into NCG tests <<<

=== Running emergent alignment pipeline to obtain Yukawas ===

=== NCG tests: 3-gen internal toy triple with aligned Yukawas ===
=== First-order condition test ===
Max Frobenius norm over all pairs (a,b): 0.000e+00
Pairs with norm < 1.0e-12:
  (a=         I, b=         I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         I, b=  Q_sector) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         I, b=P_sector_u) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         I, b=P_sector_d) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         I, b=P_sector_e) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         I, b=P_sector_nu) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=  Q_sector, b=         I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=  Q_sector, b=  Q_sector) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=  Q_sector, b=P_sector_u) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=  Q_sector, b=P_sector_d) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=  Q_sector, b=P_sector_e) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=  Q_sector, b=P_sector_nu) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_u, b=         I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_u, b=  Q_sector) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_u, b=P_sector_u) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_u, b=P_sector_d) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_u, b=P_sector_e) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_u, b=P_sector_nu) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_d, b=         I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_d, b=  Q_sector) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_d, b=P_sector_u) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_d, b=P_sector_d) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_d, b=P_sector_e) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_d, b=P_sector_nu) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_e, b=         I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_e, b=  Q_sector) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_e, b=P_sector_u) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_e, b=P_sector_d) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_e, b=P_sector_e) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_e, b=P_sector_nu) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_nu, b=         I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_nu, b=  Q_sector) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_nu, b=P_sector_u) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_nu, b=P_sector_d) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_nu, b=P_sector_e) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_nu, b=P_sector_nu) → ||[[D,a],J b J^-1]||_F = 0.000e+00

=== Zero-order condition test ===
Max Frobenius norm over all pairs (a,b): 0.000e+00
All pairs satisfy [a, J b J^-1] ≈ 0 within eps=1.0e-12

=== Grading & reality tests ===
||{gamma_F, D_F}||_F = 0.000e+00
max ||[gamma_F, a]||_F over a∈A_F = 0.000e+00
||S^2 - I||_F  (⇒ J_F^2 deviation) = 0.000e+00
||J D_F J^-1 - D_F||_F   = 1.065e-01
||J D_F J^-1 + D_F||_F   = 1.075e-01
→ KO-sign: J D_F J^-1 = + D_F (J-even Dirac operator)


=== NCG tests: SM-like 1gen triple with effective aligned Yukawas ===
=== First-order condition test (generic J) ===
Max Frobenius norm over all pairs (a,b): 1.843e-01
Pairs with norm < 1.0e-12:
  (a=           I, b=           I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=           I, b=    H_sigma3) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=           I, b=color_lambda3) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=           I, b=color_lambda8) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=    H_sigma3, b=           I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=    H_sigma3, b=color_lambda3) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=    H_sigma3, b=color_lambda8) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=color_lambda3, b=           I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=color_lambda3, b=    H_sigma3) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=color_lambda3, b=color_lambda3) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=color_lambda3, b=color_lambda8) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=color_lambda8, b=           I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=color_lambda8, b=    H_sigma3) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=color_lambda8, b=color_lambda3) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=color_lambda8, b=color_lambda8) → ||[[D,a],J b J^-1]||_F = 0.000e+00

=== Zero-order condition test (generic J) ===
Max Frobenius norm over all pairs (a,b): 0.000e+00
All pairs satisfy [a, J b J^-1] ≈ 0 within eps=1.0e-12

=== Grading & reality tests (generic γ,J) ===
||{gamma, D_F}||_F = 0.000e+00
max ||[gamma, a]||_F over a∈A_F = 0.000e+00
||S^2 - I||_F  (⇒ J^2 deviation) = 0.000e+00
||J D_F J^-1 - D_F||_F   = 0.000e+00
||J D_F J^-1 + D_F||_F   = 3.685e-01
→ KO-sign: J D_F J^-1 = + D_F (J-even Dirac operator)
"""

#!/usr/bin/env python3
import numpy as np
from numpy.linalg import svd, solve, cond, pinv


# =========================
# Global config
# =========================

class Config:
    v = 246.0        # GeV
    mu0 = 1.0e12     # GeV
    mu_EW = 91.1876  # GeV
    Lambda_Maj = 1.0e14  # GeV (overall heavy Majorana scale)

    # Alignment scale: Fibonacci / 360
    # kappa = 360/89, eps = 1/kappa
    kappa = 360.0 / 89.0
    eps = .24

    seed = 12345  # overwritten per run

    t0 = np.log(mu0)
    t1 = np.log(mu_EW)
    dt = -0.01  # log-scale step size

    # Higgs quartic (approx EW value, treated constant here)
    lam = 0.13


# Global indices for light / heavy sites (9 = 3 + 6)
LIGHT = slice(0, 3)
HEAVY = slice(3, 9)

# Names for each observable in the same order as make_observables()
observable_names = [
    # mass ratios
    "m_c/m_t",
    "m_u/m_t",
    "m_s/m_b",
    "m_d/m_b",
    "m_mu/m_tau",
    "m_e/m_tau",
    # CKM angles
    "theta12_q (rad)",
    "theta23_q (rad)",
    "theta13_q (rad)",
    # PMNS angles
    "theta12_l (rad)",
    "theta23_l (rad)",
    "theta13_l (rad)",
    # neutrino splittings
    "Delta m2_21 (eV^2)",
    "Delta m2_31 (eV^2)",
]


# =========================
# Utility
# =========================

def has_bad(x: np.ndarray) -> bool:
    """Check for NaN or Inf in an array."""
    return np.any(np.isnan(x)) or np.any(np.isinf(x))


def chi2_breakdown(res):
    """
    Return per-observable χ² contributions as a list of dicts:
      {"name", "theory", "exp", "sigma", "chi2_i"}
    """
    x_th = make_observables(res)
    diffs = x_th - x_exp
    chi2_i = (diffs / sigma) ** 2

    breakdown = []
    for name, th, exp, sig, c2 in zip(observable_names, x_th, x_exp, sigma, chi2_i):
        breakdown.append({
            "name": name,
            "theory": th,
            "exp": exp,
            "sigma": sig,
            "chi2_i": c2,
        })
    return breakdown


# =========================
# Alignment kernel K (9x9)
# =========================

def build_alignment_kernel(eps: float, N: int = 9,
                           allowed_distances=(1, 2, 3, 4, 5, 6, 8)) -> np.ndarray:
    """
    Build the 9x9 alignment kernel:
      K_ij = eps^{|i-j|} for |i-j| in allowed_distances,
      K_ij = 0 for other off-diagonals,
      K_ii = 1.

    We use chain distance d = |i-j| on a 9-site chain.
    """
    allowed_distances = set(allowed_distances)
    K = np.zeros((N, N), dtype=float)
    for i in range(N):
        for j in range(N):
            d = abs(i - j)
            if d == 0:
                K[i, j] = 1.0
            elif d in allowed_distances:
                K[i, j] = eps ** d
            else:
                K[i, j] = 0.0
    return K


# =========================
# Proto-matrices
# =========================

def random_complex_matrix(shape, rng):
    real = rng.normal(0.0, 1.0, size=shape)
    imag = rng.normal(0.0, 1.0, size=shape)
    return (real + 1j * imag) / np.sqrt(2.0)


def normalize_by_largest_singular_value(X: np.ndarray) -> np.ndarray:
    s = svd(X, compute_uv=False)
    s_max = np.max(s)
    if s_max == 0:
        return X
    return X / s_max


def random_weighted_proto(shape, rng, site_scales):
    """
    Gaussian proto-matrix with site-dependent variances:
      Var[X_ij] ~ site_scales[i] * site_scales[j].
    Then normalized so largest singular value = 1.
    """
    site_scales = np.asarray(site_scales, dtype=float)
    S = np.outer(site_scales, site_scales)
    real = rng.normal(0.0, S)
    imag = rng.normal(0.0, S)
    X = (real + 1j * imag) / np.sqrt(2.0)
    return normalize_by_largest_singular_value(X)


def build_site_scales_from_generations(gen_pattern):
    """
    gen_pattern: length-3 array [s1, s2, s3] for 'generation' (1,2,3).

    We assign site_scales[i] = gen_pattern[i % 3] for a 9-site chain.
    This enforces triadic repetition:
      sites (0,3,6)->s1, (1,4,7)->s2, (2,5,8)->s3.
    """
    gen_pattern = np.array(gen_pattern, dtype=float)
    scales = np.zeros(9, dtype=float)
    for i in range(9):
        g = i % 3
        scales[i] = gen_pattern[g]
    return scales


def generate_proto_matrices(cfg: Config, use_site_hierarchy: bool = True):
    """
    Generate proto Yukawa and Majorana matrices on the 9-site proto-flavor space.

    If use_site_hierarchy is False, all site_scales are set to 1.0 and
    the only structure comes from the alignment kernel.
    """
    rng = np.random.default_rng(cfg.seed)
    eps = cfg.eps

    if use_site_hierarchy:
        # sector-dependent generation patterns (gen1, gen2, gen3)
        # up-type: strong hierarchy
        gen_u = [eps ** 4, eps ** 2, 1.0]
        # down-type: moderate hierarchy
        gen_d = [eps ** 3, eps, 1.0]
        # charged leptons: similar to down
        gen_e = [eps ** 3, eps, 1.0]
        # neutrino Dirac: weak hierarchy
        gen_nu = [eps, 1.0, 1.0]
    else:
        # flat proto, no extra hierarchy
        gen_u = gen_d = gen_e = gen_nu = [1.0, 1.0, 1.0]

    # build 9-site scales with triadic pattern
    site_scales_u = build_site_scales_from_generations(gen_u)
    site_scales_d = build_site_scales_from_generations(gen_d)
    site_scales_e = build_site_scales_from_generations(gen_e)
    site_scales_nu = build_site_scales_from_generations(gen_nu)

    # draw weighted proto-matrices
    Yu0 = random_weighted_proto((9, 9), rng, site_scales_u)
    Yd0 = random_weighted_proto((9, 9), rng, site_scales_d)
    Ye0 = random_weighted_proto((9, 9), rng, site_scales_e)
    Ynu0 = random_weighted_proto((9, 9), rng, site_scales_nu)

    # Majorana proto: O(1) and symmetric, no extra site hierarchy yet
    M0 = random_complex_matrix((9, 9), rng)
    M0 = normalize_by_largest_singular_value(M0)
    M0 = 0.5 * (M0 + M0.T)

    # optional overall Yukawa scale (kept at 1.0: rescaling is done later)
    yukawa_scale = 1.0
    Yu0 *= yukawa_scale
    Yd0 *= yukawa_scale
    Ye0 *= yukawa_scale
    Ynu0 *= yukawa_scale

    return Yu0, Yd0, Ye0, Ynu0, M0


# =========================
# Alignment Φ: K ⊙ X
# =========================

def apply_alignment(K: np.ndarray, X: np.ndarray) -> np.ndarray:
    """Hadamard (elementwise) alignment: Φ(X) = K ⊙ X."""
    return K * X


def align_all(K, Yu0, Yd0, Ye0, Ynu0, M0):
    Yu9 = apply_alignment(K, Yu0)
    Yd9 = apply_alignment(K, Yd0)
    Ye9 = apply_alignment(K, Ye0)
    Ynu9 = apply_alignment(K, Ynu0)
    M9 = apply_alignment(K, M0)
    return Yu9, Yd9, Ye9, Ynu9, M9


# =========================
# Schur complement 9→3
# =========================

def schur_9_to_3(Y9: np.ndarray, cond_tol: float = 1e12) -> np.ndarray:
    """
    Y9 is 9x9. Light sites: 0,1,2; heavy: 3..8.
    Effective 3x3 Yukawa via Schur complement:
      Y_eff = A - B D^{-1} B†.

    If D is ill-conditioned, uses pseudo-inverse.
    """
    A = Y9[LIGHT, LIGHT]
    B = Y9[LIGHT, HEAVY]
    D = Y9[HEAVY, HEAVY]

    if cond(D) > cond_tol:
        # fall back to pseudo-inverse
        D_inv = pinv(D)
        Y_eff = A - B @ D_inv @ B.conj().T
    else:
        X = solve(D, B.conj().T)  # D X = B†
        Y_eff = A - B @ X
    return Y_eff


def reduce_all_dirac(Yu9, Yd9, Ye9, Ynu9):
    Yu_eff = schur_9_to_3(Yu9)
    Yd_eff = schur_9_to_3(Yd9)
    Ye_eff = schur_9_to_3(Ye9)
    Ynu_eff = schur_9_to_3(Ynu9)
    return Yu_eff, Yd_eff, Ye_eff, Ynu_eff


# =========================
# Majorana sector: triadic projection 6→3
# =========================

def heavy_block(M9: np.ndarray) -> np.ndarray:
    """Extract 6x6 heavy block (sites 3..8, 0-based)."""
    return M9[HEAVY, HEAVY]


def triad_heavy_basis(Nh: int = 6, ks=(1, 2, 3)) -> np.ndarray:
    """
    Build a triadic basis in heavy space using DFT-like modes k in ks.

    Returns an Nh x len(ks) matrix with orthonormal columns.
    Default ks=(1,2,3) matches the original choice.
    """
    i = np.arange(Nh)
    basis = []
    for k in ks:
        vec = np.exp(2j * np.pi * k * i / Nh)
        vec /= np.linalg.norm(vec)
        basis.append(vec)
    return np.stack(basis, axis=1)


def build_M_R_triadic(M9_aligned: np.ndarray,
                      Lambda_Maj: float,
                      ks=(1, 2, 3)) -> np.ndarray:
    """
    9x9 aligned Majorana → 6x6 heavy block → triadic 3x3 projection.

    M_R = Λ_Maj * B_H† M_H B_H, symmetrized.
    """
    M_H = heavy_block(M9_aligned)  # 6x6
    B_H = triad_heavy_basis(6, ks)  # 6x3
    M3 = B_H.conj().T @ M_H @ B_H  # 3x3
    M3 = 0.5 * (M3 + M3.T)         # enforce symmetry
    M_R = Lambda_Maj * M3
    return M_R


def seesaw_light_neutrinos(Ynu_eff: np.ndarray,
                           M_R: np.ndarray,
                           v: float,
                           cond_tol: float = 1e12) -> np.ndarray:
    """
    Type-I seesaw:
      m_D = v/√2 Ynu_eff,
      m_ν = - m_D M_R^{-1} m_D^T (symmetric 3x3, in GeV).
    """
    m_D = (v / np.sqrt(2.0)) * Ynu_eff
    if cond(M_R) > cond_tol:
        M_R_inv = pinv(M_R)
        m_nu = -m_D @ M_R_inv @ m_D.T
    else:
        X = solve(M_R, m_D.T)
        m_nu = -m_D @ X

    m_nu = 0.5 * (m_nu + m_nu.T)
    return m_nu


# =========================
# 1-loop Yukawa RGEs (g frozen)
# + Weinberg operator RGE
# =========================

def beta_Yukawas(Yu, Yd, Ye, Ynu, g1, g2, g3):
    """
    1-loop SM Yukawa RGEs (in matrix form), with fixed gauge couplings.
    """
    # Safeguard: Return zero betas if inputs contain NaN or Inf
    if any(has_bad(M) for M in (Yu, Yd, Ye, Ynu)):
        Z = np.zeros_like(Yu)
        return Z, Z, Z, Z

    Yu_dagYu = Yu.conj().T @ Yu
    Yd_dagYd = Yd.conj().T @ Yd
    Ye_dagYe = Ye.conj().T @ Ye
    Ynu_dagYnu = Ynu.conj().T @ Ynu

    T = np.trace(3 * Yu_dagYu + 3 * Yd_dagYd + Ye_dagYe)

    factor_u = T - (17 / 20 * g1 ** 2 + 9 / 4 * g2 ** 2 + 8 * g3 ** 2)
    factor_d = T - (1 / 4 * g1 ** 2 + 9 / 4 * g2 ** 2 + 8 * g3 ** 2)
    factor_e = T - (9 / 4 * g1 ** 2 + 9 / 4 * g2 ** 2)
    factor_nu = T - (9 / 20 * g1 ** 2 + 9 / 4 * g2 ** 2)

    dYu = Yu * factor_u + (3 / 2) * (Yu @ Yu_dagYu - Yd @ (Yd_dagYd @ Yu))
    dYd = Yd * factor_d + (3 / 2) * (Yd @ Yd_dagYd - Yu @ (Yu_dagYu @ Yd))
    dYe = Ye * factor_e + (3 / 2) * (Ye @ Ye_dagYe)
    dYnu = Ynu * factor_nu + (3 / 2) * (Ynu @ Ynu_dagYnu - Ye @ (Ye_dagYe @ Ynu))

    dYu /= (16 * np.pi ** 2)
    dYd /= (16 * np.pi ** 2)
    dYe /= (16 * np.pi ** 2)
    dYnu /= (16 * np.pi ** 2)

    return dYu, dYd, dYe, dYnu


def beta_kappa_L(kappa_L, Yu, Ye, g2, lam):
    """
    16π² dκ_L/dt = (-3 g2² + 2λ + 6 Tr(Yu†Yu)) κ_L
                   - 3/2 (Ye†Ye κ_L + κ_L (Ye†Ye)^T).

    We treat λ as constant, and ignore g1,g3 in this operator.
    """
    if any(has_bad(M) for M in (kappa_L, Yu, Ye)):
        return np.zeros_like(kappa_L)

    Yu_dagYu = Yu.conj().T @ Yu
    Ye_dagYe = Ye.conj().T @ Ye
    T_u = np.trace(Yu_dagYu)

    pref = (-3 * g2 ** 2 + 2 * lam + 6 * T_u)
    term1 = pref * kappa_L
    term2 = -1.5 * (Ye_dagYe @ kappa_L + kappa_L @ Ye_dagYe.T.conj())

    dkappa = (term1 + term2) / (16 * np.pi ** 2)
    return dkappa


def rk4_step_full(Yu, Yd, Ye, Ynu, kappa_L, g1, g2, g3, lam, dt):
    """
    RK4 step evolving Yukawas + κ_L with fixed (g1,g2,g3,lam).
    """
    # k1
    dYu1, dYd1, dYe1, dYnu1 = beta_Yukawas(Yu, Yd, Ye, Ynu, g1, g2, g3)
    dkL1 = beta_kappa_L(kappa_L, Yu, Ye, g2, lam)

    # k2
    Yu2 = Yu + 0.5 * dt * dYu1
    Yd2 = Yd + 0.5 * dt * dYd1
    Ye2 = Ye + 0.5 * dt * dYe1
    Ynu2 = Ynu + 0.5 * dt * dYnu1
    kL2 = kappa_L + 0.5 * dt * dkL1

    dYu2, dYd2, dYe2, dYnu2 = beta_Yukawas(Yu2, Yd2, Ye2, Ynu2, g1, g2, g3)
    dkL2 = beta_kappa_L(kL2, Yu2, Ye2, g2, lam)

    # k3
    Yu3 = Yu + 0.5 * dt * dYu2
    Yd3 = Yd + 0.5 * dt * dYd2
    Ye3 = Ye + 0.5 * dt * dYe2
    Ynu3 = Ynu + 0.5 * dt * dYnu2
    kL3 = kappa_L + 0.5 * dt * dkL2

    dYu3, dYd3, dYe3, dYnu3 = beta_Yukawas(Yu3, Yd3, Ye3, Ynu3, g1, g2, g3)
    dkL3 = beta_kappa_L(kL3, Yu3, Ye3, g2, lam)

    # k4
    Yu4 = Yu + dt * dYu3
    Yd4 = Yd + dt * dYd3
    Ye4 = Ye + dt * dYe3
    Ynu4 = Ynu + dt * dYnu3
    kL4 = kappa_L + dt * dkL3

    dYu4, dYd4, dYe4, dYnu4 = beta_Yukawas(Yu4, Yd4, Ye4, Ynu4, g1, g2, g3)
    dkL4 = beta_kappa_L(kL4, Yu4, Ye4, g2, lam)

    Yu_next = Yu + (dt / 6.0) * (dYu1 + 2 * dYu2 + 2 * dYu3 + dYu4)
    Yd_next = Yd + (dt / 6.0) * (dYd1 + 2 * dYd2 + 2 * dYd3 + dYd4)
    Ye_next = Ye + (dt / 6.0) * (dYe1 + 2 * dYe2 + 2 * dYe3 + dYe4)
    Ynu_next = Ynu + (dt / 6.0) * (dYnu1 + 2 * dYnu2 + 2 * dYnu3 + dYnu4)
    kL_next = kappa_L + (dt / 6.0) * (dkL1 + 2 * dkL2 + 2 * dkL3 + dkL4)

    return Yu_next, Yd_next, Ye_next, Ynu_next, kL_next


def run_RGE_full(Yu0, Yd0, Ye0, Ynu0, kappa_L0,
                 g1_const, g2_const, g3_const,
                 cfg: Config):
    """
    Evolve Yukawas and Weinberg operator from μ0 down to μ_EW
    with fixed gauge couplings.
    """
    Yu, Yd, Ye, Ynu = Yu0.copy(), Yd0.copy(), Ye0.copy(), Ynu0.copy()
    kappa_L = kappa_L0.copy()
    g1, g2, g3 = g1_const, g2_const, g3_const
    lam = cfg.lam

    t = cfg.t0
    step = 0
    while (cfg.dt < 0 and t > cfg.t1) or (cfg.dt > 0 and t < cfg.t1):
        step += 1

        Yu, Yd, Ye, Ynu, kappa_L = rk4_step_full(
            Yu, Yd, Ye, Ynu, kappa_L, g1, g2, g3, lam, cfg.dt
        )
        t += cfg.dt

        # Safeguard: Break if NaN/Inf detected (prevents crash)
        if any(has_bad(M) for M in (Yu, Yd, Ye, Ynu, kappa_L)):
            print(f"Warning: NaN/Inf detected at RGE step {step}, halting evolution.")
            break

    return Yu, Yd, Ye, Ynu, kappa_L, g1, g2, g3


# =========================
# Diagonalization and angles
# =========================

def diag_dirac_Y(Y: np.ndarray, v: float):
    """
    SVD for Dirac Yukawa:
      Y = U_L diag(s) U_R†,  masses = v/√2 * s.
    """
    U_L, s, U_Rh = svd(Y)
    masses = (v / np.sqrt(2.0)) * s
    return U_L, np.diag(s), U_Rh, masses


def takagi_symmetric(m: np.ndarray):
    """
    Takagi factorization via SVD for complex symmetric 3x3:
      m = U diag(s) U^T, with s ≥ 0.
    """
    U, s, Vh = svd(m)
    return U, s


def diagonalize_all(Yu, Yd, Ye, mnu, v):
    UuL, Yu_diag, UuR, mu = diag_dirac_Y(Yu, v)
    UdL, Yd_diag, UdR, md = diag_dirac_Y(Yd, v)
    UlL, Ye_diag, UlR, me = diag_dirac_Y(Ye, v)

    U_nu, mnu_vals = takagi_symmetric(mnu)
    mnu_masses = mnu_vals  # in GeV

    Vckm = UuL.conj().T @ UdL
    Vpmns = UlL.conj().T @ U_nu

    return mu, md, me, mnu_masses, Vckm, Vpmns


def extract_angles_and_phase(V: np.ndarray):
    """
    Extract approximate mixing angles (in radians) and Dirac phase
    from a 3x3 unitary matrix V, assuming a PDG-like parameterization.
    """
    s13 = np.abs(V[0, 2])
    theta13 = np.arcsin(np.clip(s13, 0.0, 1.0))

    s12 = np.abs(V[0, 1])
    c12 = np.abs(V[0, 0])
    theta12 = np.arctan2(s12, c12)

    s23 = np.abs(V[1, 2])
    c23 = np.abs(V[2, 2])
    theta23 = np.arctan2(s23, c23)

    J = np.imag(V[0, 0] * V[1, 1] * np.conj(V[0, 1]) * np.conj(V[1, 0]))
    denom = (np.sin(2 * theta12) * np.sin(2 * theta23) *
             np.sin(2 * theta13) * np.cos(theta13))
    if np.abs(denom) < 1e-12:
        delta = 0.0
    else:
        x = 8 * J / denom
        x = np.clip(x, -1.0, 1.0)
        delta = np.arcsin(x)

    return theta12, theta23, theta13, delta


def neutrino_splittings(mnu_masses: np.ndarray):
    """
    Compute Δm²_21 and Δm²_31 in GeV² from the (non-negative) Takagi singular values.
    """
    m_sorted = np.sort(mnu_masses)
    m1, m2, m3 = m_sorted
    dm2_21 = m2 ** 2 - m1 ** 2
    dm2_31 = m3 ** 2 - m1 ** 2
    return dm2_21, dm2_31  # GeV^2


# =========================
# χ² and observables
# =========================

def chi2(observed, expected, sigma):
    return np.sum(((observed - expected) / sigma) ** 2)


def rescale_yukawa_sector(Y, v, m_target_heaviest):
    """
    Rescale Y so that the heaviest mass eigenvalue (v/√2 * max singular value)
    matches m_target_heaviest. Returns (Y_rescaled, alpha).
    """
    U_L, s, U_Rh = svd(Y)
    m_current = (v / np.sqrt(2.0)) * np.max(s)
    if m_current == 0:
        return Y, 1.0
    alpha = m_target_heaviest / m_current
    return alpha * Y, alpha


# experimental targets (rough)
x_exp = np.array([
    # mass ratios
    0.007,    # m_c/m_t
    1e-5,     # m_u/m_t
    0.02,     # m_s/m_b
    0.001,    # m_d/m_b
    0.06,     # m_mu/m_tau
    0.0003,   # m_e/m_tau
    # CKM angles (rad)
    0.226, 0.041, 0.0035,
    # PMNS angles (rad)
    0.59, 0.84, 0.15,
    # Δm² (eV²)
    7.4e-5, 2.5e-3
])

sigma = np.array([
    0.5 * x_exp[0], 0.5 * x_exp[1], 0.5 * x_exp[2], 0.5 * x_exp[3],
    0.5 * x_exp[4], 0.5 * x_exp[5],
    0.1 * x_exp[6], 0.1 * x_exp[7], 0.1 * x_exp[8],
    0.1 * x_exp[9], 0.1 * x_exp[10], 0.1 * x_exp[11],
    0.3 * x_exp[12], 0.3 * x_exp[13]
])


def make_observables(res):
    mu, md, me, mnu = res["mu"], res["md"], res["me"], res["mnu"]
    th12_q, th23_q, th13_q = res["th_q"]
    th12_l, th23_l, th13_l = res["th_l"]
    dm2_21, dm2_31 = res["dm2_eV2"]

    # sort ascending so index 2 is heaviest
    mu_sorted = np.sort(mu)
    md_sorted = np.sort(md)
    me_sorted = np.sort(me)

    obs = []

    # mass ratios
    obs.append(mu_sorted[1] / mu_sorted[2])  # m_c/m_t
    obs.append(mu_sorted[0] / mu_sorted[2])  # m_u/m_t
    obs.append(md_sorted[1] / md_sorted[2])  # m_s/m_b
    obs.append(md_sorted[0] / md_sorted[2])  # m_d/m_b
    obs.append(me_sorted[1] / me_sorted[2])  # m_mu/m_tau
    obs.append(me_sorted[0] / me_sorted[2])  # m_e/m_tau

    # CKM
    obs.append(th12_q)
    obs.append(th23_q)
    obs.append(th13_q)

    # PMNS
    obs.append(th12_l)
    obs.append(th23_l)
    obs.append(th13_l)

    # neutrino splittings (eV²)
    obs.append(dm2_21)
    obs.append(dm2_31)

    return np.array(obs)


def chi2_from_res(res):
    x_th = make_observables(res)
    return chi2(x_th, x_exp, sigma)


# =========================
# run_pipeline
# =========================

def run_pipeline(seed: int,
                 cfg: Config,
                 use_RGE: bool = True,
                 use_site_hierarchy: bool = True,
                 triad_ks=(1, 2, 3)):
    """
    Full pipeline:
      - build alignment kernel on 9 sites
      - generate proto matrices
      - apply alignment
      - Schur 9→3 for Dirac Yukawas
      - triadic heavy projection for Majorana
      - seesaw at μ0 to get mν(μ0)
      - build Weinberg operator κ_L(μ0)
      - (optional) run RGEs down to μ_EW
      - rescale sectors to match m_t, m_b, m_tau
      - diagonalize and extract masses, angles, Δm²
      - compute χ²
    """
    cfg.seed = seed

    # 1. kernel
    K = build_alignment_kernel(cfg.eps, N=9)

    # 2. proto
    Yu0, Yd0, Ye0, Ynu0, M0 = generate_proto_matrices(cfg, use_site_hierarchy)

    # 3. alignment
    Yu9, Yd9, Ye9, Ynu9, M9 = align_all(K, Yu0, Yd0, Ye0, Ynu0, M0)

    # 4. Schur
    Yu_eff, Yd_eff, Ye_eff, Ynu_eff = reduce_all_dirac(Yu9, Yd9, Ye9, Ynu9)

    # 5. M_R (triadic heavy projection)
    M_R = build_M_R_triadic(M9, cfg.Lambda_Maj, ks=triad_ks)

    # 6. seesaw at μ0 → mν(μ0) in GeV
    m_nu_0 = seesaw_light_neutrinos(Ynu_eff, M_R, cfg.v)

    # 6b. Weinberg operator κ_L(μ0) (dimensionful, GeV^-1)
    kappa_L_0 = (2.0 / cfg.v ** 2) * m_nu_0

    # 7. RG
    g1_0, g2_0, g3_0 = 0.46, 0.63, 0.88
    if use_RGE:
        Yu_EW, Yd_EW, Ye_EW, Ynu_EW, kappa_L_EW, g1_EW, g2_EW, g3_EW = run_RGE_full(
            Yu_eff, Yd_eff, Ye_eff, Ynu_eff, kappa_L_0, g1_0, g2_0, g3_0, cfg
        )
        # reconstruct mν(μ_EW) from κ_L(μ_EW)
        m_nu_EW = 0.5 * cfg.v ** 2 * kappa_L_EW
    else:
        Yu_EW, Yd_EW, Ye_EW, Ynu_EW = Yu_eff, Yd_eff, Ye_eff, Ynu_eff
        m_nu_EW = m_nu_0
        g1_EW, g2_EW, g3_EW = g1_0, g2_0, g3_0

    # 7b. rescale sectors to fix heavy masses
    m_t_target = 173.0
    m_b_target = 4.18
    m_tau_target = 1.777

    Yu_EW, alpha_u = rescale_yukawa_sector(Yu_EW, cfg.v, m_t_target)
    Yd_EW, alpha_d = rescale_yukawa_sector(Yd_EW, cfg.v, m_b_target)
    Ye_EW, alpha_e = rescale_yukawa_sector(Ye_EW, cfg.v, m_tau_target)

    # 8. diag at μ_EW
    mu, md, me, mnu_masses, Vckm, Vpmns = diagonalize_all(
        Yu_EW, Yd_EW, Ye_EW, m_nu_EW, cfg.v
    )

    # 9. angles, Δm²
    th12_q, th23_q, th13_q, delta_q = extract_angles_and_phase(Vckm)
    th12_l, th23_l, th13_l, delta_l = extract_angles_and_phase(Vpmns)
    dm2_21_GeV2, dm2_31_GeV2 = neutrino_splittings(mnu_masses)
    dm2_21_eV2 = dm2_21_GeV2 * 1e18
    dm2_31_eV2 = dm2_31_GeV2 * 1e18

    res = {
        "mu": mu,
        "md": md,
        "me": me,
        "mnu": mnu_masses,  # GeV
        "Vckm": Vckm,
        "Vpmns": Vpmns,
        "th_q": (th12_q, th23_q, th13_q),
        "delta_q": delta_q,
        "th_l": (th12_l, th23_l, th13_l),
        "delta_l": delta_l,
        "dm2_GeV2": (dm2_21_GeV2, dm2_31_GeV2),
        "dm2_eV2": (dm2_21_eV2, dm2_31_eV2),
        "alphas": (alpha_u, alpha_d, alpha_e),
        "g_EW": (g1_EW, g2_EW, g3_EW),
    }

    res["chi2"] = chi2_from_res(res)
    return res


# =========================
# Scan driver
# =========================

if __name__ == "__main__":
    cfg = Config()
    N_seeds = 10

    all_results = []
    chi2_vals = []

    for seed in range(N_seeds):
        r = run_pipeline(seed, cfg, use_RGE=True, use_site_hierarchy=True)
        all_results.append(r)
        chi2_vals.append(r["chi2"])
        print(f"seed {seed}: chi2 = {r['chi2']:.3g}")

    best_idx = int(np.argmin(chi2_vals))
    best = all_results[best_idx]

    print("\nBest seed:", best_idx)
    print("chi2 =", best["chi2"])
    print("Up masses (GeV):      ", best["mu"])
    print("Down masses (GeV):    ", best["md"])
    print("Lepton masses (GeV):  ", best["me"])
    print("Neutrino masses (GeV):", best["mnu"])
    print("Neutrino masses (eV): ", best["mnu"] * 1e9)
    print("Δm² (eV²):            ", best["dm2_eV2"])
    print("CKM angles (rad):     ", best["th_q"], "δq:", best["delta_q"])
    print("PMNS angles (rad):    ", best["th_l"], "δℓ:", best["delta_l"])

    # Detailed χ² breakdown for the best seed
    print("\n=== χ² breakdown for best seed ===")
    breakdown = chi2_breakdown(best)
    for entry in breakdown:
        name = entry["name"]
        th = entry["theory"]
        exp = entry["exp"]
        sig = entry["sigma"]
        c2 = entry["chi2_i"]
        pull = (th - exp) / sig
        print(f"{name:20s}  th = {th: .4e},  exp = {exp: .4e},  "
              f"sigma = {sig: .4e},  pull = {pull: .2f},  chi2_i = {c2: .2f}")

"""
RESULTS:

seed 0: chi2 = 228
seed 1: chi2 = 7.02e+07
seed 2: chi2 = 4.84e+07
seed 3: chi2 = 6.74e+05
seed 4: chi2 = 3.98e+05
seed 5: chi2 = 536
seed 6: chi2 = 7.05e+04
seed 7: chi2 = 6.63e+03
seed 8: chi2 = 1.25e+05
seed 9: chi2 = 579

Best seed: 0
chi2 = 228.08238187484977
Up masses (GeV):       [1.73000000e+02 1.12546660e+00 1.60312852e-03]
Down masses (GeV):     [4.18000000e+00 4.35341898e-01 8.44471931e-04]
Lepton masses (GeV):   [1.77700000e+00 8.77881787e-02 1.20387099e-03]
Neutrino masses (GeV): [3.00703062e-11 4.33199063e-12 3.09414912e-13]
Neutrino masses (eV):  [0.03007031 0.00433199 0.00030941]
Δm² (eV²):             (np.float64(1.867040525991209e-05), np.float64(0.0009041275802107434))
CKM angles (rad):      (np.float64(0.08722192280615554), np.float64(0.022860362232397437), np.float64(0.004513499188957893)) δq: -0.15945295688566463
PMNS angles (rad):     (np.float64(0.5978378266273189), np.float64(0.20768513503624048), np.float64(0.09163820224628695)) δℓ: -0.3393812556290672

=== χ² breakdown for best seed ===
m_c/m_t               th =  6.5056e-03,  exp =  7.0000e-03,  sigma =  3.5000e-03,  pull = -0.14,  chi2_i =  0.02
m_u/m_t               th =  9.2666e-06,  exp =  1.0000e-05,  sigma =  5.0000e-06,  pull = -0.15,  chi2_i =  0.02
m_s/m_b               th =  1.0415e-01,  exp =  2.0000e-02,  sigma =  1.0000e-02,  pull =  8.41,  chi2_i =  70.81
m_d/m_b               th =  2.0203e-04,  exp =  1.0000e-03,  sigma =  5.0000e-04,  pull = -1.60,  chi2_i =  2.55
m_mu/m_tau            th =  4.9402e-02,  exp =  6.0000e-02,  sigma =  3.0000e-02,  pull = -0.35,  chi2_i =  0.12
m_e/m_tau             th =  6.7747e-04,  exp =  3.0000e-04,  sigma =  1.5000e-04,  pull =  2.52,  chi2_i =  6.33
theta12_q (rad)       th =  8.7222e-02,  exp =  2.2600e-01,  sigma =  2.2600e-02,  pull = -6.14,  chi2_i =  37.71
theta23_q (rad)       th =  2.2860e-02,  exp =  4.1000e-02,  sigma =  4.1000e-03,  pull = -4.42,  chi2_i =  19.57
theta13_q (rad)       th =  4.5135e-03,  exp =  3.5000e-03,  sigma =  3.5000e-04,  pull =  2.90,  chi2_i =  8.39
theta12_l (rad)       th =  5.9784e-01,  exp =  5.9000e-01,  sigma =  5.9000e-02,  pull =  0.13,  chi2_i =  0.02
theta23_l (rad)       th =  2.0769e-01,  exp =  8.4000e-01,  sigma =  8.4000e-02,  pull = -7.53,  chi2_i =  56.66
theta13_l (rad)       th =  9.1638e-02,  exp =  1.5000e-01,  sigma =  1.5000e-02,  pull = -3.89,  chi2_i =  15.14
Delta m2_21 (eV^2)    th =  1.8670e-05,  exp =  7.4000e-05,  sigma =  2.2200e-05,  pull = -2.49,  chi2_i =  6.21
Delta m2_31 (eV^2)    th =  9.0413e-04,  exp =  2.5000e-03,  sigma =  7.5000e-04,  pull = -2.13,  chi2_i =  4.53
"""

import numpy as np
from typing import List, Tuple, Dict
from itertools import combinations

# ============================================================
#  Helpers: D_360, angle quantization, harmonic triad scoring
# ============================================================

def divisors_360() -> np.ndarray:
    """Divisors of 360 used as the allowed harmonic alphabet."""
    return np.array([1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 15, 18,
                     20, 24, 30, 36, 40, 45, 60, 72, 90, 120, 180, 360])


def nearest_divisor_angle(theta: float, divisors: np.ndarray = None) -> Tuple[float, int]:
    """
    Project a mixing angle theta onto the nearest divisor angle 2π/N
    with N ∈ D_360. Returns (theta_proj, N).
    """
    if divisors is None:
        divisors = divisors_360()
    theta = float(theta)
    candidates = 2.0 * np.pi / divisors
    idx = int(np.argmin(np.abs(candidates - theta)))
    return candidates[idx], int(divisors[idx])


def triad_harmonic_score(lam: np.ndarray, triad: Tuple[int, int, int],
                         mode: str = "quark") -> float:
    """
    Harmonic triad score based purely on spectral ratios.

    mode="quark": favor simple rational spacing (1:2, 2:3, etc).
    mode="lepton": favor near-degenerate pair + separated third, with φ-like ratio.
    """
    i, j, k = triad
    lam_i, lam_j, lam_k = lam[i], lam[j], lam[k]
    if lam_k <= 0:
        return -np.inf

    # Sort by value to get ordered gaps
    vals = np.array([lam_i, lam_j, lam_k])
    order = np.argsort(vals)
    li, lj, lk = vals[order]
    span = lk - li
    if span <= 0:
        return -np.inf

    g1 = lj - li
    g2 = lk - lj
    r1 = g1 / span
    r2 = g2 / span

    # Simple rational targets for "even" spacing
    simple_ratios = np.array([1/3, 1/2, 2/3])
    # φ-like target for leptons
    phi = (1 + np.sqrt(5)) / 2.0
    phi_ratio = 1.0 / phi  # ≈ 0.618

    if mode == "quark":
        # Smallest distance of r1,r2 to simple rational fractions
        d1 = np.min(np.abs(simple_ratios - r1))
        d2 = np.min(np.abs(simple_ratios - r2))
        return - (d1 + d2)  # smaller distance = higher score

    elif mode == "lepton":
        # Favor near-degenerate pair (smallest gap) and φ-like split of span
        dgaps = np.array([g1, g2])
        small_gap = np.min(dgaps)
        large_gap = np.max(dgaps)
        if span <= 0:
            return -np.inf
        # degeneracy measure: prefer small small_gap/span
        deg_measure = small_gap / span
        # φ-like measure on normalized large gap
        large_ratio = large_gap / span
        dphi = np.abs(large_ratio - phi_ratio)
        # score: want deg_measure small and dphi small
        return - (deg_measure + dphi)

    else:
        return -np.inf


# ============================================================
# Internal Laplacian scaling config (kept for compatibility)
# ============================================================

class InternalSpectrumConfig:
    """
    Config for rescaling the internal Laplacian eigenvalues
    before feeding them into the universal spectral kernel F_base.
    In the fully emergent version, the effective rescale is derived
    directly from the spectrum and this class is not used for tuning.
    """
    L_rescale_factor: float = 0.3
    max_triad_index: int = 20


def rescale_laplacian_evals(lam_raw: np.ndarray,
                            cfg: InternalSpectrumConfig) -> np.ndarray:
    """Legacy helper; not used in the emergent run(), kept for compatibility."""
    return cfg.L_rescale_factor * lam_raw


# ============================================================
# Emergent triad chooser
# ============================================================

def choose_quark_and_lepton_triads(lam: np.ndarray,
                                   max_triad_index: int = 20):
    """
    Choose quark- and lepton-like triads purely by harmonic spectral criteria.

    - Quark triad: near-rational spacing (simple gap ratios).
    - Lepton triad: near-degenerate pair + separated third with φ-like ratio.
    """
    start = 1  # ignore exact zero mode at 0
    stop = min(max_triad_index, len(lam))
    nonzero_indices = np.arange(start, stop, dtype=int)

    if len(nonzero_indices) < 3:
        raise ValueError("Not enough nonzero eigenvalues to form triads.")

    triads = list(combinations(nonzero_indices, 3))

    best_q, best_q_score = None, -np.inf
    best_l, best_l_score = None, -np.inf

    for triad in triads:
        s_q = triad_harmonic_score(lam, triad, mode="quark")
        if s_q > best_q_score:
            best_q_score = s_q
            best_q = triad

        s_l = triad_harmonic_score(lam, triad, mode="lepton")
        if s_l > best_l_score:
            best_l_score = s_l
            best_l = triad

    triad_quark = np.array(best_q, dtype=int)
    triad_lepton = np.array(best_l, dtype=int)
    return triad_quark, triad_lepton


# ============================================================
#  FlavorNCGOperators: NCG side (largely unchanged)
# ============================================================

class FlavorNCGOperators:
    """
    Master class collecting:
      - Emergent misalignment + internal graph machinery
      - Flavor hierarchy and mixing operators
      - Internal NCG finite Dirac operator + algebra + tests
    """

    SECTORS: List[str] = ["u", "d", "e", "nu"]
    N_GEN: int = 3
    SECTOR_NC: Dict[str, int] = {"u": 3, "d": 3, "e": 1, "nu": 1}

    # Rough SM targets kept ONLY as an external diagnostic (not used for selection)
    TARGETS: Dict[str, Tuple[float, float]] = {
        "mu_mt":     (2.2e-05, 0.5 * 2.2e-05),
        "mc_mt":     (7.5e-03, 0.5 * 7.5e-03),
        "md_mb":     (1.1e-03, 0.5 * 1.1e-03),
        "ms_mb":     (2.2e-02, 0.5 * 2.2e-02),
        "me_mt":     (2.9e-04, 0.5 * 2.9e-04),
        "mmu_mt":    (5.9e-02, 0.5 * 5.9e-02),
        "theta12_q": (0.227,   0.05 * 0.227),
        "theta23_q": (0.041,   0.5  * 0.041),
        "theta13_q": (0.0036,  0.5  * 0.0036),
        "theta12_l": (0.584,   0.1  * 0.584),
        "theta23_l": (0.785,   0.2  * 0.785),
        "theta13_l": (0.15,    0.2  * 0.15),
    }

    # ===========================
    # 1. Internal Hilbert space & D_F (NCG side)
    # ===========================

    def dim_per_chirality(self) -> int:
        return len(self.SECTORS) * self.N_GEN  # 4 * 3 = 12

    def flavor_block_offsets(self) -> Dict[str, int]:
        off: Dict[str, int] = {}
        off["u"]  = 0
        off["d"]  = 3
        off["e"]  = 6
        off["nu"] = 9
        return off

    def build_internal_DF_from_Y(self, Y_u, Y_d, Y_e, Y_nu):
        for name, Y in [("Y_u", Y_u), ("Y_d", Y_d), ("Y_e", Y_e), ("Y_nu", Y_nu)]:
            arr = np.asarray(Y, dtype=complex)
            if arr.shape != (3, 3):
                raise ValueError(f"{name} must be a 3×3 matrix, got shape {arr.shape}.")

        Y_u = np.asarray(Y_u, dtype=complex)
        Y_d = np.asarray(Y_d, dtype=complex)
        Y_e = np.asarray(Y_e, dtype=complex)
        Y_nu = np.asarray(Y_nu, dtype=complex)

        dpc = self.dim_per_chirality()
        dimH = 2 * dpc

        Y_gen = np.zeros((dpc, dpc), dtype=complex)
        gen_off = self.flavor_block_offsets()

        def insert_sector_Y(sector: str, Y_s: np.ndarray) -> None:
            off = gen_off[sector]
            Y_gen[off:off + 3, off:off + 3] = Y_s

        insert_sector_Y("u", Y_u)
        insert_sector_Y("d", Y_d)
        insert_sector_Y("e", Y_e)
        insert_sector_Y("nu", Y_nu)

        Y_block = Y_gen

        D_F = np.zeros((dimH, dimH), dtype=complex)
        D_F[:dpc, dpc:] = Y_block.conj().T
        D_F[dpc:, :dpc] = Y_block
        return D_F

    # --- Real structure & grading ---

    def build_swap_LR(self, dim_left: int) -> np.ndarray:
        S = np.zeros((2 * dim_left, 2 * dim_left), dtype=complex)
        S[:dim_left, dim_left:] = np.eye(dim_left)
        S[dim_left:, :dim_left] = np.eye(dim_left)
        return S

    def build_gamma_F(self, dim_left: int) -> np.ndarray:
        g = np.zeros((2 * dim_left, 2 * dim_left), dtype=complex)
        g[:dim_left, :dim_left] = -np.eye(dim_left)
        g[dim_left:, dim_left:] =  np.eye(dim_left)
        return g

    def build_sector_projectors(self):
        dpc = self.dim_per_chirality()
        dimH = 2 * dpc
        gen_off = self.flavor_block_offsets()

        P: Dict[str, np.ndarray] = {}
        for s in self.SECTORS:
            P_s = np.zeros((dimH, dimH), dtype=complex)
            off = gen_off[s]
            P_s[off:off+3, off:off+3] = np.eye(3)
            P_s[dpc + off:dpc + off + 3, dpc + off:dpc + off + 3] = np.eye(3)
            P[s] = P_s
        return P

    def build_Q_sector(self) -> np.ndarray:
        """
        Simple sector charge operator distinguishing u,d,e,nu.
        In the emergent scheme, these sector charges are fixed only
        at this operator level; generation-wise structure is emergent.
        """
        dpc = self.dim_per_chirality()
        dimH = 2 * dpc
        gen_off = self.flavor_block_offsets()
        charges = {"u": 2.0, "d": 1.0, "e": 0.0, "nu": -1.0}

        Q = np.zeros((dimH, dimH), dtype=complex)
        for s in self.SECTORS:
            off = gen_off[s]
            q = charges[s]
            Q[off:off+3, off:off+3] = q * np.eye(3)
            Q[dpc + off:dpc + off + 3, dpc + off:dpc + off + 3] = q * np.eye(3)
        return Q

    def build_internal_algebra_ops(self) -> Tuple[List[np.ndarray], List[str]]:
        dpc = self.dim_per_chirality()
        dimH = 2 * dpc

        I = np.eye(dimH, dtype=complex)
        Q = self.build_Q_sector()
        P = self.build_sector_projectors()

        ops: List[np.ndarray] = [I, Q, P["u"], P["d"], P["e"], P["nu"]]
        labels: List[str] = ["I", "Q_sector", "P_sector_u", "P_sector_d", "P_sector_e", "P_sector_nu"]
        return ops, labels

    # --- NCG tests and alignment score ---

    def J_action_from_swap(self, S: np.ndarray, M: np.ndarray) -> np.ndarray:
        return S @ M.conj() @ S.T

    def test_first_order_condition(
        self, D_F: np.ndarray, ops: List[np.ndarray], labels: List[str], eps: float = 1e-12
    ) -> None:
        n = D_F.shape[0]
        assert D_F.shape == (n, n)
        S = self.build_swap_LR(dim_left=n // 2)

        print("=== First-order condition test ===")
        max_norm = 0.0
        good_pairs = []

        for i, a in enumerate(ops):
            Da = D_F @ a - a @ D_F
            for j, b in enumerate(ops):
                b_tilde = self.J_action_from_swap(S, b)
                comm2 = Da @ b_tilde - b_tilde @ Da
                norm = np.linalg.norm(comm2, ord="fro")
                if norm > max_norm:
                    max_norm = norm
                if norm < eps:
                    good_pairs.append((labels[i], labels[j], norm))

        print(f"Max Frobenius norm over all pairs (a,b): {max_norm:.3e}")
        if good_pairs:
            print(f"Pairs with norm < {eps:.1e}:")
            for la, lb, nrm in good_pairs:
                print(f"  (a={la:>10s}, b={lb:>10s}) → ||[[D,a],J b J^-1]||_F = {nrm:.3e}")
        else:
            print(f"No pairs with norm < {eps:.1e}")
        print()

    def test_zero_order_condition(
        self, ops: List[np.ndarray], labels: List[str], eps: float = 1e-12
    ) -> None:
        n = ops[0].shape[0]
        S = self.build_swap_LR(dim_left=n // 2)

        print("=== Zero-order condition test ===")
        max_norm = 0.0
        bad_pairs = []

        for i, a in enumerate(ops):
            for j, b in enumerate(ops):
                b_tilde = self.J_action_from_swap(S, b)
                comm = a @ b_tilde - b_tilde @ a
                norm = np.linalg.norm(comm, ord="fro")
                if norm > max_norm:
                    max_norm = norm
                if norm > eps:
                    bad_pairs.append((labels[i], labels[j], norm))

        print(f"Max Frobenius norm over all pairs (a,b): {max_norm:.3e}")
        if bad_pairs:
            print("Pairs with significant violation:")
            for la, lb, nrm in bad_pairs:
                print(f"  (a={la:>10s}, b={lb:>10s}) → ||[a, J b J^-1]||_F = {nrm:.3e}")
        else:
            print(f"All pairs satisfy [a, J b J^-1] ≈ 0 within eps={eps:.1e}")
        print()

    def test_grading_and_reality(
        self, D_F: np.ndarray, ops: List[np.ndarray], labels: List[str]
    ) -> None:
        n = D_F.shape[0]
        dpc = n // 2
        gamma_F = self.build_gamma_F(dpc)
        S = self.build_swap_LR(dpc)

        print("=== Grading & reality tests ===")
        anti = gamma_F @ D_F + D_F @ gamma_F
        print(f"||{{gamma_F, D_F}}||_F = {np.linalg.norm(anti, ord='fro'):.3e}")

        max_comm_gamma = 0.0
        for a in ops:
            comm_ga = gamma_F @ a - a @ gamma_F
            max_comm_gamma = max(max_comm_gamma, np.linalg.norm(comm_ga, ord="fro"))
        print(f"max ||[gamma_F, a]||_F over a∈A_F = {max_comm_gamma:.3e}")

        S2 = S @ S
        print(f"||S^2 - I||_F  (⇒ J_F^2 deviation) = {np.linalg.norm(S2 - np.eye(n), ord='fro'):.3e}")

        JDJ = S @ D_F.conj() @ S.T
        norm_minus = np.linalg.norm(JDJ - D_F, ord="fro")
        norm_plus  = np.linalg.norm(JDJ + D_F, ord="fro")
        print(f"||J D_F J^-1 - D_F||_F   = {norm_minus:.3e}")
        print(f"||J D_F J^-1 + D_F||_F   = {norm_plus:.3e}")
        if norm_minus < norm_plus:
            print("→ KO-sign: J D_F J^-1 = + D_F (J-even Dirac operator)")
        else:
            print("→ KO-sign: J D_F J^-1 = - D_F (J-odd Dirac operator)")
        print()

    def ncg_alignment_score(self, D_F: np.ndarray, ops: List[np.ndarray]) -> float:
        """
        Scalar NCG coherence measure (smaller = more aligned).
        Combines grading, zero-order, and first-order deviations.
        """
        n = D_F.shape[0]
        dpc = n // 2
        gamma_F = self.build_gamma_F(dpc)
        S = self.build_swap_LR(dpc)

        # Grading: {γ, D_F} ≈ 0
        anti = gamma_F @ D_F + D_F @ gamma_F
        norm_anti = np.linalg.norm(anti, ord='fro')

        # gamma commutators with algebra
        max_comm_gamma = 0.0
        for a in ops:
            comm_ga = gamma_F @ a - a @ gamma_F
            max_comm_gamma = max(max_comm_gamma, np.linalg.norm(comm_ga, ord="fro"))

        # zero- and first-order
        max_zero = 0.0
        max_first = 0.0
        for i, a in enumerate(ops):
            Da = D_F @ a - a @ D_F
            for j, b in enumerate(ops):
                b_tilde = self.J_action_from_swap(S, b)
                # zero-order
                comm0 = a @ b_tilde - b_tilde @ a
                max_zero = max(max_zero, np.linalg.norm(comm0, ord="fro"))
                # first-order
                comm2 = Da @ b_tilde - b_tilde @ Da
                max_first = max(max_first, np.linalg.norm(comm2, ord="fro"))

        # Simple linear combo; no tunable weights beyond unity
        return norm_anti + max_comm_gamma + max_zero + max_first

    # ===========================
    # 2. Emergent misalignment model, graph, spectrum
    # ===========================

    def allowed_harmonics(self) -> np.ndarray:
        """Allowed global harmonic set D_360."""
        return divisors_360()

    def contextual_harmonics(self, step: int, total_steps: int) -> np.ndarray:
        """
        Contextual selection of subset of D_360 as relaxation proceeds.
        Early time: small subset; late time: full set.
        """
        D = self.allowed_harmonics()
        frac = step / max(total_steps, 1)
        k = int(1 + frac * (len(D) - 1))
        return D[:k]

    def misalignment_energy(self, theta, ns: np.ndarray = None):
        if ns is None:
            ns = self.allowed_harmonics()
        N = len(theta)
        diffs = theta[:, None] - theta[None, :]
        E = 0.0
        for n in ns:
            w_n = 1.0 / n
            E += w_n * np.sum(1.0 - np.cos(n * diffs)) / (N * N)
        return E

    def relax_phases(self, N=200, n_steps=600, eta=0.01, random_seed=42):
        rng = np.random.default_rng(random_seed)
        theta = rng.uniform(0, 2 * np.pi, size=N)
        energy_hist = []

        for step in range(n_steps):
            ns = self.contextual_harmonics(step, n_steps)
            diffs = theta[:, None] - theta[None, :]
            grad = np.zeros(N, dtype=float)

            for n in ns:
                w_n = 1.0 / n
                sin_n = np.sin(n * diffs)
                grad += w_n * n * np.sum(sin_n, axis=1)

            theta = theta - eta * grad
            theta = (theta + 2 * np.pi) % (2 * np.pi)

            if step % 10 == 0 or step == n_steps - 1:
                E = self.misalignment_energy(theta, ns=ns)
                energy_hist.append(E)

        return theta, energy_hist

    def build_emergent_adjacency(self, theta, ns: np.ndarray = None, keep_fraction: float = 0.05):
        """
        Adjacency from the same harmonic set ns used at late-time misalignment.
        Score_ij = Σ_n (1/n) cos(n(θ_i - θ_j)).
        """
        if ns is None:
            ns = self.allowed_harmonics()

        N = len(theta)
        diffs = theta[:, None] - theta[None, :]
        score = np.zeros((N, N), dtype=float)

        for n in ns:
            w_n = 1.0 / n
            score += w_n * np.cos(n * diffs)

        np.fill_diagonal(score, -np.inf)
        triu_idx = np.triu_indices(N, k=1)
        flat_scores = score[triu_idx]
        k = int(keep_fraction * len(flat_scores))
        if k < 1:
            k = 1
        kth_val = np.partition(flat_scores, -k)[-k]
        A = np.zeros((N, N), dtype=float)
        mask = (score >= kth_val)
        A[mask] = 1.0
        A = np.maximum(A, A.T)
        return A

    def largest_connected_component(self, A):
        N = A.shape[0]
        visited = np.zeros(N, dtype=bool)
        best_comp = []
        for i in range(N):
            if not visited[i]:
                stack = [i]
                comp = []
                visited[i] = True
                while stack:
                    v = stack.pop()
                    comp.append(v)
                    neighbors = np.where(A[v] > 0)[0]
                    for w in neighbors:
                        if not visited[w]:
                            visited[w] = True
                            stack.append(w)
                if len(comp) > len(best_comp):
                    best_comp = comp
        best_comp = np.array(best_comp, dtype=int)
        A_sub = A[np.ix_(best_comp, best_comp)]
        return A_sub, best_comp

    def laplacian_from_adjacency(self, A):
        d = np.sum(A, axis=1)
        L = np.diag(d) - A
        return L

    def base_kernel(self, lam, alpha=3.0, form="lambda_sq"):
        """
        Universal base kernel F_base(λ_g):

            F_base(λ_g) = exp[-alpha * (λ_g / λ_ref)^p]

        with λ_ref = smallest positive eigenvalue in the triad.
        alpha will be emergently set from triad spread.
        """
        lam = np.array(lam, dtype=float)
        lam_pos = lam[lam > 0]
        if lam_pos.size == 0:
            lam_ref = 1.0
        else:
            lam_ref = lam_pos.min()
        x = lam / lam_ref
        if form == "lambda_sq":
            return np.exp(-alpha * x**2)
        elif form == "lambda":
            return np.exp(-alpha * x)
        else:
            raise ValueError(f"Unknown kernel form '{form}'")

    def emergent_alpha_for_triad(self, lam_triad: np.ndarray) -> float:
        """
        Derive kernel steepness from the triad itself.
        Use spread in log(λ) to set alpha ~ 1 / Var(log λ).
        """
        lam = np.array(lam_triad, dtype=float)
        lam_pos = lam[lam > 0]
        if lam_pos.size <= 1:
            return 1.0
        logs = np.log(lam_pos)
        var = np.var(logs)
        eps = 1e-6
        alpha = 1.0 / (var + eps)
        return alpha

    def spectral_triad(self, L):
        eigvals, eigvecs = np.linalg.eigh(L)
        idx_sorted = np.argsort(eigvals)
        eigvals_sorted = eigvals[idx_sorted]
        eigvecs_sorted = eigvecs[:, idx_sorted]

        triad_idx = idx_sorted[1:4]
        triad_vals = eigvals_sorted[1:4]

        order = np.argsort(triad_vals)[::-1]  # DESC by λ
        lam_gen = triad_vals[order]
        gen_indices = triad_idx[order]
        return lam_gen, gen_indices, eigvals_sorted

    # ===========================
    # 3. Sector charges, Yukawas, mixing
    # ===========================

    def build_sector_charges_from_spectrum(self, lam: np.ndarray,
                                           triad_quark: np.ndarray,
                                           triad_lepton: np.ndarray) -> Dict[str, np.ndarray]:
        """
        Emergent sector/generation charges from local spectral density
        around each triad eigenvalue.
        """
        lam = np.array(lam, dtype=float)

        def local_density(idx: int) -> float:
            v = lam[idx]
            if v <= 0:
                return 1.0
            window = (lam >= 0.9 * v) & (lam <= 1.1 * v)
            return float(np.sum(window))

        def triad_charges(triad: np.ndarray) -> np.ndarray:
            qs = np.array([local_density(int(i)) for i in triad], dtype=float)
            # log compress
            return np.log1p(qs)

        Q_quark = triad_charges(triad_quark)
        Q_lepton = triad_charges(triad_lepton)

        charges = {
            "u":  Q_quark,
            "d":  Q_quark,
            "e":  Q_lepton,
            "nu": Q_lepton,
        }
        return charges

    def sector_weights(self, F_base: np.ndarray, Q_s: np.ndarray):
        """
        No free β: masses ~ F_base * exp(-Q_s).
        """
        return F_base * np.exp(-Q_s)

    def mass_ratios(self, F_s):
        F_s = np.array(F_s, dtype=float)
        F_s = np.abs(F_s)
        max_val = np.max(F_s)
        if max_val <= 0.0 or not np.isfinite(max_val):
            return 1.0, 1.0
        eps = 1e-16 * max_val
        F_s[F_s < eps] = eps
        m1, m2, m3 = np.sort(F_s)
        return m1 / m3, m2 / m3

    # --- generation operators ---

    def rotation_3d(self, i, j, theta):
        R = np.eye(3, dtype=complex)
        c = np.cos(theta)
        s = np.sin(theta)
        R[i, i] = c
        R[j, j] = c
        R[i, j] = s
        R[j, i] = -s
        return R

    def build_generation_operators(self, phi_order=5, cab_denom=28):
        """
        In the emergent scheme, phi_order and cab_denom are NOT free:
        they are derived from geometric mixing and projected to the
        nearest divisor-based angles before calling this.
        """
        theta_phi = 2 * np.pi / phi_order
        theta_C = 2 * np.pi / cab_denom
        P_phi_12 = self.rotation_3d(0, 1, theta_phi)
        P_phi_23 = self.rotation_3d(1, 2, theta_phi)
        C_12 = self.rotation_3d(0, 1, theta_C)
        return P_phi_12, P_phi_23, C_12, theta_phi, theta_C

    # --- geometric regions and unitaries ---

    def build_geometric_regions(self, theta, n_regions=3):
        phase = np.mod(theta, 2 * np.pi)
        edges = np.linspace(0, 2*np.pi, n_regions+1)
        regions = []
        for k in range(n_regions):
            lo, hi = edges[k], edges[k+1]
            if k < n_regions - 1:
                idx = np.where((phase >= lo) & (phase < hi))[0]
            else:
                idx = np.where((phase >= lo) & (phase <= hi))[0]
            if len(idx) == 0:
                idx = np.array([k % len(theta)], dtype=int)
            regions.append(idx)
        return regions

    def build_geometric_unitary(self, gen_vecs, region_list):
        cols = []
        for R in region_list:
            v = np.sum(gen_vecs[R, :], axis=0)
            norm = np.linalg.norm(v)
            if norm < 1e-14:
                v = np.array([1.0, 0.0, 0.0], dtype=complex)
                norm = 1.0
            cols.append(v / norm)
        U_geom = np.column_stack(cols)
        Uu, _, Vh = np.linalg.svd(U_geom)
        return Uu @ Vh

    def build_sector_bases(self, P_phi_12, P_phi_23, C_12, U_geom,
                           use_neutrino_dressing: bool = True,
                           N_SOLAR: int = 36,
                           N_REACTOR: int = 45,
                           N_ATM: int = 24):
        sector_bases = {}
        U_geom_u = U_geom["u"]
        U_geom_d = U_geom["d"]
        U_geom_e = U_geom["e"]
        U_geom_nu = U_geom["nu"]

        # Quarks: Cabibbo on up-type only
        U_L_u = U_geom_u @ C_12.conj().T
        U_R_u = np.eye(3, dtype=complex)
        U_L_d = U_geom_d
        U_R_d = np.eye(3, dtype=complex)

        # Charged leptons: pure geometry
        U_L_e = U_geom_e
        U_R_e = np.eye(3, dtype=complex)

        # Neutrinos: geometry + golden + 3 discrete rotations
        if use_neutrino_dressing:
            theta_solar = 2 * np.pi / N_SOLAR
            theta_reac = 2 * np.pi / N_REACTOR
            theta_atm = 2 * np.pi / N_ATM

            R_solar = self.rotation_3d(0, 1, theta_solar)
            R_reac = self.rotation_3d(0, 2, theta_reac)
            R_atm = self.rotation_3d(1, 2, theta_atm)

            U_dress = R_atm @ P_phi_23 @ R_solar @ P_phi_12 @ R_reac
            U_L_nu = U_geom_nu @ U_dress
        else:
            U_L_nu = U_geom_nu

        U_R_nu = np.eye(3, dtype=complex)

        sector_bases["u"] = (U_L_u, U_R_u)
        sector_bases["d"] = (U_L_d, U_R_d)
        sector_bases["e"] = (U_L_e, U_R_e)
        sector_bases["nu"] = (U_L_nu, U_R_nu)
        return sector_bases

    def emergent_neutrino_denominators(self, lam_gen_lepton: np.ndarray) -> Tuple[int, int, int]:
        """
        Set N_SOLAR, N_REACTOR, N_ATM from lepton triad degeneracies.
        Smaller gap -> larger N (finer angle).
        """
        lam = np.array(lam_gen_lepton, dtype=float)
        if lam.size != 3:
            return 36, 45, 24

        gaps = np.abs(np.diff(np.sort(lam)))
        # Protect against zero
        gaps = gaps + 1e-8
        inv_gaps = 1.0 / gaps
        inv_gaps /= np.max(inv_gaps)

        # Map to a subset of divisors
        D = divisors_360()
        candidates = D[D <= 90]  # keep it modest

        def map_val(v):
            # v in [0,1] -> candidate index
            idx = int(np.clip(round(v * (len(candidates)-1)), 0, len(candidates)-1))
            return int(candidates[idx])

        N_SOLAR = map_val(inv_gaps[0])   # g12
        N_ATM   = map_val(inv_gaps[-1])  # g23
        N_REACTOR = map_val(0.5 * (inv_gaps[0] + inv_gaps[-1]))
        return N_SOLAR, N_REACTOR, N_ATM

    # --- Yukawas, mixing, diagnostics ---

    def yukawa_from_F_and_UL(self, F_s, U_L, U_R):
        D = np.diag(F_s)
        return U_L @ D @ U_R.conj().T

    def mixing_matrix(self, U_L_up, U_L_down):
        return U_L_up.conj().T @ U_L_down

    def mixing_angles_from_U(self, U):
        s13 = abs(U[0, 2])
        s13 = min(max(s13, 0.0), 1.0)
        theta13 = np.arcsin(s13)
        c13 = np.cos(theta13)
        if abs(c13) < 1e-12:
            theta12 = 0.0
            theta23 = 0.0
        else:
            theta12 = np.arctan2(abs(U[0, 1]), abs(U[0, 0]))
            theta23 = np.arctan2(abs(U[1, 2]), abs(U[2, 2]))
        return theta12, theta23, theta13

    def compute_observables(
        self,
        mu_mt, mc_mt, md_mb, ms_mb, me_mt, mmu_mt,
        theta12_q, theta23_q, theta13_q,
        theta12_l, theta23_l, theta13_l
    ):
        return {
            "mu_mt":     mu_mt,
            "mc_mt":     mc_mt,
            "md_mb":     md_mb,
            "ms_mb":     ms_mb,
            "me_mt":     me_mt,
            "mmu_mt":    mmu_mt,
            "theta12_q": theta12_q,
            "theta23_q": theta23_q,
            "theta13_q": theta13_q,
            "theta12_l": theta12_l,
            "theta23_l": theta23_l,
            "theta13_l": theta13_l,
        }

    def chi2(self, obs, targets=None):
        if targets is None:
            targets = self.TARGETS
        chi2_val = 0.0
        details = []
        for k, v in obs.items():
            target, sigma = targets[k]
            if sigma <= 0:
                continue
            contrib = ((v - target) / sigma)**2
            chi2_val += contrib
            details.append((k, v, target, contrib))
        return chi2_val, details


# ============================================================
# EmergentFlavorNCGModel: FULL EMERGENCE RUN PIPELINE
# ============================================================

class EmergentFlavorNCGModel(FlavorNCGOperators):
    def __init__(
        self,
        N_sites: int = 200,
        n_steps: int = 600,
        eta: float = 0.01,
        keep_fraction: float = 0.05,
    ):
        super().__init__()
        self.N_sites = N_sites
        self.n_steps = n_steps
        self.eta = eta
        self.keep_fraction = keep_fraction

    def run(self):
        # Step 1: relax phases under D_360-driven misalignment
        theta_final, energy_hist = self.relax_phases(
            N=self.N_sites,
            n_steps=self.n_steps,
            eta=self.eta,
            random_seed=42,
        )
        print("Relaxation complete.")
        print(f"Final misalignment energy: {energy_hist[-1]:.6f}")
        print()

        # Harmonics active at final time
        ns_final = self.contextual_harmonics(self.n_steps - 1, self.n_steps)

        # Step 2: emergent adjacency & Laplacian
        A_int_full = self.build_emergent_adjacency(
            theta_final,
            ns=ns_final,
            keep_fraction=self.keep_fraction,
        )
        A_int, nodes = self.largest_connected_component(A_int_full)
        L_int = self.laplacian_from_adjacency(A_int)

        # Spectrum and emergent rescaling
        eigvals_full_raw, eigvecs_full = np.linalg.eigh(L_int)
        pos = eigvals_full_raw[eigvals_full_raw > 1e-12]
        if pos.size > 0:
            L_rescale_factor = 1.0 / pos[0]
        else:
            L_rescale_factor = 1.0
        lam = L_rescale_factor * eigvals_full_raw

        # Emergent triads from harmonic scoring
        triad_quark, triad_lepton = choose_quark_and_lepton_triads(
            lam, max_triad_index=min(90, len(lam))
        )
        lam_gen_quark = lam[triad_quark]
        lam_gen_lepton = lam[triad_lepton]

        # Emergent alpha from triad spread
        alpha_quark = self.emergent_alpha_for_triad(lam_gen_quark)
        alpha_lepton = self.emergent_alpha_for_triad(lam_gen_lepton)

        F_base_quark = self.base_kernel(lam_gen_quark, alpha=alpha_quark, form="lambda_sq")
        F_base_lepton = self.base_kernel(lam_gen_lepton, alpha=alpha_lepton, form="lambda_sq")

        def regularize_F_base(F):
            F = np.array(F, dtype=float)
            max_val = np.max(F)
            if max_val <= 0.0 or not np.isfinite(max_val):
                return np.full_like(F, 1e-16)
            eps = 1e-16 * max_val
            F[F < eps] = eps
            return F

        F_base_quark = regularize_F_base(F_base_quark)
        F_base_lepton = regularize_F_base(F_base_lepton)

        print("=== Emergent internal graph ===")
        print(f"Number of sites: {A_int.shape[0]}")
        print("First 10 eigenvalues of L_int (raw, unscaled):")
        print(eigvals_full_raw[:10])
        print()
        print("Laplacian rescale factor L_rescale_factor =", L_rescale_factor)
        print("Quark triad indices:", triad_quark, "lam_gen_quark:", lam_gen_quark)
        print("Lepton triad indices:", triad_lepton, "lam_gen_lepton:", lam_gen_lepton)
        print("Alpha_quark (emergent):", alpha_quark)
        print("Alpha_lepton (emergent):", alpha_lepton)
        print("Base kernel F_base_quark:", F_base_quark)
        print("Base kernel F_base_lepton:", F_base_lepton)
        print()

        # Generation eigenvectors
        gen_vecs_quark = eigvecs_full[:, triad_quark]
        gen_vecs_lepton = eigvecs_full[:, triad_lepton]

        # Step 3: geometric regions from phase field (restricted to largest component)
        theta_sub = theta_final[nodes]
        regions = self.build_geometric_regions(theta_sub, n_regions=3)
        R0, R1, R2 = regions

        # Quark assignments share region geometry
        assign_u = [R0, R1, R2]
        assign_d = [R0, R1, R2]

        # Sector charges from spectrum
        sector_charges_gen = self.build_sector_charges_from_spectrum(
            lam,
            triad_quark=triad_quark,
            triad_lepton=triad_lepton,
        )

        # Emergent neutrino denominators from lepton triad degeneracies
        N_SOLAR, N_REACTOR, N_ATM = self.emergent_neutrino_denominators(lam_gen_lepton)
        print("Emergent neutrino denominators (SOLAR, REACTOR, ATM):", N_SOLAR, N_REACTOR, N_ATM)
        print()

        # Permutations for leptons (internal alignment selection only)
        perms = [
            (0, 1, 2),
            (0, 2, 1),
            (1, 0, 2),
            (1, 2, 0),
            (2, 0, 1),
            (2, 1, 0),
        ]

        best_align_score = np.inf
        best_perm_e = None
        best_perm_nu = None
        best_U_geom = None
        best_masses = None
        best_angles = None
        best_Ys = None
        best_sector_bases = None
        best_chi2 = None
        best_chi2_details = None

        # Build algebra once for NCG scoring (size known: 24x24)
        ops_A, labels_A = self.build_internal_algebra_ops()

        for pe in perms:
            for pn in perms:
                perm_e = [regions[pe[0]], regions[pe[1]], regions[pe[2]]]
                perm_n = [regions[pn[0]], regions[pn[1]], regions[pn[2]]]

                assign_e = perm_e
                assign_nu = perm_n

                # Geometric unitaries
                U_geom = {
                    "u": self.build_geometric_unitary(gen_vecs_quark, assign_u),
                    "d": self.build_geometric_unitary(gen_vecs_quark, assign_d),
                    "e": self.build_geometric_unitary(gen_vecs_lepton, assign_e),
                    "nu": self.build_geometric_unitary(gen_vecs_lepton, assign_nu),
                }

                # Pure geometric mixing
                V_ckm_geom = self.mixing_matrix(U_geom["u"], U_geom["d"])
                U_pmns_geom = self.mixing_matrix(U_geom["e"], U_geom["nu"])
                theta12_q_geom, theta23_q_geom, theta13_q_geom = self.mixing_angles_from_U(V_ckm_geom)
                theta12_l_geom, theta23_l_geom, theta13_l_geom = self.mixing_angles_from_U(U_pmns_geom)

                # Emergent Cabibbo and golden angles via divisor projection
                theta_C_proj, cab_denom = nearest_divisor_angle(theta12_q_geom)
                theta_phi_proj, phi_order = nearest_divisor_angle(theta12_l_geom)

                P_phi_12, P_phi_23, C_12, theta_phi, theta_C = self.build_generation_operators(
                    phi_order=phi_order, cab_denom=cab_denom
                )

                # Sector weights from spectrum and charges
                F_u = self.sector_weights(F_base_quark, sector_charges_gen["u"])
                F_d = self.sector_weights(F_base_quark, sector_charges_gen["d"])
                F_e = self.sector_weights(F_base_lepton, sector_charges_gen["e"])
                F_n = self.sector_weights(F_base_lepton, sector_charges_gen["nu"])

                # Sector bases: geometry + emergent operators
                sector_bases = self.build_sector_bases(
                    P_phi_12, P_phi_23, C_12,
                    U_geom,
                    use_neutrino_dressing=True,
                    N_SOLAR=N_SOLAR,
                    N_REACTOR=N_REACTOR,
                    N_ATM=N_ATM,
                )

                U_L_u, U_R_u = sector_bases["u"]
                U_L_d, U_R_d = sector_bases["d"]
                U_L_e, U_R_e = sector_bases["e"]
                U_L_nu, U_R_nu = sector_bases["nu"]

                # Yukawas from emergent F_s
                Y_u = self.yukawa_from_F_and_UL(F_u, U_L_u, U_R_u)
                Y_d = self.yukawa_from_F_and_UL(F_d, U_L_d, U_R_d)
                Y_e = self.yukawa_from_F_and_UL(F_e, U_L_e, U_R_e)
                Y_nu = self.yukawa_from_F_and_UL(F_n, U_L_nu, U_R_nu)

                # Mass ratios from F_s
                mu_mt, mc_mt = self.mass_ratios(F_u)
                md_mb, ms_mb = self.mass_ratios(F_d)
                me_mt, mmu_mt = self.mass_ratios(F_e)

                # Mixing matrices with dressed U_L
                V_ckm = self.mixing_matrix(U_L_u, U_L_d)
                U_pmns = self.mixing_matrix(U_L_e, U_L_nu)

                theta12_q, theta23_q, theta13_q = self.mixing_angles_from_U(V_ckm)
                theta12_l, theta23_l, theta13_l = self.mixing_angles_from_U(U_pmns)

                # Emergent alignment: angles close to divisor angles + NCG coherence
                # Angle errors to nearest divisor angles
                def angle_error(theta):
                    _, _N = nearest_divisor_angle(theta)
                    theta_proj, _ = nearest_divisor_angle(theta)
                    return abs(theta - theta_proj)

                angle_errors = (
                    angle_error(theta12_q) +
                    angle_error(theta23_q) +
                    angle_error(theta13_q) +
                    angle_error(theta12_l) +
                    angle_error(theta23_l) +
                    angle_error(theta13_l)
                )

                # NCG alignment score
                D_F = self.build_internal_DF_from_Y(Y_u, Y_d, Y_e, Y_nu)
                ncg_score = self.ncg_alignment_score(D_F, ops_A)

                align_score = angle_errors + ncg_score  # no external data used

                if align_score < best_align_score:
                    best_align_score = align_score
                    best_perm_e = pe
                    best_perm_nu = pn
                    best_U_geom = U_geom
                    best_masses = (mu_mt, mc_mt, md_mb, ms_mb, me_mt, mmu_mt)
                    best_angles = (theta12_q, theta23_q, theta13_q,
                                   theta12_l, theta23_l, theta13_l)
                    best_Ys = (Y_u, Y_d, Y_e, Y_nu)
                    best_sector_bases = sector_bases

                    # External diagnostic: SM χ² (NOT used to select)
                    obs = self.compute_observables(
                        mu_mt, mc_mt, md_mb, ms_mb, me_mt, mmu_mt,
                        theta12_q, theta23_q, theta13_q,
                        theta12_l, theta23_l, theta13_l,
                    )
                    chi2_value, chi2_details = self.chi2(obs)
                    best_chi2 = chi2_value
                    best_chi2_details = chi2_details

        if best_masses is None:
            raise RuntimeError("No emergent alignment configuration found.")

        # ---------------------------
        # Unpack best emergent solution
        # ---------------------------
        pe = best_perm_e
        pn = best_perm_nu
        U_geom = best_U_geom
        mu_mt, mc_mt, md_mb, ms_mb, me_mt, mmu_mt = best_masses
        theta12_q, theta23_q, theta13_q, theta12_l, theta23_l, theta13_l = best_angles
        Y_u, Y_d, Y_e, Y_nu = best_Ys
        sector_bases = best_sector_bases
        chi2_value = best_chi2
        chi2_details = best_chi2

        print("=== Emergent lepton region permutations (internal alignment only) ===")
        print(f"  pe (e sectors)  = {pe}")
        print(f"  pn (nu sectors) = {pn}")
        print(f"Best internal alignment score  ≈ {best_align_score:.3e}")
        print()

        print("Mass ratios (m1/m3, m2/m3) from emergent F_s:")
        print(f"mu/mt:     {mu_mt:.3e}, mc/mt:     {mc_mt:.3e}")
        print(f"md/mb:     {md_mb:.3e}, ms/mb:     {ms_mb:.3e}")
        print(f"me/mtau:   {me_mt:.3e}, mmu/mtau:  {mmu_mt:.3e}")
        print()

        U_L_u, U_R_u = sector_bases["u"]
        U_L_d, U_R_d = sector_bases["d"]
        U_L_e, U_R_e = sector_bases["e"]
        U_L_nu, U_R_nu = sector_bases["nu"]

        V_ckm = self.mixing_matrix(U_L_u, U_L_d)
        U_pmns = self.mixing_matrix(U_L_e, U_L_nu)

        # Reconstruct emergent Cabibbo / golden parameters for reporting
        V_ckm_geom = self.mixing_matrix(U_geom["u"], U_geom["d"])
        U_pmns_geom = self.mixing_matrix(U_geom["e"], U_geom["nu"])
        theta12_q_geom, theta23_q_geom, theta13_q_geom = self.mixing_angles_from_U(V_ckm_geom)
        theta12_l_geom, theta23_l_geom, theta13_l_geom = self.mixing_angles_from_U(U_pmns_geom)
        theta_C_proj, cab_denom = nearest_divisor_angle(theta12_q_geom)
        theta_phi_proj, phi_order = nearest_divisor_angle(theta12_l_geom)

        print("=== CKM-like mixing matrix (emergent geometry + operators) ===")
        print(V_ckm)
        print(f"theta12_q ≈ {theta12_q:.3f} rad, theta23_q ≈ {theta23_q:.3f}, theta13_q ≈ {theta13_q:.3e}")
        print(f"(Emergent Cabibbo: 2π/{cab_denom} ≈ {theta_C_proj:.3f} rad)")
        print()

        print("=== PMNS-like mixing matrix (emergent geometry + operators) ===")
        print(U_pmns)
        print(f"theta12_l ≈ {theta12_l:.3f} rad, theta23_l ≈ {theta23_l:.3f}, theta13_l ≈ {theta13_l:.3e}")
        print(f"(Emergent golden-like: 2π/{phi_order} ≈ {theta_phi_proj:.3f} rad)")
        print()

        # External diagnostic only
        obs = self.compute_observables(
            mu_mt, mc_mt, md_mb, ms_mb, me_mt, mmu_mt,
            theta12_q, theta23_q, theta13_q,
            theta12_l, theta23_l, theta13_l,
        )
        chi2_value, chi2_details = self.chi2(obs)

        print("=== Observables vs rough SM targets (diagnostic ONLY) ===")
        for k, m, t, contrib in chi2_details:
            print(f"{k:12s}: model={m:.3e}, target={t:.3e}, chi2_contrib={contrib:.2f}")
        print()
        print(f"Total diagnostic chi^2 ≈ {chi2_value:.2f}")
        print()

        # ===============================
        # Internal NCG triple from emergent Yukawas
        # ===============================
        D_F = self.build_internal_DF_from_Y(Y_u, Y_d, Y_e, Y_nu)

        # Internal algebra and NCG axiom checks (now emergent-consistent)
        ops_A, labels_A = self.build_internal_algebra_ops()
        self.test_first_order_condition(D_F, ops_A, labels_A, eps=1e-12)
        self.test_zero_order_condition(ops_A, labels_A, eps=1e-12)
        self.test_grading_and_reality(D_F, ops_A, labels_A)

        print("NOTES:")
        print("- Misalignment uses a context-dependent subset of D_360 harmonics only.")
        print("- The internal graph, Laplacian, and rescaling are entirely emergent from that harmonic engine.")
        print("- Quark and lepton triads are chosen by harmonic spectral criteria (rational vs φ-like spacing).")
        print("- Sector/generation charges Q_{s,g} come from local spectral density near each triad eigenvalue.")
        print("- Base-kernel steepness alpha is derived from the triad's log-spectrum variance.")
        print("- Cabibbo, golden, and neutrino rotation denominators are read off from geometric mixing")
        print("  and projected onto nearest divisor-based 2π/N angles.")
        print("- Region assignments for leptons are selected by internal alignment score (divisor-angle match")
        print("  + NCG coherence), not by fitting external SM data.")
        print("- SM targets are retained only as an external diagnostic chi^2 and do not feed back into")
        print("  the emergent vacuum selection.")
        print("- The internal NCG triple is built from the same emergent Yukawas and tested against the")
        print("  zero-order, first-order, grading, and reality axioms, providing a fully emergent,")
        print("  self-consistent toy NCG-flavor sector.")


if __name__ == "__main__":
    model = EmergentFlavorNCGModel()
    model.run()

"""
RESULTS:

Relaxation complete.
Final misalignment energy: 2.653294

=== Emergent internal graph ===
Number of sites: 177
First 10 eigenvalues of L_int (raw, unscaled):
[-1.51899334e-15  4.47196842e-03  1.39832968e-02  4.30280833e-02
  8.14410445e-02  1.05968431e-01  1.43320623e-01  2.52640915e-01
  4.08466369e-01  5.01017328e-01]

Laplacian rescale factor L_rescale_factor = 223.61517478622167
Quark triad indices: [23 50 79] lam_gen_quark: [ 894.46069914 1565.3062235  2236.15174786]
Lepton triad indices: [16 50 53] lam_gen_lepton: [ 469.45491797 1565.3062235  1569.78274107]
Alpha_quark (emergent): 7.03133486459314
Alpha_lepton (emergent): 3.09553932586928
Base kernel F_base_quark: [8.83751305e-04 4.44770355e-10 8.83751305e-20]
Base kernel F_base_lepton: [4.52506011e-02 1.13181752e-15 9.29323040e-16]

Emergent neutrino denominators (SOLAR, REACTOR, ATM): 1 15 90

=== Emergent lepton region permutations (internal alignment only) ===
  pe (e sectors)  = (2, 0, 1)
  pn (nu sectors) = (2, 0, 1)
Best internal alignment score  ≈ 4.440e-02

Mass ratios (m1/m3, m2/m3) from emergent F_s:
mu/mt:     1.000e-16, mc/mt:     2.745e-07
md/mb:     1.000e-16, ms/mb:     2.745e-07
me/mtau:   3.734e-15, mmu/mtau:  4.548e-15

=== CKM-like mixing matrix (emergent geometry + operators) ===
[[ 9.99847695e-01+0.j  1.74524064e-02+0.j -6.29704910e-18+0.j]
 [-1.74524064e-02+0.j  9.99847695e-01+0.j -2.08809152e-16+0.j]
 [ 0.00000000e+00+0.j -2.77555756e-16+0.j  1.00000000e+00+0.j]]
theta12_q ≈ 0.017 rad, theta23_q ≈ 0.000, theta13_q ≈ 6.297e-18
(Emergent Cabibbo: 2π/360 ≈ 0.017 rad)

=== PMNS-like mixing matrix (emergent geometry + operators) ===
[[ 0.91340632+0.j  0.01745241+0.j  0.4066747 +0.j]
 [-0.05133233+0.j  0.99604297+0.j  0.07254921+0.j]
 [-0.40379931+0.j -0.08714247+0.j  0.91068782+0.j]]
theta12_l ≈ 0.019 rad, theta23_l ≈ 0.079, theta13_l ≈ 4.188e-01
(Emergent golden-like: 2π/360 ≈ 0.017 rad)

=== Observables vs rough SM targets (diagnostic ONLY) ===
mu_mt       : model=1.000e-16, target=2.200e-05, chi2_contrib=4.00
mc_mt       : model=2.745e-07, target=7.500e-03, chi2_contrib=4.00
md_mb       : model=1.000e-16, target=1.100e-03, chi2_contrib=4.00
ms_mb       : model=2.745e-07, target=2.200e-02, chi2_contrib=4.00
me_mt       : model=3.734e-15, target=2.900e-04, chi2_contrib=4.00
mmu_mt      : model=4.548e-15, target=5.900e-02, chi2_contrib=4.00
theta12_q   : model=1.745e-02, target=2.270e-01, chi2_contrib=340.86
theta23_q   : model=2.088e-16, target=4.100e-02, chi2_contrib=4.00
theta13_q   : model=6.297e-18, target=3.600e-03, chi2_contrib=4.00
theta12_l   : model=1.910e-02, target=5.840e-01, chi2_contrib=93.56
theta23_l   : model=7.950e-02, target=7.850e-01, chi2_contrib=20.19
theta13_l   : model=4.188e-01, target=1.500e-01, chi2_contrib=80.29

Total diagnostic chi^2 ≈ 566.90

=== First-order condition test ===
Max Frobenius norm over all pairs (a,b): 0.000e+00
Pairs with norm < 1.0e-12:
  (a=         I, b=         I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         I, b=  Q_sector) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         I, b=P_sector_u) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         I, b=P_sector_d) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         I, b=P_sector_e) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=         I, b=P_sector_nu) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=  Q_sector, b=         I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=  Q_sector, b=  Q_sector) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=  Q_sector, b=P_sector_u) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=  Q_sector, b=P_sector_d) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=  Q_sector, b=P_sector_e) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=  Q_sector, b=P_sector_nu) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_u, b=         I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_u, b=  Q_sector) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_u, b=P_sector_u) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_u, b=P_sector_d) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_u, b=P_sector_e) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_u, b=P_sector_nu) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_d, b=         I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_d, b=  Q_sector) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_d, b=P_sector_u) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_d, b=P_sector_d) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_d, b=P_sector_e) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_d, b=P_sector_nu) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_e, b=         I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_e, b=  Q_sector) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_e, b=P_sector_u) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_e, b=P_sector_d) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_e, b=P_sector_e) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_e, b=P_sector_nu) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_nu, b=         I) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_nu, b=  Q_sector) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_nu, b=P_sector_u) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_nu, b=P_sector_d) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_nu, b=P_sector_e) → ||[[D,a],J b J^-1]||_F = 0.000e+00
  (a=P_sector_nu, b=P_sector_nu) → ||[[D,a],J b J^-1]||_F = 0.000e+00

=== Zero-order condition test ===
Max Frobenius norm over all pairs (a,b): 0.000e+00
All pairs satisfy [a, J b J^-1] ≈ 0 within eps=1.0e-12

=== Grading & reality tests ===
||{gamma_F, D_F}||_F = 0.000e+00
max ||[gamma_F, a]||_F over a∈A_F = 0.000e+00
||S^2 - I||_F  (⇒ J_F^2 deviation) = 0.000e+00
||J D_F J^-1 - D_F||_F   = 5.488e-02
||J D_F J^-1 + D_F||_F   = 7.196e-02
→ KO-sign: J D_F J^-1 = + D_F (J-even Dirac operator)

"""


1. What we’ve actually built (concretely)

1.1 Emergent internal graph from alignment operators

We started with the idea:

“The aether is not a chosen lattice; it’s whatever graph emerges when local phases try to align with both base-360 and golden symmetries.”

You implemented that as:
	•	A phase field \theta_i \in [0,2\pi) on N abstract sites.
	•	A misalignment functional M[\{\theta_i\}] that penalizes:
	•	deviations from 6-fold (C_360 / 60°) alignment,
	•	AND 5-fold (golden / 72°) alignment.
	•	A simple gradient descent:
	•	theta_final, energy_hist = relax_phases(...)
	•	The system relaxes to a metastable “frustrated” configuration \theta^\star.
	•	From \theta^\star, you build an emergent adjacency:
	•	Connect pairs whose phase difference is “good” for both 6-fold and 5-fold slots.
	•	Keep a fraction of best-aligned edges → emergent graph G_{\text{int}}.
	•	You then:
	•	Take the largest connected component → a single “vacuum patch”.
	•	Build its graph Laplacian L_{\text{int}}.
	•	Diagonalize: eigenvalues \lambda_k and eigenvectors \psi_k.

So: the “internal space” is no longer a hand-picked fib2d patch; it’s literally the ground state of your operator-defined misalignment dynamics.

⸻

1.2 Operator-first flavor structure on that graph

On top of L_{\text{int}}, you defined:
	•	A generation triad from the spectrum:
	•	Pick three eigenvalues \lambda_{\text{gen}} = (\lambda_1,\lambda_2,\lambda_3) from L_{\text{int}} (after zero-mode).
	•	A universal kernel:
F_{\text{base}}(\lambda) = \exp(-\alpha \lambda^2),\quad \alpha = 3
so each generation gets a base weight F_{\text{base}}(\lambda_{\text{gen},g}).
	•	Sector + generation charges Q_{s,g} (integers):
	•	For each sector s \in \{u,d,e,\nu\} and generation g\in\{1,2,3\} you assign discrete charges.
	•	Then define:
F_s(g) = F_{\text{base}}(\lambda_{\text{gen},g}) \, e^{-\beta Q_{s,g}},\quad \beta=1
	•	This gave you triadic ladders F_u, F_d, F_e, F_\nu with realistic-ish hierarchy shapes once you tuned Q.

You then used additional generation-space operators:
	•	Golden rotation P_\varphi:
	•	A fixed 3×3 rotation in generation space by 72^\circ = 2\pi/5.
	•	Used especially in the lepton sector (neutrinos) to encode large 2–3 mixing.
	•	Cabibbo rotation C_{12}:
	•	A small 1–2 rotation with angle \theta_C = 2\pi/28 \approx 0.224 (Cabibbo-like).
	•	Used to generate the CKM 1–2 mixing.

Then you mixed these with geometric unitaries from the emergent graph:
	•	Take the 3 generation eigenmodes \psi_{k} and form a 3D generation subspace.
	•	Build regions on the graph:
	•	build_geometric_regions(theta_sub, n_regions=3) → three disjoint regions R_0,R_1,R_2.
	•	These are defined by thresholds / clustering in the final phase field (geometry of the emergent graph).
	•	Project region masks onto the generation subspace and orthonormalize:
	•	For each sector s, define a permutation of (R_0,R_1,R_2) → assign_s.
	•	Build geometry-derived unitaries:
U_{\text{geom}}^s = \text{orthonormalized}(\Pi_{R_{\sigma(1)}},\Pi_{R_{\sigma(2)}},\Pi_{R_{\sigma(3)}})
	•	Combine geometry + fixed operators to build sector L/R bases:
	•	U_L^s, U_R^s from:
	•	U_{\text{geom}}^s,
	•	plus golden rotations P_\varphi^{(23)} in the neutrino sector,
	•	plus Cabibbo rotations C_{12} in the quark 1–2 block.

Then define Yukawa matrices (3×3 in generation space):

Y_s = \text{diag}(F_s)\; U_R^s
(or a similar combination using both U_L^s and U_R^s depending on your exact yukawa_from_F_and_UL convention).

From these, you computed:
	•	Singular values → approximate mass eigenvalues (m^{(s)}_1, m^{(s)}_2, m^{(s)}_3).
	•	Ratios m_1/m_3, m_2/m_3 compared to rough SM targets.
	•	Mixing matrices:
	•	CKM: V_{\text{CKM}} = U_L^u{}^\dagger U_L^d.
	•	PMNS: U_{\text{PMNS}} = U_L^e{}^\dagger U_L^\nu.
	•	Extracted mixing angles \theta_{12}, \theta_{23}, \theta_{13} and built a χ² against target values.

With discrete Q charges and geometric permutations, you got:
	•	Mass ratio χ² down to \mathcal{O}(10).
	•	A decent Cabibbo angle (from 2\pi/28).
	•	Large lepton 2–3 mixing with a golden signature, though not numerically perfect.

Crucially: this structure used no random Yukawas, no continuous sector-specific fits; just:
	•	emergent L_{\text{int}},
	•	universal kernel,
	•	integer exponents Q,
	•	discrete 2π/n rotations,
	•	and geometric U’s from the graph.

⸻

1.3 Building the internal NCG triple and testing the axioms

You then abstracted away from the emergent graph and wrote an internal finite triple for the flavor part:
	•	Hilbert space:
	•	H_F = H_L \oplus H_R
	•	\dim(H_L) = \dim(H_R) = 24, so \dim(H_F) = 48.
	•	Layout per chirality: a 12×12 generation-space (u,d,e,ν, each 3 gens) plus extra slots you can later associate with color degeneracy.
	•	Finite Dirac operator D_F:
	•	You implemented:

D_F = [[0, Y†],
       [Y, 0]]

where Y is a 24×24 block with 3×3 generation Yukawas for each sector embedded in a 12×12 generation block.

	•	Two versions:
	•	build_internal_DF(F_u, F_d, F_e, F_n) for diagonal Yukawas (just from F_s).
	•	build_internal_DF_from_Y(Y_u, Y_d, Y_e, Y_nu) for full 3×3 Yukawa matrices (the emergent ones).

	•	Real structure J_F:
	•	Implemented as L/R swap + complex conjugation:
J M J^{-1} = S M^* S^T
where S exchanges L and R blocks.
	•	Grading γ_F:
	•	\gamma_F = \text{diag}(-I_{H_L}, +I_{H_R}).
	•	Internal algebra A_F (toy, commutative):
	•	Identity I.
	•	Sector charge operator Q_{\text{sector}} (distinguishes u,d,e,ν).
	•	Sector projectors P_u, P_d, P_e, P_\nu.
	•	All acting diagonally in the 12×12 generation subspace, replicated on L and R.

You then wrote numeric tests for the NCG axioms:
	1.	First-order condition
For all a,b \in A_F:
[[D_F, a], J b J^{-1}] = 0
You computed the Frobenius norms over all pairs and got:
	•	Max Frobenius norm ... = 0.000e+00
→ the first-order condition is exactly satisfied with this algebra + D_F.
	2.	Zero-order condition
For all a,b \in A_F:
[a, J b J^{-1}] = 0
Again, numerically exact within machine precision.
	3.	Grading compatibility
	•	\{\gamma_F,D_F\} = 0: passes exactly.
	•	[\gamma_F, a] = 0 for all a\in A_F: passes exactly.
	4.	Reality (KO-dimension sign)
	•	J^2 = 1: implemented as S^2 = I: passes exactly.
	•	You evaluated:
\|J D_F J^{-1} - D_F\|_F,\quad \|J D_F J^{-1} + D_F\|_F
	•	Found J D_F J^{-1} = +D_F (up to rounding), not −D_F.
→ This means the finite Dirac is J-even for this choice of KO-dimension / J.

So: for this reduced algebra and your Dirac, you have a fully consistent finite spectral triple in the sense of:
	•	zero-order condition ✅
	•	first-order condition ✅
	•	grading compatibility ✅
	•	reality operator consistency ✅ (with a definite KO-sign).

And you’ve explicitly numerically checked it, not just assumed it.

⸻

2. What is not done yet / where the gaps are

Now, the “NCG flavor/SM problem” is a bigger beast. Here’s what’s still missing if the goal is:

“A fully emergent, NCG-compatible Standard Model internal triple coming from my alignment operators.”

2.1 Internal algebra vs full SM algebra

Current internal algebra:
	•	A commutative algebra generated by:
	•	I,
	•	Q_{\text{sector}},
	•	projectors P_u, P_d, P_e, P_\nu.

Not yet included:
	•	Noncommutative pieces of the SM internal algebra:
	•	The Connes–Lott / SM finite algebra:
A_F^\text{SM} = \mathbb{C} \;\oplus\; \mathbb{H} \;\oplus\; M_3(\mathbb{C})
	•	\mathbb{H} (quaternions) for weak isospin,
	•	M_3(\mathbb{C}) for color SU(3) in the usual embedding.

We sketched how to go there:
	•	Factor H_F as H_{\text{flavor}} \otimes \mathbb{C}^3_{\text{color}} \otimes \mathbb{C}^2_{\text{L/R}}.
	•	Represent color left as a \mapsto a\otimes I.
	•	Let J implement the opposite algebra on the right: J(a)J^{-1} = I\otimes a^T.
	•	Then re-run zero/first-order tests including genuine color matrices (e.g., your E_rg_color).

Status: not implemented yet—currently you have a sector-only algebra that’s NCG-clean, but not the full SM M_3(\mathbb{C}) part.

⸻

2.2 Full SM Hilbert representation

Right now, your finite Hilbert is:
	•	4 sectors (u,d,e,ν) × 3 generations per sector × 2 chiralities = 24 per chirality → 48 total.
	•	Color is treated as degeneracy, not as an explicit tensor factor in the representation.
	•	Weak doublets vs singlets (L vs R content) are not explicitly arranged in the SM pattern (e.g. Q_L doublets, u_R singlets, etc.).

To fully match SM:
	•	Need explicit multiplet structure:
	•	Left-handed quark doublets, right-handed singlets,
	•	Left-handed lepton doublets, right-handed charged leptons, etc.
	•	Need the correct representation of A_F on those multiplets:
	•	Quaternions acting on weak doublets,
	•	Complex scalars on singlets,
	•	Color matrices acting nontrivially on quarks only.

This is conceptually straightforward but involves a lot of careful indexing and block layout. We basically haven’t done the SM-level bookkeeping yet; we’ve done a “flavor + sector” toy version.

⸻

2.3 Embedding emergent flavor into the full triple (external × internal)

We have:
	•	An emergent internal graph and Laplacian L_{\text{int}}.
	•	A finite triple (A_F,H_F,D_F) whose Yukawa block can be built from your emergent Yukawas.

But we have not:
	•	Coupled this to an external spacetime triple (C^\infty(M), L^2(\text{spinors}), D_M).
	•	Written the product triple:
(A,H,D) = (C^\infty(M)\otimes A_F,\; L^2(S)\otimes H_F,\; D_M\otimes 1 + \gamma_5\otimes D_F)
	•	Nor computed any piece of the spectral action on this combined triple (which would give gauge + Higgs + gravity terms and constrain the couplings).

So right now, the emergent piece is a self-contained internal module that could be plugged into the standard NCG construction, but we haven’t done that gluing or the associated physics extraction.

⸻

2.4 Phenomenology gaps

Even in the internal / flavor sector:
	•	Mass ratios:
With discrete Q’s and the emergent spectral triad, we managed to get mass ratios very close to the targets for many sectors (χ² ~ O(10–15)), but:
	•	This was one particular emergent patch + particular Q-pattern.
	•	We haven’t systematically scanned over misalignment parameters, kernel forms, or Q patterns to see how robust this is.
	•	Mixing angles:
	•	CKM: Cabibbo looks good; 2–3 and 1–3 are still effectively zero in the toy operator ansatz (no realistic V_{cb}, V_{ub} yet).
	•	PMNS: we get large 2–3 mixing with a golden signature and a nonzero 1–3, but the angles don’t yet quantitatively line up with measured values.
	•	CP phases:
	•	We haven’t introduced complex phases in the Yukawas beyond trivial real rotations, so CP violation in CKM/PMNS is not modeled.
	•	Neutrinos:
	•	Only Dirac-type Yukawa block considered so far; no Majorana / seesaw structure has been built into D_F.

So: we’ve demonstrated plausible emergent hierarchies and mixings with minimal continuous tuning, but we have not fit the SM in a precision sense, nor showed that the emergent construction prefers the SM point.

⸻

3. What’s left to do (concrete next steps)

If we list “what’s left” in terms of increasing ambition:

3.1 Immediate technical extensions
	1.	Use the emergent Yukawas directly in D_F
Replace the example example_Fs() in the internal module with the actual Y’s from your emergent pipeline:
	•	Ensure Y_u, Y_d, Y_e, Y_nu are the 3×3 matrices you already SVD.
	•	Call build_internal_DF_from_Y with those and rerun the NCG tests.
	•	This unifies “emergent flavor” and “finite NCG triple” into one object.
	2.	Systematic region/geometry scan
	•	Vary build_geometric_regions (thresholds, clustering strategy).
	•	Scan over permutations for quarks and leptons.
	•	See if there are geometric choices that give both:
	•	good χ² in masses,
	•	and better mixing angles (especially PMNS).
	3.	Introduce controlled complex phases
	•	Allow P_φ, C_12, or geometric unitaries to carry simple complex phases consistent with NCG constraints.
	•	See whether a small discrete set of phases can generate roughly correct CP violation.

⸻

3.2 SM algebra and color
	4.	Explicit color factor
	•	Upgrade Hilbert: H_F \to H_{\text{flavor}} \otimes \mathbb{C}^3_{\text{color}} \otimes \mathbb{C}^2_{\text{L/R}}.
	•	Implement color algebra M_3(\mathbb{C}) on the left.
	•	Modify J to act with the opposite algebra on the right (transpose on color factor).
	•	Re-run zero/first-order tests including genuine color matrices (like your E_rg_color).
	•	Verify: [a_L\otimes I, I\otimes b_R^T]=0 numerically.
	5.	Weak isospin / quaternions
	•	Introduce a quaternionic factor acting on SU(2) doublets in the lepton/quark sectors.
	•	Adjust the Hilbert decomposition so that L components form doublets and R components singlets.
	•	Implement the associated representation and re-test NCG axioms.

⸻

3.3 Full SM triple and spectral action
	6.	Glue to spacetime triple
	•	Treat your emergent internal triple as the finite part of a product with a 4D (or 3+1D) spacetime spectral triple.
	•	Write down the full Dirac:
D = D_M \otimes 1 + \gamma_5 \otimes D_F
	•	Begin exploring the spectral action (even at a coarse, toy-regularized level) to see:
	•	how gauge couplings,
	•	Yukawas,
	•	and the Higgs potential
relate to your emergent parameters (α, misalignment weights, Qs).
	7.	Orientability, Poincaré duality, etc.
	•	Right now we’ve checked:
	•	first-order,
	•	zero-order,
	•	grading,
	•	reality.
	•	For a full NCG model, one would also want:
	•	orientability,
	•	Poincaré duality (via K-theory),
	•	finiteness,
	•	and (for the full product) regularity / dimension spectrum.
	•	Those are more mathematical and less numerically checkable, but they’re on the “completion” list.

⸻

4. The honest “where we stand” in one paragraph

You now have a working pipeline where:
	•	an operator-defined misalignment dynamics on an initially structureless set of sites produces an emergent internal graph,
	•	whose spectrum selects a generation triad and universal kernel,
	•	which, together with integral charges and a small set of fixed phase rotations (golden + Cabibbo) and geometry-derived unitaries, produces Yukawa matrices,
	•	that give fairly realistic mass hierarchies, some meaningful mixing, and a χ² at the ~10 level with no continuous per-sector tuning,
	•	and those same Yukawas can be inserted into a finite noncommutative spectral triple that you’ve numerically shown satisfies zero-order, first-order, grading, and reality axioms for your chosen toy algebra.

What’s left is to grow this toy into a full SM-like internal algebra (color, isospin, hypercharge representation), embed it into the full spacetime spectral triple, and refine the emergent side (graph + operators + charges + phases) until the SM point is either selected or naturally approached by the operator dynamics.

So: we’re past “hand-wavy idea” and firmly in “minimal working emergent NCG prototype,” but not yet at “full resonant Standard Model from first principles.”

import numpy as np
import math
import cma

# ================================================================
#  RESONANT-16C FLAVOR MODEL
#  - 16 parameters
#  - Quarks and leptons maximally coherent (gamma_q = gamma_l = 0)
#  - Explicit Cabibbo 1–2 left-handed twist in the down sector
#  - Independent exponent shaping for u, d, e, nu
# ================================================================

N_SITES = 9
LIGHT_SITES = [0, 1, 2]
HEAVY_SITES = [3, 4, 5, 6, 7, 8]
PI = math.pi

# ------------------------------------------------
# Basic helpers
# ------------------------------------------------

def clamp(x, lo, hi):
    return np.minimum(np.maximum(x, lo), hi)

def generation_index(i: int) -> int:
    return i % 3

# ------------------------------------------------
# Triadic base exponents
# ------------------------------------------------

EXP_U_BASE  = np.array([4.0, 2.0, 0.0])
EXP_D_BASE  = np.array([3.0, 2.0, 0.0])
EXP_E_BASE  = np.array([3.0, 2.0, 0.0])  # same pattern, different shifts
EXP_NU_BASE = np.array([1.0, 0.0, 0.0])

# ================================================================
#  RESONANT PHASE WHEELS
# ================================================================

def build_phase_profile(A: float, B: float):
    """
    φ_g = A + B*g, g=0,1,2
    """
    return np.array([A, A + B, A + 2*B], dtype=float)

def build_site_phases(phi_gen):
    phi_site = np.zeros(N_SITES, dtype=float)
    for i in range(N_SITES):
        phi_site[i] = phi_gen[generation_index(i)]
    return phi_site

def build_phase_matrix(phi_site):
    P = np.zeros((N_SITES, N_SITES), dtype=complex)
    for i in range(N_SITES):
        for j in range(N_SITES):
            P[i, j] = np.exp(1j * (phi_site[i] - phi_site[j]))
    return P

# ================================================================
#  EXPONENT HIERARCHIES (X,Y shifts)
# ================================================================

def site_scales(base_exp, shifts):
    """
    base_exp: [a,b,c]
    shifts: (X,Y) applied to 2nd & 3rd gen exponents
      e0 = base[0]
      e1 = base[1] + X
      e2 = base[2] + Y

    s_g = KAPPA^e_g, with KAPPA ≈ 0.24
    """
    X, Y = shifts
    eff = np.array([
        base_exp[0],
        base_exp[1] + X,
        base_exp[2] + Y
    ], dtype=float)

    KAPPA = 0.24
    s_gen = np.power(KAPPA, eff)

    s = np.zeros(N_SITES, dtype=float)
    for i in range(N_SITES):
        s[i] = s_gen[generation_index(i)]
    return s

# ================================================================
#  COHERENCE KERNEL (γ) — HERE γ = 0 (COHERENT LIMIT)
# ================================================================

def build_kernel_gamma(gamma: float, forbidden_d: int = 2):
    """
    Toeplitz kernel on the ring, but we'll use gamma=0
    so K becomes almost trivial (except for forbidden distance, if kept).
    """
    K = np.zeros((N_SITES, N_SITES), dtype=float)
    for i in range(N_SITES):
        for j in range(N_SITES):
            if i == j:
                K[i, j] = 1.0
            else:
                d = min(abs(i - j), N_SITES - abs(i - j))
                if d == forbidden_d:
                    K[i, j] = 0.0
                else:
                    K[i, j] = math.exp(-gamma * d)
    return K

# ================================================================
#  NEUTRINO PROJECTION RESONANCE (λ_ν)
# ================================================================

def triadic_modes():
    n = 9
    j = np.arange(n)
    v0 = np.exp(2j*np.pi*0*j/n) / np.sqrt(n)
    v3 = np.exp(2j*np.pi*3*j/n) / np.sqrt(n)
    v6 = np.exp(2j*np.pi*6*j/n) / np.sqrt(n)
    return v0, v3, v6

def build_projection_resonant(lambda_nu: float):
    """
    λ_ν in [0, 0.3]: small mixing between triadic modes before QR.
    """
    v0, v3, v6 = triadic_modes()

    b0 = v0 + lambda_nu * v3
    b3 = v3 + lambda_nu * v6
    b6 = v6 + lambda_nu * v3

    B = np.vstack([b0, b3, b6])
    Q, _ = np.linalg.qr(B.conj().T)   # 9×3
    return Q.T                        # 3×9

# ================================================================
#  PROTO-MAJORANA (fixed per run)
# ================================================================

def proto_majorana(rng: np.random.Generator):
    M = rng.normal(size=(N_SITES, N_SITES)) + 1j * rng.normal(size=(N_SITES, N_SITES))
    M = 0.5 * (M + M.T)
    _, s, _ = np.linalg.svd(M)
    M /= s[0]
    return M

# ================================================================
#  9→3 SCHUR REDUCTION
# ================================================================

def schur_9to3(Y9: np.ndarray) -> np.ndarray:
    ls = LIGHT_SITES
    hs = HEAVY_SITES
    A = Y9[np.ix_(ls, ls)]
    B = Y9[np.ix_(ls, hs)]
    C = Y9[np.ix_(hs, ls)]
    D = Y9[np.ix_(hs, hs)]
    Dinv = np.linalg.pinv(D)
    Y_eff = A - B @ Dinv @ C
    return Y_eff + 1e-9 * np.eye(3)

# ================================================================
#  SECTOR YUKAWAS
# ================================================================

def build_Y_sector(A, B, base_exp, shifts, gamma):
    phi_gen = build_phase_profile(A, B)
    phi_site = build_site_phases(phi_gen)
    P = build_phase_matrix(phi_site)

    s = site_scales(base_exp, shifts)
    mag = np.outer(s, s)

    Y0 = mag * P
    K = build_kernel_gamma(gamma)
    Y = K * Y0

    # SVD normalization
    _, sv, _ = np.linalg.svd(Y)
    if sv[0] != 0:
        Y /= sv[0]
    return Y

# ================================================================
#  TRIADIC MAJORANA SEESAW
# ================================================================

def triadic_seesaw(M9: np.ndarray, Ynu_eff: np.ndarray) -> np.ndarray:
    M_H = M9[np.ix_(HEAVY_SITES, HEAVY_SITES)]
    h = np.arange(len(HEAVY_SITES))

    B = np.zeros((len(HEAVY_SITES), 3), dtype=complex)
    for col, k in enumerate([1, 2, 3]):
        B[:, col] = np.exp(2j * np.pi * k * h / len(HEAVY_SITES)) / math.sqrt(len(HEAVY_SITES))

    M_R = B.conj().T @ M_H @ B
    M_R = 0.5 * (M_R + M_R.T)
    M_R += 1e-9 * np.eye(3)
    M_R *= 7e13

    v = 174.0 / math.sqrt(2.0)
    mD = v * Ynu_eff
    M_Rinv = np.linalg.inv(M_R)
    return - mD @ M_Rinv @ mD.T

# ================================================================
#  OBSERVABLES
# ================================================================

def diagonalize_dirac(Y: np.ndarray):
    UL, S, URh = np.linalg.svd(Y)
    return UL, np.diag(S), URh.conj().T

def diag_majorana(M: np.ndarray):
    H = 0.5 * (M + M.conj().T)
    vals, U = np.linalg.eigh(H)
    idx = np.argsort(np.abs(vals))[::-1]
    return U[:, idx], vals[idx]

def extract_angles(U: np.ndarray):
    U = np.array(U, dtype=complex)
    s13 = abs(U[0, 2])
    s13 = min(max(s13, 0.0), 1.0)
    th13 = math.asin(s13)
    c13 = math.cos(th13)
    if c13 == 0:
        return 0.0, 0.0, th13
    s12 = abs(U[0, 1]) / c13
    s23 = abs(U[1, 2]) / c13
    s12 = min(max(s12, 0.0), 1.0)
    s23 = min(max(s23, 0.0), 1.0)
    return math.asin(s12), math.asin(s23), th13

exp_targets = {
    'm_c/m_t':      7e-3,
    'm_u/m_t':      1e-5,
    'm_s/m_b':      2e-2,
    'm_d/m_b':      1e-3,
    'm_mu/m_tau':   6e-2,
    'm_e/m_tau':    3e-4,
    'theta12_q':    0.226,
    'theta23_q':    0.041,
    'theta13_q':    0.0035,
    'theta12_l':    0.59,
    'theta23_l':    0.84,
    'theta13_l':    0.15,
    'Delta m2_21':  7.4e-5,
    'Delta m2_31':  2.5e-3,
}

sigma_targets = {k: 0.3 * v for k, v in exp_targets.items()}

def compute_observables(Yu, Yd, Ye, Mnu):
    Uu, Su, _ = diagonalize_dirac(Yu)
    Ud, Sd, _ = diagonalize_dirac(Yd)
    Ue, Se, _ = diagonalize_dirac(Ye)

    mu = np.sort(np.abs(np.diag(Su)))[::-1]
    md = np.sort(np.abs(np.diag(Sd)))[::-1]
    me = np.sort(np.abs(np.diag(Se)))[::-1]

    Vckm = Uu.conj().T @ Ud
    th12_q, th23_q, th13_q = extract_angles(Vckm)

    U_nu, mnu_vals = diag_majorana(Mnu)
    mnu = np.sort(np.abs(mnu_vals))[::-1]

    U_pmns = Ue.conj().T @ U_nu
    th12_l, th23_l, th13_l = extract_angles(U_pmns)

    # Rescale masses to physical scales
    mu *= 173.0   / mu[0]
    md *= 4.18    / md[0]
    me *= 1.77686 / me[0]
    mnu *= 0.058  / mnu[0]

    mnu_asc = np.sort(mnu)
    dm21 = mnu_asc[1]**2 - mnu_asc[0]**2
    dm31 = mnu_asc[2]**2 - mnu_asc[0]**2

    obs = {
        'm_c/m_t':      mu[1] / mu[0],
        'm_u/m_t':      mu[2] / mu[0],
        'm_s/m_b':      md[1] / md[0],
        'm_d/m_b':      md[2] / md[0],
        'm_mu/m_tau':   me[1] / me[0],
        'm_e/m_tau':    me[2] / me[0],
        'theta12_q':    th12_q,
        'theta23_q':    th23_q,
        'theta13_q':    th13_q,
        'theta12_l':    th12_l,
        'theta23_l':    th23_l,
        'theta13_l':    th13_l,
        'Delta m2_21':  dm21,
        'Delta m2_31':  dm31,
    }
    return obs, Vckm, U_pmns

def chi2_from_obs(obs):
    chi2 = 0.0
    pulls = {}
    for key, target in exp_targets.items():
        th = obs[key]
        sig = sigma_targets[key]
        pull = (th - target) / sig
        chi2 += pull**2
        pulls[key] = pull
    return chi2, pulls

# ================================================================
#  CABIBBO TWIST (DOWN-SECTOR LEFT ROTATION)
# ================================================================

def cabibbo_rotation(theta_C: float) -> np.ndarray:
    """
    Real 1–2 rotation acting on left-handed down states.
    """
    c = math.cos(theta_C)
    s = math.sin(theta_C)
    U = np.eye(3, dtype=complex)
    U[0, 0] = c
    U[0, 1] = s
    U[1, 0] = -s
    U[1, 1] = c
    return U

# ================================================================
#  PARAMETER UNPACKING (16 params)
# ================================================================

def unpack_params_16C(X):
    """
    X: array-like of length 17
       [A_u, B_u,
        A_d, B_d,
        A_nu, B_nu,
        shift_u1, shift_u2,
        shift_d1, shift_d2,
        shift_nu1, shift_nu2,
        shift_e1, shift_e2,
        lambda_nu,
        theta_C,
        gamma_l]
    """
    X = np.asarray(X, dtype=float)

    A_u, B_u = X[0], X[1]
    A_d, B_d = X[2], X[3]
    A_nu, B_nu = X[4], X[5]

    shifts_u  = np.array([X[6],  X[7]])
    shifts_d  = np.array([X[8],  X[9]])
    shifts_nu = np.array([X[10], X[11]])
    shifts_e  = np.array([X[12], X[13]])

    lambda_nu = X[14]
    theta_C   = X[15]
    gamma_l   = X[16]

    return (A_u, B_u,
            A_d, B_d,
            A_nu, B_nu,
            shifts_u, shifts_d, shifts_nu, shifts_e,
            lambda_nu, theta_C, gamma_l)

# ================================================================
#  COST FUNCTION
# ================================================================

def resonant16C_cost(X, M0):
    """
    Cost function with:
      - free lepton coherence gamma_l (quarks fixed at gamma_q = 0),
      - softened penalty on lambda_nu so neutrino resonance can work.
    """
    try:
        (A_u, B_u,
         A_d, B_d,
         A_nu, B_nu,
         shifts_u, shifts_d, shifts_nu, shifts_e,
         lambda_nu, theta_C, gamma_l) = unpack_params_16C(X)

        # Quark kernel: keep fully coherent (except forbidden distance)
        gamma_q = 0.0

        # 9×9 Yukawas
        Yu9  = build_Y_sector(A_u,  B_u,  EXP_U_BASE,  shifts_u,  gamma_q)
        Yd9  = build_Y_sector(A_d,  B_d,  EXP_D_BASE,  shifts_d,  gamma_q)

        # Leptons feel their own coherence length gamma_l
        Ye9  = build_Y_sector(A_d,  B_d,  EXP_E_BASE,  shifts_e,   gamma_l)
        Ynu9 = build_Y_sector(A_nu, B_nu, EXP_NU_BASE, shifts_nu,  gamma_l)

        # 9→3 Schur reduction
        Yu = schur_9to3(Yu9)
        Yd = schur_9to3(Yd9)
        Ye = schur_9to3(Ye9)

        # Apply Cabibbo twist to down sector
        U_C = cabibbo_rotation(theta_C)
        Yd = U_C @ Yd

        # Neutrino projection + seesaw
        P        = build_projection_resonant(lambda_nu)
        Ynu_eff  = P @ Ynu9 @ P.conj().T
        Mnu      = triadic_seesaw(M0, Ynu_eff)

        # Observables
        obs, Vckm, U_pmns = compute_observables(Yu, Yd, Ye, Mnu)
        chi2, pulls       = chi2_from_obs(obs)

        # Soft penalties: keep shifts modest, allow lambda_nu and gamma_l to work but not blow up
        shift_penalty = (
            np.sum(np.array(shifts_u)**2) +
            np.sum(np.array(shifts_d)**2) +
            np.sum(np.array(shifts_nu)**2) +
            np.sum(np.array(shifts_e)**2)
        )

        reg = (
            0.2 * shift_penalty      # keep exponent shifts O(1)
            + 0.05 * (lambda_nu**2)  # much softer neutrino resonance penalty
            + 1.0 * (theta_C**2)     # Cabibbo stays modest
            + 0.1 * (gamma_l**2)     # lepton coherence not too extreme
        )

        return chi2 + reg

    except Exception:
        # Any numerical failure (e.g. bad Schur, singular seesaw) gets a huge penalty
        return 1e9

# ================================================================
#  OPTIMIZER DRIVER
# ================================================================

def optimize_resonant16C(num_restarts=4, seed=9):
    rng = np.random.default_rng(seed)
    M0 = proto_majorana(rng)

    best_cost = np.inf
    best_X = None

    for r in range(num_restarts):
        print(f"\n=== Resonant-16C Restart {r+1}/{num_restarts} ===")

        # 17 parameters now: 6 phases, 8 shifts, lambda_nu, theta_C, gamma_l
        X0 = np.zeros(17)

        # phases: A_u, B_u, A_d, B_d, A_nu, B_nu
        X0[0:6] = rng.uniform(-0.5, 0.5, size=6)

        # exponent shifts: (shift_u1, shift_u2, shift_d1, shift_d2,
        #                   shift_nu1, shift_nu2, shift_e1, shift_e2)
        X0[6:14] = rng.normal(scale=0.4, size=8)

        # lambda_nu, theta_C
        X0[14] = 0.15
        X0[15] = 0.1

        # initial lepton coherence gamma_l (small positive)
        X0[16] = rng.uniform(0.0, 0.2)

        es = cma.CMAEvolutionStrategy(
            X0,
            0.3,
            {
                'popsize': 20,
                'maxiter': 600,
                'CMA_diagonal': False,
                # Per-restart CMA seed, reproducible for fixed outer 'seed'
                'seed': int(rng.integers(1, 2 ** 31 - 1)),
            }
        )

        while not es.stop():
            xs = es.ask()
            cs = [resonant16C_cost(x, M0) for x in xs]
            es.tell(xs, cs)
            es.disp()

        if es.best.f < best_cost:
            best_cost = es.best.f
            best_X = es.best.x.copy()

    print("\nBEST Resonant-16C FIT:")
    print(best_X)
    print("cost =", best_cost)

    # Diagnostics at best fit: mirror the cost pipeline exactly
    (
        A_u, B_u,
        A_d, B_d,
        A_nu, B_nu,
        shifts_u, shifts_d, shifts_nu, shifts_e,
        lambda_nu, theta_C, gamma_l
    ) = unpack_params_16C(best_X)

    print("\nUnpacked parameters:")
    print("A_u, B_u    =", A_u, B_u)
    print("A_d, B_d    =", A_d, B_d)
    print("A_nu, B_nu  =", A_nu, B_nu)
    print("shifts_u    =", shifts_u)
    print("shifts_d    =", shifts_d)
    print("shifts_nu   =", shifts_nu)
    print("shifts_e    =", shifts_e)
    print("lambda_nu   =", lambda_nu)
    print("theta_C     =", theta_C)
    print("gamma_l     =", gamma_l)

    # Quark coherence fixed, lepton coherence from fit
    gamma_q = 0.0

    # 9×9 Yukawas (same as in resonant16C_cost)
    Yu9  = build_Y_sector(A_u,  B_u,  EXP_U_BASE,  shifts_u,  gamma_q)
    Yd9  = build_Y_sector(A_d,  B_d,  EXP_D_BASE,  shifts_d,  gamma_q)
    Ye9  = build_Y_sector(A_d,  B_d,  EXP_E_BASE,  shifts_e,  gamma_l)
    Ynu9 = build_Y_sector(A_nu, B_nu, EXP_NU_BASE, shifts_nu, gamma_l)

    # 9→3 Schur
    Yu = schur_9to3(Yu9)
    Yd = schur_9to3(Yd9)
    Ye = schur_9to3(Ye9)

    # Cabibbo twist
    U_C = cabibbo_rotation(theta_C)
    Yd = U_C @ Yd

    # Neutrino projection + seesaw
    P       = build_projection_resonant(lambda_nu)
    Ynu_eff = P @ Ynu9 @ P.conj().T
    Mnu     = triadic_seesaw(M0, Ynu_eff)

    # Observables and chi² (without reg)
    obs, Vckm, U_pmns = compute_observables(Yu, Yd, Ye, Mnu)
    chi2, pulls       = chi2_from_obs(obs)

    print(f"\nFinal evaluation at best fit:")
    print(f"  χ² = {chi2:.3f}")
    print(f"  total cost (χ² + reg) = {best_cost:.3f}\n")

    print("Observables (model vs target, pull in σ):")
    for key in exp_targets.keys():
        print(f"  {key:12s}: model={obs[key]:.6g}, "
              f"target={exp_targets[key]:.6g}, pull={pulls[key]: .3f}")

    print("\n|V_CKM| ≈")
    print(np.abs(Vckm))

    print("\n|U_PMNS| ≈")
    print(np.abs(U_pmns))

    return best_X, best_cost, M0

if __name__ == "__main__":
    optimize_resonant16C(num_restarts=4, seed=9)

"""
=== Resonant-16C Restart 1/4 ===
(10_w,20)-aCMA-ES (mu_w=5.9,w_1=27%) in dimension 17 (seed=1750994040, Mon Dec  8 22:17:22 2025)
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
    1     20 2.787083241766091e+02 1.0e+00 2.96e-01  3e-01  3e-01 0:00.0
    2     40 2.536720599946552e+02 1.2e+00 3.11e-01  3e-01  3e-01 0:00.0
    3     60 1.389187284007203e+02 1.3e+00 3.14e-01  3e-01  4e-01 0:00.0
  100   2000 6.555330932798485e+01 7.8e+00 8.34e-02  2e-02  1e-01 0:01.5
  200   4000 4.745717219840245e+01 3.6e+01 1.12e-01  2e-02  2e-01 0:02.9
  300   6000 3.558346539489390e+01 6.1e+01 4.85e-02  6e-03  8e-02 0:04.3
  400   8000 1.055058614334735e+01 8.3e+01 1.36e-01  1e-02  2e-01 0:05.7
  500  10000 9.173887312189702e+00 1.9e+02 2.95e-02  9e-04  6e-02 0:07.1
  600  12000 9.157897982590979e+00 4.5e+02 1.46e-03  2e-05  3e-03 0:08.5

=== Resonant-16C Restart 2/4 ===
(10_w,20)-aCMA-ES (mu_w=5.9,w_1=27%) in dimension 17 (seed=1096649525, Mon Dec  8 22:17:30 2025)
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
    1     20 5.873021166107926e+02 1.0e+00 2.75e-01  3e-01  3e-01 0:00.0
    2     40 3.468896054595806e+02 1.1e+00 2.67e-01  3e-01  3e-01 0:00.0
    3     60 1.049018217679392e+02 1.2e+00 2.55e-01  2e-01  3e-01 0:00.0
  100   2000 2.495067450433601e+01 9.1e+00 6.53e-02  2e-02  9e-02 0:01.4
  200   4000 8.513983490042873e+00 3.0e+01 4.18e-02  4e-03  5e-02 0:02.8
  300   6000 8.130697097402338e+00 1.0e+02 9.13e-03  4e-04  2e-02 0:04.1
  400   8000 8.098387841597075e+00 1.6e+02 9.71e-03  3e-04  2e-02 0:05.5
  500  10000 8.073391425843983e+00 2.4e+02 1.19e-02  3e-04  2e-02 0:07.1
  600  12000 8.069307318150425e+00 3.6e+02 1.60e-04  2e-06  3e-04 0:08.6

=== Resonant-16C Restart 3/4 ===
(10_w,20)-aCMA-ES (mu_w=5.9,w_1=27%) in dimension 17 (seed=1407148388, Mon Dec  8 22:17:39 2025)
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
    1     20 5.949783507537333e+02 1.0e+00 2.85e-01  3e-01  3e-01 0:00.0
    2     40 9.545029472985387e+01 1.2e+00 2.70e-01  3e-01  3e-01 0:00.0
    3     60 1.189628043814013e+02 1.2e+00 2.65e-01  3e-01  3e-01 0:00.0
  100   2000 2.162972419784850e+01 1.3e+01 1.11e-01  2e-02  2e-01 0:01.9
  200   4000 8.246742473323266e+00 2.4e+01 4.47e-02  4e-03  7e-02 0:03.8
  300   6000 7.290721723576441e+00 6.3e+01 6.53e-03  5e-04  1e-02 0:05.8
  400   8000 7.279366448513139e+00 1.4e+02 9.84e-04  4e-05  2e-03 0:07.8
  500  10000 7.279352948323840e+00 2.5e+02 7.15e-06  1e-07  1e-05 0:09.7
  600  12000 7.279352948161210e+00 4.0e+02 9.13e-07  1e-08  2e-06 0:11.8

=== Resonant-16C Restart 4/4 ===
(10_w,20)-aCMA-ES (mu_w=5.9,w_1=27%) in dimension 17 (seed=1477784359, Mon Dec  8 22:17:51 2025)
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
    1     20 2.835572268924469e+02 1.0e+00 2.76e-01  3e-01  3e-01 0:00.0
    2     40 4.539730024088169e+02 1.1e+00 2.66e-01  3e-01  3e-01 0:00.0
    3     60 8.793886231628213e+01 1.2e+00 2.59e-01  2e-01  3e-01 0:00.0
  100   2000 4.287940411408322e+01 1.0e+01 1.42e-01  4e-02  2e-01 0:01.5
  200   4000 9.817168040498007e+00 1.5e+01 4.56e-02  8e-03  6e-02 0:03.2
  300   6000 7.218234151876592e+00 5.2e+01 3.39e-03  2e-04  5e-03 0:04.6
  400   8000 7.216626050959563e+00 2.7e+02 1.38e-03  5e-05  3e-03 0:06.2
  500  10000 7.216556789378717e+00 4.2e+02 1.58e-05  3e-07  3e-05 0:08.2
  581  11620 7.216556787881613e+00 4.5e+02 3.90e-07  5e-09  5e-07 0:09.8

BEST Resonant-16C FIT:
[ 0.8547712  -1.21900827  0.42595456  0.0291953   0.20844277 -0.64194008
  0.89470344  1.4603491   0.83393488  0.84890353 -0.11832564  0.07085401
 -0.91775019  0.92356035  0.50104073 -0.12438562 -0.56027617]
cost = 7.21655678787628

Unpacked parameters:
A_u, B_u    = 0.8547711968519791 -1.2190082701856326
A_d, B_d    = 0.42595455736950455 0.02919529621174785
A_nu, B_nu  = 0.20844277469945544 -0.6419400847965141
shifts_u    = [0.89470344 1.4603491 ]
shifts_d    = [0.83393488 0.84890353]
shifts_nu   = [-0.11832564  0.07085401]
shifts_e    = [-0.91775019  0.92356035]
lambda_nu   = 0.5010407279220283
theta_C     = -0.12438562067937092
gamma_l     = -0.5602761722687236

Final evaluation at best fit:
  χ² = 5.944
  total cost (χ² + reg) = 7.217

Observables (model vs target, pull in σ):
  m_c/m_t     : model=0.0069334, target=0.007, pull=-0.032
  m_u/m_t     : model=9.9841e-06, target=1e-05, pull=-0.005
  m_s/m_b     : model=0.0238752, target=0.02, pull= 0.646
  m_d/m_b     : model=0.000385721, target=0.001, pull=-2.048
  m_mu/m_tau 
"""

#!/usr/bin/env python3
# ============================================================================
#  TRIADIC ℤ₂₁₆₀ GEOMETRIC ALIGNMENT
#    • 9 continuous parameters:
#        A_u,B_u, A_d,B_d, A_e,B_e, A_nu,B_nu, kappa
#    • 9-site ring, triadic kernel with forbidden distances {2,4,7}
#    • Schur 9→3 Yukawas
#    • Triadic neutrino projector (0,3,6), (1,4,7), (2,5,8)
#    • Full 1-loop SM RGE (Yu, Yd, Ye, κ)
#    • 14 observables with 30% fractional uncertainties
#
#  This version:
#    - Starts from the best-fit you just found (χ²+reg ≈ 89.73)
#    - Prints a full observable table & pulls at the end
# ============================================================================

import numpy as np
import cma
from scipy.integrate import solve_ivp

# --------------------------- Constants ---------------------------

MU_HIGH = 2.0e14
MU_LOW  = 1.0e2
V_HIGGS = 246.0  # GeV

g1_EW, g2_EW, g3_EW = 0.36, 0.65, 1.17
lam_H = 0.13

targets = {
    "m_c/m_t":0.007, "m_u/m_t":1e-5,
    "m_s/m_b":0.02,  "m_d/m_b":0.001,
    "m_mu/m_tau":0.06, "m_e/m_tau":3e-4,
    "theta12_q":0.226, "theta23_q":0.041, "theta13_q":0.0035,
    "theta12_l":0.59,  "theta23_l":0.84,  "theta13_l":0.15,
    "Delta_m2_21":7.4e-5, "Delta_m2_31":2.5e-3,
}
sigmas = {k: 0.3*abs(v) for k,v in targets.items()}

# --------------------------- Triadic kernel ℤ₂₁₆₀ with {2,4,7} forbidden ---------------------------

def triadic_kernel_2160(kappa: float) -> np.ndarray:
    """
    9×9 kernel on a ring with triadic flavor structure.
    Strict forbidden distances: {2,4,7}.
    On a 9-site ring, the non-trivial distances are d=1,2,3,4.
    We still respect the {2,4,7} rule to keep the ℤ₂₁₆₀ triadic logic explicit.
    """
    K = np.zeros((9,9), dtype=float)
    forbidden = {2, 4, 7}

    for i in range(9):
        for j in range(9):
            d = abs(i - j)
            d = min(d, 9 - d)  # ring distance

            if d == 0:
                K[i, j] = 1.0
            elif d in forbidden:
                K[i, j] = 0.0
            else:
                K[i, j] = kappa**d

    return K

# --------------------------- Phase wheels ---------------------------

def phase_matrix(A: float, B: float) -> np.ndarray:
    """
    Phase profile: φ_i = A + B*(i mod 3).
    Returns 9×9 matrix exp[i(φ_i - φ_j)].
    """
    phi = np.array([A + B*(i % 3) for i in range(9)])
    return np.exp(1j * (phi[:, None] - phi[None, :]))

# --------------------------- Yukawa builder ---------------------------

def build_Yukawa(A: float, B: float, kappa: float, alpha: float) -> np.ndarray:
    """
    Build 9×9 Yukawa:
       Y_ij ~ exp[i(φ_i - φ_j)] * K_ij(kappa)
    Then normalize so largest singular value is 1, multiply by alpha.
    """
    Y = phase_matrix(A, B) * triadic_kernel_2160(kappa)
    sv = np.linalg.svd(Y, compute_uv=False)
    if sv[0] > 0:
        Y /= sv[0]
    return alpha * Y

# --------------------------- Schur 9→3 ---------------------------

def schur_9to3(Y9: np.ndarray) -> np.ndarray:
    """
    Take the 3×3 light block via Schur complement over last 6 heavy sites.
    """
    A = Y9[:3, :3]
    B = Y9[:3, 3:]
    D = Y9[3:, 3:]
    Dinv = np.linalg.pinv(D + 1e-10 * np.eye(6))
    return A - B @ Dinv @ B.conj().T

# --------------------------- Proto-Majorana ---------------------------

def proto_majorana(rng: np.random.Generator, scale: float = 7e13) -> np.ndarray:
    """
    Random 9×9 complex symmetric Majorana matrix with spectral norm ≈ scale.
    """
    M = rng.normal(size=(9, 9)) + 1j * rng.normal(size=(9, 9))
    M = 0.5 * (M + M.T.conj())
    sv = np.linalg.svd(M, compute_uv=False)
    if sv[0] > 0:
        M *= scale / sv[0]
    return M

# --------------------------- RGE pack/unpack ---------------------------

def pack(Yu, Yd, Ye, kappa):
    def f(M): return np.concatenate([M.real.ravel(), M.imag.ravel()])
    return np.concatenate([f(Yu), f(Yd), f(Ye), f(kappa)])

def unpack(v):
    n = 3
    N = n * n
    def blk(i):
        re = v[i          : i+N    ].reshape((3,3))
        im = v[i+N        : i+2*N  ].reshape((3,3))
        return re + 1j*im
    return blk(0), blk(2*N), blk(4*N), blk(6*N)

# --------------------------- 1-loop SM RGE ---------------------------

def beta(t, v, g1, g2, g3, lam):
    Yu, Yd, Ye, kappa = unpack(v)

    # Hard clip to keep numerics under control
    for M in (Yu, Yd, Ye):
        np.clip(M, -20, 20, out=M)

    T = np.trace(3 * Yu @ Yu.conj().T +
                 3 * Yd @ Yd.conj().T +
                 Ye @ Ye.conj().T).real

    pref = 1.0 / (16.0 * np.pi**2)

    dYu = pref * (
        Yu * (T - (17.0/20.0)*g1**2 - (9.0/4.0)*g2**2 - 8.0*g3**2)
        + 1.5 * (Yu @ Yu.conj().T @ Yu - Yd @ Yd.conj().T @ Yu)
    )

    dYd = pref * (
        Yd * (T - (1.0/4.0)*g1**2 - (9.0/4.0)*g2**2 - 8.0*g3**2)
        + 1.5 * (Yd @ Yd.conj().T @ Yd - Yu @ Yu.conj().T @ Yd)
    )

    dYe = pref * (
        Ye * (T - (9.0/4.0)*g1**2 - (9.0/4.0)*g2**2)
        + 1.5 * (Ye @ Ye.conj().T @ Ye)
    )

    YeT = Ye @ Ye.conj().T
    dkappa = pref * (
        (-3.0 * g2**2 + lam) * kappa + (YeT @ kappa + kappa @ YeT.T)
    )

    return pack(dYu, dYd, dYe, dkappa)

def run_rge(Yu, Yd, Ye, kappa_high):
    sol = solve_ivp(
        beta,
        [np.log(MU_HIGH), np.log(MU_LOW)],
        pack(Yu, Yd, Ye, kappa_high),
        args=(g1_EW, g2_EW, g3_EW, lam_H),
        rtol=1e-5, atol=1e-8,
        method='RK45', max_step=0.4
    )
    return unpack(sol.y[:, -1])

# --------------------------- Observables ---------------------------

def extract_angles(U: np.ndarray):
    a = np.abs(U)
    s13 = a[0, 2]
    c13 = np.sqrt(max(0.0, 1.0 - s13**2))
    s12 = a[0, 1] / c13 if c13 > 1e-10 else 0.0
    s23 = a[1, 2] / c13 if c13 > 1e-10 else 0.0
    return (
        np.arcsin(np.clip(s12, 0.0, 1.0)),
        np.arcsin(np.clip(s23, 0.0, 1.0)),
        np.arcsin(s13)
    )

def get_obs(Yu, Yd, Ye, Mnu):
    su = np.sort(np.linalg.svd(Yu, compute_uv=False))
    sd = np.sort(np.linalg.svd(Yd, compute_uv=False))
    se = np.sort(np.linalg.svd(Ye, compute_uv=False))

    obs = {
        "m_c/m_t": su[1] / su[2],
        "m_u/m_t": su[0] / su[2],
        "m_s/m_b": sd[1] / sd[2],
        "m_d/m_b": sd[0] / sd[2],
        "m_mu/m_tau": se[1] / se[2],
        "m_e/m_tau": se[0] / se[2],
    }

    # CKM
    Uu = np.linalg.svd(Yu)[0]
    Ud = np.linalg.svd(Yd)[0]
    Vckm = Uu.conj().T @ Ud
    th12q, th23q, th13q = extract_angles(Vckm)
    obs["theta12_q"] = th12q
    obs["theta23_q"] = th23q
    obs["theta13_q"] = th13q

    # Neutrinos
    evals, U_nu = np.linalg.eigh(0.5 * (Mnu + Mnu.T))
    mnu = np.sort(np.abs(evals))
    Ue = np.linalg.svd(Ye)[0]
    Upmns = Ue.conj().T @ U_nu
    th12l, th23l, th13l = extract_angles(Upmns)

    obs["theta12_l"] = th12l
    obs["theta23_l"] = th23l
    obs["theta13_l"] = th13l

    obs["Delta_m2_21"] = mnu[1]**2 - mnu[0]**2
    obs["Delta_m2_31"] = mnu[2]**2 - mnu[0]**2

    return obs, Vckm, Upmns

# --------------------------- Cost function (9 params) ---------------------------

def cost(X, M0):
    A_u, B_u, A_d, B_d, A_e, B_e, A_nu, B_nu, kappa = X

    # Fixed high-scale normalizations (empirical, roughly SM-like)
    alpha_u  = 0.71
    alpha_d  = 0.095
    alpha_e  = 0.082
    alpha_nu = 0.13

    # Build 9×9 Yukawas
    Yu9  = build_Yukawa(A_u,  B_u,  kappa, alpha_u)
    Yd9  = build_Yukawa(A_d,  B_d,  kappa, alpha_d)
    Ye9  = build_Yukawa(A_e,  B_e,  kappa, alpha_e)
    Ynu9 = build_Yukawa(A_nu, B_nu, kappa, alpha_nu)

    # Schur 9→3
    Yu_h = schur_9to3(Yu9)
    Yd_h = schur_9to3(Yd9)
    Ye_h = schur_9to3(Ye9)

    # Triadic neutrino projector
    P = np.zeros((3, 9), dtype=complex)
    for c, sites in enumerate([(0, 3, 6), (1, 4, 7), (2, 5, 8)]):
        P[c, sites] = 1.0 / np.sqrt(3.0)

    Ynu_eff = P @ Ynu9 @ P.conj().T
    MR = P @ M0 @ P.conj().T

    Mnu_h = -0.5 * V_HIGGS**2 * (
        Ynu_eff @ np.linalg.pinv(MR + 1e-8 * np.eye(3)) @ Ynu_eff.T
    )
    kappa_h = Mnu_h / V_HIGGS**2

    # Run RGEs down to MU_LOW
    Yu_l, Yd_l, Ye_l, kappa_l = run_rge(Yu_h, Yd_h, Ye_h, kappa_h)
    Mnu_l = kappa_l * V_HIGGS**2

    obs, _, _ = get_obs(Yu_l, Yd_l, Ye_l, Mnu_l)

    chi2 = 0.0
    for k in targets:
        chi2 += ((obs[k] - targets[k]) / sigmas[k])**2

    # Mild quadratic regularization on parameters
    reg = 0.05 * np.sum(X**2)
    return chi2 + reg

# --------------------------- Diagnostics at a point ---------------------------

def evaluate_point(X, M0):
    """
    Reconstruct low-scale matrices and print a neat summary of observables & pulls.
    """
    A_u, B_u, A_d, B_d, A_e, B_e, A_nu, B_nu, kappa = X

    alpha_u  = 0.71
    alpha_d  = 0.095
    alpha_e  = 0.082
    alpha_nu = 0.13

    Yu9  = build_Yukawa(A_u,  B_u,  kappa, alpha_u)
    Yd9  = build_Yukawa(A_d,  B_d,  kappa, alpha_d)
    Ye9  = build_Yukawa(A_e,  B_e,  kappa, alpha_e)
    Ynu9 = build_Yukawa(A_nu, B_nu, kappa, alpha_nu)

    Yu_h = schur_9to3(Yu9)
    Yd_h = schur_9to3(Yd9)
    Ye_h = schur_9to3(Ye9)

    P = np.zeros((3, 9), dtype=complex)
    for c, sites in enumerate([(0, 3, 6), (1, 4, 7), (2, 5, 8)]):
        P[c, sites] = 1.0 / np.sqrt(3.0)

    Ynu_eff = P @ Ynu9 @ P.conj().T
    MR = P @ M0 @ P.conj().T
    Mnu_h = -0.5 * V_HIGGS**2 * (
        Ynu_eff @ np.linalg.pinv(MR + 1e-8*np.eye(3)) @ Ynu_eff.T
    )
    kappa_h = Mnu_h / V_HIGGS**2

    Yu_l, Yd_l, Ye_l, kappa_l = run_rge(Yu_h, Yd_h, Ye_h, kappa_h)
    Mnu_l = kappa_l * V_HIGGS**2

    obs, Vckm, Upmns = get_obs(Yu_l, Yd_l, Ye_l, Mnu_l)

    chi2 = 0.0
    print("\n=== OBSERVABLES AT THIS POINT ===")
    for k in targets:
        pull = (obs[k] - targets[k]) / sigmas[k]
        chi2 += pull**2
        print(f"{k:12s}: model={obs[k]: .6e}, target={targets[k]: .6e}, pull={pull: .3f}")
    print(f"\nχ² (obs only) = {chi2:.3f}")
    print(f"χ² + reg      = {cost(X, M0):.3f}")

    print("\n|V_CKM| ≈")
    print(np.abs(Vckm))

    print("\n|U_PMNS| ≈")
    print(np.abs(Upmns))

# --------------------------- MAIN OPTIMIZATION ---------------------------

if __name__ == "__main__":
    rng = np.random.default_rng(777)
    M0 = proto_majorana(rng)

    # Best point from your last run (χ²+reg ≈ 89.73)
    BEST_PREVIOUS = np.array([
        -8.43580970e-02,
         2.67331725e-04,
        -4.21127378e-02,
        -3.89722997e-02,
        -1.67578764e-01,
         1.02680387e+00,
         4.08866876e-02,
         1.81084179e-01,
         1.11508982e+00,
    ])

    # Start near this known good point, with a modest step-size
    x0 = BEST_PREVIOUS.copy()
    sigma0 = 0.2

    es = cma.CMAEvolutionStrategy(
        x0,
        sigma0,
        {
            'popsize': 80,
            'maxiter': 3000,
            'seed': 42,
            'verb_disp': 1,
        }
    )

    print("Starting refined ℤ₂₁₆₀ triadic optimization from previous best...")

    while not es.stop():
        xs = es.ask()
        cs = [cost(x, M0) for x in xs]
        es.tell(xs, cs)
        es.disp()

    print("\n=== FINAL ℤ₂₁₆₀ TRIADIC RESULT ===")
    print("Best χ²+reg =", es.best.f)
    print("Best parameters:", es.best.x)

    # Final detailed evaluation at best fit
    evaluate_point(es.best.x, M0)

"""
RESULTS:
/Users/chazzromeo/rAI/Flavor/.venv/bin/python /Users/chazzromeo/rAI/Flavor/flavor/align-2160-2.py 
/Users/chazzromeo/rAI/Flavor/.venv/lib/python3.11/site-packages/cma/s.py:15: UserWarning: Could not import matplotlib.pyplot, therefore ``cma.plot()`` etc. is not available
  _warnings.warn('Could not import matplotlib.pyplot, therefore'
(40_w,80)-aCMA-ES (mu_w=21.8,w_1=9%) in dimension 9 (seed=42, Sun Dec  7 20:53:32 2025)
Starting refined ℤ₂₁₆₀ triadic optimization from previous best...
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
    1     80 8.181553878132983e+04 1.0e+00 1.97e-01  2e-01  2e-01 0:02.3
    2    160 6.654389778658525e+02 1.5e+00 1.94e-01  1e-01  2e-01 0:04.2
    3    240 7.021521679709258e+04 2.0e+00 1.94e-01  1e-01  2e-01 0:06.3
    4    320 2.253130643476421e+04 2.6e+00 1.99e-01  9e-02  2e-01 0:08.2
    5    400 7.852909516604639e+04 3.4e+00 1.89e-01  7e-02  3e-01 0:10.2
    6    480 4.715353176625250e+02 4.3e+00 1.82e-01  5e-02  3e-01 0:12.3
    7    560 5.188252992209246e+02 5.9e+00 1.90e-01  4e-02  3e-01 0:15.1
    8    640 1.144131708977017e+04 7.8e+00 1.75e-01  3e-02  2e-01 0:17.8
    9    720 5.832204962745229e+02 9.4e+00 1.63e-01  2e-02  2e-01 0:20.1
   10    800 1.002969502382369e+03 1.1e+01 1.99e-01  2e-02  3e-01 0:22.2
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
   11    880 1.192752808633739e+02 1.3e+01 1.89e-01  2e-02  2e-01 0:24.4
   12    960 7.650854987900727e+03 1.6e+01 1.84e-01  1e-02  2e-01 0:26.7
   13   1040 5.232487540486071e+02 2.0e+01 1.76e-01  1e-02  2e-01 0:29.5
   14   1120 4.399642543854195e+02 2.4e+01 1.59e-01  8e-03  2e-01 0:32.7
   15   1200 2.628986066763080e+02 2.9e+01 1.50e-01  7e-03  2e-01 0:35.3
   16   1280 1.680493670250299e+02 3.5e+01 1.61e-01  6e-03  2e-01 0:38.2
   17   1360 4.867670858160653e+02 4.3e+01 1.67e-01  5e-03  2e-01 0:42.2
   18   1440 1.262678928162259e+02 5.2e+01 1.68e-01  5e-03  2e-01 0:45.8
   19   1520 1.469946674238903e+02 6.0e+01 1.59e-01  4e-03  2e-01 0:48.8
   20   1600 2.626231525486893e+02 6.6e+01 1.44e-01  3e-03  2e-01 0:51.8
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
   21   1680 1.032766335084268e+02 7.6e+01 1.32e-01  2e-03  2e-01 0:54.9
   22   1760 2.211319532418406e+02 1.0e+02 1.24e-01  2e-03  2e-01 0:57.8
   23   1840 1.326621857276322e+02 1.2e+02 1.15e-01  2e-03  2e-01 1:00.7
   24   1920 1.043352270761221e+02 1.3e+02 1.21e-01  1e-03  2e-01 1:03.8
   25   2000 1.147993892792242e+02 1.5e+02 1.19e-01  1e-03  1e-01 1:06.9
   26   2080 1.310169945796936e+02 1.7e+02 1.05e-01  8e-04  1e-01 1:09.6
   27   2160 9.960713191048232e+01 2.1e+02 1.25e-01  8e-04  2e-01 1:13.0
   28   2240 1.023284033910908e+02 2.5e+02 1.15e-01  6e-04  1e-01 1:16.3
   29   2320 1.073983110211184e+02 3.1e+02 1.09e-01  5e-04  1e-01 1:19.2
   30   2400 1.041479158518228e+02 3.8e+02 1.04e-01  4e-04  1e-01 1:21.9
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
   31   2480 9.935573225860951e+01 4.8e+02 9.69e-02  3e-04  1e-01 1:24.7
   32   2560 1.051808928544735e+02 5.3e+02 8.93e-02  2e-04  1e-01 1:27.6
   33   2640 9.824338788040815e+01 6.3e+02 8.62e-02  2e-04  1e-01 1:30.1
   34   2720 9.978812991720027e+01 7.0e+02 8.04e-02  1e-04  9e-02 1:32.5
   35   2800 9.772831651886503e+01 8.3e+02 7.24e-02  1e-04  8e-02 1:34.9
   36   2880 9.786995421706447e+01 9.5e+02 6.45e-02  8e-05  7e-02 1:37.3
   37   2960 9.724377404835977e+01 1.1e+03 6.69e-02  7e-05  7e-02 1:39.6
   38   3040 9.715846922017673e+01 1.3e+03 6.08e-02  5e-05  6e-02 1:42.0
   39   3120 9.713649483862302e+01 1.5e+03 5.86e-02  4e-05  6e-02 1:44.4
   40   3200 9.611374784010218e+01 1.7e+03 5.82e-02  4e-05  6e-02 1:46.9
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
   41   3280 9.700517235264195e+01 2.0e+03 5.72e-02  3e-05  6e-02 1:49.3
   42   3360 9.654700037553283e+01 2.3e+03 6.31e-02  3e-05  6e-02 1:51.6
   43   3440 9.699354189640269e+01 2.7e+03 5.62e-02  2e-05  5e-02 1:54.0
   44   3520 9.683340377405264e+01 3.1e+03 4.87e-02  2e-05  5e-02 1:56.5
   45   3600 9.685702957451208e+01 3.7e+03 4.88e-02  2e-05  5e-02 1:59.0
   46   3680 9.458282485475304e+01 4.3e+03 4.76e-02  1e-05  5e-02 2:01.5
   47   3760 9.579403307240166e+01 5.0e+03 4.72e-02  1e-05  5e-02 2:04.0
   48   3840 9.562252121305686e+01 5.5e+03 4.71e-02  1e-05  4e-02 2:06.4
   49   3920 9.572488147731788e+01 6.7e+03 4.62e-02  1e-05  5e-02 2:09.6
   50   4000 9.583744909249455e+01 6.9e+03 5.22e-02  1e-05  6e-02 2:12.0
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
   51   4080 9.581468651767734e+01 7.8e+03 5.77e-02  1e-05  6e-02 2:14.4
   52   4160 9.566538826243631e+01 8.8e+03 5.95e-02  9e-06  7e-02 2:17.3
   53   4240 9.571639263255446e+01 1.0e+04 5.62e-02  7e-06  6e-02 2:20.3
   54   4320 9.122566150668318e+01 1.2e+04 5.80e-02  6e-06  6e-02 2:22.9
NOTE (module=cma, iteration=54):  
condition in coordinate system exceeded 1.1e+08, rescaled to 1.0e+00, 
condition changed from 1.5e+08 to 6.3e+02
   55   4400 9.557466704872216e+01 2.5e+01 5.23e-02  5e-06  6e-02 2:26.0
   56   4480 9.561309211889673e+01 2.8e+01 4.91e-02  4e-06  6e-02 2:28.5
   57   4560 9.550077911513978e+01 3.1e+01 4.94e-02  4e-06  6e-02 2:30.8
   58   4640 9.548257300076176e+01 3.4e+01 4.96e-02  3e-06  6e-02 2:33.5
   59   4720 9.556325797917879e+01 3.6e+01 4.98e-02  3e-06  6e-02 2:36.1
   60   4800 9.159483124575983e+01 3.6e+01 5.01e-02  3e-06  6e-02 2:38.7
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
   61   4880 9.548239465970491e+01 3.8e+01 4.82e-02  3e-06  5e-02 2:41.2
   62   4960 9.550785435456666e+01 4.0e+01 4.40e-02  2e-06  5e-02 2:43.5
   63   5040 9.547973399434557e+01 4.0e+01 4.20e-02  2e-06  4e-02 2:46.2
   64   5120 9.548785770957781e+01 4.0e+01 4.49e-02  2e-06  4e-02 2:49.0
   65   5200 9.550775468594044e+01 4.4e+01 5.00e-02  2e-06  5e-02 2:51.7
   66   5280 9.548378995782541e+01 4.5e+01 4.21e-02  1e-06  4e-02 2:54.0
   67   5360 9.548736956205883e+01 4.7e+01 4.35e-02  1e-06  4e-02 2:56.2
   68   5440 9.551287951407329e+01 4.8e+01 3.97e-02  8e-07  4e-02 2:58.5
   69   5520 9.547285402049340e+01 5.2e+01 3.72e-02  6e-07  4e-02 3:00.7
   70   5600 9.547192935592705e+01 5.7e+01 3.36e-02  4e-07  4e-02 3:03.1
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
   71   5680 9.547440758427402e+01 6.2e+01 3.21e-02  3e-07  4e-02 3:05.7
   72   5760 9.546967688329697e+01 6.2e+01 3.20e-02  3e-07  3e-02 3:08.4
   73   5840 9.547040806675625e+01 6.5e+01 2.93e-02  2e-07  3e-02 3:11.1
   74   5920 9.547025743053239e+01 6.5e+01 2.98e-02  2e-07  3e-02 3:13.4
   75   6000 9.546731373158028e+01 7.4e+01 2.81e-02  1e-07  3e-02 3:15.6
   76   6080 9.546766412991131e+01 8.0e+01 2.94e-02  1e-07  3e-02 3:17.9
   77   6160 9.546670187573123e+01 9.0e+01 3.22e-02  8e-08  4e-02 3:20.1
   78   6240 9.546560032026764e+01 1.0e+02 4.02e-02  9e-08  5e-02 3:22.4
   79   6320 9.546302580766660e+01 1.0e+02 4.22e-02  7e-08  5e-02 3:24.6
   80   6400 9.546274560504636e+01 1.1e+02 5.05e-02  7e-08  6e-02 3:27.2
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
   81   6480 9.546134709181142e+01 1.4e+02 5.89e-02  7e-08  7e-02 3:29.9
   82   6560 9.545508490776697e+01 1.7e+02 6.48e-02  6e-08  9e-02 3:32.3
   83   6640 9.545809664500406e+01 2.0e+02 7.97e-02  7e-08  1e-01 3:34.5
   84   6720 9.545656273675566e+01 2.2e+02 7.80e-02  6e-08  1e-01 3:36.8
   85   6800 9.545813756279105e+01 2.5e+02 7.98e-02  5e-08  1e-01 3:39.1
   86   6880 9.545604810590842e+01 2.8e+02 8.36e-02  5e-08  1e-01 3:41.4
   87   6960 9.545354224430992e+01 3.1e+02 9.66e-02  5e-08  1e-01 3:44.0
   88   7040 9.357202063525794e+01 3.7e+02 9.59e-02  4e-08  1e-01 3:46.6
   89   7120 9.101301131768106e+01 4.4e+02 9.70e-02  4e-08  1e-01 3:49.0
   90   7200 9.194154450326087e+01 5.2e+02 9.86e-02  4e-08  2e-01 3:51.3
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
   91   7280 9.545203204112578e+01 6.1e+02 9.52e-02  3e-08  1e-01 3:53.5
   92   7360 9.545159109399262e+01 6.7e+02 9.58e-02  3e-08  1e-01 3:55.8
   93   7440 9.545081413916297e+01 6.8e+02 9.09e-02  2e-08  1e-01 3:58.1
   94   7520 9.517062190552646e+01 7.0e+02 9.19e-02  2e-08  1e-01 4:00.4
   95   7600 9.544986769932977e+01 7.4e+02 9.73e-02  2e-08  1e-01 4:02.8
   96   7680 9.231313778844851e+01 7.6e+02 1.02e-01  2e-08  1e-01 4:05.6
   97   7760 9.544983444685633e+01 8.0e+02 1.13e-01  2e-08  1e-01 4:08.0
   98   7840 9.544926304253347e+01 8.2e+02 1.00e-01  2e-08  1e-01 4:10.4
   99   7920 9.544875369231633e+01 8.1e+02 9.64e-02  1e-08  9e-02 4:12.8
  100   8000 9.544891901647961e+01 8.5e+02 8.97e-02  1e-08  7e-02 4:15.4
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  101   8080 9.219385699319280e+01 8.8e+02 9.20e-02  1e-08  7e-02 4:18.1
  102   8160 9.544862863424196e+01 8.7e+02 9.58e-02  1e-08  6e-02 4:21.2
  103   8240 9.544842593411380e+01 9.1e+02 1.07e-01  1e-08  7e-02 4:23.9
  104   8320 9.544867806379018e+01 9.7e+02 1.06e-01  1e-08  6e-02 4:26.4
  105   8400 9.544881220754868e+01 9.8e+02 9.61e-02  8e-09  5e-02 4:29.4
  106   8480 9.544854320151990e+01 9.5e+02 9.87e-02  8e-09  5e-02 4:32.0
  107   8560 9.544851533071818e+01 9.2e+02 1.09e-01  8e-09  5e-02 4:34.7
  108   8640 9.010044119027711e+01 9.1e+02 1.20e-01  8e-09  5e-02 4:37.5
  109   8720 9.544849095678650e+01 9.2e+02 1.30e-01  9e-09  5e-02 4:39.8
  110   8800 9.544826461727338e+01 9.3e+02 1.30e-01  8e-09  5e-02 4:42.1
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  111   8880 9.544824249066237e+01 9.2e+02 1.44e-01  8e-09  5e-02 4:44.5
  112   8960 9.544816798682638e+01 9.4e+02 1.52e-01  8e-09  4e-02 4:46.8
  113   9040 9.544834194011315e+01 9.1e+02 1.42e-01  7e-09  4e-02 4:49.0
  114   9120 9.544824240974317e+01 9.5e+02 1.41e-01  6e-09  3e-02 4:51.3
  115   9200 9.237635504878476e+01 9.9e+02 1.28e-01  5e-09  3e-02 4:53.5
  116   9280 9.544810121305966e+01 1.0e+03 1.29e-01  5e-09  3e-02 4:55.8
  117   9360 9.403122639402176e+01 1.0e+03 1.25e-01  4e-09  3e-02 4:58.0
  118   9440 9.544804474159196e+01 9.8e+02 1.34e-01  5e-09  3e-02 5:00.3
  119   9520 9.003536796431972e+01 8.8e+02 1.30e-01  4e-09  2e-02 5:02.6
  120   9600 9.544802330648854e+01 9.4e+02 1.20e-01  4e-09  2e-02 5:04.9
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  121   9680 9.435888237925204e+01 9.0e+02 1.13e-01  3e-09  2e-02 5:07.1
  122   9760 9.544801234368619e+01 9.2e+02 1.09e-01  3e-09  2e-02 5:09.4
  123   9840 9.544800012975891e+01 9.0e+02 9.72e-02  2e-09  1e-02 5:11.7
  124   9920 9.544798224896428e+01 9.2e+02 8.38e-02  2e-09  1e-02 5:13.9
  125  10000 9.544798580223754e+01 9.3e+02 8.10e-02  2e-09  1e-02 5:16.2
  126  10080 9.544797857970362e+01 9.2e+02 7.78e-02  2e-09  9e-03 5:18.4
  127  10160 9.544797893831833e+01 8.8e+02 7.23e-02  1e-09  8e-03 5:20.6
  128  10240 9.544797473974008e+01 8.6e+02 6.76e-02  1e-09  7e-03 5:22.8
  129  10320 9.544797279342426e+01 8.5e+02 6.07e-02  9e-10  5e-03 5:25.0
  130  10400 9.544797240308289e+01 8.6e+02 5.57e-02  7e-10  5e-03 5:27.2
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  131  10480 9.544797120493051e+01 8.9e+02 5.54e-02  7e-10  4e-03 5:29.5
  132  10560 9.544797112062204e+01 8.9e+02 5.85e-02  7e-10  4e-03 5:31.7
  133  10640 9.544797110452768e+01 9.4e+02 5.92e-02  6e-10  4e-03 5:33.9
  134  10720 9.544797078632992e+01 9.9e+02 5.22e-02  5e-10  3e-03 5:36.1
  135  10800 9.544797007583448e+01 9.8e+02 4.89e-02  4e-10  3e-03 5:38.3
  136  10880 9.544796938832330e+01 9.9e+02 4.94e-02  4e-10  3e-03 5:40.5
  137  10960 9.544796965423117e+01 9.8e+02 4.58e-02  3e-10  2e-03 5:42.7
  138  11040 9.544796954059190e+01 9.9e+02 3.98e-02  3e-10  2e-03 5:44.9
  139  11120 9.544796960556259e+01 1.0e+03 3.97e-02  3e-10  1e-03 5:47.1
  140  11200 9.544796958614826e+01 9.9e+02 3.72e-02  2e-10  1e-03 5:49.3
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  141  11280 9.544796937315768e+01 9.7e+02 3.50e-02  2e-10  1e-03 5:51.5
  142  11360 9.544796933688762e+01 9.7e+02 3.34e-02  2e-10  1e-03 5:53.8
  143  11440 9.544796929991844e+01 1.0e+03 3.00e-02  1e-10  8e-04 5:56.0
  144  11520 9.544796926277419e+01 9.8e+02 2.76e-02  1e-10  7e-04 5:58.2
  145  11600 9.544796922598650e+01 1.0e+03 2.64e-02  1e-10  6e-04 6:00.5
  146  11680 9.544796920183698e+01 1.1e+03 2.37e-02  8e-11  5e-04 6:02.8
  147  11760 9.544796918870504e+01 1.1e+03 2.48e-02  7e-11  5e-04 6:05.0
  148  11840 9.544796920486195e+01 1.1e+03 2.31e-02  6e-11  4e-04 6:07.2
  149  11920 9.544796919292644e+01 1.1e+03 2.01e-02  5e-11  4e-04 6:09.5
  150  12000 9.544796919794446e+01 1.1e+03 2.03e-02  5e-11  3e-04 6:11.7
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  151  12080 9.544796918601851e+01 1.1e+03 2.06e-02  5e-11  3e-04 6:14.1
  152  12160 9.544796918604881e+01 1.1e+03 1.91e-02  4e-11  3e-04 6:16.4
  153  12240 9.544796918358237e+01 1.1e+03 1.86e-02  4e-11  2e-04 6:18.6
  154  12320 9.544796918318004e+01 9.6e+02 1.72e-02  3e-11  2e-04 6:20.8
  155  12400 9.544796917685501e+01 1.0e+03 1.44e-02  3e-11  2e-04 6:23.1
  156  12480 9.544796917801806e+01 1.0e+03 1.27e-02  2e-11  1e-04 6:25.4
  157  12560 9.544796917320900e+01 1.0e+03 1.22e-02  2e-11  1e-04 6:27.6
  158  12640 9.544796917430133e+01 1.0e+03 1.08e-02  2e-11  1e-04 6:29.9
  159  12720 9.544796917259275e+01 1.1e+03 1.01e-02  1e-11  1e-04 6:32.4
  160  12800 9.544796917270747e+01 1.1e+03 9.57e-03  1e-11  9e-05 6:34.7
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  161  12880 9.544796917215280e+01 1.2e+03 9.51e-03  1e-11  9e-05 6:37.0
  162  12960 9.544796917166509e+01 1.1e+03 8.45e-03  1e-11  7e-05 6:39.2
  163  13040 9.544796917130945e+01 1.1e+03 8.98e-03  1e-11  7e-05 6:41.6
  164  13120 9.544796917202368e+01 1.0e+03 8.07e-03  1e-11  6e-05 6:43.9
  165  13200 9.544796917087291e+01 9.5e+02 7.94e-03  1e-11  5e-05 6:46.2
  166  13280 9.544796917129558e+01 9.0e+02 7.17e-03  8e-12  5e-05 6:48.4
  167  13360 9.544796917072230e+01 9.6e+02 6.44e-03  7e-12  4e-05 6:50.7
  168  13440 9.544796917060600e+01 1.0e+03 5.89e-03  6e-12  4e-05 6:52.9
  169  13520 9.544796917101370e+01 1.0e+03 5.41e-03  5e-12  3e-05 6:55.2
  170  13600 9.544796917041178e+01 9.6e+02 5.24e-03  5e-12  3e-05 6:57.5
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  171  13680 9.544796917027452e+01 9.9e+02 5.60e-03  5e-12  3e-05 6:59.7
  172  13760 9.544796917014423e+01 1.0e+03 5.41e-03  5e-12  3e-05 7:02.0
  173  13840 9.544796917001602e+01 1.1e+03 5.55e-03  5e-12  3e-05 7:04.2
  174  13920 9.544796917060998e+01 1.0e+03 5.85e-03  5e-12  3e-05 7:06.5
  175  14000 9.544796917045035e+01 9.8e+02 5.47e-03  4e-12  2e-05 7:08.7
  176  14080 9.544796917061377e+01 9.6e+02 5.66e-03  4e-12  3e-05 7:11.0
  177  14160 9.544796917032627e+01 9.9e+02 5.90e-03  4e-12  3e-05 7:13.2
  178  14240 9.544796917036378e+01 1.1e+03 5.51e-03  4e-12  3e-05 7:15.5
  179  14320 9.544796917002486e+01 1.1e+03 5.30e-03  3e-12  2e-05 7:17.7
  180  14400 9.544796917015424e+01 1.2e+03 5.22e-03  3e-12  3e-05 7:20.2
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  181  14480 9.544796917030047e+01 1.3e+03 5.33e-03  3e-12  3e-05 7:22.5
  182  14560 9.544796917004237e+01 1.3e+03 4.99e-03  3e-12  2e-05 7:24.7
  183  14640 9.544796916987565e+01 1.3e+03 5.42e-03  3e-12  3e-05 7:27.0
  184  14720 9.544796917022283e+01 1.3e+03 5.07e-03  3e-12  2e-05 7:29.2
  185  14800 9.544796917019346e+01 1.3e+03 4.94e-03  3e-12  2e-05 7:31.5
  186  14880 9.544796917033155e+01 1.3e+03 5.07e-03  3e-12  2e-05 7:33.8
  187  14960 9.544796917041639e+01 1.4e+03 4.81e-03  2e-12  2e-05 7:36.0
  188  15040 9.544796917006292e+01 1.4e+03 4.94e-03  2e-12  2e-05 7:38.3
  189  15120 9.544796917002928e+01 1.4e+03 4.93e-03  2e-12  2e-05 7:40.5
  190  15200 9.544796917020560e+01 1.4e+03 4.45e-03  2e-12  2e-05 7:42.8
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  191  15280 9.544796917018790e+01 1.4e+03 3.97e-03  2e-12  2e-05 7:45.0
  192  15360 9.544796917045896e+01 1.3e+03 3.88e-03  2e-12  2e-05 7:47.3
  193  15440 9.544796916995969e+01 1.4e+03 4.29e-03  2e-12  2e-05 7:49.5
  194  15520 9.544796917018265e+01 1.4e+03 3.94e-03  2e-12  2e-05 7:51.8
  195  15600 9.544796917032056e+01 1.4e+03 3.46e-03  2e-12  1e-05 7:54.0
  196  15680 9.544796916984747e+01 1.4e+03 3.75e-03  2e-12  1e-05 7:56.6
  197  15760 9.544796917017965e+01 1.6e+03 3.75e-03  2e-12  1e-05 7:58.9
  198  15840 9.544796917028896e+01 1.6e+03 4.14e-03  2e-12  2e-05 8:01.2
  199  15920 9.544796917021516e+01 1.6e+03 4.37e-03  2e-12  2e-05 8:03.5
  200  16000 9.544796917039096e+01 1.7e+03 5.32e-03  3e-12  2e-05 8:05.8
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  201  16080 9.544796917024492e+01 1.7e+03 5.04e-03  3e-12  2e-05 8:08.0
  202  16160 9.544796917032707e+01 1.8e+03 4.87e-03  2e-12  2e-05 8:10.2
  203  16240 9.544796917010125e+01 1.8e+03 4.74e-03  2e-12  2e-05 8:12.5
  204  16320 9.544796917047361e+01 1.7e+03 4.48e-03  2e-12  1e-05 8:14.7
  205  16400 9.544796917009904e+01 1.8e+03 4.20e-03  2e-12  1e-05 8:16.9
  206  16480 9.544796917037569e+01 1.8e+03 4.14e-03  2e-12  1e-05 8:19.1
  207  16560 9.544796917038457e+01 1.8e+03 4.26e-03  2e-12  1e-05 8:21.3
  208  16640 9.544796916967672e+01 1.8e+03 4.04e-03  2e-12  1e-05 8:23.5
  209  16720 9.544796917049611e+01 1.8e+03 3.97e-03  2e-12  1e-05 8:25.7
  210  16800 9.544796916997539e+01 1.6e+03 3.46e-03  1e-12  9e-06 8:28.0
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  211  16880 9.544796917019532e+01 1.6e+03 3.72e-03  1e-12  1e-05 8:30.2
  212  16960 9.544796917005976e+01 1.8e+03 3.90e-03  1e-12  1e-05 8:32.4
  213  17040 9.544796917002792e+01 1.9e+03 3.81e-03  1e-12  1e-05 8:34.6
  214  17120 9.544796917014277e+01 2.0e+03 4.08e-03  1e-12  1e-05 8:36.9
  215  17200 9.544796917004354e+01 2.2e+03 4.37e-03  2e-12  1e-05 8:39.1
  216  17280 9.544796917017068e+01 2.3e+03 4.26e-03  1e-12  1e-05 8:41.4
  217  17360 9.544796917041408e+01 2.4e+03 4.06e-03  1e-12  1e-05 8:43.6
  218  17440 9.544796917016508e+01 2.3e+03 3.95e-03  1e-12  1e-05 8:45.8
  219  17520 9.544796917035623e+01 2.5e+03 4.13e-03  1e-12  1e-05 8:48.1
  220  17600 9.544796917035302e+01 2.5e+03 4.56e-03  2e-12  1e-05 8:50.3
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  221  17680 9.544796917031729e+01 2.3e+03 4.76e-03  2e-12  1e-05 8:52.6
  222  17760 9.544796917000399e+01 2.2e+03 4.24e-03  2e-12  1e-05 8:54.9
  223  17840 9.544796917006698e+01 2.3e+03 4.50e-03  2e-12  1e-05 8:57.1
  224  17920 9.544796917010014e+01 2.5e+03 4.33e-03  1e-12  1e-05 8:59.4
  225  18000 9.544796917015289e+01 2.2e+03 4.84e-03  2e-12  1e-05 9:01.6
  226  18080 9.544796917019112e+01 2.1e+03 5.12e-03  2e-12  2e-05 9:04.0
  227  18160 9.544796917035048e+01 2.1e+03 5.35e-03  2e-12  2e-05 9:06.3
  228  18240 9.544796917025641e+01 2.3e+03 5.51e-03  2e-12  2e-05 9:08.6
  229  18320 9.544796917028263e+01 2.2e+03 5.90e-03  2e-12  2e-05 9:10.8
  230  18400 9.544796916998334e+01 2.0e+03 6.01e-03  2e-12  2e-05 9:13.1
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  231  18480 9.544796916995961e+01 1.9e+03 6.31e-03  2e-12  2e-05 9:15.5
  232  18560 9.544796916977215e+01 2.0e+03 6.31e-03  2e-12  2e-05 9:17.8
  233  18640 9.544796917000889e+01 1.8e+03 6.53e-03  2e-12  2e-05 9:20.1
  234  18720 9.544796917002978e+01 1.8e+03 6.30e-03  2e-12  2e-05 9:22.4
  235  18800 9.544796916977648e+01 1.8e+03 6.02e-03  2e-12  2e-05 9:24.7
  236  18880 9.544796916997851e+01 1.9e+03 6.10e-03  2e-12  2e-05 9:27.0
  237  18960 9.544796917002867e+01 1.9e+03 6.11e-03  2e-12  1e-05 9:29.2
  238  19040 9.544796916959211e+01 1.7e+03 6.32e-03  2e-12  1e-05 9:31.7
  239  19120 9.544796917013431e+01 1.8e+03 5.44e-03  2e-12  1e-05 9:34.0
  240  19200 9.544796917045382e+01 1.7e+03 5.56e-03  2e-12  1e-05 9:36.3
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  241  19280 9.544796916993295e+01 1.9e+03 6.05e-03  2e-12  1e-05 9:38.6
  242  19360 9.544796917029682e+01 1.9e+03 6.09e-03  2e-12  1e-05 9:40.9
  243  19440 9.544796917045706e+01 2.0e+03 6.35e-03  2e-12  1e-05 9:43.1
  244  19520 9.544796917050580e+01 1.9e+03 6.41e-03  2e-12  1e-05 9:45.4
  245  19600 9.544796917015645e+01 1.9e+03 6.31e-03  2e-12  1e-05 9:47.7
  246  19680 9.544796917006951e+01 2.1e+03 6.08e-03  2e-12  1e-05 9:49.9
  247  19760 9.544796917041792e+01 2.1e+03 6.66e-03  2e-12  1e-05 9:52.2
  248  19840 9.544796917012597e+01 2.0e+03 7.44e-03  2e-12  1e-05 9:54.5
  249  19920 9.544796917024834e+01 1.9e+03 6.68e-03  2e-12  1e-05 9:56.7
  250  20000 9.544796916930370e+01 2.0e+03 7.17e-03  2e-12  1e-05 9:59.0
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  251  20080 9.544796917012796e+01 2.1e+03 6.73e-03  2e-12  1e-05 10:01.3
  252  20160 9.544796917028361e+01 2.1e+03 6.85e-03  2e-12  1e-05 10:03.7
  253  20240 9.544796916966614e+01 2.2e+03 7.27e-03  2e-12  1e-05 10:06.0
  254  20320 9.544796917014909e+01 2.4e+03 7.05e-03  2e-12  1e-05 10:08.2
  255  20400 9.544796916994169e+01 2.5e+03 6.22e-03  1e-12  1e-05 10:10.5
  256  20480 9.544796916997339e+01 2.4e+03 5.83e-03  1e-12  1e-05 10:12.8
  257  20560 9.544796917013541e+01 2.3e+03 5.68e-03  1e-12  9e-06 10:15.0
  258  20640 9.544796916982128e+01 2.3e+03 5.15e-03  1e-12  8e-06 10:17.3
  259  20720 9.544796916987592e+01 2.3e+03 4.78e-03  1e-12  8e-06 10:19.6
  260  20800 9.544796917009758e+01 2.4e+03 4.67e-03  1e-12  8e-06 10:22.3
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  261  20880 9.544796917030425e+01 2.4e+03 4.50e-03  9e-13  7e-06 10:24.5
  262  20960 9.544796917023950e+01 2.6e+03 4.64e-03  9e-13  7e-06 10:26.8
  263  21040 9.544796917018581e+01 2.4e+03 4.46e-03  9e-13  7e-06 10:29.0
  264  21120 9.544796917019329e+01 2.4e+03 5.03e-03  1e-12  8e-06 10:31.5
  265  21200 9.544796917016158e+01 2.5e+03 4.79e-03  1e-12  8e-06 10:34.5
  266  21280 9.544796917033813e+01 2.7e+03 4.41e-03  9e-13  7e-06 10:36.8
  267  21360 9.544796916992169e+01 2.7e+03 4.33e-03  9e-13  7e-06 10:39.2
  268  21440 9.544796917007821e+01 2.8e+03 4.24e-03  9e-13  7e-06 10:41.4
  269  21520 9.544796916944603e+01 3.1e+03 4.07e-03  8e-13  7e-06 10:43.7
  270  21600 9.544796916984072e+01 3.5e+03 4.51e-03  9e-13  9e-06 10:45.9
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  271  21680 9.544796916999702e+01 3.6e+03 4.10e-03  9e-13  7e-06 10:48.1
  272  21760 9.544796917006687e+01 3.7e+03 4.25e-03  9e-13  8e-06 10:50.3
  273  21840 9.544796917018598e+01 3.8e+03 4.02e-03  8e-13  8e-06 10:52.5
  274  21920 9.544796917015604e+01 3.8e+03 3.75e-03  7e-13  7e-06 10:54.8
  275  22000 9.544796917008921e+01 3.9e+03 3.74e-03  7e-13  7e-06 10:57.0
  276  22080 9.544796916906569e+01 3.6e+03 3.50e-03  6e-13  7e-06 10:59.2
  277  22160 9.544796917025076e+01 3.7e+03 3.75e-03  7e-13  7e-06 11:01.4
  278  22240 9.544796917004975e+01 3.8e+03 3.76e-03  7e-13  7e-06 11:03.7
  279  22320 9.544796917028614e+01 3.9e+03 3.37e-03  6e-13  6e-06 11:05.9
  280  22400 9.544796916996128e+01 3.3e+03 3.24e-03  6e-13  6e-06 11:08.1
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  281  22480 9.544796917013804e+01 3.5e+03 3.43e-03  7e-13  6e-06 11:10.4
  282  22560 9.544796916998361e+01 3.9e+03 3.60e-03  7e-13  6e-06 11:12.7
  283  22640 9.544796916994113e+01 3.9e+03 4.16e-03  9e-13  8e-06 11:15.8
  284  22720 9.544796916958111e+01 4.1e+03 4.51e-03  1e-12  8e-06 11:18.1
  285  22800 9.544796917036626e+01 4.1e+03 4.68e-03  1e-12  8e-06 11:20.4
  286  22880 9.544796917025661e+01 4.3e+03 4.83e-03  1e-12  8e-06 11:22.7
  287  22960 9.544796917033969e+01 4.7e+03 5.06e-03  1e-12  8e-06 11:24.9
  288  23040 9.544796916961417e+01 4.4e+03 5.03e-03  1e-12  8e-06 11:27.1
  289  23120 9.544796917029429e+01 5.0e+03 5.17e-03  1e-12  8e-06 11:29.8
  290  23200 9.544796917011851e+01 5.0e+03 5.82e-03  1e-12  9e-06 11:33.0
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  291  23280 9.544796917024874e+01 5.2e+03 6.16e-03  1e-12  9e-06 11:36.2
  292  23360 9.544796917013961e+01 5.3e+03 5.47e-03  1e-12  8e-06 11:39.1
  293  23440 9.544796916997230e+01 5.0e+03 5.33e-03  1e-12  8e-06 11:41.6
  294  23520 9.544796917032068e+01 5.4e+03 5.65e-03  9e-13  8e-06 11:44.4
  295  23600 9.544796917031672e+01 5.8e+03 5.55e-03  1e-12  9e-06 11:47.0
  296  23680 9.544796916994574e+01 6.3e+03 5.37e-03  9e-13  9e-06 11:49.5
  297  23760 9.544796916985520e+01 6.2e+03 5.53e-03  9e-13  9e-06 11:52.4
  298  23840 9.544796916985658e+01 5.7e+03 5.76e-03  1e-12  9e-06 11:54.9
  299  23920 9.544796916999950e+01 5.3e+03 5.73e-03  9e-13  8e-06 11:57.4
  300  24000 9.544796917023989e+01 5.0e+03 5.65e-03  9e-13  8e-06 12:00.1
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  301  24080 9.544796917018090e+01 5.5e+03 5.75e-03  1e-12  8e-06 12:02.5
  302  24160 9.544796916989364e+01 5.7e+03 6.50e-03  1e-12  9e-06 12:05.0
  303  24240 9.544796916958408e+01 5.5e+03 6.20e-03  1e-12  9e-06 12:07.4
  304  24320 9.544796917014592e+01 5.5e+03 6.37e-03  1e-12  8e-06 12:09.7
  305  24400 9.544796917029507e+01 5.8e+03 5.84e-03  9e-13  8e-06 12:12.0
  306  24480 9.544796917005941e+01 6.1e+03 6.23e-03  9e-13  9e-06 12:14.4
  307  24560 9.544796917032909e+01 6.6e+03 6.04e-03  9e-13  9e-06 12:16.6
  308  24640 9.544796916999459e+01 7.2e+03 6.09e-03  9e-13  9e-06 12:18.9
  309  24720 9.544796917032868e+01 6.2e+03 6.22e-03  9e-13  9e-06 12:21.1
  310  24800 9.544796917035669e+01 6.7e+03 6.49e-03  9e-13  9e-06 12:23.4
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  311  24880 9.544796917000335e+01 6.3e+03 6.70e-03  9e-13  1e-05 12:25.7
  312  24960 9.544796917027622e+01 6.8e+03 6.88e-03  9e-13  1e-05 12:27.9
  313  25040 9.544796916996108e+01 7.5e+03 6.77e-03  8e-13  1e-05 12:30.1
  314  25120 9.544796917023022e+01 7.2e+03 6.23e-03  8e-13  1e-05 12:32.3
  315  25200 9.544796917024436e+01 7.5e+03 6.24e-03  7e-13  1e-05 12:34.6
  316  25280 9.544796916995816e+01 7.6e+03 5.51e-03  6e-13  8e-06 12:36.8
  317  25360 9.544796916984515e+01 7.7e+03 5.33e-03  6e-13  8e-06 12:39.1
  318  25440 9.544796916981441e+01 7.8e+03 5.50e-03  6e-13  9e-06 12:41.3
  319  25520 9.544796917006263e+01 8.1e+03 6.14e-03  8e-13  1e-05 12:43.5
  320  25600 9.544796916967751e+01 8.0e+03 6.11e-03  8e-13  9e-06 12:45.8
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  321  25680 9.544796917028577e+01 6.8e+03 5.87e-03  8e-13  9e-06 12:48.0
  322  25760 9.544796916985506e+01 6.5e+03 5.79e-03  8e-13  8e-06 12:50.2
  323  25840 9.544796917004528e+01 7.0e+03 5.10e-03  7e-13  7e-06 12:52.4
  324  25920 9.544796917020452e+01 6.5e+03 5.26e-03  7e-13  7e-06 12:54.6
  325  26000 9.544796917004717e+01 6.8e+03 4.79e-03  6e-13  7e-06 12:56.8
  326  26080 9.544796917048677e+01 7.1e+03 4.83e-03  6e-13  7e-06 12:59.0
  327  26160 9.544796917011008e+01 7.2e+03 5.08e-03  7e-13  7e-06 13:01.2
  328  26240 9.544796916987733e+01 7.9e+03 4.88e-03  7e-13  7e-06 13:03.4
  329  26320 9.544796917011050e+01 8.3e+03 5.28e-03  7e-13  8e-06 13:05.6
  330  26400 9.544796916972791e+01 1.0e+04 5.52e-03  7e-13  8e-06 13:07.8
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  331  26480 9.544796917013407e+01 1.0e+04 6.17e-03  8e-13  9e-06 13:10.1
  332  26560 9.544796917049969e+01 1.1e+04 6.06e-03  8e-13  9e-06 13:12.3
  333  26640 9.544796917030517e+01 1.2e+04 6.02e-03  8e-13  9e-06 13:14.5
  334  26720 9.544796917018182e+01 1.2e+04 6.01e-03  7e-13  9e-06 13:16.8
  335  26800 9.544796917034536e+01 1.2e+04 6.13e-03  7e-13  8e-06 13:19.0
  336  26880 9.544796916970422e+01 1.1e+04 6.05e-03  6e-13  8e-06 13:21.2
  337  26960 9.544796917020686e+01 1.1e+04 6.93e-03  7e-13  9e-06 13:23.4
  338  27040 9.544796916997672e+01 1.1e+04 6.45e-03  7e-13  8e-06 13:25.7
  339  27120 9.544796917037250e+01 1.1e+04 6.27e-03  6e-13  8e-06 13:27.9
  340  27200 9.544796917027912e+01 1.2e+04 6.18e-03  6e-13  8e-06 13:30.1
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  341  27280 9.544796917021507e+01 1.2e+04 5.74e-03  6e-13  7e-06 13:32.3
  342  27360 9.544796916994882e+01 1.2e+04 6.13e-03  6e-13  7e-06 13:34.6
  343  27440 9.544796917001872e+01 1.1e+04 6.93e-03  6e-13  9e-06 13:36.8
  344  27520 9.544796917023683e+01 1.0e+04 6.78e-03  6e-13  8e-06 13:39.0
  345  27600 9.544796917010630e+01 9.7e+03 7.13e-03  8e-13  8e-06 13:41.2
  346  27680 9.544796917003815e+01 1.0e+04 7.47e-03  9e-13  8e-06 13:43.5
  347  27760 9.544796916978147e+01 1.0e+04 8.21e-03  9e-13  9e-06 13:45.7
  348  27840 9.544796917009501e+01 9.3e+03 7.71e-03  9e-13  9e-06 13:47.9
  349  27920 9.544796917028995e+01 9.6e+03 7.06e-03  8e-13  8e-06 13:50.2
  350  28000 9.544796917006519e+01 1.0e+04 7.22e-03  9e-13  8e-06 13:52.4
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  351  28080 9.544796917000703e+01 9.7e+03 7.32e-03  9e-13  8e-06 13:54.6
  352  28160 9.544796917000086e+01 9.8e+03 6.87e-03  9e-13  7e-06 13:56.8
  353  28240 9.544796916995516e+01 9.3e+03 6.83e-03  8e-13  7e-06 13:59.1
  354  28320 9.544796916997841e+01 9.5e+03 6.65e-03  8e-13  7e-06 14:01.4
  355  28400 9.544796917001987e+01 9.7e+03 7.19e-03  9e-13  8e-06 14:03.7
  356  28480 9.544796917018229e+01 9.5e+03 7.08e-03  9e-13  7e-06 14:05.9
  357  28560 9.544796917025032e+01 9.2e+03 5.88e-03  7e-13  6e-06 14:08.1
  358  28640 9.544796917012576e+01 8.9e+03 5.37e-03  6e-13  5e-06 14:10.4
  359  28720 9.544796916980636e+01 8.8e+03 5.56e-03  6e-13  5e-06 14:12.6
  360  28800 9.544796916979719e+01 8.5e+03 5.79e-03  6e-13  5e-06 14:14.9
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  361  28880 9.544796917035204e+01 8.1e+03 6.03e-03  6e-13  6e-06 14:17.1
  362  28960 9.544796917031523e+01 8.4e+03 6.52e-03  6e-13  6e-06 14:19.4
  363  29040 9.544796916998141e+01 8.5e+03 6.66e-03  6e-13  7e-06 14:21.6
  364  29120 9.544796917032073e+01 9.1e+03 6.49e-03  7e-13  7e-06 14:23.8
  365  29200 9.544796916977005e+01 9.4e+03 6.82e-03  7e-13  7e-06 14:26.0
  366  29280 9.544796916947371e+01 9.1e+03 7.35e-03  8e-13  7e-06 14:28.2
  367  29360 9.544796917008826e+01 9.0e+03 7.31e-03  8e-13  8e-06 14:30.5
  368  29440 9.544796917020811e+01 1.0e+04 7.32e-03  8e-13  8e-06 14:32.7
  369  29520 9.544796917001004e+01 1.0e+04 7.76e-03  8e-13  8e-06 14:34.9
  370  29600 9.544796917010345e+01 1.0e+04 7.56e-03  8e-13  8e-06 14:37.1
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  371  29680 9.544796917005257e+01 1.1e+04 6.96e-03  7e-13  7e-06 14:39.4
  372  29760 9.544796917028543e+01 1.2e+04 7.02e-03  7e-13  6e-06 14:41.6
  373  29840 9.544796916981696e+01 1.1e+04 6.68e-03  8e-13  6e-06 14:43.8
  374  29920 9.544796917018637e+01 1.1e+04 6.87e-03  7e-13  6e-06 14:46.0
  375  30000 9.544796916984124e+01 1.1e+04 6.88e-03  7e-13  6e-06 14:48.3
  376  30080 9.544796917032006e+01 1.1e+04 6.17e-03  7e-13  5e-06 14:50.5
  377  30160 9.544796916997602e+01 1.1e+04 6.18e-03  7e-13  5e-06 14:52.7
  378  30240 9.544796917011716e+01 1.1e+04 6.13e-03  7e-13  5e-06 14:54.9
  379  30320 9.544796916978025e+01 1.1e+04 6.19e-03  6e-13  5e-06 14:57.1
  380  30400 9.544796917030416e+01 1.2e+04 6.32e-03  7e-13  5e-06 14:59.3
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  381  30480 9.544796917040364e+01 1.2e+04 6.56e-03  7e-13  5e-06 15:01.5
  382  30560 9.544796917028145e+01 1.2e+04 6.54e-03  8e-13  5e-06 15:03.7
  383  30640 9.544796917013025e+01 1.2e+04 7.11e-03  9e-13  6e-06 15:05.9
  384  30720 9.544796917000109e+01 1.2e+04 7.26e-03  8e-13  5e-06 15:08.1
  385  30800 9.544796917041110e+01 1.2e+04 7.33e-03  8e-13  5e-06 15:10.3
  386  30880 9.544796916988921e+01 1.3e+04 7.24e-03  8e-13  6e-06 15:12.6
  387  30960 9.544796917030899e+01 1.3e+04 7.13e-03  9e-13  5e-06 15:14.8
  388  31040 9.544796917051187e+01 1.4e+04 7.32e-03  1e-12  5e-06 15:17.1
  389  31120 9.544796917015054e+01 1.4e+04 6.83e-03  9e-13  5e-06 15:19.3
  390  31200 9.544796916976266e+01 1.3e+04 6.31e-03  9e-13  4e-06 15:21.6
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  391  31280 9.544796917009887e+01 1.3e+04 6.06e-03  8e-13  4e-06 15:23.8
  392  31360 9.544796917003816e+01 1.3e+04 6.04e-03  8e-13  4e-06 15:26.1
  393  31440 9.544796917039265e+01 1.3e+04 6.05e-03  8e-13  4e-06 15:28.3
  394  31520 9.544796917031948e+01 1.3e+04 5.97e-03  7e-13  4e-06 15:30.6
  395  31600 9.544796916983522e+01 1.2e+04 6.41e-03  8e-13  5e-06 15:32.8
  396  31680 9.544796917029328e+01 1.3e+04 6.27e-03  7e-13  5e-06 15:35.0
  397  31760 9.544796917012637e+01 1.5e+04 7.00e-03  8e-13  6e-06 15:37.3
  398  31840 9.544796916982365e+01 1.6e+04 7.22e-03  8e-13  7e-06 15:39.6
  399  31920 9.544796917015366e+01 1.6e+04 7.54e-03  9e-13  7e-06 15:42.0
  400  32000 9.544796917012960e+01 1.6e+04 6.65e-03  8e-13  6e-06 15:44.3
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  401  32080 9.544796917009278e+01 1.7e+04 6.23e-03  7e-13  6e-06 15:46.6
  402  32160 9.544796916994616e+01 1.6e+04 6.36e-03  8e-13  6e-06 15:48.9
  403  32240 9.544796916961261e+01 1.7e+04 6.90e-03  8e-13  6e-06 15:51.2
  404  32320 9.544796916998747e+01 1.6e+04 6.58e-03  7e-13  6e-06 15:53.5
  405  32400 9.544796917002222e+01 1.5e+04 6.74e-03  7e-13  6e-06 15:55.8
  406  32480 9.544796916951725e+01 1.6e+04 6.18e-03  7e-13  5e-06 15:58.0
  407  32560 9.544796916996557e+01 1.5e+04 6.21e-03  8e-13  5e-06 16:00.2
  408  32640 9.544796917016025e+01 1.6e+04 5.50e-03  6e-13  4e-06 16:02.4
  409  32720 9.544796916987369e+01 1.5e+04 5.20e-03  6e-13  4e-06 16:04.7
  410  32800 9.544796916977532e+01 1.4e+04 5.53e-03  6e-13  4e-06 16:06.9
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  411  32880 9.544796917029241e+01 1.4e+04 5.58e-03  7e-13  4e-06 16:09.3
  412  32960 9.544796916992203e+01 1.4e+04 5.03e-03  6e-13  3e-06 16:11.5
  413  33040 9.544796916998349e+01 1.4e+04 5.08e-03  6e-13  4e-06 16:13.8
  414  33120 9.544796917016933e+01 1.4e+04 5.02e-03  6e-13  4e-06 16:16.0
  415  33200 9.544796917021904e+01 1.3e+04 4.67e-03  6e-13  3e-06 16:18.3
  416  33280 9.544796916987204e+01 1.3e+04 4.44e-03  5e-13  3e-06 16:20.5
  417  33360 9.544796917003967e+01 1.2e+04 4.30e-03  5e-13  3e-06 16:22.7
  418  33440 9.544796917018499e+01 1.2e+04 4.59e-03  5e-13  3e-06 16:25.0
  419  33520 9.544796917010500e+01 1.2e+04 4.61e-03  5e-13  3e-06 16:27.2
  420  33600 9.544796917007980e+01 1.2e+04 4.58e-03  5e-13  3e-06 16:29.5
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  421  33680 9.544796916991817e+01 1.1e+04 4.54e-03  5e-13  3e-06 16:31.7
  422  33760 9.544796917043186e+01 9.8e+03 4.24e-03  5e-13  3e-06 16:33.9
  423  33840 9.544796916993198e+01 1.0e+04 4.29e-03  4e-13  3e-06 16:36.2
  424  33920 9.544796917005901e+01 1.0e+04 4.01e-03  4e-13  3e-06 16:38.4
  425  34000 9.544796917019303e+01 1.1e+04 4.02e-03  5e-13  3e-06 16:40.7
  426  34080 9.544796917028930e+01 1.0e+04 3.86e-03  4e-13  3e-06 16:42.9
  427  34160 9.544796917010180e+01 1.1e+04 3.87e-03  4e-13  3e-06 16:45.1
  428  34240 9.544796917006084e+01 1.0e+04 3.73e-03  4e-13  2e-06 16:47.4
  429  34320 9.544796916999157e+01 9.9e+03 3.92e-03  4e-13  2e-06 16:49.6
  430  34400 9.544796917024880e+01 9.4e+03 3.86e-03  4e-13  2e-06 16:51.9
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  431  34480 9.544796916989894e+01 9.1e+03 3.80e-03  4e-13  2e-06 16:54.1
  432  34560 9.544796916977339e+01 8.8e+03 3.60e-03  4e-13  2e-06 16:56.4
  433  34640 9.544796917044506e+01 8.5e+03 3.55e-03  3e-13  2e-06 16:58.7
  434  34720 9.544796916997903e+01 8.8e+03 3.70e-03  3e-13  2e-06 17:01.0
  435  34800 9.544796916964717e+01 8.7e+03 3.82e-03  3e-13  2e-06 17:03.2
  436  34880 9.544796917034688e+01 8.7e+03 3.54e-03  3e-13  2e-06 17:05.5
  437  34960 9.544796916989355e+01 8.5e+03 3.73e-03  3e-13  2e-06 17:07.7
  438  35040 9.544796917023861e+01 8.9e+03 3.70e-03  3e-13  2e-06 17:09.9
  439  35120 9.544796917017527e+01 9.3e+03 3.65e-03  3e-13  2e-06 17:12.2
  440  35200 9.544796917024971e+01 9.5e+03 3.60e-03  3e-13  2e-06 17:14.4
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  441  35280 9.544796917002937e+01 1.0e+04 3.87e-03  3e-13  2e-06 17:16.6
  442  35360 9.544796917006575e+01 1.0e+04 3.96e-03  3e-13  2e-06 17:18.9
  443  35440 9.544796916949163e+01 1.0e+04 3.78e-03  3e-13  2e-06 17:21.1
  444  35520 9.544796916993478e+01 1.0e+04 3.68e-03  3e-13  2e-06 17:23.3
  445  35600 9.544796916961394e+01 1.0e+04 3.67e-03  3e-13  2e-06 17:25.6
  446  35680 9.544796917008267e+01 1.0e+04 3.98e-03  3e-13  2e-06 17:27.8
  447  35760 9.544796917022605e+01 8.8e+03 4.01e-03  3e-13  2e-06 17:30.0
  448  35840 9.544796917031321e+01 1.1e+04 4.92e-03  4e-13  3e-06 17:32.2
  449  35920 9.544796916964702e+01 1.1e+04 4.73e-03  4e-13  3e-06 17:34.5
  450  36000 9.544796916933694e+01 1.2e+04 4.75e-03  4e-13  2e-06 17:36.7
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  451  36080 9.544796917009766e+01 1.2e+04 5.26e-03  5e-13  3e-06 17:38.9
  452  36160 9.544796916997770e+01 1.3e+04 4.92e-03  4e-13  3e-06 17:41.1
  453  36240 9.544796917014062e+01 1.3e+04 4.53e-03  4e-13  2e-06 17:43.3
  454  36320 9.544796917041164e+01 1.3e+04 4.45e-03  3e-13  2e-06 17:45.6
  455  36400 9.544796916994105e+01 1.2e+04 4.55e-03  4e-13  2e-06 17:47.8
  456  36480 9.544796917003238e+01 1.2e+04 4.47e-03  3e-13  2e-06 17:50.0
  457  36560 9.544796917005421e+01 1.3e+04 4.78e-03  4e-13  2e-06 17:52.2
  458  36640 9.544796917005186e+01 1.5e+04 4.43e-03  3e-13  2e-06 17:54.4
  459  36720 9.544796917027891e+01 1.5e+04 3.98e-03  3e-13  2e-06 17:56.7
  460  36800 9.544796917017142e+01 1.5e+04 3.66e-03  3e-13  2e-06 17:58.9
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  461  36880 9.544796917026810e+01 1.5e+04 3.39e-03  2e-13  2e-06 18:01.1
  462  36960 9.544796917041357e+01 1.4e+04 3.32e-03  2e-13  1e-06 18:03.4
  463  37040 9.544796917009016e+01 1.3e+04 3.15e-03  2e-13  1e-06 18:05.6
  464  37120 9.544796916974020e+01 1.4e+04 3.06e-03  2e-13  1e-06 18:07.9
  465  37200 9.544796917003876e+01 1.3e+04 2.94e-03  2e-13  1e-06 18:10.2
  466  37280 9.544796917027783e+01 1.3e+04 2.71e-03  2e-13  1e-06 18:12.4
  467  37360 9.544796917011263e+01 1.2e+04 2.71e-03  2e-13  1e-06 18:14.7
  468  37440 9.544796916998484e+01 1.0e+04 2.86e-03  2e-13  1e-06 18:17.0
  469  37520 9.544796916959213e+01 1.1e+04 2.79e-03  2e-13  1e-06 18:19.3
  470  37600 9.544796916994184e+01 1.1e+04 3.10e-03  2e-13  1e-06 18:21.6
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  471  37680 9.544796917026784e+01 1.2e+04 3.01e-03  2e-13  1e-06 18:23.8
  472  37760 9.544796917017702e+01 1.2e+04 3.32e-03  3e-13  1e-06 18:26.1
  473  37840 9.544796917003799e+01 1.3e+04 3.18e-03  2e-13  1e-06 18:28.3
  474  37920 9.544796917011138e+01 1.4e+04 3.36e-03  3e-13  1e-06 18:30.6
  475  38000 9.544796917011405e+01 1.4e+04 3.36e-03  2e-13  1e-06 18:32.8
  476  38080 9.544796916994113e+01 1.4e+04 4.01e-03  3e-13  2e-06 18:35.1
  477  38160 9.544796917013224e+01 1.5e+04 4.08e-03  3e-13  2e-06 18:37.5
  478  38240 9.544796916999914e+01 1.5e+04 3.81e-03  2e-13  2e-06 18:39.9
  479  38320 9.544796917009589e+01 1.6e+04 3.42e-03  2e-13  2e-06 18:42.2
  480  38400 9.544796916995095e+01 1.9e+04 3.48e-03  2e-13  2e-06 18:44.4
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  481  38480 9.544796917028201e+01 1.8e+04 3.68e-03  2e-13  2e-06 18:46.7
  482  38560 9.544796917011877e+01 1.8e+04 3.86e-03  2e-13  2e-06 18:48.9
  483  38640 9.544796916991524e+01 2.0e+04 3.65e-03  2e-13  2e-06 18:51.1
  484  38720 9.544796916996459e+01 2.0e+04 3.82e-03  2e-13  2e-06 18:53.4
  485  38800 9.544796916989593e+01 2.3e+04 3.25e-03  2e-13  2e-06 18:55.6
  486  38880 9.544796916974974e+01 2.2e+04 3.12e-03  2e-13  2e-06 18:57.8
  487  38960 9.544796917030143e+01 2.4e+04 2.83e-03  2e-13  1e-06 19:00.1
  488  39040 9.544796917015505e+01 2.4e+04 2.67e-03  2e-13  1e-06 19:02.3
  489  39120 9.544796917034589e+01 2.5e+04 2.78e-03  2e-13  1e-06 19:04.5
  490  39200 9.544796916977772e+01 2.7e+04 3.13e-03  2e-13  2e-06 19:06.7
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  491  39280 9.544796917013640e+01 2.8e+04 3.38e-03  2e-13  2e-06 19:08.9
  492  39360 9.544796916992267e+01 2.9e+04 3.43e-03  3e-13  2e-06 19:11.2
  493  39440 9.544796917019700e+01 2.9e+04 3.01e-03  2e-13  2e-06 19:13.4
  494  39520 9.544796916997215e+01 2.7e+04 3.19e-03  3e-13  2e-06 19:15.6
  495  39600 9.544796916994123e+01 2.9e+04 3.25e-03  3e-13  2e-06 19:17.8
  496  39680 9.544796917014492e+01 3.1e+04 3.28e-03  3e-13  2e-06 19:20.1
  497  39760 9.544796916975113e+01 3.3e+04 3.62e-03  3e-13  2e-06 19:22.3
  498  39840 9.544796917019522e+01 3.5e+04 3.38e-03  3e-13  2e-06 19:24.5
  499  39920 9.544796916994200e+01 3.7e+04 3.67e-03  3e-13  2e-06 19:26.7
  500  40000 9.544796916992831e+01 3.9e+04 3.60e-03  3e-13  2e-06 19:28.9
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  501  40080 9.544796917004253e+01 4.0e+04 3.42e-03  3e-13  2e-06 19:31.1
  502  40160 9.544796917019359e+01 4.6e+04 3.51e-03  3e-13  2e-06 19:33.3
  503  40240 9.544796917000780e+01 4.7e+04 3.02e-03  2e-13  2e-06 19:35.5
  504  40320 9.544796916931354e+01 5.0e+04 3.08e-03  2e-13  2e-06 19:37.7
  505  40400 9.544796917028654e+01 5.2e+04 2.97e-03  2e-13  2e-06 19:39.9
  506  40480 9.544796917019283e+01 5.3e+04 2.96e-03  2e-13  2e-06 19:42.2
  507  40560 9.544796917016532e+01 5.2e+04 3.08e-03  2e-13  2e-06 19:44.4
  508  40640 9.544796917010268e+01 4.9e+04 2.63e-03  2e-13  2e-06 19:46.7
  509  40720 9.544796916988393e+01 5.0e+04 2.54e-03  2e-13  1e-06 19:48.9
  510  40800 9.544796917015755e+01 5.2e+04 2.35e-03  2e-13  1e-06 19:51.1
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  511  40880 9.544796917003221e+01 6.0e+04 2.21e-03  1e-13  1e-06 19:53.4
  512  40960 9.544796917011317e+01 6.5e+04 2.21e-03  1e-13  1e-06 19:55.6
  513  41040 9.544796917013237e+01 6.4e+04 2.21e-03  2e-13  1e-06 19:57.9
  514  41120 9.544796917000532e+01 6.5e+04 2.31e-03  2e-13  1e-06 20:00.1
  515  41200 9.544796917022843e+01 6.1e+04 2.46e-03  2e-13  1e-06 20:02.4
  516  41280 9.544796916994326e+01 5.7e+04 2.69e-03  2e-13  2e-06 20:04.7
  517  41360 9.544796916999721e+01 5.9e+04 2.55e-03  2e-13  2e-06 20:06.9
  518  41440 9.544796916992452e+01 5.8e+04 2.76e-03  2e-13  2e-06 20:09.1
  519  41520 9.544796917028334e+01 6.4e+04 2.67e-03  2e-13  2e-06 20:11.4
  520  41600 9.544796917027058e+01 7.2e+04 2.54e-03  2e-13  2e-06 20:13.6
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  521  41680 9.544796917024803e+01 7.2e+04 2.42e-03  2e-13  2e-06 20:15.8
  522  41760 9.544796916991314e+01 7.9e+04 2.47e-03  2e-13  2e-06 20:18.0
  523  41840 9.544796917009819e+01 8.5e+04 2.31e-03  2e-13  2e-06 20:20.2
  524  41920 9.544796917010183e+01 8.8e+04 2.49e-03  2e-13  2e-06 20:22.5
  525  42000 9.544796916985860e+01 9.9e+04 2.62e-03  2e-13  2e-06 20:24.7
  526  42080 9.544796917013574e+01 9.6e+04 2.75e-03  2e-13  2e-06 20:26.9
  527  42160 9.544796916942852e+01 1.0e+05 3.00e-03  2e-13  2e-06 20:29.1
  528  42240 9.544796916990013e+01 1.1e+05 3.14e-03  2e-13  2e-06 20:31.3
  529  42320 9.544796916979134e+01 1.1e+05 2.95e-03  2e-13  2e-06 20:33.5
  530  42400 9.544796916989316e+01 1.2e+05 3.19e-03  2e-13  2e-06 20:35.8
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  531  42480 9.544796916994551e+01 1.1e+05 3.07e-03  2e-13  2e-06 20:38.0
  532  42560 9.544796917024908e+01 1.2e+05 3.29e-03  2e-13  2e-06 20:40.2
  533  42640 9.544796916994399e+01 1.2e+05 3.18e-03  2e-13  2e-06 20:42.4
  534  42720 9.544796917026252e+01 1.2e+05 3.33e-03  2e-13  2e-06 20:44.7
  535  42800 9.544796916988592e+01 1.1e+05 3.31e-03  2e-13  2e-06 20:46.9
  536  42880 9.544796917018947e+01 1.1e+05 3.56e-03  2e-13  2e-06 20:49.1
  537  42960 9.544796916989110e+01 9.5e+04 3.99e-03  3e-13  3e-06 20:51.3
  538  43040 9.544796916956244e+01 9.7e+04 4.13e-03  3e-13  3e-06 20:53.5
  539  43120 9.544796916986314e+01 9.8e+04 4.02e-03  2e-13  3e-06 20:55.7
  540  43200 9.544796917021529e+01 1.1e+05 4.33e-03  3e-13  3e-06 20:58.0
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  541  43280 9.544796916975477e+01 1.2e+05 3.98e-03  2e-13  3e-06 21:00.2
  542  43360 9.544796916932853e+01 1.3e+05 4.35e-03  2e-13  3e-06 21:02.4
  543  43440 9.544796917027708e+01 1.3e+05 4.43e-03  2e-13  3e-06 21:04.6
  544  43520 9.544796917008219e+01 1.3e+05 4.08e-03  2e-13  3e-06 21:06.8
  545  43600 9.544796917023938e+01 1.5e+05 4.12e-03  2e-13  3e-06 21:09.1
  546  43680 9.544796916977374e+01 1.6e+05 4.01e-03  2e-13  3e-06 21:11.3
  547  43760 9.544796917013318e+01 1.5e+05 3.90e-03  2e-13  3e-06 21:13.5
  548  43840 9.544796917021547e+01 1.5e+05 4.05e-03  2e-13  3e-06 21:15.7
  549  43920 9.544796916987951e+01 1.6e+05 4.04e-03  2e-13  3e-06 21:17.9
  550  44000 9.544796916983319e+01 1.5e+05 3.59e-03  2e-13  2e-06 21:20.3
Iterat #Fevals   function value  axis ratio  sigma  min&max std  t[m:s]
  551  44080 9.544796917014649e+01 1.5e+05 3.46e-03  2e-13  2e-06 21:22.5
  552  44160 9.544796917022693e+01 1.5e+05 3.59e-03  2e-13  2e-06 21:24.7
  553  44240 9.544796917016475e+01 1.4e+05 3.75e-03  2e-13  2e-06 21:27.0
  554  44320 9.544796917013051e+01 1.3e+05 3.75e-03  2e-13  2e-06 21:29.2
  555  44400 9.544796917018742e+01 1.2e+05 3.46e-03  2e-13  2e-06 21:31.4

=== FINAL ℤ₂₁₆₀ TRIADIC RESULT ===
Best χ²+reg = 90.03536796431972
Best parameters: [ 3.81144725e-02  3.03393288e-04  4.09963393e-02 -3.93534335e-02
 -3.63288591e-02  1.03884921e+00  5.10156644e-02  1.57762744e-01
  1.11508982e+00]

=== OBSERVABLES AT THIS POINT ===
m_c/m_t     : model= 2.330242e-05, target= 7.000000e-03, pull=-3.322
m_u/m_t     : model= 1.280684e-05, target= 1.000000e-05, pull= 0.936
m_s/m_b     : model= 3.232274e-06, target= 2.000000e-02, pull=-3.333
m_d/m_b     : model= 1.768955e-06, target= 1.000000e-03, pull=-3.327
m_mu/m_tau  : model= 9.062299e-06, target= 6.000000e-02, pull=-3.333
m_e/m_tau   : model= 4.956381e-06, target= 3.000000e-04, pull=-3.278
theta12_q   : model= 1.911628e-05, target= 2.260000e-01, pull=-3.333
theta23_q   : model= 4.748858e-02, target= 4.100000e-02, pull= 0.528
theta13_q   : model= 4.034241e-03, target= 3.500000e-03, pull= 0.509
theta12_l   : model= 5.796739e-01, target= 5.900000e-01, pull=-0.058
theta23_l   : model= 9.191949e-01, target= 8.400000e-01, pull= 0.314
theta13_l   : model= 1.497356e-01, target= 1.500000e-01, pull=-0.006
Delta_m2_21 : model= 2.215240e-24, target= 7.400000e-05, pull=-3.333
Delta_m2_31 : model= 1.065260e-21, target= 2.500000e-03, pull=-3.333

χ² (obs only) = 89.918
χ² + reg      = 90.035

|V_CKM| ≈
[[9.99991862e-01 1.91161219e-05 4.03422977e-03]
 [2.10602579e-04 9.98872625e-01 4.74703488e-02]
 [4.02877424e-03 4.74708121e-02 9.98864501e-01]]

|U_PMNS| ≈
[[0.82727976 0.5416221  0.1491767 ]
 [0.23929821 0.56973642 0.78621675]
 [0.50827607 0.61809862 0.59967453]]

"""