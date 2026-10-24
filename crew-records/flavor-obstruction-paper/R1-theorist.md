# R1 — theorist review of AUDIT-theta23-paper.md (item #4)

Read only. Paper: `theta23_obstruction_paper.tex` 147-470 (Flavor worktree). Audit: `AUDIT-theta23-paper.md`.
Every verdict below is derived on paper first; I then read the capture. I ran no numbers. "Hand check" = arithmetic done by me with the capture's own inputs.
Notation: Bloch form X = c + (δ/2) n·σ on E; a = δe/2, b = δν/2, φ = 2(θ+ψ₀), p = c_e c_ν, q = ab.

## 1. The statement
Check, claim by claim, the audit's verdicts (i), (ii), (ii-A), (iii), (iv), (v) against the paper's setup (Prop 1, Thm single, Thm nogo, Cor external). Alignment-side home of the paper's object: `content/framework/flavor-alignment.md` Thm 1–3 (𝓔(θ) = ‖[H_e, R H_ν R†]‖²; Thm 3 already restated in v0.20.1 as fit).

## 2. Assumptions (paper-internal; none are canon claims of the Alignment theory)
- A1 E is a 2-dim subspace; H_ν|_E = λ1 (Prop 1) or gap δν (Thm nogo). Paper §2 (tex 163-174). In canon: flavor-alignment.md Thm 1 premise.
- A2 R(θ) is a real SO(2) rotation on E; both Bloch vectors lie in the xz-plane (real eigenbases). Paper §2 (tex 167-171), used silently in Thm single's proof (tex 288-296). Not stated as an assumption.
- A3 Functionals are of (H_e|_E, H_ν(θ)|_E), the 2×2 compressions (Thm single, tex 266). The audit's (i) leaves this scope.
- A4 H_e, H_ν are PSD (Y_eY_e†, m_νm_ν†). Needed in (iv).

## 3. The derivation / verdicts
**(i) 3×3, E an eigenblock of H_ν — AGREE on the commutator; OVERSTATED on "Tr-type hold".**
Write H_e = [[a₁, u†],[u, H_E]] (u ∈ ℂ² the e-row), H_ν = m₁² ⊕ H_ν|_E, K = H_ν − m₁² = 0 ⊕ K_E. Then [H_e,H_ν] = [H_e,K] = [[0, u†K_E],[−K_E u, [H_E,K_E]]] (algebra). Hence
‖[H_e,H_ν]‖² = ‖[H_E,K_E]‖² + 2u†K_E²u.
K_E² = (c_ν−m₁²)² + b² + 2(c_ν−m₁²) b m(θ)·σ, so 2u†K_E²u = const + 4(c_ν−m₁²) b |u|² r_u cos(φ−φ_u), with φ_u the xz-angle of u's Bloch vector and r_u its xz-length. Total: 𝓔₃ = const + C sin²φ + B cos(φ−φ_u), B = 4(c_ν−m₁²) b |u|² r_u. A second phase φ_u (orientation of the e-row) enters ⇒ not one cosine, not reflection-symmetric about the ψ₀-point unless u ∥ an eigenvector of H_E (then B cosφ-type: still one cosine, but the stationary points move: ∂ gives cosφ* = ±B-dependent, not π/4). Audit formula u†(H_ν|_E − m₁²)²u: agreed.
- Size correction (disagree with the audit's remark "scales with the e-row coupling, ~0.15"): B ∝ |u|² and ∝ b (first order in δν); C ∝ δν² (second order). So as δν → 0 the e-row term dominates the C-term. In the paper's own regime (δν small) the 3×3 commutator norm is selected by u's orientation, not by ψ₀. Derived, not checked numerically.
- Tr(H_eH_ν) = a₁m₁² + Tr(H_E H_ν|_E): u absent, one cosine exactly. AGREE (capture floor 2e-7–7e-7).
- Tr((H_eH_ν)²): M = H_eH_ν has blocks M₁₂ = u†H_ν|_E, M₂₁ = m₁²u, so Tr M² = a₁²m₁⁴ + 2m₁² u†H_ν|_E u + Tr((H_E H_ν|_E)²). The middle term is a cos(φ−φ_u) at the k=2 frequency; the last term has k=2 and k=4 centred at ψ₀. A k=2 sinusoid re-phased by φ_u plus a k=4 piece centred elsewhere is not reflection-symmetric. So the audit's "hold" is **overstated**: capture shows asym/range 3.7e-6 (tilted default) and 3.6e-6 (CP) against a 2e-7–5e-7 floor, 300-draw max 2.9e-5 vs floor 7e-7 (`full3x3_fourier-run…txt` lines 27, 40, 78). That is a real, small effect (k=4 amp is only 3e-3 of range, so the asymmetry is suppressed), consistent with my term; it is not zero. In the θ13-only case (line 50) u is along an eigenvector of H_E, asymmetry sits at the floor: also consistent. Odd harmonics are absent either way (agree).
- Eigenframe misalignment f4 = 3 − Σ|W_ii|²: each |W_ii|² is a quadratic form in (cosθ, sinθ) of a fixed rotation, so exactly one k=2 sinusoid. AGREE (floor-level asym, k=4 amp 1e-16). Caveat: its phase is set by the 3×3 eigenframe of H_e, not by ψ₀ of the compression, so "θ₂₃ = ψ₀+θ" and its minima differ at O(tilt).
- Scope: the paper's Thm single is about H_e|_E, so (i) is a genuine extension, not a refutation of the theorem as stated.

**(ii) tilted degenerate E — AGREE.**
H_ν = m²·1 + Δ v₁v₁†, Δ = m₁²−m², v₁ = column 1 of U. Degenerate space = v₁^⊥; its flavor-basis 2-3 block is m²·1 + Δ w w†, w = (U_μ1,U_τ1): eigen-gap |Δ|(1−|U_e1|²), so ∝1 only if U_e1 = 1. Hand check with 33.4°, 8.6°: 1−c12²c13² = 0.3187 (capture 0.3186). U → U(1⊕R(t)) leaves v₁ and H_ν fixed, so P(μ→τ) is unchanged (agree). New columns: col₂' = c u₂ + s u₃, col₃' = −s u₂ + c u₃; so U_e3, U_e2 mix, θ13 and θ12 change (agree). δ=0: θ13=0 when tan t = u_e3/u_e2 = 0.1495/0.5443 → t = 15.35° (hand check, matches capture); θ12 then = arctan(√0.3187/√0.6813) = 34.36° (matches). Invariant |U_e2|²+|U_e3|² agrees.
- At m₂=m₃ all three angles are frame labels: AGREE in this sense: none of θ12, θ13, θ23 is individually invariant under the orbit; the invariants are |U_α1|. Overstatement: "any θ13 in the range" is limited to the orbit's range; and θ12 only moves by ~1°. "As a set unobservable" is fine; "θ13 is as unobservable as θ23" is true in the exact limit and wrong for any ε≠0 (the audit says this itself).
- OVERSTATED: "δ=197°: θ13 never 0". True only for real SO(2). The degenerate space carries all of U(2); a complex rotation puts col₃' ⟂ the conjugate e-row and zeroes U_e3 for any δ. The paper restricts itself to SO(2) (A2); the symmetry of H_ν|_E = λ1 is larger. This is a question for the paper, see §6 Q2.

**(ii-A) R12/R23 — AGREE, with one wording flag.**
U = R23·U13(δ)·R12 in PDG form, so U·R12(t) = R23·U13·R12(θ12+t): θ12 shifts by t and θ13, θ23, δ are exactly unchanged (the capture's 33.68→13.68 at t=20 is exactly this; hand check |U_e1|: cos13.68°·cos8.52° = 0.9609 ✓). U·R23(t) is not of that form, so θ12, θ13, θ23 all move. Flag: "near-degenerate pair" there means a small *Δm² ratio* (Δm²₂₁/|Δm²₃ℓ| = 0.0296: hand check 7.49e-5/2.534e-3 ✓), not mass-degenerate: for NO at m_lightest=0, m₂/m₁ is unbounded. The audit's own (iii) handles it correctly (relative mass splitting); the (ii-A) heading and "if the degeneracy that mattered were the data's own, the frame angle would be θ12" should say "Δm²-small", not "degenerate".

**(iii) thresholds — AGREE (hand-recomputed from `data/nufit6_central.csv`).**
NO: m₁²(1−(1−f)²) = (1−f)²Δ₃₁ − Δ₂₁. IC19 NO, Δ₃₁=2.534e-3, Δ₂₁=7.49e-5: f=10% → m₁ = 0.10202; f=5% → 0.15063; f=1% → 0.34791; Σ(10%) = 0.10202+0.10239+0.11376 = 0.31817 (capture 0.31818). IO literal (2-3), m₃²(1−(1−f)²) = (1−f)²|Δ₃₂|: f=10%, 2.510e-3 → m₃ = 0.10344, Σ = 0.33299. IO (1,2) pair: m₂² = Δ₂₁/(1−(1−f)²): f=1% → m₃ = 0.03541, Σ = 0.15750. At m_lightest=0: Σ_NO = 0.058993, Σ_IO = 0.09945, (m₂−m₁)/m₂ = 1.503%. All match captures. Conventions right: NuFIT Δm²₃ℓ with ℓ=1 (NO), ℓ=2 (IO, Δm²₃₂<0, m₃ lightest, m₁² = m₂² − Δ₂₁). Not checked: the CSV values themselves (copied from R1-literature §3.1, itself read from arXiv:2410.05380 v2) and the cosmology bounds (abstract level, R1-literature §3.2). Thresholds f = 1/5/10% are a choice, as the audit says.

**(iv) composites — AGREE on all formulas.**
With X = c+a n·σ, Y = c_ν+b m·σ: XYX = scalar + V·σ, V = 2a(c_e c_ν + ab cosφ) n + b(c_e²−a²) m (hand-derived; equals the audit's Bloch identity). ½{X,Y} = c_νX + c_eY − c_ec_ν + ab cosφ: in the linear class with θ-independent (α,β) = (c_ν,c_e). [Q₂,Y] = 2ib(V×m)·σ, V×m = 2a(c_ec_ν+ab cosφ)(n×m), so ‖[Q₂,H_ν]‖² ∝ sin²φ (p+q cosφ)². Zeros: sinφ=0 or cosφ = −p/q. PSD ⇒ c_e ≥ a, c_ν ≥ b ⇒ p/q ≥ 1 ⇒ |−p/q| ≥ 1: no interior zero (equality only at cosφ=−1, already aligned). AGREE. Interior max: x=cosφ, f=(1−x²)(p+qx)², f' = 2(p+qx)(q − px − 2qx²) ⇒ 2qx²+px−q=0 ⇒ x* = (−p+√(p²+8q²))/(4q). AGREE exactly. Hand checks: p/q=1 → x*=1/2, φ*=60° ✓; p/q=1.26 → 0.4591, 62.66° ✓; p/q=41.3 → ≈q/p=0.0242, 88.61° ✓. The other root is ≤ −1 for p ≥ q, so x* is the single interior stationary point, a max.
- Addition: Q₂'s θ-dependence is O(b): at δν → 0, V ∝ n and Q₂|_E = c_νH_e² on E, so Q₂ is inert in the paper's exact-degeneracy regime too.
- "Evades" is correct only against the paper's Prop-minimal class (fixed-coefficient linear). It does not contradict Thm single (‖[Q₂,H_ν]‖² is a G(cosφ)); it contradicts only the Remark's "stationary set frozen at {0,π/4,π/2}" (tex 415-416), which is true of one G, not of every relational functional. Thm single itself implies interior stationary points wherever G'(u)=0. That overstatement is in the paper, not the audit.
- Not run / no capture: full 3×3 composites, ‖[Q,H_e]‖² (the audit says so; analytic value ∝ sin²φ: I agree).

**(v) complex eigenbasis — AGREE.**
Hermitian H_E has Bloch vector n=(n_x,n_y,n_z); H_ν|_E real ⇒ m in xz-plane; real R(θ) rotates m only in that plane and leaves n_y fixed. n·m = r cosφ', r = √(n_x²+n_z²). Then Tr = 2c_ec_ν + ½δeδν r cosφ' and ‖[·,·]‖² = ½δe²δν²(1−r²cos²φ'), minimum ½δe²δν²(1−r²) > 0. Exactly the audit's formula. Further: if H_ν|_E is also complex (m_y≠0), n·m = n_y m_y + r r' cosφ' and the minimum is zero only if |n_y m_y|+r r' = 1. So the real-eigenbasis assumption (A2) is load-bearing for Thm monotone/nogo "minima at E=0".

**Claims with no capture:** (a) the audit's "size scales with the e-row coupling (~0.15)" — no capture varies the coupling; my derivation says B ∝ |u|² and ∝ δν (first order). (b) "Tr((HeHν)²) holds" — capture contradicts the exactness (above). (c) "δ read from |U_ij| only" — acknowledged; the δ=197° θ13-range statement therefore rests on the SO(2) restriction only. (d) NuFIT inputs are second-hand (CSV transcribed). (e) The gate's loop over `t` unused (audit self-reports). Everything else cites a capture.

## 4. Register
- Rigorous (derived on paper and reproduced numerically): 𝓔₃ = ‖[H_E,K_E]‖²+2u†K_E²u; Tr(H_eH_ν) one-cosine; f4 one-cosine; tilted-orbit θ12/θ13 change; thresholds; Q₂ Bloch identity, PSD p/q ≥ 1, x*; complex-eigenbasis formula.
- Rigorous on paper, numerics only indirect: Tr((HeHν)²) non-exactness (capture asym > floor, consistent; no direct fit of the u-term); B ∝ δν.
- Interpretive: "θ13 is as unobservable as θ23"; "obstruction is the ε→0 limit; the datum to import is the ε-term's eigenframe"; the ii-A reading that the data's relevant degeneracy is θ12.
Not flattened into each other.

## 5. EQUATIONS.md
Only the Flavor-alignment row (line 156; "rigorous · Thm 1–2 intact; Thm 3 corrected v0.20.1") touches. Current status stays correct. Proposed (no edit): none required for the registry row; a Reconciliations note could record that Thm 1–2 hold on the 2×2 compression with real Bloch vectors, and that 3×3 (e-row u) and complex-eigenbasis cases need extra input. Any change to the Thm 1–2 premise statement would be a flavor-alignment.md edit and a canon decision for Chazz.

## 6. Gaps (questions for Chazz)
- Q1 Is the 3×3 commutator (with the e-row u, E an eigenblock of H_ν) the object 𝓔 means, or only the compression 𝓔 on E? The answer decides whether the e-row term (first order in δν) is part of Thm nogo's scope.
- Q2 The symmetry of H_ν|_E = λ1 is U(2); R(θ) is only SO(2). Is the restriction to real rotations intended (A2)? It sets whether θ13 can be zeroed for δ≠0 and whether n_y components exist.
- Q3 Does the theory claim any role for E = a tilted degenerate eigenspace (rather than the flavor μτ plane) when θ12, θ13 ≠ 0? The paper's setup says "second and third generations".
- Q4 Remark (tex 415-416) "frozen at {0,π/4,π/2}": does the paper intend it for the commutator-norm family only? As written it reads for every relational functional.

## 7. Setup for the calculator
One falsifiable question: **for E an eigenblock of H_ν with a nonzero e-row u in H_e, does ‖[H_e,H_ν(θ)]‖² (3×3) equal const + C sin²φ + B cos(φ−φ_u) with B = 4(c_ν−m₁²)(δν/2)|u|² r_u, and does Tr((H_eH_ν)²) contain the extra 2m₁²u†H_ν|_E(θ)u?**
Ground truth: the 2×2 case u=0 (existing GT, B=0) must reproduce C sin²φ; and ≥ 2 values of |u| and δν.
Procedure: fit the Fourier k=2,4 content of the capture's f2 and f3 to the closed form, with u, c_ν, m₁² read from the matrices; vary δν over three decades with u fixed, and check B/C ∝ 1/δν.
Refuted if: the fitted B deviates from the closed form beyond ~1e-10 relative, B does not scale as |u|² and δν, or f3's asymmetry stays at the floor for u not ∥ an eigenvector of H_E.
Secondary (optional): U(2) orbit — confirm θ13 → 0 attainable for δ=197° with a complex rotation inside the degenerate space.
