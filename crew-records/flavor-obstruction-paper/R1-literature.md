# R1-literature — theta23 no-go paper: literature + version audit

Author: literature (read-only). Date 2026-10-23. Item #2 (R1, audit). **Revision 2** (lead's manager rulings applied: v1 = the registered docx via the extracted text; both NuFIT variants; cosmology as LCDM bounds with dataset named; Procesi/Sibirskii added; m_lightest arithmetic left to the data-analyst).
Tooling limit: no python/unzip/sed on this role. Web pages read via fetch; what I did and did not open is marked.

## 1. TL;DR
- Ask is **partly** answered inside the repo: the Thm 3 withdrawal is already in canon (`content/framework/flavor-alignment.md` §4, v0.20.1). The quasi-degenerate regime vs Δm²₃₁ and Σm_ν is nowhere in the repo, v1 or v2.
- v2 vs v1 (against the docx text): Thm 1 kept (demoted to a Prop); Thm 2 kept only as "the class we test" and its justification changed; Thm 3 **withdrawn** and replaced by a no-go; v1 §7 stability results, §3.6/§6.5/§8.3-8.4 numerics and the NuFIT-5.2 "near-maximal" numeric claim dropped; five items added (§2.4).
- NuFIT 6.0 (arXiv:2410.05380, v1 7 Oct 2024, v2 28 Oct 2025; JHEP 12 (2024) 216): both variants in §3.1 (IC19 w/o SK-atm and IC24 with SK-atm); θ₂₃ octant undetermined in both. Cosmology, LCDM: Planck 2018 TT,TE,EE+lowE+lensing+BAO <0.12 eV (see §3.2 caveat); DESI DR1+CMB <0.072 eV; DESI DR2+CMB <0.064 eV.
- m_lightest thresholds: formulas only (§3.3); data-analyst evaluates.
- Flag only: `content/vocabulary.md:166` still carries the withdrawn Thm 3 claim; lab crew fixes canon.

## 2. Inside the repo
Verdict: **partly**. (Repo = `/Users/chazzromeo/rAI/galt-crews/flavor-obstruction-paper/wt`, branch crew/flavor-obstruction-paper @ 9a1ff74; not the unmerged v0.27.0 branch.)

2.1 Canon
- `content/framework/flavor-alignment.md:33-51` Thm 1; `:53-70` Thm 2; `:72-121` Thm 3 + corrected status (lines 111-121: "the earlier 'the minimum sits at near-maximal mixing' reading of Theorem 3 is withdrawn"; E flat on exact degeneracy lines 89-94; code relaxes a circle-phase energy and scores θ₂₃ against 45°/49° target, lines 100-109).
- `content/EQUATIONS.md:156` row "Flavor alignment"; Reconciliations `:344-355`.
- `content/reference/sources/SOURCES.md:137-153` key `flavor-mixing-2026-07-02` = v1 only (docx; line 147-148 "numerical checks against NuFIT 5.2"). v2 not registered.
- `calculations/README.md:66` flavor-transport (complete, 2026-07-19), `:68` yukawa-hierarchy, `:74` triad-first-tests. Files: `calculations/flavor-transport/{README.md,CODE-AUDIT.md,run-2026-07-19.txt,test_e..test_i}`; Test G = `test_g_theta23.py`.
- Not found in repo: any Δm²₃₁ / Σm_ν / quasi-degenerate-regime check; any RG, Majorana-phase or 3+1 treatment of E.
2.2 Vocabulary rows: `content/vocabulary.md:164`, `:165`, `:166` (**stale — flag only**), `:153-154`.
2.3 Lab-crew results on the unmerged branch: only as summarised in `R1-audit.md`; not re-read by me.

2.4 v1 vs v2 — against the docx text
Sources. **v1** = the registered docx via `/Users/chazzromeo/rAI/galt-crews/flavor-obstruction/files/paper-v1-extracted.txt`; "P<n>" = paragraph/line n of that file. Caveat: the extraction **drops all display equations and several inline symbols** (e.g. P7, P263, P248-249 have blank formulas), so a v1 claim whose value was in an equation (the θ₂₃ number at P246, P248-249, P263) cannot be read from this file. **v2** = `/Users/chazzromeo/rAI/galt-crews/flavor-paper-wt/theta23_obstruction_paper.tex` ("L<n>"); **earlier draft** = `limitative_draft.tex` ("LD<n>").

| item | v1 (docx text) | v2 tex | status |
|---|---|---|---|
| Thm 1 residual SO(2) obstruction | Thm 1 P89-95: invariance of spectral triple, spectral action, and PMNS non-invariance; proof P96-105; "exact degeneracy sufficient but not necessary; approximate yields approximate SO(2)" P95; P124-125 "conditional but structural" | Prop `prop:obstruction` L188-197, proof L199-210; exact degeneracy only (L189); approximate handled via δ_ν in Thm nogo (L341) | **kept** (Thm -> Prop) |
| "θ12, θ13 unchanged, θ23 shifts by θ" | P81-83, P104, Cor 7.1 P343-345 ("first row invariant") | L195, L207-208 state "shifts θ23 by θ"; no corollary on θ12/θ13 | **kept as a statement, dropped as a corollary** (this is the item brief §4(ii) asks to test) |
| Thm 2 unique linear relational operator A=αH_e+βH_ν | Thm 2 P163-170 "No other operator satisfies these conditions"; proof P171-181; Step 2 (P175): nonlinear terms incl. H_eH_ν "modify eigenvalues"; P156 claims A commutes with π(a) | Prop `prop:minimal` L227-237; justification changed: H_eH_ν excluded as **non-Hermitian** (L234-237), and "We take Prop. as the **definition of the class we test**" (L236-237) | **weakened** (uniqueness theorem -> definition of the class tested; the reason for excluding H_eH_ν changed) |
| scope of uniqueness | P456: "claims of uniqueness ... restricted: among linear, relational operators ... simplest choice"; P132-133 "minimal relational postulate ... alternative criteria could be considered" | L138-144 "a no-go, not a prediction" | **kept the restriction, reversed the message** |
| A|_E selects a preferred μ–τ basis | P160, P177-178, P184-185: A|_E nondegenerate "selects preferred eigenvectors"; "unique orientation" | Lemma `lem:inertness` L239-256 (also LD54-77): A|_E = αH_e|_E + βλ·1; shares projectors with H_e; adopting it gives θ23 = 0 | **added; reverses v1's reading** |
| Thm 3 (unique global minimum, fixes θ23) | Thm 3 P257-264 "unique global minimum ... θ23 fixed uniquely"; proof P266-296, Step 4 P283-290: "θ* = 0 (mod π) is the unique global minimizer"; P301 "θ23 is a derived geometric quantity, not ... a phenomenological fit"; P192 "any functional vanishing iff simultaneously diagonalizable yields the same minimizer" | **withdrawn**; replaced by Thm `thm:nogo` L332-356: E = C sin² 2(θ+ψ0), C=½δ_e²δ_ν², stationary set {0, π/4, π/2}, minima at {0, π/2}, π/4 max (L347-350); abstract L74-90 states θ23 not determined | **withdrawn / replaced** |
| uniqueness of minimizer | P286 "unique ... θ*=0 (mod π)" | L347-350: minima at **both** {0, π/2} (period π/2 in θ); exact degeneracy gives E≡0 (L339, LD104-107) | v1's "mod π" uniqueness differs from v2's two minima (comparison of P286 vs L347-350; the v1 period claim depends on a formula lost in extraction) |
| numerical θ23 claim | P246 "unique minimum ... atmospheric angle in the vicinity of [value lost], well within 1σ"; P248-249 NuFIT 5.2, NO, "near-maximal" "robust under small variations" | absent; L108-111 only quotes NuFIT 6.0/PDG ranges as context; canon says the code scored against 45°/49° (flavor-alignment.md:100-109) | **withdrawn** |
| single-invariant reduction | not in v1 | Thm `thm:single` L265-298: Tr(H_e H_ν(θ)) = c + ½δ_eδ_ν cos 2(θ+ψ0) (L270) | **added** |
| aligned points always stationary | not in v1 | Cor `cor:aligned` L300-311 | **added** |
| monotone => exact alignment | not in v1 | Thm `thm:monotone` L313-330 | **added** |
| external input: Q with (a),(b) | not in v1; P433-437 and P458 only list "CP phases, RG stability, other multiplicity structures" as future directions | Cor `cor:external` L390-411 (also LD165-187); Remark L413-420; Discussion L454-466 names RG running of H_ν as the next test | **added** |
| §7 stability: spectra, hierarchies, Δm² unchanged; CKM/Jarlskog unchanged | §7 P328-377: Cor 7.1 (P343-345), Prop 7.2 (P355-359), Cor 7.3 (P368-369) | no such section | **dropped, not withdrawn** |
| numerical emergence of 1+2 block | §3.6 P115-125; §6.5 P315-327 | absent (App. A L470- is the verifier script only) | **dropped** |
| emergent internal geometry, frustrated 6-fold/5-fold phases | §8.3 P406-414; §8.4 P415-420 | absent | **dropped** |
| interpretation: "completion, not extension of NCG"; "resolves the final unfixed parameter" | P17, P382-383, P404, P446, P459 | L138-144, L426-466: "draws a boundary rather than crosses one"; near-maximal θ23 <-> approximate μ–τ symmetry conceded as standard (L429-431) | **reversed** |
| NuFIT version | NuFIT 5.2 NO (P248) | NuFit-6.0 `\bibitem{nufit6}` L537-540; PDG 2024 L533-535 | updated |
| NCG literature cited | Connes 1994, Chamseddine–Connes 1997, Devastato–Lizzi–Martinetti 2014 (P86), Barrett 2007 / van Suijlekom 2015 (P386) | v2 cites connes1994, connesCMP1996, cc1997, ccm2007, vandenDungen2012, vanSuijlekom2015 (L98-102) | changed |

Other v1 oddities relevant to a later audit (read in the text; not judged): P5-6 abstract says the functional "possesses a unique global minimum ... thereby fixing θ23"; P19-27 duplicate the paper outline (P18 and P19-27); P222 refers to "Figure X" (unresolved); P247 says the illustration is "not intended as a fit" while P246-249 give a specific value.
v2 intro L108-111 ("best-fit from ≈43° (lower octant) to ≈49° (upper octant), 3σ ≈ 41°–51°", `\cite{pdg,nufit6}`) vs NuFIT 6.0 Table 1: IC24 with SK-atm NO best fit 43.3°; IC19 without SK-atm NO 48.5°, IO 48.6°; IC24 with SK-atm IO 47.9°; 3σ ranges 41.0–50.6° across both variants (§3.1). v2 does not name a variant. v2 contains no Δm²₃₁ or Σm_ν value anywhere (lines 60-480 read).

## 3. Outside the repo
3.1 NuFIT 6.0 — OPENED (PDF pp. 1-8 of arXiv v2). Esteban, Gonzalez-Garcia, Maltoni, Martinez-Soler, Pinheiro, Schwetz, "NuFit-6.0: Updated global analysis of three-flavor neutrino oscillations", arXiv:2410.05380 [hep-ph]; v1 submitted 7 Oct 2024; **the PDF read is v2 dated 28 Oct 2025**; JHEP 12 (2024) 216; data "as of September 2024" (abstract). https://arxiv.org/abs/2410.05380
Release date of the nu-fit.org web page: **could not open** (connection refused); the arXiv v1 date is the only date cited.
Variants (the paper's own labels, p.3): "IC19 without SK atmospheric data" and "IC24 with SK atmospheric data" (tabulated χ² maps for SK-atm and IceCube/DeepCore 9.3-year). Δm²₃ℓ = Δm²₃₁ >0 for NO; = Δm²₃₂ <0 for IO. Table 1 (p.6), best fit ±1σ and 3σ range:

| quantity | IC19 w/o SK-atm, NO (Δχ²=0.6) | IC19 w/o SK-atm, IO (best) | IC24 with SK-atm, NO (best) | IC24 with SK-atm, IO (Δχ²=6.1) |
|---|---|---|---|---|
| θ23/° | 48.5 +0.7/−0.9 (3σ 41.0→50.5) | 48.6 +0.7/−0.9 (41.4→50.6) | 43.3 +1.0/−0.8 (41.3→49.9) | 47.9 +0.7/−0.9 (41.5→49.8) |
| sin²θ23 | 0.561 +0.012/−0.015 | 0.562 +0.012/−0.015 | 0.470 +0.017/−0.013 | 0.550 +0.012/−0.015 |
| δCP/° | 177 +19/−20 (96→422) | 285 +25/−28 (201→348) | 212 +26/−41 (124→364) | 274 +22/−25 (201→335) |
| Δm²₃ℓ/10⁻³ eV² | +2.534 +0.025/−0.023 (+2.463→+2.606) | −2.510 +0.024/−0.025 (−2.584→−2.438) | +2.513 +0.021/−0.019 (+2.451→+2.578) | −2.484 +0.020/−0.020 (−2.547→−2.421) |
| Δm²₂₁/10⁻⁵ eV² | 7.49 +0.19/−0.19 | 7.49 | 7.49 | 7.49 |
| θ13/° | 8.52 +0.11/−0.11 | 8.58 | 8.56 | 8.59 |
| θ12/° | 33.68 +0.73/−0.70 | 33.68 | 33.68 | 33.68 |
Octant (abstract): "θ₂₃ still suffers from the octant ambiguity, with no clear indication of whether it is larger or smaller than 45°". Fig. 1: sin²θ23 has two minima for NO (~0.47 and ~0.56) in both variants. Ordering: preference for NO with Δχ²=6.1 only when the SK-atm χ² table is added; otherwise "almost equally good fit for both orderings" (abstract). CP: NO consistent with CP conservation within 1σ; IO δCP≈270° favoured against CP conservation at >3.6σ (abstract; variant not named in the abstract).
Sec. 5 of the same paper (PDF pp.19-21: absolute-mass observables incl. cosmology) **not read**.

3.2 Cosmology, ΛCDM, 95% CL upper bound on Σm_ν (fetch/abstract-level reads; the full datasets are those the sources name):
- Planck 2018 results VI, Aghanim et al., arXiv:1807.06209, A&A 641 A6 (2020): "Planck + BAO: Σm_ν < 0.12 eV" (fetch result: "combining Planck data with BAO"). The detailed Planck likelihood combination (TT,TE,EE+lowE+lensing) is **not confirmed by me**; cite as "Planck 2018 + BAO". https://arxiv.org/abs/1807.06209
- DESI 2024 VI (DR1), Adame et al., arXiv:2404.03002: DESI BAO + CMB, flat ΛCDM, Σm_ν < 0.072 eV (the paper also gives a second, relaxed number in parentheses, 0.113 eV — its meaning is not read; do not use). https://arxiv.org/abs/2404.03002
- DESI DR2 Results II, DESI Collaboration, arXiv:2503.14738 (v1 18 Mar 2025, final 9 Oct 2025; PRD 112, 083515): DESI BAO + CMB, ΛCDM, Σm_ν < 0.064 eV. Which CMB set is **not read**. https://arxiv.org/abs/2503.14738
Caveat (one line, unquoted): non-ΛCDM extensions (e.g. w0wa dynamical dark energy) relax these bounds; I did not read those values.

3.3 m_lightest thresholds — formulas only; the data-analyst evaluates.
Conventions: NuFIT 6.0 §2: NO m1<m2<m3, Δm²₃ℓ=Δm²₃₁>0; IO m3<m1<m2, Δm²₃ℓ=Δm²₃₂<0.
- NO, lightest m1=m_l: m2=√(m_l²+Δm²₂₁), m3=√(m_l²+Δm²₃₁). Σ = m1+m2+m3.
- Exact condition (m3−m2)/m3 < f ⇔ (Δm²₃₁−Δm²₂₁)/(m3(m3+m2)) < f, with f = 0.01, 0.05, 0.10.
- Small-splitting approximation (my own algebra, not a citation): (m3−m2)/m3 ≈ (Δm²₃₁−Δm²₂₁)/(2m3²) ⇒ m3 > √((Δm²₃₁−Δm²₂₁)/(2f)); then m_l=√(m3²−Δm²₃₁). Use the exact form for evaluation.
- IO: the near-degenerate pair is (m1,m2) (splitting Δm²₂₁); m3 is lightest with m2=√(m3²+|Δm²₃₂|). Which pair plays "E" is ambiguous (§6.3).
- State survival per bound (§3.2) and per NuFIT variant (§3.1).
- Premise point (from `R1-audit.md:75-85` and §3.1 table values): the 2-3 block of H_ν is ∝𝟙 only for θ12=θ13=0 and m2=m3; NuFIT 6.0 (both variants) has θ13≈8.5°, θ12≈33.7°.

3.4 Standard anchors (what I read vs not)
- Harrison & Scott, "μ–τ reflection symmetry in lepton mixing and neutrino oscillations", Phys. Lett. B 547 (2002) 219-228, arXiv:hep-ph/0210197. Abstract read: "corresponding elements of the μ- and τ-neutrino flavour eigenstates have equal moduli", consistent with data. θ₂₃=45° and δ=±90° are **not in the abstract**; body not read; treat as to-be-verified ground truth. https://arxiv.org/abs/hep-ph/0210197
- Grimus & Lavoura, "A nonstandard CP transformation leading to maximal atmospheric neutrino mixing", Phys. Lett. B 579 (2004) 113-122, arXiv:hep-ph/0305309. Abstract-level: BMV matrix; atmospheric mixing exactly maximal, CP violation maximal in the standard phase convention; other angles and the spectrum unconstrained. https://arxiv.org/abs/hep-ph/0305309
- Casas, Espinosa, Ibarra, Navarro, "Naturalness of nearly degenerate neutrinos", Nucl. Phys. B 556 (1999) 3-22, arXiv:hep-ph/9904395. Abstract-level only: radiative corrections can generate small splittings; bimaximal mixing stable under RG in a see-saw framework; vacuum-oscillation solar solution excluded. https://arxiv.org/abs/hep-ph/9904395
- Casas, Espinosa, Ibarra, Navarro, "General RG equations for physical neutrino parameters and their phenomenological implications", Nucl. Phys. B 573 (2000) 652-684, arXiv:hep-ph/9910420. Abstract-level: RGEs for three masses, mixing angles and CP phases. The audit's "angles depend on CP parities" is **not verified**. https://arxiv.org/abs/hep-ph/9910420
- Antusch, Kersten, Lindner, Ratz, "Running neutrino masses, mixings and CP phases: analytical results and phenomenological consequences", Nucl. Phys. B 674 (2003) 401-433, arXiv:hep-ph/0305273. Abstract-level: analytic RG formulae for masses, angles, phases; compact for small θ13. Formulas for dθ23/dt not captured. https://arxiv.org/abs/hep-ph/0305273
- 3+1: Kopp, Machado, Maltoni, Schwetz, "Sterile neutrino oscillations: the global picture", JHEP 05 (2013) 050, arXiv:1303.3011. Abstract-level: 3+1/3+2/1+3+1 global fits, appearance–disappearance tension. This gives phenomenology, **not** the 4×4 parametrization; that statement is not read. https://arxiv.org/abs/1303.3011
- Invariant theory (new, per ruling): C. Procesi, "The invariant theory of n×n matrices", Adv. Math. 19 (1976) 306-381; K. S. Sibirskii, "Algebraic invariants of a set of matrices" (Russian), Sibirsk. Mat. Zh. 9 (1968) no. 1, 152-164 (English: Siberian Math. J. 9 (1968)). Read via a web search summary only: any polynomial invariant of r n×n matrices (char 0, GL(n)/O(n)/Sp(n)) is a polynomial in traces of non-commutative monomials tr(A_{i1}…A_{ik}) (Procesi). The explicit 2×2 two-matrix generating set (Tr A, Tr B, Tr A², Tr B², Tr AB) is the classical consequence that v2 uses (proof of Thm single, L278-282); **I did not open either paper**, so the n=2, r=2 list as stated is not verified from the source. https://www.ams.org/bull/1976-82-06/S0002-9904-1976-14196-1/ (Procesi's AMS Bulletin announcement, search result; not opened)
- v2's own cited reviews not opened: Xing & Zhao, Rept. Prog. Phys. 79 (2016) 076201, arXiv:1512.04207; King & Luhn, Rept. Prog. Phys. 76 (2013) 056201, arXiv:1301.1340; PDG 2024, Phys. Rev. D 110 (2024) 030001.

## 4. Standard counterparts
Ask terms -> mainstream counterpart (vocabulary.md maps rows :153-154, :164-166 only):
- E (2-D degenerate eigenspace of H_ν on μ–τ) -> degenerate-subspace freedom / quasi-degenerate m2≈m3. Mapped at :164 ("degenerate-eigenspace gauge freedom"). Lacks: "quasi-degenerate spectrum", "Majorana phases / CP parities".
- R(θ) SO(2) obstruction -> basis freedom in a degenerate eigenspace; observably, θ23 is unphysical at Δm²_atm=0. Lacks this pointer.
- A = αH_e+βH_ν -> :153-154 maps to a linear pencil.
- Tr(H_e H_ν) single invariant -> trace invariant of a pair of 2×2 Hermitians (Procesi/Sibirskii). Lacks.
- Q on E with (a),(b) -> μ–τ-breaking perturbation / RG-induced breaking (per v2's own refs). Lacks.
- θ23=π/4 forced -> μ–τ reflection / 2↔3 exchange symmetry (Harrison–Scott). Lacks.
- "RG scale", "Majorana-phase dimension", "4th-generation direction", "Cl(1,3)" (plan.md D1-D6) -> RG running; Majorana phases; 3+1 sterile; Dirac-algebra Cl(1,3): partly at :64 (sterile neutrinos as hidden alignments), rest lacking.
Proposed entries (not added): quasi-degenerate spectrum; μ–τ reflection symmetry; Majorana CP parity; RG running of θ23; 3+1 sterile mixing.

## 5. Known answers a method could be checked against
1. NuFIT 6.0 Table 1, both variants (§3.1; read).
2. v2's own identity ‖[A,B]‖_F² = ½δ_A²δ_B² sin²2(ξ_A−ξ_B) for 2×2 Hermitians (v2 L358-377; verifiable numerically; not an outside source).
3. μ–τ reflection => θ23=45°, |δ|=90° (Harrison–Scott / Grimus–Lavoura; body unread, to be verified).
4. m2=m3 exactly: θ23 unobservable (standard; no citation read).
5. RG: Antusch et al. 2003 and Casas et al. 2000 (formulas not captured; read bodies before any numeric comparison).
6. Σm_ν as a function of m_lightest from §3.3 formulas with the §3.1 Δm² values, compared to §3.2 bounds (NuFIT 6.0 Sec. 5 unread).

## 6. Open points that only Chazz can settle
1. v1 extraction loses all display equations; the exact θ₂₃ value v1 claimed (P246, P248-249, P263) and the Thm 3 minimizer formula can only be read from the docx itself.
2. v1 Thm 3 claims uniqueness "mod π" (P286); v2 finds two minima {0, π/2} (L347-350). Is v1's period a typo or a different functional?
3. IO: the near-degenerate pair is (m1,m2), not (m2,m3); the no-go's E = μ–τ plane assumes NO-like labelling. Which ordering is it meant for?
4. v2 asserts the near-degenerate regime is "relevant to Nature" (L106-116) with no Δm²₃₁ / Σm_ν numbers. Do the numbers (to be evaluated by the data-analyst) go in the paper?
5. Which NuFIT variant is the paper's reference (θ23 octant and δCP differ by variant, §3.1)?
6. `vocabulary.md:166` stale (flag only; lab crew fixes canon).
