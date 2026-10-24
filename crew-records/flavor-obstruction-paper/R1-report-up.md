# Reports for the manager (galt crew up fails: "no session 9439000a-… you may message")

## 2026-10-23 — R1 review closed
Flavor wt crew/flavor-paper: 15962c1, 6c96233, 39448b0. See AUDIT-theta23-paper.md, sections "Review" and "What the paper can still claim".
- Tr((HeHnu)^2) is not exactly one-cosine; the effect is small but real.
- The 3x3 commutator's e-row term is first order in delta_nu, so as delta_nu -> 0 the e-row orientation does the selecting, not psi0.
- A U(2) rotation can set theta13 to zero for any delta.
- Your relay, three points:
  - T0's 0.327/0.341 vs our 0.317/0.333 is a convention difference with identical inputs: T0 measures eps against the mean mass, we measure against the heavier mass.
  - The paper's invariant list is complete only at 2x2; tex 181-182 and 435 overreach.
  - In the complex case the minimiser's location is general, but theta23 there is arccos sqrt((1+-r)/2).
- The scope paragraph is written. KATRIN's m_beta < 0.45 eV does not exclude m2 ~ m3.
- Waiting for T1.

## Questions for Chazz (theorist, files/R1-theorist.md §6) — they decide the no-go's scope
Q1 Does the paper's functional mean the 3x3 commutator norm, including the e-row of H_e, or only the 2x2 compression on E?
Q2 H_nu|_E = lambda*1 has symmetry U(2), but the paper rotates only by real SO(2). Is that restriction intended?
Q3 Should E be the tilted degenerate eigenspace of the real data, rather than the mu-tau flavor plane?
Q4 Remark (tex 415-416), "stationary set frozen at {0,pi/4,pi/2}": is it meant only for the commutator-norm family?

## 2026-10-23 — NuFIT table + constructive skeleton
(1) The NuFIT 6.0 values used: analysis/data/nufit6_central.csv in the Flavor wt, committed in 4fa8f29. It holds the Table 1 central values only, with no uncertainties, for all four variant×ordering cases (IC19 w/o SK-atm, IC24 with SK-atm; NO, IO): θ12, θ13, θ23, δ, Δm²21, Δm²3ℓ. Its header records the source: arXiv:2410.05380 v2 (PDF dated 28 Oct 2025; v1 7 Oct 2024; JHEP 12 (2024) 216). The values are second-hand: the literature member transcribed them from the PDF, and nobody here has re-read them against it.
(2) The skeleton is constructive_section.tex in the Flavor wt, commit in the log after 39448b0. It is a new \input fragment and leaves the proofs untouched. It contains: §object (H_nu = H_nu^(0) + eps K, E = the nu2–nu3 pair; the object to fix is the eigenframe of K|_E, the obstruction its eps→0 limit); §conditions (a)–(e); one stub per candidate D1–D6, each with its gate (D6's gate is still to be named with the lab); a summary table. All results and verdicts read [pending]. Waiting for T1.

## 2026-10-23 — D1 drafted
The D1 subsection of constructive_section.tex is drafted from the T1 record (Flavor wt d04bcbb). It cites the record by path and commits 397b4e8/b3572ca/49fa1c6, and holds only what that record states.
Verdicts: (a) pass; (b) partial (the direction e_tau comes from data, the sign C from the model); (c) fail (it orients the solar pair 18.77° vs the atmospheric 3.57°); (d) unscored, per the challenge's finding that the null was unmatched; (e) fails outside the tuned corner. The theta13 point is labelled a consistency fact about the seed, not a D1 result.
One wording difference for you to check: the record itself scores (d) as "passes on sign". The draft says unscored, following your relay of the challenge. The challenge file is not saved yet, and its citation is a placeholder until it is.
Item #6 (data-analyst) is running. It adds an independent NumPy D1 check to the verifier: the exact flavor-diagonal transport from a mu-tau seed, the octant sign rule, and solar vs atmospheric rotation. Ready for T2.

## 2026-10-23 — D1 verifier check in
Flavor wt 7a6d549: verify_theta23_obstruction.py now has a D1 section (analysis/d1_rg_scale.py). It exits 0, ALL CHECKS PASS, and the capture is analysis/verify-run-2026-10-23-D1.txt. The old checks are unchanged.
- The record's octant rule is reproduced in all 8 cases (2 NuFIT variants × 2 orderings × SM/MSSM).
- New nuance for the lab's (c) leg: the solar pair rotates more than nu2–nu3 for IO at every mass, but for NO only above m_lightest ≈ 4.4 meV. Below that (hierarchical NO) the atmospheric pair rotates more, but the running is smallest there.
- D1 text updated to say this. The verdict is unchanged. Ready for T2.

## 2026-10-23 — D2 drafted (PROVISIONAL)
The D2 subsection of constructive_section.tex is drafted from the T2 calculator record (Flavor wt e3bcf4c). It cites the record by path and commits a139b08/b6d9ed5/71a1021, and is marked provisional in a comment and in the verdict heading.
Spot-checked against the T2 capture: 0.8414; the cos(psi/2) form (psi=150° gives 0.2178); theta13 running 0.067–34.6°; P = 0.497–0.509. The delta distances map as: IO 15°/4° = IC19 w/o SK-atm 285° / IC24 with SK-atm 274°; NO 87°/58° = 177°/212°.
Verdicts: (a) pass; (b) fail (holds only under the imposed reflection); (c) fail; (d) fail (octant at chance); (e) delta = ±90° is scorable but not supported.
No verifier check for D2 yet. I'll add one, built on the Grimus–Lavoura reflection gate, once you confirm T2. Waiting for your confirmation and T3.

## 2026-10-23 — D1 citation closed
Flavor wt 3a6d7f7: D1 now cites lab 6834d33 and the challenge at galt-crews/flavor-obstruction-lab/files/R1-T1-challenge.md. (d) reads unscored, null "unmatched, theta23-only, phases (0,0)"; I dropped the 10.3% closeness figure from the draft because it rests on the same unmatched null (your relay did not say it was struck; restore it if the record keeps it). D2 stays provisional until #17 closes.

## 2026-10-23 — D1 10.3% restored
Flavor wt 5068385: restored as a closeness figure against an unmatched null, with your exact qualification, and not as a (d) score.

## 2026-10-23 — D3 structural paragraph drafted (PROVISIONAL)
Flavor wt 89decab: Remark "flavor labels are functions of H_e" (Lagrange interpolation), plus the two consequences: in the exact mu-tau plane every f gives the flavor frame with delta_e -> f(m_tau^2)-f(m_mu^2), so the no-go covers it; with E tilted, there is a 1-parameter curve of frames and the label Q is a 1-parameter fit to 1 target. (b) is marked failed structurally (canon assigns no lepton harmonics). (c)-(e) and the numbers stay [pending] until the T3 record lands. Theorist item #7 is checking the claims independently, including the edge cases f(m_mu^2)=f(m_tau^2) and the order flip, and the definition of r. D2 is still provisional.

## 2026-10-23 — #7 closed; D3 remark corrected (PROVISIONAL)
files/R1-theorist-D3.md holds theorist-2's full note verbatim, with the author named. Flavor wt 47b8dc7 folds in your three corrections:
(1) The remark is E-restricted. The 3x3 Tr is exact; the commutator norm is not frozen (e-row second harmonic (1/2)(f_tau-f_mu)(2f_e-f_mu-f_tau)(|w_mu|^2-|w_tau|^2), phase set by H_nu's e-row).
(2) r* closed form; r*=1 at theta13=0, so the only aligning label is the inert one and all the leverage is theta13!=0.
(3) delta!=0: exact alignment is impossible, min_r|A+rB| = |A||sin(argA-argB)| > 0.
Also added: the order flip (theta23 -> 90-theta23 as a label) and d_f=0 (inert, (a) fails).
theorist-2's open gaps Q1-Q5 (in the file) go to Chazz via you: Q1 labels on flavor vs mass states; Q2 the scoring functional for a complex frame; Q3 the label form; Q4 n_e free or fixed; Q5 3x3 vs E-restricted norm. r* and the e-row coefficient are derived on paper but numerically unvalidated; its calculator checks (i)-(v) are listed in the file. Numbers stay [pending] for T3. D2 is still provisional.

## 2026-10-23 — theorist-2 Q1–Q5 resolved
Q1–Q4 are recorded in the margin of files/R1-theorist-D3.md as "manager ruling, R1-plan T3 spec (lab 285a9a8)". Q5 is the no-go scope Q1.
The list for Chazz now holds only:
- no-go scope Q1: 3×3 vs 2×2 compression;
- no-go scope Q2: U(2) vs SO(2);
- no-go scope Q3: E = tilted degenerate eigenspace vs the μ–τ plane;
- no-go scope Q4: scope of the "frozen" Remark;
- which harmonic labels, if any, e, μ, τ carry.
D2 is still provisional, pending the saved challenge.

## 2026-10-23 — D2 pre-adapted to the T2 challenge (PROVISIONAL)
Flavor wt f3d7e1d applies the six edits:
- The u(2)->o(2) result is cited as a known theorem (BRSM 1999, citation under check) and limited to the cosmology-excluded degenerate limit.
- (d) is "not applicable" (P=0.5 forced by symmetry).
- The breaking is "to reach the central values".
- (e) is "delta=+-90 with free sign, not scored", with a [pending] slot for the lab's 3-sigma scoring.
- The tau-probe leg is split SM/MSSM with eta = J_tau, values [pending the lab's record edit].
Lookups: literature item #8 is running (NuFIT 6.0 3-sigma theta23 and delta from the v2 PDF; the BRSM and Zhou citations). Results will follow here.

## 2026-10-23 — D2 FINAL; NuFIT ranges; citations
Flavor wt b44967f.
(1) NuFIT 6.0 Table 1 p.6, as printed, from the arXiv:2410.05380v2 PDF (literature-2, opened; files/R1-literature-T2.md; committed as analysis/data/nufit6_ranges.csv in 728a1d0). delta_CP bfp +1s/-1s | 3s:
- IC19 w/o SK-atm NO: 177 +19/-20 | 96 -> 422 (62)
- IC19 w/o SK-atm IO: 285 +25/-28 | 201 -> 348
- IC24 with SK-atm NO: 212 +26/-41 | 124 -> 364 (4)
- IC24 with SK-atm IO: 274 +22/-25 | 201 -> 335
theta23 3s ranges: 41.0-50.5, 41.4-50.6, 41.3-49.9, 41.5-49.8. All central values in nufit6_central.csv match the PDF.
(2) Citations:
- BRSM: confirmed (PRL 82 (1999) 683, "Degenerate and Quasi Degenerate Majorana Neutrinos"). The abstract says nothing on CP parity. The body (pp.2-4) treats THREE degenerate neutrinos: U can be rotated away only at the same CP parity, and not when one neutrino has opposite parity. Cited as such; the pair case is the lab's 2x2 sub-case.
- Zhou 1409.8600: corrected. It is a research paper (Y.-L. Zhou, "mu-tau reflection symmetry and radiative corrections"), now cited as D1 prior art. The review is Xing-Zhao 1512.04207, which was already in the paper's bibliography.
(3) D2 is final in constructive_section.tex: it cites lab 2010400 and the challenge path, and (e) is your sentence. The lab's (e)-edit commit is a placeholder until you relay it. The verifier's new D2 section (analysis/d2_reflection.py, Grimus-Lavoura gate) passes; captures analysis/d2_reflection-run-2026-10-23.txt and verify-run-2026-10-23-D2.txt. D3 is still provisional, waiting for T3.

## 2026-10-23 — D2 final cite; D3 numbers drafted (PROVISIONAL)
Flavor wt 7b982be.
- D2 now cites f1f00e3 alongside 2010400 and the challenge path. The (e) sentence already points at analysis/data/nufit6_ranges.csv (728a1d0). The D2 verifier check is in (b44967f).
- D3 numbers are from lab 3f67f4b, each traced to a PART of run-2026-10-23-t3.txt: r*_23 0.816-0.820 (delta 0) / 1.218-1.222 (180); arg(A/B) at the data delta 0.59 to -11.4 deg; q_min 9-11 (B) / 3-5 (A), null P 0.55-1.00 / 0.44-0.95; window share 4.8-6.5% (linear, A), dropping to 0.8-1.0% when theta13' is also held within 0.5 deg; reading-B octant follows sign(cos delta) (>=99.9% vs <=0.1%). The solar-pair sentence is [pending].
- DISCREPANCY for the lab: the verdict says "6 of 8 cells >= 5%", but the PART 3 table (all labels, linear, reading A) shows 7 of 8: 6.47, 5.22, 5.02, 5.63, 6.26, 5.06, 5.57, and only IO IC24 delta=180 below at 4.84. The pooled table also gives 7 of 8. The draft marks this [under check].

## 2026-10-23 — D3 pre-adapted to the T3 challenge (PROVISIONAL)
Flavor wt c301949 + 6abf12e. All five relabellings are applied:
(1) Prop. rstar: r* = -U_mu2U_mu3/(U_tau2U_tau3) = 1 + U_e2U_e3/(U_tau2U_tau3), proved in one line from column orthogonality, leading order 1 -/+ 2 tan12 sin13/sin2th23 (~0.80/1.20). The verifier's new D3 section (analysis/d3_rstar_identity.py) reproduces the lab's PART 2 r* in all 8 cells to 10 digits, and r*=1 at theta13=0 (= U_e3=0). The reading-B octant is stated as the subleading U_e3 term.
(2) delta 0/180 near-mirrors (r -> 1/r), ~4 independent cells; IO scored at delta in {0,180}, far from its data delta.
(3) CP-odd channel det[Phi,H_nu] = 2iJ Delta(n) Delta(m^2): labels only rescale J. The form and sign are [under check] (literature-3, item #9) alongside the lab's numeric check.
(4) The window share is restated as measure-dependent and at the null rate. No PART 7 ratio is quoted anywhere. The 6-vs-7 count bracket is gone with the old wording; it is moot under (4).
(5) Residual-T context, citing King-Luhn (already in the bib) [claim under check].
literature-3 (item #9) is running the four citation checks: King-Luhn, Froggatt-Nielsen, Jarlskog 1985 exact form, Feruglio-Romanino. The solar-pair sentence is still [pending].

## 2026-10-23 — D3 FINAL; citations; what the battery has shown so far
Flavor wt 602d419.
D3 is final. It cites lab 3f67f4b, ea229c1, 6ca51a3, a80e32b and ef716dd, plus the challenge path. It now contains:
- Prop. rstar (the unitarity identity, verifier D3 check);
- Prop. jarlskog (the closure, with the lab's 1e-11 verification);
- the convention-free joint-alignment sentence (no single r aligns both pairs; >=31.4 deg is left on one);
- the window share as canonical 4.8–6.5% vs signed 2.6–3.3%.
Citations (literature-3, files/R1-literature-T3.md):
- **King–Luhn:** title/journal confirmed, but the residual-T claim was NOT found in the part of the text reachable (ar5iv truncated at sec 4.2). The claim is stated verbatim in arXiv:1407.4722, "Lepton mixing, residual symmetries, and trigonometric Diophantine equations" (opened), which is now cited for it. Its authors are still to be confirmed; there is a TODO in the bib. King–Luhn and Feruglio–Romanino are cited as reviews. Feruglio–Romanino is confirmed ("Lepton Flavour Symmetries", RMP 93 (2021) 015007, 1912.06028), and no newer dedicated review was confirmed.
- **Froggatt–Nielsen:** search-only, not opened; cited as standard.
- **Jarlskog 1985:** the exact printed form was NOT seen (the PDF was unreadable). Secondary sources disagree on the overall sign, which depends on the commutator order, the difference ordering and the J convention. The proposition states the form in the lab's verified conventions (PDG J, [Phi, H_nu], the differences as written) and calls the sign convention-dependent. If the lab or anyone with text-PDF access can read PRL 55 1039 or the PDG CKM review, that would close it.

**What the battery has shown so far (D1–D3).**
- **All three fail on (b).** The orientation each supplies comes from an input the theory does not fix:
  - D1, scale: the octant follows the SM/MSSM coefficient C.
  - D2, Majorana phases: 45°/±90° come from an imposed reflection, and the breaking it needs is free.
  - D3, harmonic labels: no canon assigns lepton harmonics.
- **θ13 is the only lever for flavor-diagonal operators.** Any operator diagonal in the charged-lepton basis is a function of H_e, so it gives the flavor frame on the exact μ–τ plane and is covered by the no-go. Off that plane it acts only through U_e3: r* = 1 + U_e2U_e3/(U_τ2U_τ3), and r* = 1 at θ13 = 0.
- **The CP-odd channel is closed.** det[Φ, H_ν] is J times Vandermonde factors, so labels only rescale J.
- **The data frame needs an external input the theory does not yet supply.** Every mechanism tried can reach the measured frame only by tuning a free input, and none beats its null. In D1 the running orients the solar pair more than ν2–ν3; in D3 no single label ratio orients both pairs.
- **Still to test:** T3b (an external harmonic frame), then D4 (a 3+1 fourth direction), D5 and D6.

## 2026-10-23 — two closes
Flavor wt 17e114d. (1) The arXiv:1407.4722 author is B. Hu, "Lepton Mixing, Residual Symmetries, and Trigonometric Diophantine Equations", Phys. Rev. D 90 (2014) 073012, doi:10.1103/PhysRevD.90.073012 (lead OPENED https://arxiv.org/abs/1407.4722); the bib TODO is closed. (2) The Jarlskog sign now reads exactly as you worded it. The sterile 3+1 bounds (item #10, literature-4) are running.

## 2026-10-23 — 3+1 sterile bounds for T4 (D4)
Flavor wt b9e5ea0: analysis/data/sterile_3p1_bounds.csv. The verbatim note is files/R1-literature-T4.md (literature-4, item #10).
**CAVEAT:** "OPENED" means the page was fetched through a summarising tool. The values were NOT read off the primary tables. Re-check any value before the paper quotes it as a bound; as inputs to a lab scan they are usable with that flag.
Bounds (element: bound, CL, Dm41^2, source):
- **Gariazzo et al. 2017** (1703.00860, Table 6 PrGlo17; best fit Dm41^2 ~1.7 eV^2):
  - |Ue4|^2 0.011–0.032 (3σ allowed)
  - |Umu4|^2 0.0016–0.0051 (3σ allowed)
  - |Utau4|^2 < ~0.071 (3σ)
  - SUSPECT extra line: |Utau4|^2 < ~0.011 (90%), re-read needed
- **Dentler et al. 2018** (1803.10661):
  - |Ue4|^2 best fit ~0.009–0.010 at ~1.3 eV^2 (nu_e disappearance)
  - |Ue4|^2 < ~0.1 (99%, context to check)
  - |Umu4|^2 < ~1e-2 (qualitative, from a figure)
  - |Utau4|^2 < 0.13 (90%) / 0.17 (99%), independent of Dm41^2
- **|Utau4|^2 direct limits:**
  - IceCube DeepCore 7.5 yr (2407.01314): < 0.0574 (90%, Dm41^2 >= 1)
  - IceCube 3-param 7.6 yr (2406.00905): < 0.035 (90%, fitted 2.4–9.6 eV^2)
  - Super-K (1410.2008): < 0.18 (90%, > 0.8 eV^2, assumes Ue4 = 0)
  - MINOS/MINOS+ (1710.06488): no |Utau4| number found; most stringent sin^2 θ24
- **Tension:** Dentler 2018 excludes an LSND-type 3+1 explanation at 4.7σ (appearance vs disappearance, p < 2.6e-6). Nothing newer confirmed.
- **MicroBooNE** (2210.10216, PRL 130 011801): joint 3+1 search, no evidence, 95% exclusions.
- **KATRIN** (2201.11593): no sterile signal, 95% exclusion. Its limit near 1 eV^2 was not extracted, and a newer KATRIN result exists that was not searched.
Gaps, if the lab needs them: a newer global fit (the 2022 white paper 2203.07323 was not opened), the newer KATRIN limit, and the MINOS+ |Utau4| number.

## 2026-10-23 — process fix, hash map, D4 structural draft
(1) Understood: from now on I fix forward only. No amend or rebase once a hash is in a report-up or cited.
Hash map (the only rewrite so far): **3b55b00 → 602d419** (D3 final; same tree, commit-message fix only, "%%" → "%"). Every other hash quoted in this file is unchanged on crew/flavor-paper.
Current branch, newest first:
- 434182f D4 structural (provisional)
- b9e5ea0 sterile bounds
- 17e114d bib Hu + Jarlskog wording
- 602d419 D3 final
- 6abf12e D3 verifier
- c301949, 7b982be D3 drafts
- b44967f D2 final
- 728a1d0 NuFIT ranges
- 89a9326, f3d7e1d D2 drafts
- 47b8dc7, 89decab D3 remark
- 5068385, 3a6d7f7, d04bcbb D1
- 7a6d549 D1 verifier
- e3bcf4c D2 draft
- 9769726 skeleton
- 39448b0, 6c96233, 15962c1, 4fa8f29, 5223dc3 audit
- faf2d83 drafts in
(2) Sterile bounds: D4's text says explicitly that the csv is a sampling region only, and quotes none of its values as checked.
(3) D4 structural paragraph drafted (434182f, provisional): (a) orthogonal directions are inert; (b) the frame-label freedom persists; (c) H_e^(4) adds only the free x. It ends with the consequence: a non-orthogonal frame oriented from outside the data is needed (the T3b setting). The explicit psi_0(x) is [pending]. theorist-3 (item #11) is checking (a)–(d), including whether non-unitarity of the 3×3 block changes the picture. The note will be saved as files/R1-theorist-D4.md.

## 2026-10-23 — D4 structure verified (theorist-3), folded in (PROVISIONAL)
Flavor wt 459ed7a. Verbatim note: files/R1-theorist-D4.md. Derived on paper; nothing run.
- **(a) agree:** exact for any n and any eps. P_E Q P_E = q0 P_E. Needs E to be the full lam-eigenspace; Dm41^2 ~ 1 eV^2 ensures that.
- **(b) agree, one correction:** at eps=0 the individual |N_mu3|, |N_tau3| are not measurable (oscillations see only P_1, P_4, P_E). What changes is their label-level value. Invariants: P_E, P_1, P_4.
- **(c) agree, sharpened:** explicit psi0(x) = ½ atan2(2Re h0_23 + 2x Re(s̄2 s3), D0 + x(|s2|²−|s3|²)). It is x-independent iff (s2,s3) is an eigenvector of h0 (in particular s=0). Otherwise x sweeps a one-parameter curve of frames, like r in D3.
- **Scope (OVERSTATED as you worded it):** "independent of flavor dimension" holds for the compression class only. Full-space functionals carry a cos4 term with phase set by P_E H_e P_perp. D4 now says so.
- **(d) agree, narrowed:** orthogonal dimensions and outside-E mass operators cannot help, and a sterile with E-content adds only the free x, so (b) fails. Unitarity of N plays no role. NOT excluded: a Q non-scalar on E and oriented outside (H_e, H_nu), e.g. seesaw non-unitarity eta ~ Y_nu M^-2 Y_nu† (not a function of (H_e, H_nu) in general), or T3b's external frame.
Gaps theorist-3 raises:
- **For Chazz (canon-level):**
  - Q1 what is x: does the theory define a sterile H_e entry?
  - Q3 is a non-unitarity operator eta in the relational class, or only functions of (H_e, H_nu)?
  - Q2 (SO(2) vs U(2)) is already no-go scope Q2.
- **Not done:** Q4 the n=4 full commutator harmonics; Q5 the sizes of N_s2, N_s3 and x.
- **Calculator checks C1–C5** are listed in the file, should the lab want a D4 verifier gate: fit psi0(x); P_E Q P_E scalar for n=3..6; column/P_E invariance; FFT of the full norm for cos4; a random eta as the non-excluded case.
D4 numbers are [pending] the T4 record.

## 2026-10-23 — R2 opened (paper side)
- R2 is on this board, linked to your R2. Its plan closes after audit items #12 and #13.
- Item #12 (literature-5): sum rules, PDG CKM, QLC, anarchy, A4 modular. The result goes to files/R2-literature.md and analysis/data/r2_inputs.csv; I'll report it here.
- Item #13 (theorist-4): the octant paragraph (sigma distances from nufit6_ranges.csv) and the D11 proposition with proof (first-order zero, second-order form, other invariants).
- Skeleton: r2_section.tex (Flavor wt 848219f) with the octant slot and the D7–D11 stubs (gate, question, null), all [pending].

## 2026-10-23 — R2 items #12 and #13 in; #14 queued
**#12 inputs** (Flavor wt 46db69e): analysis/data/r2_inputs.csv; verbatim note in files/R2-literature.md.
- PDG 2024 CKM, read from the PDG page images:
  - sin θ12 0.22501±0.00068
  - sin θ13 0.003732 +0.000090/−0.000085
  - sin θ23 0.04183 +0.00079/−0.00069
  - δ 1.147±0.026 rad (65.7°)
  - Wolfenstein λ 0.22501, A 0.826, ρ̄ 0.1591, η̄ 0.3523
  - J 3.12e-5
- TM1 (King 1701.04413 Eq50–51 = Albright–Rodejohann Eq23–24, agreeing as transcribed): sin²θ12 = (1/3)(1−3s13²)/(1−s13²); cos δ = −cot2θ23 (1−5s13²)/(2√2 s13 √(1−3s13²)).
- TM2: sin²θ12 = (1/3)/(1−s13²). Its cos δ relation (A-R Eq16) came back GARBLED; re-read before gating.
- QLC: θ12+θ_C = π/4 is Minakata–Smirnov. The θ23 relation is in Raidal (Eq1), not M-S. Neither gives an explicit U_bimax V_CKM† formula as transcribed.
- Anarchy: flat in sin²θ12, sin²θ23, cos⁴θ13, δ (dGM 2003: KS 64%; dGM 2012: P_KS 44% central, >24% within 3σ).
- Modular: the Feruglio 1706.08749 best fits CONFLICT between two fetches, and both have s13² ≈ 0.045–0.049, far above data. NPPT A4 (1808.03012) gives ranges only, not tau=i or ω. NPPT S4 (1811.04933) case A best fit: τ = ±0.1045+1.01i, s23² = 0.551.
- GAP: no source retrieved for the generic first-order "θ23−45 ~ θ13 cos δ" rule.
- All non-PDG values came through a summarising fetch: usable as lab gates only after a re-read of the source.

**#13 theorist** (Flavor wt 91556dc in r2_section.tex; verbatim note in files/R2-theorist.md):
- **Octant:** sigma distances confirmed: 3.9 (IC19 NO), 1.7 (IC24 NO), 4.0 (IC19 IO), 3.2 (IC24 IO). BUT every published 3σ window contains 45, so 3.9/4.0σ is a Gaussian extrapolation the table contradicts (non-Gaussian χ²). The paper does not print them as an exclusion. Established: near-maximal (best fits within 3.6°, 3σ within [41.0, 50.6], so "≲5°" is slightly tight at 3σ).
- **Prop. Jextremum (exact):** J's critical points are exactly (45, ±90). At every fixed δ the θ23-extremum is exactly 45, and J(90−θ23) = J(θ23), so maximal CP is OCTANT-BLIND at all orders, with no δ–θ23 cross term. Tr([He,Hnu]³) ∝ J inherits this; the e-row Majorana invariants are θ23-independent.
- **NEW CANDIDATE for the lab:** maximal misalignment, max ‖[He,Hnu]‖_F. At leading order it gives |U_τ3|² = 1/2 → θ23* ≈ 44.4° (first octant, δ-independent at that order, independent of the mass scale). Theorist-4's checks C1 (true Δχ²(45) from the NuFIT χ² tables), C2 (exact argmax, gated) and C3 (Haar null) are in the file.
- Questions for Chazz: Q1 is maximal J/misalignment meant as a principle, and over which variables? Q2 which Majorana invariants count as natural? Q3 can we get the NuFIT χ² tables?

**#14** (field survey) is assigned to literature-5, which picks it up after #12.

## 2026-10-23 — your #12/#13 note crossed my report
All of it was already done before your note arrived: files/R2-literature.md and files/R2-theorist.md saved verbatim; r2_inputs.csv committed (46db69e); the octant paragraph uses theorist-4's LaTeX with no exclusion claim; Prop. Jextremum is in r2_section.tex (91556dc), now marked final-on-paper (cce4274). The follow-ups (a) TM2 exact, (b) Feruglio tables, (c) NuFIT 6.0 Δχ²(θ23=45) for all four cases are item #15, given to a fresh literature hire running in parallel with #14.

## 2026-10-23 — R2 #14 field survey and #15 follow-ups in
Notes saved verbatim: files/R2-field-survey.md (literature-5) and files/R2-literature-followups.md (literature-6). Inputs: Flavor wt 116640b, analysis/data/r2_inputs.csv.
**#15 follow-ups (read from PDF page images, OPENED-img):**
- (a) TM2, Albright–Rodejohann Eq16 RHS: cos δ tan 2θ23 = (1−2|Ue3|²)/(|Ue3|√2 √(2−3|Ue3|²)) ≈ (1/√2)(1/|Ue3|)(1 − (5/4)|Ue3|²). Its middle fraction is unverified. Eq18: (½ − sin²θ23)² ≈ (cos²δ/2)|Ue3|² ≈ (3cos²δ/2)(sin²θ12 − ⅓).
- (b) The Feruglio fits are VERIFIED and the conflict is resolved. Ex.1 (IO) gives s23² = 0.651 with rows 2↔3 permuted (0.349 unpermuted). Ex.2 (NO) gives s23² = 0.331. Ex.3 gives s23² = 0.749 ("clearly excluded"). All three have s13² far off data, so they are poor fits, not gates for D10's angles; they do gate the Y_i(τ) machinery at the printed τ.
- (c) **NuFIT Δχ²(45): NOT obtained.** nu-fit.org is unreachable from this machine (my own curl timed out too). The paper text says the octant is degenerate at Δχ² < 4 for both orderings and both variants. A new item (#16) reads Fig. 1 off the PDF image, marked DIGITISED, as the best available.
**#14 field survey, the essentials:**
- **Corrections to the brief:**
  - The "15 expressions with one charged-lepton rotation predict θ23 within 0.2°" sentence is in Denton–Gehrlein–Truelson 2606.05291 Sec 6.2. "Within 0.2° of what" is unverified (the tool said "of maximal").
  - arXiv:2603.06934 is Barger's lattice-flavon model, not an HK/DUNE paper; its two branches are (43°, 286°) and (46°, 299°).
  - arXiv:2409.11824 is Cao et al. on HK+DUNE: the band 0.48 ≤ s23² ≤ 0.54 can't separate maximal from octant.
  - The sin²θ23 = 0.473 (lower octant at 1.1σ) fit is Capozzi et al. 2503.07752. Denton's 2026 table has the upper octant (0.561), so the two recent compilations disagree on the octant.
  - The 2607.07865 range 0.435–0.583 is NOT confirmed.
- **Sharpest predictions:**
  - Littlest Seesaw CSD(n), Costa–King 2307.13895: 3 parameters (m_a, r, η) predict the 3 angles, δ and the masses. θ23 = 45.5 ± 2.4° (n=3, δ ≈ 272° or 88°); 41.4° (n=1+√6); 51.5° (n=2.5, disfavoured).
  - New LS (1906.11414): s23² = 0.599 (upper).
  - All are TM1, so cos δ ties to θ23 through the TM1 sum rule.
- **Reframings:**
  - (i) Single charged-lepton 2-3 rotation (2603.29876 Case III): θ23 ≈ θ23^ν − θ23^l, splitting θ23 into a neutrino part and a charged-lepton part.
  - (ii) τ = i residual (anti)symmetry (Kashav–Patel 2410.00565): one fixed PMNS column plus a massless/degenerate pair, with zero continuous input. This is the one place a symmetry itself produces R1's degenerate ν pair.
**Literature-5's recommendations, ranked by sharpness** (all need the source tables gated first):
1. Littlest Seesaw CSD(n) and New LS (D12): 2–3 inputs, 5+ predictions.
2. Exact μ–τ reflection, the zero-parameter baseline (done as D2).
3. Exact TM1 (the δ–θ23 correlation; D7).
4. τ = i residual (anti)symmetry (Kashav–Patel): zero continuous input; closest to R1's degenerate pair. Candidate for D10's fixed-point leg, or its own D13.
5. Single charged-lepton rotation, U = U_e† U_ν: one angle; Denton's 15 expressions.
6. A5 + gCP (Turner 1507.06224): θ23 in case-dependent ranges (38.2–45° or 45–53.3°).
7. Modular A4/S4 at generic τ: s23² > 0.54 (upper).
8. Texture zeros: m_μμ = 0 lower octant; m_μτ = 0 and m_ττ = 0 upper.
9. Anarchy and QLC as nulls (D9, D8). QLC with PDG gives θ23 ≈ 45 − 2.40 = 42.6° (lower).
**Item #16 (new hire) is running:** verify the LS tables, the New LS Table 1, the 2607.07865 range and the Denton sentence; digitise Δχ²(45) from NuFIT Fig. 1.

## 2026-10-23 — R2 #16 page-image verification in (all OPENED-img)
Verbatim note: files/R2-literature-verify.md. Flavor wt: 599b24e adds the verified rows to r2_inputs.csv; 74d44c3 adds the octant paragraph.
- **NuFIT 6.0 Δχ²(sin²θ23 = 0.5), DIGITISED from Fig. 1 (±0.3):**

  | case | relative to own-ordering min | relative to global min |
  |---|---|---|
  | IC19 NO | 3.8 | 4.4 |
  | IC19 IO | 5.7 | 5.7 |
  | IC24 NO | 2.2 | 2.2 |
  | IC24 IO | 2.8 | 8.9 |

  So maximal mixing is disfavoured at about 1.5–2.4σ, not 3.9σ. The octant paragraph now says this. The nu-fit.org tables are still unreachable.
- **Littlest Seesaw CSD(n), Costa–King Tables 3–5:**
  - The tables carry only 3σ ranges, with no masses and no χ². Each η has a normal and a flipped branch (θ23 → 90−θ23, δ → δ+180, Eq5.10).
  - n=3: only η=2.11 normal (θ23 45.5, δ 272.2) and η=4.17 flipped (θ23 44.5, δ 267.9) lie within 3σ, so CSD(3) predicts θ23 ≈ 45 ± 1 with δ ≈ 270.
  - n=1+√6: 41.4/48.6.
  - n=2.5: 51.5/38.5.
  - Structure verified: alignments (0,1,1), (1,n,n−2), the mass matrix Eq5.13, the TM1 predictions Eq5.26/5.31, and the solar sum rule Eq3.3.
- **New LS (1906.11414):**
  - Table 1 highlighted row (x=−1/2, η=±π/2): s23² 0.599, δ/π ±0.354, χ²min 2.487.
  - The paper's free fit (Eq9–10) gives s23² 0.582, χ²min 0.384.
- **Δ(96) modular LS (2607.07865):** per-pattern ranges in Table 8; the union is [0.435, 0.584], not a single printed range.
- **Denton "15 expressions":**
  - Verbatim: "All 15 possible expressions containing only one rotation in the charged lepton sector predict θ23 within 0.2°. Moreover, 13 of these predict θ23 to be less than 0.5° from maximal, in regions where DUNE and HK would be unable to provide a definitive determination of the octant."
  - The 0.2° is the WIDTH of each prediction (Table 7 caption: s23² ± 0.001), not the distance from 45.
  - The 15 values cluster at s23² ≈ 0.488, 0.511, 0.498–0.501, plus 0.462 and 0.454.
  - Whether "15 possible" equals 8+7 valid is UNRESOLVED.
- **Reading for the lab:** the single-rotation class and CSD(3) are the sharp near-maximal predictors (±0.5–1°). Exactly that band is the one HK+DUNE cannot resolve (0.48 ≤ s23² ≤ 0.54, Cao et al.).

## 2026-10-23 — D3b + D4 numbers drafted (PROVISIONAL); octant paragraph; R2 fetches item
- **1d4468a, D3b (new subsection before D4):** drafted from lab e34fc31 (gate 6e6ff74, amendment e92d5af).
  - Exact result: θ_joint = min_α arccos|U_α2| = 47.2–52.2°, label-independent.
  - The F<0.05 yardstick does not discriminate; this is flagged as a design flaw.
  - Divisor-angle look-elsewhere ≥0.133; one 1σ hit among 52 tests.
  - Verdict: (b) fails by construction; (d) unscored.
- **1d4468a, D4 numbers:** from lab dc6d4b7/a7cc8a3/e734dbc.
  - P_E Q4 P_E = 0 to 2.7e-17.
  - ψ0 moves by a median 4.1° over x = 0..10 m_τ².
  - PART B: Q_act median 16–23× the ε-term at m4² = 1 eV² (5–95%: 0.8–54×), yet no measured element changes.
  - Window rate 0.57–1.03× random, so no enrichment.
  - Verdict: (a) fails exactly, (b) fails.
  - Small wording point: the record's verdict line says "2-53x"; the 16–23× figure is the median and is the one I used, with the 5–95% band given.
- **ee7122b, octant paragraph:** now says the octant depends on dataset AND ordering (IC19: IO best, NO at Δχ² 0.6; IC24: NO best, IO at 6.1). It reports each candidate's octant as a prediction, unscored, until the octant Δχ² is in.
- **Item #17 (literature-7):**
  - (a) the NuFIT 6.0 2D (θ23, δ) χ² grids. nu-fit.org is unreachable from here, so the item tries mirrors and falls back to digitising the 2D contour figure.
  - (b) θ12 from Table 1, all four cases.
  - (c) the per-variant octant Δχ² from Fig. 1.
  - The grids, if reached, get committed under analysis/data/ and θ12 goes into nufit6_ranges.csv.

## 2026-10-23 — Valencia 2020 grids committed
Flavor wt b634f2f: analysis/data/valencia2020/ (1-dim, 2-dim, PROVENANCE.md, SHA256SUMS.txt; FILES.sha256 adds per-file hashes; 45 files). r2_section.tex now carries the rule that every 2D-grid number is labelled "Valencia 2020, not NuFIT 6.0", with 1D minima cross-checked against NuFIT Table 1 (f50ce72). #17 (literature-7) now also covers (a') the Zenodo attribution and data epoch, (b') any 2024/25 Valencia or Bari profiles, and (d') the Minakata–Smirnov Gamma_delta definition. The NuFIT grids stay "not obtained".

## 2026-10-23 — #17 in; novelty audit (#18 literature-9, #19 theorist-5) started
Flavor wt 2a25254. Verbatim notes are in files/R2-literature-T17.md.
- **(a) NuFIT 2D grids NOT obtained.** Fig. 2's (sin²θ23, δ) panel is marginalised over the ordering (IC19 vs IC24 only), so per-ordering contours can't be digitised. The NuFIT site is also at globalfit.astroparticles.es, but that fetch failed (expired certificate). D7 therefore runs on the Valencia 2020 grids, labelled as such.
- **(b) θ12 from NuFIT 6.0 Table 1**, added as columns to nufit6_ranges.csv:
  - sin²θ12 0.307 +0.012/−0.011 (IC19 NO), 0.308 for the other three; 3σ 0.275–0.345.
  - θ12 33.68 +0.73/−0.70; 3σ 31.63–35.95.
  - Correction to a past transcription: IC19 NO sin²θ23 is printed 0.561 +0.012/−0.015.
- **(c) Octant Δχ², DIGITISED from Fig. 1 (±0.3):**
  - IC19 NO ~1.0
  - IC19 IO ~3.9
  - IC24 NO ~0.5 (low confidence, ±0.5)
  - IC24 IO ~1.8
  - The text bounds them all by < 4, consistent. The octant is at most ~2σ in any case, so octants stay unscored predictions.
- **(a′) Valencia attribution:** 8 authors (de Salas, Forero, Gariazzo, Martinez-Mirave, Mena, Ternes, Tortola, Valle); the staging PROVENANCE listed 6. JHEP 02 (2021) 071; data to ~2020; Zenodo file MD5 f84740840a22ff2bbe84f34da0c1f7d2.
- **(b′)** No 2024/25 Valencia or Bari downloadable profiles were found.
- **(d′) Minakata–Smirnov:**
  - Γ_δ = diag(1, 1, e^{iδ_l}) (below Eq. 9).
  - Eq. 27: U_MNS = V_CKM† Γ_δ R23(π/4) R12(π/4).
  - Eq. 36: D23 = ½ sin²θ_C + cos²θ_C |V_cb| cos δ, with D23 = ½ − sin²θ23 (Eq. 11). Eq. 37 gives 0.025 ± 0.039.
  - These are in r2_inputs.csv.
**Novelty audit (Chazz's question):**
- #18, literature-9: verbatim quotes, read from page images, of Jarlskog, BRSM, Procesi/Sibirskii, Lam, Feruglio–Hagedorn–Ziegler, Fonseca–Grimus, Hu, Holthausen–Lim–Lindner, King's semi-direct approach, Chamseddine–Connes–Marcolli and Fritzsch–Xing.
- #19, theorist-5: waits on #18. It writes the table (our quote | closest known | verdict) for rows 1–6e, the 'Relation to known results' LaTeX, and the final verdict line.
- Deliverable: files/R2-novelty-audit.md.

## 2026-10-23 — NOVELTY AUDIT done (Chazz: "is ours authentic or already found?")
Files:
- files/R2-novelty-sources.md: two independent page-image passes (literature-7 and literature-9), verbatim quotes with page/eq.
- files/R2-novelty-audit.md: theorist-5's 10-row table, the LaTeX paragraph, the verdict and the gaps.
- Flavor wt b7820bc: relation_known_results.tex (new fragment for the Discussion) and 6 new bibitems; the unverified journal refs are marked.

**VERDICT:** the obstruction is a REFORMULATION, not new mathematics. It is three known pieces, assembled in spectral-triple language:
- the free eigenframe of a degenerate eigenspace (BRSM p.5: "in the exact degeneracy limit, the neutrino mixing matrix U0 is only defined up to an arbitrary orthogonal transformation"; Jarlskog p.1040: "if mt were equal to mc the CP phase could be removed");
- the Sibirskii 2×2 pair generators (Drensky math/0506614, Thm 2.10);
- the residual-symmetry partial pattern (Lam 0708.3665 p.2: "one column of the mixing matrix is fixed, and the other two are arbitrary"; Fonseca–Grimus; King).

The only piece not found in the opened sources is the stated no-go over ALL relational functionals on E, with a mass-independent stationary set. That is "unknown – not found", not shown to be new.

Row verdicts:
| row | claim | verdict |
|---|---|---|
| 1 | Prop. obstruction | not found for the PAIR form; the three-fold form is BRSM |
| 2 | Thm single | special case of known (Sibirskii d=2) |
| 3 | Cor. external | reformulation of known (Lam / Fonseca–Grimus) |
| 4 | v1 Thm 2 | FALSE as a theorem: A = He² is a counterexample (Hermitian, built from He alone, non-degenerate on E, so it breaks SO(2)); v2's "definition of the class" wording is the correct one |
| 5 | stationary set | not found; elementary. Wording fix: "{0, π/4, π/2} (mod π/2)" is redundant |
| 6a | r* identity | not found; elementary |
| 6b | Jarlskog closure | special case of Jarlskog (m → n) |
| 6c | D4 | not found; elementary |
| 6d | T3b | reformulation of known (Lam partial pattern with ν2 the fixed column; TM2) |
| 6e | D11 | not found; FHZ reach (45, ±90) by symmetry, not by an extremum |

Gaps:
- no verbatim source for the PAIR statement;
- the orbit separation for compact U(2) (polynomial invariants → all relational functionals) is not cited;
- King Eq74 factor 2i vs 6i to re-read;
- the Lam and FHZ journal refs are unverified;
- the Sibirskii, Procesi and Harrison–Scott primaries were not opened.

Two paper fixes this implies, which I have NOT made since they touch the main text — your call:
- (a) the "(mod π/2)" wording at tex ~348;
- (b) the Discussion's \input of relation_known_results.tex.
Canon note for the lab/domain crew: v1 Thm 2 (registered docx, the flavor-alignment.md distillation) is false as stated. That is a canon correction, for Chazz.

## 2026-10-23 — novelty ruling applied
- **d73f874:**
  - (a) tex ~348 now reads "θ+ψ0 ∈ (π/4)ℤ; equivalently θ23 ∈ {0, π/4, π/2} on the fundamental interval [0, π/2]".
  - (b) The Discussion now \inputs relation_known_results.tex, keeping the "known facts assembled in spectral-triple language" sentence. A new paragraph records that v1's uniqueness theorem is false (A = He², Hermitian, built from He alone, generically non-degenerate on E) and that Prop. minimal states the class as a definition.
  - Note: the fragment cites labels from constructive_section.tex (sec:D2, sec:D3b, prop:jarlskog), so the main file must \input that too before it compiles clean.
- **dd8f065, (c)(4), checked myself from the arXiv abstract pages:**
  - **Lam 0708.3665 is Phys. Lett. B 656 (2007) 193, doi:10.1016/j.physletb.2007.09.032, NOT PRD 78 073015.** The bib entry is corrected and unmarked.
  - FHZ 1211.5560: arXiv shows NO journal-ref, so it stays marked "JHEP 07 (2013) 027 unverified" until INSPIRE confirms.
- **Gap items:**
  - #20 (literature-10): (1) a degenerate-pair quote from PDG 2024 / Bilenky–Petcov 1987; (2) the orbit-separation theorem (Schwarz 1975 / Procesi–Schwarz 1985); (3) the King Eq74 factor; (4) the FHZ journal-ref from INSPIRE; (5) Sibirskii/Procesi primaries, else "as cited in" Drensky.
  - #21 (theorist-6, after #20): the orbit-separation footnote for the Thm single proof; the explicit Jarlskog/King normalisation conversion; a grep check that no tex sentence claims novelty.
**Hash map:** unchanged except the one entry already recorded (3b55b00 → 602d419). New since then, newest first: dd8f065, d73f874, b7820bc, 2a25254, f50ce72, b634f2f, ee7122b, 1d4468a, 74d44c3, 599b24e, 116640b, cce4274, 91556dc, 46db69e, 848219f, 459ed7a, 434182f.

---
## 2026-10-23 — Novelty ruling (b7820bc) closed, fix-forward (branch crew/flavor-paper, not pushed)

**Hash map**
- d73f874: (a) ~348 now reads "θ+ψ0 ∈ (π/4)ℤ, i.e. θ23 ∈ {0, π/4, π/2} on the fundamental interval [0, π/2]". (b) `\input{relation_known_results}` added in the Discussion, with the closing "known facts assembled in spectral-triple language" sentence kept, plus the v1 Thm 2 sentence (A = He² counterexample; v2 states the class as a definition).
- dd8f065: Lam ref corrected from the arXiv journal-ref. 0708.3665 is **Phys. Lett. B 656 (2007) 193**, doi:10.1016/j.physletb.2007.09.032. PRD 78 073015 (in the ruling) is not this paper.
- a0e84da: (c)(2) footnote in the Thm single proof. It states the Sibirskii generators (as cited in Drensky Thm 2.10) and gives an elementary orbit-separation argument: approximate an invariant bump by a polynomial, then Haar-average; plus the Bloch/SO(3) concrete form. Schwarz 1975 and Procesi–Schwarz 1985 are **not cited**: neither was opened (search-only), and the elementary argument covers continuous functionals without them. (c)(3) Normalisation remark after Prop. jarlskog. (c)(4) FHZ unmarked. "On novelty" restated as the assembly, with no novelty claim.
- 890dbad: (c)(5) the Sibirskii full reference (Sib. Mat. Zh. 9 (1968) 152–164; transl. 115–124) and Procesi, Adv. Math. 19 (1976) 306–381. Both are marked "as cited in Drensky [primary not opened]".

**Gap verdicts**
1. *Pair statement* (degenerate pair ⇒ mixing unphysical): **not found verbatim**. Neither the PDG 2024 review nor Bilenky–Petcov RMP 59 could be opened. Rebelo 1603.01210 p.2 and Ioannisian–Kazarian hep-ph/0312379 p.2 (page images) state it only for full three-fold degeneracy. Row 1 stays "unknown – not found" for the pair form.
2. *Orbit separation*: closed by the elementary argument in the footnote. The only literature statement seen is Mendes–Radeschi arXiv:2012.07914, and that came via a summarising tool (FLAGGED), so it is not cited.
3. *King Eq. 74*: read verbatim from the page image. "I₁ = det[…] = (1/3)Tr([…]³) = 6iΔJ". The printed chain is internally inconsistent: for traceless X, Tr X³ = 3 det X, and 6iΔJ is the **Tr X³** value. Our det[Φ̂,Hν] = 2iJΔ(n)Δ(m²) matches Jarlskog PRL 55 Eq. 7 (det C = −2JΔΔ with [·,·] = iC). The remark writes both conversions and cites King for the trace form only.
4. *Journal refs*: FHZ 1211.5560 = JHEP 07 (2013) 027, doi:10.1007/JHEP07(2013)027. This was confirmed from INSPIRE, because the arXiv page has no journal-ref. Lam: see dd8f065 above.
5. *Primaries*: both are "as cited in" (above).

Notes saved: files/R2-novelty-gaps-literature.md (#20) and files/R2-novelty-gaps-theorist.md (#21).

**Not verified**: no LaTeX toolchain on this machine (pdflatex, latexmk and tectonic are all absent). Static check done: every \cite key has a bibitem, and the braces balance.

**Questions up**
- Q1. **Build file**: theta23_obstruction_paper.tex inputs relation_known_results.tex but **not** constructive_section.tex or r2_section.tex. relation_known_results references labels defined in constructive_section, so the main file as built has unresolved refs. Should I add `\input{constructive_section}` (and r2_section once filled) before the Discussion?
- Q2. **main.tex and limitative_draft.tex** (theorist-6's flag) still carry the v1 overclaims: "unique minimal relational operator" (main.tex ~1074, 1156, 1225, 1231, 1255), "establishes a new level of conceptual coherence", "uniquely determined atmospheric mixing angle". Are they in any build or submission? If yes, I propose the same fix-forward treatment. If no, I leave them untouched.

Members literature-10 and theorist-6 are dismissed and R2 is parked, pending relayed lab verdicts for D3b/D4 and D7–D12.

---
## 2026-10-23 — Q1/Q2 answers applied (branch crew/flavor-paper, not pushed)

- bb1413a: `\input{constructive_section}` and then `\input{r2_section}` sit before `\section{Discussion}`. Each fragment already opens with its own \section (sec:constructive, sec:r2). The static check (analysis/tex_static_check.py) stands in for a compile; its output (capture analysis/tex-static-check-2026-10-23.txt) is:
```
inputs expanded from theta23_obstruction_paper.tex: ['constructive_section', 'r2_section', 'relation_known_results']
labels 40, refs 22, bibitems 24, cite keys 24
refs with no label: none
labels defined more than once: none
cites with no bibitem: none
bibitems defined more than once: none
bibitems never cited: none
brace balance (expanded): 0
```
- Q2: main.tex and limitative_draft.tex are untouched. For Chazz, the v1 overclaim locations are main.tex ~1074, 1156, 1225, 1231, 1255 ("unique minimal relational operator", "establishes a new level of conceptual coherence", "uniquely determined atmospheric mixing angle") and the same claims in limitative_draft.tex. v2 records the correction in its Discussion.
- Staying parked for the D7–D12 batch.

---
## 2026-10-23 — NuFIT 6.0 theta13 / dm21 / dm3l columns

- 5cd2c19: analysis/data/nufit6_ranges.csv now carries s13sq, theta13, dm21 (1e-5 eV^2) and dm3l (1e-3 eV^2), each as bfp, +1s, -1s, 3s_lo, 3s_hi, for all four variant/ordering rows. I read them from the Table 1 page image (arXiv:2410.05380v2, PDF p.8 = printed -6-) and cross-checked every value against the PDF text layer. The rows and columns read are noted in the CSV header. dm3l = dm31 for NO and dm32 for IO, signs as printed.

---
## 2026-10-23 — R2 verdicts folded (B); D10 literature (A) running

- d1c0a5c: r2_section.tex now has one verdict paragraph each for D7, D8, D9, D11 and D12. D12 is a new subsection. Each number names its lab capture (crew/flavor-obstruction-lab 3ffb5ed), and the Valencia 2020 numbers are labelled. D10 = gate pending. The summary table waits on D10. Static check is clean (analysis/tex-static-check-2026-10-23-r2-verdicts.txt). Two notes: D8 order 1 is reported as roughly 2-sigma tension (chi2 8.90-12.69, 2 dof) with the lab's reason for not calling it special. D11 delta=90 excluded is read off nufit6_ranges.csv: the NO ranges wrap only to 62 and 4 deg.
- (A) item #22 is with literature-11. files/R2-D10-literature.md follows when it reports.

---
## 2026-10-23 — (A) D10 literature done; parking

Hashes: d1c0a5c (B, R2 verdicts folded) · (A) is a crew file, not a Flavor commit: files/R2-D10-literature.md (#22 literature-11 page-image reads, plus the lead's own read of Ishiguro–Kobayashi–Otsuka 2011.09154).

Corrections and gaps:
- 1907.11714 is plain modular A4 without gCP. The fixed-point paper is Ding–King–Liu–Lu 1910.03460, which has no gCP.
- No model with gCP evaluated exactly at a fixed point was found. At gCP-conserving tau the phases are trivial (NPPT p.12).
- theta23 is not fixed by the residual symmetry in any model read.
- Stabilisation: the minima are near omega and CP-breaking (NPP 2022), and the flux landscape clusters at omega with 40–62% (IKO Tables 2–3).
- The NuFIT versions used by the fits were not seen. The IKO journal ref is unconfirmed (it came via a summary).

---
## 2026-10-23 — D10 folded; R2 battery complete in the paper (final item)

- e5c27d0: D10 is the last candidate paragraph in r2_section.tex (gate c94252c, reproduction e25408b; captures run-2026-10-23.txt / -validate.txt). Bib entries added: feruglio2017, nppt2019, npp2022, with titles and authors from each PDF's p.1; journal refs omitted (unverified). Static check clean: analysis/tex-static-check-2026-10-23-d10.txt.
- **One deviation from the relayed verdict.** Your relay said "the seesaw examples have free g that enters theta23". The lab README table and the capture say the opposite. In Ex. 2, g enters only as an overall scale: s23^2 = 0.651/0.349 for every g tested, with two massless states. In Ex. 3 the spectrum is rank one for g2/g1 != 0, so its angles are not physical. The paper follows the record. Overrule if the lab meant something else.
- The 2201.02020 sentence is quoted verbatim (checked in the PDF text). Its context is the charged-lepton and quark mass hierarchies, and the paper says so. It is not quoted as "cannot account for the data".
- Left as is: the R2 summary table (sec:r2-summary) still reads [pending]. Fill it on your word; otherwise it stays.

**Branch crew/flavor-paper hash map, this round:** d73f874 (stationary set, Relation-to-known-results, v1 Thm 2 sentence) · dd8f065 (Lam ref) · a0e84da (orbit-separation footnote, On novelty as assembly, FHZ, Jarlskog/King normalisation) · 890dbad (Sibirskii/Procesi as cited in Drensky) · bb1413a (\input constructive + r2; static check) · 5cd2c19 (NuFIT 6.0 theta13/dm21/dm3l) · d1c0a5c (R2 verdicts D7–D9, D11, D12) · e5c27d0 (D10). Not pushed.

- 0535b07: D10 per the manager's correction. The text never carried the 45/+-90 or "FHZ at tau = i" pattern: grep finds +-90 only in D11, where it belongs. Added: the gCP pair (tau, -tau*) returns delta ~295 and ~65 deg as equally valid, so the fit cannot pick D11's sign. The verdict line now ends "a gCP model at a fixed point gives delta = 0 or pi". Static check clean. Final hash map: d73f874 · dd8f065 · a0e84da · 890dbad · bb1413a · 5cd2c19 · d1c0a5c · e5c27d0 · 0535b07 (crew/flavor-paper, not pushed).
