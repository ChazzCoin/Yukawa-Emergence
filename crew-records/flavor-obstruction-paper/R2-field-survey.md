# R2-field-survey — existing flavor models and fits (item #14)

Author: literature-5 (crew flavor-obstruction-paper). Saved verbatim by the lead from the author's two crew notes; the read-only role cannot write files.

---

R2-field-survey 1/2 (#14). ALL fetches below went through a summarizing tool (so quotes/numbers = tool transcription) unless stated. Status: OPENED(summ.) = page fetched via summarizer; SEARCH = search snippet only; UNSOURCED = my background, not checked this session.

0. CORRECTIONS TO THE BRIEF
- arXiv:2409.11824 = Cao, Quyen, Hong Van, Nath, Ngoc, 'On Precision of the Leptonic Mixing Angle theta23 and its Implications for the Flavor Models' (Sep 2024, rev Oct 2024). Hyper-K+DUNE simulation: octant resolved for all values inside 1 sigma; hard to separate maximal from octant when 0.48<=sin^2 th23<=0.54; if CP conserved HK+DUNE eliminate most discrete non-Abelian models. https://arxiv.org/abs/2409.11824 OPENED(abs summ.)
- arXiv:2603.06934 is NOT a Hyper-K/DUNE discrimination paper. It is V. Barger, lattice-flavon (Froggatt-Nielsen B-lattice) lepton model: two branches (th23~43 deg, delta~286 deg) and (th23~46 deg, delta~299 deg); near-TBM neutrino sector + charged-lepton rotations + approx Z2 mu-tau symmetry; distinguishing needs DUNE/HK/IceCube/JUNO octant+delta. https://arxiv.org/abs/2603.06934 OPENED(abs summ.)
- The '15 expressions, ONE charged-lepton rotation, theta23 within 0.2 deg' claim is in Denton, Gehrlein, Truelson, 'The Future of Lepton Flavor', arXiv:2606.05291 (submitted 3 Jun 2026, exists): Sec 6.2 Discussion, verbatim as relayed: 'All 15 possible expressions containing only one rotation in the charged lepton sector predict theta23 within 0.2 deg.' No citation attached to that sentence in the tool's view; the class is 'charged-lepton corrections' = U_PMNS = U_e^dagger U_nu with U_e a single 2x2 rotation. NOT verified what 'within 0.2 deg' is relative to (I read it as 'of each other/of maximal' from the tool; ambiguous - the tool wrote 'of maximal'). https://arxiv.org/abs/2606.05291 and https://ar5iv.labs.arxiv.org/html/2606.05291 OPENED(summ.)
- Closest primary paper I found on single-rotation charged-lepton sum rules: 'Constraining the Neutrino Mixing Matrix via Single-Sector Charged-Lepton Rotations in the JUNO Precision Era', arXiv:2603.29876 (Symmetry 18(6) 954, per search). It studies 3 cases (I: (1,2); II: (1,3); III: (2,3) rotation), NOT 15 expressions; Case III Eq29: theta23^nu ~ theta23 - theta23^l; leading order small theta13: sin^2 th23^nu > sin^2 th23. Inputs NuFit 6.1 (with JUNO, Nov 2025) and NuFit 6.0. https://arxiv.org/html/2603.29876 OPENED(summ.). So the '15' count is Denton et al.'s, and the 3-case paper is a different one; link between them not shown.
- Global-fit octant: sin^2 th23 = 0.473 (+0.023 -0.013 per search snippet; fetch gave 1 sigma 0.460-0.496), 3 sigma 0.437-0.581, first (lower) octant preferred at 1.1 sigma (down from 1.6 sigma in their 2021 analysis); NO; solar+KamLAND+LBL+short-baseline reactor+SK/IC atmospheric. Source: Capozzi, Giare, Lisi, Marrone, Melchiorri, Palazzo, 'Neutrino masses and mixing: Entering the era of subpercent precision', arXiv:2503.07752 (PRD 111 (2025), journal per search). https://arxiv.org/html/2503.07752v1 OPENED(summ.). The abstract the tool returned was truncated, authors from search only.
- Denton et al Table 1 (NuFit 6.0 + JUNO 2025 compile): sin^2 th23 = 0.561 +0.012 -0.015, 3 sigma 0.430-0.596 (NO) - UPPER octant, and this is the SAME paper that says the octant is unknown. So the 2026 review table and the Capozzi 2025 fit disagree on octant. OPENED(summ.). (Matches the known IC19/IC24 split in NuFit 6.0, R1.)

1. PREDICTIVE vs DESCRIPTIVE (theta23 | octant | delta corr | free params vs observables | killed by | source)
- mu-tau reflection (exact): theta23=45 deg, delta=pi/2 or 3pi/2 | none (maximal) | delta=+-90 | 0 free for theta23,delta | theta23 != 45 at >~1 deg or delta far from +-90 | Harrison-Scott PLB547 (2002) 219, hep-ph/0210197 | SEARCH (snippet: theta23=pi/4, delta=pi/2 or 3pi/2). Breaking size/form: UNSOURCED this session. For a breaking study see 1409.8600 (radiative corrections, appeared in search) SEARCH.
- TM1: sin^2 th12 = (1/3)(1-3|Ue3|^2)/(1-|Ue3|^2)  [A-R Eq23]; cos delta = -cot2th23 (1-5s13^2)/(2 sqrt2 s13 sqrt(1-3s13^2)) [A-R Eq24, King Eq51, Costa-King Eq5.64 as relayed] | theta23 free input; delta fixed by theta23 | strong delta<->th23 correlation | 0 free (relation only: 2 observables th12, delta predicted from th13, th23) | th12 outside ~0.318, or delta vs cos-delta rule | 0812.0436, 1701.04413, 2307.13895 | OPENED(summ.)
- TM2: sin^2 th12 = 1/(3(1-|Ue3|^2)) >= 1/3 [A-R Eq15] | (1/2 - sin^2 th23)^2 ~ cos^2 delta |Ue3|^2/2 [A-R Eq18, tool-garbled] so octant not fixed, deviation bounded by delta | | sin^2 th12 > 1/3 vs data ~0.31 -> strained | 0812.0436 | OPENED(summ.)
- Littlest Seesaw CSD(n) (n=3 etc): see part 2/2 section 4. theta23 ~45.5+-2.4 (n=3), 41.4 (n=1+sqrt6), 51.5/51.0 (n=2.5); 2 inputs (r=m_b/m_a, eta) + overall m_a (3 params) vs 6 observables(3 angles, delta, 2 masses) + 1 Majorana | killed: delta vs ~272 or ~88 deg (n=3); octant | 2307.13895; 1512.07531 (King, theta23=46, delta=-87, theta12=34, theta13=8.7 per abstract) | OPENED(summ.)
- Direct discrete + gCP: A5+gCP (Turner): theta23 in ranges 38.2-45 (lower), 45-53.3 (upper) or both by case; delta and alpha21 bounded within +-90; CP-conserving delta needs 40.6<=th23<=44.2 (lower) or 45.6<=th23<=48.4 (upper) | Gnu=Z2xCP, Ge=Z2 | correlations of th12,th13,th23 with delta | 0 continuous (discrete choices; a single free angle parameter theta usually, UNSOURCED for this paper) | 1507.06224 | OPENED(summ.)
- Modular, generic tau: A4 NPPT Model I(a): sin^2 th23 > 0.54, delta=+-(50-180) deg, sum m~145 meV | upper | corr th23-delta, sum m, m_ee | complex tau (2 real) + g, g ratios etc vs 7+ observables | sum m, theta23 octant | 1808.03012 | OPENED(summ.). S4 (NPPT 1811.04933): all models predict correlations sin^2 th23 with delta, sum m, |<m>|; best case A s23^2=0.551 | OPENED(summ.). Feruglio 1706.08749 examples fit poorly on th13 (R2 note).
- Modular at fixed points: tau=i residual Z2 symmetry/antisymmetry for Weinberg operator: S_nu^T M_nu S_nu = (-1)^k M_nu; odd k -> det M=0 and Tr M=0 (one massless, possible degenerate pair); gives one fixed PMNS column (five viable: trimaximal, tri-large, ...); some solutions 'maximal CP violation for theta23=pi/4', Gamma_16 excludes max CPV | Kashav-Patel arXiv:2410.00565 | OPENED(summ.). tau=omega: texture structures (Kobayashi-Tanimoto 2307.03384 abs) OPENED(abs summ.). Delta(96) modular LS: 35 viable breaking patterns (21 NO, 14 IO), beyond CSD(n), new fixed columns; the sin^2 th23 0.435-0.583 figure NOT confirmed (ar5iv page empty; abstract has no number) 2607.07865 OPENED(abs summ.).
- Anarchy: Haar-flat; sin^2 th12, sin^2 th23, cos^4 th13 uniform; no theta23 prediction; KS 64% (2003) / 44% (2012, >24% in 3 sigma) | | | 0 | needs a precision measurement showing clustering: e.g. th23 well-constrained at 45+-1 would be 'special' | hep-ph/0301050, 1204.1249 | OPENED(summ.) (earlier note)
- QLC: theta23 + theta23^CKM ~ 45 (Raidal Eq1) => theta23 ~ 45 - 2.4 = ~42.6 (my arithmetic with PDG sin th23 0.04183 -> 2.40 deg) | lower | none stated | few | theta23 octant | hep-ph/0404046 | OPENED(summ.)
- Texture zeros: Majorana m_mumu=0 clear lower-octant preference; m_mutau=0, m_tautau=0 clear upper-octant preference (Denton et al Table, depends on absolute mass scale) | per class | | 1-2 conditions removed | octant + m_beta/ sum m | 2606.05291 | OPENED(summ.)
- GUT relations: UNSOURCED this session (not searched). Mass sum rules: 'make no specific predictions on theta23' (Denton et al.) OPENED(summ.)

---

R2-field-survey 2/2 (#14). Same status conventions (summarizing fetch tool; OPENED(summ.)).

2. OPEN PROBLEMS / HINTS (field's words as the tool relayed them; I found little verbatim text, so this table is thin)
- 'The flavor puzzle remains one of the biggest open questions in particle theory to date...' (charged fermion masses span six orders of magnitude, no obvious relations) - Denton-Gehrlein-Truelson 2606.05291 Sec 1 - OPENED(summ.)
- Octant of theta23 is one of 'four major unknowns' (with mass ordering, delta, absolute mass scale per tool wording) - same paper - OPENED(summ.)
- Octant evidence has weakened: 1.6 sigma (2021) -> 1.1 sigma first-octant preference (2025) - Capozzi et al 2503.07752 - OPENED(summ.). Denton Table 1 has upper-octant best fit 0.561 - the two 2025-26 global compilations disagree.
- Ambiguity zone: 0.48<=sin^2 th23<=0.54 cannot separate maximal from octant at HK+DUNE alone - Cao et al 2409.11824 - OPENED(summ.)
- Hints of correlations: modular A4 Model I(a) ties sin^2 th23>0.54 to delta=+-(50-180); S4 ties th23 to delta, sum m, <m>; TM1 ties cos delta to th23 exactly (Eq above). Barger: (theta23,delta) two-branch plane. None of these is a data hint; they are model correlations.
- Not found: any source stating 'nobody explains why theta23 is near maximal' verbatim. I did not locate one; do not quote it.

3. REFRAMINGS
- Single charged-lepton rotation: Denton et al: all 15 one-rotation expressions give theta23 within 0.2 deg (vs what: unverified, see 1/2). Case III (23-rotation) of 2603.29876: theta23^nu ~ theta23 - theta23^l, so theta23 is split into a neutrino-sector part and a charged-lepton angle; leading order sin^2 th23^nu > sin^2 th23. Reframe: 'what fixes theta23' becomes 'what fixes theta23^nu and theta23^l separately', the second possibly small.
- tau=i residual symmetry (Kashav-Patel 2410.00565): the symmetry fixes one column of U_PMNS and a massless/degenerate pair, so theta23 appears as a consequence of column + gCP not an input; 'most viable solutions feature large mixing angles' (tool). Reframe: theta23 near-maximal as the fixed-column geometry, not tuned.
- Exact relevance to the R1 no-go (near-degenerate m2~m3 pair, orientation undefined): Kashav-Patel's real-matrix degeneracy statement ('degeneracy between the other two is possible for real mass matrices') is the one place I saw a symmetry-driven degeneracy of an nu pair. OPENED(abs summ.).

4. LITTLEST SEESAW CSD(n), for r2_inputs.csv
Sum rules (as printed per tool, Costa-King arXiv:2307.13895, Universe 9 (2023) 472, https://ar5iv.labs.arxiv.org/html/2307.13895):
- TM1 (LS preserves first column of TBM): sin^2 th12 = 1/(3 cos^2 th13) [tool: 'first rule'] 
- Eq5.64: cos delta = -cot 2th23 (1-5 sin^2 th13) / (2 sqrt2 sin th13 sqrt(1-3 sin^2 th13)); the tool's printout dropped the sqrt on the 2 sqrt2, I added the standard form from A-R Eq24 / King Eq51; applies 'to LS models n~3'.
- Solar sum rule general form: cos delta = [tan th23 sin^2 th12 + sin^2 th13 cos^2 th12/tan th23 - ...]/[sin2th12 sin th13] (TRUNCATED by tool; model-dependent theta12^nu, TB: sin th12^nu=1/sqrt3). Original needed.
Mass matrix: m^nu = m_a (0 0 0;0 1 -1;0 -1 1) + m_b e^{i eta} (rank-1 matrix with parameter n) [tool paraphrase; n-dependent explicit vector not seen].
Best-fit/predictions (Costa-King Tables 3-5, NO, 'case' = two eta branches):
- n=3: r=m_b/m_a=0.100+-0.008; eta=2.11+-0.15 -> th12 34.32(+0.20-0.24), th23 45.5(+2.3-2.4), delta 272.2(-11.0+9.6); eta=4.17+-0.15 -> th12 34.32(+0.20-0.25), th23 45.5(+2.3-2.4), delta 87.9(+11.0-9.6).
- n=1+sqrt6~3.45: r=0.072+-0.004; eta=2.42+-0.16 -> th12 34.36(-0.21+0.18), th23 41.4+-2.6, delta 253.8(+11.7-13.8); eta=3.87+-0.16 -> th23 41.5+-2.6, delta 105.7(+13.7-11.6).
- n=2.5: r=0.15+-0.01; eta=1.5+-0.2 -> th12 34.31, th23 51.5(-2.2+1.9), delta 299.9(-9.9+9.2); eta=4.7+-0.2 -> th23 51.0(-2.3+2.0), delta 63.6(-9.3+10.1). (tool and abstract text call n=2.5 disfavoured; the disfavour statement is the tool's, not verified.)
- th13 and masses for these rows: NOT retrieved.
Other papers:
- King, 'Littlest Seesaw', arXiv:1512.07531 (JHEP 1602 (2016) 085): abstract benchmark th12=34, th13=8.7, th23=46, delta=-87 deg; NO, massless lightest; S4xZ3xZ3'. OPENED(abs summ.)
- Chen-Ding-King-Li, 'A New Littlest Seesaw Model', arXiv:1906.11414 (tri-direct CP; S4xZ4xZ9): Eq1/11: m^nu = m_a (0 0 0;0 1 -1;0 -1 1) - i (m_s/4)(4 10 -2;10 25 -5;-2 -5 1) for x=-1/2, eta=-pi/2; Eq77 Dirac columns (0,1,-1)^T and (1,-1/2,5/2)^T; Table 1 best fit: m_a=22.366 meV, r=m_s/m_a=0.145, eta=-pi/2, x=-1/2; predictions sin^2 th13=0.0220, sin^2 th12=0.318, sin^2 th23=0.599 (second octant, ~50.8 deg by my arithmetic), delta/pi~-0.354, alpha21/pi~0.317, m2~8.670 meV, m3~50.167 meV. https://ar5iv.labs.arxiv.org/html/1906.11414 OPENED(summ.). NOTE: sin^2th23=0.599 is 'final row highlighted'; check the table.

5. RECOMMENDATIONS (ranked by sharpness: few params, many predictions; all must be gated on source tables first)
1. Littlest Seesaw CSD(n=3 and 1+sqrt6, 2.5) and 1906.11414 model: 2-3 continuous inputs, predicts th12, th13, th23, delta, masses at once; th23 = 45.5 / 41.4 / 51.5 per n. [2307.13895; 1906.11414; 1512.07531]
2. mu-tau reflection exact as the zero-parameter baseline (th23=45, delta=+-90) [hep-ph/0210197, SEARCH].
3. TM1 exact (cos delta from th23, th13): delta<->th23 correlation, no new parameter [0812.0436 Eq24; 1701.04413 Eq51].
4. Tau=i residual (anti)symmetry (Kashav-Patel): zero continuous input; fixed column + massless/degenerate pair; closest to R1's degenerate-pair question [2410.00565].
5. Single charged-lepton 2x2 rotation, U_PMNS=U_e^dagger U_nu (Denton Sec 6.2; 2603.29876 Case III th23^nu ~ th23 - th23^l): one free angle.
6. A5+gCP (Turner), discrete choice + one parameter: th23 ranges by case [1507.06224].
7. A4/S4 modular at generic tau (NPPT): 2 real tau + few g's: sin^2 th23>0.54 [1808.03012, 1811.04933].
8. Texture zeros (octant by class) [2606.05291].
9. Anarchy/QLC as nulls (already D8/D9).
UNVERIFIED: whether 2607.07865's sin^2 th23 0.435-0.583 range exists (page unreadable); GUT-relations row; breaking size for mu-tau reflection; the exact 'within 0.2 deg' reference value.
