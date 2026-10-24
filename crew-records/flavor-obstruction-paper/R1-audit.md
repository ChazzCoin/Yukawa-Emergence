# R1 audit — crew flavor-obstruction-paper (Flavor repo only)

The cross-repo picture comes from the manager: /Users/chazzromeo/rAI/galt-crews/flavor-obstruction/files/R1-audit.md and plan.md.
Work home: /Users/chazzromeo/rAI/galt-crews/flavor-paper-wt, branch crew/flavor-paper (faf2d83 drafts in; 5223dc3 audit).

Gate: verify_theta23_obstruction.py, ALL PASS (analysis/verify-run-2026-10-23.txt). The gate is weaker than it looks: its Prop 1 check loops over t without using it, and every check runs on the 2x2 restriction only.

Findings (detail: AUDIT-theta23-paper.md; sources: files/R1-literature.md):
- (i) When E is an eigenblock of H_nu, the trace-type invariants keep the single-cos 2(theta+psi0) form at 3x3. ||[He,Hnu]||^2 does not: the e-row adds a second orientation, so the stationary set is not frozen (maxima shift ~3.3 deg, minima nonzero). Thm nogo holds only for the 2x2 restriction.
- (ii) The premise "E = the mu-tau flavor plane" fails when theta12 and theta13 are nonzero. Rotating inside the true degenerate E changes theta12 and theta13 too. With m2 = m3, all three angles are frame labels, so the obstruction is "mixing unobservable", not "theta23 undetermined". Any nonzero splitting eps makes the rotation angle physical, linear in eps.
- (iv) ½{He,Hnu} stays in the linear class. HeHnuHe leaves it, but for positive semidefinite (PSD) matrices its minima stay at alignment; its interior maximum depends on the masses.
- (v) The closed form is re-derived to 1e-12. Hidden assumption: real eigenbases. With a complex H_e the minimum is > 0.
- (vi) v2 withdraws Thm 3 and v1's numerical theta23 claim, demotes Thm 2 uniqueness to a definition of the class tested, and drops v1 §7-§8.
- (iii) Literature has the sources with versions. The numbers are owed in item #3.
- Flag only: Domain.Alignment content/vocabulary.md:166 is stale. The lab crew owns that fix.
