# R2-theorist — octant paragraph and D11 extremum proposition (item #13)

Author: theorist-4 (crew flavor-obstruction-paper). Saved verbatim by the lead from the author's crew note; the read-only role cannot write files. Paper only; nothing run.

---

R2-theorist (#13), for files/R2-theorist.md. Paper only; nothing run. Hand arithmetic marked (hand).

## 1. Statement
(1) The octant fact, using NuFIT 6.0 Table 1 (analysis/data/nufit6_ranges.csv). (2) D11: is J = s12 c12 s13 c13^2 s23 c23 sin(delta) extremal at theta23 = 45 and delta = +-90; is there any second-order octant preference; are the other natural invariants extremal at the data point. R2 runs outside canon by mandate (brief-R2 l.1), so none of this is Alignment canon.

## 2. Assumptions
A1 PDG parametrisation of U (standard; not stated in content/; flavor-alignment.md uses it). A2 theta12, theta13 in (0,90), so A = (1/8) sin2theta12 sin2theta13 cos theta13 > 0. A3 NuFIT 6.0 values as in the csv. A4 For (iv) only: a_i = m_l^2 with a_e, a_mu << a_tau, and b_k = m_nu,k^2.

## 3. Derivation
(1) Octant, side-appropriate 1-sigma error (algebra on csv):
IC19 NO 48.5 -0.9: (48.5-45)/0.9 = 3.9. IC19 IO 48.6 -0.9: 4.0. IC24 NO 43.3 +1.0: 1.7. IC24 IO 47.9 -0.9: 3.2. Brief values 3.9 and 1.7 confirmed.
CAVEAT, this is the main finding of (1): every printed 3-sigma window contains 45 (IC19 NO [41.0,50.5], IC19 IO [41.4,50.6], IC24 NO [41.3,49.9], IC24 IO [41.5,49.8]). So the 3.9 and 4.0 figures are linear Gaussian extrapolations that the table itself contradicts: the chi^2 profile is non-Gaussian (second-octant minimum), and 45 is excluded at less than 3 sigma in all four rows. Do not print 3.9 sigma as an exclusion. The true Delta chi^2(45) needs the NuFIT chi^2 profiles (calculator, sec. 7). Second wording fix: the 3-sigma union is [41.0, 50.6], i.e. |theta23-45| <= 5.6, so the claim needs to say: best fits within 3.6 deg, 3-sigma within about 6 deg. The claim |theta23-45| <~ 5 is slightly too tight at 3 sigma (50.6 in IC19 IO).

LaTeX:
\paragraph{The octant is not established.} NuFIT~6.0 \cite{NuFIT6} gives, for normal ordering, $\theta_{23} = 48.5^{+0.7}_{-0.9}{}^\circ$ (IC19, without SK-atm) and $43.3^{+1.0}_{-0.8}{}^\circ$ (IC24, with SK-atm); for inverted ordering, $48.6^{+0.7}_{-0.9}{}^\circ$ and $47.9^{+0.7}_{-0.9}{}^\circ$. The two data sets put the normal-ordering best fit in opposite octants. Measured in the side-appropriate $1\sigma$ error, $45^\circ$ lies $3.9\sigma$, $1.7\sigma$ (NO) and $4.0\sigma$, $3.2\sigma$ (IO) from the best fit; but the $\chi^2$ profile is not Gaussian, and every published $3\sigma$ interval ($[41.0,50.5]^\circ$, $[41.3,49.9]^\circ$ NO; $[41.4,50.6]^\circ$, $[41.5,49.8]^\circ$ IO) contains $45^\circ$. What the data establish is that $\theta_{23}$ is near-maximal: the best fits lie within $3.6^\circ$ of $45^\circ$ and the $3\sigma$ ranges within $[41.0^\circ, 50.6^\circ]$. The octant, and whether $\theta_{23}$ differs from $45^\circ$ at all, are not established. A mechanism is therefore scored on two questions: does it give near-maximal $\theta_{23}$, and which octant does it pick.

(2) D11. Algebra throughout (identities exact).
Write J = A sin(2 theta23) sin(delta).
(i) dJ/dtheta23 = 2A cos(2theta23) sin(delta); dJ/ddelta = A sin(2theta23) cos(delta). For theta23 in (0,90), sin2theta23 != 0, so both vanish iff cos delta = 0 and cos 2theta23 = 0: the only interior critical points are (45, +90) max J = +A and (45, -90) min J = -A. QED.
(ii) At those points both first derivatives vanish (shown in i), so no first-order term in theta23-45: no first-order octant preference.
(iii) Exact, all orders: x = theta23 - 45, y = delta + 90. J = -A cos(2x) cos(y). Second order: J = -A + 2A x^2 + (A/2) y^2 + O(4); mixed derivative d2J/dx dy = 2A sin2x sin y... = 0 at the point: no cross term, so a delta offset induces no theta23 shift. Stronger and exact: J is even in x for every delta, theta12, theta13 (s23 c23 = (1/2) sin 2theta23 is invariant under theta23 -> 90 - theta23). Hence at fixed delta with sin delta != 0, the theta23 extremum of J is exactly 45, independent of delta; at sin delta = 0, J = 0 identically and selects nothing. A soft principle (|J| >= (1-eps) A) gives the band cos2x |sin delta| >= 1-eps, symmetric in x. Proposition: no principle that is a function of J alone can select an octant, at any order. Only the width is delta-dependent: the theta23-curvature of J is 4A|sin delta|. At the data |sin delta| = 0.052 (IC19 NO, 177), 0.53 (IC24 NO, 212), 0.97 (IC19 IO, 285), 0.998 (IC24 IO, 274) (hand). So maximal CP sits near the IO best fits and far from the IC19 NO best fit, though 270 is inside every 3-sigma delta range.
LaTeX:
\begin{proposition} Let $J = A \sin 2\theta_{23}\sin\delta$, $A = \tfrac18\sin2\theta_{12}\sin2\theta_{13}\cos\theta_{13} > 0$. (i) On $\theta_{23}\in(0,\pi/2)$ the critical points of $J$ are exactly $(\pi/4,\pm\pi/2)$. (ii) For every fixed $\delta$ with $\sin\delta\neq0$, $\partial_{\theta_{23}}J = 0$ iff $\theta_{23} = \pi/4$. (iii) $J(\pi/2-\theta_{23},\delta) = J(\theta_{23},\delta)$ identically, so no function of $J$ distinguishes the octants. \end{proposition}
\begin{proof} $\partial_{\theta_{23}}J = 2A\cos2\theta_{23}\sin\delta$ and $\partial_\delta J = A\sin2\theta_{23}\cos\delta$; on the open interval $\sin2\theta_{23}\neq0$, which gives (i) and (ii). (iii) holds since $\sin2(\pi/2-\theta_{23}) = \sin2\theta_{23}$. \end{proof}

(iv) Other invariants.
(a) Weak-basis CP-odd Tr([He,Hnu]^3) = 6i Delta_e Delta_nu J (overall sign convention-dependent, literature-3 flag). Mass factors are angle-free, so the Proposition applies verbatim: extremal in theta23 at exactly 45, octant-blind. Rigorous.
(b) Majorana-type invariants, e-row (PDG with diag(1, e^{i a21/2}, e^{i a31/2})): S1 = Im(U_e1 U_e3*) = c12 c13 s13 sin(delta - a31/2); S2 = Im(U_e2 U_e3*) = s12 c13 s13 sin(delta + (a21-a31)/2). Both are independent of theta23 identically, so they carry no theta23 information. Extremality at the data point is undefined: a21 and a31 are unmeasured. Rigorous (algebra). Weak-basis Majorana invariants (e.g. Im Tr(He Hnu He* M M*)-type, Branco-Lavoura-Rebelo) do depend on theta23: not derived here, gap Q2.
(c) CP-even: in the e-diagonal basis, Tr([He,Hnu]^2) = -2 sum_{i<j} (a_i-a_j)^2 |Hnu_ij|^2 (exact). Invariant under Hnu -> Hnu + c 1, so independent of the absolute neutrino mass scale (exact). Leading order in a_mu/a_tau (about 3.5e-3) and Dm21/Dm31 (about 0.03): -1/2 Tr([He,Hnu]^2) = a_tau^2 Var_tau(b), and Var_tau = p3(1-p3) Dm31^2 with p3 = |U_tau3|^2 = c13^2 c23^2. That is maximal at |U_tau3|^2 = 1/2: s23^2 = (1-2 s13^2)/(2(1-s13^2)), about 0.489 at s13^2 = 0.022, theta23* about 44.4 deg (hand). That is first octant and delta-independent at leading order, with delta entering at O(Dm21/Dm31). So a maximal-misalignment principle, max ||[He,Hnu]||_F, is NOT octant-blind. It picks theta23 slightly below 45. This is algebra plus a stated approximation, not validated: calculator check C2. It is not extremal at any best fit; it lies inside all four 3-sigma windows and about 1 sigma from IC24 NO 43.3 (hand).
(d) Tr(He^m Hnu^n) = sum_i a_i^m (Hnu^n)_ii and Tr((He Hnu)^2): at the same leading order they are monotone in p3, so their only theta23 extremum is at the endpoints 0 or 90. Not extremal at data. They have interior extrema only through the dropped O(Dm21) and O(a_mu) terms; not derived.

## 4. Register
Rigorous (exact algebra, needs no numerics): Proposition (i)-(iii); J-proportionality of Tr([He,Hnu]^3) (up to sign convention); theta23-independence of S1, S2; the commutator-trace identity and its mass-scale invariance. Arithmetic, statistically approximate: the n-sigma figures (Gaussian extrapolation, contradicted by the 3-sigma windows). Leading-order and unvalidated: theta23* about 44.4 from max ||[He,Hnu]||; monotonicity in (d). Interpretive: maximal CP or maximal misalignment as a selection principle.

## 5. EQUATIONS.md
No row touched: R2 is outside canon. If maximal misalignment survives C2, it would be a new candidate, never a restatement. Nearest row: the flavor-alignment.md misalignment functional (Thm 3 withdrawn v0.20.1), which minimises misalignment; this maximises it. Propose nothing until C2.

## 6. Gaps (for Chazz)
Q1 Is maximal J or maximal misalignment meant as a principle (extremise over theta23 only, over all angles, or over delta too)? Over all angles, max ||[He,Hnu]|| gives only the curve |U_tau3|^2 = 1/2, not theta23 alone. Q2 Which Majorana invariants count as natural (e-row S1, S2 versus weak-basis)? Q3 Can we get the NuFIT chi^2 tables (needed for a true Delta chi^2(45))?

## 7. Calculator setup
C1 Delta chi^2(theta23 = 45), profiled, from the NuFIT 6.0 tables for all four rows. Refutes the 3.9 sigma figure if Delta chi^2 < 9 (expected from the windows).
C2 Exact argmax over theta23 of ||[He,Hnu]||_F^2 at NuFIT theta12, theta13, delta, Dm21, Dm31 (four rows), PDG lepton masses, m_lightest in {0, 0.01, 0.05} eV. Ground truth: invariance in m_lightest to machine precision, and at a_e = a_mu = 0, Dm21 = 0 the argmax equals arccos(1/(sqrt2 c13)). The prediction is refuted if the full argmax leaves the first octant, or moves more than about 1 deg from 44.4, or depends on delta at O(1).
C3 Null: the fraction of Haar-random U whose |U_tau3|^2 lies within the data distance of 1/2, against the same for |U_mu3|^2 = 1/2.