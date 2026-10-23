# Audit of `theta23_obstruction_paper.tex` — gate + sections (i), (ii), (iv), (v)

Crew `flavor-obstruction-paper`, item #1 (R1), data-analyst. Date 2026-10-23.
Env: python 3.13.7, numpy 2.4.4 (`/opt/homebrew/bin/python3`); sympy is not installed, so (v) is numeric to 1e-12.
Every number below is in a committed capture under `analysis/` (named per line). Scripts: `analysis/README.md`.
All inputs marked *illustrative* are magnitudes chosen to be realistic in size, **not** NuFIT values and not tied to a NuFIT version.
Sections (iii) and (vi) of the brief are not part of this item and are not covered.

## Gate — PASS
`python3 verify_theta23_obstruction.py` → every check PASS, "ALL CHECKS PASS", exit 0
(`analysis/verify-run-2026-10-23.txt`). Deviations: Thm 1 overlap 1.78e-15, Thm 3 closed form 1.24e-14.
What the gate does not show: the Prop 1 check loops over `t` but never uses it (`Hnu_deg` is constant), so it is
trivially 0; every check is on the 2×2 restriction only. The verifier does not exercise anything in (i), (ii), (iv).

## (i) Full 3×3 functionals — `analysis/full3x3_fourier-run-2026-10-23.txt`
Setup: W(θ)=1⊕R(θ) on the fixed plane E; N=720 θ-samples; "function of one cos 2(θ+ψ₀)" is tested as
(no odd harmonics) AND (reflection symmetry f(p+x)=f(p−x) for some p). Ground truth (H_e block-diagonal on E, the paper's
2×2 case) reproduces the paper: odd-harmonic content 4e-16 of range, k=2 amplitude of Tr(HeHν) = 2.491160e-02 = ½δeδν exactly, reflection
asymmetry 5e-7 (grid limit); stationary set of ‖[He,Hν]‖² = {min, max, min, max} 45° apart, minima value 2e-12.

**Verdict by case**
1. **E is an eigenblock of H_ν (the paper's premise; H_e tilted, θ12,θ13 ≠ 0, with or without CP phase 197°):**
   - Tr(HeHν), Tr((HeHν)²), eigenframe misalignment: **hold**. No odd harmonics (≤7.7e-16 of range), reflection asymmetry ≤7e-7 (Tr(HeHν), misalignment), ≤2.9e-5 (Tr((HeHν)²)); 300 random draws (seed 20261023): odd/range max 2.3e-15, frac(odd>1e-6)=0.
   - **‖[He,Hν]‖² fails** the one-cosine form: no odd harmonics, but it acquires k=2 content and is **not** reflection-symmetric: asymmetry/range 5.6e-3 (default), median 2.1e-3, max 1.07e-2 over 300 draws (noise floor 1e-6).
     Consequence for Thm "nogo": the stationary set is **not frozen** at {0,π/4,π/2}. Default tilted case: minima at θ = 62.24°, 152.09° (2×2 block-alignment at 62.17°, 152.17°), maxima at 13.89°, 110.48° (a frozen π/4 point would be at 17.17°, 107.17°: shifted by ~3.3°), minima values 4.4e-4 and 9.8e-4 (**not 0**, and unequal: the two aligned points are no longer equivalent).
     Why: the e-row of H_e (a vector u in E) enters ‖[He,Hν]‖² through u†(Hν|_E − m₁²)²u, a second orientation (u's direction) besides ψ₀ (the eigenframe of H_e|_E). Two relative angles, so one cos cannot carry it. Size scales with the e-row coupling (here s13-size, ~0.15).
   - Why Tr(HeHν) survives: H_ν(θ) has zero e-row, so H_e's frequency-1 e-row part multiplies zero. The brief's "e-row rotates at frequency 1" does not reach the trace functionals when E is an eigenblock of H_ν.
2. **E is not an eigenblock of H_ν (both matrices tilted against the fixed plane):** Thm "single" **fails**: odd harmonics appear — odd/range 0.051 (Tr(HeHν)), 0.25 (‖[He,Hν]‖²), 0.077 (misalignment) for H_ν tilt (0,5°,0); 0.093/0.33/0.14 for (3°,8.6°,0), CP 197°; reflection asymmetry 0.03–0.15. But here R(θ) of the plane is no longer a symmetry of H_ν either (Prop 1's premise is gone), so this is outside the paper's stated scope; it is what a model with E = flavor μτ plane and tilted H_ν would face.
Limits: two default parameter sets + 300 random draws in a narrow realistic window; "illustrative" spectra; single N; the reflection test is exact for the sampled series but uses harmonics ≤12.

## (ii) Tilted degenerate eigenspace — `analysis/tilted_degenerate_eigenspace-run-2026-10-23.txt`
Inputs *illustrative*: (θ12,θ13,θ23) = (33.4°, 8.6°, 49°); δ = 0 and 197°; m² = 2.5e-3 eV², Δ = m₁²−m² = −7.4e-5 eV².
Ground truth passes: θ12=θ13=0 ⇒ flavor 2-3 block = m²·1 (deviation 4.9e-20) and U→U(1⊕R(t)) moves θ23 by exactly t (49→39, 24, 9 at t = 10, 25, 40°) with θ12=θ13=0.

**Premise "E = μ–τ flavor plane": fails** for θ12,θ13 ≠ 0. The flavor 2-3 block of Hν=U diag(m₁²,m²,m²)U† has eigenvalues 2.476423e-03, 2.500000e-03 (gap 2.358e-05 eV² = |Δ|(1−|U_e1|²), 1−|U_e1|²=0.3186): not ∝ 𝟙, i.e. a "solar-size" gap of 32% of Δ. The true degenerate eigenspace (columns 2,3 of U; Hν E = m²E to 5e-20) is tilted out of the flavor plane.

**Rotating inside the true E:** Hν is unchanged (max|ΔHν| ≤ 4.3e-19) and |U_e1|,|U_μ1|,|U_τ1| are invariant (drift 0), and P(μ→τ) at an arbitrary phase is identical to 12 digits for t = 0°, 17°, 63° (0.001718686829 for δ=0). But **θ12 and θ13 do change**, so "R shifts θ23 and leaves θ12, θ13" is **false** off the θ12=θ13=0 point:
δ=0: t=0→5→10→20°: θ13 = 8.600→5.827→3.024→2.616°, θ12 = 33.400→33.927→34.248→34.277°, θ23 = 49.0→44.806→40.653→32.391°; θ13 reaches 0 (0.0068° on grid) at t=15.35°, with θ12=34.365°, θ23=36.232° — an exactly equivalent description of the same Hν with θ13=0.
δ=197°: θ13 ranges 2.423°–34.255°, never 0 (complex phases), t=20°: (29.28°, 18.85°, 31.44°) vs (33.40°, 8.60°, 49.00°). Invariant: |U_e2|²+|U_e3|²=0.318613886530 for all t.
(δ read from |U_ij| only; real-matrix cases show |δ|=0 or 180° by branch — not meaningful there.)

**Does the obstruction survive, and what does it mean:** the orbit statement survives (a 1-parameter family of PMNS matrices gives the same Hν and the same oscillation observables), but it is **not "θ23 undetermined"**: with m₂=m₃ the physical content is only the ν₁ direction (|U_e1|,|U_μ1|,|U_τ1|), and θ12, θ13, θ23 are *all* frame labels along the orbit — θ13 is as unobservable as θ23, and any given θ13 can be set to any value in the range above by choice of t. So **Δm²_atm = 0 makes the mixing angles unobservable as a set**, not θ23 alone. Lifting the degeneracy (m₃² = m²+ε) makes t physical at once: max|ΔHν(t=20°)| = 0.2908 ε (δ=0), 0.3132 ε (δ=197°), exactly linear over ε = 2.5e-3…2.5e-6 eV². So the data (ε ≠ 0) fix t; the paper's obstruction is the ε→0 limit, and "external input" is, in this regime, the orientation of the ε-term on E (consistent with its Cor. "external", condition (b), but then the object to be fixed is the ε-term's eigenframe, not θ23 in isolation).
Limits: illustrative angles; one δ=0 and one δ=197° sample; Δ chosen solar-size (the brief's tension with Δm²₃₁ is item (iii), not run here).

## (iv) Hermitian composites — `analysis/composite_operators-run-2026-10-23.txt`
2×2 on E; a=δe/2, b=δν/2, p=c_e c_ν, q=ab, φ=2(θ+ψ₀). Four mass sets (physical-like p/q=41.3 and 1.26, PSD-extreme p/q=1, indefinite toy p/q=0.4).
- **Q₁=½{He,Hν}: does NOT evade** the linear class. Identity Q₁ = c_ν X + c_e Y − c_e c_ν𝟙 + ab cosφ 𝟙 verified to 1.8e-15; the traceless part has θ-independent coefficients (α,β)=(c_ν,c_e) (fit spread 0, residual 1.7e-15); θ enters only through a scalar (inert). Its eigenframe self-consistency ([Q₁,Hν(θ)]=0) holds only at φ=180°/360° (aligned) in all four cases. **Reduces to alignment {0,π/2}.**
- **Q₂=HeHνHe: evades** the fixed-coefficient linear class — the fit coefficient α(θ) of H_e⁰ moves over [3.169, 3.326] (physical-like) and [0.652, 5.685] (wide split), residual 5e-15; Bloch identity V = (c_e²−a²)b m + 2a(c_e c_ν+ab cosφ) n verified to 4.4e-15. **But** its eigenframe still reduces to alignment for physical (PSD) Hermitians: [Q₂,Hν(θ)]=0 needs sinφ=0 or cosφ = −p/q; PSD ⇒ c ≥ δ/2 ⇒ p/q ≥ 1 (100000 random PSD draws: min p/q = 1.0069; equality only if both lowest eigenvalues are 0), so no interior zero. Only for indefinite Hν (not mm†) do interior zeros appear (toy, p/q=0.4: extra zeros at θ=36.735°, 103.158° vs aligned 69.945°).
- **Frozen-stationary-set claim does not carry over to Q₂:** ‖[Q₂,Hν(θ)]‖² has an interior **maximum** at cosφ* = (−p+√(p²+8q²))/(4q), matching the scan: φ* = 88.61° (θ=24.25°, p/q=41.3), 62.66° (11.28°, p/q=1.26), 60.00° (9.95°, p/q=1). It is a mass-dependent interior stationary point, approaching the π/4 value (φ=90°) only for p/q ≫ 1 (the physical-like case). Minima remain at alignment, so *minimising* still returns alignment {0,π/2}; the paper's "frozen {0,π/4,π/2}" is a property of the Tr/‖[He,Hν]‖ family, not of every Hermitian built from He,Hν.
Limits: 2×2 restriction only (full 3×3 composites not run); ‖[Q,He]‖² not run (analytically ∝ sin²φ, unverified here).

## (v) Independent re-derivation — `analysis/bloch_rederivation-run-2026-10-23.txt`
5000 random draws (all quadrants, independent trace parts), matrices built from eigen-projectors, not the paper's σ-decomposition. All PASS:
V1 eigenvector angle ξ−t 2.2e-16; V2 Tr(HeHν)=2c_ec_ν+½δeδν cos2(t+ψ₀) 1.8e-14; V3a [a·σ,b·σ]=2i(a×b)·σ 1.8e-15; V3b ‖[He,Hν]‖²=½δe²δν² sin²2(t+ψ₀) 3.1e-13; V4 dE/dt = 2C sin4(t+ψ₀) 2.6e-11 (relative to C); V5 Cayley–Hamilton word reductions (Tr X²Y, Tr XYXY = Tr(XY)²−2detX detY, Tr X²Y²) ≤1.7e-13. **Bloch algebra and the closed form hold.**
**Hidden assumption (beyond the paper):** the proof places both Bloch vectors in the xz-plane (real eigenbases). For a complex Hermitian H_e|_E (σ_y component, r = length of its xz projection) the identities become Tr = 2c_ec_ν + ½δeδν r cos2(t+ψ₀′) and E = ½δe²δν²(1 − r²cos²2(t+ψ₀′)) (verified 1.1e-14 / 2.0e-13, 2000 draws): the stationary set is the same family but the minimum is ½δe²δν²(1−r²) > 0 — exact alignment is unreachable (median 1−r² = 0.49 over the draws). Thm "monotone"/"nogo" statements "minima E=0 at exact alignment" need the real-eigenbasis assumption stated.

---
# Addendum (item #3): (ii-A) which pair is near-degenerate in the data; (iii) quasi-degeneracy vs Δm²₃ℓ and Σm_ν

Capture for every number below: `analysis/nufit_degeneracy-run-2026-10-23.txt` (script `analysis/nufit_degeneracy.py`; input `analysis/data/nufit6_central.csv`).
Inputs: NuFIT 6.0 Table 1 **central values only**, arXiv:2410.05380 v2, JHEP 12 (2024) 216, transcribed from `files/R1-literature.md` §3.1 (not re-read from the PDF by me). Variants named as in that paper: **IC19 w/o SK-atm** and **IC24 with SK-atm**. Cosmology bounds (ΛCDM, 95 %) as recorded in R1-literature §3.2 at abstract level: Planck 2018+BAO 0.12 eV (arXiv:1807.06209), DESI DR1+CMB 0.072 eV (arXiv:2404.03002), DESI DR2+CMB 0.064 eV (arXiv:2503.14738). No uncertainties propagated.

## (ii-A) The data's near-degenerate pair is the solar pair, not the 2–3 pair
- Δm²₂₁/|Δm²₃ℓ| = 0.02956 (IC19 w/o SK-atm, NO), 0.02984 (IC19 w/o SK-atm, IO), 0.02981 (IC24 with SK-atm, NO), 0.03015 (IC24 with SK-atm, IO). In both orderings the small splitting is Δm²₂₁, between mass states 1,2 (NO: m₁,m₂ light; IO: m₁,m₂ heavy, m₃ lightest). The 2–3 pair is split by ~Δm²₃ℓ in NO and is the far-apart pair in IO.
- U→U R₁₂(t) (rotation inside the solar pair), all four variants: the third column |U_x3| is unchanged (drift 0.00e+00 at t=10…45°), |U_e3| = 0.1482 / 0.1492 / 0.1488 / 0.1494 (IC19 NO / IC19 IO / IC24 NO / IC24 IO) constant, and θ13, θ23 read back unchanged (e.g. IC19 w/o SK-atm NO: θ13 = 8.520°, θ23 = 48.500° for all t). |U_e1|, |U_e2| and **θ12 change** (IC19 w/o SK-atm NO: |U_e1| 0.8230→0.9609 and |U_e2| 0.5484→0.2339 at t=20°; θ12 33.680°→13.680°).
- U→U R₂₃(t): the first column is unchanged (drift 0.00e+00), |U_e2| and |U_e3| change (IC19 w/o SK-atm NO, t=20°: |U_e2| 0.5484→0.4648, |U_e3| 0.1482→0.3267) and θ12, θ13, θ23 all move (θ13 8.520°→19.068°, θ23 48.500°→30.974°) — the (ii) result again.
- **Which angle:** a rotation inside the data's actually near-degenerate (solar) pair shifts **θ12** (and |U_e1|,|U_e2|), leaves θ13, θ23 and ν₃ fixed. So if the degeneracy that mattered were the data's own, the unobservable frame angle would be θ12, not θ23. The paper's obstruction on θ23 needs the 2–3 pair degenerate, which is not what the data show.
- Not shown: the capture reads angles from |U_ij|, so δ is not tracked; only central values, no ±1σ.

## (iii) What the paper's m₂≈m₃ premise requires
Ground truth (m_lightest = 0): Σm_ν = 0.05899 (IC19 w/o SK-atm, NO), 0.05878 (IC24 with SK-atm, NO) — matches the known ~0.059 eV; IO 0.09945 (IC19 w/o SK-atm), 0.09892 (IC24 with SK-atm) — matches ~0.10 eV. Ground truth PASS.
Relative splitting of the 2–3 pair at m_lightest = 0 (NO): 82.8 % (IC19 w/o SK-atm), 82.7 % (IC24 with SK-atm).

Exact condition (m₃−m₂)/m₃ < f, NO, solved by bisection (m_lightest = m₁):
| f | variant | m_lightest ≥ | Σm_ν | vs 0.12 / 0.072 / 0.064 eV |
|---|---|---|---|---|
| 1 % | IC19 w/o SK-atm, NO | 0.34791 eV | 1.04745 eV | all EXCLUDED |
| 5 % | IC19 w/o SK-atm, NO | 0.15062 eV | 0.46031 eV | all EXCLUDED |
| 10 % | IC19 w/o SK-atm, NO | 0.10202 eV | 0.31818 eV | all EXCLUDED |
| 1 % | IC24 with SK-atm, NO | 0.34642 eV | 1.04297 eV | all EXCLUDED |
| 5 % | IC24 with SK-atm, NO | 0.14998 eV | 0.45834 eV | all EXCLUDED |
| 10 % | IC24 with SK-atm, NO | 0.10158 eV | 0.31681 eV | all EXCLUDED |

IO, paper's literal "2–3" pair (m₂−m₃)/m₂ < f: relative splitting 100 % at m_lightest=0; thresholds m₃ ≥ 0.35160 / 0.15243 / 0.10344 eV (IC19 w/o SK-atm; f = 1/5/10 %), Σm_ν = 1.06179 / 0.47309 / 0.33299 eV; IC24 with SK-atm: 0.34977 / 0.15163 / 0.10291 eV, Σ = 1.05627 / 0.47063 / 0.33126 eV — all EXCLUDED by all three bounds.
IO, the data's near-degenerate pair (1,2), (m₂−m₁)/m₂ < f: already 1.503 % (IC19 w/o SK-atm) / 1.519 % (IC24 with SK-atm) at m_lightest=0, so f = 5 %, 10 % hold with Σ = 0.09945 / 0.09892 eV (allowed by Planck18+BAO 0.12, EXCLUDED by DESI DR1 0.072 and DR2 0.064); f = 1 % needs m₃ ≥ 0.03541 / 0.03577 eV, Σ = 0.15750 / 0.15786 eV (excluded by all three).

**Plain statement.** The paper's premise Hν|_E ∝ 𝟙 with E = the 2–3 pair requires m₂≈m₃: even at a loose 10 % splitting that is m_lightest ≳ 0.10 eV and Σm_ν ≳ 0.317 eV (IC24 with SK-atm, NO) – 0.318 eV (IC19 w/o SK-atm, NO), 2.6× Planck18+BAO's 0.12 eV and ~5× DESI DR2+CMB's 0.064 eV. **No ΛCDM bound listed allows it**, in either ordering, either variant, at 1, 5 or 10 %. The only pair that is near-degenerate in the data is the solar pair (Δm²₂₁/|Δm²₃ℓ| ≈ 0.03), and only for IO with its absolute scale pinned near the minimum (Σ ≈ 0.099 eV, itself excluded by the DESI ΛCDM bounds). One-line caveat: non-ΛCDM extensions relax these bounds; values not evaluated here. This tension is a statement about the premise "E = 2–3 degenerate pair"; the (i),(ii),(iv) results above (what does and does not follow once such an E is granted) are unaffected.
Limits: central values only; bounds taken from R1-literature at abstract level (dataset details for Planck's likelihood combination and DR2's CMB set unconfirmed there); the threshold f values 1/5/10 % are the manager's choice, not a physical cutoff.
