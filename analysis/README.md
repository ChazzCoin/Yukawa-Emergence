# analysis/ — theta23-paper audit scripts (crew flavor-obstruction-paper, item #1)

Run from `analysis/` with `/opt/homebrew/bin/python3` (numpy 2.4.4; system `python3` has no numpy). NumPy only.
Each script has a header (inputs, method, run line). Run captures: `<script>-run-<date>.txt` (env stamped, never edited).

| file | what it calculates | run |
|---|---|---|
| `flavor_common.py` | PDG U, angle extraction, Bloch/Hermitian builders, Fourier helper (library; self-check) | `python3 flavor_common.py` |
| `full3x3_fourier.py` | (i) Fourier/one-cosine test of Tr(HeHν), ‖[He,Hν]‖², Tr((HeHν)²), eigenframe misalignment for the full 3×3 pair; stationary set of the commutator functional | `python3 full3x3_fourier.py [seed]` (~1 min) |
| `tilted_degenerate_eigenspace.py` | (ii) E as true degenerate eigenspace with PDG θ12,θ13≠0, m2=m3: what in-E rotation does to angles, Hν, P(μ→τ); lifted-degeneracy scaling | `python3 tilted_degenerate_eigenspace.py` |
| `composite_operators.py` | (iv) ½{He,Hν} and HeHνHe on E: linear-class test, self-consistent eigenframe, stationary points | `python3 composite_operators.py` (~15 s) |
| `bloch_rederivation.py` | (v) numeric re-derivation of Bloch algebra and E closed form (+ complex-H_e case) | `python3 bloch_rederivation.py` |
| (verifier capture) `verify-run-2026-10-23.txt` | gate: output of `../verify_theta23_obstruction.py` | `python3 ../verify_theta23_obstruction.py` |
