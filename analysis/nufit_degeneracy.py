#!/usr/bin/env python3
"""
nufit_degeneracy.py -- audit (ii-A) and (iii): which neutrino pair is near-degenerate in the data,
what a rotation inside it shifts, and what m2~m3 requires of m_lightest and sum m_nu.

Calculates, for each NuFIT 6.0 variant/ordering in data/nufit6_central.csv (central values only):
  A. dm21/|dm3l|; |U_e1|,|U_e2|,|U_e3| and the third column under U -> U R12(t) and U -> U R23(t);
     PDG angles read back (flavor_common.pdg_angles) -- which angle moves.
  B. ground truth: minimal sum m_nu at m_lightest = 0 (NO ~0.059, IO ~0.10 eV);
     m_lightest and sum m_nu at which the near-degenerate pair has relative splitting < f = 1,5,10 %:
       NO: (m3-m2)/m3 < f   (the paper's m2~m3);  m_l = m1, m2^2 = m_l^2+dm21, m3^2 = m_l^2+dm31.
       IO: m3 lightest, m2^2 = m3^2+|dm32|, m1^2 = m2^2-dm21.
         pair (1,2): (m2-m1)/m2 < f   (the data's near-degenerate pair)
         pair (2,3): (m2-m3)/m2 < f   (the paper's literal "2-3" labels in IO)
     Solved by bisection on the exact condition (relative splitting is monotone decreasing in m_l).
  Compared with Sigma m_nu bounds (LCDM, 95% CL): Planck 2018+BAO 0.12 (arXiv:1807.06209),
  DESI DR1+CMB 0.072 (arXiv:2404.03002), DESI DR2+CMB 0.064 (arXiv:2503.14738) eV -- as recorded in
  files/R1-literature.md section 3.2 (abstract-level; not re-read here).
Inputs: data/nufit6_central.csv (source in its header). No fitted uncertainties are propagated.
Run: python3 nufit_degeneracy.py     Made by: crew flavor-obstruction-paper, item #3 (R1), data-analyst.
"""
import csv
import os
import numpy as np
from flavor_common import pdg_U, pdg_angles, rot2, embed_E

D = np.pi / 180
BOUNDS = [("Planck18+BAO", 0.12), ("DESI DR1+CMB", 0.072), ("DESI DR2+CMB", 0.064)]
FRACS = [0.01, 0.05, 0.10]


def load():
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "nufit6_central.csv")
    with open(path) as fh:
        rows = list(csv.DictReader(r for r in fh if not r.startswith("#")))
    return [dict(variant=r["variant"], ordering=r["ordering"], t12=float(r["theta12"]), t13=float(r["theta13"]),
                 t23=float(r["theta23"]), dl=float(r["delta"]), d21=float(r["dm21_1e-5"]) * 1e-5,
                 d3l=float(r["dm3l_1e-3"]) * 1e-3) for r in rows]


def R12(t):
    W = np.eye(3, dtype=complex); W[:2, :2] = rot2(t); return W


def masses(r, ml):
    if r["ordering"] == "NO":
        return ml, np.sqrt(ml**2 + r["d21"]), np.sqrt(ml**2 + r["d3l"])
    m3 = ml
    m2 = np.sqrt(m3**2 + abs(r["d3l"]))
    return np.sqrt(m2**2 - r["d21"]), m2, m3


def rel_split(r, ml, pair):
    m1, m2, m3 = masses(r, ml)
    if r["ordering"] == "NO":
        return (m3 - m2) / m3
    return (m2 - m1) / m2 if pair == "12" else (m2 - m3) / m2


def solve_ml(r, f, pair):
    if rel_split(r, 0.0 + 1e-18, pair) < f:
        return 0.0
    lo, hi = 1e-12, 10.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if rel_split(r, mid, pair) > f: lo = mid
        else: hi = mid
    return hi


def main():
    rows = load()
    for r in rows:
        tag = f"[{r['variant']}, {r['ordering']}]"
        print(f"\n===== {tag}: th12={r['t12']} th13={r['t13']} th23={r['t23']} delta={r['dl']} deg; "
              f"dm21={r['d21']:.4e} dm3l={r['d3l']:+.4e} eV^2")
        print(f"  A. dm21/|dm3l| = {r['d21']/abs(r['d3l']):.5f}  (solar pair is the near-degenerate one in both orderings)")
        U = pdg_U(r['t12']*D, r['t13']*D, r['t23']*D, r['dl']*D)
        a0 = pdg_angles(U)
        print("     U (|entries|) row e: %s ; third column |U_x3|: %s" % (np.round(np.abs(U[0]), 5), np.round(np.abs(U[:, 2]), 5)))
        for nm, W in (("R12", R12), ("R23", lambda t: embed_E(rot2(t)))):
            print(f"     U -> U {nm}(t):  t(deg)  |Ue1|   |Ue2|   |Ue3|   col3 drift  col1 drift   th12     th13     th23")
            c3 = c1 = 0.0
            for t in (0, 10, 20, 30, 45):
                Ut = U @ W(t * D)
                a = np.degrees(pdg_angles(Ut)[:3])
                d3 = np.max(np.abs(np.abs(Ut[:, 2]) - np.abs(U[:, 2]))); d1 = np.max(np.abs(np.abs(Ut[:, 0]) - np.abs(U[:, 0])))
                print(f"        {t:12d}  {abs(Ut[0,0]):.4f}  {abs(Ut[0,1]):.4f}  {abs(Ut[0,2]):.4f}  {d3:.2e}   {d1:.2e}   {a[0]:7.3f}  {a[1]:7.3f}  {a[2]:7.3f}")
    print("\n===== B. absolute mass scale")
    print("  ground truth: minimal sum m_nu at m_lightest = 0")
    gt_ok = True
    for r in rows:
        m = masses(r, 0.0); s = sum(m)
        ok = (0.057 < s < 0.061) if r['ordering'] == "NO" else (0.098 < s < 0.102)
        gt_ok &= ok
        print(f"    [{r['variant']}, {r['ordering']}] masses = {m[0]:.5f}, {m[1]:.5f}, {m[2]:.5f} eV  sum = {s:.5f} eV  ({'ok' if ok else 'OUT OF RANGE'}; expect ~0.059 NO / ~0.10 IO)")
    print("  ground truth:", "PASS" if gt_ok else "FAIL")
    print("  bounds (LCDM 95%): " + ", ".join(f"{n} {v} eV" for n, v in BOUNDS))
    for r in rows:
        for pair in (["23"] if r['ordering'] == "NO" else ["12", "23"]):
            lab = "pair (2,3): (m3-m2)/m3" if r['ordering'] == "NO" else ("pair (1,2): (m2-m1)/m2" if pair == "12" else "pair (2,3): (m2-m3)/m2")
            print(f"  [{r['variant']}, {r['ordering']}] {lab}; relative splitting at m_l=0: {100*rel_split(r, 1e-18, pair):.3f} %")
            for f in FRACS:
                ml = solve_ml(r, f, pair)
                m = masses(r, ml); s = sum(m)
                verdict = ", ".join(f"{n}: {'allowed' if s < v else 'EXCLUDED'}" for n, v in BOUNDS)
                print(f"     f<{100*f:4.0f}%: m_lightest >= {ml:.5f} eV, masses {m[0]:.5f}/{m[1]:.5f}/{m[2]:.5f}, sum = {s:.5f} eV -> {verdict}")


if __name__ == "__main__":
    main()
