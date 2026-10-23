#!/usr/bin/env python3
"""
d3_rstar_identity.py -- D3 check: the aligning label ratio is a column-orthogonality identity.
For real PDG U (delta = 0, 180) and every NuFIT 6.0 case in data/nufit6_central.csv:
  r* = -U_mu2 U_mu3 / (U_tau2 U_tau3)  ==  1 + U_e2 U_e3 / (U_tau2 U_tau3)   (exact),
  leading order 1 -/+ 2 tan(th12) sin(th13) / sin(2 th23);  r* = 1 when th13 = 0.
Ground truth: the lab's PART 2 values (calculations/theta23-harmonic-label/run-2026-10-23-t3.txt,
lab 3f67f4b), quoted below to 10 digits.
Made by: crew flavor-obstruction-paper lead, R1, 2026-10-23.
"""
import csv, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = {("IC19 w/o SK-atm", "NO", 0): 0.8170349147, ("IC19 w/o SK-atm", "NO", 180): 1.2179887361,
       ("IC24 with SK-atm", "NO", 0): 0.8201923095, ("IC24 with SK-atm", "NO", 180): 1.2221140769,
       ("IC19 w/o SK-atm", "IO", 0): 0.8157270425, ("IC19 w/o SK-atm", "IO", 180): 1.2196804956,
       ("IC24 with SK-atm", "IO", 0): 0.8164127100, ("IC24 with SK-atm", "IO", 180): 1.2198741390}


def pdg(t12, t13, t23, d):
    t12, t13, t23, d = np.radians([t12, t13, t23, d])
    s12, c12, s13, c13, s23, c23 = np.sin(t12), np.cos(t12), np.sin(t13), np.cos(t13), np.sin(t23), np.cos(t23)
    e = np.exp(1j * d)
    return np.array([[c12 * c13, s12 * c13, s13 / e],
                     [-s12 * c23 - c12 * s23 * s13 * e, c12 * c23 - s12 * s23 * s13 * e, s23 * c13],
                     [s12 * s23 - c12 * c23 * s13 * e, -c12 * s23 - s12 * c23 * s13 * e, c23 * c13]]).real


def main():
    rows = [r for r in open(os.path.join(HERE, "data", "nufit6_central.csv")) if not r.startswith("#")]
    ok = True
    for r in csv.DictReader(rows):
        for d in (0, 180):
            u = pdg(float(r["theta12"]), float(r["theta13"]), float(r["theta23"]), d)
            rs = -u[1, 1] * u[1, 2] / (u[2, 1] * u[2, 2])
            ri = 1 + u[0, 1] * u[0, 2] / (u[2, 1] * u[2, 2])
            lo = 1 - np.cos(np.radians(d)) * 2 * np.tan(np.radians(float(r["theta12"]))) * \
                np.sin(np.radians(float(r["theta13"]))) / np.sin(np.radians(2 * float(r["theta23"])))
            lab = LAB[(r["variant"], r["ordering"], d)]
            good = abs(rs - ri) < 1e-12 and abs(rs - lab) < 1e-9
            ok &= good
            print(f"  [{r['variant']} {r['ordering']} d={d:3d}] r*={rs:.10f} identity={ri:.10f} "
                  f"lab={lab:.10f} leading={lo:.4f} -> {'PASS' if good else 'FAIL'}")
    u0 = pdg(33.68, 0.0, 48.5, 0)
    r0 = -u0[1, 1] * u0[1, 2] / (u0[2, 1] * u0[2, 2])
    print(f"  theta13 = 0: r* = {r0:.15f} -> {'PASS' if abs(r0 - 1) < 1e-12 else 'FAIL'}")
    ok &= abs(r0 - 1) < 1e-12
    print("\n" + ("D3 CHECKS PASS" if ok else "D3 CHECK FAILED"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
