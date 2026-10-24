#!/usr/bin/env python3
"""
bloch_rederivation.py -- audit (v): independent numerical re-derivation of the Bloch-form
algebra and of E(theta) = (1/2) de^2 dnu^2 sin^2 2(theta+psi0) in theta23_obstruction_paper.tex
(Thm "single", Thm "nogo"). sympy is NOT installed, so every identity is checked numerically
to 1e-12 on random draws (sympy~=1.14 is listed in requirements.txt but absent in the env).

Independence from the paper's verifier: matrices are built from eigen-projectors
(flavor_common.bloch2), not from the sigma-decomposition used in verify_theta23_obstruction.py.
Checks, each over many random parameter draws (angles over all quadrants, gaps of both signs of
interest, trace parts independent):
  V1 eigenvector angle of R(t)^dag H R(t) is xi - t                      (paper: psi0 -> psi0 + theta)
  V2 Tr(He Hnu(t)) = 2 c_e c_nu + (1/2) de dnu cos 2(t + psi0)
  V3 [a.sigma, b.sigma] = 2i (a x b).sigma  and  ||[He,Hnu]||_F^2 = (1/2) de^2 dnu^2 sin^2 2(t+psi0)
  V4 E'(t) = 2 C sin 4(t+psi0), C = (1/2) de^2 dnu^2   (analytic derivative vs central difference)
  V5 Cayley-Hamilton word reductions behind "five generators suffice":
       Tr(X^2 Y) = tX Tr(XY) - dX tY ;  Tr(XYXY) = Tr(XY)^2 - 2 dX dY ;  Tr(X^2Y^2) = tX Tr(XY^2) - dX Tr(Y^2)
     with tX = Tr X, dX = det X
  V6 (beyond the paper) He|_E COMPLEX Hermitian (Bloch vector with a sigma_y component, r = sin chi < 1
     of its xz-projection): Tr(He Hnu(t)) = 2 c_e c_nu + (1/2) de dnu r cos 2(t+psi0') and
     E(t) = (1/2) de^2 dnu^2 (1 - r^2 cos^2 2(t+psi0')). Reports whether the minimum of E is still 0.
Run: python3 analysis/bloch_rederivation.py
Made by: crew flavor-obstruction-paper, item #1 (R1), data-analyst.
"""
import numpy as np
from flavor_common import bloch2, rot2

SY = np.array([[0, -1j], [1j, 0]])
SXm = np.array([[0, 1], [1, 0]], complex)
SZm = np.array([[1, 0], [0, -1]], complex)


def fro2(M):
    return float(np.sum(np.abs(M)**2))


def main():
    rng = np.random.default_rng(1023)
    n = 5000
    worst = dict(V1=0, V2=0, V3a=0, V3b=0, V4=0, V5a=0, V5b=0, V5c=0)
    for _ in range(n):
        ce, cn = rng.uniform(-3, 3, 2)
        de, dn = rng.uniform(0.05, 5, 2)
        xe, xn, t = rng.uniform(-np.pi, np.pi, 3)
        X, Y0 = bloch2(de, xe, ce), bloch2(dn, xn, cn)
        R = rot2(t)
        Y = R.conj().T @ Y0 @ R
        psi0 = xe - xn
        # V1: larger-eigenvalue eigenvector direction of Y is (cos(xn - t), sin(xn - t)) up to sign
        w, V = np.linalg.eigh(Y)
        v = V[:, 1]
        target = np.array([np.cos(xn - t), np.sin(xn - t)])
        worst['V1'] = max(worst['V1'], 1 - abs(np.real(np.vdot(target, v))))
        # V2
        worst['V2'] = max(worst['V2'], abs(np.real(np.trace(X @ Y)) - (2 * ce * cn + 0.5 * de * dn * np.cos(2 * (t + psi0)))))
        # V3
        a, b = rng.normal(size=3), rng.normal(size=3)
        sa = a[0] * SXm + a[1] * SY + a[2] * SZm
        sb = b[0] * SXm + b[1] * SY + b[2] * SZm
        lhs = sa @ sb - sb @ sa
        rhs = 2j * (np.cross(a, b)[0] * SXm + np.cross(a, b)[1] * SY + np.cross(a, b)[2] * SZm)
        worst['V3a'] = max(worst['V3a'], np.max(np.abs(lhs - rhs)))
        Ecomm = fro2(X @ Y - Y @ X)
        worst['V3b'] = max(worst['V3b'], abs(Ecomm - 0.5 * de**2 * dn**2 * np.sin(2 * (t + psi0))**2))
        # V4: analytic derivative vs central difference of the numerically built E
        h = 1e-3

        def Enum(tt):
            Rr = rot2(tt)
            Yy = Rr.conj().T @ Y0 @ Rr
            return fro2(X @ Yy - Yy @ X)
        d_num = (-Enum(t + 2*h) + 8*Enum(t + h) - 8*Enum(t - h) + Enum(t - 2*h)) / (12 * h)   # 5-point stencil
        d_an = 2 * (0.5 * de**2 * dn**2) * np.sin(4 * (t + psi0))
        worst['V4'] = max(worst['V4'], abs(d_num - d_an) / (0.5 * de**2 * dn**2))
        # V5
        tX, tY = np.trace(X), np.trace(Y)
        dX, dY = np.linalg.det(X), np.linalg.det(Y)
        tXY = np.trace(X @ Y)
        worst['V5a'] = max(worst['V5a'], abs(np.trace(X @ X @ Y) - (tX * tXY - dX * tY)))
        worst['V5b'] = max(worst['V5b'], abs(np.trace(X @ Y @ X @ Y) - (tXY**2 - 2 * dX * dY)))
        worst['V5c'] = max(worst['V5c'], abs(np.trace(X @ X @ Y @ Y) - (tX * np.trace(X @ Y @ Y) - dX * np.trace(Y @ Y))))
    print(f"{n} random draws (seed 1023). Worst absolute (V1: 1-|overlap|; V4: relative to C) deviations:")
    labels = dict(V1="V1 eigenvector angle xi-t", V2="V2 Tr(He Hnu(t)) closed form", V3a="V3a [a.s,b.s]=2i(axb).s",
                  V3b="V3b ||[He,Hnu]||^2 closed form", V4="V4 dE/dt = 2C sin4(t+psi0)",
                  V5a="V5a Tr(X^2Y) reduction", V5b="V5b Tr(XYXY) reduction", V5c="V5c Tr(X^2Y^2) reduction")
    allok = True
    for k, v in worst.items():
        ok = v < 1e-8 if k == 'V4' else v < 1e-12
        allok &= ok
        print(f"  {labels[k]:36s} {v:.3e}   {'PASS' if ok else 'FAIL'}")
    print("  (V4 tolerance 1e-8 (relative to C): 5-point stencil, h=1e-3; the rest 1e-12)")

    print("\nV6 (beyond the paper): He|_E complex Hermitian, Bloch vector n = (sin chi cos s, cos chi, sin chi sin s)-type with a sigma_y part")
    worst2 = dict(over=0, E=0)
    minE_ratio = []
    for _ in range(2000):
        ce, cn = rng.uniform(0, 3, 2)
        de, dn = rng.uniform(0.05, 5, 2)
        xe, xn, t = rng.uniform(-np.pi, np.pi, 3)
        chi = rng.uniform(0.05, np.pi / 2)             # sin chi = r: length of the xz-projection
        r = np.sin(chi)
        # He = ce + (de/2)[ r (cos2xe sz + sin2xe sx) + cos(chi) sy ]
        X = ce * np.eye(2) + (de / 2) * (r * (np.cos(2 * xe) * SZm + np.sin(2 * xe) * SXm) + np.cos(chi) * SY)
        Y0 = bloch2(dn, xn, cn)
        R = rot2(t)
        Y = R.conj().T @ Y0 @ R
        psi0 = xe - xn
        worst2['over'] = max(worst2['over'], abs(np.real(np.trace(X @ Y)) - (2 * ce * cn + 0.5 * de * dn * r * np.cos(2 * (t + psi0)))))
        worst2['E'] = max(worst2['E'], abs(fro2(X @ Y - Y @ X) - 0.5 * de**2 * dn**2 * (1 - r**2 * np.cos(2 * (t + psi0))**2)))
        minE_ratio.append(1 - r**2)
    ok6 = worst2['over'] < 1e-12 and worst2['E'] < 1e-12
    allok &= ok6
    print(f"  max|Tr(He Hnu) - [2 ce cnu + (1/2) de dnu r cos2(t+psi0)]| = {worst2['over']:.2e}")
    print(f"  max|E - (1/2)de^2 dnu^2 (1 - r^2 cos^2 2(t+psi0))|         = {worst2['E']:.2e}   -> {'PASS' if ok6 else 'FAIL'}")
    print(f"  consequence: E_min = (1/2)de^2 dnu^2 (1-r^2) > 0 unless r=1: over these draws min(1-r^2) = {min(minE_ratio):.3e}, median = {np.median(minE_ratio):.3e}")
    print("  stationary set of E(t): cos2(t+psi0)*sin2(t+psi0)=0 -> same {0, pi/4, pi/2} family (psi0 is the xz-projection angle)")
    print("\nALL CHECKS PASS" if allok else "\nSOME CHECK FAILED")
    return 0 if allok else 1


if __name__ == "__main__":
    raise SystemExit(main())
