#!/usr/bin/env python3
"""
verify_theta23_obstruction.py

Reproducibility appendix for
  "Spectral geometry does not fix the atmospheric mixing angle:
   a residual-symmetry obstruction in the noncommutative Standard Model"

Checks, to machine precision, every theorem in the paper:

  Prop 1  (obstruction)      : degenerate H_nu  =>  E(theta) == 0 identically
  Thm  1  (single invariant) : tr(He Hnu(t)) = 2 ce cnu + (1/2) de dnu cos 2(t+psi0),
                               and the other 4 joint invariants are t-independent
  Thm  2  (monotone no-go)   : every monotone relational functional selects
                               exact alignment theta23 in {0, 90} deg, for all masses
  Thm  3  (commutator norm)  : E(t) = (1/2) de^2 dnu^2 sin^2 2(t+psi0);
                               stationary set {0,45,90} deg; minimizer = 0 (mass-indep)

Deterministic. NumPy only. No network.

Run:  python3 verify_theta23_obstruction.py
"""

import numpy as np

# --- Pauli matrices and a 2D real rotation on the mu-tau subspace E ---------
SX = np.array([[0, 1], [1, 0]], complex)
SZ = np.array([[1, 0], [0, -1]], complex)


def rot(t):
    return np.array([[np.cos(t), -np.sin(t)],
                     [np.sin(t),  np.cos(t)]], complex)


def herm(gap, xi, c):
    """2x2 Hermitian with eigenvalue gap `gap`, eigenbasis angle `xi`, trace part `c`."""
    n = np.cos(2 * xi) * SZ + np.sin(2 * xi) * SX      # Bloch vector at angle 2*xi
    return c * np.eye(2) + (gap / 2.0) * n


def commutator_norm2(He, Hnu_t):
    C = He @ Hnu_t - Hnu_t @ He
    return float(np.sum(np.abs(C) ** 2))               # ||.||_F^2


def theta23_of(psi0, t):
    """Atmospheric angle in degrees (mod 90) at rotation parameter t."""
    return np.degrees(psi0 + t) % 90.0


def main():
    ok = True
    TH = np.linspace(-np.pi / 2, np.pi / 2, 200001)

    # reference sector operators (arbitrary O(1) choices)
    de, dnu, xi_e, xi_nu, c_e, c_nu = 3.0, 1.7, 0.35, 0.9, 1.0, 0.5
    He = herm(de, xi_e, c_e)
    Hnu0 = herm(dnu, xi_nu, c_nu)
    psi0 = xi_e - xi_nu

    def Hnu(t):
        return rot(t).conj().T @ Hnu0 @ rot(t)

    print("=" * 70)
    print("Prop 1 (obstruction): degenerate H_nu -> E(theta) identically 0")
    print("=" * 70)
    Hnu_deg = 2.5 * np.eye(2)
    emax = max(commutator_norm2(He, Hnu_deg) for t in [0, .3, .7, 1.1, 1.5])
    print(f"  max E over sampled theta = {emax:.2e}   -> {'PASS' if emax < 1e-20 else 'FAIL'}")
    ok &= emax < 1e-20

    print("\n" + "=" * 70)
    print("Thm 1 (single invariant): residual freedom lives in tr(He Hnu(theta))")
    print("=" * 70)
    overlap = np.array([np.real(np.trace(He @ Hnu(t))) for t in TH])
    pred = 2 * c_e * c_nu + 0.5 * de * dnu * np.cos(2 * (psi0 + TH))
    dev = np.max(np.abs(overlap - pred))
    print(f"  max|tr(He Hnu(t)) - [2 ce cnu + (1/2) de dnu cos2(t+psi0)]| = {dev:.2e}"
          f"   -> {'PASS' if dev < 1e-10 else 'FAIL'}")
    ok &= dev < 1e-10
    for nm, fn in [("tr He",       lambda t: np.trace(He)),
                   ("tr Hnu(t)",   lambda t: np.trace(Hnu(t))),
                   ("tr He^2",     lambda t: np.trace(He @ He)),
                   ("tr Hnu(t)^2", lambda t: np.trace(Hnu(t) @ Hnu(t)))]:
        vals = np.array([np.real(fn(t)) for t in TH])
        rng = vals.max() - vals.min()
        print(f"    {nm:12s} range over theta = {rng:.2e}   ({'t-indep' if rng < 1e-10 else 'VARIES'})")
        ok &= rng < 1e-10

    print("\n" + "=" * 70)
    print("Thm 2 (monotone no-go): monotone functionals -> exact alignment, all masses")
    print("=" * 70)
    monos = [("G=-overlap        ", lambda o: -o),
             ("G=exp(-overlap)   ", lambda o: np.exp(-o)),
             ("G=(overlap-5)^3   ", lambda o: (o - 5.0) ** 3),
             ("G=1/(overlap+10)  ", lambda o: 1.0 / (o + 10.0))]
    for name, G in monos:
        F = G(overlap)
        t23 = theta23_of(psi0, TH[np.argmin(F)])
        aligned = t23 < 0.5 or t23 > 89.5
        print(f"  {name} min theta23 = {t23:6.2f} deg   -> {'ALIGNED' if aligned else 'INTERIOR'}")
        ok &= aligned
    print("  vary the masses (de,dnu,xi_e,xi_nu):")
    for (d1, d2, x1, x2) in [(1, 1, 0.1, 0.8), (10, 0.2, 1.2, 0.3), (0.3, 5, 0.5, 1.4)]:
        A = herm(d1, x1, 2.0)
        B0 = herm(d2, x2, 1.0)
        p = x1 - x2
        ov = np.array([np.real(np.trace(A @ (rot(t).conj().T @ B0 @ rot(t)))) for t in TH])
        t23 = theta23_of(p, TH[np.argmin(-ov)])
        aligned = t23 < 0.5 or t23 > 89.5
        print(f"    ({d1},{d2},{x1},{x2}) -> min theta23 = {t23:6.2f} deg   -> {'ALIGNED' if aligned else 'INTERIOR'}")
        ok &= aligned

    print("\n" + "=" * 70)
    print("Thm 3 (commutator norm): E(t)=(1/2)de^2 dnu^2 sin^2 2(t+psi0)")
    print("=" * 70)
    E = np.array([commutator_norm2(He, Hnu(t)) for t in TH])
    Epred = 0.5 * de ** 2 * dnu ** 2 * np.sin(2 * (psi0 + TH)) ** 2
    dev = np.max(np.abs(E - Epred))
    print(f"  max|E_numeric - closed form| = {dev:.2e}   -> {'PASS' if dev < 1e-10 else 'FAIL'}")
    ok &= dev < 1e-10
    t23_min = theta23_of(psi0, TH[np.argmin(E)])
    aligned = t23_min < 0.5 or t23_min > 89.5
    print(f"  minimizer theta23 = {t23_min:.2f} deg (mod 90)   -> "
          f"{'ALIGNED (not ~49 deg)' if aligned else 'INTERIOR'}")
    ok &= aligned
    print("  stationary set {0,45,90} deg is mass-independent: C is a pure prefactor (analytic).")

    print("\n" + "=" * 70)
    print("ALL CHECKS PASS" if ok else "SOME CHECK FAILED")
    print("=" * 70)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
