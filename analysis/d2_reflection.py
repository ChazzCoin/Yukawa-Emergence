#!/usr/bin/env python3
"""
d2_reflection.py

What it calculates: the "D2: Majorana phases" ground-truth gate for the constructive section,
in the Grimus-Lavoura / Harrison-Scott form. For a complex-symmetric Majorana mass matrix m that
is mu-tau reflection symmetric, S m* S = m (S swaps mu and tau), check
  (1) |U_mu j| = |U_tau j| for every column, hence theta23 = 45 deg exactly;
  (2) delta = +-90 deg (|sin delta| = 1) whenever theta13 != 0, with both signs occurring;
  (3) the CP parities: rho = U^dagger S U* is diagonal with entries +-1 (Takagi U, positive masses);
  (4) control: a generic complex-symmetric m (no reflection) fails (1)-(3).
Inputs: none from data. m is drawn at random (seeded); m_ee and m_mutau real, m_emu = conj(m_etau),
  m_mumu = conj(m_tautau), which is the general reflection-symmetric form.
Method: Takagi factorisation U^T m U = diag(m_i >= 0) from the eigenvectors of m^dagger m, with
  column phases fixed by the diagonal of U^T m U. Angles from |U|; sin delta from the Jarlskog
  invariant J = Im(U_e1 U_mu2 conj(U_e2) conj(U_mu1)).
Run:  python3 analysis/d2_reflection.py
Made by: crew flavor-obstruction-paper lead, R1 (D2 verifier check), 2026-10-23.
Source of the result it gates: Domain.Alignment lab crew, calculations/theta23-majorana-phase/
  (gate a139b08, result 71a1021, confirmed 2010400).
"""
import numpy as np

S = np.array([[1, 0, 0], [0, 0, 1], [0, 1, 0]], dtype=complex)


def takagi(m):
    w, V = np.linalg.eigh(m.conj().T @ m)   # m = U* D U^dagger  =>  m^dagger m = U D^2 U^dagger
    d = np.diag(V.T @ m @ V)
    V = V * np.exp(-0.5j * np.angle(d))
    return V, np.sqrt(np.abs(w))


def angles(U):
    s13 = abs(U[0, 2])
    c13 = np.sqrt(1 - s13 ** 2)
    t12 = np.degrees(np.arctan2(abs(U[0, 1]), abs(U[0, 0])))
    t13 = np.degrees(np.arcsin(min(s13, 1.0)))
    t23 = np.degrees(np.arctan2(abs(U[1, 2]), abs(U[2, 2])))
    J = np.imag(U[0, 0] * U[1, 1] * np.conj(U[0, 1]) * np.conj(U[1, 0]))
    r12, r23 = np.radians(t12), np.radians(t23)
    den = np.sin(r12) * np.cos(r12) * np.sin(r23) * np.cos(r23) * c13 ** 2 * s13
    sind = J / den if den > 1e-12 else np.nan
    return t12, t13, t23, sind


def refl_matrix(rng):
    a, d = rng.normal(), rng.normal()
    b = rng.normal() + 1j * rng.normal()
    c = rng.normal() + 1j * rng.normal()
    return np.array([[a, b, np.conj(b)], [b, c, d], [np.conj(b), d, np.conj(c)]])


def generic_matrix(rng):
    X = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
    return X + X.T


def check(m):
    U, mi = takagi(m)
    diag_ok = np.max(np.abs(U.T @ m @ U - np.diag(mi))) < 1e-10 * max(mi)
    t12, t13, t23, sind = angles(U)
    rows = np.max(np.abs(np.abs(U[1]) - np.abs(U[2])))
    rho = U.conj().T @ S @ U.conj()
    offd = np.max(np.abs(rho - np.diag(np.diag(rho))))
    parity = np.max(np.abs(np.abs(np.real(np.diag(rho))) - 1)) + np.max(np.abs(np.imag(np.diag(rho))))
    return diag_ok, rows, t23, t13, sind, offd, parity, np.round(np.real(np.diag(rho))).astype(int)


def main():
    rng = np.random.default_rng(20261023)
    ok = True
    N = 2000
    worst = dict(rows=0.0, t23=0.0, sind=0.0, offd=0.0, parity=0.0)
    signs, pats, takagi_ok = set(), set(), True
    for _ in range(N):
        d, rows, t23, t13, sind, offd, parity, pat = check(refl_matrix(rng))
        takagi_ok &= d
        worst["rows"] = max(worst["rows"], rows)
        worst["t23"] = max(worst["t23"], abs(t23 - 45.0))
        if t13 > 1e-3:
            worst["sind"] = max(worst["sind"], abs(abs(sind) - 1))
            signs.add(int(np.sign(sind)))
        worst["offd"] = max(worst["offd"], offd)
        worst["parity"] = max(worst["parity"], parity)
        pats.add(tuple(int(x) for x in pat))
    print(f"  reflection-symmetric m, {N} random draws (seed 20261023):")
    print(f"    Takagi U^T m U diagonal (1e-10 rel):          {'PASS' if takagi_ok else 'FAIL'}")
    c1 = worst["rows"] < 1e-9 and worst["t23"] < 1e-7
    print(f"    max ||U_mu j|-|U_tau j||  = {worst['rows']:.2e}, max |theta23-45| = {worst['t23']:.2e} deg -> {'PASS' if c1 else 'FAIL'}")
    c2 = worst["sind"] < 1e-7 and signs == {-1, 1}
    print(f"    max ||sin delta|-1| (theta13>0) = {worst['sind']:.2e}, signs seen {sorted(signs)} -> {'PASS' if c2 else 'FAIL'}")
    c3 = worst["offd"] < 1e-9 and worst["parity"] < 1e-9
    print(f"    rho = U^+ S U* diagonal +-1: offdiag {worst['offd']:.2e}, |parity|-1 {worst['parity']:.2e} -> {'PASS' if c3 else 'FAIL'}")
    print(f"    CP-parity patterns seen: {sorted(pats)}")
    ok &= takagi_ok and c1 and c2 and c3

    gen = [check(generic_matrix(rng)) for _ in range(200)]
    frac_45 = np.mean([abs(g[2] - 45.0) < 1e-6 for g in gen])
    frac_diag = np.mean([g[5] < 1e-6 for g in gen])
    c4 = frac_45 == 0 and frac_diag == 0
    print(f"  control, generic complex-symmetric m (200 draws): fraction with theta23=45: {frac_45:.3f}, "
          f"with rho diagonal: {frac_diag:.3f} -> {'PASS' if c4 else 'FAIL'}")
    ok &= c4
    print("\n" + ("D2 CHECKS PASS" if ok else "D2 CHECK FAILED"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
