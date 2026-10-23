#!/usr/bin/env python3
"""
d1_rg_scale.py

What it calculates: the "D1: RG scale" check for the theta23 no-go paper. From a mu-tau
symmetric seed kappa0 (theta23 = 45 deg, theta13 = 0), apply the leading-log flavor-diagonal
transport kappa(M_Z) = I kappa0 I, I = diag(1, 1, 1+eps), and report
  (1) ground truth: eps = 0 gives theta23 = 45 deg, theta13 = 0 exactly, kappa0 invariant under
      the mu-tau reflection; also first-order theory vs exact diagonalisation;
  (2) the sign of theta23 - 45 deg against sign(eps) and sign(Dm2_atm), and the mapping to the
      lab record's rule  sign(theta23(M_Z)-45) = +sign(C)*sign(Dm2_atm);
  (3) which pair (solar nu1-nu2 or atmospheric nu2-nu3) the transport rotates more, vs lightest mass.

Inputs (all from analysis/data/nufit6_central.csv = NuFIT 6.0 Table 1 best fits, arXiv:2410.05380 v2;
copied 2026-10-23 from the Flavor worktree; variant named on every number):
  theta12 [deg]; dm21 [1e-5 eV^2]; dm3l [1e-3 eV^2] (dm31 for NO, dm32 for IO).
  Seed theta23 = 45 deg, theta13 = 0, Majorana phases 0 are choices fixed by the item, not data.
  m_lightest [eV]: scanned (command-line list), default 0, 0.01, 0.05, 0.1, 0.3.
  eps [dimensionless]: |eps| = 1e-4 for the sign test (either sign), 1e-8 for first-order check.
Method: kappa = U D U^T (real, phases 0), PDG U = R23 R12 at theta13 = 0. Eigenvectors of the
  real symmetric I kappa0 I matched to the seed's columns by largest overlap; angles from
  |Ue3|=sin13, tan23=|Umu3|/|Utau3|, tan12=|Ue2|/|Ue1|. Pair rotation = arcsin of |(U_seed^T U_new)_ij|.
  First order: pair (i,j) rotates by eps*U_tau,i*U_tau,j*(m_i+m_j)/(m_i-m_j)  (second-order
  perturbation of the symmetric matrix d kappa = eps(P kappa0 + kappa0 P), P = e_tau e_tau^T).
  Record: flavor-obstruction-lab/wt/calculations/theta23-rg-running/README.md (verdict, Method):
  kappa(M_Z) = E kappa(M1) E, eps = -C*J, J>0, C = -3/2 (SM), +1 (MSSM).
Run:  python3 analysis/d1_rg_scale.py [m_lightest ...]
Made by: crew flavor-obstruction-paper, data-analyst-2, item #6 (R1); first committed in the crew's Domain.Alignment worktree (aa50cef), copied here unchanged otherwise.
"""
import csv
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "data", "nufit6_central.csv")
C_REC = {"SM": -1.5, "MSSM": 1.0}          # record: Antusch et al. Eq.(4)
EPS_SIGN, EPS_FO = 1e-4, 1e-8


def rd(path):
    rows = [r for r in open(path) if not r.startswith("#")]
    return list(csv.DictReader(rows))


def masses(order, mlight, dm21, dm3l):
    """Return (m1, m2, m3) in eV. dm21, dm3l in eV^2; dm3l = dm31 (NO) or dm32 (IO, <0)."""
    if order == "NO":
        m1 = mlight
        return m1, np.sqrt(m1**2 + dm21), np.sqrt(m1**2 + dm3l)
    m3 = mlight
    m2 = np.sqrt(m3**2 - dm3l)             # dm32 < 0
    return np.sqrt(m2**2 - dm21), m2, m3


def seed_U(th12):
    c, s = np.cos(th12), np.sin(th12)
    R12 = np.array([[c, s, 0], [-s, c, 0], [0, 0, 1.0]])
    c2 = s2 = np.sqrt(0.5)
    R23 = np.array([[1, 0, 0], [0, c2, s2], [0, -s2, c2]])
    return R23 @ R12


def angles(U):
    s13 = abs(U[0, 2])
    return (np.degrees(np.arctan2(abs(U[0, 1]), abs(U[0, 0]))),
            np.degrees(np.arcsin(s13)),
            np.degrees(np.arctan2(abs(U[1, 2]), abs(U[2, 2]))))


def transport(k0, eps):
    I = np.diag([1, 1, 1 + eps])
    return I @ k0 @ I


def match(Us, k):
    """Eigenvectors of real symmetric k, columns matched to seed columns Us; return U, eigenvalues."""
    w, V = np.linalg.eigh(k)
    U = np.zeros((3, 3)); lam = np.zeros(3)
    for j in range(3):
        o = np.abs(V.T @ Us[:, j]); i = int(np.argmax(o))
        U[:, j] = V[:, i] * np.sign(V[:, i] @ Us[:, j]); lam[j] = w[i]
    return U, lam


def pair_rot(Us, U, i, j):
    R = Us.T @ U
    return np.degrees(np.arcsin(min(1.0, abs(R[i, j]))))


def main():
    mls = [float(x) for x in sys.argv[1:]] or [0.0, 0.01, 0.05, 0.1, 0.3]
    ok = True
    data = rd(CSV)
    print("D1: RG scale -- kappa(M_Z) = I kappa0 I, I = diag(1,1,1+eps); seed mu-tau symmetric")
    print(f"inputs: {os.path.relpath(CSV)} (NuFIT 6.0), m_lightest scan = {mls} eV, "
          f"|eps|(sign test) = {EPS_SIGN}, |eps|(first-order) = {EPS_FO}")
    # mu-tau reflection: nu_e -> nu_e, nu_mu -> -nu_tau, nu_tau -> -nu_mu
    Pm = np.array([[1, 0, 0], [0, 0, -1], [0, -1, 0.0]])
    rows = []                              # (variant, order, ml, sign tuples, rots)
    for r in data:
        var, order = r["variant"], r["ordering"]
        th12 = np.radians(float(r["theta12"]))
        dm21 = float(r["dm21_1e-5"]) * 1e-5
        dm3l = float(r["dm3l_1e-3"]) * 1e-3
        sgn_atm = 1 if dm3l > 0 else -1
        for ml in mls:
            m = np.array(masses(order, ml, dm21, dm3l))
            Us = seed_U(th12)
            k0 = Us @ np.diag(m) @ Us.T
            tag = f"[{var} {order} m_light={ml:g}]"
            # (1) ground truth
            a0 = angles(Us)
            refl = np.max(np.abs(Pm @ k0 @ Pm.T - k0))
            U0, lam0 = match(Us, transport(k0, 0.0))
            a00 = angles(U0)
            gt = (abs(a00[2] - 45) < 1e-10 and a00[1] < 1e-6 and refl < 1e-14
                  and np.max(np.abs(lam0 - m)) < 1e-14)
            print(f"{tag} eps=0: theta23={a00[2]:.12f} theta13={a00[1]:.3e} deg, "
                  f"reflection residual={refl:.1e}, mass residual={np.max(np.abs(lam0-m)):.1e}"
                  f" -> {'PASS' if gt else 'FAIL'}")
            ok &= gt
            # first-order check at small eps
            Ue, _ = match(Us, transport(k0, EPS_FO))
            Ut = Us[2]
            pr = {}
            for (i, j) in [(0, 1), (1, 2)]:
                if m[i] != m[j]:
                    pr[(i, j)] = np.degrees(abs(EPS_FO * Ut[i] * Ut[j] * (m[i] + m[j]) / (m[i] - m[j])))
            num = {p: pair_rot(Us, Ue, *p) for p in pr}
            relf = max(abs(num[p] - pr[p]) / pr[p] for p in pr if pr[p] > 0) if pr else 0
            fo = relf < 1e-3
            print(f"{tag} first-order pair rotation, solar/atm: theory "
                  f"{pr[(0,1)]:.3e}/{pr[(1,2)]:.3e} deg, exact {num[(0,1)]:.3e}/{num[(1,2)]:.3e}, "
                  f"max rel dev {relf:.1e} -> {'PASS' if fo else 'FAIL'}")
            ok &= fo
            # (2) sign test, both signs of eps
            sg = []
            for eps in (+EPS_SIGN, -EPS_SIGN):
                Un, _ = match(Us, transport(k0, eps))
                d23 = angles(Un)[2] - 45.0
                sg.append((eps, d23))
            s_plus = np.sign(sg[0][1]) * np.sign(sg[0][0]) * sgn_atm   # sign(d23)/(sign(eps) sign(Dm))
            s_minus = np.sign(sg[1][1]) * np.sign(sg[1][0]) * sgn_atm
            rows.append((var, order, ml, s_plus, s_minus, sg))
            # (3) solar vs atm rotation at same |eps|
            Ub, _ = match(Us, transport(k0, EPS_SIGN))
            rs, ra = pair_rot(Us, Ub, 0, 1), pair_rot(Us, Ub, 1, 2)
            print(f"{tag} eps=+{EPS_SIGN}: d23={sg[0][1]:+.4e} deg, eps=-{EPS_SIGN}: d23={sg[1][1]:+.4e} deg; "
                  f"rot solar={rs:.4e}, atm={ra:.4e} deg, solar/atm={rs/ra:.3f}")
            rows[-1] += (rs, ra, m)
    # sign summary
    print("\n(2) sign of theta23-45 relative to sign(eps)*sign(Dm2_atm):")
    svals = {r[3] for r in rows} | {r[4] for r in rows}
    print(f"    distinct values of sign(d23)/(sign(eps)*sign(Dm2_atm)) over all rows and both eps signs: {sorted(svals)}")
    conv = svals == {-1.0}
    print(f"    derived convention: sign(theta23-45) = -sign(eps)*sign(Dm2_atm) -> {'PASS' if conv else 'FAIL'}")
    ok &= conv
    print("    mapping: kappa(M_Z)=E kappa(M1) E with eps = -C*J, J>0 (record Method; C=-3/2 SM, +1 MSSM)")
    print("    => sign(eps) = -sign(C) => sign(theta23-45) = +sign(C)*sign(Dm2_atm) (the record's rule)")
    for var in sorted({r[0] for r in rows}):
        for order in ("NO", "IO"):
            sgn_atm = 1 if order == "NO" else -1
            for model, C in C_REC.items():
                eps = -C * 1.0                      # sign of eps for J>0
                rr = [x for x in rows if x[0] == var and x[1] == order]
                # d23 sign at the +/-EPS_SIGN entry matching sign(eps)
                ds = {np.sign(x[5][0][1] if eps > 0 else x[5][1][1]) for x in rr}
                oct_ = "upper" if ds == {1.0} else ("lower" if ds == {-1.0} else "MIXED")
                rule = "upper" if np.sign(C) * sgn_atm > 0 else "lower"
                good = oct_ == rule
                ok &= good
                print(f"    [{var} {order} {model}] C={C:+.1f} -> sign(eps)={int(np.sign(eps)):+d}: "
                      f"octant {oct_} (record rule: {rule}) -> {'PASS' if good else 'FAIL'}")
    # (3) summary
    print("\n(3) rotation inside the solar pair vs the atmospheric pair, same eps:")
    for var in sorted({r[0] for r in rows}):
        for order in ("NO", "IO"):
            line = []
            for x in [x for x in rows if x[0] == var and x[1] == order]:
                line.append(f"m={x[2]:g}:{x[6]/x[7]:.2f}")
            print(f"    [{var} {order}] solar/atm rotation ratio: " + "  ".join(line))
    qd = [x for x in rows if x[2] >= 0.1]
    qd_ok = all(x[6] > x[7] for x in qd)
    print(f"    quasi-degenerate rows (m_lightest >= 0.1 eV, {len(qd)} rows): solar > atm in all -> "
          f"{'PASS' if qd_ok else 'FAIL'}")
    ok &= qd_ok
    for var in sorted({r[0] for r in rows}):
        for order in ("NO", "IO"):
            below = [x[2] for x in rows if x[0] == var and x[1] == order and x[6] <= x[7]]
            print(f"    [{var} {order}] masses in the scan where solar <= atm: {below if below else 'none'}")
    print("    (qualitative; first order ratio = {[(m1+m2)^2/dm21]*s12 c12} / {[(m2+m3)^2/|dm32|]*c12}, up to O(1) factors)")
    for r in data:
        th12 = np.radians(float(r["theta12"])); dm21 = float(r["dm21_1e-5"]) * 1e-5
        dm3l = float(r["dm3l_1e-3"]) * 1e-3

        def ratio(ml):
            m = np.array(masses(r["ordering"], ml, dm21, dm3l)); Us = seed_U(th12)
            Ub, _ = match(Us, transport(Us @ np.diag(m) @ Us.T, EPS_SIGN))
            return pair_rot(Us, Ub, 0, 1) / pair_rot(Us, Ub, 1, 2)
        lo, hi = 0.0, 0.01
        if ratio(lo) >= 1:
            print(f"    [{r['variant']} {r['ordering']}] solar > atm already at m_lightest=0 (ratio {ratio(lo):.2f}); no crossover")
            continue
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if ratio(mid) < 1 else (lo, mid)
        print(f"    [{r['variant']} {r['ordering']}] crossover solar=atm at m_lightest = {hi*1e3:.3f} meV; solar > atm above it (scan to 0.3 eV)")
    print("\n" + ("D1 CHECKS PASS" if ok else "D1 CHECK FAILED"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
