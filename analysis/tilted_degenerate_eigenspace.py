#!/usr/bin/env python3
"""
tilted_degenerate_eigenspace.py -- audit (ii): E as the degenerate eigenspace of H_nu
with realistic PDG mixing (theta12, theta13 != 0), m2 = m3.

Calculates:
  A. whether the 2-3 FLAVOR block of H_nu = U diag(m1^2,m^2,m^2) U^dag is proportional to 1
     (the paper's premise "E = mu-tau flavor plane, H_nu|_E = lambda 1"), and its eigen-gap;
  B. what a rotation inside the true degenerate eigenspace, U -> U (1 (+) R(t)), does to
     H_nu (nothing) and to the PDG angles read from U (theta12, theta13, theta23, |delta|),
     |U_e2|, |U_e3|;
  C. the same rotation when the degeneracy is lifted (m3^2 = m^2 + eps): H_nu changes ~ eps,
     i.e. the rotation becomes physical.
Inputs (illustrative magnitudes, NOT fitted data and NOT tied to a NuFIT version):
  theta12=33.4, theta13=8.6, theta23=49.0 deg; delta=197 deg (or 0);
  m^2 = 2.5e-3 eV^2 (scale of the degenerate pair); Delta = m1^2 - m^2 = -7.4e-5 eV^2
  (solar-size splitting); eps in C: 2.5e-3 eV^2 down to 2.5e-6.
Method: direct numerics (numpy); PDG extraction from |U_ij| (flavor_common.pdg_angles).
Run: python3 analysis/tilted_degenerate_eigenspace.py
Made by: crew flavor-obstruction-paper, item #1 (R1), data-analyst.
"""
import numpy as np
from flavor_common import pdg_U, pdg_angles, rot2, embed_E, herm_from

D = np.pi / 180
M2 = 2.5e-3
DELTA = -7.4e-5


def ang(U):
    t12, t13, t23, dl = pdg_angles(U)
    return np.degrees(t12), np.degrees(t13), np.degrees(t23), np.degrees(dl)


def rotated(U, t):
    return U @ embed_E(rot2(t))                      # right-rotation of mass columns 2,3


def main():
    print("inputs: m^2 = %.3e eV^2, Delta = m1^2 - m^2 = %.3e eV^2" % (M2, DELTA))

    print("\n=== GROUND TRUTH: theta12 = theta13 = 0 ===")
    U0 = pdg_U(0, 0, 49 * D, 0)
    Hn = herm_from([M2 + DELTA, M2, M2], U0)
    blk = Hn[1:, 1:]
    print("  2-3 flavor block of H_nu - m^2*1 : max|.| = %.2e (expect 0 => block ∝ 1)" % np.max(np.abs(blk - M2 * np.eye(2))))
    for t in [0, 10, 25, 40]:
        a = ang(rotated(U0, t * D))
        print("  t=%2d deg -> theta12=%.3f theta13=%.3f theta23=%.3f  (expect 0, 0, 49-t=%d: R(t) shifts theta23 by exactly t, sign set by the R convention)" % (t, a[0], a[1], a[2], 49 - t))
    ok_gt = all(abs(ang(rotated(U0, t * D))[2] - (49 - t)) < 1e-9 and ang(rotated(U0, t * D))[1] < 1e-9 for t in [10, 25, 40])
    print("  GROUND TRUTH:", "PASS" if ok_gt else "FAIL")

    for dl in [0.0, 197.0]:
        print("\n=== REALISTIC tilt: theta12=33.4 theta13=8.6 theta23=49 delta=%g deg, m2=m3 ===" % dl)
        U = pdg_U(33.4 * D, 8.6 * D, 49 * D, dl * D)
        Hn = herm_from([M2 + DELTA, M2, M2], U)
        blk = Hn[1:, 1:]
        ev = np.linalg.eigvalsh(blk)
        v1 = np.abs(U[:, 0])**2
        print("  A. 2-3 flavor block of H_nu: eigenvalues = %.6e, %.6e ; gap = %.3e eV^2" % (ev[0], ev[1], ev[1] - ev[0]))
        print("     predicted gap |Delta|*(1-|U_e1|^2) = %.3e  (1-|U_e1|^2 = %.4f);  gap/|Delta| = %.4f" %
              (abs(DELTA) * (1 - v1[0]), 1 - v1[0], (ev[1] - ev[0]) / abs(DELTA)))
        print("     => flavor 2-3 block is NOT ∝ 1; true degenerate eigenspace E = span(U[:,1],U[:,2]) is tilted out of the flavor mu-tau plane")
        print("     check H_nu E = m^2 E : max|H_nu U[:,j] - m^2 U[:,j]| (j=2,3) = %.2e" %
              max(np.max(np.abs(Hn @ U[:, j] - M2 * U[:, j])) for j in (1, 2)))
        print("  B. rotate inside E: U(t) = U (1 (+) R(t))")
        print("     t(deg)   theta12   theta13   theta23   |delta|   |Ue2|     |Ue3|    |Ue2|^2+|Ue3|^2  max|dH_nu|")
        for t in [0, 5, 10, 20, 30, 45, 60, 90]:
            Ut = rotated(U, t * D)
            a = ang(Ut)
            Ht = herm_from([M2 + DELTA, M2, M2], Ut)
            print("     %5.1f  %8.3f  %8.3f  %8.3f  %8.3f  %.5f  %.5f  %.12f  %.2e" %
                  (t, a[0], a[1], a[2], abs(a[3]), abs(Ut[0, 1]), abs(Ut[0, 2]),
                   abs(Ut[0, 1])**2 + abs(Ut[0, 2])**2, np.max(np.abs(Ht - Hn))))
        ts = np.linspace(0, np.pi, 180001)
        t13 = np.array([np.degrees(np.arcsin(abs(rotated(U, t)[0, 2]))) for t in ts[::50]])
        tt = ts[::50]
        print("     theta13(t) over t in [0,180) deg: min=%.4f max=%.4f ; theta13 -> 0 at t=%.2f deg" %
              (t13.min(), t13.max(), np.degrees(tt[np.argmin(t13)])))
        i0 = np.argmin(t13)
        a = ang(rotated(U, tt[i0]))
        print("     at that t: theta12=%.3f theta13=%.4f theta23=%.3f (an orbit point with theta13=0, same H_nu)" % (a[0], a[1], a[2]))
        # frame-independent content: first column of U (nu_1 direction) is invariant
        print("     invariant under t: |U_e1|,|U_mu1|,|U_tau1| = %s ; max drift over t = %.2e" %
              (np.round(np.abs(U[:, 0]), 6),
               max(np.max(np.abs(np.abs(rotated(U, t)[:, 0]) - np.abs(U[:, 0]))) for t in ts[::5000])))
        # observable: transition probability from H_nu itself, any baseline phase
        x = 3.7e3
        def Pmt(Hm):
            w, V = np.linalg.eigh(Hm)
            S = V @ np.diag(np.exp(-1j * w * x)) @ V.conj().T
            return abs(S[2, 1])**2
        print("     P(mu->tau) at an arbitrary phase x=%.1f eV^-2: t=0: %.12f ; t=17: %.12f ; t=63: %.12f" %
              (x, Pmt(herm_from([M2 + DELTA, M2, M2], U)), Pmt(herm_from([M2 + DELTA, M2, M2], rotated(U, 17 * D))),
               Pmt(herm_from([M2 + DELTA, M2, M2], rotated(U, 63 * D)))))

        print("  C. lift the degeneracy: m3^2 = m^2 + eps, rotate by t=20 deg inside E")
        print("     eps (eV^2)   max|H_nu(t=20)-H_nu(0)|   ratio to eps   theta12  theta13  theta23 (t=0 -> t=20)")
        for eps in [2.5e-3, 2.5e-4, 2.5e-5, 2.5e-6]:
            H0 = herm_from([M2 + DELTA, M2, M2 + eps], U)
            H1 = herm_from([M2 + DELTA, M2, M2 + eps], rotated(U, 20 * D))
            a0, a1 = ang(U), ang(rotated(U, 20 * D))
            print("     %.1e      %.3e                 %.4f        (%.2f,%.2f,%.2f)->(%.2f,%.2f,%.2f)" %
                  (eps, np.max(np.abs(H1 - H0)), np.max(np.abs(H1 - H0)) / eps, a0[0], a0[1], a0[2], a1[0], a1[1], a1[2]))


if __name__ == "__main__":
    main()
