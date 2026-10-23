#!/usr/bin/env python3
"""
full3x3_fourier.py -- audit (i): does theta enter functionals of the FULL 3x3 pair
(H_e, H_nu) only through cos 2(theta+psi0)?

Calculates: for H_nu(theta) = W(theta)^dag H_nu W(theta), W = 1 (+) R(theta) (rotation in
the fixed 2-3 plane E), the theta-dependence of
  f1 = Tr(H_e H_nu)          f2 = ||[H_e,H_nu]||_F^2
  f3 = Tr((H_e H_nu)^2)      f4 = 3 - sum_i |(V_e^dag V_nu)_ii|^2   (eigenframe misalignment)
as a Fourier series over theta in [0,2pi): amplitude per harmonic k, odd-k (frequency-1
family) vs even-k content, and the asymmetry from "function of cos 2(theta+psi0)"
(= period pi AND reflection symmetry f(p+x)=f(p-x) for some p).

Inputs (all dimensionless, spectra normalised to the largest charged-lepton eigenvalue):
  H_e eigenvalues  (m_e^2,m_mu^2,m_tau^2)/m_tau^2 = (8.3e-8, 3.5e-3, 1)
      from m_e=0.511 MeV, m_mu=105.66 MeV, m_tau=1776.86 MeV (PDG; values rounded, version not pinned)
  H_nu: block form diag(m1^2) (+) [E block, eigenvalues L2 < L3, gap d_nu, orientation xi_nu]
      defaults (m1^2, L2, L3) = (0.90, 1.00, 1.05): illustrative quasi-degenerate, NOT fitted data
  H_e mixing relative to E: PDG-form P(t12e,t13e,t23e,dl), H_e = P diag P^dag.
      default (33.4, 8.6, 45) deg, dl = 0: illustrative angle sizes, NOT sourced NuFIT values.
Method: direct numerics (numpy), N=720 samples of theta, FFT; ground-truth case first.
Run: python3 analysis/full3x3_fourier.py [seed]
Made by: crew flavor-obstruction-paper, item #1 (R1), data-analyst.
"""
import sys
import numpy as np
from flavor_common import pdg_U, rot2, embed_E, herm_from, fourier_amps

N = 720
TH = np.arange(N) * 2 * np.pi / N
EIG_E = np.array([0.511e-3**2, 0.10566**2, 1.77686**2]) / 1.77686**2
M1SQ, L2, L3 = 0.90, 1.00, 1.05
DEG = np.pi / 180


def build(t12e, t13e, t23e, dl, xi_nu, eig_e=EIG_E, m1sq=M1SQ, l2=L2, l3=L3, tilt_nu=None):
    """Return (He, Hnu0) in the E-adapted basis: Hnu0 block-diagonal, He generic."""
    P = pdg_U(t12e, t13e, t23e, dl)
    He = herm_from(eig_e, P)
    c, s = np.cos(xi_nu), np.sin(xi_nu)
    Vnu = np.eye(3, dtype=complex)
    Vnu[1:, 1:] = [[-s, c], [c, s]]                     # columns: l2 vec (smaller), l3 vec (larger)
    if tilt_nu is not None:                             # H_nu NOT block-diagonal wrt E: Vnu <- PDG(tilt_nu) . Vnu
        Vnu = pdg_U(*tilt_nu) @ Vnu
    Hnu0 = herm_from([m1sq, l2, l3], Vnu)
    return He, Hnu0


def functionals(He, Hnu0):
    f1, f2, f3, f4 = [], [], [], []
    we, Ve = np.linalg.eigh(He)
    for t in TH:
        W = embed_E(rot2(t))
        H = W.conj().T @ Hnu0 @ W
        f1.append(np.real(np.trace(He @ H)))
        C = He @ H - H @ He
        f2.append(np.real(np.sum(np.abs(C)**2)))
        f3.append(np.real(np.trace(He @ H @ He @ H)))
        wn, Vn = np.linalg.eigh(H)
        Wm = Ve.conj().T @ Vn
        f4.append(3 - np.sum(np.abs(np.diag(Wm))**2))
    return np.array([f1, f2, f3, f4])


NAMES = ["Tr(He Hnu)", "||[He,Hnu]||^2", "Tr((He Hnu)^2)", "misalign 3-sum|W_ii|^2"]


def summarize(F, label, verbose=True):
    out = []
    for nm, f in zip(NAMES, F):
        a = fourier_amps(f, 8)
        odd = np.sqrt(np.sum(a[1::2]**2))
        even = np.sqrt(np.sum(a[2::2]**2))
        # reflection asymmetry, exact from the Fourier series f = sum a_k cos(kx+phi_k):
        # rms[f(p+x)-f(p-x)]^2 = sum_k 2 a_k^2 sin^2(k p + phi_k); minimised over p on a fine grid
        c = np.fft.rfft(f)[1:13] / N
        k = np.arange(1, 13)
        def asym2(p):
            return np.sum(2 * (2 * np.abs(c))**2 * np.sin(k * p[:, None] + np.angle(c))**2, axis=1)
        pp = np.linspace(0, np.pi, 3601)
        p0 = pp[np.argmin(asym2(pp))]
        pf = np.linspace(p0 - 2e-3, p0 + 2e-3, 4001)       # local refinement
        best = np.sqrt(np.min(asym2(pf)))
        scale = f.max() - f.min()
        out.append((odd / scale if scale > 0 else 0.0, even / scale if scale > 0 else 0.0,
                    best / scale if scale > 0 else 0.0, a))
        if verbose:
            print(f"  {nm:24s} range={scale:.3e}  odd-k amp/range={out[-1][0]:.2e}  "
                  f"even-k amp/range={out[-1][1]:.2e}  reflection-asym/range={out[-1][2]:.2e}")
            print(f"     amps k=1..4: " + "  ".join(f"{x:.3e}" for x in a[1:5]))
    return out


def stationary_report(He, Hn, label):
    """Location of the stationary points / minima of f2 in theta vs the block-alignment angles
    (theta where [He|_E, Hnu(theta)|_E] = 0, i.e. the paper's 2x2 alignment)."""
    fine = np.linspace(0, np.pi, 18001)
    f2, blk = [], []
    for t in fine:
        W = embed_E(rot2(t))
        H = W.conj().T @ Hn @ W
        C = He @ H - H @ He
        f2.append(np.real(np.sum(np.abs(C)**2)))
        Cb = He[1:, 1:] @ H[1:, 1:] - H[1:, 1:] @ He[1:, 1:]
        blk.append(np.real(np.sum(np.abs(Cb)**2)))
    f2, blk = np.array(f2), np.array(blk)
    d = np.diff(f2)
    st = [fine[i + 1] for i in range(len(d) - 1) if d[i] * d[i + 1] < 0]
    kinds = ["min" if f2[np.argmin(abs(fine - t))] < f2[(np.argmin(abs(fine - t)) + 300) % len(fine)] and
             f2[np.argmin(abs(fine - t))] < f2[(np.argmin(abs(fine - t)) - 300) % len(fine)] else "max" for t in st]
    al = [fine[i] for i in range(1, len(blk) - 1) if blk[i] < blk[i - 1] and blk[i] <= blk[i + 1] and blk[i] < 1e-8]
    print(f"  [{label}] f2 stationary theta (deg, over [0,180)): "
          + ", ".join(f"{np.degrees(t):.2f}({k})" for t, k in zip(st, kinds)))
    print(f"  [{label}] f2 value at minima: " + ", ".join(f"{f2[np.argmin(abs(fine - t))]:.3e}" for t, k in zip(st, kinds) if k == 'min'))
    print(f"  [{label}] block-alignment theta (2x2 [He_E,Hnu_E]=0), deg: " + ", ".join(f"{np.degrees(t):.2f}" for t in al))


def main():
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 20261023
    print("N samples =", N, " seed =", seed)
    print("H_e eigenvalues =", EIG_E)
    print("H_nu block: (m1^2, L2, L3) =", (M1SQ, L2, L3))

    print("\n=== GROUND TRUTH: H_e block-diagonal wrt E (t12e=t13e=0), paper's 2x2 situation ===")
    He, Hn = build(0, 0, 25 * DEG, 0, 0.3)
    F = functionals(He, Hn)
    gt = summarize(F, "GT")
    ok = all(r[0] < 1e-10 and r[2] < 1e-5 for r in gt[:4])
    # analytic: Tr(He Hnu)= const + (1/2) de dnu cos2(th+psi0) on E; check amplitude of k=2
    de = EIG_E[2] - EIG_E[1]; dn = L3 - L2
    print(f"  paper amplitude (1/2) de dnu = {0.5*de*dn:.6e};  measured k=2 amp of f1 = {fourier_amps(F[0])[2]:.6e}")
    ok &= abs(fourier_amps(F[0])[2] - 0.5 * de * dn) < 1e-10
    print("  GROUND TRUTH:", "PASS (odd-k absent, f1 amplitude = paper's)" if ok else "FAIL")

    stationary_report(He, Hn, "GT")
    print("\n=== NOVEL: H_e tilted wrt E, default (33.4, 8.6, 45) deg, dl=0, xi_nu=0.3 ===")
    He, Hn = build(33.4 * DEG, 8.6 * DEG, 45 * DEG, 0, 0.3)
    F = functionals(He, Hn)
    summarize(F, "tilted")
    stationary_report(He, Hn, "tilted")

    print("\n=== NOVEL: same with CP phase dl = 197 deg (complex H_e) ===")
    He, Hn = build(33.4 * DEG, 8.6 * DEG, 45 * DEG, 197 * DEG, 0.3)
    summarize(functionals(He, Hn), "tilted+CP")

    print("\n=== NOVEL: tilt only through theta13 (t12e=0, t13e=8.6 deg) ===")
    He, Hn = build(0, 8.6 * DEG, 45 * DEG, 0, 0.3)
    summarize(functionals(He, Hn), "t13 only")

    print("\n=== NOVEL: BOTH tilted wrt the fixed plane E (E is NOT an eigenblock of H_nu): H_e (33.4,8.6,45) deg, H_nu extra tilt (0,5,0) deg ===")
    He, Hn = build(33.4 * DEG, 8.6 * DEG, 45 * DEG, 0, 0.3, tilt_nu=(0, 5 * DEG, 0, 0))
    summarize(functionals(He, Hn), "both tilted")
    print("=== same, H_nu extra tilt (3,8.6,0) deg, dl_e=197 deg ===")
    He, Hn = build(33.4 * DEG, 8.6 * DEG, 45 * DEG, 197 * DEG, 0.3, tilt_nu=(3 * DEG, 8.6 * DEG, 0, 0))
    summarize(functionals(He, Hn), "both tilted 2")

    print("\n=== RANDOM SWEEP: 300 draws, t12e in [30,36], t13e in [7,10], t23e in [40,50] deg,")
    print("    dl in [0,2pi), xi_nu in [0,pi/2), dnu in [0.01,0.3], m1sq in [0.5,1.5]; odd-k/range and reflection-asym/range ===")
    rng = np.random.default_rng(seed)
    res = np.zeros((300, 4, 3))
    for i in range(300):
        dnu_ = rng.uniform(0.01, 0.3)
        He, Hn = build(rng.uniform(30, 36) * DEG, rng.uniform(7, 10) * DEG, rng.uniform(40, 50) * DEG,
                       rng.uniform(0, 2 * np.pi), rng.uniform(0, np.pi / 2),
                       m1sq=rng.uniform(0.5, 1.5), l2=1.0, l3=1.0 + dnu_)
        for j, r in enumerate(summarize(functionals(He, Hn), "r", verbose=False)):
            res[i, j] = r[:3]
    for j, nm in enumerate(NAMES):
        o, e, a = res[:, j, 0], res[:, j, 1], res[:, j, 2]
        print(f"  {nm:24s} odd/range: median={np.median(o):.2e} max={o.max():.2e} | "
              f"reflection-asym/range: median={np.median(a):.2e} max={a.max():.2e} | frac(odd>1e-6)={np.mean(o>1e-6):.2f}")


if __name__ == "__main__":
    main()
