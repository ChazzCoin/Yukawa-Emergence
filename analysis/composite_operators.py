#!/usr/bin/env python3
"""
composite_operators.py -- audit (iv): do the Hermitian composites Q1 = (1/2){H_e,H_nu} and
Q2 = H_e H_nu H_e evade Prop. "minimal" (the linear class A = alpha H_e + beta H_nu), and does
their eigenframe on E still reduce to alignment {0, pi/4, pi/2}?

Setup: 2x2 restriction to E, as in the paper. X = H_e|_E (c_e, delta_e, xi_e),
Y = H_nu(theta)|_E (c_nu, delta_nu, xi_nu - theta); phi = 2(theta+psi0) is the Bloch angle between
them. a = delta_e/2, b = delta_nu/2.
Calculates (all by direct 2x2 matrix products, numpy):
  1. identity check  Q1 = c_nu X + c_e Y - c_e c_nu 1 + (a b) cos(phi) 1        (ground truth, algebra)
  2. identity check  traceless(Q2) = (c_e^2 - a^2) b m + 2a(c_e c_nu + a b cos phi) n   (Bloch vectors)
  3. does the traceless part of Q lie in span{H_e^0, H_nu^0} with theta-INDEPENDENT coefficients?
     (fit per theta; report the spread of the fitted coefficients over theta)
  4. self-consistent eigenframe: zeros of ||[Q, H_nu(theta)]||_F over theta
  5. stationary points of ||[Q2, H_nu(theta)]||^2 vs the closed form
     cos phi* = (-p + sqrt(p^2 + 8 q^2)) / (4 q), p = c_e c_nu, q = a b
  6. physical positivity: H = Y Y^dag >= 0 implies c >= delta/2, hence p/q >= 1.
Inputs (illustrative; masses in units of m_tau^2 for the charged leptons, quasi-degenerate nu pair):
  m_mu = 0.10566, m_tau = 1.77686 (GeV, PDG, rounded; version not pinned) -> c_e, delta_e
  nu pair eigenvalues (1.00, 1.05) (illustrative, not data); xi_e - xi_nu = psi0 = 0.35 rad.
Run: python3 analysis/composite_operators.py
Made by: crew flavor-obstruction-paper, item #1 (R1), data-analyst.
"""
import numpy as np
from flavor_common import bloch2, rot2, SX, SZ

NG = 60001
TH = np.linspace(0, np.pi, NG)


def fro2(M):
    return float(np.sum(np.abs(M)**2))


def comm(A, B):
    return A @ B - B @ A


def params(c_e, d_e, c_n, d_n, psi0):
    X = bloch2(d_e, psi0, c_e)                       # xi_e = psi0, xi_nu = 0
    Y0 = bloch2(d_n, 0.0, c_n)
    return X, Y0


def Yt(Y0, t):
    R = rot2(t)
    return R.conj().T @ Y0 @ R


def bloch_vec(M):
    """(x, z) Bloch components of the traceless part, sigma_x and sigma_z coefficients."""
    return np.array([np.real(np.trace(M @ SX)) / 2, np.real(np.trace(M @ SZ)) / 2])


def run_case(label, c_e, d_e, c_n, d_n, psi0=0.35, verbose=True):
    X, Y0 = params(c_e, d_e, c_n, d_n, psi0)
    a, b = d_e / 2, d_n / 2
    p, q = c_e * c_n, a * b
    print(f"\n--- {label}: c_e={c_e:.6g} d_e={d_e:.6g} c_nu={c_n:.6g} d_nu={d_n:.6g} psi0={psi0}   p/q = {p/q:.6g} ---")
    m1 = m2 = m3 = 0.0
    maxerr1 = maxerr2 = 0.0
    coef1, coef2 = [], []
    z1, z2 = [], []
    f1n = np.empty(NG)
    f2 = np.empty(NG)
    for i, t in enumerate(TH):
        Y = Yt(Y0, t)
        phi = 2 * (t + psi0)
        Q1 = 0.5 * (X @ Y + Y @ X)
        Q2 = X @ Y @ X
        # 1. identity for Q1
        pred1 = c_n * X + c_e * Y - c_e * c_n * np.eye(2) + a * b * np.cos(phi) * np.eye(2)
        maxerr1 = max(maxerr1, np.max(np.abs(Q1 - pred1)))
        # 2. Bloch-vector identity for Q2
        n = bloch_vec(X) / a
        m = bloch_vec(Y) / b
        V = (c_e**2 - a**2) * b * m + 2 * a * (c_e * c_n + a * b * np.cos(phi)) * n
        maxerr2 = max(maxerr2, np.max(np.abs(bloch_vec(Q2) - V)))
        f2[i] = fro2(comm(Q2, Y))
        f1n[i] = fro2(comm(Q1, Y))
        if i % 2000 == 0:
            # 3. fit traceless(Q) = al * X0 + be * Y0t   (X0, Yt0 traceless parts)
            X0 = X - np.trace(X) / 2 * np.eye(2)
            Y00 = Y - np.trace(Y) / 2 * np.eye(2)
            B = np.stack([bloch_vec(X0), bloch_vec(Y00)], axis=1)
            for Q, store in ((Q1, coef1), (Q2, coef2)):
                tl = bloch_vec(Q)
                sol, res, rk, sv = np.linalg.lstsq(B, tl, rcond=None)
                store.append((sol[0], sol[1], np.linalg.norm(B @ sol - tl)))
    print(f"  1. max|Q1 - [c_nu X + c_e Y - c_e c_nu 1 + a b cos(phi) 1]| = {maxerr1:.2e}")
    print(f"  2. max|Bloch(Q2) - closed form|                             = {maxerr2:.2e}")
    for nm, co in (("Q1", coef1), ("Q2", coef2)):
        co = np.array(co)
        print(f"  3. {nm}: traceless part = alpha(theta) X0 + beta(theta) Y0(theta): fit residual max={co[:,2].max():.2e};"
              f" alpha in [{co[:,0].min():.6g},{co[:,0].max():.6g}], beta in [{co[:,1].min():.6g},{co[:,1].max():.6g}]"
              f"  -> {'theta-INDEPENDENT coefficients (inside linear class + scalar)' if (np.ptp(co[:,0]) < 1e-9 and np.ptp(co[:,1]) < 1e-9) else 'theta-DEPENDENT coefficients (outside the fixed-coefficient linear class)'}")
    print(f"     Q1 coefficients predicted (alpha,beta)=(c_nu,c_e)=({c_n:.6g},{c_e:.6g})")
    # 4. self-consistent eigenframe: zeros of ||[Q,Y(theta)]||^2 (relative threshold 1e-6 of the max)
    for nm, g in (("Q1", f1n), ("Q2", f2)):
        zeros = [TH[i] for i in range(1, NG - 1) if g[i] < g[i - 1] and g[i] <= g[i + 1] and g[i] < 1e-6 * g.max()]
        print(f"  4. {nm}: theta of zeros of ||[{nm},H_nu(theta)]||^2 : "
              + (", ".join(f"{np.degrees(z):.3f} deg (phi={np.degrees(2*(z+psi0))%360:.2f} deg)" for z in zeros) or "none")
              + f"   [aligned = theta = -psi0 mod 90 deg = {np.degrees(-psi0)%90:.3f}]")
    # 5. stationary points of f2
    d = np.diff(f2)
    idx = [i + 1 for i in range(len(d) - 1) if d[i] * d[i + 1] < 0]
    print("  5. Q2: stationary theta of ||[Q2,H_nu]||^2 over [0,180) deg: "
          + ", ".join(f"{np.degrees(TH[i]):.3f}({'min' if f2[i] < f2[max(i-2000,0)] and f2[i] < f2[min(i+2000,NG-1)] else 'max'},f={f2[i]:.3e})" for i in idx))
    if q > 0:
        disc = np.sqrt(p**2 + 8 * q**2)
        cstar = (-p + disc) / (4 * q)
        phistar = np.arccos(np.clip(cstar, -1, 1))
        print(f"     closed form cos phi* = {cstar:.9f} -> phi* = {np.degrees(phistar):.4f} deg  -> theta* = {np.degrees(phistar/2 - psi0)%90:.4f} deg"
              f" (and mod 90 siblings); aligned = phi in {{0,180}} deg; pi/4 point = phi = 90 deg")


def main():
    mu, tau = 0.10566, 1.77686
    c_e, d_e = (mu**2 + tau**2) / 2, tau**2 - mu**2
    # ground-truth structural check lives inside run_case (identities 1,2 at 1e-14)
    run_case("physical-like: charged mu-tau block, quasi-degenerate nu pair (1.00,1.05)", c_e, d_e, 1.025, 0.05)
    run_case("physical-like, wide nu splitting (nu pair (0.2, 1.8))", c_e, d_e, 1.0, 1.6)
    run_case("PSD extreme: lowest eigenvalue 0 for both (p/q = 1)", 1.0, 2.0, 0.5, 1.0)
    run_case("INDEFINITE toy (not m m^dag): c_nu small so p/q < 1", 1.0, 1.0, 0.1, 1.0)

    print("\n--- 6. positivity: p/q = (c_e/a)(c_nu/b) over random PSD 2x2 pairs (eigenvalues >= 0) ---")
    rng = np.random.default_rng(7)
    lam = rng.uniform(0, 5, size=(100000, 4))
    c_e_, a_ = (lam[:, 0] + lam[:, 1]) / 2, np.abs(lam[:, 0] - lam[:, 1]) / 2
    c_n_, b_ = (lam[:, 2] + lam[:, 3]) / 2, np.abs(lam[:, 2] - lam[:, 3]) / 2
    r = (c_e_ * c_n_) / (a_ * b_)
    print(f"  100000 random draws: min p/q = {r.min():.6f}  (>=1 always: {bool(np.all(r >= 1 - 1e-12))})")
    print("  analytic: c >= a  for PSD (c = (l1+l2)/2 >= |l1-l2|/2), so p/q >= 1, equality iff both smaller eigenvalues are 0.")


if __name__ == "__main__":
    main()
