#!/usr/bin/env python3
"""
full_flavor_pipeline.py

A single-script “full flavor” pipeline in the spectral-triple spirit:

1) Choose flavor sites on a cycle (default: 360-cycle sites [1,2,5]).
2) Build a κ-kernel K_ij = κ^{d(i,j)} (or exp(-κ d)).
3) Build sector Yukawas (Yu,Yd,Ye,Ynu) from a polynomial in K, dressed by
   left/right diagonal phase projectors.
4) Extract masses via SVD and mixings via left-unitary misalignments:
      V_CKM  = UuL† UdL
      U_PMNS = UeL† Uν
5) Optional Type-I seesaw with a Majorana MR and Takagi factorization.

Run:
  python full_flavor_pipeline.py
  python full_flavor_pipeline.py --config my_config.json --save out.json
"""

from __future__ import annotations
import argparse, json, math, cmath
from dataclasses import dataclass
from typing import Dict, Any, Tuple, List, Optional

import numpy as np

# -----------------------------
# Utilities
# -----------------------------

def wrap_pi(x: float) -> float:
    """Wrap angle to (-pi, pi]."""
    y = (x + math.pi) % (2 * math.pi) - math.pi
    return y

def deg(x: float) -> float:
    return 180.0 * x / math.pi

def complex_diag_phases(phases: List[float]) -> np.ndarray:
    """Diagonal phase matrix diag(exp(i*phase_k)). phases in radians."""
    return np.diag([np.exp(1j * p) for p in phases]).astype(np.complex128)

def cyclic_distance(a: int, b: int, cycle: int) -> int:
    """Distance on Z_cycle."""
    d = abs(a - b) % cycle
    return min(d, cycle - d)

def build_kappa_kernel_from_sites(
    sites: List[int],
    kappa: float,
    cycle: int = 360,
    kind: str = "power",
) -> np.ndarray:
    """
    Build NxN kernel K from flavor sites on a cycle:
      power: K_ij = kappa^{d(i,j)}
      exp  : K_ij = exp(-kappa * d(i,j))
    """
    N = len(sites)
    K = np.zeros((N, N), dtype=np.complex128)
    for i in range(N):
        for j in range(N):
            d = cyclic_distance(int(sites[i]), int(sites[j]), int(cycle))
            if kind == "power":
                K[i, j] = (kappa ** d)
            elif kind == "exp":
                K[i, j] = np.exp(-kappa * d)
            else:
                raise ValueError(f"Unknown kernel kind: {kind}")
    return K

def poly_in_kernel(K: np.ndarray, coeffs: List[complex]) -> np.ndarray:
    """
    Compute P(K) = c0 I + c1 K + c2 K^2 + ...
    """
    N = K.shape[0]
    out = np.zeros_like(K, dtype=np.complex128)
    Ki = np.eye(N, dtype=np.complex128)
    for c in coeffs:
        out += c * Ki
        Ki = Ki @ K
    return out

def make_yukawa(K: np.ndarray, sector_cfg: Dict[str, Any]) -> np.ndarray:
    """
    Sector texture:
      Y = P_L * (c0 I + c1 K + c2 K^2 + ...) * P_R†

    where P_L,P_R are diagonal phase matrices (or identity).
    """
    coeffs = sector_cfg.get("poly_coeffs", [1.0, 0.0, 0.0])
    coeffs_c = [complex(c[0], c[1]) if isinstance(c, list) else complex(c) for c in coeffs]

    phases_L = sector_cfg.get("phases_L", None)
    phases_R = sector_cfg.get("phases_R", None)

    P_L = complex_diag_phases(phases_L) if phases_L is not None else np.eye(K.shape[0], dtype=np.complex128)
    P_R = complex_diag_phases(phases_R) if phases_R is not None else np.eye(K.shape[0], dtype=np.complex128)

    core = poly_in_kernel(K, coeffs_c)
    Y = P_L @ core @ P_R.conj().T
    return Y

def svd_diagonalize(Y: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Y = U_L diag(s) U_R†  via SVD.
    Returns (singular_values_desc, U_L, U_R).
    """
    U, s, Vh = np.linalg.svd(Y)
    # Ensure descending order (numpy already does, but keep explicit)
    idx = np.argsort(s)[::-1]
    s = s[idx]
    U = U[:, idx]
    Vh = Vh[idx, :]
    UR = Vh.conj().T
    return s, U, UR

def pdg_angles_and_delta(U: np.ndarray) -> Dict[str, float]:
    """
    Extract PDG-like mixing angles (θ12,θ23,θ13) and δ from a unitary matrix U.
    Uses absolute values for angles; δ via a standard rephasing-invariant combination.

    Conventions:
      s13 = |U_{e3}|
      s12 = |U_{e2}| / c13
      s23 = |U_{μ3}| / c13
      δ   = arg(-U_{e1} U_{μ3} U*_{e3} U*_{μ1})
    """
    U = np.array(U, dtype=np.complex128)
    s13 = abs(U[0, 2])
    s13 = min(max(s13, 0.0), 1.0)
    c13 = math.sqrt(max(0.0, 1.0 - s13 * s13))

    if c13 < 1e-15:
        # pathological
        theta13 = math.asin(s13)
        return {"theta12": float("nan"), "theta23": float("nan"), "theta13": theta13, "delta": float("nan")}

    s12 = abs(U[0, 1]) / c13
    s23 = abs(U[1, 2]) / c13
    s12 = min(max(s12, 0.0), 1.0)
    s23 = min(max(s23, 0.0), 1.0)

    theta13 = math.asin(s13)
    theta12 = math.asin(s12)
    theta23 = math.asin(s23)

    # Rephasing-invariant δ estimate
    inv = -U[0, 0] * U[1, 2] * np.conj(U[0, 2]) * np.conj(U[1, 0])
    delta = wrap_pi(cmath.phase(inv))

    return {"theta12": theta12, "theta23": theta23, "theta13": theta13, "delta": delta}

def jarlskog(U: np.ndarray) -> float:
    """J = Im(U11 U22 U12* U21*) with 0-based indices (0,0)(1,1)(0,1)(1,0)."""
    return float(np.imag(U[0, 0] * U[1, 1] * np.conj(U[0, 1]) * np.conj(U[1, 0])))

def takagi_factorization(M: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Takagi factorization for (approximately) complex symmetric M:
      M = U diag(m) U^T with m >= 0.

    Practical approach:
      Use SVD: M = U Σ V†. For symmetric M, V ≈ U*.
      We symmetrize M slightly, then take U from SVD.
    Returns (masses_desc, U).
    """
    Ms = 0.5 * (M + M.T)  # enforce symmetry
    U, s, Vh = np.linalg.svd(Ms)
    # s are nonnegative singular values; for symmetric Ms, these are Takagi masses.
    idx = np.argsort(s)[::-1]
    s = s[idx]
    U = U[:, idx]
    # Fix arbitrary phases so that U^T Ms U is ~ diagonal real-positive
    D = U.T @ Ms @ U
    phase_fix = np.exp(-0.5j * np.angle(np.diag(D) + 1e-30))
    U = U @ np.diag(phase_fix)
    return s, U

def typeI_seesaw_mnu(Ynu: np.ndarray, MR: np.ndarray, v: float = 1.0) -> np.ndarray:
    """
    mν = - v^2 * Yν * MR^{-1} * Yν^T
    """
    MRi = np.linalg.inv(MR)
    return -(v * v) * (Ynu @ MRi @ Ynu.T)

# -----------------------------
# Default configuration
# -----------------------------

DEFAULT_CONFIG: Dict[str, Any] = {
    "cycle": 360,
    "sites": [1, 2, 5],          # 3 flavor sites on the 360-cycle
    "kappa": 0.24,
    "kernel_kind": "power",      # "power" or "exp"
    "v": 1.0,                    # set 174.0 if you want absolute GeV scaling
    "sectors": {
        # poly_coeffs can be real numbers or [re,im]
        # You can tune these to match your target ratios/angles.
        "Yu": {
            "poly_coeffs": [ [0.00,0.00], [1.00,0.00], [0.20,0.00] ],
            "phases_L": [0.00, 0.15, -0.10],
            "phases_R": [0.00, -0.20, 0.25],
        },
        "Yd": {
            "poly_coeffs": [ [0.00,0.00], [0.75,0.00], [0.35,0.00] ],
            "phases_L": [0.00, 0.05, 0.12],
            "phases_R": [0.00, -0.10, 0.08],
        },
        "Ye": {
            "poly_coeffs": [ [0.00,0.00], [0.80,0.00], [0.30,0.00] ],
            "phases_L": [0.00, -0.08, 0.18],
            "phases_R": [0.00, 0.12, -0.06],
        },
        "Ynu": {
            "poly_coeffs": [ [0.00,0.00], [1.00,0.00], [0.10,0.00] ],
            "phases_L": [0.00, 0.30, -0.22],
            "phases_R": [0.00, -0.18, 0.10],
        },
        # Optional Majorana: if omitted, PMNS is computed from Dirac-like SVD (less physical).
        "MR": {
            # symmetric MR in same basis; simplest: diagonal (real-positive)
            "matrix": [
                [1.0, 0.0, 0.0],
                [0.0, 0.2, 0.0],
                [0.0, 0.0, 0.05]
            ],
            "scale": 1.0  # overall scale factor (e.g. 2e14 if you want eV/GeV bookkeeping elsewhere)
        }
    }
}

# -----------------------------
# I/O helpers
# -----------------------------

def load_config(path: Optional[str]) -> Dict[str, Any]:
    if path is None:
        return json.loads(json.dumps(DEFAULT_CONFIG))
    with open(path, "r", encoding="utf-8") as f:
        cfg = json.load(f)
    return cfg

def mat_from_cfg(Mcfg: Any) -> np.ndarray:
    M = np.array(Mcfg, dtype=np.complex128)
    return M

# -----------------------------
# Main pipeline
# -----------------------------

def run_pipeline(cfg: Dict[str, Any]) -> Dict[str, Any]:
    cycle = int(cfg.get("cycle", 360))
    sites = [int(x) for x in cfg.get("sites", [1,2,5])]
    kappa = float(cfg.get("kappa", 0.24))
    kind  = str(cfg.get("kernel_kind", "power"))
    v     = float(cfg.get("v", 1.0))

    if len(sites) != 3:
        raise ValueError("This script assumes 3 flavors (3 sites). Provide exactly 3 sites.")

    K = build_kappa_kernel_from_sites(sites=sites, kappa=kappa, cycle=cycle, kind=kind)

    sectors = cfg.get("sectors", {})
    Yu = make_yukawa(K, sectors.get("Yu", {}))
    Yd = make_yukawa(K, sectors.get("Yd", {}))
    Ye = make_yukawa(K, sectors.get("Ye", {}))
    Ynu = make_yukawa(K, sectors.get("Ynu", {}))

    su, UuL, _ = svd_diagonalize(Yu)
    sd, UdL, _ = svd_diagonalize(Yd)
    se, UeL, _ = svd_diagonalize(Ye)

    Vckm = UuL.conj().T @ UdL
    ckm = pdg_angles_and_delta(Vckm)

    # Neutrinos
    MR_cfg = sectors.get("MR", None)
    if MR_cfg is not None and "matrix" in MR_cfg:
        MR = mat_from_cfg(MR_cfg["matrix"])
        MR = 0.5 * (MR + MR.T)  # enforce symmetry
        MR *= float(MR_cfg.get("scale", 1.0))
        mnu = typeI_seesaw_mnu(Ynu, MR, v=v)
        mnu_masses, Unu = takagi_factorization(mnu)
        pmns = UeL.conj().T @ Unu
        pmns_angles = pdg_angles_and_delta(pmns)

        # Mass-squared diffs (ordering by descending mass here; interpret carefully)
        m1, m2, m3 = mnu_masses[::-1]  # ascending (m1<=m2<=m3) by convention
        dm21 = float(m2*m2 - m1*m1)
        dm31 = float(m3*m3 - m1*m1)
    else:
        # fallback: treat Ynu like Dirac and use SVD-left unitary (not a Majorana PMNS)
        snu, UnuL, _ = svd_diagonalize(Ynu)
        mnu_masses = snu
        pmns = UeL.conj().T @ UnuL
        pmns_angles = pdg_angles_and_delta(pmns)
        dm21 = float("nan"); dm31 = float("nan")

    out: Dict[str, Any] = {
        "inputs": {"cycle": cycle, "sites": sites, "kappa": kappa, "kernel_kind": kind, "v": v},
        "kernel_K": K.tolist(),
        "singular_values": {
            "Yu": su.tolist(),
            "Yd": sd.tolist(),
            "Ye": se.tolist(),
            "Ynu": mnu_masses.tolist(),
        },
        "mass_ratios": {
            "up":   [float(su[0]/su[0]), float(su[1]/su[0]), float(su[2]/su[0])],
            "down": [float(sd[0]/sd[0]), float(sd[1]/sd[0]), float(sd[2]/sd[0])],
            "lep":  [float(se[0]/se[0]), float(se[1]/se[0]), float(se[2]/se[0])],
        },
        "CKM": {
            "V": Vckm.tolist(),
            "angles_deg": {k: deg(vv) for k, vv in ckm.items() if k != "delta"} | {"delta_deg": deg(ckm["delta"])},
            "J": jarlskog(Vckm),
        },
        "PMNS": {
            "U": pmns.tolist(),
            "angles_deg": {k: deg(vv) for k, vv in pmns_angles.items() if k != "delta"} | {"delta_deg": deg(pmns_angles["delta"])},
            "J": jarlskog(pmns),
            "dm21": dm21,
            "dm31": dm31,
        },
    }
    return out

def pretty_print(report: Dict[str, Any]) -> None:
    inp = report["inputs"]
    print("\nFlavor pipeline")
    print(f"  cycle={inp['cycle']}  sites={inp['sites']}  kappa={inp['kappa']}  kind={inp['kernel_kind']}\n")

    su = report["singular_values"]["Yu"]
    sd = report["singular_values"]["Yd"]
    se = report["singular_values"]["Ye"]
    print("Yukawa singular values (descending):")
    print(f"  Yu: {np.array(su)}")
    print(f"  Yd: {np.array(sd)}")
    print(f"  Ye: {np.array(se)}\n")

    print("Mass ratios (largest normalized to 1):")
    print(f"  up   (t,c,u): {report['mass_ratios']['up']}")
    print(f"  down (b,s,d): {report['mass_ratios']['down']}")
    print(f"  lep  (τ,μ,e): {report['mass_ratios']['lep']}\n")

    ckm = report["CKM"]["angles_deg"]
    print("CKM (PDG-like):")
    print(f"  θ12={ckm['theta12']:.3f}°  θ23={ckm['theta23']:.3f}°  θ13={ckm['theta13']:.3f}°  δ={ckm['delta_deg']:.3f}°")
    print(f"  J={report['CKM']['J']:.6e}\n")

    pmns = report["PMNS"]["angles_deg"]
    print("PMNS (PDG-like):")
    print(f"  θ12={pmns['theta12']:.3f}°  θ23={pmns['theta23']:.3f}°  θ13={pmns['theta13']:.3f}°  δ={pmns['delta_deg']:.3f}°")
    print(f"  J={report['PMNS']['J']:.6e}")
    print(f"  Δm^2_21={report['PMNS']['dm21']:.6e}  Δm^2_31={report['PMNS']['dm31']:.6e}\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=str, default=None, help="Path to JSON config (optional).")
    ap.add_argument("--save", type=str, default=None, help="Save report JSON to this file.")
    args = ap.parse_args()

    cfg = load_config(args.config)
    report = run_pipeline(cfg)
    pretty_print(report)

    if args.save:
        with open(args.save, "w", encoding="utf-8") as f:
            json.dump(report, f, indent=2)
        print(f"Saved report to: {args.save}")

if __name__ == "__main__":
    main()

"""
/Users/chazzromeo/rAI/Flavor/.venv/bin/python /Users/chazzromeo/rAI/Flavor/alignment_v5.py 

Flavor pipeline
  cycle=360  sites=[1, 2, 5]  kappa=0.24  kind=power

Yukawa singular values (descending):
  Yu: [1.54843399 1.19946681 0.87522005]
  Yd: [1.46914857 1.09944779 0.77186512]
  Ye: [1.45422335 1.09946682 0.78099109]

Mass ratios (largest normalized to 1):
  up   (t,c,u): [1.0, 0.7746321900619537, 0.5652291635309713]
  down (b,s,d): [1.0, 0.748357115655774, 0.5253826127401872]
  lep  (τ,μ,e): [1.0, 0.7560508684955204, 0.5370503061808038]

CKM (PDG-like):
  θ12=0.688°  θ23=0.334°  θ13=2.884°  δ=92.228°
  J=3.508112e-06

PMNS (PDG-like):
  θ12=85.641°  θ23=3.123°  θ13=27.587°  δ=26.500°
  J=6.704310e-04
  Δm^2_21=4.624190e+01  Δm^2_31=5.854426e+02

"""