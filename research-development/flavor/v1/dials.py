from dataclasses import dataclass, field
from typing import Tuple, Dict
import numpy as np


# -----------------------------
# Small helper dataclasses
# -----------------------------

@dataclass
class SectorPhases:
    """Phase seed and N_eff for one sector."""
    n0: int
    delta: int
    N_eff: int


@dataclass
class FNCharges:
    """Froggatt–Nielsen left-handed charges + small parameter."""
    QL_u: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    QL_d: Tuple[float, float, float] = (1.0, 0.5, 0.0)
    QL_e: Tuple[float, float, float] = (0.5, 0.0, 0.0)
    QL_nu: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    eps_L: float = 0.3


@dataclass
class RGESettings:
    """Toy RGE settings (for stub_rge_run)."""
    alpha: float = 0.1
    mu_high: float = 1e14
    mu_EW: float = 173.0


@dataclass
class SeesawSettings:
    """Overall scales for the seesaw."""
    v: float = 174.0          # EW vev (GeV)
    Lambda_Maj: float = 7e13  # heavy Majorana scale (GeV)


@dataclass
class CMASettings:
    """Optimizer hyperparameters."""
    num_restarts: int = 8
    sigma_init: float = 0.3
    popsize: int = 20
    maxiter: int = 200
    seed: int = 9


# -----------------------------
# Main config dataclass
# -----------------------------

@dataclass
class FlavorModelConfig:
    """
    One-stop shop for all tunable knobs in the 9-site flavor model.
    Everything that is NOT purely structural should live here.
    """

    # ---- geometric coherence & kernels ----
    kappa: float = 0.242          # your discovered coherence scale
    eps_align: float = 4.13       # current "effective" EPS (decoupled from 1/kappa)
    forbidden_d: int = 2          # forbidden distance on the 9-site ring

    # sector-dependent decoherence exponents γ_x
    gamma_u: float = 0.00
    gamma_d: float = 0.03
    gamma_e: float = 0.05
    gamma_nu: float = 0.08
    gamma_maj: float = 0.08

    # ---- baseline triadic exponents per sector (FN-like hierarchy priors) ----
    exp_up_base:  Tuple[float, float, float] = (4.0, 2.0, 0.0)
    exp_down_base: Tuple[float, float, float] = (3.0, 2.0, 0.0)
    exp_lep_base:  Tuple[float, float, float] = (3.0, 2.0, 0.0)
    exp_nu_base:   Tuple[float, float, float] = (1.0, 0.0, 0.0)

    # ---- harmonic phase seeds per sector ----
    # (n0, delta, N_eff) for each sector
    phases_u: SectorPhases = field(default_factory=lambda: SectorPhases(0, 6,   360))
    phases_d: SectorPhases = field(default_factory=lambda: SectorPhases(0, 3,   180))
    phases_e: SectorPhases = field(default_factory=lambda: SectorPhases(0, 8,   120))
    phases_nu: SectorPhases = field(default_factory=lambda: SectorPhases(0, 24,  90))

    # ---- FN charges & eps_L ----
    fn: FNCharges = field(default_factory=FNCharges)

    # ---- Seesaw & RGE ----
    seesaw: SeesawSettings = field(default_factory=SeesawSettings)
    rge:    RGESettings    = field(default_factory=RGESettings)

    # ---- projection & leakages (default starting values) ----
    a_tq_init: float = 0.0     # initial Cabibbo twist amplitude
    phi_tq_init: float = 0.0   # initial Cabibbo twist phase

    proj_eps_03_init: float = 0.0
    proj_eps_30_init: float = 0.0
    proj_eps_36_init: float = 0.0

    # you also had eps12, eps21 as small regularized knobs:
    eps12_init: float = 0.0
    eps21_init: float = 0.0

    # ---- penalty weights ----
    lambda_geom: float = 0.05
    lambda_exp:  float = 0.10
    lambda_proj: float = 1.00

    # ---- experimental targets & errors ----
    exp_targets: Dict[str, float] = field(default_factory=lambda: {
        'm_c/m_t':      7e-3,
        'm_u/m_t':      1e-5,
        'm_s/m_b':      2e-2,
        'm_d/m_b':      1e-3,
        'm_mu/m_tau':   6e-2,
        'm_e/m_tau':    3e-4,
        'theta12_q':    0.226,
        'theta23_q':    0.041,
        'theta13_q':    0.0035,
        'theta12_l':    0.59,
        'theta23_l':    0.84,
        'theta13_l':    0.15,
        'Delta m2_21':  7.4e-5,
        'Delta m2_31':  2.5e-3,
    })

    sigma_targets: Dict[str, float] = field(default_factory=dict)

    # ---- optimizer settings ----
    cma: CMASettings = field(default_factory=CMASettings)

    # ---- misc / seeds ----
    M0_seed: int = 9   # seed for proto-Majorana generation

    def __post_init__(self):
        """
        Fill in derived pieces:
          - default sigma_targets if not provided
        """
        if not self.sigma_targets:
            self.sigma_targets = {
                'm_c/m_t':     0.5 * self.exp_targets['m_c/m_t'],
                'm_u/m_t':     0.5 * self.exp_targets['m_u/m_t'],
                'm_s/m_b':     0.5 * self.exp_targets['m_s/m_b'],
                'm_d/m_b':     0.5 * self.exp_targets['m_d/m_b'],
                'm_mu/m_tau':  0.5 * self.exp_targets['m_mu/m_tau'],
                'm_e/m_tau':   0.5 * self.exp_targets['m_e/m_tau'],
                'theta12_q':   0.1 * self.exp_targets['theta12_q'],
                'theta23_q':   0.1 * self.exp_targets['theta23_q'],
                'theta13_q':   0.1 * self.exp_targets['theta13_q'],
                'theta12_l':   0.1 * self.exp_targets['theta12_l'],
                'theta23_l':   0.1 * self.exp_targets['theta23_l'],
                'theta13_l':   0.1 * self.exp_targets['theta13_l'],
                'Delta m2_21': 0.3 * self.exp_targets['Delta m2_21'],
                'Delta m2_31': 0.3 * self.exp_targets['Delta m2_31'],
            }

    # --- optional helpers ---

    @property
    def exp_up_base_arr(self) -> np.ndarray:
        return np.array(self.exp_up_base, dtype=float)

    @property
    def exp_down_base_arr(self) -> np.ndarray:
        return np.array(self.exp_down_base, dtype=float)

    @property
    def exp_lep_base_arr(self) -> np.ndarray:
        return np.array(self.exp_lep_base, dtype=float)

    @property
    def exp_nu_base_arr(self) -> np.ndarray:
        return np.array(self.exp_nu_base, dtype=float)

    def as_gamma_dict(self) -> Dict[str, float]:
        """Convenient mapping for building sector kernels."""
        return {
            "u":   self.gamma_u,
            "d":   self.gamma_d,
            "e":   self.gamma_e,
            "nu":  self.gamma_nu,
            "maj": self.gamma_maj,
        }