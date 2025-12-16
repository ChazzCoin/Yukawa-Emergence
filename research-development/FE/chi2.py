# FE/chi2.py

import numpy as np


def mixing_angles_from_matrix(U):
    """
    Extract approximate mixing angles from a 3x3 unitary.
    Fully emergent: no Euler angle assumptions.
    """
    U = np.asarray(U)
    t12 = np.abs(U[0,1]) / max(1e-15, np.abs(U[0,0]))
    θ12 = np.arctan(t12)

    θ13 = np.abs(U[0,2])

    θ23 = np.abs(U[1,2]) / max(1e-15, np.abs(U[2,2]))
    θ23 = np.arctan(θ23)

    return θ12, θ23, θ13


def reduced_mass_ratios(masses):
    """
    Given array [m3, m2, m1], return {m2/m3, m1/m3}.
    """
    m3, m2, m1 = masses
    if m3 < 1e-15:
        return 0.0, 0.0
    return m2/m3, m1/m3


def chi2_term(model, target, sigma):
    """
    Standard chi^2 = ((model - target)/sigma)^2.
    """
    if sigma <= 0:
        sigma = max(1e-6, 0.1 * abs(target))
    return ((model - target) / sigma)**2


def compute_global_chi2(spectra, mixings):
    """
    spectra: { 'u': [...], 'd': [...], 'e': [...], 'nu': [...] }
    mixings: { 'CKM': M, 'PMNS': M }

    Returns:
        chi2_total, contributions (dict)
    """
    contributions = {}

    # --- 1) MASS RATIOS ---
    # SM reference ratios (dimensionless, approximate)
    SM_ratios = {
        "u_m2/m3":   0.005,
        "u_m1/m3":   0.00004,
        "d_m2/m3":   0.05,
        "d_m1/m3":   0.003,
        "e_m2/m3":   0.006,
        "e_m1/m3":   0.0000048,
        "nu_m2/m3":  0.17,
        "nu_m1/m3":  0.01
    }

    sig_ratio = 0.20  # 20% tolerance (dimensionless)

    for sector in ["u", "d", "e", "nu"]:
        m2m3, m1m3 = reduced_mass_ratios(spectra[sector])

        contributions[f"{sector}_m2/m3"] = chi2_term(
            m2m3, SM_ratios[f"{sector}_m2/m3"], sig_ratio * SM_ratios[f"{sector}_m2/m3"]
        )
        contributions[f"{sector}_m1/m3"] = chi2_term(
            m1m3, SM_ratios[f"{sector}_m1/m3"], sig_ratio * SM_ratios[f"{sector}_m1/m3"]
        )

    # --- 2) MIXING ANGLES ---
    # SM central values
    SM_mix = {
        "CKM_θ12": 0.226,
        "CKM_θ23": 0.041,
        "CKM_θ13": 0.0036,
        "PMNS_θ12": 0.59,
        "PMNS_θ23": 0.78,
        "PMNS_θ13": 0.15,
    }

    sig_ckm = 0.20  # 20% fractional tolerance
    sig_pmns = 0.20

    # CKM
    θ12, θ23, θ13 = mixing_angles_from_matrix(mixings["CKM"])
    contributions["CKM_θ12"] = chi2_term(θ12, SM_mix["CKM_θ12"], sig_ckm*SM_mix["CKM_θ12"])
    contributions["CKM_θ23"] = chi2_term(θ23, SM_mix["CKM_θ23"], sig_ckm*SM_mix["CKM_θ23"])
    contributions["CKM_θ13"] = chi2_term(θ13, SM_mix["CKM_θ13"], sig_ckm*SM_mix["CKM_θ13"])

    # PMNS
    θ12, θ23, θ13 = mixing_angles_from_matrix(mixings["PMNS"])
    contributions["PMNS_θ12"] = chi2_term(θ12, SM_mix["PMNS_θ12"], sig_pmns*SM_mix["PMNS_θ12"])
    contributions["PMNS_θ23"] = chi2_term(θ23, SM_mix["PMNS_θ23"], sig_pmns*SM_mix["PMNS_θ23"])
    contributions["PMNS_θ13"] = chi2_term(θ13, SM_mix["PMNS_θ13"], sig_pmns*SM_mix["PMNS_θ13"])

    # --- 3) TOTAL ---
    chi2_total = float(sum(contributions.values()))
    return chi2_total, contributions