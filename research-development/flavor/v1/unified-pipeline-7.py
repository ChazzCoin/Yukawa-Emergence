import numpy as np
import math

# Fix the best phase configuration found so far
PHASE_U = (0, 2)
PHASE_D = (0, 3)
PHASE_E = (0, 10)
PHASE_NU = (0, 25)
NOISE_LEVEL = 0.05

# Rebind run_alignment_high_scale again explicitly to lock these in
def run_alignment_high_scale(
    seed=0,
    exponents_u=(4, 2, 0),
    exponents_d=(3, 1, 0),
    exponents_e=(4, 2, 0),
    exponents_nu=(1, 0, 0),
    triad_ks=(1, 2, 3),
    use_site_hierarchy=True,
):
    K = build_alignment_kernel(eps, N=9, d_star=7)
    Yu0, Yd0, Ye0, Ynu0, M0 = generate_aligned_proto_matrices(
        seed,
        use_site_hierarchy=use_site_hierarchy,
        exponents_u=exponents_u,
        exponents_d=exponents_d,
        exponents_e=exponents_e,
        exponents_nu=exponents_nu,
        phase_u=PHASE_U,
        phase_d=PHASE_D,
        phase_e=PHASE_E,
        phase_nu=PHASE_NU,
        noise_level=NOISE_LEVEL,
    )
    Yu9, Yd9, Ye9, Ynu9, M9 = align_all(K, Yu0, Yd0, Ye0, Ynu0, M0)
    Yu_eff, Yd_eff, Ye_eff, Ynu_eff = reduce_all_dirac(Yu9, Yd9, Ye9, Ynu9)
    M_R = build_M_R_triadic(M9, Lambda_Maj, ks=triad_ks)
    mnu = seesaw_light_neutrinos(Ynu_eff, M_R, v_HIGGS)
    return Yu_eff, Yd_eff, Ye_eff, mnu

# Use the existing run_full_pipeline_with_RGE_and_rescaling and chi2_from_res infrastructure
def run_pipeline_for_seed_fixed_phases(seed, mu_high=1e14):
    base = run_full_pipeline_with_RGE_and_rescaling(seed=seed, mu_high=mu_high)
    mu_vals = base["mu"]
    md_vals = base["md"]
    me_vals = base["me"]
    mnu_vals = base["mnu"]
    Vckm = base["Vckm"]
    Vpmns = base["Vpmns"]

    th12_q, th23_q, th13_q, delta_q = extract_angles_and_phase(Vckm)
    th12_l, th23_l, th13_l, delta_l = extract_angles_and_phase(Vpmns)

    dm2_21_GeV2, dm2_31_GeV2 = neutrino_splittings(mnu_vals)
    dm2_21_eV2 = dm2_21_GeV2 * 1e18
    dm2_31_eV2 = dm2_31_GeV2 * 1e18

    res = {
        "mu": mu_vals,
        "md": md_vals,
        "me": me_vals,
        "mnu": mnu_vals,
        "th_q": (th12_q, th23_q, th13_q, delta_q),
        "th_l": (th12_l, th23_l, th13_l, delta_l),
        "dm2_eV2": (dm2_21_eV2, dm2_31_eV2),
    }
    res["chi2"] = chi2_from_res(res)
    return res

if __name__ == "__main__":
    # Analyze the best seed (seed = 10) in detail: observables, pulls, etc.
    best_seed = 10
    best_res = run_pipeline_for_seed_fixed_phases(best_seed)
    obs_vec = make_observables(best_res)
    pulls = (obs_vec - x_exp) / sigma

    names = [
        "m_c/m_t",
        "m_u/m_t",
        "m_s/m_b",
        "m_d/m_b",
        "m_mu/m_tau",
        "m_e/m_tau",
        "theta12^q (CKM)",
        "theta23^q (CKM)",
        "theta13^q (CKM)",
        "theta12^l (PMNS)",
        "theta23^l (PMNS)",
        "theta13^l (PMNS)",
        "Δm21^2 (eV^2)",
        "Δm31^2 (eV^2)",
    ]

    rows = []
    for i, name in enumerate(names):
        rows.append({
            "observable": name,
            "theory": obs_vec[i],
            "exp": x_exp[i],
            "sigma": sigma[i],
            "pull_sigma": pulls[i]
        })

    rows_sorted = sorted(rows, key=lambda r: abs(r["pull_sigma"]), reverse=True)

    # Also extract masses and angles for this seed
    mu_vals, md_vals, me_vals, mnu_vals = best_res["mu"], best_res["md"], best_res["me"], best_res["mnu"]
    th12_q, th23_q, th13_q, delta_q = best_res["th_q"]
    th12_l, th23_l, th13_l, delta_l = best_res["th_l"]

    best_seed, best_res["chi2"], rows_sorted[:6], {
        "mu": mu_vals,
        "md": md_vals,
        "me": me_vals,
        "mnu": mnu_vals,
        "th_q_deg": [math.degrees(th12_q), math.degrees(th23_q), math.degrees(th13_q), math.degrees(delta_q)],
        "th_l_deg": [math.degrees(th12_l), math.degrees(th23_l), math.degrees(th13_l), math.degrees(delta_l)],
    }

if __name__ == "__main__":
    # Scan seeds 0..10 for this fixed phase configuration
    chi2_list = []
    res_list = []
    for seed in range(11):
        res = run_pipeline_for_seed_fixed_phases(seed)
        chi2_list.append(res["chi2"])
        res_list.append(res)
        print(f"seed {seed:2d}: chi2 = {res['chi2']:.3g}")

    chi2_arr = np.array(chi2_list)
    best_idx = int(np.argmin(chi2_arr))
    best_seed = best_idx
    best_res = res_list[best_idx]

    best_seed, chi2_arr[best_idx], chi2_arr.min(), chi2_arr.max()
