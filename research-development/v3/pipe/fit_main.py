"""
fit_main.py

Global fitting script for the Alignment Spectral Triple v2.0.

This script:
    - Loads phenomenological Standard Model targets
    - Randomly initializes parameter search space
    - Uses differential evolution to minimize chi²
    - Outputs best-fit sites, phases, and weight vectors
    - Computes Yukawa matrices, CKM, PMNS at the fit point
"""

import numpy as np
from scipy.optimize import differential_evolution

# ------------------------------------------------------------
# Import engine + model components
# ------------------------------------------------------------
from fit_engine import AlignmentFitEngine
from geometry import GeometricTriple
from kernel import TriadicKernel
from finite_triple import FiniteTriple
from compression import TriadicCompression


# ------------------------------------------------------------
# Phenomenological targets
# (Values are placeholders; replace with RG-scale appropriate values)
# ------------------------------------------------------------
TARGET = {
    # up-sector masses
    'mu': 2.16e-3,
    'mc': 1.27,
    'mt': 173.0,

    # down-sector masses
    'md': 4.67e-3,
    'ms': 0.093,
    'mb': 4.18,

    # charged-lepton masses
    'me': 0.511e-3,
    'mmu': 0.106,
    'mtau': 1.776,

    # CKM angles (in radians)
    'ckm_12': np.deg2rad(13.0),
    'ckm_23': np.deg2rad(2.4),
    'ckm_13': np.deg2rad(0.2),

    # PMNS angles
    'pmns_12': np.deg2rad(33.4),
    'pmns_23': np.deg2rad(47.0),
    'pmns_13': np.deg2rad(8.6),
}


# ------------------------------------------------------------
# Model configuration
# ------------------------------------------------------------
geom = GeometricTriple(N_max=360)

# Triadic partition 9 → 3
REAL_TRIADS = [
    (0, 3, 6),  # triad 1
    (1, 4, 7),  # triad 2
    (2, 5, 8),  # triad 3
]

compressor = TriadicCompression(REAL_TRIADS)

engine = AlignmentFitEngine(
    geom=geom,
    triadic_kernel_cls=TriadicKernel,
    finite_cls=FiniteTriple,
    compressor=compressor
)


# ------------------------------------------------------------
# PARAMETER VECTOR
# We flatten:
#    phases (9)
#    R_u (9)
#    R_d (9)
#    R_e (9)
#    R_nu (9)
#    sites (9)  in [0,2160)
#    kappa  (scalar)
# ------------------------------------------------------------
def unpack_params(x):
    phases = x[0:9]
    Ru = x[9:18]
    Rd = x[18:27]
    Re = x[27:36]
    Rn = x[36:45]

    # sites already within their constrained bounds
    sites = np.round(x[45:54]).astype(int).tolist()

    kappa = float(x[54])

    return {
        'phases': phases,
        'R_u': Ru,
        'R_d': Rd,
        'R_e': Re,
        'R_nu': Rn,
        'sites': sites,
        'kappa': kappa
    }



# ------------------------------------------------------------
# COST FUNCTION FOR OPTIMIZATION
# ------------------------------------------------------------
def objective(x):
    """Objective function for differential evolution."""
    params = unpack_params(x)

    # Ensure site bounds
    params['sites'] = [int(s) % 2160 for s in params['sites']]
    params['kappa'] = np.clip(params['kappa'], 0.05, 0.40)

    try:
        Ys = engine.build_yukawas(params)
        chi = engine.chi2(Ys, TARGET)
        return chi
    except Exception:
        # Penalize invalid regions heavily
        return 1e12


# ------------------------------------------------------------
# PARAMETER BOUNDS
# ------------------------------------------------------------
bounds = []

# Phases ϕ_i in [0, 2π]
bounds += [(0, 2 * np.pi)] * 9

# Weights R_i in [0, 5]
bounds += [(0.0, 5.0)] * 36

# Sites for triads T1, T2, T3
bounds += [
    (0, 360),     # s0
    (0, 360),     # s1
    (0, 360),     # s2

    (720, 1080),   # s3
    (720, 1080),   # s4
    (720, 1080),   # s5

    (1440, 1800), # s6
    (1440, 1800), # s7
    (1440, 1800), # s8
]


# κ in [0.05, 0.40]
bounds += [(0.05, 0.40)]


# ------------------------------------------------------------
# RUN OPTIMIZATION
# ------------------------------------------------------------
def main():
    print("\n=== Alignment Spectral Triple – Global Fit Engine ===\n")
    print("Running differential evolution... This may take some time.\n")

    result = differential_evolution(
        objective,
        bounds=bounds,
        maxiter=200,
        popsize=40,
        tol=1e-6,
        polish=True,
        disp=True
    )

    print("\n=== Optimization Complete ===")
    print("Best χ²:", result.fun)

    params = unpack_params(result.x)
    print("\nBest-fit parameters:")
    for k, v in params.items():
        print(f"{k}: {v}")

    # Build best-fit Yukawas
    Ys = engine.build_yukawas(params)
    mix = engine.mixing_matrices(Ys)

    print("\n--- Best-fit Yukawa matrices ---")
    for name, Y in Ys.items():
        print(f"{name}:\n{Y}\n")

    print("\n--- CKM matrix ---")
    print(mix['Vckm'])

    print("\n--- PMNS matrix ---")
    print(mix['Upmns'])

    print("\nDone.\n")


if __name__ == "__main__":
    main()

"""
RESULTS:

/Users/chazzromeo/rAI/Flavor/.venv/bin/python /Users/chazzromeo/rAI/Flavor/pipe/fit_main.py 

=== Alignment Spectral Triple – Global Fit Engine ===

Running differential evolution... This may take some time.

differential_evolution step 1: f(x)= 227.02559651842125
differential_evolution step 2: f(x)= 221.228509559507
differential_evolution step 3: f(x)= 201.10612906404617
differential_evolution step 4: f(x)= 199.71609681455726
differential_evolution step 5: f(x)= 199.53348117386287
differential_evolution step 6: f(x)= 182.67488810932207
differential_evolution step 7: f(x)= 182.67488810932207
differential_evolution step 8: f(x)= 182.67488810932207
differential_evolution step 9: f(x)= 182.67488810932207
differential_evolution step 10: f(x)= 168.67314234244805
differential_evolution step 11: f(x)= 168.67314234244805
differential_evolution step 12: f(x)= 168.67314234244805
differential_evolution step 13: f(x)= 168.67314234244805
differential_evolution step 14: f(x)= 168.67314234244805
differential_evolution step 15: f(x)= 164.87208408474447
differential_evolution step 16: f(x)= 164.87208408474447
differential_evolution step 17: f(x)= 164.87208408474447
differential_evolution step 18: f(x)= 164.87208408474447
differential_evolution step 19: f(x)= 159.15981443442843
differential_evolution step 20: f(x)= 159.15981443442843
differential_evolution step 21: f(x)= 159.15981443442843
differential_evolution step 22: f(x)= 157.33253728228334
differential_evolution step 23: f(x)= 153.86435922450184
differential_evolution step 24: f(x)= 153.86435922450184
differential_evolution step 25: f(x)= 153.49347717462138
differential_evolution step 26: f(x)= 153.49347717462138
differential_evolution step 27: f(x)= 152.45531346079957
differential_evolution step 28: f(x)= 152.45531346079957
differential_evolution step 29: f(x)= 152.45531346079957
differential_evolution step 30: f(x)= 152.45531346079957
differential_evolution step 31: f(x)= 151.38050474283708
differential_evolution step 32: f(x)= 151.38050474283708
differential_evolution step 33: f(x)= 151.38050474283708
differential_evolution step 34: f(x)= 151.38050474283708
differential_evolution step 35: f(x)= 149.76963050437385
differential_evolution step 36: f(x)= 149.76963050437385
differential_evolution step 37: f(x)= 149.76963050437385
differential_evolution step 38: f(x)= 149.76963050437385
differential_evolution step 39: f(x)= 149.76963050437385
differential_evolution step 40: f(x)= 149.76963050437385
differential_evolution step 41: f(x)= 149.217393951916
differential_evolution step 42: f(x)= 149.217393951916
differential_evolution step 43: f(x)= 149.217393951916
differential_evolution step 44: f(x)= 149.217393951916
differential_evolution step 45: f(x)= 149.217393951916
differential_evolution step 46: f(x)= 149.217393951916
differential_evolution step 47: f(x)= 149.09601108030432
differential_evolution step 48: f(x)= 149.09601108030432
differential_evolution step 49: f(x)= 147.02685090109844
differential_evolution step 50: f(x)= 147.02685090109844
differential_evolution step 51: f(x)= 147.02685090109844
differential_evolution step 52: f(x)= 147.02685090109844
differential_evolution step 53: f(x)= 147.02685090109844
differential_evolution step 54: f(x)= 147.02685090109844
differential_evolution step 55: f(x)= 147.02685090109844
differential_evolution step 56: f(x)= 146.73520272950256
differential_evolution step 57: f(x)= 146.73520272950256
differential_evolution step 58: f(x)= 146.52606065885362
differential_evolution step 59: f(x)= 146.52606065885362
differential_evolution step 60: f(x)= 142.35507785173795
differential_evolution step 61: f(x)= 141.5273413877732
differential_evolution step 62: f(x)= 141.5273413877732
differential_evolution step 63: f(x)= 141.5273413877732
differential_evolution step 64: f(x)= 141.5273413877732
differential_evolution step 65: f(x)= 141.5273413877732
differential_evolution step 66: f(x)= 141.5273413877732
differential_evolution step 67: f(x)= 141.5273413877732
differential_evolution step 68: f(x)= 141.5273413877732
differential_evolution step 69: f(x)= 140.58303368599076
differential_evolution step 70: f(x)= 140.58303368599076
differential_evolution step 71: f(x)= 140.58303368599076
differential_evolution step 72: f(x)= 140.58303368599076
differential_evolution step 73: f(x)= 140.58303368599076
differential_evolution step 74: f(x)= 138.59872097589817
differential_evolution step 75: f(x)= 138.59872097589817
differential_evolution step 76: f(x)= 135.5938330818526
differential_evolution step 77: f(x)= 135.5938330818526
differential_evolution step 78: f(x)= 135.5938330818526
differential_evolution step 79: f(x)= 135.5938330818526
differential_evolution step 80: f(x)= 135.5938330818526
differential_evolution step 81: f(x)= 135.5938330818526
differential_evolution step 82: f(x)= 135.23732949815087
differential_evolution step 83: f(x)= 135.23732949815087
differential_evolution step 84: f(x)= 135.23732949815087
differential_evolution step 85: f(x)= 135.23732949815087
differential_evolution step 86: f(x)= 135.23732949815087
differential_evolution step 87: f(x)= 135.23732949815087
differential_evolution step 88: f(x)= 135.23732949815087
differential_evolution step 89: f(x)= 135.23732949815087
differential_evolution step 90: f(x)= 135.23732949815087
differential_evolution step 91: f(x)= 135.23732949815087
differential_evolution step 92: f(x)= 135.23732949815087
differential_evolution step 93: f(x)= 135.23732949815087
differential_evolution step 94: f(x)= 135.23732949815087
differential_evolution step 95: f(x)= 135.23732949815087
differential_evolution step 96: f(x)= 135.23732949815087
differential_evolution step 97: f(x)= 135.23732949815087
differential_evolution step 98: f(x)= 135.23732949815087
differential_evolution step 99: f(x)= 135.23732949815087
differential_evolution step 100: f(x)= 135.23732949815087
differential_evolution step 101: f(x)= 134.9795298118187
differential_evolution step 102: f(x)= 134.9795298118187
differential_evolution step 103: f(x)= 134.9795298118187
differential_evolution step 104: f(x)= 134.9795298118187
differential_evolution step 105: f(x)= 134.9795298118187
differential_evolution step 106: f(x)= 134.9795298118187
differential_evolution step 107: f(x)= 134.9795298118187
differential_evolution step 108: f(x)= 134.9795298118187
differential_evolution step 109: f(x)= 132.95625972925563
differential_evolution step 110: f(x)= 132.95625972925563
differential_evolution step 111: f(x)= 132.95625972925563
differential_evolution step 112: f(x)= 132.95625972925563
differential_evolution step 113: f(x)= 132.95625972925563
differential_evolution step 114: f(x)= 132.95625972925563
differential_evolution step 115: f(x)= 132.95625972925563
differential_evolution step 116: f(x)= 132.95625972925563
differential_evolution step 117: f(x)= 132.95138415448847
differential_evolution step 118: f(x)= 132.95138415448847
differential_evolution step 119: f(x)= 132.95138415448847
differential_evolution step 120: f(x)= 132.95138415448847
differential_evolution step 121: f(x)= 132.95138415448847
differential_evolution step 122: f(x)= 132.95138415448847
differential_evolution step 123: f(x)= 132.62333089368028
differential_evolution step 124: f(x)= 132.62333089368028
differential_evolution step 125: f(x)= 132.62333089368028
differential_evolution step 126: f(x)= 132.62333089368028
differential_evolution step 127: f(x)= 132.62333089368028
differential_evolution step 128: f(x)= 132.62333089368028
differential_evolution step 129: f(x)= 132.62333089368028
differential_evolution step 130: f(x)= 132.62333089368028
differential_evolution step 131: f(x)= 132.62333089368028
differential_evolution step 132: f(x)= 132.62333089368028
differential_evolution step 133: f(x)= 132.62333089368028
differential_evolution step 134: f(x)= 132.62333089368028
differential_evolution step 135: f(x)= 132.62333089368028
differential_evolution step 136: f(x)= 132.62333089368028
differential_evolution step 137: f(x)= 132.62333089368028
differential_evolution step 138: f(x)= 132.62333089368028
differential_evolution step 139: f(x)= 132.38307113018664
differential_evolution step 140: f(x)= 132.38307113018664
differential_evolution step 141: f(x)= 132.38307113018664
differential_evolution step 142: f(x)= 132.38307113018664
differential_evolution step 143: f(x)= 132.38307113018664
differential_evolution step 144: f(x)= 132.38307113018664
differential_evolution step 145: f(x)= 132.38307113018664
differential_evolution step 146: f(x)= 132.38307113018664
differential_evolution step 147: f(x)= 132.38307113018664
differential_evolution step 148: f(x)= 132.38307113018664
differential_evolution step 149: f(x)= 132.38307113018664
differential_evolution step 150: f(x)= 132.38307113018664
differential_evolution step 151: f(x)= 132.38307113018664
differential_evolution step 152: f(x)= 131.16070123323118
differential_evolution step 153: f(x)= 131.16070123323118
differential_evolution step 154: f(x)= 131.16070123323118
differential_evolution step 155: f(x)= 131.16070123323118
differential_evolution step 156: f(x)= 131.16070123323118
differential_evolution step 157: f(x)= 131.16070123323118
differential_evolution step 158: f(x)= 131.16070123323118
differential_evolution step 159: f(x)= 131.16070123323118
differential_evolution step 160: f(x)= 129.6351061463733
differential_evolution step 161: f(x)= 129.6351061463733
differential_evolution step 162: f(x)= 129.6351061463733
differential_evolution step 163: f(x)= 129.6351061463733
differential_evolution step 164: f(x)= 129.6351061463733
differential_evolution step 165: f(x)= 129.6351061463733
differential_evolution step 166: f(x)= 129.6351061463733
differential_evolution step 167: f(x)= 129.6351061463733
differential_evolution step 168: f(x)= 129.6351061463733
differential_evolution step 169: f(x)= 129.6351061463733
differential_evolution step 170: f(x)= 129.6351061463733
differential_evolution step 171: f(x)= 129.6351061463733
differential_evolution step 172: f(x)= 129.6351061463733
differential_evolution step 173: f(x)= 129.6351061463733
differential_evolution step 174: f(x)= 129.6351061463733
differential_evolution step 175: f(x)= 129.6351061463733
differential_evolution step 176: f(x)= 129.6351061463733
differential_evolution step 177: f(x)= 129.6351061463733
differential_evolution step 178: f(x)= 129.6351061463733
differential_evolution step 179: f(x)= 129.6351061463733
differential_evolution step 180: f(x)= 129.6351061463733
differential_evolution step 181: f(x)= 129.6351061463733
differential_evolution step 182: f(x)= 129.6351061463733
differential_evolution step 183: f(x)= 129.6351061463733
differential_evolution step 184: f(x)= 129.6351061463733
differential_evolution step 185: f(x)= 129.6351061463733
differential_evolution step 186: f(x)= 129.6351061463733
differential_evolution step 187: f(x)= 129.4769326550521
differential_evolution step 188: f(x)= 129.4769326550521
differential_evolution step 189: f(x)= 129.4769326550521
differential_evolution step 190: f(x)= 129.4769326550521
differential_evolution step 191: f(x)= 129.4769326550521
differential_evolution step 192: f(x)= 129.4769326550521
differential_evolution step 193: f(x)= 129.4769326550521
differential_evolution step 194: f(x)= 129.4769326550521
differential_evolution step 195: f(x)= 129.4769326550521
differential_evolution step 196: f(x)= 129.4769326550521
differential_evolution step 197: f(x)= 129.4769326550521
differential_evolution step 198: f(x)= 129.4769326550521
differential_evolution step 199: f(x)= 129.4769326550521
differential_evolution step 200: f(x)= 129.4769326550521
Polishing solution with 'L-BFGS-B'

=== Optimization Complete ===
Best χ²: 129.4769326550521

Best-fit parameters:
phases: [3.635595   0.2610597  2.36416248 5.66657328 6.18764673 1.78530843
 2.07157794 4.54991969 4.45837379]
R_u: [1.20829299 1.23372736 1.89332317 1.47187618 1.57712759 0.65250414
 0.75085864 0.79832153 0.71304985]
R_d: [0.25540073 0.56017475 0.09469115 0.65720903 0.54612726 0.15053069
 0.43332006 0.36341715 0.82150549]
R_e: [0.38782779 0.28820986 0.04965427 0.24621418 0.13229251 0.22619418
 0.33592897 0.48586046 0.51975447]
R_nu: [2.47554227 3.98611597 2.17925899 3.32632087 2.65905672 1.90644176
 1.89897239 3.11560845 4.18560559]
sites: [265, 138, 324, 1068, 927, 950, 1537, 1512, 1626]
kappa: 0.24707066954201284

--- Best-fit Yukawa matrices ---
Yu:
[[ 1.39672672e+00-4.42775550e-26j -2.96605582e-11-5.60234750e-10j
  -1.57000958e-10+4.06167863e-10j]
 [-2.96605582e-11+5.60234750e-10j  1.54891063e+00-1.00380474e-25j
  -1.34764931e-09-6.04645752e-10j]
 [-1.57000958e-10-4.06167863e-10j -1.34764931e-09+6.04645752e-10j
   1.50629145e+00+2.58493941e-26j]]

Yd:
[[ 2.28306506e-01+1.59114195e-27j  3.59437324e-11-2.51554673e-11j
   1.00934461e-11+2.52589998e-11j]
 [ 3.59437324e-11+2.51554673e-11j  2.48040922e-01-1.86823598e-26j
  -4.13971899e-11+1.76429064e-10j]
 [ 1.00934461e-11-2.52589998e-11j -4.13971899e-11-1.76429064e-10j
   2.35499056e-01+3.23117427e-27j]]

Ye:
[[ 1.07960030e-01+6.18803391e-27j -2.26026186e-11-3.55246710e-11j
   1.88482160e-11-1.76549061e-11j]
 [-2.26026186e-11+3.55246710e-11j  1.12208873e-01+1.04172954e-27j
   8.39428186e-12+4.45917449e-11j]
 [ 1.88482160e-11+1.76549061e-11j  8.39428186e-12-4.45917449e-11j
   1.07924688e-01+3.23117427e-27j]]

Yn:
[[ 6.93293873e+00-5.36726545e-26j  2.67593526e-10-2.45791692e-09j
   5.41434432e-10+1.05638533e-09j]
 [ 2.67593526e-10+2.45791692e-09j  1.08889064e+01-4.45439291e-26j
  -4.29990809e-09+2.82989160e-09j]
 [ 5.41434432e-10-1.05638533e-09j -4.29990809e-09-2.82989160e-09j
   8.63432803e+00+0.00000000e+00j]]


--- CKM matrix ---
[[-5.29266804e-01+8.48455450e-01j -1.77343880e-08-3.58584006e-08j
   2.50983291e-09+1.88621780e-09j]
 [ 7.76099295e-09-3.92441133e-08j -7.32362707e-01-6.80914727e-01j
   1.20474255e-09-2.57510421e-09j]
 [-2.72000635e-10-3.12779385e-09j -8.71117752e-10+2.70623729e-09j
   1.00000000e+00+3.28107111e-18j]]

--- PMNS matrix ---
[[-7.80649325e-01+6.24969304e-01j  4.70008449e-09+8.75454420e-09j
  -9.42199814e-09-3.90599611e-10j]
 [-7.11117023e-09+6.19338417e-09j -2.00604774e-07-7.02445958e-07j
   1.00000000e+00+4.86224243e-16j]
 [-8.89704354e-09+4.42440083e-09j -2.75484766e-01-9.61305437e-01j
  -7.30528678e-07+6.70700750e-10j]]

Done.

"""