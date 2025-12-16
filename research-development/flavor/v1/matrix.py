import numpy as np
from scipy.linalg import eigh

def num_divisors(n):
    if n <= 0: return 1.0
    count = 0
    for i in range(1, int(np.sqrt(n))+1):
        if n % i == 0:
            count += 1 if i*i == n else 2
    return float(count)

allowed_distances = {1,2,3,4,5,6,8}   # divisors of 360 ≤8

kappa = 0.24                          # your single geometric parameter
heavy_scale = 2e12                    # perfectly allowed intermediate/GUT scale

Phi = np.zeros((9,9))
for i in range(9):
    for j in range(9):
        d = abs(i-j)
        if d == 0:
            Phi[i,j] = 1.0
        elif d in allowed_distances:
            Phi[i,j] = (kappa ** d) / num_divisors(d)

Y0 = np.ones((9,9))                   # democratic proto-Yukawa

Y = Phi * Y0

light = [3,4,5]
heavy = [0,1,2,6,7,8]

A = Y[np.ix_(light,light)]
B = Y[np.ix_(light,heavy)]
D = Y[np.ix_(heavy,heavy)] + heavy_scale * np.eye(6)

Meff = A - B @ np.linalg.inv(D) @ B.T
Meff = (Meff + Meff.T.conj()) / 2

# Masses and mixing
evals, evecs = eigh(Meff @ Meff.T.conj())
masses = np.sqrt(np.maximum(evals, 0))
masses = np.sort(masses)[::-1]

print("Mass ratios (heaviest : medium : light) =", masses / masses[0])
print("Cabibbo angle θ12 ≈", np.degrees(np.arcsin(np.abs(evecs[0,1]))), "°")
print("\nMeff ≈\n", np.round(Meff, 4))