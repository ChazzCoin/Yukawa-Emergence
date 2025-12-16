import numpy as np

# Placeholder for harmonic pipeline
# Define divisor set
D360 = [1,2,3,4,5,6,8,9,10,12]

phi = 97/60

def triadic_magnitudes(n0, phi=phi):
    # returns s1,s2,s3 as complex magnitudes
    s1 = phi**-1 * np.exp(1j*n0)
    s2 = phi**-2 * np.exp(1j*2*n0)
    s3 = phi**-3 * np.exp(1j*3*n0)
    return [s1, s2, s3]

def assign_sites():
    # map 9 sites into 3 generations
    gens = {0:0,1:1,2:2,3:0,4:1,5:2,6:0,7:1,8:2}
    return gens

def build_proto_Y(s_vals, phases):
    Y = np.zeros((9,9), dtype=complex)
    gens = assign_sites()
    for i in range(9):
        for j in range(9):
            si = s_vals[gens[i]]
            sj = s_vals[gens[j]]
            Y[i,j] = si * np.conj(sj) * np.exp(1j*(phases[i]-phases[j]))
    return Y

# Example usage projections
if __name__ == "__main__":
    n0 = 3
    s_vals = triadic_magnitudes(n0)
    phases = np.zeros(9)
    Y0 = build_proto_Y(s_vals, phases)
    print(Y0)
