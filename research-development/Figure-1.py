import numpy as np
import matplotlib.pyplot as plt

# # Charged-lepton masses (GeV, arbitrary normalization ok)
# m_mu, m_tau = 0.106, 1.777
# He = np.diag([m_mu**2, m_tau**2])
#
# # Degenerate neutrino block
# Hnu = np.eye(2)

def R(theta):
    return np.array([[np.cos(theta), -np.sin(theta)],
                     [np.sin(theta),  np.cos(theta)]])

def alignment_energy(theta):
    comm = He @ R(theta).T @ Hnu @ R(theta) - R(theta).T @ Hnu @ R(theta) @ He
    return np.linalg.norm(comm, 'fro')**2

# theta = np.linspace(0, np.pi, 500)
# E = np.array([alignment_energy(t) for t in theta])
#
# plt.plot(theta, E)
# plt.xlabel(r'$\theta$')
# plt.ylabel(r'$E(\theta)$')
# plt.title('Alignment Energy in the $\mu$–$\tau$ Plane')
# plt.tight_layout()
# plt.show()
#

def theta23_from_U(U):
    return np.arctan2(abs(U[1,2]), abs(U[2,2]))

# mock Ue, Uν for demonstration
Ue = np.eye(3)
Unu = np.eye(3)

theta_vals = np.linspace(0, np.pi, 200)
theta23_vals = []

for t in theta_vals:
    R3 = np.eye(3)
    R3[1:3,1:3] = R(t)
    U = Ue.T @ Unu @ R3
    theta23_vals.append(theta23_from_U(U))

plt.plot(theta_vals, theta23_vals)
plt.xlabel(r'$\theta$')
plt.ylabel(r'$\theta_{23}$')
plt.show()
