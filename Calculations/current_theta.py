import numpy as np

# alpha = 0.383
alpha = 0.383 * np.sqrt(2)

u0 = 4.*np.pi * 1e-7
qe = 1.602e-19
mp = 1.6726219e-27

factor = (alpha * u0 * qe) / (np.sqrt(2*mp))

# theta = np.radians(10)  # Convert degrees to radians
# Ep = 10.0e6 * qe  # Convert MeV to Joules

# theta = np.radians(8)  # Convert degrees to radians
# Ep = 25.0e6 * qe  # Convert MeV to Joules

theta = np.radians(5)  # Convert degrees to radians
Ep = 10.0e6 * qe  # Convert MeV to Joules

current = np.sqrt(Ep)* theta / factor / 1e3

print(f"Alpha: {alpha}")
print(f"Theta: {np.degrees(theta):.2f} degrees")
print(f"Current: {current:.2e} kA")
