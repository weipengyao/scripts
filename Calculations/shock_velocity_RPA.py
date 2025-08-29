#!/usr/bin/env python3
"""
Impulsive radiation-pressure acceleration model (non-HB).

v_max = P_rad * tau / (rho * delta),
with delta = c / omega_pe and omega_pe = sqrt(ne * q_e^2 / (eps0 * m_e)).

You can choose either:
  - give ne_direct [m^-3], or
  - give overdensity = ne/nc and wavelength lambda0 to compute nc and ne.

Edit the 'USER INPUTS' section below and run.
"""

import math

# ---------------- Physical constants (SI) ----------------
c     = 2.99792458e8         # m/s
me    = 9.1093837015e-31     # kg
qe    = 1.602176634e-19      # C
eps0  = 8.8541878128e-12     # F/m
pi    = math.pi

# ---------------- Helper functions ----------------
def critical_density(lambda_m: float) -> float:
    """n_c [m^-3] for wavelength lambda_m [m]."""
    omega = 2.0 * pi * c / lambda_m
    return eps0 * me * omega**2 / qe**2

def omega_pe_from_ne(ne: float) -> float:
    """Electron plasma frequency [rad/s] from ne [m^-3]."""
    return math.sqrt(ne * qe**2 / (eps0 * me))

def skin_depth_from_ne(ne: float) -> float:
    """Collisionless electron skin depth delta [m] from ne [m^-3]."""
    return c / omega_pe_from_ne(ne)

def radiation_pressure(I_Wm2: float, R: float = 1.0) -> float:
    """P_rad = (1+R) * I / c  [Pa]."""
    return (1.0 + R) * I_Wm2 / c

# ---------------- USER INPUTS ----------------
rho         = 2.65e3          # kg/m^3  (e.g., SiO2 ~ 2650 kg/m^3)
I_Wcm2      = 3.0e17          # W/cm^2  (peak intensity)
lambda_nm   = 400.0           # nm      (used only if overdensity is provided)
R           = 1.0             # reflectivity (0..1), 1 = perfect mirror
tau_fs      = 186.0           # fs      (pulse duration or effective push time)

# Choose *one* of the following two ways to set ne:
ne_direct   = None            # set to a value in m^-3 to force ne, e.g., 6.9e29
overdensity = 100.0           # ne/nc if ne_direct is None; set to None to disable

# ---------------- Calculations ----------------
# Intensity to SI
I_Wm2 = I_Wcm2 * 1e4

# Determine ne
if ne_direct is not None:
    ne = float(ne_direct)
elif overdensity is not None:
    lambda_m = lambda_nm * 1e-9
    nc = critical_density(lambda_m)
    ne = overdensity * nc
else:
    raise ValueError("Please set either `ne_direct` (m^-3) or `overdensity` (ne/nc).")

# Skin depth, radiation pressure, and velocity
delta = skin_depth_from_ne(ne)                 # m
P_rad = radiation_pressure(I_Wm2, R=R)         # Pa
tau   = tau_fs * 1e-15                         # s
v_max = P_rad * tau / (rho * delta)            # m/s

# ---------------- Output ----------------
print("=== Inputs ===")
print(f"rho        = {rho:.3e} kg/m^3")
print(f"I          = {I_Wcm2:.3e} W/cm^2  ({I_Wm2:.3e} W/m^2)")
print(f"lambda     = {lambda_nm:.1f} nm")
print(f"reflectivity R = {R:.2f}")
print(f"tau        = {tau_fs:.3f} fs")

if ne_direct is not None:
    print(f"ne (direct)= {ne:.3e} m^-3")
else:
    print(f"nc(λ)      = {nc:.3e} m^-3")
    print(f"overdensity= {overdensity:g}")
    print(f"ne         = {ne:.3e} m^-3")

print("\n=== Results ===")
print(f"omega_pe   = {omega_pe_from_ne(ne):.3e} rad/s")
print(f"delta      = {delta*1e9:.3f} nm")
print(f"P_rad      = {P_rad:.3e} Pa  ({P_rad/1e5:.3e} bar)")
print(f"v_max      = {v_max:.3e} m/s")