#!/usr/bin/env python3
"""
Radiation-pressure hole-boring velocity calculator (simple version).

You can directly set the input parameters below (rho, intensity, lambda0, reflectivity).
"""

import math

# ---------------- Input parameters ----------------
rho      = 2.65e3        # mass density [kg/m^3]
I_Wcm2   = 3.0e17        # laser intensity [W/cm^2]
lambda0  = 400e-9        # laser wavelength [m]
reflectivity = 1.0       # 0 = full absorption, 1 = perfect reflection

# ---------------- Physical constants ----------------
c   = 2.998e8             # m/s
me  = 9.109e-31           # kg
qe  = 1.602e-19           # C
eps0 = 8.854e-12          # F/m

# ---------------- Calculations ----------------
# Convert intensity to SI
I_Wm2 = I_Wcm2 * 1e4      # W/cm^2 -> W/m^2

# Critical density at given wavelength
omega = 2.0 * math.pi * c / lambda0
nc = eps0 * me * omega**2 / qe**2   # [m^-3]

# Radiation pressure
P = (1.0 + reflectivity) * I_Wm2 / c   # Pa

# Hole-boring velocity
v_HB = math.sqrt(P / rho)              # m/s

# ---------------- Output ----------------
print(f"Inputs:")
print(f"  rho = {rho:.3e} kg/m^3")
print(f"  I   = {I_Wcm2:.3e} W/cm^2  ({I_Wm2:.3e} W/m^2)")
print(f"  λ   = {lambda0*1e9:.1f} nm")
print(f"  Reflectivity = {reflectivity}")

print("\nResults:")
print(f"  Critical density n_c = {nc:.3e} m^-3")
print(f"  Radiation pressure P = {P:.3e} Pa")
print(f"  Hole-boring velocity v_HB = {v_HB:.3e} m/s")
