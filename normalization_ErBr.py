#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 08:49:52 2022

@author: yz
"""
me = 9.1e-31 # [kg]
c  = 3.0e8   # [m/s]
e  = 1.6e-19 # [C]
# ne = 1.0e18  # [cm^-3]
nc = 1.0e21  # [cm^-3]
ne = 0.02*nc  # [cm^-3]

wp = 5.64e4 * ne**0.5 # [rad/s]
Er = me*c*wp/e
Br = me*wp/e

print("Er = {:.0e} V/m".format(Er))
print("Br = {:.0e} T".format(Br))

B1 = 0.0021 * Br

wpi = 1.32e3 * ne**0.5 # [rad/s], assuming H+
print("1/wpi = {:.0e} ns".format(6.28/wpi*1e9))

wc = 0.014 * wp

B = wc*me/e
print("B = {:.1f} T".format(B))