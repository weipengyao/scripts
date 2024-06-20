import numpy as np

me = 9.1e-31
mp = 1836.*me
qe = 1.6e-19

E = 150e3 # eV
B = 20.0 # T

v = np.sqrt(2.*E*qe/mp)
r = mp * v / qe / B
t = r/v

print('v = {:.1e} m/s'.format(v))
print('r = {:.1e} m'.format(r))
print('t = {:.1e} s'.format(t))

