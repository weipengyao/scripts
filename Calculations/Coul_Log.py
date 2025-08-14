import numpy as np

def constants():
    # constants in CGS
    c    = 3e10            # cm/s, light speed
    qe   = 4.8e-10         # elementary charge
    hbar = 1.0e-27         # Planck constant
    me   = 9.1e-28         # electron mass, g
    return c, qe, hbar, me

c, qe, hbar, me = constants()

# setup initial condition
def initial_condition(Te, Z0, ne, t_end):
    # Te = 0.1            # initial temperature
    # Z0 = 1.               # plasma material
    # ne = 2.0e17           # electron number density cm^-3
    # t_end = 10e-9         # estimation time
    omg_pe  = 5.64e4 * ne**0.5
    # print('plasma electron angular frequency = {:.2e} rad/s'.format(omg_pe))
    return Te, Z0, ne, t_end, omg_pe

Te, Z0, ne, t_end, omg_pe = initial_condition(10., 1., 1.0e18, 10e-9)

# define the laser
def setup_laser(E_ts, tau_ts, d_ts, lmd_ts):
    '''
    E_ts   = 110.0   # Laser energy (J)
    tau_ts = 5.0e-9  # Laser duration (s) 
    d_ts   = 0.1     # Laser focal spot diameter (cm) 
    lmd_ts = 1.0*1e-4 # cm, 1 omega laser (1 um)
    '''
    S_ts   = np.pi * (d_ts/2.)**2 # cm^2
    I_ts   = E_ts / tau_ts / S_ts # W/cm2
    # print('Laser intensity for TS = {:.2e} W/cm2'.format(I_ts))
    omg_ts = 2.*np.pi*c/lmd_ts
    # print('Laser angular frequency for TS = {:.2e} rad/s'.format(omg_ts))
    return E_ts, tau_ts, d_ts, lmd_ts, S_ts, I_ts, omg_ts

E_ts, tau_ts, d_ts, lmd_ts, S_ts, I_ts, omg_ts = setup_laser(110.0, 2.0e-9, 0.1, 1.0*1e-4)

# calculation of the Coulomb logarithm
def LMD2(temp, ne, Z0, omg_ts):
    Ld = 7.43e2 * temp**0.5 * ne**(-0.5)  # cm
    ni = ne/Z0  # cm^-3
    Ri = ni**(-1./3.)  # cm
    vte  = 4.19e7 * temp**0.5              # cm/s
    bmax = np.minimum(np.maximum(Ld, Ri), vte/omg_ts)
    bmin = np.minimum(np.maximum(Z0*qe**2 / (temp*1.6e-12), hbar/np.sqrt(me*temp*1.6e-12)), Ri) 
    return np.log(1 + 0.7 * bmax/bmin)

print(LMD2(Te, ne, Z0, omg_ts))