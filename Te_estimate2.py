## Estimate of the electron temperature evolution heating by nanosecond lasers
##
## for low temperature and low density experiment (Te < 300 eV, ne/nc = 1e-2)
## neglecting the electron thermal diffusion
## accounting only for the inverse Bremsstrahlung laser absorption
##
## Equation:
## (3/2) n_e dT_e / dt = \nu_B I_0 / qe
##
## Parameters:
##     n_e   : number density of electron [cm^-3]
##     T_e0  : you need an initial value for T_e [eV]
##     \nu_B : bremsstrahlung laser absorption coefficient [cm^-1]
##     I_0   : laser intensity, usually assume a Gaussian time profile [W/cm^2]
##     qe    : elementry charge, [eV = J / 1.6e-19]

import numpy as np

def constants():
    # constants in CGS
    c    = 3e10            # cm/s, light speed
    qe   = 4.8e-10         # elementary charge
    hbar = 1.0e-27         # Planck constant
    me   = 9.1e-28         # electron mass, g
    return c, qe, hbar, me

c, qe, hbar, me = constants()

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
    print('Laser intensity for TS = {:.2e} W/cm2'.format(I_ts))
    omg_ts = 2.*np.pi*c/lmd_ts
    print('Laser angular frequency for TS = {:.2e} rad/s'.format(omg_ts))
    return E_ts, tau_ts, d_ts, lmd_ts, S_ts, I_ts, omg_ts

E_ts, tau_ts, d_ts, lmd_ts, S_ts, I_ts, omg_ts = setup_laser(110.0, 5.0e-9, 0.1, 1.0*1e-4)



# setup initial condition
def initial_condition(Te, Z0, ne, t_end):
    # Te = [0.1]            # initial temperature
    # Z0 = 1.               # plasma material
    # ne = 2.0e17           # electron number density cm^-3
    # t_end = 10e-9         # estimation time
    omg_pe  = 5.64e4 * ne**0.5
    print('plasma electron angular frequency = {:.2e} rad/s'.format(omg_pe))
    return Te, Z0, ne, t_end, omg_pe

Te, Z0, ne, t_end, omg_pe = initial_condition([0.1], 1., 2.0e17, 10e-9)


# create the array for time evolution
def time_evo(num, t_end):
    time = np.linspace(0,t_end,num)
    d_time = time[1] - time[0]
    return time, d_time

time, d_time = time_evo(1000, t_end)


# Gaussian profile for the Laser time-evolution
def tgaussian(fwhm,order,center):
    '''
    fwhm = tau_ts         # equals to tau_ts above
    order = 4.            # order of gaussian profile
    center = t_end / 2.0  # center of gaussian profile
    '''
    sigma = (0.5*fwhm)**order/np.log(2.)
    return np.exp(-(time-center)**order / sigma)

factor = 0.25
I_laser = I_ts*tgaussian(tau_ts,4.,t_end / 2.0)*factor

# calculation of the Coulomb logarithm
def LMD2(temp, ne, Z0, omg_ts):
    Ld = 7.43e2 * temp**0.5 * ne**(-0.5)  # cm
    ni = ne/Z0  # cm^-3
    Ri = ni**(-1./3.)  # cm
    vte  = 4.19e7 * temp**0.5              # cm/s
    bmax = np.minimum(np.maximum(Ld, Ri), vte/omg_ts)
    bmin = np.minimum(np.maximum(Z0*qe**2 / (temp*1.6e-12), hbar/np.sqrt(me*temp*1.6e-12)), Ri)  # bug from SI->CGS ???
    return np.log(1 + 0.7 * bmax/bmin)

# Inverse bremsstrahlung absorption coefficient for radiation of angular frequency 𝜔:
# Considered the temperature dependence
def nu_B(temp,coul_alog, Z0, ne, omg_ts, omg_pe):
    return 3.1e-7 * Z0 * ne**2 * coul_alog(temp, ne, Z0, omg_ts) * temp**(-1.5) * omg_ts**-2. * (1.-omg_pe**2/omg_ts**2)**(-0.5)   # cm^-1

def integrate_Te_evolution(
    Te0,             # a list with initial temperature [eV]
    ne,              # electron number density [cm^-3]
    I_laser,         # laser intensity array [W/cm^2]
    nu_B_func,       # function: nu_B(temp, coul_alog_func)
    coul_alog_func,  # function: coul_alog(temp)
    d_time           # time step [s]
):
    """
    Integrate electron temperature over time using inverse Bremsstrahlung heating.

    Returns:
        Te_list: list of electron temperatures [eV] for each time step
    """
    Te_list = Te0
    num_steps = len(I_laser)

    for t in range(num_steps):
        # (3/2) n_e dT_e/dt = ν_B * I / q_e
        dt = (nu_B_func(Te_list[t], coul_alog_func, Z0, ne, omg_ts, omg_pe) * I_laser[t] / (1.5 * ne)) \
             / 1.6e-19 * d_time
        Te_list.append(Te_list[t] + dt)

    return Te_list

Te = integrate_Te_evolution(
    Te0=Te,
    ne=ne,
    I_laser=I_laser,
    nu_B_func=nu_B,
    coul_alog_func=LMD2,
    d_time=d_time
)

print('final Te = {:.0f} eV'.format(Te[-1]))
# print(Te[-1])
