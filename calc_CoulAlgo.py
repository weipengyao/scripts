#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt


# In[2]:


## Estimate of the electron temperature evolution

## for low temperature and low density experiment (Te < 300 eV, ne/nc = 1e-2)
## neglecting the electron thermal diffusion
## accounting only for the inverse Bremsstrahlung laser absorption

## (3/2) n_e dT_e / dt = \nu_B I_0 / qe

## n_e: number density of electron [cm^-3]
## T_e0: you need an initial value for T_e [eV]
## \nu_B: bremsstrahlung laser absorption coefficient [cm^-1]
## I_0: laser intensity, usually assume a Gaussian time profile [W/cm^2]
## qe: elementry charge, [eV = J / 1.6e-19]


# In[29]:


# constants in CGS
c  = 3e10              # cm/s, light speed
qe   = 4.8e-10         # elementary charge
hbar = 1.0e-27         # Planck constant
me   = 9.1e-28         # electron mass, g

## for TS, we use 1ns, 15J, 200um (in plasma) in 2w:
E_ts   = 15.0  # J
tau_ts = 3e-9  # s [?]
d_ts   = 200e-4 # cm
S_ts   = np.pi * (d_ts/2.)**2 # cm^2
# I_ts   = E_ts / tau_ts / S_ts # W/cm2
I_ts   = 3.8e13          # W/cm2,  Here I just took the parameter in the paper
print('Laser intensity for TS = {:.2e} W/cm2'.format(I_ts))


lmd_ts = 1.0*1e-4 # cm, 1 omega laser (1 um)
omg_ts = 2.*np.pi*c/lmd_ts
print('Laser angular frequency for TS = {:.2e} rad/s'.format(omg_ts))


# In[51]:


Te = [10.]            # initial temperature
Te2 = [10.]
Z0 = 2.               # material
fwhm = 6.0e-10 * 1    # equals to tau_ts above
ne = 1e-2*1.0e21      # cm^-3

order = 2.            # order of gaussian profile
t_end = 3.0e-9        # estimation time
center = t_end / 2.0  # center of gaussian profile

omg_pe  = 5.64e4 * ne**0.5
print('plasma electron angular frequency = {:.2e} rad/s'.format(omg_pe))


# create the array for time evolution
num = 100
time = np.linspace(0,t_end,num)
d_time = time[1] - time[0]


# Gaussian profile for the Laser time-evolution
def tgaussian(fwhm,order,center):
    import numpy as np
    sigma = (0.5*fwhm)**order/np.log(2.)    # correction for gaussian profile expressed by fwhm
#     sigma = fwhm**order 
    return np.exp(-(time-center)**order / sigma)

# Coulomb alogrithm, two different methods

def LMD1(temp):
    import numpy as np
    vte  = 4.19e7 * temp**0.5              # cm/s
    V_d1 = hbar/np.sqrt(me*temp*1.6e-12)   # cm, erg s / (g erg)^(1/2), erg -> g cm^2 s^-2, Note here eV -> erg by a factor of 1.6e-12
    V_d2 = Z0*qe**2 / (temp*1.6e-12)       # cm, qe -> m^1/2 l^3/2 / t, 
    V    = np.maximum(omg_ts,omg_pe) * np.maximum(V_d1, V_d2)   # cm/s, from NRL (2019) P58 Eq.32
    return np.log(vte/V)

def LMD2(temp):
    import numpy as np
    Ld = 7.43e2 * temp**0.5 * ne**(-0.5)  # cm
    ni = ne/Z0  # cm^-3
    Ri = ni**(-1./3.)  # cm
    vte  = 4.19e7 * temp**0.5              # cm/s
    bmax = np.minimum(np.maximum(Ld, Ri), vte/omg_ts)
    bmin = np.minimum(np.maximum(Z0*qe**2 / (temp*1.6e-12), hbar/np.sqrt(me*temp*1.6e-12)), Ri)  # bug from SI->CGS ???
    return np.log(1 + 0.7 * bmax/bmin)

# Inverse bremsstrahlung absorption coefficient for radiation of angular frequency 𝜔:
# Considered the temperature dependence
def nu_B(temp,coul_alog):
    return 3.1e-7 * Z0 * ne**2 * coul_alog(temp) * temp**(-1.5) * omg_ts**-2. * (1.-omg_pe**2/omg_ts**2)**(-0.5)   # cm^-1, how to get cm^-1 from this equation????

I_laser = I_ts*tgaussian(fwhm,order,center)

# here's the integration, a typo in the formula (2 or e)?
for t in range(num): 
    dt = (nu_B(Te[t],LMD1)*I_laser[t]/1.5/ne)/1.6e-19 * d_time
    Te.append(Te[t]+dt)
    dt2 = (nu_B(Te2[t],LMD2)*I_laser[t]/1.5/ne)/1.6e-19 * d_time
    Te2.append(Te2[t]+dt2)


# In[55]:


import matplotlib as mpl
import matplotlib.pyplot as plt
plt.rc('text', usetex=False)
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['axes.labelsize'] = 18
mpl.rcParams['legend.fontsize'] = 16#20
# mpl.rcParams['legend.fontname'] = 'Comic Sans MS'
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
mpl.rcParams['lines.linewidth'] = 1.2

width  = 3.14*2 # single column, 8cm
height = width
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

p1,  = ax1.plot(time,Te[:-1],'-b',label='Te, with LMD1')
p12, = ax1.plot(time,Te2[:-1],'--b',label='Te, with LMD2')
p2,  = ax2.semilogy(time,I_laser,'-r',label='I')

lines = [p1, p12, p2]
ax1.legend(lines, [l.get_label() for l in lines])
ax1.yaxis.label.set_color(p1.get_color())
ax2.yaxis.label.set_color(p2.get_color())

# ax1.spines["left"].set_edgecolor(p1.get_color())
ax2.spines["right"].set_edgecolor(p2.get_color())
ax2.spines["left"].set_edgecolor(p1.get_color())

ax1.set_xlabel('time [s]')
ax1.set_ylabel('Te [eV]')
ax2.set_ylabel('I_ts [W/cm2]')
ax2.set_ylim([5e11,4e13])
ax1.set_ylim([0,360])

ax1.tick_params(axis='y', colors=p1.get_color())
ax2.tick_params(axis='y', which='both', colors=p2.get_color())
# fig.legend(loc="upper left", bbox_to_anchor=(0,1), bbox_transform=ax1.transAxes)
fig.set_size_inches(width, height)


# In[54]:


fig.savefig('~/Desktop/'+'TS_heating_estimation'+'.pdf',bbox_inches='tight',dpi=300,transparent=True)


# In[ ]:





# In[47]:


# Z0   = 2
# qe   = 4.8e-10         # CGS
# hbar = 1.0e-27         # Planck constant
# me   = 9.1e-28         # electron mass, g

# temp = np.linspace(10,500,1000)
# V_d1 = hbar/np.sqrt(me*temp*1.6e-12)   # cm, erg s / (g erg)^(1/2), erg -> g cm^2 s^-2, Note here eV -> erg by a factor of 1.6e-12
# V_d2 = Z0*qe**2 / (temp*1.6e-12)       # cm, qe -> m^1/2 l^3/2 / t, 
# plt.plot(temp, V_d2, '-r', label='Ze^2/kT')
# plt.plot(temp, V_d1, '-b', label='hbar/(mkT)')
# plt.xlabel('temperature [eV]')
# plt.legend()
# plt.grid()
# plt.show()


# In[51]:


# plt.plot(temp, nu_B(temp))
# plt.xlabel('temp [eV]')
# plt.ylabel('nu_B [cm^-1 ?]')


# In[61]:


temp = np.linspace(10,350,1000)
lmd1 = LMD1(temp)
lmd2 = LMD2(temp)
plt.plot(temp, lmd1, '-r', label='mine')
plt.plot(temp, lmd2, '-b', label='yours')
plt.xlabel('temperature [eV]')
plt.ylabel('Coulomb Algorithm')
plt.legend()
plt.grid()
plt.show()


# In[58]:


LMD2(10)


# In[5]:


# qe   = 4.8e-10         # CGS
# hbar = 1.0e-27         # Planck constant
# me   = 9.1e-28         # electron mass, g

# ne   = 1e-2*1.0e21     # cm^-3,  in our case, 1e18
# Te0  = 10              # eV,     in our case, ?
# tau  = 600e-12         # s,      in our case, ?
# I0   = 3.8e13          # W/cm2,  in our case, ?
# Z0   = 1.0             # Hydrogen
# # nu_B calculated separately below according to other parameters
# c   = 3e10            # cm/s

# omg_pe  = 5.64e4 * ne**0.5
# print('plasma electron angular frequency = {:.2e} rad/s'.format(omg_pe))


# In[136]:


# ## for \nu_B, in NRL:
# # Te   = np.linspace(Te0,300,100)
# Te0 = 100.
# vte  = 4.19e7 * Te0**0.5   # cm/s
# V    = np.maximum(omg_ts,omg_pe) * (hbar/np.sqrt(me*Te0))
# nu_B = 3.1e-7 * Z0 * ne**2 * np.log(vte/V) * Te0**(-1.5) * omg_ts**-2. * (1.-omg_pe**2/omg_ts**2)**(-0.5)   # cm^-1
# print(
# 'Inverse bremsstrahlung absorption coefficient \n \
# for radiation of angular frequency {:.2e} rad/s \n \
# at initial temperature Te0 = {:.1f} eV is {:.2e} cm^-1'.format(omg_ts,Te0, nu_B)
# )


# In[ ]:




