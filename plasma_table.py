import numpy as np

###############################################################################
### Fundamental parameters
e_charge = 1.6022e-19 # [C] electron charge
e_charge_cgs = 4.8e-10  # [statcoul] 
mass_proton = 1.6726e-24 # [g] proton mass 
mass_electron = 9.109e-28  # [g] electron mass
c_cgs         = 3.0e10     # [cm/s]
gamma = 5./3 # adiabatic index
K_Bol = 1.38e-16 # Boltsmann constant in CGS units
###############################################################################

### Input parameters for the plasma system in CGS units, except T_e, T_i are mostly [eV] (1.16e4 for K conversion)
# V [cm/c] - speed 
# B [G] - magnetic field strength (1 T = 1e4 G)
# T_e [eV] - electron temperature
# T_i [ev] - ion temperature
# n_e [cm-3] - electron density 
# space_scale [cm] - scale of the object
# A averaged mass number
# Z averaged charge state
###############################################################################

###############################################################################
### Input parameters for specified systems:

## magnetic field compression
# V = 3.0e7; 
# B = 60.e4; 
# T_e = 300.; 
# T_i = 300.; 
# n_e = 1.5e20; 
# space_scale = 0.1; # this we don't care for now
# A = 64.0; 
# Z = 19.0


V = 3.0e7; 
B = 60.e4; 
T_e = 300.; 
T_i = 300.; 
n_e = 1.5e20; 
space_scale = 0.1; # this we don't care for now
A = 1.0; 
Z = 1.0


###############################################################################

m_i = A * mass_proton   # [g]
n_i = n_e/Z             # [cm-3]
rho = n_i*m_i           # [g/cm3]

### Printing input parameters: ### 
print ('')
print ('Parameters:------------------')
print ('')
print ('Flow velocity V           = {:.1e}'.format(V),           '[cm/s]')
print ('B-field                   = {:.1e}'.format(B),           '[G]')
print ('Electron temperature T_e  = {:.1f}'.format(T_e),         '[eV]')
print ('Ion temperature T_i       = {:.1f}'.format(T_i),         '[eV]')
print ('Electron density n_e      = {:.1e}'.format(n_e),    	   '[cm-3]')
print ('Ion density n_i           = {:.1e}'.format(n_i),         '[cm-3]')   
# print 'Density rho               = {:.1e}'.format(rho),         '[g/cm-3]'
print ('Charge state Z            = {:.1f}'.format(Z))
print ('Mass number A             = {:.1f}'.format(A))
print ('Spacial scale             = {:.1e}'.format(space_scale), '[cm]')

print ('')

print ('Calculations:----------------')
print ('')

### Calculation of Coulomb logarithm (formula from Braginskii 1965 p.215)
if T_e < 50:
    cou_log_brag = (23.4-1.15*np.log10(n_e)+3.45*np.log10(T_e))
else:
    cou_log_brag = (25.3-1.15*np.log10(n_e)+2.3*np.log10(T_e))
print ('Coulomb Log. (Braginskii) = {:.1f}'.format(cou_log_brag))
print ('')

# ### Calculation of Thermal equilibration (formula from NRL p.34)
Ther_equilib_ie_cou_Brag = 1.0/(1.8e-19*n_e*cou_log_brag*Z**2*(mass_electron*m_i)**(1./2)/(mass_electron*T_i + m_i*T_e)**(3./2))
Ther_equilib_ei_cou_Brag = 1.0/(1.8e-19*n_i*cou_log_brag*Z**2*(mass_electron*m_i)**(1./2)/(mass_electron*T_i + m_i*T_e)**(3./2))
print ('Thermal equilibration i-e = {:.2f}'.format(Ther_equilib_ie_cou_Brag*1e+9), '[ns]')
print ('Thermal equilibration e-i = {:.2f}'.format(Ther_equilib_ei_cou_Brag*1e+9), '[ns]')

print ('')

### Sound velocity and Alfven velocity
c_sound = np.sqrt(5.0*(n_i*T_i*K_Bol*1.16e4+n_e*T_e*K_Bol*1.16e4)/(3.0*rho)) # P_ther*gamma = rho*V_s^2 (1 eV = 1.16e4 K)
print ('Sound velocity            = {:.1e}'.format(c_sound),  '[cm/s]')
c_alfven = B/(4*np.pi*n_i*m_i)**(1./2) # Formula from NRL p.29
print ('Alfven velocity           = {:.1e}'.format(c_alfven), '[cm/s]')
c_magnetosonic = np.sqrt(c_sound**2 + c_alfven**2)
print ('c_magnetosonic velocity   = {:.1e}'.format(c_magnetosonic), '[cm/s]')

print ('')

### electron collision time
f_e = 2.91e-6*Z*(n_e)*cou_log_brag*(T_e)**(-3./2) # [s-1]  Formula from NRL p.28
tau_e=1.0/f_e # [s]
#thermal electron velocity
vth_e= 4.19e7*(T_e)**(1./2) # [cm/sec] Formula from NRL p.29
### electron mean free path
mfp_e = max(vth_e, V)*tau_e # [cm]
print ('Electron thermal velocity = {:.2f}'.format(vth_e*1e-5), '[km/s]')
print ('Electron mean free path   = {:.2e}'.format(mfp_e),      '[cm]')
print ('Electron collision time   = {:.2e}'.format(tau_e*1e+9), '[ns]')
# print ('Electron collisionality   = {:.2e}'.format(mfp_e / space_scale))

print ('')

### ion collision time
f_i = 4.8e-8*Z**4*(A)**(-1./2)*n_i*cou_log_brag*(T_i)**(-3./2) # [s-1]  Formula from NRL p.28
tau_i= 1.0/f_i # [s]
#thermal ion velocity
vth_i= 9.79e5*(A)**(-1./2)*(T_i)**(1./2) # [cm/s] Formula from NRL p.29
 ### ion mean free path
mfp_i = max(vth_i, V)*tau_i # [cm]
print ('Ion thermal velocity      = {:.1e}'.format(vth_i*1e-5), '[km/s]')
print ('Ion mean free path        = {:.1e}'.format(mfp_i),      '[cm]')
print ('Ion collision time        = {:.2e}'.format(tau_i*1e+9), '[ns]')
# print ('Ion collisionality        = {:.2e}'.format(mfp_i / space_scale))

print ('')

### electron gyroradius and gyrofrequency
gyro_e = mass_electron * c_cgs/(e_charge_cgs*B) * max(vth_e, V) # [cm]
print ('Electron Larmor radius    = {:.2e}'.format(gyro_e),     '[cm]')
gyrofreq_e = vth_e/gyro_e/(2.*np.pi) # [s-1]
print ('Electron gyrofrequency    = {:.2e}'.format(gyrofreq_e), '[s-1]')
print ('Electron Larmor period    = {:.1e}'.format(1./gyrofreq_e),     '[s]')
#electron magnetization => collision_time/gyrotime
omega_tau_e = e_charge_cgs*B*tau_e/mass_electron/c_cgs
print ('Electron magnetization    = {:.2f}'.format(omega_tau_e))
print ('Hall parameter for eon He = {:.2f}'.format(mfp_e/gyro_e))
   
print ('') 
 
### ions gyroradius and gyrofrequency
gyro_freq_i = (Z*e_charge_cgs*B)/m_i/c_cgs/(2.*np.pi)  # [s^-1]
print ('gyro_freq_i               = {:.1e}'.format(gyro_freq_i))
print ('Ion Larmor period         = {:.1e}'.format(1e9/gyro_freq_i),     '[ns]')
gyro_i = m_i*c_cgs/(Z*e_charge_cgs*B) * max(vth_i, V) # [cm]
print ('Ion Larmor radius         = {:.2e}'.format(gyro_i*1e4),     '[um]')
gyrofreq_i = vth_i/gyro_i/(2.*np.pi) # [s-1]
print ('Ion gyrofrequency         = {:.2e}'.format(gyrofreq_i), '[s-1]')
print ('Ion gyro-period           = {:.2e}'.format(1.e9/gyrofreq_i), '[ns]')
#ion magnetization => collision_time/gyrotime
omega_tau_i = e_charge_cgs*Z*B*tau_i/m_i/c_cgs
print ('Ion magnetization         = {:.6f}'.format(omega_tau_i))
print ('Hall parameter for ion Hi = {:.6f}'.format(mfp_i/gyro_i))
# print ('Ion Larmor radius / L_sys = {:.2e}'.format(gyro_i / space_scale))

print ('')

### ion inertial length
di = 2.28e7 * Z**(-1) * (A/n_i)**(0.5)  #[cm]
print ('Ion inertial length di    = {:.2e}'.format(di),     '[cm]')

### electron inertial length
de = 5.31e5 * (1/n_e)**(0.5)  #[cm]
print ('Eon inertial length de    = {:.2e}'.format(de),     '[cm]')

### Cooling time  
cool_rate = 1.7*10e-25*(Z**2)*T_e**(1./2) # normalized cooling rate
cooling_time_thin = 4.0e-36*A*(1+Z)*T_e/(Z*rho*cool_rate)
print ('Cooling rate              = {:.2e}'.format(cool_rate))
print ('Cooling time - thin       = {:.2e}'.format(cooling_time_thin*1e+9), '[ns]')

print ('')

### beta thermic and beta dynamic
P_ther = (n_i*T_i+n_e*T_e)*K_Bol*1.16e4
P_dyn = rho*V**2

beta_ther_CGS = P_ther/(B**2.0/(8*np.pi))
print ('Beta thermic              = {:.1e}'.format(beta_ther_CGS)) # [P_ther/(B^2/(8*pi))]
beta_dyn_CGS = rho*V**2.0/(B**2.0/(8*np.pi)) 
print ('Beta dynamic              = {:.1e}'.format(beta_dyn_CGS)) # [P_dyn/(B^2/(8*pi))]

P_ther_i = (n_i*T_i)*K_Bol*1.16e4
P_ther_e = (n_e*T_e)*K_Bol*1.16e4
beta_ther_i_CGS = P_ther_i/(B**2.0/(8*np.pi))
beta_ther_e_CGS = P_ther_e/(B**2.0/(8*np.pi))
print ('Beta thermic of ion       = {:.1e}'.format(beta_ther_i_CGS)) # [P_ther/(B^2/(8*pi))]
print ('Beta thermic of eon       = {:.1e}'.format(beta_ther_e_CGS)) # [P_ther/(B^2/(8*pi))]


mach_sonore = V/c_sound
print ('Mach number               = {:.2e}'.format(mach_sonore))
mach_alfvenique = V/c_alfven
print ('Alfven Mach number        = {:.2e}'.format(mach_alfvenique))    
mach_magnetosonic = V/c_magnetosonic
print ('Magnetosonic Mach number  = {:.2e}'.format(mach_magnetosonic))    

print ('')

### Reunolds, Magnetic Reynolds and Pecklet - crutial parameters for MHD, they all should be >1 (even >>1)

#ion kinematic viscosity with B, alpha = 1
visc_i_B = 2.0e8*T_e/(Z*B)
visc_i_B = visc_i_B #(cm2/s)
#Ion Reynolds number with magnetic field
Re_i_B = space_scale*V/visc_i_B
print ('Reynolds number (with B)  = {:.1e}'.format(Re_i_B))

# ion kinematic viscosity without B
Mu = 2.0e19 * T_i**(2.5) / cou_log_brag / A**0.5 / Z**4 / n_i
Re = space_scale * V / Mu
print ('Reynolds number (w/o  B)  = {:.1e}'.format(Re))
print ('')

#sptizer diffusivity
Spitzer_condu = n_e*e_charge**2.0*tau_e/mass_electron
print ('Spitzer conductivity      = {:.1e}'.format(Spitzer_condu))
Spitzer_resi = 1.0/Spitzer_condu
print ('Spitzer resistivity       = {:.1e}'.format(Spitzer_resi))
mag_diff = Spitzer_resi/(4.*np.pi)
print ('Magnetic diffusivity      = {:.2e}'.format(mag_diff*1e-4), '[m^2/s]')
mag_diff_time = space_scale**2/mag_diff
print ('Magnetic diffusion time   = {:.2e}'.format(mag_diff_time*1e+9), '[ns]')     
# Magnetic Reynolds ### advection (induction) of B by plasma motion / diffusion of B in the plasma
Rm = V*space_scale/mag_diff  
print ('Magnetic Reynolds number  = {:.1e}'.format(Rm))

print ('')
# Double Check
# Dm = 8.2e5 * cou_log_brag * Z / T_e**(1.5)
# Rm1 = space_scale * V / Dm
# print 'Magnetic Reynolds number 1  = {:.1f}'.format(Rm1)





### electron thermal diffusivity with magnetic field, alpha = 1
kappa_e_B = 8.6e9*np.sqrt(A)*T_e/(Z*B)
kappa_e_B = kappa_e_B
### electronic peclet number with B
Pe_e_B = space_scale*V/kappa_e_B
print ('Peclet number (with B)    = {:.1e}'.format(Pe_e_B))

# electron thermal diffusivity without magnetic field
Oe = 2.21e21 * T_e**2.5 / cou_log_brag / (Z+1) / n_e
Pe = space_scale * V / Oe
print ('Peclet number (w/o  B)    = {:.1e}'.format(Pe))

print ('')



### Euler and Alfven numbers
Euler_CGS = (P_dyn/P_ther)**(1./2)
# print 'Euler number              = {:.2f}'.format(Euler_CGS) # [sqrt(P_dyn/P_ther)]

Alfven_CGS = B/(4*np.pi*P_ther)**(1./2)
# print 'Alfven number             = {:.2f}'.format(Alfven_CGS) # [B/sqrt(4pi*P)] 

print ('')

## optical depth
g_bar = 1.2

od = 5e-38 * n_e * n_i * Z**2 * g_bar * space_scale * T_e**-3.5
print ('Bremsstrahlung opt. dep.  = {:.1e}'.format(od))

