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
# A averaged mass number (for (C_2F_4)_n = (2*12.0096*m_c+4*18.998*m_f*2)/(2*12.0096+4*18.998*2)) = 17.32*m_p)
# Z averaged charge state
###############################################################################

###############################################################################
### Input parameters for specified systems:
### uncomment only one line of the system of the interest to have the output values 
### *** LABORATORY *** ###
### ELFIE Lab Parameters (10, 20 and 30 T)
#V = 6e7; B = 10e4; T_e = 40; T_i = 40; n_e = 1.0e19; space_scale = 0.1; A = 17.32; Z = 5.92
#V = 6e7; B = 20e4; T_e = 40; T_i = 40; n_e = 1.0e19; space_scale = 0.1; A = 17.32; Z = 5.92
#V = 6e7; B = 30e4; T_e = 40; T_i = 40; n_e = 1.0e19; space_scale = 0.1; A = 17.32; Z = 5.92
### PEARL Lab line-focusing Parameters (2, 5 and 10 eV)
# V = 0.25e7; B = 13.5e4; T_e = 50; T_i = 30; n_e = 1.0e19; space_scale = 0.1; A = 17.32; Z = 6
# at the interface
# V = 0.25e7; B = 13.5e4; T_e = 50; T_i = 25; n_e = 1.0e18; space_scale = 0.01; A = 17.32; Z = 6
#V = 1e7; B = 13.5e4; T_e = 5; T_i = 5; n_e = 1.0e18; space_scale = 0.1; A = 17.32; Z = 2.2
# V = 1e7; B = 13.5e4; T_e = 100; T_i = 10; n_e = 1.0e18; space_scale = 0.1; A = 17.32; Z = 8
### PEARL Lab fingers Parameters
# V = 1.e6; B = 13.5e4; T_e = 1.; T_i = 1.; n_e = 1.e18; space_scale = 0.1; A = 12.0; Z = 1.0

# V = 0.8e7; B = 13.5e4; T_e = 30; T_i = 55; n_e = 1.5e18*4; space_scale = 0.01; A = 17.32; Z = 4

### KROT Lab Parameters (plasma facility, no lasers)
# V = 0.25e7; B = 450; T_e = 4; T_i = 1; n_e = 3.0e13; space_scale = 3; A = 10.43; Z = 0.272
# V = 0.25e7; B = 450; T_e = 4; T_i = 1; n_e = 1.0e13; space_scale = 3; A = 10.43; Z = 2.67
# V = 0.35e7; B = 100; T_e = 1; T_i = 0.5; n_e = 3.0e13; space_scale = 3; A = 10.43; Z = 2.67
### Plasma Slab (from Benjamin's)
# V = 5.e7; B = 10.e4; T_e = 100.; T_i = 100.; n_e = 1.e19; space_scale = 0.1; A = 17.32; Z = 5.92

### Magnetized shock (Xavier's FLASH simulation)
# V = 0.5e8; B = 20.e4; T_e = 100.; T_i = 200.; n_e = 1.e18; space_scale = 0.1; A = 1.0; Z = 1.0
# V = 8.0e7; B = 20.e4; T_e = 60.; T_i = 20.; n_e = 3.e18; space_scale = 0.1; A = 1.0; Z = 1.0
# V = 8.6e7; B = 77.e4; T_e = 150.; T_i = 200.; n_e = 4.e18; space_scale = 0.1; A = 1.0; Z = 1.0

# Piston
# V = 1.2e8; B = 5.e4; T_e = 80.; T_i = 40.; n_e = 1.e19; space_scale = 0.1; A = 17.3; Z = 8.0
# Shock Local
# V = 1.5e8; B = 60.e4; T_e = 100.; T_i = 200.; n_e = 1.e18; space_scale = 0.1; A = 1.0; Z = 1.0
# Shock upstream
# V = 1.5e8; B = 20.e4; T_e = 50.; T_i = 20.; n_e = 1.e18; space_scale = 0.1; A = 1.0; Z = 1.0

# LULI shock proposal 2022
# V = 1.5e8; B = 20.e4; T_e = 10.; T_i = 0.1; n_e = 1.e18; space_scale = 0.1; A = 1.0; Z = 1.0
# V = 1.5e8; B = 20.e4; T_e = 100.; T_i = 200.; n_e = 1.e18; space_scale = 0.1; A = 1.0; Z = 1.0
# V = 5.5e7; B = 40.e4; T_e = 10.; T_i = 0.1; n_e = 4.e19; space_scale = 0.1; A = 1.0; Z = 1.0
# V = 5.5e7; B = 40.e4; T_e = 100.; T_i = 200.; n_e = 2.5e19; space_scale = 0.1; A = 1.0; Z = 1.0
# V = 1.5e8; B = 60.e4; T_e = 10.; T_i = 0.1; n_e = 2.e19; space_scale = 0.1; A = 1.0; Z = 1.0
# V = 1.5e8; B = 60.e4; T_e = 100.; T_i = 200.; n_e = 2.e19; space_scale = 0.1; A = 1.0; Z = 1.0
# V = 1.5e8; B = 60.e4; T_e = 100.; T_i = 200.; n_e = 1.e18; space_scale = 0.1; A = 1.0; Z = 1.0

## Apollon fast shock 2022 Oct.
# V = 7.5e8; B = 60.e4; T_e = 100.; T_i = 200.; n_e = 1.e19; space_scale = 0.1; A = 4.0; Z = 2.0
# V = 8.2e8; B = 60.e4; T_e = 100.; T_i = 200.; n_e = 5.e19; space_scale = 0.1; A = 1.0; Z = 1.0

# B compression
# normalization
# V = 0.5*9.1e+06; B = 400.e4; T_e = 300.; T_i = 300.; n_e = 1.3e21; space_scale = 0.1; A = 64.0; Z = 19.0
# foot
# V = 0.5*1.3e+07; B = 40.e4; T_e = 300.; T_i = 300.; n_e = 1.3e21; space_scale = 0.1; A = 64.0; Z = 19.0
# cs
V = 0.5*9.1e+06; B = 400.e4; T_e = 30.; T_i = 30.; n_e = 1.3e21; space_scale = 0.03; A = 64.0; Z = 9.0

# B table for astro
# V = 2.5e6; B = 20e-6; T_e = 0.8e3; T_i = 0.8e3; n_e = 7.0; space_scale = 4.6e19; A = 1; Z = 1
# Scollini in the CS
# V = 2.5e7; B = 30e-9*1e4; T_e = 1e6/1e4; T_i = 1e6/1e4; n_e = 2.0; space_scale = 1.5e7*1e5; A = 1; Z = 1
# Scollini in the foot
# V = 2.5e7; B = 10e-9*1e4; T_e = 2.5e4/1e4; T_i = 2.5e4/1e4; n_e = 2.0; space_scale = 1.5e7*1e5; A = 1; Z = 1



# CME at 10 G
# V = 4.6e+08; B = 100; T_e = 1e7/1.16e4; T_i = 1e7/1.16e4; n_e = 1e12; space_scale = 0.1; A = 1.28; Z = 1
# V = 0.8e+08; B = 10; T_e = 2.3e5/1.16e4; T_i = 2.3e5/1.16e4; n_e = 1e10; space_scale = 2e11; A = 1.28; Z = 1
# CME at 100 G
# V = 2.5e+08; B = 100; T_e = 2.3e6/1.16e4; T_i = 2.3e6/1.16e4; n_e = 1e11; space_scale = 2e11; A = 1.28; Z = 1
# CME at 1 G
# V = 0.24e+08; B = 1; T_e = 2.3e4/1.16e4; T_i = 2.3e4/1.16e4; n_e = 1e9; space_scale = 2e11; A = 1.28; Z = 1

# CME in lab with B = 3e5 G
# V = 0.6e+08; B = 3e5; T_e = 5.0e5/1.16e4; T_i = 5.0e5/1.16e4; n_e = 7.5e18; space_scale = 1.0; A = 17.3; Z = 5.9



# MLPI 2022
# V = 1.0e5; B = 20.e4; T_e = 800.; T_i = 10.; n_e = 3.e17; space_scale = 0.6; A = 1.0; Z = 1.0
# V = 1.0e7; B = 20.e4; T_e = 100.; T_i = 60.; n_e = 6.e18; space_scale = 0.6; A = 1.0; Z = 1.0
# V = 6e9; B = 20.e4; T_e = 100.; T_i = 100.; n_e = 0.02*1.1e21; space_scale = 0.6; A = 1.0; Z = 1.0

# MRTI
# V = 3e7*2; B = 20.e4*1; T_e = 100.; T_i = 100.; n_e = 1e19; space_scale = 0.6; A = 1.0; Z = 1.0
# V = 0.0; B = 20.e4*1; T_e = 700.; T_i = 700.; n_e = 6e18; space_scale = 0.6; A = 1.0; Z = 1.0


# TNSA-B propagation
# V = 4.2e8; B = 20.e4; T_e = 140.; T_i = 140; n_e = 2.e18; space_scale = 0.1; A = 1.0; Z = 1.0 

# RAL-2022-TAW before RCF analysis using PROBLEM
# V = 1.0e7; B = 20.e4; T_e = 150.; T_i = 100; n_e = 1.44e19; space_scale = 0.1; A = 20.0; Z = 10.0  
## NOTE here V is for the turbulence

# Weakly Collisional shock
# V = 5.0e7; B = 0.01; T_e = 50.; T_i = 50; n_e = 1.0e19; space_scale = 0.1; A = 14.0; Z = 3.0  


### Earth's bow shock
# V = 5.0e7; B = 2.5e-4; T_e = 5.0; T_i = 15.0; n_e = 10.0; space_scale = 1e7; A = 1.0; Z = 1.0
### Non-relativistic SNR W28
# V = 2.0e6; B = 7.0e-4; T_e = 0.1; T_i = T_e; n_e = 1.0e3; space_scale = 6e19; A = 1.0; Z = 1.0
### Kes 78
# V = 1.0e7; B = 1.5e-3; T_e = 1.0; T_i = T_e; n_e = 1.0e5; space_scale = 6e19; A = 1.0; Z = 1.0
### W44
# V = 0.5e7; B = 0.4e-3; T_e = 1.0; T_i = T_e; n_e = 1.0e4; space_scale = 6e19; A = 1.0; Z = 1.0
### IC443
# V = 0.6e7; B = 0.3e-3; T_e = 1.0; T_i = T_e; n_e = 1.0e4; space_scale = 6e19; A = 1.0; Z = 1.0
### W49B
# V = 0.3e7; B = 0.12e-3; T_e = 1.0; T_i = T_e; n_e = 3.0e3; space_scale = 6e19; A = 1.0; Z = 1.0
### termination shock
# V = 3.5e7; B = 1.0e-6; T_e = 1.0; T_i = T_e; n_e = 1.0e-3; space_scale = 6e8; A = 1.0; Z = 1.0

# V = 4.5e7; B = 32.0e4; T_e = 300.0; T_i = 670.0; n_e = 3.0e+18; space_scale = 0.01; A = 17.3; Z = 8.0

# V = 4.5e7; B = 30.0e4; T_e = 70.0; T_i = 80.0; n_e = 5.0e+19; space_scale = 0.01; A = 17.3; Z = 8.0


### Nat. Phys., A. Bondarenko, et al., (2017)
# V = 6.0e7; B = 7.1e2; T_e = 4.3; T_i = 0.5; n_e = 7.2e12; space_scale = 50.0; A = 12.0; Z = 4.0



### Pancake -> three stages
### 1st stage @ 5ns at shock layer
# V = 3.5e7;  B = 25.e4; T_e = 150.; T_i = 400.; n_e = 1.e19; space_scale = 0.01; A = 17.32; Z = 6.5
### 2st stage @ 30ns at flutes
# V = 1.5e7; B = 29.e4; T_e = 120.; T_i = 90.; n_e = 5.e17; space_scale = 0.01; A = 17.32; Z = 7.0
### 3st stage @ 50ns at pancake
# V = 4.0e7; B = 30.e4; T_e = 100.; T_i = 100.; n_e = 1.e18; space_scale = 0.01; A = 17.32; Z = 7.0

###############################################################################
### *** ASTRO *** ###
### Solar CME front structure Parameters
#V = 5e7; B = 1; T_e = 172; T_i = 172; n_e = 1e8; space_scale = 1e+11; A = 1; Z = 1;
### Solar CME cavity Parameters
#V = 5e7; B = 1; T_e = 172; T_i = 172; n_e = 1e7; space_scale = 1e+11; A = 1; Z = 1;
### Solar CME prominence core Parameters
#V = 5e7; B = 10; T_e = 6.9; T_i = 6.9; n_e = 1e11; space_scale = 1e+11; A = 1; Z = 1;
### Stellar CME from Argiroffi 2019 Parameters
#V = 1e7; B = 3; T_e = 345; T_i = 345; n_e = 5.5e8; space_scale = 1e+11; A = 1; Z = 1;
### Hypothetic Stellar CME with 75 G (potentially suppressed by B-field)
#V = 38e7; B = 75; T_e = 345; T_i = 345; n_e = 4.0e10; space_scale = 1e+11; A = 1; Z = 1;
### Hypothetic Stellar CME with 300 km/s speed and 1e9 m-3 density (exercise to match lab parameters)
#V = 3e7; B = 1; T_e = 2; T_i = 2; n_e = 1.0e9; space_scale = 1e+11; A = 1; Z = 1;
### 'Tongues' of equatorial accretion (V: 100s km/s, B: 10-100 G, T: 1e3-1e4 K (0.1-1 eV), n_e: 1e11-1e14 cm-3)
# V = 1.0e7; B = 100; T_e = 0.1; T_i = 0.1; n_e = 1e13; space_scale = 1e+10; A = 1; Z = 1;
###############################################################################

# m_i = A*(mass_proton+mass_electron) # [gg]
m_i = A * mass_proton   # [g]
# n_e = n_i * Z           # [cm-3]
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
print ('Spatial scale             = {:.1e}'.format(space_scale), '[cm]')

print ('')

print ('Calculations:----------------')
print ('')
### Calculation of Coulomb logarithm (formula from Braginskii 1965 p.215)
if T_e < 50:
    cou_log_brag = (23.4-1.15*np.log10(n_e)+3.45*np.log10(T_e))
else:
    # cou_log_brag = (25.3-1.15*np.log10(n_e)+2.3*np.log10(T_e))*(T_e>50.0)
    cou_log_brag = (25.3-1.15*np.log10(n_e)+2.3*np.log10(T_e))
# cou_log_brag = 10.
# print ('Coulomb Log. (Braginskii) = {:.1f}'.format(cou_log_brag))
# print ('')

# ### Calculation of Thermal equilibration (formula from NRL p.34)
Ther_equilib_ie_cou_Brag = 1.0/(1.8e-19*n_e*cou_log_brag*Z**2*(mass_electron*m_i)**(1./2)/(mass_electron*T_i + m_i*T_e)**(3./2))
Ther_equilib_ei_cou_Brag = 1.0/(1.8e-19*n_i*cou_log_brag*Z**2*(mass_electron*m_i)**(1./2)/(mass_electron*T_i + m_i*T_e)**(3./2))
# print ('Thermal equilibration i-e = {:.2f}'.format(Ther_equilib_ie_cou_Brag*1e+9), '[ns]')
# print ('Thermal equilibration e-i = {:.2f}'.format(Ther_equilib_ei_cou_Brag*1e+9), '[ns]')

# print ('')

### Sound velocity and Alfven velocity
c_sound = np.sqrt(5.0*(n_i*T_i*K_Bol*1.16e4+n_e*T_e*K_Bol*1.16e4)/(3.0*rho)) # P_ther*gamma = rho*V_s^2 (1 eV = 1.16e4 K)
print ('Sound velocity            = {:.1e}'.format(c_sound),  '[cm/s]')
c_alfven = B/(4*np.pi*n_i*m_i)**(1./2) # Formula from NRL p.29
print ('Alfven velocity           = {:.1e}'.format(c_alfven), '[cm/s]')
c_magnetosonic = np.sqrt(c_sound**2 + c_alfven**2)
print ('c_magnetosonic velocity   = {:.1e}'.format(c_magnetosonic), '[cm/s]')

# print ('')

### electron collision time
f_e = 2.91e-6*Z*(n_e)*cou_log_brag*(T_e)**(-3./2) # [s-1]  Formula from NRL p.28
tau_e=1.0/f_e # [s]
#thermal electron velocity
vth_e= 4.19e7*(T_e)**(1./2) # [cm/sec] Formula from NRL p.29
### electron mean free path
mfp_e = max(vth_e, V)*tau_e # [cm]
# print ('Electron thermal velocity = {:.2f}'.format(vth_e*1e-5), '[km/s]')
print ('Electron mean free path   = {:.2e}'.format(mfp_e),      '[cm]')
print ('Electron collision time   = {:.2e}'.format(tau_e*1e+9), '[ns]')
# print ('Electron collisionality   = {:.2e}'.format(mfp_e / space_scale))

# print ('')

### ion collision time
f_i = 4.8e-8*Z**4*(A)**(-1./2)*n_i*cou_log_brag*(T_i)**(-3./2) # [s-1]  Formula from NRL p.28
tau_i= 1.0/f_i # [s]
#thermal ion velocity
vth_i= 9.79e5*(A)**(-1./2)*(T_i)**(1./2) # [cm/s] Formula from NRL p.29
 ### ion mean free path
mfp_i = max(vth_i, V)*tau_i # [cm]
# print ('Ion thermal velocity      = {:.1e}'.format(vth_i*1e-5), '[km/s]')
print ('Ion mean free path        = {:.1e}'.format(mfp_i),      '[cm]')
print ('Ion collision time        = {:.2e}'.format(tau_i*1e+9), '[ns]')
# print ('Ion collisionality        = {:.2e}'.format(mfp_i / space_scale))

# print ('')

### electron gyroradius and gyrofrequency
# gyro_e = mass_electron/(e_charge*B) * max(vth_e, V) # [cm]
gyro_e = mass_electron * c_cgs/(e_charge_cgs*B) * max(vth_e, V) # [cm]
print ('Electron Larmor radius    = {:.2e}'.format(gyro_e),     '[cm]')
gyrofreq_e = vth_e/gyro_e/(2.*np.pi) # [s-1]
# print ('Electron gyrofrequency    = {:.2e}'.format(gyrofreq_e), '[s-1]')
print ('Electron Larmor period    = {:.1e}'.format(1./gyrofreq_e),     '[s]')
#electron magnetization => collision_time/gyrotime
omega_tau_e = e_charge_cgs*B*tau_e/mass_electron/c_cgs
# print ('Electron magnetization    = {:.2f}'.format(omega_tau_e))
print ('Hall parameter for eon He = {:.2f}'.format(mfp_e/gyro_e))
   
# print ('') 
 
### ions gyroradius and gyrofrequency
gyro_freq_i = (Z*e_charge_cgs*B)/m_i/c_cgs/(2.*np.pi)  # [s^-1]
gyro_i = m_i*c_cgs/(Z*e_charge_cgs*B) * max(vth_i, V) # [cm]

# print ('gyro_freq_i               = {:.1e}'.format(gyro_freq_i))
print ('Ion Larmor radius         = {:.2e}'.format(gyro_i*1e4),     '[um]')
print ('Ion Larmor period         = {:.1e}'.format(1e9/gyro_freq_i),     '[ns]')
# gyrofreq_i = vth_i/gyro_i/(2.*np.pi) # [s-1]
# print ('Ion gyrofrequency         = {:.2e}'.format(gyrofreq_i), '[s-1]')
# print ('Ion gyro-period           = {:.2e}'.format(1.e9/gyrofreq_i), '[ns]')
#ion magnetization => collision_time/gyrotime
omega_tau_i = e_charge_cgs*Z*B*tau_i/m_i/c_cgs
# print ('Ion magnetization         = {:.2f}'.format(omega_tau_i))
print ('Hall parameter for ion Hi = {:.8f}'.format(mfp_i/gyro_i))
# print ('Ion Larmor radius / L_sys = {:.2e}'.format(gyro_i / space_scale))

# print ('')

### ion inertial length
di = 2.28e7 * Z**(-1) * (A/n_i)**(0.5)  #[cm]
# print ('Ion inertial length di    = {:.2e}'.format(di),     '[cm]')

### electron inertial length
de = 5.31e5 * (1/n_e)**(0.5)  #[cm]
# print ('Electron inertial length de    = {:.2e}'.format(de),     '[cm]')

### Cooling time  
cool_rate = 1.7*10e-25*(Z**2)*T_e**(1./2) # normalized cooling rate
cooling_time_thin = 4.0e-36*A*(1+Z)*T_e/(Z*rho*cool_rate)
# print ('Cooling rate              = {:.2e}'.format(cool_rate))
# print ('Cooling time - thin       = {:.2e}'.format(cooling_time_thin*1e+9), '[ns]')

# print ('')

### beta thermic and beta dynamic
P_ther = (n_i*T_i+n_e*T_e)*K_Bol*1.16e4
P_dyn = rho*V**2

beta_ther_CGS = P_ther/(B**2.0/(8*np.pi))
print ('Beta thermic              = {:.1e}'.format(beta_ther_CGS)) # [P_ther/(B^2/(8*pi))]
beta_dyn_CGS = rho*V**2.0/(B**2.0/(8*np.pi)) 
print ('Beta dynamic              = {:.1e}'.format(beta_dyn_CGS)) # [P_dyn/(B^2/(8*pi))]

print ('Beta dynamic + thermal    = {:.1e}'.format(beta_dyn_CGS+beta_ther_CGS)) # [P_dyn/(B^2/(8*pi))]

P_ther_i = (n_i*T_i)*K_Bol*1.16e4
P_ther_e = (n_e*T_e)*K_Bol*1.16e4
beta_ther_i_CGS = P_ther_i/(B**2.0/(8*np.pi))
beta_ther_e_CGS = P_ther_e/(B**2.0/(8*np.pi))
# print ('Beta thermic of ion       = {:.1e}'.format(beta_ther_i_CGS)) # [P_ther/(B^2/(8*pi))]
# print ('Beta thermic of eon       = {:.1e}'.format(beta_ther_e_CGS)) # [P_ther/(B^2/(8*pi))]


mach_sonore = V/c_sound
print ('Mach number               = {:.2e}'.format(mach_sonore))
mach_alfvenique = V/c_alfven
print ('Alfven Mach number        = {:.2e}'.format(mach_alfvenique))    
mach_magnetosonic = V/c_magnetosonic
print ('Magnetosonic Mach number  = {:.2e}'.format(mach_magnetosonic))    

# print ('')

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
# print ('')

#sptizer diffusivity
Spitzer_condu = n_e*e_charge**2.0*tau_e/mass_electron
# print ('Spitzer conductivity = {:.1e}'.format(Spitzer_condu))
Spitzer_resi = 1.0/Spitzer_condu
# print ('Spitzer resistivity = {:.1e}'.format(Spitzer_resi))
mag_diff = Spitzer_resi/(4.*np.pi)
# print ('Magnetic diffusivity      = {:.2e}'.format(mag_diff*1e-4), '[m^2/s]')
mag_diff_time = space_scale**2/mag_diff
# print ('Magnetic diffusion time   = {:.2e}'.format(mag_diff_time*1e+9), '[ns]')     
# Magnetic Reynolds ### advection (induction) of B by plasma motion / diffusion of B in the plasma
Rm = V*space_scale/mag_diff  
print ('Magnetic Reynolds number  = {:.1e}'.format(Rm))

# print ('')
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

# print ('')



### Euler and Alfven numbers
Euler_CGS = (P_dyn/P_ther)**(1./2)
print ('Euler number              = {:.1e}'.format(Euler_CGS)) # [sqrt(P_dyn/P_ther)]

Alfven_CGS = B/(4*np.pi*P_ther)**(1./2)
# print 'Alfven number             = {:.2f}'.format(Alfven_CGS) # [B/sqrt(4pi*P)] 

# print ('')

## optical depth
g_bar = 1.2

od = 5e-38 * n_e * n_i * Z**2 * g_bar * space_scale * T_e**-3.5
# print ('Bremsstrahlung optical depth = {:.1e}'.format(od))




###############################################################################
###############################################################################
###############################################################################
### PRINT PRINT PRINT ### Just to round the values

# ### Printing input parameters: ### 
# print 'Flow velocity V  =', round(V*1e-5,2), '[km/s]'
# print 'B-field =', round(B), '[G]'
# print 'Electron temperature T_e =', round(T_e,1), '[eV]'
# print 'Ion temperature T_i =', round(T_i,1),  '[eV]'
# print 'Electron density n_e =', round(n_e,-10), '[cm-3]'
# print 'Ion density n_i =', round(n_i,-10), '[cm-3]'   
# print 'Density rho =', round(rho,13), '[g/cm-3]'
# print 'Charge state Z =', round(Z,3)
# print 'Mass number A =', round(A,2)
# print 'Spacial scale =', round(space_scale), '[cm]'
# print ''
# print 'Coulomb Logarithm (from Braginskii) =', round(cou_log_brag,2)
# print 'Thermal equilibration i-e =', round(Ther_equilib_ie_cou_Brag *1e+9,0), '[ns]'
# print 'Thermal equilibration e-i =', round(Ther_equilib_ei_cou_Brag*1e+9,0), '[ns]'
# print ''
# print 'Sound velocity =', round(c_sound*1e-5,2), '[km/s]'
# print 'Alfven velocity =', round(c_alfven*1e-5,2), '[km/s]'
# print ''
# print 'Electron thermal velocity =', round(vth_e*1e-5,2), '[km/s]'
# print 'Electron mean free path =', round(mfp_e,4), '[cm]'
# print 'Electron collision time =', round(tau_e*1e+9,2), '[ns]'
# print ''
# print 'Ion thermal velocity =', round(vth_i*1e-5,0), '[km/s]'
# print 'Ion mean free path =', round(mfp_i,4), '[cm]'
# print 'Ion collision time =', round(tau_i*1e+9,0), '[ns]'
# print ''
# print 'Electron Larmor radius =', round(gyro_e,5), '[cm]'
# print 'Electron gyrofrequency =', round(gyrofreq_e,-8), '[s-1]'
# print 'Electron magnetization =', round(omega_tau_e,2)
# print ''   
# print 'Ion Larmor radius =', round(gyro_i,3), '[cm]'
# print 'Ion gyrofrequency =', round(gyrofreq_i,-4), '[s-1]'
# print 'Ion magnetization =', round(omega_tau_i,2)
# print ''
# #print 'Cooling rate =', round(cool_rate,26)
# #print 'Cooling time - thin =', round(cooling_time_thin*1e+9,-5), '[ns]'
# #print ''
# print 'Mach number =', round(mach_sonore,2)
# print 'Alfven Mach number =', round(mach_sup_alfvenique,2)    
# print ''
# print 'Magnetic diffusion time =', round(mag_diff_time*1e+9,-22), '[ns]'     
# print 'Magnetic Reynolds number =', round(Rm,-10)
# print 'Reynolds number =', round(Re_i_B,2)
# print 'Peclet number =', round(Pe_e_B,-8)
# print ''
# print 'Beta thermic =', round(beta_ther_CGS,3)
# print 'Beta dynamic =', round(beta_dyn_CGS,3)
# print 'Euler number =', round(Euler_CGS,2)
# print 'Alfven number =', round(Alfven_CGS,2)

###############################################################################
###############################################################################
###############################################################################
