import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from copy import copy


plt.ion()

# # Mancic, A., et al. "Isochoric heating of solids by laser-accelerated protons: Experimental characterization and self-consistent hydrodynamic modeling." High Energy Density Physics 6.1 (2010): 21-28.
def func(x, par):
	[a,b,c,d,e] = par
	return a + b*x + c*x**2 + d*x**3 +e*x**4
	
# There were some suspicious values in the tail of the spectra
df3=pd.read_csv('spectrum_noB_noGas.csv', sep=',')
x = np.array(df3['Energy [MeV]'])
y = np.array(df3['Nb protons [/MeV/sr]'])
z = copy(y)		# just intitialisation of the output array
	
maxE=20.4
maxE3=max(x)

# There is an error in coefficients in Mancic et al.. so use these 
# It is necessary to rescale it according to our maximum energy.
# Fitting the data from Mancic, Fig 3 in the folder "fitting"
par3 =[ 7.08445636e+00,  4.94911721e+00*maxE/maxE3, -6.00922396e-01*(maxE/maxE3)**2,  2.99448952e-02*(maxE/maxE3)**3, -6.47401643e-04*(maxE/maxE3)**4]


# own calculation
for i in range(len(z)):
	ang = func(x[i],par3)*np.pi/180.0
	z[i] = y[i]*(2*np.pi*(1-np.cos(ang)))

# Plot original input	/MeV/sr
plt.figure(5)
plt.clf()
plt.plot(x,y, '-', label='input')
plt.yscale('log')
plt.xlabel(r'$E$ [MeV]')
plt.xlim(0,)
plt.ylabel(r'Nb protons [/MeV/sr]')
plt.legend()
plt.savefig('energy_spectrum_per_Omega.png', dpi=200, transparent=False)	

# Plot output /MeV
plt.figure(6)
plt.clf()
plt.plot(x,z, '-', label='output')
plt.yscale('log')
plt.xlabel(r'$E$ [MeV]')
plt.ylabel(r'Nb protons [/MeV]')
plt.xlim(0,)
plt.legend()
plt.savefig('energy_spectrum.png', dpi=200, transparent=False)
print('Total number of protons above %0.1f MeV: %e'%(x[3], np.trapz(z[3:],x[3:])))

# write output
np.savetxt('output.csv', np.column_stack((x,z)), delimiter=",", header='"Energy [MeV]","Nb protons [/MeV]"', comments='')

