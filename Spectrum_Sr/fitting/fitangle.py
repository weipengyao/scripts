import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from copy import copy
from scipy.optimize import curve_fit

def func(x, a,b,c,d,e):
	return a + b*x + c*x**2 + d*x**3 +e*x**4

plt.ion()
plt.figure(4)
plt.clf()
df=pd.read_csv('datafit.csv', sep=',')
plt.plot(df['E'], df['a'], 'o')
plt.xlabel('E [MeV]')
plt.ylabel('half-angle [deg]')
popt, pcov = curve_fit(func, df['E'], df['a'])

xfit = np.linspace(0,np.max(df['E']))
yfit = copy(xfit)
for i in range(len(xfit)):
	yfit[i] = func(xfit[i], popt[0], popt[1], popt[2], popt[3], popt[4])
plt.plot(xfit,yfit)
plt.xlim(0,)
plt.ylim(0,)


