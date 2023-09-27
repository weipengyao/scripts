#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 23:17:58 2023

@author: yao
"""

import happi

wkdir = [
        # '/Volumes/LaCie/collisional_shock/t4_ni1e19_v150_T50_dx02_N3_n_ionize_n_coll/',
         # '/Volumes/LaCie/collisional_shock/t4_ni1e19_v150_T50_dx02_N3_n_ionize_w_ii_coll/',
          # '/Volumes/LaCie/collisional_shock/t5/t5_ni1e19_v150_T50_dx02_N3_w_coll3_r0/',
          # '/Users/yao/Desktop/t5/t5_ni1e19_v150_T50_dx02_N3_n_coll4/',
          # '/Users/yao/Desktop/t5/t5_ni1e19_v150_T50_dx02_N3_w_coll4_r0/'
          # '/Users/yao/Desktop/t6/t6_ni1e19_v150_T50_dx01_N3_w_coll/',
          # '/Users/yao/Desktop/t6/t6_ni1e19_v150_T50_dx02_N3_w_coll/',
          # '/Users/yao/Desktop/t6/t6_ni1e19_v150_T50_dx01_N7_w_coll_r0/',
          # '/Users/yao/Desktop/t6/t6_ni4.3e18_v150_T50_dx01_N7_w_coll_r0/',
          '/Volumes/LaCie/collisional_shock/t7_2/t7_ne1e19_v150_T50_dx01_N5_n_coll_r1/',
          '/Volumes/LaCie/collisional_shock/t7_2/t7_ne1e19_v150_T50_dx01_N5_w_coll_r3/',
         ]


S0 = happi.Open(wkdir[0], reference_angular_frequency_SI=177788752730887.47)
S1 = happi.Open(wkdir[1], reference_angular_frequency_SI=177788752730887.47)
# S2 = happi.Open(wkdir[2], reference_angular_frequency_SI=56375055300167.87)

#%%

# prepare constants, units

me = 9.1e-31
mp = 1836.*me
qe = 1.6e-19
ep = 8.9e-12  # epsilon_0
c  = 3.0e8
wr = S0.namelist.w_r
de = c / wr
Lx = S0.namelist.L_x.real * de * 1e3      # in mm
dx = S0.namelist.d_x * de * 1e3           # in mm

Te = S0.namelist.T_e * 511.e3             # in eV
ne = 7.0e19                              # in cm-3
ld = 7.43e2 * Te**0.5 * ne**(-0.5) * 10. # in mm
dt = S0.namelist.d_t

B0 = S0.namelist.B_z * (me * wr / qe)
wc = qe * B0 / me

#%%

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

jetcmap = plt.cm.get_cmap("jet", 9) #generate a jet map with 10 values "rainbow", "jet"
jet_vals = jetcmap(np.arange(9)) #extract those values as an array 
jet_vals[0] = [1.0, 1, 1.0, 1] #change the first value 
newcmap = mpl.colors.LinearSegmentedColormap.from_list("mine", jet_vals) 

#%%

Ex0 = S0.Field(0,'Ex', units=['mm','ns','cm^-3','km/s','V/m'],
               label='w/o. coll.')

Ex1 = S1.Field(0,'Ex', units=['mm','ns','cm^-3','km/s','V/m'],
               label='with coll.')
                              
# Ex2 = S2.Field(0,'Ex', units=['mm','ns','cm^-3','km/s','V/m'],
               # label='3')

happi.multiSlide(Ex0,Ex1)#,Ex2)                              
                             
                              

#%%

xpx00 = S0.ParticleBinning("#0+#1",units=['mm','ns','cm^-3','km/s'],
                          data_log=True,
                           vmin=-3,
                           vmax=3,
                           # xmin=4,
                           # xmax=7,
                          cmap=newcmap)

xpx01 = S1.ParticleBinning("#0+#1",units=['mm','ns','cm^-3','km/s'],
                          data_log=True,
                            vmin=-3,
                            vmax=3,
                           # xmin=4,
                           # xmax=7,
                          cmap=newcmap
                          )

# xpx02 = S2.ParticleBinning("#0+#1",units=['mm','ns','cm^-3','km/s'],
#                           data_log=True,
#                            vmin=2,
#                            vmax=5,
#                            # xmin=4,
#                            # xmax=7,
#                           cmap=newcmap)#.slide()

happi.multiSlide(xpx00,xpx01,
                 # xpx02, 
                 shape=[2,1])

#%% eon x-px

xpx00 = S0.ParticleBinning("#2+#3",units=['mm','ns','cm^-3','km/s'],
                          data_log=True,
                           vmin=-3,
                           vmax=3,
                           # xmin=4,
                           # xmax=7,
                          cmap=newcmap)

xpx01 = S1.ParticleBinning("#2+#3",units=['mm','ns','cm^-3','km/s'],
                          data_log=True,
                           vmin=-3,
                           vmax=3,
                           # xmin=4,
                           # xmax=7,
                          cmap=newcmap)#.slide()

# xpx02 = S2.ParticleBinning("#1",units=['mm','ns','cm^-3','km/s'],
#                           data_log=True,
#                            vmin=2,
#                            vmax=5,
#                            # xmin=4,
#                            # xmax=7,
#                           cmap=newcmap)#.slide()

happi.multiSlide(xpx00,xpx01,
                 # xpx02, 
                  shape=[2,1]
                 )


#%% eon number density

xpx00 = S0.ParticleBinning("#6+#7",
                           units=['mm','ns','cm^-3','km/s'],
                           label='w/o coll.',
                           xmin=0.4,
                           xmax=1.4,
                          # data_log=True,
                            vmin=1.0e19,
                            vmax=3.5e19,
                           # xmin=4,
                           # xmax=7,
                          # cmap=newcmap
                          )

xpx01 = S1.ParticleBinning("#6+#7",
                           units=['mm','ns','cm^-3','km/s'],
                           label='with coll.',
                           xmin=0.4,
                           xmax=1.4,
                            vmin=1.0e19,
                            vmax=3.5e19,
                          # data_log=True,
                            # vmin=-3,
                           # vmax=3,
                           # xmin=4,
                           # xmax=7,
                          # cmap=newcmap
                          )#.slide()

# xpx02 = S2.ParticleBinning("#1",units=['mm','ns','cm^-3','km/s'],
#                           data_log=True,
#                            vmin=2,
#                            vmax=5,
#                            # xmin=4,
#                            # xmax=7,
#                           cmap=newcmap)#.slide()

happi.multiSlide(xpx00,xpx01,
                 # xpx02, 
                 # shape=[2,1]
                 )
#%% ion number density

xpx00 = S0.ParticleBinning("#10+#11",
                           units=['mm','ns','cm^-3','km/s'],
                           label='w/o coll.',
                           xmin=0.4,
                           xmax=1.4,
                          # data_log=True,
                            # vmin=1.0e19,
                            # vmax=2.5e19,
                           # xmin=4,
                           # xmax=7,
                          # cmap=newcmap
                          )

xpx01 = S1.ParticleBinning("#10+#11",
                           units=['mm','ns','cm^-3','km/s'],
                           label='with coll.',
                           xmin=0.4,
                           xmax=1.4,
                            # vmin=1.0e19,
                            # vmax=2.5e19,
                          # data_log=True,
                            # vmin=-3,
                           # vmax=3,
                           # xmin=4,
                           # xmax=7,
                          # cmap=newcmap
                          )#.slide()

# xpx02 = S2.ParticleBinning("#1",units=['mm','ns','cm^-3','km/s'],
#                           data_log=True,
#                            vmin=2,
#                            vmax=5,
#                            # xmin=4,
#                            # xmax=7,
#                           cmap=newcmap)#.slide()

happi.multiSlide(xpx00,xpx01,
                 # xpx02, 
                 # shape=[2,1]
                 )
#%%

xpx00 = S0.ParticleBinning("#2",units=['mm','ns','cm^-3','km/s'],
                          data_log=True,
                           vmin=2,
                           vmax=5,
                           # xmin=4,
                           # xmax=7,
                          cmap=newcmap)

xpx01 = S1.ParticleBinning("#2",units=['mm','ns','cm^-3','km/s'],
                          data_log=True,
                           vmin=2,
                           vmax=5,
                           # xmin=4,
                           # xmax=7,
                          cmap=newcmap)#.slide()

xpx02 = S2.ParticleBinning("#2",units=['mm','ns','cm^-3','km/s'],
                          data_log=True,
                           vmin=2,
                           vmax=5,
                           # xmin=4,
                           # xmax=7,
                          cmap=newcmap)#.slide()

happi.multiSlide(xpx00,xpx01,xpx02, shape=[3,1])

#%%

S2.ParticleBinning("#4+#5",units=['mm','ns','cm^-3','km/s'],
                          # data_log=True,
                           # vmin=2,
                           # vmax=5,
                           # xmin=4,
                           # xmax=7,
                          cmap=newcmap).slide()

#%%

Te00 = S0.ParticleBinning('(#4+#5)/(#6+#7)',
                         units=['mm','ns','cm^-3','km/s','eV'],
                         color='blue',
                         linestyle='-',
                          vmin=40,
                          vmax=90,
                            xmin=0.4,
                            xmax=1.4,
                          # vmax=1200.,
                         label='w/o coll.')

Te01 = S1.ParticleBinning('(#4+#5)/(#6+#7)',
                         units=['mm','ns','cm^-3','km/s','eV'],
                         color='black',
                         linestyle='-.',
                          vmin=40,
                          vmax=90,
                            xmin=0.4,
                            xmax=1.4,
                          # vmax=1200.,
                         label='with coll.')

# Te02 = S2.ParticleBinning('(#4+#5)/(#6+#7)',
#                          units=['mm','ns','cm^-3','km/s','eV'],
#                          color='red',
#                          linestyle='--',
#                           vmin=50,
#                           vmax=120,
#                             # xmin=0,
#                             # xmax=2.66,
#                           # vmax=1200.,
#                          label='ni=4.3e18, Z=7')

happi.multiSlide(Te00, Te01)#,Te02)
#%%

Te00 = S0.ParticleBinning('(#5)',
                         units=['mm','ns','cm^-3','km/s','eV'],
                         color='blue',
                         linestyle='-',
                          # vmin=50,
                          # vmax=120,
                            # xmin=0,
                            # xmax=2.66,
                          # vmax=1200.,
                         label='ni=1e19, Z=3')

Te01 = S1.ParticleBinning('(#5)',
                         units=['mm','ns','cm^-3','km/s','eV'],
                         color='black',
                         linestyle='-.',
                          # vmin=50,
                          # vmax=120,
                            # xmin=0,
                            # xmax=2.66,
                          # vmax=1200.,
                         label='ni=1e19, Z=7')

Te02 = S2.ParticleBinning('(#5)',
                         units=['mm','ns','cm^-3','km/s','eV'],
                         color='red',
                         linestyle='--',
                          # vmin=50,
                          # vmax=120,
                            # xmin=0,
                            # xmax=2.66,
                          # vmax=1200.,
                         label='ni=4.3e18, Z=7')

happi.multiSlide(Te00, Te01,Te02)

#%%

S2.Field(0,'Ex',units=["V/m"]).slide()

#%%

Te00 = S0.ParticleBinning('(#7)',
                         units=['mm','ns','cm^-3','km/s','eV'],
                         color='blue',
                         linestyle='-',
                          # vmin=50,
                          # vmax=120,
                            # xmin=0,
                            # xmax=2.66,
                          # vmax=1200.,
                         label='ni=1e19, Z=3')

Te01 = S1.ParticleBinning('(#7)',
                         units=['mm','ns','cm^-3','km/s','eV'],
                         color='black',
                         linestyle='-.',
                          # vmin=50,
                          # vmax=120,
                            # xmin=0,
                            # xmax=2.66,
                          # vmax=1200.,
                         label='ni=1e19, Z=7')

Te02 = S2.ParticleBinning('(#7)',
                         units=['mm','ns','cm^-3','km/s','eV'],
                         color='red',
                         linestyle='--',
                          # vmin=50,
                          # vmax=120,
                            # xmin=0,
                            # xmax=2.66,
                          # vmax=1200.,
                         label='ni=4.3e18, Z=7')

happi.multiSlide(Te00, Te01,Te02)

#%%

Te00 = S0.ParticleBinning('(#4)/(#6)',
                         units=['mm','ns','cm^-3','km/s','eV'],
                         color='blue',
                         linestyle='-',
                          vmin=50,
                          vmax=120,
                            # xmin=0,
                            # xmax=2.66,
                          # vmax=1200.,
                         label='ni=1e19, Z=3')

Te01 = S1.ParticleBinning('(#4)/(#6)',
                         units=['mm','ns','cm^-3','km/s','eV'],
                         color='black',
                         linestyle='-.',
                          vmin=50,
                          vmax=120,
                            # xmin=0,
                            # xmax=2.66,
                          # vmax=1200.,
                         label='ni=1e19, Z=7')

Te02 = S2.ParticleBinning('(#4)/(#6)',
                         units=['mm','ns','cm^-3','km/s','eV'],
                         color='red',
                         linestyle='--',
                          vmin=50,
                          vmax=120,
                            # xmin=0,
                            # xmax=2.66,
                          # vmax=1200.,
                         label='ni=4.3e18, Z=7')

happi.multiSlide(Te00, Te01,Te02)

#%%


xpx0 = S0.ParticleBinning(0,units=['mm','ns','cm^-3','km/s'],
                          data_log=True,
                           vmin=2,
                           vmax=5,
                           # xmin=4,
                           # xmax=7,
                          cmap=newcmap)#.slide()

xpx1 = S1.ParticleBinning(0,units=['mm','ns','cm^-3','km/s'],
                          data_log=True,
                            vmin=2,
                            vmax=5,
                            # xmin=4,
                            # xmax=7,
                          cmap=newcmap)


happi.multiSlide(xpx0,xpx1,shape=[2,1])
# happi.multiPlot(xpx0,xpx1,shape=[2,1],
#                 skipAnimation=True,
#                 # movie='/Users/yao/Desktop/test_xpx.gif'
#                 )

#%%

ni0 = S0.ParticleBinning("#5*3",
                       units=['ns','km/s','mm','cm^-3'],
                       color='blue',
                       linestyle='-',
                       label='ni w/o coll.',
                       # xmin=4,
                       # xmax=7,
                       # vmin=0,
                       # vmax=3.5e19,
                       )#.slide()
#%%
happi.multiSlide(ni0,ne0)

#%%
ni1 = S0.ParticleBinning(5,
                       units=['mm','ns','cm^-3','km/s'],
                       color='blue',
                       linestyle='--',
                       label='ni no coll.',
                       # xmin=4,
                       # xmax=7,
                       # vmin=0,
                       # vmax=3.5e19,
                       )#.slide()

ni2 = S1.ParticleBinning(5,
                        units=['mm','ns','cm^-3','km/s'],
                        color='red',
                        linestyle='-',
                        label='ni with coll.',
                        # xmin=4,
                        # xmax=7,
                        # vmin=0,
                        # vmax=3.5e19,
                        )#.slide()

happi.multiSlide(ni1,ni2)
# happi.multiPlot(ni1,ni2,
#                 # ne0,ne1,
#                 skipAnimation=True,)

#%%

te0 = S0.ParticleBinning('(#2)/(#3)',
                         units=['mm','ns','cm^-3','km/s','eV'],
                         color='cyan',
                         linestyle=':',
                           # xmin=1,
                           # xmax=4,
                         label='Te w/o coll.')#.slide()

ti0 = S0.ParticleBinning('(#4/#6)+(#5/#7)',
                         units=['mm','ns','cm^-3','km/s','eV'],
                         color='blue',
                         linestyle='-',
                           # xmin=1,
                           # xmax=4,
                          # vmax=1200.,
                         label='Ti w/o coll.')#.slide()

te1 = S1.ParticleBinning('(#2)/(#3)',
                         units=['mm','ns','cm^-3','km/s','eV'],
                         color='purple',
                         linestyle=':',
                           # xmin=1,
                           # xmax=4,
                         label='Te with coll.')#.slide()

ti1 = S1.ParticleBinning('(#4)/(#5)',
                         units=['mm','ns','cm^-3','km/s','eV'],
                         color='red',
                         linestyle='-',
                           # xmin=1,
                           # xmax=4,
                         label='Ti with coll.')#.slide()

happi.multiSlide(ti0, ti1, te0, te1)

# happi.multiPlot(ti0, ti1, te0, te1,skipAnimation=True)

#%%

te0 = S0.ParticleBinning('511*1e3*(#9)/(#3)',
                          units=['mm','ns'],
                         color='cyan',
                         linestyle=':',
                           xmin=1,
                           xmax=4,
                         label='Te w/o coll.')#.slide()

ti0 = S0.ParticleBinning('14*1836*511*1e3*(#12)/(#5)',
                         units=['mm','ns'],
                         color='blue',
                         linestyle='-',
                           xmin=1,
                           xmax=4,
                          # vmax=1200.,
                         label='Ti w/o coll.')#.slide()

te1 = S1.ParticleBinning('511*1e3*(#9)/(#3)',
                         units=['mm','ns'],
                         color='purple',
                         linestyle=':',
                           xmin=1,
                           xmax=4,
                         label='Te with coll.')#.slide()

ti1 = S1.ParticleBinning('14*1836*511*1e3*(#12)/(#5)',
                         units=['mm','ns'],
                         color='red',
                         linestyle='-',
                           xmin=1,
                           xmax=4,
                         label='Ti with coll.')#.slide()

happi.multiSlide(ti0, ti1, te0, te1)

#%%

te0 = S0.ParticleBinning('511*1e3*(#10)/(#3)',
                          units=['mm','ns'],
                         color='cyan',
                         linestyle=':',
                           xmin=1,
                           xmax=4,
                         label='Te w/o coll.')#.slide()

ti0 = S0.ParticleBinning('14*1836*511*1e3*(#13)/(#5)',
                         units=['mm','ns'],
                         color='blue',
                         linestyle='-',
                           xmin=1,
                           xmax=4,
                          # vmax=1200.,
                         label='Ti w/o coll.')#.slide()

te1 = S1.ParticleBinning('511*1e3*(#10)/(#3)',
                         units=['mm','ns'],
                         color='purple',
                         linestyle=':',
                           xmin=1,
                           xmax=4,
                         label='Te with coll.')#.slide()

ti1 = S1.ParticleBinning('14*1836*511*1e3*(#13)/(#5)',
                         units=['mm','ns'],
                         color='red',
                         linestyle='-',
                           xmin=1,
                           xmax=4,
                         label='Ti with coll.')#.slide()

happi.multiSlide(ti0, ti1, te0, te1)

#%%

DATASMILEI_Eon_Ty_ = S0.ParticleBinning(9, units=['mm','ns'])
DATASMILEI_Eon_Ty1 = np.array(DATASMILEI_Eon_Ty_.getData())

DATASMILEI_Ion_Ty_ = S0.ParticleBinning(12, units=['mm','ns'])
DATASMILEI_Ion_Ty1 = np.array(DATASMILEI_Ion_Ty_.getData())

DATASMILEI_Eon_Ne_ = S0.ParticleBinning(3, units=['mm','ns'])
DATASMILEI_Eon_Ne1 = np.array(DATASMILEI_Eon_Ne_.getData())

DATASMILEI_Ion_Ni_ = S0.ParticleBinning(5, units=['mm','ns'])
DATASMILEI_Ion_Ni1 = np.array(DATASMILEI_Ion_Ni_.getData())

Nx = (np.shape(DATASMILEI_Eon_Ty1))[1] #the number of discrete x-values
Nout = (np.shape(DATASMILEI_Eon_Ty1))[0]# the number of output

SMILEI_Eon_Ty1 = np.zeros((Nout,Nx))
SMILEI_Ion_Ty1 = np.zeros((Nout,Nx))
#x_dat = np.zeros(Nx)
x_dat = np.array(DATASMILEI_Eon_Ty_.getAxis("x"))



#%%

it =65
every = 1

SMILEI_Eon_Ty1[it,:] = 511*1e3*(DATASMILEI_Eon_Ty1[it*every,:])/DATASMILEI_Eon_Ne1[it,:] #returns the expectation value for the temperature
SMILEI_Ion_Ty1[it,:] = 14*1836*511*1e3*(DATASMILEI_Ion_Ty1[it*every,:])/DATASMILEI_Ion_Ni1[it,:]

plt.figure(it)
plt.xlim(0., 1)
#plt.ylim(0., 200)

plt.plot(x_dat , SMILEI_Eon_Ty1[it,:], label = "Electrons")
plt.plot(x_dat , SMILEI_Ion_Ty1[it,:], label = "Ions")
plt.xlabel('Length (mm)')
plt.legend(loc=1)
plt.ylabel('$\\rm Temperature \, (eV)$')
plt.xlim([1,4])
plt.tight_layout()
#plt.savefig("Temperature_vs_length_"+str(i)+".png", dpi = 300)
plt.show()


#%%

xpx0 = S0.ParticleBinning(6,units=['mm','ns','cm^-3','km/s'],
                          data_log=True,
                            vmin=5,
                            vmax=8,
                          cmap=newcmap,
                          aspect='equal')

xpx1 = S1.ParticleBinning(6,units=['mm','ns','cm^-3','km/s'],
                          data_log=True,
                            vmin=5,
                            vmax=8,
                          cmap=newcmap,
                          aspect='equal')

happi.multiSlide(xpx0,xpx1,shape=[2,1])
# happi.multiPlot(xpx0,xpx1,shape=[2,1],
#                 skipAnimation=True,
#                 figsize=[4,7],
#                 # movie='/Users/yao/Desktop/test_xpx.gif'
#                 )


#%%
S0.Field(0,'Ex',units=['V/m'],
         vsym=True).slide()





#%%

ne0 = S0.ParticleBinning(3,
                       units=['mm','ns','cm^-3','km/s'],
                       color='cyan',
                       linestyle='-',
                       linewidth=2.0,
                       label='ne w/o coll.'
                       )#.slide()
ne1 = S1.ParticleBinning(3,
                       units=['mm','ns','cm^-3','km/s'],
                       color='purple',
                       linestyle='--',
                       linewidth=1.0,
                       label='ne with coll.'
                       )#.slide()

happi.multiSlide(ne0,ni0)
# happi.multiPlot(ne0,ne1,
#                 skipAnimation=True,)

#%%

es0 = S0.ParticleBinning(1,
                   units=['mm','ns','cm^-3','km/s','keV'],
                   sum={'x':'all'},
                   data_log=True,
                   )#.slide()

es1 = S1.ParticleBinning(1,
                   units=['mm','ns','cm^-3','km/s','keV'],
                   sum={'x':'all'},
                   data_log=True,
                   )#.slide()

happi.multiSlide(es0,es1)

#%%

tel0 = S0.ParticleBinning('#6/#7',
                         units=['mm','ns','cm^-3','km/s','eV'],
                         color='blue',
                         linestyle='-',
                         label='Te left w/o coll.')#.slide()

ter0 = S0.ParticleBinning('#8/#9',
                         units=['mm','ns','cm^-3','km/s','eV'],
                         color='cyan',
                         linestyle='-',
                         label='Te right w/o coll.')#.slide()

til0 = S0.ParticleBinning('#10/#11',
                         units=['mm','ns','cm^-3','km/s','eV'],
                         color='red',
                         linestyle='-',
                         label='Ti left w/o coll.')#.slide()

tir0 = S0.ParticleBinning('#12/#13',
                         units=['mm','ns','cm^-3','km/s','eV'],
                         color='purple',
                         linestyle='-',
                         label='Ti right w/o coll.')#.slide()

happi.multiSlide(tel0,ter0)


#%%

te0 = S0.ParticleBinning('(#6+#8)/(#7+#9)',
                         units=['mm','ns','cm^-3','km/s','eV'],
                         color='blue',
                         linestyle='-',
                         label='Te w/o coll.')#.slide()

ti0 = S0.ParticleBinning('(#10+#12)/(#11+#13)',
                         units=['mm','ns','cm^-3','km/s','eV'],
                         color='red',
                         linestyle='-',
                         label='Ti w/o coll.').slide()

#%%

te0 = S0.ParticleBinning('(#8)/(#9)',
                         units=['mm','ns','cm^-3','km/s','eV'],
                         color='cyan',
                         linestyle=':',
                          xmin=Lx/2.,
                          vmax=1200.,
                         label='Te w/o coll.')#.slide()

ti0 = S0.ParticleBinning('(#12)/(#13)',
                         units=['mm','ns','cm^-3','km/s','eV'],
                         color='blue',
                         linestyle='-',
                          xmin=Lx/2.,
                          vmax=1200.,
                         label='Ti w/o coll.')#.slide()

te1 = S1.ParticleBinning('(#8)/(#9)',
                         units=['mm','ns','cm^-3','km/s','eV'],
                         color='purple',
                         linestyle=':',
                          xmin=Lx/2.,
                          vmax=1200.,
                         label='Te with coll.')#.slide()

ti1 = S1.ParticleBinning('(#12)/(#13)',
                         units=['mm','ns','cm^-3','km/s','eV'],
                         color='red',
                         linestyle='-',
                          xmin=Lx/2.,
                          vmax=1200.,
                         label='Ti with coll.')#.slide()

happi.multiSlide(ti0, ti1, te0, te1)
# happi.multiPlot(ti0, ti1, te0, te1,
#                 skipAnimation=True,)


#%%


import h5py


f = h5py.File(wkdir[1]+"BinaryProcesses0.h5", "r")

for key in f.keys():
    print(key)
    
every = 101010 # adjust from the above results


for it in range(30): # I put 10, but this should be the number of output that you have 
    it+=1
    data = f['t'+str(it*every).zfill(8)] # You might need to change this 
    s_param = np.array(data['s'])
    print(it*every,"  ",s_param.max())
    
    #%%