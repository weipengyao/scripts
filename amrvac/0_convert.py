import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec

fs=22

sim_name='2P_2GeV_dt2_260_iso_xmax100'

my_N_step=np.linspace(0,95,96)#np.array([40])#[10, 20, 30, 40, 50, 60, 70, 80, 90, 100]) #, 80])

b1_array=[]
b2_array=[]
b3_array=[]

for step in (my_N_step):

    x, y, rho, b1, b2, b3, jp1, up1=np.loadtxt('XYZ_%s/visit_ex_db_%04d.xyz' % (sim_name,step),unpack=True,skiprows=2,usecols=[1,2,4,5,6,7,8,9])

    # Number of grid points along x and y
    Nx=np.unique(x).size
    Ny=np.unique(y).size

    # Size of simulation domains in unit of rg of simulated particles
    xmin=np.min(x)
    xmax=np.max(x)
    ymin=np.min(y)
    ymax=np.max(y)

    # Reformat the data file, by increasing x, then by y within each block of constant x
    xyz_combined=np.vstack((x, y, rho, b1, b2, b3, jp1, up1)).T
    xyz_sorted=xyz_combined[np.lexsort((xyz_combined[:, 1], xyz_combined[:, 0]))]
    unique_xyz=np.unique(xyz_sorted, axis=0)
    x, y, rho, b1, b2, b3, jp1, up1=unique_xyz.T
    x=x.reshape(Nx,Ny)
    y=y.reshape(Nx,Ny)
    rho=rho.reshape(Nx,Ny)
    b1=b1.reshape(Nx,Ny)
    b2=b2.reshape(Nx,Ny)
    b3=b3.reshape(Nx,Ny)
    jp1=jp1.reshape(Nx,Ny)
    up1=up1.reshape(Nx,Ny)

    # Append current step B1, B2, B3
    b1_array.append(b1)
    b2_array.append(b2)
    b3_array.append(b3)

    x_array=x
    y_array=y

# Convert lists to 3D arrays (time, Nx, Ny)
b1_array = np.array(b1_array)
b2_array = np.array(b2_array)
b3_array = np.array(b3_array)

# Save to .npz file
np.savez("Results_%s/sim_%s.npz" % (sim_name, sim_name), x=x_array, y=y_array, b1=b1_array, b2=b2_array, b3=b3_array)