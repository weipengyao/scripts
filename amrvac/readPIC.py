import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec

fs=10 # -> Fontsize for plots

# Converting simulation units to physical units
mpCGS=1.67e-24                          # g
vA=0.05*3.0e10                          # cm/s
ni=0.01                                 # cm^-3
B0=vA*np.sqrt(4.0*np.pi*ni*mpCGS)*1.0e6 # uG
print('Field strength of the ordered field = %.2e uG' % B0)

# Function to load the simulation data
def load_data(sim_name):
    data=np.load("Results_%s/sim_%s.npz" % (sim_name, sim_name))
    x=data["x"]
    y=data["y"]
    b1_sim=data["b1"]
    b2_sim=data["b2"]
    b3_sim=data["b3"]

    b0=b1_sim[0,0,0]
    B1=b1_sim*B0/b0
    B2=b2_sim*B0/b0
    B3=b3_sim*B0/b0

    return x, y, B1, B2, B3

# Plot components of B-field 
def plot_compB(x, y, B1, B2, B3, step, i_comp):

    # Extent of simulation domains in unit of rg of simulated particles
    xmin=np.min(x)
    xmax=np.max(x)
    ymin=np.min(y)
    ymax=np.max(y)


    # Plot components of B-field
    fig, ax=plt.subplots(figsize=(10,4.5))
    if(i_comp==1):
        im=ax.imshow(B1[step,:,:].T, extent=[xmin, xmax, ymin, ymax], origin='lower', cmap='magma')#, vmin=-200, vmax=200)
        ax.set_title(r'$B_1(x,y)$', fontsize=fs)
    if(i_comp==2):
        im=ax.imshow(B2[step,:,:].T, extent=[xmin, xmax, ymin, ymax], origin='lower', cmap='magma')#, vmin=-200, vmax=200)
        ax.set_title(r'$B_2(x,y)$', fontsize=fs)
    if(i_comp==3):
        im=ax.imshow(B3[step,:,:].T, extent=[xmin, xmax, ymin, ymax], origin='lower', cmap='magma')#, vmin=-200, vmax=200)
        ax.set_title(r'$B_3(x,y)$', fontsize=fs)

    ax.set_xlabel(r'$x/r_{\rm CR}$', fontsize=fs)
    ax.set_ylabel(r'$y/r_{\rm CR}$', fontsize=fs)
    ax.tick_params(axis='x', labelsize=fs)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_xlim(-40,40)
    ax.set_ylim(-20,20)
    ax.set_title(r'$\Omega_0 t=%d$, %s' % (step*10, sim_name), fontsize=fs+2, x=0.01, y=0.98, ha='left')

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    plt.colorbar(im, cax=cax).ax.tick_params(labelsize=fs)    
    plt.savefig('Results_%s/fg_compB%d_%s_%d.png' % (sim_name, i_comp, sim_name, step))
    plt.close()

# Plot the turbulence spectrum
def plot_spectrum(x, y, B1, B2, B3, step):

    # Extent of simulation domains in unit of rg of simulated particles
    xmin=np.min(x)
    xmax=np.max(x)
    ymin=np.min(y)
    ymax=np.max(y)

    # Extract the B-field components at the time step considered
    b1=B1[step,:,:]
    b2=B2[step,:,:]
    b3=B3[step,:,:]

    Nx=np.unique(x).size
    Ny=np.unique(y).size

    kx=(np.fft.fftshift(np.fft.fftfreq(Nx,(xmax-xmin)/Nx)))*2.0*np.pi
    ky=(np.fft.fftshift(np.fft.fftfreq(Ny,(ymax-ymin)/Ny)))*2.0*np.pi
    kxmin=np.min(kx)
    kxmax=np.max(kx)
    kymin=np.min(ky)
    kymax=np.max(ky)

    # Apply 2D FFT
    b1_fft=np.fft.fft2(b1)
    b2_fft=np.fft.fft2(b2)
    b3_fft=np.fft.fft2(b3)

    # Shift the zero frequency component to the center
    b1_fft=np.fft.fftshift(b1_fft)
    b2_fft=np.fft.fftshift(b2_fft)
    b3_fft=np.fft.fftshift(b3_fft)

    # Compute the magnitude 
    P_11_mag=np.abs(b1_fft)**2
    P_22_mag=np.abs(b2_fft)**2
    P_33_mag=np.abs(b3_fft)**2

    # Create plot with the power spectra
    fig = plt.figure(figsize=(24, 18))
    gs = gridspec.GridSpec(3, 3, height_ratios=[1, 1, 0.8], width_ratios=[1, 1, 1])

    ax0 = fig.add_subplot(gs[0, 0])
    im = ax0.imshow(b1.T, extent=[xmin, xmax, ymin, ymax], origin='lower', cmap='magma')
    ax0.set_title(r'$B_1(x,y)$', fontsize=fs)
    ax0.set_xlabel(r'$x/r_{\rm CR}$', fontsize=fs)
    ax0.set_ylabel(r'$y/r_{\rm CR}$', fontsize=fs)
    ax0.tick_params(axis='x', labelsize=fs)
    ax0.tick_params(axis='y', labelsize=fs)
    divider = make_axes_locatable(ax0)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    plt.colorbar(im, cax=cax).ax.tick_params(labelsize=fs)

    ax1 = fig.add_subplot(gs[0, 1])
    im = ax1.imshow(b2.T, extent=[xmin, xmax, ymin, ymax], origin='lower', cmap='magma')
    ax1.set_title(r'$B_2(x,y)$', fontsize=fs)
    ax1.set_xlabel(r'$x/r_{\rm CR}$', fontsize=fs)
    ax1.set_ylabel(r'$y/r_{\rm CR}$', fontsize=fs)
    ax1.tick_params(axis='x', labelsize=fs)
    ax1.tick_params(axis='y', labelsize=fs)
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    plt.colorbar(im, cax=cax).ax.tick_params(labelsize=fs)

    ax2 = fig.add_subplot(gs[0, 2])
    im = ax2.imshow(b3.T, extent=[xmin, xmax, ymin, ymax], origin='lower', cmap='magma')
    ax2.set_title(r'$B_3(x,y)$', fontsize=fs)
    ax2.set_xlabel(r'$x$', fontsize=fs)
    ax2.set_ylabel(r'$y$', fontsize=fs)
    ax2.tick_params(axis='x', labelsize=fs)
    ax2.tick_params(axis='y', labelsize=fs)
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    plt.colorbar(im, cax=cax).ax.tick_params(labelsize=fs)

    ax3 = fig.add_subplot(gs[1, 0])
    im = ax3.imshow(np.log10(P_11_mag).T, extent=[kxmin, kxmax, kymin, kymax], origin='lower', cmap='magma')
    ax3.set_title(r'$P_{11}(k_x,k_y)$', fontsize=fs)
    ax3.set_xlabel(r'$k_x$', fontsize=fs)
    ax3.set_ylabel(r'$k_y$', fontsize=fs)
    ax3.tick_params(axis='x', labelsize=fs)
    ax3.tick_params(axis='y', labelsize=fs)
    divider = make_axes_locatable(ax3)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    plt.colorbar(im, cax=cax).ax.tick_params(labelsize=fs)

    ax4 = fig.add_subplot(gs[1, 1])
    im = ax4.imshow(np.log10(P_22_mag).T, extent=[kxmin, kxmax, kymin, kymax], origin='lower', cmap='magma')
    ax4.set_title(r'$P_{22}(k_x,k_y)$', fontsize=fs)
    ax4.set_xlabel(r'$k_x$', fontsize=fs)
    ax4.set_ylabel(r'$k_y$', fontsize=fs)
    ax4.tick_params(axis='x', labelsize=fs)
    ax4.tick_params(axis='y', labelsize=fs)
    divider = make_axes_locatable(ax4)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    plt.colorbar(im, cax=cax).ax.tick_params(labelsize=fs)

    ax5 = fig.add_subplot(gs[1, 2])
    im = ax5.imshow(np.log10(P_33_mag).T, extent=[kxmin, kxmax, kymin, kymax], origin='lower', cmap='magma')
    ax5.set_title(r'$P_{33}(k_x,k_y)$', fontsize=fs)
    ax5.set_xlabel(r'$k_x$', fontsize=fs)
    ax5.set_ylabel(r'$k_y$', fontsize=fs)
    ax5.tick_params(axis='x', labelsize=fs)
    ax5.tick_params(axis='y', labelsize=fs)
    divider = make_axes_locatable(ax5)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    plt.colorbar(im, cax=cax).ax.tick_params(labelsize=fs)

    ax6 = fig.add_subplot(gs[2, 1])  # Center-left
    for i in range((Ny // 2) - 10, (Ny // 2) + 10):
        ax6.plot(kx, P_22_mag[:, i], linestyle='-', color='grey', linewidth=1)
    ax6.plot(kx, np.mean(P_22_mag, axis=1), 'g-', linewidth=3, label=r'$P_{22}(k_x,k_y={\rm const})$')
    ax6.set_xscale('log')
    ax6.set_yscale('log')
    ax6.set_ylim(1.0e-4 * np.max(np.mean(P_22_mag, axis=1)), 1.0e1 * np.max(np.mean(P_22_mag, axis=1)))
    ax6.set_xlabel(r'$k_x$', fontsize=fs)
    ax6.set_ylabel(r'$P_{22}(k_x,k_y={\rm const})$', fontsize=fs)
    ax6.tick_params(axis='x', labelsize=fs)
    ax6.tick_params(axis='y', labelsize=fs)
    ax6.legend(loc='lower right', prop={"size": fs})
    ax6.grid(linestyle='--')
    ax6.set_title(r'$P_{22}$', fontsize=fs)

    ax7 = fig.add_subplot(gs[2, 2])  # Center-right
    for i in range((Ny // 2) - 10, (Ny // 2) + 10):
        ax7.plot(kx, P_33_mag[:, i], linestyle='-', color='grey', linewidth=1)
    ax7.plot(kx, np.mean(P_33_mag, axis=1), 'g-', linewidth=3, label=r'$P_{33}(k_x,k_y={\rm const})$')
    ax7.set_xscale('log')
    ax7.set_yscale('log')
    ax7.set_ylim(1.0e-4 * np.max(np.mean(P_22_mag, axis=1)), 1.0e1 * np.max(np.mean(P_22_mag, axis=1)))  # Same limits
    ax7.set_xlabel(r'$k_x$', fontsize=fs)
    ax7.set_ylabel(r'$P_{33}(k_x,k_y={\rm const})$', fontsize=fs)
    ax7.tick_params(axis='x', labelsize=fs)
    ax7.tick_params(axis='y', labelsize=fs)
    ax7.legend(loc='lower right', prop={"size": fs})
    ax7.grid(linestyle='--')
    ax7.set_title(r'$P_{33}$', fontsize=fs)

    fig.suptitle(r'$n_{\rm CR}/n_i=%s, \Omega_0 t=%.1f$' % (sim_name, (step * 10.0)), fontsize=fs+2, x=0.01, y=0.98, ha='left')
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig('Results_%s/fg_spectrum_%s_%04d.png' % (sim_name, sim_name, step))

sim_name='2P_2GeV_dt2_260_iso_xmax100'

x, y, B1, B2, B3=load_data(sim_name)

Nt, Nx, Ny=B1.shape
print("Shape of the B-field cubes for 3 components: number of time steps = %d, size of 2D_grid = (%d, %d)" % (Nt, Nx, Ny))

plot_compB(x, y, B1, B2, B3, 31, 2)

