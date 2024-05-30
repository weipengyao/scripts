import yt
import matplotlib.pyplot  as plt
import numpy as np
import astropy

def get_yt_data_1d(filename):

    j_idx = 0
    k_idx = 0
    fields_list = ["x", "dens", "tele", "tion", "trad"]

    data_yt = yt.load(filename)
    data_yt_map = data_yt.covering_grid(
        level=0, left_edge=[0, 0.0, 0.0], dims=data_yt.domain_dimensions
    )

    data_dict = {}
    for f in fields_list:
        data_dict[f] = data_yt_map[f][:,j_idx,k_idx]
    return data_dict

def plot_dens_temps(filedir, time_now):
	fig, ax1 = plt.subplots()

	ax1.semilogy(data_dict['x'], data_dict['dens'], '-k', lw=1.0)
	ax1.set_ylabel(r'Density (g/cm$^3$)', color='k', fontsize=16)
	ax1.tick_params(axis='y',color='k', labelcolor='k')
	ax1.set_ylim([1e-7,1e-3])
	ax1.set_xlim([0.0,1.5])
	ax1.set_xlabel('x (cm)', fontsize=16)
	ax1.grid(which='both', alpha=0.5,linestyle=':', lw=1.0)

	ax2 = ax1.twinx()
	ax2.semilogy(data_dict['x'], data_dict['tele']/11600, '-.r', label='Te', lw=1.5)
	ax2.semilogy(data_dict['x'], data_dict['tion']/11600, '--b', label='Ti', lw=1.5)
	# ax2.semilogy(data_dict['x'], data_dict['trad']/11600, ':m', label='Tr', lw=1.0)
	ax2.set_ylabel('Temperatures (eV)', color='b', fontsize=16)
	ax2.tick_params(axis='y', color='b', labelcolor='b')
	ax2.spines['right'].set_color('b')
	ax2.spines['left'].set_color('k')
	ax2.set_ylim([1e0, 1e4])


	fig.legend(['dens','Te','Ti', 'Tr'], bbox_to_anchor=(0.85, 0.9),fontsize=12, fancybox=False, frameon=False)
	ax1.set_title('Time = {:.1f} ns'.format(time_now))

	fig.tight_layout()
	fig.savefig(filedir+'dens_temp_{:.1f}.png'.format(time_now),dpi=300)

# here, choose the data folder.
filedir = '/Users/yao/sshfs_mesopsl_data/job.65067/'

for i in range(101):
	i = i + 0
	filename = 'lasslab_hdf5_chk_'+str(i).zfill(4)
	data_dict = get_yt_data_1d(filedir+filename)
	time_now = i*0.5
	plot_dens_temps(filedir,time_now)

