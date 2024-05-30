import yt
import matplotlib.pyplot  as plt
import numpy as np
import astropy

def get_yt_data_1d(filename):

    j_idx = 0
    k_idx = 0
    fields_list = ["x", "dens", "tele", "tion", "trad", "velx", "ye"]

    data_yt = yt.load(filename)
    data_yt_map = data_yt.covering_grid(
        level=0, left_edge=[0, 0.0, 0.0], dims=data_yt.domain_dimensions
    )

    data_dict = {}
    for f in fields_list:
        data_dict[f] = data_yt_map[f][:,j_idx,k_idx]
    return data_dict

def plot_variables(filedir, time_now):
	fig, axs = plt.subplots(2, 1, sharex=True)

	axs[0].semilogy(data_dict['x'], data_dict['dens'], '-k', lw=1.0)
	axs[0].set_ylabel(r'Density (g/cm$^3$)', color='k', fontsize=16)
	axs[0].tick_params(axis='y',color='k', labelcolor='k')
	axs[0].set_ylim([2e-6,2e-4])
	axs[0].set_xlim([0.0,1.5])
	# axs[0].set_xlabel('x (cm)', fontsize=16)
	axs[0].grid(which='both', alpha=0.5,linestyle=':', lw=0.5)

	ax2 = axs[0].twinx()
	ax2.semilogy(data_dict['x'], data_dict['tele']/11600, '-.r', label='Te', lw=1.5)
	ax2.semilogy(data_dict['x'], data_dict['tion']/11600, '--b', label='Ti', lw=1.5)
	ax2.set_ylabel('Temperatures (eV)', color='b', fontsize=16)
	ax2.tick_params(axis='y', color='b', labelcolor='b')
	ax2.spines['right'].set_color('b')
	ax2.spines['left'].set_color('k')
	ax2.set_ylim([2e1, 2e3])


	fig.legend(['dens','Te','Ti'], bbox_to_anchor=(0.8, 0.9),fontsize=12, fancybox=False, frameon=False)
	axs[0].set_title('Time = {:.1f} ns'.format(time_now))

	axs[1].plot(data_dict['x'], data_dict['velx'], '-c', lw=1.0)
	axs[1].set_ylabel(r'Velocity (cm/s)', color='c', fontsize=16)
	axs[1].tick_params(axis='y',color='c', labelcolor='c')
	axs[1].set_ylim([0,3.0e7])
	axs[1].set_xlim([0.0,1.5])
	axs[1].set_xlabel('x (cm)', fontsize=16)
	axs[1].grid(which='both', alpha=0.5,linestyle=':', lw=1.0)

	ax3 = axs[1].twinx()
	ax3.plot(data_dict['x'], data_dict['ye']*14, '-m', lw=1.5)
	ax3.tick_params(axis='y', color='m', labelcolor='m')
	ax3.set_ylabel('Zeff', color='m', fontsize=16)
	ax3.spines['right'].set_color('m')
	ax3.spines['left'].set_color('c')
	ax3.set_ylim([1, 7])
	
	fig.tight_layout()
	fig.savefig(filedir+'Measure_{:.1f}.png'.format(time_now),dpi=300)



# here, choose the data folder.
# 10 mbar
# filedir = '/Users/yao/sshfs_mesopsl_data/job.64935/'  # fe=0.06
# filedir = '/Users/yao/sshfs_mesopsl_data/job.65387/'  # fe=0.1
# filedir = '/Users/yao/sshfs_mesopsl_data/job.65096/'  # fe=0.2
# filedir = '/Users/yao/sshfs_mesopsl_data/job.64960/'  # fe=1.0

# 5 mbar
# filedir = '/Users/yao/sshfs_mesopsl_data/job.65064/'  # fe=0.06
filedir = '/Users/yao/sshfs_mesopsl_data/job.65388/'  # fe=0.1
# filedir = '/Users/yao/sshfs_mesopsl_data/job.65095/'  # fe=0.2
# filedir = '/Users/yao/sshfs_mesopsl_data/job.65063/'  # fe=1.0

# 1 mbar
# filedir = '/Users/yao/sshfs_mesopsl_data/job.65066/'  # fe=0.06
# filedir = '/Users/yao/sshfs_mesopsl_data/job.65389/'  # fe=0.1
# filedir = '/Users/yao/sshfs_mesopsl_data/job.65094/'  # fe=0.2
# filedir = '/Users/yao/sshfs_mesopsl_data/job.65067/'  # fe=1.0

for i in range(81):
	i = i + 0
	filename = 'lasslab_hdf5_chk_'+str(i).zfill(4)
	data_dict = get_yt_data_1d(filedir+filename)
	time_now = i*0.5
	plot_variables(filedir,time_now)

