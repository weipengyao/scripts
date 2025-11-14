import sys

import matplotlib as mpl

sys.path.append("/Users/yao/Smilei")
import happi
import matplotlib.pyplot as plt
import numpy as np

jetcmap = plt.cm.get_cmap(
    "jet", 9
)  # generate a jet map with 10 values "rainbow", "jet", YlOrRd
jet_vals = jetcmap(np.arange(9))  # extract those values as an array
jet_vals[0] = [1.0, 1, 1.0, 1]  # change the first value
jet_vals[8] = [0.0, 0, 0.0, 1]  # change the first value
newcmap = mpl.colors.LinearSegmentedColormap.from_list("mine", jet_vals)

from matplotlib import font_manager

font_dirs = ["/Users/yao/Documents/Calibri and Cambria Fonts/"]
font_files = font_manager.findSystemFonts(fontpaths=font_dirs)

for font_file in font_files:
    font_manager.fontManager.addfont(font_file)

# set font
plt.rcParams["font.family"] = "Calibri"

plt.rc("text", usetex=False)
plt.rc("xtick", labelsize=14)
plt.rc("ytick", labelsize=14)
plt.rc("axes", labelsize=14)
plt.rc("legend", fontsize=12)


# Load simulation data
#
## prepare simulation data

wkdir = [
    #  '/Users/yao/Documents/Data/IFE/test0/laser_propagation_3d', # with only test particles
    #  '/Users/yao/Documents/Data/IFE/test0/laser_propagation_3d_target/', # with a thin foil target
    #  '/Users/yao/Desktop/test1/',
    #  '/Users/yao/Desktop/IFE_FlatTop_3D/',
    #  '/Users/yao/Desktop/ife_yao0/',
    #  '/Users/yao/Desktop/data/ife_yao2_cp2/',
    #  '/Users/yao/Desktop/data/ife_yao2_cp3/',
    #  '/Users/yao/Desktop/ife_yao3/',
    #  '/Users/yao/Documents/Data/IFE/ife_yao1_Ly25.6/',
    "/Users/yao/Documents/Data/IFE/ife_yao4_benchmark_probe_reso/",  # a0=225, tp=5t0, ne=60nc
    "/Users/yao/Documents/Data/IFE/ife_yao6/",  # a0=300, tp=5t0, ne=60nc
    "/Users/yao/Documents/Data/IFE/ife_yao9/",  # a0=325, tp=5t0, ne=60nc
    "/Users/yao/Documents/Data/IFE/ife_yao7/",  # a0=350, tp=5t0, ne=60nc
    "/Users/yao/Documents/Data/IFE/ife_yao10/",  # a0=375, tp=5t0, ne=60nc
    "/Users/yao/Documents/Data/IFE/ife_yao8/",  # a0=400, tp=5t0, ne=60nc
    "/Users/yao/Documents/Data/IFE/ife_yao5/",  # a0=450, tp=5t0, ne=60nc
    "/Users/yao/Documents/Data/IFE/ife_yao11/",  # a0=225, tp=50t0, ne=60nc
    "/Users/yao/Documents/Data/IFE/ife_yao12/",  # a0=225, tp=5t0, ne=15nc
]

S0 = happi.Open(
    wkdir[0], reference_angular_frequency_SI=2.0 * np.pi * 3e8 / (1.0 * 1e-6)
)
S1 = happi.Open(
    wkdir[1], reference_angular_frequency_SI=2.0 * np.pi * 3e8 / (1.0 * 1e-6)
)
S2 = happi.Open(
    wkdir[2], reference_angular_frequency_SI=2.0 * np.pi * 3e8 / (1.0 * 1e-6)
)
S3 = happi.Open(
    wkdir[3], reference_angular_frequency_SI=2.0 * np.pi * 3e8 / (1.0 * 1e-6)
)
S4 = happi.Open(
    wkdir[4], reference_angular_frequency_SI=2.0 * np.pi * 3e8 / (1.0 * 1e-6)
)
S5 = happi.Open(
    wkdir[5], reference_angular_frequency_SI=2.0 * np.pi * 3e8 / (1.0 * 1e-6)
)
S6 = happi.Open(
    wkdir[6], reference_angular_frequency_SI=2.0 * np.pi * 3e8 / (1.0 * 1e-6)
)
S7 = happi.Open(
    wkdir[7], reference_angular_frequency_SI=2.0 * np.pi * 3e8 / (1.0 * 1e-6)
)
S8 = happi.Open(
    wkdir[8], reference_angular_frequency_SI=2.0 * np.pi * 3e8 / (1.0 * 1e-6)
)


# Get simulation box size
Lx = S0.namelist.Lx / 2 / np.pi  # in um
print("Lx = ", Lx)
Ly = S0.namelist.Ly / 2 / np.pi  # in um
print("Ly = ", Ly)
Lz = S0.namelist.Lz / 2 / np.pi  # in um
print("Lz = ", Lz)


def fft_field_analysis(
    S,
    diag_num,
    field_name,
    slice_direction,
    slice_position,
    shift,
    timestep,
    k0x,
    k0y,
    wx,
    wy,
):
    """
    Perform FFT analysis on a given field from the simulation data.

    Parameters:
    S : happi.Simulation
        The simulation object.
    diag_num : int
        The diagnostic number of the field.
    field_name : str
        The name of the field to analyze.
    slice_direction : str
        e.g., 'x', 'y', 'z', etc.
    slice_position : float
        The position to slice the field (e.g., Lz/2).
    shift : bool
        Whether to shift the sample region to avoid boundary condition issues.
    timestep : int
        The timestep at which to perform the analysis.
    k_filter_x, k_filter_y : float
        The mask to filter out the signals. Laser is 1.

    Returns:
    kx : np.ndarray
        The wave numbers in the x-direction.
    ky : np.ndarray
        The wave numbers in the y-direction.
    spectrum : np.ndarray
        The FFT spectrum of the field.
    """

    twopi = 2.0 * np.pi
    var0_time = S.Probe(diag_num, field_name, units=["fs"]).getTimesteps()
    print("timesteps of var : ", var0_time)
    print("num. of timesteps: ", var0_time.size)

    var0 = S.Probe(diag_num, field_name, units=["fs"]).getData(
        timestep=var0_time[timestep]
    )[0]
    print("shape of var: ", var0.shape)

    var0_time_real = S.Probe(diag_num, field_name, units=["fs"]).getTimes()
    print("time of var : {:.2f} fs".format(var0_time_real[timestep]))

    field = S.Probe(
        diag_num,
        field_name,
        subset={slice_direction: slice_position},
        units=["um", "fs", "1e9 Gauss"],
        vsym=True,
    )

    data = np.array(field.getData()[timestep])

    if shift == True:
        nx = data.shape[0]
        ny = data.shape[1]
        sf = 0  # shift of cells to get rid of boundary condition
        print("note that the sample region is shifted to get rid of the BC issue")
        print("No. of cells being shifted, ", sf)
        #         print(nx, ny)
        data = data[sf : (nx - sf), sf : (ny - sf)]

    field_x = np.array(field.getAxis("axis1"))
    field_y = np.array(field.getAxis("axis2"))

    Nx = data.shape[0]
    Ny = data.shape[1]

    dx = (field_x[1] - field_x[0])[0]
    dy = (field_y[1] - field_y[0])[1]

    kx = np.fft.fftfreq(Nx, d=dx) * 2 * np.pi  # spatial freq extents
    ky = np.fft.fftfreq(Ny, d=dy) * 2 * np.pi
    KX, KY = np.meshgrid(kx, ky, indexing="ij")

    # This is another way to construct kx, ky grids -- but somehow doesn't work as expected
    # kx=np.linspace(-Nx/2,Nx/2-1,Nx)/((Nx)*dx)  ## -1? min&max??
    # ky=np.linspace(-Ny/2,Ny/2-1,Ny)/((Ny)*dy)
    # KX, KY = np.meshgrid(twopi*kx, twopi*ky, indexing='ij')

    # construct mask
    # super gaussian mask
    pwr = 12  # super gaussian power
    mask = np.exp(-(np.sqrt(((KX - k0x) / wx) ** 2 + ((KY - k0y) / wy) ** 2) ** pwr))
    # a quick check of the mask
    plt.imshow(mask.T, extent=(kx.min(), kx.max(), ky.min(), ky.max()), origin="lower")

    # another notch filter mask, which is not used for now
    def notch_mask(KX, KY, k0x, k0y, wx, wy, p=12):
        # four symmetric notches around ±k0x and ±k0y
        notch_1 = np.exp(
            -((((KX - k0x) / wx) ** 2 + ((KY - k0y) / wy) ** 2) ** (p / 2))
        )
        notch_2 = np.exp(
            -((((KX + k0x) / wx) ** 2 + ((KY - k0y) / wy) ** 2) ** (p / 2))
        )
        notch_3 = np.exp(
            -((((KX - k0x) / wx) ** 2 + ((KY + k0y) / wy) ** 2) ** (p / 2))
        )
        notch_4 = np.exp(
            -((((KX + k0x) / wx) ** 2 + ((KY + k0y) / wy) ** 2) ** (p / 2))
        )
        return (1 - notch_1) * (1 - notch_2) * (1 - notch_3) * (1 - notch_4)

    # mask = notch_mask(KX, KY, k0x, k0y, wx, wy, p=12)
    # plt.imshow(mask.T, extent=(kx.min(), kx.max(), ky.min(), ky.max()), origin='lower')
    # plt.xlim([0,1])
    # plt.ylim([0,1])

    data_fft = np.fft.fft2(data, axes=(0, 1))  # shape (Nx, Ny)
    Diag_data = np.real(np.fft.ifft2(data_fft * mask[:, :]))

    # data_fft  = np.fft.fft2(data, axes=(0,1))/Nx/Ny*2   # shape (Nx, Ny)
    data_fft_shift = np.fft.fftshift(data_fft)
    data_fft_shift_masked = data_fft_shift * mask[:, :]
    Diag_data_shift_back = np.fft.ifftshift(data_fft_shift * mask[:, :])
    # Diag_data = np.real(np.fft.ifft2(Diag_data_shift_back))
    # Diag_data = np.real(np.fft.ifft2(data_fft_shift* mask[:,:]))

    return kx, ky, data, data_fft, data_fft_shift, data_fft_shift_masked, Diag_data


k0x = 0.0
k0y = 0.0
wx = 0.5
wy = 2.0

# def fft_field_analysis(S, diag_num, field_name, slice_direction, slice_position, shift, timestep, k0x, k0y, wx, wy):

# kx0, ky0, Bx0, Bx0_fft, Bx0_fft_shift, Bx0_masked = fft_field_analysis(S0, 0, "Bx", 'z', Lz/2., -1, wx, wy)
# kx1, ky1, Bx1, Bx1_fft, Bx1_masked = fft_field_analysis(S1, 0, "Bx", 'z', Lz/2., -1, wx, wy)
# kx2, ky2, Bx2, Bx2_fft, Bx2_masked = fft_field_analysis(S2, 0, "Bx", 'z', Lz/2., -1, wx, wy)
# kx3, ky3, Bx3, Bx3_fft, Bx3_masked = fft_field_analysis(S3, 0, "Bx", 'z', Lz/2., -1, wx, wy)
# kx4, ky4, Bx4, Bx4_fft, Bx4_masked = fft_field_analysis(S4, 0, "Bx", 'z', Lz/2., -1, wx, wy)
# kx5, ky5, Bx5, Bx5_fft, Bx5_masked = fft_field_analysis(S5, 0, "Bx", 'z', Lz/2., -1, wx, wy)
# kx6, ky6, Bx6, Bx6_fft, Bx6_masked = fft_field_analysis(S6, 0, "Bx", 'z', Lz/2., -1, wx, wy)
# kx7, ky7, Bx7, Bx7_fft, Bx7_masked = fft_field_analysis(S7, 0, "Bx", 'z', Lz/2., -1, wx, wy)
# kx8, ky8, Bx8, Bx8_fft, Bx8_masked = fft_field_analysis(S8, 0, "Bx", 'z', Lz/2., -1, wx, wy)
#
kx3, ky3, Bx3, Bx3_fft, Bx3_fft_shifted, Bx3_fft_shifted_masked, Bx3_masked = (
    fft_field_analysis(S3, 0, "Bx", "z", Lz / 2.0, False, -1, k0x, k0y, wx, wy)
)

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

fig, ax = plt.subplots(figsize=(6, 5))

im = ax.imshow(
    # Bx3_masked.T,
    Bx3.T,
    origin="lower",
    extent=[0, Lx, 0, Ly],
    vmin=-20,
    vmax=20,
    aspect="equal",
    cmap="bwr",
)

# Make a colorbar axis that exactly matches the height of `ax`
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)  # width & gap
cb = fig.colorbar(im, cax=cax, label="Bx (GG)")

ax.set_xlabel("x (um)")
ax.set_ylabel("y (um)")
fig.tight_layout()
# plt.show()

# plt.tight_layout()
plt.savefig("/Users/yao/Desktop/Bx_a0_350.png", dpi=300, bbox_inches="tight")
#
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

fig, ax = plt.subplots(figsize=(6, 5))

im = ax.imshow(
    Bx3_masked.T,
    # Bx3.T,
    origin="lower",
    extent=[0, Lx, 0, Ly],
    vmin=-3,
    vmax=3,
    aspect="equal",
    cmap="bwr",
)

# Make a colorbar axis that exactly matches the height of `ax`
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)  # width & gap
cb = fig.colorbar(im, cax=cax, label="Bx (GG)")

ax.set_xlabel("x (um)")
ax.set_ylabel("y (um)")
fig.tight_layout()
plt.show()

# plt.tight_layout()
plt.savefig("/Users/yao/Desktop/Bx_a0_350_fft.png", dpi=300, bbox_inches="tight")


plt.figure()
plt.imshow(
    # np.log10(np.abs(Bx3_fft_shifted.T)),
    np.abs(Bx3_fft_shifted_masked.T),
    extent=(kx3.min(), kx3.max(), ky3.min(), ky3.max()),
    origin="lower",
    #    vmin=-5,
    #    vmax=-1,
    aspect="auto",
)
plt.colorbar(label="|Bx FFT|")
plt.xlim([-1, 1])
plt.ylim([-1, 1])
plt.xlabel("kx (1/um)")
plt.ylabel("ky (1/um)")
# plt.title('FFT of Bx field at z=Lz/2, time=133.10 fs')
plt.tight_layout()
plt.savefig("/Users/yao/Desktop/Bx_a0_350_mask.png", dpi=300, bbox_inches="tight")
