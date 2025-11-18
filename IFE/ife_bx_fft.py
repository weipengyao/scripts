from turtle import pos

import matplotlib

matplotlib.use("MacOSX")  # for mac users

# print(matplotlib.get_backend())

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

import scienceplots as splt

plt.style.use(["nature", "notebook", "grid", "high-vis"])

wkdir = [
    ## different benchmark cases
    #  '/Users/yao/Documents/Data/IFE/test0/laser_propagation_3d', # with only test particles
    #  '/Users/yao/Documents/Data/IFE/test0/laser_propagation_3d_target/', # with a thin foil target
    #  '/Users/yao/Desktop/test1/',
    #  '/Users/yao/Desktop/IFE_FlatTop_3D/',
    #  '/Users/yao/Desktop/ife_yao0/',
    #  '/Users/yao/Desktop/data/ife_yao2_cp2/',
    #  '/Users/yao/Desktop/data/ife_yao2_cp3/',
    #  '/Users/yao/Desktop/ife_yao3/',
    #  '/Users/yao/Documents/Data/IFE/ife_yao1_Ly25.6/',
    ## start to scan a0 from here
    # "/Users/yao/Documents/Data/IFE/ife_yao4_benchmark_probe_reso/",  # a0=225, tp=5t0, ne=60nc
    # "/Users/yao/Documents/Data/IFE/ife_yao6/",  # a0=300, tp=5t0, ne=60nc
    # "/Users/yao/Documents/Data/IFE/ife_yao9/",  # a0=325, tp=5t0, ne=60nc
    # "/Users/yao/Documents/Data/IFE/ife_yao7/",  # a0=350, tp=5t0, ne=60nc
    # "/Users/yao/Documents/Data/IFE/ife_yao10/",  # a0=375, tp=5t0, ne=60nc
    # "/Users/yao/Documents/Data/IFE/ife_yao8/",  # a0=400, tp=5t0, ne=60nc
    # "/Users/yao/Documents/Data/IFE/ife_yao5/",  # a0=450, tp=5t0, ne=60nc
    # "/Users/yao/Documents/Data/IFE/ife_yao11/",  # a0=225, tp=50t0, ne=60nc
    # "/Users/yao/Documents/Data/IFE/ife_yao12/",  # a0=225, tp=5t0, ne=15nc
    ## start to check the density effect
    "/Users/yao/Documents/Data/IFE/ife_yao13/",  # a0=350, tp=5t0, ne=240nc
]

S0 = happi.Open(
    wkdir[0], reference_angular_frequency_SI=2.0 * np.pi * 3e8 / (1.0 * 1e-6)
)
# S1 = happi.Open(
#     wkdir[1], reference_angular_frequency_SI=2.0 * np.pi * 3e8 / (1.0 * 1e-6)
# )
# S2 = happi.Open(
#     wkdir[2], reference_angular_frequency_SI=2.0 * np.pi * 3e8 / (1.0 * 1e-6)
# )
# S3 = happi.Open(
#     wkdir[3], reference_angular_frequency_SI=2.0 * np.pi * 3e8 / (1.0 * 1e-6)
# )
# S4 = happi.Open(
#     wkdir[4], reference_angular_frequency_SI=2.0 * np.pi * 3e8 / (1.0 * 1e-6)
# )
# S5 = happi.Open(
#     wkdir[5], reference_angular_frequency_SI=2.0 * np.pi * 3e8 / (1.0 * 1e-6)
# )
# S6 = happi.Open(
#     wkdir[6], reference_angular_frequency_SI=2.0 * np.pi * 3e8 / (1.0 * 1e-6)
# )
# S7 = happi.Open(
#     wkdir[7], reference_angular_frequency_SI=2.0 * np.pi * 3e8 / (1.0 * 1e-6)
# )
# S8 = happi.Open(
#     wkdir[8], reference_angular_frequency_SI=2.0 * np.pi * 3e8 / (1.0 * 1e-6)
# )

# Get simulation box size
Lx = S0.namelist.Lx / 2 / np.pi  # in um
print("Lx = ", Lx)
Ly = S0.namelist.Ly / 2 / np.pi  # in um
print("Ly = ", Ly)
Lz = S0.namelist.Lz / 2 / np.pi  # in um
print("Lz = ", Lz)


def get_fft(S, time, vmin0, vmax0):
    Bx = S.Probe(0, "Bx", units=["um", "fs", "1e5 T"])
    data = np.array(Bx.getData()[time])
    field_x = np.array(Bx.getAxis("axis1"))
    field_y = np.array(Bx.getAxis("axis2"))
    print("time = ", Bx.getTimes()[time], " fs")

    plt.figure(figsize=(6, 5))
    plt.imshow(
        data.T,
        extent=(field_x.min(), field_x.max(), field_y.min(), field_y.max()),
        origin="lower",
        cmap="bwr",
        aspect="auto",
        vmin=vmin0,
        vmax=vmax0,
    )
    plt.colorbar()
    plt.xlabel("x (um)")
    plt.ylabel("y (um)")
    plt.title("Bx field in real space (before FFT)")
    plt.tight_layout()
    plt.savefig("data_before_FFT.png", dpi=300)
    # plt.show()

    Nx = data.shape[0]
    Ny = data.shape[1]
    dx = (field_x[1] - field_x[0])[0]
    dy = (field_y[1] - field_y[0])[1]
    kx = np.fft.fftfreq(Nx, d=dx)  # / (2 * np.pi)  # spatial freq extents
    ky = np.fft.fftfreq(Ny, d=dy)  # / (2 * np.pi)
    KX, KY = np.meshgrid(kx, ky, indexing="ij")

    data_fft = np.fft.fft2(data, axes=(0, 1))  # shape (Nx, Ny)

    # plt.figure(figsize=(6, 5))
    # plt.imshow(
    #     np.log10(np.abs(data_fft).T),
    #     extent=(kx.min(), kx.max(), ky.min(), ky.max()),
    #     origin="lower",
    # )
    # plt.colorbar()
    # plt.xlabel("kx (1/um)")
    # plt.ylabel("ky (1/um)")
    # plt.title("Fourier Transform of Bx in log10 scale")
    # plt.tight_layout()
    # plt.show()

    data_fft_shift = np.fft.fftshift(data_fft)
    kx_shift = np.fft.fftshift(kx)
    ky_shift = np.fft.fftshift(ky)

    plt.figure(figsize=(6, 5))
    plt.imshow(
        np.log10(np.abs(data_fft_shift).T),
        extent=(kx_shift.min(), kx_shift.max(), ky_shift.min(), ky_shift.max()),
        origin="lower",
    )
    plt.colorbar()
    plt.xlabel("kx (1/um)")
    plt.ylabel("ky (1/um)")
    plt.title("Shifted (& Zoom) Fourier Transform of Bx in log10 scale")
    plt.xlim(-2, 2)
    plt.ylim(-2, 2)
    plt.tight_layout()
    plt.savefig("FFT_shift_zoom.png", dpi=300)
    # plt.show(block=True)

    print(
        "With these plots, you can identify the dominant spatial frequencies in the Bx field."
    )
    print(
        "You can now design your mask to filter out unwanted frequencies accordingly."
    )

    return data, data_fft, kx, ky, KX, KY, field_x, field_y


def apply_mask_and_ifft(data_fft, kx, ky, KX, KY, k0x, k0y, wx, wy):
    pwr = 12  # super gaussian power

    # # adjust the mask parameters here according to the FFT plot
    # k0x = 0.4
    # k0y = 0.4
    # wx = 0.3
    # wy = 0.4

    mask = 1 - np.exp(
        -(
            np.sqrt(((np.abs(KX) - k0x) / wx) ** 2 + ((np.abs(KY) - k0y) / wy) ** 2)
            ** pwr
        )
    )
    # a quick check of the mask
    # plt.figure(figsize=(6, 5))
    # plt.imshow(mask.T, extent=(kx.min(), kx.max(), ky.min(), ky.max()), origin="lower")
    # plt.colorbar()
    # plt.xlabel("kx (1/um)")
    # plt.ylabel("ky (1/um)")
    # plt.title("Fourier Mask")
    # plt.tight_layout()
    # plt.show()

    mask_shift = np.fft.fftshift(mask)

    # plt.figure(figsize=(6, 5))
    # plt.imshow(
    #     mask_shift.T, extent=(kx.min(), kx.max(), ky.min(), ky.max()), origin="lower"
    # )
    # plt.colorbar()
    # plt.xlabel("kx (1/um)")
    # plt.ylabel("ky (1/um)")
    # plt.title("Shifted & Zoomed Fourier Mask")
    # plt.xlim(-2, 2)
    # plt.ylim(-2, 2)
    # plt.tight_layout()
    # plt.show()

    data_fft_masked = data_fft * mask
    data_fft_masked_shift = np.fft.fftshift(data_fft_masked)

    plt.figure(figsize=(6, 5))
    plt.imshow(
        np.log10(np.abs(data_fft_masked_shift).T),
        extent=(kx.min(), kx.max(), ky.min(), ky.max()),
        origin="lower",
    )
    plt.colorbar()
    plt.xlabel("kx (1/um)")
    plt.ylabel("ky (1/um)")
    plt.title("Masked & Shifted Fourier Transform of Bx in log10 scale")
    plt.xlim(-2, 2)
    plt.ylim(-2, 2)
    plt.tight_layout()
    plt.savefig("FFT_masked_shift_zoom.png", dpi=300)
    # plt.show(block=True)

    return data_fft_masked


def compute_ifft(data_fft_masked, field_x, field_y, vmin0, vmax0):
    data_ifft = np.real(np.fft.ifft2(data_fft_masked))

    plt.figure(figsize=(6, 5))
    plt.imshow(
        data_ifft.T,
        extent=(field_x.min(), field_x.max(), field_y.min(), field_y.max()),
        origin="lower",
        cmap="bwr",
        aspect="auto",
        vmin=vmin0,
        vmax=vmax0,
    )
    plt.colorbar()
    plt.xlabel("x (um)")
    plt.ylabel("y (um)")
    plt.title("Filtered Bx field in real space (after IFFT)")
    plt.tight_layout()
    plt.savefig("data_after_FFT.png", dpi=300)
    # plt.show(block=True)

    return data_ifft


def plot_comparison(Bx_original, Bx_filtered):
    Nx = Bx_original.shape[0]
    Ny = Bx_original.shape[1]
    xx = np.linspace(0, Lx, Bx_original.shape[0])
    yy = np.linspace(0, Ly, Bx_original.shape[1])

    plt.figure(figsize=(6, 5))
    plt.plot(xx, Bx_original[:, Ny // 2], label="Original Bx")
    plt.plot(xx, Bx_filtered[:, Ny // 2], label="Filtered Bx")
    plt.xlabel("x (um)")
    plt.ylabel("Bx (1e5 T)")
    plt.title("Lineout of Bx at y=0")
    plt.xlim(0, Lx)
    plt.legend()
    plt.tight_layout()
    plt.savefig("data_lineout.png", dpi=300)
    # plt.show(block=True)
    #
    pos = 450
    # print(xx[pos])
    #
    print("max Bx_ori = ", np.max(np.abs(Bx_original[pos:, Ny // 2])), " GG")
    print("max Bx_filt = ", np.max(np.abs(Bx_filtered[pos:, Ny // 2])), " GG")
    print(
        "error bar = ",
        np.max(np.abs(Bx_original[pos:, Ny // 2]))
        - np.max(np.abs(Bx_filtered[pos:, Ny // 2])),
        " GG",
    )

    cell_y = 48  # 1 cell = 0.04 um along y
    center = int(Bx_original.shape[1] / 2)

    plt.figure(figsize=(6, 5))
    plt.plot(
        xx,
        np.average(Bx_original[:, (center - cell_y) : (center + cell_y)], axis=1),
        label="Original Bx",
        color="blue",
        linestyle="-",
    )
    plt.plot(
        xx,
        np.average(Bx_filtered[:, (center - cell_y) : (center + cell_y)], axis=1),
        label="Filtered Bx",
        color="red",
        linestyle="--",
    )
    plt.xlabel("x (um)")
    plt.ylabel("Bx (1e5 T)")
    plt.title(
        "Averaged Bx field at y = {:.1f} um with a width of {:.1f} um".format(
            yy[center], 2 * cell_y * (Ly / Bx_original.shape[1])
        )
    )
    plt.legend()
    plt.xlim(0, Lx)
    plt.tight_layout()
    plt.savefig("data_lineout_averaged.png", dpi=300)
    # plt.show(block=True)


def post_process_Bx_field(S, k0x, k0y, wx, wy, timestep, vmin0, vmax0):
    Bx0, Bx0_fft, kx0, ky0, KX0, KY0, field_x0, field_y0 = get_fft(
        S, timestep, vmin0, vmax0
    )
    Bx0_fft_masked = apply_mask_and_ifft(Bx0_fft, kx0, ky0, KX0, KY0, k0x, k0y, wx, wy)
    Bx0_ifft = compute_ifft(Bx0_fft_masked, field_x0, field_y0, vmin0, vmax0)
    plot_comparison(Bx0, Bx0_ifft)


post_process_Bx_field(
    S0, k0x=0.4, k0y=0.4, wx=0.3, wy=0.4, timestep=14, vmin0=-2, vmax0=2
)
# post_process_Bx_field(
# S2, k0x=0.4, k0y=0.4, wx=0.3, wy=0.4, timestep=-1, vmin0=-2, vmax0=2
# )
# post_process_Bx_field(S3, k0x=0.4, k0y=0.4, wx=0.3, wy=0.4, timestep=-1, vmin0=-2, vmax0=2)
