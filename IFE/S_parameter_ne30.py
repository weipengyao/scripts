import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager

sys.path.append("/Users/yao/Smilei")
import happi

jetcmap = mpl.colormaps.get_cmap("jet").resampled(9)
jet_vals = jetcmap(np.arange(9))
jet_vals[0] = [1.0, 1.0, 1.0, 1.0]
jet_vals[8] = [0.0, 0.0, 0.0, 1.0]
newcmap = mpl.colors.LinearSegmentedColormap.from_list("mine", jet_vals)

font_dirs = ["/Users/yao/Documents/Calibri and Cambria Fonts/"]
for font_file in font_manager.findSystemFonts(fontpaths=font_dirs):
    font_manager.fontManager.addfont(font_file)

plt.rcParams["font.family"] = "Calibri"
plt.rc("text", usetex=False)
plt.rc("xtick", labelsize=14)
plt.rc("ytick", labelsize=14)
plt.rc("axes", labelsize=14)
plt.rc("legend", fontsize=12)

import scienceplots  # noqa: F401

plt.style.use(["nature", "notebook", "grid", "high-vis"])


DATA_ROOT = "/Users/yao/Documents/Data/IFE"
l0_SI = 0.8e-6
wr = 2.0 * np.pi * 3e8 / l0_SI

cases = [
    ("a0=200",         f"{DATA_ROOT}/ife_ne30_a200_lx2/"),
    ("a0=400",         f"{DATA_ROOT}/ife_ne30_a400_lx2/"),
    ("a0=600 (b)",     f"{DATA_ROOT}/ife_ne30_a600_lx2_b/"),
    ("a0=600 (r0+r1)", [f"{DATA_ROOT}/ife_ne30_a600_lx2_r0/",
                        f"{DATA_ROOT}/ife_ne30_a600_lx2_r1/"]),
    ("a0=800",         f"{DATA_ROOT}/ife_ne30_a800_lx2/"),
    ("a0=1000",        f"{DATA_ROOT}/ife_ne30_a1000_lx2/"),
]


def get_S_para_t(S, label):
    Ey_all = np.array(S.Probe(0, "Ey", units=["um", "fs"]).getData())
    ne_all = np.array(S.ParticleBinning(3, units=["um", "fs"]).getData())
    nt = min(Ey_all.shape[0], ne_all.shape[0])
    yc = Ey_all.shape[2] // 2
    zc = ne_all.shape[3] // 2
    print(f"[{label}] Ey {Ey_all.shape}, ne {ne_all.shape}, yc={yc}, zc={zc}, nt={nt}")

    S_para_t = np.zeros(nt)
    for i in range(nt):
        S_para_t[i] = np.sum(ne_all[i][:, yc, zc] * np.abs(Ey_all[i][:, yc]) ** 2)

    time = np.array(S.Probe(0, "Ey", units=["um", "fs"]).getTimes())[:nt]
    return time, S_para_t


fig, ax = plt.subplots()
for label, path in cases:
    S = happi.Open(path, reference_angular_frequency_SI=wr)
    t, S_t = get_S_para_t(S, label)
    ax.plot(t, S_t, label=label)

ax.set_xlabel("time (fs)")
ax.set_ylabel("S_para (NOT normalized)")
ax.set_yscale("log")
ax.legend()
ax.set_title(r"$S=\int_0^{L_x} dx\, n_e |E_y|^2$")

out = "/Users/yao/Desktop/S_para_t_ne30_scan.png"
fig.savefig(out, dpi=300, bbox_inches="tight")
print(f"saved: {out}")
