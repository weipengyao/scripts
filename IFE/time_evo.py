import matplotlib.pyplot as plt
import numpy as np
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
plt.rc("axes", labelsize=16)
plt.rc("legend", fontsize=12)

t = np.array([216, 233, 250, 266, 282, 300], dtype=float)

cases_ordered = [
    (
        250,
        "#21 a0=250 LL",
        np.array([0.38, 0.43, 0.43, 0.32, 0.35, 0.31]),
        np.array([0.36, 0.33, 0.05, 0.06, 0.04, 0.04]),
        "LL",
    ),
    (
        300,
        "#19 a0=300 LL",
        np.array([0.67, 0.64, 0.73, 0.81, 0.94, 0.79]),
        np.array([0.39, 0.22, 0.21, 0.13, 0.05, 0.02]),
        "LL",
    ),
    (
        350,
        "#16 a0=350 LL",
        np.array([1.00, 1.18, 0.95, 1.03, 1.04, 0.89]),
        np.array([0.48, 0.20, 0.33, 0.20, 0.18, 0.12]),
        "LL",
    ),
    (
        350,
        "#17 a0=350 cLL",
        np.array([1.07, 0.75, 0.71, 0.59, 0.56, 0.62]),
        np.array([0.01, 0.50, 0.30, 0.39, 0.44, 0.33]),
        "cLL",
    ),
    (
        400,
        "#18 a0=400 LL",
        np.array([1.75, 1.48, 1.41, 1.34, 1.20, 1.06]),
        np.array([0.66, 0.70, 0.18, 0.08, 0.01, 0.07]),
        "LL",
    ),
    (
        450,
        "#20 a0=450 LL",
        np.array([2.52, 2.22, 1.81, 1.67, 1.51, 1.34]),
        np.array([0.11, 0.06, 0.16, 0.04, 0.02, 0.04]),
        "LL",
    ),
    (
        450,
        "#23 a0=450 cLL",
        np.array([2.68, 1.60, 1.71, 1.40, 1.54, 1.23]),
        np.array([0.42, 0.60, 0.06, 0.02, 0.02, 0.01]),
        "cLL",
    ),
    (
        500,
        "#22 a0=500 LL",
        np.array([2.22, 2.01, 1.54, 1.40, 1.28, 1.24]),
        np.array([1.08, 0.28, 0.08, 0.06, 0.01, 0.06]),
        "LL",
    ),
]

# ---------------------------------------------
# Plot — uses error bars and publication style
# ---------------------------------------------
plt.figure(figsize=(10, 6))

for a0, label, vals, errs, model_type in cases_ordered:
    linestyle = "-" if model_type == "LL" else "--"
    marker = "o" if model_type == "LL" else "s"

    plt.errorbar(
        t,
        vals,
        yerr=errs,
        linestyle=linestyle,
        marker=marker,
        capsize=4,
        label=label,
        markersize=8,
    )

plt.xlabel("Time (fs)", fontsize=16)
plt.ylabel("Magnetic Field Strength (GG)", fontsize=16)
# plt.title("Time Evolution of B-field (ordered by $a_0$)")
plt.grid(True, which="both", linestyle=":", linewidth=0.6, alpha=0.5)
plt.legend(ncol=2, fontsize=14)
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)
plt.xlim(210, 310)
plt.ylim(0, 3.5)
plt.tight_layout()
# plt.show()
plt.savefig("/Users/yao/Desktop/time_evo_Bfield_ordered_a0.png", dpi=300)
