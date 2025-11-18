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
plt.rc("axes", labelsize=14)
plt.rc("legend", fontsize=12)

import scienceplots as splt

plt.style.use(["nature", "notebook", "grid", "high-vis"])

# ---- Fine-tuned data (from the image) ----
# a0 values (x-axis) and corresponding Bx in GGauss (y-axis)
a0 = np.array([280, 300, 325, 350, 375, 400, 500], dtype=float)
Bx = np.array([0.28, 0.53, 0.78, 1.27, 1.60, 2.05, 3.30], dtype=float)

a1 = np.array([225, 300, 325, 350, 375, 400, 450], dtype=float)
a2 = np.array([300, 325, 350, 375, 400], dtype=float)
Bx2 = np.array([0.45, 0.75, 0.82, 0.85, 0.70], dtype=float)
yerr2 = np.array([0.14, 0.10, 0.16, 0.09, 0.08], dtype=float)
# Bx1 = np.array([0.46, 0.38, 3.14, 1.88, 2.21], dtype=float) # filter: 0.5, 15
# Bx1 = np.array([0.17, 0.37, 1.77, 1.60, 1.85], dtype=float) # filter: 0.3, 1.5
Bx1_low = np.array(
    [0.07, 0.08, 0.22, 0.43, 0.37, 0.2, 0.3], dtype=float
)  # updated FFT, k0x, k0y, wx, wy: 0.4, 0.4, 0.3, 0.4; lineout at y=0 with average over 5 um
Bx1_up = np.array(
    [0.3, 0.6, 0.8, 0.8, 0.8, 0.6, 0.45], dtype=float
)  # updated FFT, k0x, k0y, wx, wy: 0.4, 0.4, 0.3, 0.4; lineout at y=0

# Central value as the midpoint
Bx1 = 0.5 * (Bx1_low + Bx1_up)

# Asymmetric error bars
yerr_lower = Bx1 - Bx1_low
yerr_upper = Bx1_up - Bx1
yerr = np.vstack([yerr_lower, yerr_upper])

# ---- Fit a power law: Bx = C * a0^p ----
log_a0 = np.log(a0)
log_Bx = np.log(Bx)
p, logC = np.polyfit(log_a0, log_Bx, 1)
C = np.exp(logC)

# Smooth curve for the fitted relation
a0_fit = np.linspace(200, 520, 300)
Bx_fit = C * a0_fit**p

# fit a power law to the new data points as well (optional)
log_a1 = np.log(a1)
log_Bx1 = np.log(Bx1)
p1, logC1 = np.polyfit(log_a1, log_Bx1, 1)
C1 = np.exp(logC1)

# Smooth curve for the fitted relation of new data points (optional)
a1_fit = np.linspace(200, 520, 300)
Bx1_fit = C1 * a1_fit**p1

log_a2 = np.log(a2)
log_Bx2 = np.log(Bx2)
p2, logC2 = np.polyfit(log_a2, log_Bx2, 1)
C2 = np.exp(logC2)

# ---- Plot ----
plt.figure(figsize=(6, 5))
plt.scatter(
    a0, Bx, s=120, color="b", label="former 3D PIC data"
)  # default color, similar to the figure
plt.plot(
    a0_fit, Bx_fit, "--", lw=1.6, color="b", label="fitting former"
)  # dashed fit line
# plt.scatter(
#     a1, Bx1, s=100, color="red", marker="x", label="new 3D PIC data"
# )  # new data points
plt.errorbar(
    a1,
    Bx1,
    yerr=yerr,
    fmt="x",
    capsize=4,
    linewidth=1.5,
    markersize=6,
    color="red",
    label="new 3D PIC data -- 300 fs",
)
# plt.plot(a1_fit, Bx1_fit, ":", lw=1.6, color="r", label="fitting new 300 fs")

plt.errorbar(
    a2,
    Bx2,
    yerr=yerr2,
    fmt="o",
    capsize=6,
    linewidth=1.5,
    markersize=8,
    color="green",
    label="new 3D PIC data -- < 300 fs",
)
Bx2_fit = C2 * a1_fit**p2
# plt.plot(a1_fit, Bx2_fit, "-.", lw=1.6, color="g", label="fitting new < 300 fs")


# Axes style to match the original figure
plt.xlabel(r"$a_0$", fontsize=14)
plt.ylabel(r"$B_x$ [GGauss]", fontsize=14)
plt.xlim(200, 520)
plt.ylim(0.0, 6.0)
plt.legend()
plt.title(f"Fitted power law: Bx = {C:.3e} * a0^{p:.2f}")
plt.tight_layout()
# plt.show()

# Print the fitted law for reference
# print(f"Fitted power law: Bx = {C:.3e} * a0^{p:.2f}")
plt.savefig("/Users/yao/Desktop/scaling.png", dpi=300, bbox_inches="tight")
