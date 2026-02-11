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

# a1 = np.array([225, 300, 325, 350, 375, 400, 450], dtype=float)
# Bx1 = np.array([0.46, 0.38, 3.14, 1.88, 2.21], dtype=float) # filter: 0.5, 15
# Bx1 = np.array([0.17, 0.37, 1.77, 1.60, 1.85], dtype=float) # filter: 0.3, 1.5
# Bx1_low = np.array(
#     [0.07, 0.08, 0.22, 0.43, 0.37, 0.2, 0.3], dtype=float
# )  # updated FFT, k0x, k0y, wx, wy: 0.4, 0.4, 0.3, 0.4; lineout at y=0 with average over 5 um
# Bx1_up = np.array(
#     [0.3, 0.6, 0.8, 0.8, 0.8, 0.6, 0.45], dtype=float
# )  # updated FFT, k0x, k0y, wx, wy: 0.4, 0.4, 0.3, 0.4; lineout at y=0
# Central value as the midpoint
# Bx1 = 0.5 * (Bx1_low + Bx1_up)
#
# # Asymmetric error bars
# yerr_lower = Bx1 - Bx1_low
# yerr_upper = Bx1_up - Bx1
# yerr = np.vstack([yerr_lower, yerr_upper])

a2 = np.array([300, 325, 350, 375, 400, 450], dtype=float)
Bx2 = np.array([0.43, 0.60, 0.68, 0.90, 0.73, 0.92], dtype=float)
yerr2 = np.array([0.01, 0.14, 0.16, 0.01, 0.01, 0.10], dtype=float)

a3 = np.array([250, 300, 350, 400, 450, 500], dtype=float)
# Bx3 = np.array([0.31, 0.79, 1.18, 1.41, 2.52, 2.01], dtype=float)
# yerr3 = np.array([0.04, 0.02, 0.20, 0.18, 0.11, 0.28], dtype=float)
# at t=216 fs
# Bx3 = np.array([0.38, 0.67, 1.00, 1.75, 2.52, 2.22], dtype=float)
# yerr3 = np.array([0.36, 0.39, 0.48, 0.66, 0.11, 1.08], dtype=float)
# at t=233 fs, and optimized the ave and err
Bx3 = np.array([0.59, 0.75, 1.28, 1.83, 2.19, 2.15], dtype=float)
yerr3 = np.array([0.16, 0.11, 0.10, 0.35, 0.03, 0.14], dtype=float)
# at t=250 fs
# Bx3 = np.array([0.45, 0.73, 0.95, 1.41, 1.81, 1.54], dtype=float)
# yerr3 = np.array([0.05, 0.21, 0.33, 0.18, 0.16, 0.08], dtype=float)

# ne=120nc, res=80, ppc=4
a4    = np.array([  300,  400,  500,  600,  700], dtype=float)
Bx4   = np.array([ 0.55, 1.41, 2.16, 2.25, 2.92], dtype=float)
yerr4 = np.array([0.004, 0.16, 0.17, 0.26, 0.28], dtype=float)

# ne=120nc, res=80, ppc=1
a5    = np.array([  200,  700], dtype=float)
Bx5   = np.array([ 0.39, 3.20], dtype=float)
yerr5 = np.array([0.001, 0.06], dtype=float)

# ne=60nc, res=80, ppc=4
a6    = np.array([  500,  600], dtype=float)
Bx6   = np.array([ 1.93, 3.13], dtype=float)
yerr6 = np.array([ 0.12, 0.31], dtype=float)


# ---- Fit a power law: Bx = C * a0^p ----
log_a0 = np.log(a0)
log_Bx = np.log(Bx)
p, logC = np.polyfit(log_a0, log_Bx, 1)
C = np.exp(logC)

# Smooth curve for the fitted relation
a0_fit = np.linspace(180, 720, 500)
Bx_fit = C * a0_fit**p

# fit a power law to the new data points as well (optional)
# log_a1 = np.log(a1)
# log_Bx1 = np.log(Bx1)
# p1, logC1 = np.polyfit(log_a1, log_Bx1, 1)
# C1 = np.exp(logC1)

# Smooth curve for the fitted relation of new data points (optional)
# a1_fit = np.linspace(200, 520, 300)
# Bx1_fit = C1 * a1_fit**p1

log_a2 = np.log(a2)
log_Bx2 = np.log(Bx2)
p2, logC2 = np.polyfit(log_a2, log_Bx2, 1)
C2 = np.exp(logC2)

log_a3 = np.log(a3[:-1])
log_Bx3 = np.log(Bx3[:-1])
p3, logC3 = np.polyfit(log_a3, log_Bx3, 1)
C3 = np.exp(logC3)

Bx3_fit = C3 * a0_fit**p3

# ---- Plot ----
plt.figure(figsize=(6, 5))
plt.scatter(
    a0, Bx, s=120, color="b", label="former: with RR",
    facecolors="none",
)  # default color, similar to the figure
# plt.plot(
#     a0_fit, Bx_fit, "--", lw=1.6, color="b", label=f"fit Bx = {C:.2e} * a0^{p:.2f}"
# )  # dashed fit line
# plt.scatter(
#     a1, Bx1, s=100, color="red", marker="x", label="new 3D PIC data"
# )  # new data points
# plt.errorbar(
#     a1,
#     Bx1,
#     yerr=yerr,
#     fmt="x",
#     capsize=4,
#     linewidth=1.5,
#     markersize=6,
#     color="red",
#     label="new 3D PIC data -- 300 fs",
# )
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
    label="w/o RR",
    markerfacecolor="none",
    
)
# Bx2_fit = C2 * a1_fit**p2
# plt.plot(a1_fit, Bx2_fit, "-.", lw=1.6, color="g", label="fitting new < 300 fs")
#
plt.errorbar(
    a3,
    Bx3,
    yerr=yerr3,
    fmt="*",
    capsize=6,
    linewidth=1.5,
    markersize=12,
    color="red",
    label="w. RR=LL",
    markerfacecolor="none",
)
# plt.plot(
#     a0_fit, Bx3_fit, "-.", lw=1.6, color="r", label=f"fit Bx = {C3:.2e} * a0^{p3:.2f}"
# )

plt.errorbar(
    a4,
    Bx4,
    yerr=yerr4,
    fmt="s",
    capsize=6,
    linewidth=1.5,
    markersize=8,
    color="purple",
    label="new: w. RR=LL, ne=120",
    markerfacecolor="none",
)

# plt.errorbar(
#     a5,
#     Bx5,
#     yerr=yerr5,
#     fmt="s",
#     # linestyle='--',
#     capsize=6,
#     linewidth=1.5,
#     markersize=8,
#     color="purple",
#     label="new: w. RR=LL, ppc=1",
#     markerfacecolor="none",
# )

plt.errorbar(
    a6,
    Bx6,
    yerr=yerr6,
    fmt="s",
    # linestyle='--',
    capsize=6,
    linewidth=1.5,
    markersize=8,
    color="cyan",
    label="new: w. RR=LL, ne=60",
    markerfacecolor="none",
)

# Axes style to match the original figure
plt.xlabel(r"$a_0$", fontsize=16)
plt.ylabel(r"$B_x$ [GGauss]", fontsize=16)
plt.xlim(150, 750)
plt.ylim(0.0, 4.0)

# handles, labels = plt.gca().get_legend_handles_labels()
# # Reorder them manually in the sequence you want:
# # order = [3, 4, 0, 1, 2]  # example: new w/o RR → new w. RR → former → fit1 → fit2
# order = [0, 1, 4, 2, 3]
# plt.legend(
#     [handles[i] for i in order], [labels[i] for i in order], fontsize=14, frameon=True
# )
plt.legend(fontsize=12, frameon=True)
# plt.title(f"Fitted power law: Bx = {C:.3e} * a0^{p:.2f}")
# plt.title("at an early time=216 fs")
plt.tight_layout()
# plt.show()

# Print the fitted law for reference
# print(f"Fitted power law: Bx = {C:.3e} * a0^{p:.2f}")
plt.savefig("/Users/yao/Desktop/scaling.png", dpi=300, bbox_inches="tight")
