#!/usr/bin/env python3
"""
Schematic: co-moving frame field decomposition.
θ' = angle between E and centripetal (-r̂).
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
import matplotlib.patheffects as pe

R = 1.0
phi_e = 55.0           # electron at ~2 o'clock
theta_lag = 35.0       # phase lag θ'

phi_r = np.radians(phi_e)
theta_r = np.radians(theta_lag)

ex, ey = R * np.cos(phi_r), R * np.sin(phi_r)
r_hat = np.array([np.cos(phi_r), np.sin(phi_r)])
t_hat = np.array([-np.sin(phi_r), np.cos(phi_r)])
E_dir = np.cos(theta_r) * (-r_hat) + np.sin(theta_r) * t_hat

v_len  = 0.62
E_len  = 1.18
comp_c = E_len * np.cos(theta_r)
comp_t = E_len * np.sin(theta_r)

def ang_deg(v): return np.degrees(np.arctan2(v[1], v[0])) % 360
neg_r_deg = ang_deg(-r_hat)
E_deg     = ang_deg(E_dir)

fig, ax = plt.subplots(figsize=(5.2, 5.8))
ax.set_aspect("equal")
ax.set_xlim(-1.90, 1.70)
ax.set_ylim(-1.60, 2.35)
ax.axis("off")

bg  = [pe.withStroke(linewidth=3.5, foreground="white")]
bg2 = [pe.withStroke(linewidth=2.5, foreground="white")]

# ── Orbit circle ──
tc = np.linspace(0, 2 * np.pi, 300)
ax.plot(R * np.cos(tc), R * np.sin(tc),
        color="0.78", ls="--", lw=0.8, zorder=1)
ax.plot(0, 0, "+", color="0.68", ms=7, mew=0.9, zorder=1)

# Orbit direction arrow
ap = np.radians(phi_e - 60)
ax.annotate("",
    xy=(R * np.cos(ap + 0.14), R * np.sin(ap + 0.14)),
    xytext=(R * np.cos(ap), R * np.sin(ap)),
    arrowprops=dict(arrowstyle="->,head_width=0.04,head_length=0.04",
                    color="0.62", lw=0.8))

# Faint -r̂ reference line
ref_end = np.array([ex, ey]) + (comp_c + 0.22) * (-r_hat)
ax.plot([ex, ref_end[0]], [ey, ref_end[1]],
        color="0.82", ls="-", lw=0.7, zorder=1)
# tiny label
ax.text(ref_end[0] - 0.06, ref_end[1] - 0.10,
        r"$-\hat{r}$", fontsize=10, color="0.58",
        ha="center", va="top", path_effects=bg2, zorder=3)

# ── Electron ──
ax.plot(ex, ey, "ko", ms=9, zorder=10)
ax.text(ex + 0.14, ey + 0.08, r"$e^-$", fontsize=14,
        fontweight="bold", va="center", zorder=10, path_effects=bg)

# ── v (blue, tangent) ──
v_end = np.array([ex, ey]) + v_len * t_hat
ax.annotate("", xy=v_end, xytext=(ex, ey),
            arrowprops=dict(arrowstyle="->,head_width=0.10,head_length=0.12",
                            color="C0", lw=2.4, shrinkA=5, shrinkB=0),
            zorder=8)
ax.text(v_end[0] + 0.02, v_end[1] + 0.13,
        r"$\mathbf{v}$", fontsize=17, color="C0", fontweight="bold",
        ha="center", path_effects=bg, zorder=9)

# ── eE_L (red) ──
E_end = np.array([ex, ey]) + E_len * E_dir
ax.annotate("", xy=E_end, xytext=(ex, ey),
            arrowprops=dict(arrowstyle="->,head_width=0.10,head_length=0.12",
                            color="C3", lw=2.4, shrinkA=5, shrinkB=0),
            zorder=7)
# Place label on left side of arrow tip
E_label_offset = 0.08 * t_hat + 0.06 * (-r_hat)
ax.text(E_end[0] + E_label_offset[0] - 0.16, E_end[1] + E_label_offset[1] + 0.06,
        r"$e\mathbf{E}_L$", fontsize=17, color="C3", fontweight="bold",
        ha="right", path_effects=bg, zorder=9)

# ── Centripetal component (green dashed, along -r̂) ──
c_end = np.array([ex, ey]) + comp_c * (-r_hat)
ax.annotate("", xy=c_end, xytext=(ex, ey),
            arrowprops=dict(arrowstyle="->,head_width=0.07,head_length=0.09",
                            color="C2", lw=1.7, linestyle=(0, (5, 3)),
                            shrinkA=5, shrinkB=0),
            zorder=6)

# ── Tangential component (orange dashed, along t̂) ──
t_end = np.array([ex, ey]) + comp_t * t_hat
ax.annotate("", xy=t_end, xytext=(ex, ey),
            arrowprops=dict(arrowstyle="->,head_width=0.07,head_length=0.09",
                            color="C1", lw=1.7, linestyle=(0, (5, 3)),
                            shrinkA=5, shrinkB=0),
            zorder=6)

# Parallelogram completion
ax.plot([c_end[0], E_end[0]], [c_end[1], E_end[1]],
        color="0.65", ls=":", lw=0.7, zorder=2)
ax.plot([t_end[0], E_end[0]], [t_end[1], E_end[1]],
        color="0.65", ls=":", lw=0.7, zorder=2)

# ── θ' arc between E and -r̂ ──
arc_r = 0.38
# Arc CCW from E_deg to neg_r_deg
arc = Arc((ex, ey), 2 * arc_r, 2 * arc_r,
          angle=0, theta1=E_deg, theta2=neg_r_deg,
          color="0.20", lw=1.4, zorder=5)
ax.add_patch(arc)
mid_ang = np.radians(0.5 * (E_deg + neg_r_deg))
ax.text(ex + (arc_r + 0.16) * np.cos(mid_ang),
        ey + (arc_r + 0.16) * np.sin(mid_ang),
        r"$\theta'$", fontsize=15, ha="center", va="center",
        color="0.20", zorder=5, path_effects=bg)

# ── Labels on components ──
# cos θ' label: to the right of the centripetal arrow
c_mid = np.array([ex, ey]) + 0.45 * comp_c * (-r_hat)
c_perp = t_hat  # perpendicular direction for offset
ax.text(c_mid[0] - 0.22 * c_perp[0], c_mid[1] - 0.22 * c_perp[1],
        r"$\cos\theta'$", fontsize=14, color="C2",
        ha="center", va="center", path_effects=bg, zorder=8)

# sin θ' label: above the tangential arrow
t_mid = np.array([ex, ey]) + 0.55 * comp_t * t_hat
ax.text(t_mid[0] + 0.16 * r_hat[0], t_mid[1] + 0.16 * r_hat[1],
        r"$\sin\theta'$", fontsize=14, color="C1",
        ha="center", va="center", path_effects=bg, zorder=8)

# ── Equations (top left) ──
eq_x, eq_y = -1.78, 2.22
ax.text(eq_x, eq_y,
        r"centripetal:   $eE_L\cos\theta' = \gamma^\prime m\omega^\prime c$",
        fontsize=11, color="C2", path_effects=bg, va="center")
ax.text(eq_x, eq_y - 0.30,
        r"radiation:      $eE_L\sin\theta' = P^\prime_{\rm rad}/c$",
        fontsize=11, color="C1", path_effects=bg, va="center")

# ── Bottom label ──
ax.text(-0.10, -1.48, "circular orbit  (co-moving frame)",
        fontsize=11, ha="center", color="0.50", style="italic")

fig.tight_layout(pad=0.3)
fig.savefig("schematic.pdf", bbox_inches="tight", dpi=300)
fig.savefig("schematic.png", bbox_inches="tight", dpi=200)
plt.close(fig)
print("Done: schematic.pdf / .png")
