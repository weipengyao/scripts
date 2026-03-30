#!/usr/bin/env python3
"""
Generate Figure 2: η_rad vs a₀ for the radiative orbit note.

Three curves:
  1. Exact self-consistent solution (Eq. 26): η = 2√(1-η/2) · ξ'a₀³ f
  2. Approximate formula (Eq. 27):            η ≈ 2 ξ'a₀³ f
  3. Low-field asymptote (Eq. 28):            η ≈ (2/5) ξ'a₀³

Usage:
    python3 plot_eta_rad.py

Output:
    eta_rad.pdf
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from scipy.integrate import quad

# ── Physical parameters ──────────────────────────────────────────────
r_e   = 2.8179e-15          # m, classical electron radius
lam   = 1.0e-6              # m, laser wavelength
c     = 2.9979e8            # m/s
omega = 2 * np.pi * c / lam
xi    = 2 * r_e * omega / (3 * c)
a_cr  = xi ** (-1.0 / 3.0)

print(f"ξ    = {xi:.4e}")
print(f"a_cr = {a_cr:.1f}")


# ── Solve Eq. (14): γ'²(1 + ξ'² γ'⁶) = ā² ──────────────────────────
def solve_gamma(a, xi_val):
    """Return γ' for given local amplitude ā."""
    def res(g):
        return g**2 * (1.0 + xi_val**2 * g**6) - a**2
    if a < 1e-10:
        return 0.0
    return brentq(res, 1e-6, a)


# ── Skin-depth averaging function f(ξ', a₀) (Eq. 25) ───────────────
def compute_f(a0, xi_val):
    """f = (1/a₀⁵) ∫₀^{a₀} γ'^4(ξ', ā) dā"""
    def integrand(abar):
        g = solve_gamma(abar, xi_val)
        return g**4
    val, _ = quad(integrand, 0, a0, limit=200)
    return val / a0**5


# ── Exact η_rad from Eq. (26) ───────────────────────────────────────
def eta_exact(a0, xi_val):
    """Solve η = 2√(1 - η/2) · ξ'a₀³ f  for η."""
    f_val = compute_f(a0, xi_val)
    C = xi_val * a0**3 * f_val   # the combination ξ'a₀³f

    def res(eta):
        return eta - 2.0 * np.sqrt(max(1.0 - 0.5 * eta, 0.0)) * C
    # η lies in [0, 2) (since 1-η/2 ≥ 0)
    try:
        return brentq(res, 0, 1.999)
    except ValueError:
        return np.nan


# ── Approximate η_rad from Eq. (27) ─────────────────────────────────
def eta_approx(a0, xi_val):
    """η ≈ 2 ξ'a₀³ f"""
    return 2.0 * xi_val * a0**3 * compute_f(a0, xi_val)


# ── Compute curves ──────────────────────────────────────────────────
a0_arr = np.linspace(300, 900, 200)

print("Computing exact solution (this takes ~30s)...")
eta_ex  = np.array([eta_exact(a0, xi)  for a0 in a0_arr])
eta_app = np.array([eta_approx(a0, xi) for a0 in a0_arr])
eta_low = (2.0 / 5.0) * xi * a0_arr**3   # Eq. 28

# Saturation level
eta_sat = (np.sqrt(17) - 1) / 4.0
print(f"Saturation level: η_∞ = (√17 − 1)/4 ≈ {eta_sat:.4f}")

# Print a few values
for a0_check in [300, 440, 600, 900]:
    idx = np.argmin(np.abs(a0_arr - a0_check))
    print(f"  a₀={a0_check:4d}:  η_exact={eta_ex[idx]:.4f},"
          f"  η_approx={eta_app[idx]:.4f},"
          f"  η_low={eta_low[idx]:.4f}")


# ── Figure ───────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(6, 4.2))

ax.plot(a0_arr, eta_ex,  "C3-",  lw=2.2,
        label=r"exact (Eq.26)")
ax.plot(a0_arr, eta_app, "C0--", lw=1.8,
        label=r"approx (Eq.27)")
ax.plot(a0_arr, eta_low, "C2:",  lw=1.8,
        label=r"low-field: $\frac{2}{5}\xi^\prime a_0^3$  (Eq.28)")

# Saturation line at 0.78
ax.axhline(eta_sat, color="0.50", ls="-", lw=0.8, alpha=0.6)
ax.text(870, eta_sat + 0.02,
        rf"$\eta_{{\rm rad}}^{{(\infty)}} \approx {eta_sat:.2f}$",
        fontsize=10, color="0.40", ha="right", va="bottom")

# a_cr vertical line
ax.axvline(a_cr, color="0.60", ls=":", lw=1.0)
ax.text(a_cr + 8, 0.03,
        rf"$a_{{\rm cr}} \approx {a_cr:.0f}$",
        fontsize=10, color="0.50")

ax.set_xlabel(r"Normalised laser amplitude $a_0$", fontsize=12)
ax.set_ylabel(r"Radiative efficiency $\eta_{\rm rad}$", fontsize=12)
ax.set_xlim(300, 900)
ax.set_ylim(0, 1.05)
ax.legend(fontsize=10, loc="best")
ax.tick_params(labelsize=10)

fig.tight_layout()
fig.savefig("eta_rad.pdf", bbox_inches="tight")
fig.savefig("eta_rad.png", bbox_inches="tight", dpi=150)
plt.close(fig)
print("Wrote eta_rad.pdf / .png")
