"""
Estimate electron temperature evolution heated by nanosecond lasers

Assumptions:
- Low temperature and low density (Te < 300 eV, ne/nc ~ 1e-2)
- Neglect electron thermal diffusion
- Heating only by inverse Bremsstrahlung (IB) absorption

Governing equation (CGS units mixed with eV for Te):
    (3/2) n_e dT_e/dt = ν_B * I(t) / q_e
where:
    n_e  : electron number density [cm^-3]
    T_e  : electron temperature [eV]
    ν_B  : IB absorption coefficient [cm^-1]
    I(t) : laser intensity [W/cm^2]
    q_e  : conversion between J and eV (handled explicitly below)
"""
from __future__ import annotations
import numpy as np
from dataclasses import dataclass

# ---------------------------
# Physical constants (CGS)
# ---------------------------
C_LIGHT = 3.0e10        # cm/s
QE_CGS  = 4.8e-10       # statcoulomb (esu)
HBAR    = 1.0e-27       # erg*s
ME_G    = 9.1e-28       # g
J_PER_EV = 1.602e-19    # J per eV (for explicit J↔eV conversion)

# ---------------------------
# Helper functions
# ---------------------------

def laser_intensity(E_J: float, tau_s: float, diameter_cm: float) -> tuple[float, float]:
    """Return (I, S) where I is average intensity [W/cm^2] and S is spot area [cm^2]."""
    area = np.pi * (diameter_cm / 2.0) ** 2
    intensity = E_J / tau_s / area
    return intensity, area


def angular_frequency(lambda_cm: float) -> float:
    """Laser angular frequency [rad/s] for wavelength in cm."""
    return 2.0 * np.pi * C_LIGHT / lambda_cm


def plasma_frequency(ne_cm3: float) -> float:
    """Electron plasma angular frequency [rad/s] (CGS shortcut)."""
    return 5.64e4 * np.sqrt(ne_cm3)


def tgaussian(time: np.ndarray, fwhm: float, order: float, center: float) -> np.ndarray:
    """Generalized time-Gaussian profile normalized to 1 at its peak."""
    sigma = (0.5 * fwhm) ** order / np.log(2.0)
    return np.exp(-((time - center) ** order) / sigma)


def coulomb_logarithm(temp_eV: float, ne_cm3: float, Z0: float, omg_ts: float) -> float:
    """Coulomb logarithm with simple bounding choices.

    Notes:
        - Uses standard CGS shortcuts for Debye length and thermal speed.
        - Keeps the original code's bounding logic, including the comment about
          possible SI→CGS confusion for b_min; retained for consistency.
    """
    # Debye length [cm]
    Ld = 7.43e2 * np.sqrt(temp_eV) * ne_cm3 ** (-0.5)
    # Ion spacing [cm]
    ni = ne_cm3 / Z0
    Ri = ni ** (-1.0 / 3.0)
    # Electron thermal speed [cm/s]
    vte = 4.19e7 * np.sqrt(temp_eV)

    bmax = np.minimum(np.maximum(Ld, Ri), vte / omg_ts)
    bmin = np.minimum(
        np.maximum(Z0 * QE_CGS ** 2 / (temp_eV * 1.6e-12), HBAR / np.sqrt(ME_G * temp_eV * 1.6e-12)),
        Ri,
    )
    return np.log(1.0 + 0.7 * bmax / bmin)


def nu_ib(temp_eV: float, ne_cm3: float, Z0: float, omg_ts: float, omg_pe: float) -> float:
    """Inverse Bremsstrahlung absorption coefficient ν_B [cm^-1].

    Implements the temperature dependence:
        ν_B ∝ Z0 * n_e^2 * Λ(T) * T^{-3/2} * ω^{-2} * (1 - ω_pe^2/ω^2)^{-1/2}

    The Coulomb logarithm Λ(T) is evaluated internally using `coulomb_logarithm`.
    """
    Lambda = coulomb_logarithm(temp_eV, ne_cm3, Z0, omg_ts)
    return (
        3.1e-7
        * Z0
        * ne_cm3 ** 2
        * Lambda
        * temp_eV ** (-1.5)
        * omg_ts ** -2.0
        * (1.0 - (omg_pe ** 2) / (omg_ts ** 2)) ** (-0.5)
    )


@dataclass
class SimulationParams:
    # Laser
    E_J: float = 110.0      # J
    tau_s: float = 5.0e-9   # s (pulse duration, also used as FWHM)
    diameter_cm: float = 0.1  # cm
    lambda_cm: float = 1.0e-4  # 1 µm

    # Plasma & material
    Z0: float = 1.0
    ne_cm3: float = 2.0e17
    Te0_eV: float = 0.1

    # Time grid
    t_end_s: float = 10e-9
    n_steps: int = 1000

    # Temporal laser profile
    order: float = 4.0


@dataclass
class SimulationOutputs:
    time_s: np.ndarray
    Te_eV: np.ndarray
    I_Wcm2: np.ndarray
    I0_Wcm2: float
    omg_ts: float
    omg_pe: float


# ---------------------------
# Core integrator
# ---------------------------

def evolve_temperature(params: SimulationParams) -> SimulationOutputs:
    """Integrate temperature in time for given parameters.

    Returns arrays for time, Te(t), and I(t), and key derived scalars.
    """
    # Laser base quantities
    I0, area = laser_intensity(params.E_J, params.tau_s, params.diameter_cm)
    omg_ts = angular_frequency(params.lambda_cm)

    # Plasma frequency
    omg_pe = plasma_frequency(params.ne_cm3)

    # Time grid
    time = np.linspace(0.0, params.t_end_s, params.n_steps)
    dt = time[1] - time[0]

    # Laser temporal profile
    center = params.t_end_s / 2.0
    profile = tgaussian(time, fwhm=params.tau_s, order=params.order, center=center)
    I_t = I0 * profile

    # Integrate Te(t)
    Te = np.empty(params.n_steps + 1, dtype=float)
    Te[0] = params.Te0_eV

    for i in range(params.n_steps):
        nuB = nu_ib(Te[i], params.ne_cm3, params.Z0, omg_ts, omg_pe)
        # RHS: dTe/dt = [ν_B * I(t) / (1.5 * n_e)] / (J/eV)
        dTe_dt = (nuB * I_t[i] / (1.5 * params.ne_cm3)) / J_PER_EV
        Te[i + 1] = Te[i] + dTe_dt * dt

    return SimulationOutputs(
        time_s=time,
        Te_eV=Te,  # length n_steps + 1
        I_Wcm2=I_t,
        I0_Wcm2=I0,
        omg_ts=omg_ts,
        omg_pe=omg_pe,
    )


# ---------------------------
# Main (script-style execution)
# ---------------------------

def main() -> None:
    params = SimulationParams()
    out = evolve_temperature(params)

    print(f"Laser intensity for TS = {out.I0_Wcm2:.2e} W/cm2")
    print(f"Laser angular frequency for TS = {out.omg_ts:.2e} rad/s")
    print(f"plasma electron angular frequency = {out.omg_pe:.2e} rad/s")
    print(f"final Te = {out.Te_eV[-1]:.1f} eV")


if __name__ == "__main__":
    main()
