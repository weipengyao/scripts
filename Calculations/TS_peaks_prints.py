import math
from typing import Tuple

def _alpha_SI(lambda_i_m: float, theta_deg: float, ne_m3: float, Te_eV: float) -> float:
    """Compute α = 1/(k λ_De) using SI (Te in eV)."""
    if lambda_i_m <= 0 or ne_m3 <= 0 or Te_eV <= 0:
        raise ValueError("lambda_i_m, ne_m3, and Te_eV must be positive.")
    eps0 = 8.854_187_8128e-12  # F/m
    e    = 1.602_176_634e-19   # C
    s = math.sin(math.radians(theta_deg) / 2.0)
    if s <= 0:
        raise ValueError("theta_deg must be in (0, 180].")
    # k = (4.0 * math.pi / lambda_i_m) * s                # m^-1
    k = 2. * math.pi / lambda_i_m   # m^-1, equivalent and more stable
    print('k={:.4e} m^-1'.format(k))
    lam_De = math.sqrt(eps0 * e * Te_eV / (ne_m3 * e**2))  # m
    print('lam_De={:.4e} m'.format(lam_De))
    return 1.0 / (k * lam_De)

def ts_regime_SI(lambda_i_m: float,
                 theta_deg: float,
                 ne_m3: float,
                 Te_eV: float) -> str:
    """
    Determine the Thomson-scattering regime from α=1/(k λ_De).
      α < 0.1 -> 'non-collective'
      α ≥ 1   -> 'collective'
      else    -> 'mixed'
    """
    alpha = _alpha_SI(lambda_i_m, theta_deg, ne_m3, Te_eV)
    if alpha < 0.1:
        return "non-collective"
    elif alpha >= 1.0:
        return "collective"
    else:
        return "mixed"

def estimate_Te_from_ion_feature(delta_lambda_m: float,
                                 lambda_i_m: float,
                                 theta_deg: float,
                                 ne_m3: float,
                                 Z: float,
                                 m_i_kg: float,
                                 Ti_over_Te: float = 0.1,
                                 tol: float = 1e-8,
                                 max_iter: int = 200) -> float:
    """
    Estimate electron temperature Te (in eV) from the ion-feature peak separation Δλ,
    with fixed Ti/Te (default 0.1). SI for all inputs except temperatures (in eV).
    """
    c    = 299_792_458.0
    e    = 1.602_176_634e-19
    eps0 = 8.854_187_8128e-12
    pi   = math.pi

    if delta_lambda_m <= 0 or lambda_i_m <= 0:
        raise ValueError("delta_lambda_m and lambda_i_m must be positive.")
    if ne_m3 <= 0 or Z <= 0 or m_i_kg <= 0:
        raise ValueError("ne_m3, Z, and m_i_kg must be positive.")

    s = math.sin(math.radians(theta_deg) / 2.0)
    if s <= 0.0:
        raise ValueError("theta_deg must be in (0, 180].")

    frac = delta_lambda_m / lambda_i_m
    pref = (4.0 / c) * s
    k_iaw = (4.0 * pi / lambda_i_m) * s
    print('k_iaw={:.4e} m^-1'.format(k_iaw))

    def model_frac(Te_eV: float) -> float:
        kBTe_J = e * Te_eV
        lam_De = math.sqrt(eps0 * kBTe_J / (ne_m3 * e**2))
        denom  = 1.0 + (k_iaw * lam_De) ** 2
        shape  = (Z / denom) + 3.0 * Ti_over_Te
        return pref * math.sqrt((kBTe_J / m_i_kg) * shape)

    # Bisection solve model_frac(Te) = frac
    lo, hi = 1e-3, 1e5  # eV
    f_lo, f_hi = model_frac(lo) - frac, model_frac(hi) - frac
    for _ in range(20):
        if f_lo * f_hi <= 0:
            break
        if f_lo > 0:
            lo *= 0.1; f_lo = model_frac(lo) - frac
        else:
            hi *= 10.0; f_hi = model_frac(hi) - frac
    if f_lo * f_hi > 0:
        raise RuntimeError("Failed to bracket a solution for Te. Check inputs.")

    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        f_mid = model_frac(mid) - frac
        if abs(f_mid) < tol:
            return mid
        if f_lo * f_mid <= 0:
            hi, f_hi = mid, f_mid
        else:
            lo, f_lo = mid, f_mid
    return 0.5 * (lo + hi)


def estimate_Te_and_regime(delta_lambda_m: float,
                           lambda_i_m: float,
                           theta_deg: float,
                           ne_m3: float,
                           Z: float,
                           m_i_kg: float,
                           Ti_over_Te: float = 0.1) -> Tuple[float, str]:
    """
    Parameters
    ----------
    delta_lambda_m : float
        Peak-to-peak wavelength separation Δλ of the ion-acoustic feature in meters.
        This is the measured distance between the two ion-feature peaks in the TS
        spectrum.
    lambda_i_m : float
        Wavelength of the incident probe beam in meters.
    theta_deg : float
        Scattering angle θ in degrees, defined between incident (k_i) and scattered
        (k_s) wave vectors.
    ne_m3 : float
        Electron number density in m^-3.
    Z : float
        Mean ion charge state (dimensionless).
    m_i_kg : float
        Ion mass in kilograms.
    Ti_over_Te : float, optional
        Fixed ion-to-electron temperature ratio Ti/Te (default = 0.1, i.e., Te/Ti = 10).

    Returns
    -------
    Te_eV : float
        Inferred electron temperature in eV from Eq. (2.18).
    regime : str
        String `"collective"` if the regime check passes.
    """
    alpha_highT = _alpha_SI(lambda_i_m, theta_deg, ne_m3, Te_eV=1e4)   # 10 keV
    alpha_lowT  = _alpha_SI(lambda_i_m, theta_deg, ne_m3, Te_eV=1.0)   # 1 eV

    if alpha_lowT < 0.1:
        raise RuntimeError("TS is definitively non-collective for Te≈1 eV; "
                           "ion-feature formula (Eq. 2.18) is not applicable.")

    Te_eV = estimate_Te_from_ion_feature(delta_lambda_m, lambda_i_m, theta_deg,
                                         ne_m3, Z, m_i_kg, Ti_over_Te=Ti_over_Te)

    regime = ts_regime_SI(lambda_i_m, theta_deg, ne_m3, Te_eV)
    if regime == "non-collective":
        raise RuntimeError(f"Regime is '{regime}', not 'mixed' or 'collective'; "
                           "We can't extract Te.")

    return Te_eV, regime


attempt_Te, attempt_regime = estimate_Te_and_regime(
    delta_lambda_m = 0.2421e-9,   # example Δλ in meters
    lambda_i_m     = 526.5e-9,    # probe wavelength in meters
    theta_deg      = 90.0,      # scattering angle
    ne_m3          = 1.0e24,    # electron density in m^-3
    Z              = 1.0,       # mean ion charge
    m_i_kg         = 1.6726219e-27,  # proton mass
    Ti_over_Te     = 0.1
)
print(attempt_regime)
print('Te = {:.2e}'.format(attempt_Te))

# def test1():
#     # Expected values from a known calculation
#     expect_Te      = 111.473605544
#     expect_regime  = "collective"

#     attempt_Te, attempt_regime = estimate_Te_and_regime(
#         delta_lambda_m = 5.53e-10,   # example Δλ in meters
#         lambda_i_m     = 532e-9,    # probe wavelength in meters
#         theta_deg      = 90.0,      # scattering angle
#         ne_m3          = 1.0e25,    # electron density in m^-3
#         Z              = 1.0,       # mean ion charge
#         m_i_kg         = 1.6726219e-27,  # proton mass
#         Ti_over_Te     = 0.1
#     )
#     print(expect_regime)
#     print('Te = {:.2e}'.format(attempt_Te))

# test1()
