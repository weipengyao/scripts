def shock_characterization(laser_intensity, pulse_duration, mass_density, skin_depth, phi_deg, M0, beta, reflectivity=1.0, c=3e8, gamma=5/3):
    """
    Calculate the shock pressure and velocity in a material subjected to laser irradiation, as well as the shock angle sigma baesed on the given flow deflection angle phi, Mach number M, plasma beta beta, and adiabatic index gamma.

    Parameters:
    laser_intensity (float): Laser intensity in W/cm^2.
    pulse_duration (float): Pulse duration in fs.
    mass_density (float): Mass density of the material in kg/m^3.
    skin_depth (float): Skin depth of the material in nm.
    reflectivity (float): Reflectivity of the material (default is 1.0).
    c (float): Speed of light in m/s (default is 3e8 m/s).
    phi (float): Flow deflection angle in degrees.
    M (float): Mach number.
    beta (float): Plasma beta.
    gamma (float): Adiabatic index (default is 5/3).

    Returns:
    tuple: Shock pressure in bar and shock velocity in m/s.
    """

    def shock_velocity(laser_intensity, pulse_duration, mass_density, skin_depth,
                       reflectivity=1.0, c=299_792_458.0):
        """
        Compute the shock velocity driven by radiation pressure of an ultrashort laser.
        """
        # unit conversions
        I0_W_per_m2 = laser_intensity * 1e4
        tau_s = pulse_duration * 1e-15
        delta_m = skin_depth * 1e-9

        P_rad_Pa = (1.0 + reflectivity) * I0_W_per_m2 / c          # radiation pressure (Pa = N/m^2)
        v = P_rad_Pa * tau_s / (mass_density * delta_m)
        return v

    def sigma_from_phi(phi_deg, M0, beta, gamma=5/3):
        import math

        phi = math.radians(phi_deg)

        def r_from_sigma(sigma):
            """
            Compute density compression ratio r given shock angle σ [radians].
            """
            Mn0 = M0 * math.sin(sigma)  # normal Mach number
            A = 2*(2 - gamma)
            B = 2*gamma*(beta + 1) + beta*gamma*(gamma - 1)*Mn0**2
            C = -beta*gamma*(gamma + 1)*Mn0**2
            disc = B*B - 4*A*C
            r1 = (-B + math.sqrt(disc)) / (2*A)
            r2 = (-B - math.sqrt(disc)) / (2*A)
            return max(r1, r2)  # pick the physically relevant branch (r > 1)

        def f(sigma):
            """
            Root of this function corresponds to the correct σ.
            """
            r = r_from_sigma(sigma)
            t = math.tan(sigma)
            return math.tan(phi) - t*(r - 1)/(r + t*t)

        # Bracket the physical solution within typical oblique shock range
        a, b = math.radians(20.0), math.radians(60.0)
        fa, fb = f(a), f(b)

        # Solve f(σ) = 0 by bisection
        for _ in range(80):
            m = 0.5*(a + b)
            fm = f(m)
            if fa*fm <= 0:
                b, fb = m, fm
            else:
                a, fa = m, fm

        sigma = 0.5*(a + b)
        r = r_from_sigma(sigma)
        return math.degrees(sigma), r

    shock_velocity = shock_velocity(laser_intensity, pulse_duration, mass_density, skin_depth, reflectivity, c)
    shock_angle, compression_ratio = sigma_from_phi(phi_deg, M0, beta, gamma=5/3)

    return shock_velocity, shock_angle

# if __name__ == "__main__":
#     # Example parameters
#     laser_intensity = 3e17  # W/cm^2
#     pulse_duration = 186     # fs
#     mass_density = 2.65e3    # kg/m^3 (Aluminum)
#     skin_depth = 6.4         # nm
#     phi_deg = 20            # degrees
#     M0 = 2.81                 # Mach number
#     beta = 2                # Plasma beta

#     shock_velocity, shock_angle = shock_characterization(laser_intensity, pulse_duration, mass_density, skin_depth, phi_deg, M0, beta)

#     print(f"Shock Velocity: {shock_velocity:.2f} m/s")
#     print(f"Shock Angle: {shock_angle:.2f} degrees")


# import pytest

# # keep the parameterized test function definition in the notebook
# @pytest.mark.parametrize(
#     "phi_deg, laser_intensity, pulse_duration, mass_density, skin_depth, reflectivity, expect_v, expect_sigma",
#     [
#         (20.0, 3e17, 186, 2.65e3, 6.4, 1.0, 219_339.6226415094, 55.1561928762),
#         (15.0, 1e17, 100, 2200.0, 10.0, 1.0, 30_303.030303030308, 44.3585025807),
#         (10.0, 5e16, 250, 2700.0, 5.0, 0.2, 37_037.03703703704, 37.2436435962),
#     ]
# )
# def test_shock_characterization(phi_deg, laser_intensity, pulse_duration,
#                                 mass_density, skin_depth, reflectivity,
#                                 expect_v, expect_sigma):
#     v, sigma = shock_characterization(
#         laser_intensity=laser_intensity,
#         pulse_duration=pulse_duration,
#         mass_density=mass_density,
#         skin_depth=skin_depth,
#         phi_deg=phi_deg,
#         M0=2.81,
#         beta=2.0,
#         reflectivity=reflectivity
#     )

#     assert abs(v - expect_v) <= 1e-3 * max(1.0, abs(expect_v))
#     assert abs(sigma - expect_sigma) <= 0.5
#     print(f"pass: φ={phi_deg}° → σ≈{expect_sigma:.1f}°, v≈{expect_v:.1f} m/s")

# # Run pytest inside Jupyter
# pytest.main(["-q", "--tb=short", "-s"])

# import math

def assert_close_rel(got, expect, rtol=1e-3, msg=""):
    assert abs(got - expect) <= rtol * max(1.0, abs(expect)), \
        f"{msg} got {got:.6g}, expect {expect:.6g} (rtol={rtol})"

def test_shock_phi_20_deg_example():
    # Example from the snippet
    v, sigma = shock_characterization(
        laser_intensity=3e17,   # W/cm^2
        pulse_duration=186,     # fs
        mass_density=2.65e3,    # kg/m^3
        skin_depth=6.4,         # nm
        phi_deg=20.0,           # deg
        M0=2.81,
        beta=2.0
    )
    expect_v = 219_339.6226415094   # m/s
    expect_sigma = 55.1561928762    # deg
    assert_close_rel(v, expect_v, 1e-3, "velocity:")
    assert abs(sigma - expect_sigma) <= 0.5, f"sigma got {sigma:.3f}, expect ~{expect_sigma:.3f}"

def test_shock_phi_15_deg_moderate():
    v, sigma = shock_characterization(
        laser_intensity=1e17,   # W/cm^2
        pulse_duration=100,     # fs
        mass_density=2200.0,    # kg/m^3
        skin_depth=10.0,        # nm
        phi_deg=15.0,           # deg
        M0=2.81,
        beta=2.0
    )
    expect_v = 30_303.030303030308  # m/s
    expect_sigma = 44.3585025807    # deg (~45° blue-curve target)
    assert_close_rel(v, expect_v, 1e-3, "velocity:")
    assert abs(sigma - expect_sigma) <= 0.5

def test_shock_phi_10_deg_lowR():
    v, sigma = shock_characterization(
        laser_intensity=5e16,   # W/cm^2
        pulse_duration=250,     # fs
        mass_density=2700.0,    # kg/m^3
        skin_depth=5.0,         # nm
        phi_deg=10.0,           # deg
        M0=2.81,
        beta=2.0,
        reflectivity=0.2
    )
    expect_v = 37_037.03703703704  # m/s
    expect_sigma = 37.2436435962   # deg (~37° blue-curve target)
    assert_close_rel(v, expect_v, 1e-3, "velocity:")
    assert abs(sigma - expect_sigma) <= 0.5

test_shock_phi_20_deg_example()
test_shock_phi_15_deg_moderate()
test_shock_phi_10_deg_lowR()
