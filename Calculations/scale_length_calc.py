import math

def compute_scale_length(x1_um, n1_cc, x2_um, n2_cc):
    """
    Compute the electron density scale length L (in microns)
    assuming an exponential density profile: n_e(x) = n0 * exp(x / L)

    Parameters:
        x1_um (float): Position x1 in microns
        n1_cc (float): Electron density at x1 in cm^-3
        x2_um (float): Position x2 in microns
        n2_cc (float): Electron density at x2 in cm^-3

    Returns:
        float: Scale length L in microns
    """
    if n1_cc <= 0 or n2_cc <= 0:
        raise ValueError("Electron densities must be positive numbers.")

    delta_x = x2_um - x1_um
    log_ratio = math.log(n2_cc / n1_cc)
    L_um = delta_x / log_ratio
    return L_um

# Example usage:
if __name__ == "__main__":
    x1 = -20.0       # in microns
    n1 = 1e20        # in cm^-3
    x2 = 0.0       # in microns
    n2 = 5e22        # in cm^-3

    L = compute_scale_length(x1, n1, x2, n2)
    print(f"Scale length L = {L:.4f} microns")