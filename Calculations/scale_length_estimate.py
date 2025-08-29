import numpy as np

def density_scale_length(x1_um, n1_cc, x2_um, n2_cc):
    """
    Compute the exponential density scale length L_n from two points.

    Parameters
    ----------
    x1_um : float
        Position 1 in microns
    n1_cc : float
        Electron density at x1 in 1/cc
    x2_um : float
        Position 2 in microns
    n2_cc : float
        Electron density at x2 in 1/cc

    Returns
    -------
    L_n : float
        Exponential density scale length in microns
    """
    if n1_cc <= 0 or n2_cc <= 0:
        raise ValueError("Densities must be positive values.")

    return (x2_um - x1_um) / np.log(n2_cc / n1_cc)

# Example usage:
x1, n1 = -146.0, 1e19   # micron, 1/cc
x2, n2 = -115.0, 1.0e22   # micron, 1/cc

L_n = density_scale_length(x1, n1, x2, n2)
print(f"Exponential scale length L_n = {L_n:.2f} micron")
