import math

def sigma_from_phi(phi_deg, M0, beta, gamma=5/3):
    """
    Solve for shock angle σ (deg) given deflection φ (deg),
    total upstream Mach M0, plasma beta β, and γ.
    Uses Eq. (6) for r with M_{n,0}=M0*sinσ and Eq. (8) for geometry.
    """
    phi = math.radians(phi_deg)

    def r_from_sigma(sigma):
        Mn0 = M0 * math.sin(sigma)
        A = 2*(2 - gamma)
        B = 2*gamma*(beta + 1) + beta*gamma*(gamma - 1)*Mn0**2
        C = -beta*gamma*(gamma + 1)*Mn0**2
        disc = B*B - 4*A*C
        r1 = (-B + math.sqrt(disc)) / (2*A)
        r2 = (-B - math.sqrt(disc)) / (2*A)
        return max(r1, r2)     # pick the physical (>1) branch

    def f(sigma):
        r = r_from_sigma(sigma)
        t = math.tan(sigma)
        return math.tan(phi) - t*(r - 1)/(r + t*t)  # Eq. (8)

    # bracket the physical oblique-shock range (Mn0>1 ⇒ σ not too small)
    a, b = math.radians(30.0), math.radians(80.0)
    fa, fb = f(a), f(b)
    # bisection
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

# Your case: φ=20°, M0=2.81, β=2, γ=5/3
sigma_deg, r = sigma_from_phi(20.0, 2.81, 2.0, gamma=5/3)
print(f"sigma ≈ {sigma_deg:.2f}°,  r ≈ {r:.3f}")