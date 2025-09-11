import numpy as np
from astropy import units as u
from astropy.constants import G, M_sun, R_earth

# White dwarf: 0.6 Msun, 1.0 R_earth
M_star = 0.6 * M_sun
R_star = 1.0 * R_earth

def _qty(x, assumed_unit):
    """Attach assumed_unit if x is unitless; pass through if already a Quantity."""
    q = u.Quantity(x)
    return q if q.unit != u.dimensionless_unscaled else (q.value * assumed_unit)

def n_pts_in_transit(P, k, inc, W=30*u.day, cadence=200*u.s):
    """
    Circular orbit. Vectorized.
    - P: period (Quantity or number; if number -> days)
    - k: Rp/R* (float/array)
    - inc: inclination (Quantity or number; if number -> degrees)
    - W: observing window (Quantity), default 30 d
    - cadence: sampling cadence (Quantity), default 200 s

    Returns: number of in-transit points (dimensionless)
    """
    # Attach default units if needed
    P = _qty(P,  u.day).to(u.s)
    i = _qty(inc, u.deg).to(u.rad)

    a = (G * M_star * P**2 / (4*np.pi**2))**(1/3)   # semi-major axis

    # impact parameter (circular)
    b = (a / R_star) * np.cos(i.value)

    # geometry factor; guard small sin(i)
    sin_i = np.maximum(np.sin(i.value), 1e-12)
    inside = (1.0 + k)**2 - b**2
    geom = np.sqrt(np.maximum(0.0, inside)) / sin_i

    # total in-transit time across window W
    T_tot = (W.to(u.s) / np.pi) * (R_star / a) * geom

    # zero if it misses
    T_tot = T_tot * (b <= (1.0 + k))

    return (T_tot / cadence).to(u.dimensionless_unscaled)

if __name__ == "__main__":
    points = n_pts_in_transit(1, 5, 89.5)
    print(points)