"""
Shared dielectric parameter calculations.

This module centralizes dielectric formulas so instrument drivers do not
maintain duplicate implementations that can drift over time.
"""

VACUUM_PERMITTIVITY_F_PER_M = 8.854e-12


def compute_eps_r_from_area(capacitance, thickness_nm, area_um2):
    """Compute relative permittivity (epsilon_r) from capacitance.

    epsilon_r = C * t / (epsilon_0 * A)

    Args:
        capacitance: Scalar or array-like capacitance in F.
        thickness_nm: Dielectric thickness in nm.
        area_um2: Electrode area in um^2.
    """
    t_m = thickness_nm * 1e-9
    area_m2 = area_um2 * 1e-12
    return capacitance * t_m / (VACUUM_PERMITTIVITY_F_PER_M * area_m2)
