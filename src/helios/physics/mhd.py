"""Characteristic speeds of a magnetised plasma, in practical units.

Assumption (stated once, used everywhere): the plasma is protons plus
electrons, quasi-neutral, so mass density is n * m_p and thermal pressure is
n k (T_p + T_e). Alpha particles are neglected; where they carry ~4% of the
number density this underestimates the mass density by ~15%.
"""

from __future__ import annotations

import numpy as np

from ..constants import CC, GAMMA, KB, KMS, MP, MU0, NT


def mass_density(n_cc):
    """Mass density in kg m^-3 from proton number density in cm^-3."""
    return np.asarray(n_cc, dtype=float) * CC * MP


def thermal_pressure(n_cc, tp_k, te_k=None):
    """Thermal pressure n k (T_p + T_e) in Pa. T_e defaults to T_p."""
    tp = np.asarray(tp_k, dtype=float)
    te = tp if te_k is None else np.asarray(te_k, dtype=float)
    return np.asarray(n_cc, dtype=float) * CC * KB * (tp + te)


def alfven_speed(b_nt, n_cc):
    """Alfven speed |B| / sqrt(mu0 rho) in km/s."""
    b = np.asarray(b_nt, dtype=float) * NT
    return b / np.sqrt(MU0 * mass_density(n_cc)) / KMS


def sound_speed(tp_k, te_k=None, gamma=GAMMA):
    """Sound speed sqrt(gamma k (T_p + T_e) / m_p) in km/s. T_e defaults to T_p."""
    tp = np.asarray(tp_k, dtype=float)
    te = tp if te_k is None else np.asarray(te_k, dtype=float)
    return np.sqrt(gamma * KB * (tp + te) / MP) / KMS


def magnetosonic_speed(v_a, c_s):
    """Fast speed for perpendicular propagation, sqrt(vA^2 + cs^2)."""
    return np.sqrt(np.asarray(v_a) ** 2 + np.asarray(c_s) ** 2)


def fast_speed(v_a, c_s, theta_rad):
    """Fast magnetosonic speed for propagation at angle theta to B.

    v_f^2 = 1/2 [ (vA^2 + cs^2) + sqrt((vA^2 + cs^2)^2 - 4 vA^2 cs^2 cos^2 theta) ]
    """
    s = np.asarray(v_a) ** 2 + np.asarray(c_s) ** 2
    disc = np.maximum(s**2 - 4 * np.asarray(v_a) ** 2 * np.asarray(c_s) ** 2 * np.cos(theta_rad) ** 2, 0.0)
    return np.sqrt(0.5 * (s + np.sqrt(disc)))


def plasma_beta(n_cc, tp_k, b_nt, te_k=None):
    """Ratio of thermal to magnetic pressure."""
    b = np.asarray(b_nt, dtype=float) * NT
    return thermal_pressure(n_cc, tp_k, te_k) / (b**2 / (2 * MU0))
