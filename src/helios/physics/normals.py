"""Shock normal estimates from upstream and downstream averages.

All are consequences of the Rankine-Hugoniot conditions and need only the
mean fields on either side (Paschmann & Daly, "Analysis Methods for
Multi-Spacecraft Data", ISSI SR-001, ch. 10; Schwartz 1998). With
dB = B_d - B_u and dV = V_d - V_u:

    magnetic coplanarity   MC   n ∝ (B_d × B_u) × dB
    velocity coplanarity   VC   n ∝ dV                 (approximate; best at high Mach number)
    mixed mode 1           MX1  n ∝ (B_u × dV) × dB
    mixed mode 2           MX2  n ∝ (B_d × dV) × dB
    mixed mode 3           MX3  n ∝ (dB × dV) × dB     (Abraham-Shrauner & Yun 1976)

MC fails for exactly parallel or perpendicular shocks, where B_u, B_d and dB
are collinear or the cross product vanishes. The mixed modes stay well
conditioned there, which is why MX3 is the usual default.

Orientation: every normal is returned pointing into the upstream medium.
Mass conservation gives V_u.n - V_sh = N_d (V_u - V_d).n / (N_d - N_u), which
must be negative (plasma flows from upstream into the shock), so for a
compressive shock n is chosen with dV.n > 0.
"""

from __future__ import annotations

import numpy as np

METHODS = ("MC", "VC", "MX1", "MX2", "MX3")


def _unit(v):
    norm = np.linalg.norm(v)
    if norm == 0 or not np.isfinite(norm):
        return np.full(3, np.nan)
    return v / norm


def shock_normal(b_up, b_down, v_up, v_down, method: str = "MX3") -> np.ndarray:
    """Unit shock normal pointing upstream, by the named method."""
    bu, bd = np.asarray(b_up, float), np.asarray(b_down, float)
    vu, vd = np.asarray(v_up, float), np.asarray(v_down, float)
    db, dv = bd - bu, vd - vu
    if method == "MC":
        n = np.cross(np.cross(bd, bu), db)
    elif method == "VC":
        n = dv
    elif method == "MX1":
        n = np.cross(np.cross(bu, dv), db)
    elif method == "MX2":
        n = np.cross(np.cross(bd, dv), db)
    elif method == "MX3":
        n = np.cross(np.cross(db, dv), db)
    else:
        raise ValueError(f"Unknown method {method!r}; choose from {METHODS}")
    n = _unit(n)
    if np.all(np.isfinite(n)) and n @ dv < 0:
        n = -n
    return n


def angle_between_deg(a, b) -> float:
    """Angle between two directions, ignoring sign, in degrees (0-90)."""
    a, b = _unit(np.asarray(a, float)), _unit(np.asarray(b, float))
    return float(np.degrees(np.arccos(np.clip(abs(a @ b), 0, 1))))
