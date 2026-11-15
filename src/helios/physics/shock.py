"""Shock parameters from upstream and downstream averages.

Selection criteria and conventions follow the Heliospheric Shock Database
(University of Helsinki; Kilpua et al. 2015; method documentation of
2024-09-06), so that results can be compared with that catalogue like for like:

    B_d / B_u >= 1.2,   N_d / N_u >= 1.2,   T_d / T_u >= 1/1.2,
    |V| jump >= 20 km/s  (V_d > V_u for fast forward, V_u > V_d for fast reverse)

    normal: mixed mode MX3, falling back to magnetic coplanarity
    shock speed: V_sh = (N_d V_d - N_u V_u) . n / (N_d - N_u)      (mass flux)
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

from ..constants import GAMMA
from . import mhd
from .normals import shock_normal
from .rankine_hugoniot import PlasmaState

B_RATIO_MIN = 1.2
N_RATIO_MIN = 1.2
T_RATIO_MIN = 1 / 1.2
DV_MIN_KMS = 20.0


@dataclass(frozen=True)
class ShockParameters:
    shock_type: str  # "FF", "FR" or "none"
    meets_criteria: bool
    b_ratio: float
    n_ratio: float
    t_ratio: float
    dv_kms: float  # |V_d| - |V_u|
    normal: tuple[float, float, float]
    normal_method: str
    theta_bn_deg: float
    shock_speed_kms: float  # along the normal, observer frame
    upstream_normal_speed_kms: float  # |V_u.n - V_sh|: inflow speed in the shock frame
    alfven_speed_kms: float
    sound_speed_kms: float
    mach_alfven: float
    mach_magnetosonic: float  # against sqrt(vA^2 + cs^2)
    mach_fast: float  # against the fast speed along the normal
    plasma_beta: float

    def as_dict(self) -> dict:
        return asdict(self)


def analyze_shock(up: PlasmaState, down: PlasmaState, te_over_tp: float = 1.0,
                  gamma: float = GAMMA, method: str = "MX3") -> ShockParameters:
    """Classify a jump between two plasma states and derive its shock parameters."""
    bu, bd = np.linalg.norm(up.b), np.linalg.norm(down.b)
    vu, vd = np.linalg.norm(up.v), np.linalg.norm(down.v)
    b_ratio, n_ratio, t_ratio = bd / bu, down.n / up.n, down.tp / up.tp
    dv = vd - vu

    jumps_ok = b_ratio >= B_RATIO_MIN and n_ratio >= N_RATIO_MIN and t_ratio >= T_RATIO_MIN
    if dv >= DV_MIN_KMS:
        shock_type = "FF"
    elif -dv >= DV_MIN_KMS:
        shock_type = "FR"
    else:
        shock_type = "none"
    meets = bool(jumps_ok and shock_type != "none")

    n = shock_normal(up.b, down.b, up.v, down.v, method)
    used = method
    if not np.all(np.isfinite(n)):
        n = shock_normal(up.b, down.b, up.v, down.v, "MC")
        used = "MC"

    te = up.tp * te_over_tp
    va = float(mhd.alfven_speed(bu, up.n))
    cs = float(mhd.sound_speed(up.tp, te, gamma))
    beta = float(mhd.plasma_beta(up.n, up.tp, bu, te))

    if np.all(np.isfinite(n)) and down.n != up.n:
        theta = np.degrees(np.arccos(np.clip(abs(up.b @ n) / bu, 0, 1)))
        v_sh = float((down.n * down.v - up.n * up.v) @ n / (down.n - up.n))
        un = abs(float(up.v @ n) - v_sh)
        vf = float(mhd.fast_speed(va, cs, np.radians(theta)))
        vms = float(mhd.magnetosonic_speed(va, cs))
        machs = (un / va, un / vms, un / vf)
    else:
        theta, v_sh, un, machs = np.nan, np.nan, np.nan, (np.nan, np.nan, np.nan)

    return ShockParameters(
        shock_type=shock_type if meets else "none", meets_criteria=meets,
        b_ratio=float(b_ratio), n_ratio=float(n_ratio), t_ratio=float(t_ratio), dv_kms=float(dv),
        normal=tuple(float(x) for x in n), normal_method=used, theta_bn_deg=float(theta),
        shock_speed_kms=v_sh, upstream_normal_speed_kms=float(un),
        alfven_speed_kms=va, sound_speed_kms=cs,
        mach_alfven=float(machs[0]), mach_magnetosonic=float(machs[1]), mach_fast=float(machs[2]),
        plasma_beta=beta,
    )
