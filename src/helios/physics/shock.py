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


HSD_TE_1AU_K = 146_277.0
HSD_TE_INDEX = -0.664


def hsd_electron_temperature(time=None, r_au: float = 1.0) -> dict:
    """Electron temperature assumption of the Heliospheric Shock Database, as keyword arguments.

    Electron temperatures are rarely measured well, so the database assumes
    T_e = 146 277 K * R^-0.664 (R in au; the 1 au value for near-Earth
    spacecraft). The sound speed, magnetosonic speed, beta and Mach numbers
    all depend on this choice; jump ratios, normals and shock speed do not.

    The method documentation (2024-09-06) also mentions T_e = 2 T_p for
    shocks before August 2015. The catalogue export does not follow that:
    for all 43 analysable Wind shocks of 2012 its sound speeds are reproduced
    to 0.0% with the fixed value and miss by 16% with 2 T_p. The fixed value
    is therefore used for every date. `time` is accepted for call-site
    compatibility and ignored.
    """
    return {"te_k": HSD_TE_1AU_K * r_au**HSD_TE_INDEX}


def analyze_shock(up: PlasmaState, down: PlasmaState, te_over_tp: float = 1.0, te_k: float | None = None,
                  gamma: float = GAMMA, method: str = "MX3") -> ShockParameters:
    """Classify a jump between two plasma states and derive its shock parameters.

    The electron temperature is either a fixed value `te_k` (kelvin) or, if
    that is None, `te_over_tp` times the upstream proton temperature.
    """
    # Field strengths use <|B|> (as the catalogue does); directions use the mean vector.
    bu, bd = up.field_strength, down.field_strength
    bu_vec = np.linalg.norm(up.b)
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

    te = up.tp * te_over_tp if te_k is None else float(te_k)
    va = float(mhd.alfven_speed(bu, up.n))
    cs = float(mhd.sound_speed(up.tp, te, gamma))
    beta = float(mhd.plasma_beta(up.n, up.tp, bu, te))

    if np.all(np.isfinite(n)) and down.n != up.n:
        theta = np.degrees(np.arccos(np.clip(abs(up.b @ n) / bu_vec, 0, 1)))
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
