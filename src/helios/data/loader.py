"""Turn raw HAPI tables into the two standard frames the analysis uses.

    mag     columns bx, by, bz           (nT, GSE)
    plasma  columns vx, vy, vz, n, tp    (km/s GSE, cm^-3, K)

Both are indexed by UTC time at each instrument's native cadence. They are
deliberately NOT merged onto one time base: the field is sampled 4 to 30
times faster than the plasma, and window means are taken per instrument, as
the Heliospheric Shock Database does.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from ..constants import KB, KMS, MP
from ..physics.rankine_hugoniot import PlasmaState


@dataclass(frozen=True)
class EventData:
    spacecraft: str
    mag: pd.DataFrame
    plasma: pd.DataFrame


def thermal_speed_to_temperature(w_kms):
    """Proton temperature from the most probable thermal speed: T = m_p w^2 / (2 k)."""
    return MP * (np.asarray(w_kms, dtype=float) * KMS) ** 2 / (2 * KB)


def standardize_mag(raw: pd.DataFrame, vector: str) -> pd.DataFrame:
    out = raw[[f"{vector}_0", f"{vector}_1", f"{vector}_2"]].copy()
    out.columns = ["bx", "by", "bz"]
    return out.dropna()


def standardize_plasma(raw: pd.DataFrame, spacecraft: str) -> pd.DataFrame:
    out = pd.DataFrame(index=raw.index)
    out[["vx", "vy", "vz"]] = raw[["V_GSE_0", "V_GSE_1", "V_GSE_2"]]
    out["n"] = raw["Np"]
    if spacecraft == "wind":
        out["tp"] = thermal_speed_to_temperature(raw["THERMAL_SPD"])  # WI_K0_SWE gives a thermal speed
    elif spacecraft == "ace":
        out["tp"] = raw["Tpr"]  # AC_H0_SWE gives the radial proton temperature in K
    elif spacecraft == "dscovr":
        out["tp"] = raw["THERMAL_TEMP"]
    else:
        raise ValueError(f"no plasma mapping for {spacecraft!r}")
    return out


def window_state(data: EventData, start: pd.Timestamp, stop: pd.Timestamp, min_points: int = 3) -> PlasmaState | None:
    """Mean state over [start, stop], each instrument averaged at its own cadence.

    Plasma moments are averaged quantity by quantity, ignoring gaps in each
    separately (a missing temperature does not discard a valid density).
    Returns None if any quantity has fewer than `min_points` samples.
    """
    m = data.mag.loc[start:stop]
    p = data.plasma.loc[start:stop]
    if len(m) < min_points or p.notna().sum().min() < min_points:
        return None
    pm = p.mean()
    xyz = m[["bx", "by", "bz"]]
    # Binned field data carry the mean of |B| per bin; native data need it computed.
    b_mag = float(m["bmag"].mean()) if "bmag" in m else float(np.linalg.norm(xyz.to_numpy(), axis=1).mean())
    return PlasmaState(b=xyz.mean().to_numpy(), v=pm[["vx", "vy", "vz"]].to_numpy(), n=float(pm["n"]), tp=float(pm["tp"]), b_mag=b_mag)


def shock_windows(data: EventData, t: pd.Timestamp, shock_type: str):
    """Upstream and downstream states using the Heliospheric Shock Database default windows."""
    minute = pd.Timedelta(minutes=1)
    if shock_type == "FF":
        up = window_state(data, t - 9 * minute, t - 1 * minute)
        down = window_state(data, t + 2 * minute, t + 10 * minute)
    else:
        up = window_state(data, t + 1 * minute, t + 9 * minute)
        down = window_state(data, t - 10 * minute, t - 2 * minute)
    return up, down
