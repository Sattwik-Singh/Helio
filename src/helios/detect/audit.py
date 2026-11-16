"""Physics audit: tests a candidate must pass to be called a shock.

The selection criteria ask only whether the mean state changed between two
8-minute windows. Two further questions separate shocks from look-alikes:

1. Is the upstream flow super-fast?  M_fast > 1.
   A fast shock cannot exist otherwise: the Rankine-Hugoniot relations have no
   compressive fast solution for sub-fast inflow. This is a requirement of
   the physics with no tunable threshold. It does depend on the assumed
   electron temperature, through the sound speed.

2. Is the change abrupt?  sharpness >= SHARPNESS_MIN.
   Sharpness is the fraction of the upstream-to-downstream change in |B|
   that happens within +-45 s of the reported time. A shock ramp is seconds
   wide, so the value is near 1; a compression that builds over ten minutes
   gives a value near 0. The threshold was chosen on Wind 2023 and is the one
   tuned number here, so it must be judged on other years.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from ..data.loader import EventData, shock_windows

SHARPNESS_MIN = 0.3
EDGE_S = 45.0
SIDE_S = 60.0


@dataclass(frozen=True)
class Audit:
    mach_fast: float
    sharpness: float
    super_fast: bool
    abrupt: bool

    @property
    def passed(self) -> bool:
        return self.super_fast and self.abrupt


def sharpness(data: EventData, t: pd.Timestamp, shock_type: str) -> float:
    """Fraction of the window-to-window |B| change that occurs within +-EDGE_S of t."""
    up, down = shock_windows(data, t, shock_type)
    if up is None or down is None:
        return float("nan")
    if "bmag" in data.mag:
        b = data.mag["bmag"]
    else:
        b = pd.Series(np.linalg.norm(data.mag[["bx", "by", "bz"]].to_numpy(), axis=1), index=data.mag.index)
    edge, side = pd.Timedelta(seconds=EDGE_S), pd.Timedelta(seconds=SIDE_S)
    before = b.loc[t - edge - side: t - edge].mean()
    after = b.loc[t + edge: t + edge + side].mean()
    total = down.field_strength - up.field_strength
    if not np.isfinite(before) or not np.isfinite(after) or total == 0:
        return float("nan")
    return float((after - before) / total if shock_type == "FF" else (before - after) / total)


def audit(data: EventData, detection) -> Audit:
    s = sharpness(data, detection.time, detection.params.shock_type)
    m = detection.params.mach_fast
    return Audit(mach_fast=m, sharpness=s, super_fast=bool(m > 1.0), abrupt=bool(np.isfinite(s) and s >= SHARPNESS_MIN))
