"""Rule-based shock detector: the baseline every later method must beat.

It applies the Heliospheric Shock Database selection criteria at every time
step, comparing means over that database's default analysis windows
(upstream 9 to 1 min before, downstream 2 to 10 min after, mirrored for
reverse shocks), then keeps one detection per cluster of triggers.

Known weakness, by design: the windows are 8 minutes long and separated by a
3-minute gap, so a gradual compression over ~10 minutes (a stream interaction
region, a density pile-up) satisfies the criteria as well as a discontinuity
does. A sharpness score is reported but not thresholded here; separating
real shocks from such look-alikes is the job of the later stages.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from ..physics.rankine_hugoniot import PlasmaState
from ..physics.shock import ShockParameters, analyze_shock

# Analysis windows in seconds, relative to the candidate time (Helsinki defaults).
NEAR_GAP_FIRST = 60.0  # excluded on the upstream side
NEAR_GAP_SECOND = 120.0  # excluded on the downstream side
WINDOW = 480.0  # length of each averaging window


@dataclass(frozen=True)
class Detection:
    time: float  # same units as the input index (seconds)
    params: ShockParameters
    sharpness: float  # fraction of the |B| jump that occurs within +-1 sample of the time
    upstream: PlasmaState
    downstream: PlasmaState


def _state(df: pd.DataFrame) -> PlasmaState | None:
    if len(df) < 3 or df.isna().any().any():
        df = df.dropna()
        if len(df) < 3:
            return None
    m = df.mean()
    return PlasmaState(b=m[["bx", "by", "bz"]].to_numpy(), v=m[["vx", "vy", "vz"]].to_numpy(), n=float(m["n"]), tp=float(m["tp"]))


def windows(df: pd.DataFrame, t: float, shock_type: str) -> tuple[PlasmaState | None, PlasmaState | None]:
    """Upstream and downstream mean states for a candidate at time t."""
    before_up = df.loc[t - NEAR_GAP_FIRST - WINDOW: t - NEAR_GAP_FIRST]
    after_down = df.loc[t + NEAR_GAP_SECOND: t + NEAR_GAP_SECOND + WINDOW]
    after_up = df.loc[t + NEAR_GAP_FIRST: t + NEAR_GAP_FIRST + WINDOW]
    before_down = df.loc[t - NEAR_GAP_SECOND - WINDOW: t - NEAR_GAP_SECOND]
    if shock_type == "FF":
        return _state(before_up), _state(after_down)
    return _state(after_up), _state(before_down)


def detect(df: pd.DataFrame, te_over_tp: float = 1.0, merge_s: float = 900.0) -> list[Detection]:
    """Detections in a time series indexed by seconds with the standard columns.

    Triggers closer together than `merge_s` are treated as one event and the
    time of the steepest |B| change among them is reported.
    """
    t = df.index.to_numpy(dtype=float)
    bmag = np.linalg.norm(df[["bx", "by", "bz"]].to_numpy(), axis=1)
    grad = np.abs(np.gradient(bmag))
    triggers: list[tuple[float, str, float]] = []
    for i, ti in enumerate(t):
        for kind in ("FF", "FR"):
            up, down = windows(df, ti, kind)
            if up is None or down is None:
                continue
            p = analyze_shock(up, down, te_over_tp=te_over_tp)
            if p.meets_criteria and p.shock_type == kind:
                triggers.append((ti, kind, grad[i]))

    detections: list[Detection] = []
    cluster: list[tuple[float, str, float]] = []

    def flush():
        if not cluster:
            return
        ti, kind, _ = max(cluster, key=lambda c: c[2])
        up, down = windows(df, ti, kind)
        params = analyze_shock(up, down, te_over_tp=te_over_tp)
        i = int(np.argmin(np.abs(t - ti)))
        lo, hi = max(i - 1, 0), min(i + 1, len(t) - 1)
        total = abs(np.linalg.norm(down.b) - np.linalg.norm(up.b))
        sharp = float(abs(bmag[hi] - bmag[lo]) / total) if total > 0 else 0.0
        detections.append(Detection(time=float(ti), params=params, sharpness=min(sharp, 1.0), upstream=up, downstream=down))

    for trig in triggers:
        if cluster and trig[0] - cluster[-1][0] > merge_s:
            flush()
            cluster = []
        cluster.append(trig)
    flush()
    return detections
