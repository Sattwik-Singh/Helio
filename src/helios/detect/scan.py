"""Baseline detector for real spacecraft data (field and plasma at different cadences).

Same rule as detect/baseline.py — the Heliospheric Shock Database criteria on
its default windows — but evaluated on a time grid over an EventData pair,
with the catalogue's electron temperature convention. The reported time is
refined to the steepest change in |B| at the magnetometer's native cadence.

This version slices pandas frames at every grid point: fine for hours of
data, too slow for years. A vectorised version is the first task of Phase 2.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from ..data.loader import EventData, shock_windows
from ..physics.shock import ShockParameters, analyze_shock, hsd_electron_temperature


@dataclass(frozen=True)
class RealDetection:
    time: pd.Timestamp
    params: ShockParameters
    n_triggers: int  # grid points in the cluster that satisfied the criteria


def _steepest_change(mag: pd.DataFrame, start: pd.Timestamp, stop: pd.Timestamp, half_s: float = 30.0) -> pd.Timestamp | None:
    """Time in [start, stop] where mean |B| just after differs most from mean |B| just before."""
    seg = mag.loc[start - pd.Timedelta(seconds=half_s): stop + pd.Timedelta(seconds=half_s)]
    if len(seg) < 8:
        return None
    b = pd.Series(np.linalg.norm(seg.to_numpy(), axis=1), index=seg.index)
    win = f"{int(half_s)}s"
    before = b.rolling(win, closed="left").mean()
    after = b[::-1].rolling(win, closed="left").mean()[::-1]
    jump = (after - before).abs().loc[start:stop].dropna()
    return jump.idxmax() if len(jump) else None


def scan(data: EventData, step_s: float = 30.0, merge_s: float = 900.0) -> list[RealDetection]:
    start = max(data.mag.index[0], data.plasma.index[0]) + pd.Timedelta(minutes=10)
    stop = min(data.mag.index[-1], data.plasma.index[-1]) - pd.Timedelta(minutes=10)
    grid = pd.date_range(start, stop, freq=pd.Timedelta(seconds=step_s))
    triggers: list[tuple[pd.Timestamp, str]] = []
    for t in grid:
        te = hsd_electron_temperature(t)
        for kind in ("FF", "FR"):
            up, down = shock_windows(data, t, kind)
            if up is None or down is None:
                continue
            p = analyze_shock(up, down, **te)
            if p.meets_criteria and p.shock_type == kind:
                triggers.append((t, kind))

    out: list[RealDetection] = []
    cluster: list[tuple[pd.Timestamp, str]] = []

    def flush():
        if not cluster:
            return
        kinds = [k for _, k in cluster]
        kind = max(set(kinds), key=kinds.count)
        t = _steepest_change(data.mag, cluster[0][0], cluster[-1][0]) or cluster[len(cluster) // 2][0]
        up, down = shock_windows(data, t, kind)
        if up is None or down is None:
            return
        out.append(RealDetection(time=t, params=analyze_shock(up, down, **hsd_electron_temperature(t)), n_triggers=len(cluster)))

    for trig in triggers:
        if cluster and (trig[0] - cluster[-1][0]).total_seconds() > merge_s:
            flush()
            cluster = []
        cluster.append(trig)
    flush()
    return out
