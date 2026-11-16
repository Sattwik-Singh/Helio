"""Vectorised baseline detector for long intervals.

Same rule as detect/scan.py (the Heliospheric Shock Database criteria on its
default windows), but every window mean for every candidate time is computed
at once from cumulative sums, so a year of data takes seconds instead of hours.

Window mean over [a, b] of samples v at sorted times t:
    i = first index with t >= a,  j = first index with t > b
    mean = (C[j] - C[i]) / (K[j] - K[i])
where C is the cumulative sum of v with gaps as zero and K the cumulative
count of valid samples. Both i and j come from one binary search per edge.

Candidates that pass are clustered, timed at the steepest change in |B| and
then re-analysed with the ordinary (scalar) code path, so the reported
parameters are identical to what analyze_shock gives for that time.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..data.loader import EventData, shock_windows
from ..physics.shock import (B_RATIO_MIN, DV_MIN_KMS, N_RATIO_MIN, T_RATIO_MIN, analyze_shock,
                             hsd_electron_temperature)
from .scan import RealDetection, _steepest_change

MIN_POINTS = 3
MINUTE_NS = 60_000_000_000


class _Cumulative:
    def __init__(self, df: pd.DataFrame):
        self.t = df.index.tz_convert("UTC").tz_localize(None).to_numpy("datetime64[ns]").astype("int64")
        v = df.to_numpy(dtype=float)
        ok = np.isfinite(v)
        zero = np.zeros((1, v.shape[1]))
        self.c = np.vstack([zero, np.cumsum(np.where(ok, v, 0.0), axis=0)])
        self.k = np.vstack([zero, np.cumsum(ok, axis=0)])
        self.cols = list(df.columns)

    def mean(self, a_ns: np.ndarray, b_ns: np.ndarray):
        """Per-column means and valid-sample counts over [a, b] for arrays of window edges."""
        i = np.searchsorted(self.t, a_ns, side="left")
        j = np.searchsorted(self.t, b_ns, side="right")
        count = self.k[j] - self.k[i]
        with np.errstate(invalid="ignore", divide="ignore"):
            return (self.c[j] - self.c[i]) / count, count


def candidate_mask(data: EventData, grid_ns: np.ndarray):
    """Boolean arrays (FF, FR): grid times at which the shock criteria hold."""
    mag_df = data.mag.copy()
    if "bmag" not in mag_df:
        mag_df["bmag"] = np.linalg.norm(mag_df[["bx", "by", "bz"]].to_numpy(), axis=1)
    mag, pls = _Cumulative(mag_df[["bx", "by", "bz", "bmag"]]), _Cumulative(data.plasma[["vx", "vy", "vz", "n", "tp"]])
    m = MINUTE_NS
    # windows relative to the candidate time: (start, stop) in minutes
    sides = {"before_near": (-9, -1), "before_far": (-10, -2), "after_near": (1, 9), "after_far": (2, 10)}
    stats = {}
    for name, (a, b) in sides.items():
        bm, bk = mag.mean(grid_ns + a * m, grid_ns + b * m)
        pm, pk = pls.mean(grid_ns + a * m, grid_ns + b * m)
        ok = (bk.min(axis=1) >= MIN_POINTS) & (pk.min(axis=1) >= MIN_POINTS)
        stats[name] = {"bmag": bm[:, 3], "speed": np.linalg.norm(pm[:, 0:3], axis=1), "n": pm[:, 3], "tp": pm[:, 4], "ok": ok}

    def passes(up, down, forward):
        with np.errstate(invalid="ignore", divide="ignore"):
            jumps = ((down["bmag"] / up["bmag"] >= B_RATIO_MIN) & (down["n"] / up["n"] >= N_RATIO_MIN)
                     & (down["tp"] / up["tp"] >= T_RATIO_MIN))
            dv = down["speed"] - up["speed"]
            speed = dv >= DV_MIN_KMS if forward else -dv >= DV_MIN_KMS
        return up["ok"] & down["ok"] & jumps & speed

    ff = passes(stats["before_near"], stats["after_far"], True)  # upstream first
    fr = passes(stats["after_near"], stats["before_far"], False)  # downstream first
    return ff, fr


def scan_fast(data: EventData, step_s: float = 30.0, merge_s: float = 900.0) -> list[RealDetection]:
    start = max(data.mag.index[0], data.plasma.index[0]) + pd.Timedelta(minutes=10)
    stop = min(data.mag.index[-1], data.plasma.index[-1]) - pd.Timedelta(minutes=10)
    grid = pd.date_range(start, stop, freq=pd.Timedelta(seconds=step_s))
    grid_ns = grid.tz_convert("UTC").tz_localize(None).to_numpy("datetime64[ns]").astype("int64")
    ff, fr = candidate_mask(data, grid_ns)
    hit = np.flatnonzero(ff | fr)
    if len(hit) == 0:
        return []

    mag_xyz = data.mag[["bx", "by", "bz"]]
    out: list[RealDetection] = []
    breaks = np.flatnonzero(np.diff(grid_ns[hit]) > merge_s * 1e9) + 1
    for cluster in np.split(hit, breaks):
        kind = "FF" if ff[cluster].sum() >= fr[cluster].sum() else "FR"
        t = _steepest_change(mag_xyz, grid[cluster[0]], grid[cluster[-1]]) or grid[cluster[len(cluster) // 2]]
        up, down = shock_windows(data, t, kind)
        if up is None or down is None:
            t = grid[cluster[len(cluster) // 2]]
            up, down = shock_windows(data, t, kind)
            if up is None or down is None:
                continue
        out.append(RealDetection(time=t, params=analyze_shock(up, down, **hsd_electron_temperature(t)), n_triggers=len(cluster)))
    return out
