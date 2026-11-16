"""Read the compact .npz cache files written by scripts/fetch_range.py."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .loader import EventData


def load_npz(path) -> pd.DataFrame:
    z = np.load(path, allow_pickle=False)
    index = pd.DatetimeIndex(z["time_ns"].astype("datetime64[ns]"), tz="UTC", name="time")
    return pd.DataFrame(z["values"].astype("float64"), index=index, columns=[str(c) for c in z["columns"]])


def load_cached(spacecraft: str, directory, tag: str | None = None) -> EventData:
    """EventData from every matching cache file in `directory`, concatenated in time.

    The field frame keeps its `bmag` column (mean of |B| per 15-second bin),
    which the window averaging uses for field strengths.
    """
    directory = Path(directory)
    pattern = f"{spacecraft}_{tag or '*'}"
    mags = sorted(directory.glob(f"{pattern}_mag15s.npz"))
    plss = sorted(directory.glob(f"{pattern}_plasma.npz"))
    if not mags or not plss:
        raise FileNotFoundError(f"no cache files for {spacecraft!r} in {directory}")
    mag = pd.concat([load_npz(p) for p in mags]).sort_index()
    pls = pd.concat([load_npz(p) for p in plss]).sort_index()
    mag = mag[~mag.index.duplicated()].drop(columns=["count"], errors="ignore")
    pls = pls[~pls.index.duplicated()]
    return EventData(spacecraft, mag, pls)
