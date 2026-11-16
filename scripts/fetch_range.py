"""Download a date range of field and plasma data for one spacecraft and cache it compactly.

    python scripts/fetch_range.py wind 2023-01-01 2024-01-01 --out data/cache

Requests go to CDAWeb one day at a time. The magnetic field arrives at native
cadence (3 s for Wind: ten million rows a year) and is reduced on the spot to
15-second bins holding the mean vector, the mean of |B| and the sample count.
Both are kept because they differ whenever the field direction fluctuates,
and the catalogue's field strengths are means of |B|. Plasma moments are kept
at native cadence. Output: two .npz files with float32 columns and int64
epoch-nanosecond times.
"""

from __future__ import annotations

import argparse
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from helios.data import hapi
from helios.data.datasets import SPACECRAFT
from helios.data.loader import standardize_mag, standardize_plasma

BIN = "15s"


def fetch_day(dataset, params, day, meta, tries=4):
    start = day.strftime("%Y-%m-%dT00:00:00Z")
    stop = (day + pd.Timedelta(days=1)).strftime("%Y-%m-%dT00:00:00Z")
    for attempt in range(tries):
        try:
            return hapi.fetch(dataset, params, start, stop, meta=meta)
        except hapi.HapiError as err:
            if "1201" in str(err) or "no data" in str(err).lower():  # HAPI 1201: OK, no data for the interval
                return None
            last = err
        except Exception as err:  # noqa: BLE001 - network hiccups are retried
            last = err
        time.sleep(3 * (attempt + 1))
    print(f"  giving up on {dataset} {start[:10]}: {last!r}", flush=True)
    return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("spacecraft")
    ap.add_argument("start")
    ap.add_argument("stop")
    ap.add_argument("--out", default="data/cache")
    ap.add_argument("--workers", type=int, default=6)
    args = ap.parse_args()

    craft = args.spacecraft
    spec = SPACECRAFT[craft]
    days = pd.date_range(args.start, args.stop, freq="D", inclusive="left", tz="UTC")
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    tag = f"{craft}_{days[0].strftime('%Y%m%d')}_{(days[-1] + pd.Timedelta(days=1)).strftime('%Y%m%d')}"

    mag_meta = hapi.info(spec["mag"]["dataset"])
    pls_meta = hapi.info(spec["plasma"]["dataset"])
    vector = spec["mag"]["candidates"][0]

    def one_day(day):
        raw_m = fetch_day(spec["mag"]["dataset"], spec["mag"]["candidates"], day, mag_meta)
        raw_p = fetch_day(spec["plasma"]["dataset"], spec["plasma"]["candidates"], day, pls_meta)
        mag = pls = None
        if raw_m is not None and len(raw_m):
            m = standardize_mag(raw_m, vector)
            m["bmag"] = np.linalg.norm(m[["bx", "by", "bz"]].to_numpy(), axis=1)
            g = m.resample(BIN)
            mag = g.mean()
            mag["count"] = g["bx"].count()
            mag = mag[mag["count"] > 0]
        if raw_p is not None and len(raw_p):
            pls = standardize_plasma(raw_p, craft)
            pls = pls[pls.notna().any(axis=1)]
        return day, mag, pls

    t0 = time.time()
    mags, plss, missing = [], [], []
    with ThreadPoolExecutor(args.workers) as pool:
        for k, (day, mag, pls) in enumerate(pool.map(one_day, days), 1):
            if mag is None or pls is None:
                missing.append((day.strftime("%Y-%m-%d"), mag is None, pls is None))
            if mag is not None:
                mags.append(mag)
            if pls is not None:
                plss.append(pls)
            if k % 20 == 0 or k == len(days):
                print(f"  {k}/{len(days)} days, {time.time() - t0:.0f} s", flush=True)

    def save(frames, name, cols):
        df = pd.concat(frames).sort_index()
        df = df[~df.index.duplicated()]
        t = df.index.tz_convert("UTC").tz_localize(None).to_numpy("datetime64[ns]").astype("int64")
        np.savez_compressed(out / f"{tag}_{name}.npz", time_ns=t, columns=np.array(cols),
                            values=df[cols].to_numpy(dtype="float32"))
        print(f"saved {tag}_{name}.npz: {len(df)} rows", flush=True)

    save(mags, "mag15s", ["bx", "by", "bz", "bmag", "count"])
    save(plss, "plasma", ["vx", "vy", "vz", "n", "tp"])
    print(f"days with missing data (date, no field, no plasma): {missing}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
