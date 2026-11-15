"""First contact with the data: check dataset names, list their real parameters,
and download a two-hour window around a few catalogued shocks.

Run where cdaweb.gsfc.nasa.gov is reachable (your machine, or the "Probe data"
GitHub Actions workflow):

    python scripts/probe_hapi.py

Writes data/probe/datasets.json (parameter lists as the server reports them)
and one CSV per event and instrument. Nothing here is analysis; it exists to
replace assumptions in helios/data/datasets.py with facts.
"""

from __future__ import annotations

import json
import sys
import traceback
from pathlib import Path

import pandas as pd

from helios.data import hapi
from helios.data.datasets import PROBE_EVENTS, SPACECRAFT

OUT = Path(__file__).resolve().parents[1] / "data" / "probe"


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    report: dict = {"server": hapi.CDAWEB, "datasets": {}, "events": [], "errors": []}
    metas: dict[str, dict] = {}

    for craft, products in SPACECRAFT.items():
        for kind, spec in products.items():
            ds = spec["dataset"]
            try:
                meta = hapi.info(ds)
                metas[ds] = meta
                params = [{"name": p["name"], "units": p.get("units"), "size": p.get("size"), "fill": p.get("fill")} for p in meta["parameters"]]
                names = {p["name"] for p in params}
                report["datasets"][ds] = {
                    "spacecraft": craft, "kind": kind, "start": meta.get("startDate"), "stop": meta.get("stopDate"),
                    "cadence": meta.get("cadence"), "parameters": params,
                    "candidates_found": [c for c in spec["candidates"] if c in names],
                    "candidates_missing": [c for c in spec["candidates"] if c not in names],
                }
                print(f"ok    {ds}: {len(params)} parameters; missing candidates: {report['datasets'][ds]['candidates_missing']}")
            except Exception as err:  # noqa: BLE001 - report every failure, keep going
                report["errors"].append({"dataset": ds, "error": repr(err)})
                print(f"FAIL  {ds}: {err!r}")

    for ev in PROBE_EVENTS:
        t = pd.Timestamp(ev["time"])
        start = (t - pd.Timedelta(hours=1)).strftime("%Y-%m-%dT%H:%M:%SZ")
        stop = (t + pd.Timedelta(hours=1)).strftime("%Y-%m-%dT%H:%M:%SZ")
        for kind, spec in SPACECRAFT[ev["spacecraft"]].items():
            ds = spec["dataset"]
            entry = {"event": ev, "dataset": ds}
            try:
                found = report["datasets"][ds]["candidates_found"]
                if not found:
                    raise RuntimeError("none of the candidate parameters exist in this dataset")
                df = hapi.fetch(ds, found, start, stop, meta=metas[ds])
                name = f"{ev['spacecraft']}_{kind}_{t.strftime('%Y%m%dT%H%M%S')}.csv"
                df.to_csv(OUT / name)
                entry.update({"file": name, "rows": int(len(df)), "columns": list(df.columns),
                              "nan_fraction": float(df.isna().mean().mean()) if len(df) else None})
                print(f"ok    {name}: {len(df)} rows")
            except Exception as err:  # noqa: BLE001
                entry["error"] = repr(err)
                report["errors"].append(entry)
                print(f"FAIL  {ds} @ {ev['time']}: {err!r}")
                traceback.print_exc(limit=1)
            report["events"].append(entry)

    (OUT / "datasets.json").write_text(json.dumps(report, indent=2, default=str) + "\n")
    print(f"\nWrote {OUT / 'datasets.json'} with {len(report['errors'])} error(s).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
