"""Benchmark the baseline detector against the Heliospheric Shock Database.

    python scripts/benchmark.py wind data/cache data/catalog/ipshocks_export.csv [--tolerance 120]

Scans every cached interval for the spacecraft, then matches detections to
the catalogue's shocks for that spacecraft inside the scanned time span.
Catalogue shocks that fall in a data gap (no analysable windows) are counted
separately: the detector cannot be blamed for data that is not there, but
they are reported, not hidden.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from helios.data.cache import load_cached
from helios.data.catalog import load_catalog
from helios.data.loader import shock_windows
from helios.detect.audit import audit
from helios.detect.fast import scan_fast
from helios.evaluation.matching import match_events
from helios.physics.normals import angle_between_deg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("spacecraft")
    ap.add_argument("cache")
    ap.add_argument("catalog")
    ap.add_argument("--tolerance", type=float, default=120.0, help="matching tolerance in seconds")
    ap.add_argument("--out", default="results")
    args = ap.parse_args()

    data = load_cached(args.spacecraft, args.cache)
    t_start = max(data.mag.index[0], data.plasma.index[0])
    t_stop = min(data.mag.index[-1], data.plasma.index[-1])
    days = (t_stop - t_start).total_seconds() / 86400

    t0 = time.time()
    detections = scan_fast(data)
    scan_s = time.time() - t0

    cat = load_catalog(args.catalog, args.spacecraft)
    cat = cat[(cat.time >= t_start) & (cat.time <= t_stop)].reset_index(drop=True)
    analysable = np.array([all(s is not None for s in shock_windows(data, t, k)) for t, k in zip(cat.time, cat.type)])

    audits = [audit(data, d) for d in detections]
    det_s = np.array([d.time.timestamp() for d in detections])
    cat_s = np.array([t.timestamp() for t in cat.time])
    m = match_events(det_s, cat_s, tolerance_s=args.tolerance)
    matched_cat = {j for _, j in m.pairs}
    missed_with_data = [j for j in m.misses if analysable[j]]
    missed_in_gap = [j for j in m.misses if not analysable[j]]
    n_analysable = int(analysable.sum())
    recall_analysable = len([j for j in matched_cat if analysable[j]]) / n_analysable if n_analysable else float("nan")

    print(f"{args.spacecraft.upper()}  {t_start:%Y-%m-%d} to {t_stop:%Y-%m-%d}  ({days:.0f} days, scanned in {scan_s:.1f} s)")
    print(f"catalogue shocks in span: {len(cat)}  (with analysable data: {n_analysable}; in data gaps: {len(cat) - n_analysable})")
    print(f"detections: {len(detections)}   matched within {args.tolerance:.0f} s: {len(m.pairs)}")
    print(f"precision {m.precision:.3f}   recall {m.recall:.3f} (all catalogue shocks)   recall {recall_analysable:.3f} (shocks with analysable data)")
    print(f"false positives: {len(m.false_positives)}  ({len(m.false_positives) / days * 365:.0f} per year)   median |time error| {m.median_abs_time_error_s:.0f} s")

    rows = []
    for i, d in enumerate(detections):
        j = next((j for a, j in m.pairs if a == i), None)
        row = {"time": d.time.isoformat(), "type": d.params.shock_type, "matched": j is not None, "triggers": d.n_triggers,
               "b_ratio": d.params.b_ratio, "n_ratio": d.params.n_ratio, "t_ratio": d.params.t_ratio, "dv": d.params.dv_kms,
               "theta_bn": d.params.theta_bn_deg, "v_shock": d.params.shock_speed_kms, "mach_ms": d.params.mach_magnetosonic,
               "mach_fast": d.params.mach_fast, "sharpness": audits[i].sharpness, "super_fast": audits[i].super_fast,
               "abrupt": audits[i].abrupt, "audit_passed": audits[i].passed, "audit_v2_passed": audits[i].passed_v2,
               "normal_sigma_deg": audits[i].normal_sigma_deg, "mach_fast_low": audits[i].mach_fast_low,
               "theta_sigma_deg": audits[i].theta_sigma_deg, "v_shock_sigma": audits[i].v_shock_sigma}
        if j is not None:
            c = cat.iloc[j]
            row.update({"catalog_time": c.time.isoformat(), "dt_s": d.time.timestamp() - c.time.timestamp(), "catalog_type": c.type,
                        "normal_angle_deg": angle_between_deg(d.params.normal, [c.nx, c.ny, c.nz]),
                        "cat_b_ratio": c.b_ratio, "cat_n_ratio": c.n_ratio, "cat_mach_ms": c.mach_ms, "cat_v_shock": c.v_shock, "cat_theta_bn": c.theta_bn})
        rows.append(row)
    table = pd.DataFrame(rows)

    # The same benchmark after the physics audit, and with each test alone.
    stages = {}
    for name, keep in (("baseline", np.ones(len(table), bool)), ("super-fast only", table.super_fast.to_numpy()),
                       ("abrupt only", table.abrupt.to_numpy()), ("audit v1", table.audit_passed.to_numpy()),
                       ("audit v2", table.audit_v2_passed.to_numpy())):
        mm = match_events(det_s[keep], cat_s, tolerance_s=args.tolerance)
        stages[name] = {"detections": int(keep.sum()), "matched": len(mm.pairs), "false_positives": len(mm.false_positives),
                        "precision": mm.precision, "recall_all": mm.recall, "f1": mm.f1}
    print("\nstage              detections  matched  false+  precision  recall    F1")
    for name, r in stages.items():
        print(f"{name:<18} {r['detections']:>10}  {r['matched']:>7}  {r['false_positives']:>6}  {r['precision']:>9.3f}  {r['recall_all']:>6.3f}  {r['f1']:>5.3f}")
    out = Path(args.out)
    out.mkdir(exist_ok=True)
    tag = f"{args.spacecraft}_{t_start:%Y%m%d}_{t_stop:%Y%m%d}"
    table.to_csv(out / f"benchmark_{tag}_detections.csv", index=False)
    misses = cat.iloc[m.misses][["time", "type", "b_ratio", "n_ratio", "t_ratio", "dv", "mach_ms", "interval_min"]].copy()
    misses["analysable"] = analysable[m.misses]
    misses.to_csv(out / f"benchmark_{tag}_misses.csv", index=False)
    summary = {"spacecraft": args.spacecraft, "start": t_start.isoformat(), "stop": t_stop.isoformat(), "days": days,
               "tolerance_s": args.tolerance, "catalogue_shocks": len(cat), "catalogue_analysable": n_analysable,
               "detections": len(detections), "matched": len(m.pairs), "false_positives": len(m.false_positives),
               "missed_with_data": len(missed_with_data), "missed_in_gap": len(missed_in_gap),
               "precision": m.precision, "recall_all": m.recall, "recall_analysable": recall_analysable,
               "median_abs_time_error_s": m.median_abs_time_error_s, "stages": stages}
    (out / f"benchmark_{tag}_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    if len(m.pairs):
        t = table[table.matched]
        print("\nmatched events vs catalogue (median absolute difference):")
        print(f"  normal {t.normal_angle_deg.median():.1f} deg   B ratio {(t.b_ratio / t.cat_b_ratio - 1).abs().median() * 100:.1f}%"
              f"   N ratio {(t.n_ratio / t.cat_n_ratio - 1).abs().median() * 100:.1f}%   shock speed {(t.v_shock.abs() / t.cat_v_shock - 1).abs().median() * 100:.1f}%"
              f"   M_ms {(t.mach_ms / t.cat_mach_ms - 1).abs().median() * 100:.1f}%   type agreement {(t.type == t.catalog_type).mean() * 100:.0f}%")
    if len(misses):
        print("\nmissed catalogue shocks:")
        print(misses.to_string(index=False))
    return summary, table, misses


if __name__ == "__main__":
    main()
