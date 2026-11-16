"""Run the baseline detector blind over each two-hour probe window and compare
what it finds with the catalogued shock time.

    python scripts/scan_probe_events.py
"""

import pandas as pd

from compare_probe_events import load
from helios.data.reference import REFERENCE
from helios.detect.scan import scan
from helios.physics.normals import angle_between_deg

for ref in REFERENCE:
    found = scan(load(ref))
    t_ref = pd.Timestamp(ref["time"])
    print(f"\n{ref['time']}  {ref['spacecraft'].upper()}  catalogue: {ref['type']}  ->  {len(found)} detection(s) in the 2-hour window")
    for d in found:
        dt = (d.time - t_ref).total_seconds()
        print(f"  {d.time.strftime('%H:%M:%S')}  ({dt:+.0f} s from catalogue)  type {d.params.shock_type}  criteria met at refined time: {d.params.meets_criteria}"
              f"  B ratio {d.params.b_ratio:.2f}  N ratio {d.params.n_ratio:.2f}  normal within {angle_between_deg(d.params.normal, ref['normal']):.1f} deg"
              f"  triggers {d.n_triggers}")
