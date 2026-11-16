"""Go/no-go check: analyse the probe shocks at their catalogued times and
compare every derived parameter with the Heliospheric Shock Database.

    python scripts/compare_probe_events.py

Uses the catalogue's own electron temperature assumption, so that sound speed,
beta and Mach numbers are compared like for like.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from helios.data.datasets import SPACECRAFT
from helios.data.loader import EventData, shock_windows, standardize_mag, standardize_plasma
from helios.data.reference import REFERENCE
from helios.physics import mhd
from helios.physics.normals import METHODS, angle_between_deg, shock_normal
from helios.physics.shock import analyze_shock, hsd_electron_temperature

PROBE = Path(__file__).resolve().parents[1] / "data" / "probe"


def load(ref) -> EventData:
    stamp = pd.Timestamp(ref["time"]).strftime("%Y%m%dT%H%M%S")
    craft = ref["spacecraft"]
    read = lambda kind: pd.read_csv(PROBE / f"{craft}_{kind}_{stamp}.csv", index_col=0, parse_dates=True)  # noqa: E731
    vector = SPACECRAFT[craft]["mag"]["candidates"][0]
    return EventData(craft, standardize_mag(read("mag"), vector), standardize_plasma(read("plasma"), craft))


def main():
    rows = []
    for ref in REFERENCE:
        data = load(ref)
        t = pd.Timestamp(ref["time"])
        up, down = shock_windows(data, t, ref["type"])
        if up is None or down is None:
            print(f"{ref['time']} {ref['spacecraft']}: not enough data in the analysis windows")
            continue
        p = analyze_shock(up, down, **hsd_electron_temperature(t))
        vms = float(mhd.magnetosonic_speed(p.alfven_speed_kms, p.sound_speed_kms))
        print(f"\n{ref['time']}  {ref['spacecraft'].upper()}  catalogue type {ref['type']}  ->  Helios type {p.shock_type}, criteria met: {p.meets_criteria}")
        pairs = [
            ("B ratio", p.b_ratio, ref["b_ratio"]), ("N ratio", p.n_ratio, ref["n_ratio"]), ("T ratio", p.t_ratio, ref["t_ratio"]),
            ("|dV| km/s", abs(p.dv_kms), ref["dv"]), ("v_A km/s", p.alfven_speed_kms, ref["v_alfven"]),
            ("v_ms km/s", vms, ref["v_ms"]), ("beta", p.plasma_beta, ref["beta"]), ("theta_Bn deg", p.theta_bn_deg, ref["theta_bn"]),
            ("V_shock km/s", abs(p.shock_speed_kms), ref["v_shock"]), ("M_ms", p.mach_magnetosonic, ref["mach_ms"]),
        ]
        for name, ours, theirs in pairs:
            print(f"  {name:<13} Helios {ours:8.2f}   catalogue {theirs:8.2f}   diff {100 * (ours - theirs) / theirs:+6.1f}%")
            rows.append((name, ours, theirs))
        angles = {m: angle_between_deg(shock_normal(up.b, down.b, up.v, down.v, m), ref["normal"]) for m in METHODS}
        print(f"  normal        Helios {np.round(p.normal, 2)}   catalogue {ref['normal']}   angle between: {angles['MX3']:.1f} deg")
        print("  normal angle to catalogue by method: " + ", ".join(f"{m} {a:.1f}" for m, a in angles.items()))
        rows.append(("normal angle", angles["MX3"], 0.0))

    print("\nMedian absolute difference from the catalogue, over the events analysed:")
    for name in dict.fromkeys(r[0] for r in rows):
        sel = [r for r in rows if r[0] == name]
        if name == "normal angle":
            print(f"  {name:<13} {np.median([r[1] for r in sel]):.1f} deg")
        else:
            print(f"  {name:<13} {np.median([abs(r[1] - r[2]) / r[2] for r in sel]) * 100:.1f}%")


if __name__ == "__main__":
    main()
