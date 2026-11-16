"""Read an export of the Heliospheric Shock Database (ipshocks.helsinki.fi).

The export is the database's own download: one row per shock with upstream
and downstream values, derived parameters, their error estimates and the
spacecraft position. The file is not redistributed with this project; fetch
it from the database and pass its path.

    "This work uses data from the Heliospheric Shock Database, generated and
    maintained at the University of Helsinki."
"""

from __future__ import annotations

import pandas as pd

RENAME = {
    "Spacecraft": "spacecraft", "Shock Type": "type",
    "Magnetic Field Ratio": "b_ratio", "Proton Density Ratio": "n_ratio", "Proton Temperature Ratio": "t_ratio",
    "Speed Jump (km/s)": "dv", "Sound Speed Upstream (km/s)": "v_sound", "Alfven Speed Upstream (km/s)": "v_alfven",
    "Magnetosonic Speed Upstream (km/s)": "v_ms", "Plasma Beta Upstream": "beta",
    "Shock Normal X": "nx", "Shock Normal Y": "ny", "Shock Normal Z": "nz",
    "Shock Normal Coordinate System": "normal_frame", "Shock Theta (degrees)": "theta_bn",
    "Shock Speed (km/s)": "v_shock", "Alfven Mach Number": "mach_alfven", "Magnetosonic Mach Number": "mach_ms",
    "Spacecraft X": "sc_x", "Spacecraft Y": "sc_y", "Spacecraft Z": "sc_z", "Spacecraft Coordinate System": "sc_frame",
    "Length of the Analysis Interval (minutes)": "interval_min",
}


def load_catalog(path, spacecraft: str | None = None) -> pd.DataFrame:
    """Catalogue rows sorted by time, with short column names and a UTC `time` column."""
    df = pd.read_csv(path, index_col=False)
    df.columns = [c.strip() for c in df.columns]
    df["time"] = pd.to_datetime(
        dict(year=df["Year"], month=df["Month (1-12)"], day=df["Day (1-31)"], hour=df["Hour (0-23)"],
             minute=df["Minute (0-59)"], second=df["Second (0-59)"]), utc=True)
    df = df.rename(columns=RENAME)
    keep = ["time"] + [c for c in RENAME.values() if c in df.columns]
    df = df[keep].sort_values("time").reset_index(drop=True)
    if spacecraft:
        df = df[df["spacecraft"].str.lower() == spacecraft.lower()].reset_index(drop=True)
    return df
