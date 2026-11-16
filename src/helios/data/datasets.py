"""CDAWeb datasets for each spacecraft, as served by its HAPI server.

Dataset identifiers and parameter names were checked against the server by
scripts/probe_hapi.py (results in data/probe/datasets.json). Two things that
probe established:

  * WI_H0_MFI is not a dataset on the HAPI server; it is split into
    WI_H0_MFI@0 (1-minute), @1 (3-second) and @2 (1-hour).
  * Coverage is uneven: DSCOVR_H1_FC plasma spans only 2016-06 to 2019-06,
    and AC_H0_SWE ends in 2024-07.

Instruments and native cadences follow the Heliospheric Shock Database
method documentation (2024-09-06), Appendix A.1.
"""

SPACECRAFT = {
    "wind": {
        "mag": {"dataset": "WI_H0_MFI@1", "instrument": "MFI", "cadence_s": 3, "candidates": ["B3GSE"]},
        "plasma": {"dataset": "WI_K0_SWE", "instrument": "SWE", "cadence_s": 90, "candidates": ["V_GSE", "Np", "THERMAL_SPD"]},
    },
    "ace": {
        "mag": {"dataset": "AC_H0_MFI", "instrument": "MAG", "cadence_s": 16, "candidates": ["BGSEc"]},
        "plasma": {"dataset": "AC_H0_SWE", "instrument": "SWEPAM", "cadence_s": 64, "candidates": ["V_GSE", "Np", "Tpr"]},
    },
    "dscovr": {
        "mag": {"dataset": "DSCOVR_H0_MAG", "instrument": "PlasMag", "cadence_s": 1, "candidates": ["B1GSE"]},
        "plasma": {"dataset": "DSCOVR_H1_FC", "instrument": "PlasMag Faraday cup", "cadence_s": 60, "candidates": ["V_GSE", "Np", "THERMAL_TEMP"]},
    },
    "omni": {
        "merged": {"dataset": "OMNI_HRO2_1MIN", "instrument": "multi-spacecraft, shifted to the bow shock nose", "cadence_s": 60,
                   "candidates": ["BX_GSE", "BY_GSE", "BZ_GSE", "Vx", "Vy", "Vz", "proton_density", "T"]},
    },
}

# A few shocks listed in the Heliospheric Shock Database (ipshocks.helsinki.fi), used only to
# choose time windows for the first data probe.
PROBE_EVENTS = [
    {"time": "2018-04-19T23:46:34Z", "spacecraft": "wind", "type": "FF"},
    {"time": "2018-03-09T17:08:21Z", "spacecraft": "wind", "type": "FF"},
    {"time": "2018-03-09T17:22:40Z", "spacecraft": "ace", "type": "FF"},
    {"time": "2018-01-09T00:32:12Z", "spacecraft": "wind", "type": "FR"},
]
