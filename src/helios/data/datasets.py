"""Candidate CDAWeb datasets for each spacecraft.

PROVISIONAL. Dataset identifiers are the CDAWeb names for these products;
parameter names are from memory of the CDAWeb variable lists and have NOT yet
been checked against the server from this project (the development sandbox
cannot reach cdaweb.gsfc.nasa.gov). scripts/probe_hapi.py asks the server for
each dataset's real parameter list; nothing downstream should rely on this
table until that probe has run and this note has been removed.

Instruments and native cadences follow the Heliospheric Shock Database
method documentation (2024-09-06), Appendix A.1.
"""

SPACECRAFT = {
    "wind": {
        "mag": {"dataset": "WI_H0_MFI", "instrument": "MFI", "cadence_s": 3, "candidates": ["B3GSE", "BGSE"]},
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
