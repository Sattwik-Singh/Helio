"""Reference values for the four probe shocks, from the Heliospheric Shock Database.

Transcribed from the public listing at https://ipshocks.helsinki.fi/database
(generated and maintained at the University of Helsinki). Column meanings
were inferred from the listing's order and checked for self-consistency
(e.g. beta against the two speeds); they will be replaced by the catalogue's
own export, with its column headers, as soon as that is available.
"""

REFERENCE = [
    {"time": "2018-04-19T23:46:34Z", "spacecraft": "wind", "type": "FF", "b_ratio": 2.01, "dv": 42, "n_ratio": 3.00,
     "t_ratio": 3.18, "v_alfven": 17, "v_ms": 50, "beta": 9.1, "normal": [-0.77, -0.61, -0.17], "theta_bn": 63,
     "v_shock": 310, "mach_ms": 1.6},
    {"time": "2018-03-09T17:08:21Z", "spacecraft": "wind", "type": "FF", "b_ratio": 2.01, "dv": 28, "n_ratio": 2.51,
     "t_ratio": 1.71, "v_alfven": 27, "v_ms": 56, "beta": 4.0, "normal": [-0.64, 0.04, -0.77], "theta_bn": 52,
     "v_shock": 306, "mach_ms": 1.4},
    {"time": "2018-03-09T17:22:40Z", "spacecraft": "ace", "type": "FF", "b_ratio": 2.05, "dv": 30, "n_ratio": 2.51,
     "t_ratio": 2.12, "v_alfven": 25, "v_ms": 55, "beta": 4.4, "normal": [-0.53, -0.23, -0.82], "theta_bn": 89,
     "v_shock": 285, "mach_ms": 1.5},
    {"time": "2018-01-09T00:32:12Z", "spacecraft": "wind", "type": "FR", "b_ratio": 1.54, "dv": 69, "n_ratio": 2.01,
     "t_ratio": 1.95, "v_alfven": 66, "v_ms": 91, "beta": 1.1, "normal": [0.99, 0.06, -0.16], "theta_bn": 46,
     "v_shock": 436, "mach_ms": 1.4},
]
