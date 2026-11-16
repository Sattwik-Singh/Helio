# Helios

Physics-audited detection and analysis of interplanetary shocks in solar wind data.

**Status: benchmarked on three years of Wind data** against the Heliospheric
Shock Database. See [PLAN.md](PLAN.md) for the roadmap and what is verified.

## Benchmark: Wind, 2012, 2022 and 2023

A rule-based detector applies the database's own selection criteria across
each year; a physics audit then removes candidates that fail further tests.
Detections are matched to the catalogue within 120 seconds. Thresholds were
chosen on 2023 only; 2022 and 2012 were run afterwards, unchanged.

| Year | Stage | Detections | Matched | False positives | Precision | Recall |
|---|---|---|---|---|---|---|
| 2023 (tuning) | criteria only | 64 | 35 of 36 | 29 | 0.55 | 0.97 |
| | audit v1 | 49 | 35 | 14 | 0.71 | 0.97 |
| | audit v2 | 46 | 35 | 11 | 0.76 | 0.97 |
| 2022 (held out) | criteria only | 52 | 21 of 21 | 31 | 0.40 | 1.00 |
| | audit v1 | 40 | 21 | 19 | 0.53 | 1.00 |
| | audit v2 | 34 | 20 | 14 | 0.59 | 0.95 |
| 2012 (held out) | criteria only | 76 | 42 of 45 | 34 | 0.55 | 0.93 |
| | audit v1 | 55 | 41 | 14 | 0.75 | 0.91 |
| | audit v2 | 51 | 40 | 11 | 0.78 | 0.89 |

- **Recall of the criteria alone:** 98 of 102 catalogued shocks across the
  three years, with a median timing error of 4–9 s. Two of the four misses
  fall beside data gaps.
- **Audit v1** requires super-fast upstream flow (M_fast > 1, which the jump
  conditions demand and which needs no threshold) and an abrupt field change
  (sharpness ≥ 0.3). **Audit v2** adds a shock normal that stays within 20°
  under resampling of the analysis windows.
- **The audit is a trade, not a free gain.** On the held-out years it removes
  39% to 68% of the false positives but also loses up to two catalogued
  shocks per year.
- **What the remaining false positives are is not yet known.** They are
  weaker on average than catalogued shocks (median field ratio 1.4 against
  1.9 in 2023), and several look like sharp, weak shocks in the data. They
  may be real events absent from the catalogue, or not shocks at all. Until
  that is reviewed, precision is measured against the catalogue as it stands.
- **For matched events the derived parameters agree with the catalogue**
  (medians, all three years): shock normal within 0.6°, shock speed within
  0.3%, magnetosonic Mach number within 0.7%.

### Other findings

- **A learned ranking** (logistic regression on six audit features, trained
  on 2023's 64 candidates: `scripts/rank_candidates.py`) separates catalogued
  from other candidates with AUC 0.93 on 2022 and 0.87 on 2012. At its 0.5
  threshold it reaches precision 0.81 and 0.87 at recall 0.81 and 0.73.
- **Uncertainties** from resampling the windows are available for every
  parameter. They are 3–6 times smaller than the catalogue's error estimates
  (they are standard errors under an independence assumption, so lower
  bounds) but rank the events the same way (Spearman 0.84–0.89).
- **A negative result:** residuals of the tangential electric field and
  tangential momentum jump conditions barely separate catalogued shocks from
  other candidates (AUC about 0.6). With about five plasma samples per
  window, their noise is as large as the signal.
- **A documentation discrepancy:** the catalogue's method document describes
  a different electron temperature assumption before August 2015, but its
  2012 sound speeds are reproduced to 0.0% with the fixed 146 277 K and miss
  by 16% with the documented alternative.

Reproduce with `python scripts/benchmark.py wind <cache dir> <catalogue export>`.
The catalogue export is downloaded from the database and is not redistributed here.

## Parameter check on four shocks

Four shocks from the Heliospheric Shock Database (three at Wind, one at ACE,
2018), analysed at the catalogued time with the catalogue's own windows and
conventions (`python scripts/compare_probe_events.py`):

| Parameter | Median difference from the catalogue |
|---|---|
| B, N, T downstream/upstream ratios | 0.2%, 0.3%, 0.1% |
| Speed jump | 0.5% |
| Alfvén and magnetosonic speeds, plasma beta | 0.3%, 0.4%, 0.7% |
| Shock angle θ_Bn | 0.1% |
| Shock speed | 0.1% |
| Magnetosonic Mach number | 2.1% |
| Shock normal direction | 0.2° |

Run blind over the two-hour window around each shock
(`python scripts/scan_probe_events.py`), the baseline detector reported
exactly one event per window, 1 to 8 seconds from the catalogued time.

On these four real shocks the choice of normal method mattered a great deal:
mixed mode MX3 agreed with the catalogue (which also uses it) to within 1.4°,
while magnetic coplanarity differed by 20° to 88°.

## What works today

- **Rankine–Hugoniot solver.** Exact downstream state of an oblique ideal-MHD
  fast shock, verified against all six conservation laws to 1e-10.
- **Shock analysis.** Normals by magnetic coplanarity, velocity coplanarity and
  mixed modes MX1–MX3; shock speed from mass flux; θ_Bn; Alfvén, magnetosonic
  and fast Mach numbers. Conventions follow the Heliospheric Shock Database.
- **Synthetic events.** Forward and reverse shocks with known parameters,
  a finite ramp and correlated noise.
- **Baseline detector** using the Heliospheric Shock Database criteria, and
  catalogue matching (precision, recall, timing error).

## A first measurement

Error of each shock-normal method on 400 synthetic shocks per noise level
(median / 90th percentile, degrees; `python scripts/normal_noise_study.py`):

| Noise | MC | VC | MX1 | MX2 | MX3 |
|---|---|---|---|---|---|
| 0% | 0.0 / 0.0 | 3.4 / 12.0 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 |
| 2% | 3.4 / 19.3 | 3.6 / 12.4 | 1.3 / 2.8 | 1.3 / 2.9 | 1.3 / 2.9 |
| 5% | 8.6 / 37.7 | 3.9 / 13.0 | 3.1 / 7.1 | 3.1 / 7.1 | 3.3 / 7.2 |
| 10% | 17.4 / 54.8 | 5.2 / 14.3 | 6.3 / 14.0 | 6.3 / 14.0 | 6.4 / 14.1 |

Magnetic coplanarity is exact without noise but degrades fastest. Velocity
coplanarity is biased even without noise yet barely worsens. The mixed modes
are the best choice until noise reaches about 10%. This holds for the
synthetic noise model only.

## Install and test

```bash
pip install -e ".[dev]"
pytest
```

## Data and acknowledgements

Solar wind data come from NASA's CDAWeb through its HAPI server (Wind MFI and
SWE, ACE MAG and SWEPAM). The
reference catalogue is the Heliospheric Shock Database, generated and
maintained at the University of Helsinki (Kilpua et al. 2015).

## License

MIT
