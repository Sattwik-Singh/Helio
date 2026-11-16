# Helios

Physics-audited detection and analysis of interplanetary shocks in solar wind data.

**Status: first benchmark done.** Two years of Wind data have been scanned
and scored against the Heliospheric Shock Database. See [PLAN.md](PLAN.md)
for the roadmap and what is verified.

## Benchmark: Wind, 2022 and 2023

A rule-based detector applies the database's own selection criteria across
the whole year, then a physics audit removes candidates that cannot be fast
shocks. Detections are matched to the catalogue within 120 seconds.

| Year | Stage | Detections | Matched | False positives | Precision | Recall |
|---|---|---|---|---|---|---|
| 2023 (tuning) | criteria only | 64 | 35 of 36 | 29 | 0.55 | 0.97 |
| 2023 (tuning) | + physics audit | 49 | 35 of 36 | 14 | 0.71 | 0.97 |
| 2022 (held out) | criteria only | 52 | 21 of 21 | 31 | 0.40 | 1.00 |
| 2022 (held out) | + physics audit | 40 | 21 of 21 | 19 | 0.53 | 1.00 |

- **Recall is high:** 56 of 57 catalogued shocks are found, with a median timing
  error of 4–6 s. The one miss was detected 131 s from the catalogued time,
  beside a data gap.
- **The audit has two tests.** The upstream flow must be super-fast
  (M_fast > 1), which the jump conditions require and which needs no
  threshold; and the field change must be abrupt (sharpness ≥ 0.3), a
  threshold chosen on 2023. On the held-out year the audit still loses no
  catalogued shock, but improves precision less.
- **What the remaining false positives are is not yet known.** Several look
  like sharp, weak shocks in the data. They may be real events absent from the
  catalogue, or not shocks at all; that needs review, and until then
  precision is measured against the catalogue as it stands.
- **For matched events the derived parameters agree with the catalogue**
  (medians): shock normal within 0.5°, shock speed within 0.2%, magnetosonic
  Mach number within 0.7%.

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
