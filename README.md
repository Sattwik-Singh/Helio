# Helios

Physics-audited detection and analysis of interplanetary shocks in solar wind data.

**Status: early, but running on real data.** On four catalogued shocks the
pipeline reproduces the reference catalogue's parameters and finds each shock
blind. That is a check that the method is implemented correctly, not a
performance benchmark; the benchmark over a full year is the next milestone.
See [PLAN.md](PLAN.md) for the roadmap and what is verified.

## First results on spacecraft data

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
