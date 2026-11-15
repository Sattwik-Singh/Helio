# Helios

Physics-audited detection and analysis of interplanetary shocks in solar wind data.

**Status: early.** The physics core, a baseline detector and the evaluation
code are written and tested on synthetic shocks. Nothing has been run on
spacecraft data yet, and no detection performance is claimed. See
[PLAN.md](PLAN.md) for the roadmap and what is verified.

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
synthetic noise model only; whether it holds for real solar wind is a Phase 4
question.

## Install and test

```bash
pip install -e ".[dev]"
pytest
```

## Data and acknowledgements

Solar wind data will come from NASA's CDAWeb through its HAPI server. The
reference catalogue is the Heliospheric Shock Database, generated and
maintained at the University of Helsinki (Kilpua et al. 2015).

## License

MIT
