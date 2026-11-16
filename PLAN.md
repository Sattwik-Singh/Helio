# Helios — project plan and checklist

**Goal.** An open, reproducible pipeline that finds interplanetary shocks in
solar wind data, derives their physical parameters, checks each one against
physics, and is scored against an expert catalogue — with an interface that
shows the evidence for every event.

**The claim the project will be judged on.** Precision and recall against the
Heliospheric Shock Database (University of Helsinki, ipshocks.helsinki.fi),
plus how closely the derived shock parameters agree with that catalogue's.
Everything else supports that number.

**What is and is not new.** Automatic shock detection exists: the Helsinki
database itself uses a machine-learning candidate search (IPSVM) followed by
expert review. Helios does not claim to be the first detector. Its intended
contributions are (1) a fully open pipeline with a published benchmark,
(2) a physics audit that scores every candidate against the Rankine–Hugoniot
relations, (3) a cross-spacecraft test: predict when a shock seen at one
spacecraft should reach another and compare with what was observed, and
(4) an agent that decides what further evidence to fetch for ambiguous
candidates and writes an auditable justification.

Status: `[x]` done and tested, `[~]` written but unverified, `[ ]` to do.

## Phase 0 — Foundations (no data needed)

- [x] Package scaffold, tests, continuous integration
- [x] MHD characteristic speeds (Alfvén, sound, fast, plasma beta)
- [x] Rankine–Hugoniot solver: exact downstream state of an oblique fast shock
      (all six conservation laws hold to 1e-10 on 300 random shocks; reduces
      to gas dynamics for parallel shocks; field is frozen-in for perpendicular)
- [x] Shock normal methods: magnetic coplanarity, velocity coplanarity, mixed
      modes MX1–MX3; exact on exact shocks, always oriented upstream
- [x] Shock analysis: type, jump ratios, normal, θ_Bn, shock speed from mass
      flux, Mach numbers — following the Helsinki conventions and criteria
- [x] Synthetic event generator (forward and reverse shocks, ramp, correlated noise)
- [x] Baseline rule-based detector using the Helsinki criteria and windows
- [x] Catalogue matching: precision, recall, F1, timing error
- [x] First measurement: how each normal method degrades with noise
      (`scripts/normal_noise_study.py`)

## Phase 1 — First contact with data

- [x] HAPI client, corrected against the live server (parameters must be sent
      in dataset order; CDAWeb uses the HAPI 2 request names)
- [x] Dataset table for Wind, ACE, DSCOVR, OMNI verified by the probe
      (Wind's 3-second field is `WI_H0_MFI@1`; DSCOVR plasma covers only 2016-06 to 2019-06)
- [x] Data probe script and GitHub Actions workflow; four shock windows downloaded
- [x] Loader: field and plasma kept at native cadence, averaged per instrument
- [x] **Go/no-go check passed on 4 shocks** (3 Wind, 1 ACE; `scripts/compare_probe_events.py`):
      every catalogue parameter reproduced — median difference 0.1–0.7% for
      jump ratios, speeds, beta, θ_Bn and shock speed, 2.1% for the
      magnetosonic Mach number, 0.2° for the shock normal. Getting there
      required two conventions from the catalogue's documentation: its
      electron temperature (146 277 K at 1 au from August 2015, 2 T_p before)
      and field strengths as the mean of |B|, not the magnitude of the mean vector.
- [x] Baseline detector run blind over the four 2-hour windows
      (`scripts/scan_probe_events.py`): one detection in each, 1–8 s from the
      catalogued time, no false alarms. Eight hours of data centred on known
      shocks — encouraging, not a benchmark.
- [x] Helsinki catalogue export obtained (2,328 shocks; kept out of the repository, see `helios/data/catalog.py`)
- [ ] Extend the go/no-go check from 4 shocks to 20+, including pre-2015 events

## Phase 2 — Baseline benchmark

- [x] Vectorised scanner (`detect/fast.py`): a year in about 2 seconds, identical to the reference scanner on the probe events
- [x] Range downloader and `Fetch data` workflow; Wind 2022 and 2023 cached on the `data` branch (15-second field bins, native plasma)
- [x] Baseline benchmark against the catalogue (`scripts/benchmark.py`, 120 s matching tolerance):

  | Year | Catalogue shocks | Detections | Matched | False positives | Precision | Recall |
  |---|---|---|---|---|---|---|
  | 2023 | 36 | 64 | 35 | 29 | 0.55 | 0.97 |
  | 2022 | 21 | 52 | 21 | 31 | 0.40 | 1.00 |

  The one miss (2023-05-12) was found 131 s from the catalogued time, next to a data gap.
  For matched events the parameters agree with the catalogue: normal 0.4–0.5°, shock speed 0.2%, M_ms 0.6–0.7% (medians).
- [x] First look at the false positives: about a third have sub-fast upstream flow, others are gradual; the rest look shock-like
- [ ] Use native 3-second field for timing (the cache is 15-second bins)
- [ ] Review the unmatched candidates that pass the audit: weak shocks missing from the catalogue, or something else?

## Phase 3 — Physics audit and learned detector

- [x] Physics audit, first version (`detect/audit.py`): M_fast > 1 (required by the jump conditions, no threshold)
      and sharpness ≥ 0.3 (threshold chosen on 2023). Precision 0.55 → 0.71 on 2023 and, **held out**, 0.40 → 0.53 on
      2022, with no catalogued shock lost in either year. The gain is smaller on the held-out year.
- [ ] Extend the audit: Rankine–Hugoniot residuals with uncertainties, agreement between normal methods
- [ ] Parameter uncertainties by resampling the upstream/downstream windows
- [ ] Learned detector (1D CNN or gradient-boosted trees on window features),
      trained with a time-based split: train on early years, test on later ones
- [ ] Ablation: baseline vs + audit vs + learned detector, same test years
- [ ] Sensitivity to the assumed electron temperature

## Phase 4 — Multiple spacecraft

- [ ] Spacecraft positions (the plasma datasets already carry them: `SC_pos_gse`, `SC_pos_GSE`)
- [ ] First test case is already in hand: the 2018-03-09 shock reached Wind at 17:08:21 and ACE at 17:22:40
- [ ] Associate the same shock at Wind, ACE and DSCOVR
- [ ] Arrival-time test: from one spacecraft's normal and speed, predict the
      arrival at the others; compare with observation; report the error distribution
- [ ] Which normal method predicts arrival times best on real shocks?

## Phase 5 — Agent

- [ ] Tool interface: fetch data, analyse a window, audit, query the catalogue
- [ ] Agent loop for ambiguous candidates: chooses which instrument or
      spacecraft to consult next, ends with accept / reject and stated reasons
- [ ] Evaluate the agent as a classifier on held-out candidates, against the
      audit score alone. Keep it only if it adds measurable value.
- [ ] Every decision logged with the evidence it used

## Phase 6 — Interface and release

- [ ] Event catalogue file with parameters, uncertainties and audit scores
- [ ] Web interface: event list, time-series evidence plots, 3D view of
      spacecraft positions and the shock plane sweeping past them
- [ ] README with the benchmark, methods, limitations and data acknowledgements
- [ ] Deployed demo

## Known risks

| Risk | Why it matters | How it is retired |
|---|---|---|
| Data access from the development sandbox | It cannot reach NASA servers | Probe and downloads run on GitHub Actions or locally |
| Parameter conventions differ from the catalogue | Electron temperature, window choice and averaging change Mach numbers | Phase 1 go/no-go on 20 known shocks |
| Baseline floods with false positives | Gradual compressions pass the criteria (shown on synthetic data) | Sharpness and audit score in Phase 3 |
| Catalogue is not complete ground truth | Real shocks missing from it count as false positives | Manually review a sample of "false positives" |
| Agent adds cost but not accuracy | An LLM in the loop must earn its place | Phase 5 evaluation against the audit score alone |
