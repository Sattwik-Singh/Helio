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

## Phase 1 — First contact with data  ← current

- [~] HAPI client (CSV parsing tested offline; never yet run against the server)
- [~] Dataset table for Wind, ACE, DSCOVR, OMNI (identifiers and parameter
      names provisional until the probe runs)
- [~] Data probe script and GitHub Actions workflow
- [ ] Run the probe; replace provisional parameter names with the server's
- [ ] Obtain the Helsinki catalogue export for Wind, ACE, DSCOVR
- [ ] Loader: merge magnetic field and plasma onto one time base, flag gaps,
      convert thermal speed to temperature where needed
- [ ] **Go/no-go check:** analyse 20 catalogued shocks by hand-picked time.
      Do Helios's B, N, T ratios, normal, θ_Bn, speed and Mach numbers agree
      with the catalogue? If not, find out why before building further.

## Phase 2 — Baseline benchmark

- [ ] Download one full year for one spacecraft (Wind), cached locally
- [ ] Run the baseline detector over the year; score against the catalogue
- [ ] Characterise failures: what do false positives look like? what is missed?
- [ ] Use native cadence (3 s field, 90 s plasma) instead of one-minute means
- [ ] Publish the baseline numbers. This is the bar every later stage must clear.

## Phase 3 — Physics audit and learned detector

- [ ] Physics audit score per candidate: sharpness, Rankine–Hugoniot residuals
      with uncertainties, agreement between normal methods, M_fast > 1
- [ ] Parameter uncertainties by resampling the upstream/downstream windows
- [ ] Learned detector (1D CNN or gradient-boosted trees on window features),
      trained with a time-based split: train on early years, test on later ones
- [ ] Ablation: baseline vs + audit vs + learned detector, same test years
- [ ] Sensitivity to the assumed electron temperature

## Phase 4 — Multiple spacecraft

- [ ] Spacecraft positions (SSCWeb or CDAWeb orbit datasets)
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
