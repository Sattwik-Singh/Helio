"""Synthetic spacecraft time series containing a shock with known parameters.

The shock is an exact Rankine-Hugoniot solution (physics/rankine_hugoniot.py)
swept past a stationary observer, with a finite ramp and additive noise. Real
solar wind is far messier; these series exist to verify that the analysis
code recovers known answers, not to stand in for data.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from ..physics.rankine_hugoniot import PlasmaState, ShockSolution, solve_fast_shock

COLUMNS = ["bx", "by", "bz", "vx", "vy", "vz", "n", "tp"]


@dataclass(frozen=True)
class SyntheticEvent:
    data: pd.DataFrame  # indexed by time in seconds; columns bx by bz (nT), vx vy vz (km/s), n (cm^-3), tp (K)
    shock_time: float  # seconds
    shock_type: str  # "FF" or "FR"
    truth: ShockSolution


def _vector(state: PlasmaState) -> np.ndarray:
    return np.concatenate([state.b, state.v, [state.n, state.tp]])


def _red_noise(rng, n, corr):
    """Unit-variance AR(1) noise: a crude stand-in for correlated solar wind fluctuations."""
    x = np.empty(n)
    x[0] = rng.standard_normal()
    eps = rng.standard_normal(n) * np.sqrt(1 - corr**2)
    for i in range(1, n):
        x[i] = corr * x[i - 1] + eps[i]
    return x


def make_event(upstream: PlasmaState, normal, mach_fast: float, *, duration_s: float = 3600.0,
               cadence_s: float = 60.0, ramp_s: float = 5.0, noise: float = 0.0,
               noise_correlation: float = 0.7, seed: int = 0, te_over_tp: float = 1.0) -> SyntheticEvent:
    """One shock crossing at the middle of the interval.

    `noise` is the fluctuation amplitude as a fraction of each quantity's
    local mean magnitude (of |B| and of the thermal speed scale for velocity).
    The shock type follows from the geometry: if the shock moves along its
    normal in the observer frame the upstream side arrives first (FF);
    otherwise it is carried backwards over the observer (FR).
    """
    sol = solve_fast_shock(upstream, normal, mach_fast, te_over_tp=te_over_tp)
    rng = np.random.default_rng(seed)
    t = np.arange(0.0, duration_s, cadence_s)
    t0 = duration_s / 2
    forward = sol.shock_speed > 0
    # s goes 0 -> 1 from the first-seen side to the second
    s = 0.5 * (1 + np.tanh((t - t0) / max(ramp_s, 1e-9)))
    first, second = (sol.upstream, sol.downstream) if forward else (sol.downstream, sol.upstream)
    a, b = _vector(first), _vector(second)
    values = a[None, :] + s[:, None] * (b - a)[None, :]

    if noise > 0:
        bmag = np.linalg.norm(values[:, 0:3], axis=1)
        vscale = 0.1 * np.linalg.norm(values[:, 3:6], axis=1)
        scales = np.column_stack([bmag, bmag, bmag, vscale, vscale, vscale, values[:, 6], values[:, 7]])
        fluct = np.column_stack([_red_noise(rng, len(t), noise_correlation) for _ in COLUMNS])
        values = values + noise * scales * fluct
        values[:, 6] = np.maximum(values[:, 6], 1e-3)
        values[:, 7] = np.maximum(values[:, 7], 1.0)

    df = pd.DataFrame(values, index=pd.Index(t, name="t"), columns=COLUMNS)
    return SyntheticEvent(data=df, shock_time=t0, shock_type="FF" if forward else "FR", truth=sol)


def make_quiet(upstream: PlasmaState, *, duration_s: float = 3600.0, cadence_s: float = 60.0,
               noise: float = 0.05, noise_correlation: float = 0.7, seed: int = 0) -> pd.DataFrame:
    """A shock-free interval with the same noise model, for false-alarm tests."""
    rng = np.random.default_rng(seed)
    t = np.arange(0.0, duration_s, cadence_s)
    base = np.tile(_vector(upstream), (len(t), 1))
    bmag = np.linalg.norm(upstream.b)
    vscale = 0.1 * np.linalg.norm(upstream.v)
    scales = np.array([bmag, bmag, bmag, vscale, vscale, vscale, upstream.n, upstream.tp])
    fluct = np.column_stack([_red_noise(rng, len(t), noise_correlation) for _ in COLUMNS])
    values = base + noise * scales[None, :] * fluct
    values[:, 6] = np.maximum(values[:, 6], 1e-3)
    values[:, 7] = np.maximum(values[:, 7], 1.0)
    return pd.DataFrame(values, index=pd.Index(t, name="t"), columns=COLUMNS)
