"""Uncertainties and consistency checks from resampling the analysis windows.

Every shock parameter is a function of four window means. Resampling the
samples inside each window (bootstrap) and recomputing gives a distribution
for each parameter, hence an uncertainty that reflects how much the plasma
fluctuated during the windows.

Caveat: samples are drawn independently, but solar wind fluctuations are
correlated in time, so these uncertainties are lower bounds. Field samples
(15 s bins, ~32 per window) are far more numerous than plasma samples (~5 per
window at Wind), so plasma noise dominates anything involving velocity.

Consistency checks. A real shock satisfies all the jump conditions at once,
so quantities the analysis did not use to build the normal are genuine tests:

  * tangential electric field: [u_n B_t - B_n u_t] = 0
  * tangential momentum:       [rho u_n u_t - B_n B_t / mu0] = 0
  * the three mixed-mode normals (MX1, MX2, MX3) should coincide

(Continuity of B_n is not a test here: the MX3 normal is perpendicular to
the field change by construction. Normal momentum and energy involve the
electron pressure, which is assumed, so they are reported but not used.)
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from ..data.loader import EventData
from .normals import angle_between_deg, shock_normal
from .rankine_hugoniot import PlasmaState, conservation_residuals
from .shock import analyze_shock, hsd_electron_temperature

MINUTE = pd.Timedelta(minutes=1)


def window_samples(data: EventData, t: pd.Timestamp, shock_type: str):
    """Raw samples (field array, plasma array) for the upstream and downstream windows."""
    if shock_type == "FF":
        up, down = (t - 9 * MINUTE, t - 1 * MINUTE), (t + 2 * MINUTE, t + 10 * MINUTE)
    else:
        up, down = (t + 1 * MINUTE, t + 9 * MINUTE), (t - 10 * MINUTE, t - 2 * MINUTE)
    out = []
    for a, b in (up, down):
        m = data.mag.loc[a:b]
        cols = ["bx", "by", "bz"] + (["bmag"] if "bmag" in m else [])
        mm = m[cols].to_numpy()
        if "bmag" not in m:
            mm = np.column_stack([mm, np.linalg.norm(mm, axis=1)])
        pp = data.plasma.loc[a:b][["vx", "vy", "vz", "n", "tp"]].dropna().to_numpy()
        out.append((mm, pp))
    return out


def _state(mm: np.ndarray, pp: np.ndarray) -> PlasmaState:
    m, p = mm.mean(axis=0), pp.mean(axis=0)
    return PlasmaState(b=m[0:3], v=p[0:3], n=float(p[3]), tp=float(p[4]), b_mag=float(m[3]))


@dataclass(frozen=True)
class Consistency:
    e_field: float  # tangential electric field residual (relative)
    momentum: float  # tangential momentum residual (relative)
    normal_spread_deg: float  # largest angle between the MX1, MX2 and MX3 normals
    normal_sigma_deg: float  # bootstrap scatter of the MX3 normal about its mean direction
    theta_sigma_deg: float
    v_shock_sigma: float
    mach_fast_sigma: float
    mach_fast_low: float  # 16th percentile of M_fast: is the flow super-fast within uncertainty?
    n_plasma: int  # plasma samples in the smaller window

    def as_dict(self):
        return self.__dict__.copy()


def _residuals(up: PlasmaState, down: PlasmaState, te: dict):
    p = analyze_shock(up, down, **te)
    n = np.array(p.normal)
    if not np.all(np.isfinite(n)) or not np.isfinite(p.shock_speed_kms):
        return p, np.nan, np.nan, np.nan
    r = conservation_residuals(up, down, n, p.shock_speed_kms, te_over_tp=te.get("te_over_tp", 1.0), te_k=te.get("te_k"))
    normals = [shock_normal(up.b, down.b, up.v, down.v, m) for m in ("MX1", "MX2", "MX3")]
    spread = max(angle_between_deg(a, b) for i, a in enumerate(normals) for b in normals[i + 1:]) if all(
        np.all(np.isfinite(x)) for x in normals) else np.nan
    return p, r["tangential_electric_field"], r["tangential_momentum"], spread


def consistency(data: EventData, t: pd.Timestamp, shock_type: str, n_boot: int = 200, seed: int = 0) -> Consistency | None:
    (mu, pu), (md, pd_) = window_samples(data, t, shock_type)
    if min(len(mu), len(md)) < 3 or min(len(pu), len(pd_)) < 3:
        return None
    te = hsd_electron_temperature(t)
    up, down = _state(mu, pu), _state(md, pd_)
    p0, e0, m0, spread0 = _residuals(up, down, te)
    n0 = np.array(p0.normal)
    rng = np.random.default_rng(seed)
    angles, thetas, speeds, machs = [], [], [], []
    for _ in range(n_boot):
        pick = lambda a: a[rng.integers(0, len(a), len(a))]  # noqa: E731
        p = analyze_shock(_state(pick(mu), pick(pu)), _state(pick(md), pick(pd_)), **te)
        if not np.all(np.isfinite(p.normal)):
            continue
        angles.append(angle_between_deg(p.normal, n0))
        thetas.append(p.theta_bn_deg)
        speeds.append(p.shock_speed_kms)
        machs.append(p.mach_fast)
    if len(angles) < n_boot // 2:
        return None
    return Consistency(
        e_field=float(e0), momentum=float(m0), normal_spread_deg=float(spread0),
        normal_sigma_deg=float(np.sqrt(np.mean(np.square(angles)))), theta_sigma_deg=float(np.std(thetas)),
        v_shock_sigma=float(np.std(speeds)), mach_fast_sigma=float(np.std(machs)),
        mach_fast_low=float(np.percentile(machs, 16)), n_plasma=int(min(len(pu), len(pd_))),
    )
