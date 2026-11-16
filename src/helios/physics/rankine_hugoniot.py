"""Rankine-Hugoniot jump conditions for an ideal-MHD fast shock.

Given the upstream state and the shock geometry, solve for the downstream
state. This is used to manufacture synthetic shocks whose parameters are known
exactly, so every analysis routine can be tested against a ground truth before
it ever sees spacecraft data.

Conservation laws across the shock, in the shock rest frame, with n the unit
normal, subscript n the normal component and t the tangential part
(brackets denote downstream minus upstream):

    [rho v_n] = 0                                         mass
    [B_n] = 0                                             div B = 0
    [rho v_n^2 + p + B_t^2 / 2 mu0] = 0                   normal momentum
    [rho v_n v_t - B_n B_t / mu0] = 0                     tangential momentum
    [v_n B_t - B_n v_t] = 0                               tangential electric field
    [v_n (rho v^2 / 2 + gamma p / (gamma-1) + B^2/mu0) - B_n (v . B) / mu0] = 0   energy

Method: for a trial compression ratio r = rho2/rho1 the mass, tangential
momentum and induction equations give the downstream field in closed form,

    B_t2 = B_t1 * r (M_An^2 - 1) / (M_An^2 - r),      M_An = |v_n1| / v_An,

and normal momentum then gives p2. The energy equation is left as one scalar
residual in r, solved by bracketing. r = 1 (no shock) is always a root, so the
residual is divided by (r - 1) to remove it.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from ..constants import CC, GAMMA, KB, KMS, MP, MU0, NT
from . import mhd


@dataclass(frozen=True)
class PlasmaState:
    """Bulk plasma state in practical units (nT, km/s, cm^-3, K)."""

    b: np.ndarray  # magnetic field vector, nT
    v: np.ndarray  # bulk velocity vector, km/s
    n: float  # proton number density, cm^-3
    tp: float  # proton temperature, K
    # Mean field magnitude <|B|> over the averaging window, when the state is a window average.
    # It exceeds |<B>| whenever the field direction fluctuates, so the two must not be confused.
    b_mag: float | None = None

    @property
    def field_strength(self) -> float:
        """<|B|> if known, otherwise the magnitude of the mean vector."""
        return float(np.linalg.norm(self.b)) if self.b_mag is None else float(self.b_mag)

    def __post_init__(self):
        object.__setattr__(self, "b", np.asarray(self.b, dtype=float))
        object.__setattr__(self, "v", np.asarray(self.v, dtype=float))


@dataclass(frozen=True)
class ShockSolution:
    upstream: PlasmaState
    downstream: PlasmaState
    normal: np.ndarray  # unit normal, pointing into the upstream medium
    shock_speed: float  # speed of the shock along `normal` in the observer frame, km/s
    compression: float  # rho2 / rho1
    theta_bn_deg: float  # acute angle between upstream B and the normal
    mach_alfven: float  # upstream normal flow speed / Alfven speed
    mach_fast: float  # upstream normal flow speed / fast speed along the normal
    te_over_tp: float


def _si(state: PlasmaState, te_over_tp: float, te_k: float | None = None):
    rho = state.n * CC * MP
    te = state.tp * te_over_tp if te_k is None else te_k
    p = state.n * CC * KB * (state.tp + te)
    return rho, state.v * KMS, state.b * NT, p


def conservation_residuals(up: PlasmaState, down: PlasmaState, normal, shock_speed: float,
                           te_over_tp: float = 1.0, gamma: float = GAMMA, te_k: float | None = None) -> dict[str, float]:
    """Relative violation of each conservation law between two states.

    Velocities are transformed to the shock rest frame using `shock_speed`
    along `normal`. Each residual is normalised by the size of the upstream
    flux, so values near machine precision mean the law is satisfied.
    """
    n = np.asarray(normal, dtype=float)
    n = n / np.linalg.norm(n)

    def fluxes(state):
        rho, v, b, p = _si(state, te_over_tp, te_k)
        u = v - shock_speed * KMS * n
        un, bn = u @ n, b @ n
        ut, bt = u - un * n, b - bn * n
        return {
            "mass": np.array([rho * un]),
            "normal_field": np.array([bn]),
            "normal_momentum": np.array([rho * un**2 + p + (bt @ bt) / (2 * MU0)]),
            "tangential_momentum": rho * un * ut - bn * bt / MU0,
            "tangential_electric_field": un * bt - bn * ut,
            "energy": np.array([un * (0.5 * rho * (u @ u) + gamma * p / (gamma - 1) + (b @ b) / MU0) - bn * (u @ b) / MU0]),
        }

    f1, f2 = fluxes(up), fluxes(down)
    rho1, v1, b1, p1 = _si(up, te_over_tp, te_k)
    un1 = (v1 - shock_speed * KMS * n) @ n
    scales = {
        "mass": abs(rho1 * un1),
        "normal_field": np.linalg.norm(b1),
        "normal_momentum": rho1 * un1**2 + p1 + (b1 @ b1) / (2 * MU0),
        "tangential_momentum": rho1 * un1**2 + (b1 @ b1) / MU0,
        "tangential_electric_field": abs(un1) * np.linalg.norm(b1),
        "energy": abs(un1) * (0.5 * rho1 * un1**2 + gamma * p1 / (gamma - 1) + (b1 @ b1) / MU0),
    }
    return {k: float(np.linalg.norm(f2[k] - f1[k]) / scales[k]) for k in f1}


def solve_fast_shock(upstream: PlasmaState, normal, mach_fast: float,
                     te_over_tp: float = 1.0, gamma: float = GAMMA) -> ShockSolution:
    """Downstream state of a fast shock with the given fast Mach number.

    `normal` points into the upstream medium (the direction the shock travels
    relative to the upstream plasma). `mach_fast` is the upstream normal flow
    speed in the shock frame divided by the fast speed along the normal, and
    must exceed 1 for a shock to exist. The electron-to-proton temperature
    ratio is held fixed across the shock.
    """
    if mach_fast <= 1:
        raise ValueError("A fast shock requires a fast Mach number greater than 1.")
    n = np.asarray(normal, dtype=float)
    n = n / np.linalg.norm(n)
    rho1, v1, b1, p1 = _si(upstream, te_over_tp)

    bn = b1 @ n
    bt1 = b1 - bn * n
    bt1_mag = np.linalg.norm(bt1)
    b_mag = np.linalg.norm(b1)
    theta = np.arccos(np.clip(abs(bn) / b_mag, 0, 1))

    va = b_mag / np.sqrt(MU0 * rho1)
    cs = np.sqrt(gamma * p1 / rho1)
    vf = float(mhd.fast_speed(va, cs, theta))
    vn1 = mach_fast * vf  # upstream normal speed in the shock frame (positive)
    van2 = bn**2 / (MU0 * rho1)
    man2 = vn1**2 / van2 if van2 > 0 else np.inf  # normal Alfven Mach number squared

    def downstream(r):
        """B_t2 magnitude and p2 for compression ratio r."""
        bt2 = bt1_mag * r * (man2 - 1) / (man2 - r) if np.isfinite(man2) else bt1_mag * r
        p2 = p1 + rho1 * vn1**2 * (1 - 1 / r) + (bt1_mag**2 - bt2**2) / (2 * MU0)
        return bt2, p2

    def energy_residual(r):
        bt2, p2 = downstream(r)
        # tangential velocity change in the shock frame, along the B_t direction
        dvt = bn * (bt2 - bt1_mag) / (MU0 * rho1 * (-vn1)) if bn != 0 else 0.0
        # choose the tangential frame where v_t1 = 0 (allowed: the laws are Galilean invariant along the shock)
        vt2 = dvt
        vn2 = vn1 / r
        # energy flux: u_n (rho u^2/2 + gamma p/(gamma-1) + B^2/mu0) - B_n (u.B)/mu0.
        # The B_n^2 parts of the two field terms cancel, leaving B_t^2 and u_t.B_t.
        e1 = -vn1 * (0.5 * rho1 * vn1**2 + gamma * p1 / (gamma - 1) + bt1_mag**2 / MU0)
        e2 = -vn2 * (0.5 * rho1 * r * (vn2**2 + vt2**2) + gamma * p2 / (gamma - 1) + bt2**2 / MU0) \
            - bn * (vt2 * bt2) / MU0
        scale = vn1 * (0.5 * rho1 * vn1**2 + gamma * p1 / (gamma - 1) + b_mag**2 / MU0)
        return (e2 - e1) / scale / (r - 1)

    r_max = (gamma + 1) / (gamma - 1)
    upper = min(r_max, man2) if np.isfinite(man2) else r_max
    grid = np.linspace(1 + 1e-6, upper - 1e-9 * upper, 4000)
    vals = np.array([energy_residual(r) for r in grid])
    crossings = np.where(np.sign(vals[:-1]) * np.sign(vals[1:]) < 0)[0]
    root = None
    for i in crossings:  # take the first physical root: positive pressure, compressive
        r = brentq(energy_residual, grid[i], grid[i + 1], xtol=1e-14, rtol=1e-14)
        if downstream(r)[1] > 0:
            root = r
            break
    if root is None:
        raise RuntimeError("No fast-shock solution found for this upstream state and Mach number.")

    r = root
    bt2, p2 = downstream(r)
    t_hat = bt1 / bt1_mag if bt1_mag > 0 else np.zeros(3)
    b2 = bn * n + bt2 * t_hat
    dvt = bn * (bt2 - bt1_mag) / (MU0 * rho1 * (-vn1)) if bn != 0 else 0.0
    # Observer frame: the shock moves along n at V_sh = v1.n + vn1.
    shock_speed = (v1 @ n + vn1) / KMS
    v2 = v1 + (vn1 - vn1 / r) * n + dvt * t_hat
    n2 = upstream.n * r
    tp2 = p2 / (n2 * CC * KB * (1 + te_over_tp))
    down = PlasmaState(b=b2 / NT, v=v2 / KMS, n=n2, tp=tp2)
    return ShockSolution(
        upstream=upstream, downstream=down, normal=n, shock_speed=float(shock_speed), compression=float(r),
        theta_bn_deg=float(np.degrees(theta)), mach_alfven=float(vn1 / va), mach_fast=float(mach_fast),
        te_over_tp=te_over_tp,
    )
