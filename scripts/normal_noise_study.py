"""How well does each shock-normal method survive noise?

Synthetic shocks with exactly known normals are given correlated noise, then
analysed with the same 8-minute upstream and downstream windows the detector
uses. Prints the median and 90th-percentile angular error of each method.

    python scripts/normal_noise_study.py
"""

import numpy as np

from helios.detect.baseline import windows
from helios.physics.normals import METHODS, angle_between_deg, shock_normal
from helios.physics.rankine_hugoniot import PlasmaState
from helios.synth.generator import make_event


def random_case(rng):
    b = rng.normal(0, 4, 3)
    v = np.array([-rng.uniform(300, 600), rng.normal(0, 20), rng.normal(0, 20)])
    up = PlasmaState(b=b, v=v, n=rng.uniform(2, 15), tp=rng.uniform(2e4, 3e5))
    normal = np.array([-1.0, 0, 0]) + rng.normal(0, 0.4, 3)
    return up, normal / np.linalg.norm(normal), rng.uniform(1.3, 3.0)


def main(n_events: int = 400):
    print(f"{n_events} synthetic fast-forward shocks per noise level, 60 s cadence, 8-minute windows\n")
    print(f"{'noise':>6} | " + " | ".join(f"{m:>13}" for m in METHODS) + "   (median / 90th percentile error, degrees)")
    for noise in (0.0, 0.02, 0.05, 0.10):
        rng = np.random.default_rng(11)
        errs = {m: [] for m in METHODS}
        for k in range(n_events):
            up, normal, mach = random_case(rng)
            ev = make_event(up, normal, mach, noise=noise, seed=k)
            u, d = windows(ev.data, ev.shock_time, ev.shock_type)
            for m in METHODS:
                n = shock_normal(u.b, d.b, u.v, d.v, m)
                if np.all(np.isfinite(n)):
                    errs[m].append(angle_between_deg(n, ev.truth.normal))
        cells = [f"{np.median(errs[m]):5.1f} / {np.percentile(errs[m], 90):5.1f}" for m in METHODS]
        print(f"{noise:>6.0%} | " + " | ".join(cells))


if __name__ == "__main__":
    main()
