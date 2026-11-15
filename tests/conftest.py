import numpy as np
import pytest

from helios.physics.rankine_hugoniot import PlasmaState


@pytest.fixture
def upstream():
    """A typical slow solar wind state at 1 au (GSE: x points to the Sun, so the wind flows along -x)."""
    return PlasmaState(b=[3.0, -4.0, 2.0], v=[-400.0, 10.0, -5.0], n=6.0, tp=8e4)


def random_case(rng):
    b = rng.normal(0, 4, 3)
    v = np.array([-rng.uniform(300, 600), rng.normal(0, 20), rng.normal(0, 20)])
    up = PlasmaState(b=b, v=v, n=rng.uniform(2, 15), tp=rng.uniform(2e4, 3e5))
    normal = np.array([-1.0, 0, 0]) + rng.normal(0, 0.5, 3)
    if rng.random() < 0.3:
        normal[0] = abs(normal[0])  # some reverse shocks
    return up, normal / np.linalg.norm(normal), rng.uniform(1.15, 4.0)
