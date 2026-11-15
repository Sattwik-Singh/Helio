import numpy as np
import pytest

from conftest import random_case
from helios.constants import GAMMA
from helios.physics import mhd
from helios.physics.rankine_hugoniot import PlasmaState, conservation_residuals, solve_fast_shock


def test_all_conservation_laws_hold_for_random_shocks():
    """The downstream state must satisfy every jump condition, checked with an independent formula."""
    rng = np.random.default_rng(1)
    worst = 0.0
    for _ in range(300):
        up, normal, mach = random_case(rng)
        sol = solve_fast_shock(up, normal, mach)
        res = conservation_residuals(sol.upstream, sol.downstream, sol.normal, sol.shock_speed)
        worst = max(worst, max(res.values()))
    assert worst < 1e-10


def test_parallel_shock_reduces_to_gas_dynamics():
    """With B along the normal and c_s > v_A the field drops out: r = (g+1) M^2 / ((g-1) M^2 + 2)."""
    up = PlasmaState(b=[-5.0, 0, 0], v=[-400.0, 0, 0], n=5.0, tp=1e5)
    assert mhd.sound_speed(up.tp) > mhd.alfven_speed(5.0, up.n)
    for mach in (1.2, 2.0, 5.0):
        sol = solve_fast_shock(up, [-1, 0, 0], mach)
        expected = (GAMMA + 1) * mach**2 / ((GAMMA - 1) * mach**2 + 2)
        assert np.isclose(sol.compression, expected, rtol=1e-9)
        assert np.allclose(sol.downstream.b, up.b, atol=1e-9)  # field unchanged


def test_perpendicular_shock_compresses_field_with_density():
    """B perpendicular to the normal is frozen in: B2/B1 = rho2/rho1."""
    up = PlasmaState(b=[0.0, 5.0, 0.0], v=[-420.0, 0, 0], n=6.0, tp=6e4)
    sol = solve_fast_shock(up, [-1, 0, 0], 2.5)
    assert np.isclose(np.linalg.norm(sol.downstream.b) / 5.0, sol.compression, rtol=1e-10)
    assert sol.theta_bn_deg == pytest.approx(90.0)


def test_compression_is_bounded_and_grows_with_mach(upstream):
    rs = [solve_fast_shock(upstream, [-0.9, 0.3, 0.2], m).compression for m in (1.05, 1.5, 2, 4, 10, 50)]
    assert all(1 < r < (GAMMA + 1) / (GAMMA - 1) for r in rs)
    assert all(b > a for a, b in zip(rs, rs[1:]))
    assert rs[-1] > 3.9  # strong-shock limit is 4 for gamma = 5/3


def test_entropy_increases(upstream):
    """A shock is irreversible: p / rho^gamma must rise across it."""
    rng = np.random.default_rng(2)
    for _ in range(100):
        up, normal, mach = random_case(rng)
        sol = solve_fast_shock(up, normal, mach)
        s_up = up.n * up.tp / up.n**GAMMA
        s_down = sol.downstream.n * sol.downstream.tp / sol.downstream.n**GAMMA
        assert s_down > s_up


def test_weak_shock_limit(upstream):
    sol = solve_fast_shock(upstream, [-0.9, 0.3, 0.2], 1.0005)
    assert sol.compression - 1 < 5e-3


def test_sub_fast_flow_is_rejected(upstream):
    with pytest.raises(ValueError):
        solve_fast_shock(upstream, [-1, 0, 0], 0.9)
