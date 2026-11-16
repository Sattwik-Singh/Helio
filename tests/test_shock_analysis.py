import numpy as np

from conftest import random_case
from helios.physics.normals import METHODS, angle_between_deg, shock_normal
from helios.physics.rankine_hugoniot import PlasmaState, solve_fast_shock
from helios.physics.shock import analyze_shock, hsd_electron_temperature


def test_coplanarity_normals_are_exact_on_exact_shocks():
    """MC and the three mixed modes follow from the jump conditions, so they must recover n exactly."""
    rng = np.random.default_rng(3)
    worst = {m: 0.0 for m in METHODS if m != "VC"}
    for _ in range(300):
        up, normal, mach = random_case(rng)
        sol = solve_fast_shock(up, normal, mach)
        if not 5 < sol.theta_bn_deg < 85:  # MC is singular at 0 and 90 degrees
            continue
        for m in worst:
            n = shock_normal(sol.upstream.b, sol.downstream.b, sol.upstream.v, sol.downstream.v, m)
            worst[m] = max(worst[m], angle_between_deg(n, sol.normal))
    assert max(worst.values()) < 1e-4


def test_normal_points_upstream_for_forward_and_reverse_shocks():
    rng = np.random.default_rng(4)
    for _ in range(200):
        up, normal, mach = random_case(rng)
        sol = solve_fast_shock(up, normal, mach)
        n = shock_normal(sol.upstream.b, sol.downstream.b, sol.upstream.v, sol.downstream.v, "MX3")
        assert n @ sol.normal > 0.999


def test_velocity_coplanarity_is_only_approximate(upstream):
    """VC ignores the tangential velocity jump: useful, but biased, and better at high Mach number."""
    errs = []
    for mach in (1.5, 3.0, 8.0):
        sol = solve_fast_shock(upstream, [-0.9, 0.3, 0.2], mach)
        n = shock_normal(sol.upstream.b, sol.downstream.b, sol.upstream.v, sol.downstream.v, "VC")
        errs.append(angle_between_deg(n, sol.normal))
    assert errs[0] > 1.0 and errs[0] > errs[1] > errs[2]


def test_analysis_recovers_the_parameters_the_shock_was_built_with():
    rng = np.random.default_rng(5)
    for _ in range(200):
        up, normal, mach = random_case(rng)
        sol = solve_fast_shock(up, normal, mach)
        p = analyze_shock(sol.upstream, sol.downstream)
        assert np.isclose(p.shock_speed_kms, sol.shock_speed, rtol=1e-8, atol=1e-6)
        assert np.isclose(p.mach_fast, sol.mach_fast, rtol=1e-8)
        assert np.isclose(p.mach_alfven, sol.mach_alfven, rtol=1e-8)
        assert np.isclose(p.n_ratio, sol.compression, rtol=1e-12)
        assert np.isclose(p.theta_bn_deg, sol.theta_bn_deg, atol=1e-5)


def test_forward_and_reverse_classification(upstream):
    ff = solve_fast_shock(upstream, [-0.95, 0.2, 0.1], 2.0)
    fr = solve_fast_shock(upstream, [0.95, 0.2, -0.1], 2.0)
    assert ff.shock_speed > 0 and analyze_shock(ff.upstream, ff.downstream).shock_type == "FF"
    assert fr.shock_speed < 0 and analyze_shock(fr.upstream, fr.downstream).shock_type == "FR"


def test_criteria_reject_weak_and_non_shock_jumps(upstream):
    weak = solve_fast_shock(upstream, [-0.95, 0.2, 0.1], 1.03)
    assert not analyze_shock(weak.upstream, weak.downstream).meets_criteria
    # density and field rise but the speed does not change: a pressure-balanced structure, not a shock
    fake = PlasmaState(b=np.asarray(upstream.b) * 1.5, v=upstream.v, n=upstream.n * 1.5, tp=upstream.tp)
    p = analyze_shock(upstream, fake)
    assert not p.meets_criteria and p.shock_type == "none"


def test_catalogue_electron_temperature_convention(upstream):
    assert hsd_electron_temperature("2018-03-09T17:08:21Z") == {"te_k": 146277.0}
    assert hsd_electron_temperature("2003-10-29T06:00:00Z") == {"te_k": 146277.0}  # same before 2015: see docstring
    assert np.isclose(hsd_electron_temperature("2020-01-01", r_au=0.5)["te_k"], 146277.0 * 0.5**-0.664)
    # hotter electrons raise the sound speed and beta, and lower the Mach number; the geometry is untouched
    sol = solve_fast_shock(upstream, [-0.9, 0.3, 0.2], 2.0)
    cold = analyze_shock(sol.upstream, sol.downstream, te_over_tp=1.0)
    hot = analyze_shock(sol.upstream, sol.downstream, te_k=146277.0)
    assert hot.sound_speed_kms > cold.sound_speed_kms and hot.plasma_beta > cold.plasma_beta
    assert hot.mach_magnetosonic < cold.mach_magnetosonic
    assert hot.normal == cold.normal and hot.shock_speed_kms == cold.shock_speed_kms
