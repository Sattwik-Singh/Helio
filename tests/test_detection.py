import numpy as np

from helios.detect.baseline import detect
from helios.evaluation.matching import match_events
from helios.physics.normals import angle_between_deg
from helios.synth.generator import make_event, make_quiet


def test_clean_shock_is_found_at_the_right_time(upstream):
    for normal, kind in (([-0.95, 0.2, 0.1], "FF"), ([0.95, 0.2, -0.1], "FR")):
        ev = make_event(upstream, normal, 2.0, cadence_s=60.0)
        assert ev.shock_type == kind
        found = detect(ev.data)
        assert len(found) == 1
        assert abs(found[0].time - ev.shock_time) <= 60.0
        assert found[0].params.shock_type == kind
        assert angle_between_deg(found[0].params.normal, ev.truth.normal) < 0.5
        assert np.isclose(found[0].params.shock_speed_kms, ev.truth.shock_speed, rtol=0.01)


def test_noisy_shocks_are_recovered_and_quiet_wind_is_left_alone(upstream):
    """5% correlated noise: every one of 40 shocks found within 2 minutes, and no alarms in 40 quiet hours."""
    hits = 0
    for seed in range(40):
        ev = make_event(upstream, [-0.9, 0.3, 0.2], 2.0, noise=0.05, seed=seed)
        m = match_events([d.time for d in detect(ev.data)], [ev.shock_time], tolerance_s=120)
        hits += len(m.pairs)
    assert hits == 40
    false_alarms = sum(len(detect(make_quiet(upstream, noise=0.05, seed=seed))) for seed in range(40))
    assert false_alarms == 0


def test_gradual_compression_fools_the_baseline(upstream):
    """Documented weakness: a 10-minute ramp passes the criteria. Its sharpness score is what gives it away."""
    sharp = detect(make_event(upstream, [-0.9, 0.3, 0.2], 2.0, ramp_s=5.0).data)
    slow = detect(make_event(upstream, [-0.9, 0.3, 0.2], 2.0, ramp_s=300.0).data)
    assert len(sharp) == 1 and len(slow) == 1
    assert sharp[0].sharpness > 0.9
    assert slow[0].sharpness < 0.5


def test_matching_scores():
    m = match_events([100, 1000, 5000, 5100], [130, 5050, 9000], tolerance_s=300)
    assert m.pairs == [(0, 0), (2, 1)] or m.pairs == [(0, 0), (3, 1)]
    assert len(m.false_positives) == 2 and m.misses == [2]
    assert np.isclose(m.precision, 0.5) and np.isclose(m.recall, 2 / 3)
