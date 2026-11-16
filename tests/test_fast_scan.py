"""The vectorised scanner must agree with the slow reference scanner on real data."""

from pathlib import Path

import pandas as pd
import pytest

from helios.data.datasets import SPACECRAFT
from helios.data.loader import EventData, standardize_mag, standardize_plasma
from helios.detect.fast import scan_fast
from helios.detect.scan import scan

PROBE = Path(__file__).resolve().parents[1] / "data" / "probe"
EVENTS = [("wind", "20180419T234634", "FF"), ("ace", "20180309T172240", "FF"), ("wind", "20180109T003212", "FR")]


def _load(craft, stamp):
    read = lambda kind: pd.read_csv(PROBE / f"{craft}_{kind}_{stamp}.csv", index_col=0, parse_dates=True)  # noqa: E731
    return EventData(craft, standardize_mag(read("mag"), SPACECRAFT[craft]["mag"]["candidates"][0]),
                     standardize_plasma(read("plasma"), craft))


@pytest.mark.parametrize("craft,stamp,kind", EVENTS)
def test_fast_scan_matches_reference_scan(craft, stamp, kind):
    data = _load(craft, stamp)
    slow, fast = scan(data), scan_fast(data)
    assert len(slow) == len(fast) == 1
    assert fast[0].params.shock_type == slow[0].params.shock_type == kind
    assert abs((fast[0].time - slow[0].time).total_seconds()) <= 1.0
    assert fast[0].n_triggers == slow[0].n_triggers
    assert fast[0].params.b_ratio == pytest.approx(slow[0].params.b_ratio, rel=1e-9)
    assert fast[0].params.shock_speed_kms == pytest.approx(slow[0].params.shock_speed_kms, rel=1e-9)


@pytest.mark.parametrize("craft,stamp,kind", EVENTS)
def test_catalogued_shocks_pass_the_physics_audit(craft, stamp, kind):
    from helios.detect.audit import audit

    data = _load(craft, stamp)
    result = audit(data, scan_fast(data)[0])
    assert result.super_fast and result.mach_fast > 1
    assert result.abrupt and result.sharpness > 0.6
    assert result.passed
