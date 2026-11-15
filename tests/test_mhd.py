import numpy as np

from helios.physics import mhd


def test_alfven_speed_reference_value():
    # v_A [km/s] = 21.812 * B[nT] / sqrt(n[cm^-3]) for a proton plasma
    assert np.isclose(mhd.alfven_speed(5.0, 5.0), 21.812 * 5 / np.sqrt(5), rtol=1e-4)


def test_sound_speed_reference_value():
    # c_s = sqrt(gamma k (Tp + Te) / m_p); Tp = Te = 1e5 K, gamma = 5/3 -> 52.5 km/s
    assert np.isclose(mhd.sound_speed(1e5), 52.46, rtol=1e-3)


def test_fast_speed_limits():
    va, cs = 60.0, 45.0
    assert np.isclose(mhd.fast_speed(va, cs, np.pi / 2), np.hypot(va, cs))  # perpendicular: sqrt(vA^2 + cs^2)
    assert np.isclose(mhd.fast_speed(va, cs, 0.0), max(va, cs))  # parallel: the larger of the two
    thetas = np.linspace(0, np.pi / 2, 50)
    assert np.all(np.diff(mhd.fast_speed(va, cs, thetas)) >= -1e-12)  # grows monotonically with angle


def test_plasma_beta_is_cs2_over_va2_scaled():
    n, t, b = 7.0, 9e4, 6.0
    beta = mhd.plasma_beta(n, t, b)
    cs, va = mhd.sound_speed(t), mhd.alfven_speed(b, n)
    assert np.isclose(beta, (2 / (5 / 3)) * cs**2 / va**2, rtol=1e-12)
