import math

import pytest

from schwarzschild import analytic, deflection_angle, integrate, photon_state

R0 = 50.0


def run_photon(b, **kwargs):
    return integrate(photon_state(R0, b), 1.0, b, eps=0.0, h0=0.05, r_max=60.0, **kwargs)


def test_critical_impact_parameter_value():
    assert analytic.critical_impact_parameter() == pytest.approx(5.196152, abs=1e-6)


def test_photon_just_below_critical_impact_parameter_is_captured():
    traj = run_photon(5.19)
    assert traj.captured
    assert traj.relative_constraint_error().max() < 1e-6


def test_photon_just_above_critical_impact_parameter_escapes():
    traj = run_photon(5.20)
    assert traj.reason == "escaped"
    # It winds around the photon sphere before leaving.
    assert traj.r.min() == pytest.approx(3.0687, abs=1e-3)
    assert traj.phi[-1] > 2 * math.pi
    assert traj.relative_constraint_error().max() < 1e-10


def test_weak_field_deflection_matches_expansion():
    b = 40.0
    numerical = deflection_angle(b)
    assert numerical == pytest.approx(analytic.weak_field_deflection(b, order=3), rel=2e-3)
    # The leading-order Einstein value 4M/b is clearly worse.
    assert abs(numerical - analytic.weak_field_deflection(b, order=1)) > 0.05 * numerical


def test_deflection_diverges_near_critical_impact_parameter():
    b_c = analytic.critical_impact_parameter()
    assert deflection_angle(b_c + 1e-3) > 2 * math.pi
    assert math.isnan(deflection_angle(b_c - 0.1))


def test_deflection_series_orders_converge():
    b = 20.0
    values = [analytic.weak_field_deflection(b, order=n) for n in (1, 2, 3, 4)]
    steps = [abs(values[i + 1] - values[i]) for i in range(3)]
    assert steps[0] > steps[1] > steps[2]
