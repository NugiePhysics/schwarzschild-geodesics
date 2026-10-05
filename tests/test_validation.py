"""Invalid input must fail loudly instead of returning meaningless numbers."""

import numpy as np
import pytest

from schwarzschild import (
    Schwarzschild,
    deflection_angle,
    energy_at_turning_point,
    initial_radial_momentum,
    integrate,
    integrate_fixed,
    rk4_step,
)


@pytest.mark.parametrize("r0", [1.5, 2.0])
def test_integrate_rejects_start_inside_horizon(r0):
    with pytest.raises(ValueError, match="outside the horizon"):
        integrate([0.0, r0, 0.0, 0.0], 1.0, 1.0, 1.0)


def test_horizon_check_uses_the_metric_mass():
    with pytest.raises(ValueError, match="outside the horizon"):
        integrate([0.0, 5.0, 0.0, 0.0], 1.0, 1.0, 1.0, metric=Schwarzschild(M=3.0))


@pytest.mark.parametrize("kwargs", [{"h0": 0.0}, {"h0": -0.1}, {"alpha": 0.0}])
def test_integrate_rejects_non_positive_steps(kwargs):
    with pytest.raises(ValueError, match="positive"):
        integrate([0.0, 10.0, 0.0, 0.0], 0.96, 3.8, 1.0, **kwargs)


@pytest.mark.parametrize("func", [initial_radial_momentum, energy_at_turning_point])
def test_initial_conditions_reject_radius_inside_horizon(func):
    args = (1.0, 1.0, 1.0, 1.0) if func is initial_radial_momentum else (1.0, 1.0, 1.0)
    with pytest.raises(ValueError, match="outside the horizon"):
        func(*args)


def test_integrate_fixed_ends_exactly_at_lam_end():
    y0, E, L = [0.0, 20.0, 0.0, 0.0], 0.97, 4.2
    # 1.0 / 0.3 is not an integer: the step must be adjusted to 1/3, not truncated to 3 * 0.3.
    expected = y0
    for _ in range(3):
        expected = rk4_step(np.asarray(expected, dtype=float), 1.0 / 3.0, E, L)
    assert integrate_fixed(y0, E, L, 0.3, 1.0) == pytest.approx(expected, rel=1e-14)


def test_deflection_angle_raises_when_integration_is_cut_short():
    with pytest.raises(RuntimeError, match="stopped early"):
        deflection_angle(40.0, max_steps=10)


def test_trajectories_compare_by_identity():
    run = integrate([0.0, 10.0, 0.0, 0.0], 0.96, 3.8, 1.0, h0=0.1, max_steps=5)
    other = integrate([0.0, 10.0, 0.0, 0.0], 0.96, 3.8, 1.0, h0=0.1, max_steps=5)
    assert run == run
    assert run != other  # no element-wise array comparison (which would raise)
