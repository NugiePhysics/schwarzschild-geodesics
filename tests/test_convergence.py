import pytest

from schwarzschild import convergence_study, energy_at_turning_point


def test_rk4_global_error_scales_as_h_to_the_fourth():
    r0, L = 20.0, 4.2
    E = energy_at_turning_point(r0, L)
    errors, order = convergence_study(
        [0.0, r0, 0.0, 0.0], E, L, step_sizes=[0.8, 0.4, 0.2, 0.1], lam_end=400.0
    )
    assert order == pytest.approx(4.0, abs=0.3)
    assert errors[0] > errors[1] > errors[2] > errors[3]
