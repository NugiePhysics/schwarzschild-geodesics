import numpy as np
import pytest

from schwarzschild import (
    Schwarzschild,
    effective_potential,
    hamiltonian_constraint,
    initial_radial_momentum,
    rhs,
)


def dp_r_explicit(r, p_r, E, L, M):
    """Explicit Schwarzschild form of dp_r/dlambda, eq. (24) of docs/derivation.md."""
    return -M * E**2 / (r - 2 * M) ** 2 - M * p_r**2 / r**2 + L**2 / r**3


@pytest.mark.parametrize("M", [1.0, 2.5])
def test_general_rhs_matches_explicit_schwarzschild(M):
    rng = np.random.default_rng(0)
    metric = Schwarzschild(M)
    E, L = 0.97, 4.0
    for r, p_r in zip(rng.uniform(2.1 * M, 50 * M, 200), rng.uniform(-5, 5, 200), strict=True):
        general = rhs([0.0, r, 0.0, p_r], E, L, metric)[3]
        assert general == pytest.approx(dp_r_explicit(r, p_r, E, L, M), rel=1e-10)


def test_rhs_components():
    M, r, p_r, E, L = 1.0, 10.0, -0.3, 0.95, 3.7
    dt, dr, dphi, _ = rhs([0.0, r, 0.0, p_r], E, L, Schwarzschild(M))
    f = 1 - 2 * M / r
    assert dt == pytest.approx(E / f)
    assert dr == pytest.approx(f * p_r)
    assert dphi == pytest.approx(L / r**2)


def test_effective_potential_expanded_form():
    M, L, r = 1.0, 4.2, np.linspace(2.5, 40.0, 50)
    expected = 1 - 2 * M / r + L**2 / r**2 - 2 * M * L**2 / r**3
    assert effective_potential(r, L, 1.0, Schwarzschild(M)) == pytest.approx(expected)


@pytest.mark.parametrize("eps", [0.0, 1.0])
@pytest.mark.parametrize("sign", [-1.0, 1.0])
def test_initial_momentum_satisfies_constraint(eps, sign):
    r0, E, L = 30.0, 1.0, 4.0
    p_r = initial_radial_momentum(r0, E, L, eps, sign=sign)
    assert np.sign(p_r) == sign
    assert hamiltonian_constraint([0.0, r0, 0.0, p_r], E, L, eps) == pytest.approx(0, abs=1e-14)


def test_forbidden_starting_position_raises():
    with pytest.raises(ValueError, match="not allowed"):
        initial_radial_momentum(5.0, 0.5, 6.0, 1.0)


def test_custom_metric_follows_the_protocol():
    class SameAsSchwarzschild:
        horizon = 2.0

        def f(self, r):
            return 1.0 - 2.0 / np.asarray(r)

        def df(self, r):
            return 2.0 / np.asarray(r) ** 2

    y = [0.0, 8.0, 0.0, -0.2]
    assert rhs(y, 0.96, 3.9, SameAsSchwarzschild()) == pytest.approx(rhs(y, 0.96, 3.9))
