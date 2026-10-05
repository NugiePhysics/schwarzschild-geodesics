import numpy as np
import pytest

from schwarzschild import (
    analytic,
    effective_potential,
    energy_at_turning_point,
    integrate,
    precession_per_orbit,
)


def test_circular_orbit_stays_circular_with_keplerian_frequency():
    r_c = 10.0
    E, L = analytic.circular_orbit_constants(r_c)
    traj = integrate([0.0, r_c, 0.0, 0.0], E, L, eps=1.0, h0=0.1, lam_max=400.0)
    assert traj.reason == "lam_max"
    assert np.abs(traj.r - r_c).max() < 1e-8
    omega = (traj.phi[-1] - traj.phi[0]) / (traj.t[-1] - traj.t[0])
    assert omega == pytest.approx(analytic.circular_orbit_angular_frequency(r_c), rel=1e-9)


def test_circular_orbit_constants_at_isco():
    E, L = analytic.circular_orbit_constants(analytic.isco_radius())
    assert E == pytest.approx(np.sqrt(8 / 9))
    assert L == pytest.approx(np.sqrt(12))


def test_no_circular_orbit_inside_photon_sphere():
    with pytest.raises(ValueError):
        analytic.circular_orbit_constants(3.0)


@pytest.fixture(scope="module")
def eccentric_orbit():
    """Bound orbit started at apoapsis, r_0 = 20M and L = 4.2M, over about 4 radial periods."""
    r0, L = 20.0, 4.2
    E = energy_at_turning_point(r0, L)
    return integrate([0.0, r0, 0.0, 0.0], E, L, eps=1.0, h0=0.1, lam_max=1800.0)


def test_eccentric_orbit_turning_points(eccentric_orbit):
    r_peri = eccentric_orbit.r.min()
    assert eccentric_orbit.r.max() == pytest.approx(20.0, abs=1e-3)
    assert r_peri == pytest.approx(10.3308, abs=1e-3)
    E, L = eccentric_orbit.E, eccentric_orbit.L
    assert effective_potential(r_peri, L, 1.0) == pytest.approx(E**2, abs=1e-6)


def test_hamiltonian_constraint_is_conserved(eccentric_orbit):
    assert np.abs(eccentric_orbit.constraint).max() < 1e-12


def test_numerical_precession_matches_elliptic_integral_result(eccentric_orbit):
    numerical = precession_per_orbit(eccentric_orbit)
    assert len(numerical) >= 3
    expected = analytic.precession_per_orbit(r_apo=20.0, r_peri=eccentric_orbit.r.min())
    assert numerical.mean() == pytest.approx(expected, rel=1e-5)


def test_precession_vanishes_in_the_newtonian_limit():
    # Weak-field estimate: 6 pi M / (a (1 - e^2)); here r_apo = 2e4 M and r_peri = 1e4 M.
    r_apo, r_peri = 2e4, 1e4
    a, e = 0.5 * (r_apo + r_peri), (r_apo - r_peri) / (r_apo + r_peri)
    weak_field = 6 * np.pi / (a * (1 - e**2))
    assert analytic.precession_per_orbit(r_apo, r_peri) == pytest.approx(weak_field, rel=1e-3)


def test_unbound_orbit_is_rejected_by_analytic_precession():
    with pytest.raises(ValueError):
        analytic.precession_per_orbit(r_apo=20.0, r_peri=4.0)
