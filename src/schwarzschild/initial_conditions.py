"""Initial conditions that satisfy the Hamiltonian constraint ``H = -eps/2``."""

from __future__ import annotations

import numpy as np

from .metric import Metric, effective_potential, resolve_metric


def _require_outside_horizon(r0: float, metric: Metric) -> None:
    if r0 <= metric.horizon:
        raise ValueError(f"r0 = {r0} must lie outside the horizon r = {metric.horizon}")


def initial_radial_momentum(
    r0: float,
    E: float,
    L: float,
    eps: float,
    sign: float = -1.0,
    metric: Metric | None = None,
) -> float:
    """Radial momentum ``p_r(0) = sign * sqrt(E^2 - V_eff(r0)) / f(r0)``.

    ``sign = -1`` for a particle moving inwards, ``+1`` for outwards. The square root is
    needed only here; during the integration the sign of ``p_r`` follows from the ODEs.

    Raises
    ------
    ValueError
        If ``r0`` is not outside the horizon, or ``E^2 < V_eff(r0)`` (a forbidden region).
    """
    metric = resolve_metric(metric)
    _require_outside_horizon(r0, metric)
    D = E**2 - effective_potential(r0, L, eps, metric)
    if D < 0:
        raise ValueError("E^2 < V_eff(r0): the starting position is not allowed")
    return float(sign * np.sqrt(D) / metric.f(r0))


def energy_at_turning_point(
    r0: float, L: float, eps: float = 1.0, metric: Metric | None = None
) -> float:
    """Energy ``E = sqrt(V_eff(r0))`` for which ``r0`` is a turning point (``p_r = 0``)."""
    metric = resolve_metric(metric)
    _require_outside_horizon(r0, metric)
    return float(np.sqrt(effective_potential(r0, L, eps, metric)))


def photon_state(r0: float, b: float, metric: Metric | None = None) -> np.ndarray:
    """State ``(t, r, phi, p_r)`` of an ingoing photon with impact parameter ``b = L/E``.

    Photon paths depend only on ``b``, so use ``E = 1`` and ``L = b`` with this state.
    """
    p_r0 = initial_radial_momentum(r0, 1.0, b, 0.0, sign=-1.0, metric=metric)
    return np.array([0.0, r0, 0.0, p_r0])
