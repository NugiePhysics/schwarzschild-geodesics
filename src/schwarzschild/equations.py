"""Hamilton's equations for equatorial geodesics.

With the conserved energy ``E = -p_t`` and angular momentum ``L = p_phi`` and the equatorial
conditions ``theta = pi/2``, ``p_theta = 0``, the Hamiltonian

    H = 1/2 [ -E^2 / f + f p_r^2 + L^2 / r^2 ]

yields four first-order ODEs for the state ``y = (t, r, phi, p_r)``:

    dt/dlambda   = E / f
    dr/dlambda   = f p_r
    dphi/dlambda = L / r^2
    dp_r/dlambda = -f'/2 (E^2 / f^2 + p_r^2) + L^2 / r^3

The parameter ``eps`` (massive vs. massless) never appears in these equations; it only fixes
the initial radial momentum through the constraint ``H = -eps/2``.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .metric import Metric, resolve_metric


def rhs(y: ArrayLike, E: float, L: float, metric: Metric | None = None) -> NDArray[np.float64]:
    """Right-hand side ``dy/dlambda`` of the equatorial Hamilton equations."""
    metric = resolve_metric(metric)
    _, r, _, p_r = y
    fr = metric.f(r)
    return np.array(
        [
            E / fr,
            fr * p_r,
            L / r**2,
            -0.5 * metric.df(r) * (E**2 / fr**2 + p_r**2) + L**2 / r**3,
        ]
    )


def hamiltonian_constraint(
    state: ArrayLike, E: float, L: float, eps: float, metric: Metric | None = None
) -> NDArray[np.float64]:
    """Constraint ``C = H + eps/2``, which is exactly zero along a geodesic.

    ``state`` may be a single state of shape ``(4,)`` or a trajectory of shape ``(N, 4)``.
    """
    metric = resolve_metric(metric)
    state = np.asarray(state, dtype=float)
    r, p_r = state[..., 1], state[..., 3]
    fr = metric.f(r)
    return 0.5 * (-(E**2) / fr + fr * p_r**2 + L**2 / r**2) + 0.5 * eps


def relative_constraint_error(
    state: ArrayLike, E: float, L: float, eps: float, metric: Metric | None = None
) -> NDArray[np.float64]:
    """``|C| / (E^2 / 2f)``: the constraint violation relative to its largest term.

    Close to the horizon ``E^2 / f`` blows up, so this is a fairer accuracy indicator than
    the absolute value of ``C``.
    """
    metric = resolve_metric(metric)
    state = np.asarray(state, dtype=float)
    scale = E**2 / (2.0 * metric.f(state[..., 1]))
    return np.abs(hamiltonian_constraint(state, E, L, eps, metric)) / scale
