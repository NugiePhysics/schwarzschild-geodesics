"""Quantities extracted from numerical trajectories."""

from __future__ import annotations

import math

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .initial_conditions import photon_state
from .integrators import Trajectory, integrate, integrate_fixed
from .metric import Metric


def periapsis_angles(traj: Trajectory) -> NDArray[np.float64]:
    """Azimuth ``phi`` at each periapsis passage (``p_r`` going from negative to positive).

    The crossing is located by linear interpolation of ``phi`` at ``p_r = 0``.
    """
    p_r, phi = traj.p_r, traj.phi
    idx = np.where((p_r[:-1] < 0) & (p_r[1:] >= 0))[0]
    w = -p_r[idx] / (p_r[idx + 1] - p_r[idx])
    return phi[idx] + w * (phi[idx + 1] - phi[idx])


def precession_per_orbit(traj: Trajectory) -> NDArray[np.float64]:
    """Periapsis advance ``Delta phi - 2 pi`` between consecutive periapsis passages."""
    return np.diff(periapsis_angles(traj)) - 2.0 * math.pi


def deflection_angle(
    b: float,
    *,
    metric: Metric | None = None,
    r_far: float = 1000.0,
    h0: float = 0.5,
    alpha: float = 0.02,
    max_steps: int = 500_000,
) -> float:
    """Total deflection ``alpha = Delta phi - pi`` of a photon with impact parameter ``b``.

    The photon starts and ends at ``r ~ r_far``; the (straight-line) angle ``arcsin(b/r)``
    covered outside that radius is added back at both ends. The remaining error is of order
    ``M b / r_far^2``. Returns ``nan`` if the photon is captured.

    Raises
    ------
    RuntimeError
        If the integration stops before the photon is either captured or back at ``r_far``.
    """
    traj = integrate(
        photon_state(r_far, b, metric),
        1.0,
        b,
        0.0,
        metric=metric,
        h0=h0,
        alpha=alpha,
        r_max=r_far + 1e-6,
        max_steps=max_steps,
    )
    if traj.captured:
        return math.nan
    if traj.reason != "escaped":
        raise RuntimeError(
            f"integration stopped early ({traj.reason}) for b = {b}; increase max_steps or h0"
        )
    dphi = traj.phi[-1] + math.asin(b / r_far) + math.asin(b / traj.r[-1])
    return float(dphi - math.pi)


def convergence_study(
    y0: ArrayLike,
    E: float,
    L: float,
    step_sizes: ArrayLike,
    lam_end: float,
    metric: Metric | None = None,
    reference_factor: int = 16,
) -> tuple[NDArray[np.float64], float]:
    """Global error of fixed-step RK4 for several step sizes, and the measured order.

    The reference solution uses ``min(step_sizes) / reference_factor``. The error is the
    largest difference in ``r`` and ``phi`` at ``lam_end``; the order is the slope of
    ``log(error)`` against ``log(h)`` (4 for RK4).
    """
    step_sizes = np.asarray(step_sizes, dtype=float)
    reference = integrate_fixed(y0, E, L, step_sizes.min() / reference_factor, lam_end, metric)
    errors = np.array(
        [
            np.abs(integrate_fixed(y0, E, L, h, lam_end, metric)[[1, 2]] - reference[[1, 2]]).max()
            for h in step_sizes
        ]
    )
    order = float(np.polyfit(np.log(step_sizes), np.log(errors), 1)[0])
    return errors, order
