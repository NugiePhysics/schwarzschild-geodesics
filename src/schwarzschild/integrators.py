"""Fourth-order Runge-Kutta integration of the equatorial geodesic equations."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .equations import hamiltonian_constraint, relative_constraint_error, rhs
from .metric import Metric, resolve_metric


@dataclass(frozen=True, eq=False)
class Trajectory:
    """Result of :func:`integrate`.

    Attributes
    ----------
    lam
        Affine parameter (proper time for massive particles), shape ``(N,)``.
    state
        States ``(t, r, phi, p_r)``, shape ``(N, 4)``.
    constraint
        Hamiltonian constraint ``C = H + eps/2`` at each step, shape ``(N,)``.
    reason
        Why the integration stopped: ``"horizon"`` (captured), ``"escaped"`` (reached
        ``r_max``), ``"lam_max"`` or ``"max_steps"``.
    E, L, eps, metric
        The parameters the trajectory was computed with.
    """

    lam: NDArray[np.float64]
    state: NDArray[np.float64]
    constraint: NDArray[np.float64]
    reason: str
    E: float
    L: float
    eps: float
    metric: Metric

    @property
    def t(self) -> NDArray[np.float64]:
        return self.state[:, 0]

    @property
    def r(self) -> NDArray[np.float64]:
        return self.state[:, 1]

    @property
    def phi(self) -> NDArray[np.float64]:
        return self.state[:, 2]

    @property
    def p_r(self) -> NDArray[np.float64]:
        return self.state[:, 3]

    @property
    def captured(self) -> bool:
        """Whether the particle fell into the horizon."""
        return self.reason == "horizon"

    def xy(self) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Cartesian coordinates ``(x, y) = (r cos(phi), r sin(phi))`` in the orbital plane."""
        return self.r * np.cos(self.phi), self.r * np.sin(self.phi)

    def relative_constraint_error(self) -> NDArray[np.float64]:
        """``|C| / (E^2 / 2f)`` at each step."""
        return relative_constraint_error(self.state, self.E, self.L, self.eps, self.metric)


def rk4_step(
    y: NDArray[np.float64], h: float, E: float, L: float, metric: Metric | None = None
) -> NDArray[np.float64]:
    """One classical RK4 step of size ``h``."""
    k1 = rhs(y, E, L, metric)
    k2 = rhs(y + 0.5 * h * k1, E, L, metric)
    k3 = rhs(y + 0.5 * h * k2, E, L, metric)
    k4 = rhs(y + h * k3, E, L, metric)
    return y + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def integrate(
    y0: ArrayLike,
    E: float,
    L: float,
    eps: float,
    *,
    metric: Metric | None = None,
    h0: float = 0.01,
    max_steps: int = 500_000,
    alpha: float = 0.02,
    delta: float = 1e-3,
    r_max: float = 1e3,
    lam_max: float = np.inf,
) -> Trajectory:
    """Integrate a geodesic with RK4 and an adaptive step near the horizon.

    The step is ``h = min(h0, alpha * (r - r_h))``, so it shrinks in proportion to the
    distance from the horizon where ``p_r`` and ``dp_r/dlambda`` blow up. Because the
    equations are autonomous, varying ``h`` does not affect the scheme.

    Parameters
    ----------
    y0
        Initial state ``(t, r, phi, p_r)``.
    E, L
        Conserved energy and angular momentum per unit mass.
    eps
        ``1`` for massive particles, ``0`` for photons. Used only for the constraint.
    h0, alpha
        Maximum step and the horizon-proximity factor.
    max_steps, lam_max
        Stop after this many steps, or once the affine parameter reaches ``lam_max``.
    delta
        Stop when ``r <= r_h (1 + delta)`` (the particle is considered captured).
    r_max
        Stop when ``r >= r_max`` (the particle has escaped).

    Raises
    ------
    ValueError
        If ``h0`` or ``alpha`` is not positive, or ``y0`` is not outside the horizon.
    """
    if h0 <= 0 or alpha <= 0:
        raise ValueError("h0 and alpha must be positive")
    metric = resolve_metric(metric)
    r_h = metric.horizon
    y = np.array(y0, dtype=float)
    if y[1] <= r_h:
        raise ValueError(f"the initial radius r = {y[1]} must lie outside the horizon r = {r_h}")
    lam = 0.0
    lams, states = [lam], [y.copy()]
    reason = "max_steps"
    for _ in range(max_steps):
        h = min(h0, alpha * (y[1] - r_h))
        y = rk4_step(y, h, E, L, metric)
        lam += h
        lams.append(lam)
        states.append(y.copy())
        if y[1] <= r_h * (1.0 + delta):
            reason = "horizon"
            break
        if y[1] >= r_max:
            reason = "escaped"
            break
        if lam >= lam_max:
            reason = "lam_max"
            break
    state = np.array(states)
    return Trajectory(
        lam=np.array(lams),
        state=state,
        constraint=hamiltonian_constraint(state, E, L, eps, metric),
        reason=reason,
        E=E,
        L=L,
        eps=eps,
        metric=metric,
    )


def integrate_fixed(
    y0: ArrayLike,
    E: float,
    L: float,
    h: float,
    lam_end: float,
    metric: Metric | None = None,
) -> NDArray[np.float64]:
    """Fixed-step RK4 up to ``lam_end``; returns only the final state.

    ``h`` is adjusted slightly, if needed, so that a whole number of steps ends exactly at
    ``lam_end``. Meant for convergence studies far from the horizon, where the adaptive step is
    not needed.
    """
    if h <= 0 or lam_end <= 0:
        raise ValueError("h and lam_end must be positive")
    n_steps = max(1, round(lam_end / h))
    h = lam_end / n_steps
    y = np.array(y0, dtype=float)
    for _ in range(n_steps):
        y = rk4_step(y, h, E, L, metric)
    return y
