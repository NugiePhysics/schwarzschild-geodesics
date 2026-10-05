"""Static, spherically symmetric metrics.

The line element is

    ds^2 = -f(r) dt^2 + dr^2 / f(r) + r^2 (dtheta^2 + sin^2(theta) dphi^2),

and everything else in the package only needs ``f(r)``, ``f'(r)`` and the horizon radius.
Supporting another metric (Reissner-Nordstrom, de Sitter-Schwarzschild, ...) therefore means
writing a small class that implements the :class:`Metric` protocol.

Units are geometrized, ``G = c = 1``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import numpy as np
from numpy.typing import ArrayLike, NDArray


class Metric(Protocol):
    """What the integrator needs to know about a spacetime."""

    @property
    def horizon(self) -> float:
        """Coordinate radius of the event horizon."""
        ...

    def f(self, r: ArrayLike) -> NDArray[np.float64]:
        """Metric function, ``g_tt = -f``."""
        ...

    def df(self, r: ArrayLike) -> NDArray[np.float64]:
        """Radial derivative ``f'(r)``."""
        ...


@dataclass(frozen=True)
class Schwarzschild:
    """Schwarzschild black hole of mass ``M``: ``f(r) = 1 - 2M/r``."""

    M: float = 1.0

    @property
    def horizon(self) -> float:
        return 2.0 * self.M

    def f(self, r: ArrayLike) -> NDArray[np.float64]:
        return 1.0 - 2.0 * self.M / np.asarray(r, dtype=float)

    def df(self, r: ArrayLike) -> NDArray[np.float64]:
        return 2.0 * self.M / np.asarray(r, dtype=float) ** 2


def resolve_metric(metric: Metric | None) -> Metric:
    """Return ``metric``, or the default ``Schwarzschild(M=1)`` if it is ``None``."""
    return Schwarzschild() if metric is None else metric


def effective_potential(
    r: ArrayLike, L: float, eps: float, metric: Metric | None = None
) -> NDArray[np.float64]:
    """Effective potential ``V_eff = f (eps + L^2 / r^2)``.

    ``eps = 1`` for massive particles and ``eps = 0`` for photons, so that the radial motion
    obeys ``(dr/dlambda)^2 + V_eff = E^2``.
    """
    metric = resolve_metric(metric)
    r = np.asarray(r, dtype=float)
    return metric.f(r) * (eps + L**2 / r**2)
