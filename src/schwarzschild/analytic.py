"""Closed-form Schwarzschild results used to validate the numerics.

All functions use geometrized units (``G = c = 1``) and take the mass ``M`` as an argument.
"""

from __future__ import annotations

import math


def critical_impact_parameter(M: float = 1.0) -> float:
    """Photon capture threshold ``b_c = 3 sqrt(3) M``."""
    return 3.0 * math.sqrt(3.0) * M


def isco_radius(M: float = 1.0) -> float:
    """Innermost stable circular orbit, ``r = 6M``."""
    return 6.0 * M


def circular_orbit_constants(r_c: float, M: float = 1.0) -> tuple[float, float]:
    """Energy and angular momentum ``(E, L)`` of a timelike circular orbit of radius ``r_c``.

    The orbit is stable only for ``r_c >= 6M``; circular orbits with ``3M < r_c < 6M`` exist
    but are unstable.
    """
    if r_c <= 3.0 * M:
        raise ValueError("timelike circular orbits require r_c > 3M")
    L = math.sqrt(M * r_c**2 / (r_c - 3.0 * M))
    E = (1.0 - 2.0 * M / r_c) / math.sqrt(1.0 - 3.0 * M / r_c)
    return E, L


def circular_orbit_angular_frequency(r_c: float, M: float = 1.0) -> float:
    """Coordinate angular velocity ``dphi/dt = sqrt(M / r_c^3)`` (Kepler's third law)."""
    return math.sqrt(M / r_c**3)


def weak_field_deflection(b: float, M: float = 1.0, order: int = 3) -> float:
    """Light deflection angle as a series in ``u = M/b`` (valid for ``b >> M``).

    ``alpha = 4u + (15 pi / 4) u^2 + (128 / 3) u^3 + (3465 pi / 64) u^4 + ...``;
    ``order`` (1 to 4) is the number of terms kept. ``order = 1`` is Einstein's ``4M/b``.
    """
    if not 1 <= order <= 4:
        raise ValueError("order must be between 1 and 4")
    u = M / b
    coefficients = [4.0, 15.0 * math.pi / 4.0, 128.0 / 3.0, 3465.0 * math.pi / 64.0]
    return sum(c * u ** (n + 1) for n, c in enumerate(coefficients[:order]))


def _agm(a: float, b: float) -> float:
    """Arithmetic-geometric mean."""
    while abs(a - b) > 1e-15 * a:
        a, b = 0.5 * (a + b), math.sqrt(a * b)
    return a


def _ellipk(m: float) -> float:
    """Complete elliptic integral of the first kind, ``K(m)`` with parameter ``m = k^2``."""
    return 0.5 * math.pi / _agm(1.0, math.sqrt(1.0 - m))


def precession_per_orbit(r_apo: float, r_peri: float, M: float = 1.0) -> float:
    """Periapsis advance per radial period of a bound timelike orbit, in radians.

    With ``u = 1/r`` the orbit equation is ``(du/dphi)^2 = 2M (u - u_1)(u - u_2)(u - u_3)``.
    The turning points ``u_1 = 1/r_apo`` and ``u_2 = 1/r_peri`` fix the third root through
    ``u_1 + u_2 + u_3 = 1/(2M)``, and the angle swept per radial period is

        Delta phi = 4 K(k^2) / sqrt(2M (u_3 - u_1)),   k^2 = (u_2 - u_1) / (u_3 - u_1).

    The precession is ``Delta phi - 2 pi``.
    """
    u_a, u_p = 1.0 / r_apo, 1.0 / r_peri
    u_3 = 1.0 / (2.0 * M) - u_a - u_p
    if not u_a < u_p < u_3:
        raise ValueError("(r_apo, r_peri) do not describe a stable bound orbit")
    m = (u_p - u_a) / (u_3 - u_a)
    return 4.0 * _ellipk(m) / math.sqrt(2.0 * M * (u_3 - u_a)) - 2.0 * math.pi
