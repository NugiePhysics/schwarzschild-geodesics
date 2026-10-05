"""Timelike and null geodesics around a Schwarzschild black hole.

Hamilton's equations are integrated with RK4 in the equatorial plane. Plotting helpers live in
:mod:`schwarzschild.plotting` and closed-form reference results in :mod:`schwarzschild.analytic`.
"""

from . import analytic
from .diagnostics import (
    convergence_study,
    deflection_angle,
    periapsis_angles,
    precession_per_orbit,
)
from .equations import hamiltonian_constraint, relative_constraint_error, rhs
from .initial_conditions import energy_at_turning_point, initial_radial_momentum, photon_state
from .integrators import Trajectory, integrate, integrate_fixed, rk4_step
from .metric import Metric, Schwarzschild, effective_potential

__version__ = "0.1.0"

__all__ = [
    "Metric",
    "Schwarzschild",
    "Trajectory",
    "__version__",
    "analytic",
    "convergence_study",
    "deflection_angle",
    "effective_potential",
    "energy_at_turning_point",
    "hamiltonian_constraint",
    "initial_radial_momentum",
    "integrate",
    "integrate_fixed",
    "periapsis_angles",
    "photon_state",
    "precession_per_orbit",
    "relative_constraint_error",
    "rhs",
    "rk4_step",
]
