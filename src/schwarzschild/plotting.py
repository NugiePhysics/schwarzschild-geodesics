"""Matplotlib helpers shared by the notebook and the figure script."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from numpy.typing import ArrayLike, NDArray

from .integrators import Trajectory

#: Categorical colors in fixed order (blue, orange, aqua, yellow, magenta, green, violet, red).
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]


def apply_style() -> None:
    """Thin lines, faint grid, no top/right spines, fixed color cycle."""
    plt.rcParams.update(
        {
            "figure.dpi": 110,
            "axes.prop_cycle": plt.cycler(color=COLORS),
            "axes.grid": True,
            "grid.alpha": 0.25,
            "grid.linewidth": 0.6,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "lines.linewidth": 1.6,
            "legend.frameon": False,
        }
    )


def to_cartesian(r: ArrayLike, phi: ArrayLike) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Orbital-plane Cartesian coordinates ``(r cos(phi), r sin(phi))``."""
    r, phi = np.asarray(r, dtype=float), np.asarray(phi, dtype=float)
    return r * np.cos(phi), r * np.sin(phi)


def draw_black_hole(
    ax: Axes, M: float = 1.0, photon_sphere: bool = True, isco: bool = False
) -> None:
    """Draw the horizon (black disk), photon sphere (dashed) and ISCO (dotted) on ``ax``."""
    ax.add_patch(plt.Circle((0, 0), 2.0 * M, color="black", zorder=5))
    theta = np.linspace(0, 2 * np.pi, 400)
    if photon_sphere:
        ax.plot(
            3 * M * np.cos(theta), 3 * M * np.sin(theta), ls="--", lw=0.8, color="0.5", zorder=4
        )
    if isco:
        ax.plot(6 * M * np.cos(theta), 6 * M * np.sin(theta), ls=":", lw=0.8, color="0.5", zorder=4)
    ax.set_aspect("equal")
    ax.set_xlabel("x / M")
    ax.set_ylabel("y / M")


def plot_photon_paths(
    ax: Axes, photons: dict[float, Trajectory], r0: float, M: float = 1.0, half_width: float = 30.0
) -> None:
    """Plot photon trajectories keyed by impact parameter ``b``, all started at radius ``r0``.

    Each path is rotated so the photon arrives parallel to the x axis, at ``y = -b``.
    """
    for color, (b, traj) in zip(COLORS, photons.items(), strict=False):
        x, y = traj.xy()
        angle = np.pi + np.arcsin(b / r0)
        xr = x * np.cos(angle) - y * np.sin(angle)
        yr = x * np.sin(angle) + y * np.cos(angle)
        outcome = "captured" if traj.captured else "escapes"
        ax.plot(xr, yr, color=color, label=rf"$b = {b:.2f}M$ ({outcome})")
    draw_black_hole(ax, M)
    ax.set(xlim=(-half_width, half_width), ylim=(-half_width, half_width))
    ax.legend(loc="lower left", fontsize=8)
