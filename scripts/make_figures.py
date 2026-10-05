"""Regenerate the figures in ``docs/figures`` (used by the README).

Run from the repository root:

    python scripts/make_figures.py
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from schwarzschild import (  # noqa: E402
    analytic,
    convergence_study,
    deflection_angle,
    effective_potential,
    energy_at_turning_point,
    integrate,
    photon_state,
)
from schwarzschild.plotting import (  # noqa: E402
    COLORS,
    apply_style,
    draw_black_hole,
    plot_photon_paths,
)

OUT = Path(__file__).resolve().parents[1] / "docs" / "figures"
M = 1.0
DPI = 150

# Bound orbit started at apoapsis (r_0 = 20M, L = 4.2M).
R0, L_BOUND = 20.0, 4.2
E_BOUND = energy_at_turning_point(R0, L_BOUND)
B_C = analytic.critical_impact_parameter(M)


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / name, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print("wrote", OUT / name)


def photon_runs(impact_parameters, r0=50.0):
    return {
        b: integrate(photon_state(r0, b), 1.0, b, eps=0.0, h0=0.02, r_max=r0 + 10.0)
        for b in impact_parameters
    }


def bound_orbit(lam_max=2600.0):
    return integrate([0.0, R0, 0.0, 0.0], E_BOUND, L_BOUND, eps=1.0, h0=0.05, lam_max=lam_max)


def hero():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5.2))
    plot_photon_paths(ax1, photon_runs([4.0, 5.19, 5.20, 6.0, 8.0]), 50.0, M)
    ax1.set_title(r"Photons: $b_c = 3\sqrt{3}\,M \approx 5.196\,M$")

    ax2.plot(*bound_orbit().xy(), color=COLORS[0], lw=1.0)
    draw_black_hole(ax2, M, isco=True)
    ax2.set_title(r"Bound orbit with periapsis precession ($L = 4.2\,M$)")
    save(fig, "hero.png")


def potentials():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))
    r = np.linspace(2.0 * M, 40.0 * M, 2000)
    for L in [3.0, np.sqrt(12.0), 4.0, 4.2, 5.0]:
        label = r"$L=\sqrt{12}M$ (ISCO)" if np.isclose(L, np.sqrt(12)) else rf"$L={L:.1f}M$"
        ax1.plot(r, effective_potential(r, L, 1.0), label=label)
    ax1.axhline(1.0, color="0.5", lw=0.8, ls=":")
    ax1.set(ylim=(0.7, 1.15), xlabel="r / M", ylabel=r"$V_{\mathrm{eff}}$")
    ax1.set_title(r"Massive particles ($\epsilon = 1$)")
    ax1.legend()

    r = np.linspace(2.0 * M, 15.0 * M, 2000)
    for b in [4.0, B_C, 6.0]:
        label = r"$b=b_c=3\sqrt{3}M$" if b == B_C else rf"$b={b:.1f}M$"
        ax2.plot(r, effective_potential(r, b, 0.0), label=label)
    ax2.axhline(1.0, color="0.5", lw=0.8, ls=":")
    ax2.axvline(3.0 * M, color="0.5", lw=0.8, ls="--")
    ax2.set(xlabel="r / M", ylabel=r"$V_{\mathrm{eff}}$  ($E = 1$)")
    ax2.set_title(r"Photons ($\epsilon = 0$)")
    ax2.legend()
    save(fig, "effective_potential.png")


def deflection():
    b_values = np.concatenate([np.linspace(B_C + 1e-3, 6.0, 15), np.linspace(6.2, 40.0, 25)])
    alpha = [deflection_angle(b) for b in b_values]
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.semilogy(b_values, alpha, "o-", ms=4, color=COLORS[0], label="numerical (RK4)")
    for order, style, color, label in [
        (1, "--", COLORS[1], r"weak field, $4M/b$"),
        (3, ":", COLORS[2], "third-order expansion"),
    ]:
        values = [analytic.weak_field_deflection(b, M, order) for b in b_values]
        ax.semilogy(b_values, values, ls=style, color=color, label=label)
    ax.axvline(B_C, color="0.5", lw=0.8, ls=":")
    ax.set(xlabel="b / M", ylabel=r"deflection angle $\alpha$ (rad)")
    ax.set_title("Light bending by a Schwarzschild black hole")
    ax.legend()
    save(fig, "deflection.png")


def convergence():
    step_sizes = np.array([0.8, 0.4, 0.2, 0.1])
    errors, order = convergence_study([0.0, R0, 0.0, 0.0], E_BOUND, L_BOUND, step_sizes, 400.0)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.loglog(
        step_sizes, errors, "o-", ms=6, color=COLORS[0], label=f"numerical (order {order:.2f})"
    )
    ax.loglog(
        step_sizes,
        errors[0] * (step_sizes / step_sizes[0]) ** 4,
        ls="--",
        color="0.5",
        label=r"$\propto h^4$",
    )
    ax.set(xlabel="h / M", ylabel="global error")
    ax.set_title("RK4 convergence")
    ax.legend()
    save(fig, "convergence.png")


if __name__ == "__main__":
    apply_style()
    hero()
    potentials()
    deflection()
    convergence()
