"""Phase 2 validation: Hess-Smith panel method vs exact and published results.

Run with ``uv run python validation/panel_method_validation.py``. Writes figures to
``docs/figures/`` and a results table to ``validation/results/panel_method.md``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from reference import section_value  # noqa: E402

from aero import plotting as P  # noqa: E402
from aero.geometry import Joukowski  # noqa: E402
from aero.panel_method_2d import solve_airfoil, solve_panels  # noqa: E402
from aero.thin_airfoil import ThinAirfoil  # noqa: E402

DEG = np.pi / 180
RESULTS = Path(__file__).resolve().parent / "results"
N_PANELS = 300


def section_summary(code: str) -> dict[str, float]:
    s0, s1 = solve_airfoil(code, 0.0, N_PANELS), solve_airfoil(code, 4 * DEG, N_PANELS)
    slope = (s1.cl - s0.cl) / (4 * DEG)
    ta = ThinAirfoil.from_naca4(code)
    return {
        "panel_cl_alpha": slope * DEG,
        "panel_alpha_L0": -s0.cl / slope / DEG,
        "panel_cm_c4": s1.cm,
        "thin_cl_alpha": 2 * np.pi * DEG,
        "thin_alpha_L0": ta.alpha_zero_lift / DEG,
        "thin_cm_c4": ta.cm_quarter_chord,
    }


def fig_joukowski():
    jk = Joukowski(0.1, 0.1)
    alpha = 5 * DEG
    fig, ax = plt.subplots(figsize=(7, 4.2))
    exact_pts = jk.coordinates(800)[1:-1]  # TE is 0/0 in the map
    ax.plot(
        exact_pts[:, 0],
        jk.surface_cp(exact_pts, alpha),
        color=P.REFERENCE,
        lw=1.2,
        ls="--",
        label="Exact (conformal map)",
        zorder=5,
    )
    for i, n in enumerate((50, 200)):
        sol = solve_panels(jk.coordinates(n), alpha)
        ax.plot(
            sol.control_points[:, 0],
            sol.cp,
            "o",
            color=P.SERIES[i],
            ms=3.5,
            label=f"Panel method, N = {n}  (cl = {sol.cl:.3f})",
        )
    P.plot_cp(ax, [], [], P.SERIES[0])
    ax.set_ylim(1.2, -2.6)
    ax.set_title(f"Joukowski airfoil at α = 5°: exact cl = {jk.cl(alpha):.3f}")
    ax.legend(loc="lower right")
    return P.save(fig, "joukowski_cp.png")


def fig_cp(code: str, alphas=(0, 4, 8)):
    fig, ax = plt.subplots(figsize=(7, 4.2))
    for i, a in enumerate(alphas):
        sol = solve_airfoil(code, a * DEG, N_PANELS)
        xl, cpl, xu, cpu = sol.surfaces()
        ax.plot(xu, cpu, color=P.SERIES[i], label=f"α = {a}°  (cl = {sol.cl:.2f})")
        ax.plot(xl, cpl, color=P.SERIES[i], lw=1.2, ls=(0, (4, 2)))
    P.plot_cp(ax, [], [], P.SERIES[0])
    ax.set_title(f"{code}: inviscid surface pressure (solid upper, dashed lower)")
    ax.legend(loc="upper right")
    return P.save(fig, f"cp_{code.lower().replace(' ', '')}.png")


def fig_lift_curves():
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    alphas = np.linspace(-4, 12, 9)
    for ax, code in zip(axes, ("NACA 0012", "NACA 2412"), strict=True):
        ta = ThinAirfoil.from_naca4(code)
        cl_panel = [solve_airfoil(code, a * DEG, N_PANELS).cl for a in alphas]
        ax.plot(alphas, [ta.cl(a * DEG) for a in alphas], color=P.SERIES[0], label="Thin airfoil")
        ax.plot(alphas, cl_panel, "o-", color=P.SERIES[1], ms=5, label="Panel method")
        a0, _, _, _ = section_value(code, "alpha_L0_deg")
        slope, lo, hi, _ = section_value(code, "cl_alpha_per_deg")
        ax.fill_between(
            alphas,
            lo * (alphas - a0),
            hi * (alphas - a0),
            color=P.GRID,
            label="Measured slope band",
        )
        ax.plot(
            alphas,
            slope * (alphas - a0),
            color=P.REFERENCE,
            lw=1.2,
            ls="--",
            label="Measured (AvD, Re 6e6)",
        )
        ax.set_title(code)
        ax.set_xlabel("Angle of attack α [deg]")
    axes[0].set_ylabel("Lift coefficient cl")
    axes[0].legend(loc="upper left")
    return P.save(fig, "lift_curves.png")


def fig_streamlines():
    sol = solve_airfoil("NACA 2412", 8 * DEG, 200)
    fig, ax = plt.subplots(figsize=(8, 4))
    P.streamlines(ax, sol, extent=(-0.4, 1.4, -0.45, 0.45))
    ax.set_title("NACA 2412 at α = 8°: potential-flow streamlines (shade = speed)")
    ax.set_xlabel("x / c")
    return P.save(fig, "streamlines_naca2412.png")


def results_table() -> str:
    lines = [
        "| Airfoil | Quantity | Thin airfoil | Panel (N=300) | Measured | Source |",
        "|---|---|---|---|---|---|",
    ]
    for code in ("NACA 0012", "NACA 2412"):
        r = section_summary(code)
        for q, key, fmt in (
            ("alpha_L0_deg", "alpha_L0", "{:+.2f}°"),
            ("cl_alpha_per_deg", "cl_alpha", "{:.4f} /deg"),
            ("cm_c4", "cm_c4", "{:+.4f}"),
        ):
            v, _, _, src = section_value(code, q)
            lines.append(
                f"| {code} | {q} | {fmt.format(r['thin_' + key])} | "
                f"{fmt.format(r['panel_' + key])} | {fmt.format(v)} | {src} |"
            )
    jk = Joukowski(0.1, 0.1)
    lines += [
        "",
        "| Joukowski (mx=0.1, my=0.1), α = 5° | N | cl panel | cl exact | error |",
        "|---|---|---|---|---|",
    ]
    for n in (100, 200, 400, 800):
        cl = solve_panels(jk.coordinates(n), 5 * DEG).cl
        ex = jk.cl(5 * DEG)
        lines.append(f"| | {n} | {cl:.4f} | {ex:.4f} | {100 * (cl - ex) / ex:+.2f}% |")
    return "\n".join(lines) + "\n"


def main() -> None:
    P.use_style()
    for path in (
        fig_joukowski(),
        fig_cp("NACA 0012"),
        fig_cp("NACA 2412"),
        fig_lift_curves(),
        fig_streamlines(),
    ):
        print("wrote", path.relative_to(Path.cwd()) if path.is_relative_to(Path.cwd()) else path)
    RESULTS.mkdir(exist_ok=True)
    table = results_table()
    (RESULTS / "panel_method.md").write_text(table)
    print(table)


if __name__ == "__main__":
    main()
