"""Phase 3 validation: integral boundary layer and profile drag.

Run with ``uv run python validation/boundary_layer_validation.py``. Writes figures to
``docs/figures/`` and a results table to ``validation/results/boundary_layer.md``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from reference import section_value  # noqa: E402

from aero import plotting as P  # noqa: E402
from aero.boundary_layer import head, thwaites, viscous_drag  # noqa: E402
from aero.panel_method_2d import solve_airfoil  # noqa: E402

DEG = np.pi / 180
RESULTS = Path(__file__).resolve().parent / "results"


def fig_flat_plate():
    re = 1e7
    x = np.linspace(1e-4, 1, 2000)
    lam = thwaites(x, np.ones_like(x), re, theta0=0.0)
    turb = head(x, np.ones_like(x), re, theta0=0.036 * x[0] * (re * x[0]) ** -0.2)
    rex = re * x
    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.loglog(rex, lam.cf, color=P.SERIES[0], label="Thwaites (laminar)")
    ax.loglog(rex, 0.664 / np.sqrt(rex), color=P.REFERENCE, ls="--", lw=1.2, label="Blasius")
    ax.loglog(rex, turb.cf, color=P.SERIES[1], label="Head + Ludwieg–Tillmann (turbulent)")
    ax.loglog(rex, 0.0576 * rex**-0.2, color=P.REFERENCE, ls=":", lw=1.2, label="1/7-power law")
    ax.set_xlim(1e4, 1e7)
    ax.set_xlabel("Reynolds number Re_x")
    ax.set_ylabel("Skin friction coefficient cf")
    ax.set_title("Flat plate skin friction")
    ax.legend(loc="lower left")
    return P.save(fig, "bl_flat_plate.png")


def fig_airfoil_bl(code="NACA 0012", alpha_deg=4.0, re=6e6):
    r = viscous_drag(solve_airfoil(code, alpha_deg * DEG, 240), re)
    fig, axes = plt.subplots(3, 1, figsize=(7, 7.5), sharex=True)
    for i, surf in enumerate((r.upper, r.lower)):
        c = P.SERIES[i]
        m = surf.x >= 0
        axes[0].plot(surf.x[m], 1e3 * surf.theta[m], color=c, label=f"{surf.name} surface")
        axes[1].plot(surf.x[m], surf.H[m], color=c)
        axes[2].plot(surf.x[m], 1e3 * surf.cf[m], color=c)
        if surf.x_transition is not None:
            for ax in axes:
                ax.axvline(surf.x_transition, color=c, lw=0.8, ls=":")
    axes[0].set_ylabel("θ / c  [×10⁻³]")
    axes[1].set_ylabel("Shape factor H")
    axes[2].set_ylabel("cf  [×10⁻³]")
    axes[2].set_ylim(0, 8)  # clip the leading-edge peak (cf -> large as theta -> 0)
    axes[2].set_xlabel("x / c   (dotted: predicted transition)")
    axes[0].set_title(f"{code}, α = {alpha_deg:g}°, Re = {re:.0e}: cd = {r.cd:.5f}")
    axes[0].legend(loc="upper left")
    return P.save(fig, "bl_naca0012.png")


def polar(code, re, alphas):
    out = []
    for a in alphas:
        r = viscous_drag(solve_airfoil(code, a * DEG, 240), re)
        out.append((r.inviscid.cl, r.cd))
    return np.array(out)


def fig_drag_polar(code="NACA 0012"):
    alphas = np.arange(-8, 8.1, 1.0)
    fig, ax = plt.subplots(figsize=(7, 4.2))
    for i, re in enumerate((3e6, 6e6, 9e6)):
        p = polar(code, re, alphas)
        ax.plot(1e4 * p[:, 1], p[:, 0], "o-", ms=4, color=P.SERIES[i], label=f"Re = {re:.0e}")
    v, lo, hi, _ = section_value(code, "cd_min")
    ax.axvspan(1e4 * lo, 1e4 * hi, color=P.GRID, label="Measured cd_min band (AvD, Re 6e6)")
    ax.set_xlabel("Profile drag cd  [counts, ×10⁻⁴]")
    ax.set_ylabel("Lift coefficient cl (inviscid)")
    ax.set_title(f"{code} drag polar: Thwaites / Michel / Head / Squire–Young")
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5))
    return P.save(fig, "drag_polar_naca0012.png")


def results_table() -> str:
    v, lo, hi, src = section_value("NACA 0012", "cd_min")
    lines = [
        "| Case | Computed | Reference | Source |",
        "|---|---|---|---|",
    ]
    x = np.linspace(0, 1, 2001)
    lam = thwaites(x, np.ones_like(x), 1e6)
    coeff = lam.theta[-1] * np.sqrt(1e6)
    lines.append(f"| Flat plate θ√Re_x / x (laminar) | {coeff:.4f} | 0.664 | Blasius |")
    xh = np.linspace(0, 0.2, 20001)
    sep = thwaites(xh, 1 - xh, 1e5).s_sep
    lines.append(f"| Howarth Ue = 1 − x, separation x | {sep:.4f} | 0.1199 | Howarth (exact) |")
    for re in (3e6, 6e6, 9e6):
        r = viscous_drag(solve_airfoil("NACA 0012", 0.0, 240), re)
        ref = f"{v:.4f} ({lo:.4f}–{hi:.4f})" if re == 6e6 else "—"
        lines.append(
            f"| NACA 0012 α = 0°, Re = {re:.0e}: cd (x_tr = {r.upper.x_transition:.2f}) "
            f"| {r.cd:.5f} | {ref} | {src if re == 6e6 else ''} |"
        )
    return "\n".join(lines) + "\n"


def main() -> None:
    P.use_style()
    for path in (fig_flat_plate(), fig_airfoil_bl(), fig_drag_polar()):
        print("wrote", path)
    RESULTS.mkdir(exist_ok=True)
    table = results_table()
    (RESULTS / "boundary_layer.md").write_text(table)
    print(table)


if __name__ == "__main__":
    main()
