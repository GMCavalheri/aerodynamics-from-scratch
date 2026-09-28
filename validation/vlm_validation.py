"""Phase 4 validation: vortex lattice method vs classical wing theory.

Run with ``uv run python validation/vlm_validation.py``. Writes figures to ``docs/figures/``
and a results table to ``validation/results/vlm.md``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from reference import load  # noqa: E402

from aero import plotting as P  # noqa: E402
from aero.vortex_lattice_3d import Wing, build_lattice, solve_vlm  # noqa: E402

DEG = np.pi / 180
ALPHA = 4 * DEG
RESULTS = Path(__file__).resolve().parent / "results"


def helmholtz(ar):
    return 2 * np.pi * ar / (2 + np.sqrt(ar**2 + 4))


def lifting_line(ar):
    return 2 * np.pi / (1 + 2 / ar)


def fig_lattice():
    wing = Wing.tapered(6, 0.4, sweep=30 * DEG)
    lat = build_lattice(wing, 24, 6)
    fig, ax = plt.subplots(figsize=(8, 4))
    c = lat.corners
    for i in range(c.shape[0]):
        ax.plot(c[i, :, 1], c[i, :, 0], color=P.GRID, lw=0.8)
    for j in range(c.shape[1]):
        ax.plot(c[:, j, 1], c[:, j, 0], color=P.GRID, lw=0.8)
    a, b = lat.bound_a.reshape(-1, 3), lat.bound_b.reshape(-1, 3)
    for p, q in zip(a, b, strict=True):
        ax.plot([p[1], q[1]], [p[0], q[0]], color=P.SERIES[0], lw=1.2)
    cp = lat.control.reshape(-1, 3)
    ax.plot(cp[:, 1], cp[:, 0], "o", color=P.SERIES[1], ms=2.5, label="Control points (3/4 c)")
    ax.plot([], [], color=P.SERIES[0], label="Bound vortices (1/4 c)")
    ax.invert_yaxis()
    ax.set_aspect("equal")
    ax.set_xlabel("y  (span)")
    ax.set_ylabel("x  (downstream)")
    ax.set_title("Vortex lattice: AR 6, λ = 0.4, Λc/4 = 30°, 24 × 6 panels, cosine spacing")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.32), ncol=2)
    ax.grid(False)
    return P.save(fig, "vlm_lattice.png")


def fig_elliptic_loading():
    fig, ax = plt.subplots(figsize=(7, 4.2))
    eta = np.linspace(-1, 1, 400)
    ax.plot(
        eta, np.sqrt(1 - eta**2), color=P.REFERENCE, ls="--", lw=1.2, label="Elliptic (Prandtl)"
    )
    for i, ar in enumerate((4, 8, 20)):
        s = solve_vlm(Wing.elliptic(ar), ALPHA, 40, 4)
        g = s.strip_gamma
        ax.plot(
            2 * s.y_mid / s.wing.span,
            g / g.max(),
            "o",
            ms=3.5,
            color=P.SERIES[i],
            label=f"VLM, AR {ar}  (e = {s.span_efficiency:.4f})",
        )
    ax.set_xlabel("Spanwise station η = 2y / b")
    ax.set_ylabel("Γ / Γ_root")
    ax.set_title("Elliptic planforms: VLM span loading vs Prandtl")
    ax.legend(loc="lower center")
    return P.save(fig, "vlm_elliptic_loading.png")


def fig_lift_slope():
    ars = np.array([1, 1.5, 2, 3, 4, 6, 8, 12, 16, 20, 30])
    cla = [solve_vlm(Wing.elliptic(ar), ALPHA, 40, 8).CL / ALPHA for ar in ars]
    fine = np.linspace(0.8, 30, 300)
    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.plot(fine, lifting_line(fine), color=P.SERIES[0], label="Lifting line  2π / (1 + 2/AR)")
    ax.plot(fine, helmholtz(fine), color=P.SERIES[2], label="Helmholtz  2πAR / (2 + √(AR² + 4))")
    ax.plot(ars, cla, "o", color=P.SERIES[1], ms=6, label="VLM (elliptic planform)")
    kin = next(r for r in load("vlm_benchmarks.csv") if r["source"] == "[Kinner]")
    ax.plot(
        float(kin["aspect_ratio"]),
        float(kin["cl_alpha_per_rad"]),
        "D",
        color=P.REFERENCE,
        ms=7,
        label="Kinner, exact circular wing",
    )
    ax.axhline(2 * np.pi, color=P.TEXT_MUTED, lw=0.8, ls=":")
    ax.text(29.5, 2 * np.pi - 0.08, "2π (2D)", ha="right", va="top", color=P.TEXT_MUTED, fontsize=9)
    ax.set_xlabel("Aspect ratio AR")
    ax.set_ylabel("Lift-curve slope CLα [1/rad]")
    ax.set_title("Finite-wing lift-curve slope")
    ax.legend(loc="lower right")
    return P.save(fig, "vlm_lift_slope.png")


def fig_convergence():
    ns = np.array([6, 10, 16, 24, 40, 64, 100])
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
    for i, (label, wing) in enumerate(
        (("Elliptic AR 8", Wing.elliptic(8)), ("Rectangular AR 8", Wing.rectangular(8)))
    ):
        for k, spacing in enumerate(("cosine", "uniform")):
            e = [solve_vlm(wing, ALPHA, n, 4, spacing).span_efficiency for n in ns]
            axes[i].semilogx(
                ns,
                e,
                "o-" if k == 0 else "s--",
                ms=4,
                color=P.SERIES[k],
                label=f"{spacing} spacing",
            )
        axes[i].set_title(label)
        axes[i].set_xlabel("Spanwise panels")
    axes[0].axhline(1.0, color=P.REFERENCE, lw=1.0, ls=":")
    axes[0].set_ylabel("Span efficiency e")
    axes[0].legend(loc="upper right")
    return P.save(fig, "vlm_convergence.png")


def results_table() -> str:
    lines = ["| Case | VLM | Reference | Source |", "|---|---|---|---|"]
    for ar in (4, 8, 20):
        s = solve_vlm(Wing.elliptic(ar), ALPHA, 40, 8)
        lines.append(f"| Elliptic AR {ar}: e | {s.span_efficiency:.4f} | 1 | Prandtl |")
        lines.append(
            f"| Elliptic AR {ar}: CLα [1/rad] | {s.CL / ALPHA:.4f} | "
            f"{lifting_line(ar):.4f} (LL) / {helmholtz(ar):.4f} (Helmholtz) | closed form |"
        )
    for row in load("vlm_benchmarks.csv"):
        ar = float(row["aspect_ratio"])
        if row["taper"] == "elliptic":
            s = solve_vlm(Wing.elliptic(ar), 2 * DEG, 40, 20)
        else:
            w = Wing.tapered(ar, float(row["taper"]), sweep=float(row["sweep_c4_deg"]) * DEG)
            s = solve_vlm(w, 2 * DEG, 8, 1)
        lines.append(
            f"| {row['case']} (AR {ar:g}): CLα [1/rad] | {s.CL / (2 * DEG):.4f} | "
            f"{float(row['cl_alpha_per_rad']):.4f} | {row['source']} |"
        )
    return "\n".join(lines) + "\n"


def main() -> None:
    P.use_style()
    for path in (fig_lattice(), fig_elliptic_loading(), fig_lift_slope(), fig_convergence()):
        print("wrote", path)
    RESULTS.mkdir(exist_ok=True)
    table = results_table()
    (RESULTS / "vlm.md").write_text(table)
    print(table)


if __name__ == "__main__":
    main()
