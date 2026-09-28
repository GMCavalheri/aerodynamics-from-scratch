"""Phase 5: planform trade studies with the VLM.

Run with ``uv run python validation/planform_studies.py``. Writes figures to
``docs/figures/`` and a results table to ``validation/results/planform_studies.md``.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from aero import plotting as P
from aero.vortex_lattice_3d import Wing, solve_vlm
from aero.vortex_lattice_3d.studies import (
    ALPHA,
    MESH,
    induced_drag_vs_aspect_ratio,
    sweep_study,
    taper_study,
)

DEG = np.pi / 180
RESULTS = Path(__file__).resolve().parent / "results"
AR = 8


def fig_taper():
    st = taper_study(AR, np.linspace(0, 1, 41))
    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.plot(st.taper, st.e_vlm, color=P.SERIES[0], label="VLM")
    ax.plot(st.taper, st.e_lifting_line, color=P.SERIES[1], label="Lifting line (Glauert)")
    ax.axhline(1.0, color=P.REFERENCE, lw=1.0, ls="--", label="Elliptic planform")
    for e, lam, c in (
        (st.e_vlm, st.optimum_vlm, P.SERIES[0]),
        (st.e_lifting_line, st.optimum_lifting_line, P.SERIES[1]),
    ):
        ax.plot(lam, e.max(), "o", color=c, ms=7)
        ax.annotate(
            f"λ = {lam:.3g}",
            (lam, e.max()),
            textcoords="offset points",
            xytext=(6, -14),
            color=P.TEXT_MUTED,
            fontsize=9,
        )
    ax.set_xlabel("Taper ratio λ = c_tip / c_root")
    ax.set_ylabel("Span efficiency e")
    ax.set_title(f"Unswept tapered wings, AR {AR}: induced-drag efficiency vs taper")
    ax.set_ylim(0.85, 1.01)
    ax.legend(loc="lower right")
    return P.save(fig, "planform_taper.png"), st


def fig_loading():
    wings = [
        Wing.elliptic(AR, name="Elliptic"),
        Wing.rectangular(AR, name="Rectangular"),
        Wing.tapered(AR, 0.4, name="Tapered λ = 0.4"),
        Wing.tapered(AR, 0.0, name="Pointed λ = 0"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4))
    for i, w in enumerate(wings):
        s = solve_vlm(w, ALPHA, **MESH)
        eta = 2 * s.y_mid / w.span
        half = eta >= 0
        axes[0].plot(
            eta[half],
            s.span_loading[half] / s.CL,
            color=P.SERIES[i],
            label=f"{w.name}  (e = {s.span_efficiency:.3f})",
        )
        axes[1].plot(eta[half], s.cl_span[half] / s.CL, color=P.SERIES[i])
    axes[0].set_title("Span loading  c·cl / (c̄·CL)")
    axes[1].set_title("Section lift  cl / CL  (where stall starts)")
    for ax in axes:
        ax.set_xlabel("Semispan station η = 2y / b")
        ax.set_xlim(0, 1)
    axes[1].set_ylim(0, 1.8)
    axes[0].legend(loc="lower left")
    return P.save(fig, "planform_loading.png")


def fig_aspect_ratio():
    ars = np.array([3, 4, 5, 6, 8, 10, 12, 16, 20])
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4))
    rows = []
    # colours follow the planform (same slots as planform_loading.png)
    for slot, name, make in (
        (0, "Elliptic", Wing.elliptic),
        (1, "Rectangular", Wing.rectangular),
        (2, "Tapered λ = 0.4", lambda ar: Wing.tapered(ar, 0.4)),
    ):
        cdi, e = induced_drag_vs_aspect_ratio(make, ars, cl=0.5)
        axes[0].plot(ars, 1e4 * cdi, "o-", ms=4, color=P.SERIES[slot], label=name)
        axes[1].plot(ars, 100 * (1 / e - 1), "o-", ms=4, color=P.SERIES[slot], label=name)
        rows.append((name, cdi, e))
    axes[0].set_ylabel("CDi at CL = 0.5  [counts]")
    axes[0].set_title("Induced drag falls as 1 / AR")
    axes[1].set_ylabel("Penalty vs ideal, 1/e − 1  [%]")
    axes[1].set_title("Excess induced drag over elliptic loading")
    for ax in axes:
        ax.set_xlabel("Aspect ratio AR")
    axes[0].legend(loc="upper right")
    return P.save(fig, "planform_aspect_ratio.png"), ars, rows


def fig_sweep():
    sweeps_deg = np.arange(-45, 61, 5)
    sols = sweep_study(AR, 0.4, sweeps_deg * DEG)
    cla = np.array([s.CL / ALPHA for s in sols])
    e = np.array([s.span_efficiency for s in sols])
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.9))
    axes[0].plot(sweeps_deg, cla, color=P.SERIES[0], label="VLM")
    axes[0].plot(
        sweeps_deg,
        cla[sweeps_deg == 0] * np.cos(sweeps_deg * DEG),
        color=P.REFERENCE,
        lw=1.0,
        ls="--",
        label="CLα(0) · cos Λ",
    )
    axes[0].set_title("Lift-curve slope CLα [1/rad]")
    axes[0].legend(loc="lower center")
    axes[1].plot(sweeps_deg, e, color=P.SERIES[0])
    axes[1].set_title("Span efficiency e")
    for ax in axes[:2]:
        ax.set_xlabel("Quarter-chord sweep Λ [deg]")
    for i, sw in enumerate((-30, 0, 30, 60)):
        s = sols[int(np.nonzero(sweeps_deg == sw)[0][0])]
        eta = 2 * s.y_mid / s.wing.span
        half = eta >= 0
        axes[2].plot(eta[half], s.cl_span[half] / s.CL, color=P.SERIES[i], label=f"Λ = {sw}°")
    axes[2].set_title("Section lift cl / CL")
    axes[2].set_xlabel("Semispan station η")
    axes[2].set_xlim(0, 1)
    axes[2].legend(loc="lower left")
    fig.suptitle(f"Sweep effects, AR {AR}, λ = 0.4", x=0.01, y=1.06, ha="left", fontweight="bold")
    return P.save(fig, "planform_sweep.png"), sweeps_deg, cla, e


def main() -> None:
    P.use_style()
    p1, st = fig_taper()
    p2 = fig_loading()
    p3, ars, rows = fig_aspect_ratio()
    p4, sweeps, cla, e = fig_sweep()
    for p in (p1, p2, p3, p4):
        print("wrote", p)
    lines = [
        f"Optimum taper ratio (AR {AR}, unswept): VLM λ = {st.optimum_vlm:.3g} "
        f"(e = {st.e_vlm.max():.4f}); lifting line λ = {st.optimum_lifting_line:.3g} "
        f"(e = {st.e_lifting_line.max():.4f}).",
        "",
        "| Planform | " + " | ".join(f"AR {a}" for a in ars) + " |",
        "|---|" + "---|" * len(ars),
    ]
    for name, cdi, eff in rows:
        lines.append(f"| {name}: e | " + " | ".join(f"{x:.3f}" for x in eff) + " |")
        lines.append(f"| {name}: CDi (CL 0.5) | " + " | ".join(f"{x:.4f}" for x in cdi) + " |")
    lines += ["", "| Sweep Λc/4 | CLα [1/rad] | e |", "|---|---|---|"]
    for sw, c, ef in zip(sweeps, cla, e, strict=True):
        if sw % 15 == 0:
            lines.append(f"| {sw}° | {c:.3f} | {ef:.4f} |")
    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "planform_studies.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
