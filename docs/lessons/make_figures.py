"""Teaching figures for the lessons in docs/lessons/.

Run with ``uv run python docs/lessons/make_figures.py``; writes to ``docs/figures/lessons/``.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from aero import plotting as P
from aero.boundary_layer import viscous_drag
from aero.boundary_layer.thwaites import LAMBDA_SEPARATION, shape_factor, shear_correlation
from aero.boundary_layer.transition import michel_re_theta
from aero.geometry import NACA4, Joukowski
from aero.panel_method_2d import solve_airfoil, solve_panels
from aero.thin_airfoil import ThinAirfoil
from aero.vortex_lattice_3d import Wing, build_lattice, solve_vlm

DEG = np.pi / 180
OUT = P.FIGURES / "lessons"


def save(fig, name):
    return P.save(fig, name, OUT)


# --- Lesson 1 --------------------------------------------------------------------------------
def l1_glauert_transform():
    theta = np.linspace(0, np.pi, 13)
    x = 0.5 * (1 - np.cos(theta))
    fig, ax = plt.subplots(figsize=(7, 4.2))
    t = np.linspace(0, np.pi, 200)
    ax.plot(0.5 * (1 - np.cos(t)), 0.5 * np.sin(t), color=P.TEXT_MUTED, lw=1.0)
    ax.plot([0, 1], [0, 0], color=P.TEXT, lw=2.5)
    for th, xi in zip(theta, x, strict=True):
        yc = 0.5 * np.sin(th)
        ax.plot([xi, xi], [0, yc], color=P.SERIES[0], lw=0.8, ls=":")
        ax.plot(xi, yc, "o", color=P.SERIES[0], ms=5)
        ax.plot(xi, 0, "o", color=P.SERIES[1], ms=5, zorder=5)
    ax.text(0.0, -0.1, "θ = 0 (leading edge)", fontsize=9, color=P.TEXT_MUTED)
    ax.text(1.0, -0.1, "θ = π (trailing edge)", fontsize=9, color=P.TEXT_MUTED, ha="right")
    ax.text(0.5, 0.57, "equal steps in θ on the circle …", ha="center", color=P.SERIES[0])
    ax.text(
        0.5, 0.08, "… cluster points at both edges of the chord", ha="center", color=P.SERIES[1]
    )
    ax.set_title("Glauert's substitution  x = (1 − cos θ) / 2", pad=14)
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.15, 0.65)
    ax.set_aspect("equal")
    ax.axis("off")
    return save(fig, "l1_glauert_transform.png")


def l1_loading_decomposition():
    ta = ThinAirfoil.from_naca4("2412")
    alpha = 4 * DEG
    alpha_ideal = -ta.A0(0.0)  # A0 = alpha - B0 vanishes at the ideal angle alpha = B0
    x = np.linspace(0.004, 0.999, 500)
    theta = np.arccos(1 - 2 * x)
    flat = 4 * (alpha - alpha_ideal) * (1 + np.cos(theta)) / np.sin(theta)
    total = ta.delta_cp(x, alpha)
    camber = total - flat
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(x, total, color=P.SERIES[0], label="Total loading, α = 4°")
    ax.plot(x, flat, color=P.SERIES[1], label="Flat-plate part  4(α − α_ideal)(1 + cos θ)/sin θ")
    ax.plot(
        x,
        camber,
        color=P.SERIES[2],
        label=f"Camber part at α_ideal = {alpha_ideal / DEG:.2f}°  (4 Σ An sin nθ)",
    )
    ax.set_ylim(0, 5)
    ax.set_xlabel("x / c")
    ax.set_ylabel("ΔCp = Cp,lower − Cp,upper")
    ax.set_title("NACA 2412: loading = angle-of-attack part + camber part")
    ax.legend(loc="upper right")
    return save(fig, "l1_loading_decomposition.png")


# --- Lesson 2 --------------------------------------------------------------------------------
def l2_panel_geometry():
    nodes = NACA4.from_designation("0012").coordinates(16)
    mid = 0.5 * (nodes[:-1] + nodes[1:])
    d = np.diff(nodes, axis=0)
    t = d / np.linalg.norm(d, axis=1, keepdims=True)
    n = np.column_stack([-t[:, 1], t[:, 0]])
    fig, ax = plt.subplots(figsize=(9, 3.8))
    ax.fill(nodes[:, 0], nodes[:, 1], color=P.GRID)
    ax.plot(nodes[:, 0], nodes[:, 1], "-", color=P.TEXT_MUTED, lw=1.2)
    ax.plot(nodes[:, 0], nodes[:, 1], "o", color=P.TEXT, ms=4, label="Nodes (panel end points)")
    ax.plot(mid[:, 0], mid[:, 1], "o", color=P.SERIES[1], ms=5, label="Control points (midpoints)")
    ax.quiver(
        mid[:, 0],
        mid[:, 1],
        n[:, 0],
        n[:, 1],
        color=P.SERIES[0],
        scale=18,
        width=0.004,
        label="Outward normals n̂",
    )
    for i in (0, len(mid) - 1):
        ax.plot(nodes[i : i + 2, 0], nodes[i : i + 2, 1], color=P.SERIES[3], lw=3)
    for i in (0, 3, 7, 12, 15):
        lx = mid[i, 0] - (0.035 if i in (0, 15) else 0.0)
        ly = mid[i, 1] + (0.022 if mid[i, 1] < 0 else -0.022)
        ax.text(lx, ly, str(i + 1), fontsize=8, color=P.TEXT, ha="center", va="center")
    ax.text(1.03, 0.06, "Kutta condition:\npanels 1 and N", fontsize=9, color=P.TEXT)
    ax.annotate(
        "",
        (0.35, -0.125),
        xytext=(0.65, -0.125),
        arrowprops={"arrowstyle": "->", "color": P.TEXT_MUTED},
    )
    ax.text(
        0.5,
        -0.155,
        "numbering runs clockwise: TE → lower → LE → upper → TE",
        ha="center",
        fontsize=9,
        color=P.TEXT_MUTED,
    )
    ax.set_aspect("equal")
    ax.set_xlim(-0.08, 1.25)
    ax.set_ylim(-0.17, 0.14)
    ax.axis("off")
    ax.set_title("Discretising an airfoil into N = 16 straight panels")
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, -0.02), ncol=3, fontsize=8)
    return save(fig, "l2_panel_geometry.png")


def l2_cylinder():
    fig, ax = plt.subplots(figsize=(7, 4))
    th = np.linspace(0, 2 * np.pi, 400)
    ax.plot(
        np.degrees(th),
        1 - 4 * np.sin(th) ** 2,
        color=P.REFERENCE,
        ls="--",
        lw=1.2,
        label="Exact  1 − 4 sin²θ",
    )
    for i, n in enumerate((8, 32)):
        ang = -np.linspace(0, 2 * np.pi, n + 1)
        nodes = np.column_stack([np.cos(ang), np.sin(ang)])
        sol = solve_panels(nodes, 0.0, lifting=False)
        a = np.mod(np.arctan2(*sol.control_points[:, ::-1].T), 2 * np.pi)
        order = np.argsort(a)
        ax.plot(
            np.degrees(a[order]),
            sol.cp[order],
            "o",
            color=P.SERIES[i],
            ms=5,
            label=f"Source panels, N = {n}",
        )
    ax.set_xlabel("Angle around the cylinder θ [deg]")
    ax.set_ylabel("Pressure coefficient Cp")
    ax.set_xticks([0, 90, 180, 270, 360])
    ax.set_title("First test of any panel code: the circular cylinder")
    ax.legend(loc="lower center")
    return save(fig, "l2_cylinder.png")


# --- Lesson 3 --------------------------------------------------------------------------------
def l3_edge_velocity():
    sol = solve_airfoil("0012", 4 * DEG, 200)
    r = viscous_drag(sol, 3e6)
    fig, ax = plt.subplots(figsize=(7, 4))
    for i, s in enumerate((r.upper, r.lower)):
        ax.plot(s.x, s.ue, color=P.SERIES[i], label=f"{s.name} surface")
    ax.plot(r.upper.x[0], 0, "o", color=P.TEXT, ms=6)
    ax.annotate(
        "stagnation point: Ue = 0,\nboth boundary layers start here",
        (r.upper.x[0], 0),
        xytext=(0.12, 0.15),
        fontsize=9,
        color=P.TEXT,
        arrowprops={"arrowstyle": "->", "color": P.TEXT_MUTED},
    )
    ax.annotate(
        "suction peak, then\nadverse pressure gradient\n(Ue falling)",
        (0.05, 1.55),
        xytext=(0.3, 1.38),
        fontsize=9,
        color=P.TEXT,
        arrowprops={"arrowstyle": "->", "color": P.TEXT_MUTED},
    )
    ax.set_xlabel("x / c")
    ax.set_ylabel("Edge velocity Ue / V∞")
    ax.set_title("NACA 0012, α = 4°: what the boundary layer sees")
    ax.legend(loc="upper right")
    return save(fig, "l3_edge_velocity.png")


def l3_thwaites_correlations():
    lam = np.linspace(-0.09, 0.25, 400)
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
    for ax, f, name in (
        (axes[0], shape_factor, "Shape factor H(λ) = δ*/θ"),
        (axes[1], shear_correlation, "Shear parameter ℓ(λ) = τw θ / (μ Ue)"),
    ):
        ax.axvspan(-0.1, 0, color=P.GRID, alpha=0.6)
        ax.plot(lam, f(lam), color=P.SERIES[0])
        ax.axvline(LAMBDA_SEPARATION, color=P.SERIES[7], lw=1.0, ls="--")
        ax.set_xlabel("Pressure-gradient parameter λ = (θ²/ν) dUe/ds")
        ax.set_title(name)
        ax.set_xlim(-0.1, 0.25)
    axes[0].text(-0.087, 3.3, "separation\nλ = −0.09", color=P.SERIES[7], fontsize=9)
    axes[0].text(-0.05, 2.1, "adverse", ha="center", color=P.TEXT_MUTED, fontsize=9)
    axes[0].text(0.12, 2.1, "favourable", ha="center", color=P.TEXT_MUTED, fontsize=9)
    return save(fig, "l3_thwaites_correlations.png")


def l3_transition():
    re = 3e6
    r = viscous_drag(solve_airfoil("0012", 0.0, 240), re)
    s = r.upper
    lam = s.laminar
    rex = re * s.s
    fig, ax = plt.subplots(figsize=(7, 4))
    m = rex > 1e4
    ax.loglog(rex[m], lam.re_theta[m], color=P.SERIES[0], label="Laminar Reθ (Thwaites)")
    ax.loglog(
        rex[m],
        michel_re_theta(rex[m]),
        color=P.REFERENCE,
        ls="--",
        lw=1.2,
        label="Michel's criterion",
    )
    k = s.i_tr
    ax.plot(
        rex[k],
        lam.re_theta[k],
        "o",
        color=P.SERIES[1],
        ms=8,
        label=f"Transition at x/c = {s.x_transition:.2f}",
    )
    ax.set_xlabel("Reynolds number based on arc length Re_s")
    ax.set_ylabel("Momentum-thickness Reynolds number Reθ")
    ax.set_title("NACA 0012, α = 0°, Re = 3×10⁶: predicting transition")
    ax.legend(loc="upper left")
    return save(fig, "l3_transition.png")


# --- Lesson 4 --------------------------------------------------------------------------------
def l4_horseshoe():
    wing = Wing.rectangular(4, chord=1.0)
    lat = build_lattice(wing, 4, 1, spacing="uniform")
    fig = plt.figure(figsize=(8, 5))
    ax = fig.add_subplot(projection="3d")
    c = lat.corners
    for i in range(c.shape[0]):
        ax.plot(c[i, :, 1], c[i, :, 0], c[i, :, 2], color=P.TEXT_MUTED, lw=1)
    for j in range(c.shape[1]):
        ax.plot(c[:, j, 1], c[:, j, 0], c[:, j, 2], color=P.TEXT_MUTED, lw=1)
    far = 3.5
    for j, (a, b) in enumerate(zip(lat.bound_a[0], lat.bound_b[0], strict=True)):
        col = P.SERIES[j % 4]
        ax.plot([a[1], a[1]], [far, a[0]], [0, 0], color=col, lw=1.6)
        ax.plot([a[1], b[1]], [a[0], b[0]], [0, 0], color=col, lw=2.4)
        ax.plot([b[1], b[1]], [b[0], far], [0, 0], color=col, lw=1.6)
        cp = lat.control[0, j]
        ax.scatter(cp[1], cp[0], cp[2], color=P.TEXT, s=18)
    ax.text(-2, 0.1, 0.25, "bound vortex at c/4", color=P.TEXT, fontsize=9)
    ax.text(-2, 0.8, 0.25, "control point at 3c/4", color=P.TEXT, fontsize=9)
    ax.text(1.2, far, 0.1, "trailing legs → ∞", color=P.TEXT, fontsize=9)
    ax.set_xlabel("y (span)")
    ax.set_ylabel("x (downstream)")
    ax.set_zticks([])
    ax.set_zlim(-0.5, 0.5)
    ax.view_init(elev=35, azim=-60)
    ax.set_box_aspect((4, 3.5, 0.6))
    ax.set_title("Four horseshoe vortices on a rectangular wing")
    return save(fig, "l4_horseshoe.png")


def l4_trefftz_downwash():
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.9), sharex=True)
    for i, (name, wing) in enumerate(
        (("Elliptic", Wing.elliptic(8)), ("Rectangular", Wing.rectangular(8)))
    ):
        s = solve_vlm(wing, 5 * DEG, 40, 4)
        eta = 2 * s.y_mid / wing.span
        wn, _ = s.trefftz()
        g = s.strip_gamma
        inner = np.abs(eta) < 0.95  # discrete tip vortices make the outermost strips noisy
        axes[0].plot(eta, g / g.max(), color=P.SERIES[i], label=name)
        axes[1].plot(eta[inner], wn[inner] / s.CL, color=P.SERIES[i], label=name)
    axes[0].set_title("Circulation Γ(y) / Γ_max")
    axes[1].set_title("Trefftz-plane downwash w / CL (normalised)")
    for ax in axes:
        ax.set_xlabel("η = 2y / b")
    axes[1].legend(loc="lower center")
    return save(fig, "l4_trefftz_downwash.png")


# --- Lesson 6 --------------------------------------------------------------------------------
def l6_convergence_orders():
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4))
    jk = Joukowski(0.1, 0.1)
    ns = np.array([50, 100, 200, 400, 800])
    err = [abs(solve_panels(jk.coordinates(n), 5 * DEG).cl / jk.cl(5 * DEG) - 1) for n in ns]
    axes[0].loglog(ns, err, "o-", color=P.SERIES[0], label="Joukowski cl vs exact")
    af = NACA4.from_designation("0012")
    ref = solve_airfoil(af, 5 * DEG, 1600).cl
    err2 = [abs(solve_airfoil(af, 5 * DEG, n).cl / ref - 1) for n in ns]
    axes[0].loglog(ns, err2, "s-", color=P.SERIES[1], label="NACA 0012 cl vs N = 1600")
    axes[0].loglog(ns, err[0] * ns[0] / ns, color=P.TEXT_MUTED, ls=":", lw=1, label="slope −1")
    axes[0].set_xticks(ns, [str(n) for n in ns])
    axes[0].minorticks_off()
    axes[0].set_xlabel("Panels N")
    axes[0].set_ylabel("Relative error in cl")
    axes[0].set_title("Panel method")
    axes[0].legend(loc="lower left")
    nsp = np.array([6, 10, 16, 24, 40, 64, 100])
    e_ref = solve_vlm(Wing.elliptic(8), 4 * DEG, 200, 2).span_efficiency
    for k, sp in enumerate(("cosine", "uniform")):
        e = np.array([solve_vlm(Wing.elliptic(8), 4 * DEG, n, 2, sp).span_efficiency for n in nsp])
        axes[1].loglog(
            nsp,
            np.abs(e - e_ref),
            "o-" if k == 0 else "s--",
            color=P.SERIES[k],
            label=f"{sp} spacing",
        )
    axes[1].loglog(nsp, 1.0 / nsp, color=P.TEXT_MUTED, ls=":", lw=1, label="slope −1")
    axes[1].set_xticks(nsp, [str(n) for n in nsp])
    axes[1].minorticks_off()
    axes[1].set_xlabel("Spanwise panels")
    axes[1].set_ylabel("|e − e_converged|, elliptic AR 8")
    axes[1].set_title("Vortex lattice")
    axes[1].legend(loc="lower left")
    return save(fig, "l6_convergence_orders.png")


def main() -> None:
    P.use_style()
    for f in (
        l1_glauert_transform,
        l1_loading_decomposition,
        l2_panel_geometry,
        l2_cylinder,
        l3_edge_velocity,
        l3_thwaites_correlations,
        l3_transition,
        l4_horseshoe,
        l4_trefftz_downwash,
        l6_convergence_orders,
    ):
        print("wrote", f().relative_to(Path.cwd()))


if __name__ == "__main__":
    main()
