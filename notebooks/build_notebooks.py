"""Generate and execute the demo notebooks from the cell sources below.

Run with ``uv run python notebooks/build_notebooks.py``. Keeping the sources here (instead of
hand-editing .ipynb JSON) makes the notebooks reviewable in diffs and reproducible.
"""

from __future__ import annotations

from pathlib import Path

import nbformat
from nbclient import NotebookClient

HERE = Path(__file__).resolve().parent

SETUP = """\
import numpy as np
import matplotlib.pyplot as plt
from aero import plotting as P

P.use_style()
DEG = np.pi / 180"""

NOTEBOOKS: dict[str, list[tuple[str, str]]] = {
    "01_thin_airfoil.ipynb": [
        (
            "md",
            """\
# 1 · Thin airfoil theory

The camber line is replaced by a vortex sheet whose strength is expanded in a Glauert series
in $\\theta$, with $x = (1 - \\cos\\theta)/2$. Only $A_0, A_1, A_2$ matter for the integrated
loads:

$$c_l = 2\\pi(\\alpha - \\alpha_{L0}), \\qquad c_{m,c/4} = \\frac{\\pi}{4}(A_2 - A_1).$$

See [docs/thin_airfoil.md](../docs/thin_airfoil.md) for the derivation.""",
        ),
        (
            "code",
            SETUP
            + """
from aero.geometry import NACA4
from aero.thin_airfoil import ThinAirfoil""",
        ),
        (
            "code",
            """\
for code in ["0012", "2412", "4412", "6409"]:
    ta = ThinAirfoil.from_naca4(code)
    print(f"NACA {code}:  alpha_L0 = {ta.alpha_zero_lift / DEG:+.3f} deg,  "
          f"cm_c/4 = {ta.cm_quarter_chord:+.4f},  A1 = {ta.A1:.4f}, A2 = {ta.A2:.4f}")""",
        ),
        (
            "md",
            "Anderson's worked example gives $\\alpha_{L0} = -2.077°$ and "
            "$c_{m,c/4} = -0.053$ for the NACA 2412, matching the output above.",
        ),
        (
            "code",
            """\
ta = ThinAirfoil.from_naca4("2412")
x = np.linspace(0.005, 1, 400)
fig, ax = plt.subplots(figsize=(7, 4))
for i, a in enumerate([0, 4, 8]):
    ax.plot(x, ta.delta_cp(x, a * DEG), color=P.SERIES[i], label=f"α = {a}°, cl = {ta.cl(a * DEG):.2f}")
ax.set_ylim(0, 6)
ax.set_xlabel("x / c"); ax.set_ylabel("ΔCp = Cp,lower − Cp,upper")
ax.set_title("NACA 2412 thin-airfoil loading"); ax.legend();""",
        ),
        (
            "md",
            "The leading-edge singularity of the $A_0$ term is the thin-airfoil model of the "
            "suction peak; the camber term ($A_1, A_2$) gives the loading at $\\alpha = \\alpha_{ideal}$.",
        ),
    ],
    "02_panel_method.ipynb": [
        (
            "md",
            """\
# 2 · Hess–Smith panel method

Constant-strength source panels plus one vortex strength shared by all panels. $N$ tangency
conditions at the panel midpoints and the Kutta condition close the $(N+1)$-unknown system.
Details in [docs/panel_method_2d.md](../docs/panel_method_2d.md).""",
        ),
        (
            "code",
            SETUP
            + """
from aero.geometry import NACA4, Joukowski
from aero.panel_method_2d import solve_airfoil, solve_panels""",
        ),
        (
            "code",
            """\
sol = solve_airfoil("2412", 6 * DEG, n_panels=200)
print(f"cl (Kutta-Joukowski) = {sol.cl:.4f}")
print(f"cl (pressure)        = {sol.cl_pressure:.4f}")
print(f"cd (pressure)        = {sol.cd_pressure:.1e}   <- d'Alembert: should be ~0")
print(f"cm_c/4               = {sol.cm:.4f}")""",
        ),
        (
            "code",
            """\
fig, ax = plt.subplots(figsize=(8, 4))
P.streamlines(ax, sol, extent=(-0.4, 1.4, -0.45, 0.45))
ax.set_title("NACA 2412, α = 6°");""",
        ),
        (
            "md",
            "## Exact benchmark: Joukowski airfoil\n\nThe conformal map gives the exact potential "
            "flow, so the panel-method error can be measured directly.",
        ),
        (
            "code",
            """\
jk = Joukowski(0.1, 0.1)
for n in (50, 100, 200, 400, 800):
    s = solve_panels(jk.coordinates(n), 5 * DEG)
    print(f"N = {n:4d}: cl = {s.cl:.4f}  (exact {jk.cl(5 * DEG):.4f}, error {100 * (s.cl / jk.cl(5 * DEG) - 1):+.2f}%)")""",
        ),
        (
            "md",
            "The error falls by less than half per doubling (observed order about 0.7–0.8): the "
            "cusped trailing edge is the hardest case for constant-strength panels.",
        ),
    ],
    "03_boundary_layer.ipynb": [
        (
            "md",
            """\
# 3 · Boundary layer and profile drag

Thwaites (laminar) → Michel (transition) → Head (turbulent) → Squire–Young (drag), all driven
by the inviscid edge velocity from the panel method. See
[docs/boundary_layer.md](../docs/boundary_layer.md).""",
        ),
        (
            "code",
            SETUP
            + """
from aero.panel_method_2d import solve_airfoil
from aero.boundary_layer import thwaites, viscous_drag""",
        ),
        (
            "code",
            """\
x = np.linspace(0, 0.2, 20001)
print("Howarth Ue = 1 - x: Thwaites separation at x =", round(thwaites(x, 1 - x, 1e5).s_sep, 4),
      "(exact 0.1199)")""",
        ),
        (
            "code",
            """\
for re in (1e6, 3e6, 6e6, 9e6):
    r = viscous_drag(solve_airfoil("0012", 0.0, 240), re)
    print(f"Re = {re:.0e}: cd = {r.cd:.5f} (friction {r.cd_friction:.5f}), "
          f"transition at x/c = {r.upper.x_transition:.2f} ({r.upper.transition_cause})")""",
        ),
        (
            "code",
            """\
r = viscous_drag(solve_airfoil("2412", 4 * DEG, 240), 3e6)
fig, ax = plt.subplots(figsize=(7, 3.5))
for i, s in enumerate((r.upper, r.lower)):
    ax.plot(s.x, s.H, color=P.SERIES[i], label=f"{s.name} (transition x = {s.x_transition:.2f})")
ax.set_xlabel("x / c"); ax.set_ylabel("Shape factor H")
ax.set_title(f"NACA 2412, α = 4°, Re = 3e6: cd = {r.cd:.5f}"); ax.legend();""",
        ),
    ],
    "04_vortex_lattice.ipynb": [
        (
            "md",
            """\
# 4 · Vortex lattice method

Horseshoe vortices on the quarter-chord line of each panel, tangency at the three-quarter
chord, induced drag in the Trefftz plane. See [docs/vortex_lattice_3d.md](../docs/vortex_lattice_3d.md).""",
        ),
        (
            "code",
            SETUP
            + """
from aero.vortex_lattice_3d import Wing, solve_vlm, solve_lifting_line""",
        ),
        (
            "code",
            """\
for ar in (4, 8, 20):
    s = solve_vlm(Wing.elliptic(ar), 5 * DEG, 40, 6)
    print(f"Elliptic AR {ar:2d}: CL_alpha = {s.CL / (5 * DEG):.3f}/rad "
          f"(LL {2 * np.pi / (1 + 2 / ar):.3f}),  e = {s.span_efficiency:.4f}")""",
        ),
        (
            "code",
            """\
w = Wing.tapered(6, 0.4, sweep=35 * DEG, twist=lambda eta: -3 * DEG * eta)
s = solve_vlm(w, 6 * DEG, 40, 8)
print(f"CL = {s.CL:.3f}, CDi = {s.CDi:.5f}, e = {s.span_efficiency:.3f}, Cm (root c/4) = {s.Cm:.3f}")
eta = 2 * s.y_mid / w.span
fig, ax = plt.subplots(figsize=(7, 3.5))
ax.plot(eta, s.cl_span, color=P.SERIES[0])
ax.set_xlabel("η"); ax.set_ylabel("section cl")
ax.set_title("Swept, tapered, washed-out wing at α = 6°");""",
        ),
    ],
    "05_planform_studies.ipynb": [
        (
            "md",
            """\
# 5 · Planform trade studies

Reproduces the classical result that elliptic loading minimises induced drag and shows how
taper, aspect ratio and sweep move a wing towards or away from it.""",
        ),
        (
            "code",
            SETUP
            + """
from aero.vortex_lattice_3d.studies import taper_study, sweep_study""",
        ),
        (
            "code",
            """\
st = taper_study(8, np.linspace(0, 1, 21))
print(f"optimum taper: VLM {st.optimum_vlm:.2f} (e = {st.e_vlm.max():.4f}), "
      f"lifting line {st.optimum_lifting_line:.2f} (e = {st.e_lifting_line.max():.4f})")
fig, ax = plt.subplots(figsize=(7, 3.5))
ax.plot(st.taper, st.e_vlm, color=P.SERIES[0], label="VLM")
ax.plot(st.taper, st.e_lifting_line, color=P.SERIES[1], label="Lifting line")
ax.set_xlabel("taper ratio λ"); ax.set_ylabel("e"); ax.legend();""",
        ),
        (
            "code",
            """\
sweeps = np.array([-30, 0, 30, 45])
for sw, s in zip(sweeps, sweep_study(8, 0.4, sweeps * DEG), strict=True):
    print(f"sweep {sw:+3d} deg: CL_alpha = {s.CL / (4 * DEG):.3f}, e = {s.span_efficiency:.4f}, "
          f"outboard cl/CL = {s.cl_span[5] / s.CL:.3f}")""",
        ),
        (
            "md",
            "Aft sweep loads the tips (tip-stall tendency), forward sweep unloads them; both "
            "reduce the lift-curve slope roughly like $\\cos\\Lambda$.",
        ),
    ],
}


def _merge_streams(outputs):
    """Join consecutive stdout/stderr chunks; how the kernel splits them depends on timing."""
    merged = []
    for out in outputs:
        prev = merged[-1] if merged else None
        if (
            prev is not None
            and out.get("output_type") == "stream"
            and prev.get("output_type") == "stream"
            and prev.get("name") == out.get("name")
        ):
            prev["text"] += out["text"]
        else:
            merged.append(out)
    return merged


def build(name: str, cells: list[tuple[str, str]]) -> Path:
    nb = nbformat.v4.new_notebook()
    nb.metadata["kernelspec"] = {
        "name": "python3",
        "display_name": "Python 3",
        "language": "python",
    }
    nb.cells = [
        nbformat.v4.new_markdown_cell(src) if kind == "md" else nbformat.v4.new_code_cell(src)
        for kind, src in cells
    ]
    for i, cell in enumerate(nb.cells):
        cell.id = f"cell-{i:02d}"  # stable ids keep rebuilds diff-free
    NotebookClient(
        nb,
        timeout=600,
        kernel_name="python3",
        record_timing=False,
        resources={"metadata": {"path": str(HERE)}},
    ).execute()
    for cell in nb.cells:
        cell.get("outputs", [])[:] = _merge_streams(cell.get("outputs", []))
    path = HERE / name
    nbformat.write(nb, path)
    return path


def main() -> None:
    for name, cells in NOTEBOOKS.items():
        print("built", build(name, cells).name)


if __name__ == "__main__":
    main()
