import numpy as np
import pytest

from aero.geometry import Joukowski
from aero.panel_method_2d import solve_airfoil, solve_panels
from aero.thin_airfoil import ThinAirfoil

DEG = np.pi / 180


def cylinder_nodes(n, radius=1.0):
    # clockwise, starting at theta = 0
    th = -np.linspace(0, 2 * np.pi, n + 1)
    return np.column_stack([radius * np.cos(th), radius * np.sin(th)])


def test_cylinder_pressure_distribution():
    # Exact potential flow: Cp = 1 - 4 sin^2(theta) (Anderson Sec. 3.13). On a regular polygon
    # constant source panels reproduce it exactly at the midpoints, so demand round-off accuracy.
    sol = solve_panels(cylinder_nodes(128), 0.0, lifting=False)
    xc, yc = sol.control_points.T
    th = np.arctan2(yc, xc)
    assert np.allclose(sol.cp, 1 - 4 * np.sin(th) ** 2, atol=1e-12)
    assert np.sum(sol.sigma * sol.lengths) == pytest.approx(0.0, abs=1e-10)


def test_cylinder_lifting_symmetric_gives_no_lift():
    sol = solve_panels(cylinder_nodes(64), 0.0)
    assert sol.cl == pytest.approx(0.0, abs=1e-10)


def test_symmetric_airfoil_zero_lift_at_zero_alpha():
    sol = solve_airfoil("0012", 0.0)
    assert sol.cl == pytest.approx(0.0, abs=1e-10)
    assert sol.cm == pytest.approx(0.0, abs=1e-10)


def test_net_source_vanishes_with_refinement():
    # A closed body has zero net source strength; the discrete Hess-Smith solution only
    # satisfies this to O(1/N), so check first-order convergence to zero.
    net = [
        abs(np.sum(s.sigma * s.lengths))
        for s in (solve_airfoil("2412", 5 * DEG, n) for n in (100, 200, 400))
    ]
    assert net[2] < 1e-3
    assert net[0] / net[1] == pytest.approx(2.0, rel=0.1)
    assert net[1] / net[2] == pytest.approx(2.0, rel=0.1)


def test_kutta_joukowski_matches_pressure_integration():
    sol = solve_airfoil("2412", 6 * DEG, n_panels=300)
    assert sol.cl_pressure == pytest.approx(sol.cl, rel=5e-3)
    assert abs(sol.cd_pressure) < 2e-3  # d'Alembert: zero in exact potential flow


def test_flow_tangency_in_field():
    sol = solve_airfoil("2412", 4 * DEG)
    # velocity evaluated just outside the surface is tangent to it
    xc, yc = (sol.control_points + 1e-9 * sol.normals).T
    u, v = sol.velocity(xc, yc)
    assert np.allclose(u * sol.normals[:, 0] + v * sol.normals[:, 1], 0, atol=1e-6)
    tx, ty = sol.normals[:, 1], -sol.normals[:, 0]
    assert np.allclose(u * tx + v * ty, sol.vt, atol=1e-6)


def test_grid_convergence():
    cls = [solve_airfoil("0012", 5 * DEG, n).cl for n in (50, 100, 200, 400)]
    diffs = np.abs(np.diff(cls))
    assert np.all(diffs[1:] < diffs[:-1])
    assert diffs[-1] < 1e-3


def test_lift_slope_includes_thickness_effect():
    # Thickness raises the inviscid lift slope above 2 pi; a classical estimate is
    # 2 pi (1 + 0.77 t/c) (Abbott & von Doenhoff, Sec. 4.3 / Moran Sec. 4.9).
    a1, a2 = -2 * DEG, 4 * DEG
    s1, s2 = solve_airfoil("0012", a1), solve_airfoil("0012", a2)
    slope = (s2.cl - s1.cl) / (a2 - a1)
    assert slope == pytest.approx(2 * np.pi * (1 + 0.77 * 0.12), rel=0.02)


def test_thin_section_approaches_thin_airfoil_theory():
    ta = ThinAirfoil.from_naca4("2402")
    s0, s1 = solve_airfoil("2402", 0.0, 300), solve_airfoil("2402", 4 * DEG, 300)
    alpha_l0 = -s0.cl / ((s1.cl - s0.cl) / (4 * DEG))
    assert alpha_l0 / DEG == pytest.approx(ta.alpha_zero_lift / DEG, abs=0.1)
    assert s1.cm == pytest.approx(ta.cm_quarter_chord, abs=0.005)


def test_joukowski_exact_solution():
    # Exact conformal-map solution (Anderson Sec. 4.14). The cusped trailing edge holds
    # Hess-Smith below first order (observed ~0.7-0.8), so require the error to drop by more
    # than 1.5x per doubling (order > 0.58).
    jk = Joukowski(0.1, 0.1)
    alpha = 5 * DEG
    errs, cp_errs = [], []
    for n in (200, 400, 800):
        sol = solve_panels(jk.coordinates(n), alpha)
        errs.append(abs(sol.cl - jk.cl(alpha)) / jk.cl(alpha))
        away_from_cusp = sol.control_points[:, 0] < 0.9
        cp_err = np.abs(sol.cp - jk.surface_cp(sol.control_points, alpha))[away_from_cusp]
        cp_errs.append(np.median(cp_err))
    assert errs[-1] < 0.01
    assert cp_errs[-1] < 0.01
    assert errs[0] / errs[1] > 1.5 and errs[1] / errs[2] > 1.5
    assert cp_errs[0] / cp_errs[1] > 1.5 and cp_errs[1] / cp_errs[2] > 1.5


def test_joukowski_symmetric_section():
    # Symmetric section: no lift at alpha = 0, and thickness lifts cl above the flat plate.
    jk = Joukowski(0.08, 0.0)
    assert jk.cl(0.0) == pytest.approx(0.0)
    assert jk.cl(4 * DEG) > 2 * np.pi * np.sin(4 * DEG)
