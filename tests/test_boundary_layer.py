import numpy as np
import pytest

from aero.boundary_layer import head, michel_re_theta, thwaites, viscous_drag
from aero.boundary_layer.head import h1_of_h, h_of_h1
from aero.panel_method_2d import solve_airfoil

DEG = np.pi / 180


def test_thwaites_flat_plate_vs_blasius():
    # Blasius: theta = 0.664 x / sqrt(Re_x), cf = 0.664 / sqrt(Re_x) (White Sec. 4-3).
    # Thwaites gives 0.671 and 0.656, i.e. ~1% off.
    re = 1e6
    x = np.linspace(0, 1, 2001)
    bl = thwaites(x, np.ones_like(x), re, theta0=0.0)
    xs = x[100:]
    rex = re * xs
    assert np.allclose(bl.theta[100:], 0.664 * xs / np.sqrt(rex), rtol=0.015)
    assert np.allclose(bl.cf[100:], 0.664 / np.sqrt(rex), rtol=0.015)
    # Thwaites H(0) = 2.61 (Blasius 2.59); the two correlation branches differ by 1e-4 at 0
    assert np.allclose(bl.H[100:], 2.61, atol=1e-3)
    assert bl.i_sep is None


def test_thwaites_howarth_retarded_flow_separation():
    # Howarth's linearly retarded flow Ue = 1 - x: exact separation at x = 0.1199;
    # Thwaites' method predicts x ~ 0.123 (White Sec. 4-6).
    x = np.linspace(0, 0.2, 20001)
    bl = thwaites(x, 1 - x, 1e5)
    assert bl.s_sep == pytest.approx(0.1199, rel=0.03)


def test_thwaites_stagnation_point():
    # Ue = K s: theta^2 = 0.075 nu / K everywhere (Hiemenz flow, constant thickness).
    re, K = 1e6, 2.0
    s = np.linspace(0, 0.05, 501)
    bl = thwaites(s, K * s, re)
    assert np.allclose(bl.theta, np.sqrt(0.075 / (re * K)), rtol=1e-3)
    assert np.allclose(bl.lam, 0.075, rtol=1e-3)


def test_michel_criterion_value():
    # 1.174 (1 + 22400/1e6) (1e6)^0.46
    assert michel_re_theta(1e6) == pytest.approx(1.174 * 1.0224 * 1e6**0.46)


def test_head_shape_factor_correlation_roundtrip():
    H = np.linspace(1.2, 2.8, 50)
    assert np.allclose(h_of_h1(h1_of_h(H)), H, rtol=1e-10)


def test_head_flat_plate_vs_power_law():
    # Turbulent flat plate from x = 0 (1/7-power law, White Sec. 6-6):
    # theta ~ 0.036 x Re_x^-0.2, cf ~ 0.0576 Re_x^-0.2. Head + Ludwieg-Tillmann should land
    # within ~10% and settle at H ~ 1.3-1.4.
    re = 1e7
    x = np.linspace(1e-4, 1.0, 400)
    tb = head(x, np.ones_like(x), re, theta0=0.036 * x[0] * (re * x[0]) ** -0.2, H0=1.4)
    rex = re * x[-1]
    assert tb.theta[-1] == pytest.approx(0.036 * rex**-0.2, rel=0.1)
    assert tb.cf[-1] == pytest.approx(0.0576 * rex**-0.2, rel=0.1)
    assert 1.25 < tb.H[-1] < 1.45
    assert tb.i_sep is None


def test_symmetric_airfoil_zero_alpha_is_symmetric():
    r = viscous_drag(solve_airfoil("0012", 0.0, 200), 6e6)
    assert r.upper.x_transition == pytest.approx(r.lower.x_transition)
    assert r.upper.squire_young == pytest.approx(r.lower.squire_young, rel=1e-6)


def test_naca0012_profile_drag_magnitude():
    # Measured / XFOIL NACA 0012 profile drag near alpha = 0 at Re ~ 6e6 is roughly
    # 0.005-0.006 (Abbott & von Doenhoff App. IV). This uncoupled Thwaites/Michel/Head/
    # Squire-Young chain is expected to land within ~30% (transition is predicted early).
    r = viscous_drag(solve_airfoil("0012", 0.0, 200), 6e6)
    assert 0.005 < r.cd < 0.0075
    assert 0.7 < r.cd_friction / r.cd < 1.0  # skin friction dominates at low alpha


def test_drag_trends():
    sol = solve_airfoil("0012", 0.0, 200)
    cds = [viscous_drag(sol, re).cd for re in (2e6, 4e6, 8e6)]
    assert cds[0] > cds[1] > cds[2]  # thinner boundary layers at higher Re
    lo = viscous_drag(solve_airfoil("0012", 0.0, 200), 6e6)
    hi = viscous_drag(solve_airfoil("0012", 6 * DEG, 200), 6e6)
    assert hi.upper.x_transition < lo.upper.x_transition  # suction peak moves transition fwd
    assert hi.lower.x_transition > lo.lower.x_transition
    assert hi.cd > lo.cd


def test_forced_transition():
    sol = solve_airfoil("0012", 0.0, 200)
    free = viscous_drag(sol, 3e6)
    tripped = viscous_drag(sol, 3e6, x_trip=0.05)
    assert tripped.upper.transition_cause == "forced"
    assert tripped.upper.x_transition == pytest.approx(0.05, abs=0.01)
    assert tripped.cd > free.cd
