import numpy as np
import pytest
from scipy.integrate import quad

from aero import singularities as s

RNG = np.random.default_rng(0)
PTS = RNG.uniform(-3, 3, size=(2, 50))


def test_point_source_is_radial_with_1_over_2pi_r():
    x, y = PTS
    u, v = s.source_velocity(2.0, 0.3, -0.2, x, y)
    dx, dy = x - 0.3, y + 0.2
    r = np.hypot(dx, dy)
    assert np.allclose(np.hypot(u, v), 2.0 / (2 * np.pi * r))
    assert np.allclose(u * dy - v * dx, 0)  # no tangential component
    assert np.all(u * dx + v * dy > 0)  # outflow


def test_point_vortex_is_tangential_ccw():
    x, y = PTS
    u, v = s.vortex_velocity(3.0, 0.0, 0.0, x, y)
    r = np.hypot(x, y)
    assert np.allclose(np.hypot(u, v), 3.0 / (2 * np.pi * r))
    assert np.allclose(u * x + v * y, 0)
    assert np.all(x * v - y * u > 0)  # counter-clockwise


def test_vortex_circulation_around_loop():
    th = np.linspace(0, 2 * np.pi, 2001)
    x, y = 0.2 + 1.5 * np.cos(th), -0.1 + 1.5 * np.sin(th)
    u, v = s.vortex_velocity(1.7, 0.2, -0.1, x, y)
    circ = np.trapezoid(u * np.gradient(x) + v * np.gradient(y))
    assert circ == pytest.approx(1.7, rel=1e-4)


def test_doublet_velocity_is_gradient_of_potential():
    x, y = PTS
    h = 1e-6
    a = 0.7
    phi = lambda xx, yy: s.doublet_potential(1.3, 0.1, 0.2, xx, yy, a)  # noqa: E731
    u, v = s.doublet_velocity(1.3, 0.1, 0.2, x, y, a)
    assert np.allclose(u, (phi(x + h, y) - phi(x - h, y)) / (2 * h), rtol=1e-5, atol=1e-8)
    assert np.allclose(v, (phi(x, y + h) - phi(x, y - h)) / (2 * h), rtol=1e-5, atol=1e-8)


def test_doublet_is_limit_of_source_sink_pair():
    # source at -eps/2, sink at +eps/2 on the x axis, mu = sigma * eps
    x, y = PTS
    eps, mu = 1e-5, 2.0
    us, vs = s.source_velocity(mu / eps, -eps / 2, 0, x, y)
    uk, vk = s.source_velocity(-mu / eps, eps / 2, 0, x, y)
    u, v = s.doublet_velocity(mu, 0, 0, x, y)
    assert np.allclose(us + uk, u, rtol=1e-4, atol=1e-6)
    assert np.allclose(vs + vk, v, rtol=1e-4, atol=1e-6)


def test_uniform_flow_plus_doublet_is_cylinder():
    # Anderson, Fundamentals of Aerodynamics, Sec. 3.13: V_theta = -2 U sin(theta) on r = R
    U, R = 1.5, 0.8
    th = np.linspace(0.01, 2 * np.pi - 0.01, 100)
    x, y = R * np.cos(th), R * np.sin(th)
    u, v = s.doublet_velocity(2 * np.pi * U * R**2, 0, 0, x, y)
    u = u + U
    assert np.allclose(u * np.cos(th) + v * np.sin(th), 0, atol=1e-12)  # no penetration
    assert np.allclose(-u * np.sin(th) + v * np.cos(th), -2 * U * np.sin(th))


PANEL = (0.2, -0.1, 1.1, 0.4)  # x1, y1, x2, y2


@pytest.mark.parametrize("kind", ["source", "vortex"])
def test_panel_matches_quadrature_of_point_singularities(kind):
    x1, y1, x2, y2 = PANEL
    L = np.hypot(x2 - x1, y2 - y1)
    point = s.source_velocity if kind == "source" else s.vortex_velocity
    panel = s.source_panel_velocity if kind == "source" else s.vortex_panel_velocity
    for x, y in PTS.T[:15]:

        def comp(i, x=x, y=y):
            def f(t):
                return point(1.0, x1 + t * (x2 - x1), y1 + t * (y2 - y1), x, y)[i] * L

            return quad(f, 0, 1, limit=200)[0]

        u, v = panel(x, y, *PANEL)
        assert u == pytest.approx(comp(0), rel=1e-7, abs=1e-10)
        assert v == pytest.approx(comp(1), rel=1e-7, abs=1e-10)


def test_panel_far_field_is_point_singularity():
    x1, y1, x2, y2 = PANEL
    L = np.hypot(x2 - x1, y2 - y1)
    xm, ym = (x1 + x2) / 2, (y1 + y2) / 2
    x, y = 150.0, -90.0
    assert np.allclose(
        s.source_panel_velocity(x, y, *PANEL), s.source_velocity(L, xm, ym, x, y), rtol=1e-4
    )
    assert np.allclose(
        s.vortex_panel_velocity(x, y, *PANEL), s.vortex_velocity(L, xm, ym, x, y), rtol=1e-4
    )


def test_panel_self_induced_jump():
    # Just off the panel on the normal side: source gives v_n = sigma / 2,
    # a CCW vortex sheet gives u_t = -gamma / 2.
    x1, y1, x2, y2 = PANEL
    tx, ty, nx, ny, _ = s.panel_frame(*PANEL)
    xm, ym = (x1 + x2) / 2, (y1 + y2) / 2
    for on_panel, (x, y) in [(False, (xm + 1e-9 * nx, ym + 1e-9 * ny)), (True, (xm, ym))]:
        mask = np.array(on_panel)
        u, v = s.source_panel_velocity(x, y, *PANEL, on_panel=mask)
        assert u * nx + v * ny == pytest.approx(0.5, abs=1e-6)
        assert u * tx + v * ty == pytest.approx(0.0, abs=1e-6)
        u, v = s.vortex_panel_velocity(x, y, *PANEL, on_panel=mask)
        assert u * tx + v * ty == pytest.approx(-0.5, abs=1e-6)
        assert u * nx + v * ny == pytest.approx(0.0, abs=1e-6)
