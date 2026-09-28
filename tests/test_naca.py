import numpy as np
import pytest

from aero.geometry import NACA4


def test_parse_designation():
    af = NACA4.from_designation("NACA 2412")
    assert (af.m, af.p, af.t) == pytest.approx((0.02, 0.4, 0.12))
    assert af.name == "NACA 2412"
    with pytest.raises(ValueError):
        NACA4.from_designation("23012")


@pytest.mark.parametrize("closed_te", [True, False])
def test_max_thickness_is_t(closed_te):
    # Abbott & von Doenhoff: the thickness polynomial peaks at 2 * y_t = t at x ~ 0.30.
    af = NACA4.from_designation("0012", closed_te=closed_te)
    x = np.linspace(0, 1, 20001)
    full = 2 * af.thickness(x)
    assert full.max() == pytest.approx(0.12, rel=2e-3)
    assert x[full.argmax()] == pytest.approx(0.30, abs=0.01)


def test_trailing_edge_thickness():
    assert NACA4(0, 0, 0.12, closed_te=True).thickness(1.0) == pytest.approx(0.0, abs=1e-12)
    assert NACA4(0, 0, 0.12, closed_te=False).thickness(1.0) > 1e-3


def test_camber_peak_and_continuity():
    af = NACA4.from_designation("4412")
    x = np.linspace(0, 1, 10001)
    yc = af.camber(x)
    assert yc.max() == pytest.approx(0.04, rel=1e-6)
    assert x[yc.argmax()] == pytest.approx(0.4, abs=1e-3)
    assert af.camber(0.0) == pytest.approx(0.0) and af.camber(1.0) == pytest.approx(0.0)
    # slope is continuous and zero at x = p
    assert af.camber_slope(0.4) == pytest.approx(0.0)
    # slope matches finite differences of the camber line
    h = 1e-6
    xs = np.array([0.1, 0.3, 0.6, 0.9])
    fd = (af.camber(xs + h) - af.camber(xs - h)) / (2 * h)
    assert np.allclose(fd, af.camber_slope(xs), atol=1e-6)


def test_coordinates_ordering_and_symmetry():
    af = NACA4.from_designation("0012")
    pts = af.coordinates(100)
    assert pts.shape == (101, 2)
    assert np.allclose(pts[0], [1, 0]) and np.allclose(pts[-1], [1, 0])
    le = pts[50]
    assert np.allclose(le, [0, 0])
    # clockwise: first half is the lower surface
    assert np.all(pts[1:50, 1] < 0) and np.all(pts[51:-1, 1] > 0)
    # symmetric section: upper surface mirrors the lower one
    assert np.allclose(pts[::-1, 0], pts[:, 0])
    assert np.allclose(pts[::-1, 1], -pts[:, 1])
    # clockwise contour has negative signed area (shoelace)
    x, y = pts[:, 0], pts[:, 1]
    area = 0.5 * np.sum(x[:-1] * y[1:] - x[1:] * y[:-1])
    assert area < 0
    # NACA 00xx area ~ 0.685 t (Abbott & von Doenhoff)
    assert -area == pytest.approx(0.685 * 0.12, rel=0.01)


def test_coordinates_rejects_odd_panel_count():
    with pytest.raises(ValueError):
        NACA4.from_designation("0012").coordinates(99)
