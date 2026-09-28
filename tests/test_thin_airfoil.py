import numpy as np
import pytest
from scipy.integrate import quad

from aero.thin_airfoil import ThinAirfoil

DEG = np.pi / 180


def test_symmetric_section():
    # Flat plate / symmetric: alpha_L0 = 0, cl = 2 pi alpha, cm_c/4 = 0, x_cp = c/4.
    ta = ThinAirfoil.from_naca4("0012")
    assert ta.alpha_zero_lift == pytest.approx(0.0, abs=1e-14)
    assert ta.cl(5 * DEG) == pytest.approx(2 * np.pi * 5 * DEG)
    assert ta.cm_quarter_chord == pytest.approx(0.0, abs=1e-14)
    assert ta.x_cp(3 * DEG) == pytest.approx(0.25)
    assert ta.cm_le(4 * DEG) == pytest.approx(-ta.cl(4 * DEG) / 4)


def test_parabolic_arc_closed_form():
    # z = 4h x (1 - x) -> dz/dx = 4h cos(theta): A1 = 4h, An>1 = 0,
    # alpha_L0 = -2h, cm_c/4 = -pi h.
    h = 0.03
    ta = ThinAirfoil(lambda x: 4 * h * (1 - 2 * x))
    assert ta.A1 == pytest.approx(4 * h)
    assert np.allclose(ta.An[1:], 0, atol=1e-12)
    assert ta.alpha_zero_lift == pytest.approx(-2 * h)
    assert ta.cm_quarter_chord == pytest.approx(-np.pi * h)


def test_naca2412_anderson_example():
    # Anderson, Fundamentals of Aerodynamics, Example 4.6/4.7:
    # alpha_L0 = -2.077 deg, cm_c/4 = -0.053 for NACA 2412.
    ta = ThinAirfoil.from_naca4("2412")
    assert ta.alpha_zero_lift / DEG == pytest.approx(-2.077, abs=2e-3)
    assert ta.cm_quarter_chord == pytest.approx(-0.053, abs=5e-4)


def test_results_scale_linearly_with_camber():
    a = ThinAirfoil.from_naca4("2412")
    b = ThinAirfoil.from_naca4("4412")
    assert b.alpha_zero_lift == pytest.approx(2 * a.alpha_zero_lift)
    assert b.cm_quarter_chord == pytest.approx(2 * a.cm_quarter_chord)


def test_loading_integrates_to_cl_and_cm():
    ta = ThinAirfoil.from_naca4("2412")
    alpha = 4 * DEG

    # integrate in theta to handle the LE square-root singularity
    def dcp(th):
        return ta.delta_cp(0.5 * (1 - np.cos(th)), alpha) * 0.5 * np.sin(th)

    cl = quad(dcp, 0, np.pi, limit=200)[0]
    cm_le = -quad(lambda th: dcp(th) * 0.5 * (1 - np.cos(th)), 0, np.pi, limit=200)[0]
    assert cl == pytest.approx(ta.cl(alpha), rel=1e-6)
    assert cm_le == pytest.approx(ta.cm_le(alpha), rel=1e-6)
    assert ta.cm_le(alpha) + ta.cl(alpha) / 4 == pytest.approx(ta.cm_quarter_chord)
