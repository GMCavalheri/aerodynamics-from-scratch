"""Panel method vs the cited measurements in validation/reference_data.

Inviscid theory is expected to *over*-predict lift slope and |cm| compared with experiment
(the boundary layer de-cambers the section), so these checks assert the physically expected
direction and a bounded gap, not a match.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "validation"))

from reference import section_value  # noqa: E402

from aero.panel_method_2d import solve_airfoil  # noqa: E402

DEG = np.pi / 180


def summary(code):
    s0, s1 = solve_airfoil(code, 0.0, 300), solve_airfoil(code, 4 * DEG, 300)
    slope = (s1.cl - s0.cl) / (4 * DEG)
    return slope * DEG, -s0.cl / slope / DEG, s1.cm


@pytest.mark.parametrize("code", ["NACA 0012", "NACA 2412"])
def test_zero_lift_angle_within_measured_band(code):
    _, a0, _ = summary(code)
    _, lo, hi, _ = section_value(code, "alpha_L0_deg")
    assert lo - 0.1 <= a0 <= hi + 0.1


@pytest.mark.parametrize("code", ["NACA 0012", "NACA 2412"])
def test_inviscid_lift_slope_exceeds_measured_by_bounded_margin(code):
    slope, _, _ = summary(code)
    measured, _, _, _ = section_value(code, "cl_alpha_per_deg")
    assert 1.0 < slope / measured < 1.2


def test_cambered_moment_sign_and_magnitude():
    _, _, cm = summary("NACA 2412")
    measured, _, _, _ = section_value("NACA 2412", "cm_c4")
    assert cm < measured < 0  # nose-down, inviscid magnitude larger
    assert abs(cm - measured) < 0.02
