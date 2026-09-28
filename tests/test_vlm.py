import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import quad

from aero.vortex_lattice_3d import Wing, solve_lifting_line, solve_vlm
from aero.vortex_lattice_3d.biot_savart import (
    horseshoe_velocity,
    segment_velocity,
    semi_infinite_velocity,
)

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "validation"))
from reference import load  # noqa: E402

DEG = np.pi / 180


# --- Biot-Savart -----------------------------------------------------------------------------
def test_long_segment_tends_to_infinite_line():
    # Infinite straight vortex: |V| = Gamma / (2 pi h), direction by the right-hand rule.
    a, b = np.array([0.0, -1e5, 0.0]), np.array([0.0, 1e5, 0.0])
    v = segment_velocity(np.array([0.0, 0.3, 0.5]), a, b, gamma=2.0)
    # filament along +y, point above it: dl x r = y x z = +x
    assert np.allclose(v, [2.0 / (2 * np.pi * 0.5), 0, 0], rtol=1e-6)


def test_segment_matches_quadrature():
    a, b = np.array([0.1, -0.4, 0.2]), np.array([0.5, 0.7, -0.1])
    p = np.array([0.9, 0.2, 0.6])
    t = b - a

    def dv(s, k):
        r = p - (a + s * t)
        return np.cross(t, r)[k] / (4 * np.pi * np.linalg.norm(r) ** 3)

    ref = [quad(dv, 0, 1, args=(k,))[0] for k in range(3)]
    assert np.allclose(segment_velocity(p, a, b), ref, rtol=1e-8)


def test_semi_infinite_is_half_infinite_at_its_foot():
    v = semi_infinite_velocity(np.array([0.0, 0.0, 1.0]), np.zeros(3), [1.0, 0, 0])
    assert np.allclose(v, [0, -1 / (4 * np.pi), 0])  # x cross z = -y


def test_horseshoe_on_axis_is_downwash():
    # Behind the bound vortex, between the legs, positive Gamma induces downwash.
    v = horseshoe_velocity(np.array([0.5, 0, 0]), np.array([0, -1.0, 0]), np.array([0, 1.0, 0]))
    assert v[2] < 0 and abs(v[0]) < 1e-12 and abs(v[1]) < 1e-12


def test_core_cutoff_on_filament():
    a, b = np.zeros(3), np.array([0, 1.0, 0])
    assert np.allclose(segment_velocity(np.array([0, 0.5, 0]), a, b), 0)


# --- lifting line ----------------------------------------------------------------------------
@pytest.mark.parametrize("ar", [4, 8, 20])
def test_lifting_line_elliptic_is_exact(ar):
    ll = solve_lifting_line(Wing.elliptic(ar), 5 * DEG)
    assert ll.CL / (5 * DEG) == pytest.approx(2 * np.pi / (1 + 2 / ar), rel=1e-6)
    assert ll.span_efficiency == pytest.approx(1.0, abs=1e-8)


# --- VLM -------------------------------------------------------------------------------------
@pytest.mark.parametrize(("ar", "shape_tol"), [(4, 0.06), (8, 0.04), (20, 0.02)])
def test_elliptic_wing_has_elliptic_loading(ar, shape_tol):
    # Prandtl: elliptic planform -> elliptic loading -> CDi = CL^2 / (pi AR), e = 1.
    # Lifting-surface effects make low-AR loading only nearly elliptic (the deviation is
    # mostly a small sin(3 theta) term, which barely changes e).
    s = solve_vlm(Wing.elliptic(ar), 5 * DEG, 40, 4)
    assert s.span_efficiency == pytest.approx(1.0, abs=5e-3)
    assert s.CDi == pytest.approx(s.CL**2 / (np.pi * ar), rel=5e-3)
    eta = 2 * s.y_mid / s.wing.span
    g = s.strip_gamma
    assert np.allclose(g / g.max(), np.sqrt(1 - eta**2), atol=shape_tol)


def test_high_aspect_ratio_matches_lifting_line():
    # Lifting-line theory is the AR -> infinity limit of lifting-surface theory.
    ar = 100
    s = solve_vlm(Wing.elliptic(ar), 5 * DEG, 40, 2)
    assert s.CL / (5 * DEG) == pytest.approx(2 * np.pi / (1 + 2 / ar), rel=5e-3)


def test_lift_slope_below_lifting_line_and_helmholtz_bracket():
    # Lifting-surface CL_alpha of a finite wing lies below lifting-line's 2 pi / (1 + 2/AR);
    # Helmholtz' 2 pi AR / (2 + sqrt(AR^2 + 4)) is a closer estimate (within ~3% for AR 8).
    ar = 8
    cla = solve_vlm(Wing.elliptic(ar), 5 * DEG, 40, 4).CL / (5 * DEG)
    assert cla < 2 * np.pi / (1 + 2 / ar)
    assert cla == pytest.approx(2 * np.pi * ar / (2 + np.sqrt(ar**2 + 4)), rel=0.03)


def test_benchmarks_from_reference_data():
    for row in load("vlm_benchmarks.csv"):
        ar, sweep = float(row["aspect_ratio"]), float(row["sweep_c4_deg"]) * DEG
        if row["taper"] == "elliptic":
            wing = Wing.elliptic(ar)
            s = solve_vlm(wing, 2 * DEG, 40, 20)
        else:
            wing = Wing.tapered(ar, float(row["taper"]), sweep=sweep)
            s = solve_vlm(wing, 2 * DEG, 8, 1)  # Bertin: 4 horseshoes per semispan
        cla = s.CL / (2 * DEG)
        assert cla == pytest.approx(float(row["cl_alpha_per_rad"]), rel=float(row["tolerance_rel"]))


def test_symmetry_and_linearity():
    w = Wing.tapered(6, 0.5, sweep=30 * DEG, dihedral=5 * DEG)
    s1, s2 = solve_vlm(w, 2 * DEG), solve_vlm(w, 4 * DEG)
    assert np.allclose(s1.strip_gamma, s1.strip_gamma[::-1], rtol=1e-9)
    assert s2.CL == pytest.approx(2 * s1.CL, rel=2e-3)  # sin(alpha) vs alpha
    assert solve_vlm(w, 0.0).CL == pytest.approx(0.0, abs=1e-12)


def test_mesh_convergence():
    w = Wing.tapered(8, 0.4)
    cl = [solve_vlm(w, 5 * DEG, ns, 4).CL for ns in (10, 20, 40, 80)]
    d = np.abs(np.diff(cl))
    assert np.all(d[1:] < d[:-1])
    assert d[-1] / cl[-1] < 1e-3


def test_washout_unloads_the_tips():
    base = Wing.rectangular(8)
    washed = Wing(base.span, base.chord, twist=lambda eta: -4 * DEG * eta)
    s0, s1 = solve_vlm(base, 5 * DEG), solve_vlm(washed, 5 * DEG)
    outboard = np.abs(2 * s0.y_mid / base.span) > 0.8
    assert np.all(s1.cl_span[outboard] < s0.cl_span[outboard])
    assert s1.CL < s0.CL
    assert solve_vlm(washed, 0.0).CL < 0


def test_sweep_reduces_lift_slope():
    cl = [solve_vlm(Wing.rectangular(6, sweep=sw * DEG), 4 * DEG).CL for sw in (0, 30, 45)]
    assert cl[0] > cl[1] > cl[2]


def test_trefftz_and_near_field_lift_agree():
    s = solve_vlm(Wing.tapered(8, 0.4), 5 * DEG, 40, 4)
    cl_trefftz = 2 * np.sum(s.strip_gamma * np.diff(s.lattice.y_nodes)) / s.area
    assert s.CL == pytest.approx(cl_trefftz, rel=5e-3)
