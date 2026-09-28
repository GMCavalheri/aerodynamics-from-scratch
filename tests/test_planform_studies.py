import numpy as np
import pytest

from aero.vortex_lattice_3d import Wing
from aero.vortex_lattice_3d.studies import (
    induced_drag_vs_aspect_ratio,
    sweep_study,
    taper_study,
)

DEG = np.pi / 180


def test_elliptic_minimises_induced_drag():
    # Munk / Prandtl: among planar wings of given span and lift, elliptic loading has the
    # least induced drag, so no other planform beats e = 1.
    others = [Wing.rectangular(8), Wing.tapered(8, 0.4), Wing.tapered(8, 0.0)]
    cdi_e, e_ell = induced_drag_vs_aspect_ratio(Wing.elliptic, [8])
    for w in others:
        cdi, e = induced_drag_vs_aspect_ratio(lambda ar, w=w: w, [8])
        assert e[0] < e_ell[0]
        assert cdi[0] > cdi_e[0]


def test_optimum_taper_ratio():
    # Classical lifting-line result: an unswept tapered wing is closest to elliptic for
    # lambda ~ 0.3-0.4 (Anderson Sec. 5.3.2, Fig. 5.20); the VLM optimum is similar.
    st = taper_study(8, np.linspace(0, 1, 21))
    assert 0.25 <= st.optimum_lifting_line <= 0.45
    assert 0.3 <= st.optimum_vlm <= 0.55
    assert st.e_vlm.max() > 0.99
    assert st.e_vlm[0] < 0.9  # pointed tip overloads the tips
    assert st.e_vlm[-1] < st.e_vlm.max()


def test_induced_drag_inverse_with_aspect_ratio():
    ars = np.array([4, 8, 16])
    cdi, e = induced_drag_vs_aspect_ratio(Wing.elliptic, ars)
    assert np.allclose(cdi * ars, cdi[0] * ars[0], rtol=0.01)
    assert np.allclose(e, 1, atol=0.005)


def test_aft_sweep_moves_loading_outboard_and_cuts_lift_slope():
    sweeps = np.array([0, 30, 45]) * DEG
    sols = sweep_study(8, 0.4, sweeps)
    tip = [s.cl_span[5] / s.CL for s in sols]  # an outboard strip (eta ~ -0.86)
    cla = [s.CL for s in sols]
    assert tip[0] < tip[1] < tip[2]
    assert cla[0] > cla[1] > cla[2]
    # CL_alpha roughly follows the simple-sweep-theory trend (cos Lambda scaling of the
    # normal component) better than being independent of sweep
    assert cla[2] / cla[0] == pytest.approx(np.cos(45 * DEG), abs=0.15)


def test_forward_sweep_moves_loading_inboard():
    fwd, straight = sweep_study(8, 0.4, [-30 * DEG, 0.0])
    assert fwd.cl_span[5] / fwd.CL < straight.cl_span[5] / straight.CL
