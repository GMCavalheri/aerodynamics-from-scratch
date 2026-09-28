"""Viscous drag of an airfoil from the inviscid panel solution (no interaction).

1. Locate the stagnation point where the panel tangential velocity changes sign and split
   the surface into upper and lower boundary layers, each starting there with Ue = 0.
2. March Thwaites' method until Michel's criterion (or laminar separation, which is assumed
   to close as a short bubble) triggers transition.
3. Continue with Head's turbulent method to the trailing edge.
4. Apply the Squire-Young formula at the trailing edge,
   cd = 2 theta_TE Ue_TE^((H_TE + 5) / 2), summed over both surfaces
   (Squire & Young 1937; Moran Sec. 7.9).

The trailing edge of a finite-angle section is an inviscid stagnation point, so the last
few percent of chord see an unphysical deceleration. Following common practice for
uncoupled methods the march stops at ``x_stop`` (default 0.98 c).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from aero.boundary_layer.head import TurbulentBL, head
from aero.boundary_layer.thwaites import LaminarBL, thwaites
from aero.boundary_layer.transition import michel_transition
from aero.panel_method_2d import PanelSolution


@dataclass
class SurfaceBL:
    """Boundary layer on one surface, stations ordered from the stagnation point."""

    name: str
    x: np.ndarray
    s: np.ndarray
    ue: np.ndarray
    laminar: LaminarBL
    turbulent: TurbulentBL | None
    i_tr: int  # index where the turbulent part starts (len(s) if laminar throughout)
    transition_cause: str  # "michel", "laminar separation", "forced" or "none"

    @property
    def theta(self) -> np.ndarray:
        th = self.laminar.theta.copy()
        if self.turbulent is not None:
            th[self.i_tr :] = self.turbulent.theta
        return th

    @property
    def H(self) -> np.ndarray:
        H = self.laminar.H.copy()
        if self.turbulent is not None:
            H[self.i_tr :] = self.turbulent.H
        return H

    @property
    def cf(self) -> np.ndarray:
        cf = self.laminar.cf.copy()
        if self.turbulent is not None:
            cf[self.i_tr :] = self.turbulent.cf
        return cf

    @property
    def x_transition(self) -> float | None:
        return float(self.x[self.i_tr]) if self.i_tr < len(self.x) else None

    @property
    def x_separation(self) -> float | None:
        """Turbulent separation location (H > 2.4), if it occurs before the TE."""
        if self.turbulent is None or self.turbulent.i_sep is None:
            return None
        return float(self.x[self.i_tr + self.turbulent.i_sep])

    @property
    def squire_young(self) -> float:
        th, H, ue = self.theta[-1], self.H[-1], self.ue[-1]
        return float(2 * th * ue ** ((H + 5) / 2))

    @property
    def friction_drag(self) -> float:
        # integral of tau_w / (0.5 rho V^2 c) along the surface, with tau_w = cf 0.5 rho Ue^2
        return float(np.trapezoid(self.cf * self.ue**2, self.s))


@dataclass
class ViscousResult:
    inviscid: PanelSolution
    re: float
    upper: SurfaceBL
    lower: SurfaceBL

    @property
    def cd(self) -> float:
        """Profile drag coefficient (Squire-Young, both surfaces)."""
        return self.upper.squire_young + self.lower.squire_young

    @property
    def cd_friction(self) -> float:
        return self.upper.friction_drag + self.lower.friction_drag

    @property
    def cd_pressure(self) -> float:
        return self.cd - self.cd_friction


def split_at_stagnation(sol: PanelSolution):
    """Return ((x, s, ue) upper, (x, s, ue) lower) starting at the stagnation point."""
    vt = sol.vt
    xc = sol.control_points
    # panels run TE -> lower -> LE -> upper -> TE; vt < 0 on the lower surface
    k = int(np.nonzero((vt[:-1] < 0) & (vt[1:] >= 0))[0][0])
    f = -vt[k] / (vt[k + 1] - vt[k])
    stag = xc[k] + f * (xc[k + 1] - xc[k])

    def march(points, speeds):
        pts = np.vstack([stag, points])
        s = np.concatenate([[0.0], np.cumsum(np.hypot(*np.diff(pts, axis=0).T))])
        return pts[:, 0], s, np.concatenate([[0.0], np.abs(speeds)])

    upper = march(xc[k + 1 :], vt[k + 1 :])
    lower = march(xc[k::-1], vt[k::-1])
    return upper, lower


def _surface(name, x, s, ue, re, x_stop, x_trip):
    keep = x <= x_stop
    keep[: np.argmin(x) + 1] = True  # never cut the nose region
    x, s, ue = x[keep], s[keep], ue[keep]
    lam = thwaites(s, ue, re)
    i_tr = michel_transition(s, lam.re_theta, re)
    cause = "michel"
    if lam.i_sep is not None and (i_tr is None or lam.i_sep < i_tr):
        i_tr, cause = lam.i_sep, "laminar separation"
    if x_trip is not None:
        after_le = np.arange(len(x)) > np.argmin(x)
        cand = np.nonzero(after_le & (x >= x_trip))[0]
        if len(cand) and (i_tr is None or cand[0] < i_tr):
            i_tr, cause = int(cand[0]), "forced"
    if i_tr is None:
        return SurfaceBL(name, x, s, ue, lam, None, len(s), "none")
    i_tr = max(i_tr, 1)
    turb = head(s[i_tr:], ue[i_tr:], re, theta0=float(lam.theta[i_tr]))
    return SurfaceBL(name, x, s, ue, lam, turb, i_tr, cause)


def viscous_drag(
    sol: PanelSolution, re: float, x_stop: float = 0.98, x_trip: float | None = None
) -> ViscousResult:
    """Boundary layers on both surfaces and the Squire-Young profile drag.

    ``x_trip`` forces transition at that chord fraction (like a trip strip) if free
    transition has not already happened upstream.
    """
    (xu, su, uu), (xl, sl, ul) = split_at_stagnation(sol)
    upper = _surface("upper", xu, su, uu, re, x_stop, x_trip)
    lower = _surface("lower", xl, sl, ul, re, x_stop, x_trip)
    return ViscousResult(sol, re, upper, lower)
