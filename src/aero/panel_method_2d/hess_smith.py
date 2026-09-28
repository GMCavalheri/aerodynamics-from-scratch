"""Hess-Smith panel method for 2D lifting bodies.

The body contour is split into N straight panels, each carrying its own constant source
strength sigma_j, plus a single vortex strength gamma shared by every panel. The N + 1
unknowns follow from

* no penetration at every panel midpoint (N equations), and
* the Kutta condition: equal and opposite tangential velocity on the two trailing-edge
  panels, i.e. the flow leaves the trailing edge smoothly (1 equation).

References: Hess & Smith (1967), Prog. Aero. Sci. 8; Katz & Plotkin, *Low-Speed
Aerodynamics*, Sec. 11.3; Moran, *An Introduction to Theoretical and Computational
Aerodynamics*, Sec. 4.9.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from aero.geometry import NACA4
from aero.singularities import panel_frame, source_panel_velocity, vortex_panel_velocity


@dataclass
class PanelSolution:
    """Solution of the panel method at one angle of attack (chord = 1, V_inf = 1)."""

    nodes: np.ndarray  # (N+1, 2) panel end points, clockwise from the trailing edge
    alpha: float  # angle of attack [rad]
    sigma: np.ndarray  # (N,) source strengths
    gamma: float  # vortex strength per unit length (CCW positive)
    vt: np.ndarray  # (N,) tangential velocity at control points (along the panel direction)
    x_ref: float = 0.25  # moment reference point on the chord

    @property
    def n_panels(self) -> int:
        return len(self.sigma)

    @property
    def control_points(self) -> np.ndarray:
        return 0.5 * (self.nodes[:-1] + self.nodes[1:])

    @property
    def lengths(self) -> np.ndarray:
        return np.hypot(*np.diff(self.nodes, axis=0).T)

    @property
    def normals(self) -> np.ndarray:
        n = self.nodes
        _, _, nx, ny, _ = panel_frame(n[:-1, 0], n[:-1, 1], n[1:, 0], n[1:, 1])
        return np.column_stack([nx, ny])

    @property
    def chord(self) -> float:
        x = self.nodes[:, 0]
        return float(x.max() - x.min())

    @property
    def cp(self) -> np.ndarray:
        """Pressure coefficient at the control points."""
        return 1.0 - self.vt**2

    @property
    def circulation(self) -> float:
        """Total clockwise circulation (positive for positive lift)."""
        return float(-self.gamma * self.lengths.sum())

    @property
    def cl(self) -> float:
        """Lift coefficient from the Kutta-Joukowski theorem, cl = 2 Gamma / (V c)."""
        return 2.0 * self.circulation / self.chord

    def _force_coefficients(self) -> tuple[float, float, float]:
        n, ds, cp = self.normals, self.lengths, self.cp
        fx = -np.sum(cp * n[:, 0] * ds) / self.chord
        fy = -np.sum(cp * n[:, 1] * ds) / self.chord
        r = self.control_points - np.array([self.x_ref * self.chord, 0.0])
        mz = np.sum(r[:, 0] * (-cp * n[:, 1] * ds) - r[:, 1] * (-cp * n[:, 0] * ds))
        return fx, fy, -mz / self.chord**2  # nose-up positive = clockwise

    @property
    def cl_pressure(self) -> float:
        """Lift coefficient from integrating the surface pressure."""
        fx, fy, _ = self._force_coefficients()
        return fy * np.cos(self.alpha) - fx * np.sin(self.alpha)

    @property
    def cd_pressure(self) -> float:
        """Pressure drag; zero in exact potential flow (d'Alembert), so a discretisation check."""
        fx, fy, _ = self._force_coefficients()
        return fx * np.cos(self.alpha) + fy * np.sin(self.alpha)

    @property
    def cm(self) -> float:
        """Pitching moment coefficient about (x_ref, 0), nose-up positive."""
        return self._force_coefficients()[2]

    def velocity(self, x, y) -> tuple[np.ndarray, np.ndarray]:
        """Velocity (u, v) anywhere in the field, e.g. for streamlines."""
        x = np.asarray(x, dtype=float)[..., None]
        y = np.asarray(y, dtype=float)[..., None]
        n = self.nodes
        seg = (n[:-1, 0], n[:-1, 1], n[1:, 0], n[1:, 1])
        us, vs = source_panel_velocity(x, y, *seg)
        uv, vv = vortex_panel_velocity(x, y, *seg)
        u = np.cos(self.alpha) + us @ self.sigma + self.gamma * uv.sum(-1)
        v = np.sin(self.alpha) + vs @ self.sigma + self.gamma * vv.sum(-1)
        return u, v

    def surfaces(self) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Cp split at the leading-edge node: (x_lower, cp_lower, x_upper, cp_upper)."""
        i_le = int(np.argmin(self.nodes[:, 0]))
        xc = self.control_points[:, 0]
        return xc[:i_le], self.cp[:i_le], xc[i_le:], self.cp[i_le:]


def _influence_matrices(nodes: np.ndarray):
    x1, y1 = nodes[:-1, 0], nodes[:-1, 1]
    x2, y2 = nodes[1:, 0], nodes[1:, 1]
    xc, yc = 0.5 * (x1 + x2), 0.5 * (y1 + y2)
    tx, ty, nx, ny, _ = panel_frame(x1, y1, x2, y2)
    on_panel = np.eye(len(xc), dtype=bool)
    args = (xc[:, None], yc[:, None], x1, y1, x2, y2, on_panel)
    us, vs = source_panel_velocity(*args)
    uv, vv = vortex_panel_velocity(*args)
    # rows: control point i, columns: panel j
    a_n = us * nx[:, None] + vs * ny[:, None]
    a_t = us * tx[:, None] + vs * ty[:, None]
    b_n = (uv * nx[:, None] + vv * ny[:, None]).sum(axis=1)
    b_t = (uv * tx[:, None] + vv * ty[:, None]).sum(axis=1)
    return a_n, a_t, b_n, b_t, (tx, ty, nx, ny)


def solve_panels(
    nodes: np.ndarray, alpha: float, lifting: bool = True, x_ref: float = 0.25
) -> PanelSolution:
    """Solve the Hess-Smith system for a closed contour at angle of attack ``alpha`` [rad].

    ``nodes`` must run clockwise from the trailing edge (lower surface first), so that each
    panel normal (tangent rotated +90 deg) points into the fluid. ``lifting=False`` drops the
    vortex and the Kutta condition, giving the pure source solution (e.g. for a cylinder).
    """
    nodes = np.asarray(nodes, dtype=float)
    n = len(nodes) - 1
    a_n, a_t, b_n, b_t, (tx, ty, nx, ny) = _influence_matrices(nodes)
    ca, sa = np.cos(alpha), np.sin(alpha)
    vn_inf = ca * nx + sa * ny
    vt_inf = ca * tx + sa * ty

    if lifting:
        m = np.zeros((n + 1, n + 1))
        m[:n, :n] = a_n
        m[:n, n] = b_n
        m[n, :n] = a_t[0] + a_t[-1]
        m[n, n] = b_t[0] + b_t[-1]
        rhs = np.concatenate([-vn_inf, [-(vt_inf[0] + vt_inf[-1])]])
        sol = np.linalg.solve(m, rhs)
        sigma, gamma = sol[:n], float(sol[n])
    else:
        sigma, gamma = np.linalg.solve(a_n, -vn_inf), 0.0

    vt = a_t @ sigma + gamma * b_t + vt_inf
    return PanelSolution(nodes, float(alpha), sigma, gamma, vt, x_ref)


def solve_airfoil(
    airfoil: NACA4 | str, alpha: float, n_panels: int = 200, x_ref: float = 0.25
) -> PanelSolution:
    """Convenience wrapper: panel a NACA 4-digit section and solve at ``alpha`` [rad]."""
    if isinstance(airfoil, str):
        airfoil = NACA4.from_designation(airfoil)
    return solve_panels(airfoil.coordinates(n_panels), alpha, x_ref=x_ref)
