"""Vortex lattice method (VLM) for finite wings.

Each panel carries a horseshoe vortex whose bound segment lies on the panel quarter-chord
line and whose trailing legs run to infinity downstream (+x). Enforcing flow tangency at the
three-quarter-chord control points gives a dense linear system for the circulations
(Katz & Plotkin, *Low-Speed Aerodynamics*, Sec. 12.3; Bertin & Cummings, Ch. 7).

Lift and moment come from the Kutta-Joukowski theorem on each bound segment. Induced drag
is computed in the Trefftz plane far downstream, where the wake reduces to 2D point
vortices at the spanwise panel edges:

    D_i = -(rho / 2) * sum_j Gamma_j (w_j . n_j) ds_j

which is much less sensitive to discretisation than near-field force integration.

Nondimensional: V_inf = 1, rho = 1.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from aero.vortex_lattice_3d.biot_savart import horseshoe_velocity
from aero.vortex_lattice_3d.geometry import Lattice, Wing, build_lattice


def influence_matrix(lat: Lattice) -> np.ndarray:
    """Normal velocity at every control point due to a unit horseshoe on every panel."""
    cp = lat.control.reshape(-1, 3)
    a = lat.bound_a.reshape(-1, 3)
    b = lat.bound_b.reshape(-1, 3)
    n = lat.normal.reshape(-1, 3)
    vel = horseshoe_velocity(cp[:, None, :], a[None], b[None])
    return np.einsum("ijk,ik->ij", vel, n)


def freestream(alpha: float) -> np.ndarray:
    return np.array([np.cos(alpha), 0.0, np.sin(alpha)])


@dataclass
class VLMSolution:
    wing: Wing
    lattice: Lattice
    alpha: float
    gamma: np.ndarray  # (nc, ns) horseshoe strengths
    x_ref: float = 0.0  # moment reference x (e.g. root quarter chord = 0)

    @property
    def strip_gamma(self) -> np.ndarray:
        """Total bound circulation of each spanwise strip (what the wake sees)."""
        return self.gamma.sum(axis=0)

    @property
    def y_mid(self) -> np.ndarray:
        """Spanwise stations of the strips (where the control points sit)."""
        y = self.lattice.y_nodes
        return y[:-1] + self.lattice.span_frac * np.diff(y)

    @property
    def chord_mid(self) -> np.ndarray:
        eta = np.abs(self.y_mid) / (self.wing.span / 2)
        return self.wing.chord(eta)

    @property
    def area(self) -> float:
        """Projected area of the lattice (used as reference area)."""
        return float(self.lattice.area.sum())

    def _panel_forces(self) -> np.ndarray:
        v = freestream(self.alpha)
        seg = self.lattice.bound_b - self.lattice.bound_a
        return self.gamma[..., None] * np.cross(v, seg)  # rho = V = 1

    @property
    def force_coefficients(self) -> np.ndarray:
        """Body-axis force coefficients (CX, CY, CZ) from Kutta-Joukowski."""
        return 2 * self._panel_forces().sum(axis=(0, 1)) / self.area

    @property
    def CL(self) -> float:
        cx, _, cz = self.force_coefficients
        return float(cz * np.cos(self.alpha) - cx * np.sin(self.alpha))

    @property
    def cl_span(self) -> np.ndarray:
        """Section lift coefficient at strip mid-points, cl = 2 Gamma / (V c)."""
        return 2 * self.strip_gamma / np.maximum(self.chord_mid, 1e-12)

    @property
    def span_loading(self) -> np.ndarray:
        """c * cl / c_avg: normalised spanwise loading (area under it / b equals CL)."""
        c_avg = self.area / self.wing.span
        return 2 * self.strip_gamma / c_avg

    @property
    def Cm(self) -> float:
        """Pitching moment about (x_ref, 0, 0), nose-up positive, referenced to the MAC."""
        f = self._panel_forces()
        mid = 0.5 * (self.lattice.bound_a + self.lattice.bound_b)
        r = mid - np.array([self.x_ref, 0.0, 0.0])
        m = np.cross(r, f).sum(axis=(0, 1))
        return float(2 * m[1] / (self.area * self.wing.mean_aerodynamic_chord))

    def trefftz(self) -> tuple[np.ndarray, np.ndarray]:
        """Normal wash w_n and segment lengths ds in the Trefftz plane (per strip)."""
        te = self.lattice.corners[-1]  # trailing-edge node row, (ns + 1, 3)
        yz = te[:, 1:]
        g = self.strip_gamma
        gs = np.concatenate([[0.0], g, [0.0]])
        shed = gs[:-1] - gs[1:]  # (ns + 1,) net +x filament strength at each node
        seg = np.diff(yz, axis=0)
        ds = np.linalg.norm(seg, axis=1)
        n = np.column_stack([-seg[:, 1], seg[:, 0]]) / ds[:, None]  # up for +y segments
        mid = yz[:-1] + self.lattice.span_frac[:, None] * seg
        d = mid[:, None, :] - yz[None, :, :]
        r2 = np.sum(d**2, axis=-1)
        vel = shed[None, :, None] * np.stack([-d[..., 1], d[..., 0]], axis=-1)
        vel = (vel / (2 * np.pi * r2[..., None])).sum(axis=1)
        return np.sum(vel * n, axis=1), ds

    @property
    def CDi(self) -> float:
        wn, ds = self.trefftz()
        return float(-np.sum(self.strip_gamma * wn * ds) / self.area)

    @property
    def span_efficiency(self) -> float:
        return float(self.CL**2 / (np.pi * self.wing.aspect_ratio * self.CDi))


def solve_vlm(
    wing: Wing,
    alpha: float,
    n_span: int = 40,
    n_chord: int = 6,
    spacing: str = "cosine",
    x_ref: float = 0.0,
) -> VLMSolution:
    """Solve the VLM for ``wing`` at angle of attack ``alpha`` [rad]."""
    lat = build_lattice(wing, n_span, n_chord, spacing)
    aic = influence_matrix(lat)
    rhs = -lat.normal.reshape(-1, 3) @ freestream(alpha)
    gamma = np.linalg.solve(aic, rhs).reshape(lat.control.shape[:2])
    return VLMSolution(wing, lat, alpha, gamma, x_ref)
