"""Joukowski airfoils and their exact potential-flow solution.

The map z = zeta + b^2 / zeta sends a circle of radius a centred at zeta0 = (-mx, my),
passing through zeta = b, onto a cusped, cambered airfoil with its trailing edge at z = 2b.
The flow around the circle with the Kutta condition at zeta = b is known in closed form,
so this gives an exact benchmark for panel methods (Anderson, *Fundamentals of
Aerodynamics*, Sec. 4.14; Katz & Plotkin, Sec. 6.6).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Joukowski:
    """Joukowski airfoil from thickness (``mx``) and camber (``my``) offsets of the circle."""

    mx: float = 0.1
    my: float = 0.0
    b: float = 1.0

    @property
    def zeta0(self) -> complex:
        return complex(-self.mx, self.my)

    @property
    def radius(self) -> float:
        return abs(self.b - self.zeta0)

    @property
    def beta(self) -> float:
        """Angle of the trailing-edge point seen from the circle centre (camber angle)."""
        return float(np.arctan2(self.my, self.b + self.mx))

    def _theta(self, s: np.ndarray) -> np.ndarray:
        # s in [0, 2 pi] runs clockwise from the trailing edge (lower surface first)
        return -self.beta - np.asarray(s, dtype=float)

    def _circle(self, s):
        return self.zeta0 + self.radius * np.exp(1j * self._theta(s))

    def _map(self, zeta):
        return zeta + self.b**2 / zeta

    @property
    def chord(self) -> float:
        z = self._map(self._circle(np.linspace(0, 2 * np.pi, 20001)))
        return float(z.real.max() - z.real.min())

    @property
    def x_le(self) -> float:
        z = self._map(self._circle(np.linspace(0, 2 * np.pi, 20001)))
        return float(z.real.min())

    def coordinates(self, n_panels: int = 200) -> np.ndarray:
        """Contour normalised to unit chord with LE at x = 0, clockwise from the TE.

        Nodes are equally spaced in the circle angle, which clusters them at LE and TE.
        """
        s = np.linspace(0, 2 * np.pi, n_panels + 1)
        z = self._map(self._circle(s))
        z = (z - self.x_le) / self.chord
        pts = np.column_stack([z.real, z.imag])
        pts[-1] = pts[0]
        return pts

    def circulation(self, alpha: float, v_inf: float = 1.0) -> float:
        """Clockwise circulation from the Kutta condition, 4 pi a V sin(alpha + beta)."""
        return 4 * np.pi * self.radius * v_inf * np.sin(alpha + self.beta)

    def cl(self, alpha: float) -> float:
        return 2 * self.circulation(alpha) / self.chord

    def surface_cp(self, points: np.ndarray, alpha: float) -> np.ndarray:
        """Exact Cp at airfoil surface points given in normalised coordinates.

        Each point is mapped back to the circle by inverting z = zeta + b^2/zeta and picking
        the root outside the ``b``-circle, then the circle-plane velocity is transformed.
        """
        z = (points[:, 0] + 1j * points[:, 1]) * self.chord + self.x_le
        disc = np.sqrt(z**2 - 4 * self.b**2)
        r1, r2 = (z + disc) / 2, (z - disc) / 2
        zeta = np.where(np.abs(r1 - self.zeta0) >= np.abs(r2 - self.zeta0), r1, r2)
        # project onto the circle (the points may be panel midpoints slightly inside it)
        zp = zeta - self.zeta0
        zp = self.radius * zp / np.abs(zp)
        zeta = self.zeta0 + zp
        a, g = self.radius, self.circulation(alpha)
        dfdzeta = (
            np.exp(-1j * alpha) - a**2 * np.exp(1j * alpha) / zp**2 + 1j * g / (2 * np.pi * zp)
        )
        dzdzeta = 1 - self.b**2 / zeta**2
        q = np.abs(dfdzeta / dzdzeta)
        return 1 - q**2
