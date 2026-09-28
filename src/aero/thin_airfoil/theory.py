"""Thin airfoil theory for an arbitrary mean camber line.

The vortex sheet on the camber line is written as a Glauert series in the transformed
chordwise coordinate x = (1 - cos theta) / 2:

    gamma(theta) / (2 V) = A0 (1 + cos theta) / sin theta + sum_n An sin(n theta)

    A0 = alpha - 1/pi   int_0^pi dz/dx dtheta
    An =        2/pi    int_0^pi dz/dx cos(n theta) dtheta

giving cl = pi (2 A0 + A1) = 2 pi (alpha - alpha_L0) and cm_c/4 = pi/4 (A2 - A1)
(Anderson, *Fundamentals of Aerodynamics*, Secs. 4.7-4.8).
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
from scipy.integrate import quad

from aero.geometry import NACA4


def _x_of_theta(theta):
    return 0.5 * (1.0 - np.cos(theta))


class ThinAirfoil:
    """Thin airfoil theory solution for a camber line given by its slope dz/dx(x).

    Args:
        camber_slope: callable returning dz/dx at chordwise positions 0 <= x <= 1.
        breakpoints: x locations where the slope has kinks (improves the quadrature).
        n_terms: number of Glauert coefficients A1..An to keep.
    """

    def __init__(
        self,
        camber_slope: Callable[[float], float],
        breakpoints: tuple[float, ...] = (),
        n_terms: int = 20,
    ):
        self.camber_slope = camber_slope
        pts = [float(np.arccos(1 - 2 * x)) for x in breakpoints if 0 < x < 1]

        def integral(weight):
            def f(th):
                return float(camber_slope(_x_of_theta(th))) * weight(th)

            return quad(f, 0.0, np.pi, points=pts or None, limit=200, epsabs=1e-13)[0]

        # B0 = 1/pi int dz/dx dtheta, so A0 = alpha - B0
        self._b0 = integral(lambda th: 1.0) / np.pi
        self.An = np.array(
            [2 / np.pi * integral(lambda th, n=n: np.cos(n * th)) for n in range(1, n_terms + 1)]
        )

    @classmethod
    def from_naca4(cls, airfoil: NACA4 | str, n_terms: int = 20) -> ThinAirfoil:
        if isinstance(airfoil, str):
            airfoil = NACA4.from_designation(airfoil)
        bps = () if airfoil.is_symmetric else (airfoil.p,)
        return cls(airfoil.camber_slope, bps, n_terms)

    # --- Glauert coefficients -------------------------------------------------------------
    def A0(self, alpha: float) -> float:
        return alpha - self._b0

    @property
    def A1(self) -> float:
        return float(self.An[0])

    @property
    def A2(self) -> float:
        return float(self.An[1])

    # --- integrated quantities -----------------------------------------------------------
    lift_slope = 2 * np.pi  # per radian, independent of camber

    @property
    def alpha_zero_lift(self) -> float:
        """Zero-lift angle of attack [rad]: alpha_L0 = B0 - A1 / 2."""
        return self._b0 - 0.5 * self.A1

    def cl(self, alpha: float) -> float:
        return np.pi * (2 * self.A0(alpha) + self.A1)

    @property
    def cm_quarter_chord(self) -> float:
        """Moment coefficient about c/4 (independent of alpha: c/4 is the aerodynamic centre)."""
        return np.pi / 4 * (self.A2 - self.A1)

    def cm_le(self, alpha: float) -> float:
        """Moment coefficient about the leading edge (nose-up positive)."""
        return -(self.cl(alpha) / 4 + np.pi / 4 * (self.A1 - self.A2))

    def x_cp(self, alpha: float) -> float:
        """Centre of pressure location x/c (diverges as cl -> 0 for cambered sections)."""
        return 0.25 * (1 + np.pi / self.cl(alpha) * (self.A1 - self.A2))

    def delta_cp(self, x: np.ndarray, alpha: float) -> np.ndarray:
        """Pressure loading Cp_lower - Cp_upper = 2 gamma / V along the chord."""
        theta = np.arccos(1 - 2 * np.asarray(x, dtype=float))
        n = np.arange(1, len(self.An) + 1)
        series = np.sin(np.multiply.outer(theta, n)) @ self.An
        with np.errstate(divide="ignore"):
            lead = self.A0(alpha) * (1 + np.cos(theta)) / np.sin(theta)
        return 4 * (lead + series)
