"""Thwaites' integral method for laminar boundary layers.

Thwaites (1949) integrated the momentum-integral equation using one-parameter correlations
in the pressure-gradient parameter lambda = (theta^2 / nu) dUe/ds:

    theta^2 = 0.45 nu / Ue^6  int_0^s Ue^5 ds

Shear and shape-factor correlations l(lambda), H(lambda) are the Cebeci & Bradshaw (1977)
fits quoted in White, *Viscous Fluid Flow*, 3rd ed., Sec. 4-6. Separation is predicted at
lambda = -0.09.

Everything is nondimensional: lengths by chord c, velocities by V_inf, Re = V_inf c / nu.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

LAMBDA_SEPARATION = -0.09


def shear_correlation(lam):
    """l(lambda) = tau_w theta / (mu Ue)."""
    lam = np.asarray(lam, dtype=float)
    fav = 0.22 + 1.57 * lam - 1.8 * lam**2
    adv = 0.22 + 1.402 * lam + 0.018 * lam / (lam + 0.107)
    return np.where(lam >= 0, fav, adv)


def shape_factor(lam):
    """H(lambda) = delta* / theta."""
    lam = np.asarray(lam, dtype=float)
    fav = 2.61 - 3.75 * lam + 5.24 * lam**2
    adv = 2.088 + 0.0731 / (lam + 0.14)
    return np.where(lam >= 0, fav, adv)


@dataclass
class LaminarBL:
    s: np.ndarray
    ue: np.ndarray
    theta: np.ndarray
    lam: np.ndarray
    H: np.ndarray
    cf: np.ndarray
    re: float
    i_sep: int | None  # first station with lambda <= -0.09, if any

    @property
    def delta_star(self) -> np.ndarray:
        return self.H * self.theta

    @property
    def re_theta(self) -> np.ndarray:
        return self.re * self.ue * self.theta

    @property
    def s_sep(self) -> float | None:
        return None if self.i_sep is None else float(self.s[self.i_sep])


def thwaites(s, ue, re: float, theta0: float = 0.0) -> LaminarBL:
    """March Thwaites' method along stations ``s`` with edge velocity ``ue``.

    If ``ue[0] == 0`` (stagnation point) the Hiemenz limit theta^2 = 0.075 / (Re dUe/ds)
    is used there. ``theta0`` sets a nonzero initial momentum thickness otherwise.
    """
    s = np.asarray(s, dtype=float)
    ue = np.asarray(ue, dtype=float)
    due = np.gradient(ue, s)
    # int Ue^5 ds, exact for piecewise-linear Ue (plain trapezoids are badly off near a
    # stagnation point, where Ue ~ s and the integrand ~ s^5)
    u1, u2 = ue[:-1], ue[1:]
    seg = np.diff(s) / 6 * sum(u1**k * u2 ** (5 - k) for k in range(6))
    integral = np.concatenate([[0.0], np.cumsum(seg)])
    with np.errstate(divide="ignore", invalid="ignore"):
        theta2 = (0.45 * integral / re + theta0**2 * ue[0] ** 6) / ue**6
    if ue[0] == 0.0:
        theta2[0] = 0.075 / (re * due[0])
    theta = np.sqrt(theta2)
    lam = np.clip(theta2 * re * due, -0.1, 0.25)
    l = shear_correlation(lam)
    H = shape_factor(lam)
    with np.errstate(divide="ignore", invalid="ignore"):
        cf = np.where(ue > 0, 2 * l / (re * ue * theta), 0.0)
    below = np.nonzero(theta2 * re * due <= LAMBDA_SEPARATION)[0]
    return LaminarBL(s, ue, theta, lam, H, cf, re, int(below[0]) if len(below) else None)
