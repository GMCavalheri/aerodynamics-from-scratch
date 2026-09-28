"""Head's entrainment method for turbulent boundary layers.

Momentum integral plus Head's (1958) entrainment equation, closed with the
Ludwieg-Tillmann skin-friction law (Cebeci & Bradshaw, *Momentum Transfer in Boundary
Layers*, Sec. 8.2; Moran, *Theoretical and Computational Aerodynamics*, Sec. 7.8):

    d theta / ds       = Cf/2 - (H + 2) (theta / Ue) dUe/ds
    d(Ue theta H1)/ds  = Ue * 0.0306 (H1 - 3)^-0.6169
    Cf                 = 0.246 * 10^(-0.678 H) Re_theta^-0.268

Separation is flagged when H exceeds ``H_SEPARATION``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator

H_SEPARATION = 2.4


def h1_of_h(H):
    """Head's mass-flow shape factor H1 = (delta - delta*) / theta."""
    H = np.asarray(H, dtype=float)
    lo = 3.3 + 0.8234 * np.maximum(H - 1.1, 1e-9) ** -1.287
    hi = 3.3 + 1.5501 * (H - 0.6778) ** -3.064
    return np.where(H <= 1.6, lo, hi)


def h_of_h1(H1):
    """Inverse of :func:`h1_of_h`."""
    H1 = np.maximum(np.asarray(H1, dtype=float), 3.3 + 1e-6)
    lo = 1.1 + ((H1 - 3.3) / 0.8234) ** (-1 / 1.287)
    hi = 0.6778 + ((H1 - 3.3) / 1.5501) ** (-1 / 3.064)
    return np.where(H1 >= 5.3, lo, hi)  # H1(1.6) = 5.3 joins the two branches


def entrainment(H1):
    return 0.0306 * np.maximum(np.asarray(H1, dtype=float) - 3.0, 1e-6) ** -0.6169


def ludwieg_tillmann(H, re_theta):
    return 0.246 * 10 ** (-0.678 * np.asarray(H)) * np.asarray(re_theta) ** -0.268


@dataclass
class TurbulentBL:
    s: np.ndarray
    ue: np.ndarray
    theta: np.ndarray
    H: np.ndarray
    cf: np.ndarray
    re: float
    i_sep: int | None

    @property
    def s_sep(self) -> float | None:
        return None if self.i_sep is None else float(self.s[self.i_sep])


def head(s, ue, re: float, theta0: float, H0: float = 1.4) -> TurbulentBL:
    """Integrate Head's method from ``s[0]`` with initial theta0 and shape factor H0."""
    s = np.asarray(s, dtype=float)
    ue = np.asarray(ue, dtype=float)
    ue_f = PchipInterpolator(s, ue)
    due_f = ue_f.derivative()

    def rhs(x, y):
        theta, q = y  # q = Ue * theta * H1
        u, du = ue_f(x), due_f(x)
        theta = max(theta, 1e-12)
        H1 = q / (u * theta)
        H = h_of_h1(H1)
        cf = ludwieg_tillmann(H, re * u * theta)
        return [0.5 * cf - (H + 2) * theta / u * du, u * entrainment(H1)]

    y0 = [theta0, ue[0] * theta0 * float(h1_of_h(H0))]
    if len(s) == 1:
        theta, q = np.array([y0[0]]), np.array([y0[1]])
    else:
        sol = solve_ivp(rhs, (s[0], s[-1]), y0, t_eval=s, method="LSODA", rtol=1e-8, atol=1e-12)
        theta, q = sol.y
    H = h_of_h1(q / (ue * theta))
    cf = ludwieg_tillmann(H, re * ue * theta)
    sep = np.nonzero(H >= H_SEPARATION)[0]
    return TurbulentBL(s, ue, theta, H, cf, re, int(sep[0]) if len(sep) else None)
