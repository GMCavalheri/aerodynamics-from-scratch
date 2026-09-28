"""Laminar-turbulent transition prediction.

Michel's (1951) criterion as refitted by Cebeci & Smith (1974):

    Re_theta,tr = 1.174 (1 + 22400 / Re_x) Re_x^0.46

valid for 0.1e6 < Re_x < 40e6 (White, *Viscous Fluid Flow*, Sec. 5-5; Cebeci & Cousteix,
*Modeling and Computation of Boundary-Layer Flows*).
"""

from __future__ import annotations

import numpy as np


def michel_re_theta(re_x):
    re_x = np.asarray(re_x, dtype=float)
    with np.errstate(divide="ignore"):
        return 1.174 * (1 + 22400 / re_x) * re_x**0.46


def michel_transition(s, re_theta, re: float) -> int | None:
    """Index of the first station where Re_theta exceeds Michel's curve, or ``None``."""
    re_x = re * np.asarray(s, dtype=float)
    crit = michel_re_theta(np.maximum(re_x, 1.0))
    hit = np.nonzero((re_theta >= crit) & (re_x > 0))[0]
    return int(hit[0]) if len(hit) else None
