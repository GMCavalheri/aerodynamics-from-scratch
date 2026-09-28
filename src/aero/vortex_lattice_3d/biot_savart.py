"""Biot-Savart law for straight vortex filaments (Katz & Plotkin, Sec. 10.4.5).

Conventions: x downstream, y to the right (starboard), z up. A filament of circulation
Gamma running from A to B induces, by the right-hand rule,

    V = Gamma / (4 pi) * (r1 x r2) / |r1 x r2|^2 * r0 . (r1/|r1| - r2/|r2|)

with r1 = P - A, r2 = P - B, r0 = B - A. Points closer than ``core`` to the filament axis
get zero velocity (a simple cut-off core that removes the line singularity).
"""

from __future__ import annotations

import numpy as np

FOUR_PI = 4.0 * np.pi


def segment_velocity(p, a, b, gamma=1.0, core=1e-10):
    """Velocity at points ``p`` (..., 3) from filaments a -> b (..., 3), broadcasting."""
    r1 = p - a
    r2 = p - b
    r0 = b - a
    cross = np.cross(r1, r2)
    cross2 = np.sum(cross**2, axis=-1)
    n1 = np.linalg.norm(r1, axis=-1)
    n2 = np.linalg.norm(r2, axis=-1)
    len0 = np.linalg.norm(r0, axis=-1)
    safe = (cross2 > (core * len0) ** 2) & (n1 > core) & (n2 > core)
    with np.errstate(divide="ignore", invalid="ignore"):
        k = np.sum(r0 * (r1 / n1[..., None] - r2 / n2[..., None]), axis=-1) / cross2
    k = np.where(safe, k, 0.0)
    return (gamma * k / FOUR_PI)[..., None] * cross


def semi_infinite_velocity(p, a, direction, gamma=1.0, core=1e-10):
    """Velocity from a filament starting at ``a`` and running to infinity along ``direction``.

    V = Gamma / (4 pi) * (d x r) / |d x r|^2 * (1 + d . r / |r|), r = P - A.
    Reverse the sign of ``gamma`` for a filament coming in from infinity to ``a``.
    """
    d = np.asarray(direction, dtype=float)
    d = d / np.linalg.norm(d, axis=-1, keepdims=True)
    r = p - a
    cross = np.cross(d, r)
    cross2 = np.sum(cross**2, axis=-1)
    nr = np.linalg.norm(r, axis=-1)
    safe = (cross2 > core**2) & (nr > core)
    with np.errstate(divide="ignore", invalid="ignore"):
        k = (1 + np.sum(d * r, axis=-1) / nr) / cross2
    k = np.where(safe, k, 0.0)
    return (gamma * k / FOUR_PI)[..., None] * cross


def horseshoe_velocity(p, a, b, direction=(1.0, 0.0, 0.0), gamma=1.0, core=1e-10):
    """Horseshoe vortex: from infinity to A, bound segment A -> B, then B to infinity.

    With A to port of B (A.y < B.y) and the legs trailing downstream, positive ``gamma``
    produces positive lift.
    """
    return (
        -semi_infinite_velocity(p, a, direction, gamma, core)
        + segment_velocity(p, a, b, gamma, core)
        + semi_infinite_velocity(p, b, direction, gamma, core)
    )
