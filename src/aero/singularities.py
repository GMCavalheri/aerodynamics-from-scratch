"""Elementary 2D potential-flow singularities.

Conventions: vortices are positive counter-clockwise; sources are positive outflow.
All functions broadcast over their array arguments and return velocity components
``(u, v)`` in the global frame.

Panel formulas (constant strength per unit length over a straight panel) follow
Katz & Plotkin, *Low-Speed Aerodynamics*, 2nd ed., Sec. 10.2, rewritten for a
counter-clockwise-positive vortex.
"""

from __future__ import annotations

import numpy as np

TWO_PI = 2.0 * np.pi


def source_velocity(sigma, x0, y0, x, y):
    """Velocity of a point source of strength ``sigma`` (volume flux) at (x0, y0)."""
    dx, dy = np.subtract(x, x0), np.subtract(y, y0)
    r2 = dx**2 + dy**2
    k = np.asarray(sigma) / (TWO_PI * r2)
    return k * dx, k * dy


def vortex_velocity(gamma, x0, y0, x, y):
    """Velocity of a counter-clockwise point vortex of circulation ``gamma`` at (x0, y0)."""
    dx, dy = np.subtract(x, x0), np.subtract(y, y0)
    r2 = dx**2 + dy**2
    k = np.asarray(gamma) / (TWO_PI * r2)
    return -k * dy, k * dx


def doublet_potential(mu, x0, y0, x, y, angle=0.0):
    """Potential of a point doublet, phi = mu/(2 pi) * (d . r) / r^2, d = (cos a, sin a).

    With ``angle = 0`` and ``mu = 2 pi U R^2`` it superposes with a uniform stream U in +x
    to give the flow around a circular cylinder of radius R.
    """
    dx, dy = np.subtract(x, x0), np.subtract(y, y0)
    r2 = dx**2 + dy**2
    return np.asarray(mu) / TWO_PI * (dx * np.cos(angle) + dy * np.sin(angle)) / r2


def doublet_velocity(mu, x0, y0, x, y, angle=0.0):
    """Velocity (gradient of :func:`doublet_potential`)."""
    dx, dy = np.subtract(x, x0), np.subtract(y, y0)
    r2 = dx**2 + dy**2
    c, s = np.cos(angle), np.sin(angle)
    proj = dx * c + dy * s
    k = np.asarray(mu) / TWO_PI
    return k * (c / r2 - 2 * dx * proj / r2**2), k * (s / r2 - 2 * dy * proj / r2**2)


def panel_frame(x1, y1, x2, y2):
    """Unit tangent (node 1 -> node 2), unit normal (tangent rotated +90 deg), and length."""
    dx, dy = np.subtract(x2, x1), np.subtract(y2, y1)
    length = np.hypot(dx, dy)
    tx, ty = dx / length, dy / length
    return tx, ty, -ty, tx, length


def panel_integrals(x, y, x1, y1, x2, y2, on_panel=None):
    """Geometric integrals for a straight panel seen from field point (x, y).

    Returns ``(log_ratio, beta)`` with ``log_ratio = ln(r1 / r2)`` and ``beta`` the angle
    subtended by the panel (theta2 - theta1), positive on the normal (+90 deg) side.
    ``on_panel`` is an optional boolean mask marking field points lying on the panel itself
    (e.g. its own control point); there ``beta`` is set to +pi, the limit approached from
    the normal side, which is ambiguous in floating point.
    """
    r1x, r1y = np.subtract(x, x1), np.subtract(y, y1)
    r2x, r2y = np.subtract(x, x2), np.subtract(y, y2)
    log_ratio = 0.5 * np.log((r1x**2 + r1y**2) / (r2x**2 + r2y**2))
    beta = np.arctan2(r1x * r2y - r1y * r2x, r1x * r2x + r1y * r2y)
    if on_panel is not None:
        beta = np.where(on_panel, np.pi, beta)
    return log_ratio, beta


def _to_global(ul, vl, x1, y1, x2, y2):
    tx, ty, nx, ny, _ = panel_frame(x1, y1, x2, y2)
    return ul * tx + vl * nx, ul * ty + vl * ny


def source_panel_velocity(x, y, x1, y1, x2, y2, on_panel=None):
    """Velocity induced at (x, y) by a unit-strength constant source panel."""
    log_ratio, beta = panel_integrals(x, y, x1, y1, x2, y2, on_panel)
    return _to_global(log_ratio / TWO_PI, beta / TWO_PI, x1, y1, x2, y2)


def vortex_panel_velocity(x, y, x1, y1, x2, y2, on_panel=None):
    """Velocity induced at (x, y) by a unit-strength constant CCW vortex panel."""
    log_ratio, beta = panel_integrals(x, y, x1, y1, x2, y2, on_panel)
    return _to_global(-beta / TWO_PI, log_ratio / TWO_PI, x1, y1, x2, y2)
