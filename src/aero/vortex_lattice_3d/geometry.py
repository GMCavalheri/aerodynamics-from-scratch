"""Wing planform description and vortex-lattice discretisation."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Wing:
    """A symmetric wing defined by its semi-span chord distribution.

    Attributes:
        span: tip-to-tip span b.
        chord: callable c(eta) for 0 <= eta = 2|y|/b <= 1.
        sweep: quarter-chord sweep angle [rad] (positive = aft).
        dihedral: dihedral angle [rad].
        twist: callable geometric twist epsilon(eta) [rad], positive nose-up (washout < 0).
        name: label for plots.
    """

    span: float
    chord: Callable[[np.ndarray], np.ndarray]
    sweep: float = 0.0
    dihedral: float = 0.0
    twist: Callable[[np.ndarray], np.ndarray] | None = None
    name: str = "wing"

    # --- factories -------------------------------------------------------------------------
    @classmethod
    def rectangular(cls, aspect_ratio: float, chord: float = 1.0, **kw) -> Wing:
        return cls(
            aspect_ratio * chord,
            lambda eta: np.full_like(np.asarray(eta, float), chord),
            name=kw.pop("name", f"Rectangular AR {aspect_ratio:g}"),
            **kw,
        )

    @classmethod
    def tapered(cls, aspect_ratio: float, taper: float, root_chord: float = 1.0, **kw) -> Wing:
        # S = b c_r (1 + lambda) / 2 and AR = b^2 / S  ->  b = AR c_r (1 + lambda) / 2
        span = aspect_ratio * root_chord * (1 + taper) / 2
        return cls(
            span,
            lambda eta: root_chord * (1 - (1 - taper) * np.asarray(eta, float)),
            name=kw.pop("name", f"Tapered AR {aspect_ratio:g}, λ = {taper:g}"),
            **kw,
        )

    @classmethod
    def elliptic(cls, aspect_ratio: float, root_chord: float = 1.0, **kw) -> Wing:
        # S = pi b c_r / 4  ->  b = AR pi c_r / 4
        span = aspect_ratio * np.pi * root_chord / 4
        return cls(
            span,
            lambda eta: root_chord * np.sqrt(np.clip(1 - np.asarray(eta, float) ** 2, 0, None)),
            name=kw.pop("name", f"Elliptic AR {aspect_ratio:g}"),
            **kw,
        )

    # --- integral properties -----------------------------------------------------------------
    @property
    def area(self) -> float:
        # Gauss-Chebyshev-friendly substitution eta = sin(t) handles the elliptic tip
        t = np.linspace(0, np.pi / 2, 4001)
        return float(self.span * np.trapezoid(self.chord(np.sin(t)) * np.cos(t), t))

    @property
    def aspect_ratio(self) -> float:
        return self.span**2 / self.area

    @property
    def mean_aerodynamic_chord(self) -> float:
        t = np.linspace(0, np.pi / 2, 4001)
        c = self.chord(np.sin(t))
        return float(self.span * np.trapezoid(c**2 * np.cos(t), t) / self.area)

    def twist_at(self, eta):
        eta = np.asarray(eta, float)
        return np.zeros_like(eta) if self.twist is None else self.twist(eta)


@dataclass
class Lattice:
    """Vortex lattice: nc x ns panels, each with a horseshoe vortex and a control point."""

    corners: np.ndarray  # (nc + 1, ns + 1, 3) panel corner grid (LE row first)
    bound_a: np.ndarray  # (nc, ns, 3) port end of bound vortex (1/4 panel chord)
    bound_b: np.ndarray  # (nc, ns, 3) starboard end
    control: np.ndarray  # (nc, ns, 3) control point (3/4 panel chord, mid-span)
    normal: np.ndarray  # (nc, ns, 3) unit normal at the control point
    area: np.ndarray  # (nc, ns) panel areas
    y_nodes: np.ndarray  # (ns + 1,) spanwise node positions
    chord_nodes: np.ndarray  # (ns + 1,) local chords at the nodes
    span_frac: np.ndarray  # (ns,) spanwise position of the control points within each strip


def spanwise_nodes(n_span: int, spacing: str = "cosine") -> tuple[np.ndarray, np.ndarray]:
    """Spanwise nodes and control stations on eta_signed in [-1, 1].

    Cosine spacing clusters nodes at the tips, y = -cos(phi) with uniform phi. Its control
    stations sit at the phi mid-points rather than the y mid-points: with that placement
    the discrete Trefftz-plane drag of an elliptic loading is exact (the same property
    that underlies Lan's quasi-vortex-lattice method, J. Aircraft 11(9), 1974), and the
    lattice converges much faster near the tips.
    """
    if spacing == "cosine":
        phi = np.linspace(0, np.pi, n_span + 1)
        return -np.cos(phi), -np.cos(0.5 * (phi[:-1] + phi[1:]))
    if spacing == "uniform":
        nodes = np.linspace(-1, 1, n_span + 1)
        return nodes, 0.5 * (nodes[:-1] + nodes[1:])
    raise ValueError(f"unknown spacing {spacing!r}")


def build_lattice(
    wing: Wing, n_span: int = 40, n_chord: int = 6, spacing: str = "cosine"
) -> Lattice:
    """Discretise the planform into a flat (camber-free) vortex lattice.

    The quarter-chord line is swept by ``wing.sweep`` and raised by the dihedral; chordwise
    divisions are uniform. As usual in linear VLM, twist is not built into the geometry
    (which would warp the panels and put the trailing legs off the control-point plane) but
    enters the boundary condition: each control-point normal is rotated nose-up by the local
    twist epsilon, n = n_flat cos(epsilon) + x sin(epsilon).
    """
    eta_s, eta_c = spanwise_nodes(n_span, spacing)
    span_frac = (eta_c - eta_s[:-1]) / np.diff(eta_s)
    eta = np.abs(eta_s)
    y = eta_s * wing.span / 2
    c = wing.chord(eta)
    x_qc = np.abs(y) * np.tan(wing.sweep)
    z_qc = np.abs(y) * np.tan(wing.dihedral)

    frac = np.linspace(0, 1, n_chord + 1)[:, None]  # chordwise fraction of local chord
    dx = (frac - 0.25) * c  # relative to the quarter chord
    corners = np.stack(
        [x_qc + dx, np.broadcast_to(y, dx.shape), np.broadcast_to(z_qc, dx.shape)], axis=-1
    )

    def lerp(f):
        # point at chordwise fraction f (within each panel) along both side edges
        return corners[:-1] + f * (corners[1:] - corners[:-1])

    qa = lerp(0.25)  # (nc, ns+1, 3)
    tq = lerp(0.75)
    bound_a, bound_b = qa[:, :-1], qa[:, 1:]
    control = tq[:, :-1] + span_frac[None, :, None] * (tq[:, 1:] - tq[:, :-1])

    d1 = corners[1:, 1:] - corners[:-1, :-1]
    d2 = corners[:-1, 1:] - corners[1:, :-1]
    cr = np.cross(d1, d2)
    area = 0.5 * np.linalg.norm(cr, axis=-1)
    normal = cr / np.maximum(np.linalg.norm(cr, axis=-1, keepdims=True), 1e-300)
    # orient normals upward, then tilt them by the local twist
    normal *= np.sign(normal[..., 2:3] + 1e-300)
    eps = wing.twist_at(np.abs(eta_c))[None, :, None]
    normal = normal * np.cos(eps) + np.array([1.0, 0.0, 0.0]) * np.sin(eps)
    return Lattice(corners, bound_a, bound_b, control, normal, area, y, c, span_frac)
