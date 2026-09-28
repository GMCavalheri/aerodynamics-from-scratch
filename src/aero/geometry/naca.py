"""NACA 4-digit airfoil geometry.

Equations from Abbott & von Doenhoff, *Theory of Wing Sections* (1959), Sec. 6.4.
All lengths are normalised by the chord (0 <= x <= 1).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# Thickness polynomial coefficients. The last one is -0.1015 for the original open
# trailing edge and -0.1036 for a closed (zero-thickness) trailing edge.
_A = (0.2969, -0.1260, -0.3516, 0.2843)
_A4_OPEN = -0.1015
_A4_CLOSED = -0.1036


@dataclass(frozen=True)
class NACA4:
    """A NACA 4-digit section, e.g. ``NACA4.from_designation("2412")``.

    Attributes:
        m: maximum camber as a fraction of chord (first digit / 100).
        p: chordwise location of maximum camber (second digit / 10).
        t: maximum thickness as a fraction of chord (last two digits / 100).
        closed_te: use the modified coefficient that closes the trailing edge.
    """

    m: float
    p: float
    t: float
    closed_te: bool = True

    @classmethod
    def from_designation(cls, code: str, closed_te: bool = True) -> NACA4:
        code = code.strip().upper().removeprefix("NACA").strip()
        if len(code) != 4 or not code.isdigit():
            raise ValueError(f"expected a 4-digit NACA designation, got {code!r}")
        return cls(int(code[0]) / 100, int(code[1]) / 10, int(code[2:]) / 100, closed_te)

    @property
    def name(self) -> str:
        return f"NACA {round(self.m * 100)}{round(self.p * 10)}{round(self.t * 100):02d}"

    @property
    def is_symmetric(self) -> bool:
        return self.m == 0 or self.p == 0

    def thickness(self, x: np.ndarray) -> np.ndarray:
        """Half-thickness y_t(x) measured perpendicular to the camber line."""
        x = np.asarray(x, dtype=float)
        a4 = _A4_CLOSED if self.closed_te else _A4_OPEN
        a0, a1, a2, a3 = _A
        return 5 * self.t * (a0 * np.sqrt(x) + a1 * x + a2 * x**2 + a3 * x**3 + a4 * x**4)

    def camber(self, x: np.ndarray) -> np.ndarray:
        """Mean camber line y_c(x)."""
        x = np.asarray(x, dtype=float)
        if self.is_symmetric:
            return np.zeros_like(x)
        m, p = self.m, self.p
        fwd = m / p**2 * (2 * p * x - x**2)
        aft = m / (1 - p) ** 2 * ((1 - 2 * p) + 2 * p * x - x**2)
        return np.where(x < p, fwd, aft)

    def camber_slope(self, x: np.ndarray) -> np.ndarray:
        """dy_c/dx."""
        x = np.asarray(x, dtype=float)
        if self.is_symmetric:
            return np.zeros_like(x)
        m, p = self.m, self.p
        fwd = 2 * m / p**2 * (p - x)
        aft = 2 * m / (1 - p) ** 2 * (p - x)
        return np.where(x < p, fwd, aft)

    def surfaces(self, x: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Upper and lower surface points (xu, yu, xl, yl) at camber-line stations x."""
        x = np.asarray(x, dtype=float)
        yt = self.thickness(x)
        yc = self.camber(x)
        theta = np.arctan(self.camber_slope(x))
        xu = x - yt * np.sin(theta)
        yu = yc + yt * np.cos(theta)
        xl = x + yt * np.sin(theta)
        yl = yc - yt * np.cos(theta)
        return xu, yu, xl, yl

    def coordinates(self, n_panels: int = 160) -> np.ndarray:
        """Closed contour for a panel method, shape (n_panels + 1, 2).

        Points use cosine (full-circle) spacing, which clusters nodes at the leading and
        trailing edges, and run clockwise: TE -> lower surface -> LE -> upper surface -> TE.
        The first and last points coincide at the trailing edge when ``closed_te`` is set.
        ``n_panels`` must be even so that a node sits exactly on the leading edge.
        """
        if n_panels < 4 or n_panels % 2:
            raise ValueError("n_panels must be an even integer >= 4")
        n_side = n_panels // 2
        beta = np.linspace(0.0, np.pi, n_side + 1)
        x = 0.5 * (1 - np.cos(beta))  # 0 (LE) -> 1 (TE)
        xu, yu, xl, yl = self.surfaces(x)
        lower = np.column_stack([xl[::-1], yl[::-1]])  # TE -> LE
        upper = np.column_stack([xu[1:], yu[1:]])  # LE (excluded) -> TE
        return np.vstack([lower, upper])
