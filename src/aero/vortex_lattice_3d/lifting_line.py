"""Prandtl's lifting-line theory solved with Glauert's Fourier-series method.

Gamma(theta) = 2 b V sum_n A_n sin(n theta),  y = -(b/2) cos(theta)

Collocation of the monoplane equation at N stations gives A_n, then
CL = pi AR A1, CDi = pi AR sum n A_n^2, e = 1 / (1 + delta), delta = sum_{n>1} n (A_n/A1)^2
(Anderson, *Fundamentals of Aerodynamics*, Sec. 5.3). Used as the classical reference for
the VLM on unswept, high-aspect-ratio wings.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from aero.vortex_lattice_3d.geometry import Wing


@dataclass
class LiftingLineSolution:
    wing: Wing
    alpha: float
    A: np.ndarray  # A_1..A_N

    @property
    def CL(self) -> float:
        return float(np.pi * self.wing.aspect_ratio * self.A[0])

    @property
    def CDi(self) -> float:
        n = np.arange(1, len(self.A) + 1)
        return float(np.pi * self.wing.aspect_ratio * np.sum(n * self.A**2))

    @property
    def span_efficiency(self) -> float:
        n = np.arange(1, len(self.A) + 1)
        return float(1 / (1 + np.sum(n[1:] * (self.A[1:] / self.A[0]) ** 2)))

    def gamma(self, y: np.ndarray) -> np.ndarray:
        theta = np.arccos(np.clip(-2 * np.asarray(y) / self.wing.span, -1, 1))
        n = np.arange(1, len(self.A) + 1)
        return 2 * self.wing.span * np.sin(np.multiply.outer(theta, n)) @ self.A


def solve_lifting_line(
    wing: Wing,
    alpha: float,
    n: int = 40,
    a0: float = 2 * np.pi,
    alpha_l0: float = 0.0,
) -> LiftingLineSolution:
    """Solve lifting-line theory with ``n`` Fourier terms (section lift slope ``a0``)."""
    theta = np.linspace(0, np.pi, n + 2)[1:-1]
    eta = np.abs(np.cos(theta))
    c = wing.chord(eta)
    k = np.arange(1, n + 1)
    s = np.sin(np.outer(theta, k))
    m = s * (4 * wing.span / (a0 * c))[:, None] + s * k[None, :] / np.sin(theta)[:, None]
    rhs = alpha - alpha_l0 + wing.twist_at(eta)
    return LiftingLineSolution(wing, alpha, np.linalg.solve(m, rhs))
