"""Planform trade studies built on the VLM (and lifting line for comparison).

All studies run at a small angle of attack; the quantities returned (span efficiency,
lift-curve slope, normalised loadings) are independent of alpha in linear theory.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from aero.vortex_lattice_3d.geometry import Wing
from aero.vortex_lattice_3d.lifting_line import solve_lifting_line
from aero.vortex_lattice_3d.vlm import VLMSolution, solve_vlm

ALPHA = 4 * np.pi / 180
MESH = {"n_span": 40, "n_chord": 6}


@dataclass
class TaperStudy:
    taper: np.ndarray
    e_vlm: np.ndarray
    e_lifting_line: np.ndarray

    @property
    def optimum_vlm(self) -> float:
        return float(self.taper[np.argmax(self.e_vlm)])

    @property
    def optimum_lifting_line(self) -> float:
        return float(self.taper[np.argmax(self.e_lifting_line)])


def taper_study(aspect_ratio: float, tapers, sweep: float = 0.0) -> TaperStudy:
    """Span efficiency vs taper ratio (Anderson Fig. 5.20 / Glauert's delta curves)."""
    tapers = np.asarray(tapers, dtype=float)
    e_vlm, e_ll = [], []
    for lam in tapers:
        w = Wing.tapered(aspect_ratio, lam, sweep=sweep)
        e_vlm.append(solve_vlm(w, ALPHA, **MESH).span_efficiency)
        e_ll.append(solve_lifting_line(w, ALPHA, 60).span_efficiency)
    return TaperStudy(tapers, np.array(e_vlm), np.array(e_ll))


def induced_drag_vs_aspect_ratio(make_wing, aspect_ratios, cl: float = 0.5):
    """CDi at a fixed CL for wings built by ``make_wing(AR)``; returns (CDi, e)."""
    cdi, e = [], []
    for ar in aspect_ratios:
        s = solve_vlm(make_wing(ar), ALPHA, **MESH)
        e.append(s.span_efficiency)
        cdi.append(cl**2 / (np.pi * ar * s.span_efficiency))
    return np.array(cdi), np.array(e)


def sweep_study(aspect_ratio: float, taper: float, sweeps) -> list[VLMSolution]:
    return [solve_vlm(Wing.tapered(aspect_ratio, taper, sweep=sw), ALPHA, **MESH) for sw in sweeps]
