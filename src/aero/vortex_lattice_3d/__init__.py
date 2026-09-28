"""3D vortex lattice method and lifting-line reference solver."""

from aero.vortex_lattice_3d.geometry import Lattice, Wing, build_lattice
from aero.vortex_lattice_3d.lifting_line import LiftingLineSolution, solve_lifting_line
from aero.vortex_lattice_3d.vlm import VLMSolution, solve_vlm

__all__ = [
    "Lattice",
    "LiftingLineSolution",
    "VLMSolution",
    "Wing",
    "build_lattice",
    "solve_lifting_line",
    "solve_vlm",
]
