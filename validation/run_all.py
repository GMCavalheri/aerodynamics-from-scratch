"""Regenerate every validation figure and results table.

Run with ``uv run python validation/run_all.py``.
"""

from __future__ import annotations

import runpy
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCRIPTS = [
    "panel_method_validation.py",
    "boundary_layer_validation.py",
    "vlm_validation.py",
    "planform_studies.py",
]

if __name__ == "__main__":
    for name in SCRIPTS:
        print(f"== {name}")
        runpy.run_path(str(HERE / name), run_name="__main__")
