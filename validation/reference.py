"""Load the cited reference values in ``validation/reference_data``."""

from __future__ import annotations

import csv
from pathlib import Path

DATA = Path(__file__).resolve().parent / "reference_data"


def load(name: str) -> list[dict[str, str]]:
    with open(DATA / name) as f:
        rows = [line for line in f if not line.startswith("#")]
    return list(csv.DictReader(rows))


def section_value(airfoil: str, quantity: str) -> tuple[float, float, float, str]:
    """(value, low, high, source) for one row of ``naca_section_data.csv``."""
    for row in load("naca_section_data.csv"):
        if row["airfoil"] == airfoil and row["quantity"] == quantity:
            return float(row["value"]), float(row["low"]), float(row["high"]), row["source"]
    raise KeyError((airfoil, quantity))
