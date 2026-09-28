# CLAUDE.md

From-scratch aerodynamics toolkit (thin airfoil → 2D panel → boundary layer → 3D VLM).
Plan: `AERO-01-aerodynamics-from-scratch-PLAN.md` (tick the phase checkboxes as phases land).

## Commands
- `uv sync` — install deps (numpy, scipy, matplotlib; dev: pytest, ruff)
- `uv run pytest` — tests
- `uv run ruff check . && uv run ruff format --check .` — lint/format (CI runs both)
- `uv run python validation/<script>.py` — regenerate validation figures into `docs/figures/`
- `uv run python notebooks/build_notebooks.py` — regenerate and execute the notebooks (edit cell
  sources there, never the .ipynb JSON)

## Layout
- `src/aero/` — single package; subpackages mirror the plan: `geometry`, `thin_airfoil`,
  `panel_method_2d`, `boundary_layer`, `vortex_lattice_3d`; `singularities.py` shared.
- `tests/` — pytest, one file per module. `validation/` — scripts + cited reference data.
- `docs/` — theory notes per module, validation report, figures. `notebooks/` — executed demos.

## Conventions
- Everything nondimensional: chord = 1, freestream speed = 1 unless an argument says otherwise.
- Angles are **radians** internally; public helpers that take degrees say so in the name (`*_deg`).
- Airfoil coordinates run clockwise: TE → lower surface → LE → upper surface → TE.
- Panel `i` goes from node `i` to node `i+1`; normals point out of the body.
- Pure NumPy/SciPy — no aerodynamics libraries (no XFOIL wrappers, no AeroSandbox, etc.).

## Testing rules
- Every test asserts against a closed-form result or a cited published value; name the source
  in a comment. Never loosen a tolerance to make a failing solver pass without explaining why.
- Reference data in `validation/reference_data/` must carry its citation. Do not invent data.

## Workflow
- Commit small logical changes directly to `main`, push after `ruff` + `pytest` pass.
