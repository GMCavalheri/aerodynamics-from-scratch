# PLAN.md — AERO-01: Aerodynamics From Scratch

## Overview
A from-scratch implementation of classical and computational aerodynamics, moving from thin airfoil theory through 2D panel methods to the 3D vortex lattice method (VLM), built entirely in Python/NumPy with no black-box CFD libraries. This is the first of five projects covering the classical pillars of aerospace engineering.

## Objective
Build a working aerodynamic analysis toolkit capable of predicting lift, drag (induced), and pitching moment for airfoils and finite wings, validated against known analytical and experimental results (e.g., NACA airfoil data).

## Why This Project
Aerodynamics is the pillar most candidates approach only through commercial software (XFOIL, ANSYS Fluent) without understanding the underlying numerics. Implementing panel methods and VLM from scratch demonstrates genuine command of potential flow theory and numerical methods — the same mathematical machinery (linear systems, singularity distributions, boundary value problems) that appears throughout aerospace engineering.

## Knowledge Covered
1. **Thin Airfoil Theory**: analytical solution for cambered airfoils, lift curve slope, moment coefficient about the quarter-chord.
2. **2D Panel Method (Source/Vortex Panel)**: discretize an arbitrary airfoil geometry (e.g., NACA 4-digit series) into panels, solve for singularity strengths satisfying the no-penetration boundary condition, compute surface pressure distribution and Cl.
3. **Viscous Corrections**: integrate a simple boundary layer model (e.g., Thwaites' method) to estimate skin friction drag and predict separation onset — enough to discuss viscous drag without a full Navier-Stokes solve.
4. **3D Vortex Lattice Method (VLM)**: extend to finite wings — discretize the wing planform into horseshoe vortices, solve for circulation distribution, compute total lift, induced drag (via Trefftz plane analysis), and spanwise loading.
5. **Wing Planform Studies**: compare elliptical, tapered, and swept planforms; reproduce the classical result that elliptical loading minimizes induced drag.

## Prerequisites / Background
- Potential flow theory: sources, sinks, vortices, doublets, superposition
- Linear algebra (solving systems for unknown singularity strengths)
- Basic numerical integration

## Data / Dataset
- **Airfoil geometry**: NACA 4-digit airfoils generated analytically; additional coordinates from the public UIUC Airfoil Coordinates Database.
- **Reference results**: published NACA experimental airfoil data (e.g., NACA 0012, NACA 2412) and XFOIL runs (free, open-source) for Cl, Cp, and Cd comparison.
- **Closed-form benchmarks**: thin airfoil theory and Prandtl lifting-line solutions (elliptical wing).
- No proprietary data required.

## Tech Stack
- Python, NumPy, SciPy (linear solvers)
- Matplotlib/Plotly for pressure distributions, streamlines, spanwise loading plots
- Pytest for validation against known analytical cases

## Repository Structure
```
aerodynamics-from-scratch/
  src/
    thin_airfoil/
    panel_method_2d/
    boundary_layer/
    vortex_lattice_3d/
  validation/         # comparison against NACA/XFOIL reference data
  notebooks/          # exploratory analysis and visualizations
  docs/                # written explanations of the theory per module
  tests/
```

## Execution Plan
- [x] **Phase 1 — Thin airfoil theory**: implement and validate against NACA 4-digit analytical results.
- [x] **Phase 2 — 2D panel method**: implement source/vortex panel method, validate Cl and Cp distribution against published NACA airfoil data (e.g., NACA 0012, NACA 2412).
- [ ] **Phase 3 — Boundary layer / viscous drag**: implement Thwaites' method, estimate Cd, compare against XFOIL reference values.
- [ ] **Phase 4 — 3D Vortex Lattice Method**: implement horseshoe vortex discretization, solve for circulation, validate against Prandtl's lifting-line results for elliptical wings.
- [ ] **Phase 5 — Planform trade studies**: compare taper ratio, aspect ratio, and sweep effects on induced drag and spanwise loading.
- [ ] **Phase 6 — Documentation and publishing**: write theory notes per module, publish repository with visualizations (pressure distributions, streamlines, loading diagrams).

## Validation Strategy
- Compare 2D panel method output against published NACA airfoil experimental/XFOIL data.
- Compare VLM total lift and induced drag against Prandtl's classical lifting-line theory for an elliptical wing (closed-form solution exists).
- Unit tests for each singularity type (source, vortex, doublet) against known closed-form velocity fields.

## Deliverables
- Public repository with modular, documented aerodynamic solvers
- Validation report comparing results against classical/published data
- Visualizations: pressure distributions, spanwise lift distribution, induced drag vs. aspect ratio

## Portfolio Differentiators
- Most candidates use aerodynamics software as a black box; this project proves you understand what's inside it.
- Sets up the numerical/physics foundation directly reused in the Aircraft Conceptual Design capstone project.

## Possible Extensions
- Add a simple 2D Euler solver (finite volume) as a bridge toward full CFD
- Extend VLM to unsteady aerodynamics (time-varying circulation, basic flutter precursor analysis)
