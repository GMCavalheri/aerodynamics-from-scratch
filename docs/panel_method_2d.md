# 2D panel method (Hess–Smith)

Code: [`src/aero/panel_method_2d/hess_smith.py`](../src/aero/panel_method_2d/hess_smith.py),
[`src/aero/singularities.py`](../src/aero/singularities.py),
[`src/aero/geometry/`](../src/aero/geometry) ·
Notebook: [`02_panel_method.ipynb`](../notebooks/02_panel_method.ipynb)

## Singularities

Incompressible, irrotational flow satisfies Laplace's equation, so solutions superpose. The
building blocks (counter-clockwise vortices positive):

| Element | Velocity at distance $r$ |
|---|---|
| Source $\sigma$ | $\mathbf{V} = \frac{\sigma}{2\pi r}\,\hat{\mathbf{e}}_r$ |
| Vortex $\Gamma$ | $\mathbf{V} = \frac{\Gamma}{2\pi r}\,\hat{\mathbf{e}}_\theta$ |
| Doublet $\mu$ (axis $\hat{\mathbf d}$) | $\phi = \frac{\mu}{2\pi}\frac{\hat{\mathbf d}\cdot\mathbf r}{r^2}$ |

A straight panel of length $L$ with a constant source (or vortex) density induces, in panel
coordinates, velocities built from two geometric integrals:

$$
\ln\frac{r_1}{r_2}, \qquad \beta = \theta_2 - \theta_1 ,
$$

the log of the distance ratio to the end points and the angle subtended by the panel. For a
source panel $u = \frac{\sigma}{2\pi}\ln\frac{r_1}{r_2}$, $v = \frac{\sigma}{2\pi}\beta$; a vortex
panel is the same field rotated by 90°. On the panel itself $\beta \to \pi$ from the fluid side,
which gives the self-induced normal velocity $\sigma/2$ (a jump the code sets explicitly, since
$\beta$ is ambiguous in floating point there).

## Hess–Smith formulation

The contour is split into $N$ straight panels, ordered clockwise from the trailing edge so
every panel normal points into the fluid. Each panel carries its own source strength
$\sigma_j$; a single vortex strength $\gamma$ is shared by all panels (so the circulation is
$\Gamma = \gamma \sum L_j$). The $N+1$ unknowns follow from

1. **No penetration** at every panel midpoint:
   $\sum_j A^n_{ij}\sigma_j + \gamma\sum_j B^n_{ij} + \mathbf V_\infty\cdot\hat{\mathbf n}_i = 0$;
2. **Kutta condition**: equal tangential speeds on the two trailing-edge panels,
   $V_{t,1} + V_{t,N} = 0$ (their tangents point in opposite directions).

Outputs: surface speed $V_t$, $C_p = 1 - (V_t/V_\infty)^2$, lift from Kutta–Joukowski
($c_l = 2\Gamma/(V_\infty c)$) and from integrating $C_p$, pitching moment, and the velocity
anywhere in the field (used for streamlines).

## Geometry

- **NACA 4-digit** ([`naca.py`](../src/aero/geometry/naca.py)): Abbott & von Doenhoff thickness
  polynomial (closed-TE coefficient $-0.1036$ by default), camber line, cosine spacing so nodes
  cluster at both edges.
- **Joukowski** ([`joukowski.py`](../src/aero/geometry/joukowski.py)): $z = \zeta + b^2/\zeta$
  maps a circle through $\zeta = b$ onto a cusped airfoil. The circle flow with the Kutta
  condition, $\Gamma = 4\pi a V_\infty\sin(\alpha + \beta)$, is exact, giving an exact $C_p$
  benchmark.

## Validation

| Check | Result |
|---|---|
| Circular cylinder (non-lifting), $C_p = 1 - 4\sin^2\theta$ | exact to round-off at the control points, any N |
| Joukowski airfoil, exact $c_l$ and $C_p$ | observed order 0.66 → 0.78; −0.9 % $c_l$ at N = 800 |
| Pressure-integrated vs Kutta–Joukowski $c_l$ | agree within 0.5 % |
| d'Alembert: pressure drag | $|c_d| < 2\times10^{-3}$ (discretisation error) |
| Thickness effect on lift slope | $2\pi(1 + 0.77\,t/c)$ within 2 % |
| Thin section (NACA 2402) | $\alpha_{L0}$, $c_{m,c/4}$ match thin airfoil theory |
| NACA 0012 / 2412 vs measurements | see [validation report](validation-report.md) |

The cusped Joukowski trailing edge is the worst case for constant-strength panels (the Kutta
condition is applied between two nearly parallel panels), which is why convergence there is
below first order (observed order 0.66 → 0.78 between N = 100 and 800). On a smooth NACA 0012
the observed order approaches 1. On a regular polygon the cylinder solution is exact at the
control points, so the cylinder checks assembly and signs but cannot measure discretisation
error.

## Limits

Inviscid, so lift slope and $|c_m|$ are over-predicted compared with experiment (no
boundary-layer de-cambering) and there is no stall. Constant-strength panels give at
most a first-order method; linear-vorticity panels would converge faster.
