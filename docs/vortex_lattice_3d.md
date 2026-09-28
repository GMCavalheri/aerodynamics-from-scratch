# 3D vortex lattice method

Code: [`src/aero/vortex_lattice_3d/`](../src/aero/vortex_lattice_3d) ·
Notebooks: [`04_vortex_lattice.ipynb`](../notebooks/04_vortex_lattice.ipynb),
[`05_planform_studies.ipynb`](../notebooks/05_planform_studies.ipynb)

## Biot–Savart

A straight vortex filament of circulation $\Gamma$ from $A$ to $B$ induces at $P$

```math
\mathbf V = \frac{\Gamma}{4\pi}\,\frac{\mathbf r_1\times\mathbf r_2}{|\mathbf r_1\times\mathbf r_2|^2}\,
\mathbf r_0\cdot\left(\frac{\mathbf r_1}{r_1} - \frac{\mathbf r_2}{r_2}\right),
\qquad \mathbf r_1 = P - A,\ \mathbf r_2 = P - B,\ \mathbf r_0 = B - A,
```

and a semi-infinite filament starting at $A$ along $\hat{\mathbf d}$ induces
$\frac{\Gamma}{4\pi}\frac{\hat{\mathbf d}\times\mathbf r}{|\hat{\mathbf d}\times\mathbf r|^2}(1 + \hat{\mathbf d}\cdot\hat{\mathbf r})$.
A horseshoe vortex is: in from $+\infty$ to $A$, bound segment $A\to B$, out from $B$ to
$+\infty$. Points on a filament axis get zero velocity (cut-off core).

## Discretisation

The planform (chord law $c(\eta)$, quarter-chord sweep, dihedral) is split into
$n_{span}\times n_{chord}$ panels. Each panel carries a horseshoe with its bound segment on
the panel's quarter-chord line and a control point at three-quarter chord — the classic
"1/4–3/4 rule", which reproduces the 2D flat-plate result $c_l = 2\pi\alpha$ exactly with one
panel.

**Spanwise spacing.** Nodes use cosine spacing, $y = -\tfrac b2\cos\phi$ with uniform $\phi$.
The control points are placed at the $\phi$ mid-points, not the $y$ mid-points. With that
placement the discrete Trefftz-plane drag of an elliptic loading is exact (the property behind
Lan's quasi-vortex-lattice method), and $e$ settles to within 0.1 % of its converged value by about 24
spanwise panels. With $y$-midpoints or uniform spacing, $e$ is over-predicted by several percent
and converges only like $1/N$ (see `docs/figures/vlm_convergence.png`).

**Twist** enters the boundary condition, not the geometry: each control-point normal is
rotated by the local twist, $\hat{\mathbf n} = \hat{\mathbf n}_{flat}\cos\varepsilon + \hat{\mathbf x}\sin\varepsilon$.
Building twist into the panel geometry warps the panels and lifts the control points off the
plane of the trailing legs, which produced spurious tip loads.

## Solution and loads

Flow tangency at every control point,
$\left(\mathbf V_\infty + \sum_j \Gamma_j \mathbf V_{ij}\right)\cdot\hat{\mathbf n}_i = 0$, is a
dense linear system for the $\Gamma_j$.

- **Lift and moment:** Kutta–Joukowski on each bound segment,
  $\mathbf F_j = \rho\ \Gamma_j\ \mathbf V_\infty\times\mathbf l_j$.
- **Induced drag:** in the Trefftz plane far downstream the wake is a row of 2D point vortices
  of strength $\Delta\Gamma_k$ at the strip edges;
  $D_i = -\tfrac{\rho}{2}\sum_j\Gamma_j\ (\mathbf w_j\cdot\hat{\mathbf n}_j)\ \Delta s_j$.
  This far-field evaluation is much less sensitive to the lattice than integrating near-field
  forces.
- **Span efficiency:** $e = C_L^2/(\pi AR\ C_{D_i})$.

## Lifting-line reference

[`lifting_line.py`](../src/aero/vortex_lattice_3d/lifting_line.py) solves Prandtl's monoplane
equation with Glauert's Fourier series,
$\Gamma = 2bV_\infty\sum A_n\sin n\theta$: $C_L = \pi AR A_1$,
$C_{D_i} = \pi AR\sum nA_n^2$. For an elliptic planform it gives $A_{n>1} = 0$ and
$C_{L\alpha} = 2\pi/(1 + 2/AR)$ exactly.

## Validation

| Check | VLM | Reference |
|---|---|---|
| Elliptic wing, AR 4 / 8 / 20: $e$ | 0.998 / 0.999 / 1.000 | 1 (Prandtl) |
| Elliptic AR 100: $C_{L\alpha}$ | within 0.5 % | lifting line |
| Circular wing ($AR = 4/\pi$): $C_{L\alpha}$ | 1.781 | 1.790 (Kinner, exact lifting surface) |
| AR 5, λ 0.5, Λ 45°: $C_{L\alpha}$ | 3.415 /rad | 3.44 /rad (Bertin & Smith, Ex. 7.2) |
| Biot–Savart: infinite-line limit, quadrature, semi-infinite foot | to 1e-6 | closed form |

Lifting-surface $C_{L\alpha}$ sits below lifting-line theory at low aspect ratio (the chord is
not small compared with the span); Helmholtz's $2\pi AR/(2 + \sqrt{AR^2 + 4})$ is within 3 %
at AR 8.

**Rectangular wings.** Lifting line gives $e = 0.937$ at AR 8, the VLM 0.972. The two
distributions agree to about 2 % inboard and differ mainly within a couple of chords of the tip, where lifting-line's
strip assumption over-loads a wide, square tip. A Fourier fit of the VLM circulation gives the
same $e$ as the Trefftz-plane calculation, so the difference is in the model, not the drag
evaluation.

## Limits

Linear (small angles, planar wake along $+x$), inviscid (no profile drag, no stall), no
thickness; camber could be added through the normals the same way as twist.
