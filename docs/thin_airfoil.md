# Thin airfoil theory

Code: [`src/aero/thin_airfoil/theory.py`](../src/aero/thin_airfoil/theory.py) ·
Notebook: [`01_thin_airfoil.ipynb`](../notebooks/01_thin_airfoil.ipynb)

## Model

For a thin, slightly cambered section at small angle of attack, the airfoil is replaced by a
vortex sheet $\gamma(x)$ lying on the chord line. The sheet must make the camber line
$z(x)$ a streamline. Linearising the tangency condition gives the fundamental equation of
thin airfoil theory:

```math
\frac{1}{2\pi}\int_0^c \frac{\gamma(\xi)\,d\xi}{x-\xi} = V_\infty\left(\alpha - \frac{dz}{dx}\right),
\qquad \gamma(c) = 0 \;\text{(Kutta)}.
```

## Glauert's solution

With $x = \tfrac{c}{2}(1-\cos\theta)$ the solution is

```math
\gamma(\theta) = 2V_\infty\left[A_0\,\frac{1+\cos\theta}{\sin\theta} + \sum_{n\ge1} A_n \sin n\theta\right],
```

```math
A_0 = \alpha - \frac{1}{\pi}\int_0^\pi \frac{dz}{dx}\,d\theta_0, \qquad
A_n = \frac{2}{\pi}\int_0^\pi \frac{dz}{dx}\cos n\theta_0\,d\theta_0 .
```

Integrating the loading gives

| Quantity | Result |
|---|---|
| Lift | $c_l = \pi(2A_0 + A_1) = 2\pi(\alpha - \alpha_{L0})$ |
| Zero-lift angle | $\alpha_{L0} = -\frac{1}{\pi}\int_0^\pi \frac{dz}{dx}(\cos\theta_0 - 1)\ d\theta_0$ |
| Moment about c/4 | $c_{m,c/4} = \frac{\pi}{4}(A_2 - A_1)$, independent of $\alpha$ → c/4 is the aerodynamic centre |
| Moment about LE | $c_{m,LE} = -\left[\frac{c_l}{4} + \frac{\pi}{4}(A_1 - A_2)\right]$ |
| Centre of pressure | $x_{cp}/c = \frac14\left[1 + \frac{\pi}{c_l}(A_1 - A_2)\right]$ |
| Loading | $\Delta C_p = 2\gamma/V_\infty$ |

## Implementation notes

- Any camber line can be used: `ThinAirfoil(camber_slope, breakpoints)` takes a callable
  $dz/dx$. `ThinAirfoil.from_naca4("2412")` wires in the NACA 4-digit camber line.
- The Fourier integrals are evaluated with adaptive quadrature (`scipy.integrate.quad`) in
  $\theta$. The NACA 4-digit camber line has a curvature jump at the maximum-camber location
  $p$, so $\theta_p = \arccos(1-2p)$ is passed as a breakpoint.

## Validation

| Check | Expected | Source |
|---|---|---|
| Symmetric section | $\alpha_{L0} = 0$, $c_{l\alpha} = 2\pi$, $c_{m,c/4} = 0$ | theory |
| Parabolic arc $z = 4hx(1-x)$ | $A_1 = 4h$, $\alpha_{L0} = -2h$, $c_{m,c/4} = -\pi h$ | closed form |
| NACA 2412 | $\alpha_{L0} = -2.077°$, $c_{m,c/4} = -0.053$ | Anderson, *Fundamentals of Aerodynamics*, Ex. 4.6–4.7 |
| Loading integrates back to $c_l$, $c_{m,LE}$ | to 1e-6 | consistency |

Measured NACA 2412 data give $\alpha_{L0} \approx -2.1°$ and $c_{m,c/4} \approx -0.045$: thin
airfoil theory gets the zero-lift angle right and slightly overestimates the nose-down moment
(the boundary layer de-cambers the real section).

## Limits

No thickness (so no leading-edge suction peak shape, no thickness effect on lift slope),
linearised boundary condition, inviscid. The [panel method](panel_method_2d.md) removes the
first two.
