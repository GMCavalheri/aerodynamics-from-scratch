# Boundary layer and profile drag

Code: [`src/aero/boundary_layer/`](../src/aero/boundary_layer) ·
Notebook: [`03_boundary_layer.ipynb`](../notebooks/03_boundary_layer.ipynb)

The panel method supplies the inviscid edge velocity $U_e(s)$ along each surface. Integral
boundary-layer methods march the momentum-integral equation

```math
\frac{d\theta}{ds} = \frac{c_f}{2} - (H + 2)\frac{\theta}{U_e}\frac{dU_e}{ds}
```

downstream from the stagnation point, closed with empirical correlations. There is no
viscous–inviscid interaction: the boundary layer does not feed back on $U_e$.

## 1. Laminar: Thwaites' method

Thwaites noticed that $\frac{U_e}{\nu}\frac{d\theta^2}{ds} \approx 0.45 - 6\lambda$ with
$\lambda = \frac{\theta^2}{\nu}\frac{dU_e}{ds}$, which integrates in closed form:

```math
\theta^2(s) = \frac{0.45\,\nu}{U_e^6}\int_0^s U_e^5\,ds .
```

$H(\lambda)$ and the shear parameter $\ell(\lambda) = \tau_w\theta/(\mu U_e)$ come from the
Cebeci–Bradshaw fits; laminar separation is predicted at $\lambda = -0.09$.

At a stagnation point $U_e \approx Ks$ and the integrand behaves like $s^5$, where plain
trapezoids are badly wrong (a factor 3 on the first interval). The code integrates $U_e^5$
exactly for piecewise-linear $U_e$ — $\int = \frac{h}{6}\sum_{k=0}^5 u_1^k u_2^{5-k}$ — which
reproduces the Hiemenz limit $\theta^2 = 0.075\ \nu/K$.

## 2. Transition: Michel's criterion

Transition is declared where

```math
Re_\theta \ge 1.174\left(1 + \frac{22400}{Re_x}\right)Re_x^{0.46},
```

or at laminar separation if that comes first (assumed to close as a short bubble). A
`x_trip` argument forces transition, like a trip strip.

## 3. Turbulent: Head's entrainment method

Head's entrainment equation with the mass-flow shape factor $H_1 = (\delta - \delta^{\ast})/\theta$,

```math
\frac{d}{ds}(U_e\theta H_1) = U_e\,0.0306\,(H_1 - 3)^{-0.6169},
```

plus the momentum integral and the Ludwieg–Tillmann friction law
$c_f = 0.246\cdot10^{-0.678H}Re_\theta^{-0.268}$, is integrated with `solve_ivp` (LSODA). The
$H_1(H)$ correlation is inverted in closed form. Turbulent separation is flagged at $H > 2.4$.

## 4. Drag: Squire–Young

The wake momentum deficit far downstream follows from the trailing-edge state:

```math
c_d = 2\,\theta_{TE}\left(\frac{U_{e,TE}}{V_\infty}\right)^{(H_{TE}+5)/2}
```

summed over both surfaces. The trailing edge of a finite-angle section is an inviscid
stagnation point, so the last few percent of chord see an unphysical deceleration; the march
stops at $x/c = 0.98$ (results change by about 4 % between 0.95 and 0.995).

## Validation

| Check | Result | Reference |
|---|---|---|
| Flat plate, laminar $\theta\sqrt{Re_x}/x$ | 0.671 | 0.664 (Blasius) |
| Howarth $U_e = 1 - x$ separation | $x = 0.123$ | 0.1199 (exact) |
| Hiemenz stagnation flow | $\theta^2 = 0.075\nu/K$ | exact |
| Turbulent flat plate $\theta$, $c_f$ | within 10 % | 1/7-power law |
| NACA 0012, $\alpha = 0$, Re = 6e6 | $c_d = 0.0068$ | ≈ 0.0058 (Abbott & von Doenhoff) |

The airfoil drag is about 15 % high: Michel's criterion puts transition near 0.2c on a NACA
0012 at Re = 6e6, earlier than an $e^N$ method would, so more of the surface is turbulent.
Trends (drag falls with Re, rises with $|\alpha|$, transition moves forward on the suction
side) are all reproduced.

## Limits

No viscous–inviscid coupling (no displacement effect on lift, no stall prediction, no
separation bubbles), empirical transition, and the trailing-edge treatment above. XFOIL-class
accuracy needs a coupled solver with an $e^N$ transition model.
