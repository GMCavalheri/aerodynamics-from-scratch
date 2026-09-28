# Lesson 5 · Planform design: taper, aspect ratio and sweep

> We now have a tool that computes the lift and induced drag of any wing in a fraction of a
> second. Engineers use exactly this kind of tool in conceptual design. In this lesson we act as
> the designer: we ask *which shape* and answer with numbers.

**Code:** [`src/aero/vortex_lattice_3d/studies.py`](../../src/aero/vortex_lattice_3d/studies.py),
[`validation/planform_studies.py`](../../validation/planform_studies.py) ·
**Notebook:** [`05_planform_studies.ipynb`](../../notebooks/05_planform_studies.ipynb) ·
**Previous:** [Lesson 4](04-vortex-lattice-method.md) · **Next:** [Lesson 6](06-verification-and-validation.md)

## Learning objectives

1. Write induced drag in terms of span loading and explain why **span**, not area, is the lever.
2. Relate planform shape to span loading, section lift and stall behaviour.
3. Find the best taper ratio for an unswept wing and explain why it moves with sweep.
4. Explain what sweep does to lift slope and span loading, and why swept wings need washout.
5. Connect these trends to the proportions of real aircraft.

---

## 1. The design equation

From Lesson 4, $C_{D_i} = C_L^2/(\pi AR\,e)$. Write it in dimensional form with $L = qSC_L$ and
$AR = b^2/S$:

$$
\boxed{\;D_i = \frac{L^2}{q\,\pi\,b^2\,e}\;}
\qquad q = \tfrac12\rho V_\infty^2 .
$$

Read it carefully:

- For a given weight ($L = W$) and speed, induced drag depends only on the **span $b$** and the
  **span efficiency $e$**. Wing area does not appear.
- $b$ enters **squared**. Span is the big lever; planform shape (through $e$) is the fine
  adjustment, typically worth a few percent.
- Induced drag scales like $1/q$, so it dominates at **low speed** (take-off, climb, loiter).

That is why gliders have huge spans and why winglets, which increase effective span, pay off.

## 2. Planform → loading → stall

Four planforms of the same aspect ratio (8), analysed with the VLM:

![Planform loading](../figures/planform_loading.png)

**Left: span loading $c\,c_l$.** This is what the wake sees and what sets induced drag. The
elliptic wing is elliptic ($e = 0.999$). The rectangular wing carries too much load near the
tips ($e = 0.972$). The pointed wing ($\lambda = 0$) piles the load inboard but ends with a
steep drop; its $e$ of 0.879 is the worst.

**Right: section lift $c_l/C_L$.** This is what each airfoil section feels, and **where it
peaks is where stall starts**:

| Planform | $c_l$ peak | Stall starts | Consequence |
|---|---|---|---|
| Rectangular | at the root | root | benign: ailerons keep working; the classic trainer wing |
| Elliptic | uniform | everywhere at once | abrupt stall, no warning |
| Tapered λ = 0.4 | mid-span | mid/outboard | acceptable with a little washout |
| Pointed λ = 0 | at the tip | tip | dangerous: roll-off, loss of aileron control |

So the best-drag planform (elliptic) and the safest planform (rectangular) are different
shapes. Design is the compromise between them.

## 3. Taper ratio

Tapering a rectangular wing moves area, and therefore load, inboard. Too little taper leaves the
tips overloaded; too much starves them. There is an optimum:

![Taper study](../figures/planform_taper.png)

| Method | Optimum λ | $e_{max}$ |
|---|---|---|
| VLM (this code) | 0.45 | 0.996 |
| Lifting line (Glauert) | 0.375 | 0.988 |
| Classical textbook value | ≈ 0.35–0.4 | — |

Three lessons from this plot:

1. **A tapered wing gets within half a percent of elliptic.** Elliptic planforms (the Spitfire)
   are expensive to build for almost no gain. Straight taper is the practical answer.
2. **The optimum is flat.** Any λ between 0.3 and 0.65 is within 0.5 % of the best $e$, so taper
   is chosen for structure (a deeper root spar) and stall behaviour as much as for drag.
3. **The two theories agree on the trend, not the numbers.** Lifting line treats each section as
   a 2D airfoil. Near a wide tip that over-loads the tip, so it penalises low taper more.

### Sweep moves the optimum

| Sweep $\Lambda_{c/4}$ | Optimum λ (VLM, AR 8) |
|---|---|
| 0° | 0.45 |
| 15° | 0.30 |
| 30° | 0.23 |
| 45° | 0.15 |

Aft sweep shifts load outboard (section 5). Designers counter it with **more taper**, which is
why swept transport wings have taper ratios around 0.25–0.3 instead of 0.4.

## 4. Aspect ratio

![Aspect ratio study](../figures/planform_aspect_ratio.png)

- **Left:** at fixed $C_L$, induced drag falls like $1/AR$. Going from AR 6 to AR 12 halves it.
- **Right:** the penalty *relative to elliptic*, $1/e - 1$. It grows with AR for the rectangular
  wing (0.3 % at AR 3 → 10.5 % at AR 20). High-aspect-ratio wings are more sensitive to getting
  the loading right, which is why sailplanes use carefully tailored multi-taper planforms.

Why not AR 30 on an airliner? **Structure.** Wing root bending moment grows with span, and so
does wing weight. The optimum aspect ratio balances induced drag against structural weight, and
it has crept up (≈ 8–9 on older jets, 9.5–11 on current composite wings) as composites made long,
thin wings affordable.

## 5. Sweep

Transport aircraft sweep their wings for a **compressibility** reason that potential flow does
not model: by **simple sweep theory**, only the velocity component normal to the leading edge,
$V_\infty\cos\Lambda$, drives the pressure distribution, which delays the drag rise to a higher
Mach number. Sweep has low-speed costs, and the VLM shows them clearly:

![Sweep study](../figures/planform_sweep.png)

- **Lift slope falls** roughly like $\cos\Lambda$ (left). A swept wing needs a higher angle of
  attack, which is one reason airliners land nose-high with large flaps.
- **Span efficiency falls** for large sweep in either direction (centre).
- **Aft sweep moves load outboard** (right): the tip $c_l/C_L$ goes from 0.85 at 0° to 1.07 at
  45°. The result is **tip stall**, and because the tips of an aft-swept wing are behind the centre of gravity, losing tip
  lift pitches the nose **up**, deeper into the stall. This is the swept-wing **pitch-up**
  problem, fixed with washout, more taper, wing fences, vortex generators and leading-edge slats.
- **Forward sweep moves load inboard** (right, Λ = −30°). Root stall first is benign, which is why
  the X-29 and the HFB 320 Hansa Jet flew forward-swept wings. The catch is aeroelastic
  divergence: bending twists the tips nose-up, which needs stiff composite skins.

## 6. Putting it together: why wings look the way they do

| Aircraft type | Typical choices | Driven by |
|---|---|---|
| Trainer (e.g. Cessna 150/152 class) | rectangular or mildly tapered, AR ≈ 7, washout | cost, benign root-first stall |
| Airliner | AR ≈ 9–11, λ ≈ 0.25–0.3, Λ ≈ 25–35°, washout | Mach, induced drag vs weight, tip-stall control |
| Sailplane | AR 20–30+, multi-taper near-elliptic | induced drag dominates at low speed |
| Fighter | low AR, high sweep or delta | supersonic drag, manoeuvre, not induced drag |

Every row is an application of the equations and figures in this lesson.

## 7. Worked example: span beats shape

Two wings carry the same lift at the same speed. Wing A is rectangular, AR 8. Wing B is a
perfect ellipse with 5 % **less span**. Which has less induced drag?

$$
\frac{D_{i,B}}{D_{i,A}} = \frac{b_A^2\,e_A}{b_B^2\,e_B} = \frac{0.972}{0.95^2 \times 1.000} = 1.077 .
$$

The "perfect" elliptic wing has **7.7 % more** induced drag. Five percent of span outweighs the
whole planform effect.

## 8. Exercises

1. Using the VLM, find the taper ratio that maximises $e$ for AR 6 and AR 12. Does the optimum
   depend on aspect ratio?
2. A designer adds winglets that increase the *effective* span by 4 % at the same $e$. By how much
   does induced drag drop?
3. For an AR 8, λ = 0.4 wing with 30° of aft sweep, how much linear washout is needed to bring the
   tip $c_l/C_L$ back to the unswept value? What happens to $e$?
4. Why does the elliptic-wing Spitfire not appear in the table in section 6 as "optimal"? List
   two reasons its planform was chosen that are not about induced drag.

<details>
<summary>Answers</summary>

1. The optimum moves only weakly: λ = 0.475 at AR 6, 0.45 at AR 8, 0.425 at AR 12. The
   within-0.5 % band is 0.3–0.75 at AR 6 but narrows to 0.3–0.55 at AR 12, since high-AR wings
   are more sensitive to the loading.
2. $D_i \propto 1/b^2$: $1/1.04^2 = 0.925$, so about **7.5 % less**. That is why winglets are worth
   their weight on long-range aircraft.
3. About 2° of linear washout. At α = 6° the outboard $c_l/C_L$ drops from 0.98 to 0.86 (the
   unswept wing has 0.85), and $e$ at that $C_L$ even rises, from 0.984 to 0.995. At other lift
   coefficients the same twist costs $e$ (Lesson 4, exercise 4). Try it with
   `Wing(..., twist=lambda eta: -eps * eta)`.
4. Thin root with room for retracted undercarriage and eight wing guns; good handling (the
   Spitfire also had washout). Its elliptic shape gave a small aerodynamic benefit at high cost.

</details>

## Key takeaways

- $D_i = L^2/(q\pi b^2 e)$: span is the lever, $e$ is the fine tuning.
- Planform sets span loading (drag) *and* section lift distribution (where it stalls).
- Straight taper λ ≈ 0.35–0.45 gets within 0.5 % of elliptic; the optimum is flat.
- Aft sweep loads the tips, lowers the optimum taper and causes pitch-up; washout and taper fix it.
- Real wings trade induced drag against structure, cost, Mach number and stall safety.

## Further reading

- Raymer, *Aircraft Design: A Conceptual Approach*, Ch. 4.
- Anderson, *Fundamentals of Aerodynamics*, Sec. 5.3.
- Küchemann, *The Aerodynamic Design of Aircraft*, Pergamon, 1978.
- Kroo, "Drag due to lift: concepts for prediction and reduction", *Annu. Rev. Fluid Mech.* 33, 2001.
