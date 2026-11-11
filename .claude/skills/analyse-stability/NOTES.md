# Project notes: analyse-stability

- The cod–sprat state is a stable fixed point only in a narrow window. With
  cod `z0 = 0.8`, `z_ext = 0.8`, it is stable for cod `gamma` of about 0.13–0.17
  times mizer's default (1.07e-11 is used). It cycles above that, and cod can
  invade the sprat-only state above about 0.25 times the default. With adult
  cod mortality around 0.3/yr the cod state is a deep predator–prey cycle
  (period about 16 yr) instead of a fixed point.
- A cod Beverton–Holt ceiling also stabilises the cycle, but only at
  reproduction levels of about 0.8 or more, so it was not used.
- Raising the cod `beta` to 60–100 destroyed the bistable window (cycles, then
  the cod state is lost before the fold).
- Other stabilisers that were tried and rejected:
  - Cod cannibalism (interaction Cod→Cod = 1) gives a stable cod state at a
    cod gamma of 0.2× default, but cod can then invade the sprat-only state
    (R0 ≈ 4.6), so there is no bistability. At 0.12× default the cod state is
    lost, with cannibalism strength 0.3 or 1.
  - Feeding cod on the large end of the plankton spectrum
    (`interaction_resource` 0.5–2 for cod) as a depletable alternative food
    changes almost nothing. The stunted sprat deplete those plankton sizes too.
  - Less alternative food (`E_ext` 0 or 1 instead of 2.5): cod cannot persist
    at all.
  - A higher cod critical feeding level (`fc = 0.3`): the cod state cycles and
    is lost before any stable window opens.
- Imposing cod's predation-mortality profile on sprat as fixed `ext_mort`
  (no cod present) gives a stable sprat state for fractions up to about 0.3
  and at 1, but limit cycles for fractions of about 0.45–0.8.
- Under `second_order_w(params) <- c(flux = "van_leer")` the cod state is still
  stable (max real part about −0.014/yr, versus −0.13/yr upwind). Check it with
  `findSteadyState(solver = "newton")` + `getStability()`, not by projection;
  see the run-simulation notes for the van Leer NaN bug.

## Combined model with cod larvae (`R/combined_model.R`)

The combined model now uses the numerics of the cultivation-depensation
repository: second-order scheme (centred flux, bin averages), 320 size classes
and size diffusion with c = 0.1. Findings for that version:

- The cod state is much more strongly damped than under upwind: max real part
  −0.095/yr at the baseline (`w_shift` 0.2 g, θ 0) against −0.060/yr, and
  −0.04 to −0.17/yr up to θ = 0.02 for `w_shift` 0.1–0.6 g. Under upwind it
  cycled for `w_shift` ≥ 0.65 g.
- Size diffusion carries cod across the food bottleneck: the food barrier
  is never absolute (R0 in the cod-free state ≥ 0.27 at θ = 0, instead of 0),
  and it fails at `w_shift` ≈ 0.39 g instead of 0.46 g.
- Started from the anchor (θ = 0), Newton's method finds a **saddle** at
  θ ≥ 0.025 for some `w_shift`: a real eigenvalue of about +0.1/yr and up to
  30 % less cod (at 0.45 g, θ = 0.03: +0.109/yr). Starting from the solution at
  the previous θ instead finds the stable cod state (−0.13/yr). Isolated
  "unstable" points between stable neighbours are the giveaway. The driver
  therefore continues along θ for each `w_shift`.
- Cost of one point, single-threaded: Newton about 30 s, `getStability()`
  about 36 s, scheme-consistent R0 under 1 s.

Findings for the first, upwind version (150 classes, no physical diffusion):

- Modelling cod eggs and larvae explicitly (1 mg instead of 0.5 g, larval
  food up to `w_shift` = 0.2 g) halves the damping of the cod state: max real
  part −0.06/yr against −0.13/yr in the original model. The basin shrinks with
  it. Introductions of 30 % and 100 % of the cod state into the cod-free state
  overshoot and die out, whereas they establish in the original model.
- At `w_shift` = 0.2 g, sprat predation on cod larvae destabilises the cod
  state (θ = 0.02: −0.011/yr; θ = 0.05: cod collapse). At `w_shift` ≥ 0.45 g
  it stabilises it, and at `w_shift` 0.7–0.8 g it turns a cycling cod state
  (θ = 0) into a stable fixed point (θ ≥ 0.01–0.015).
- A later diet shift (`w_shift` 0.25–0.6 g, with `erepro` recalibrated) is a
  better food-leak knob than more adult alternative food. `E_ext` ≥ 4.2 makes
  the cod state cycle before the food barrier lets cod invade.
- In the combined model, higher cod `z0` (1.0–1.5) and lower cod `gamma`
  (0.85–0.9×) destabilise the cod state rather than stabilise it.
- `findSteadyState(solver = "newton")` occasionally stalls on this model
  (residual about 1e-2 instead of 1e-7, at isolated grid points). The
  eigenvalues reported there are meaningless. A second Newton pass from the
  stalled result fixes the isolated cases. Do not fall back to projection:
  at θ ≥ 0.035 the copied-in starting state lies outside the basin of a cod
  state that Newton finds to be stable, and cod die out.
- Near θ ≈ 0.04 (at `w_shift` = 0.6 g) the cod state reaches a fold. A real
  eigenvalue rises from −0.085/yr at θ = 0.038 to −0.0001/yr at θ = 0.040,
  and Newton stalls with residuals of 1e-4 to 1e-2. The prototype scan stops
  at θ = 0.035.

## Background model (`R/background_model.R`)

Cod and sprat in a power-law background, built with the mizerEcopath
workflow (`newAllometricParams()` → `setFeedingLevels()` →
`makeInteracting()` → `setResourceInteraction()`) and then
`set_background_response()`. The background mortality is fixed. The
background prey is depletable (semichemostat, `background_level` 0.5).
Upwind, 200 size classes.

- **The background prey is in the resource slot** (`n_pp` above
  `w_pp_cutoff`, the plankton below it), with custom `Encounter` and
  `ResourceMort` rate functions. That is what makes `getStability()`, the
  Newton solver and `cod_R0()` treat it as part of the state. A
  `setComponent()` component would be held fixed in the Jacobian.
- **`getStability()` cannot show how strongly a stable state is damped
  here.** `resource_level(p) <-` sets the replenishment rate proportional to
  the grazing, so plankton and background-prey sizes that no fish eats have
  rates down to zero. Each contributes a real eigenvalue −(rate + mortality):
  a ladder from about −1e-3/yr up to rounding error (±1e-13/yr). The largest
  can come out as +6e-14, so `$stable` is `FALSE` for a stable state. Count a
  state as unstable only if `max_real_part` > 1e-6, and measure damping from
  a projection (`return_rate()`). Matching the ladder against
  −(`resource_rate` + `getResourceMort`) and dropping those eigenvalues does
  not isolate the fish modes: weakly grazed sizes are shifted just enough to
  escape a 1e-3 relative tolerance.
- **Where cod can invade the cod-free state, the coexistence state is
  stable; behind a barrier it is mostly a saddle.** The scan (81 points per
  variant) has a barrier at 40 points with cod larvae eating plankton and at
  24 without. The coexistence state is stable at only 3 and 5 of those, in a
  narrow band next to R0 = 1. Elsewhere behind the barrier it has one real
  positive eigenvalue, up to 0.2/yr.
- Cod started 20 % above a saddle settle with 1.1–1.6 times the cod and
  fewer sprat, cycle (strong links), or die out (strongest links). They never
  drive sprat extinct or grow without bound: the background prey limits them.
- Two slow points, both with cod larvae eating plankton:
  - `phi_prey` 0.36, `phi_larv` 0.72: R0 = 1.03 and the coexistence state is
    unstable (0.005/yr). Cod cycle between 0.55 and 1.4 times their
    coexistence biomass with a period of about 2300 years.
  - `phi_prey` 0.6, `phi_larv` 0.24: a complex leading pair
    (0.0010 ± 0.0013i /yr). Cod started above pass the state and die out,
    which takes more than 1000 years. The scan counts cod below 1 % after
    1000 years as dying.
- **The abundance of the background prey (`background_abundance`) changes
  eigenvalues and basins but not steady states or R0.** Ten times the
  abundance makes the prey respond ten times as slowly: at `phi_prey` 0.72,
  `phi_larv` 0.24 the instability grows at 0.046/yr instead of 0.031/yr, and
  cod started 20 % above the state die out instead of settling above it.
- Background prey that only large cod eat replenishes slowly (0.006/yr at
  5 kg). After removing cod it is still recovering after 300 years, so
  `background_cod_free()` sets the resource to its steady state for the fish
  at the end of the projection.
- A cod Beverton–Holt relation with reproduction level L multiplies R0 in
  the cod-free state by 1/(1 − L). It does not stabilise the coexistence
  state at any of the four points tried where it is unstable.
- Newton's method stalls on the cod-free state (residual about 2e-4, cod
  exactly zero) and mostly on the states above a saddle (residual 1e-5 to
  1e-4). Those states are found by projection instead.
- Cost, single-threaded: about 0.1 s per projected year (three times the
  cost with a fixed background), so `scan_point()` takes 100–160 s.
