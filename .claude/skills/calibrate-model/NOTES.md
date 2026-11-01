# Project notes: calibrate-model

- The model (`R/build_model.R`) is deliberately **not** tuned with
  `tuneSteadyState()`. That function holds births fixed, which cuts the very
  feedback the model is about: cod predation raising sprat reproduction. Both
  steady states are found with fixed parameters by `projectUntilSettled()`.
  Tested: from the cod-free state plus 3e13 g of cod, `tuneSteadyState()` kept
  sprat recruitment at its cod-free value, settled on only 1.25e11 g of cod and
  set cod `erepro` to 1.7e34. Projection reaches 8.67e12 g. Reported upstream
  for the skills as https://github.com/sizespectrum/mizerAgents/issues/12 and
  for the help page as https://github.com/sizespectrum/mizer/issues/631.
- The cod-free state is reached by setting cod abundance to zero. The cod state
  is reached by adding cod to the cod-free state at 3e13 g (mizer's default
  initial cod spectrum shape, rescaled). 1e13 g dies out. Introductions far
  above equilibrium can overshoot, eat out the sprat and also die out.
- `findSteadyState(solver = "newton")` takes about 15 s per call on this
  model and agrees with the projected states to well under 1 %.
- `projectUntilSettled()` with its default `distance_tol` stops about 0.5 %
  short of the fixed point. Use `distance_tol = 0.001` (as `scanModel()` does)
  when points must lie on a smooth curve.
- To hold the cod state fixed while changing cod egg/larval parameters (the
  combined model, `R/combined_model.R`): with `R_max = Inf` every steady state
  has R0 = 1, so `erepro` selects the state. Copy the anchor state in,
  rescale `erepro` by 1/R0 there, then polish with
  `findSteadyState(solver = "newton")`. Cod biomass stays within about 2 %.
  The cod-free state needs no recalculation: it does not depend on any cod
  parameter. With size diffusion, R0 must come from the discretised
  transport equation (`cod_R0()`), not from the upwind recursion of `getR0()`.
- Carrying the original states over to the second-order grid
  (`combined_anchor()`): interpolate, then calibrate `erepro` and run Newton
  *twice*. The first Newton run from the interpolated cod state stalls
  (residual 0.06). Newton from the interpolated cod-free state stalls too
  (residual 0.02, and a second pass does not help); 200 years of
  `project(method = "tr_bdf2")` first, then Newton, gives 5e-9.
