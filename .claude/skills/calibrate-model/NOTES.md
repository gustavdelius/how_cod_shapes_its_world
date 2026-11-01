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
