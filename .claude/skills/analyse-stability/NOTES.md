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
- Under `second_order_w(params) <- c(flux = "van_leer")` the cod state is still
  stable (max real part about −0.014/yr, versus −0.13/yr upwind). Check it with
  `findSteadyState(solver = "newton")` + `getStability()`, not by projection;
  see the run-simulation notes for the van Leer NaN bug.
