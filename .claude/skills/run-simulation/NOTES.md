# Project notes: run-simulation

- **van Leer NaN bug (mizer 3.4.0.9000 master and 3.4.0.9001).** Projecting the
  cod state with `second_order_w(params) <- c(flux = "van_leer")` gives `NaN`
  for the whole Sprat spectrum after about 240–256 simulated years, with Euler
  or second-order time stepping. Cause: in `flux_limiter_chi()` (R/transport.R)
  the guard `ifelse(den == 0, 0, num / den)` misses subnormal denominators, so
  `r = ±Inf` and van Leer gives `Inf/Inf = NaN`. Restarting from
  `finalParams()` does not reproduce it, because the limiter is averaged with
  `attr(n, "chi")` from the previous step. Reported as
  https://github.com/sizespectrum/mizer/issues/630; check whether it is fixed
  before relying on van Leer projections here.
- Fishing ramp: F rising to 0.6/yr over 40 years is too brief above the
  persistence limit (about 0.39/yr) to collapse cod, which recover when fishing
  stops. 0.8/yr collapses them for good.
- `scanEffort()` with no gear sets every gear. The model has two gears, so use
  `scanEffort("Trawl")`, and give `effort` as a named vector or a time × gear
  array.
- Combined model (second order, 320 classes, size diffusion): a 200-year
  `project(method = "tr_bdf2", dt = 0.1)` takes about 20 s. One grid point
  (two Newton solves, `getStability()`, four scheme-consistent R0 solves)
  takes a few minutes. Run parallel workers with `OPENBLAS_NUM_THREADS=1
  OMP_NUM_THREADS=1`, or the BLAS threads of 8 workers compete.
