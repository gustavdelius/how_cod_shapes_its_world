# Project notes: set-up-fishing

- The model has two gears, both at effort 0 in the saved states: "Trawl" on
  cod (sigmoid length, l50 = 35 cm) and "Sprat trawl" on sprat (knife edge at
  1 g). Give effort as a named vector or a time × gear array.
- Sprat trawl and cod invasion fitness R0 in the cod-free state:
  - Knife edge 1 g: F = 0.25 → R0 1.18, F = 0.5 → R0 1.55 (fixed points).
    F ≥ 1 sends sprat into limit cycles, and R0 at one moment means nothing.
  - Knife edge 3 g: F = 0.25 → R0 0.2 only. F ≥ 0.5 → sprat cycles.
- Recovery by prey thinning (cod-free state with 1 % cod immigration, sprat
  trawl at F = 0.5 from year 10): 10, 15 or 20 years are not enough, 30 years
  are. F = 0.3 for 50 years also works.

## Background model (`R/background_model.R`, `background/run_fishing.R`)

- The model carries the two gears of `R/build_model.R` (`gear_params(p) <-
  cod_sprat_gears` at the end of `new_background_params()`), both at effort 0.
- The cod states are far more fragile than in the combined model (fold at
  0.27/yr there):
  - Where the coexistence state is stable (`phi_prey` 0.96, `phi_larv` 0.36,
    cod larvae not eating plankton) it persists at F_cod = 0.001/yr and
    collapses at 0.002/yr. Use a fishing grid that is fine near zero.
  - The cod state above the threshold at `phi_prey` 0.72, `phi_larv` 0.24
    (cod larvae eating plankton) persists at 0.025/yr and collapses at
    0.05/yr.
  - Cod never return on the downward sweep.
- A sprat fishery (knife edge at 1 g) cuts both ways: it removes the sprat
  that eat cod larvae and the adult sprat that cod eat. R0 of cod in the
  cod-free state peaks at 0.99 (F = 0.25) at the first point and never
  reaches 1. At the second it is 1.20 at F = 0.5 and 1.06 at F = 1. Cod
  biomass in the cod state falls with sprat fishing in both cases, and cod
  die out at F = 0.5 at the first point.
