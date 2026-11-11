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
- The cod states withstand very little cod fishing (the combined model's fold
  is at 0.27/yr). Use a fishing grid that is fine near zero.
  - Where the coexistence state is stable behind a barrier (`phi_prey` 0.96,
    `phi_larv` 0.36, cod larvae not eating plankton) it persists at
    F_cod = 0.002/yr and collapses at 0.003/yr.
  - The cod state above the threshold at `phi_prey` 0.72, `phi_larv` 0.24
    (cod larvae eating plankton) persists at 0.003/yr and collapses at
    0.004/yr. With a fixed background it held to 0.025/yr, because it lay
    further above the threshold (1.63 instead of 1.23 times the coexistence
    biomass).
  - Cod never return on the downward sweep in these two cases.
  - Without a barrier (`phi_prey` 0.36, `phi_larv` 0.48, cod larvae eating
    plankton) there is no hysteresis. Cod fall steadily (0.74 at 0.01/yr,
    0.23 at 0.025/yr) and are gone at 0.04/yr. At 0.03 and 0.035/yr
    `settle_fishing()` does not converge within 3000 years.
  - Reason for the fragility: food is the only density dependence of cod,
    and the background prey at level 0.5 can at most double. R0 of a rare
    cod is at most 3.06 in the scan.
- A sprat fishery (knife edge at 1 g) opens the cod-free state: R0 of cod
  there exceeds 1 from F = 0.25 at the first point (1.96 at F = 0.5) and
  from F = 0.1 at the second (5.3 at F = 1). It also takes the adult sprat
  that cod eat: cod biomass in the cod state falls to 0.16 and 0.49 at F = 1.
  Without a barrier cod gain (1.43 times at F = 0.25).
- `background/run_fishing.R` takes about 20 minutes on 4 cores.
