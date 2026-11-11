# Project notes: build-model

## Background model (`R/background_model.R`)

- The allometric workflow lives in two packages: `newAllometricParams()` and
  `setFeedingLevels()` in mizerEcopath, `makeInteracting()`,
  `setResourceInteraction()` and `alignResource()` in mizerFromEcopath, which
  `library(mizerEcopath)` attaches.
- `newAllometricParams()` sets the search volume coefficient `gamma` to 1,
  so the interaction matrix carries the absolute predation strength (values
  around 2e-12 for this model). `max_interaction()` computes the largest
  admissible value for each link; `phi_prey` and `phi_larv` are fractions of
  it.
- At exactly the largest cod-on-sprat strength the external encounter of cod
  is zero at one size, and `setResourceInteraction()` (default `tol = 1e-10`)
  then lets cod absorb no plankton at all. The scan therefore stops at
  φ = 0.96.
- With plankton up to 1 g, how much plankton cod can absorb depends on
  `phi_prey` (100 % of larval food at 0, 46 % at 0.95). With a cutoff of
  0.1 g it does not (99.7 % throughout), because the cod sizes that depend on
  sprat no longer eat plankton. The model uses 0.1 g.
- `setResourceInteraction()` clamps the resource level to just below 1 when
  it starts from constant resource dynamics. Set `resource_level(p) <-`
  afterwards.
- The pipeline does not work with size diffusion as it stands:
  `steadySingleSpecies()` with `D_ext` leaves a diffusive flux out of the
  largest size class of each species (residual 0.16/yr for sprat). The model
  uses the upwind scheme without diffusion.
- The built model is recalculation-safe: `species_params(p) <-
  species_params(p)` changes no rate. (`all.equal()` on the rate arrays
  reports a difference, but only in the `params` attribute they carry.)
- The background prey (`set_background_response()`) is stored in the resource
  spectrum above `w_pp_cutoff`, with a kernel matrix (rows species × size,
  columns `w_full`) in `other_params(p)$background` and the rate functions
  `background_encounter()` and `background_resource_mort()`.
  `resource_level(p) <-` takes a vector over `w_full`, so the plankton and
  the background prey can sit at different levels. `background_level = 1`
  skips the step and gives the fixed background.
- The model needs mizer ≥ 3.4.0.9003 and mizerEcopath 0.4.0 (branch
  `extract-mizerFromEcopath`) with mizerFromEcopath. With mizer 3.4.0.9001
  and mizerEcopath 0.3.2.9000, `setResourceInteraction()` fails with
  "non-conformable arrays".
- `setResourceInteraction()` warns that "the resource capacity has been set
  manually and so it was not rebalanced". It is harmless here:
  `resource_level(p) <-` sets the capacity afterwards.
- `readParams()` of mizer 3.4.0.9003 fails on `model/*.rds`, which were saved
  by a 3.4.0.9001 build that predates the `defaults_edition` slot:
  "defaults_edition must be a single number of at least 1". The background
  model does not read those files; `combined/run_prototype.R` and
  `combined/fishing_branches.R` do.
