# Cod and sprat in a power-law background (see R/background_model.R and
# background_model.qmd): the invasion fitness of cod in the cod-free state and
# the stability of the coexistence state, over the strengths of the two
# predation links, with and without cod feeding on the plankton. The
# background prey is at its default level.
#
# Run from the repository root, with BLAS limited to one thread per worker:
#   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript background/run_scan.R
# It takes about an hour and a half on 4 cores and writes
# background/scan.rds, which background_model.qmd reads.

library(mizer)
library(parallel)
options(mizer_info_level = 0)
source("R/build_model.R")
source("R/combined_model.R")
source("R/background_model.R")

# Set the environment variable MC_CORES to use fewer workers
cores <- as.integer(Sys.getenv("MC_CORES", min(8, detectCores())))

#' Abundance of the background prey in `state` relative to `params`: the
#' smallest and the largest ratio over its sizes (NA if there is none)
background_prey_range <- function(state, params) {
    prey <- !is_plankton(params) & initialNResource(params) > 0
    if (!any(prey)) return(c(NA, NA))
    range((initialNResource(state) / initialNResource(params))[prey])
}

#' Everything the regime map needs at one point
#'
#' The coexistence state is a steady state by construction. Its stability is
#' read from the leading eigenvalue, and cod are started 20 % below and 20 %
#' above it to see where they go. Where the state is unstable the run from
#' above lasts 1000 years, and a state it reaches that still has both species
#' is polished with Newton's method.
scan_point <- function(phi_prey, phi_larv, cod_resource, ...) {
    p <- new_background_params(phi_prey, phi_larv, cod_resource, ...)
    B <- getBiomass(p)
    lead <- leading_eigenvalue(p)
    unstable <- Re(lead) > 1e-6
    no_cod <- background_cod_free(p)
    below <- background_project(p, 0.8)
    above <- background_project(p, 1.2, t_max = if (unstable) 1000 else 300)
    b_above <- attr(above, "biomass")
    last <- tail(seq_len(nrow(b_above)), 100)
    prey_free <- background_prey_range(no_cod, p)
    row <- cbind(
        data.frame(phi_prey = phi_prey, phi_larv = phi_larv,
                   cod_resource = cod_resource,
                   lead_re = Re(lead), lead_im = Im(lead),
                   sprat_free = getBiomass(no_cod)[["Sprat"]] / B[["Sprat"]],
                   sprat_free_drift = biomass_drift(no_cod)[["Sprat"]],
                   background_free_min = prey_free[1],
                   background_free_max = prey_free[2]),
        background_decomposition(p, no_cod),
        below_cod = getBiomass(below)[["Cod"]] / B[["Cod"]],
        below_sprat = getBiomass(below)[["Sprat"]] / B[["Sprat"]],
        below_return = return_rate(below, B[["Cod"]]),
        above_cod = getBiomass(above)[["Cod"]] / B[["Cod"]],
        above_sprat = getBiomass(above)[["Sprat"]] / B[["Sprat"]],
        above_return = return_rate(above, B[["Cod"]]),
        above_cod_min = min(b_above[last, "Cod"]) / B[["Cod"]],
        above_cod_max = max(b_above[last, "Cod"]) / B[["Cod"]],
        above_residual = NA, above_lead_re = NA, above_lead_im = NA)
    if (unstable && all(getBiomass(above) / B > 1e-6)) {
        polished <- newton_state(above)
        row$above_residual <- attr(polished, "convergence")$residual
        if (row$above_residual < 1e-5) {
            lead_above <- leading_eigenvalue(polished)
            row$above_cod <- getBiomass(polished)[["Cod"]] / B[["Cod"]]
            row$above_sprat <- getBiomass(polished)[["Sprat"]] / B[["Sprat"]]
            row$above_lead_re <- Re(lead_above)
            row$above_lead_im <- Im(lead_above)
        }
    }
    row
}

if (sys.nframe() == 0) {
    phis <- seq(0, 0.96, by = 0.12)
    grid <- expand.grid(phi_prey = phis, phi_larv = phis,
                        cod_resource = c(TRUE, FALSE))
    scan <- do.call(rbind, mclapply(seq_len(nrow(grid)), function(i) {
        scan_point(grid$phi_prey[i], grid$phi_larv[i], grid$cod_resource[i])
    }, mc.cores = cores, mc.preschedule = FALSE))
    saveRDS(list(scan = scan,
                 mizer_version = as.character(packageVersion("mizer")),
                 mizerEcopath_version =
                     as.character(packageVersion("mizerEcopath"))),
            "background/scan.rds")
    print(format(scan, digits = 3), row.names = FALSE)
}
