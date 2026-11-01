# Prototype: the emergent Allee effect and cultivation-depensation in one
# mizer model (see R/combined_model.R and combined_model.qmd).
#
# Scans the plane of diet-shift size and larval predation, and simulates cod
# introductions. Run from the repository root, with BLAS limited to one thread
# per worker:
#   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript combined/run_prototype.R
# It takes about an hour on 8 cores and writes combined/results.rds, which
# combined_model.qmd reads.

library(mizer)
library(parallel)
options(mizer_info_level = 0)
source("R/build_model.R")
source("R/combined_model.R")

cores <- min(8, detectCores())

# The two states of the original model, carried over to the second-order
# numerics with size diffusion
anchor <- combined_anchor(readParams("model/cod_sprat_with_cod.rds"),
                          readParams("model/cod_sprat_no_cod.rds"))

states_at <- function(w_shift, theta) {
    combined_states(new_combined_params(theta = theta, w_shift = w_shift),
                    anchor)
}

# Regime map ---------------------------------------------------------------

#' Follow the cod state along increasing theta at fixed w_shift
#'
#' Each step starts Newton's method from the cod state of the previous step.
#' Started from the anchor instead, it finds a saddle at theta >= 0.025 for
#' some w_shift, with a positive real eigenvalue and up to 30 % less cod.
scan_chain <- function(w_shift, thetas) {
    ref <- anchor
    rows <- list()
    for (theta in thetas) {
        s <- combined_states(
            new_combined_params(theta = theta, w_shift = w_shift), ref)
        rows[[length(rows) + 1]] <- cbind(data.frame(
            w_shift = w_shift, theta = theta,
            erepro = species_params(s$cod)["Cod", "erepro"],
            cod_biomass = getBiomass(s$cod)[["Cod"]],
            residual = attr(s$cod, "convergence")$residual,
            max_real_part = getStability(s$cod)$max_real_part),
            invasion_decomposition(s$cod, s$no_cod))
        if (attr(s$cod, "convergence")$residual < 1e-4) ref$cod <- s$cod
    }
    do.call(rbind, rows)
}
# The scan stops at theta = 0.035. In the upwind version the cod state reached
# a fold near theta = 0.04.
scan <- do.call(rbind, mclapply(seq(0.1, 0.6, by = 0.05), scan_chain,
                                thetas = seq(0, 0.035, by = 0.005),
                                mc.cores = cores))
# Where Newton's method did not converge there is no cod state to report
scan$converged <- scan$residual < 1e-4
scan[!scan$converged, c("max_real_part", "R0", "food", "larval_predation",
                        "interaction")] <- NA

# Introductions of cod into the cod-free state ------------------------------

intro_points <- data.frame(w_shift = c(0.2, 0.3, 0.5, 0.5, 0.5),
                           theta = c(0, 0.01, 0, 0.01, 0.02))
intro_factors <- c(0.01, 0.3, 1)
intro_runs <- expand.grid(point = seq_len(nrow(intro_points)),
                          factor = intro_factors)
intro <- do.call(rbind, mclapply(seq_len(nrow(intro_runs)), function(i) {
    pt <- intro_points[intro_runs$point[i], ]
    fac <- intro_runs$factor[i]
    s <- states_at(pt$w_shift, pt$theta)
    p <- s$no_cod
    n <- initialN(p)
    n["Cod", ] <- initialN(s$cod)["Cod", ] * fac
    initialN(p) <- n
    sim <- project(p, t_max = 300, t_save = 1, method = "tr_bdf2",
                   progress_bar = FALSE)
    data.frame(w_shift = pt$w_shift, theta = pt$theta, factor = fac,
               time = getTimes(sim),
               cod = getBiomass(sim)[, "Cod"] / getBiomass(s$cod)[["Cod"]])
}, mc.cores = cores))

saveRDS(list(scan = scan, intro = intro,
             mizer_version = as.character(packageVersion("mizer"))),
        "combined/results.rds")

# Summary ------------------------------------------------------------------

print(format(scan[scan$w_shift %in% c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6) &
                      scan$theta %in% c(0, 0.01, 0.02, 0.035), ],
             digits = 3), row.names = FALSE)
