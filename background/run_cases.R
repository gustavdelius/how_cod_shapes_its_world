# Cod and sprat in a power-law background (see R/background_model.R and
# background_model.qmd): the computations at single points of the regime map.
#   * the size spectra with and without cod;
#   * cod started 20 % below and above the coexistence state;
#   * the same over 3000 years at two points of the scan that are slow;
#   * cod introductions into the cod-free state;
#   * the sensitivity to the parameters that set how the species respond;
#   * density dependence in the reproduction of cod.
#
# Run from the repository root, with BLAS limited to one thread per worker:
#   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript background/run_cases.R
# It takes about 35 minutes on 4 cores and writes background/cases.rds, which
# background_model.qmd reads.

library(mizer)
library(parallel)
options(mizer_info_level = 0)
source("R/build_model.R")
source("R/combined_model.R")
source("R/background_model.R")
source("background/run_scan.R")  # scan_point(); does not run the scan

# The points -----------------------------------------------------------------

# Four kinds of point on the regime map. `cod_resource` says whether cod feed
# on the plankton.
points <- data.frame(
    label = c("Coexistence state stable",
              "Threshold, cod settle above it",
              "No barrier",
              "Threshold, cod cycle above it"),
    phi_prey = c(0.96, 0.72, 0.36, 0.96),
    phi_larv = c(0.36, 0.24, 0.48, 0.72),
    cod_resource = c(FALSE, TRUE, TRUE, TRUE))

build <- function(i, ...) {
    new_background_params(points$phi_prey[i], points$phi_larv[i],
                          points$cod_resource[i], ...)
}

# Spectra with and without cod -------------------------------------------------

spectra <- do.call(rbind, mclapply(seq_len(nrow(points)), function(i) {
    p <- build(i)
    no_cod <- background_cod_free(p)
    present <- initialN(p)["Sprat", ] > 0
    plankton <- is_plankton(p) & initialNResource(p) > 0
    prey <- !is_plankton(p) & initialNResource(p) > 0
    rbind(
        data.frame(point = points$label[i], what = "Sprat",
                   w = w(p)[present],
                   ratio = (initialN(no_cod)["Sprat", ] /
                                initialN(p)["Sprat", ])[present]),
        data.frame(point = points$label[i], what = "Plankton",
                   w = w_full(p)[plankton],
                   ratio = (initialNResource(no_cod) /
                                initialNResource(p))[plankton]),
        data.frame(point = points$label[i], what = "Background prey",
                   w = w_full(p)[prey],
                   ratio = (initialNResource(no_cod) /
                                initialNResource(p))[prey]),
        data.frame(point = points$label[i], what = "Sprat feeding level",
                   w = w(p)[present],
                   ratio = getFeedingLevel(no_cod)["Sprat", present]))
}, mc.cores = cores))

# Cod started below and above the coexistence state ---------------------------

departure_runs <- expand.grid(point = seq_len(nrow(points)),
                              factor = c(0.8, 1.2))
departures <- do.call(rbind, mclapply(seq_len(nrow(departure_runs)), function(k) {
    i <- departure_runs$point[k]
    p <- build(i)
    B <- getBiomass(p)
    b <- attr(background_project(p, departure_runs$factor[k], t_max = 1000),
              "biomass")
    data.frame(point = points$label[i], factor = departure_runs$factor[k],
               time = as.numeric(rownames(b)), cod = b[, "Cod"] / B[["Cod"]],
               sprat = b[, "Sprat"] / B[["Sprat"]])
}, mc.cores = cores, mc.preschedule = FALSE))

# Two slow points of the scan --------------------------------------------------

# Both with cod feeding on the plankton. At the first, cod can invade and the
# coexistence state is weakly unstable. At the second, its leading eigenvalues
# are a complex pair with a small positive real part.
slow_runs <- expand.grid(point = 1:2, factor = c(0.8, 1.2))
slow_runs$phi_prey <- c(0.36, 0.6)[slow_runs$point]
slow_runs$phi_larv <- c(0.72, 0.24)[slow_runs$point]
slow <- do.call(rbind, mclapply(seq_len(nrow(slow_runs)), function(k) {
    g <- slow_runs[k, ]
    p <- new_background_params(g$phi_prey, g$phi_larv, TRUE)
    B <- getBiomass(p)
    b <- attr(background_project(p, g$factor, t_max = 3000), "biomass")
    data.frame(phi_prey = g$phi_prey, phi_larv = g$phi_larv,
               factor = g$factor, time = as.numeric(rownames(b)),
               cod = b[, "Cod"] / B[["Cod"]], sprat = b[, "Sprat"] / B[["Sprat"]])
}, mc.cores = cores, mc.preschedule = FALSE))

# Introductions of cod into the cod-free state --------------------------------

intro_factors <- c(0.01, 0.3, 1, 3)
intro_runs <- expand.grid(point = seq_len(nrow(points)),
                          factor = intro_factors)
intro <- do.call(rbind, mclapply(seq_len(nrow(intro_runs)), function(k) {
    i <- intro_runs$point[k]
    p <- build(i)
    B <- getBiomass(p)
    q <- background_cod_free(p)
    n <- initialN(q)
    n["Cod", ] <- initialN(p)["Cod", ] * intro_runs$factor[k]
    initialN(q) <- n
    sim <- project(q, t_max = 400, t_save = 1, progress_bar = FALSE)
    b <- getBiomass(sim)
    data.frame(point = points$label[i], factor = intro_runs$factor[k],
               time = getTimes(sim), cod = b[, "Cod"] / B[["Cod"]],
               sprat = b[, "Sprat"] / B[["Sprat"]])
}, mc.cores = cores, mc.preschedule = FALSE))

# Sensitivity ------------------------------------------------------------------

# At the second point: what sets how strongly the species respond
variants <- list(
    Baseline = list(),
    `Background level 0.2` = list(background_level = 0.2),
    `Background level 0.8` = list(background_level = 0.8),
    `Fixed background` = list(background_level = 1),
    `Background prey 10 times as abundant` = list(background_abundance = 10),
    `Plankton level 0.2` = list(resource_level = 0.2),
    `Plankton level 0.8` = list(resource_level = 0.8),
    `Feeding level 0.4` = list(f0 = 0.4),
    `Feeding level 0.8` = list(f0 = 0.8),
    `Plankton up to 1 g` = list(w_pp_cutoff = 1),
    `320 size classes` = list(no_w = 320))
sensitivity <- do.call(rbind, mclapply(names(variants), function(v) {
    cbind(variant = v,
          do.call(scan_point, c(list(points$phi_prey[2], points$phi_larv[2],
                                     points$cod_resource[2]),
                                variants[[v]])))
}, mc.cores = cores, mc.preschedule = FALSE))

# Density dependence in the reproduction of cod ---------------------------------

# A Beverton-Holt relation for cod with reproduction level L in the
# coexistence state. Cod are started 50 % below and above that state.
bh_runs <- expand.grid(phi_prey = c(0.48, 0.96), phi_larv = c(0.24, 0.72),
                       cod_resource = c(FALSE, TRUE), level = c(0.5, 0.8))
reproduction <- do.call(rbind, mclapply(seq_len(nrow(bh_runs)), function(k) {
    g <- bh_runs[k, ]
    p <- new_background_params(g$phi_prey, g$phi_larv, g$cod_resource,
                               cod_reproduction_level = g$level)
    B <- getBiomass(p)[["Cod"]]
    below <- background_project(p, 0.5)
    above <- background_project(p, 1.5)
    cbind(g, lead_re = Re(leading_eigenvalue(p)),
          R0 = cod_R0(background_cod_free(p)),
          below_cod = getBiomass(below)[["Cod"]] / B,
          above_cod = getBiomass(above)[["Cod"]] / B)
}, mc.cores = cores, mc.preschedule = FALSE))

saveRDS(list(points = points, spectra = spectra, departures = departures,
             slow = slow, intro = intro,
             sensitivity = sensitivity, reproduction = reproduction,
             mizer_version = as.character(packageVersion("mizer"))),
        "background/cases.rds")

print(format(reproduction, digits = 3), row.names = FALSE)
