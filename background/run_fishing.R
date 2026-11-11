# Cod and sprat in a power-law background (see R/background_model.R and
# background_model.qmd): equilibrium branches under cod fishing, and the
# effect of a sprat fishery. This follows combined/fishing_branches.R.
#
# Run from the repository root, with BLAS limited to one thread per worker:
#   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript background/run_fishing.R
# It takes about 20 minutes on 4 cores and writes background/fishing.rds,
# which background_model.qmd reads.

library(mizer)
library(parallel)
options(mizer_info_level = 0)
source("R/build_model.R")
source("R/combined_model.R")
source("R/background_model.R")

cores <- as.integer(Sys.getenv("MC_CORES", min(8, detectCores())))

# Three cod states, at the first three points of background/run_cases.R:
#   * the coexistence state where it is stable and cod cannot invade;
#   * where it is the threshold, the state that cod reach from 20 % above it.
#     That state does not have power-law rates;
#   * the coexistence state where cod can invade.
cases <- data.frame(
    label = c("Coexistence state stable",
              "Threshold, cod settle above it",
              "No barrier"),
    phi_prey = c(0.96, 0.72, 0.36),
    phi_larv = c(0.36, 0.24, 0.48),
    cod_resource = c(FALSE, TRUE, TRUE),
    from_above = c(FALSE, TRUE, FALSE))
# All three are lost at a very low fishing mortality, so the sequences are
# fine near zero
F_cod_seq <- list(c(seq(0, 0.01, by = 0.001), 0.02, 0.03),
                  c(seq(0, 0.01, by = 0.001), 0.02, 0.03),
                  seq(0, 0.05, by = 0.005))
F_sprat_values <- c(0, 0.1, 0.25, 0.5, 1)

effort <- function(F_cod, F_sprat) c(Trawl = F_cod, `Sprat trawl` = F_sprat)

#' Project until every biomass changes by less than `tol` per year
#'
#' Cod below 1e-6 of `B_ref` and declining count as extinct.
settle_fishing <- function(params, effort, B_ref, t_round = 100,
                           max_rounds = 30, tol = 1e-5) {
    p <- params
    for (k in seq_len(max_rounds)) {
        sim <- project(p, effort = effort, t_max = t_round,
                       t_save = t_round / 10, progress_bar = FALSE)
        p <- finalParams(sim)
        b <- getBiomass(sim)
        nt <- nrow(b)
        rate <- log(b[nt, ] / b[nt - 1, ]) / (t_round / 10)
        extinct <- b[nt, ] < 1e-6 * B_ref & rate < 0
        if (all(is.na(rate) | abs(rate) < tol | extinct)) break
    }
    attr(p, "converged") <- k < max_rounds
    p
}

#' Follow the stable states along a sequence of cod fishing mortalities
#'
#' Each step starts from the state reached at the previous one. With `seed`
#' cod are topped up to at least that abundance before each step. Biomasses
#' are reported relative to `B_ref`.
trace_branch <- function(start, F_cod_seq, B_ref, seed = NULL) {
    p <- start
    rows <- vector("list", length(F_cod_seq))
    for (i in seq_along(F_cod_seq)) {
        if (!is.null(seed)) {
            n <- initialN(p)
            n["Cod", ] <- pmax(n["Cod", ], seed)
            initialN(p) <- n
        }
        p <- settle_fishing(p, effort(F_cod_seq[i], 0), B_ref)
        b <- getBiomass(p) / B_ref
        rows[[i]] <- data.frame(F_cod = F_cod_seq[i], Sprat = b[["Sprat"]],
                                Cod = b[["Cod"]],
                                converged = attr(p, "converged"))
    }
    do.call(rbind, rows)
}

sweeps <- c("Increasing F_cod (from the cod state)",
            "Decreasing F_cod (from the cod-free state + seed)")

#' The coexistence state, the cod state and the cod-free state of case i
case_states <- function(i) {
    coexistence <- new_background_params(cases$phi_prey[i], cases$phi_larv[i],
                                         cases$cod_resource[i])
    B_ref <- getBiomass(coexistence)
    cod <- coexistence
    if (cases$from_above[i]) {
        cod <- settle_fishing(background_project(coexistence, 1.2, 100),
                              effort(0, 0), B_ref)
    }
    list(coexistence = coexistence, cod = cod, B_ref = B_ref,
         no_cod = background_cod_free(coexistence))
}

jobs <- expand.grid(case = seq_len(nrow(cases)), dir = c("up", "down"),
                    stringsAsFactors = FALSE)
branches <- do.call(rbind, mclapply(seq_len(nrow(jobs)), function(k) {
    i <- jobs$case[k]
    s <- case_states(i)
    if (jobs$dir[k] == "up") {
        d <- trace_branch(s$cod, F_cod_seq[[i]], s$B_ref)
        d$sweep <- sweeps[1]
    } else {
        seed <- initialN(s$coexistence)["Cod", ] * 1e-6
        d <- trace_branch(s$no_cod, rev(F_cod_seq[[i]]), s$B_ref, seed = seed)
        d$sweep <- sweeps[2]
    }
    cbind(case = cases$label[i], d)
}, mc.cores = min(cores, nrow(jobs))))

# A sprat fishery, without cod fishing: the invasion fitness of cod in the
# cod-free state, and what becomes of the cod state
sprat_runs <- expand.grid(case = seq_len(nrow(cases)), F_sprat = F_sprat_values)
sprat_fishery <- do.call(rbind, mclapply(seq_len(nrow(sprat_runs)), function(k) {
    i <- sprat_runs$case[k]
    F_sprat <- sprat_runs$F_sprat[k]
    s <- case_states(i)
    no_cod <- settle_fishing(s$no_cod, effort(0, F_sprat), s$B_ref)
    cod <- settle_fishing(s$cod, effort(0, F_sprat), s$B_ref)
    initial_effort(no_cod) <- effort(0, F_sprat)
    data.frame(case = cases$label[i], F_sprat = F_sprat,
               R0 = cod_R0(no_cod) / cod_R0(s$coexistence),
               sprat_free = getBiomass(no_cod)[["Sprat"]] / s$B_ref[["Sprat"]],
               cod = getBiomass(cod)[["Cod"]] / s$B_ref[["Cod"]],
               sprat = getBiomass(cod)[["Sprat"]] / s$B_ref[["Sprat"]],
               converged = attr(cod, "converged"))
}, mc.cores = cores, mc.preschedule = FALSE))

saveRDS(list(cases = cases, branches = branches, sweeps = sweeps,
             sprat_fishery = sprat_fishery,
             mizer_version = as.character(packageVersion("mizer"))),
        "background/fishing.rds")

print(format(branches, digits = 3), row.names = FALSE)
print(format(sprat_fishery, digits = 3), row.names = FALSE)
