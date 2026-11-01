# Equilibrium branches of the combined model under cod fishing, with and
# without a sprat fishery: the analogue of Figure 3 of paper.Rmd in the
# cultivation-depensation repository (see combined_model.qmd).
#
# Run from the repository root:
#   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript combined/fishing_branches.R
# It takes about 15 minutes on 4 cores and writes
# combined/fishing_branches.rds, which combined_model.qmd reads.

library(mizer)
library(parallel)
options(mizer_info_level = 0)
source("R/build_model.R")
source("R/combined_model.R")

# The model point: both barriers act (food factor 0.69, larval-predation
# factor 0.53 without fishing)
w_shift <- 0.3
theta <- 0.01
# Sprat fishing in the right panel. Above about 0.28 it lets cod invade the
# cod-free state.
F_sprat_values <- c(0, 0.5)
F_cod_seq <- seq(0, 0.4, by = 0.01)

effort <- function(F_cod, F_sprat) c(Trawl = F_cod, `Sprat trawl` = F_sprat)

#' Project until every biomass changes by less than `tol` per year
#'
#' Uses mizer's first-order (Euler) time step with dt = 0.1, which leaves the
#' fixed points unchanged. Larger steps are unstable in this model: Euler with
#' dt = 0.5 has a one-step spectral radius of 1.08 in the cod-free state, and
#' TR-BDF2 with dt = 0.5 drifts away from a stable cod state after about 150
#' years. Cod below 1e-6 of `B_ref` and declining count as extinct.
settle <- function(params, effort, B_ref, t_round = 100, max_rounds = 30,
                   tol = 1e-5) {
    p <- params
    for (k in seq_len(max_rounds)) {
        sim <- project(p, effort = effort, t_max = t_round, dt = 0.1,
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
#' cod are topped up to at least that abundance before each step, so that
#' they can invade wherever their invasion fitness exceeds 1. Biomasses are
#' reported relative to `B_ref`.
trace_branch <- function(start, F_sprat, F_cod_seq, B_ref, seed = NULL) {
    p <- start
    rows <- vector("list", length(F_cod_seq))
    for (i in seq_along(F_cod_seq)) {
        if (!is.null(seed)) {
            n <- initialN(p)
            n["Cod", ] <- pmax(n["Cod", ], seed)
            initialN(p) <- n
        }
        p <- settle(p, effort(F_cod_seq[i], F_sprat), B_ref)
        b <- getBiomass(p) / B_ref
        rows[[i]] <- data.frame(F_sprat = F_sprat, F_cod = F_cod_seq[i],
                                Sprat = b[["Sprat"]], Cod = b[["Cod"]],
                                converged = attr(p, "converged"))
    }
    do.call(rbind, rows)
}

compute_branches <- function() {
    anchor <- combined_anchor(readParams("model/cod_sprat_with_cod.rds"),
                              readParams("model/cod_sprat_no_cod.rds"))
    states <- combined_states(
        new_combined_params(theta = theta, w_shift = w_shift), anchor)
    B_ref <- getBiomass(states$cod)
    seed <- initialN(states$cod)["Cod", ] * 1e-6
    jobs <- expand.grid(F_sprat = F_sprat_values, dir = c("up", "down"),
                        stringsAsFactors = FALSE)
    branches <- do.call(rbind, mclapply(seq_len(nrow(jobs)), function(i) {
        F_sprat <- jobs$F_sprat[i]
        if (jobs$dir[i] == "up") {
            d <- trace_branch(states$cod, F_sprat, F_cod_seq, B_ref)
            d$sweep <- "Increasing F_cod (from the cod state)"
        } else {
            start <- settle(states$no_cod, effort(0, F_sprat), B_ref)
            d <- trace_branch(start, F_sprat, rev(F_cod_seq), B_ref,
                              seed = seed)
            d$sweep <- "Decreasing F_cod (from the cod-free state + seed)"
        }
        d
    }, mc.cores = 4))

    # Invasion fitness of a rare cod in the cod-free state. That state does
    # not depend on cod fishing, so it is settled once per sprat fishing level.
    invasion <- do.call(rbind, lapply(F_sprat_values, function(F_sprat) {
        p <- settle(states$no_cod, effort(0, F_sprat), B_ref)
        R0 <- function(F_cod) {
            initial_effort(p) <- effort(F_cod, F_sprat)
            cod_R0(p)
        }
        R0_0 <- R0(0)
        data.frame(F_sprat = F_sprat, R0_unfished_cod = R0_0,
                   F_inv = if (R0_0 > 1) uniroot(function(f) R0(f) - 1,
                                                 c(0, 2))$root else NA)
    }))

    res <- list(branches = branches, invasion = invasion,
                w_shift = w_shift, theta = theta, B_ref = B_ref,
                mizer_version = as.character(packageVersion("mizer")))
    saveRDS(res, "combined/fishing_branches.rds")
    res
}

res <- compute_branches()
branches <- res$branches
invasion <- res$invasion

# Summary --------------------------------------------------------------------

sweeps <- c("Increasing F_cod (from the cod state)",
            "Decreasing F_cod (from the cod-free state + seed)")
# Cod above this fraction of their unfished biomass count as abundant
abundant <- 0.1
print(invasion)
for (F_sprat in F_sprat_values) {
    up <- branches[branches$F_sprat == F_sprat & branches$sweep == sweeps[1], ]
    last <- max(up$F_cod[up$Cod > abundant])
    cat(sprintf(paste("F_sprat = %.2f: cod state persists up to F_cod = %.2f",
                      "(cod at %.2f of unfished), collapses at %.2f\n"),
                F_sprat, last, up$Cod[up$F_cod == last],
                min(up$F_cod[up$F_cod > last])))
}
cat("Steps that did not settle:\n")
print(branches[!branches$converged, c("F_sprat", "F_cod", "sweep", "Cod")],
      row.names = FALSE)
