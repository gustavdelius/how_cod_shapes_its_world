# Helper functions for analysing the cod-sprat model. Source after
# R/build_model.R.

# Life stages, following the three stages of the stage-structured model in the
# paper: juveniles, small adults that still grow, and large adults.
stage_breaks <- list(
    Sprat = c(Juveniles = 0, `Small adults` = 5, `Large adults` = 7),
    Cod = c(Juveniles = 0, `Small adults` = 500, `Large adults` = 4000)
)

#' Biomass of each life stage of a species
#'
#' @param object A MizerParams (one value per stage) or MizerSim (a time x
#'   stage matrix).
stage_biomass <- function(object, species) {
    lower <- stage_breaks[[species]]
    upper <- c(lower[-1], Inf)
    out <- mapply(function(lo, up) {
        b <- getBiomass(object, min_w = lo, max_w = up)
        # A MizerSim gives a time x species array, a MizerParams a vector
        if (is.null(dim(b))) b[species] else b[, species]
    }, lower, upper)
    if (is.null(dim(out))) {
        out <- matrix(out, nrow = 1, dimnames = list(NULL, names(lower)))
    }
    out
}

#' Invasion fitness of a species: its lifetime reproductive output R0
#'
#' The expected number of recruits produced over its lifetime by one recruit
#' of `species` that experiences the growth and mortality rates of the current
#' state of `params`. When `species` is rare this is the criterion for
#' invasion: it can increase from low density only if R0 > 1. The steady-state
#' abundance per recruit uses the same first-order upwind discretisation as
#' mizer's projection, so R0 is 1 at a steady state of `project()` up to a
#' discretisation error. Density-dependent reproduction (R_max) is ignored,
#' which is exact at low density.
getR0 <- function(params, species = "Cod") {
    g <- getEGrowth(params)[species, ]
    mu <- getMort(params)[species, ]
    e_repro <- getERepro(params)[species, ]
    dw <- dw(params)
    i0 <- params@w_min_idx[[species]]
    N <- numeric(length(g))
    N[i0] <- 1 / (g[i0] + mu[i0] * dw[i0])
    for (j in (i0 + 1):length(g)) {
        N[j] <- g[j - 1] * N[j - 1] / (g[j] + mu[j] * dw[j])
    }
    erepro <- species_params(params)[species, "erepro"]
    0.5 * erepro * sum(e_repro * N * dw) / w(params)[i0]
}

#' Recruitment with a constant immigration term
#'
#' Beverton-Holt recruitment plus a constant immigration of recruits given in
#' `other_params(params)$immigration` (a named vector, recruits per year).
#' Used only in the fishing scenario, to give cod every chance to recover.
immigrationRDD <- function(rdi, species_params, params, ...) {
    BevertonHoltRDD(rdi = rdi, species_params = species_params, ...) +
        other_params(params)$immigration
}

#' Set a small constant immigration of cod recruits
set_cod_immigration <- function(params, rate) {
    other_params(params)$immigration <-
        setNames(c(0, rate), c("Sprat", "Cod"))[species_params(params)$species]
    setRateFunction(params, "RDD", "immigrationRDD")
}

#' Fishing effort vector with only the cod trawl active
cod_effort <- function(effort) c(Trawl = effort, `Sprat trawl` = 0)

#' Follow an attractor while the cod fishing effort is stepped through
#' `efforts`
#'
#' Each step starts from the state reached at the previous one. Returns a data
#' frame with the stage biomasses, the cod yield and the attractor type.
effort_branch <- function(params, efforts, t_max = 500) {
    rows <- list()
    for (effort in efforts) {
        initial_effort(params) <- cod_effort(effort)
        # A tighter tolerance than projectUntilSettled()'s default, as in
        # scanModel(), so that the points lie on a smooth curve
        sim <- projectUntilSettled(params, t_max = t_max, distance_tol = 0.001,
                                   progress_bar = FALSE, info_level = 0)
        conv <- attr(sim, "convergence")
        params <- finalParams(sim)
        initial_effort(params) <- cod_effort(effort)
        cod_stages <- stage_biomass(params, "Cod")
        sprat_stages <- stage_biomass(params, "Sprat")
        colnames(cod_stages) <- paste("Cod", colnames(cod_stages))
        colnames(sprat_stages) <- paste("Sprat", colnames(sprat_stages))
        rows[[length(rows) + 1]] <- data.frame(
            effort = effort,
            attractor = if (is.null(conv$attractor)) NA else conv$attractor,
            cod = getBiomass(params)[["Cod"]],
            sprat = getBiomass(params)[["Sprat"]],
            yield = getYield(params)[["Cod"]],
            cod_stages, sprat_stages,
            check.names = FALSE)
    }
    do.call(rbind, rows)
}
