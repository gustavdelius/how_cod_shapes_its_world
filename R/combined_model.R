# A prototype that combines the two mechanisms that can stop cod recovering:
#   * the emergent Allee effect of Van Leeuwen, De Roos & Persson (2008),
#     which R/build_model.R reproduces: without cod, the stunted sprat
#     population produces too few small sprat for juvenile cod to grow on;
#   * cultivation-depensation (Walters & Kitchell 2001), as in the
#     cultivation-depensation repository: without cod, the abundant adult
#     sprat eat cod eggs and larvae.
#
# The combined model is the cod-sprat model of R/build_model.R with three
# changes:
#   * cod enter the model as 1 mg eggs instead of 0.5 g settled juveniles.
#     The cod energetics (h, ks) are those of the original model;
#   * cod larvae and small juveniles feed on an unmodelled larval food: a
#     fixed external encounter rate that falls off around the diet-shift size
#     `w_shift`. Above that size cod depend on sprat, as before. The larval
#     food is not depleted, so cod larvae do not compete with sprat;
#   * sprat eat cod eggs and larvae, with interaction strength `theta`.
#
# The original model folded cod egg and larval survival into erepro. That
# survival is now modelled, so cod erepro is recalibrated for every parameter
# set, to keep the cod state that of the original model. The cod-free state
# does not depend on any cod parameter, so it is also the original one.
#
# Two parameters move the model between the two mechanisms:
#   w_shift  how far cod can grow on larval food before they need sprat. The
#            later the shift, the weaker the food bottleneck.
#   theta    how readily sprat eat cod eggs and larvae.
#
# Source after R/build_model.R and R/analysis_functions.R.

#' External encounter rate of cod: the alternative food E_ext w^n of the
#' original model plus a larval food E_larv w^n that falls off around
#' `w_shift` with steepness `k`
combined_ext_encounter <- function(params, w_shift, E_larv, k = 6) {
    ee <- ext_encounter(params)
    E_ext <- species_params(params)["Cod", "E_ext"]
    n <- species_params(params)["Cod", "n"]
    ww <- w(params)
    ee["Cod", ] <- (E_ext + E_larv / (1 + (ww / w_shift)^k)) * ww^n
    ee
}

#' Create the combined MizerParams object (state not yet set)
#'
#' @param theta Interaction strength of sprat on cod: predation on cod eggs
#'   and larvae.
#' @param w_shift Size (g) around which the larval food of cod runs out.
#' @param f_larv Feeding level that the larval food alone gives.
#' @param k Steepness of the fall-off of the larval food.
new_combined_params <- function(theta = 0, w_shift = 0.2, f_larv = 0.6,
                                k = 6) {
    original <- species_params(new_cod_sprat_params())
    sp <- cod_sprat_species
    sp$w_min[sp$species == "Cod"] <- 0.001
    # Keep the energetics of the original model. Mizer would otherwise
    # re-derive h and ks of cod from the new w_min.
    sp$h <- original[sp$species, "h"]
    sp$ks <- original[sp$species, "ks"]
    inter <- cod_sprat_interaction
    inter["Sprat", "Cod"] <- theta
    params <- newMultispeciesParams(
        sp,
        interaction = inter,
        gear_params = cod_sprat_gears,
        initial_effort = 0,
        no_w = 150,
        info_level = 0)
    h <- original["Cod", "h"]
    ext_encounter(params) <- combined_ext_encounter(
        params, w_shift, E_larv = h * f_larv / (1 - f_larv), k = k)
    other_params(params)$larval_food <- list(w_shift = w_shift,
                                             f_larv = f_larv, k = k)
    setMetadata(
        params,
        title = "Cod-sprat model with cod larvae",
        description = paste(
            "Combines the emergent Allee effect of Van Leeuwen et al. (2008)",
            "with predation by sprat on cod eggs and larvae",
            "(cultivation-depensation)."))
}

#' Copy the abundances of `from` into `params`, which must have the same
#' size grid
copy_state <- function(params, from) {
    initialN(params) <- initialN(from)
    initialNResource(params) <- initialNResource(from)
    params
}

#' Fill in cod below 0.5 g, the recruitment size of the original model, with
#' the steady-state profile that matches the abundance at 0.5 g
fill_cod_larvae <- function(params) {
    i_settle <- sum(w(params) < 0.5) + 1
    N <- per_recruit_abundance(getEGrowth(params)["Cod", ],
                               getMort(params)["Cod", ], dw(params),
                               params@w_min_idx[["Cod"]])
    n <- initialN(params)
    larvae <- seq_len(i_settle - 1)
    n["Cod", larvae] <- N[larvae] * n["Cod", i_settle] / N[i_settle]
    initialN(params) <- n
    params
}

#' Rescale cod erepro so that cod are in steady state (R0 = 1) in the state
#' stored in `params`
calibrate_cod_erepro <- function(params) {
    sp <- given_species_params(params)
    sp["Cod", "erepro"] <- sp["Cod", "erepro"] / getR0(params)
    given_species_params(params) <- sp
    params
}

#' The cod state and the cod-free state of the combined model
#'
#' The cod state starts from the cod state of the original model, with cod
#' eggs and larvae filled in. Cod erepro is rescaled to make it a steady
#' state, and Newton's method then finds the fixed point. That converges even
#' where the fixed point is unstable. When Newton's method stalls, it is
#' restarted once from where it stopped. Check
#' `attr(cod, "convergence")$residual` before using the result. Projection is
#' no fallback: at strong larval predation the starting state lies outside the
#' basin of a cod state that is locally stable, and cod die out.
#'
#' @param original_cod,original_no_cod The two states of the original model.
combined_states <- function(params, original_cod, original_no_cod) {
    p <- copy_state(params, original_cod) |>
        fill_cod_larvae() |>
        calibrate_cod_erepro()
    newton <- function(p) suppressWarnings(
        findSteadyState(p, solver = "newton", info_level = 0))
    cod <- newton(p)
    if (attr(cod, "convergence")$residual > 1e-5) {
        cod <- newton(cod)
    }
    list(cod = cod, no_cod = copy_state(cod, original_no_cod))
}

#' Decompose the invasion fitness of cod in the cod-free state
#'
#' R0 of a rare cod in the cod-free state, relative to the cod state, where
#' R0 = 1. "food" takes only the growth and reproduction of cod from the
#' cod-free state, "larval_predation" only their mortality. "interaction" is
#' the factor by which R0 differs from the product of the two; it is below 1
#' when hungry juvenile cod linger in the sizes at which sprat eat them.
invasion_decomposition <- function(cod, no_cod) {
    r0_cod <- getR0(cod)
    out <- data.frame(
        R0 = getR0(no_cod) / r0_cod,
        food = getR0(no_cod, mortality_from = cod) / r0_cod,
        larval_predation = getR0(cod, mortality_from = no_cod) / r0_cod)
    out$interaction <- out$R0 / (out$food * out$larval_predation)
    out
}
