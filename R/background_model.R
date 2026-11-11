# Cod and sprat embedded in a fixed power-law background.
#
# The combined model (R/combined_model.R) has only cod, sprat and plankton, so
# each species is the other's whole world. Here the two species live in a
# community that is not modelled: in the coexistence state every individual
# meets food at a rate proportional to w^n and dies at a rate proportional to
# w^(n - 1), as it would in a community with a Sheldon spectrum. Only part of
# that food and that mortality is attributed to the modelled species and to
# the plankton resource. The rest is a fixed background that does not respond
# when cod or sprat change. The question is whether the two mechanisms that
# can stop cod recovering survive this dilution:
#   * the emergent Allee effect of Van Leeuwen, De Roos & Persson (2008):
#     without cod there is less sprat food for cod;
#   * cultivation-depensation (Walters & Kitchell 2001): without cod the
#     more abundant sprat eat more cod eggs and larvae.
#
# The model is built with the allometric workflow of mizerEcopath:
#   newAllometricParams()     non-interacting species with power-law
#                             encounter and mortality, at steady state;
#   setFeedingLevels()        satiation and metabolism (h, ks) at unchanged
#                             growth, so that the species respond to food;
#   makeInteracting()         cod eat sprat and sprat eat cod eggs and larvae.
#                             The predation that this creates is subtracted
#                             from the external encounter and mortality;
#   setResourceInteraction()  feeding on a dynamic plankton resource replaces
#                             as much of the external encounter as it can.
# None of these steps changes a rate in the coexistence state, so that state
# is the same steady state for every strength of the interactions. No
# recalibration of erepro is needed to hold it fixed, unlike in the combined
# model.
#
# The strength of each predation link is given as a fraction `phi` of the
# largest strength for which the background stays non-negative:
#   phi_prey  cod eating sprat. At 1 the cod size that depends most on sprat
#             (about 300 g) gets all its food from sprat.
#   phi_larv  sprat eating cod eggs and larvae. At 1 the cod size that sprat
#             eat most (about 0.03 g) dies only from sprat predation.
# At 0 the species see only the background.
#
# The numerics are those the workflow produces: mizer's first-order upwind
# scheme without size diffusion, as in R/build_model.R. With size diffusion
# steadySingleSpecies() leaves an unbalanced diffusive flux out of the
# largest size class of each species, which the combined model avoids with a
# senescence mortality that would break the power law.
#
# Source after R/build_model.R (for the fishing gears) and
# R/combined_model.R (for cod_R0()). Needs the mizerEcopath package.

library(mizerEcopath)

# Life histories as in R/build_model.R, with both species entering the model
# as 1 mg eggs. The biomasses are those of the cod state of that model.
background_species <- data.frame(
    species = c("Sprat", "Cod"),
    w_max = c(15, 10000),
    w_mat = c(5, 500),
    w_min = 0.001,
    age_mat = c(1.5, 4),
    beta = c(1000, 30),
    sigma = c(2, 1.5),
    alpha = 0.6,
    biomass_observed = c(1.57e13, 8.67e12),
    a = 0.01,
    b = 3
)

#' Largest interaction strengths that leave the background non-negative
#'
#' For a non-interacting model with external encounter and mortality only.
#' Entry (i, j) is the strength of predation by i on j at which, at some
#' size, either all food of i comes from j or all mortality of j comes from i.
#' Predation is linear in the interaction strength while the feeding levels
#' are fixed, so it is measured at a negligible strength and scaled up.
max_interaction <- function(params) {
    eps <- 1e-20
    present <- initialN(params) > 0
    out <- interaction_matrix(params)
    for (i in rownames(out)) for (j in colnames(out)) {
        p <- params
        p@interaction[] <- 0
        p@interaction[i, j] <- eps
        mort <- getPredMort(p)[j, ] / eps
        p@ext_encounter[] <- 0
        encounter <- getEncounter(p)[i, ] / eps
        food_share <- (encounter / ext_encounter(params)[i, ])[present[i, ]]
        mort_share <- (mort / ext_mort(params)[j, ])[present[j, ]]
        out[i, j] <- 1 / max(food_share, mort_share)
    }
    out
}

#' Create the model in its coexistence state
#'
#' @param phi_prey Strength of cod predation on sprat, as a fraction of the
#'   largest admissible strength.
#' @param phi_larv Strength of sprat predation on cod eggs and larvae, as a
#'   fraction of the largest admissible strength.
#' @param cod_resource Whether cod feed on the plankton resource. If not,
#'   that part of their food stays in the fixed background, so that cod
#'   larvae do not compete with sprat.
#' @param resource_level Plankton abundance relative to its carrying capacity
#'   in the coexistence state. The lower it is, the more the plankton
#'   responds to grazing.
#' @param f0,fc Feeding level in the coexistence state and critical feeding
#'   level.
#' @param w_pp_cutoff Largest plankton size (g).
#' @param no_w Number of size classes.
#' @param cod_reproduction_level Reproduction level of cod in the coexistence
#'   state: the ratio of their recruitment to the maximum R_max of a
#'   Beverton-Holt relation. At the default of 0 there is no such relation
#'   (R_max = Inf) and all density dependence is emergent.
new_background_params <- function(phi_prey = 0.5, phi_larv = 0.5,
                                  cod_resource = TRUE, resource_level = 0.5,
                                  f0 = 0.6, fc = 0.2, w_pp_cutoff = 0.1,
                                  no_w = 200, cod_reproduction_level = 0) {
    p <- newAllometricParams(background_species, no_w = no_w)
    p <- setFeedingLevels(p, f = f0, f_c = fc)

    # Just short of the largest strength, where rounding could make the
    # background negative
    theta <- max_interaction(p) * (1 - 1e-9)
    inter <- theta
    inter[] <- 0
    inter["Cod", "Sprat"] <- phi_prey * theta["Cod", "Sprat"]
    inter["Sprat", "Cod"] <- phi_larv * theta["Sprat", "Cod"]
    p <- makeInteracting(p, inter)

    # alignResource() sets the capacity by hand and mizer warns that it has
    # not rebalanced it. setResourceInteraction() does that.
    p <- suppressWarnings(alignResource(p, w_pp_cutoff = w_pp_cutoff))
    background_food <- ext_encounter(p)
    p <- setResourceInteraction(p,
                                resource_dynamics = "resource_semichemostat")
    if (!cod_resource) {
        # Hand the plankton that cod eat back to the background
        sp <- species_params(p)
        sp["Cod", "interaction_resource"] <- 0
        species_params(p, recalculate = FALSE) <- sp
        food <- ext_encounter(p)
        food["Cod", ] <- background_food["Cod", ]
        ext_encounter(p) <- food
    }
    # Sets the carrying capacity so that the plankton replenishes at the
    # rate at which it is now eaten
    resource_level(p) <- resource_level

    if (cod_reproduction_level > 0) {
        # Raises erepro so that the coexistence state stays a steady state
        p <- setBevertonHolt(p, reproduction_level =
                                 c(Sprat = 0, Cod = cod_reproduction_level))
    }
    gear_params(p) <- cod_sprat_gears
    other_params(p)$background <- list(
        phi_prey = phi_prey, phi_larv = phi_larv,
        cod_resource = cod_resource, theta_max = theta)
    setMetadata(
        p,
        title = "Cod and sprat in a power-law background",
        description = paste(
            "Cod and sprat with allometric encounter and mortality in the",
            "coexistence state, part of which is predation between them."))
}

#' Shares of food and mortality that the modelled world accounts for
#'
#' For the state stored in `params`, by species and size: the share of the
#' encounter rate that comes from the other species ("prey") and from the
#' plankton ("resource"), and the share of the mortality that is predation.
background_shares <- function(params) {
    total <- getEncounter(params)
    p <- params
    p@ext_encounter[] <- 0
    with_resource <- getEncounter(p)
    p@species_params$interaction_resource <- 0
    prey <- getEncounter(p)
    list(prey = prey / total,
         resource = (with_resource - prey) / total,
         predation = getPredMort(params) / getMort(params))
}

#' The state with the abundance of cod multiplied by `factor`, projected for
#' `t_max` years
#'
#' The yearly biomasses are kept as attribute "biomass".
background_project <- function(params, factor = 1, t_max = 300) {
    n <- initialN(params)
    n["Cod", ] <- n["Cod", ] * factor
    initialN(params) <- n
    sim <- project(params, t_max = t_max, t_save = 1, progress_bar = FALSE)
    structure(finalParams(sim), biomass = getBiomass(sim))
}

#' The cod-free state: cod are removed and the sprat and the plankton settle
background_cod_free <- function(params, t_max = 300) {
    background_project(params, factor = 0, t_max = t_max)
}

#' Relative rate of change (1/yr) of the biomasses over the last ten years of
#' a state returned by background_project()
biomass_drift <- function(state) {
    b <- attr(state, "biomass")
    nt <- nrow(b)
    log(b[nt, ] / b[nt - 10, ]) / 10
}

#' Leading real part of the eigenvalues of the Jacobian at a steady state
#'
#' Returned as a complex number. A positive real part is the growth rate of
#' the instability. A stable state cannot be told from this number how
#' strongly it is damped: setting the resource level makes the plankton
#' replenish at a rate proportional to the rate at which it is eaten, so
#' plankton too small for any fish to eat hardly changes at all. Those size
#' classes contribute a ladder of real eigenvalues from about -1e-3/yr up to
#' rounding error (+-1e-13/yr), which hides the leading fish mode and makes
#' getStability() report `stable = FALSE`. Count a state as unstable only if
#' the real part is well above rounding error, and measure the damping of a
#' stable state from a projection (return_rate()).
leading_eigenvalue <- function(params) {
    ev <- getStability(params)$eigenvalues
    ev[which.max(Re(ev))]
}

#' Rate (1/yr) at which cod return to their biomass `B_cod` in a state
#' returned by background_project(), between years `from` and `to`
#'
#' Negative when cod move away.
return_rate <- function(state, B_cod, from = 50, to = 150) {
    deviation <- abs(attr(state, "biomass")[, "Cod"] / B_cod - 1)
    -log(deviation[to + 1] / deviation[from + 1]) / (to - from)
}

#' Decompose the invasion fitness of cod in the cod-free state
#'
#' R0 of a rare cod in the cod-free state, relative to the coexistence state,
#' where R0 = 1. Each factor takes one thing from the cod-free state and the
#' rest from the coexistence state:
#'   prey              the sprat that cod eat (the route of Van Leeuwen et al.)
#'   competition       the plankton that cod eat
#'   food              both: growth and reproduction of cod
#'   larval_predation  the mortality of cod (cultivation-depensation)
#' "interaction" is the factor by which R0 differs from food x larval
#' predation. It is below 1 when hungry cod linger at the sizes at which
#' sprat eat them.
background_decomposition <- function(cod, no_cod) {
    r0 <- cod_R0(cod)
    sprat_only <- cod
    n <- initialN(cod)
    n["Sprat", ] <- initialN(no_cod)["Sprat", ]
    initialN(sprat_only) <- n
    plankton_only <- cod
    initialNResource(plankton_only) <- initialNResource(no_cod)
    out <- data.frame(
        R0 = cod_R0(no_cod) / r0,
        prey = cod_R0(sprat_only, mortality_from = cod) / r0,
        competition = cod_R0(plankton_only, mortality_from = cod) / r0,
        food = cod_R0(no_cod, mortality_from = cod) / r0,
        larval_predation = cod_R0(cod, mortality_from = no_cod) / r0)
    out$interaction <- out$R0 / (out$food * out$larval_predation)
    out
}
