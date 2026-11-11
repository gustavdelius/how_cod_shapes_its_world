# Cod and sprat embedded in a power-law background.
#
# The combined model (R/combined_model.R) has only cod, sprat and plankton, so
# each species is the other's whole world. Here the two species live in a
# community that is not modelled: in the coexistence state every individual
# meets food at a rate proportional to w^n and dies at a rate proportional to
# w^(n - 1), as it would in a community with a Sheldon spectrum. Only part of
# that food and that mortality is attributed to the modelled species and to
# the plankton resource. The rest is the background. Its mortality is fixed.
# Its food is a spectrum of other prey that is depleted by what cod and sprat
# eat of it. The question is whether the two mechanisms that can stop cod
# recovering survive this dilution:
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
# A last step, set_background_response(), turns the external encounter that
# is left into feeding on a spectrum of background prey above the largest
# plankton size.
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
# The background prey shares mizer's resource spectrum with the plankton: the
# plankton lies below `w_pp_cutoff` and the background prey above it. Both are
# semichemostats, each with its own level (abundance relative to carrying
# capacity) in the coexistence state. The lower the level, the more the
# abundance gives way under predation. Because the background prey is part of
# the resource, getStability(), the Newton solver and cod_R0() cover it. With
# a background level of 1 the background food is fixed.
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
#' @param background_level The same for the background prey. At 1 the
#'   background food is fixed.
#' @param background_abundance Abundance of the background prey as a multiple
#'   of the power law of the plankton. It does not change how far the
#'   background prey gives way, only how fast: the more of it there is, the
#'   lower the mortality that cod and sprat cause, and the rate at which it
#'   replenishes is set by that mortality.
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
                                  background_level = 0.5,
                                  background_abundance = 1,
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
    other_params(p)$background <- list(
        phi_prey = phi_prey, phi_larv = phi_larv,
        cod_resource = cod_resource, background_level = background_level,
        theta_max = theta)
    level <- resource_level
    if (background_level < 1) {
        p <- set_background_response(p, background_abundance)
        level <- ifelse(is_plankton(p), resource_level, background_level)
    }
    # Sets the carrying capacities so that the plankton and the background
    # prey replenish at the rate at which they are now eaten
    resource_level(p) <- level

    if (cod_reproduction_level > 0) {
        # Raises erepro so that the coexistence state stays a steady state
        p <- setBevertonHolt(p, reproduction_level =
                                 c(Sprat = 0, Cod = cod_reproduction_level))
    }
    gear_params(p) <- cod_sprat_gears
    setMetadata(
        p,
        title = "Cod and sprat in a power-law background",
        description = paste(
            "Cod and sprat with allometric encounter and mortality in the",
            "coexistence state, part of which is predation between them."))
}

#' Which sizes of the resource spectrum are plankton
#'
#' The rest, above the largest plankton size, is background prey.
is_plankton <- function(params) {
    w_full(params) <= resource_params(params)$w_pp_cutoff
}

#' Let the background food give way under predation
#'
#' Turns external encounter into feeding on a spectrum of background prey:
#' the power law of the plankton, continued above the largest plankton size.
#' The background prey is held in the resource spectrum, above the plankton.
#'
#' Of the external encounter of a predator, the part that responds is at most
#' what a Sheldon spectrum would supply from prey above the largest plankton
#' size. So small larvae, which could reach such prey only with the tail of
#' their predation kernel, keep a fixed food. The predator feeds on the
#' background prey with an availability that depends on its size and is
#' chosen so that its encounter rate does not change.
#'
#' Sets the rate functions background_encounter() and
#' background_resource_mort(). The carrying capacity of the background prey
#' still has to be set, with `resource_level<-`.
#'
#' @param abundance Abundance of the background prey as a multiple of the
#'   power law of the plankton.
set_background_response <- function(params, abundance = 1) {
    plankton <- is_plankton(params)
    rp <- resource_params(params)
    sheldon <- rp$kappa * w_full(params)^(-rp$lambda)
    energy <- w_full(params) * dw_full(params)

    # Encounter per unit of availability, from background prey and from a
    # whole Sheldon spectrum
    kernel <- getPredKernel(params)
    per_availability <- function(n_prey) {
        search_vol(params) * apply(kernel, c(1, 2), function(k)
            sum(k * n_prey * energy))
    }
    from_background <- per_availability(sheldon * !plankton)
    share <- from_background / per_availability(sheldon)
    share[!is.finite(share)] <- 0
    responsive <- pmin(ext_encounter(params), getEncounter(params) * share)
    availability <- responsive / (abundance * from_background)
    availability[!is.finite(availability)] <- 0

    # Predation kernel on the background prey, with a row for each species
    # and predator size in the order of a species x size array
    kernel <- sweep(kernel, c(1, 2), availability * search_vol(params), "*")
    kernel[, , plankton] <- 0
    dim(kernel) <- c(prod(dim(kernel)[1:2]), dim(kernel)[3])

    initialNResource(params) <- initialNResource(params) * plankton +
        abundance * sheldon * !plankton
    ext_encounter(params) <- ext_encounter(params) - responsive
    other_params(params)$background <- c(
        other_params(params)$background,
        list(plankton = plankton, prey_energy = energy * !plankton,
             kernel = kernel, availability = availability))
    params <- setRateFunction(params, "Encounter", "background_encounter")
    setRateFunction(params, "ResourceMort", "background_resource_mort")
}

#' Encounter rate of the background prey
#'
#' A species x size array.
background_prey_encounter <- function(params, n_pp) {
    bg <- other_params(params)$background
    encounter <- bg$kernel %*% (as.numeric(n_pp) * bg$prey_energy)
    dim(encounter) <- dim(bg$availability)
    encounter
}

#' Encounter rate, with the plankton and the background prey in the resource
#' spectrum
#'
#' Replaces mizerEncounter(), which is left to deal with the fish prey, the
#' plankton and the fixed external encounter.
background_encounter <- function(params, n, n_pp, n_other, t = 0, ...) {
    bg <- other_params(params)$background
    mizerEncounter(params, n = n, n_pp = n_pp * bg$plankton,
                   n_other = n_other, t = t, ...) +
        background_prey_encounter(params, n_pp)
}

#' Mortality of the plankton and of the background prey
#'
#' Replaces mizerResourceMort(), which is left to deal with the plankton. The
#' mortality of the background prey removes what the predators eat of it.
background_resource_mort <- function(params, n, n_pp, n_other, t = 0,
                                     pred_rate, ...) {
    bg <- other_params(params)$background
    encounter <- background_encounter(params, n = n, n_pp = n_pp,
                                      n_other = n_other, t = t, ...)
    feeding_level <- mizerFeedingLevel(params, n = n, n_pp = n_pp,
                                       n_other = n_other, t = t,
                                       encounter = encounter, ...)
    predators <- sweep((1 - feeding_level) * n, 2, dw(params), "*")
    mizerResourceMort(params, n = n, n_pp = n_pp, n_other = n_other, t = t,
                      pred_rate = pred_rate, ...) * bg$plankton +
        as.vector(crossprod(bg$kernel, as.vector(predators)))
}

#' Mortality of the background prey caused by each species
#'
#' A species x size array over the sizes of the resource spectrum, for the
#' state stored in `params`.
background_prey_mort <- function(params) {
    predators <- sweep((1 - getFeedingLevel(params)) * initialN(params), 2,
                       dw(params), "*")
    kernel <- other_params(params)$background$kernel
    dim(kernel) <- c(dim(predators), ncol(kernel))
    mort <- apply(kernel, 3, function(k) rowSums(k * predators))
    dimnames(mort) <- list(sp = rownames(predators), w = w_full(params))
    mort
}

#' Shares of food and mortality that the modelled world accounts for
#'
#' For the state stored in `params`, by species and size: the share of the
#' encounter rate that comes from the other species ("prey"), from the
#' plankton ("resource") and from the background prey ("background"), and the
#' share of the mortality that is predation. The remainder is fixed.
background_shares <- function(params) {
    total <- getEncounter(params)
    responsive <- !is.null(other_params(params)$background$kernel)
    background <- 0 * total
    if (responsive) {
        background[] <- background_prey_encounter(params,
                                                  initialNResource(params))
    }
    p <- params
    p@ext_encounter[] <- 0
    prey <- getEncounter(p, n_pp = 0 * initialNResource(p))
    list(prey = prey / total,
         resource = (total - prey - background - ext_encounter(params)) /
             total,
         background = background / total,
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

#' The cod-free state: cod are removed and the sprat, the plankton and the
#' background prey settle
#'
#' Background prey that only large cod eat replenishes very slowly once the
#' cod are gone, and is still on its way after `t_max` years. The resource is
#' therefore put at the steady state that goes with the fish at the end of
#' the projection.
background_cod_free <- function(params, t_max = 300) {
    state <- background_project(params, factor = 0, t_max = t_max)
    biomass <- attr(state, "biomass")
    rate <- resource_rate(state)
    steady <- rate * resource_capacity(state) /
        (rate + getResourceMort(state))
    n_pp <- initialNResource(state)
    n_pp[is.finite(steady)] <- steady[is.finite(steady)]
    initialNResource(state) <- n_pp
    structure(state, biomass = biomass)
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
#' and the background prey replenish at a rate proportional to the rate at
#' which they are eaten, so plankton too small and background prey too large
#' for any fish to eat hardly change at all. Those size classes contribute a ladder of real eigenvalues from about -1e-3/yr up to
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
#'   background        the background prey that cod eat
#'   food              all three: growth and reproduction of cod
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
    plankton <- is_plankton(cod)
    plankton_only <- cod
    initialNResource(plankton_only) <- ifelse(
        plankton, initialNResource(no_cod), initialNResource(cod))
    background_only <- cod
    initialNResource(background_only) <- ifelse(
        plankton, initialNResource(cod), initialNResource(no_cod))
    out <- data.frame(
        R0 = cod_R0(no_cod) / r0,
        prey = cod_R0(sprat_only, mortality_from = cod) / r0,
        competition = cod_R0(plankton_only, mortality_from = cod) / r0,
        background = cod_R0(background_only, mortality_from = cod) / r0,
        food = cod_R0(no_cod, mortality_from = cod) / r0,
        larval_predation = cod_R0(cod, mortality_from = no_cod) / r0)
    out$interaction <- out$R0 / (out$food * out$larval_predation)
    out
}
