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
# The numerics follow the cultivation-depensation repository rather than the
# original model:
#   * mizer's second-order scheme: centred flux and bin-averaged rates, on 320
#     size classes;
#   * size diffusion D(w) = D_ext w^(n + 1) for individual variability in
#     growth, with D_ext = c (alpha f0 h - ks), so that a juvenile cohort at
#     feeding level f0 = 0.6 gains variance of ln w at c = 0.1 per unit of ln w
#     grown;
#   * size grids up to twice w_inf, and a senescence mortality
#     z_sen (w / w_inf)^k_sen that stops diffusion carrying adults beyond w_inf;
#   * R0 from the steady state of the discretised transport equation itself
#     (cod_R0()), not from the upwind recursion of getR0(), which is exact
#     only for the first-order upwind scheme without diffusion.
#
# The original model folded cod egg and larval survival into erepro. That
# survival is now modelled, so cod erepro is recalibrated for every parameter
# set to keep the cod state of the original model, transferred to the new
# numerics (combined_anchor()). The cod-free state does not depend on any cod
# parameter, so it is the same throughout.
#
# Two parameters move the model between the two mechanisms:
#   w_shift  how far cod can grow on larval food before they need sprat. The
#            later the shift, the weaker the food bottleneck.
#   theta    how readily sprat eat cod eggs and larvae.
#
# cod_R0() uses mizer internals, as the cultivation-depensation repository
# does; it was written against mizer 3.4.0.9001. An exported replacement is
# requested in https://github.com/sizespectrum/mizer/issues/632.
#
# Source after R/build_model.R.

#' Bin averages of f(w) over the size classes of `params`, as the
#' second-order scheme needs them (midpoint rule in log w, `m` points per bin)
bin_average <- function(params, f, m = 16) {
    wl <- w(params)
    dlog <- log((wl + dw(params)) / wl)
    s <- (seq_len(m) - 0.5) / m
    vapply(seq_along(wl), function(i) {
        x <- wl[i] * exp(s * dlog[i])
        sum(f(x) * x) * dlog[i] / m / dw(params)[i]
    }, numeric(1))
}

#' External encounter rate of cod: the alternative food E_ext w^n of the
#' original model plus a larval food E_larv w^n that falls off around
#' `w_shift` with steepness `k`
combined_ext_encounter <- function(params, w_shift, E_larv, k = 6) {
    ee <- ext_encounter(params)
    E_ext <- species_params(params)["Cod", "E_ext"]
    n <- species_params(params)["Cod", "n"]
    ee["Cod", ] <- bin_average(params, function(x)
        (E_ext + E_larv / (1 + (x / w_shift)^k)) * x^n)
    ee
}

#' Create the combined MizerParams object (state not yet set)
#'
#' @param theta Interaction strength of sprat on cod: predation on cod eggs
#'   and larvae.
#' @param w_shift Size (g) around which the larval food of cod runs out.
#' @param f_larv Feeding level that the larval food alone gives.
#' @param k Steepness of the fall-off of the larval food.
#' @param no_w Number of size classes.
#' @param c_diff Size diffusion: growth in the variance of ln w per unit of
#'   ln w grown by a juvenile at feeding level 0.6.
#' @param z_sen,k_sen Senescence mortality z_sen (w / w_inf)^k_sen.
new_combined_params <- function(theta = 0, w_shift = 0.2, f_larv = 0.6,
                                k = 6, no_w = 320, c_diff = 0.1,
                                z_sen = 0.5, k_sen = 10) {
    original <- species_params(new_cod_sprat_params())
    sp <- cod_sprat_species
    sp$w_min[sp$species == "Cod"] <- 0.001
    sp$w_max <- 2 * sp$w_inf
    # Keep the energetics of the original model. Mizer would otherwise
    # re-derive h and ks of cod from the new w_min.
    sp$h <- original[sp$species, "h"]
    sp$ks <- original[sp$species, "ks"]
    # Size diffusion as in the cultivation-depensation repository, with
    # mizer's default reference feeding level f0 = 0.6
    sp$D_ext <- c_diff * (0.6 * 0.6 * sp$h - sp$ks)  # alpha = 0.6
    inter <- cod_sprat_interaction
    inter["Sprat", "Cod"] <- theta
    params <- newMultispeciesParams(
        sp,
        interaction = inter,
        gear_params = cod_sprat_gears,
        initial_effort = 0,
        no_w = no_w,
        second_order_w = TRUE,
        info_level = 0)
    second_order_w(params) <- "centred"
    w_inf <- species_params(params)$w_inf
    ext_mort(params) <- ext_mort(params) +
        t(vapply(w_inf, function(wi)
            bin_average(params, function(x) z_sen * (x / wi)^k_sen),
            numeric(length(w(params)))))
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

#' Interpolate the abundances of `from` onto the size grid of `params`
#' (linearly in log abundance and log size)
interpolate_state <- function(params, from) {
    interp <- function(x, y, x_new) {
        pos <- y > 0
        out <- exp(approx(log(x[pos]), log(y[pos]), log(x_new))$y)
        out[is.na(out)] <- 0
        out
    }
    n <- initialN(params)
    for (s in rownames(n)) {
        n[s, ] <- 0
        if (any(initialN(from)[s, ] > 0)) {
            n[s, ] <- interp(w(from), initialN(from)[s, ], w(params))
        }
    }
    initialN(params) <- n
    initialNResource(params) <- interp(w_full(from), initialNResource(from),
                                       w_full(params))
    params
}

#' Lifetime reproductive output R0 of a cod recruit
#'
#' Holds growth, mortality, reproductive investment and diffusion at their
#' values in the state stored in `params`, with the mortality taken from
#' `mortality_from` instead if given. Taking growth and mortality from
#' different states separates the effect of food from that of predation.
#' Solves the steady-state transport equation of the model's own size
#' discretisation for the abundance produced by one recruit and returns its
#' egg production, so the result is consistent with the second-order scheme
#' and with size diffusion. The recruitment flux needed to sustain that
#' abundance is read from the balance of the first size class. Ported from
#' R0_from_rates() in the cultivation-depensation repository.
#'
#' @return R0, with the steady-state abundance per recruit as attribute
#'   "per_recruit".
cod_R0 <- function(params, mortality_from = params, species = "Cod") {
    g <- getEGrowth(params)
    mort <- getMort(mortality_from)
    D <- getDiffusion(params)
    sp_names <- species_params(params)$species
    N0 <- setNames(numeric(length(sp_names)), sp_names)
    N0[species] <- 1
    n <- mizer:::get_steady_state_n(params, g, mort, D, N0,
                                    max_iterations = 20000)
    coefs <- mizer:::get_transport_coefs(
        params, n, g, mort, dt = 1,
        recruitment_flux = rep(0, length(sp_names)), d = D,
        flux_limiter = mizer:::flux_limiter_scheme(params))
    j <- params@w_min_idx[[species]]
    R_req <- dw(params)[j] * ((coefs$b[species, j] - 1) * n[species, j] +
                                  coefs$c[species, j] * n[species, j + 1])
    rdi <- mizerRDI(params, n = n, n_pp = initialNResource(params),
                    n_other = initialNOther(params), e_growth = g,
                    mort = mort, e_repro = getERepro(params), diffusion = D)
    structure(unname(rdi[species] / R_req),
              per_recruit = n[species, ] / R_req)
}

#' Fill in cod below 0.5 g, the recruitment size of the original model, with
#' the steady-state profile that matches the abundance at 0.5 g
fill_cod_larvae <- function(params) {
    N <- attr(cod_R0(params), "per_recruit")
    i_settle <- sum(w(params) < 0.5) + 1
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
    sp["Cod", "erepro"] <- sp["Cod", "erepro"] / cod_R0(params)
    given_species_params(params) <- sp
    params
}

#' Polish a state with Newton's method, restarting once if it stalls
#'
#' Check `attr(result, "convergence")$residual` before using the result.
newton_state <- function(params) {
    newton <- function(p) suppressWarnings(
        findSteadyState(p, solver = "newton", info_level = 0))
    out <- newton(params)
    if (attr(out, "convergence")$residual > 1e-5) out <- newton(out)
    out
}

#' The two states of the original model, carried over to the combined model
#'
#' The cod state of the original model (upwind scheme, cod from 0.5 g) is
#' interpolated onto the grid of the combined model at its baseline
#' (`w_shift` = 0.2 g, no larval predation), cod eggs and larvae are filled
#' in, and cod erepro is calibrated so that it is a steady state. The
#' cod-free state is interpolated and settled by projection.
#'
#' @param original_cod,original_no_cod The two states of the original model.
combined_anchor <- function(original_cod, original_no_cod, ...) {
    p <- new_combined_params(...)
    cod <- interpolate_state(p, original_cod) |>
        fill_cod_larvae() |>
        calibrate_cod_erepro() |>
        newton_state()
    # The first calibration is made in the copied-in state, which is not quite
    # a steady state of the new numerics. Repeat it at the solution.
    cod <- calibrate_cod_erepro(cod) |> newton_state()
    sim <- project(interpolate_state(cod, original_no_cod), t_max = 200,
                   t_save = 200, method = "tr_bdf2", progress_bar = FALSE)
    no_cod <- newton_state(finalParams(sim))
    list(cod = cod, no_cod = no_cod)
}

#' The cod state and the cod-free state of the combined model
#'
#' Starts from the anchor states, rescales cod erepro so that the anchor cod
#' state is a steady state of `params`, and finds the fixed point with
#' Newton's method. That converges even where the fixed point is unstable.
#' Projection is no substitute: at strong larval predation the anchor lies
#' outside the basin of a cod state that is locally stable, and cod die out.
#'
#' @param anchor The result of combined_anchor().
combined_states <- function(params, anchor) {
    cod <- copy_state(params, anchor$cod) |>
        calibrate_cod_erepro() |>
        newton_state()
    list(cod = cod, no_cod = copy_state(cod, anchor$no_cod))
}

#' Decompose the invasion fitness of cod in the cod-free state
#'
#' R0 of a rare cod in the cod-free state, relative to the cod state, where
#' R0 = 1. "food" takes only the growth and reproduction of cod from the
#' cod-free state, "larval_predation" only their mortality. "interaction" is
#' the factor by which R0 differs from the product of the two; it is below 1
#' when hungry juvenile cod linger in the sizes at which sprat eat them.
invasion_decomposition <- function(cod, no_cod) {
    r0_cod <- cod_R0(cod)
    out <- data.frame(
        R0 = cod_R0(no_cod) / r0_cod,
        food = cod_R0(no_cod, mortality_from = cod) / r0_cod,
        larval_predation = cod_R0(cod, mortality_from = no_cod) / r0_cod)
    out$interaction <- out$R0 / (out$food * out$larval_predation)
    out
}
