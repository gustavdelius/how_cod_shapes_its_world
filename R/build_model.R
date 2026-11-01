# Build a two-species mizer model (cod and sprat) that reproduces the emergent
# Allee effect of Van Leeuwen, De Roos & Persson (2008) "How cod shapes its
# world", J. Sea Res. 60: 89-104.
#
# The model has the same set of feeding links as the stage-structured
# biomass model in the paper:
#   * sprat feed only on the plankton resource, with juveniles and adults
#     competing for it;
#   * cod feed only on sprat and on an unmodelled "alternative food", which on
#     its own covers only part of their maintenance needs;
#   * there is no cannibalism and no feeding of sprat on cod;
#   * there is no imposed density dependence in recruitment (R_max = Inf), so
#     all density dependence is emergent.
#
# Running this script writes two MizerParams objects with identical
# parameters but different initial states:
#   model/cod_sprat_with_cod.rds   the stable cod-sprat equilibrium
#   model/cod_sprat_no_cod.rds     the stable sprat-only equilibrium

library(mizer)

cod_sprat_species <- data.frame(
    species = c("Sprat", "Cod"),
    # Sizes (g). Cod enter the model as settled juveniles of 0.5 g; the
    # pelagic egg and larval phase is not modelled (as in the paper) and its
    # survival is folded into erepro.
    w_inf = c(15, 10000),
    w_mat = c(5, 500),
    w_min = c(0.001, 0.5),
    # Maximum intake rate h is derived by mizer from the age at maturity.
    age_mat = c(1.5, 4),
    # Feeding kernel: sprat eat zooplankton ~1000 times smaller than
    # themselves, cod eat prey ~30 times smaller.
    beta = c(1000, 30),
    sigma = c(2, 1.5),
    interaction_resource = c(1, 0),
    # Search volume of cod (sprat uses the mizer default)
    gamma = c(NA, 1.07e-11),
    # Alternative food for cod: an external encounter rate E_ext * w^n.
    # On its own it gives a feeding level of about 0.1, roughly half the
    # critical feeding level, so cod cannot persist without sprat.
    E_ext = c(0, 2.5),
    # Background mortality z0 + z_ext * w^(-1/3) (1/year)
    z0 = c(0.25, 0.8),
    z_ext = c(0, 0.8),
    # Reproductive efficiency; R_max is left at its default of Inf
    erepro = c(0.1, 0.1),
    # Length-weight relation w = a l^b (g, cm), used by the length-based
    # fishing selectivity
    a = c(0.01, 0.01),
    b = c(3, 3)
)

# Predator (rows) x prey (columns): only cod eat sprat
cod_sprat_interaction <- matrix(
    c(0, 0,
      1, 0),
    nrow = 2, byrow = TRUE,
    dimnames = list(predator = c("Sprat", "Cod"), prey = c("Sprat", "Cod")))

# Two fisheries, both switched off initially. With catchability 1 the effort
# equals the fishing mortality on fully selected fish (1/year).
#   * "Trawl" fishes cod, selecting by length around the length at maturity.
#   * "Sprat trawl" fishes sprat above 1 g; it is used only to explore the
#     paper's suggestion that thinning the prey can let cod recover.
cod_sprat_gears <- data.frame(
    gear = c("Trawl", "Sprat trawl"),
    species = c("Cod", "Sprat"),
    sel_func = c("sigmoid_length", "knife_edge"),
    l50 = c(35, NA),  # cm, roughly the length of a 500 g cod
    l25 = c(30, NA),
    knife_edge_size = c(NA, 1),  # g
    catchability = 1
)

#' Create the cod-sprat MizerParams object (initial state not yet settled)
new_cod_sprat_params <- function() {
    params <- newMultispeciesParams(
        cod_sprat_species,
        interaction = cod_sprat_interaction,
        gear_params = cod_sprat_gears,
        initial_effort = 0,
        no_w = 150,
        info_level = 0)
    setMetadata(
        params,
        title = "Cod-sprat emergent Allee effect",
        description = paste(
            "Two-species mizer model demonstrating the emergent Allee effect",
            "of Van Leeuwen, De Roos & Persson (2008), J. Sea Res. 60: 89-104."))
}

#' Settle a model onto the attractor reached from its current state
settle <- function(params, t_max = 500) {
    projectUntilSettled(params, t_max = t_max, progress_bar = FALSE,
                        info_level = 0) |>
        finalParams()
}

#' The cod-free steady state
cod_free_state <- function(params) {
    n <- initialN(params)
    n["Cod", ] <- 0
    initialN(params) <- n
    settle(params)
}

#' The cod-sprat steady state, reached from the cod-free state by adding cod
#' at high abundance
cod_state <- function(params, cod_biomass = 3e13) {
    params <- cod_free_state(params)
    n <- initialN(params)
    template <- initialN(new_cod_sprat_params())["Cod", ]
    n["Cod", ] <- template * cod_biomass /
        sum(template * w(params) * dw(params))
    initialN(params) <- n
    settle(params, t_max = 1000)
}

if (sys.nframe() == 0) {
    params <- new_cod_sprat_params()
    params_no_cod <- cod_free_state(params)
    params_cod <- cod_state(params)
    dir.create("model", showWarnings = FALSE)
    saveParams(params_no_cod, "model/cod_sprat_no_cod.rds")
    saveParams(params_cod, "model/cod_sprat_with_cod.rds")
    print(getBiomass(params_no_cod))
    print(getBiomass(params_cod))
}
