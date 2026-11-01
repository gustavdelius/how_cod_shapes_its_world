# Prototype: the emergent Allee effect and cultivation-depensation in one
# mizer model (see R/combined_model.R and combined/README.md).
#
# Run from the repository root:
#   Rscript combined/run_prototype.R
# It takes about 5 minutes on 8 cores and writes combined/results.rds and the
# figures in combined/.

library(mizer)
library(ggplot2)
library(patchwork)
library(parallel)
options(mizer_info_level = 0)
source("R/build_model.R")
source("R/analysis_functions.R")
source("R/combined_model.R")

original_cod <- readParams("model/cod_sprat_with_cod.rds")
original_no_cod <- readParams("model/cod_sprat_no_cod.rds")
cores <- min(8, detectCores())

states_at <- function(w_shift, theta) {
    combined_states(new_combined_params(theta = theta, w_shift = w_shift),
                    original_cod, original_no_cod)
}

# Regime map ---------------------------------------------------------------

evaluate <- function(w_shift, theta) {
    s <- states_at(w_shift, theta)
    cbind(data.frame(
        w_shift = w_shift, theta = theta,
        erepro = species_params(s$cod)["Cod", "erepro"],
        cod_biomass = getBiomass(s$cod)[["Cod"]],
        residual = attr(s$cod, "convergence")$residual,
        max_real_part = getStability(s$cod)$max_real_part),
        invasion_decomposition(s$cod, s$no_cod))
}
# Near theta = 0.04 the cod state reaches a fold (a real eigenvalue goes to
# zero), so the scan stops at 0.035
grid <- expand.grid(w_shift = seq(0.2, 0.8, by = 0.05),
                    theta = seq(0, 0.035, by = 0.005))
scan <- do.call(rbind, mclapply(seq_len(nrow(grid)), function(i)
    evaluate(grid$w_shift[i], grid$theta[i]), mc.cores = cores))
# Where Newton's method did not converge there is no cod state to report
scan$converged <- scan$residual < 1e-4
scan[!scan$converged, c("max_real_part", "R0", "food", "larval_predation",
                        "interaction")] <- NA

# Introductions of cod into the cod-free state ------------------------------

intro_points <- data.frame(w_shift = c(0.2, 0.45, 0.6, 0.6, 0.6),
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
    sim <- project(p, t_max = 300, t_save = 1, progress_bar = FALSE)
    data.frame(w_shift = pt$w_shift, theta = pt$theta, factor = fac,
               time = getTimes(sim),
               cod = getBiomass(sim)[, "Cod"] / getBiomass(s$cod)[["Cod"]])
}, mc.cores = cores))

saveRDS(list(scan = scan, intro = intro,
             mizer_version = as.character(packageVersion("mizer"))),
        "combined/results.rds")

# Figures ------------------------------------------------------------------

theme_set(theme_minimal(base_size = 12) +
              theme(panel.grid.minor = element_blank(),
                    legend.position = "top",
                    strip.text = element_text(face = "bold", hjust = 0)))

# R0 is 0 where the food barrier is absolute; cap the log scale there
map <- scan
map$log_R0 <- pmin(pmax(log10(map$R0), -2), 0.5)
map$log_food <- pmin(pmax(log10(map$food), -2), 0.5)
p_map <- ggplot(map, aes(w_shift, theta)) +
    geom_raster(aes(fill = log_R0)) +
    geom_contour(aes(z = log_R0), breaks = 0, colour = "black",
                 linewidth = 0.9) +
    geom_contour(aes(z = log_food), breaks = 0, colour = "black",
                 linewidth = 0.6, linetype = "dashed") +
    geom_point(data = map[map$converged & map$max_real_part > 0, ],
               shape = 4, size = 2) +
    annotate("label", x = 0.29, y = 0.003, size = 3.3,
             label = "Food barrier\n(Van Leeuwen et al.)") +
    annotate("label", x = 0.38, y = 0.024, size = 3.3,
             label = "Both barriers") +
    annotate("label", x = 0.68, y = 0.027, size = 3.3,
             label = "Larval-predation barrier\n(cultivation–depensation)") +
    annotate("label", x = 0.71, y = 0.0025, size = 3.3,
             label = "Cod invade:\nno bistability") +
    scale_fill_gradient2(low = "#2a78d6", mid = "white", high = "#eb6834",
                         midpoint = 0, limits = c(-2, 0.5), na.value = "grey80",
                         breaks = c(-2, -1, 0, 0.5),
                         labels = c("≤ 0.01", "0.1", "1", "3"),
                         name = expression(R[0] ~ "of cod in the cod-free state")) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = "Size at which cod larval food runs out, w_shift (g)",
         y = "Sprat predation on cod larvae, θ") +
    theme(legend.key.width = unit(1.5, "cm"))
ggsave("combined/regime_map.png", p_map, width = 8, height = 5.5, dpi = 110,
       bg = "white")

dec <- scan[scan$theta %in% c(0.01, 0.02) & scan$converged, ]
dec_long <- rbind(
    data.frame(dec[c("w_shift", "theta")], factor = "Food (growth only)",
               value = dec$food),
    data.frame(dec[c("w_shift", "theta")],
               factor = "Larval predation (mortality only)",
               value = dec$larval_predation),
    data.frame(dec[c("w_shift", "theta")], factor = "Both (R0)",
               value = dec$R0))
dec_long$factor <- factor(dec_long$factor, levels = unique(dec_long$factor))
dec_long$theta <- paste("θ =", dec_long$theta)
p_dec <- ggplot(dec_long[dec_long$value > 0, ],
                aes(w_shift, value, colour = factor)) +
    geom_hline(yintercept = 1, colour = "grey45", linetype = "dashed") +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.4) +
    facet_wrap(~theta) +
    scale_y_log10() +
    scale_colour_manual(values = c("#eb6834", "#2a78d6", "grey20"),
                        name = NULL) +
    labs(x = "Size at which cod larval food runs out, w_shift (g)",
         y = "R0 of a rare cod in the cod-free state")
ggsave("combined/decomposition.png", p_dec, width = 8, height = 4, dpi = 110,
       bg = "white")

intro$point <- factor(sprintf("w_shift = %.2f g, θ = %.2f", intro$w_shift,
                              intro$theta),
                      levels = sprintf("w_shift = %.2f g, θ = %.2f",
                                       intro_points$w_shift,
                                       intro_points$theta))
intro$introduced <- factor(paste0(100 * intro$factor, " %"),
                           levels = paste0(100 * intro_factors, " %"))
p_intro <- ggplot(intro, aes(time, cod, colour = introduced)) +
    geom_hline(yintercept = 1, colour = "grey45", linetype = "dashed") +
    geom_line(linewidth = 0.9) +
    facet_wrap(~point, ncol = 3) +
    scale_y_log10() +
    coord_cartesian(ylim = c(1e-6, 3)) +
    scale_colour_manual(values = c("#86b6ef", "#2a78d6", "#0d366b"),
                        name = "Cod introduced (fraction of cod state)") +
    labs(x = "Years after introduction",
         y = "Cod biomass / cod-state biomass")
ggsave("combined/introductions.png", p_intro, width = 9, height = 5.5,
       dpi = 110, bg = "white")

# Summary ------------------------------------------------------------------

print(format(scan[scan$w_shift %in% c(0.2, 0.3, 0.45, 0.6, 0.8) &
                      scan$theta %in% c(0, 0.01, 0.02, 0.04), ],
             digits = 3), row.names = FALSE)
