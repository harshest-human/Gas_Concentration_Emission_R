######################################################################
# 02_regression_pairs.R   -- Tier A.1
# ---------------------------------------------------------------------
# OLS (y ~ x) regression for every analyzer pair on every variable
# of interest. Reports slope, intercept, R^2, SEE (residual SD), p
# of the slope, and 95 % CI on the slope.
#
# Reviewer linkage:
#   R2-L249, R2-L353, R1-L301  -- "Bland-Altman alone is not enough,
#   please present pairwise linear regression as well".
#
# Inputs:  v6 long-format datasets `concentration_reshaped` and
#          `emission_reshaped` (must be in scope).
# Outputs: result_data/tables/Version_7/02_regression_pairs.csv
#          result_data/plots/Version_7/02_regression_<var>.png
######################################################################

stopifnot(exists("concentration_reshaped"), exists("emission_reshaped"))

v7_regression_pair <- function(df_long, var_name, pair) {
        wide <- v7_to_wide(df_long, var_name)
        v   <- v7_pair_vectors(wide, pair[1], pair[2])
        if (is.null(v)) return(NULL)
        fit <- stats::lm(v$y ~ v$x)
        s   <- summary(fit)
        ci  <- tryCatch(stats::confint(fit, level = 0.95), error = function(e) matrix(NA, 2, 2))
        tibble::tibble(
                variable    = var_name,
                analyzer_x  = pair[1],
                analyzer_y  = pair[2],
                n           = v$n,
                slope       = unname(coef(fit)[2]),
                slope_lo95  = ci[2, 1],
                slope_hi95  = ci[2, 2],
                intercept   = unname(coef(fit)[1]),
                R2          = s$r.squared,
                adj_R2      = s$adj.r.squared,
                SEE         = s$sigma,
                slope_p     = s$coefficients[2, 4],
                # Diagnostic: does the 95 % slope CI cover 1?
                slope_eq1   = (ci[2, 1] <= 1 & ci[2, 2] >= 1),
                # Diagnostic: does the intercept differ from 0?
                intercept_p = s$coefficients[1, 4]
        )
}

v7_regression_var <- function(df_long, var_name, pairs = v7_pairs_all) {
        purrr::map_dfr(pairs, ~ v7_regression_pair(df_long, var_name, .x))
}

# ---- run on key variables ------------------------------------------
v7_reg_concentration <- dplyr::bind_rows(
        v7_regression_var(concentration_reshaped, "CO2_mgm3"),
        v7_regression_var(concentration_reshaped, "CH4_mgm3"),
        v7_regression_var(concentration_reshaped, "NH3_mgm3"),
        v7_regression_var(concentration_reshaped, "delta_CO2"),
        v7_regression_var(concentration_reshaped, "delta_CH4"),
        v7_regression_var(concentration_reshaped, "delta_NH3")
)
v7_reg_emission <- dplyr::bind_rows(
        v7_regression_var(emission_reshaped, "Q_vent"),
        v7_regression_var(emission_reshaped, "e_CH4_ghLU"),
        v7_regression_var(emission_reshaped, "e_NH3_ghLU")
)
v7_regression_table <- dplyr::bind_rows(v7_reg_concentration, v7_reg_emission) %>%
        dplyr::mutate(across(where(is.numeric), ~ v7_round(.x, 4)))

v7_write_csv(v7_regression_table, "02_regression_pairs.csv")

# ---- diagnostic scatter plots --------------------------------------
v7_plot_regression <- function(df_long, var_name, pairs = v7_pairs_cross) {
        plots <- lapply(pairs, function(pr) {
                wide <- v7_to_wide(df_long, var_name)
                v <- v7_pair_vectors(wide, pr[1], pr[2])
                if (is.null(v)) return(ggplot2::ggplot() + ggplot2::theme_void())
                d <- tibble::tibble(x = v$x, y = v$y)
                ggplot2::ggplot(d, ggplot2::aes(x = x, y = y)) +
                        ggplot2::geom_point(alpha = 0.4, size = 0.7) +
                        ggplot2::geom_abline(slope = 1, intercept = 0,
                                             linetype = "dashed", colour = "grey40") +
                        ggplot2::geom_smooth(method = "lm", formula = y ~ x,
                                             se = TRUE, colour = "red", linewidth = 0.6) +
                        ggplot2::labs(x = pr[1], y = pr[2],
                                      title = paste0(var_name, "  ", pr[1], " vs ", pr[2])) +
                        ggplot2::theme_classic(base_size = 9)
        })
        patchwork::wrap_plots(plots, ncol = 3)
}

for (vn in c("delta_CO2","delta_CH4","delta_NH3","Q_vent","e_CH4_ghLU","e_NH3_ghLU")) {
        df_use <- if (vn %in% c("Q_vent","e_CH4_ghLU","e_NH3_ghLU")) emission_reshaped else concentration_reshaped
        p <- v7_plot_regression(df_use, vn, pairs = v7_pairs_cross)
        v7_save_plot(p, paste0("02_regression_", vn, ".png"), width = 11, height = 11)
}

message("v7 module 02 (regression pairs): done.")
