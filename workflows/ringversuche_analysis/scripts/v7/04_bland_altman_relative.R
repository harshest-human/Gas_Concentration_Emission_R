######################################################################
# 04_bland_altman_relative.R   -- Tier A.3
# ---------------------------------------------------------------------
# Bland-Altman analysis on the RELATIVE difference scale (% of pair
# mean) instead of the absolute scale, and reporting numerical
# bias and 95 % limits of agreement.
#
# Reviewer linkage: R2-L353, R1-L327
#   "Use relative differences in Bland-Altman; the absolute LoA grow
#    with the magnitude and obscure the agreement pattern".
#
# For each pair x_i, y_i the formulas are:
#   diff_pct_i = 100 * (y_i - x_i) / mean(x_i, y_i)
#   bias       = mean(diff_pct)
#   loa_low    = bias - 1.96 * sd(diff_pct)
#   loa_high   = bias + 1.96 * sd(diff_pct)
#   95 % CI on bias and on each LoA via the standard
#   normal-approximation formulas (Bland & Altman 1999).
######################################################################

stopifnot(exists("emission_reshaped"), exists("concentration_reshaped"))

v7_ba_relative <- function(x, y) {
        ok <- is.finite(x) & is.finite(y) & ((x + y) != 0)
        x <- x[ok]; y <- y[ok]
        n <- length(x)
        if (n < 5L) return(NULL)
        m <- (x + y) / 2
        d <- 100 * (y - x) / m
        bias  <- mean(d)
        sd_d  <- sd(d)
        loa_l <- bias - 1.96 * sd_d
        loa_h <- bias + 1.96 * sd_d
        # Bland & Altman (1999) standard errors
        se_bias <- sd_d / sqrt(n)
        se_loa  <- sd_d * sqrt(3 / n)
        tibble::tibble(
                n           = n,
                bias_pct    = bias,
                bias_lo95   = bias - 1.96 * se_bias,
                bias_hi95   = bias + 1.96 * se_bias,
                loa_low     = loa_l,
                loa_low_lo95= loa_l - 1.96 * se_loa,
                loa_low_hi95= loa_l + 1.96 * se_loa,
                loa_high    = loa_h,
                loa_high_lo95 = loa_h - 1.96 * se_loa,
                loa_high_hi95 = loa_h + 1.96 * se_loa,
                sd_pct      = sd_d
        )
}

v7_ba_pair <- function(df_long, var_name, pair, locs = c("North background","South background")) {
        wide <- v7_to_wide(df_long, var_name, locs)
        v <- v7_pair_vectors(wide, pair[1], pair[2])
        if (is.null(v)) return(NULL)
        out <- v7_ba_relative(v$x, v$y)
        if (is.null(out)) return(NULL)
        out %>% dplyr::mutate(variable = var_name, analyzer_x = pair[1], analyzer_y = pair[2],
                              .before = 1L)
}

v7_ba_var <- function(df_long, var_name, pairs = v7_pairs_all) {
        purrr::map_dfr(pairs, ~ v7_ba_pair(df_long, var_name, .x))
}

v7_ba_table <- dplyr::bind_rows(
        v7_ba_var(concentration_reshaped, "CO2_mgm3"),
        v7_ba_var(concentration_reshaped, "CH4_mgm3"),
        v7_ba_var(concentration_reshaped, "NH3_mgm3"),
        v7_ba_var(concentration_reshaped, "delta_CO2"),
        v7_ba_var(concentration_reshaped, "delta_CH4"),
        v7_ba_var(concentration_reshaped, "delta_NH3"),
        v7_ba_var(emission_reshaped,      "Q_vent"),
        v7_ba_var(emission_reshaped,      "e_CH4_ghLU"),
        v7_ba_var(emission_reshaped,      "e_NH3_ghLU")
) %>%
        dplyr::mutate(across(where(is.numeric), ~ v7_round(.x, 3)))

v7_write_csv(v7_ba_table, "04_bland_altman_relative.csv")

# ---- relative B-A plot for the three lab-internal pairs ------------
v7_ba_plot <- function(df_long, var_name, pair, location_filter = NULL) {
        wide <- v7_to_wide(df_long, var_name)
        if (!is.null(location_filter)) wide <- dplyr::filter(wide, location == location_filter)
        v <- v7_pair_vectors(wide, pair[1], pair[2])
        if (is.null(v)) return(ggplot2::ggplot() + ggplot2::theme_void())
        m <- (v$x + v$y) / 2
        d <- 100 * (v$y - v$x) / m
        ok <- is.finite(d)
        df <- tibble::tibble(mean_v = m[ok], diff_pct = d[ok])
        bias  <- mean(df$diff_pct);  sd_d <- sd(df$diff_pct)
        loa_l <- bias - 1.96 * sd_d; loa_h <- bias + 1.96 * sd_d
        ggplot2::ggplot(df, ggplot2::aes(x = mean_v, y = diff_pct)) +
                ggplot2::geom_point(alpha = 0.4, size = 0.7) +
                ggplot2::geom_hline(yintercept = bias,  colour = "blue", linetype = "dashed") +
                ggplot2::geom_hline(yintercept = loa_l, colour = "red",  linetype = "dotted") +
                ggplot2::geom_hline(yintercept = loa_h, colour = "red",  linetype = "dotted") +
                ggplot2::annotate("text", x = max(df$mean_v, na.rm = TRUE), y = bias,
                                  hjust = 1.05, vjust = -0.3, size = 3,
                                  label = sprintf("bias = %.1f%%", bias)) +
                ggplot2::labs(title = paste0(var_name, "  ", pair[1], " vs ", pair[2],
                                             ifelse(is.null(location_filter), "", paste0("  | ", location_filter))),
                              x = "Mean of pair", y = "Relative diff (%)") +
                ggplot2::theme_classic(base_size = 9)
}

# Render lab-internal ATB / LUFA / UB pairs for emission variables
for (pair_name in names(v7_pairs_lab)) {
        pair <- v7_pairs_lab[[pair_name]]
        plots <- list(
                v7_ba_plot(emission_reshaped, "Q_vent",     pair),
                v7_ba_plot(emission_reshaped, "e_CH4_ghLU", pair),
                v7_ba_plot(emission_reshaped, "e_NH3_ghLU", pair)
        )
        v7_save_plot(patchwork::wrap_plots(plots, ncol = 3),
                     paste0("04_BA_relative_lab_", pair_name, ".png"),
                     width = 14, height = 4)
}

message("v7 module 04 (relative Bland-Altman): done.")
