######################################################################
# 09_wind_conditional_ba.R   -- Tier B.2
# ---------------------------------------------------------------------
# Bland-Altman conditioned on wind direction sector (N / E / S / W).
# This addresses the reviewer's intuition that some lab pairs may
# agree well when wind is from one direction but disagree when wind
# is from the perpendicular direction (because background gradient
# changes character).
#
# We compute % bias and 95 % LoA for each (pair, variable, sector)
# triplet, then plot a 4-panel BA grid for the lab-internal pairs
# on the emission variables.
#
# Reviewer linkage: R1-L335.
######################################################################

stopifnot(exists("emission_reshaped"), exists("emission_result"))

# ---- attach wind sector to long-format emission data ---------------
v7_attach_sector <- function(df_long) {
        met <- emission_result %>%
                dplyr::select(DATE.TIME, wd_mst, ws_mst) %>%
                dplyr::distinct(DATE.TIME, .keep_all = TRUE) %>%
                dplyr::mutate(sector = v7_wd_sector(wd_mst))
        df_long %>% dplyr::left_join(met, by = "DATE.TIME")
}

# ---- BA stats for one (pair, var, sector) --------------------------
v7_ba_sector_pair <- function(df_long, var_name, pair, sector) {
        df <- v7_attach_sector(df_long) %>%
                dplyr::filter(.data$sector == !!sector,
                              .data$var == var_name,
                              .data$analyzer %in% pair) %>%
                dplyr::select(DATE.TIME, location, analyzer, value) %>%
                tidyr::pivot_wider(names_from = analyzer, values_from = value,
                                   values_fn = ~ mean(.x, na.rm = TRUE))
        if (!all(c(pair[1], pair[2]) %in% names(df))) return(NULL)
        df <- df %>%
                dplyr::rename(x = !!pair[1], y = !!pair[2]) %>%
                dplyr::filter(is.finite(x), is.finite(y), (x + y) != 0)
        if (nrow(df) < 5L) return(NULL)
        mn <- (df$x + df$y) / 2
        d  <- 100 * (df$y - df$x) / mn
        bias  <- mean(d); sd_d <- sd(d)
        tibble::tibble(
                variable    = var_name,
                analyzer_x  = pair[1],
                analyzer_y  = pair[2],
                sector      = sector,
                n           = nrow(df),
                bias_pct    = bias,
                loa_low     = bias - 1.96 * sd_d,
                loa_high    = bias + 1.96 * sd_d,
                sd_pct      = sd_d
        )
}

v7_ba_sector_var <- function(df_long, var_name, pairs = v7_pairs_lab,
                             sectors = c("N","E","S","W")) {
        out <- list()
        for (pn in names(pairs)) {
                for (sec in sectors) {
                        out[[paste(pn, sec, sep = "|")]] <-
                                v7_ba_sector_pair(df_long, var_name, pairs[[pn]], sec)
                }
        }
        dplyr::bind_rows(out)
}

v7_ba_sector_table <- dplyr::bind_rows(
        v7_ba_sector_var(emission_reshaped, "Q_vent"),
        v7_ba_sector_var(emission_reshaped, "e_CH4_ghLU"),
        v7_ba_sector_var(emission_reshaped, "e_NH3_ghLU")
) %>%
        dplyr::mutate(across(where(is.numeric), ~ v7_round(.x, 3)))

v7_write_csv(v7_ba_sector_table, "09_wind_conditional_BA.csv")

# ---- summary plot: bias by sector for each pair, faceted by variable
p_sector <- v7_ba_sector_table %>%
        dplyr::mutate(
                pair   = paste0(analyzer_x, " / ", analyzer_y),
                sector = factor(sector, levels = c("N","E","S","W"))
        ) %>%
        ggplot2::ggplot(ggplot2::aes(x = sector, y = bias_pct,
                                     fill = pair, ymin = loa_low, ymax = loa_high)) +
        ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
        ggplot2::geom_errorbar(position = ggplot2::position_dodge(width = 0.8),
                               width = 0.2, alpha = 0.8) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
        ggplot2::facet_wrap(~ variable, scales = "free_y", ncol = 1) +
        ggplot2::labs(x = "Wind direction sector", y = "Relative bias (%) [error bars = 95 % LoA]",
                      fill = NULL) +
        ggplot2::theme_classic(base_size = 10) +
        ggplot2::theme(legend.position = "bottom")
v7_save_plot(p_sector, "09_wind_conditional_BA_summary.png", width = 10, height = 9)

message("v7 module 09 (wind-conditional Bland-Altman): done.")
