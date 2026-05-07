######################################################################
# 08_meteo_regression.R   -- Tier B.1
# ---------------------------------------------------------------------
# Multivariate regression to characterise how meteorological drivers
# influence the inter-analyzer disagreement.
#
# We model:
#   |Diff_pair| ~ wd_mst + ws_mst + temp_in + RH_in
# where |Diff_pair| is the absolute pairwise difference between two
# analyzers for delta-concentration / Q_vent / emissions.
#
# Wind direction is encoded as sin(wd) and cos(wd) so the linear
# model can pick up directional sensitivity without forcing a
# break at 0/360.
#
# Reviewer linkage: R1-L335.
######################################################################

stopifnot(exists("emission_reshaped"), exists("concentration_reshaped"),
          exists("emission_result"))   # emission_result has wide met cols

# ---- assemble a wide long-format table with met joined ------------
v7_join_meteo <- function(df_long) {
        # emission_result has DATE.TIME, met columns, and analyzer-level wide cols.
        # We extract just DATE.TIME, wd_mst, ws_mst, temp_in, RH_in (one row per timestamp).
        met <- emission_result %>%
                dplyr::select(DATE.TIME, dplyr::any_of(c("wd_mst","ws_mst","temp_in","RH_in"))) %>%
                dplyr::distinct(DATE.TIME, .keep_all = TRUE)
        df_long %>% dplyr::left_join(met, by = "DATE.TIME")
}

v7_meteo_pair <- function(df_long, var_name, pair) {
        df <- v7_join_meteo(df_long)
        wide <- df %>%
                dplyr::filter(.data$var == var_name, .data$analyzer %in% pair) %>%
                dplyr::select(DATE.TIME, location, analyzer, value, wd_mst, ws_mst, temp_in, RH_in) %>%
                tidyr::pivot_wider(names_from = analyzer, values_from = value,
                                   values_fn = ~ mean(.x, na.rm = TRUE))
        if (!all(c(pair[1], pair[2]) %in% names(wide))) return(NULL)
        wide <- wide %>%
                dplyr::mutate(diff_abs = abs(.data[[pair[2]]] - .data[[pair[1]]]),
                              wd_sin   = sin(wd_mst * pi / 180),
                              wd_cos   = cos(wd_mst * pi / 180)) %>%
                dplyr::filter(complete.cases(.))
        if (nrow(wide) < 12L) return(NULL)
        fit <- stats::lm(diff_abs ~ wd_sin + wd_cos + ws_mst + temp_in + RH_in, data = wide)
        s   <- summary(fit)
        broom::tidy(fit) %>%
                dplyr::mutate(
                        variable     = var_name,
                        analyzer_x   = pair[1],
                        analyzer_y   = pair[2],
                        n            = nrow(wide),
                        R2_overall   = s$r.squared,
                        adj_R2       = s$adj.r.squared,
                        F_p          = stats::pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3],
                                                 lower.tail = FALSE)
                ) %>%
                dplyr::select(variable, analyzer_x, analyzer_y, n, term, estimate, std.error,
                              statistic, p.value, R2_overall, adj_R2, F_p)
}

v7_meteo_var <- function(df_long, var_name, pairs = v7_pairs_lab) {
        purrr::map_dfr(pairs, ~ v7_meteo_pair(df_long, var_name, .x))
}

v7_meteo_table <- dplyr::bind_rows(
        v7_meteo_var(concentration_reshaped, "delta_CO2"),
        v7_meteo_var(concentration_reshaped, "delta_CH4"),
        v7_meteo_var(concentration_reshaped, "delta_NH3"),
        v7_meteo_var(emission_reshaped,      "Q_vent"),
        v7_meteo_var(emission_reshaped,      "e_CH4_ghLU"),
        v7_meteo_var(emission_reshaped,      "e_NH3_ghLU")
) %>%
        dplyr::mutate(across(where(is.numeric), ~ v7_round(.x, 4)))

v7_write_csv(v7_meteo_table, "08_meteo_regression.csv")
message("v7 module 08 (multivariate meteo regression): done.")
