######################################################################
# 11_cross_correlation.R   -- Tier B.4
# ---------------------------------------------------------------------
# Cross-correlation function (CCF) between analyzer time series.
# Identifies whether one analyzer leads or lags another.
#
# For every analyzer pair on a given variable / location, we:
#   1. Pivot to a regular 7.5-min grid (the cleaning step's native
#      cadence).
#   2. Run stats::ccf() up to +- 8 lags (so +- 60 minutes).
#   3. Report lag of maximum |correlation|, that correlation, and
#      the lag-0 correlation.
#
# Reviewer linkage: R2-L353 ("are the analyzers in phase?").
######################################################################

stopifnot(exists("concentration_reshaped"), exists("emission_reshaped"))

v7_ccf_pair <- function(df_long, var_name, pair, location_filter, max_lag = 8L) {
        df <- df_long %>%
                dplyr::filter(.data$var == var_name,
                              .data$location == location_filter,
                              .data$analyzer %in% pair) %>%
                dplyr::group_by(DATE.TIME, analyzer) %>%
                dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
                tidyr::pivot_wider(names_from = analyzer, values_from = value)
        if (!all(c(pair[1], pair[2]) %in% names(df))) return(NULL)
        df <- df %>%
                dplyr::arrange(DATE.TIME) %>%
                dplyr::filter(is.finite(.data[[pair[1]]]),
                              is.finite(.data[[pair[2]]]))
        if (nrow(df) < 3L * max_lag) return(NULL)
        ccf_out <- tryCatch(
                stats::ccf(df[[pair[1]]], df[[pair[2]]],
                           lag.max = max_lag, plot = FALSE,
                           na.action = stats::na.exclude),
                error = function(e) NULL
        )
        if (is.null(ccf_out)) return(NULL)
        rs <- as.numeric(ccf_out$acf)
        ls <- as.numeric(ccf_out$lag)
        idx_max <- which.max(abs(rs))
        tibble::tibble(
                variable        = var_name,
                location        = location_filter,
                analyzer_x      = pair[1],
                analyzer_y      = pair[2],
                n               = nrow(df),
                lag_at_max      = ls[idx_max],
                r_at_max        = rs[idx_max],
                r_at_lag_0      = rs[ls == 0L],
                # Conventional sign reading: positive lag means y leads x
                interpretation  = dplyr::case_when(
                        ls[idx_max] >  0 ~ paste0(pair[2], " leads ", pair[1]),
                        ls[idx_max] <  0 ~ paste0(pair[1], " leads ", pair[2]),
                        TRUE             ~ "in phase"
                )
        )
}

v7_ccf_var <- function(df_long, var_name,
                       pairs = v7_pairs_lab,
                       locs  = c("North background","South background")) {
        out <- list()
        for (pn in names(pairs)) {
                for (loc in locs) {
                        out[[paste(pn, loc, sep = "|")]] <-
                                v7_ccf_pair(df_long, var_name, pairs[[pn]], loc)
                }
        }
        dplyr::bind_rows(out)
}

v7_ccf_table <- dplyr::bind_rows(
        v7_ccf_var(concentration_reshaped, "delta_CO2"),
        v7_ccf_var(concentration_reshaped, "delta_CH4"),
        v7_ccf_var(concentration_reshaped, "delta_NH3"),
        v7_ccf_var(emission_reshaped,      "Q_vent"),
        v7_ccf_var(emission_reshaped,      "e_CH4_ghLU"),
        v7_ccf_var(emission_reshaped,      "e_NH3_ghLU")
) %>%
        dplyr::mutate(across(where(is.numeric), ~ v7_round(.x, 3)))

v7_write_csv(v7_ccf_table, "11_cross_correlation.csv")
message("v7 module 11 (cross-correlation lag): done.")
