######################################################################
# 06_extended_pcc.R   -- Tier A.5
# ---------------------------------------------------------------------
# Extended Pearson correlation analysis. v6 already produces a
# correlogram (function `emicorrgram` lives in the v6 main script).
# Module 06 enriches it with:
#   - per-location correlation matrices (3 locations)
#   - per-time-resolution: raw 7.5-min, hourly mean, daily mean
#   - significance threshold + Bonferroni-corrected p
#   - tabular output for the response letter
#
# Reviewer linkage: R1-L301, R2-L391
#   "Demonstrate analyzer co-movement at multiple aggregation levels".
######################################################################

stopifnot(exists("concentration_reshaped"), exists("emission_reshaped"))

v7_pcc_at_resolution <- function(df_long, var_name, resolution = c("raw","hour","day")) {
        resolution <- match.arg(resolution)

        df <- df_long %>%
                dplyr::filter(.data$var == var_name, .data$analyzer != "baseline")

        if (resolution == "hour") {
                df <- df %>%
                        dplyr::mutate(grp = lubridate::floor_date(DATE.TIME, "hour")) %>%
                        dplyr::group_by(grp, location, analyzer) %>%
                        dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
                        dplyr::rename(DATE.TIME = grp)
        } else if (resolution == "day") {
                df <- df %>%
                        dplyr::mutate(grp = as.Date(DATE.TIME)) %>%
                        dplyr::group_by(grp, location, analyzer) %>%
                        dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
                        dplyr::rename(DATE.TIME = grp)
        }

        out <- list()
        for (loc in unique(df$location)) {
                wide <- df %>%
                        dplyr::filter(.data$location == loc) %>%
                        dplyr::select(DATE.TIME, analyzer, value) %>%
                        tidyr::pivot_wider(names_from = analyzer, values_from = value,
                                           values_fn = ~ mean(.x, na.rm = TRUE)) %>%
                        dplyr::select(-DATE.TIME)
                if (ncol(wide) < 2L || nrow(wide) < 4L) next
                rc <- Hmisc::rcorr(as.matrix(wide), type = "pearson")
                cm <- rc$r; pm <- rc$P
                pairs_idx <- which(upper.tri(cm), arr.ind = TRUE)
                if (nrow(pairs_idx) == 0L) next
                tab <- tibble::tibble(
                        variable    = var_name,
                        location    = loc,
                        resolution  = resolution,
                        analyzer_1  = colnames(cm)[pairs_idx[, 1]],
                        analyzer_2  = colnames(cm)[pairs_idx[, 2]],
                        r           = cm[pairs_idx],
                        p           = pm[pairs_idx]
                )
                out[[paste(var_name, loc, resolution, sep = "|")]] <- tab
        }
        dplyr::bind_rows(out)
}

v7_pcc_var <- function(df_long, var_name) {
        dplyr::bind_rows(
                v7_pcc_at_resolution(df_long, var_name, "raw"),
                v7_pcc_at_resolution(df_long, var_name, "hour"),
                v7_pcc_at_resolution(df_long, var_name, "day")
        )
}

v7_pcc_table <- dplyr::bind_rows(
        v7_pcc_var(concentration_reshaped, "CO2_mgm3"),
        v7_pcc_var(concentration_reshaped, "CH4_mgm3"),
        v7_pcc_var(concentration_reshaped, "NH3_mgm3"),
        v7_pcc_var(concentration_reshaped, "delta_CO2"),
        v7_pcc_var(concentration_reshaped, "delta_CH4"),
        v7_pcc_var(concentration_reshaped, "delta_NH3"),
        v7_pcc_var(emission_reshaped,      "Q_vent"),
        v7_pcc_var(emission_reshaped,      "e_CH4_ghLU"),
        v7_pcc_var(emission_reshaped,      "e_NH3_ghLU")
) %>%
        dplyr::mutate(
                p_bonf  = pmin(1, p * dplyr::n()),
                sig_05  = p < 0.05,
                sig_bonf = p_bonf < 0.05,
                across(where(is.numeric), ~ v7_round(.x, 4))
        )

v7_write_csv(v7_pcc_table, "06_extended_pcc.csv")
message("v7 module 06 (extended PCC): done.")
