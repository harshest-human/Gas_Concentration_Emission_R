######################################################################
# 07_ftir2_excluded.R   -- Tier A.6
# ---------------------------------------------------------------------
# The editor and Reviewer 2 noted that FTIR.2 (LUFA) had a known
# NH3 leak / drift during part of the campaign. v7 quantifies how
# much the conclusions move when FTIR.2's NH3 column (and FTIR.2 in
# general for sensitivity tests) is removed.
#
# This module re-computes the core agreement statistics from
# modules 02-04 with FTIR.2 excluded, side by side with the full
# dataset, so the response letter can present a single table of
# "with FTIR.2" vs "without FTIR.2" deltas.
#
# Reviewer linkage: Editor's framing letter, R2-L391.
######################################################################

stopifnot(exists("concentration_reshaped"), exists("emission_reshaped"))

# ---- subset helper --------------------------------------------------
v7_drop_ftir2_for_var <- function(df_long, var_name) {
        if (var_name == "NH3_mgm3" || var_name == "delta_NH3" || var_name == "e_NH3_ghLU") {
                df_long %>% dplyr::filter(!(analyzer == "FTIR.2"))
        } else {
                df_long
        }
}

# ---- regression and CCC and B-A on the reduced dataset --------------
v7_pairs_no_ftir2 <- Filter(function(p) !any(p == "FTIR.2"), v7_pairs_all)

# ---- regression sensitivity ---------------------------------------
v7_reg_no_ftir2 <- dplyr::bind_rows(
        v7_regression_var(concentration_reshaped, "delta_CO2", v7_pairs_no_ftir2),
        v7_regression_var(concentration_reshaped, "delta_CH4", v7_pairs_no_ftir2),
        v7_regression_var(v7_drop_ftir2_for_var(concentration_reshaped, "delta_NH3"),
                          "delta_NH3", v7_pairs_no_ftir2),
        v7_regression_var(emission_reshaped,      "Q_vent",     v7_pairs_no_ftir2),
        v7_regression_var(emission_reshaped,      "e_CH4_ghLU", v7_pairs_no_ftir2),
        v7_regression_var(v7_drop_ftir2_for_var(emission_reshaped, "e_NH3_ghLU"),
                          "e_NH3_ghLU", v7_pairs_no_ftir2)
) %>% dplyr::mutate(scenario = "FTIR.2 excluded", .before = 1L) %>%
        dplyr::mutate(across(where(is.numeric), ~ v7_round(.x, 4)))

v7_reg_full <- dplyr::bind_rows(
        v7_regression_var(concentration_reshaped, "delta_CO2"),
        v7_regression_var(concentration_reshaped, "delta_CH4"),
        v7_regression_var(concentration_reshaped, "delta_NH3"),
        v7_regression_var(emission_reshaped,      "Q_vent"),
        v7_regression_var(emission_reshaped,      "e_CH4_ghLU"),
        v7_regression_var(emission_reshaped,      "e_NH3_ghLU")
) %>% dplyr::mutate(scenario = "full dataset", .before = 1L) %>%
        dplyr::mutate(across(where(is.numeric), ~ v7_round(.x, 4)))

v7_sens_regression <- dplyr::bind_rows(v7_reg_full, v7_reg_no_ftir2)
v7_write_csv(v7_sens_regression, "07_ftir2_sensitivity_regression.csv")

# ---- CCC sensitivity -----------------------------------------------
v7_sens_ccc <- dplyr::bind_rows(
        # full
        purrr::map_dfr(c("delta_CH4","delta_NH3","Q_vent","e_CH4_ghLU","e_NH3_ghLU"), function(vn) {
                df_use <- if (vn %in% c("Q_vent","e_CH4_ghLU","e_NH3_ghLU")) emission_reshaped else concentration_reshaped
                v7_ccc_var(df_use, vn) %>% dplyr::mutate(scenario = "full dataset", .before = 1L)
        }),
        # excluded
        purrr::map_dfr(c("delta_CH4","delta_NH3","Q_vent","e_CH4_ghLU","e_NH3_ghLU"), function(vn) {
                df_base <- if (vn %in% c("Q_vent","e_CH4_ghLU","e_NH3_ghLU")) emission_reshaped else concentration_reshaped
                df_use  <- v7_drop_ftir2_for_var(df_base, vn) %>%
                        dplyr::filter(analyzer != "FTIR.2")
                v7_ccc_var(df_use, vn, pairs = v7_pairs_no_ftir2) %>%
                        dplyr::mutate(scenario = "FTIR.2 excluded", .before = 1L)
        })
) %>% dplyr::mutate(across(where(is.numeric), ~ v7_round(.x, 4)))

v7_write_csv(v7_sens_ccc, "07_ftir2_sensitivity_ccc.csv")

# ---- compact summary "delta" table -------------------------------
v7_sens_summary <- v7_sens_ccc %>%
        dplyr::filter(grepl("^FTIR\\.|CRDS", analyzer_x), grepl("^FTIR\\.|CRDS", analyzer_y)) %>%
        dplyr::group_by(scenario, variable) %>%
        dplyr::summarise(
                n_pairs       = dplyr::n(),
                mean_rho_c    = mean(rho_c, na.rm = TRUE),
                median_rho_c  = stats::median(rho_c, na.rm = TRUE),
                .groups = "drop"
        ) %>%
        tidyr::pivot_wider(names_from = scenario, values_from = c(n_pairs, mean_rho_c, median_rho_c)) %>%
        dplyr::mutate(
                d_mean_rho_c   = .data[["mean_rho_c_FTIR.2 excluded"]] - .data[["mean_rho_c_full dataset"]],
                d_median_rho_c = .data[["median_rho_c_FTIR.2 excluded"]] - .data[["median_rho_c_full dataset"]]
        ) %>%
        dplyr::mutate(across(where(is.numeric), ~ v7_round(.x, 4)))

v7_write_csv(v7_sens_summary, "07_ftir2_sensitivity_summary.csv")
message("v7 module 07 (FTIR.2 excluded sensitivity): done.")
