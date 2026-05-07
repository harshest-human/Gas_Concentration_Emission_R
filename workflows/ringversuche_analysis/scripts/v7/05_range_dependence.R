######################################################################
# 05_range_dependence.R   -- Tier A.4
# ---------------------------------------------------------------------
# Heteroscedasticity check: does the magnitude of the difference grow
# with the magnitude of the measurement? If yes, fixed LoA in B-A
# under-estimate disagreement at high values.
#
# We regress |Diff_i| ~ Mean_i for every analyzer pair on every
# variable. A slope significantly > 0 indicates range-dependent error.
#
# We also report the relative dependence (slope expressed in % of
# the mean), so the manuscript can say e.g. "the residual
# inter-analyzer error grows by 3 % per unit of |delta_CO2|".
#
# Reviewer linkage: R1-L327, R2-L391.
######################################################################

stopifnot(exists("concentration_reshaped"), exists("emission_reshaped"))

v7_range_pair <- function(df_long, var_name, pair) {
        wide <- v7_to_wide(df_long, var_name)
        v <- v7_pair_vectors(wide, pair[1], pair[2])
        if (is.null(v)) return(NULL)
        m <- (v$x + v$y) / 2
        d <- abs(v$y - v$x)
        ok <- is.finite(m) & is.finite(d)
        if (sum(ok) < 5L) return(NULL)
        fit <- stats::lm(d[ok] ~ m[ok])
        s   <- summary(fit)
        # Spearman rank correlation as a non-parametric sanity check
        sp  <- suppressWarnings(stats::cor.test(m[ok], d[ok], method = "spearman"))
        tibble::tibble(
                variable    = var_name,
                analyzer_x  = pair[1],
                analyzer_y  = pair[2],
                n           = sum(ok),
                slope       = unname(coef(fit)[2]),
                slope_p     = s$coefficients[2, 4],
                R2          = s$r.squared,
                spearman_rho = unname(sp$estimate),
                spearman_p   = sp$p.value,
                # interpretation: relative slope = slope * 100 / mean(m)
                rel_slope_pct_per_unit_mean = 100 * unname(coef(fit)[2]) / mean(m[ok], na.rm = TRUE)
        )
}

v7_range_var <- function(df_long, var_name, pairs = v7_pairs_all) {
        purrr::map_dfr(pairs, ~ v7_range_pair(df_long, var_name, .x))
}

v7_range_table <- dplyr::bind_rows(
        v7_range_var(concentration_reshaped, "delta_CO2"),
        v7_range_var(concentration_reshaped, "delta_CH4"),
        v7_range_var(concentration_reshaped, "delta_NH3"),
        v7_range_var(emission_reshaped,      "Q_vent"),
        v7_range_var(emission_reshaped,      "e_CH4_ghLU"),
        v7_range_var(emission_reshaped,      "e_NH3_ghLU")
) %>%
        dplyr::mutate(
                hetero_sig = slope_p < 0.05,
                across(where(is.numeric), ~ v7_round(.x, 4))
        )

v7_write_csv(v7_range_table, "05_range_dependence.csv")
message("v7 module 05 (range dependence): done.")
