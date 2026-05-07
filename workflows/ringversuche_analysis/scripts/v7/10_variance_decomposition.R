######################################################################
# 10_variance_decomposition.R   -- Tier B.3
# ---------------------------------------------------------------------
# Random-effects variance decomposition: how much of the spread in
# the long-format reshaped values is between analyzers, between
# days, between hour-of-day, vs residual within (DATE.TIME analyzer)?
#
# This answers the reviewer's question (R1-L301): "is the
# inter-analyzer variability the dominant source of uncertainty,
# or is it day-to-day temporal variability?".
#
# We fit lme4::lmer with random intercepts:
#   value ~ 1 + (1 | analyzer) + (1 | day) + (1 | hour) + (1 | location)
# per variable, on the de-meaned long table, and report variance
# components and ICC for each random factor.
######################################################################

stopifnot(exists("concentration_reshaped"), exists("emission_reshaped"))

v7_var_decomp <- function(df_long, var_name) {
        df <- df_long %>%
                dplyr::filter(.data$var == var_name, .data$analyzer != "baseline") %>%
                dplyr::mutate(day  = factor(as.Date(DATE.TIME)),
                              hour = factor(format(DATE.TIME, "%H"))) %>%
                dplyr::filter(is.finite(value))
        if (nrow(df) < 50L) return(NULL)
        fit <- tryCatch(
                lme4::lmer(value ~ 1 + (1 | analyzer) + (1 | day) + (1 | hour) + (1 | location),
                           data = df, REML = TRUE,
                           control = lme4::lmerControl(check.conv.singular = "ignore")),
                error = function(e) NULL
        )
        if (is.null(fit)) return(NULL)
        vc <- as.data.frame(lme4::VarCorr(fit))
        # vc has columns: grp, var1, var2, vcov, sdcor
        total <- sum(vc$vcov, na.rm = TRUE)
        out <- vc %>%
                dplyr::transmute(
                        variable  = var_name,
                        component = grp,
                        variance  = vcov,
                        sd        = sdcor,
                        share_pct = 100 * vcov / total
                )
        out
}

v7_vc_table <- dplyr::bind_rows(
        v7_var_decomp(concentration_reshaped, "delta_CO2"),
        v7_var_decomp(concentration_reshaped, "delta_CH4"),
        v7_var_decomp(concentration_reshaped, "delta_NH3"),
        v7_var_decomp(emission_reshaped,      "Q_vent"),
        v7_var_decomp(emission_reshaped,      "e_CH4_ghLU"),
        v7_var_decomp(emission_reshaped,      "e_NH3_ghLU")
) %>%
        dplyr::mutate(across(where(is.numeric), ~ v7_round(.x, 4)))

v7_write_csv(v7_vc_table, "10_variance_decomposition.csv")

# ---- bar chart of share by variable -------------------------------
p_vc <- v7_vc_table %>%
        ggplot2::ggplot(ggplot2::aes(x = variable, y = share_pct, fill = component)) +
        ggplot2::geom_col(position = "stack", width = 0.7) +
        ggplot2::scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 100)) +
        ggplot2::labs(x = NULL, y = "Variance share", fill = NULL) +
        ggplot2::theme_classic(base_size = 11) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
                       legend.position = "bottom")
v7_save_plot(p_vc, "10_variance_decomposition.png", width = 10, height = 6)

message("v7 module 10 (variance decomposition): done.")
