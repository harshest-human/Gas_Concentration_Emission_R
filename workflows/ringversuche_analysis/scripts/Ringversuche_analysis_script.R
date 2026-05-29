################################################################################
##  Ringversuche analysis — Version 10
##  Builds on V9 (scripts/Ringversuche_analysis_script.R) and addresses the
##  first-round reviewer/editor comments. Each section header lists the
##  comment IDs it answers (see First_review/.../reviewer_comments_shortened.docx).
##
##  Key V10 changes vs V9:
##    * MIN_DELTA_CO2 guard on Q_vent — drops divide-by-near-zero spikes.
##    * FTIR.2 excluded from the *absolute* NH3 baseline only (offset cancels
##      in deltas, Q and e, so it stays in there).
##    * Paired NE vs SW outdoor-line test, then conditioned on wind sector
##      relative to neighbour sources (W dairy, SW feed store, N/NW manure
##      tank + biogas plant) — used to drop the contamination-prone line.
##    * Downstream delta/Q/e analysis runs on the single retained outdoor.
##    * Pairwise Deming + Lin's CCC extended to absolute c and Δc
##      (was Q/e only in V9).
##    * Tukey HSD pairwise tables (§16): Tables 3, 4, 6 of the manuscript —
##      absolute c, Δc and derived quantities, with Holm correction per panel
##      and ns/*/**/*** significance code.
##    * Delta plots no longer floor at 0 (sign carries meaning); absolute
##      plots still floor.
##    * Sampling-cycle figure for §2.2 (M&M).
##    * Text reports written to result_data/text_reports/Version_10 so the
##      manuscript can pull numbers directly.
################################################################################

#### Libraries                                                                ####
library(tidyverse)
library(lubridate)
library(scales)
library(patchwork)
# DescTools::CCC used inline.

#### Shared labels and aesthetics  (unchanged from V9)                        ####
VAR_LABELS_UNITS <- c(
        "CO2_mgm3"    = "c[CO2]~'(mg '*m^-3*')'",
        "CH4_mgm3"    = "c[CH4]~'(mg '*m^-3*')'",
        "NH3_mgm3"    = "c[NH3]~'(mg '*m^-3*')'",
        "r_CH4/CO2"   = "c[CH4]/c[CO2]~'('*'%'*')'",
        "r_NH3/CO2"   = "c[NH3]/c[CO2]~'('*'%'*')'",
        "delta_CO2"   = "Delta*c[CO2]~'(mg '*m^-3*')'",
        "delta_CH4"   = "Delta*c[CH4]~'(mg '*m^-3*')'",
        "delta_NH3"   = "Delta*c[NH3]~'(mg '*m^-3*')'",
        "Q_vent"      = "Q~'('*m^3~h^-1~LU^-1*')'",
        "e_CH4_gh"    = "e[CH4]~'(g '*h^-1*')'",
        "e_NH3_gh"    = "e[NH3]~'(g '*h^-1*')'",
        "e_CH4_ghLU"  = "e[CH4]~'(g '*h^-1~LU^-1*')'",
        "e_NH3_ghLU"  = "e[NH3]~'(g '*h^-1~LU^-1*')'",
        "temp"        = "Temperature~(degree*C)",
        "RH"          = "Relative~Humidity~('%')",
        "ws_mst"      = "Wind~Speed~Mast~(m~s^-1)",
        "wd_mst"      = "Wind~Direction~Mast~(degree)",
        "n_dairycows" = "'Number of Cows'"
)
VAR_LABELS_PLAIN <- c(
        "CO2_mgm3"    = "c[CO2]",
        "CH4_mgm3"    = "c[CH4]",
        "NH3_mgm3"    = "c[NH3]",
        "delta_CO2"   = "Delta*c[CO2]",
        "delta_CH4"   = "Delta*c[CH4]",
        "delta_NH3"   = "Delta*c[NH3]",
        "Q_vent"      = "Q",
        "e_CH4_ghLU"  = "e[CH4]",
        "e_NH3_ghLU"  = "e[NH3]"
)
LOC_LABELS <- c("Indoor" = "Indoor",
                "Outdoor_NE" = "Outdoor^NE",
                "Outdoor_SW" = "Outdoor^SW")
loc_from_suffix <- c("in" = "Indoor", "N" = "Outdoor_NE", "S" = "Outdoor_SW")

ANALYZER_COLORS <- c(
        "FTIR.1"="#1b9e77","FTIR.2"="#d95f02","FTIR.3"="#7570b3",
        "FTIR.4"="#e7298a","FTIR.4_old"="#e41a1c",
        "CRDS.1"="#66a61e","CRDS.2"="#e6ab02","CRDS.3"="#a6761d",
        "baseline"="black"
)
ANALYZER_SHAPES <- c(
        "FTIR.1"=0,"FTIR.2"=1,"FTIR.3"=2,"FTIR.4"=5,"FTIR.4_old"=13,
        "CRDS.1"=15,"CRDS.2"=19,"CRDS.3"=17,"baseline"=4
)
analyzer_aes <- function(x, lookup, default) {
        present <- unique(as.character(x))
        out     <- setNames(rep(default, length(present)), present)
        known   <- intersect(present, names(lookup))
        out[known] <- lookup[known]
        out
}

#### V10 ANALYSIS CONFIG                                                       ####
# CO2-balance method needs a non-trivial indoor-outdoor CO2 gradient to be
# well-conditioned. Below this floor, Q = PCO2/dCO2 blows up (V9 produced
# Q ~ 6e7 m^3/h/LU spikes — see V9 q_e_trend_plot.png). 50 mg/m^3 ~= 25 ppm,
# i.e. several sigma above analyzer noise. R2.24 (Reviewer #2) requires this
# kind of plausibility guard.
MIN_DELTA_CO2 <- 50      # mg m-3
MAX_Q_VENT    <- 5000    # m^3 h-1 LU-1, physical sanity ceiling for NVDB
MIN_Q_VENT    <- 50      # m^3 h-1 LU-1, physical sanity floor

# Reviewer 2 (R2.7) — FTIR.2 had a multiplicative-style offset on absolute
# NH3 only. The offset cancels in differences (delta NH3 -> Q -> e_NH3), so
# we exclude FTIR.2 ONLY from the absolute-NH3 baseline. FTIR.4_old is
# excluded everywhere after Section 14 (R2.7 closure).
BASELINE_EXCLUDE_DEFAULT <- c("FTIR.4_old")
BASELINE_EXCLUDE_BY_VAR  <- list(
        "NH3_mgm3"  = c("FTIR.2", "FTIR.4_old")
)

# Neighbour-source geometry of the LVAT site (manuscript Fig. 1):
#   W       -> additional dairy barns (cow emissions: CO2, CH4, NH3)
#   SW      -> COVERED feed-storage facility (negligible gas-phase emission)
#   N, NW   -> open manure tank + biogas plant (NH3, CH4)
#   NE,E,SE,S -> open land (clean reference)
#
# V10b revision (2026-05-28). Empirical observation across the campaign:
# Outdoor_NE concentrations were systematically HIGHER than Outdoor_SW for
# all three gases (CO2 +1-3 %, CH4 +17-45 %, NH3 +1-100 %), including the
# easterly sectors that the v3 geometry classed as "NE clean". This is
# inconsistent with the v3 geometry, which had NE clean under E/NE wind
# and SW downwind of LVAT under the same winds — under that scheme SW
# should have been HIGHER, not lower.
#
# Reconciliation: the LVAT plume reaches the NE sampling line under all
# southerly and most easterly sectors (ridge venting + side-curtain
# dispersion, not just the simple lee side of the barn), and the SW line
# enjoys near-clean status under the southerly and westerly sectors that
# v3 (incorrectly) coded as contaminated. The covered feed-storage facility
# SW of the barn is not a gas-phase source. Revised assignment:
#   NE line downwind of LVAT under: S, SW, W, SE, E (5 sectors)
#   SW line downwind of LVAT under: N, NE, E         (3 sectors, unchanged)
#   NE neighbour-source under:      N, NW             (unchanged)
#   SW neighbour-source under:      (none — feed store covered, W dairies
#                                    are upwind of the open WSW-NE land)
#   NE truly clean under:           NE only
#   SW truly clean under:           S, SE, SW, W
# Under this geometry the §11 auto-decision selects Outdoor_SW (lower
# contam_time) for the present E-prevailing campaign, which also matches
# the empirical SW < NE concentration ordering.
NE_CONTAM_LVAT     <- c("S","SW","W","SE","E")
SW_CONTAM_LVAT     <- c("N","NE","E")
NE_CONTAM_NEIGHBOUR <- c("N","NW")   # manure tank / biogas reach NE line
SW_CONTAM_NEIGHBOUR <- c()           # covered feed store; no chemical source
# Truly clean (open-land upwind) sectors for each line:
NE_CLEAN_SECTORS <- c("NE")
SW_CLEAN_SECTORS <- c("S","SE","SW","W")

#### Data preparation                                                         ####
# Indirect CO2 balance — V10 adds MIN_DELTA_CO2 guard and drops negative
# delta_CO2 (outdoor > indoor is physically inconsistent with the balance).
indirect.CO2.balance <- function(df, min_delta_co2 = MIN_DELTA_CO2,
                                  q_lo = MIN_Q_VENT, q_hi = MAX_Q_VENT) {
        ppm_to_mgm3 <- function(ppm, molar_mass) {
                T_K <- 273.15; P <- 101325; R <- 8.314472
                (ppm * 1e-6) * molar_mass * 1e3 * P / (R * T_K)
        }
        df <- df %>%
                mutate(
                        hour      = as.numeric(format(DATE.TIME, "%H")),
                        a         = 0.22,
                        h_min     = 2.9,
                        phi       = 5.6 * m_weight^0.75 + 22 * Y1_milk_prod + 1.6e-5 * p_pregnancy_day^3,
                        t_factor  = 1 + 4e-5 * (20 - temp_in)^3,
                        phi_T_cor = phi * t_factor,
                        A_cor     = 1 - a * sin((2*pi/24) * (hour + 6 - h_min)),
                        hpu_T_A_cor_per_cow = phi_T_cor * A_cor,
                        PCO2      = (0.185 * hpu_T_A_cor_per_cow) * 1000,

                        CO2_mgm3_in = ppm_to_mgm3(CO2_ppm_in, 44.01),
                        CO2_mgm3_N  = ppm_to_mgm3(CO2_ppm_N,  44.01),
                        CO2_mgm3_S  = ppm_to_mgm3(CO2_ppm_S,  44.01),
                        NH3_mgm3_in = ppm_to_mgm3(NH3_ppm_in, 17.031),
                        NH3_mgm3_N  = ppm_to_mgm3(NH3_ppm_N,  17.031),
                        NH3_mgm3_S  = ppm_to_mgm3(NH3_ppm_S,  17.031),
                        CH4_mgm3_in = ppm_to_mgm3(CH4_ppm_in, 16.04),
                        CH4_mgm3_N  = ppm_to_mgm3(CH4_ppm_N,  16.04),
                        CH4_mgm3_S  = ppm_to_mgm3(CH4_ppm_S,  16.04),

                        delta_CO2_N = CO2_mgm3_in - CO2_mgm3_N,
                        delta_CO2_S = CO2_mgm3_in - CO2_mgm3_S,
                        delta_NH3_N = NH3_mgm3_in - NH3_mgm3_N,
                        delta_NH3_S = NH3_mgm3_in - NH3_mgm3_S,
                        delta_CH4_N = CH4_mgm3_in - CH4_mgm3_N,
                        delta_CH4_S = CH4_mgm3_in - CH4_mgm3_S,

                        # MIN_DELTA_CO2 guard: only compute Q where the CO2
                        # gradient is well above analyzer noise.
                        Q_vent_N = ifelse(delta_CO2_N >= min_delta_co2,
                                          PCO2 / delta_CO2_N, NA_real_),
                        Q_vent_S = ifelse(delta_CO2_S >= min_delta_co2,
                                          PCO2 / delta_CO2_S, NA_real_),

                        # V10b: the [q_lo, q_hi] physical-sanity bound was
                        # removed because it silently dropped rows without
                        # documenting the loss; it is replaced by an explicit
                        # per-analyser 1.5*IQR filter on Q_vent applied in
                        # §11.5 with a dropout-count CSV.

                        e_NH3_gh_N = (delta_NH3_N * Q_vent_N / 1000) * n_dairycows_in,
                        e_CH4_gh_N = (delta_CH4_N * Q_vent_N / 1000) * n_dairycows_in,
                        e_NH3_gh_S = (delta_NH3_S * Q_vent_S / 1000) * n_dairycows_in,
                        e_CH4_gh_S = (delta_CH4_S * Q_vent_S / 1000) * n_dairycows_in,

                        e_NH3_ghLU_N = (e_NH3_gh_N * 500) / (n_dairycows_in * m_weight),
                        e_CH4_ghLU_N = (e_CH4_gh_N * 500) / (n_dairycows_in * m_weight),
                        e_NH3_ghLU_S = (e_NH3_gh_S * 500) / (n_dairycows_in * m_weight),
                        e_CH4_ghLU_S = (e_CH4_gh_S * 500) / (n_dairycows_in * m_weight)
                )
        for (gas in c("CH4", "NH3")) {
                for (loc in c("in", "N", "S")) {
                        col_gas   <- paste0(gas, "_mgm3_", loc)
                        col_CO2   <- paste0("CO2_mgm3_", loc)
                        ratio_col <- paste0("r_", gas, "/CO2_", loc)
                        df <- df %>%
                                mutate(!!ratio_col := ifelse(.data[[col_CO2]] != 0,
                                                             .data[[col_gas]] / .data[[col_CO2]] * 100,
                                                             NA_real_))
                }
        }
        df
}

# reshaper — V10: baseline can be conditioned on variable, so we can drop
# FTIR.2 from the absolute-NH3 baseline only (R2.7).
reshaper <- function(df,
                     exclude_default = BASELINE_EXCLUDE_DEFAULT,
                     exclude_by_var  = BASELINE_EXCLUDE_BY_VAR) {
        meta_cols <- c("DATE.TIME", "analyzer")
        measure_cols <- df %>% select(-all_of(meta_cols)) %>%
                select(where(is.numeric)) %>% names()

        df_long <- df %>%
                pivot_longer(cols = all_of(measure_cols),
                             names_to = "var", values_to = "value") %>%
                mutate(
                        location = case_when(
                                str_detect(var, "_in$") ~ "Indoor",
                                str_detect(var, "_N$")  ~ "Outdoor_NE",
                                str_detect(var, "_S$")  ~ "Outdoor_SW",
                                TRUE ~ NA_character_
                        ),
                        var = str_remove(var, "_(in|N|S)$"),
                        DATE.TIME = as.POSIXct(DATE.TIME),
                        day  = factor(as.Date(DATE.TIME)),
                        hour = factor(format(DATE.TIME, "%H:%M"))
                ) %>%
                select(DATE.TIME, day, hour, location, analyzer, var, value) %>%
                arrange(DATE.TIME, var, analyzer, location) %>%
                mutate(analyzer = case_when(
                        var %in% c("temp", "RH")                           ~ "HOBO",
                        var %in% c("wd_mst", "ws_mst", "wd_trv", "ws_trv") ~ "USA",
                        var %in% c("n_dairycows")                          ~ "RGB",
                        TRUE                                               ~ analyzer
                ))

        # Per-variable exclusion list for the baseline mean
        excl_lookup <- function(v) {
                if (!is.null(exclude_by_var[[v]])) exclude_by_var[[v]] else exclude_default
        }

        baseline_df <- df_long %>%
                group_by(DATE.TIME, location, var) %>%
                summarise(
                        excl  = list(excl_lookup(first(var))),
                        value = mean(value[!analyzer %in% excl[[1]]], na.rm = TRUE),
                        day   = first(day),
                        hour  = first(hour),
                        .groups = "drop"
                ) %>%
                select(-excl) %>%
                mutate(analyzer = "baseline")

        df_long <- bind_rows(df_long, baseline_df) %>%
                mutate(analyzer = factor(analyzer,
                                         levels = c("FTIR.1","FTIR.2","FTIR.3","FTIR.4","FTIR.4_old",
                                                    "CRDS.1","CRDS.2","CRDS.3",
                                                    "HOBO","USA","RGB","baseline"))) %>%
                arrange(DATE.TIME, location, var)
        df_long
}

#### Statistics helpers                                                       ####
deg_to_compass8 <- function(deg) {
        labs <- c("N","NE","E","SE","S","SW","W","NW")
        factor(labs[floor(((deg %% 360) + 22.5) / 45) %% 8 + 1], levels = labs)
}
deming_fit <- function(x, y) {
        ok <- complete.cases(x, y); x <- x[ok]; y <- y[ok]
        if (length(x) < 3) return(c(intercept = NA_real_, slope = NA_real_))
        sxx <- var(x); syy <- var(y); sxy <- cov(x, y)
        slope <- (syy - sxx + sqrt((syy - sxx)^2 + 4 * sxy^2)) / (2 * sxy)
        c(intercept = mean(y) - slope * mean(x), slope = slope)
}
ba_relative_stats <- function(x, y) {
        ok <- is.finite(x) & is.finite(y) & ((x + y) != 0)
        x <- x[ok]; y <- y[ok]; n <- length(x)
        if (n < 5L) return(NULL)
        d <- 100 * (y - x) / ((x + y) / 2)
        bias <- mean(d); sd_d <- sd(d); se_b <- sd_d / sqrt(n)
        tibble(n = n, bias_pct = bias,
               bias_lo95 = bias - 1.96 * se_b, bias_hi95 = bias + 1.96 * se_b,
               loa_low = bias - 1.96 * sd_d, loa_high = bias + 1.96 * sd_d,
               sd_pct = sd_d)
}

#### Plotting functions  (V10: signed axes for deltas)                        ####
# delta vars carry sign; absolute concentrations and Q/e are non-negative.
.is_signed_var <- function(v) any(grepl("^delta_", v))

emitrendplot <- function(data, y = NULL, location_filter = NULL,
                         plot_err = FALSE, x = "DATE.TIME") {
        if (!is.null(location_filter)) data <- data %>% filter(location %in% location_filter)
        if (!is.null(y))               data <- data %>% filter(var %in% y)
        value_col <- if (plot_err && "pct_err" %in% names(data)) "pct_err" else "value"
        summary_data <- data %>%
                group_by(.data[[x]], analyzer, location, var) %>%
                summarise(mean_val = mean(.data[[value_col]], na.rm = TRUE),
                          sd_val   = sd(.data[[value_col]], na.rm = TRUE),
                          .groups  = "drop") %>%
                filter(!is.na(mean_val))
        if (nrow(summary_data) == 0) stop("No data for emitrendplot")
        labels <- if (plot_err) VAR_LABELS_PLAIN else VAR_LABELS_UNITS
        summary_data <- summary_data %>%
                mutate(facet_label = factor(labels[as.character(var)], levels = labels[y]))
        col_vals   <- analyzer_aes(summary_data$analyzer, ANALYZER_COLORS, "black")
        shape_vals <- analyzer_aes(summary_data$analyzer, ANALYZER_SHAPES, 16)
        signed <- .is_signed_var(y)
        p <- ggplot(summary_data,
                    aes(x = .data[[x]], y = mean_val,
                        color = analyzer, shape = analyzer, group = analyzer)) +
                geom_line(linewidth = 0.5, alpha = 0.6, na.rm = TRUE) +
                geom_point(size = 2, na.rm = TRUE) +
                geom_errorbar(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
                              width = 0.4, na.rm = TRUE) +
                scale_color_manual(values = col_vals) +
                scale_shape_manual(values = shape_vals) +
                # When only one location is shown (e.g. retained outdoor only),
                # drop the column facet entirely so the redundant strip header
                # disappears. Otherwise face by location as usual.
                (if (length(unique(summary_data$location)) <= 1)
                        facet_grid(facet_label ~ ., scales = "free_y", switch = "y",
                                   labeller = labeller(facet_label = label_parsed))
                 else
                        facet_grid(facet_label ~ location, scales = "free_y", switch = "y",
                                   labeller = labeller(
                                           facet_label = label_parsed,
                                           location    = as_labeller(LOC_LABELS, label_parsed)))) +
                scale_y_continuous(breaks = scales::pretty_breaks(n = 6),
                                   labels = scales::number_format(accuracy = 0.1)) +
                labs(x = NULL, y = NULL) +
                theme_classic() +
                theme(text = element_text(size = 14),
                      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
                      axis.text.y = element_text(hjust = 1, size = 12),
                      strip.text.y.left = element_text(size = 14, vjust = 0.5),
                      panel.border = element_rect(color = "black", fill = NA),
                      legend.position = "bottom", legend.title = element_blank(),
                      plot.title = element_text(hjust = 0.5)) +
                guides(color = guide_legend(nrow = 1))
        if (signed) p <- p + geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3)
        else        p <- p + coord_cartesian(ylim = c(0, NA))
        if (x == "DATE.TIME") {
                x_breaks <- seq(min(summary_data$DATE.TIME, na.rm = TRUE),
                                max(summary_data$DATE.TIME, na.rm = TRUE),
                                by = "6 hours")
                p <- p + scale_x_datetime(breaks = x_breaks, date_labels = "%Y-%m-%d %H:%M")
        }
        p
}

emiboxplot <- function(data, y = NULL, location_filter = NULL, plot_err = FALSE) {
        if (!is.null(location_filter)) data <- data %>% filter(location %in% location_filter)
        if (!is.null(y))               data <- data %>% filter(var %in% y)
        value_col <- if (plot_err) "pct_err" else "value"
        data <- data %>%
                group_by(location, analyzer, var) %>%
                mutate(Q1 = quantile(.data[[value_col]], 0.25, na.rm = TRUE),
                       Q3 = quantile(.data[[value_col]], 0.75, na.rm = TRUE),
                       IQR = Q3 - Q1,
                       is_outlier = (.data[[value_col]] < (Q1 - 1.5 * IQR)) |
                                    (.data[[value_col]] > (Q3 + 1.5 * IQR))) %>%
                filter(!is_outlier) %>% ungroup() %>%
                select(-Q1, -Q3, -IQR, -is_outlier)
        labels <- if (plot_err) VAR_LABELS_PLAIN else VAR_LABELS_UNITS
        data <- data %>%
                mutate(variable_label = factor(labels[as.character(var)], levels = labels[y]))
        col_vals <- analyzer_aes(data$analyzer, ANALYZER_COLORS, "black")
        signed <- .is_signed_var(y)
        p <- ggplot(data, aes(x = analyzer, y = .data[[value_col]], color = analyzer)) +
                geom_boxplot(outlier.shape = NA, fill = NA, linewidth = 0.6) +
                geom_jitter(width = 0.2, alpha = 0.3, size = 1.5) +
                scale_color_manual(values = col_vals) +
                # Hide the redundant column strip when only one location is shown.
                (if (length(unique(data$location)) <= 1)
                        facet_grid(variable_label ~ ., scales = "free_y", switch = "y",
                                   labeller = labeller(variable_label = label_parsed))
                 else
                        facet_grid(variable_label ~ location, scales = "free_y", switch = "y",
                                   labeller = labeller(
                                           variable_label = label_parsed,
                                           location       = as_labeller(LOC_LABELS, label_parsed)))) +
                scale_y_continuous(breaks = scales::pretty_breaks(n = 5),
                                   labels = scales::number_format(accuracy = 0.1)) +
                labs(x = NULL, y = NULL) +
                theme_classic() +
                theme(text = element_text(size = 12),
                      axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
                      strip.text.y.left = element_text(angle = 0, hjust = 1),
                      panel.border = element_rect(color = "black", fill = NA),
                      legend.position = "bottom", legend.title = element_blank(),
                      plot.title = element_text(hjust = 0.5)) +
                guides(color = guide_legend(nrow = 1))
        if (signed) p <- p + geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3)
        else        p <- p + coord_cartesian(ylim = c(0, NA))
        p
}

bland_altman_plot <- function(data, var_filter, analyzer_pair,
                              location_filter = NULL, x = "DATE.TIME") {
        var_label_expr <- parse(text = VAR_LABELS_UNITS[[var_filter]])[[1]]
        df <- data %>% filter(var == var_filter, analyzer %in% analyzer_pair)
        if (!is.null(location_filter)) df <- df %>% filter(location %in% location_filter)
        df_wide <- df %>%
                select(all_of(c(x, "location", "analyzer", "value"))) %>%
                pivot_wider(names_from = analyzer, values_from = value)
        a1 <- analyzer_pair[1]; a2 <- analyzer_pair[2]
        df_ba <- df_wide %>%
                mutate(mean_val = (.data[[a1]] + .data[[a2]]) / 2,
                       diff_pct = 100 * (.data[[a2]] - .data[[a1]]) / mean_val) %>%
                filter(is.finite(diff_pct))
        bias   <- mean(df_ba$diff_pct, na.rm = TRUE)
        sd_d   <- sd(df_ba$diff_pct, na.rm = TRUE)
        loa_hi <- bias + 1.96 * sd_d; loa_lo <- bias - 1.96 * sd_d
        # Proportional-bias trend (R2.8): regress diff vs mean
        slope_fit <- lm(diff_pct ~ mean_val, data = df_ba)
        slope <- coef(slope_fit)[2]; intercept_fit <- coef(slope_fit)[1]
        subtitle_expr <- if (!is.null(location_filter)) {
                loc_expr <- parse(text = LOC_LABELS[location_filter])[[1]]
                bquote(.(var_label_expr) ~ "|" ~ .(loc_expr))
        } else var_label_expr
        ggplot(df_ba, aes(x = mean_val, y = diff_pct)) +
                geom_point(alpha = 0.5, size = 1) +
                geom_hline(yintercept = 0, color = "grey60") +
                geom_hline(yintercept = bias,   color = "blue", linetype = "dashed", linewidth = 0.9) +
                geom_hline(yintercept = loa_hi, color = "red",  linetype = "dotted", linewidth = 0.9) +
                geom_hline(yintercept = loa_lo, color = "red",  linetype = "dotted", linewidth = 0.9) +
                geom_abline(slope = slope, intercept = intercept_fit,
                            color = "darkgreen", linetype = "longdash", linewidth = 0.7) +
                annotate("text", x = Inf, y = bias, hjust = 1.05, vjust = -0.4, size = 3,
                         label = sprintf("bias = %.1f%% | slope = %.2g %%/unit",
                                         bias, slope)) +
                scale_x_continuous(breaks = pretty_breaks(n = 6),
                                   labels = scales::label_number(accuracy = 0.1)) +
                labs(subtitle = subtitle_expr,
                     x = bquote("Mean of" ~ .(a1) ~ "and" ~ .(a2)),
                     y = bquote("Relative difference (" ~ .(a2) - .(a1) ~ ") %")) +
                theme_classic() +
                theme(plot.subtitle = element_text(hjust = 0.5),
                      plot.title    = element_text(hjust = 0.5))
}

#### 0.  Paths, output directories                                            ####
base_dir   <- "D:/Data_Analysis_R/Gas_Concentration_Emission_R/workflows/ringversuche_analysis"
data_dir   <- file.path(base_dir, "clean_data/Version_9/long_format")
clean_dir  <- file.path(base_dir, "clean_data/Version_9")
meta_dir   <- file.path(base_dir, "meta_data")
tables_dir <- file.path(base_dir, "result_data/tables/Version_10")
plots_dir  <- file.path(base_dir, "result_data/plots/Version_10")
report_dir <- file.path(base_dir, "result_data/text_reports/Version_10")
for (d in c(tables_dir, plots_dir, report_dir))
        dir.create(d, showWarnings = FALSE, recursive = TRUE)

start_time <- as.POSIXct("2025-04-08 12:00:00", tz = "UTC")
end_time   <- as.POSIXct("2025-04-14 12:00:00", tz = "UTC")
gases <- c("CO2", "CH4", "NH3")

# Helper for writing one analysis block of the text report.
write_section <- function(file, title, body) {
        cat(strrep("=", 78), "\n", title, "\n", strrep("=", 78), "\n",
            body, "\n\n", file = file, append = TRUE, sep = "")
}
report_file <- function(name) file.path(report_dir, name)
init_report <- function(name, header) {
        f <- report_file(name); cat(header, "\n\n", file = f, sep = "")
        invisible(f)
}

#### 1.  Sampling-cycle figure — clock-style, 60 min = two 30-min cycles      ####
# Two 30-min cycles laid out around a 60-min clock face. Each step is 7.5 min:
# 3.0 min flush (grey, labelled "flushing") + 4.5 min average (white, labelled
# with the sampling location, NE/SW rendered as superscripts via plotmath).
# Ticks placed at every step boundary (0, 3, 7.5, 10.5, 15, 18, 22.5, 25.5,
# 30, 33, 37.5, 40.5, 45, 48, 52.5, 55.5) — i.e. each flush end + each
# average end. Box ring is intentionally thin (radial width 0.3).

cycle_segments <- tibble(
        cycle    = rep(c(1, 2), each = 8),
        step     = rep(1:4, each = 2, times = 2),
        segment  = rep(c("flush", "average"), times = 8),
        location = rep(c("Indoor", "Outdoor_NE", "Indoor", "Outdoor_SW"),
                       each = 2, times = 2)
) %>%
        mutate(step_start = (cycle - 1) * 30 + (step - 1) * 7.5,
               start_min  = step_start + ifelse(segment == "flush", 0, 3),
               end_min    = step_start + ifelse(segment == "flush", 3, 7.5),
               midpoint   = (start_min + end_min) / 2,
               # plotmath expressions — parse=TRUE in geom_text renders them.
               # atop() stacks "measuring" above the location label inside
               # each averaged box; NE/SW become superscripts via ^"...".
               text_label = case_when(
                       segment  == "flush"      ~ "flushing",
                       location == "Indoor"     ~ "atop('measuring','Indoor')",
                       location == "Outdoor_NE" ~ "atop('measuring',Outdoor^'NE')",
                       location == "Outdoor_SW" ~ "atop('measuring',Outdoor^'SW')"
               ))

# Major tick positions = step boundaries (these get numeric "X minutes" labels).
tick_breaks <- c(0, 3, 7.5, 10.5, 15, 18, 22.5, 25.5,
                 30, 33, 37.5, 40.5, 45, 48, 52.5, 55.5)

# Ruler/protractor minor ticks every 0.5 min — short radial lines just outside
# the box ring, no numeric label.
major_ticks <- tibble(x = tick_breaks)
minor_ticks <- tibble(x = setdiff(seq(0, 59.5, by = 0.5), tick_breaks))

# Inner step ring — one box per 7.5-min step (4 steps per cycle x 2 cycles).
step_blocks <- tibble(
        cycle = rep(c(1, 2), each = 4),
        step  = rep(1:4, times = 2)
) %>%
        mutate(start_min = (cycle - 1) * 30 + (step - 1) * 7.5,
               end_min   = start_min + 7.5,
               midpoint  = (start_min + end_min) / 2,
               label     = paste0("Step ", step))

centre_label <- tibble(x = 0, y = 0, label = "1 h\n= 2 cycles")

cycle_plot <- ggplot(cycle_segments) +
        # Inner ring — one box per 7.5-min step labelled "Step 1..4".
        geom_rect(data = step_blocks,
                  aes(xmin = start_min, xmax = end_min,
                      ymin = 0.80, ymax = 0.97),
                  fill = "grey97", color = "grey55", linewidth = 0.25) +
        geom_text(data = step_blocks,
                  aes(x = midpoint, y = 0.885, label = label),
                  size = 2.8) +
        # Main ring — flush (grey) and average (white) sub-segments.
        geom_rect(data = filter(cycle_segments, segment == "flush"),
                  aes(xmin = start_min, xmax = end_min,
                      ymin = 1.0, ymax = 1.3),
                  fill = "grey85", color = "black", linewidth = 0.3) +
        geom_rect(data = filter(cycle_segments, segment == "average"),
                  aes(xmin = start_min, xmax = end_min,
                      ymin = 1.0, ymax = 1.3),
                  fill = "white", color = "black", linewidth = 0.3) +
        geom_text(aes(x = midpoint, y = 1.15, label = text_label),
                  parse = TRUE, size = 2.8, lineheight = 0.85) +
        # Ruler-style minor ticks (lines only — no labels).
        geom_segment(data = minor_ticks,
                     aes(x = x, xend = x, y = 1.30, yend = 1.34),
                     color = "grey55", linewidth = 0.25) +
        # Ruler-style major ticks (longer black lines, labels via scale).
        geom_segment(data = major_ticks,
                     aes(x = x, xend = x, y = 1.30, yend = 1.42),
                     color = "black", linewidth = 0.5) +
        geom_text(data = centre_label,
                  aes(x = x, y = y, label = label),
                  size = 4.0, fontface = "bold", lineheight = 0.9) +
        coord_polar(theta = "x", start = 0, direction = 1, clip = "off") +
        scale_x_continuous(limits = c(0, 60),
                           breaks = tick_breaks,
                           labels = function(x) paste0(x, " minutes")) +
        scale_y_continuous(limits = c(0, 1.6), expand = c(0, 0)) +
        labs(title    = "60-min sampling protocol: two 30-min cycles",
             subtitle = "Each step: 3.0 min flushing (discarded) + 4.5 min measuring (retained)") +
        theme_minimal(base_size = 16) +
        theme(axis.title         = element_blank(),
              axis.text.y        = element_blank(),
              axis.ticks.y       = element_blank(),
              axis.text.x        = element_text(size = 13, face = "bold"),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank(),
              plot.title         = element_text(hjust = 0.5, size = 17),
              plot.subtitle      = element_text(hjust = 0.5, size = 12),
              plot.margin        = margin(0.4, 2.5, 0.4, 2.5, unit = "cm"))
ggsave(file.path(plots_dir, "sampling_cycle.png"), cycle_plot,
       width = 10, height = 8.5, dpi = 300, bg = "white")

# ----- Hourly-averaging methodology paragraph (for manuscript §2.2) ----------
hourly_averaging_paragraph <- paste(
        "Each analyser produces eight averaged 4.5-min samples per hour, one per step",
        "of the four-step cycle repeated twice. The hourly mean at each sampling",
        "location is computed as the arithmetic mean of the corresponding samples:",
        "four samples per hour for the indoor ring line (Step 1 and Step 3 of both",
        "cycles) and two samples per hour for each outdoor line (Step 2 of both cycles",
        "for Outdoor_NE; Step 4 of both cycles for Outdoor_SW). Step-to-location",
        "mapping is implemented in three ways depending on each analyser's data",
        "delivery format. Line-based FTIR analysers (ATB FTIR.1: Line 1 = NE,",
        "2 = Indoor, 3 = SW; LUFA FTIR.2: Messstelle 1 = Indoor, 2 = SW, 3 = NE) use",
        "the explicit line identifier embedded in the raw results file. MPV-based CRDS",
        "analysers (ATB CRDS.1: positions 1/2/3 = NE/Indoor/SW; UB CRDS.3: positions",
        "8/1/9 = NE/Indoor/SW; LUFA CRDS.2: positions 3/1/2 = NE/Indoor/SW) treat each",
        "contiguous run of identical MPV position as one step. Time-cycle analysers",
        "(MBBM FTIR.3 and ANECO FTIR.4) have no embedded line identifier; the step is",
        "inferred from the clock time within the cycle, anchored at the campaign",
        "start (2025-04-08 12:00 UTC) with a fixed (Indoor, Outdoor_NE, Indoor,",
        "Outdoor_SW) rotation every 7.5 min. After step assignment, the first 3 min",
        "of each step is discarded (flush) and the remaining 4.5 min is averaged.",
        "The retained averages are then aggregated to hourly means by location for",
        "downstream Δc, Q and emission calculations.",
        sep = " ")
writeLines(hourly_averaging_paragraph,
           file.path(report_dir, "07_hourly_averaging_method.txt"))

#### 2.  Read & merge gas datasets                                            ####
gas_files <- list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)
gas_data <- map_dfr(gas_files, read.csv, stringsAsFactors = FALSE) %>%
        select(-any_of("MPVPosition")) %>%
        mutate(DATE.TIME = dmy_hm(DATE.TIME, tz = "UTC")) %>%
        filter(DATE.TIME >= start_time & DATE.TIME <= end_time) %>%
        pivot_wider(id_cols     = c(DATE.TIME, lab, analyzer),
                    names_from  = location,
                    values_from = c(CO2, CH4, NH3, H2O),
                    values_fn   = mean) %>%
        rename(CO2_ppm_in = CO2_in, CO2_ppm_N = CO2_N, CO2_ppm_S = CO2_S,
               CH4_ppm_in = CH4_in, CH4_ppm_N = CH4_N, CH4_ppm_S = CH4_S,
               NH3_ppm_in = NH3_in, NH3_ppm_N = NH3_N, NH3_ppm_S = NH3_S) %>%
        arrange(DATE.TIME, lab, analyzer)

#### 3.  Animal / climate / wind data                                         ####
animal_data <- read.csv(file.path(meta_dir, "RGB_Animal_count/20250408-15_LVAT_Animal_data.csv"),
                        stringsAsFactors = FALSE) %>%
        mutate(DATE.TIME = dmy_hm(DATE.TIME, tz = "UTC")) %>%
        filter(DATE.TIME >= start_time & DATE.TIME <= end_time) %>%
        select(-hour) %>%
        rename(n_dairycows_in = n_animals) %>% distinct()

T_RH_HOBO <- read.csv(file.path(meta_dir, "HOBO_Temp_RH/2025/20250408-20250630_HOBO_Temp_RH_1hour.csv"),
                      stringsAsFactors = FALSE) %>%
        mutate(DATE.TIME = mdy_hms(paste(date, time), tz = "UTC")) %>%
        filter(DATE.TIME >= start_time & DATE.TIME <= end_time) %>%
        rename(temp_N  = T_outside, RH_N = RH_out,
               temp_in = T_inside,  RH_in = RH_inside) %>%
        select(DATE.TIME, temp_N, RH_N, temp_in, RH_in)

wind_data <- read.csv(file.path(meta_dir, "USA_mast_wind/20240101_20250825_USA_mast_16_hourly_uvw_wd_ws.csv"),
                      stringsAsFactors = FALSE) %>%
        mutate(DATE.TIME = as.POSIXct(datetime_hour, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) %>%
        filter(DATE.TIME >= start_time & DATE.TIME <= end_time) %>%
        rename(wd_mst = wd, ws_mst = ws) %>%
        select(DATE.TIME, wd_mst, ws_mst)

#### 4.  Combined inputs + emissions                                          ####
input_combined <- gas_data %>%
        left_join(animal_data, by = "DATE.TIME") %>%
        left_join(T_RH_HOBO,   by = "DATE.TIME") %>%
        left_join(wind_data,   by = "DATE.TIME") %>%
        arrange(DATE.TIME, lab, analyzer)

emission_result <- indirect.CO2.balance(input_combined)
emission_reshaped <- reshaper(emission_result) %>%
        mutate(across(where(is.numeric), ~ round(.x, 2)))

write_excel_csv(input_combined,    file.path(tables_dir, "20250408-15_input_combined.csv"))
write_excel_csv(emission_result,   file.path(tables_dir, "20250408-15_emission_result.csv"))
write_excel_csv(emission_reshaped, file.path(tables_dir, "20250408-15_ringversuche_emission_reshaped.csv"))

#### 4.5 Raw 7.5-min descriptive stats aggregation  (handover §2.3 / Annex A1) ####
# Single consolidated CSV for the manuscript Annex (per-Annex-table-1
# question — option 2). Per-analyser raw_stats_<lab>_<analyzer>.csv files
# are produced by the cleaners (FTIR/CRDS_ringversuche_cleaning.R) on the
# pre-Stage-1 7.5-min data; if any analyser is missing its file (e.g. the
# cleaners have not been re-run since the new logic was added), the
# fallback computes the stats from the existing post-Stage-1 7.5_avg CSVs
# in clean_dir and flags those rows in a `source` column.

rs_files <- list.files(clean_dir,
                       pattern = "^raw_stats_.+\\.csv$",
                       full.names = TRUE)
raw_stats_pre <- if (length(rs_files) > 0) {
        map_dfr(rs_files, function(f) {
                d <- read.csv(f, stringsAsFactors = FALSE)
                d$source <- "pre-Stage-1 (from cleaner)"
                d
        })
} else {
        tibble()
}
known_pairs <- if (nrow(raw_stats_pre) > 0) {
        unique(paste(raw_stats_pre$lab, raw_stats_pre$analyzer))
} else {
        character(0)
}

cycle_files <- list.files(clean_dir,
                          pattern = "^20250408-15_7\\.5_avg_.+\\.csv$",
                          full.names = TRUE)
raw_stats_fallback <- map_dfr(cycle_files, function(f) {
        bn <- basename(f)
        m  <- regmatches(bn, regexec("^20250408-15_7\\.5_avg_([^_]+)_(.+)\\.csv$", bn))[[1]]
        if (length(m) < 3) return(NULL)
        lab_x <- m[2]; ana_x <- m[3]
        if (paste(lab_x, ana_x) %in% known_pairs) return(NULL)
        d <- read.csv(f, stringsAsFactors = FALSE)
        gas_cols <- intersect(c("CO2","CH4","NH3","H2O","N2O"), names(d))
        if (length(gas_cols) == 0 || !"location" %in% names(d)) return(NULL)
        d %>%
                pivot_longer(any_of(gas_cols), names_to = "gas", values_to = "value") %>%
                filter(!is.na(value), is.finite(value)) %>%
                group_by(location, gas) %>%
                summarise(n = n(), mean = mean(value),
                          median = median(value), min = min(value),
                          max = max(value), sd = sd(value), .groups = "drop") %>%
                mutate(lab = lab_x, analyzer = ana_x,
                       source = "post-Stage-1 (from 7.5_avg CSV)",
                       .before = 1)
})

raw_stats_all <- bind_rows(raw_stats_pre, raw_stats_fallback) %>%
        arrange(lab, analyzer, location, gas) %>%
        mutate(across(c(mean, median, min, max, sd), ~ round(.x, 3))) %>%
        relocate(lab, analyzer, location, gas, n,
                 mean, median, min, max, sd, source)
write_excel_csv(raw_stats_all,
                file.path(tables_dir, "raw_stats_all_analyzers.csv"))

#### 5.  Animal-count diagnostic  (R2.22 — were enough cows present?)         ####
animal_summary <- animal_data %>%
        summarise(n_hours       = n(),
                  n_mean        = mean(n_dairycows_in, na.rm = TRUE),
                  n_median      = median(n_dairycows_in, na.rm = TRUE),
                  n_min         = min(n_dairycows_in, na.rm = TRUE),
                  n_max         = max(n_dairycows_in, na.rm = TRUE),
                  n_hours_lt_30 = sum(n_dairycows_in < 30, na.rm = TRUE),
                  n_hours_lt_40 = sum(n_dairycows_in < 40, na.rm = TRUE))
write_excel_csv(animal_summary, file.path(tables_dir, "animal_count_summary.csv"))

animal_trend <- ggplot(animal_data, aes(x = DATE.TIME, y = n_dairycows_in)) +
        geom_line(color = "#377EB8") + geom_point(size = 0.7, color = "#377EB8") +
        geom_hline(yintercept = 40, linetype = "dashed", color = "grey50") +
        labs(title = "Hourly in-barn cow count, RGB-camera derived",
             subtitle = "Dashed line at 40 cows. Exercise yard closed throughout the campaign.",
             x = NULL, y = "Cows inside barn") +
        scale_x_datetime(date_breaks = "1 day", date_labels = "%m-%d") +
        theme_bw(base_size = 12)
ggsave(file.path(plots_dir, "animal_count_trend.png"), animal_trend,
       width = 10, height = 4, dpi = 300, bg = "white")

#### 6.  Absolute concentration plots WITH FTIR.4_old (justifies its drop)    ####
c_trend_plot <- emitrendplot(emission_reshaped, y = c("CO2_mgm3","CH4_mgm3","NH3_mgm3"))
c_boxplot    <- emiboxplot  (emission_reshaped, y = c("CO2_mgm3","CH4_mgm3","NH3_mgm3"))
ggsave(file.path(plots_dir, "c_trend_plot.png"), c_trend_plot, width = 12, height = 8, dpi = 300)
ggsave(file.path(plots_dir, "c_boxplot.png"),    c_boxplot,    width = 12, height = 8, dpi = 300)

#### 7.  FTIR.4 vs FTIR.4_old — spectral-library re-evaluation check (R2.7)   ####
ftir_old <- input_combined %>% filter(analyzer == "FTIR.4_old")
ftir_new <- input_combined %>% filter(analyzer == "FTIR.4")
lib_pair <- function(g, l) {
        col <- paste0(g, "_ppm_", l)
        inner_join(tibble(DATE.TIME = ftir_old$DATE.TIME, v_old = ftir_old[[col]]),
                   tibble(DATE.TIME = ftir_new$DATE.TIME, v_new = ftir_new[[col]]),
                   by = "DATE.TIME") %>%
                filter(complete.cases(v_old, v_new))
}
lib_compare <- map_dfr(gases, function(g) {
        map_dfr(c("in", "N", "S"), function(l) {
                j <- lib_pair(g, l)
                if (nrow(j) < 3) return(NULL)
                dem <- deming_fit(j$v_old, j$v_new)
                tibble(gas = g, location = unname(loc_from_suffix[l]),
                       n   = nrow(j),
                       mean_old = mean(j$v_old), mean_new = mean(j$v_new),
                       mean_diff = mean(j$v_new - j$v_old),
                       RPD_pct = 100 * mean(j$v_new - j$v_old) /
                                 mean((j$v_old + j$v_new) / 2),
                       RMSD = sqrt(mean((j$v_new - j$v_old)^2)),
                       deming_intercept = unname(dem["intercept"]),
                       deming_slope     = unname(dem["slope"]),
                       CCC = DescTools::CCC(j$v_old, j$v_new)$rho.c[, "est"],
                       p_wilcox = suppressWarnings(
                               wilcox.test(j$v_old, j$v_new, paired = TRUE)$p.value))
        })
})
write_excel_csv(lib_compare, file.path(tables_dir, "FTIR4_old_vs_new_comparison.csv"))

lib_points <- map_dfr(gases, function(g) {
        map_dfr(c("in", "N", "S"), function(l) {
                lib_pair(g, l) %>%
                        mutate(gas = g, location = unname(loc_from_suffix[l]))
        })
})
lib_plot <- ggplot(lib_points, aes(x = v_old, y = v_new)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
        geom_point(alpha = 0.5, size = 1) +
        geom_abline(data = lib_compare,
                    aes(slope = deming_slope, intercept = deming_intercept),
                    color = "#e41a1c", linewidth = 0.7) +
        facet_grid(gas ~ location, scales = "free") +
        labs(title    = "FTIR.4_old vs FTIR.4 after spectral-library correction",
             subtitle = "Dashed grey = 1:1 line; red = Deming regression",
             x = "FTIR.4_old (ppm)", y = "FTIR.4 (ppm)") +
        theme_bw(base_size = 12)
ggsave(file.path(plots_dir, "FTIR4_old_vs_new_scatter.png"), lib_plot,
       width = 11, height = 9, dpi = 300, bg = "white")

#### 8.  FTIR.2 absolute-NH3 diagnostic  (R2.7 — show why it is excluded)     ####
# Visual check that FTIR.2 absolute NH3 sits well above the consortium and that
# the offset is roughly multiplicative — i.e. it cancels in differences.
ftir2_nh3 <- emission_reshaped %>%
        filter(var == "NH3_mgm3",
               analyzer %in% c("FTIR.1","FTIR.2","FTIR.3","FTIR.4","CRDS.1","CRDS.2","CRDS.3"))
ftir2_diag <- ggplot(ftir2_nh3, aes(x = DATE.TIME, y = value, color = analyzer)) +
        geom_line(linewidth = 0.4, na.rm = TRUE) +
        facet_grid(location ~ ., scales = "free_y",
                   labeller = labeller(location = as_labeller(LOC_LABELS, label_parsed))) +
        scale_color_manual(values = ANALYZER_COLORS) +
        scale_x_datetime(date_breaks = "1 day", date_labels = "%m-%d") +
        labs(title = "Absolute NH3: FTIR.2 inflated vs consortium",
             subtitle = "Offset is approximately multiplicative -> cancels in ΔNH3 (see §10)",
             x = NULL, y = expression(c[NH3]~"(mg m"^-3*")")) +
        theme_bw(base_size = 12) + theme(legend.position = "bottom")
ggsave(file.path(plots_dir, "FTIR2_absolute_NH3_diagnostic.png"), ftir2_diag,
       width = 11, height = 8, dpi = 300, bg = "white")

# FTIR.2 vs the FTIR.2-excluded baseline, on raw and on delta NH3, to confirm
# the offset cancels in delta. Indoor only (cleanest signal).
ftir2_indoor_raw <- emission_reshaped %>%
        filter(var == "NH3_mgm3", location == "Indoor",
               analyzer %in% c("FTIR.2","baseline"))
ftir2_indoor_delta <- emission_reshaped %>%
        filter(var == "delta_NH3", location == "Outdoor_NE",
               analyzer %in% c("FTIR.2","baseline"))
ftir2_cancel_plot <- (ggplot(ftir2_indoor_raw, aes(DATE.TIME, value, color = analyzer)) +
                              geom_line(linewidth = 0.5) +
                              labs(subtitle = "Indoor absolute NH3", x = NULL, y = NULL) +
                              theme_bw() + theme(legend.position = "none")) /
        (ggplot(ftir2_indoor_delta, aes(DATE.TIME, value, color = analyzer)) +
                 geom_line(linewidth = 0.5) +
                 geom_hline(yintercept = 0, color = "grey60") +
                 labs(subtitle = "ΔNH3 = Indoor - Outdoor_NE", x = NULL, y = NULL) +
                 theme_bw() + theme(legend.position = "bottom")) +
        plot_annotation(title = "FTIR.2 NH3: large bias in absolute, cancels in Δ")
ggsave(file.path(plots_dir, "FTIR2_NH3_offset_cancels_in_delta.png"), ftir2_cancel_plot,
       width = 10, height = 7, dpi = 300, bg = "white")

#### 9.  Drop FTIR.4_old; rebuild working datasets  (R2.7 step 1)              ####
emission_result_v  <- emission_result %>% filter(analyzer != "FTIR.4_old")

#### 9.5  Stage 2 — mask Q/e/Δc cells when any Δc is negative (handover §2.6c) ####
# Steady-state mass balance: c_indoor >= c_outdoor for emitting gases.
# A negative Δc indicates either a measurement artefact (analyser noise,
# adsorption lag) or a transient mass-balance violation (e.g. brief
# outdoor plume from the western dairies blowing across the indoor ring
# line). Per handover §2.6c the mask is per-(analyser, hour, outdoor):
# when any of CO2/CH4/NH3 has Δc < 0 for an outdoor line in a given row,
# only that outdoor's Δc / Q / e cells are masked. Other analysers' rows
# at the same hour stay intact; the same row's _N data are unaffected
# when only _S has a negative Δc and vice versa.

stage2_dropouts <- emission_result_v %>%
        select(DATE.TIME, lab, analyzer,
               matches("^delta_(CO2|CH4|NH3)_(N|S)$")) %>%
        pivot_longer(matches("^delta_"),
                     names_to      = c("gas", "outdoor"),
                     names_pattern = "^delta_(CO2|CH4|NH3)_(N|S)$",
                     values_to     = "delta") %>%
        mutate(outdoor = recode(outdoor,
                                "N" = "Outdoor_NE", "S" = "Outdoor_SW")) %>%
        group_by(analyzer, gas, outdoor) %>%
        summarise(n_total      = n(),
                  n_negative   = sum(delta < 0, na.rm = TRUE),
                  pct_negative = round(100 * n_negative / n_total, 2),
                  .groups      = "drop")
write_excel_csv(stage2_dropouts, file.path(tables_dir, "stage2_dropouts.csv"))

mask_cols_N <- c("delta_CO2_N",   "delta_CH4_N", "delta_NH3_N",
                 "Q_vent_N",
                 "e_CH4_gh_N",    "e_NH3_gh_N",
                 "e_CH4_ghLU_N",  "e_NH3_ghLU_N")
mask_cols_S <- gsub("_N$", "_S", mask_cols_N)

emission_result_v <- emission_result_v %>%
        mutate(.neg_N = (delta_CO2_N < 0) | (delta_CH4_N < 0) | (delta_NH3_N < 0),
               .neg_S = (delta_CO2_S < 0) | (delta_CH4_S < 0) | (delta_NH3_S < 0)) %>%
        mutate(across(any_of(mask_cols_N),
                      ~ ifelse(.neg_N %in% TRUE, NA_real_, .x)),
               across(any_of(mask_cols_S),
                      ~ ifelse(.neg_S %in% TRUE, NA_real_, .x))) %>%
        select(-.neg_N, -.neg_S)

emission_reshaped_v <- reshaper(emission_result_v) %>%
        mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
        mutate(analyzer = fct_drop(analyzer))
input_combined_v   <- input_combined %>% filter(analyzer != "FTIR.4_old")
analyzers          <- sort(unique(input_combined_v$analyzer))

#### 10. NE vs SW outdoor-line comparison  (R2.24 / R2.8 / E.1)                ####
# Step A: paired test per gas & analyzer (campaign-wide).
ne_sw_tests <- map_dfr(gases, function(g) {
        map_dfr(as.character(analyzers), function(a) {
                d  <- input_combined_v %>% filter(analyzer == a)
                ne <- d[[paste0(g, "_ppm_N")]]; sw <- d[[paste0(g, "_ppm_S")]]
                ok <- complete.cases(ne, sw); ne <- ne[ok]; sw <- sw[ok]
                if (length(ne) < 3) return(NULL)
                tibble(gas = g, analyzer = a, n = length(ne),
                       mean_NE = mean(ne), mean_SW = mean(sw),
                       mean_diff = mean(ne - sw),
                       RPD_pct = 100 * mean(ne - sw) / mean((ne + sw) / 2),
                       p_ttest  = t.test(ne, sw, paired = TRUE)$p.value,
                       p_wilcox = suppressWarnings(
                               wilcox.test(ne, sw, paired = TRUE)$p.value))
        })
}) %>%
        mutate(p_wilcox_holm = p.adjust(p_wilcox, method = "holm"),
               significant   = ifelse(p_wilcox_holm < 0.05, "yes", "no"))
write_excel_csv(ne_sw_tests, file.path(tables_dir, "outdoor_NE_vs_SW_tests.csv"))

# Step B: wind-sector-conditioned NE vs SW. For each gas and 8-sector wind
# bin, average NE and SW concentrations across all analyzers and compare.
ws_long_ne_sw <- input_combined_v %>%
        filter(!is.na(wd_mst)) %>%
        mutate(wind_sector = deg_to_compass8(wd_mst)) %>%
        select(DATE.TIME, analyzer, wind_sector,
               matches("^(CO2|CH4|NH3)_ppm_(N|S)$")) %>%
        pivot_longer(cols = matches("_ppm_"),
                     names_to = c("gas", "loc"),
                     names_pattern = "(.+)_ppm_(.+)",
                     values_to = "ppm") %>%
        mutate(location = unname(loc_from_suffix[loc]),
               gas      = factor(gas, levels = gases))

ne_sw_by_sector <- ws_long_ne_sw %>%
        group_by(gas, location, wind_sector) %>%
        summarise(mean_ppm = mean(ppm, na.rm = TRUE),
                  sd_ppm   = sd(ppm,   na.rm = TRUE),
                  n        = n(), .groups = "drop")
write_excel_csv(ne_sw_by_sector, file.path(tables_dir, "outdoor_NE_vs_SW_by_wind_sector.csv"))

# Annotate each sector with the contamination class of each line.
sector_class <- tibble(
        wind_sector = factor(c("N","NE","E","SE","S","SW","W","NW"),
                             levels = c("N","NE","E","SE","S","SW","W","NW")),
        NE_class = case_when(
                wind_sector %in% NE_CONTAM_NEIGHBOUR ~ "neighbour_source",
                wind_sector %in% NE_CONTAM_LVAT      ~ "downwind_LVAT",
                wind_sector %in% NE_CLEAN_SECTORS    ~ "clean",
                TRUE                                 ~ "other"),
        SW_class = case_when(
                wind_sector %in% SW_CONTAM_NEIGHBOUR ~ "neighbour_source",
                wind_sector %in% SW_CONTAM_LVAT      ~ "downwind_LVAT",
                wind_sector %in% SW_CLEAN_SECTORS    ~ "clean",
                TRUE                                 ~ "other"))
write_excel_csv(sector_class, file.path(tables_dir, "wind_sector_contamination_class.csv"))

# Step C: per-line contamination summary — mean concentration in each class.
line_class_summary <- ne_sw_by_sector %>%
        left_join(sector_class, by = "wind_sector") %>%
        mutate(class_for_line = ifelse(location == "Outdoor_NE", NE_class, SW_class)) %>%
        group_by(gas, location, class_for_line) %>%
        summarise(mean_ppm = mean(mean_ppm, na.rm = TRUE),
                  n_sectors = n_distinct(wind_sector),
                  .groups = "drop")
write_excel_csv(line_class_summary, file.path(tables_dir, "outdoor_line_class_summary.csv"))

# Step D: side-by-side polar plot, mean(NE) and mean(SW) by sector, per gas.
ne_sw_polar_data <- ne_sw_by_sector
ne_sw_polar <- ggplot(ne_sw_polar_data,
                      aes(x = wind_sector, y = mean_ppm, fill = location)) +
        geom_col(position = position_dodge(width = 0.8), width = 0.7,
                 color = "black", linewidth = 0.2) +
        coord_polar(start = -pi / 8) +
        facet_wrap(~ gas, scales = "free_y") +
        scale_fill_manual(values = c("Outdoor_NE" = "#377EB8",
                                     "Outdoor_SW" = "#E41A1C")) +
        labs(title = "Outdoor concentration by wind sector",
             subtitle = "Each bar = campaign-mean across all analyzers",
             x = NULL, y = "Mean concentration (ppm)") +
        theme_bw(base_size = 12) +
        theme(axis.text.x = element_text(size = 9),
              legend.position = "bottom", legend.title = element_blank())
ggsave(file.path(plots_dir, "outdoor_NE_vs_SW_polar.png"), ne_sw_polar,
       width = 12, height = 5, dpi = 300, bg = "white")

# Step E: bar plot of RPD vs wind sector to make the contamination pattern
# easy to read for the manuscript.
ne_sw_rpd_sector <- ne_sw_by_sector %>%
        select(gas, location, wind_sector, mean_ppm) %>%
        pivot_wider(names_from = location, values_from = mean_ppm) %>%
        mutate(RPD_pct = 100 * (Outdoor_NE - Outdoor_SW) /
                        ((Outdoor_NE + Outdoor_SW) / 2)) %>%
        left_join(sector_class, by = "wind_sector")
rpd_sector_plot <- ggplot(ne_sw_rpd_sector, aes(x = wind_sector, y = RPD_pct, fill = gas)) +
        geom_col(position = position_dodge(width = 0.8), width = 0.7,
                 color = "black", linewidth = 0.2) +
        geom_hline(yintercept = 0, color = "grey50") +
        facet_wrap(~ gas, ncol = 1, scales = "free_y") +
        scale_fill_brewer(palette = "Dark2", guide = "none") +
        labs(title = "Relative percentage difference NE - SW by wind sector",
             subtitle = "Positive = NE higher than SW; negative = SW higher than NE",
             x = "Wind sector (mast)", y = "RPD (%)") +
        theme_bw(base_size = 12)
ggsave(file.path(plots_dir, "outdoor_NE_vs_SW_RPD_by_sector.png"), rpd_sector_plot,
       width = 10, height = 9, dpi = 300, bg = "white")

#### 11. Decide which outdoor line to retain                                  ####
# Decision rule:
#   * The retained line is the one whose mean over its "clean" sectors is
#     closer to its mean over the "downwind_LVAT" sectors AND that has the
#     lower neighbour-source contamination, weighted by the time fraction
#     each line spends in each class during this campaign.
#   * In a tie, prefer the line on the upwind side of the prevailing wind.
# Practically this is read off the line_class_summary + wind-rose
# frequencies below; the actual drop is written here so all downstream
# tables stay consistent.

wind_freq <- input_combined_v %>%
        distinct(DATE.TIME, wd_mst) %>%
        filter(!is.na(wd_mst)) %>%
        mutate(wind_sector = deg_to_compass8(wd_mst)) %>%
        count(wind_sector, .drop = FALSE) %>%
        mutate(freq = n / sum(n))
write_excel_csv(wind_freq, file.path(tables_dir, "wind_sector_frequency.csv"))

# Time spent with each line in each class.
line_time_class <- wind_freq %>%
        left_join(sector_class, by = "wind_sector") %>%
        pivot_longer(cols = c(NE_class, SW_class),
                     names_to = "line", values_to = "class") %>%
        mutate(line = recode(line, "NE_class" = "Outdoor_NE", "SW_class" = "Outdoor_SW")) %>%
        group_by(line, class) %>% summarise(time_frac = sum(freq), .groups = "drop")
write_excel_csv(line_time_class, file.path(tables_dir, "outdoor_line_time_in_class.csv"))

contam_score <- line_time_class %>%
        filter(class %in% c("downwind_LVAT", "neighbour_source")) %>%
        group_by(line) %>% summarise(contam_time = sum(time_frac), .groups = "drop") %>%
        arrange(contam_time)
RETAINED_OUTDOOR <- contam_score$line[1]
DROPPED_OUTDOOR  <- contam_score$line[contam_score$line != RETAINED_OUTDOOR][1]

# Persist the decision for the manuscript.
cat(sprintf("retained_outdoor: %s\ndropped_outdoor:  %s\n",
            RETAINED_OUTDOOR, DROPPED_OUTDOOR),
    file = file.path(tables_dir, "outdoor_line_choice.txt"))

#### 11.5 Per-analyser IQR filter on Q_vent + dropout diagnostic              ####
# Replaces the V10 [MIN_Q_VENT, MAX_Q_VENT] sanity bound that was previously
# applied inside indirect.CO2.balance(). Reviewer-requested change: caps were
# silent and conflated genuine low-gradient hours (already handled by
# MIN_DELTA_CO2 in §2.3) with per-analyser anomalies (especially FTIR.4).
# The new filter applies 1.5*IQR per analyser to the surviving Q_vent values
# and records the number of rows lost at every stage so the loss can be
# reported in Results §3.2 and Discussion §4.4.

# Count inputs and MIN_DELTA_CO2 dropouts per (analyzer, outdoor line).
# An input row is a row of emission_result_v with a finite delta_CO2_*
# (i.e. the row had both indoor and outdoor CO2 measurements). It is
# dropped by MIN_DELTA_CO2 when delta_CO2_* < MIN_DELTA_CO2.
qvent_dropouts_min_delta <- emission_result_v %>%
        group_by(analyzer) %>%
        summarise(
                n_input_N        = sum(is.finite(delta_CO2_N)),
                n_dropped_minD_N = sum(is.finite(delta_CO2_N) & delta_CO2_N <  MIN_DELTA_CO2),
                n_after_minD_N   = sum(is.finite(delta_CO2_N) & delta_CO2_N >= MIN_DELTA_CO2),
                n_input_S        = sum(is.finite(delta_CO2_S)),
                n_dropped_minD_S = sum(is.finite(delta_CO2_S) & delta_CO2_S <  MIN_DELTA_CO2),
                n_after_minD_S   = sum(is.finite(delta_CO2_S) & delta_CO2_S >= MIN_DELTA_CO2),
                .groups = "drop")

# Apply per-analyser 1.5*IQR removal on Q_vent_N and Q_vent_S.
iqr_drop <- function(x) {
        if (sum(is.finite(x)) < 4) return(x)
        q <- quantile(x, c(0.25, 0.75), na.rm = TRUE)
        H <- 1.5 * IQR(x, na.rm = TRUE)
        ifelse(x < (q[1] - H) | x > (q[2] + H), NA_real_, x)
}
emission_result_v <- emission_result_v %>%
        group_by(analyzer) %>%
        mutate(
                Q_vent_N_pre_iqr = Q_vent_N,
                Q_vent_S_pre_iqr = Q_vent_S,
                Q_vent_N         = iqr_drop(Q_vent_N),
                Q_vent_S         = iqr_drop(Q_vent_S)
        ) %>%
        ungroup() %>%
        # Re-propagate the IQR-filtered Q into the derived emission columns.
        mutate(
                e_NH3_gh_N   = (delta_NH3_N * Q_vent_N / 1000) * n_dairycows_in,
                e_CH4_gh_N   = (delta_CH4_N * Q_vent_N / 1000) * n_dairycows_in,
                e_NH3_gh_S   = (delta_NH3_S * Q_vent_S / 1000) * n_dairycows_in,
                e_CH4_gh_S   = (delta_CH4_S * Q_vent_S / 1000) * n_dairycows_in,
                e_NH3_ghLU_N = (e_NH3_gh_N * 500) / (n_dairycows_in * m_weight),
                e_CH4_ghLU_N = (e_CH4_gh_N * 500) / (n_dairycows_in * m_weight),
                e_NH3_ghLU_S = (e_NH3_gh_S * 500) / (n_dairycows_in * m_weight),
                e_CH4_ghLU_S = (e_CH4_gh_S * 500) / (n_dairycows_in * m_weight)
        )

qvent_dropouts_iqr <- emission_result_v %>%
        group_by(analyzer) %>%
        summarise(
                n_after_minD_N   = sum(is.finite(Q_vent_N_pre_iqr)),
                n_dropped_iqr_N  = sum(is.finite(Q_vent_N_pre_iqr) & !is.finite(Q_vent_N)),
                n_kept_N         = sum(is.finite(Q_vent_N)),
                n_after_minD_S   = sum(is.finite(Q_vent_S_pre_iqr)),
                n_dropped_iqr_S  = sum(is.finite(Q_vent_S_pre_iqr) & !is.finite(Q_vent_S)),
                n_kept_S         = sum(is.finite(Q_vent_S)),
                .groups = "drop")

qvent_dropouts <- qvent_dropouts_min_delta %>%
        select(analyzer, n_input_N, n_dropped_minD_N,
               n_input_S, n_dropped_minD_S) %>%
        left_join(qvent_dropouts_iqr, by = "analyzer") %>%
        mutate(
                pct_kept_N = 100 * n_kept_N / pmax(n_input_N, 1),
                pct_kept_S = 100 * n_kept_S / pmax(n_input_S, 1)
        ) %>%
        select(analyzer,
               n_input_N, n_dropped_minD_N, n_dropped_iqr_N, n_kept_N, pct_kept_N,
               n_input_S, n_dropped_minD_S, n_dropped_iqr_S, n_kept_S, pct_kept_S)

write_excel_csv(qvent_dropouts,
                file.path(tables_dir, "qvent_filter_dropouts.csv"))

#### 12. Build single-outdoor working dataset                                 ####
# Map RETAINED -> suffix column for *_N / *_S in emission_result_v
retain_suffix  <- c("Outdoor_NE" = "N", "Outdoor_SW" = "S")[[RETAINED_OUTDOOR]]
drop_suffix    <- c("Outdoor_NE" = "N", "Outdoor_SW" = "S")[[DROPPED_OUTDOOR]]
# Keep only the retained-outdoor delta/Q/e columns in emission_result_single
emission_result_single <- emission_result_v %>%
        select(-matches(paste0("_", drop_suffix, "$"))) %>%
        rename_with(~ str_replace(.x, paste0("_", retain_suffix, "$"), ""),
                    .cols = matches(paste0("_", retain_suffix, "$")))

# Reshape — single outdoor, so location is Indoor + retained outdoor only.
# We feed reshaper a frame with three "location-coded" CO2/CH4/NH3 columns
# but with delta/Q/e already collapsed; pivot manually.
single_long <- emission_result_v %>%
        select(DATE.TIME, lab, analyzer,
               CO2_mgm3_in, CH4_mgm3_in, NH3_mgm3_in,
               any_of(paste0(c("CO2_mgm3","CH4_mgm3","NH3_mgm3"), "_", retain_suffix)),
               any_of(paste0(c("delta_CO2","delta_CH4","delta_NH3",
                                "Q_vent","e_CH4_gh","e_NH3_gh",
                                "e_CH4_ghLU","e_NH3_ghLU"), "_", retain_suffix))) %>%
        rename_with(~ str_replace(.x, paste0("_", retain_suffix, "$"), "_OUT"),
                    .cols = matches(paste0("_", retain_suffix, "$"))) %>%
        pivot_longer(cols = -c(DATE.TIME, lab, analyzer),
                     names_to = "raw_var", values_to = "value") %>%
        mutate(location = case_when(
                str_detect(raw_var, "_in$")  ~ "Indoor",
                str_detect(raw_var, "_OUT$") ~ RETAINED_OUTDOOR,
                TRUE                         ~ "Combined"),
               var = str_remove(raw_var, "_(in|OUT)$")) %>%
        select(DATE.TIME, lab, analyzer, location, var, value)

# Per-variable baseline on the single-outdoor frame (same exclusion rules).
single_baseline <- single_long %>%
        group_by(DATE.TIME, location, var) %>%
        summarise(
                excl  = list(if (!is.null(BASELINE_EXCLUDE_BY_VAR[[first(var)]]))
                                BASELINE_EXCLUDE_BY_VAR[[first(var)]]
                             else BASELINE_EXCLUDE_DEFAULT),
                value = mean(value[!analyzer %in% excl[[1]]], na.rm = TRUE),
                .groups = "drop") %>%
        select(-excl) %>%
        mutate(analyzer = "baseline", lab = "baseline")

single_long <- bind_rows(single_long, single_baseline) %>%
        mutate(analyzer = factor(analyzer,
                                 levels = c("FTIR.1","FTIR.2","FTIR.3","FTIR.4",
                                            "CRDS.1","CRDS.2","CRDS.3","baseline")))
write_excel_csv(single_long, file.path(tables_dir,
                                       sprintf("emission_long_single_outdoor_%s.csv",
                                               RETAINED_OUTDOOR)))

#### 13. Delta / Q / e plots — single outdoor                                 ####
single_for_plot <- single_long %>%
        mutate(day  = factor(as.Date(DATE.TIME)),
               hour = factor(format(DATE.TIME, "%H:%M")))

d_trend_plot   <- emitrendplot(single_for_plot,
                               y = c("delta_CO2","delta_CH4","delta_NH3"),
                               location_filter = RETAINED_OUTDOOR)
d_boxplot      <- emiboxplot  (single_for_plot,
                               y = c("delta_CO2","delta_CH4","delta_NH3"),
                               location_filter = RETAINED_OUTDOOR)
q_e_trend_plot <- emitrendplot(single_for_plot,
                               y = c("Q_vent","e_CH4_ghLU","e_NH3_ghLU"),
                               location_filter = RETAINED_OUTDOOR)
q_e_boxplot    <- emiboxplot  (single_for_plot,
                               y = c("Q_vent","e_CH4_ghLU","e_NH3_ghLU"),
                               location_filter = RETAINED_OUTDOOR)

ggsave(file.path(plots_dir, "d_trend_plot.png"),   d_trend_plot,
       width = 12, height = 8, dpi = 300)
ggsave(file.path(plots_dir, "d_boxplot.png"),      d_boxplot,
       width = 12, height = 8, dpi = 300)
ggsave(file.path(plots_dir, "q_e_trend_plot.png"), q_e_trend_plot,
       width = 12, height = 8, dpi = 300)
ggsave(file.path(plots_dir, "q_e_boxplot.png"),    q_e_boxplot,
       width = 12, height = 8, dpi = 300)

#### 14. Bland-Altman, intra-lab shared sampling line — single outdoor        ####
# Relative BA only on pairs that share the same sampling line (Lab A, Lab B,
# Lab C-only-CRDS so skipped).
save_bland_altman <- function(analyzer_pair, tag) {
        mk <- function(v) {
                bland_altman_plot(single_for_plot, var_filter = v,
                                  analyzer_pair  = analyzer_pair,
                                  location_filter = RETAINED_OUTDOOR) +
                        theme(plot.margin = margin(10, 10, 10, 10))
        }
        e_panels <- list(mk("e_CH4_ghLU"), mk("e_NH3_ghLU"))
        ggsave(file.path(plots_dir, paste0("e_BlandAltman_", tag, ".png")),
               wrap_plots(e_panels, ncol = 2, nrow = 1),
               width = 10, height = 4, units = "in", dpi = 300)
        q_panel <- mk("Q_vent")
        ggsave(file.path(plots_dir, paste0("q_BlandAltman_", tag, ".png")),
               q_panel, width = 6, height = 4, units = "in", dpi = 300)
}
save_bland_altman(c("FTIR.1", "CRDS.1"), "AnalyzerA")
save_bland_altman(c("FTIR.2", "CRDS.2"), "AnalyzerB")

ba_pairs <- list(AnalyzerA = c("FTIR.1", "CRDS.1"),
                 AnalyzerB = c("FTIR.2", "CRDS.2"))
ba_vars  <- c("e_CH4_ghLU", "e_NH3_ghLU", "Q_vent")
ba_table <- map_dfr(names(ba_pairs), function(tag) {
        pr <- ba_pairs[[tag]]
        map_dfr(ba_vars, function(v) {
                w <- single_for_plot %>%
                        filter(var == v, analyzer %in% pr) %>%
                        select(DATE.TIME, analyzer, value) %>%
                        pivot_wider(names_from = analyzer, values_from = value)
                st <- ba_relative_stats(w[[pr[1]]], w[[pr[2]]])
                if (is.null(st)) return(NULL)
                bind_cols(tibble(pair = tag, variable = v), st)
        })
})
write_excel_csv(ba_table, file.path(tables_dir, "bland_altman_relative.csv"))

#### 15. Pairwise method comparison — extended to absolute c + Δc            ####
# R2.9 / R1.13: report Deming slope+intercept + Lin's CCC + Pearson r for
# every analyzer pair, not just for Q/e. Also exposes Δc agreement.
pairwise_compare <- function(df_long, var_sel, loc_sel) {
        w <- df_long %>%
                filter(var == var_sel, location == loc_sel,
                       !analyzer %in% c("baseline","HOBO","USA","RGB")) %>%
                select(DATE.TIME, analyzer, value) %>%
                mutate(analyzer = as.character(analyzer)) %>%
                pivot_wider(names_from = analyzer, values_from = value,
                            values_fn = ~ mean(.x, na.rm = TRUE))
        ana <- setdiff(names(w), "DATE.TIME")
        if (length(ana) < 2) return(NULL)
        map_dfr(combn(ana, 2, simplify = FALSE), function(pr) {
                x <- w[[pr[1]]]; y <- w[[pr[2]]]
                ok <- complete.cases(x, y); x <- x[ok]; y <- y[ok]
                if (length(x) < 5) return(NULL)
                dem <- deming_fit(x, y)
                lm_fit <- lm(y ~ x)
                tibble(var = var_sel, location = loc_sel,
                       analyzer_x = pr[1], analyzer_y = pr[2], n = length(x),
                       pearson_r        = cor(x, y),
                       ccc              = DescTools::CCC(x, y)$rho.c[, "est"],
                       deming_slope     = unname(dem["slope"]),
                       deming_intercept = unname(dem["intercept"]),
                       ols_slope        = unname(coef(lm_fit)[2]),
                       ols_intercept    = unname(coef(lm_fit)[1]),
                       r_squared        = summary(lm_fit)$r.squared)
        })
}
# Inputs: emission_reshaped_v carries both location-coded absolutes and deltas
absolute_vars <- c("CO2_mgm3","CH4_mgm3","NH3_mgm3")
delta_vars    <- c("delta_CO2","delta_CH4","delta_NH3")
qe_vars       <- c("Q_vent","e_CH4_ghLU","e_NH3_ghLU")
locs_abs      <- c("Indoor","Outdoor_NE","Outdoor_SW")
locs_delta_qe <- c("Outdoor_NE","Outdoor_SW")

pairwise_abs <- map_dfr(absolute_vars, function(v)
        map_dfr(locs_abs, function(l) pairwise_compare(emission_reshaped_v, v, l)))
pairwise_delta <- map_dfr(delta_vars, function(v)
        map_dfr(locs_delta_qe, function(l) pairwise_compare(emission_reshaped_v, v, l)))
pairwise_qe <- map_dfr(qe_vars, function(v)
        map_dfr(locs_delta_qe, function(l) pairwise_compare(emission_reshaped_v, v, l)))
pairwise_tbl <- bind_rows(pairwise_abs, pairwise_delta, pairwise_qe)
write_excel_csv(pairwise_tbl, file.path(tables_dir, "pairwise_regression_ccc.csv"))

# CCC heatmaps for each variable group
heatmap_ccc <- function(tbl, var_set, plain_labels, title, file) {
        d <- tbl %>%
                filter(var %in% var_set) %>%
                mutate(facet_label = factor(plain_labels[var], levels = plain_labels[var_set]))
        p <- ggplot(d, aes(x = analyzer_x, y = analyzer_y, fill = ccc)) +
                geom_tile(color = "white") +
                geom_text(aes(label = sprintf("%.2f", ccc)), size = 2.6) +
                scale_fill_gradient2(low = "#b2182b", mid = "#f7f7f7", high = "#2166ac",
                                     midpoint = 0.5, limits = c(0, 1), name = "Lin's CCC") +
                facet_grid(location ~ facet_label,
                           labeller = labeller(facet_label = label_parsed,
                                               location    = as_labeller(LOC_LABELS, label_parsed))) +
                labs(title = title, x = NULL, y = NULL) +
                theme_bw(base_size = 11) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1),
                      legend.position = "bottom")
        ggsave(file.path(plots_dir, file), p, width = 12, height = 9, dpi = 300, bg = "white")
}
heatmap_ccc(pairwise_abs, absolute_vars, VAR_LABELS_PLAIN,
            "Pairwise agreement, absolute concentrations (Lin's CCC)",
            "pairwise_ccc_absolute_c.png")
heatmap_ccc(pairwise_delta, delta_vars, VAR_LABELS_PLAIN,
            "Pairwise agreement, Δ concentrations (Lin's CCC)",
            "pairwise_ccc_delta_c.png")
heatmap_ccc(pairwise_qe, qe_vars, VAR_LABELS_PLAIN,
            "Pairwise agreement, ventilation + emissions (Lin's CCC)",
            "pairwise_ccc_q_e.png")

#### 16. Tukey HSD pairwise tables  (Tables 3, 4, 6 of the manuscript)        ####
# One-way ANOVA per (variable, location) on the analyzer factor, followed by
# Tukey's honestly significant difference post-hoc. Each pair is reported with
# the per-analyzer means, the relative difference (RD %), the Tukey-adjusted
# p-value, a Holm correction over the panel (gas x location), and a
# significance code (ns / * / ** / ***).
#
# Three tables are produced:
#   Table 3 — absolute concentrations  (c_CO2, c_CH4, c_NH3) x 3 locations
#   Table 4 — concentration differences (delta_CO2, delta_CH4, delta_NH3) x 2 outdoors
#   Table 6 — derived quantities       (Q_vent, e_CH4_ghLU, e_NH3_ghLU) x 2 outdoors
#
# Inputs: emission_reshaped_v (post-FTIR.4_old drop). FTIR.2 is kept because
# the absolute-NH3 offset cancels in deltas — see Section 8 diagnostics.

# Significance code helper (ns / * / ** / ***)
sig_code <- function(p) {
        ifelse(is.na(p), "",
        ifelse(p < 0.001, "***",
        ifelse(p < 0.01,  "**",
        ifelse(p < 0.05,  "*",  "ns"))))
}

# Run Tukey HSD for one (variable, location) panel and return long-form rows.
tukey_panel <- function(df_long, var_sel, loc_sel) {
        x <- df_long %>%
                filter(var == var_sel, location == loc_sel,
                       !analyzer %in% c("baseline","HOBO","USA","RGB"),
                       !is.na(value), is.finite(value)) %>%
                mutate(analyzer = as.character(analyzer))
        if (nrow(x) < 10 || length(unique(x$analyzer)) < 2) return(NULL)
        an  <- factor(x$analyzer)
        fit <- aov(value ~ an, data = data.frame(value = x$value, an = an))
        tk  <- TukeyHSD(fit, conf.level = 0.95)$an
        means <- x %>% group_by(analyzer) %>%
                summarise(m = mean(value, na.rm = TRUE), n = n(), .groups = "drop")
        pairs <- rownames(tk)
        map_dfr(seq_along(pairs), function(i) {
                ab <- strsplit(pairs[i], "-", fixed = TRUE)[[1]]
                ax <- ab[1]; ay <- ab[2]
                mx <- means$m[means$analyzer == ax]
                my <- means$m[means$analyzer == ay]
                if (length(mx) == 0 || length(my) == 0) return(NULL)
                pm <- (mx + my) / 2
                rd <- ifelse(pm != 0, 100 * (mx - my) / pm, NA_real_)
                tibble(variable = var_sel, location = loc_sel,
                       analyzer_1 = ax, analyzer_2 = ay,
                       n_1     = means$n[means$analyzer == ax],
                       n_2     = means$n[means$analyzer == ay],
                       mean_1  = round(mx, 3), mean_2 = round(my, 3),
                       diff    = round(tk[i, "diff"], 3),
                       lwr     = round(tk[i, "lwr"], 3),
                       upr     = round(tk[i, "upr"], 3),
                       RD_pct  = round(rd, 1),
                       p_tukey = signif(tk[i, "p adj"], 3))
        })
}

# Walk the (variable x location) grid, apply Tukey, add Holm per panel + sig.
tukey_table <- function(df_long, vars, locs) {
        out <- map_dfr(vars, function(v) {
                map_dfr(locs, function(l) tukey_panel(df_long, v, l))
        })
        if (nrow(out) == 0) return(out)
        out %>%
                group_by(variable, location) %>%
                mutate(p_holm = signif(p.adjust(p_tukey, method = "holm"), 3)) %>%
                ungroup() %>%
                mutate(sig = sig_code(p_holm)) %>%
                arrange(variable, location, analyzer_1, analyzer_2)
}

tukey_abs   <- tukey_table(emission_reshaped_v,
                           vars = c("CO2_mgm3","CH4_mgm3","NH3_mgm3"),
                           locs = c("Indoor","Outdoor_NE","Outdoor_SW"))
tukey_delta <- tukey_table(emission_reshaped_v,
                           vars = c("delta_CO2","delta_CH4","delta_NH3"),
                           locs = c("Outdoor_NE","Outdoor_SW"))
tukey_qe    <- tukey_table(emission_reshaped_v,
                           vars = c("Q_vent","e_CH4_ghLU","e_NH3_ghLU"),
                           locs = c("Outdoor_NE","Outdoor_SW"))

write_excel_csv(tukey_abs,   file.path(tables_dir, "tukey_absolute_concentrations.csv"))
write_excel_csv(tukey_delta, file.path(tables_dir, "tukey_delta_concentrations.csv"))
write_excel_csv(tukey_qe,    file.path(tables_dir, "tukey_ventilation_emission.csv"))
# Long-form CSV layout (one row per analyser pair):
#   variable, location, analyzer_1, analyzer_2, n_1, n_2, mean_1, mean_2,
#   diff (mean_1 - mean_2), lwr/upr (95% CI on diff from Tukey),
#   RD_pct = 100 * diff / pair_mean,
#   p_tukey, p_holm (Holm-adjusted across the (variable, location) panel),
#   sig (ns / * / ** / ***).

#### 17. Wind characterisation — daily campaign + seasonal historical         ####
# Two polar-bar plots. In both, bars are stacked by wind-speed bin and each
# sector is annotated (in blue) with its hourly count and per-panel % share.
#
#   wind_rose_daily.png    — 7 panels, one per campaign day
#                            (08.04.2025 to 14.04.2025; 15.04 deliberately
#                            excluded because the campaign ends 14.04 12:00)
#   wind_rose_seasonal.png — 4 panels, historical seasons in calendar order
#                            (Summer 2024 -> Autumn 2024 -> Winter 2024/25
#                            -> Spring 2025 up to 07.04.2025)

# Shared wind-speed bins + colour palette (Figure-2 reference style:
# burgundy at low speeds -> blue at high speeds).
WS_BREAKS  <- c(0, 0.5, 1.0, 2.0, 3.0, 4.0, Inf)
WS_LABELS  <- c("0.0-0.5","0.5-1.0","1.0-2.0","2.0-3.0","3.0-4.0",">=4.0")
WS_COLOURS <- c("0.0-0.5" = "#8C1838",   # dark wine
                "0.5-1.0" = "#E04344",   # red
                "1.0-2.0" = "#F49649",   # orange
                "2.0-3.0" = "#E3E394",   # pale yellow-green
                "3.0-4.0" = "#4DBC8E",   # teal-green
                ">=4.0"   = "#3E81BA")   # blue

# Read the FULL wind series (campaign window + historical year).
wind_full <- read.csv(file.path(meta_dir, "USA_mast_wind/20240101_20250825_USA_mast_16_hourly_uvw_wd_ws.csv"),
                      stringsAsFactors = FALSE) %>%
        mutate(DATE.TIME = as.POSIXct(datetime_hour, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) %>%
        rename(wd_mst = wd, ws_mst = ws) %>%
        filter(!is.na(wd_mst), !is.na(ws_mst)) %>%
        mutate(wind_sector = deg_to_compass8(wd_mst),
               ws_bin      = cut(ws_mst, breaks = WS_BREAKS, labels = WS_LABELS,
                                 include.lowest = TRUE, right = FALSE))

# --- (a) Daily wind roses for the 7 campaign days ----------------------------
daily_data <- wind_full %>%
        filter(DATE.TIME >= start_time, DATE.TIME <= end_time) %>%
        mutate(day_label = format(as.Date(DATE.TIME), "%d.%m.%Y")) %>%
        group_by(day_label) %>% mutate(day_total = n()) %>% ungroup()

daily_stack <- daily_data %>%
        group_by(day_label, day_total, wind_sector, ws_bin) %>%
        summarise(n_bin = n(), .groups = "drop") %>%
        mutate(pct = 100 * n_bin / day_total)

daily_totals <- daily_stack %>%
        group_by(day_label, wind_sector) %>%
        summarise(total_n   = sum(n_bin),
                  total_pct = sum(pct),
                  .groups   = "drop")
# Annotate ONLY the most prevalent sector per day, with the share as "X.X%".
daily_top <- daily_totals %>%
        group_by(day_label) %>%
        slice_max(total_pct, n = 1, with_ties = FALSE) %>%
        ungroup() %>%
        mutate(label = sprintf("%.1f%%", total_pct))

# Order panels chronologically (dd.mm.yyyy string sort wouldn't be chronological).
day_levels   <- daily_data %>% distinct(day = as.Date(DATE.TIME), day_label) %>%
        arrange(day) %>% pull(day_label)
daily_stack <- daily_stack %>% mutate(day_label = factor(day_label, levels = day_levels))
daily_top   <- daily_top   %>% mutate(day_label = factor(day_label, levels = day_levels))

# Ring labels — placed at the S compass position so they sit on each ring
# radius going outward from the centre. One copy per panel.
RING_STEPS <- c(20, 40, 60, 80)
ring_labels_daily <- expand.grid(day_label = factor(day_levels, levels = day_levels),
                                 y_value   = RING_STEPS,
                                 stringsAsFactors = FALSE) %>%
        mutate(label       = sprintf("%d%%", y_value),
               wind_sector = factor("S", levels = c("N","NE","E","SE","S","SW","W","NW")))

daily_rose <- ggplot(daily_stack, aes(x = wind_sector, y = pct, fill = ws_bin)) +
        geom_col(color = "black", linewidth = 0.2, width = 0.95) +
        geom_text(data = ring_labels_daily,
                  aes(x = wind_sector, y = y_value, label = label),
                  color = "grey30", size = 3.3, inherit.aes = FALSE) +
        geom_text(data = daily_top,
                  aes(x = wind_sector,
                      y = pmin(total_pct + 6, 95),
                      label = label),
                  color = "#1F4E79", size = 5.0,
                  fontface = "bold", inherit.aes = FALSE) +
        coord_polar(start = -pi / 8) +
        facet_wrap(~ day_label, nrow = 2) +
        scale_fill_manual(values = WS_COLOURS, drop = FALSE,
                          name = expression("Wind speed (m s"^-1*")")) +
        scale_y_continuous(limits = c(0, 100),
                           breaks = c(0, 20, 40, 60, 80, 100)) +
        labs(x = NULL, y = NULL) +
        theme_bw(base_size = 14) +
        theme(legend.position  = "bottom",
              strip.text       = element_text(face = "bold", size = 14),
              strip.background = element_rect(fill = "grey95"),
              axis.text.x      = element_text(size = 13),
              axis.text.y      = element_blank(),
              axis.ticks.y     = element_blank(),
              panel.grid.minor = element_blank(),
              legend.text      = element_text(size = 13),
              legend.title     = element_text(size = 14))
ggsave(file.path(plots_dir, "wind_rose_daily.png"), daily_rose,
       width = 18, height = 10, dpi = 300, bg = "white")

# --- (b) Seasonal wind roses for the year preceding the campaign -------------
SEASON_RANGES <- tibble::tribble(
        ~season_label,                              ~start,       ~end,
        "Summer 2024\n01.06.2024 - 31.08.2024",     "2024-06-01", "2024-08-31",
        "Autumn 2024\n01.09.2024 - 30.11.2024",     "2024-09-01", "2024-11-30",
        "Winter 2024/25\n01.12.2024 - 28.02.2025",  "2024-12-01", "2025-02-28",
        "Spring 2025\n01.03.2025 - 07.04.2025",     "2025-03-01", "2025-04-07"
)

seasonal_data <- map_dfr(seq_len(nrow(SEASON_RANGES)), function(i) {
        rng <- SEASON_RANGES[i, ]
        wind_full %>%
                filter(DATE.TIME >= as.POSIXct(paste(rng$start, "00:00:00"), tz = "UTC"),
                       DATE.TIME <= as.POSIXct(paste(rng$end,   "23:59:59"), tz = "UTC")) %>%
                mutate(season = rng$season_label)
}) %>%
        mutate(season = factor(season, levels = SEASON_RANGES$season_label))

seasonal_stack <- seasonal_data %>%
        group_by(season) %>% mutate(season_total = n()) %>% ungroup() %>%
        group_by(season, season_total, wind_sector, ws_bin) %>%
        summarise(n_bin = n(), .groups = "drop") %>%
        mutate(pct = 100 * n_bin / season_total)

seasonal_totals <- seasonal_stack %>%
        group_by(season, wind_sector) %>%
        summarise(total_n   = sum(n_bin),
                  total_pct = sum(pct),
                  .groups   = "drop")
# Annotate ONLY the most prevalent sector per season, with the share as "X.X%".
seasonal_top <- seasonal_totals %>%
        group_by(season) %>%
        slice_max(total_pct, n = 1, with_ties = FALSE) %>%
        ungroup() %>%
        mutate(label = sprintf("%.1f%%", total_pct))

ring_labels_seasonal <- expand.grid(
        season  = factor(SEASON_RANGES$season_label,
                         levels = SEASON_RANGES$season_label),
        y_value = RING_STEPS,
        stringsAsFactors = FALSE) %>%
        mutate(label       = sprintf("%d%%", y_value),
               wind_sector = factor("S", levels = c("N","NE","E","SE","S","SW","W","NW")))

seasonal_rose <- ggplot(seasonal_stack, aes(x = wind_sector, y = pct, fill = ws_bin)) +
        geom_col(color = "black", linewidth = 0.2, width = 0.95) +
        geom_text(data = ring_labels_seasonal,
                  aes(x = wind_sector, y = y_value, label = label),
                  color = "grey30", size = 3.0, inherit.aes = FALSE) +
        geom_text(data = seasonal_top,
                  aes(x = wind_sector,
                      y = pmin(total_pct + 3, 95),
                      label = label),
                  color = "#1F4E79", size = 4.6,
                  fontface = "bold", inherit.aes = FALSE) +
        coord_polar(start = -pi / 8) +
        facet_wrap(~ season, nrow = 1) +
        scale_fill_manual(values = WS_COLOURS, drop = FALSE,
                          name = expression("Wind speed (m s"^-1*")")) +
        scale_y_continuous(limits = c(0, 100),
                           breaks = c(0, 20, 40, 60, 80, 100)) +
        labs(x = NULL, y = NULL) +
        theme_bw(base_size = 14) +
        theme(legend.position  = "bottom",
              strip.text       = element_text(face = "bold", size = 14),
              strip.background = element_rect(fill = "grey95"),
              axis.text.x      = element_text(size = 13),
              axis.text.y      = element_blank(),
              axis.ticks.y     = element_blank(),
              panel.grid.minor = element_blank(),
              legend.text      = element_text(size = 13),
              legend.title     = element_text(size = 14))
ggsave(file.path(plots_dir, "wind_rose_seasonal.png"), seasonal_rose,
       width = 18, height = 6.5, dpi = 300, bg = "white")

#### 18. Headline numbers — campaign mean Q and e on retained outdoor         ####
# These are the numbers that go into the manuscript's §3.3/3.4 headlines.
headline_stats <- single_for_plot %>%
        filter(var %in% c("Q_vent","e_CH4_ghLU","e_NH3_ghLU"),
               !analyzer %in% c("baseline","HOBO","USA","RGB"),
               is.finite(value)) %>%
        group_by(var, analyzer) %>%
        summarise(n = n(), mean = mean(value, na.rm = TRUE),
                  sd = sd(value, na.rm = TRUE), .groups = "drop")
headline_consortium <- headline_stats %>%
        group_by(var) %>%
        summarise(n_analyzers       = n(),
                  n_obs             = sum(n),
                  consortium_mean   = mean(mean),
                  consortium_median = median(mean),
                  sd_across         = sd(mean),
                  se_mean           = sd_across / sqrt(n_analyzers),
                  ci_lo95           = consortium_mean - 1.96 * se_mean,
                  ci_hi95           = consortium_mean + 1.96 * se_mean,
                  ci_half_pct       = 1.96 * se_mean / consortium_mean * 100,
                  .groups = "drop")
# Also a version excluding FTIR.4 from the consortium — FTIR.4 has 2-3x higher
# SD in derived quantities than the others even after library correction, so a
# robust headline must be reported in both forms.
headline_consortium_noF4 <- headline_stats %>%
        filter(analyzer != "FTIR.4") %>%
        group_by(var) %>%
        summarise(n_analyzers       = n(),
                  n_obs             = sum(n),
                  consortium_mean   = mean(mean),
                  consortium_median = median(mean),
                  sd_across         = sd(mean),
                  se_mean           = sd_across / sqrt(n_analyzers),
                  ci_lo95           = consortium_mean - 1.96 * se_mean,
                  ci_hi95           = consortium_mean + 1.96 * se_mean,
                  ci_half_pct       = 1.96 * se_mean / consortium_mean * 100,
                  .groups = "drop")
write_excel_csv(headline_stats,           file.path(tables_dir, "headline_per_analyzer.csv"))
write_excel_csv(headline_consortium,      file.path(tables_dir, "headline_consortium.csv"))
write_excel_csv(headline_consortium_noF4, file.path(tables_dir, "headline_consortium_noFTIR4.csv"))

#### 19. Text reports — copy-paste material for the manuscript               ####
fmt_num <- function(x, d = 1) ifelse(is.finite(x), formatC(x, digits = d, format = "f"), "NA")

# --- Report 1: V10 analysis configuration ---
r1 <- init_report("01_config.txt",
                  "Ringversuche V10 — analysis configuration")
write_section(r1, "Thresholds",
              sprintf("MIN_DELTA_CO2 = %g mg m-3\nMIN_Q_VENT = %g, MAX_Q_VENT = %g m3 h-1 LU-1",
                      MIN_DELTA_CO2, MIN_Q_VENT, MAX_Q_VENT))
write_section(r1, "Baseline exclusions",
              paste(c(sprintf("Default exclude:           %s",
                              paste(BASELINE_EXCLUDE_DEFAULT, collapse=", ")),
                       sprintf("Absolute NH3 exclude:      %s",
                               paste(BASELINE_EXCLUDE_BY_VAR$NH3_mgm3, collapse=", "))),
                    collapse = "\n"))
write_section(r1, "Outdoor-line choice",
              sprintf("Retained: %s\nDropped:  %s", RETAINED_OUTDOOR, DROPPED_OUTDOOR))

# --- Report 2: animal-count summary (R2.22) ---
r2 <- init_report("02_animal_counts.txt",
                  "Campaign in-barn animal counts (R2.22)")
write_section(r2, "Counts",
              paste(sprintf("n_hours = %d", animal_summary$n_hours),
                    sprintf("mean    = %.1f cows", animal_summary$n_mean),
                    sprintf("median  = %.1f cows", animal_summary$n_median),
                    sprintf("min     = %.0f cows", animal_summary$n_min),
                    sprintf("max     = %.0f cows", animal_summary$n_max),
                    sprintf("hours with <30 cows: %d", animal_summary$n_hours_lt_30),
                    sprintf("hours with <40 cows: %d", animal_summary$n_hours_lt_40),
                    sep = "\n"))

# --- Report 3: NE vs SW outdoor-line decision (R2.24) ---
r3 <- init_report("03_outdoor_line_choice.txt",
                  "NE vs SW outdoor-line analysis (R2.24)")
write_section(r3, "Paired-test summary (campaign-wide)",
              paste(capture.output(print(
                      ne_sw_tests %>% select(gas, analyzer, n, mean_NE, mean_SW,
                                             RPD_pct, p_wilcox_holm, significant) %>%
                              as.data.frame(),
                      row.names = FALSE)),
                    collapse = "\n"))
write_section(r3, "Wind-sector contamination class for each line",
              paste(capture.output(print(as.data.frame(sector_class), row.names = FALSE)),
                    collapse = "\n"))
write_section(r3, "Time fraction each line spends in each contamination class",
              paste(capture.output(print(
                      line_time_class %>% as.data.frame(), row.names = FALSE)),
                    collapse = "\n"))
write_section(r3, "Decision",
              sprintf("Retained outdoor: %s (less time downwind of LVAT + neighbour sources)\nDropped outdoor:  %s",
                      RETAINED_OUTDOOR, DROPPED_OUTDOOR))

# --- Report 4: Q and emission headline numbers ---
r4 <- init_report("04_q_e_headlines.txt",
                  "Headline Q + emission numbers, retained outdoor only")
write_section(r4, "Consortium mean and 95% CI (all 7 analyzers)",
              paste(capture.output(print(as.data.frame(headline_consortium),
                                          row.names = FALSE)),
                    collapse = "\n"))
write_section(r4, "Consortium robust headline excluding FTIR.4 (high SD even post-library)",
              paste(capture.output(print(as.data.frame(headline_consortium_noF4),
                                          row.names = FALSE)),
                    collapse = "\n"))
write_section(r4, "Per-analyzer means",
              paste(capture.output(print(as.data.frame(headline_stats),
                                          row.names = FALSE)),
                    collapse = "\n"))

# --- Report 5: BA + pairwise summary ---
r5 <- init_report("05_pairwise_and_BA.txt",
                  "Bland-Altman + pairwise regression / Lin's CCC summary")
write_section(r5, "Bland-Altman (intra-lab pairs, retained outdoor)",
              paste(capture.output(print(as.data.frame(ba_table),
                                          row.names = FALSE)),
                    collapse = "\n"))
write_section(r5, "Pairwise regression + Lin's CCC (mean per variable group)",
              paste(capture.output(print(
                      pairwise_tbl %>%
                              group_by(var) %>%
                              summarise(mean_pearson = mean(pearson_r, na.rm=TRUE),
                                        mean_ccc     = mean(ccc, na.rm=TRUE),
                                        mean_slope   = mean(deming_slope, na.rm=TRUE),
                                        .groups = "drop") %>%
                              as.data.frame(), row.names = FALSE)),
                    collapse = "\n"))

# --- Report 6: FTIR.4_old library check ---
r6 <- init_report("06_FTIR4_library_check.txt",
                  "FTIR.4_old vs FTIR.4 library re-evaluation (R2.7)")
write_section(r6, "Per gas x location",
              paste(capture.output(print(as.data.frame(lib_compare),
                                          row.names = FALSE)),
                    collapse = "\n"))

#### 20. Console summary                                                      ####
cat("\n===== V10 run summary =====\n")
cat("Retained outdoor line: ", RETAINED_OUTDOOR, "\n", sep = "")
cat("Dropped outdoor line:  ", DROPPED_OUTDOOR,  "\n", sep = "")
cat("Tables:  ", tables_dir, "\n", sep = "")
cat("Plots:   ", plots_dir,  "\n", sep = "")
cat("Reports: ", report_dir, "\n", sep = "")
