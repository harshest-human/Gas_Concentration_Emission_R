#### Load libraries                                                             ####
library(tidyverse)   # ggplot2, dplyr, tidyr, stringr, purrr, readr, tibble
library(lubridate)   # date-time parsing
library(scales)      # axis breaks and number formatting
library(patchwork)   # wrap_plots()
library(dplyr)       # for %>% and data manipulation
# DescTools is used via DescTools::CCC() for Lin's concordance and need not be attached.

#### Shared configuration                                                       ####
# Variable labels (plotmath)
# With measurement units - used by trend plots, boxplots and Bland-Altman.
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
        "wd_trv"      = "Wind~Direction~Traverse~(degree)",
        "ws_trv"      = "Wind~Speed~Traverse~(m~s^-1)",
        "n_dairycows" = "'Number of Cows'"
)
# Without units - used by percent-error plots.
VAR_LABELS_PLAIN <- c(
        "CO2_mgm3"    = "c[CO2]",
        "CH4_mgm3"    = "c[CH4]",
        "NH3_mgm3"    = "c[NH3]",
        "r_CH4/CO2"   = "c[CH4]/c[CO2]",
        "r_NH3/CO2"   = "c[NH3]/c[CO2]",
        "delta_CO2"   = "Delta*c[CO2]",
        "delta_CH4"   = "Delta*c[CH4]",
        "delta_NH3"   = "Delta*c[NH3]",
        "Q_vent"      = "Q",
        "e_CH4_gh"    = "e[CH4]",
        "e_NH3_gh"    = "e[NH3]",
        "e_CH4_ghLU"  = "e[CH4]",
        "e_NH3_ghLU"  = "e[NH3]",
        "temp"        = "Temperature",
        "RH"          = "Relative~Humidity",
        "ws_mst"      = "Wind~Speed~Mast",
        "wd_mst"      = "Wind~Direction~Mast",
        "wd_trv"      = "Wind~Direction~Traverse",
        "ws_trv"      = "Wind~Speed~Traverse",
        "n_dairycows" = "'Number of Cows'"
)
# Location labels (plotmath) - superscript compass for the two backgrounds
LOC_LABELS <- c(
        "Indoor"     = "Indoor",
        "Outdoor_NE" = "Outdoor^NE",
        "Outdoor_SW" = "Outdoor^SW"
)
# Maps the _in/_N/_S column suffix to the renamed location labels
loc_from_suffix <- c("in" = "Indoor", "N" = "Outdoor_NE", "S" = "Outdoor_SW")

# Analyzer aesthetics. FTIR.4 = new spectral library (canonical);
# FTIR.4_old = superseded library run, drawn in red and dropped after Section 9.
ANALYZER_COLORS <- c(
        "FTIR.1" = "#1b9e77", "FTIR.2" = "#d95f02", "FTIR.3" = "#7570b3",
        "FTIR.4" = "#e7298a", "FTIR.4_old" = "#e41a1c",
        "CRDS.1" = "#66a61e", "CRDS.2" = "#e6ab02", "CRDS.3" = "#a6761d",
        "baseline" = "black"
)
ANALYZER_SHAPES <- c(
        "FTIR.1" = 0, "FTIR.2" = 1, "FTIR.3" = 2, "FTIR.4" = 5, "FTIR.4_old" = 13,
        "CRDS.1" = 15, "CRDS.2" = 19, "CRDS.3" = 17, "baseline" = 4
)
# Expand a lookup so every analyzer present in `x` has a value, filling
# gaps with `default`. Replaces the repeated *_full construction.
analyzer_aes <- function(x, lookup, default) {
        present <- unique(as.character(x))
        out     <- setNames(rep(default, length(present)), present)
        known   <- intersect(present, names(lookup))
        out[known] <- lookup[known]
        out
}

#### Data preparation functions                                                 ####
# indirect.CO2.balance(): CO2-balance ventilation rate and gas emissions
indirect.CO2.balance <- function(df) {
        # ppm -> mg/m^3 at 0 degC (273.15 K) and 1 atm
        ppm_to_mgm3 <- function(ppm, molar_mass) {
                T_K <- 273.15
                P   <- 101325
                R   <- 8.314472
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

                        Q_vent_N = ifelse(delta_CO2_N != 0, PCO2 / delta_CO2_N, NA_real_),
                        Q_vent_S = ifelse(delta_CO2_S != 0, PCO2 / delta_CO2_S, NA_real_),

                        e_NH3_gh_N = (delta_NH3_N * Q_vent_N / 1000) * n_dairycows_in,
                        e_CH4_gh_N = (delta_CH4_N * Q_vent_N / 1000) * n_dairycows_in,
                        e_NH3_gh_S = (delta_NH3_S * Q_vent_S / 1000) * n_dairycows_in,
                        e_CH4_gh_S = (delta_CH4_S * Q_vent_S / 1000) * n_dairycows_in,

                        e_NH3_ghLU_N = (e_NH3_gh_N * 500) / (n_dairycows_in * m_weight),
                        e_CH4_ghLU_N = (e_CH4_gh_N * 500) / (n_dairycows_in * m_weight),
                        e_NH3_ghLU_S = (e_NH3_gh_S * 500) / (n_dairycows_in * m_weight),
                        e_CH4_ghLU_S = (e_CH4_gh_S * 500) / (n_dairycows_in * m_weight)
                )

        # CH4/CO2 and NH3/CO2 ratios (%) per location
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

# reshaper(): wide -> long format, with a per-(time, location, var) baseline row
reshaper <- function(df) {
        meta_cols <- c("DATE.TIME", "analyzer")

        # Identify numeric measurement columns
        measure_cols <- df %>%
                select(-all_of(meta_cols)) %>%
                select(where(is.numeric)) %>%
                names()

        # Pivot to long format
        df_long <- df %>%
                pivot_longer(
                        cols = all_of(measure_cols),
                        names_to = "var",
                        values_to = "value"
                ) %>%
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
                # Map special analyzers
                mutate(analyzer = case_when(
                        var %in% c("temp", "RH")                           ~ "HOBO",
                        var %in% c("wd_mst", "ws_mst", "wd_trv", "ws_trv") ~ "USA",
                        var %in% c("n_dairycows")                          ~ "RGB",
                        TRUE                                               ~ analyzer
                ))

        # ---- Add baseline per DATE.TIME, location, var ----
        baseline_df <- df_long %>%
                group_by(DATE.TIME, location, var) %>%
                summarise(
                        value = mean(value, na.rm = TRUE),
                        day   = first(day),    # copy day
                        hour  = first(hour),   # copy hour
                        .groups = "drop"
                ) %>%
                mutate(analyzer = "baseline")

        df_long <- bind_rows(df_long, baseline_df) %>%
                mutate(analyzer = factor(analyzer,
                                         levels = c("FTIR.1","FTIR.2","FTIR.3","FTIR.4","FTIR.4_old",
                                                    "CRDS.1","CRDS.2","CRDS.3",
                                                    "HOBO","USA","RGB","baseline"))) %>%
                arrange(DATE.TIME, location, var)

        return(df_long)
}

#### Statistics helper functions                                                ####
# 8-point compass sector from wind direction (degrees)
deg_to_compass8 <- function(deg) {
        labs <- c("N", "NE", "E", "SE", "S", "SW", "W", "NW")
        factor(labs[floor(((deg %% 360) + 22.5) / 45) %% 8 + 1], levels = labs)
}

# Deming regression (orthogonal fit, equal error variance lambda = 1)
deming_fit <- function(x, y) {
        ok <- complete.cases(x, y); x <- x[ok]; y <- y[ok]
        sxx <- var(x); syy <- var(y); sxy <- cov(x, y)
        slope <- (syy - sxx + sqrt((syy - sxx)^2 + 4 * sxy^2)) / (2 * sxy)
        c(intercept = mean(y) - slope * mean(x), slope = slope)
}

# Relative Bland-Altman statistics (v7 logic, Bland & Altman 1999):
# differences expressed as % of the pair mean, with 95% CIs on bias.
ba_relative_stats <- function(x, y) {
        ok <- is.finite(x) & is.finite(y) & ((x + y) != 0)
        x <- x[ok]; y <- y[ok]; n <- length(x)
        if (n < 5L) return(NULL)
        d     <- 100 * (y - x) / ((x + y) / 2)
        bias  <- mean(d); sd_d <- sd(d)
        se_b  <- sd_d / sqrt(n)
        tibble(
                n         = n,
                bias_pct  = bias,
                bias_lo95 = bias - 1.96 * se_b,
                bias_hi95 = bias + 1.96 * se_b,
                loa_low   = bias - 1.96 * sd_d,
                loa_high  = bias + 1.96 * sd_d,
                sd_pct    = sd_d
        )
}

#### Plotting functions                                                         ####
# emitrendplot(): mean +/- SD over time/hour/day, faceted by variable x location
emitrendplot <- function(data, y = NULL, location_filter = NULL, plot_err = FALSE, x = "DATE.TIME") {
        # Filter by location and variable
        if (!is.null(location_filter)) {
                data <- data %>% filter(location %in% location_filter)
        }
        if (!is.null(y)) {
                data <- data %>% filter(var %in% y)
        }

        # Choose value column
        value_col <- if (plot_err && "pct_err" %in% names(data)) "pct_err" else "value"
        if (!value_col %in% names(data)) {
                stop("Column '", value_col, "' not found in data.")
        }

        # Compute summary (mean +/- SD)
        summary_data <- data %>%
                group_by(.data[[x]], analyzer, location, var) %>%
                summarise(
                        mean_val = mean(.data[[value_col]], na.rm = TRUE),
                        sd_val   = sd(.data[[value_col]], na.rm = TRUE),
                        .groups = "drop"
                ) %>%
                filter(!is.na(mean_val))

        if (nrow(summary_data) == 0) stop("No data available for the selected variables and plot_err setting.")

        # Facet labels (shared config)
        labels <- if (plot_err) VAR_LABELS_PLAIN else VAR_LABELS_UNITS
        summary_data <- summary_data %>%
                mutate(facet_label = factor(labels[as.character(var)], levels = labels[y]))

        # Analyzer aesthetics (shared config)
        col_vals   <- analyzer_aes(summary_data$analyzer, ANALYZER_COLORS, "black")
        shape_vals <- analyzer_aes(summary_data$analyzer, ANALYZER_SHAPES, 16)

        # Base plot
        p <- ggplot(summary_data, aes(x = .data[[x]], y = mean_val,
                                      color = analyzer, shape = analyzer, group = analyzer))

        if (x %in% c("day", "hour")) {
                p <- p +
                        geom_point(size = 2, position = position_dodge(width = 0.5)) +
                        geom_errorbar(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
                                      width = 0.4, position = position_dodge(width = 0.5))
        } else {
                p <- p +
                        geom_line(linewidth = 0.5, alpha = 0.6, na.rm = TRUE) +
                        geom_point(size = 2, na.rm = TRUE) +
                        geom_errorbar(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
                                      width = 0.4, na.rm = TRUE)
        }

        # Add styling/scales
        p <- p +
                scale_color_manual(values = col_vals) +
                scale_shape_manual(values = shape_vals) +
                facet_grid(facet_label ~ location, scales = "free_y", switch = "y",
                           labeller = labeller(
                                   facet_label = label_parsed,
                                   location = as_labeller(LOC_LABELS, label_parsed)
                           )) +
                scale_y_continuous(
                        breaks = scales::pretty_breaks(n = 6),
                        labels = scales::label_number(accuracy = 1, big.mark = "")
                ) +
                labs(x = NULL, y = NULL) +
                theme_classic() +
                theme(
                        text = element_text(size = 14),
                        axis.text = element_text(size = 14),
                        axis.title = element_text(size = 14),
                        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
                        axis.text.y = element_text(hjust = 1, size = 12),
                        strip.text.y.left = element_text(size = 14, vjust = 0.5),
                        panel.border = element_rect(color = "black", fill = NA),
                        legend.position = "bottom",
                        legend.title = element_blank(),
                        plot.title = element_text(hjust = 0.5)
                ) +
                guides(color = guide_legend(nrow = 1))

        # Time axis formatting
        if (x == "DATE.TIME") {
                x_breaks <- seq(from = min(summary_data$DATE.TIME, na.rm = TRUE),
                                to   = max(summary_data$DATE.TIME, na.rm = TRUE),
                                by   = "6 hours")
                p <- p + scale_x_datetime(breaks = x_breaks, date_labels = "%Y-%m-%d %H:%M")
        }

        return(p)
}

# emiboxplot(): distribution per analyzer, faceted by variable x location.
# Boxes are drawn as coloured outlines (no fill) over coloured jittered points.
emiboxplot <- function(data, y = NULL, location_filter = NULL, plot_err = FALSE) {
        if (!is.null(location_filter)) {
                data <- data %>% filter(location %in% location_filter)
        }
        if (!is.null(y)) {
                data <- data %>% filter(var %in% y)
        }

        value_col <- if (plot_err) "pct_err" else "value"
        if (!value_col %in% names(data)) stop("Column '", value_col, "' not found in data.")

        # Remove outliers using IQR method per location, analyzer, variable (across all days)
        data <- data %>%
                group_by(location, analyzer, var) %>%
                mutate(
                        Q1  = quantile(.data[[value_col]], 0.25, na.rm = TRUE),
                        Q3  = quantile(.data[[value_col]], 0.75, na.rm = TRUE),
                        IQR = Q3 - Q1,
                        is_outlier = (.data[[value_col]] < (Q1 - 1.5 * IQR)) |
                                     (.data[[value_col]] > (Q3 + 1.5 * IQR))
                ) %>%
                filter(!is_outlier) %>%
                ungroup() %>%
                select(-Q1, -Q3, -IQR, -is_outlier)

        # Facet labels (shared config)
        labels <- if (plot_err) VAR_LABELS_PLAIN else VAR_LABELS_UNITS
        data <- data %>%
                mutate(variable_label = factor(labels[as.character(var)], levels = labels[y]))

        col_vals <- analyzer_aes(data$analyzer, ANALYZER_COLORS, "black")

        p <- ggplot(data, aes(x = analyzer, y = .data[[value_col]], color = analyzer)) +
                geom_boxplot(outlier.shape = NA, fill = NA, linewidth = 0.6) +
                geom_jitter(width = 0.2, alpha = 0.3, size = 1.5) +
                scale_color_manual(values = col_vals) +
                facet_grid(variable_label ~ location, scales = "free_y", switch = "y",
                           labeller = labeller(
                                   variable_label = label_parsed,
                                   location = as_labeller(LOC_LABELS, label_parsed)
                           )) +
                scale_y_continuous(
                        breaks = scales::pretty_breaks(n = 5),
                        labels = scales::label_number(accuracy = 0.1, big.mark = "")
                ) +
                labs(x = NULL, y = NULL) +
                theme_classic() +
                theme(
                        text = element_text(size = 12),
                        axis.text = element_text(size = 11),
                        axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
                        strip.text = element_text(size = 12),
                        strip.text.y.left = element_text(angle = 0, hjust = 1),
                        panel.border = element_rect(color = "black", fill = NA),
                        legend.position = "bottom",
                        legend.title = element_blank(),
                        plot.title = element_text(hjust = 0.5)
                ) +
                guides(color = guide_legend(nrow = 1))

        return(p)
}

# bland_altman_plot(): RELATIVE Bland-Altman (% of pair mean) for an analyzer
# pair and one variable. Relative scale keeps the limits of agreement
# comparable across the wide concentration range (v7 module 04 logic).
bland_altman_plot <- function(data, var_filter, analyzer_pair, location_filter = NULL, x = "DATE.TIME") {
        var_label_expr <- parse(text = VAR_LABELS_UNITS[[var_filter]])[[1]]

        df <- data %>% filter(var == var_filter)
        if (!is.null(location_filter)) {
                df <- df %>% filter(location %in% location_filter)
        }
        df <- df %>% filter(analyzer %in% analyzer_pair)

        df_wide <- df %>%
                select(all_of(c(x, "location", "analyzer", "value"))) %>%
                pivot_wider(names_from = analyzer, values_from = value)

        a1 <- analyzer_pair[1]
        a2 <- analyzer_pair[2]

        df_ba <- df_wide %>%
                mutate(
                        mean_val = (.data[[a1]] + .data[[a2]]) / 2,
                        diff_pct = 100 * (.data[[a2]] - .data[[a1]]) / mean_val
                ) %>%
                filter(is.finite(diff_pct))

        bias   <- mean(df_ba$diff_pct, na.rm = TRUE)
        sd_d   <- sd(df_ba$diff_pct, na.rm = TRUE)
        loa_hi <- bias + 1.96 * sd_d
        loa_lo <- bias - 1.96 * sd_d

        subtitle_expr <- if (!is.null(location_filter)) {
                loc_expr <- parse(text = LOC_LABELS[location_filter])[[1]]
                bquote(.(var_label_expr) ~ "|" ~ .(loc_expr))
        } else {
                var_label_expr
        }

        ggplot(df_ba, aes(x = mean_val, y = diff_pct)) +
                geom_point(alpha = 0.5, size = 1) +
                geom_hline(yintercept = 0, color = "grey60") +
                geom_hline(yintercept = bias,   color = "blue", linetype = "dashed", linewidth = 0.9) +
                geom_hline(yintercept = loa_hi, color = "red",  linetype = "dotted", linewidth = 0.9) +
                geom_hline(yintercept = loa_lo, color = "red",  linetype = "dotted", linewidth = 0.9) +
                annotate("text", x = Inf, y = bias, hjust = 1.05, vjust = -0.4, size = 3,
                         label = sprintf("bias = %.1f%%", bias)) +
                scale_x_continuous(breaks = pretty_breaks(n = 6),
                                   labels = scales::label_number(accuracy = 0.1)) +
                labs(
                        subtitle = subtitle_expr,
                        x = bquote("Mean of" ~ .(a1) ~ "and" ~ .(a2)),
                        y = bquote("Relative difference (" ~ .(a2) - .(a1) ~ ") %")
                ) +
                theme_classic() +
                theme(
                        plot.subtitle = element_text(hjust = 0.5),
                        plot.title = element_text(hjust = 0.5)
                )
}

#### 1. Paths and time range                                                    ####
base_dir   <- "D:/Data_Analysis_R/Gas_Concentration_Emission_R/workflows/ringversuche_analysis"
data_dir   <- file.path(base_dir, "clean_data/Version_9/long_format")
meta_dir   <- file.path(base_dir, "meta_data")
tables_dir <- file.path(base_dir, "result_data/tables/Version_9")
plots_dir  <- file.path(base_dir, "result_data/plots/Version_9")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

start_time <- as.POSIXct("2025-04-08 12:00:00", tz = "UTC")
end_time   <- as.POSIXct("2025-04-14 12:00:00", tz = "UTC")

gases <- c("CO2", "CH4", "NH3")

#### 2. Read & merge gas datasets                                               ####
gas_files <- list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)

gas_data <- map_dfr(gas_files, read.csv, stringsAsFactors = FALSE) %>%
        select(-any_of("MPVPosition")) %>%  # CRDS files carry an extra column; drop so pivot id_cols stay clean
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

#### 3. Read animal, temperature, wind datasets                                 ####
animal_data <- read.csv(file.path(meta_dir, "RGB_Animal_count/20250408-15_LVAT_Animal_data.csv"),
                        stringsAsFactors = FALSE) %>%
        mutate(DATE.TIME = dmy_hm(DATE.TIME, tz = "UTC")) %>%
        filter(DATE.TIME >= start_time & DATE.TIME <= end_time) %>%
        select(-hour) %>%
        rename(n_dairycows_in = n_animals) %>%
        distinct()  # source CSV has each timestamp duplicated; collapse exact dups

T_RH_HOBO <- read.csv(file.path(meta_dir, "HOBO_Temp_RH/2025/20250408-20250630_HOBO_Temp_RH_1hour.csv"),
                      stringsAsFactors = FALSE) %>%
        # HOBO file mixes 4-digit (MM/DD/YYYY) and 2-digit (MM/DD/YY) year formats;
        # mdy_hms() handles both transparently.
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

#### 4. Combine all input parameters                                            ####
input_combined <- gas_data %>%
        left_join(animal_data, by = "DATE.TIME") %>%
        left_join(T_RH_HOBO,   by = "DATE.TIME") %>%
        left_join(wind_data,   by = "DATE.TIME") %>%
        arrange(DATE.TIME, lab, analyzer)

#### 5. Emissions per (lab, analyzer)                                           ####
emission_result <- indirect.CO2.balance(input_combined)
#### 6. Emissions reshaped                                                      ####
emission_reshaped <- reshaper(emission_result) %>%
        mutate(across(where(is.numeric), ~ round(.x, 2)))

#### 7. Write csv tables                                                        ####
write_excel_csv(input_combined,  file.path(tables_dir, "20250408-15_input_combined.csv"))
write_excel_csv(emission_result, file.path(tables_dir, "20250408-15_emission_result.csv"))
write_excel_csv(emission_reshaped, file.path(tables_dir, "20250408-15_ringversuche_emission_reshaped.csv"))

#### 8. Absolute concentration plots (incl. FTIR.4_old)                         ####
# Plotted FIRST, with FTIR.4_old still in, so the old vs new spectral-library
# difference is visible before it is dropped from the rest of the analysis.
c_trend_plot <- emitrendplot(emission_reshaped, y = c("CO2_mgm3", "CH4_mgm3", "NH3_mgm3"))
c_boxplot    <- emiboxplot(emission_reshaped,   y = c("CO2_mgm3", "CH4_mgm3", "NH3_mgm3"))
ggsave(file.path(plots_dir, "c_trend_plot.png"), c_trend_plot, width = 12, height = 8, dpi = 300)
ggsave(file.path(plots_dir, "c_boxplot.png"),    c_boxplot,    width = 12, height = 8, dpi = 300)

#### 9. FTIR.4 vs FTIR.4_old (spectral-library check)                           ####
# Justifies dropping FTIR.4_old: the two are the SAME instrument, FTIR.4_old
# evaluated with the superseded spectral library, FTIR.4 re-evaluated from the
# spectra with the corrected library. Fully paired (identical timestamps).
ftir_old <- input_combined %>% filter(analyzer == "FTIR.4_old")
ftir_new <- input_combined %>% filter(analyzer == "FTIR.4")

# Paired FTIR.4_old (v_old) vs FTIR.4 (v_new) values for one gas x location
lib_pair <- function(g, l) {
        col <- paste0(g, "_ppm_", l)
        inner_join(
                tibble(DATE.TIME = ftir_old$DATE.TIME, v_old = ftir_old[[col]]),
                tibble(DATE.TIME = ftir_new$DATE.TIME, v_new = ftir_new[[col]]),
                by = "DATE.TIME") %>%
                filter(complete.cases(v_old, v_new))
}

lib_compare <- map_dfr(gases, function(g) {
        map_dfr(c("in", "N", "S"), function(l) {
                j <- lib_pair(g, l)
                if (nrow(j) < 3) return(NULL)
                dem <- deming_fit(j$v_old, j$v_new)
                tibble(
                        gas       = g,
                        location  = unname(loc_from_suffix[l]),
                        n         = nrow(j),
                        mean_old  = mean(j$v_old),
                        mean_new  = mean(j$v_new),
                        mean_diff = mean(j$v_new - j$v_old),
                        RPD_pct   = 100 * mean(j$v_new - j$v_old) / mean((j$v_old + j$v_new) / 2),
                        RMSD      = sqrt(mean((j$v_new - j$v_old)^2)),
                        deming_intercept = unname(dem["intercept"]),
                        deming_slope     = unname(dem["slope"]),
                        CCC       = DescTools::CCC(j$v_old, j$v_new)$rho.c[, "est"],
                        p_wilcox  = suppressWarnings(
                                wilcox.test(j$v_old, j$v_new, paired = TRUE)$p.value)
                )
        })
})
write_excel_csv(lib_compare, file.path(tables_dir, "FTIR4_old_vs_new_comparison.csv"))

lib_points <- map_dfr(gases, function(g) {
        map_dfr(c("in", "N", "S"), function(l) {
                lib_pair(g, l) %>% mutate(gas = g, location = unname(loc_from_suffix[l]))
        })
})

lib_plot <- ggplot(lib_points, aes(x = v_old, y = v_new)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
        geom_point(alpha = 0.5, size = 1) +
        geom_abline(data = lib_compare,
                    aes(slope = deming_slope, intercept = deming_intercept),
                    color = "#e41a1c", linewidth = 0.7) +
        facet_grid(gas ~ location, scales = "free") +
        labs(title = "FTIR.4_old vs FTIR.4 after spectral-library correction",
             subtitle = "Dashed grey = 1:1 line; red = Deming regression",
             x = "FTIR.4_old (ppm)", y = "FTIR.4 (ppm)") +
        theme_bw(base_size = 12)
ggsave(file.path(plots_dir, "FTIR4_old_vs_new_scatter.png"), lib_plot,
       width = 11, height = 9, dpi = 300, bg = "white")

#### 10. Drop FTIR.4_old from all further analysis                              ####
# FTIR.4_old (old spectral library) showed large, systematic errors against
# its own re-evaluation FTIR.4 (see Section 9: CO2 ~ -6 to -7 %, CH4/NH3
# off by > 80 %). The deviation is a library artefact, not a real measurement
# difference, so FTIR.4 (corrected library) is kept and FTIR.4_old is dropped
# from every plot and statistic below. Datasets are rebuilt without it so the
# baseline mean is also recomputed cleanly.
emission_result_v <- emission_result %>% filter(analyzer != "FTIR.4_old")
emission_reshaped_v <- reshaper(emission_result_v) %>%
        mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
        mutate(analyzer = fct_drop(analyzer))
input_combined_v <- input_combined %>% filter(analyzer != "FTIR.4_old")
analyzers <- sort(unique(input_combined_v$analyzer))

#### 11. Delta, ventilation and emission plots (no FTIR.4_old)                  ####
all_plots <- list(
        d_trend_plot   = emitrendplot(emission_reshaped_v, y = c("delta_CO2", "delta_CH4", "delta_NH3")),
        d_boxplot      = emiboxplot(emission_reshaped_v,   y = c("delta_CO2", "delta_CH4", "delta_NH3")),
        q_e_trend_plot = emitrendplot(emission_reshaped_v, y = c("Q_vent", "e_CH4_ghLU", "e_NH3_ghLU")),
        q_e_boxplot    = emiboxplot(emission_reshaped_v,   y = c("Q_vent", "e_CH4_ghLU", "e_NH3_ghLU"))
)
iwalk(all_plots, function(p, nm) {
        ggsave(file.path(plots_dir, paste0(nm, ".png")), p, width = 12, height = 8, dpi = 300)
})

#### 12. Relative Bland-Altman plots (no FTIR.4_old)                            ####
# One row of relative Bland-Altman panels per lab-internal analyzer pair.
save_bland_altman <- function(analyzer_pair, tag) {
        locs <- c("Outdoor_NE", "Outdoor_SW")
        mk <- function(v, loc) {
                bland_altman_plot(emission_reshaped_v, var_filter = v,
                                  analyzer_pair = analyzer_pair, location_filter = loc) +
                        theme(plot.margin = margin(10, 10, 10, 10))
        }

        e_panels <- list()
        for (loc in locs) {
                for (v in c("e_CH4_ghLU", "e_NH3_ghLU")) {
                        e_panels <- c(e_panels, list(mk(v, loc)))
                }
        }
        ggsave(file.path(plots_dir, paste0("e_BlandAltman_", tag, ".png")),
               wrap_plots(e_panels, ncol = 4, nrow = 1),
               width = 16, height = 4, units = "in", dpi = 300)

        q_panels <- lapply(locs, function(loc) mk("Q_vent", loc))
        ggsave(file.path(plots_dir, paste0("q_BlandAltman_", tag, ".png")),
               wrap_plots(q_panels, ncol = 2, nrow = 1),
               width = 10, height = 4, units = "in", dpi = 300)
}
save_bland_altman(c("FTIR.1", "CRDS.1"), "AnalyzerA")
save_bland_altman(c("FTIR.2", "CRDS.2"), "AnalyzerB")

# Numeric relative Bland-Altman table for the same pairs
ba_pairs <- list(AnalyzerA = c("FTIR.1", "CRDS.1"),
                 AnalyzerB = c("FTIR.2", "CRDS.2"))
ba_vars  <- c("e_CH4_ghLU", "e_NH3_ghLU", "Q_vent")
ba_table <- map_dfr(names(ba_pairs), function(tag) {
        pr <- ba_pairs[[tag]]
        map_dfr(ba_vars, function(v) {
                map_dfr(c("Outdoor_NE", "Outdoor_SW"), function(loc) {
                        w <- emission_reshaped_v %>%
                                filter(var == v, location == loc, analyzer %in% pr) %>%
                                select(DATE.TIME, analyzer, value) %>%
                                pivot_wider(names_from = analyzer, values_from = value)
                        st <- ba_relative_stats(w[[pr[1]]], w[[pr[2]]])
                        if (is.null(st)) return(NULL)
                        bind_cols(tibble(pair = tag, variable = v, location = loc), st)
                })
        })
})
write_excel_csv(ba_table, file.path(tables_dir, "bland_altman_relative.csv"))

#### 13. Pairwise method comparison: regression + Lin's CCC (no FTIR.4_old)     ####
# For every analyzer pair, regress and score agreement on absolute
# concentrations. Lin's CCC captures correlation AND bias in one number,
# so it replaces the old Pearson correlogram.
pairwise_compare <- function(df_long, var_sel, loc_sel) {
        w <- df_long %>%
                filter(var == var_sel, location == loc_sel,
                       !analyzer %in% c("baseline", "HOBO", "USA", "RGB")) %>%
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
                tibble(
                        var = var_sel, location = loc_sel,
                        analyzer_x = pr[1], analyzer_y = pr[2], n = length(x),
                        pearson_r        = cor(x, y),
                        ccc              = DescTools::CCC(x, y)$rho.c[, "est"],
                        deming_slope     = unname(dem["slope"]),
                        deming_intercept = unname(dem["intercept"])
                )
        })
}

conc_vars <- c("CO2_mgm3", "CH4_mgm3", "NH3_mgm3")
conc_locs <- c("Indoor", "Outdoor_NE", "Outdoor_SW")
pairwise_tbl <- map_dfr(conc_vars, function(v) {
        map_dfr(conc_locs, function(l) pairwise_compare(emission_reshaped_v, v, l))
})
write_excel_csv(pairwise_tbl, file.path(tables_dir, "pairwise_regression_ccc.csv"))

# CCC heatmap (analyzer x analyzer), faceted by gas x location - the CCC
# equivalent of the dropped correlogram.
ccc_plot <- pairwise_tbl %>%
        mutate(facet_label = factor(VAR_LABELS_PLAIN[var], levels = VAR_LABELS_PLAIN[conc_vars])) %>%
        ggplot(aes(x = analyzer_x, y = analyzer_y, fill = ccc)) +
        geom_tile(color = "white") +
        geom_text(aes(label = sprintf("%.2f", ccc)), size = 2.6) +
        scale_fill_gradient2(low = "#b2182b", mid = "#f7f7f7", high = "#2166ac",
                             midpoint = 0.5, limits = c(0, 1), name = "Lin's CCC") +
        facet_grid(location ~ facet_label,
                   labeller = labeller(facet_label = label_parsed,
                                       location = as_labeller(LOC_LABELS, label_parsed))) +
        labs(title = "Pairwise agreement (Lin's concordance correlation)",
             x = NULL, y = NULL) +
        theme_bw(base_size = 11) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "bottom")
ggsave(file.path(plots_dir, "pairwise_ccc_heatmap.png"), ccc_plot,
       width = 12, height = 9, dpi = 300, bg = "white")

#### 14. Wind-sector and wind-speed analysis (no FTIR.4_old)                    ####
# Long format: analyzer x timestamp x gas x location, with wind sector + speed.
ws_long <- input_combined_v %>%
        filter(!is.na(wd_mst), !is.na(ws_mst)) %>%
        mutate(wind_sector = deg_to_compass8(wd_mst),
               speed_class = cut_number(ws_mst, 4)) %>%
        select(DATE.TIME, analyzer, wind_sector, speed_class,
               matches("^(CO2|CH4|NH3)_ppm_(in|N|S)$")) %>%
        pivot_longer(cols = matches("_ppm_"),
                     names_to = c("gas", "loc"), names_pattern = "(.+)_ppm_(.+)",
                     values_to = "ppm") %>%
        mutate(location = unname(loc_from_suffix[loc]),
               gas      = factor(gas, levels = gases))

# Campaign-average concentration per analyzer x gas x location x wind sector
ws_summary <- ws_long %>%
        group_by(analyzer, gas, location, wind_sector) %>%
        summarise(mean_ppm = mean(ppm, na.rm = TRUE), n = n(), .groups = "drop")
write_excel_csv(ws_summary, file.path(tables_dir, "wind_sector_concentration.csv"))

# Wind rose - campaign wind-sector frequency (one count per timestamp)
wind_rose <- input_combined_v %>%
        distinct(DATE.TIME, wd_mst) %>%
        filter(!is.na(wd_mst)) %>%
        mutate(wind_sector = deg_to_compass8(wd_mst)) %>%
        count(wind_sector, .drop = FALSE) %>%
        ggplot(aes(x = wind_sector, y = n)) +
        geom_col(fill = "#377EB8", color = "black", width = 1) +
        geom_text(aes(label = n), vjust = -0.3, size = 3.5) +
        coord_polar(start = -pi / 8) +
        labs(title = "Campaign wind-sector frequency", x = NULL, y = "Timestamps") +
        theme_minimal(base_size = 13)
ggsave(file.path(plots_dir, "wind_rose.png"), wind_rose,
       width = 6, height = 6, dpi = 300, bg = "white")

# Mean concentration by wind sector, faceted gas (rows) x analyzer (cols)
ws_plot <- ggplot(ws_summary, aes(x = wind_sector, y = mean_ppm, fill = location)) +
        geom_col(position = position_dodge(width = 0.8), width = 0.7,
                 color = "black", linewidth = 0.2) +
        facet_grid(gas ~ analyzer, scales = "free_y") +
        scale_fill_manual(values = c("Indoor"     = "#4DAF4A",
                                     "Outdoor_NE" = "#377EB8",
                                     "Outdoor_SW" = "#E41A1C")) +
        labs(title = "Average gas concentration by wind sector and analyzer",
             x = "Wind sector", y = "Mean concentration (ppm)") +
        theme_bw(base_size = 12) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
              legend.position = "bottom", legend.title = element_blank())
ggsave(file.path(plots_dir, "wind_sector_concentration.png"), ws_plot,
       width = 22, height = 9, dpi = 300, bg = "white")

# Bivariate polar plot: wind direction (angle) x wind-speed class (radius),
# fill = mean concentration. Reveals whether highs are direction-driven
# (advected source) or low-speed-driven (local accumulation). One plot per
# gas (own colour scale), faceted by location, averaged over analyzers.
ws_speed_summary <- ws_long %>%
        group_by(gas, location, wind_sector, speed_class) %>%
        summarise(mean_ppm = mean(ppm, na.rm = TRUE), n = n(), .groups = "drop")
write_excel_csv(ws_speed_summary, file.path(tables_dir, "wind_sector_speed_concentration.csv"))

make_polar <- function(g) {
        ggplot(filter(ws_speed_summary, gas == g),
               aes(x = wind_sector, y = speed_class, fill = mean_ppm)) +
                geom_tile(color = "white") +
                coord_polar(start = -pi / 8) +
                facet_wrap(~ location) +
                scale_fill_viridis_c(name = paste0(g, " (ppm)"), option = "C") +
                labs(title = g, x = NULL, y = "Wind-speed class (m/s)") +
                theme_minimal(base_size = 11) +
                theme(axis.text.x = element_text(size = 8),
                      legend.position = "right")
}
polar_fig <- wrap_plots(lapply(gases, make_polar), ncol = 1)
ggsave(file.path(plots_dir, "wind_direction_speed_polar.png"), polar_fig,
       width = 12, height = 15, dpi = 300, bg = "white")

#### 15. Outdoor_NE vs Outdoor_SW comparison (no FTIR.4_old)                    ####
# Paired comparison (same analyzer & timestamp) for each gas x analyzer.
ne_sw_tests <- map_dfr(gases, function(g) {
        map_dfr(as.character(analyzers), function(a) {
                d  <- input_combined_v %>% filter(analyzer == a)
                ne <- d[[paste0(g, "_ppm_N")]]
                sw <- d[[paste0(g, "_ppm_S")]]
                ok <- complete.cases(ne, sw)
                ne <- ne[ok]; sw <- sw[ok]
                if (length(ne) < 3) return(NULL)
                tibble(
                        gas       = g,
                        analyzer  = a,
                        n         = length(ne),
                        mean_NE   = mean(ne),
                        mean_SW   = mean(sw),
                        mean_diff = mean(ne - sw),
                        RPD_pct   = 100 * mean(ne - sw) / mean((ne + sw) / 2),
                        p_ttest   = t.test(ne, sw, paired = TRUE)$p.value,
                        p_wilcox  = suppressWarnings(
                                wilcox.test(ne, sw, paired = TRUE)$p.value)
                )
        })
}) %>%
        mutate(p_wilcox_holm = p.adjust(p_wilcox, method = "holm"),
               significant   = ifelse(p_wilcox_holm < 0.05, "yes", "no"))
write_excel_csv(ne_sw_tests, file.path(tables_dir, "outdoor_NE_vs_SW_tests.csv"))

ne_sw_plot <- ggplot(ne_sw_tests, aes(x = analyzer, y = RPD_pct, fill = significant)) +
        geom_col(color = "black", linewidth = 0.2) +
        geom_hline(yintercept = 0) +
        facet_wrap(~ gas, ncol = 1, scales = "free_y") +
        scale_fill_manual(values = c(yes = "#E41A1C", no = "grey70")) +
        labs(title = "Outdoor NE vs SW: relative percentage difference",
             subtitle = "RPD = 100 x mean(NE - SW) / mean concentration; red = Holm-adjusted Wilcoxon p < 0.05",
             x = NULL, y = "RPD (%)") +
        theme_bw(base_size = 12) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "bottom")
ggsave(file.path(plots_dir, "outdoor_NE_vs_SW_RPD.png"), ne_sw_plot,
       width = 9, height = 9, dpi = 300, bg = "white")

#### 16. Background-choice summary                                              ####
bg_summary <- ne_sw_tests %>%
        group_by(gas) %>%
        summarise(
                analyzers_tested = n(),
                analyzers_signif = sum(significant == "yes"),
                median_RPD_pct   = median(RPD_pct),
                mean_NE          = mean(mean_NE),
                mean_SW          = mean(mean_SW),
                .groups = "drop"
        ) %>%
        mutate(lower_background = ifelse(mean_NE < mean_SW, "Outdoor_NE", "Outdoor_SW"))
write_excel_csv(bg_summary, file.path(tables_dir, "background_choice_summary.csv"))
print(bg_summary)
