#### Load libraries                             ####
library(tidyverse)   # ggplot2, dplyr, tidyr, stringr, purrr, readr, tibble
library(lubridate)   # date-time parsing
library(scales)      # axis breaks and number formatting
library(patchwork)   # wrap_plots()
# Hmisc is used via Hmisc::rcorr() inside emicorrgram() and need not be attached.

#### Shared configuration                       ####
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
# Without units - used by correlograms and percent-error plots.
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
# Analyzer aesthetics 
ANALYZER_COLORS <- c(
        "FTIR.1" = "#1b9e77", "FTIR.2" = "#d95f02", "FTIR.3" = "#7570b3",
        "FTIR.4" = "#e7298a", "FTIR.4_v2" = "#c71c6d",
        "CRDS.1" = "#66a61e", "CRDS.2" = "#e6ab02", "CRDS.3" = "#a6761d",
        "baseline" = "black"
)
ANALYZER_SHAPES <- c(
        "FTIR.1" = 0, "FTIR.2" = 1, "FTIR.3" = 2, "FTIR.4" = 5, "FTIR.4_v2" = 5,
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

#### Data preparation functions                 ####
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
                                str_detect(var, "_in$") ~ "Barn inside",
                                str_detect(var, "_N$")  ~ "North background",
                                str_detect(var, "_S$")  ~ "South background",
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
                                         levels = c("FTIR.1","FTIR.2","FTIR.3","FTIR.4","FTIR.4_v2",
                                                    "CRDS.1","CRDS.2","CRDS.3",
                                                    "HOBO","USA","RGB","baseline"))) %>%
                arrange(DATE.TIME, location, var)

        return(df_long)
}

#### Plotting functions                         ####
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
                                   location = label_value
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

# emicorrgram(): pairwise Pearson correlation between analyzers, per variable x location
emicorrgram <- function(data, target_variables, locations = NULL) {
        # Helper: correlation matrix for one variable & location
        compute_corr <- function(var_sel, loc_sel) {
                filtered <- data %>%
                        filter(var == var_sel) %>%
                        select(DATE.TIME, location, analyzer, value)

                if (!is.null(loc_sel)) {
                        filtered <- filtered %>% filter(location %in% loc_sel)
                }

                filtered <- filtered %>%
                        group_by(DATE.TIME, location, analyzer) %>%
                        summarise(value = mean(value, na.rm = TRUE), .groups = "drop")

                pivoted <- filtered %>%
                        pivot_wider(names_from = analyzer, values_from = value,
                                    values_fn = function(x, ...) mean(x, na.rm = TRUE))

                subdata_num <- pivoted %>%
                        select(-DATE.TIME, -location) %>%
                        mutate(across(everything(), as.numeric))

                if (ncol(subdata_num) < 2) return(NULL)

                cor_res <- Hmisc::rcorr(as.matrix(subdata_num), type = "pearson")

                expand.grid(
                        Var1 = colnames(cor_res$r),
                        Var2 = colnames(cor_res$r),
                        stringsAsFactors = FALSE
                ) %>%
                        mutate(
                                correlation = as.vector(cor_res$r),
                                pvalue = as.vector(cor_res$P),
                                target_variable = var_sel,
                                location = loc_sel
                        ) %>%
                        filter(as.numeric(factor(Var1)) > as.numeric(factor(Var2)))
        }

        # Loop over variables and locations
        corr_df <- map_df(target_variables, function(tv) {
                map_df(locations, function(loc) compute_corr(tv, loc))
        }) %>%
                mutate(facet_label = factor(VAR_LABELS_PLAIN[target_variable],
                                            levels = VAR_LABELS_PLAIN[target_variables]))

        # Plot
        p <- ggplot(corr_df, aes(x = Var1, y = Var2, fill = correlation)) +
                geom_tile(color = "white") +
                geom_text(aes(label = paste0(
                        round(correlation, 2),
                        ifelse(pvalue <= 0.001, "\n***",
                               ifelse(pvalue <= 0.01, "\n**",
                                      ifelse(pvalue <= 0.05, "\n*", "\nns")))
                )), size = 3, color = "white") +
                scale_y_discrete(position = "right") +
                scale_fill_gradientn(
                        colors = c("darkred","red","orange2","gold2","yellow",
                                   "greenyellow","green1","green3","darkgreen"),
                        limits = c(0,1),
                        name = "PCC"
                ) +
                facet_grid(location ~ facet_label, scales = "free", space = "free", switch = "y",
                           labeller = labeller(facet_label = label_parsed)) +
                theme_classic() +
                theme(
                        axis.title = element_blank(),
                        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                        axis.text.y = element_text(angle = 45, hjust = 1, size = 10),
                        legend.position = "bottom",
                        strip.text = element_text(size = 12),
                        panel.border = element_rect(colour = "black", fill = NA)
                )

        return(p)
}

# emiboxplot(): distribution per analyzer, faceted by variable x location (IQR outliers removed)
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

        fill_vals <- analyzer_aes(data$analyzer, ANALYZER_COLORS, "black")

        p <- ggplot(data, aes(x = analyzer, y = .data[[value_col]], fill = analyzer)) +
                geom_boxplot(outlier.shape = NA, alpha = 0.8) +
                geom_jitter(width = 0.2, alpha = 0.3, size = 1.5) +
                scale_fill_manual(values = fill_vals) +
                facet_grid(variable_label ~ location, scales = "free_y", switch = "y",
                           labeller = labeller(
                                   variable_label = label_parsed,
                                   location = label_value
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
                guides(fill = guide_legend(nrow = 1))

        return(p)
}

# bland_altman_plot(): agreement between a pair of analyzers for one variable
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
                        mean_val = ( .data[[a1]] + .data[[a2]] ) / 2,
                        diff_val = .data[[a1]] - .data[[a2]]
                )

        bias      <- mean(df_ba$diff_val, na.rm = TRUE)
        sd_diff   <- sd(df_ba$diff_val, na.rm = TRUE)
        loa_upper <- bias + 1.96 * sd_diff
        loa_lower <- bias - 1.96 * sd_diff

        subtitle_expr <- if (!is.null(location_filter)) {
                bquote(.(var_label_expr) ~ "|" ~ .(location_filter))
        } else {
                var_label_expr
        }

        p <- ggplot(df_ba, aes(x = mean_val, y = diff_val)) +
                geom_point(alpha = 0.6) +
                geom_hline(yintercept = bias,      color = "blue", linetype = "dashed", linewidth = 1) +
                geom_hline(yintercept = loa_upper, color = "red",  linetype = "dotted", linewidth = 1) +
                geom_hline(yintercept = loa_lower, color = "red",  linetype = "dotted", linewidth = 1) +
                scale_x_continuous(
                        breaks = pretty_breaks(n = 8),
                        labels = scales::label_number(accuracy = 0.1)
                ) +
                scale_y_continuous(
                        breaks = pretty_breaks(n = 8),
                        labels = scales::label_number(accuracy = 0.1)
                ) +
                labs(
                        subtitle = subtitle_expr,
                        x = bquote("Mean of" ~ .(a1) ~ "and" ~ .(a2)),
                        y = bquote("Difference (" ~ .(a1) - .(a2) ~ ")")
                ) +
                theme_classic() +
                theme(
                        plot.subtitle = element_text(hjust = 0.5),
                        plot.title = element_text(hjust = 0.5)
                )

        return(p)
}

#### 1. Paths and time range                    ####
base_dir   <- "D:/Data_Analysis_R/Gas_Concentration_Emission_R/workflows/ringversuche_analysis"
data_dir   <- file.path(base_dir, "clean_data/Version_6/long_format")
meta_dir   <- file.path(base_dir, "meta_data")
tables_dir <- file.path(base_dir, "result_data/tables/Version_6")
plots_dir  <- file.path(base_dir, "result_data/plots/Version_6")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

start_time <- as.POSIXct("2025-04-08 12:00:00", tz = "UTC")
end_time   <- as.POSIXct("2025-04-14 12:00:00", tz = "UTC")

#### 2. Read & merge gas datasets###
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

#### 3. Read animal, temperature, wind datasets ####
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

#### 4. Combine all input paramters             ####
input_combined <- gas_data %>%
        left_join(animal_data, by = "DATE.TIME") %>%
        left_join(T_RH_HOBO,   by = "DATE.TIME") %>%
        left_join(wind_data,   by = "DATE.TIME") %>%
        arrange(DATE.TIME, lab, analyzer)

#### 5. Emissions per (lab, analyzer)           ####
emission_result <- indirect.CO2.balance(input_combined)
#### 6. Emissions reshaped                      ####
emission_reshaped <-  reshaper(emission_result) %>%
        mutate(across(where(is.numeric), ~ round(.x, 2)))

#### 7. Write csv tables                        ####
write_excel_csv(input_combined,  file.path(tables_dir, "20250408-15_input_combined.csv"))
write_excel_csv(emission_result, file.path(tables_dir, "20250408-15_emission_result.csv"))
write_excel_csv(emission_reshaped, file.path(tables_dir, "20250408-15_ringversuche_emission_reshaped.csv"))

#### 8. Boxplots and time series                ####
# Concentration, delta and ventilation/emission plots - all saved at a common 12 x 8 in
all_plots <- list(
        c_trend_plot   = emitrendplot(emission_reshaped, y = c("CO2_mgm3", "CH4_mgm3", "NH3_mgm3")),
        c_boxplot      = emiboxplot(emission_reshaped,   y = c("CO2_mgm3", "CH4_mgm3", "NH3_mgm3")),
        d_trend_plot   = emitrendplot(emission_reshaped, y = c("delta_CO2", "delta_CH4", "delta_NH3")),
        d_boxplot      = emiboxplot(emission_reshaped,   y = c("delta_CO2", "delta_CH4", "delta_NH3")),
        q_e_trend_plot = emitrendplot(emission_reshaped, y = c("Q_vent", "e_CH4_ghLU", "e_NH3_ghLU")),
        q_e_boxplot    = emiboxplot(emission_reshaped,   y = c("Q_vent", "e_CH4_ghLU", "e_NH3_ghLU"))
)
iwalk(all_plots, function(p, nm) {
        ggsave(file.path(plots_dir, paste0(nm, ".png")), p, width = 12, height = 8, dpi = 300)
})

#### 9. Bland-Altman plots                      ####
# Save one row of Bland-Altman panels per analyzer pair:
save_bland_altman <- function(analyzer_pair, tag) {
        locs <- c("North background", "South background")
        mk <- function(v, loc) {
                bland_altman_plot(emission_reshaped, var_filter = v,
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

#### 10. Correlograms                           ####
d_corrgram <- emicorrgram(emission_reshaped,
                          target_variables = c("delta_CO2", "delta_CH4", "delta_NH3"),
                          locations = c("North background", "South background"))

q_e_corrgram <- emicorrgram(emission_reshaped,
                            target_variables = c("e_CH4_ghLU", "e_NH3_ghLU", "Q_vent"),
                            locations = c("North background", "South background"))

ggsave(file.path(plots_dir, "d_corrgram.png"), plot = d_corrgram,
       width = 10, height = 6, units = "in", dpi = 300)

ggsave(file.path(plots_dir, "q_e_corrgram.png"), plot = q_e_corrgram,
       width = 10, height = 6, units = "in", dpi = 300)
