
# Development of indirect.CO2.balance function
indirect.CO2.balance <- function(df) {
        library(dplyr)
        
        ppm_to_mgm3 <- function(ppm, molar_mass) {
                T_K <- 273.15      # Kelvin
                P   <- 101325      # Pa
                R   <- 8.314472    # J/mol/K
                (ppm * 1e-6) * molar_mass * 1e3 * P / (R * T_K) # mg/m³
        }
        
        df %>%
                mutate(
                        hour = as.numeric(format(DATE.TIME, "%H")),
                        a = 0.22,
                        h_min = 2.9, # hour of minimal activity
                        
                        phi = 5.6 * m_weight^0.75 + 22 * Y1_milk_prod + 1.6e-5 * p_pregnancy_day^3,
                        t_factor = 1 + 4e-5 * (20 - temp_in)^3,
                        phi_T_cor = phi * t_factor,
                        A_cor = 1 - a * sin((2 * pi / 24) * (hour + 6 - h_min)),
                        hpu_T_A_cor_per_cow = phi_T_cor * A_cor,
                        
                        PCO2 = (0.185 * hpu_T_A_cor_per_cow) * 1000,
                        
                        delta_CO2_mgm3 = ppm_to_mgm3(delta_CO2, 44.01),
                        delta_NH3_mgm3 = ppm_to_mgm3(delta_NH3, 17.031),
                        delta_CH4_mgm3 = ppm_to_mgm3(delta_CH4, 16.04),
                        
                        Q_vent = PCO2 / delta_CO2_mgm3,
                        
                        # Instantaneous emissions (g/h) divided by 1000 to convert mg to g
                        e_NH3_gh = (delta_NH3_mgm3 * Q_vent / 1000) * n_dairycows,
                        e_CH4_gh = (delta_CH4_mgm3 * Q_vent / 1000) * n_dairycows,
                        
                        # Annual emissions (kg/year) divided by 1000 to convert g to kg
                        e_NH3_ghLU = (e_NH3_gh * 500) / (n_dairycows * m_weight),
                        e_CH4_ghLU = (e_CH4_gh * 500) / (n_dairycows * m_weight)
                )
}

# Standard-parameter CO2 balance for multi-campaign comparisons.
#
# The function deliberately does not apply an indoor-temperature correction or
# a diurnal animal-activity correction to animal CO2 production. DWD
# temperature is used only in the ideal-gas conversion from ppm to mg m-3.
# The default animal constants reproduce the agreed full-occupancy scenario:
# 58 cows, 500 kg cow-1, median milk production of 31.21864407 kg cow-1 d-1,
# and a constant mean pregnancy day of 119.5510204 d.
indirect.CO2.balance.st <- function(
        df,
        inside_suffix = "in",
        outside_suffix = "out",
        dwd_temp_col = "temp_dwd",
        pressure_col = NULL,
        n_dairy_cows = 58,
        cow_weight_kg = 500,
        milk_kg_cow_d = 31.21864407,
        pregnancy_day = 119.5510204,
        pco2_per_hpu_m3_h = 0.185,
        assumed_pressure_pa = 101325,
        min_delta_co2_ppm = 0,
        low_delta_co2_ppm = 200,
        annualisation_hours = 8760
) {
        # Keep the function self-contained when this file is sourced into a
        # clean R session; do not require callers to attach dplyr first.
        if (!requireNamespace("dplyr", quietly = TRUE)) {
                stop("Package 'dplyr' is required by indirect.CO2.balance.st().")
        }
        `%>%` <- dplyr::`%>%`

        required_gas_cols <- unlist(lapply(
                c("CO2", "CH4", "NH3"),
                function(gas) paste0(gas, "_ppm_", c(inside_suffix, outside_suffix))
        ))
        required_cols <- c(required_gas_cols, dwd_temp_col)
        if (!is.null(pressure_col)) required_cols <- c(required_cols, pressure_col)

        missing_cols <- setdiff(required_cols, names(df))
        if (length(missing_cols) > 0) {
                stop(
                        "indirect.CO2.balance.st(): missing required column(s): ",
                        paste(missing_cols, collapse = ", ")
                )
        }
        if (!is.numeric(n_dairy_cows) || length(n_dairy_cows) != 1L ||
            !is.finite(n_dairy_cows) || n_dairy_cows <= 0) {
                stop("n_dairy_cows must be one positive finite number.")
        }
        if (!is.numeric(cow_weight_kg) || length(cow_weight_kg) != 1L ||
            !is.finite(cow_weight_kg) || cow_weight_kg <= 0) {
                stop("cow_weight_kg must be one positive finite number.")
        }
        if (!is.numeric(annualisation_hours) || length(annualisation_hours) != 1L ||
            !is.finite(annualisation_hours) || annualisation_hours <= 0) {
                stop("annualisation_hours must be one positive finite number.")
        }

        gas_col <- function(gas, suffix) paste0(gas, "_ppm_", suffix)
        mass_col <- function(gas, suffix) paste0(gas, "_mgm3_", suffix)

        # Constants used by the ideal-gas conversion.
        gas_constant <- 8.314462618 # Pa m3 mol-1 K-1
        molar_mass <- c(CO2 = 44.0095, CH4 = 16.04246, NH3 = 17.03052)

        # Uncorrected animal heat production at the agreed standard inputs.
        phi_w_per_cow <- 5.6 * cow_weight_kg^0.75 +
                22 * milk_kg_cow_d +
                1.6e-5 * pregnancy_day^3
        hpu_per_cow <- phi_w_per_cow / 1000
        pco2_cow_m3_h <- pco2_per_hpu_m3_h * hpu_per_cow
        pco2_barn_m3_h <- pco2_cow_m3_h * n_dairy_cows
        livestock_units <- n_dairy_cows * cow_weight_kg / 500

        result <- df %>%
                dplyr::mutate(
                        n_dairycows_std = n_dairy_cows,
                        m_weight_std_kg = cow_weight_kg,
                        Y1_milk_prod_std = milk_kg_cow_d,
                        p_pregnancy_day_std = pregnancy_day,
                        LU_std = livestock_units,
                        phi_std_W_cow = phi_w_per_cow,
                        t_factor_std = 1,
                        A_cor_std = 1,
                        hpu_std_cow = hpu_per_cow,
                        PCO2_std_m3_h_cow = pco2_cow_m3_h,
                        PCO2_std_m3_h_barn = pco2_barn_m3_h,
                        gas_temperature_K = .data[[dwd_temp_col]] + 273.15,
                        gas_pressure_Pa = if (is.null(pressure_col)) {
                                assumed_pressure_pa
                        } else {
                                .data[[pressure_col]]
                        },
                        pressure_assumed = is.null(pressure_col),
                        gas_conversion_valid = is.finite(gas_temperature_K) &
                                gas_temperature_K > 0 &
                                is.finite(gas_pressure_Pa) &
                                gas_pressure_Pa > 0
                )

        for (gas in names(molar_mass)) {
                for (suffix in c(inside_suffix, outside_suffix)) {
                        ppm_name <- gas_col(gas, suffix)
                        mgm3_name <- mass_col(gas, suffix)
                        result <- result %>%
                                dplyr::mutate(
                                        !!mgm3_name := dplyr::if_else(
                                                gas_conversion_valid,
                                                .data[[ppm_name]] * 1e-6 *
                                                        gas_pressure_Pa * molar_mass[[gas]] * 1e3 /
                                                        (gas_constant * gas_temperature_K),
                                                NA_real_
                                        )
                                )
                }
        }

        co2_in_ppm <- gas_col("CO2", inside_suffix)
        co2_out_ppm <- gas_col("CO2", outside_suffix)
        co2_in_mgm3 <- mass_col("CO2", inside_suffix)
        co2_out_mgm3 <- mass_col("CO2", outside_suffix)

        result <- result %>%
                dplyr::mutate(
                        delta_CO2_ppm = .data[[co2_in_ppm]] - .data[[co2_out_ppm]],
                        delta_CO2_mgm3 = .data[[co2_in_mgm3]] - .data[[co2_out_mgm3]],
                        delta_co2_positive = is.finite(delta_CO2_ppm) & delta_CO2_ppm > 0,
                        delta_co2_lt_200ppm = delta_co2_positive &
                                delta_CO2_ppm < low_delta_co2_ppm,
                        Q_vent_m3_h_barn = dplyr::if_else(
                                delta_co2_positive & delta_CO2_ppm >= min_delta_co2_ppm,
                                PCO2_std_m3_h_barn * 1e6 / delta_CO2_ppm,
                                NA_real_
                        ),
                        Q_vent_m3_h_cow = Q_vent_m3_h_barn / n_dairy_cows,
                        Q_vent_m3_h_LU = Q_vent_m3_h_barn / LU_std
                )

        for (gas in c("CH4", "NH3")) {
                gas_in_mgm3 <- mass_col(gas, inside_suffix)
                gas_out_mgm3 <- mass_col(gas, outside_suffix)
                delta_name <- paste0("delta_", gas, "_mgm3")
                emission_name <- paste0("e_", gas, "_gh_barn")
                emission_cow_name <- paste0("e_", gas, "_gh_cow")
                emission_lu_name <- paste0("e_", gas, "_ghLU")
                annual_name <- paste0("e_", gas, "_kg_year_LU_rate_equivalent")
                annual_short_name <- paste0("e_", gas, "_kg_year_LU")
                negative_name <- paste0("delta_", tolower(gas), "_negative")

                result <- result %>%
                        dplyr::mutate(
                                !!delta_name := .data[[gas_in_mgm3]] - .data[[gas_out_mgm3]],
                                !!negative_name := is.finite(.data[[delta_name]]) &
                                        .data[[delta_name]] < 0,
                                !!emission_name := .data[[delta_name]] *
                                        Q_vent_m3_h_barn / 1000,
                                !!emission_cow_name := .data[[emission_name]] / n_dairy_cows,
                                !!emission_lu_name := .data[[emission_name]] / LU_std,
                                # Annualised rate equivalent. This assumes the
                                # instantaneous g h-1 LU-1 rate persists for
                                # annualisation_hours; it is not a measured
                                # annual integral unless representative periods
                                # have first been weighted appropriately.
                                !!annual_name := .data[[emission_lu_name]] *
                                        annualisation_hours / 1000,
                                !!annual_short_name := .data[[annual_name]]
                        )
        }

        for (gas in c("CH4", "NH3")) {
                for (suffix in c(inside_suffix, outside_suffix)) {
                        ratio_name <- paste0("r_", gas, "_CO2_", suffix, "_pct")
                        gas_mgm3 <- mass_col(gas, suffix)
                        co2_mgm3 <- mass_col("CO2", suffix)
                        result <- result %>%
                                dplyr::mutate(
                                        !!ratio_name := dplyr::if_else(
                                                is.finite(.data[[co2_mgm3]]) &
                                                        .data[[co2_mgm3]] != 0,
                                                .data[[gas_mgm3]] / .data[[co2_mgm3]] * 100,
                                                NA_real_
                                        )
                                )
                }
        }

        result
}

# Backward-compatible name retained for any existing scripts that already use
# the earlier `.std` spelling.
indirect.CO2.balance.std <- indirect.CO2.balance.st

# Development of pivot longer function
reshaper <- function(df) {
        library(dplyr)
        library(tidyr)
        library(stringr)
        library(lubridate)
        
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
                select(DATE.TIME, analyzer, var, value) %>%
                arrange(DATE.TIME, var, analyzer) %>%
                arrange(DATE.TIME)
        
        return(df_long)
}

# Errorbarplot Function
emierrorbarplot <- function(data, y = NULL) {
        library(dplyr)
        library(ggplot2)
        library(scales)

        analyzer_labels <- c(
                "crds" = "CRDS",
                "logas_ndir" = "PRONOVA",
                "logas_tdlas" = "CUBIC",
                "otice" = "OTICE"
        )
        
        if ("var" %in% names(data) && !"variable" %in% names(data)) {
                data <- data %>% rename(variable = var)
        }
        
        if (!is.null(y)) {
                data <- data %>% filter(variable %in% y)
        }
        
        if (!"value" %in% names(data)) stop("Column 'value' not found in data.")
        
        facet_labels <- c(
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
        
        data <- data %>%
                mutate(
                        analyzer = as.character(analyzer),
                        analyzer = recode(analyzer, !!!analyzer_labels, .default = analyzer),
                        facet_label = factor(
                                facet_labels[as.character(variable)],
                                levels = facet_labels[y]
                        )
                )

        analyzer_colors <- c(
                "CUBIC" = "#1b9e77",
                "PRONOVA" = "#7570b3",
                "CRDS" = "darkgray",
                "OTICE" = "#C45A11"
        )

        analyzer_shapes <- c(
                "CUBIC" = 0,
                "PRONOVA" = 1,
                "CRDS" = 2,
                "OTICE" = 15
        )
        
        all_analyzers <- unique(data$analyzer)
        analyzer_colors_full <- setNames(rep("black", length(all_analyzers)), all_analyzers)
        analyzer_shapes_full <- setNames(rep(16, length(all_analyzers)), all_analyzers)
        analyzer_colors_full[names(analyzer_colors)] <- analyzer_colors
        analyzer_shapes_full[names(analyzer_shapes)] <- analyzer_shapes
        
        summary_data <- data %>%
                group_by(analyzer, variable, facet_label) %>%
                summarise(
                        mean = mean(value, na.rm = TRUE),
                        sd   = sd(value, na.rm = TRUE),
                        n    = n(),
                        se   = sd / sqrt(n),
                        .groups = "drop"
                )
        
        p <- ggplot(summary_data, aes(x = analyzer, y = mean, color = analyzer, shape = analyzer)) +
                geom_point(position = position_dodge(width = 0.6), size = 3) +
                geom_errorbar(
                        aes(ymin = mean - se, ymax = mean + se),
                        width = 0.2,
                        position = position_dodge(width = 0.6)
                ) +
                scale_color_manual(values = analyzer_colors_full) +
                scale_shape_manual(values = analyzer_shapes_full) +
                facet_grid(facet_label ~ ., scales = "free_y", switch = "y",
                           labeller = labeller(facet_label = label_parsed)) +
                scale_y_continuous(
                        breaks = function(limits) seq(limits[1], limits[2], length.out = 6),
                        labels = scales::label_number(accuracy = 0.1, big.mark = "")
                ) +
                labs(x = NULL, y = NULL) +
                theme_classic() +
                theme(
                        text = element_text(size = 14),
                        axis.text = element_text(size = 14),
                        axis.title = element_text(size = 14),
                        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
                        axis.text.y = element_text(hjust = 1, size = 12),
                        strip.text.y.left = element_text(size = 11, vjust = 0.5),
                        panel.border = element_rect(color = "black", fill = NA),
                        legend.position = "bottom",
                        legend.title = element_blank(),
                        plot.title = element_text(hjust = 0.5)
                ) +
                guides(color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))
        
        print(p)
        return(p)
}

# Trend Plot Function 
emitrendplot <- function(data, y = NULL) {
        library(dplyr)
        library(ggplot2)
        library(scales)

        analyzer_labels <- c(
                "crds" = "CRDS",
                "logas_ndir" = "PRONOVA",
                "logas_tdlas" = "CUBIC",
                "otice" = "OTICE"
        )
        
        if ("var" %in% names(data) && !"variable" %in% names(data)) {
                data <- data %>% rename(variable = var)
        }
        
        if (!is.null(y)) {
                data <- data %>% filter(variable %in% y)
        }
        
        if (!"value" %in% names(data)) stop("Column 'value' not found in data.")
        if (!"DATE.TIME" %in% names(data)) stop("Column 'DATE.TIME' not found in data.")
        
        summary_data <- data %>%
                group_by(DATE.TIME, analyzer, variable) %>%
                summarise(
                        mean_val = mean(value, na.rm = TRUE),
                        sd_val   = sd(value, na.rm = TRUE),
                        .groups = "drop"
                ) %>%
                filter(!is.na(mean_val))
        
        if (nrow(summary_data) == 0) {
                stop("No data available for the selected variables.")
        }
        
        facet_labels <- c(
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
        
        summary_data <- summary_data %>%
                mutate(
                        analyzer = as.character(analyzer),
                        analyzer = recode(analyzer, !!!analyzer_labels, .default = analyzer),
                        facet_label = factor(
                                facet_labels[as.character(variable)],
                                levels = facet_labels[y]
                        )
                )

        analyzer_colors <- c(
                "CUBIC" = "#1b9e77",
                "PRONOVA" = "#7570b3",
                "CRDS" = "darkgray",
                "OTICE" = "#C45A11"
        )

        analyzer_shapes <- c(
                "CUBIC" = 0,
                "PRONOVA" = 1,
                "CRDS" = 2,
                "OTICE" = 15
        )
        
        all_analyzers <- unique(summary_data$analyzer)
        analyzer_colors_full <- setNames(rep("black", length(all_analyzers)), all_analyzers)
        analyzer_shapes_full <- setNames(rep(16, length(all_analyzers)), all_analyzers)
        analyzer_colors_full[names(analyzer_colors)] <- analyzer_colors
        analyzer_shapes_full[names(analyzer_shapes)] <- analyzer_shapes
        
        x_breaks <- seq(
                from = min(summary_data$DATE.TIME, na.rm = TRUE),
                to   = max(summary_data$DATE.TIME, na.rm = TRUE),
                by   = "24 hours"
        )
        
        p <- ggplot(
                summary_data,
                aes(x = DATE.TIME, y = mean_val, color = analyzer, shape = analyzer, group = analyzer)
        ) +
                geom_line(linewidth = 0.5, alpha = 0.6, na.rm = TRUE) +
                geom_point(size = 2, na.rm = TRUE) +
                geom_errorbar(
                        aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
                        width = 0.4,
                        na.rm = TRUE
                ) +
                scale_color_manual(values = analyzer_colors_full) +
                scale_shape_manual(values = analyzer_shapes_full) +
                facet_grid(facet_label ~ ., scales = "free_y", switch = "y",
                           labeller = labeller(facet_label = label_parsed)) +
                scale_y_continuous(
                        breaks = scales::pretty_breaks(n = 6),
                        labels = scales::label_number(accuracy = 1, big.mark = "")
                ) +
                scale_x_datetime(breaks = x_breaks, date_labels = "%Y-%m-%d %H:%M") +
                labs(x = NULL, y = NULL) +
                theme_classic() +
                theme(
                        text = element_text(size = 14),
                        axis.text = element_text(size = 14),
                        axis.title = element_text(size = 14),
                        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
                        axis.text.y = element_text(hjust = 1, size = 12),
                        strip.text.y.left = element_text(size = 11, vjust = 0.5),
                        panel.border = element_rect(color = "black", fill = NA),
                        legend.position = "bottom",
                        legend.title = element_blank(),
                        plot.title = element_text(hjust = 0.5)
                ) +
                guides(color = guide_legend(nrow = 1))
        
        print(p)
        return(p)
}
