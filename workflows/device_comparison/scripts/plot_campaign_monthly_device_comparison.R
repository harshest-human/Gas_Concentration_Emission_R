library(dplyr)
library(lubridate)
library(readr)
library(tidyr)
library(ggplot2)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args_all[grepl(file_arg, args_all)])
helper_dir <- if (length(script_path) > 0) dirname(normalizePath(script_path[1], winslash = "/", mustWork = FALSE)) else getwd()
source(file.path(helper_dir, "device_comparison_helpers.R"))

read_crds_delta <- function(project_dir, campaign_start) {
  crds_dir <- file.path(project_dir, "workflows", "crds_routine_cleaning", "clean_data", "crds_clean", "hourly_in_out_avg")
  files <- list.files(crds_dir, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
  bind_rows(lapply(files, function(path) read_csv(path, col_types = cols(.default = col_character()), show_col_types = FALSE))) %>%
    mutate(
      DATE.HOUR = parse_mixed_datetime(DATE.HOUR),
      CO2_in = as.numeric(CO2_in), CO2_S = as.numeric(CO2_S),
      CH4_in = as.numeric(CH4_in), CH4_S = as.numeric(CH4_S),
      NH3_in = as.numeric(NH3_in), NH3_S = as.numeric(NH3_S)
    ) %>%
    filter(!is.na(DATE.HOUR), DATE.HOUR >= campaign_start) %>%
    group_by(DATE.HOUR) %>%
    summarise(
      delta_CO2 = mean(CO2_in - CO2_S, na.rm = TRUE),
      delta_CH4 = mean(CH4_in - CH4_S, na.rm = TRUE),
      delta_NH3 = mean(NH3_in - NH3_S, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(analyzer = "crds")
}

build_delta_plot <- function(month_df, month_label) {
  plot_df <- month_df %>%
    pivot_longer(cols = c(delta_CO2, delta_CH4, delta_NH3), names_to = "gas", values_to = "value") %>%
    filter(!is.na(value)) %>%
    mutate(
      gas = recode(gas, delta_CO2 = "Delta CO2 (ppm)", delta_CH4 = "Delta CH4 (ppm)", delta_NH3 = "Delta NH3 (ppm)"),
      analyzer = factor(analyzer, levels = names(device_labels), labels = unname(device_labels))
    )

  ggplot(plot_df, aes(x = DATE.HOUR, y = value, color = analyzer)) +
    geom_line(linewidth = 0.6, na.rm = TRUE) +
    facet_grid(gas ~ ., scales = "free_y", switch = "y") +
    scale_color_manual(values = setNames(unname(device_colors), unname(device_labels[names(device_colors)]))) +
    labs(title = paste("Device delta comparison,", month_label), x = NULL, y = NULL, color = NULL) +
    device_plot_theme_classic()
}

project_dir <- resolve_project_dir()
campaign_start <- campaign_start_default
plot_dir <- file.path(project_dir, "workflows", "device_comparison", "result_data", "plots", "campaign_monthly_compare")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

crds_delta <- read_crds_delta(project_dir, campaign_start)

ndir_path <- file.path(project_dir, "workflows", "device_comparison", "clean_data", "logas_ndir_campaign", "logas_ndir_campaign_hourly.csv")
tdlas_path <- file.path(project_dir, "workflows", "device_comparison", "clean_data", "logas_tdlas_campaign", "logas_tdlas_campaign_hourly.csv")
otice_path <- file.path(project_dir, "workflows", "device_comparison", "clean_data", "otice_campaign", "otice_campaign_hourly.csv")

device_frames <- list(crds_delta)

if (file.exists(ndir_path)) {
  device_frames[[length(device_frames) + 1]] <- read_csv(ndir_path, show_col_types = FALSE) %>%
    mutate(DATE.HOUR = parse_mixed_datetime(DATE.HOUR), analyzer = "logas_ndir") %>%
    select(DATE.HOUR, analyzer, delta_CO2, delta_CH4, delta_NH3)
}

if (file.exists(tdlas_path)) {
  device_frames[[length(device_frames) + 1]] <- read_csv(tdlas_path, show_col_types = FALSE) %>%
    mutate(DATE.HOUR = parse_mixed_datetime(DATE.HOUR), analyzer = "logas_tdlas") %>%
    select(DATE.HOUR, analyzer, delta_CO2, delta_CH4, delta_NH3)
}

if (file.exists(otice_path)) {
  device_frames[[length(device_frames) + 1]] <- read_csv(otice_path, show_col_types = FALSE) %>%
    mutate(DATE.HOUR = parse_mixed_datetime(DATE.HOUR)) %>%
    group_by(DATE.HOUR) %>%
    summarise(
      delta_CO2 = mean(delta_CO2_fitted, na.rm = TRUE),
      delta_CH4 = NA_real_,
      delta_NH3 = mean(delta_NH3_fitted, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(analyzer = "otice")
}

device_delta <- bind_rows(device_frames) %>%
  filter(!is.na(DATE.HOUR), DATE.HOUR >= campaign_start) %>%
  mutate(month_id = format(DATE.HOUR, "%Y-%m"))

write_csv(device_delta, file.path(project_dir, "workflows", "device_comparison", "result_data", "tables", "device_delta_campaign_all.csv"))

months <- sort(unique(device_delta$month_id))

for (month_id in months) {
  month_df <- device_delta %>% filter(month_id == !!month_id)
  if (nrow(month_df) == 0) next
  month_plot <- build_delta_plot(month_df, month_id)
  ggsave(file.path(plot_dir, paste0("device_delta_", gsub("-", "_", month_id), ".png")), month_plot, width = 11, height = 7, dpi = 150)
}

cat("Wrote monthly device comparison plots to:\n", normalizePath(plot_dir, winslash = "/", mustWork = FALSE), "\n")
