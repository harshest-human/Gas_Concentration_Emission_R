library(dplyr)
library(lubridate)
library(readr)
library(tidyr)
library(ggplot2)

resolve_project_dir <- function() {
  project_marker <- function(path) {
    dir.exists(file.path(path, "workflows", "device_comparison")) &&
      dir.exists(file.path(path, "workflows", "crds_routine_cleaning"))
  }
  candidate_paths <- c(
    normalizePath(getwd(), winslash = "/", mustWork = FALSE),
    normalizePath(file.path(getwd(), "Gas_Concentration_Emission_R"), winslash = "/", mustWork = FALSE)
  )
  candidate_paths <- unique(candidate_paths[nzchar(candidate_paths)])
  for (path in candidate_paths) if (project_marker(path)) return(path)
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

parse_mixed_datetime <- function(x, time_zone = "Europe/Berlin") {
  x_chr <- as.character(x)
  parsed <- as.POSIXct(rep(NA_character_, length(x_chr)), tz = time_zone)
  formats <- c("%Y-%m-%dT%H:%M:%SZ", "%Y-%m-%d %H:%M:%S", "%Y/%m/%d %H:%M:%S")
  for (fmt in formats) {
    idx <- which(is.na(parsed) & !is.na(x_chr))
    if (length(idx) == 0) break
    trial <- as.POSIXct(x_chr[idx], format = fmt, tz = time_zone)
    parsed[idx[!is.na(trial)]] <- trial[!is.na(trial)]
  }
  parsed
}

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
      gas = recode(gas, delta_CO2 = "Delta CO2 (ppm)", delta_CH4 = "Delta CH4 (ppm)", delta_NH3 = "Delta NH3 (ppm)")
    )

  ggplot(plot_df, aes(x = DATE.HOUR, y = value, color = analyzer)) +
    geom_line(linewidth = 0.6, na.rm = TRUE) +
    facet_grid(gas ~ ., scales = "free_y", switch = "y") +
    scale_color_manual(values = c(crds = "#4B5563", logas_ndir = "#7570B3", logas_tdlas = "#1B9E77", otice = "#C45A11")) +
    labs(title = paste("Device delta comparison,", month_label), x = NULL, y = NULL, color = NULL) +
    theme_classic(base_size = 15) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", hjust = 0.5, size = 17),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      strip.text.y.left = element_text(size = 13),
      panel.border = element_rect(color = "black", fill = NA)
    )
}

project_dir <- resolve_project_dir()
campaign_start <- as.POSIXct("2025-09-29 00:00:00", tz = "Europe/Berlin")
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
