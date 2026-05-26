library(dplyr)
library(lubridate)
library(readr)
library(tidyr)

resolve_project_dir <- function() {
  project_marker <- function(path) {
    dir.exists(file.path(path, "workflows", "device_comparison")) &&
      dir.exists(file.path(path, "workflows", "crds_routine_cleaning"))
  }

  candidate_paths <- character()
  args_all <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- sub(file_arg, "", args_all[grepl(file_arg, args_all)])

  if (length(script_path) > 0) {
    candidate_paths <- c(candidate_paths, normalizePath(file.path(dirname(script_path[1]), "..", "..", ".."), winslash = "/", mustWork = FALSE))
  }

  candidate_paths <- c(
    candidate_paths,
    normalizePath(getwd(), winslash = "/", mustWork = FALSE),
    normalizePath(file.path(getwd(), "Gas_Concentration_Emission_R"), winslash = "/", mustWork = FALSE)
  )

  candidate_paths <- unique(candidate_paths[nzchar(candidate_paths)])

  for (path in candidate_paths) {
    if (project_marker(path)) {
      return(path)
    }
  }

  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

parse_mixed_datetime <- function(x, time_zone = "Europe/Berlin") {
  x_chr <- as.character(x)
  parsed <- as.POSIXct(rep(NA_character_, length(x_chr)), tz = time_zone)
  formats <- c(
    "%Y-%m-%dT%H:%M:%SZ",
    "%Y-%m-%d %H:%M:%S",
    "%Y/%m/%d %H:%M:%S"
  )

  for (fmt in formats) {
    idx <- which(is.na(parsed) & !is.na(x_chr))
    if (length(idx) == 0) break
    trial <- as.POSIXct(x_chr[idx], format = fmt, tz = time_zone)
    parsed[idx[!is.na(trial)]] <- trial[!is.na(trial)]
  }

  parsed
}

safe_mean <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

campaign_start <- as.POSIXct("2025-09-29 00:00:00", tz = "Europe/Berlin")
project_dir <- resolve_project_dir()

input_dir <- file.path(project_dir, "workflows", "device_comparison", "clean_data", "otice_clean")
output_dir <- file.path(project_dir, "workflows", "device_comparison", "clean_data", "otice_campaign")
crds_hourly_dir <- file.path(project_dir, "workflows", "crds_routine_cleaning", "clean_data", "crds_clean", "hourly_in_out_avg")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

otice_raw <- read_csv(file.path(input_dir, "otice_node_hourly_raw_all.csv"), show_col_types = FALSE) %>%
  mutate(DATE.HOUR = parse_mixed_datetime(DATE.HOUR)) %>%
  filter(!is.na(DATE.HOUR), DATE.HOUR >= campaign_start)

otice_fit <- read_csv(file.path(input_dir, "otice_node_hourly_calibrated_all.csv"), show_col_types = FALSE) %>%
  mutate(DATE.HOUR = parse_mixed_datetime(DATE.HOUR)) %>%
  filter(!is.na(DATE.HOUR), DATE.HOUR >= campaign_start) %>%
  select(DATE.HOUR, node, analyzer, crds_location, OTICE_CO2_fitted, OTICE_NH3_fitted)

crds_files <- list.files(crds_hourly_dir, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
crds_s <- bind_rows(lapply(crds_files, function(path) read_csv(path, col_types = cols(.default = col_character()), show_col_types = FALSE))) %>%
  mutate(
    DATE.HOUR = parse_mixed_datetime(DATE.HOUR),
    CO2_S = as.numeric(CO2_S),
    NH3_S = as.numeric(NH3_S)
  ) %>%
  filter(!is.na(DATE.HOUR), DATE.HOUR >= campaign_start) %>%
  group_by(DATE.HOUR) %>%
  summarise(
    CRDS_CO2_S = safe_mean(CO2_S),
    CRDS_NH3_S = safe_mean(NH3_S),
    .groups = "drop"
  )

otice_campaign <- otice_raw %>%
  left_join(otice_fit, by = c("DATE.HOUR", "node")) %>%
  left_join(crds_s, by = "DATE.HOUR") %>%
  mutate(
    delta_CO2_raw = OTICE_CO2_raw - CRDS_CO2_S,
    delta_NH3_raw = OTICE_NH3_raw - CRDS_NH3_S,
    delta_CO2_fitted = OTICE_CO2_fitted - CRDS_CO2_S,
    delta_NH3_fitted = OTICE_NH3_fitted - CRDS_NH3_S,
    month_id = format(DATE.HOUR, "%Y-%m")
  ) %>%
  select(
    DATE.HOUR, month_id, node, analyzer, crds_location,
    OTICE_CO2_raw, OTICE_NH3_raw,
    OTICE_CO2_fitted, OTICE_NH3_fitted,
    CRDS_CO2_S, CRDS_NH3_S,
    delta_CO2_raw, delta_NH3_raw,
    delta_CO2_fitted, delta_NH3_fitted
  ) %>%
  arrange(DATE.HOUR, node)

write_csv(otice_campaign, file.path(output_dir, "otice_campaign_hourly.csv"))

cat("Wrote OTICE campaign file to:\n", normalizePath(file.path(output_dir, "otice_campaign_hourly.csv"), winslash = "/", mustWork = FALSE), "\n")
