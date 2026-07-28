library(dplyr)
library(lubridate)
library(readr)
library(tidyr)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args_all[grepl(file_arg, args_all)])
helper_dir <- if (length(script_path) > 0) dirname(normalizePath(script_path[1], winslash = "/", mustWork = FALSE)) else getwd()
source(file.path(helper_dir, "device_comparison_helpers.R"))

campaign_start <- campaign_start_default
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

output_path <- file.path(output_dir, "otice_campaign_hourly.csv")
write_csv(otice_campaign, output_path)

cat("Wrote OTICE campaign file to:\n", normalizePath(output_path, winslash = "/", mustWork = FALSE), "\n")
