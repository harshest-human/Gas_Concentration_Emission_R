library(dplyr)
library(lubridate)
library(readr)
library(readxl)
library(tidyr)

resolve_project_dir <- function() {
  project_marker <- function(path) dir.exists(file.path(path, "workflows", "device_comparison"))
  candidate_paths <- c(
    normalizePath(getwd(), winslash = "/", mustWork = FALSE),
    normalizePath(file.path(getwd(), "Gas_Concentration_Emission_R"), winslash = "/", mustWork = FALSE)
  )
  candidate_paths <- unique(candidate_paths[nzchar(candidate_paths)])
  for (path in candidate_paths) if (project_marker(path)) return(path)
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

safe_mean <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

campaign_start <- as.POSIXct("2025-09-29 00:00:00", tz = "Europe/Berlin")
project_dir <- resolve_project_dir()
raw_dir <- "D:/Data_Analysis_R/owncloud_sync_data/Cubic_raw"
output_dir <- file.path(project_dir, "workflows", "device_comparison", "clean_data", "logas_tdlas_campaign")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

tdlas_files <- list.files(raw_dir, pattern = "\\.xlsx$", full.names = TRUE)

tdlas_hourly <- bind_rows(lapply(tdlas_files, function(path) {
  read_excel(path, col_types = c("text", "numeric", "numeric", "numeric", "numeric", "numeric"))
})) %>%
  mutate(
    Time = ymd_hms(Time, tz = "Europe/Berlin", quiet = TRUE),
    DATE.HOUR = floor_date(Time, "hour"),
    stream = if_else(Type == 1, "in", "S")
  ) %>%
  filter(!is.na(DATE.HOUR), DATE.HOUR >= campaign_start) %>%
  pivot_longer(cols = any_of(c("CO2", "CH4", "NH3")), names_to = "gas", values_to = "value") %>%
  mutate(value = as.numeric(value), variable = paste0(gas, "_", stream)) %>%
  group_by(DATE.HOUR, variable) %>%
  summarise(value = safe_mean(value), .groups = "drop") %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  mutate(
    delta_CO2 = CO2_in - CO2_S,
    delta_CH4 = CH4_in - CH4_S,
    delta_NH3 = NH3_in - NH3_S,
    month_id = format(DATE.HOUR, "%Y-%m")
  ) %>%
  arrange(DATE.HOUR)

write_csv(tdlas_hourly, file.path(output_dir, "logas_tdlas_campaign_hourly.csv"))

cat("Wrote LoGAS TDLAS campaign file to:\n", normalizePath(file.path(output_dir, "logas_tdlas_campaign_hourly.csv"), winslash = "/", mustWork = FALSE), "\n")
