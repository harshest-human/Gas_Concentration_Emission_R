library(dplyr)
library(lubridate)
library(readr)
library(readxl)
library(tidyr)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args_all[grepl(file_arg, args_all)])
helper_dir <- if (length(script_path) > 0) dirname(normalizePath(script_path[1], winslash = "/", mustWork = FALSE)) else getwd()
source(file.path(helper_dir, "device_comparison_helpers.R"))

project_dir <- resolve_project_dir()
campaign_start <- campaign_start_default
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

output_path <- file.path(output_dir, "logas_tdlas_campaign_hourly.csv")
write_csv(tdlas_hourly, output_path)

cat("Wrote CUBIC campaign file to:\n", normalizePath(output_path, winslash = "/", mustWork = FALSE), "\n")
