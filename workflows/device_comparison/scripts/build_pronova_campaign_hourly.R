library(dplyr)
library(lubridate)
library(readr)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args_all[grepl(file_arg, args_all)])
helper_dir <- if (length(script_path) > 0) dirname(normalizePath(script_path[1], winslash = "/", mustWork = FALSE)) else getwd()
source(file.path(helper_dir, "device_comparison_helpers.R"))

read_log_file <- function(file_path) {
  tryCatch(
    read.table(
      file_path,
      header = TRUE,
      sep = "\t",
      dec = ",",
      check.names = FALSE,
      stringsAsFactors = FALSE,
      fill = TRUE,
      comment.char = "",
      colClasses = "character"
    ),
    error = function(e) NULL
  )
}

project_dir <- resolve_project_dir()
campaign_start <- campaign_start_default
raw_dir <- "D:/Data_Analysis_R/owncloud_sync_data/Pronova_LoGAS_raw/Messdaten"
output_dir <- file.path(project_dir, "workflows", "device_comparison", "clean_data", "logas_ndir_campaign")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ndir_files <- list.files(raw_dir, pattern = "^differenzmessung_.*\\.txt$", full.names = TRUE)

ndir_hourly <- bind_rows(lapply(ndir_files, read_log_file)) %>%
  mutate(
    DATE.TIME = parse_date_time(`Datum Uhrzeit`, orders = c("dmy HMS", "dmy HM"), tz = "Europe/Berlin", quiet = TRUE),
    DATE.HOUR = floor_date(DATE.TIME, "hour")
  ) %>%
  filter(!is.na(DATE.HOUR), DATE.HOUR >= campaign_start) %>%
  group_by(DATE.HOUR) %>%
  summarise(
    CO2_raw = safe_mean(as.numeric(gsub(",", ".", `CO2 in ppm`))),
    CH4_raw = safe_mean(as.numeric(gsub(",", ".", `CH4 in ppm`))),
    NH3_raw = safe_mean(as.numeric(gsub(",", ".", `NH3 in ppm`))),
    .groups = "drop"
  ) %>%
  mutate(
    delta_CO2 = CO2_raw,
    delta_CH4 = CH4_raw,
    delta_NH3 = NH3_raw,
    month_id = format(DATE.HOUR, "%Y-%m")
  ) %>%
  arrange(DATE.HOUR)

output_path <- file.path(output_dir, "logas_ndir_campaign_hourly.csv")
write_csv(ndir_hourly, output_path)

cat("Wrote PRONOVA campaign file to:\n", normalizePath(output_path, winslash = "/", mustWork = FALSE), "\n")
