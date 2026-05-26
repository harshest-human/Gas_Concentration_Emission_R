library(dplyr)
library(lubridate)
library(readr)

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

safe_mean <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

campaign_start <- as.POSIXct("2025-09-29 00:00:00", tz = "Europe/Berlin")
project_dir <- resolve_project_dir()
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

write_csv(ndir_hourly, file.path(output_dir, "logas_ndir_campaign_hourly.csv"))

cat("Wrote LoGAS NDIR campaign file to:\n", normalizePath(file.path(output_dir, "logas_ndir_campaign_hourly.csv"), winslash = "/", mustWork = FALSE), "\n")
