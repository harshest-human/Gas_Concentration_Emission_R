# =============================================================================
# Export calibrated OTICE hourly data for multi-device emission comparison
#
# OTICE has calibrated inside/node concentrations. It does not have an outside
# sampling line, so this export uses CRDS location S as the shared outside
# background for OTICE delta-concentration calculations.
# =============================================================================

library(dplyr)
library(lubridate)
library(readr)

timezone_local <- "Europe/Berlin"

period_start <- ymd_hms("2025-12-09 00:00:00", tz = timezone_local)
period_end <- ymd_hms("2025-12-22 23:00:00", tz = timezone_local)

parse_hour <- function(x) {
  parsed <- suppressWarnings(ymd_hms(as.character(x), tz = timezone_local))
  missing <- is.na(parsed)
  parsed[missing] <- suppressWarnings(parse_date_time(
    as.character(x),
    orders = c("ymd HMS", "ymd HM", "Ymd HMS", "Ymd HM"),
    tz = timezone_local
  ))[missing]
  parsed
}

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args_all[grepl(file_arg, args_all)])
base_dir <- if (length(script_path) > 0) dirname(normalizePath(script_path[1])) else getwd()
tables_dir <- file.path(base_dir, "output", "OTICE_versus_CRDS_concentration", "tables")

picarro_project_dir <- "D:/Data_Analysis_R/Picarro-G2508_CRDS_gas_measurement"
out_dir <- file.path(picarro_project_dir, "data_processed", "device_hourly")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

otice_node_hourly_path <- file.path(tables_dir, "dec_node_hourly_comparison.csv")
crds_hourly_path <- file.path(picarro_project_dir, "crds_data", "crds_hour", "2025", "H_CRDS8_20251209_20251231.csv")

otice_node_hourly <- read_csv(otice_node_hourly_path, col_types = cols(.default = col_character())) |>
  mutate(datetime_hour = parse_hour(datetime_hour)) |>
  mutate(across(any_of(c("OTICE_CO2_fitted", "OTICE_NH3_fitted")), as.numeric)) |>
  filter(datetime_hour >= period_start, datetime_hour <= period_end)

crds_hourly <- read_csv(crds_hourly_path, col_types = cols(.default = col_character())) |>
  mutate(DATE.HOUR = parse_hour(DATE.HOUR)) |>
  mutate(across(any_of(c("CO2_S", "NH3_S")), as.numeric)) |>
  filter(DATE.HOUR >= period_start, DATE.HOUR <= period_end)

otice_inside <- otice_node_hourly |>
  group_by(datetime_hour) |>
  summarise(
    CO2_in = mean(OTICE_CO2_fitted, na.rm = TRUE),
    NH3_in = mean(OTICE_NH3_fitted, na.rm = TRUE),
    n_otice_nodes = n_distinct(node[!is.na(OTICE_CO2_fitted) | !is.na(OTICE_NH3_fitted)]),
    .groups = "drop"
  )

otice_for_emission <- otice_inside |>
  left_join(
    crds_hourly |>
      select(DATE.HOUR, CO2_S, NH3_S),
    by = c("datetime_hour" = "DATE.HOUR")
  ) |>
  transmute(
    DATE.HOUR = format(datetime_hour, "%Y-%m-%d %H:%M:%S"),
    analyzer = "OTICE",
    CO2_in,
    CO2_S,
    CH4_in = NA_real_,
    CH4_S = NA_real_,
    NH3_in,
    NH3_S,
    delta_CO2 = CO2_in - CO2_S,
    delta_CH4 = NA_real_,
    delta_NH3 = NH3_in - NH3_S,
    n_otice_nodes,
    outside_reference = "CRDS8 location S",
    note = "OTICE inside is 48-hour calibrated node mean; outside background is CRDS S."
  )

out_file <- file.path(out_dir, "OTICE_hourly_20251209_20251222.csv")
write_csv(otice_for_emission, out_file)

cat("Wrote OTICE hourly comparison file:\n", out_file, "\n")
cat("Rows:", nrow(otice_for_emission), "\n")
cat("Time range:", format(min(otice_for_emission$DATE.HOUR, na.rm = TRUE)), "to",
    format(max(otice_for_emission$DATE.HOUR, na.rm = TRUE)), "\n")
