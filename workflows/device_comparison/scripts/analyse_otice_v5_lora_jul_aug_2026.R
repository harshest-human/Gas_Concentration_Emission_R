#!/usr/bin/env Rscript

# Compile and compare selected OTICE_v5_LoRA analysers, 2026-07-05--2026-08-03.
# Run from the repository root with:
# Rscript workflows/device_comparison/scripts/analyse_otice_v5_lora_jul_aug_2026.R

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(lubridate)
  library(readr)
})

find_repo_root <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(path, ".git"))) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) stop("Could not find repository root.")
    path <- parent
  }
}

repo_root <- find_repo_root()
raw_dir <- file.path(repo_root, "workflows", "device_comparison", "raw_data", "otice_raw", "2026")
output_root <- file.path(repo_root, "workflows", "device_comparison", "result_data", "otice_v5_lora_2026-07-05_2026-08-03")
clean_dir <- file.path(output_root, "clean_data")
table_dir <- file.path(output_root, "tables")
plot_dir <- file.path(output_root, "plots")
dir.create(clean_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

selected_ids <- c("46", "47", "48", "49", "54", "55", "59")
analyser_levels <- c("46_out", "47", "48", "49", "54", "55", "59")
selected_files <- file.path(raw_dir, paste0("ATB_", selected_ids, "_2026-08-05.csv"))
missing_files <- selected_files[!file.exists(selected_files)]
if (length(missing_files) > 0) stop("Missing input file(s): ", paste(missing_files, collapse = ", "))

read_analyser <- function(path, id) {
  read_csv(path, show_col_types = FALSE, na = c("", "NA", "NaN")) %>%
    transmute(
      Date.Time = ymd_hms(created_at, tz = "UTC", quiet = TRUE),
      analyser = if_else(id == "46", "46_out", id),
      analyser_model = "OTICE_v5_LoRA",
      NH3 = as.numeric(field1) / 100,
      TB600_Air_Temperature = as.numeric(field2),
      TB600_Air_Relative_Humidity = as.numeric(field3),
      EE894_CO2_Average = as.numeric(field5),
      EE894_CO2_Raw = as.numeric(field6),
      EE894_Air_Temperature = as.numeric(field7),
      EE894_Air_Pressure = as.numeric(field8),
      source_file = basename(path)
    )
}

start_time <- as.POSIXct("2026-07-05 00:00:00", tz = "UTC")
end_time_exclusive <- as.POSIXct("2026-08-04 00:00:00", tz = "UTC")

clean_data <- bind_rows(Map(read_analyser, selected_files, selected_ids)) %>%
  filter(!is.na(Date.Time), Date.Time >= start_time, Date.Time < end_time_exclusive) %>%
  arrange(Date.Time, factor(analyser, levels = c("46_out", "47", "48", "49", "54", "55", "59")))

if (nrow(clean_data) == 0) stop("No observations remained after date filtering.")

measurement_columns <- c(
  "NH3", "TB600_Air_Temperature", "TB600_Air_Relative_Humidity",
  "EE894_CO2_Average", "EE894_CO2_Raw", "EE894_Air_Temperature", "EE894_Air_Pressure"
)

aggregate_period <- function(data, period_name) {
  data %>%
    mutate(Date.Time = floor_date(Date.Time, unit = period_name)) %>%
    group_by(Date.Time, analyser, analyser_model) %>%
    summarise(
      across(all_of(measurement_columns),
             list(mean = ~ mean(.x, na.rm = TRUE), sd = ~ sd(.x, na.rm = TRUE), n = ~ sum(!is.na(.x))),
             .names = "{.col}_{.fn}"),
      .groups = "drop"
    ) %>%
    mutate(across(ends_with(c("_mean", "_sd")), ~ ifelse(is.nan(.x), NA_real_, .x))) %>%
    arrange(Date.Time, analyser)
}

hourly <- aggregate_period(clean_data, "hour")

write_csv(clean_data, file.path(clean_dir, "otice_v5_lora_rowwise_clean.csv"), na = "")
write_csv(hourly, file.path(clean_dir, "otice_v5_lora_hourly.csv"), na = "")

# Match each analyser's hourly raw CO2 mean to analyser 47 before calculating RPE.
reference_hourly <- hourly %>%
  filter(analyser == "47") %>%
  select(Date.Time, reference_mean = EE894_CO2_Raw_mean)

rpe_hourly <- hourly %>%
  select(Date.Time, analyser, device_mean = EE894_CO2_Raw_mean) %>%
  inner_join(reference_hourly, by = "Date.Time") %>%
  mutate(relative_percentage_error = if_else(
    is.na(reference_mean) | reference_mean == 0,
    NA_real_,
    100 * (device_mean - reference_mean) / reference_mean
  ))

stats <- clean_data %>%
  group_by(analyser) %>%
  summarise(
    n = sum(!is.na(EE894_CO2_Raw)),
    mean = mean(EE894_CO2_Raw, na.rm = TRUE),
    sd = sd(EE894_CO2_Raw, na.rm = TRUE),
    min = min(EE894_CO2_Raw, na.rm = TRUE),
    max = max(EE894_CO2_Raw, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    rpe_hourly %>%
      group_by(analyser) %>%
      summarise(
        matched_hours = sum(!is.na(relative_percentage_error)),
        relative_percentage_error = mean(relative_percentage_error, na.rm = TRUE),
        mean_absolute_percentage_error = mean(abs(relative_percentage_error), na.rm = TRUE),
        .groups = "drop"
      ),
    by = "analyser"
  ) %>%
  right_join(tibble(analyser = analyser_levels), by = "analyser") %>%
  mutate(n = coalesce(n, 0L), matched_hours = coalesce(matched_hours, 0L)) %>%
  arrange(factor(analyser, levels = analyser_levels))

write_csv(stats, file.path(table_dir, "otice_v5_lora_co2_raw_statistics.csv"), na = "")
write_csv(rpe_hourly, file.path(table_dir, "otice_v5_lora_hourly_relative_percentage_error.csv"), na = "")

plot_raw_series <- function(data, value_column, title, y_label, filename) {
  plot_data <- data %>%
    filter(is.finite(.data[[value_column]])) %>%
    mutate(analyser = factor(analyser, levels = analyser_levels))
  p <- ggplot(plot_data, aes(Date.Time, .data[[value_column]], colour = analyser, group = analyser)) +
    geom_line(linewidth = 0.3, alpha = 0.65) +
    labs(
      title = title, subtitle = "Unaggregated cleaned observations",
      x = "Date/time (UTC)", y = y_label, colour = "Analyser"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", plot.title.position = "plot")
  ggsave(file.path(plot_dir, filename), p, width = 12, height = 6.5, dpi = 300, bg = "white")
}

plot_raw_series(
  clean_data, "EE894_CO2_Raw", "Raw CO2 time series by analyser",
  expression("EE894 raw " * CO[2] * " (ppm)"), "otice_v5_lora_co2_raw_timeseries.png"
)
plot_raw_series(
  clean_data, "NH3", "NH3 time series by analyser",
  expression(NH[3]), "otice_v5_lora_nh3_raw_timeseries.png"
)

# Smooth hourly period means and SDs independently so the ribbon remains a
# descriptive mean +/- 1 SD band (rather than a model confidence interval).
smooth_hourly_series <- function(data, mean_column, sd_column, span = 0.45) {
  data %>%
    select(Date.Time, analyser, mean = all_of(mean_column), sd = all_of(sd_column)) %>%
    filter(is.finite(mean)) %>%
    group_by(analyser) %>%
    group_modify(~ {
      series <- arrange(.x, Date.Time)
      x <- as.numeric(difftime(series$Date.Time, min(series$Date.Time), units = "hours"))
      if (length(unique(x)) < 4) {
        return(mutate(series, smooth_mean = mean, smooth_sd = sd))
      }
      mean_fit <- loess(mean ~ x, data = series, span = span, degree = 1,
                        control = loess.control(surface = "direct"))
      smooth_mean <- as.numeric(predict(mean_fit, newdata = data.frame(x = x)))
      valid_sd <- is.finite(series$sd)
      if (sum(valid_sd) >= 4) {
        sd_series <- series[valid_sd, ]
        sd_series$x_sd <- x[valid_sd]
        sd_fit <- loess(sd ~ x_sd, data = sd_series, span = span, degree = 1,
                        control = loess.control(surface = "direct"))
        smooth_sd <- as.numeric(predict(sd_fit, newdata = data.frame(x_sd = x)))
      } else {
        smooth_sd <- series$sd
      }
      mutate(series, smooth_mean = smooth_mean, smooth_sd = pmax(smooth_sd, 0))
    }) %>%
    ungroup() %>%
    mutate(analyser = factor(analyser, levels = analyser_levels))
}

plot_smoothed_hourly <- function(data, title, y_label, filename) {
  ribbon_data <- data %>% filter(is.finite(smooth_mean), is.finite(smooth_sd))
  p <- ggplot(data, aes(Date.Time, smooth_mean, colour = analyser, group = analyser)) +
    geom_ribbon(
      data = ribbon_data,
      aes(ymin = smooth_mean - smooth_sd, ymax = smooth_mean + smooth_sd, fill = analyser),
      alpha = 0.13, colour = NA
    ) +
    geom_line(linewidth = 0.85) +
    labs(
      title = title,
      subtitle = "LOESS-smoothed hourly means; shaded bands are smoothed mean +/- 1 SD",
      x = "Date/time (UTC)", y = y_label, colour = "Analyser", fill = "Analyser"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", plot.title.position = "plot")
  ggsave(file.path(plot_dir, filename), p, width = 12, height = 6.5, dpi = 300, bg = "white")
}

co2_smoothed <- smooth_hourly_series(hourly, "EE894_CO2_Raw_mean", "EE894_CO2_Raw_sd")
nh3_smoothed <- smooth_hourly_series(hourly, "NH3_mean", "NH3_sd")

plot_smoothed_hourly(
  co2_smoothed,
  "Smoothed hourly raw CO2 time series by analyser",
  expression("EE894 raw " * CO[2] * " (ppm)"),
  "otice_v5_lora_co2_raw_hourly_smoothed_mean_sd.png"
)
plot_smoothed_hourly(
  nh3_smoothed,
  "Smoothed hourly NH3 time series by analyser",
  expression(NH[3]),
  "otice_v5_lora_nh3_hourly_smoothed_mean_sd.png"
)

message("Wrote ", format(nrow(clean_data), big.mark = ","), " clean observations.")
message("Outputs: ", output_root)
