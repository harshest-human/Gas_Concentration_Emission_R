# =============================================================================
# Clean logas_tdlas data
#
# Outputs:
# - hourly clean concentration and delta csv files
# - monthly concentration and delta plots
# =============================================================================

library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(readxl)
library(readr)
library(purrr)

default_device_old <- getOption("device")
options(device = function(...) {
  pdf(file = tempfile(fileext = ".pdf"), ...)
})
on.exit(options(device = default_device_old), add = TRUE)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args_all[grepl(file_arg, args_all)])

project_dir <- if (length(script_path) > 0) {
  normalizePath(file.path(dirname(script_path[1]), "..", ".."), winslash = "/", mustWork = FALSE)
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

timezone_local <- "Europe/Berlin"

clean_dir <- file.path(project_dir, "workflows", "device_comparison", "clean_data", "logas_tdlas_clean")
plot_dir <- file.path(project_dir, "workflows", "device_comparison", "result_data", "plots")
dir.create(clean_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

period_lookup <- tibble(
  period_id = c("sep_oct", "nov_2025", "dec_2025"),
  folder = c("october2025", "nov_2025", "dec_2025"),
  file_tag = c("sep_oct", "nov_2025", "dec_2025"),
  start = as.POSIXct(c("2025-09-29 00:00:00", "2025-11-01 00:00:00", "2025-12-01 00:00:00"), tz = timezone_local),
  end = as.POSIXct(c("2025-11-01 00:00:00", "2025-12-01 00:00:00", "2026-01-01 00:00:00"), tz = timezone_local)
)

for (folder_name in period_lookup$folder) {
  dir.create(file.path(plot_dir, folder_name), recursive = TRUE, showWarnings = FALSE)
}

parse_datetime_local <- function(x) {
  parse_date_time(
    as.character(x),
    orders = c("ymd HMS", "ymd HM", "Ymd HMS", "Ymd HM"),
    tz = timezone_local,
    quiet = TRUE
  )
}

safe_mean <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  mean(x, na.rm = TRUE)
}

add_period_id <- function(df, datetime_col = "DATE.HOUR") {
  df |>
    mutate(
      period_id = case_when(
        .data[[datetime_col]] >= period_lookup$start[1] & .data[[datetime_col]] < period_lookup$end[1] ~ period_lookup$period_id[1],
        .data[[datetime_col]] >= period_lookup$start[2] & .data[[datetime_col]] < period_lookup$end[2] ~ period_lookup$period_id[2],
        .data[[datetime_col]] >= period_lookup$start[3] & .data[[datetime_col]] < period_lookup$end[3] ~ period_lookup$period_id[3],
        TRUE ~ NA_character_
      )
    ) |>
    filter(!is.na(period_id))
}

plot_device_concentrations <- function(df, title_text) {
  label_map <- c(
    "CO2_in" = "CO2 inside (ppm)",
    "CO2_S" = "CO2 outside S (ppm)",
    "CH4_in" = "CH4 inside (ppm)",
    "CH4_S" = "CH4 outside S (ppm)",
    "NH3_in" = "NH3 inside (ppm)",
    "NH3_S" = "NH3 outside S (ppm)"
  )

  plot_df <- df |>
    pivot_longer(cols = all_of(names(label_map)), names_to = "variable", values_to = "value") |>
    mutate(variable = factor(variable, levels = names(label_map), labels = label_map))

  ggplot(plot_df, aes(x = DATE.HOUR, y = value)) +
    geom_line(linewidth = 0.6, color = "#145A96", na.rm = TRUE) +
    facet_grid(variable ~ ., scales = "free_y", switch = "y") +
    labs(title = title_text, x = NULL, y = NULL) +
    theme_classic(base_size = 15) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 17),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      strip.text.y.left = element_text(size = 13),
      panel.border = element_rect(color = "black", fill = NA)
    )
}

plot_device_delta <- function(df, title_text) {
  label_map <- c(
    "delta_CO2" = "Delta CO2 (ppm)",
    "delta_CH4" = "Delta CH4 (ppm)",
    "delta_NH3" = "Delta NH3 (ppm)"
  )

  plot_df <- df |>
    pivot_longer(cols = all_of(names(label_map)), names_to = "variable", values_to = "value") |>
    mutate(variable = factor(variable, levels = names(label_map), labels = label_map))

  ggplot(plot_df, aes(x = DATE.HOUR, y = value)) +
    geom_line(linewidth = 0.6, color = "#1B9E77", na.rm = TRUE) +
    facet_grid(variable ~ ., scales = "free_y", switch = "y") +
    labs(title = title_text, x = NULL, y = NULL) +
    theme_classic(base_size = 15) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 17),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      strip.text.y.left = element_text(size = 13),
      panel.border = element_rect(color = "black", fill = NA)
    )
}

logas_tdlas_files <- list.files(
  path = file.path(project_dir, "workflows", "device_comparison", "raw_data", "logas_tdlas_raw"),
  pattern = "^[^~].*\\.xlsx$",
  full.names = TRUE
)

logas_tdlas_hourly <- map_dfr(logas_tdlas_files, read_excel) |>
  mutate(
    Time = parse_datetime_local(Time),
    DATE.HOUR = floor_date(Time, "hour")
  ) |>
  filter(!is.na(DATE.HOUR)) |>
  pivot_longer(
    cols = any_of(c("CH4", "NH3", "CO2")),
    names_to = "gas",
    values_to = "value"
  ) |>
  mutate(
    value = as.numeric(value),
    stream = ifelse(Type == 1, "in", "S"),
    variable = paste0(gas, "_", stream)
  ) |>
  group_by(DATE.HOUR, variable) |>
  summarise(value = safe_mean(value), .groups = "drop") |>
  pivot_wider(names_from = variable, values_from = value) |>
  mutate(
    delta_CO2 = CO2_in - CO2_S,
    delta_CH4 = CH4_in - CH4_S,
    delta_NH3 = NH3_in - NH3_S
  ) |>
  add_period_id() |>
  arrange(DATE.HOUR)

write_csv(logas_tdlas_hourly, file.path(clean_dir, "logas_tdlas_hourly_all.csv"))

for (i in seq_len(nrow(period_lookup))) {
  period_row <- period_lookup[i, ]
  period_df <- logas_tdlas_hourly |>
    filter(period_id == period_row$period_id)

  if (nrow(period_df) == 0) {
    next
  }

  write_csv(period_df, file.path(clean_dir, paste0("logas_tdlas_hourly_", period_row$file_tag, ".csv")))

  concentration_plot <- plot_device_concentrations(
    period_df,
    title_text = paste("logas_tdlas concentrations", period_row$file_tag)
  )

  delta_plot <- plot_device_delta(
    period_df,
    title_text = paste("logas_tdlas delta concentrations", period_row$file_tag)
  )

  ggsave(
    file.path(plot_dir, period_row$folder, paste0("logas_tdlas_concentrations_", period_row$file_tag, ".png")),
    concentration_plot,
    width = 11,
    height = 7,
    dpi = 150
  )

  ggsave(
    file.path(plot_dir, period_row$folder, paste0("logas_tdlas_delta_", period_row$file_tag, ".png")),
    delta_plot,
    width = 10,
    height = 6.5,
    dpi = 150
  )
}

cat("Wrote cleaned logas_tdlas outputs to:\n", normalizePath(clean_dir), "\n")
