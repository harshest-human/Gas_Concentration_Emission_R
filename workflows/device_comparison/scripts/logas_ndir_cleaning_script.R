# =============================================================================
# Clean logas_ndir data
#
# Outputs:
# - hourly clean delta csv files
# - monthly delta plots
# =============================================================================

library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(readr)

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

clean_dir <- file.path(project_dir, "workflows", "device_comparison", "clean_data", "logas_ndir_clean")
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
    orders = c("dmy HMS", "dmy HM"),
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
    geom_line(linewidth = 0.6, color = "#7570B3", na.rm = TRUE) +
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

logas_ndir_files <- list.files(
  path = file.path(project_dir, "workflows", "device_comparison", "raw_data", "logas_ndir_raw", "Messdaten"),
  pattern = "^differenzmessung_.*\\.txt$",
  full.names = TRUE
)

logas_ndir_hourly <- logas_ndir_files |>
  lapply(read_log_file) |>
  bind_rows() |>
  mutate(
    DATE.TIME = parse_datetime_local(`Datum Uhrzeit`),
    DATE.HOUR = floor_date(DATE.TIME, "hour")
  ) |>
  filter(!is.na(DATE.HOUR)) |>
  group_by(DATE.HOUR) |>
  summarise(
    delta_CO2 = safe_mean(as.numeric(gsub(",", ".", `CO2 in ppm`))),
    delta_CH4 = safe_mean(as.numeric(gsub(",", ".", `CH4 in ppm`))),
    delta_NH3 = safe_mean(as.numeric(gsub(",", ".", `NH3 in ppm`))),
    .groups = "drop"
  ) |>
  add_period_id() |>
  arrange(DATE.HOUR)

write_csv(logas_ndir_hourly, file.path(clean_dir, "logas_ndir_hourly_all.csv"))

for (i in seq_len(nrow(period_lookup))) {
  period_row <- period_lookup[i, ]
  period_df <- logas_ndir_hourly |>
    filter(period_id == period_row$period_id)

  if (nrow(period_df) == 0) {
    next
  }

  write_csv(period_df, file.path(clean_dir, paste0("logas_ndir_hourly_", period_row$file_tag, ".csv")))

  delta_plot <- plot_device_delta(
    period_df,
    title_text = paste("logas_ndir delta concentrations", period_row$file_tag)
  )

  ggsave(
    file.path(plot_dir, period_row$folder, paste0("logas_ndir_delta_", period_row$file_tag, ".png")),
    delta_plot,
    width = 10,
    height = 6.5,
    dpi = 150
  )
}

cat("Wrote cleaned logas_ndir outputs to:\n", normalizePath(clean_dir), "\n")
