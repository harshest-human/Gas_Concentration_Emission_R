library(dplyr)
library(ggplot2)
library(lubridate)
library(readr)
library(scales)
library(tidyr)

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
site_id <- "LVAT"
site_label <- "LVAT, Gross Kreutz (Havel), Brandenburg, Germany"

input_dir <- file.path(
  project_dir,
  "workflows",
  "ringversuche_analysis",
  "meta_data",
  "USA_mast_wd_ws"
)

result_dir <- file.path(
  project_dir,
  "workflows",
  "ringversuche_analysis",
  "results",
  "lvat_historical_wind"
)

plot_dir <- file.path(result_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

rplots_pdf_path <- file.path(project_dir, "Rplots.pdf")
if (file.exists(rplots_pdf_path)) {
  unlink(rplots_pdf_path, force = TRUE)
}

compass16 <- c(
  "N", "NNE", "NE", "ENE",
  "E", "ESE", "SE", "SSE",
  "S", "SSW", "SW", "WSW",
  "W", "WNW", "NW", "NNW"
)

compass8 <- c("N", "NE", "E", "SE", "S", "SW", "W", "NW")

speed_palette <- c(
  "0.0-0.5" = "#D8EEF2",
  "0.5-1.0" = "#A6D7E3",
  "1.0-2.0" = "#5BA9C9",
  "2.0-3.0" = "#2F78A8",
  "3.0-4.0" = "#1F4E79",
  ">=4.0" = "#0B2239"
)

to_compass <- function(direction_deg, n = 16) {
  labels <- if (n == 16) compass16 else compass8
  width <- 360 / n
  idx <- floor(((direction_deg %% 360) + width / 2) / width) %% n + 1
  labels[idx]
}

weighted_direction_deg <- function(direction_deg, speed = NULL) {
  if (length(direction_deg) == 0 || all(is.na(direction_deg))) {
    return(NA_real_)
  }

  if (is.null(speed)) {
    speed <- rep(1, length(direction_deg))
  }

  valid <- !is.na(direction_deg) & !is.na(speed)
  if (!any(valid)) {
    return(NA_real_)
  }

  direction_rad <- direction_deg[valid] * pi / 180
  weights <- speed[valid]

  mean_u <- weighted.mean(sin(direction_rad), w = weights)
  mean_v <- weighted.mean(cos(direction_rad), w = weights)

  if (isTRUE(all.equal(mean_u, 0)) && isTRUE(all.equal(mean_v, 0))) {
    return(NA_real_)
  }

  (atan2(mean_u, mean_v) * 180 / pi) %% 360
}

input_files <- list.files(
  input_dir,
  pattern = "_1hour\\.csv$",
  full.names = TRUE
)

if (length(input_files) == 0) {
  stop("No 1-hour wind files found in: ", input_dir)
}

wind_raw <- lapply(input_files, function(path) {
  read_csv(path, show_col_types = FALSE) |>
    mutate(source_file = basename(path))
}) |>
  bind_rows()

wind_data <- wind_raw |>
  mutate(
    datetime = suppressWarnings(dmy_hms(paste(date, time), tz = timezone_local)),
    wind_speed = as.numeric(wind_speed),
    wind_direction = as.numeric(wind_direction)
  ) |>
  filter(!is.na(datetime), !is.na(wind_speed), !is.na(wind_direction)) |>
  arrange(datetime) |>
  mutate(
    month = floor_date(datetime, "month"),
    date_only = as.Date(datetime, tz = timezone_local),
    analysis_period = case_when(
      datetime < ymd_hms("2025-04-01 00:00:00", tz = timezone_local) ~ "Before April 2025",
      datetime >= ymd_hms("2025-04-01 00:00:00", tz = timezone_local) &
        datetime < ymd_hms("2025-05-01 00:00:00", tz = timezone_local) ~ "April 2025",
      TRUE ~ "After April 2025"
    ),
    direction_16 = factor(to_compass(wind_direction, n = 16), levels = compass16),
    direction_8 = factor(to_compass(wind_direction, n = 8), levels = compass8),
    speed_bin = cut(
      wind_speed,
      breaks = c(0, 0.5, 1, 2, 3, 4, Inf),
      labels = c("0.0-0.5", "0.5-1.0", "1.0-2.0", "2.0-3.0", "3.0-4.0", ">=4.0"),
      include.lowest = TRUE,
      right = FALSE
    )
  )

if (nrow(wind_data) == 0) {
  stop("No valid wind observations remained after parsing and NA filtering.")
}

comparison_data <- wind_data |>
  filter(analysis_period %in% c("Before April 2025", "April 2025")) |>
  mutate(
    analysis_period = factor(
      analysis_period,
      levels = c("Before April 2025", "April 2025")
    )
  )

overall_stats <- wind_data |>
  summarise(
    observations = n(),
    start_time = min(datetime),
    end_time = max(datetime),
    mean_speed = mean(wind_speed, na.rm = TRUE),
    median_speed = median(wind_speed, na.rm = TRUE),
    p90_speed = quantile(wind_speed, probs = 0.9, na.rm = TRUE),
    max_speed = max(wind_speed, na.rm = TRUE),
    calm_lt_05_pct = mean(wind_speed < 0.5, na.rm = TRUE) * 100
  )

sector_summary <- wind_data |>
  count(direction_8, sort = TRUE) |>
  mutate(percentage = n / sum(n) * 100)

top_three_sectors <- sector_summary |>
  slice_head(n = 3)

monthly_summary <- wind_data |>
  group_by(month) |>
  summarise(
    observations = n(),
    mean_speed = mean(wind_speed, na.rm = TRUE),
    median_speed = median(wind_speed, na.rm = TRUE),
    p90_speed = quantile(wind_speed, probs = 0.9, na.rm = TRUE),
    max_speed = max(wind_speed, na.rm = TRUE),
    prevailing_direction_deg = weighted_direction_deg(wind_direction, speed = wind_speed),
    .groups = "drop"
  ) |>
  mutate(
    prevailing_direction = factor(
      to_compass(prevailing_direction_deg, n = 16),
      levels = compass16
    ),
    month_label = format(month, "%Y-%m")
  )

monthly_direction_segments <- monthly_summary |>
  arrange(month) |>
  mutate(
    wrap_break = if_else(
      row_number() == 1,
      FALSE,
      abs(prevailing_direction_deg - lag(prevailing_direction_deg)) > 180
    ),
    segment_id = cumsum(wrap_break)
  )

windrose_data <- wind_data |>
  count(direction_16, speed_bin) |>
  complete(direction_16, speed_bin, fill = list(n = 0)) |>
  mutate(percentage = n / sum(n) * 100)

period_summary <- comparison_data |>
  group_by(analysis_period) |>
  summarise(
    observations = n(),
    period_start = min(datetime),
    period_end = max(datetime),
    mean_speed = mean(wind_speed, na.rm = TRUE),
    median_speed = median(wind_speed, na.rm = TRUE),
    p90_speed = quantile(wind_speed, probs = 0.9, na.rm = TRUE),
    max_speed = max(wind_speed, na.rm = TRUE),
    calm_lt_05_pct = mean(wind_speed < 0.5, na.rm = TRUE) * 100,
    prevailing_direction_deg = weighted_direction_deg(wind_direction, speed = wind_speed),
    .groups = "drop"
  ) |>
  mutate(
    prevailing_direction = factor(
      to_compass(prevailing_direction_deg, n = 16),
      levels = compass16
    )
  )

period_sector_summary <- comparison_data |>
  count(analysis_period, direction_8) |>
  group_by(analysis_period) |>
  mutate(percentage = n / sum(n) * 100) |>
  arrange(analysis_period, desc(percentage), direction_8) |>
  ungroup()

daily_period_speed <- comparison_data |>
  group_by(analysis_period, date_only) |>
  summarise(
    mean_speed = mean(wind_speed, na.rm = TRUE),
    p90_speed = quantile(wind_speed, probs = 0.9, na.rm = TRUE),
    .groups = "drop"
  )

period_windrose_data <- comparison_data |>
  count(analysis_period, direction_16, speed_bin) |>
  complete(analysis_period, direction_16, speed_bin, fill = list(n = 0)) |>
  group_by(analysis_period) |>
  mutate(percentage = n / sum(n) * 100) |>
  ungroup()

make_period_windrose <- function(period_name) {
  period_plot_data <- period_windrose_data |>
    filter(analysis_period == period_name)

  ggplot(period_plot_data, aes(x = direction_16, y = percentage, fill = speed_bin)) +
    geom_col(color = "white", linewidth = 0.25, width = 1) +
    coord_polar(start = -pi / 16) +
    scale_fill_manual(values = speed_palette, drop = FALSE) +
    scale_y_continuous(
      labels = label_number(accuracy = 0.1),
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(
      title = paste("Windrose at LVAT:", period_name),
      subtitle = "Hourly wind direction frequency stratified by wind speed class",
      x = NULL,
      y = "Frequency (%)",
      fill = "Wind speed\n(m/s)",
      caption = "Source: USA 16 mast hourly data"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 18),
      plot.subtitle = element_text(size = 12),
      legend.position = "right"
    )
}

make_period_speed_plot <- function(period_name) {
  period_plot_data <- daily_period_speed |>
    filter(analysis_period == period_name)

  ggplot(period_plot_data, aes(x = date_only, y = mean_speed)) +
    geom_col(fill = "#2F78A8", width = 0.85) +
    geom_line(aes(y = p90_speed, group = 1), color = "#B24C2F", linewidth = 0.9) +
    geom_point(aes(y = p90_speed), color = "#B24C2F", size = 1.8) +
    scale_x_date(date_breaks = "7 days", date_labels = "%Y-%m-%d") +
    scale_y_continuous(
      name = "Wind speed (m/s)",
      limits = c(0, max(period_plot_data$p90_speed, na.rm = TRUE) * 1.15),
      labels = label_number(accuracy = 0.1)
    ) +
    labs(
      title = paste("Daily Wind Speed at LVAT:", period_name),
      subtitle = "Bars show daily mean hourly wind speed; red line shows daily 90th percentile",
      x = NULL,
      caption = "Source: USA 16 mast hourly data"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold", size = 18),
      plot.subtitle = element_text(size = 12)
    )
}

monthly_speed_plot <- ggplot(monthly_summary, aes(x = month, y = mean_speed)) +
  geom_col(fill = "#2F78A8", width = 25) +
  geom_line(aes(y = p90_speed, group = 1), color = "#B24C2F", linewidth = 1) +
  geom_point(aes(y = p90_speed), color = "#B24C2F", size = 2.5) +
  scale_x_datetime(
    date_breaks = "1 month",
    date_labels = "%Y-%m",
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(
    name = "Wind speed (m/s)",
    limits = c(0, max(monthly_summary$p90_speed) * 1.15),
    labels = label_number(accuracy = 0.1)
  ) +
  labs(
    title = "Historical Monthly Wind Speed at LVAT",
    subtitle = "Bars show monthly mean hourly wind speed; red line shows monthly 90th percentile",
    x = NULL,
    caption = "Source: USA 16 mast hourly data"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 12)
  )

monthly_direction_plot <- ggplot(
  monthly_direction_segments,
  aes(x = month, y = prevailing_direction_deg)
) +
  geom_line(aes(group = segment_id), color = "#1F4E79", linewidth = 1) +
  geom_point(aes(size = mean_speed), color = "#D98C2B") +
  geom_text(aes(label = prevailing_direction), nudge_y = 12, size = 3.8, color = "#1F4E79") +
  scale_size_continuous(name = "Mean speed (m/s)", range = c(2.5, 7)) +
  scale_x_datetime(
    date_breaks = "1 month",
    date_labels = "%Y-%m",
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(
    breaks = c(0, 45, 90, 135, 180, 225, 270, 315),
    labels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW"),
    limits = c(0, 360)
  ) +
  labs(
    title = "Monthly Prevailing Wind Direction at LVAT",
    subtitle = "Point size reflects mean monthly wind speed",
    x = NULL,
    y = "Prevailing wind direction",
    caption = "Source: USA 16 mast hourly data"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 12)
  )

windrose_plot <- ggplot(
  windrose_data,
  aes(x = direction_16, y = percentage, fill = speed_bin)
) +
  geom_col(color = "white", linewidth = 0.25, width = 1) +
  coord_polar(start = -pi / 16) +
  scale_fill_manual(values = speed_palette, drop = FALSE) +
  scale_y_continuous(
    labels = label_number(accuracy = 0.1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "Historical Windrose at LVAT",
    subtitle = "Hourly wind direction frequency stratified by wind speed class",
    x = NULL,
    y = "Frequency (%)",
    fill = "Wind speed\n(m/s)",
    caption = "Source: USA 16 mast hourly data"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 12),
    legend.position = "right"
  )

summary_table_path <- file.path(result_dir, "lvat_historical_wind_monthly_summary.csv")
write_csv(monthly_summary, summary_table_path)

windrose_path <- file.path(plot_dir, "lvat_historical_windrose.png")
speed_trend_path <- file.path(plot_dir, "lvat_historical_wind_speed_trend.png")
direction_trend_path <- file.path(plot_dir, "lvat_historical_prevailing_direction_trend.png")
before_april_windrose_path <- file.path(plot_dir, "lvat_before_april_2025_windrose.png")
april_windrose_path <- file.path(plot_dir, "lvat_april_2025_windrose.png")
before_april_speed_path <- file.path(plot_dir, "lvat_before_april_2025_daily_wind_speed.png")
april_speed_path <- file.path(plot_dir, "lvat_april_2025_daily_wind_speed.png")

ggsave(
  filename = windrose_path,
  plot = windrose_plot,
  width = 12,
  height = 12,
  units = "in",
  dpi = 400,
  bg = "white"
)

ggsave(
  filename = speed_trend_path,
  plot = monthly_speed_plot,
  width = 14,
  height = 8,
  units = "in",
  dpi = 400,
  bg = "white"
)

ggsave(
  filename = direction_trend_path,
  plot = monthly_direction_plot,
  width = 14,
  height = 8,
  units = "in",
  dpi = 400,
  bg = "white"
)

ggsave(
  filename = before_april_windrose_path,
  plot = make_period_windrose("Before April 2025"),
  width = 12,
  height = 12,
  units = "in",
  dpi = 400,
  bg = "white"
)

ggsave(
  filename = april_windrose_path,
  plot = make_period_windrose("April 2025"),
  width = 12,
  height = 12,
  units = "in",
  dpi = 400,
  bg = "white"
)

ggsave(
  filename = before_april_speed_path,
  plot = make_period_speed_plot("Before April 2025"),
  width = 14,
  height = 8,
  units = "in",
  dpi = 400,
  bg = "white"
)

ggsave(
  filename = april_speed_path,
  plot = make_period_speed_plot("April 2025"),
  width = 14,
  height = 8,
  units = "in",
  dpi = 400,
  bg = "white"
)

top_sector_text <- paste0(
  top_three_sectors$direction_8,
  " (",
  sprintf("%.1f", top_three_sectors$percentage),
  "%)"
) |>
  paste(collapse = ", ")

highest_month <- monthly_summary |>
  slice_max(order_by = mean_speed, n = 1, with_ties = FALSE)

lowest_month <- monthly_summary |>
  slice_min(order_by = mean_speed, n = 1, with_ties = FALSE)

period_summary_path <- file.path(result_dir, "lvat_before_april_vs_april_2025_summary.csv")
period_sector_summary_path <- file.path(result_dir, "lvat_before_april_vs_april_2025_direction_frequency.csv")
write_csv(period_summary, period_summary_path)
write_csv(period_sector_summary, period_sector_summary_path)

before_april_stats <- period_summary |>
  filter(analysis_period == "Before April 2025")

april_stats <- period_summary |>
  filter(analysis_period == "April 2025")

top_directions_for_text <- function(period_name) {
  period_sector_summary |>
    filter(analysis_period == period_name) |>
    slice_head(n = 3) |>
    transmute(text = paste0(direction_8, " (", sprintf("%.1f", percentage), "%)")) |>
    pull(text) |>
    paste(collapse = ", ")
}

before_april_top_text <- top_directions_for_text("Before April 2025")
april_top_text <- top_directions_for_text("April 2025")

methods_results_text <- paste0(
  "Historical local wind conditions for ", site_label,
  " were characterized from the available hourly USA 16 mast dataset using all 1-hour observations from ",
  format(overall_stats$start_time, "%Y-%m-%d %H:%M"),
  " to ",
  format(overall_stats$end_time, "%Y-%m-%d %H:%M"),
  " local time (n = ",
  format(overall_stats$observations, big.mark = ","),
  " valid records; Europe/Berlin). The two hourly source files in ",
  basename(input_dir),
  " were merged, timestamps were parsed to local time, and the single daylight-saving gap record with missing wind speed and direction was excluded. Wind direction was classified into 16-sector and 8-sector compass classes, and wind speed was grouped into 0.0-0.5, 0.5-1.0, 1.0-2.0, 2.0-3.0, 3.0-4.0, and >=4.0 m/s classes for windrose visualization. Across the full period, hourly wind speed averaged ",
  sprintf("%.2f", overall_stats$mean_speed),
  " m/s (median ",
  sprintf("%.2f", overall_stats$median_speed),
  " m/s; 90th percentile ",
  sprintf("%.2f", overall_stats$p90_speed),
  " m/s; maximum ",
  sprintf("%.2f", overall_stats$max_speed),
  " m/s), while ",
  sprintf("%.1f", overall_stats$calm_lt_05_pct),
  "% of hours remained below 0.5 m/s. The prevailing wind regime was dominated by ",
  top_sector_text,
  ", indicating mainly westerly to northwesterly flow during the historical record. Monthly mean wind speed was highest in ",
  format(highest_month$month, "%Y-%m"),
  " (",
  sprintf("%.2f", highest_month$mean_speed),
  " m/s) and lowest in ",
  format(lowest_month$month, "%Y-%m"),
  " (",
  sprintf("%.2f", lowest_month$mean_speed),
  " m/s). Monthly prevailing direction was most often westerly to northwesterly, with short periods of northeasterly flow in 2025-02 and southwesterly flow in 2025-04. These outputs provide the local wind context for interpreting air-flow exposure and potential building-influenced transport conditions around the LVAT experimental area."
)

text_path <- file.path(result_dir, "lvat_historical_wind_methods_results.txt")
writeLines(methods_results_text, con = text_path, useBytes = TRUE)

comparison_text <- paste0(
  "For interpretation of the April 2025 experiment period at LVAT, the available mast record was separated into a baseline period before April 2025 and the month of April 2025 itself. ",
  "Before April 2025 (",
  format(before_april_stats$period_start, "%Y-%m-%d %H:%M"),
  " to ",
  format(before_april_stats$period_end, "%Y-%m-%d %H:%M"),
  "; n = ",
  format(before_april_stats$observations, big.mark = ","),
  " hourly observations), mean wind speed was ",
  sprintf("%.2f", before_april_stats$mean_speed),
  " m/s (median ",
  sprintf("%.2f", before_april_stats$median_speed),
  " m/s; 90th percentile ",
  sprintf("%.2f", before_april_stats$p90_speed),
  " m/s), and the dominant directions were ",
  before_april_top_text,
  ", with a speed-weighted prevailing direction of ",
  as.character(before_april_stats$prevailing_direction),
  ". This indicates a baseline regime characterized mainly by westerly to northwesterly flow. ",
  "In April 2025 (",
  format(april_stats$period_start, "%Y-%m-%d %H:%M"),
  " to ",
  format(april_stats$period_end, "%Y-%m-%d %H:%M"),
  "; n = ",
  format(april_stats$observations, big.mark = ","),
  " hourly observations), mean wind speed was ",
  sprintf("%.2f", april_stats$mean_speed),
  " m/s (median ",
  sprintf("%.2f", april_stats$median_speed),
  " m/s; 90th percentile ",
  sprintf("%.2f", april_stats$p90_speed),
  " m/s), and the dominant directions were ",
  april_top_text,
  ", while the speed-weighted prevailing direction was ",
  as.character(april_stats$prevailing_direction),
  ". Relative to the pre-April baseline, April 2025 was slightly calmer overall and showed a broader directional spread, with substantial easterly, southeasterly, and westerly contributions rather than a strongly northwesterly regime. This makes April 2025 meteorologically distinct enough to justify showing separately when discussing experimental transport and exposure conditions around the barn and mast."
)

comparison_text_path <- file.path(result_dir, "lvat_before_april_vs_april_2025_methods_results.txt")
writeLines(comparison_text, con = comparison_text_path, useBytes = TRUE)

message("Analysis complete.")
message("Result directory: ", result_dir)
