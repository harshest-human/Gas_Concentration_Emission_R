# =============================================================================
# OTICE sensor aging analysis
# December-focused node-level comparison against mapped CRDS reference
# using hourly averages. Calibration is learned from the first 24 hours
# of the December campaign and applied across the rest of December.
# =============================================================================

library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(patchwork)

# -----------------------------------------------------------------------------
# 1. Paths and settings
# -----------------------------------------------------------------------------
timezone_local <- "Europe/Berlin"
args_trailing <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args_trailing[grepl(file_arg, args_trailing)])
base_dir <- if (length(script_path) > 0) dirname(normalizePath(script_path[1])) else getwd()

otice_dir <- file.path(base_dir, "processed_data", "OTICE_processed")
crds_dir  <- file.path(base_dir, "processed_data", "CRDS8_processed")
out_dir   <- file.path(base_dir, "output", "OTICE_sensor_aging_december")
plots_hourly_dir <- file.path(out_dir, "plots", "hourly")
plots_daily_dir  <- file.path(out_dir, "plots", "daily")
tables_hourly_dir <- file.path(out_dir, "tables", "hourly")
tables_daily_dir  <- file.path(out_dir, "tables", "daily")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_hourly_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_daily_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_hourly_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_daily_dir, recursive = TRUE, showWarnings = FALSE)

unlink(list.files(plots_hourly_dir, full.names = TRUE), force = TRUE)
unlink(list.files(plots_daily_dir, full.names = TRUE), force = TRUE)
unlink(list.files(tables_hourly_dir, full.names = TRUE), force = TRUE)
unlink(list.files(tables_daily_dir, full.names = TRUE), force = TRUE)

# -----------------------------------------------------------------------------
# 2. Helper functions
# -----------------------------------------------------------------------------
safe_mean <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

clean_otice_location <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("^O", "", x)
  x <- sub("^0+", "", x)
  x
}

clean_crds_location <- function(x) {
  x <- tolower(trimws(as.character(x)))
  numeric_mask <- grepl("^[0-9]+$", x)
  x[numeric_mask] <- as.character(as.integer(x[numeric_mask]))
  x
}

build_campaign_reference <- function(period_id, period_label, start_time, end_time,
                                     analyzer, reference_rows) {
  bind_rows(lapply(reference_rows, function(row) {
    tibble(
      period_id       = period_id,
      period_label    = period_label,
      start_time      = as.POSIXct(start_time, tz = timezone_local),
      end_time        = as.POSIXct(end_time,   tz = timezone_local),
      node            = row$node,
      analyzer        = analyzer,
      crds_location   = as.character(row$crds_location),
      nh3_sensor_id   = as.character(row$nh3_sensor_id),
      co2_sensor_id   = as.character(row$co2_sensor_id)
    )
  }))
}

save_plot <- function(plot_obj, filename, width = 14, height = 8) {
  ggsave(filename = filename, plot = plot_obj, width = width, height = height, dpi = 150)
}

get_plot_style <- function(gas_label) {
  if (gas_label == "CO2") {
    return(list(
      raw_label = "OTICE raw",
      cal_label = "OTICE after first-24h campaign calibration",
      raw_color = "#79BCE8",
      cal_color = "#0B4F8A",
      ref_color = "#3F3F46"
    ))
  }

  list(
    raw_label = "OTICE raw",
    cal_label = "OTICE after first-24h campaign calibration",
    raw_color = "#7BC67B",
    cal_color = "#1B7F3A",
    ref_color = "#3F3F46"
  )
}

build_first24h_baseline <- function(df, sensor_episode_col) {
  episode_starts <- df |>
    group_by(.data[[sensor_episode_col]]) |>
    summarise(baseline_start = min(datetime_hour, na.rm = TRUE), .groups = "drop")

  df |>
    left_join(episode_starts, by = setNames(sensor_episode_col, sensor_episode_col)) |>
    filter(datetime_hour >= baseline_start,
           datetime_hour < baseline_start + days(1))
}

make_node_timeseries_plot <- function(df, node_id, gas_label, raw_col, ref_col, pred_col, y_label,
                                      baseline_start, baseline_end, period_caption) {
  style <- get_plot_style(gas_label)
  color_values <- setNames(
    c(style$ref_color, style$raw_color, style$cal_color),
    c("CRDS reference", style$raw_label, style$cal_label)
  )
  linetype_values <- setNames(
    c("solid", "dotted", "solid"),
    c("CRDS reference", style$raw_label, style$cal_label)
  )

  node_df <- df |>
    filter(node == node_id) |>
    select(
      datetime_hour,
      crds_location,
      raw_value = all_of(raw_col),
      ref_value = all_of(ref_col),
      calibrated_value = all_of(pred_col)
    ) |>
    pivot_longer(
      cols = c(ref_value, raw_value, calibrated_value),
      names_to = "source",
      values_to = "value"
    ) |>
    mutate(
      source = recode(
        source,
        ref_value = "CRDS reference",
        raw_value = style$raw_label,
        calibrated_value = style$cal_label
      )
    )

  crds_location_label <- node_df$crds_location |>
    unique() |>
    (\(x) x[!is.na(x)])()
  if (length(crds_location_label) == 0) {
    crds_location_label <- "unknown"
  } else {
    crds_location_label <- paste(crds_location_label, collapse = ", ")
  }

  ggplot(node_df, aes(x = datetime_hour, y = value, color = source, linetype = source)) +
    annotate(
      "rect",
      xmin = baseline_start,
      xmax = baseline_end,
      ymin = -Inf,
      ymax = Inf,
      alpha = 0.08,
      fill = "goldenrod"
    ) +
    geom_line(linewidth = 0.9, alpha = 0.9, na.rm = TRUE) +
    geom_point(size = 1.0, alpha = 0.5, na.rm = TRUE) +
    scale_color_manual(values = color_values) +
    scale_linetype_manual(values = linetype_values) +
    scale_x_datetime(
      date_breaks = "1 day",
      date_labels = "%d %b",
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    labs(
      title = paste0(
        gas_label, ": OTICE node ", node_id,
        " compared with CRDS sampling point ", crds_location_label
      ),
      subtitle = paste(
        paste0(period_caption, " hourly comparison."),
        "Yellow band marks the first 24 hours used for calibration."
      ),
      x = NULL,
      y = y_label,
      color = NULL,
      linetype = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold")
    )
}

make_node_daily_plot <- function(df, node_id, gas_label, raw_col, ref_col, pred_col, y_label,
                                 baseline_start, baseline_end, period_caption) {
  style <- get_plot_style(gas_label)
  color_values <- setNames(
    c(style$ref_color, style$raw_color, style$cal_color),
    c("CRDS reference", style$raw_label, style$cal_label)
  )
  linetype_values <- setNames(
    c("solid", "dotted", "solid"),
    c("CRDS reference", style$raw_label, style$cal_label)
  )

  node_df <- df |>
    filter(node == node_id) |>
    select(
      date_day,
      crds_location,
      raw_value = all_of(raw_col),
      ref_value = all_of(ref_col),
      calibrated_value = all_of(pred_col)
    ) |>
    pivot_longer(
      cols = c(ref_value, raw_value, calibrated_value),
      names_to = "source",
      values_to = "value"
    ) |>
    mutate(
      source = recode(
        source,
        ref_value = "CRDS reference",
        raw_value = style$raw_label,
        calibrated_value = style$cal_label
      )
    )

  crds_location_label <- node_df$crds_location |>
    unique() |>
    (\(x) x[!is.na(x)])()
  if (length(crds_location_label) == 0) {
    crds_location_label <- "unknown"
  } else {
    crds_location_label <- paste(crds_location_label, collapse = ", ")
  }

  ggplot(node_df, aes(x = date_day, y = value, color = source, linetype = source)) +
    annotate(
      "rect",
      xmin = as.Date(baseline_start, tz = timezone_local),
      xmax = as.Date(baseline_end, tz = timezone_local),
      ymin = -Inf,
      ymax = Inf,
      alpha = 0.08,
      fill = "goldenrod"
    ) +
    geom_line(linewidth = 0.9, alpha = 0.9, na.rm = TRUE) +
    geom_point(size = 1.6, alpha = 0.75, na.rm = TRUE) +
    scale_color_manual(values = color_values) +
    scale_linetype_manual(values = linetype_values) +
    scale_x_date(
      date_breaks = "1 day",
      date_labels = "%d %b",
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    labs(
      title = paste0(
        gas_label, ": OTICE node ", node_id,
        " compared with CRDS sampling point ", crds_location_label
      ),
      subtitle = paste(
        paste0(period_caption, " daily averages."),
        "Yellow band marks the first 24 hours used for calibration."
      ),
      x = NULL,
      y = y_label,
      color = NULL,
      linetype = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold")
    )
}

make_overall_timeseries_plot <- function(df, gas_label, raw_col, ref_col, pred_col, y_label,
                                         baseline_start, baseline_end, period_caption) {
  style <- get_plot_style(gas_label)
  color_values <- setNames(
    c(style$ref_color, style$raw_color, style$cal_color),
    c("CRDS reference mean", style$raw_label, style$cal_label)
  )
  linetype_values <- setNames(
    c("solid", "dotted", "solid"),
    c("CRDS reference mean", style$raw_label, style$cal_label)
  )

  plot_df <- df |>
    select(
      datetime_hour,
      raw_value = all_of(raw_col),
      ref_value = all_of(ref_col),
      calibrated_value = all_of(pred_col)
    ) |>
    pivot_longer(
      cols = c(ref_value, raw_value, calibrated_value),
      names_to = "source",
      values_to = "value"
    ) |>
    mutate(
      source = recode(
        source,
        ref_value = "CRDS reference mean",
        raw_value = style$raw_label,
        calibrated_value = style$cal_label
      )
    )

  ggplot(plot_df, aes(x = datetime_hour, y = value, color = source, linetype = source)) +
    annotate(
      "rect",
      xmin = baseline_start,
      xmax = baseline_end,
      ymin = -Inf,
      ymax = Inf,
      alpha = 0.08,
      fill = "goldenrod"
    ) +
    geom_line(linewidth = 1.0, alpha = 0.9, na.rm = TRUE) +
    geom_point(size = 1.2, alpha = 0.5, na.rm = TRUE) +
    scale_color_manual(values = color_values) +
    scale_linetype_manual(values = linetype_values) +
    scale_x_datetime(
      date_breaks = "1 day",
      date_labels = "%d %b",
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    labs(
      title = paste0(gas_label, ": all OTICE nodes average compared with CRDS average"),
      subtitle = paste(
        paste0(period_caption, " hourly averages across matched nodes."),
        "Yellow band marks the first 24 hours used for calibration."
      ),
      x = NULL,
      y = y_label,
      color = NULL,
      linetype = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold")
    )
}

make_overall_daily_plot <- function(df, gas_label, raw_col, ref_col, pred_col, y_label,
                                    baseline_start, baseline_end, period_caption) {
  style <- get_plot_style(gas_label)
  color_values <- setNames(
    c(style$ref_color, style$raw_color, style$cal_color),
    c("CRDS reference mean", style$raw_label, style$cal_label)
  )
  linetype_values <- setNames(
    c("solid", "dotted", "solid"),
    c("CRDS reference mean", style$raw_label, style$cal_label)
  )

  plot_df <- df |>
    select(
      date_day,
      raw_value = all_of(raw_col),
      ref_value = all_of(ref_col),
      calibrated_value = all_of(pred_col)
    ) |>
    pivot_longer(
      cols = c(ref_value, raw_value, calibrated_value),
      names_to = "source",
      values_to = "value"
    ) |>
    mutate(
      source = recode(
        source,
        ref_value = "CRDS reference mean",
        raw_value = style$raw_label,
        calibrated_value = style$cal_label
      )
    )

  ggplot(plot_df, aes(x = date_day, y = value, color = source, linetype = source)) +
    annotate(
      "rect",
      xmin = as.Date(baseline_start, tz = timezone_local),
      xmax = as.Date(baseline_end, tz = timezone_local),
      ymin = -Inf,
      ymax = Inf,
      alpha = 0.08,
      fill = "goldenrod"
    ) +
    geom_line(linewidth = 1.0, alpha = 0.9, na.rm = TRUE) +
    geom_point(size = 1.8, alpha = 0.75, na.rm = TRUE) +
    scale_color_manual(values = color_values) +
    scale_linetype_manual(values = linetype_values) +
    scale_x_date(
      date_breaks = "1 day",
      date_labels = "%d %b",
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    labs(
      title = paste0(gas_label, ": all OTICE nodes daily average compared with CRDS daily average"),
      subtitle = paste(
        paste0(period_caption, " daily averages derived from hourly all-node averages."),
        "Yellow band marks the first 24 hours used for calibration."
      ),
      x = NULL,
      y = y_label,
      color = NULL,
      linetype = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold")
    )
}

# -----------------------------------------------------------------------------
# 3. Campaign metadata with node-specific sensor IDs
# -----------------------------------------------------------------------------
campaign_reference_map <- bind_rows(
  build_campaign_reference(
    "C1", "2025-09-29 to 2025-10-06",
    "2025-09-29 00:00:00", "2025-10-06 00:00:00", "CRDS9",
    list(
      list(node = "2",  nh3_sensor_id = "9",  co2_sensor_id = "5",  crds_location = "24"),
      list(node = "6",  nh3_sensor_id = "4",  co2_sensor_id = "3",  crds_location = "24"),
      list(node = "7",  nh3_sensor_id = "20", co2_sensor_id = "4",  crds_location = "24"),
      list(node = "11", nh3_sensor_id = "13", co2_sensor_id = "10", crds_location = "24"),
      list(node = "17", nh3_sensor_id = "17", co2_sensor_id = "7",  crds_location = "24"),
      list(node = "18", nh3_sensor_id = "14", co2_sensor_id = "6",  crds_location = "24"),
      list(node = "14", nh3_sensor_id = "16", co2_sensor_id = "14", crds_location = "24")
    )
  ),
  build_campaign_reference(
    "C2", "2025-10-06 to 2025-10-21",
    "2025-10-06 00:00:00", "2025-10-21 00:00:00", "CRDS8",
    list(
      list(node = "2",  nh3_sensor_id = "9",  co2_sensor_id = "5",  crds_location = "30"),
      list(node = "6",  nh3_sensor_id = "4",  co2_sensor_id = "3",  crds_location = "30"),
      list(node = "7",  nh3_sensor_id = "20", co2_sensor_id = "4",  crds_location = "30"),
      list(node = "11", nh3_sensor_id = "13", co2_sensor_id = "10", crds_location = "30"),
      list(node = "17", nh3_sensor_id = "17", co2_sensor_id = "7",  crds_location = "30"),
      list(node = "18", nh3_sensor_id = "14", co2_sensor_id = "6",  crds_location = "30"),
      list(node = "14", nh3_sensor_id = "16", co2_sensor_id = "14", crds_location = "30"),
      list(node = "5",  nh3_sensor_id = "2",  co2_sensor_id = "15", crds_location = "30")
    )
  ),
  build_campaign_reference(
    "C3", "2025-10-21 to 2025-11-10",
    "2025-10-21 00:00:00", "2025-11-10 00:00:00", "CRDS8",
    list(
      list(node = "2",  nh3_sensor_id = "9",  co2_sensor_id = "5",  crds_location = "30"),
      list(node = "6",  nh3_sensor_id = "4",  co2_sensor_id = "3",  crds_location = "30"),
      list(node = "7",  nh3_sensor_id = "12", co2_sensor_id = "4",  crds_location = "30"),
      list(node = "11", nh3_sensor_id = "13", co2_sensor_id = "10", crds_location = "30"),
      list(node = "17", nh3_sensor_id = "17", co2_sensor_id = "7",  crds_location = "30"),
      list(node = "18", nh3_sensor_id = "14", co2_sensor_id = "6",  crds_location = "30"),
      list(node = "14", nh3_sensor_id = "16", co2_sensor_id = "14", crds_location = "30"),
      list(node = "5",  nh3_sensor_id = "2",  co2_sensor_id = "15", crds_location = "30")
    )
  ),
  build_campaign_reference(
    "C4", "2025-11-10 to 2025-11-12",
    "2025-11-10 00:00:00", "2025-11-12 00:00:00", "CRDS8",
    list(
      list(node = "2",  nh3_sensor_id = "9",  co2_sensor_id = "5",  crds_location = "30"),
      list(node = "6",  nh3_sensor_id = "4",  co2_sensor_id = "3",  crds_location = "30"),
      list(node = "7",  nh3_sensor_id = "12", co2_sensor_id = "4",  crds_location = "30"),
      list(node = "11", nh3_sensor_id = "13", co2_sensor_id = "10", crds_location = "30"),
      list(node = "17", nh3_sensor_id = "5",  co2_sensor_id = "12", crds_location = "30"),
      list(node = "18", nh3_sensor_id = "14", co2_sensor_id = "6",  crds_location = "30"),
      list(node = "14", nh3_sensor_id = "16", co2_sensor_id = "14", crds_location = "30"),
      list(node = "3",  nh3_sensor_id = "2",  co2_sensor_id = "15", crds_location = "30")
    )
  ),
  build_campaign_reference(
    "C5", "2025-11-12 to 2025-11-19",
    "2025-11-12 00:00:00", "2025-11-19 00:00:00", "CRDS8",
    list(
      list(node = "2",  nh3_sensor_id = "9",  co2_sensor_id = "5",  crds_location = "30"),
      list(node = "6",  nh3_sensor_id = "4",  co2_sensor_id = "3",  crds_location = "30"),
      list(node = "7",  nh3_sensor_id = "12", co2_sensor_id = "4",  crds_location = "30"),
      list(node = "11", nh3_sensor_id = "13", co2_sensor_id = "10", crds_location = "30"),
      list(node = "17", nh3_sensor_id = "5",  co2_sensor_id = "12", crds_location = "30"),
      list(node = "18", nh3_sensor_id = "14", co2_sensor_id = "6",  crds_location = "30"),
      list(node = "14", nh3_sensor_id = "16", co2_sensor_id = "13", crds_location = "30"),
      list(node = "3",  nh3_sensor_id = "2",  co2_sensor_id = "15", crds_location = "30")
    )
  ),
  build_campaign_reference(
    "C6", "2025-11-26 to 2025-12-01",
    "2025-11-26 00:00:00", "2025-12-01 00:00:00", "CRDS8",
    list(
      list(node = "2",  nh3_sensor_id = "9",  co2_sensor_id = "5",  crds_location = "36"),
      list(node = "6",  nh3_sensor_id = "4",  co2_sensor_id = "3",  crds_location = "30"),
      list(node = "7",  nh3_sensor_id = "12", co2_sensor_id = "4",  crds_location = "15"),
      list(node = "11", nh3_sensor_id = "13", co2_sensor_id = "10", crds_location = "45"),
      list(node = "17", nh3_sensor_id = "5",  co2_sensor_id = "12", crds_location = "42"),
      list(node = "18", nh3_sensor_id = "14", co2_sensor_id = "6",  crds_location = "12"),
      list(node = "3",  nh3_sensor_id = "2",  co2_sensor_id = "15", crds_location = "18")
    )
  ),
  build_campaign_reference(
    "C7", "2025-12-01 to 2025-12-09",
    "2025-12-01 00:00:00", "2025-12-09 00:00:00", "CRDS8",
    list(
      list(node = "2",  nh3_sensor_id = "9",  co2_sensor_id = "11", crds_location = "36"),
      list(node = "6",  nh3_sensor_id = "4",  co2_sensor_id = "2",  crds_location = "30"),
      list(node = "7",  nh3_sensor_id = "12", co2_sensor_id = "4",  crds_location = "15"),
      list(node = "11", nh3_sensor_id = "13", co2_sensor_id = "10", crds_location = "45"),
      list(node = "17", nh3_sensor_id = "5",  co2_sensor_id = "12", crds_location = "42"),
      list(node = "18", nh3_sensor_id = "14", co2_sensor_id = "14", crds_location = "12"),
      list(node = "3",  nh3_sensor_id = "2",  co2_sensor_id = "15", crds_location = "18")
    )
  ),
  build_campaign_reference(
    "C8", "2025-12-09 to 2025-12-31",
    "2025-12-09 00:00:00", "2026-01-01 00:00:00", "CRDS8",
    list(
      list(node = "2",  nh3_sensor_id = "9",  co2_sensor_id = "11", crds_location = "5"),
      list(node = "6",  nh3_sensor_id = "4",  co2_sensor_id = "8",  crds_location = "4"),
      list(node = "7",  nh3_sensor_id = "12", co2_sensor_id = "4",  crds_location = "3"),
      list(node = "11", nh3_sensor_id = "13", co2_sensor_id = "10", crds_location = "7"),
      list(node = "17", nh3_sensor_id = "5",  co2_sensor_id = "12", crds_location = "6"),
      list(node = "18", nh3_sensor_id = "16", co2_sensor_id = "13", crds_location = "1"),
      list(node = "3",  nh3_sensor_id = "2",  co2_sensor_id = "15", crds_location = "2")
    )
  )
)

write.csv(campaign_reference_map,
          file.path(tables_hourly_dir, "campaign_reference_map.csv"),
          row.names = FALSE)

c1_campaign_map <- campaign_reference_map |>
  filter(period_id == "C1")

c2_campaign_map <- campaign_reference_map |>
  filter(period_id == "C2")

december_campaign_map <- campaign_reference_map |>
  filter(period_id == "C8")

c1_start <- as.POSIXct("2025-09-29 00:00:00", tz = timezone_local)
c1_end <- as.POSIXct("2025-10-06 00:00:00", tz = timezone_local)
c2_start <- as.POSIXct("2025-10-06 00:00:00", tz = timezone_local)
c2_end <- as.POSIXct("2025-10-21 00:00:00", tz = timezone_local)
october_start <- c1_start
october_end <- c2_end

december_start <- as.POSIXct("2025-12-09 00:00:00", tz = timezone_local)
december_end <- as.POSIXct("2026-01-01 00:00:00", tz = timezone_local)
december_calibration_end <- december_start + days(1)

# -----------------------------------------------------------------------------
# 4. Read and standardize OTICE minute data
# -----------------------------------------------------------------------------
otice_files <- list.files(otice_dir, pattern = "^min_calibrated.*\\.csv$", full.names = TRUE)
otice_files <- otice_files[file.info(otice_files)$size > 0]

OTICE_dataset <- lapply(otice_files, read.csv, stringsAsFactors = FALSE) |>
  bind_rows() |>
  transmute(
    DATE.TIME = as.POSIXct(.data[["Datetime_Berlin"]], format = "%Y-%m-%d %H:%M:%S", tz = timezone_local),
    node      = clean_otice_location(.data[["Node"]]),
    analyzer  = as.character(.data[["Type"]]),
    NH3_raw   = suppressWarnings(as.numeric(.data[["NH3_ppm"]])),
    CO2_raw   = suppressWarnings(as.numeric(.data[["CO2.RAW"]])),
    NH3_cal   = suppressWarnings(as.numeric(.data[["NH3_ppm_barn"]])),
    CO2_cal   = suppressWarnings(as.numeric(.data[["CO2.AVG_barn"]]))
  ) |>
  filter(!is.na(DATE.TIME)) |>
  distinct(node, DATE.TIME, .keep_all = TRUE)

otice_hourly_node <- OTICE_dataset |>
  mutate(datetime_hour = floor_date(DATE.TIME, "hour")) |>
  group_by(node, datetime_hour) |>
  summarise(
    OTICE_NH3_raw = safe_mean(NH3_raw),
    OTICE_CO2_raw = safe_mean(CO2_raw),
    OTICE_NH3_cal = safe_mean(NH3_cal),
    OTICE_CO2_cal = safe_mean(CO2_cal),
    .groups = "drop"
  )

# -----------------------------------------------------------------------------
# 5. Read and standardize CRDS data
# -----------------------------------------------------------------------------
crds_files <- list.files(crds_dir, pattern = "\\.csv$", full.names = TRUE)
crds_files <- crds_files[file.info(crds_files)$size > 0]

CRDS_dataset <- lapply(crds_files, function(file_path) {
  df <- read.csv(file_path, stringsAsFactors = FALSE)
  df$location <- as.character(df$location)
  df$analyzer <- as.character(df$analyzer)
  df
}) |>
  bind_rows() |>
  transmute(
    DATE.TIME = as.POSIXct(.data[["DATE.TIME"]], format = "%Y-%m-%d %H:%M:%S", tz = timezone_local),
    crds_location = clean_crds_location(.data[["location"]]),
    analyzer      = toupper(trimws(as.character(.data[["analyzer"]]))),
    CRDS_NH3      = suppressWarnings(as.numeric(.data[["NH3"]])),
    CRDS_CO2      = suppressWarnings(as.numeric(.data[["CO2"]]))
  ) |>
  filter(!is.na(DATE.TIME)) |>
  distinct(analyzer, crds_location, DATE.TIME, .keep_all = TRUE)

crds_hourly_reference <- CRDS_dataset |>
  mutate(datetime_hour = floor_date(DATE.TIME, "hour")) |>
  group_by(analyzer, crds_location, datetime_hour) |>
  summarise(
    CRDS_NH3 = safe_mean(CRDS_NH3),
    CRDS_CO2 = safe_mean(CRDS_CO2),
    .groups = "drop"
  )

# -----------------------------------------------------------------------------
# 6. Build node-level matched hourly dataset
# -----------------------------------------------------------------------------
node_hourly_matched <- december_campaign_map |>
  inner_join(otice_hourly_node, by = "node", relationship = "many-to-many") |>
  filter(datetime_hour >= start_time, datetime_hour < end_time) |>
  inner_join(
    crds_hourly_reference,
    by = c("analyzer", "crds_location", "datetime_hour"),
    relationship = "many-to-many"
  ) |>
  arrange(node, datetime_hour)

node_hourly_matched <- node_hourly_matched |>
  mutate(
    nh3_episode_id = paste0("node", node, "_NH3_sensor", nh3_sensor_id),
    co2_episode_id = paste0("node", node, "_CO2_sensor", co2_sensor_id)
  )

write.csv(node_hourly_matched,
          file.path(tables_hourly_dir, "node_hourly_matched_december.csv"),
          row.names = FALSE)

c1_hourly_matched <- c1_campaign_map |>
  inner_join(otice_hourly_node, by = "node", relationship = "many-to-many") |>
  filter(datetime_hour >= start_time, datetime_hour < end_time) |>
  inner_join(
    crds_hourly_reference,
    by = c("analyzer", "crds_location", "datetime_hour"),
    relationship = "many-to-many"
  ) |>
  arrange(node, datetime_hour) |>
  mutate(
    nh3_episode_id = paste0("node", node, "_", period_id, "_NH3_sensor", nh3_sensor_id),
    co2_episode_id = paste0("node", node, "_", period_id, "_CO2_sensor", co2_sensor_id)
  )

c2_hourly_matched <- c2_campaign_map |>
  inner_join(otice_hourly_node, by = "node", relationship = "many-to-many") |>
  filter(datetime_hour >= start_time, datetime_hour < end_time) |>
  inner_join(
    crds_hourly_reference,
    by = c("analyzer", "crds_location", "datetime_hour"),
    relationship = "many-to-many"
  ) |>
  arrange(node, datetime_hour) |>
  mutate(
    nh3_episode_id = paste0("node", node, "_", period_id, "_NH3_sensor", nh3_sensor_id),
    co2_episode_id = paste0("node", node, "_", period_id, "_CO2_sensor", co2_sensor_id)
  )

october_hourly_matched <- bind_rows(c1_hourly_matched, c2_hourly_matched) |>
  arrange(node, datetime_hour)

write.csv(october_hourly_matched,
          file.path(tables_hourly_dir, "node_hourly_matched_october.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 7. First 24 hours of December calibration baseline
# -----------------------------------------------------------------------------
december_baseline <- node_hourly_matched |>
  filter(datetime_hour >= december_start,
         datetime_hour < december_calibration_end)

fit_calibration_model <- function(df, sensor_episode_col, raw_col, ref_col, gas_label,
                                  min_points = 8) {
  split_df <- split(df, df[[sensor_episode_col]])

  bind_rows(lapply(split_df, function(d) {
    ok <- !is.na(d[[raw_col]]) & !is.na(d[[ref_col]])
    if (sum(ok) < min_points) {
      return(NULL)
    }

    fit <- lm(d[[ref_col]][ok] ~ d[[raw_col]][ok])
    tibble(
      gas = gas_label,
      sensor_episode_id = d[[sensor_episode_col]][1],
      node = d$node[1],
      n_baseline = sum(ok),
      intercept = unname(coef(fit)[1]),
      slope = unname(coef(fit)[2]),
      baseline_R2 = summary(fit)$r.squared
    )
  }))
}

nh3_calibration_models <- fit_calibration_model(
  december_baseline,
  sensor_episode_col = "nh3_episode_id",
  raw_col = "OTICE_NH3_raw",
  ref_col = "CRDS_NH3",
  gas_label = "NH3"
)

co2_calibration_models <- fit_calibration_model(
  december_baseline,
  sensor_episode_col = "co2_episode_id",
  raw_col = "OTICE_CO2_raw",
  ref_col = "CRDS_CO2",
  gas_label = "CO2"
)

write.csv(nh3_calibration_models,
          file.path(tables_hourly_dir, "nh3_calibration_models_first24h_december.csv"),
          row.names = FALSE)

write.csv(co2_calibration_models,
          file.path(tables_hourly_dir, "co2_calibration_models_first24h_december.csv"),
          row.names = FALSE)

october_nh3_baseline <- bind_rows(
  build_first24h_baseline(c1_hourly_matched, "nh3_episode_id"),
  build_first24h_baseline(c2_hourly_matched, "nh3_episode_id")
)
october_co2_baseline <- bind_rows(
  build_first24h_baseline(c1_hourly_matched, "co2_episode_id"),
  build_first24h_baseline(c2_hourly_matched, "co2_episode_id")
)

october_nh3_calibration_models <- fit_calibration_model(
  october_nh3_baseline,
  sensor_episode_col = "nh3_episode_id",
  raw_col = "OTICE_NH3_raw",
  ref_col = "CRDS_NH3",
  gas_label = "NH3"
)

october_co2_calibration_models <- fit_calibration_model(
  october_co2_baseline,
  sensor_episode_col = "co2_episode_id",
  raw_col = "OTICE_CO2_raw",
  ref_col = "CRDS_CO2",
  gas_label = "CO2"
)

write.csv(october_nh3_calibration_models,
          file.path(tables_hourly_dir, "nh3_calibration_models_first24h_october.csv"),
          row.names = FALSE)

write.csv(october_co2_calibration_models,
          file.path(tables_hourly_dir, "co2_calibration_models_first24h_october.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 9. Apply first-24h December calibration across the December campaign
# -----------------------------------------------------------------------------
nh3_predictions <- node_hourly_matched |>
  left_join(
    nh3_calibration_models |>
      select(sensor_episode_id, nh3_intercept = intercept, nh3_slope = slope),
    by = c("nh3_episode_id" = "sensor_episode_id")
  ) |>
  mutate(
    NH3_predicted_from_december24h = nh3_intercept + nh3_slope * OTICE_NH3_raw,
    NH3_bias_after_december24h_cal = NH3_predicted_from_december24h - CRDS_NH3
  ) |>
  select(
    period_id, period_label, datetime_hour, node, nh3_sensor_id, nh3_episode_id,
    crds_location, OTICE_NH3_raw, OTICE_NH3_cal, CRDS_NH3,
    NH3_predicted_from_december24h, NH3_bias_after_december24h_cal
  )

co2_predictions <- node_hourly_matched |>
  left_join(
    co2_calibration_models |>
      select(sensor_episode_id, co2_intercept = intercept, co2_slope = slope),
    by = c("co2_episode_id" = "sensor_episode_id")
  ) |>
  mutate(
    CO2_predicted_from_december24h = co2_intercept + co2_slope * OTICE_CO2_raw,
    CO2_bias_after_december24h_cal = CO2_predicted_from_december24h - CRDS_CO2
  ) |>
  select(
    period_id, period_label, datetime_hour, node, co2_sensor_id, co2_episode_id,
    crds_location, OTICE_CO2_raw, OTICE_CO2_cal, CRDS_CO2,
    CO2_predicted_from_december24h, CO2_bias_after_december24h_cal
  )

write.csv(nh3_predictions,
          file.path(tables_hourly_dir, "nh3_predictions_after_first24h_december_calibration.csv"),
          row.names = FALSE)

write.csv(co2_predictions,
          file.path(tables_hourly_dir, "co2_predictions_after_first24h_december_calibration.csv"),
          row.names = FALSE)

october_nh3_predictions <- october_hourly_matched |>
  left_join(
    october_nh3_calibration_models |>
      select(sensor_episode_id, nh3_intercept = intercept, nh3_slope = slope),
    by = c("nh3_episode_id" = "sensor_episode_id")
  ) |>
  mutate(
    NH3_predicted_from_october24h = nh3_intercept + nh3_slope * OTICE_NH3_raw,
    NH3_bias_after_october24h_cal = NH3_predicted_from_october24h - CRDS_NH3
  ) |>
  select(
    period_id, period_label, datetime_hour, node, nh3_sensor_id, nh3_episode_id,
    crds_location, OTICE_NH3_raw, OTICE_NH3_cal, CRDS_NH3,
    NH3_predicted_from_october24h, NH3_bias_after_october24h_cal
  )

october_co2_predictions <- october_hourly_matched |>
  left_join(
    october_co2_calibration_models |>
      select(sensor_episode_id, co2_intercept = intercept, co2_slope = slope),
    by = c("co2_episode_id" = "sensor_episode_id")
  ) |>
  mutate(
    CO2_predicted_from_october24h = co2_intercept + co2_slope * OTICE_CO2_raw,
    CO2_bias_after_october24h_cal = CO2_predicted_from_october24h - CRDS_CO2
  ) |>
  select(
    period_id, period_label, datetime_hour, node, co2_sensor_id, co2_episode_id,
    crds_location, OTICE_CO2_raw, OTICE_CO2_cal, CRDS_CO2,
    CO2_predicted_from_october24h, CO2_bias_after_october24h_cal
  )

write.csv(october_nh3_predictions,
          file.path(tables_hourly_dir, "nh3_predictions_after_first24h_october_calibration.csv"),
          row.names = FALSE)

write.csv(october_co2_predictions,
          file.path(tables_hourly_dir, "co2_predictions_after_first24h_october_calibration.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 10. December-only node-level visual comparison
# -----------------------------------------------------------------------------
december_nh3_comparison <- nh3_predictions |>
  filter(period_id == "C8",
         datetime_hour >= december_start,
         datetime_hour < december_end) |>
  select(period_id, period_label, datetime_hour, node, nh3_sensor_id, crds_location,
         OTICE_NH3_raw, OTICE_NH3_cal, CRDS_NH3, NH3_predicted_from_december24h)

december_co2_comparison <- co2_predictions |>
  filter(period_id == "C8",
         datetime_hour >= december_start,
         datetime_hour < december_end) |>
  select(period_id, period_label, datetime_hour, node, co2_sensor_id, crds_location,
         OTICE_CO2_raw, OTICE_CO2_cal, CRDS_CO2, CO2_predicted_from_december24h)

write.csv(december_nh3_comparison,
          file.path(tables_hourly_dir, "december_nh3_node_vs_crds.csv"),
          row.names = FALSE)

write.csv(december_co2_comparison,
          file.path(tables_hourly_dir, "december_co2_node_vs_crds.csv"),
          row.names = FALSE)

december_nh3_daily <- december_nh3_comparison |>
  mutate(date_day = as.Date(datetime_hour, tz = timezone_local)) |>
  group_by(period_id, period_label, date_day, node, nh3_sensor_id, crds_location) |>
  summarise(
    OTICE_NH3_raw = safe_mean(OTICE_NH3_raw),
    OTICE_NH3_cal = safe_mean(OTICE_NH3_cal),
    CRDS_NH3 = safe_mean(CRDS_NH3),
    NH3_predicted_from_december24h = safe_mean(NH3_predicted_from_december24h),
    n_hourly_pairs = sum(!is.na(OTICE_NH3_raw) & !is.na(CRDS_NH3)),
    .groups = "drop"
  )

december_co2_daily <- december_co2_comparison |>
  mutate(date_day = as.Date(datetime_hour, tz = timezone_local)) |>
  group_by(period_id, period_label, date_day, node, co2_sensor_id, crds_location) |>
  summarise(
    OTICE_CO2_raw = safe_mean(OTICE_CO2_raw),
    OTICE_CO2_cal = safe_mean(OTICE_CO2_cal),
    CRDS_CO2 = safe_mean(CRDS_CO2),
    CO2_predicted_from_december24h = safe_mean(CO2_predicted_from_december24h),
    n_hourly_pairs = sum(!is.na(OTICE_CO2_raw) & !is.na(CRDS_CO2)),
    .groups = "drop"
  )

write.csv(december_nh3_daily,
          file.path(tables_daily_dir, "december_nh3_node_vs_crds_daily.csv"),
          row.names = FALSE)

write.csv(december_co2_daily,
          file.path(tables_daily_dir, "december_co2_node_vs_crds_daily.csv"),
          row.names = FALSE)

december_nh3_all_nodes_hourly <- december_nh3_comparison |>
  group_by(period_id, period_label, datetime_hour) |>
  summarise(
    OTICE_NH3_raw = safe_mean(OTICE_NH3_raw),
    OTICE_NH3_cal = safe_mean(OTICE_NH3_cal),
    CRDS_NH3 = safe_mean(CRDS_NH3),
    NH3_predicted_from_december24h = safe_mean(NH3_predicted_from_december24h),
    n_nodes = sum(!is.na(OTICE_NH3_raw) & !is.na(CRDS_NH3)),
    .groups = "drop"
  )

december_co2_all_nodes_hourly <- december_co2_comparison |>
  group_by(period_id, period_label, datetime_hour) |>
  summarise(
    OTICE_CO2_raw = safe_mean(OTICE_CO2_raw),
    OTICE_CO2_cal = safe_mean(OTICE_CO2_cal),
    CRDS_CO2 = safe_mean(CRDS_CO2),
    CO2_predicted_from_december24h = safe_mean(CO2_predicted_from_december24h),
    n_nodes = sum(!is.na(OTICE_CO2_raw) & !is.na(CRDS_CO2)),
    .groups = "drop"
  )

write.csv(december_nh3_all_nodes_hourly,
          file.path(tables_hourly_dir, "december_nh3_all_nodes_vs_crds_hourly.csv"),
          row.names = FALSE)

write.csv(december_co2_all_nodes_hourly,
          file.path(tables_hourly_dir, "december_co2_all_nodes_vs_crds_hourly.csv"),
          row.names = FALSE)

december_nh3_all_nodes_daily <- december_nh3_all_nodes_hourly |>
  mutate(date_day = as.Date(datetime_hour, tz = timezone_local)) |>
  group_by(period_id, period_label, date_day) |>
  summarise(
    OTICE_NH3_raw = safe_mean(OTICE_NH3_raw),
    OTICE_NH3_cal = safe_mean(OTICE_NH3_cal),
    CRDS_NH3 = safe_mean(CRDS_NH3),
    NH3_predicted_from_december24h = safe_mean(NH3_predicted_from_december24h),
    n_hourly_points = sum(!is.na(OTICE_NH3_raw) & !is.na(CRDS_NH3)),
    mean_nodes_per_hour = safe_mean(n_nodes),
    .groups = "drop"
  )

december_co2_all_nodes_daily <- december_co2_all_nodes_hourly |>
  mutate(date_day = as.Date(datetime_hour, tz = timezone_local)) |>
  group_by(period_id, period_label, date_day) |>
  summarise(
    OTICE_CO2_raw = safe_mean(OTICE_CO2_raw),
    OTICE_CO2_cal = safe_mean(OTICE_CO2_cal),
    CRDS_CO2 = safe_mean(CRDS_CO2),
    CO2_predicted_from_december24h = safe_mean(CO2_predicted_from_december24h),
    n_hourly_points = sum(!is.na(OTICE_CO2_raw) & !is.na(CRDS_CO2)),
    mean_nodes_per_hour = safe_mean(n_nodes),
    .groups = "drop"
  )

write.csv(december_nh3_all_nodes_daily,
          file.path(tables_daily_dir, "december_nh3_all_nodes_vs_crds_daily.csv"),
          row.names = FALSE)

write.csv(december_co2_all_nodes_daily,
          file.path(tables_daily_dir, "december_co2_all_nodes_vs_crds_daily.csv"),
          row.names = FALSE)

october_nh3_comparison <- october_nh3_predictions |>
  filter(period_id %in% c("C1", "C2"),
         datetime_hour >= october_start,
         datetime_hour < october_end) |>
  select(period_id, period_label, datetime_hour, node, nh3_sensor_id, crds_location,
         OTICE_NH3_raw, OTICE_NH3_cal, CRDS_NH3, NH3_predicted_from_october24h)

october_co2_comparison <- october_co2_predictions |>
  filter(period_id %in% c("C1", "C2"),
         datetime_hour >= october_start,
         datetime_hour < october_end) |>
  select(period_id, period_label, datetime_hour, node, co2_sensor_id, crds_location,
         OTICE_CO2_raw, OTICE_CO2_cal, CRDS_CO2, CO2_predicted_from_october24h)

write.csv(october_nh3_comparison,
          file.path(tables_hourly_dir, "october_nh3_node_vs_crds.csv"),
          row.names = FALSE)

write.csv(october_co2_comparison,
          file.path(tables_hourly_dir, "october_co2_node_vs_crds.csv"),
          row.names = FALSE)

october_nh3_daily <- october_nh3_comparison |>
  mutate(date_day = as.Date(datetime_hour, tz = timezone_local)) |>
  group_by(period_id, period_label, date_day, node, nh3_sensor_id, crds_location) |>
  summarise(
    OTICE_NH3_raw = safe_mean(OTICE_NH3_raw),
    OTICE_NH3_cal = safe_mean(OTICE_NH3_cal),
    CRDS_NH3 = safe_mean(CRDS_NH3),
    NH3_predicted_from_october24h = safe_mean(NH3_predicted_from_october24h),
    n_hourly_pairs = sum(!is.na(OTICE_NH3_raw) & !is.na(CRDS_NH3)),
    .groups = "drop"
  )

october_co2_daily <- october_co2_comparison |>
  mutate(date_day = as.Date(datetime_hour, tz = timezone_local)) |>
  group_by(period_id, period_label, date_day, node, co2_sensor_id, crds_location) |>
  summarise(
    OTICE_CO2_raw = safe_mean(OTICE_CO2_raw),
    OTICE_CO2_cal = safe_mean(OTICE_CO2_cal),
    CRDS_CO2 = safe_mean(CRDS_CO2),
    CO2_predicted_from_october24h = safe_mean(CO2_predicted_from_october24h),
    n_hourly_pairs = sum(!is.na(OTICE_CO2_raw) & !is.na(CRDS_CO2)),
    .groups = "drop"
  )

write.csv(october_nh3_daily,
          file.path(tables_daily_dir, "october_nh3_node_vs_crds_daily.csv"),
          row.names = FALSE)

write.csv(october_co2_daily,
          file.path(tables_daily_dir, "october_co2_node_vs_crds_daily.csv"),
          row.names = FALSE)

october_nh3_all_nodes_hourly <- october_nh3_comparison |>
  group_by(period_id, period_label, datetime_hour) |>
  summarise(
    OTICE_NH3_raw = safe_mean(OTICE_NH3_raw),
    OTICE_NH3_cal = safe_mean(OTICE_NH3_cal),
    CRDS_NH3 = safe_mean(CRDS_NH3),
    NH3_predicted_from_october24h = safe_mean(NH3_predicted_from_october24h),
    n_nodes = sum(!is.na(OTICE_NH3_raw) & !is.na(CRDS_NH3)),
    .groups = "drop"
  )

october_co2_all_nodes_hourly <- october_co2_comparison |>
  group_by(period_id, period_label, datetime_hour) |>
  summarise(
    OTICE_CO2_raw = safe_mean(OTICE_CO2_raw),
    OTICE_CO2_cal = safe_mean(OTICE_CO2_cal),
    CRDS_CO2 = safe_mean(CRDS_CO2),
    CO2_predicted_from_october24h = safe_mean(CO2_predicted_from_october24h),
    n_nodes = sum(!is.na(OTICE_CO2_raw) & !is.na(CRDS_CO2)),
    .groups = "drop"
  )

write.csv(october_nh3_all_nodes_hourly,
          file.path(tables_hourly_dir, "october_nh3_all_nodes_vs_crds_hourly.csv"),
          row.names = FALSE)

write.csv(october_co2_all_nodes_hourly,
          file.path(tables_hourly_dir, "october_co2_all_nodes_vs_crds_hourly.csv"),
          row.names = FALSE)

october_nh3_all_nodes_daily <- october_nh3_all_nodes_hourly |>
  mutate(date_day = as.Date(datetime_hour, tz = timezone_local)) |>
  group_by(period_id, period_label, date_day) |>
  summarise(
    OTICE_NH3_raw = safe_mean(OTICE_NH3_raw),
    OTICE_NH3_cal = safe_mean(OTICE_NH3_cal),
    CRDS_NH3 = safe_mean(CRDS_NH3),
    NH3_predicted_from_october24h = safe_mean(NH3_predicted_from_october24h),
    n_hourly_points = sum(!is.na(OTICE_NH3_raw) & !is.na(CRDS_NH3)),
    mean_nodes_per_hour = safe_mean(n_nodes),
    .groups = "drop"
  )

october_co2_all_nodes_daily <- october_co2_all_nodes_hourly |>
  mutate(date_day = as.Date(datetime_hour, tz = timezone_local)) |>
  group_by(period_id, period_label, date_day) |>
  summarise(
    OTICE_CO2_raw = safe_mean(OTICE_CO2_raw),
    OTICE_CO2_cal = safe_mean(OTICE_CO2_cal),
    CRDS_CO2 = safe_mean(CRDS_CO2),
    CO2_predicted_from_october24h = safe_mean(CO2_predicted_from_october24h),
    n_hourly_points = sum(!is.na(OTICE_CO2_raw) & !is.na(CRDS_CO2)),
    mean_nodes_per_hour = safe_mean(n_nodes),
    .groups = "drop"
  )

write.csv(october_nh3_all_nodes_daily,
          file.path(tables_daily_dir, "october_nh3_all_nodes_vs_crds_daily.csv"),
          row.names = FALSE)

write.csv(october_co2_all_nodes_daily,
          file.path(tables_daily_dir, "october_co2_all_nodes_vs_crds_daily.csv"),
          row.names = FALSE)

if (nrow(december_nh3_comparison) > 0) {
  for (node_id in sort(unique(december_nh3_comparison$node))) {
    save_plot(
      make_node_timeseries_plot(
        december_nh3_comparison,
        node_id = node_id,
        gas_label = "NH3",
        raw_col = "OTICE_NH3_raw",
        ref_col = "CRDS_NH3",
        pred_col = "NH3_predicted_from_december24h",
        y_label = "NH3 (ppm)",
        baseline_start = december_start,
        baseline_end = december_calibration_end,
        period_caption = "December 2025"
      ),
      file.path(plots_hourly_dir, paste0("December_NH3_node_", node_id, ".png")),
      width = 14,
      height = 8
    )
  }
}

if (nrow(december_nh3_daily) > 0) {
  for (node_id in sort(unique(december_nh3_daily$node))) {
    save_plot(
      make_node_daily_plot(
        december_nh3_daily,
        node_id = node_id,
        gas_label = "NH3",
        raw_col = "OTICE_NH3_raw",
        ref_col = "CRDS_NH3",
        pred_col = "NH3_predicted_from_december24h",
        y_label = "NH3 (ppm)",
        baseline_start = december_start,
        baseline_end = december_calibration_end,
        period_caption = "December 2025"
      ),
      file.path(plots_daily_dir, paste0("December_NH3_node_", node_id, "_daily.png")),
      width = 14,
      height = 8
    )
  }
}

if (nrow(december_nh3_all_nodes_hourly) > 0) {
  save_plot(
    make_overall_timeseries_plot(
      december_nh3_all_nodes_hourly,
      gas_label = "NH3",
      raw_col = "OTICE_NH3_raw",
      ref_col = "CRDS_NH3",
      pred_col = "NH3_predicted_from_december24h",
      y_label = "NH3 (ppm)",
      baseline_start = december_start,
      baseline_end = december_calibration_end,
      period_caption = "December 2025"
    ),
    file.path(plots_hourly_dir, "December_NH3_all_nodes_hourly.png"),
    width = 14,
    height = 8
  )
}

if (nrow(december_nh3_all_nodes_daily) > 0) {
  save_plot(
    make_overall_daily_plot(
      december_nh3_all_nodes_daily,
      gas_label = "NH3",
      raw_col = "OTICE_NH3_raw",
      ref_col = "CRDS_NH3",
      pred_col = "NH3_predicted_from_december24h",
      y_label = "NH3 (ppm)",
      baseline_start = december_start,
      baseline_end = december_calibration_end,
      period_caption = "December 2025"
    ),
    file.path(plots_daily_dir, "December_NH3_all_nodes_daily.png"),
    width = 14,
    height = 8
  )
}

if (nrow(december_co2_comparison) > 0) {
  for (node_id in sort(unique(december_co2_comparison$node))) {
    save_plot(
      make_node_timeseries_plot(
        december_co2_comparison,
        node_id = node_id,
        gas_label = "CO2",
        raw_col = "OTICE_CO2_raw",
        ref_col = "CRDS_CO2",
        pred_col = "CO2_predicted_from_december24h",
        y_label = "CO2 (ppm)",
        baseline_start = december_start,
        baseline_end = december_calibration_end,
        period_caption = "December 2025"
      ),
      file.path(plots_hourly_dir, paste0("December_CO2_node_", node_id, ".png")),
      width = 14,
      height = 8
    )
  }
}

if (nrow(december_co2_daily) > 0) {
  for (node_id in sort(unique(december_co2_daily$node))) {
    save_plot(
      make_node_daily_plot(
        december_co2_daily,
        node_id = node_id,
        gas_label = "CO2",
        raw_col = "OTICE_CO2_raw",
        ref_col = "CRDS_CO2",
        pred_col = "CO2_predicted_from_december24h",
        y_label = "CO2 (ppm)",
        baseline_start = december_start,
        baseline_end = december_calibration_end,
        period_caption = "December 2025"
      ),
      file.path(plots_daily_dir, paste0("December_CO2_node_", node_id, "_daily.png")),
      width = 14,
      height = 8
    )
  }
}

if (nrow(december_co2_all_nodes_hourly) > 0) {
  save_plot(
    make_overall_timeseries_plot(
      december_co2_all_nodes_hourly,
      gas_label = "CO2",
      raw_col = "OTICE_CO2_raw",
      ref_col = "CRDS_CO2",
      pred_col = "CO2_predicted_from_december24h",
      y_label = "CO2 (ppm)",
      baseline_start = december_start,
      baseline_end = december_calibration_end,
      period_caption = "December 2025"
    ),
    file.path(plots_hourly_dir, "December_CO2_all_nodes_hourly.png"),
    width = 14,
    height = 8
  )
}

if (nrow(december_co2_all_nodes_daily) > 0) {
  save_plot(
    make_overall_daily_plot(
      december_co2_all_nodes_daily,
      gas_label = "CO2",
      raw_col = "OTICE_CO2_raw",
      ref_col = "CRDS_CO2",
      pred_col = "CO2_predicted_from_december24h",
      y_label = "CO2 (ppm)",
      baseline_start = december_start,
      baseline_end = december_calibration_end,
      period_caption = "December 2025"
    ),
    file.path(plots_daily_dir, "December_CO2_all_nodes_daily.png"),
    width = 14,
    height = 8
  )
}

if (nrow(october_nh3_comparison) > 0) {
  for (node_id in sort(unique(october_nh3_comparison$node))) {
    node_baseline_start <- min(
      october_nh3_comparison$datetime_hour[october_nh3_comparison$node == node_id],
      na.rm = TRUE
    )
    save_plot(
      make_node_timeseries_plot(
        october_nh3_comparison,
        node_id = node_id,
        gas_label = "NH3",
        raw_col = "OTICE_NH3_raw",
        ref_col = "CRDS_NH3",
        pred_col = "NH3_predicted_from_october24h",
        y_label = "NH3 (ppm)",
        baseline_start = node_baseline_start,
        baseline_end = node_baseline_start + days(1),
        period_caption = "September-October 2025"
      ),
      file.path(plots_hourly_dir, paste0("October_NH3_node_", node_id, ".png")),
      width = 14,
      height = 8
    )
  }
}

if (nrow(october_nh3_daily) > 0) {
  for (node_id in sort(unique(october_nh3_daily$node))) {
    node_baseline_start <- min(
      october_nh3_comparison$datetime_hour[october_nh3_comparison$node == node_id],
      na.rm = TRUE
    )
    save_plot(
      make_node_daily_plot(
        october_nh3_daily,
        node_id = node_id,
        gas_label = "NH3",
        raw_col = "OTICE_NH3_raw",
        ref_col = "CRDS_NH3",
        pred_col = "NH3_predicted_from_october24h",
        y_label = "NH3 (ppm)",
        baseline_start = node_baseline_start,
        baseline_end = node_baseline_start + days(1),
        period_caption = "September-October 2025"
      ),
      file.path(plots_daily_dir, paste0("October_NH3_node_", node_id, "_daily.png")),
      width = 14,
      height = 8
    )
  }
}

if (nrow(october_nh3_all_nodes_hourly) > 0) {
  save_plot(
    make_overall_timeseries_plot(
      october_nh3_all_nodes_hourly,
      gas_label = "NH3",
      raw_col = "OTICE_NH3_raw",
      ref_col = "CRDS_NH3",
      pred_col = "NH3_predicted_from_october24h",
      y_label = "NH3 (ppm)",
      baseline_start = october_start,
      baseline_end = october_start + days(1),
      period_caption = "September-October 2025"
    ),
    file.path(plots_hourly_dir, "October_NH3_all_nodes_hourly.png"),
    width = 14,
    height = 8
  )
}

if (nrow(october_nh3_all_nodes_daily) > 0) {
  save_plot(
    make_overall_daily_plot(
      october_nh3_all_nodes_daily,
      gas_label = "NH3",
      raw_col = "OTICE_NH3_raw",
      ref_col = "CRDS_NH3",
      pred_col = "NH3_predicted_from_october24h",
      y_label = "NH3 (ppm)",
      baseline_start = october_start,
      baseline_end = october_start + days(1),
      period_caption = "September-October 2025"
    ),
    file.path(plots_daily_dir, "October_NH3_all_nodes_daily.png"),
    width = 14,
    height = 8
  )
}

if (nrow(october_co2_comparison) > 0) {
  for (node_id in sort(unique(october_co2_comparison$node))) {
    node_baseline_start <- min(
      october_co2_comparison$datetime_hour[october_co2_comparison$node == node_id],
      na.rm = TRUE
    )
    save_plot(
      make_node_timeseries_plot(
        october_co2_comparison,
        node_id = node_id,
        gas_label = "CO2",
        raw_col = "OTICE_CO2_raw",
        ref_col = "CRDS_CO2",
        pred_col = "CO2_predicted_from_october24h",
        y_label = "CO2 (ppm)",
        baseline_start = node_baseline_start,
        baseline_end = node_baseline_start + days(1),
        period_caption = "September-October 2025"
      ),
      file.path(plots_hourly_dir, paste0("October_CO2_node_", node_id, ".png")),
      width = 14,
      height = 8
    )
  }
}

if (nrow(october_co2_daily) > 0) {
  for (node_id in sort(unique(october_co2_daily$node))) {
    node_baseline_start <- min(
      october_co2_comparison$datetime_hour[october_co2_comparison$node == node_id],
      na.rm = TRUE
    )
    save_plot(
      make_node_daily_plot(
        october_co2_daily,
        node_id = node_id,
        gas_label = "CO2",
        raw_col = "OTICE_CO2_raw",
        ref_col = "CRDS_CO2",
        pred_col = "CO2_predicted_from_october24h",
        y_label = "CO2 (ppm)",
        baseline_start = node_baseline_start,
        baseline_end = node_baseline_start + days(1),
        period_caption = "September-October 2025"
      ),
      file.path(plots_daily_dir, paste0("October_CO2_node_", node_id, "_daily.png")),
      width = 14,
      height = 8
    )
  }
}

if (nrow(october_co2_all_nodes_hourly) > 0) {
  save_plot(
    make_overall_timeseries_plot(
      october_co2_all_nodes_hourly,
      gas_label = "CO2",
      raw_col = "OTICE_CO2_raw",
      ref_col = "CRDS_CO2",
      pred_col = "CO2_predicted_from_october24h",
      y_label = "CO2 (ppm)",
      baseline_start = october_start,
      baseline_end = october_start + days(1),
      period_caption = "September-October 2025"
    ),
    file.path(plots_hourly_dir, "October_CO2_all_nodes_hourly.png"),
    width = 14,
    height = 8
  )
}

if (nrow(october_co2_all_nodes_daily) > 0) {
  save_plot(
    make_overall_daily_plot(
      october_co2_all_nodes_daily,
      gas_label = "CO2",
      raw_col = "OTICE_CO2_raw",
      ref_col = "CRDS_CO2",
      pred_col = "CO2_predicted_from_october24h",
      y_label = "CO2 (ppm)",
      baseline_start = october_start,
      baseline_end = october_start + days(1),
      period_caption = "September-October 2025"
    ),
    file.path(plots_daily_dir, "October_CO2_all_nodes_daily.png"),
    width = 14,
    height = 8
  )
}

# -----------------------------------------------------------------------------
# 11. Console summary
# -----------------------------------------------------------------------------
cat("\n============================================================\n")
cat("OTICE December node-versus-CRDS comparison\n")
cat("============================================================\n")
cat("Focus campaign: 2025-12-09 to 2025-12-31 (C8)\n")
cat("Calibration baseline: first 24 hours of 2025-12-09\n")
cat("September-October campaign: 2025-09-29 to 2025-10-21 (C1 + C2)\n")
cat("September-October calibration baseline: first 24 hours from each node's initial matched timestamp\n")
cat("Matched hourly rows in October:", nrow(october_hourly_matched), "\n")
cat("October NH3 rows:", nrow(october_nh3_comparison), "\n")
cat("October CO2 rows:", nrow(october_co2_comparison), "\n")
cat("October NH3 daily rows:", nrow(october_nh3_daily), "\n")
cat("October CO2 daily rows:", nrow(october_co2_daily), "\n")
cat("October NH3 all-node hourly rows:", nrow(october_nh3_all_nodes_hourly), "\n")
cat("October CO2 all-node hourly rows:", nrow(october_co2_all_nodes_hourly), "\n")
cat("October NH3 all-node daily rows:", nrow(october_nh3_all_nodes_daily), "\n")
cat("October CO2 all-node daily rows:", nrow(october_co2_all_nodes_daily), "\n")
cat("Matched hourly rows in December:", nrow(node_hourly_matched), "\n")
cat("December NH3 rows:", nrow(december_nh3_comparison), "\n")
cat("December CO2 rows:", nrow(december_co2_comparison), "\n")
cat("December NH3 daily rows:", nrow(december_nh3_daily), "\n")
cat("December CO2 daily rows:", nrow(december_co2_daily), "\n")
cat("December NH3 all-node hourly rows:", nrow(december_nh3_all_nodes_hourly), "\n")
cat("December CO2 all-node hourly rows:", nrow(december_co2_all_nodes_hourly), "\n")
cat("December NH3 all-node daily rows:", nrow(december_nh3_all_nodes_daily), "\n")
cat("December CO2 all-node daily rows:", nrow(december_co2_all_nodes_daily), "\n")
cat("\nFiles written to:", out_dir, "\n")
cat("  tables/hourly/campaign_reference_map.csv\n")
cat("  tables/hourly/node_hourly_matched_december.csv\n")
cat("  tables/hourly/node_hourly_matched_october.csv\n")
cat("  tables/hourly/nh3_calibration_models_first24h_december.csv\n")
cat("  tables/hourly/co2_calibration_models_first24h_december.csv\n")
cat("  tables/hourly/nh3_calibration_models_first24h_october.csv\n")
cat("  tables/hourly/co2_calibration_models_first24h_october.csv\n")
cat("  tables/hourly/nh3_predictions_after_first24h_december_calibration.csv\n")
cat("  tables/hourly/co2_predictions_after_first24h_december_calibration.csv\n")
cat("  tables/hourly/nh3_predictions_after_first24h_october_calibration.csv\n")
cat("  tables/hourly/co2_predictions_after_first24h_october_calibration.csv\n")
cat("  tables/hourly/october_nh3_node_vs_crds.csv\n")
cat("  tables/hourly/october_co2_node_vs_crds.csv\n")
cat("  tables/hourly/october_nh3_all_nodes_vs_crds_hourly.csv\n")
cat("  tables/hourly/october_co2_all_nodes_vs_crds_hourly.csv\n")
cat("  tables/hourly/december_nh3_node_vs_crds.csv\n")
cat("  tables/hourly/december_co2_node_vs_crds.csv\n")
cat("  tables/hourly/december_nh3_all_nodes_vs_crds_hourly.csv\n")
cat("  tables/hourly/december_co2_all_nodes_vs_crds_hourly.csv\n")
cat("  tables/daily/october_nh3_node_vs_crds_daily.csv\n")
cat("  tables/daily/october_co2_node_vs_crds_daily.csv\n")
cat("  tables/daily/october_nh3_all_nodes_vs_crds_daily.csv\n")
cat("  tables/daily/october_co2_all_nodes_vs_crds_daily.csv\n")
cat("  tables/daily/december_nh3_node_vs_crds_daily.csv\n")
cat("  tables/daily/december_co2_node_vs_crds_daily.csv\n")
cat("  tables/daily/december_nh3_all_nodes_vs_crds_daily.csv\n")
cat("  tables/daily/december_co2_all_nodes_vs_crds_daily.csv\n")
cat("  plots/hourly/*.png\n")
cat("  plots/daily/*.png\n")
