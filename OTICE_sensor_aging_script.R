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

otice_dir <- "processed_data/OTICE_processed"
crds_dir  <- "processed_data/CRDS8_processed"
out_dir   <- file.path("output", "OTICE_sensor_aging_december")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "plots"), showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), showWarnings = FALSE)

unlink(list.files(file.path(out_dir, "plots"), full.names = TRUE), force = TRUE)
unlink(list.files(file.path(out_dir, "tables"), full.names = TRUE), force = TRUE)

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

make_timeseries_plot <- function(df, gas_label, raw_col, ref_col, pred_col, y_label,
                                 baseline_start, baseline_end) {

  plot_df <- df |>
    select(node, period_id, period_label, datetime_hour,
           raw_value = all_of(raw_col),
           ref_value = all_of(ref_col),
           calibrated_value = all_of(pred_col)) |>
    pivot_longer(
      cols = c(ref_value, raw_value, calibrated_value),
      names_to = "source",
      values_to = "value"
    ) |>
    mutate(
      source = recode(
        source,
        ref_value = "CRDS reference",
        raw_value = "OTICE raw",
        calibrated_value = "OTICE after first-24h December calibration"
      )
    )

  ggplot(plot_df, aes(x = datetime_hour, y = value, color = source)) +
    annotate(
      "rect",
      xmin = baseline_start,
      xmax = baseline_end,
      ymin = -Inf,
      ymax = Inf,
      alpha = 0.08,
      fill = "goldenrod"
    ) +
    geom_line(linewidth = 0.6, na.rm = TRUE) +
    facet_wrap(~ node, scales = "free_y", ncol = 2) +
    scale_color_manual(
      values = c(
        "CRDS reference" = "#1a1a2e",
        "OTICE raw" = "#e76f51",
        "OTICE after first-24h December calibration" = "#2a9d8f"
      )
    ) +
    scale_x_datetime(
      date_breaks = "1 day",
      date_labels = "%d %b",
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    labs(
      title = paste0(gas_label, ": December node comparison through time"),
      subtitle = paste(
        "Yellow band = first 24 hours of the December campaign used for calibration.",
        "Compare OTICE raw and December-calibrated lines against CRDS over time."
      ),
      x = NULL,
      y = y_label,
      color = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold")
    )
}

make_node_timeseries_plot <- function(df, node_id, gas_label, raw_col, ref_col, pred_col, y_label,
                                      baseline_start, baseline_end) {
  node_df <- df |>
    filter(node == node_id) |>
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
        ref_value = "CRDS reference",
        raw_value = "OTICE raw",
        calibrated_value = "OTICE after first-24h December calibration"
      )
    )

  ggplot(node_df, aes(x = datetime_hour, y = value, color = source)) +
    annotate(
      "rect",
      xmin = baseline_start,
      xmax = baseline_end,
      ymin = -Inf,
      ymax = Inf,
      alpha = 0.08,
      fill = "goldenrod"
    ) +
    geom_line(linewidth = 0.7, na.rm = TRUE) +
    geom_point(size = 1.1, alpha = 0.65, na.rm = TRUE) +
    scale_color_manual(
      values = c(
        "CRDS reference" = "#1a1a2e",
        "OTICE raw" = "#e76f51",
        "OTICE after first-24h December calibration" = "#2a9d8f"
      )
    ) +
    scale_x_datetime(
      date_breaks = "1 day",
      date_labels = "%d %b",
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    labs(
      title = paste0(gas_label, " - December comparison for OTICE node ", node_id),
      subtitle = "Hourly OTICE raw and first-24h-December calibrated values versus mapped CRDS reference",
      x = NULL,
      y = y_label,
      color = NULL
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
      list(node = "2",  nh3_sensor_id = "9",  co2_sensor_id = "5",  crds_location = c("24", "27")),
      list(node = "6",  nh3_sensor_id = "4",  co2_sensor_id = "3",  crds_location = c("24", "27")),
      list(node = "7",  nh3_sensor_id = "20", co2_sensor_id = "4",  crds_location = c("24", "27")),
      list(node = "11", nh3_sensor_id = "13", co2_sensor_id = "10", crds_location = c("24", "27")),
      list(node = "17", nh3_sensor_id = "17", co2_sensor_id = "7",  crds_location = c("24", "27")),
      list(node = "18", nh3_sensor_id = "14", co2_sensor_id = "6",  crds_location = c("24", "27")),
      list(node = "14", nh3_sensor_id = "16", co2_sensor_id = "14", crds_location = c("24", "27"))
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
          file.path(out_dir, "tables", "campaign_reference_map.csv"),
          row.names = FALSE)

december_campaign_map <- campaign_reference_map |>
  filter(period_id == "C8")

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
    DATE.TIME = as.POSIXct(Datetime_Berlin, format = "%Y-%m-%d %H:%M:%S", tz = timezone_local),
    node      = clean_otice_location(Node),
    analyzer  = as.character(Type),
    NH3_raw   = suppressWarnings(as.numeric(NH3_ppm)),
    CO2_raw   = suppressWarnings(as.numeric(CO2.RAW)),
    NH3_cal   = suppressWarnings(as.numeric(NH3_ppm_barn)),
    CO2_cal   = suppressWarnings(as.numeric(CO2.AVG_barn))
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
    DATE.TIME = as.POSIXct(DATE.TIME, format = "%Y-%m-%d %H:%M:%S", tz = timezone_local),
    crds_location = clean_crds_location(location),
    analyzer      = toupper(trimws(as.character(analyzer))),
    CRDS_NH3      = suppressWarnings(as.numeric(NH3)),
    CRDS_CO2      = suppressWarnings(as.numeric(CO2))
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
          file.path(out_dir, "tables", "node_hourly_matched_december.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 7. First 24 hours of December calibration baseline
# -----------------------------------------------------------------------------
december_baseline <- node_hourly_matched |>
  filter(datetime_hour >= december_start,
         datetime_hour < december_calibration_end)

fit_calibration_model <- function(df, sensor_episode_col, raw_col, ref_col, gas_label) {
  split_df <- split(df, df[[sensor_episode_col]])

  bind_rows(lapply(split_df, function(d) {
    ok <- !is.na(d[[raw_col]]) & !is.na(d[[ref_col]])
    if (sum(ok) < 10) {
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
          file.path(out_dir, "tables", "nh3_calibration_models_first24h_december.csv"),
          row.names = FALSE)

write.csv(co2_calibration_models,
          file.path(out_dir, "tables", "co2_calibration_models_first24h_december.csv"),
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
    OTICE_NH3_raw, OTICE_NH3_cal, CRDS_NH3,
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
    OTICE_CO2_raw, OTICE_CO2_cal, CRDS_CO2,
    CO2_predicted_from_december24h, CO2_bias_after_december24h_cal
  )

write.csv(nh3_predictions,
          file.path(out_dir, "tables", "nh3_predictions_after_first24h_december_calibration.csv"),
          row.names = FALSE)

write.csv(co2_predictions,
          file.path(out_dir, "tables", "co2_predictions_after_first24h_december_calibration.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 10. December-only node-level visual comparison
# -----------------------------------------------------------------------------
december_nh3_comparison <- nh3_predictions |>
  filter(period_id == "C8",
         datetime_hour >= december_start,
         datetime_hour < december_end) |>
  select(period_id, period_label, datetime_hour, node, nh3_sensor_id,
         OTICE_NH3_raw, OTICE_NH3_cal, CRDS_NH3, NH3_predicted_from_december24h)

december_co2_comparison <- co2_predictions |>
  filter(period_id == "C8",
         datetime_hour >= december_start,
         datetime_hour < december_end) |>
  select(period_id, period_label, datetime_hour, node, co2_sensor_id,
         OTICE_CO2_raw, OTICE_CO2_cal, CRDS_CO2, CO2_predicted_from_december24h)

write.csv(december_nh3_comparison,
          file.path(out_dir, "tables", "december_nh3_node_vs_crds.csv"),
          row.names = FALSE)

write.csv(december_co2_comparison,
          file.path(out_dir, "tables", "december_co2_node_vs_crds.csv"),
          row.names = FALSE)

if (nrow(december_nh3_comparison) > 0) {
  save_plot(
    make_timeseries_plot(
      december_nh3_comparison,
      gas_label = "NH3",
      raw_col = "OTICE_NH3_raw",
      ref_col = "CRDS_NH3",
      pred_col = "NH3_predicted_from_december24h",
      y_label = "NH3 (ppm)",
      baseline_start = december_start,
      baseline_end = december_calibration_end
    ),
    file.path(out_dir, "plots", "December_NH3_all_nodes_vs_CRDS.png"),
    width = 16,
    height = 12
  )

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
        baseline_end = december_calibration_end
      ),
      file.path(out_dir, "plots", paste0("December_NH3_node_", node_id, ".png")),
      width = 14,
      height = 8
    )
  }
}

if (nrow(december_co2_comparison) > 0) {
  save_plot(
    make_timeseries_plot(
      december_co2_comparison,
      gas_label = "CO2",
      raw_col = "OTICE_CO2_raw",
      ref_col = "CRDS_CO2",
      pred_col = "CO2_predicted_from_december24h",
      y_label = "CO2 (ppm)",
      baseline_start = december_start,
      baseline_end = december_calibration_end
    ),
    file.path(out_dir, "plots", "December_CO2_all_nodes_vs_CRDS.png"),
    width = 16,
    height = 12
  )

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
        baseline_end = december_calibration_end
      ),
      file.path(out_dir, "plots", paste0("December_CO2_node_", node_id, ".png")),
      width = 14,
      height = 8
    )
  }
}

# -----------------------------------------------------------------------------
# 11. Console summary
# -----------------------------------------------------------------------------
cat("\n============================================================\n")
cat("OTICE December node-versus-CRDS comparison\n")
cat("============================================================\n")
cat("Focus campaign: 2025-12-09 to 2025-12-31 (C8)\n")
cat("Calibration baseline: first 24 hours of 2025-12-09\n")
cat("Matched hourly rows in December:", nrow(node_hourly_matched), "\n")
cat("December NH3 rows:", nrow(december_nh3_comparison), "\n")
cat("December CO2 rows:", nrow(december_co2_comparison), "\n")
cat("\nFiles written to:", out_dir, "\n")
cat("  tables/campaign_reference_map.csv\n")
cat("  tables/node_hourly_matched_december.csv\n")
cat("  tables/nh3_calibration_models_first24h_december.csv\n")
cat("  tables/co2_calibration_models_first24h_december.csv\n")
cat("  tables/nh3_predictions_after_first24h_december_calibration.csv\n")
cat("  tables/co2_predictions_after_first24h_december_calibration.csv\n")
cat("  tables/december_nh3_node_vs_crds.csv\n")
cat("  tables/december_co2_node_vs_crds.csv\n")
cat("  plots/December_NH3_all_nodes_vs_CRDS.png\n")
cat("  plots/December_CO2_all_nodes_vs_CRDS.png\n")
cat("  plots/*.png\n")
