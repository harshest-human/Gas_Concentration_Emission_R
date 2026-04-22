# =============================================================================
# OTICE versus CRDS concentration
# Fresh hourly matching workflow for aging analysis
#
# Main idea:
# - import OTICE and CRDS directly from processed_data
# - standardize timestamps to hourly bins
# - map each OTICE node to the best CRDS reference within each plot period
# - calibrate with the first 12 matched hourly values only
# - apply the fixed fit forward to the rest of the data
# - summarise to daily values only for tables and plots
# =============================================================================

library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(scales)

# -----------------------------------------------------------------------------
# 1. Paths and settings
# -----------------------------------------------------------------------------
timezone_local <- "Europe/Berlin"
calibration_hours <- 48

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args_all[grepl(file_arg, args_all)])
base_dir <- if (length(script_path) > 0) dirname(normalizePath(script_path[1])) else getwd()

otice_dir <- file.path(base_dir, "processed_data", "OTICE_processed")
crds_dir <- file.path(base_dir, "processed_data", "CRDS8_processed")
out_dir <- file.path(base_dir, "output", "OTICE_versus_CRDS_concentration")
tables_dir <- file.path(out_dir, "tables")
plots_daily_dir <- file.path(out_dir, "plots", "daily")

if (dir.exists(out_dir)) {
  unlink(out_dir, recursive = TRUE, force = TRUE)
}

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_daily_dir, recursive = TRUE, showWarnings = FALSE)

period_lookup <- tibble(
  period_id = c("SepOct", "November", "December"),
  period_label = c("September-October 2025", "November 2025", "December 2025"),
  period_start = as.POSIXct(
    c("2025-09-29 00:00:00", "2025-11-01 00:00:00", "2025-12-01 00:00:00"),
    tz = timezone_local
  ),
  period_end = as.POSIXct(
    c("2025-11-01 00:00:00", "2025-12-01 00:00:00", "2026-01-01 00:00:00"),
    tz = timezone_local
  )
)

# -----------------------------------------------------------------------------
# 2. Helper functions
# -----------------------------------------------------------------------------
safe_mean <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

safe_rmse <- function(x, y) {
  ok <- complete.cases(x, y)
  if (!any(ok)) return(NA_real_)
  sqrt(mean((x[ok] - y[ok])^2))
}

safe_mae <- function(x, y) {
  ok <- complete.cases(x, y)
  if (!any(ok)) return(NA_real_)
  mean(abs(x[ok] - y[ok]))
}

safe_cor <- function(x, y, min_points = 6) {
  ok <- complete.cases(x, y)
  if (sum(ok) < min_points) return(NA_real_)
  if (sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  cor(x[ok], y[ok])
}

safe_zrmse <- function(x, y, min_points = 6) {
  ok <- complete.cases(x, y)
  if (sum(ok) < min_points) return(NA_real_)
  if (sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  x_scaled <- as.numeric(scale(x[ok]))
  y_scaled <- as.numeric(scale(y[ok]))
  sqrt(mean((x_scaled - y_scaled)^2))
}

clean_otice_node <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("^O", "", x)
  x <- sub("^0+", "", x)
  x
}

clean_crds_location <- function(x) {
  x <- trimws(as.character(x))
  numeric_mask <- grepl("^[0-9]+$", x)
  x[numeric_mask] <- as.character(as.integer(x[numeric_mask]))
  x
}

assign_period_id <- function(datetime_value) {
  case_when(
    datetime_value >= as.POSIXct("2025-09-29 00:00:00", tz = timezone_local) &
      datetime_value < as.POSIXct("2025-11-01 00:00:00", tz = timezone_local) ~ "SepOct",
    datetime_value >= as.POSIXct("2025-11-01 00:00:00", tz = timezone_local) &
      datetime_value < as.POSIXct("2025-12-01 00:00:00", tz = timezone_local) ~ "November",
    datetime_value >= as.POSIXct("2025-12-01 00:00:00", tz = timezone_local) &
      datetime_value < as.POSIXct("2026-01-01 00:00:00", tz = timezone_local) ~ "December",
    TRUE ~ NA_character_
  )
}

period_prefix <- function(period_id_value) {
  case_when(
    period_id_value == "SepOct" ~ "sep_oct",
    period_id_value == "November" ~ "nov",
    period_id_value == "December" ~ "dec",
    TRUE ~ period_id_value
  )
}

node_file_label <- function(node_value) {
  paste0("O", suppressWarnings(as.integer(as.character(node_value))))
}

fmt_num <- function(x, digits = 2) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}

gas_style <- function(gas_label) {
  if (gas_label == "CO2") {
    list(
      y_label = "CO2 (ppm)",
      raw_color = "#7DB7E8",
      fit_color = "#145A96",
      ref_color = "#3F3F46",
      accuracy = 1
    )
  } else {
    list(
      y_label = "NH3 (ppm)",
      raw_color = "#8BD38B",
      fit_color = "#1F8A4C",
      ref_color = "#3F3F46",
      accuracy = 0.01
    )
  }
}

build_mapping_score <- function(joined_hourly) {
  overlap_hours <- nrow(joined_hourly)
  cor_co2 <- safe_cor(joined_hourly$OTICE_CO2_raw, joined_hourly$CRDS_CO2)
  cor_nh3 <- safe_cor(joined_hourly$OTICE_NH3_raw, joined_hourly$CRDS_NH3)
  zrmse_co2 <- safe_zrmse(joined_hourly$OTICE_CO2_raw, joined_hourly$CRDS_CO2)
  zrmse_nh3 <- safe_zrmse(joined_hourly$OTICE_NH3_raw, joined_hourly$CRDS_NH3)

  mean_cor <- if (all(is.na(c(cor_co2, cor_nh3)))) NA_real_ else mean(c(cor_co2, cor_nh3), na.rm = TRUE)
  mean_zrmse <- if (all(is.na(c(zrmse_co2, zrmse_nh3)))) NA_real_ else mean(c(zrmse_co2, zrmse_nh3), na.rm = TRUE)
  gas_count <- sum(!is.na(c(cor_co2, cor_nh3)))

  score <- if (is.na(mean_cor)) {
    -Inf
  } else {
    mean_cor - 0.20 * ifelse(is.na(mean_zrmse), 2, mean_zrmse) + 0.005 * overlap_hours + 0.05 * gas_count
  }

  tibble(
    overlap_hours = overlap_hours,
    gas_count = gas_count,
    cor_co2 = cor_co2,
    cor_nh3 = cor_nh3,
    zrmse_co2 = zrmse_co2,
    zrmse_nh3 = zrmse_nh3,
    mean_cor = mean_cor,
    mean_zrmse = mean_zrmse,
    mapping_score = score
  )
}

fit_initial_window_model <- function(df, raw_col, ref_col, calibration_hours_value = 12) {
  usable <- df |>
    filter(complete.cases(.data[[raw_col]], .data[[ref_col]])) |>
    arrange(datetime_hour)

  calibration_df <- usable |>
    slice_head(n = calibration_hours_value)

  predicted_all <- rep(NA_real_, nrow(df))
  intercept <- NA_real_
  slope <- NA_real_
  r_squared <- NA_real_
  model_type <- "none"

  if (nrow(calibration_df) >= 4 &&
      sd(calibration_df[[raw_col]]) > 0 &&
      sd(calibration_df[[ref_col]]) > 0) {
    fit <- lm(calibration_df[[ref_col]] ~ calibration_df[[raw_col]])
    intercept <- unname(coef(fit)[1])
    slope <- unname(coef(fit)[2])
    r_squared <- summary(fit)$r.squared
    model_type <- "linear"
    predicted_all <- intercept + slope * df[[raw_col]]
  } else if (nrow(calibration_df) >= 2) {
    intercept <- mean(calibration_df[[ref_col]] - calibration_df[[raw_col]])
    slope <- 1
    model_type <- "offset_only"
    predicted_all <- df[[raw_col]] + intercept
  } else if (nrow(calibration_df) == 1) {
    intercept <- calibration_df[[ref_col]][1] - calibration_df[[raw_col]][1]
    slope <- 1
    model_type <- "single_point_offset"
    predicted_all <- df[[raw_col]] + intercept
  }

  calibration_start <- if (nrow(calibration_df) > 0) min(calibration_df$datetime_hour) else as.POSIXct(NA, tz = timezone_local)
  calibration_end <- if (nrow(calibration_df) > 0) max(calibration_df$datetime_hour) + hours(1) else as.POSIXct(NA, tz = timezone_local)

  raw_rmse <- safe_rmse(df[[raw_col]], df[[ref_col]])
  fitted_rmse <- safe_rmse(predicted_all, df[[ref_col]])
  raw_mae <- safe_mae(df[[raw_col]], df[[ref_col]])
  fitted_mae <- safe_mae(predicted_all, df[[ref_col]])

  list(
    predicted = predicted_all,
    calibration_rows = calibration_df,
    stats = tibble(
      model_type = model_type,
      n_calibration_hours = nrow(calibration_df),
      calibration_start = calibration_start,
      calibration_end = calibration_end,
      intercept = intercept,
      slope = slope,
      r_squared = r_squared,
      raw_rmse = raw_rmse,
      fitted_rmse = fitted_rmse,
      raw_mae = raw_mae,
      fitted_mae = fitted_mae,
      rmse_improvement = raw_rmse - fitted_rmse,
      mae_improvement = raw_mae - fitted_mae
    )
  )
}

build_node_stats_text <- function(model_row) {
  if (nrow(model_row) == 0) return("No fit statistics available")

  eq_text <- if (is.na(model_row$slope[1]) || is.na(model_row$intercept[1])) {
    "fit: not available"
  } else if (model_row$model_type[1] %in% c("offset_only", "single_point_offset")) {
    paste0("fit: y = x + ", fmt_num(model_row$intercept[1], 3))
  } else {
    paste0("fit: y = ", fmt_num(model_row$slope[1], 3), "x + ", fmt_num(model_row$intercept[1], 3))
  }

  paste(
    paste0("reference: ", model_row$analyzer[1], " location ", model_row$crds_location[1]),
    paste0("model: ", model_row$model_type[1]),
    paste0("cal hrs = ", model_row$n_calibration_hours[1]),
    paste0("R2 = ", fmt_num(model_row$r_squared[1], 3)),
    paste0("RMSE raw ", fmt_num(model_row$raw_rmse[1], 3), " -> ", fmt_num(model_row$fitted_rmse[1], 3)),
    eq_text,
    sep = " | "
  )
}

build_average_stats_text <- function(stats_row) {
  paste(
    paste0("days = ", stats_row$n_days[1]),
    paste0("nodes = ", stats_row$n_nodes[1]),
    paste0("RMSE raw ", fmt_num(stats_row$raw_rmse[1], 3), " -> ", fmt_num(stats_row$fitted_rmse[1], 3)),
    paste0("MAE raw ", fmt_num(stats_row$raw_mae[1], 3), " -> ", fmt_num(stats_row$fitted_mae[1], 3)),
    sep = " | "
  )
}

make_node_daily_plot <- function(df, model_row, period_label_value, node_value, gas_label) {
  style <- gas_style(gas_label)

  plot_df <- if (gas_label == "CO2") {
    df |>
      select(date_day, CRDS = CRDS_CO2, OTICE_raw = OTICE_CO2_raw, OTICE_fitted = OTICE_CO2_fitted) |>
      pivot_longer(cols = c(CRDS, OTICE_raw, OTICE_fitted), names_to = "series", values_to = "value")
  } else {
    df |>
      select(date_day, CRDS = CRDS_NH3, OTICE_raw = OTICE_NH3_raw, OTICE_fitted = OTICE_NH3_fitted) |>
      pivot_longer(cols = c(CRDS, OTICE_raw, OTICE_fitted), names_to = "series", values_to = "value")
  }

  color_values <- c(CRDS = style$ref_color, OTICE_raw = style$raw_color, OTICE_fitted = style$fit_color)
  legend_labels <- c(CRDS = "CRDS reference", OTICE_raw = "OTICE raw", OTICE_fitted = "OTICE fitted")
  linetype_values <- c(CRDS = "solid", OTICE_raw = "dotted", OTICE_fitted = "solid")

  ggplot(plot_df, aes(x = date_day, y = value, color = series, linetype = series)) +
    annotate(
      "rect",
      xmin = as.Date(model_row$calibration_start[1], tz = timezone_local),
      xmax = as.Date(model_row$calibration_end[1], tz = timezone_local),
      ymin = -Inf,
      ymax = Inf,
      fill = "goldenrod",
      alpha = 0.10
    ) +
    geom_line(linewidth = 0.9, na.rm = TRUE) +
    geom_point(size = 2, alpha = 0.75, na.rm = TRUE) +
    scale_color_manual(values = color_values, labels = legend_labels) +
    scale_linetype_manual(values = linetype_values, labels = legend_labels) +
    scale_x_date(date_breaks = "1 day", date_labels = "%d %b", expand = expansion(mult = c(0.01, 0.02))) +
    scale_y_continuous(labels = label_number(accuracy = style$accuracy), n.breaks = 12) +
    labs(
      title = paste0(period_label_value, " ", gas_label, ": OTICE node ", node_value),
      subtitle = build_node_stats_text(model_row),
      x = NULL,
      y = style$y_label,
      color = NULL,
      linetype = NULL,
      caption = "Yellow band marks the actual time span of the first matched 48 hourly values used for calibration."
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

make_average_daily_plot <- function(df, avg_stats_row, period_label_value, gas_label, reference_mix_text) {
  style <- gas_style(gas_label)

  plot_df <- if (gas_label == "CO2") {
    df |>
      select(date_day, CRDS = CRDS_CO2_mean, OTICE_raw = OTICE_CO2_raw_mean, OTICE_fitted = OTICE_CO2_fitted_mean) |>
      pivot_longer(cols = c(CRDS, OTICE_raw, OTICE_fitted), names_to = "series", values_to = "value")
  } else {
    df |>
      select(date_day, CRDS = CRDS_NH3_mean, OTICE_raw = OTICE_NH3_raw_mean, OTICE_fitted = OTICE_NH3_fitted_mean) |>
      pivot_longer(cols = c(CRDS, OTICE_raw, OTICE_fitted), names_to = "series", values_to = "value")
  }

  color_values <- c(CRDS = style$ref_color, OTICE_raw = style$raw_color, OTICE_fitted = style$fit_color)
  legend_labels <- c(CRDS = "CRDS mean", OTICE_raw = "OTICE raw mean", OTICE_fitted = "OTICE fitted mean")
  linetype_values <- c(CRDS = "solid", OTICE_raw = "dotted", OTICE_fitted = "solid")

  ggplot(plot_df, aes(x = date_day, y = value, color = series, linetype = series)) +
    geom_line(linewidth = 0.9, na.rm = TRUE) +
    geom_point(size = 2, alpha = 0.75, na.rm = TRUE) +
    scale_color_manual(values = color_values, labels = legend_labels) +
    scale_linetype_manual(values = linetype_values, labels = legend_labels) +
    scale_x_date(date_breaks = "1 day", date_labels = "%d %b", expand = expansion(mult = c(0.01, 0.02))) +
    scale_y_continuous(labels = label_number(accuracy = style$accuracy), n.breaks = 12) +
    labs(
      title = paste0(period_label_value, " ", gas_label, ": average of OTICE nodes versus average CRDS"),
      subtitle = build_average_stats_text(avg_stats_row),
      x = NULL,
      y = style$y_label,
      color = NULL,
      linetype = NULL,
      caption = paste0("Reference mix: ", reference_mix_text)
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

# -----------------------------------------------------------------------------
# 3. Import OTICE data and build hourly data
# -----------------------------------------------------------------------------
otice_files <- list.files(otice_dir, pattern = "^min_calibrated.*\\.csv$", full.names = TRUE)
otice_files <- otice_files[file.info(otice_files)$size > 0]

otice_raw <- bind_rows(lapply(otice_files, function(file_path) {
  read.csv(file_path, stringsAsFactors = FALSE)
})) |>
  transmute(
    datetime = as.POSIXct(.data[["Datetime_Berlin"]], format = "%Y-%m-%d %H:%M:%S", tz = timezone_local),
    datetime_hour = floor_date(datetime, "hour"),
    period_id = assign_period_id(datetime_hour),
    node = clean_otice_node(.data[["Node"]]),
    OTICE_NH3_raw = suppressWarnings(as.numeric(.data[["NH3_ppm"]])),
    OTICE_CO2_raw = suppressWarnings(as.numeric(.data[["CO2.RAW"]])),
    OTICE_NH3_internal = suppressWarnings(as.numeric(.data[["NH3_ppm_barn"]])),
    OTICE_CO2_internal = suppressWarnings(as.numeric(.data[["CO2.AVG_barn"]]))
  ) |>
  filter(!is.na(datetime_hour), !is.na(period_id), !is.na(node))

otice_hourly <- otice_raw |>
  group_by(period_id, datetime_hour, node) |>
  summarise(
    OTICE_NH3_raw = safe_mean(OTICE_NH3_raw),
    OTICE_CO2_raw = safe_mean(OTICE_CO2_raw),
    OTICE_NH3_internal = safe_mean(OTICE_NH3_internal),
    OTICE_CO2_internal = safe_mean(OTICE_CO2_internal),
    otice_points = n(),
    .groups = "drop"
  ) |>
  left_join(period_lookup |> select(period_id, period_label), by = "period_id")

# -----------------------------------------------------------------------------
# 4. Import CRDS data and build hourly data
# -----------------------------------------------------------------------------
crds_files <- list.files(crds_dir, pattern = "\\.csv$", full.names = TRUE)
crds_files <- crds_files[file.info(crds_files)$size > 0]

crds_raw <- bind_rows(lapply(crds_files, function(file_path) {
  df <- read.csv(file_path, stringsAsFactors = FALSE)
  df$location <- as.character(df$location)
  df$analyzer <- as.character(df$analyzer)
  df
})) |>
  transmute(
    datetime = as.POSIXct(.data[["DATE.TIME"]], format = "%Y-%m-%d %H:%M:%S", tz = timezone_local),
    datetime_hour = floor_date(datetime, "hour"),
    period_id = assign_period_id(datetime_hour),
    analyzer = toupper(trimws(as.character(.data[["analyzer"]]))),
    crds_location = clean_crds_location(.data[["location"]]),
    CRDS_NH3 = suppressWarnings(as.numeric(.data[["NH3"]])),
    CRDS_CO2 = suppressWarnings(as.numeric(.data[["CO2"]]))
  ) |>
  filter(!is.na(datetime_hour), !is.na(period_id), !is.na(analyzer), !is.na(crds_location))

crds_hourly <- crds_raw |>
  group_by(period_id, datetime_hour, analyzer, crds_location) |>
  summarise(
    CRDS_NH3 = safe_mean(CRDS_NH3),
    CRDS_CO2 = safe_mean(CRDS_CO2),
    crds_points = n(),
    .groups = "drop"
  ) |>
  left_join(period_lookup |> select(period_id, period_label), by = "period_id")

# -----------------------------------------------------------------------------
# 5. Fresh hourly mapping for each node within each plot period
# -----------------------------------------------------------------------------
mapping_best_rows <- list()
mapping_score_rows <- list()
map_index <- 1

for (period_id_value in period_lookup$period_id) {
  period_otice <- otice_hourly |> filter(period_id == period_id_value)
  period_crds <- crds_hourly |> filter(period_id == period_id_value)

  if (nrow(period_otice) == 0 || nrow(period_crds) == 0) {
    next
  }

  candidate_refs <- period_crds |>
    distinct(analyzer, crds_location)

  for (node_value in sort(unique(period_otice$node))) {
    node_hourly <- period_otice |> filter(node == node_value)

    score_table <- bind_rows(lapply(seq_len(nrow(candidate_refs)), function(i) {
      candidate <- candidate_refs[i, ]

      joined_hourly <- node_hourly |>
        select(period_id, period_label, datetime_hour, node, OTICE_NH3_raw, OTICE_CO2_raw, OTICE_NH3_internal, OTICE_CO2_internal) |>
        inner_join(
          period_crds |>
            filter(analyzer == candidate$analyzer, crds_location == candidate$crds_location) |>
            select(period_id, datetime_hour, analyzer, crds_location, CRDS_NH3, CRDS_CO2),
          by = c("period_id", "datetime_hour")
        )

      if (nrow(joined_hourly) < 2) {
        return(NULL)
      }

      bind_cols(
        tibble(
          period_id = period_id_value,
          node = node_value,
          analyzer = candidate$analyzer,
          crds_location = candidate$crds_location
        ),
        build_mapping_score(joined_hourly)
      )
    }))

    if (nrow(score_table) == 0) {
      next
    }

    best_row <- score_table |>
      arrange(desc(mapping_score), desc(overlap_hours), desc(gas_count), desc(mean_cor), mean_zrmse) |>
      slice(1) |>
      left_join(period_lookup |> select(period_id, period_label), by = "period_id") |>
      relocate(period_id, period_label, node, analyzer, crds_location)

    mapping_best_rows[[map_index]] <- best_row
    mapping_score_rows[[map_index]] <- score_table |>
      left_join(period_lookup |> select(period_id, period_label), by = "period_id") |>
      relocate(period_id, period_label, node, analyzer, crds_location)
    map_index <- map_index + 1
  }
}

reference_mapping <- bind_rows(mapping_best_rows) |>
  arrange(factor(period_id, levels = period_lookup$period_id), as.numeric(node), analyzer, crds_location)

reference_mapping_scores <- bind_rows(mapping_score_rows) |>
  arrange(factor(period_id, levels = period_lookup$period_id), as.numeric(node), desc(mapping_score))

# -----------------------------------------------------------------------------
# 6. Build matched hourly node-versus-reference data
# -----------------------------------------------------------------------------
matched_hourly <- otice_hourly |>
  inner_join(
    reference_mapping |>
      select(period_id, node, analyzer, crds_location, overlap_hours, gas_count, cor_co2, cor_nh3, mean_cor, mean_zrmse, mapping_score),
    by = c("period_id", "node")
  ) |>
  inner_join(
    crds_hourly |>
      select(period_id, datetime_hour, analyzer, crds_location, CRDS_NH3, CRDS_CO2, crds_points),
    by = c("period_id", "datetime_hour", "analyzer", "crds_location")
  ) |>
  arrange(period_id, node, datetime_hour)

# -----------------------------------------------------------------------------
# 7. Fit fixed first-12-hour calibrations per node and period
# -----------------------------------------------------------------------------
group_keys <- matched_hourly |>
  distinct(period_id, period_label, node, analyzer, crds_location) |>
  arrange(factor(period_id, levels = period_lookup$period_id), as.numeric(node))

hourly_prediction_groups <- list()
node_model_stat_rows <- list()

for (i in seq_len(nrow(group_keys))) {
  group_info <- group_keys[i, ]

  group_df <- matched_hourly |>
    filter(
      period_id == group_info$period_id,
      node == group_info$node,
      analyzer == group_info$analyzer,
      crds_location == group_info$crds_location
    ) |>
    arrange(datetime_hour)

  co2_fit <- fit_initial_window_model(group_df, "OTICE_CO2_raw", "CRDS_CO2", calibration_hours)
  nh3_fit <- fit_initial_window_model(group_df, "OTICE_NH3_raw", "CRDS_NH3", calibration_hours)

  group_df$OTICE_CO2_fitted <- co2_fit$predicted
  group_df$OTICE_NH3_fitted <- nh3_fit$predicted

  hourly_prediction_groups[[i]] <- group_df

  node_model_stat_rows[[length(node_model_stat_rows) + 1]] <- bind_cols(
    group_info,
    tibble(gas = "CO2"),
    co2_fit$stats
  )

  node_model_stat_rows[[length(node_model_stat_rows) + 1]] <- bind_cols(
    group_info,
    tibble(gas = "NH3"),
    nh3_fit$stats
  )
}

node_hourly_comparison <- bind_rows(hourly_prediction_groups)
node_model_stats <- bind_rows(node_model_stat_rows) |>
  arrange(factor(period_id, levels = period_lookup$period_id), gas, as.numeric(node))

# -----------------------------------------------------------------------------
# 8. Summarise hourly outputs to daily node-level outputs
# -----------------------------------------------------------------------------
node_daily_comparison <- node_hourly_comparison |>
  mutate(date_day = as.Date(datetime_hour, tz = timezone_local)) |>
  group_by(period_id, period_label, date_day, node, analyzer, crds_location) |>
  summarise(
    OTICE_NH3_raw = safe_mean(OTICE_NH3_raw),
    OTICE_CO2_raw = safe_mean(OTICE_CO2_raw),
    OTICE_NH3_fitted = safe_mean(OTICE_NH3_fitted),
    OTICE_CO2_fitted = safe_mean(OTICE_CO2_fitted),
    CRDS_NH3 = safe_mean(CRDS_NH3),
    CRDS_CO2 = safe_mean(CRDS_CO2),
    matched_hours = n(),
    .groups = "drop"
  ) |>
  arrange(factor(period_id, levels = period_lookup$period_id), as.numeric(node), date_day)

# -----------------------------------------------------------------------------
# 9. Build daily average tables from node-level fitted values
# -----------------------------------------------------------------------------
average_daily_comparison <- node_daily_comparison |>
  group_by(period_id, period_label, date_day) |>
  summarise(
    OTICE_NH3_raw_mean = safe_mean(OTICE_NH3_raw),
    OTICE_CO2_raw_mean = safe_mean(OTICE_CO2_raw),
    OTICE_NH3_fitted_mean = safe_mean(OTICE_NH3_fitted),
    OTICE_CO2_fitted_mean = safe_mean(OTICE_CO2_fitted),
    CRDS_NH3_mean = safe_mean(CRDS_NH3),
    CRDS_CO2_mean = safe_mean(CRDS_CO2),
    n_nodes = n_distinct(node),
    .groups = "drop"
  ) |>
  arrange(factor(period_id, levels = period_lookup$period_id), date_day)

average_stat_rows <- list()

for (period_id_value in period_lookup$period_id) {
  avg_df <- average_daily_comparison |> filter(period_id == period_id_value)
  if (nrow(avg_df) == 0) next

  average_stat_rows[[length(average_stat_rows) + 1]] <- tibble(
    period_id = period_id_value,
    period_label = unique(avg_df$period_label)[1],
    gas = "CO2",
    n_days = nrow(avg_df),
    n_nodes = max(avg_df$n_nodes, na.rm = TRUE),
    raw_rmse = safe_rmse(avg_df$OTICE_CO2_raw_mean, avg_df$CRDS_CO2_mean),
    fitted_rmse = safe_rmse(avg_df$OTICE_CO2_fitted_mean, avg_df$CRDS_CO2_mean),
    raw_mae = safe_mae(avg_df$OTICE_CO2_raw_mean, avg_df$CRDS_CO2_mean),
    fitted_mae = safe_mae(avg_df$OTICE_CO2_fitted_mean, avg_df$CRDS_CO2_mean)
  )

  average_stat_rows[[length(average_stat_rows) + 1]] <- tibble(
    period_id = period_id_value,
    period_label = unique(avg_df$period_label)[1],
    gas = "NH3",
    n_days = nrow(avg_df),
    n_nodes = max(avg_df$n_nodes, na.rm = TRUE),
    raw_rmse = safe_rmse(avg_df$OTICE_NH3_raw_mean, avg_df$CRDS_NH3_mean),
    fitted_rmse = safe_rmse(avg_df$OTICE_NH3_fitted_mean, avg_df$CRDS_NH3_mean),
    raw_mae = safe_mae(avg_df$OTICE_NH3_raw_mean, avg_df$CRDS_NH3_mean),
    fitted_mae = safe_mae(avg_df$OTICE_NH3_fitted_mean, avg_df$CRDS_NH3_mean)
  )
}

average_model_stats <- bind_rows(average_stat_rows)

# -----------------------------------------------------------------------------
# 10. Save tables
# -----------------------------------------------------------------------------
write.csv(otice_hourly, file.path(tables_dir, "otice_hourly.csv"), row.names = FALSE)
write.csv(crds_hourly, file.path(tables_dir, "crds_hourly.csv"), row.names = FALSE)
write.csv(reference_mapping_scores, file.path(tables_dir, "reference_mapping_scores.csv"), row.names = FALSE)
write.csv(reference_mapping, file.path(tables_dir, "reference_mapping_by_period.csv"), row.names = FALSE)
write.csv(node_model_stats, file.path(tables_dir, "node_model_stats.csv"), row.names = FALSE)
write.csv(node_hourly_comparison, file.path(tables_dir, "node_hourly_comparison.csv"), row.names = FALSE)
write.csv(node_daily_comparison, file.path(tables_dir, "node_daily_comparison.csv"), row.names = FALSE)
write.csv(average_model_stats, file.path(tables_dir, "average_model_stats.csv"), row.names = FALSE)
write.csv(average_daily_comparison, file.path(tables_dir, "average_daily_comparison.csv"), row.names = FALSE)

for (period_id_value in period_lookup$period_id) {
  prefix_value <- period_prefix(period_id_value)

  write.csv(
    reference_mapping |> filter(period_id == period_id_value),
    file.path(tables_dir, paste0(prefix_value, "_reference_mapping.csv")),
    row.names = FALSE
  )

  write.csv(
    node_model_stats |> filter(period_id == period_id_value),
    file.path(tables_dir, paste0(prefix_value, "_node_model_stats.csv")),
    row.names = FALSE
  )

  write.csv(
    node_hourly_comparison |> filter(period_id == period_id_value),
    file.path(tables_dir, paste0(prefix_value, "_node_hourly_comparison.csv")),
    row.names = FALSE
  )

  write.csv(
    node_daily_comparison |> filter(period_id == period_id_value),
    file.path(tables_dir, paste0(prefix_value, "_node_daily_comparison.csv")),
    row.names = FALSE
  )

  write.csv(
    average_daily_comparison |> filter(period_id == period_id_value),
    file.path(tables_dir, paste0(prefix_value, "_average_daily_comparison.csv")),
    row.names = FALSE
  )
}

# -----------------------------------------------------------------------------
# 11. Make daily plots with yellow calibration band
# -----------------------------------------------------------------------------
for (period_id_value in period_lookup$period_id) {
  period_row <- period_lookup |> filter(period_id == period_id_value)
  prefix_value <- period_prefix(period_id_value)

  period_node_daily <- node_daily_comparison |> filter(period_id == period_id_value)
  period_average_daily <- average_daily_comparison |> filter(period_id == period_id_value)
  period_mapping <- reference_mapping |> filter(period_id == period_id_value)

  if (nrow(period_node_daily) == 0) next

  reference_mix_text <- period_mapping |>
    count(analyzer, crds_location, name = "nodes_here") |>
    arrange(desc(nodes_here), analyzer, crds_location) |>
    transmute(label = paste0(nodes_here, "x ", analyzer, "-", crds_location)) |>
    pull(label) |>
    paste(collapse = "; ")

  for (gas_label in c("CO2", "NH3")) {
    avg_stats_row <- average_model_stats |>
      filter(period_id == period_id_value, gas == gas_label)

    if (nrow(period_average_daily) > 0 && nrow(avg_stats_row) > 0) {
      avg_plot <- make_average_daily_plot(
        df = period_average_daily,
        avg_stats_row = avg_stats_row,
        period_label_value = period_row$period_label[1],
        gas_label = gas_label,
        reference_mix_text = reference_mix_text
      )

      ggsave(
        filename = file.path(plots_daily_dir, paste0(prefix_value, "_", gas_label, "_all_nodes.png")),
        plot = avg_plot,
        width = 15,
        height = 8,
        dpi = 150
      )
    }

    for (node_value in sort(unique(period_node_daily$node))) {
      node_df <- period_node_daily |>
        filter(node == node_value)

      model_row <- node_model_stats |>
        filter(period_id == period_id_value, gas == gas_label, node == node_value)

      if (nrow(node_df) == 0 || nrow(model_row) == 0) next

      node_plot <- make_node_daily_plot(
        df = node_df,
        model_row = model_row,
        period_label_value = period_row$period_label[1],
        node_value = node_value,
        gas_label = gas_label
      )

      ggsave(
        filename = file.path(plots_daily_dir, paste0(prefix_value, "_", gas_label, "_", node_file_label(node_value), ".png")),
        plot = node_plot,
        width = 15,
        height = 8,
        dpi = 150
      )
    }
  }
}

# -----------------------------------------------------------------------------
# 12. Console summary
# -----------------------------------------------------------------------------
cat("\n============================================================\n")
cat("OTICE versus CRDS aging workflow\n")
cat("============================================================\n")
cat("Base directory:", base_dir, "\n")
cat("Calibration window:", calibration_hours, "matched hourly values\n")
cat("Output directory:", out_dir, "\n")
cat("OTICE hourly rows:", nrow(otice_hourly), "\n")
cat("CRDS hourly rows:", nrow(crds_hourly), "\n")
cat("Mapped node-period references:", nrow(reference_mapping), "\n")
cat("Matched node hourly rows:", nrow(node_hourly_comparison), "\n")
cat("Matched node daily rows:", nrow(node_daily_comparison), "\n")
cat("Average daily rows:", nrow(average_daily_comparison), "\n")
cat("\nFiles written to:\n")
cat("  ", tables_dir, "\n")
cat("  ", plots_daily_dir, "\n")
