# =============================================================================
# Clean otice data
#
# Outputs:
# - node-level hourly raw and 48-hour fitted values
# - hourly averaged raw and fitted inside concentrations
# - hourly fitted delta values using CRDS S as outside background
# - reference mapping and model statistics
# - monthly node plots
# =============================================================================

library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
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
calibration_hours <- 48

clean_dir <- file.path(project_dir, "workflows", "device_comparison", "clean_data", "otice_clean")
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
    orders = c(
      "ymd HMS", "ymd HM",
      "Ymd HMS", "Ymd HM",
      "Y/m/d HMS", "Y/m/d HM"
    ),
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

safe_cor <- function(x, y, min_points = 6) {
  ok <- complete.cases(x, y)
  if (sum(ok) < min_points) {
    return(NA_real_)
  }
  if (sd(x[ok]) == 0 || sd(y[ok]) == 0) {
    return(NA_real_)
  }
  cor(x[ok], y[ok])
}

safe_zrmse <- function(x, y, min_points = 6) {
  ok <- complete.cases(x, y)
  if (sum(ok) < min_points) {
    return(NA_real_)
  }
  if (sd(x[ok]) == 0 || sd(y[ok]) == 0) {
    return(NA_real_)
  }
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
    mean_cor = mean_cor,
    mean_zrmse = mean_zrmse,
    mapping_score = score
  )
}

fit_initial_window_model <- function(df, raw_col, ref_col, calibration_hours_value = 48) {
  usable <- df |>
    filter(complete.cases(.data[[raw_col]], .data[[ref_col]])) |>
    arrange(DATE.HOUR)

  calibration_df <- usable |>
    slice_head(n = calibration_hours_value)

  predicted_all <- rep(NA_real_, nrow(df))
  model_type <- "none"
  intercept <- NA_real_
  slope <- NA_real_

  if (nrow(calibration_df) >= 4 &&
      sd(calibration_df[[raw_col]]) > 0 &&
      sd(calibration_df[[ref_col]]) > 0) {
    fit <- lm(calibration_df[[ref_col]] ~ calibration_df[[raw_col]])
    intercept <- unname(coef(fit)[1])
    slope <- unname(coef(fit)[2])
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

  list(
    predicted = predicted_all,
    stats = tibble(
      model_type = model_type,
      n_calibration_hours = nrow(calibration_df),
      calibration_start = if (nrow(calibration_df) > 0) min(calibration_df$DATE.HOUR) else as.POSIXct(NA, tz = timezone_local),
      calibration_end = if (nrow(calibration_df) > 0) max(calibration_df$DATE.HOUR) + hours(1) else as.POSIXct(NA, tz = timezone_local),
      intercept = intercept,
      slope = slope
    )
  )
}

plot_otice_node <- function(df, gas_name, title_text) {
  if (gas_name == "CO2") {
    plot_df <- df |>
      select(DATE.HOUR, OTICE_raw = OTICE_CO2_raw, OTICE_fitted = OTICE_CO2_fitted, CRDS = CRDS_CO2) |>
      pivot_longer(cols = -DATE.HOUR, names_to = "series", values_to = "value")
    line_colors <- c("OTICE_raw" = "#7DB7E8", "OTICE_fitted" = "#145A96", "CRDS" = "#3F3F46")
    y_label <- "CO2 (ppm)"
  } else {
    plot_df <- df |>
      select(DATE.HOUR, OTICE_raw = OTICE_NH3_raw, OTICE_fitted = OTICE_NH3_fitted, CRDS = CRDS_NH3) |>
      pivot_longer(cols = -DATE.HOUR, names_to = "series", values_to = "value")
    line_colors <- c("OTICE_raw" = "#8BD38B", "OTICE_fitted" = "#1F8A4C", "CRDS" = "#3F3F46")
    y_label <- "NH3 (ppm)"
  }

  ggplot(plot_df, aes(x = DATE.HOUR, y = value, color = series)) +
    geom_line(linewidth = 0.6, na.rm = TRUE) +
    scale_color_manual(values = line_colors) +
    labs(title = title_text, x = NULL, y = y_label, color = NULL) +
    theme_classic(base_size = 15) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      legend.position = "bottom",
      panel.border = element_rect(color = "black", fill = NA)
    )
}

otice_files <- list.files(
  path = clean_dir,
  pattern = "^min_calibrated_otice_data_.*\\.csv$",
  full.names = TRUE
)

otice_hourly_nodes <- map_dfr(
  otice_files,
  ~ read_csv(.x, col_types = cols(.default = col_character()))
) |>
  mutate(
    Datetime_Berlin = parse_datetime_local(Datetime_Berlin),
    DATE.HOUR = floor_date(Datetime_Berlin, "hour"),
    node = clean_otice_node(Node),
    OTICE_CO2_raw = as.numeric(CO2.AVG),
    OTICE_NH3_raw = as.numeric(NH3_ppm)
  ) |>
  filter(!is.na(DATE.HOUR), !is.na(node), node != "") |>
  group_by(DATE.HOUR, node) |>
  summarise(
    OTICE_CO2_raw = safe_mean(OTICE_CO2_raw),
    OTICE_NH3_raw = safe_mean(OTICE_NH3_raw),
    otice_points = n(),
    .groups = "drop"
  ) |>
  add_period_id() |>
  arrange(period_id, node, DATE.HOUR)

write_csv(otice_hourly_nodes, file.path(clean_dir, "otice_node_hourly_raw_all.csv"))

crds_step_files <- list.files(
  path = file.path(project_dir, "workflows", "crds_routine_cleaning", "clean_data", "crds_clean", "step_avg"),
  pattern = "\\.csv$",
  full.names = TRUE
)

crds_hourly_refs <- bind_rows(lapply(crds_step_files, function(file_path) {
  read_csv(file_path, col_types = cols(.default = col_character()))
})) |>
  transmute(
    DATE.TIME = parse_datetime_local(DATE.TIME),
    DATE.HOUR = floor_date(DATE.TIME, "hour"),
    analyzer = toupper(trimws(as.character(analyzer))),
    crds_location = clean_crds_location(location),
    CRDS_NH3 = as.numeric(NH3),
    CRDS_CO2 = as.numeric(CO2)
  ) |>
  filter(!is.na(DATE.HOUR), !is.na(analyzer), !is.na(crds_location)) |>
  group_by(DATE.HOUR, analyzer, crds_location) |>
  summarise(
    CRDS_NH3 = safe_mean(CRDS_NH3),
    CRDS_CO2 = safe_mean(CRDS_CO2),
    crds_points = n(),
    .groups = "drop"
  ) |>
  add_period_id()

mapping_rows <- list()
mapping_index <- 1

for (period_id_value in period_lookup$period_id) {
  otice_period <- otice_hourly_nodes |>
    filter(period_id == period_id_value)
  crds_period <- crds_hourly_refs |>
    filter(period_id == period_id_value)

  if (nrow(otice_period) == 0 || nrow(crds_period) == 0) {
    next
  }

  candidate_refs <- crds_period |>
    distinct(analyzer, crds_location)

  for (node_value in sort(unique(otice_period$node))) {
    node_hourly <- otice_period |>
      filter(node == node_value)

    score_table <- bind_rows(lapply(seq_len(nrow(candidate_refs)), function(i) {
      candidate <- candidate_refs[i, ]

      joined_hourly <- node_hourly |>
        inner_join(
          crds_period |>
            filter(analyzer == candidate$analyzer, crds_location == candidate$crds_location) |>
            select(period_id, DATE.HOUR, analyzer, crds_location, CRDS_NH3, CRDS_CO2),
          by = c("period_id", "DATE.HOUR")
        )

      if (nrow(joined_hourly) < 2) {
        return(NULL)
      }

      bind_cols(
        tibble(period_id = period_id_value, node = node_value, analyzer = candidate$analyzer, crds_location = candidate$crds_location),
        build_mapping_score(joined_hourly)
      )
    }))

    if (nrow(score_table) == 0) {
      next
    }

    mapping_rows[[mapping_index]] <- score_table |>
      arrange(desc(mapping_score), desc(overlap_hours), desc(gas_count), desc(mean_cor), mean_zrmse) |>
      slice(1)
    mapping_index <- mapping_index + 1
  }
}

reference_mapping <- bind_rows(mapping_rows) |>
  arrange(period_id, as.numeric(node))

write_csv(reference_mapping, file.path(clean_dir, "otice_reference_mapping.csv"))

otice_matched_hourly <- otice_hourly_nodes |>
  inner_join(reference_mapping, by = c("period_id", "node")) |>
  inner_join(
    crds_hourly_refs |>
      select(period_id, DATE.HOUR, analyzer, crds_location, CRDS_NH3, CRDS_CO2, crds_points),
    by = c("period_id", "DATE.HOUR", "analyzer", "crds_location")
  ) |>
  arrange(period_id, node, DATE.HOUR)

node_groups <- otice_matched_hourly |>
  distinct(period_id, node, analyzer, crds_location) |>
  arrange(period_id, as.numeric(node))

calibrated_groups <- list()
model_stats <- list()

for (i in seq_len(nrow(node_groups))) {
  group_key <- node_groups[i, ]

  group_df <- otice_matched_hourly |>
    filter(
      period_id == group_key$period_id,
      node == group_key$node,
      analyzer == group_key$analyzer,
      crds_location == group_key$crds_location
    ) |>
    arrange(DATE.HOUR)

  co2_fit <- fit_initial_window_model(group_df, "OTICE_CO2_raw", "CRDS_CO2", calibration_hours)
  nh3_fit <- fit_initial_window_model(group_df, "OTICE_NH3_raw", "CRDS_NH3", calibration_hours)

  group_df$OTICE_CO2_fitted <- co2_fit$predicted
  group_df$OTICE_NH3_fitted <- nh3_fit$predicted

  calibrated_groups[[i]] <- group_df
  model_stats[[length(model_stats) + 1]] <- bind_cols(group_key, tibble(gas = "CO2"), co2_fit$stats)
  model_stats[[length(model_stats) + 1]] <- bind_cols(group_key, tibble(gas = "NH3"), nh3_fit$stats)
}

otice_node_hourly_calibrated <- bind_rows(calibrated_groups) |>
  arrange(period_id, as.numeric(node), DATE.HOUR)

otice_node_model_stats <- bind_rows(model_stats) |>
  arrange(period_id, gas, as.numeric(node))

write_csv(otice_node_hourly_calibrated, file.path(clean_dir, "otice_node_hourly_calibrated_all.csv"))
write_csv(otice_node_model_stats, file.path(clean_dir, "otice_node_model_stats.csv"))

for (i in seq_len(nrow(period_lookup))) {
  period_row <- period_lookup[i, ]
  period_nodes_raw <- otice_hourly_nodes |>
    filter(period_id == period_row$period_id)
  period_nodes_fit <- otice_node_hourly_calibrated |>
    filter(period_id == period_row$period_id)

  if (nrow(period_nodes_raw) > 0) {
    write_csv(period_nodes_raw, file.path(clean_dir, paste0("otice_node_hourly_raw_", period_row$file_tag, ".csv")))
  }

  if (nrow(period_nodes_fit) > 0) {
    write_csv(period_nodes_fit, file.path(clean_dir, paste0("otice_node_hourly_calibrated_", period_row$file_tag, ".csv")))
  }
}

otice_hourly_inside <- otice_node_hourly_calibrated |>
  group_by(period_id, DATE.HOUR) |>
  summarise(
    CO2_in_raw = safe_mean(OTICE_CO2_raw),
    CO2_in_fitted = safe_mean(OTICE_CO2_fitted),
    NH3_in_raw = safe_mean(OTICE_NH3_raw),
    NH3_in_fitted = safe_mean(OTICE_NH3_fitted),
    n_otice_nodes = n_distinct(node[!is.na(OTICE_CO2_fitted) | !is.na(OTICE_NH3_fitted)]),
    .groups = "drop"
  ) |>
  arrange(period_id, DATE.HOUR)

write_csv(otice_hourly_inside, file.path(clean_dir, "otice_hourly_inside_all.csv"))

crds_outside_files <- list.files(
  file.path(project_dir, "workflows", "crds_routine_cleaning", "clean_data", "crds_clean", "hourly_in_out_avg"),
  pattern = "\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)

crds_outside_hourly <- map_dfr(
  crds_outside_files,
  ~ read_csv(.x, col_types = cols(.default = col_character()))
) |>
  mutate(
    DATE.HOUR = parse_datetime_local(DATE.HOUR),
    CO2_S = as.numeric(CO2_S),
    NH3_S = as.numeric(NH3_S)
  ) |>
  filter(!is.na(DATE.HOUR)) |>
  group_by(DATE.HOUR) |>
  summarise(
    CO2_S = safe_mean(CO2_S),
    NH3_S = safe_mean(NH3_S),
    .groups = "drop"
  ) |>
  add_period_id()

otice_hourly_for_comparison <- otice_hourly_inside |>
  left_join(crds_outside_hourly, by = c("period_id", "DATE.HOUR")) |>
  mutate(
    analyzer = "otice",
    CH4_in = NA_real_,
    CH4_S = NA_real_,
    CO2_in = CO2_in_fitted,
    NH3_in = NH3_in_fitted,
    delta_CO2_raw = CO2_in_raw - CO2_S,
    delta_NH3_raw = NH3_in_raw - NH3_S,
    delta_CO2_fitted = CO2_in_fitted - CO2_S,
    delta_NH3_fitted = NH3_in_fitted - NH3_S,
    delta_CO2 = delta_CO2_fitted,
    delta_CH4 = NA_real_,
    delta_NH3 = delta_NH3_fitted
  ) |>
  arrange(period_id, DATE.HOUR)

write_csv(otice_hourly_for_comparison, file.path(clean_dir, "otice_hourly_for_comparison_all.csv"))

for (i in seq_len(nrow(period_lookup))) {
  period_row <- period_lookup[i, ]
  inside_df <- otice_hourly_inside |>
    filter(period_id == period_row$period_id)
  compare_df <- otice_hourly_for_comparison |>
    filter(period_id == period_row$period_id)

  if (nrow(inside_df) > 0) {
    write_csv(inside_df, file.path(clean_dir, paste0("otice_hourly_inside_", period_row$file_tag, ".csv")))
  }

  if (nrow(compare_df) > 0) {
    write_csv(compare_df, file.path(clean_dir, paste0("otice_hourly_for_comparison_", period_row$file_tag, ".csv")))
  }
}

for (i in seq_len(nrow(period_lookup))) {
  period_row <- period_lookup[i, ]
  period_df <- otice_node_hourly_calibrated |>
    filter(period_id == period_row$period_id)

  if (nrow(period_df) == 0) {
    next
  }

  for (node_value in sort(unique(period_df$node))) {
    node_df <- period_df |>
      filter(node == node_value)

    co2_plot <- plot_otice_node(
      node_df,
      gas_name = "CO2",
      title_text = paste("OTICE node O", node_value, "CO2", period_row$file_tag)
    )

    nh3_plot <- plot_otice_node(
      node_df,
      gas_name = "NH3",
      title_text = paste("OTICE node O", node_value, "NH3", period_row$file_tag)
    )

    ggsave(
      file.path(plot_dir, period_row$folder, paste0("otice_node_O", node_value, "_CO2_", period_row$file_tag, ".png")),
      co2_plot,
      width = 10,
      height = 5.8,
      dpi = 150
    )

    ggsave(
      file.path(plot_dir, period_row$folder, paste0("otice_node_O", node_value, "_NH3_", period_row$file_tag, ".png")),
      nh3_plot,
      width = 10,
      height = 5.8,
      dpi = 150
    )
  }
}

cat("Wrote cleaned otice outputs to:\n", normalizePath(clean_dir), "\n")
