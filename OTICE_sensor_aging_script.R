# =============================================================================
# OTICE versus CRDS concentration comparison
# Fresh daily-resolution workflow for September to December 2025
# 1. import raw processed data
# 2. build daily averages
# 3. map each OTICE node to the best CRDS reference for each month
# 4. fit simple monthly calibrations
# 5. save tables and daily plots
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

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args_all[grepl(file_arg, args_all)])
base_dir <- if (length(script_path) > 0) dirname(normalizePath(script_path[1])) else getwd()

otice_dir <- file.path(base_dir, "processed_data", "OTICE_processed")
crds_dir <- file.path(base_dir, "processed_data", "CRDS8_processed")
out_dir <- file.path(base_dir, "output", "OTICE_versus_CRDS_concentration")
plots_daily_dir <- file.path(out_dir, "plots", "daily")
tables_dir <- file.path(out_dir, "tables")

if (dir.exists(out_dir)) {
  unlink(out_dir, recursive = TRUE, force = TRUE)
}

dir.create(plots_daily_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

month_lookup <- tibble(
  month_id = c("2025-09", "2025-10", "2025-11", "2025-12"),
  month_name = c("September", "October", "November", "December"),
  month_label = c("September 2025", "October 2025", "November 2025", "December 2025"),
  month_start = as.Date(c("2025-09-01", "2025-10-01", "2025-11-01", "2025-12-01")),
  month_end = as.Date(c("2025-10-01", "2025-11-01", "2025-12-01", "2026-01-01")),
  make_plot = c(FALSE, TRUE, TRUE, TRUE)
)

# -----------------------------------------------------------------------------
# 2. Small helper functions
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

safe_cor <- function(x, y, min_points = 4) {
  ok <- complete.cases(x, y)
  if (sum(ok) < min_points) return(NA_real_)
  if (sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  cor(x[ok], y[ok])
}

safe_zrmse <- function(x, y, min_points = 4) {
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

assign_month_id <- function(date_value) {
  case_when(
    date_value >= as.Date("2025-09-01") & date_value < as.Date("2025-10-01") ~ "2025-09",
    date_value >= as.Date("2025-10-01") & date_value < as.Date("2025-11-01") ~ "2025-10",
    date_value >= as.Date("2025-11-01") & date_value < as.Date("2025-12-01") ~ "2025-11",
    date_value >= as.Date("2025-12-01") & date_value < as.Date("2026-01-01") ~ "2025-12",
    TRUE ~ NA_character_
  )
}

fmt_num <- function(x, digits = 2) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
}

month_file_prefix <- function(month_name_value) {
  gsub(" ", "_", month_name_value)
}

fit_month_model <- function(df, raw_col, ref_col) {
  valid <- complete.cases(df[[raw_col]], df[[ref_col]])
  n_points <- sum(valid)

  predicted_all <- rep(NA_real_, nrow(df))
  model_type <- "none"
  intercept <- NA_real_
  slope <- NA_real_
  r_squared <- NA_real_

  if (n_points >= 4 &&
      sd(df[[raw_col]][valid]) > 0 &&
      sd(df[[ref_col]][valid]) > 0) {
    fit <- lm(df[[ref_col]][valid] ~ df[[raw_col]][valid])
    intercept <- unname(coef(fit)[1])
    slope <- unname(coef(fit)[2])
    r_squared <- summary(fit)$r.squared
    predicted_all <- intercept + slope * df[[raw_col]]
    model_type <- "linear"
  } else if (n_points >= 2) {
    intercept <- mean(df[[ref_col]][valid] - df[[raw_col]][valid])
    slope <- 1
    predicted_all <- df[[raw_col]] + intercept
    model_type <- "offset_only"
  } else if (n_points >= 1) {
    intercept <- 0
    slope <- 1
    predicted_all <- df[[raw_col]]
    model_type <- "identity"
  }

  raw_rmse <- safe_rmse(df[[raw_col]], df[[ref_col]])
  fitted_rmse <- safe_rmse(predicted_all, df[[ref_col]])
  raw_mae <- safe_mae(df[[raw_col]], df[[ref_col]])
  fitted_mae <- safe_mae(predicted_all, df[[ref_col]])

  list(
    predicted = predicted_all,
    stats = tibble(
      model_type = model_type,
      n_points = n_points,
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

build_mapping_score <- function(joined_daily) {
  overlap_days <- nrow(joined_daily)
  cor_co2 <- safe_cor(joined_daily$OTICE_CO2_raw, joined_daily$CRDS_CO2)
  cor_nh3 <- safe_cor(joined_daily$OTICE_NH3_raw, joined_daily$CRDS_NH3)
  zrmse_co2 <- safe_zrmse(joined_daily$OTICE_CO2_raw, joined_daily$CRDS_CO2)
  zrmse_nh3 <- safe_zrmse(joined_daily$OTICE_NH3_raw, joined_daily$CRDS_NH3)

  mean_cor <- if (all(is.na(c(cor_co2, cor_nh3)))) NA_real_ else mean(c(cor_co2, cor_nh3), na.rm = TRUE)
  mean_zrmse <- if (all(is.na(c(zrmse_co2, zrmse_nh3)))) NA_real_ else mean(c(zrmse_co2, zrmse_nh3), na.rm = TRUE)
  gas_count <- sum(!is.na(c(cor_co2, cor_nh3)))

  score <- if (is.na(mean_cor)) {
    -Inf
  } else {
    mean_cor - 0.20 * ifelse(is.na(mean_zrmse), 2, mean_zrmse) + 0.02 * overlap_days + 0.05 * gas_count
  }

  tibble(
    overlap_days = overlap_days,
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

gas_style <- function(gas_label) {
  if (gas_label == "CO2") {
    list(
      y_label = "CO2 (ppm)",
      ref_color = "#3F3F46",
      raw_color = "#7DB7E8",
      fit_color = "#145A96",
      accuracy = 1
    )
  } else {
    list(
      y_label = "NH3 (ppm)",
      ref_color = "#3F3F46",
      raw_color = "#8BD38B",
      fit_color = "#1F8A4C",
      accuracy = 0.01
    )
  }
}

build_model_text <- function(model_row) {
  if (nrow(model_row) == 0) return("No model statistics available.")

  eq_text <- if (is.na(model_row$slope[1]) || is.na(model_row$intercept[1])) {
    "fit: not available"
  } else if (model_row$model_type[1] == "offset_only") {
    paste0("fit: y = x + ", fmt_num(model_row$intercept[1], 3))
  } else {
    paste0(
      "fit: y = ",
      fmt_num(model_row$slope[1], 3),
      "x + ",
      fmt_num(model_row$intercept[1], 3)
    )
  }

  paste(
    paste0("reference: ", model_row$analyzer[1], " location ", model_row$crds_location[1]),
    paste0("model: ", model_row$model_type[1]),
    paste0("n = ", model_row$n_points[1]),
    paste0("R2 = ", fmt_num(model_row$r_squared[1], 3)),
    paste0("RMSE raw ", fmt_num(model_row$raw_rmse[1], 3), " -> ", fmt_num(model_row$fitted_rmse[1], 3)),
    eq_text,
    sep = " | "
  )
}

make_node_daily_plot <- function(df, model_row, month_label_value, month_name_value, node_value, gas_label) {
  style <- gas_style(gas_label)

  if (gas_label == "CO2") {
    plot_df <- df |>
      select(date_day, CRDS = CRDS_CO2, OTICE_raw = OTICE_CO2_raw, OTICE_fitted = OTICE_CO2_fitted) |>
      pivot_longer(cols = c(CRDS, OTICE_raw, OTICE_fitted), names_to = "series", values_to = "value")
  } else {
    plot_df <- df |>
      select(date_day, CRDS = CRDS_NH3, OTICE_raw = OTICE_NH3_raw, OTICE_fitted = OTICE_NH3_fitted) |>
      pivot_longer(cols = c(CRDS, OTICE_raw, OTICE_fitted), names_to = "series", values_to = "value")
  }

  color_values <- c(
    CRDS = style$ref_color,
    OTICE_raw = style$raw_color,
    OTICE_fitted = style$fit_color
  )

  label_values <- c(
    CRDS = "CRDS reference",
    OTICE_raw = "OTICE raw",
    OTICE_fitted = "OTICE fitted"
  )

  ggplot(plot_df, aes(x = date_day, y = value, color = series)) +
    geom_line(linewidth = 0.9, na.rm = TRUE) +
    geom_point(size = 2, alpha = 0.7, na.rm = TRUE) +
    scale_color_manual(values = color_values, labels = label_values) +
    scale_x_date(date_breaks = "2 days", date_labels = "%d %b", expand = expansion(mult = c(0.01, 0.02))) +
    scale_y_continuous(labels = label_number(accuracy = style$accuracy), n.breaks = 12) +
    labs(
      title = paste0(month_label_value, " ", gas_label, ": OTICE node ", node_value),
      subtitle = build_model_text(model_row),
      x = NULL,
      y = style$y_label,
      color = NULL,
      caption = paste0(
        "Fresh monthly mapping was chosen from CRDS analyzer-location candidates using daily CO2/NH3 similarity. ",
        "Plot file month: ", month_name_value, "."
      )
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

make_average_daily_plot <- function(df, average_stats_row, reference_mix_text, month_label_value, month_name_value, gas_label) {
  style <- gas_style(gas_label)

  if (gas_label == "CO2") {
    plot_df <- df |>
      select(date_day, CRDS = CRDS_CO2_mean, OTICE_raw = OTICE_CO2_raw_mean, OTICE_fitted = OTICE_CO2_fitted_mean) |>
      pivot_longer(cols = c(CRDS, OTICE_raw, OTICE_fitted), names_to = "series", values_to = "value")
  } else {
    plot_df <- df |>
      select(date_day, CRDS = CRDS_NH3_mean, OTICE_raw = OTICE_NH3_raw_mean, OTICE_fitted = OTICE_NH3_fitted_mean) |>
      pivot_longer(cols = c(CRDS, OTICE_raw, OTICE_fitted), names_to = "series", values_to = "value")
  }

  color_values <- c(
    CRDS = style$ref_color,
    OTICE_raw = style$raw_color,
    OTICE_fitted = style$fit_color
  )

  label_values <- c(
    CRDS = "CRDS reference mean",
    OTICE_raw = "OTICE raw mean",
    OTICE_fitted = "OTICE fitted mean"
  )

  stats_text <- paste(
    paste0("days = ", average_stats_row$n_points[1]),
    paste0("RMSE raw ", fmt_num(average_stats_row$raw_rmse[1], 3), " -> ", fmt_num(average_stats_row$fitted_rmse[1], 3)),
    paste0("MAE raw ", fmt_num(average_stats_row$raw_mae[1], 3), " -> ", fmt_num(average_stats_row$fitted_mae[1], 3)),
    paste0("mapped nodes = ", average_stats_row$n_nodes[1]),
    sep = " | "
  )

  ggplot(plot_df, aes(x = date_day, y = value, color = series)) +
    geom_line(linewidth = 0.9, na.rm = TRUE) +
    geom_point(size = 2, alpha = 0.7, na.rm = TRUE) +
    scale_color_manual(values = color_values, labels = label_values) +
    scale_x_date(date_breaks = "2 days", date_labels = "%d %b", expand = expansion(mult = c(0.01, 0.02))) +
    scale_y_continuous(labels = label_number(accuracy = style$accuracy), n.breaks = 12) +
    labs(
      title = paste0(month_label_value, " ", gas_label, ": average of mapped OTICE nodes versus average CRDS"),
      subtitle = stats_text,
      x = NULL,
      y = style$y_label,
      color = NULL,
      caption = paste0(
        "Displayed fitted mean is the average of node-level fitted values. ",
        "Reference mix: ", reference_mix_text, ". Plot file month: ", month_name_value, "."
      )
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

# -----------------------------------------------------------------------------
# 3. Import OTICE data and build daily averages
# -----------------------------------------------------------------------------
otice_files <- list.files(otice_dir, pattern = "^min_calibrated.*\\.csv$", full.names = TRUE)
otice_files <- otice_files[file.info(otice_files)$size > 0]

OTICE_dataset <- bind_rows(lapply(otice_files, function(file_path) {
  read.csv(file_path, stringsAsFactors = FALSE)
})) |>
  transmute(
    datetime = as.POSIXct(.data[["Datetime_Berlin"]], format = "%Y-%m-%d %H:%M:%S", tz = timezone_local),
    date_day = as.Date(datetime, tz = timezone_local),
    month_id = assign_month_id(date_day),
    node = clean_otice_node(.data[["Node"]]),
    OTICE_NH3_raw = suppressWarnings(as.numeric(.data[["NH3_ppm"]])),
    OTICE_CO2_raw = suppressWarnings(as.numeric(.data[["CO2.RAW"]])),
    OTICE_NH3_internal = suppressWarnings(as.numeric(.data[["NH3_ppm_barn"]])),
    OTICE_CO2_internal = suppressWarnings(as.numeric(.data[["CO2.AVG_barn"]]))
  ) |>
  filter(!is.na(datetime), !is.na(month_id), !is.na(node))

otice_daily <- OTICE_dataset |>
  group_by(month_id, date_day, node) |>
  summarise(
    OTICE_NH3_raw = safe_mean(OTICE_NH3_raw),
    OTICE_CO2_raw = safe_mean(OTICE_CO2_raw),
    OTICE_NH3_internal = safe_mean(OTICE_NH3_internal),
    OTICE_CO2_internal = safe_mean(OTICE_CO2_internal),
    otice_points = n(),
    .groups = "drop"
  ) |>
  left_join(month_lookup |> select(month_id, month_name, month_label), by = "month_id")

# -----------------------------------------------------------------------------
# 4. Import CRDS data and build daily averages
# -----------------------------------------------------------------------------
crds_files <- list.files(crds_dir, pattern = "\\.csv$", full.names = TRUE)
crds_files <- crds_files[file.info(crds_files)$size > 0]

CRDS_dataset <- bind_rows(lapply(crds_files, function(file_path) {
  df <- read.csv(file_path, stringsAsFactors = FALSE)
  df$location <- as.character(df$location)
  df$analyzer <- as.character(df$analyzer)
  df
})) |>
  transmute(
    datetime = as.POSIXct(.data[["DATE.TIME"]], format = "%Y-%m-%d %H:%M:%S", tz = timezone_local),
    date_day = as.Date(datetime, tz = timezone_local),
    month_id = assign_month_id(date_day),
    analyzer = toupper(trimws(as.character(.data[["analyzer"]]))),
    crds_location = clean_crds_location(.data[["location"]]),
    CRDS_NH3 = suppressWarnings(as.numeric(.data[["NH3"]])),
    CRDS_CO2 = suppressWarnings(as.numeric(.data[["CO2"]]))
  ) |>
  filter(!is.na(datetime), !is.na(month_id), !is.na(analyzer), !is.na(crds_location))

crds_daily <- CRDS_dataset |>
  group_by(month_id, date_day, analyzer, crds_location) |>
  summarise(
    CRDS_NH3 = safe_mean(CRDS_NH3),
    CRDS_CO2 = safe_mean(CRDS_CO2),
    crds_points = n(),
    .groups = "drop"
  ) |>
  left_join(month_lookup |> select(month_id, month_name, month_label), by = "month_id")

# -----------------------------------------------------------------------------
# 5. Fresh monthly mapping: choose the best CRDS reference for each node
# -----------------------------------------------------------------------------
mapping_rows <- list()
mapping_score_rows <- list()
mapping_index <- 1

for (month_id_value in month_lookup$month_id) {
  month_otice <- otice_daily |>
    filter(month_id == month_id_value)

  month_crds <- crds_daily |>
    filter(month_id == month_id_value)

  if (nrow(month_otice) == 0 || nrow(month_crds) == 0) {
    next
  }

  candidate_refs <- month_crds |>
    distinct(analyzer, crds_location)

  for (node_value in sort(unique(month_otice$node))) {
    node_daily <- month_otice |>
      filter(node == node_value)

    score_table <- bind_rows(lapply(seq_len(nrow(candidate_refs)), function(i) {
      candidate <- candidate_refs[i, ]

      joined_daily <- node_daily |>
        select(date_day, month_id, node, OTICE_NH3_raw, OTICE_CO2_raw, OTICE_NH3_internal, OTICE_CO2_internal) |>
        inner_join(
          month_crds |>
            filter(analyzer == candidate$analyzer, crds_location == candidate$crds_location) |>
            select(date_day, analyzer, crds_location, CRDS_NH3, CRDS_CO2),
          by = "date_day"
        )

      if (nrow(joined_daily) < 2) {
        return(NULL)
      }

      score_row <- build_mapping_score(joined_daily)

      bind_cols(
        tibble(
          month_id = month_id_value,
          node = node_value,
          analyzer = candidate$analyzer,
          crds_location = candidate$crds_location
        ),
        score_row
      )
    }))

    if (nrow(score_table) == 0) {
      next
    }

    best_row <- score_table |>
      arrange(desc(mapping_score), desc(overlap_days), desc(gas_count), desc(mean_cor), mean_zrmse) |>
      slice(1)

    best_row <- best_row |>
      left_join(month_lookup |> select(month_id, month_name, month_label), by = "month_id")

    mapping_rows[[mapping_index]] <- best_row
    mapping_score_rows[[mapping_index]] <- score_table |>
      left_join(month_lookup |> select(month_id, month_name, month_label), by = "month_id")
    mapping_index <- mapping_index + 1
  }
}

reference_mapping <- bind_rows(mapping_rows) |>
  relocate(month_id, month_name, month_label, node, analyzer, crds_location)

reference_mapping_scores <- bind_rows(mapping_score_rows) |>
  relocate(month_id, month_name, month_label, node, analyzer, crds_location)

# -----------------------------------------------------------------------------
# 6. Build matched daily node-versus-reference table
# -----------------------------------------------------------------------------
matched_daily <- otice_daily |>
  inner_join(
    reference_mapping |>
      select(month_id, node, analyzer, crds_location, overlap_days, gas_count, cor_co2, cor_nh3, mean_cor, mean_zrmse, mapping_score),
    by = c("month_id", "node")
  ) |>
  inner_join(
    crds_daily |>
      select(month_id, date_day, analyzer, crds_location, CRDS_NH3, CRDS_CO2, crds_points),
    by = c("month_id", "date_day", "analyzer", "crds_location")
  ) |>
  arrange(month_id, node, date_day)

# -----------------------------------------------------------------------------
# 7. Fit fresh monthly calibrations for each node
# -----------------------------------------------------------------------------
group_keys <- matched_daily |>
  distinct(month_id, month_name, month_label, node, analyzer, crds_location) |>
  arrange(month_id, as.numeric(node))

predicted_groups <- list()
co2_model_rows <- list()
nh3_model_rows <- list()

for (i in seq_len(nrow(group_keys))) {
  group_info <- group_keys[i, ]

  group_df <- matched_daily |>
    filter(
      month_id == group_info$month_id,
      node == group_info$node,
      analyzer == group_info$analyzer,
      crds_location == group_info$crds_location
    ) |>
    arrange(date_day)

  co2_fit <- fit_month_model(group_df, "OTICE_CO2_raw", "CRDS_CO2")
  nh3_fit <- fit_month_model(group_df, "OTICE_NH3_raw", "CRDS_NH3")

  group_df$OTICE_CO2_fitted <- co2_fit$predicted
  group_df$OTICE_NH3_fitted <- nh3_fit$predicted

  predicted_groups[[i]] <- group_df

  co2_model_rows[[i]] <- bind_cols(
    group_info,
    tibble(gas = "CO2"),
    co2_fit$stats
  )

  nh3_model_rows[[i]] <- bind_cols(
    group_info,
    tibble(gas = "NH3"),
    nh3_fit$stats
  )
}

node_daily_comparison <- bind_rows(predicted_groups)
node_model_stats <- bind_rows(co2_model_rows, nh3_model_rows) |>
  arrange(month_id, gas, as.numeric(node))

# -----------------------------------------------------------------------------
# 8. Build monthly node-average daily tables
# -----------------------------------------------------------------------------
average_daily_comparison <- node_daily_comparison |>
  group_by(month_id, month_name, month_label, date_day) |>
  summarise(
    OTICE_NH3_raw_mean = safe_mean(OTICE_NH3_raw),
    OTICE_CO2_raw_mean = safe_mean(OTICE_CO2_raw),
    OTICE_NH3_fitted_mean = safe_mean(OTICE_NH3_fitted),
    OTICE_CO2_fitted_mean = safe_mean(OTICE_CO2_fitted),
    CRDS_NH3_mean = safe_mean(CRDS_NH3),
    CRDS_CO2_mean = safe_mean(CRDS_CO2),
    n_nodes = n_distinct(node),
    .groups = "drop"
  )

average_stats_rows <- list()

for (i in seq_len(nrow(month_lookup))) {
  month_id_value <- month_lookup$month_id[i]
  month_name_value <- month_lookup$month_name[i]
  month_label_value <- month_lookup$month_label[i]

  month_avg <- average_daily_comparison |>
    filter(month_id == month_id_value) |>
    arrange(date_day)

  if (nrow(month_avg) == 0) {
    next
  }

  co2_avg_fit <- fit_month_model(month_avg, "OTICE_CO2_raw_mean", "CRDS_CO2_mean")
  nh3_avg_fit <- fit_month_model(month_avg, "OTICE_NH3_raw_mean", "CRDS_NH3_mean")

  average_stats_rows[[length(average_stats_rows) + 1]] <- bind_cols(
    tibble(
      month_id = month_id_value,
      month_name = month_name_value,
      month_label = month_label_value,
      gas = "CO2",
      n_nodes = max(month_avg$n_nodes, na.rm = TRUE)
    ),
    co2_avg_fit$stats
  )

  average_stats_rows[[length(average_stats_rows) + 1]] <- bind_cols(
    tibble(
      month_id = month_id_value,
      month_name = month_name_value,
      month_label = month_label_value,
      gas = "NH3",
      n_nodes = max(month_avg$n_nodes, na.rm = TRUE)
    ),
    nh3_avg_fit$stats
  )
}

average_model_stats <- bind_rows(average_stats_rows) |>
  arrange(month_id, gas)

# -----------------------------------------------------------------------------
# 9. Save tables
# -----------------------------------------------------------------------------
write.csv(otice_daily, file.path(tables_dir, "otice_daily.csv"), row.names = FALSE)
write.csv(crds_daily, file.path(tables_dir, "crds_daily.csv"), row.names = FALSE)
write.csv(reference_mapping_scores, file.path(tables_dir, "reference_mapping_scores.csv"), row.names = FALSE)
write.csv(reference_mapping, file.path(tables_dir, "reference_mapping_by_month.csv"), row.names = FALSE)
write.csv(node_model_stats, file.path(tables_dir, "node_model_stats.csv"), row.names = FALSE)
write.csv(node_daily_comparison, file.path(tables_dir, "node_daily_comparison.csv"), row.names = FALSE)
write.csv(average_model_stats, file.path(tables_dir, "average_model_stats.csv"), row.names = FALSE)
write.csv(average_daily_comparison, file.path(tables_dir, "average_daily_comparison.csv"), row.names = FALSE)

for (i in seq_len(nrow(month_lookup))) {
  month_id_value <- month_lookup$month_id[i]
  month_name_value <- month_lookup$month_name[i]
  prefix_value <- month_file_prefix(month_name_value)

  write.csv(
    reference_mapping |> filter(month_id == month_id_value),
    file.path(tables_dir, paste0(prefix_value, "_reference_mapping.csv")),
    row.names = FALSE
  )

  write.csv(
    node_model_stats |> filter(month_id == month_id_value),
    file.path(tables_dir, paste0(prefix_value, "_node_model_stats.csv")),
    row.names = FALSE
  )

  write.csv(
    node_daily_comparison |> filter(month_id == month_id_value),
    file.path(tables_dir, paste0(prefix_value, "_node_daily_comparison.csv")),
    row.names = FALSE
  )

  write.csv(
    average_daily_comparison |> filter(month_id == month_id_value),
    file.path(tables_dir, paste0(prefix_value, "_average_daily_comparison.csv")),
    row.names = FALSE
  )
}

# -----------------------------------------------------------------------------
# 10. Make daily plots for October, November and December
# -----------------------------------------------------------------------------
for (i in seq_len(nrow(month_lookup))) {
  month_row <- month_lookup[i, ]

  if (!month_row$make_plot) {
    next
  }

  month_nodes <- node_daily_comparison |>
    filter(month_id == month_row$month_id)

  month_average <- average_daily_comparison |>
    filter(month_id == month_row$month_id)

  month_mapping <- reference_mapping |>
    filter(month_id == month_row$month_id)

  if (nrow(month_nodes) == 0 || nrow(month_average) == 0 || nrow(month_mapping) == 0) {
    next
  }

  reference_mix_text <- month_mapping |>
    count(analyzer, crds_location, name = "nodes_here") |>
    arrange(desc(nodes_here), analyzer, crds_location) |>
    transmute(label = paste0(nodes_here, "x ", analyzer, "-", crds_location)) |>
    pull(label) |>
    paste(collapse = "; ")

  for (gas_label in c("CO2", "NH3")) {
    avg_stats_row <- average_model_stats |>
      filter(month_id == month_row$month_id, gas == gas_label)

    avg_plot <- make_average_daily_plot(
      df = month_average,
      average_stats_row = avg_stats_row,
      reference_mix_text = reference_mix_text,
      month_label_value = month_row$month_label,
      month_name_value = month_row$month_name,
      gas_label = gas_label
    )

    ggsave(
      filename = file.path(
        plots_daily_dir,
        paste0(month_file_prefix(month_row$month_name), "_", gas_label, "_all_nodes_daily.png")
      ),
      plot = avg_plot,
      width = 14,
      height = 8,
      dpi = 150
    )

    for (node_value in sort(unique(month_nodes$node))) {
      node_df <- month_nodes |>
        filter(node == node_value)

      node_model_row <- node_model_stats |>
        filter(month_id == month_row$month_id, gas == gas_label, node == node_value)

      if (nrow(node_df) == 0 || nrow(node_model_row) == 0) {
        next
      }

      node_plot <- make_node_daily_plot(
        df = node_df,
        model_row = node_model_row,
        month_label_value = month_row$month_label,
        month_name_value = month_row$month_name,
        node_value = node_value,
        gas_label = gas_label
      )

      ggsave(
        filename = file.path(
          plots_daily_dir,
          paste0(month_file_prefix(month_row$month_name), "_", gas_label, "_node_", node_value, "_daily.png")
        ),
        plot = node_plot,
        width = 14,
        height = 8,
        dpi = 150
      )
    }
  }
}

# -----------------------------------------------------------------------------
# 11. Short console summary
# -----------------------------------------------------------------------------
cat("\n============================================================\n")
cat("OTICE versus CRDS daily concentration workflow\n")
cat("============================================================\n")
cat("Base directory:", base_dir, "\n")
cat("Output directory:", out_dir, "\n")
cat("OTICE daily rows:", nrow(otice_daily), "\n")
cat("CRDS daily rows:", nrow(crds_daily), "\n")
cat("Mapped node-month references:", nrow(reference_mapping), "\n")
cat("Matched node daily rows:", nrow(node_daily_comparison), "\n")
cat("Average daily rows:", nrow(average_daily_comparison), "\n")
cat("\nFiles written to:\n")
cat("  ", tables_dir, "\n")
cat("  ", plots_daily_dir, "\n")
