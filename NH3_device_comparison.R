# =============================================================================
# NH3 & CO2 Device Performance Comparison: OTICE (per node) vs CRDS
# =============================================================================
# Reads output files from NH3_comparison.R (run that first).
#
# Compares each OTICE node against the CRDS reference for:
#   - NH3_OTICE      vs NH3_CRDS_in  (raw calibrated)
#   - NH3_OTICE_corr vs NH3_CRDS_in  (barn-corrected)
#
# Outputs per comparison type:
#   - Summary statistics table (R², RMSE, MAE, Bias, slope, intercept)
#   - Scatter plots per node
#   - Time series overlay (CRDS + OTICE + OTICE_corr)
#   - Bland-Altman agreement plots
# =============================================================================

# install.packages(c("ggplot2", "patchwork", "broom"))

library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(patchwork)
library(broom)

# -----------------------------------------------------------------------------
# 1. Paths
# -----------------------------------------------------------------------------
project_dir <- "D:/Data_Analysis_R/Gas_Concentration_Emission_R"
output_dir  <- file.path(project_dir, "output")
dir.create(output_dir, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 2. Load paired hourly data (CRDS node-average) from NH3_comparison.R output
# -----------------------------------------------------------------------------
paired_file <- file.path(output_dir, "comparison_hourly.csv")
node_file   <- file.path(output_dir, "OTICE_node_hourly.csv")

for (f in c(paired_file, node_file)) {
  if (!file.exists(f))
    stop("File not found: ", f, "\nRun NH3_comparison.R first.")
}

# CRDS hourly reference (one row per hour)
# Note: distinct() removes the DST clock-change duplicate that arises when
# POSIXct timestamps are serialised to CSV strings and read back.
crds <- read.csv(paired_file, stringsAsFactors = FALSE) |>
  mutate(datetime_hour = as.POSIXct(datetime_hour, format = "%Y-%m-%d %H:%M:%S",
                                    tz = "Europe/Berlin")) |>
  select(datetime_hour, NH3_CRDS_in, CO2_CRDS_in) |>
  distinct(datetime_hour, .keep_all = TRUE)

message("CRDS reference rows: ", nrow(crds))

# Per-node OTICE hourly (one row per node per hour)
otice <- read.csv(node_file, stringsAsFactors = FALSE) |>
  mutate(
    datetime_hour  = as.POSIXct(datetime_hour, format = "%Y-%m-%d %H:%M:%S",
                                tz = "Europe/Berlin"),
    NH3_OTICE      = as.numeric(NH3_OTICE),
    NH3_OTICE_corr = as.numeric(NH3_OTICE_corr),
    CO2_OTICE      = as.numeric(CO2_OTICE),
    CO2_OTICE_corr = as.numeric(CO2_OTICE_corr)
  ) |>
  select(datetime_hour, Node, NH3_OTICE, NH3_OTICE_corr, CO2_OTICE, CO2_OTICE_corr)

message("OTICE per-node rows: ", nrow(otice))

# Join each node to the CRDS reference at the same hour
paired <- otice |>
  left_join(crds, by = "datetime_hour", relationship = "many-to-many") |>
  filter(!is.na(NH3_CRDS_in)) |>
  arrange(Node, datetime_hour)

message("Paired rows: ", nrow(paired))
message("Nodes: ", paste(sort(unique(paired$Node)), collapse = ", "))

# -----------------------------------------------------------------------------
# 3. Per-node statistics
#    calc_stats_node(): generic — sensor_col vs ref_col for one node's data
# -----------------------------------------------------------------------------
calc_stats_node <- function(df, sensor_col, ref_col) {
  y  <- df[[sensor_col]]
  x  <- df[[ref_col]]
  ok <- !is.na(y) & !is.na(x)
  if (sum(ok) < 5) return(NULL)
  fit <- lm(y[ok] ~ x[ok])
  res <- y[ok] - x[ok]
  tibble(
    comparison = paste0(sensor_col, " vs ", ref_col),
    n          = sum(ok),
    r2         = round(summary(fit)$r.squared,  4),
    slope      = round(unname(coef(fit)[2]),     4),
    intercept  = round(unname(coef(fit)[1]),     4),
    RMSE_ppm   = round(sqrt(mean(res^2)),        4),
    MAE_ppm    = round(mean(abs(res)),           4),
    Bias_ppm   = round(mean(res),                4)   # positive = sensor > CRDS
  )
}

stats_table <- paired |>
  group_by(Node) |>
  group_modify(~ bind_rows(
    calc_stats_node(.x, "NH3_OTICE",      "NH3_CRDS_in"),
    calc_stats_node(.x, "NH3_OTICE_corr", "NH3_CRDS_in")
  )) |>
  ungroup() |>
  arrange(comparison, desc(r2))

cat("\n===== NH3 Device Comparison Statistics (per node, OTICE vs CRDS) =====\n")
print(stats_table, n = Inf, width = 120)

write.csv(stats_table,
          file.path(output_dir, "device_comparison_stats.csv"),
          row.names = FALSE)
message("\nStats saved to: ", file.path(output_dir, "device_comparison_stats.csv"))

# -----------------------------------------------------------------------------
# 4. Helper: build scatter plot for one comparison pair
# -----------------------------------------------------------------------------
make_scatter <- function(df, sensor_col, ref_col, stats_tbl) {
  st <- stats_tbl |> filter(comparison == paste0(sensor_col, " vs ", ref_col))
  max_val <- max(c(df[[ref_col]], df[[sensor_col]]), na.rm = TRUE) * 1.05
  min_val <- min(c(df[[ref_col]], df[[sensor_col]]), na.rm = TRUE)

  labels <- st |>
    mutate(label = paste0(
      "R\u00b2=", r2,
      "\nSlope=", slope,
      "\nRMSE=",  RMSE_ppm, " ppm",
      "\nBias=",  Bias_ppm, " ppm"
    ))

  df |>
    left_join(labels |> select(Node, label), by = "Node") |>
    filter(!is.na(.data[[ref_col]]), !is.na(.data[[sensor_col]])) |>
    ggplot(aes(x = .data[[ref_col]], y = .data[[sensor_col]])) +
    geom_point(alpha = 0.35, size = 0.8, color = "steelblue") +
    geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "black", linewidth = 0.5) +
    geom_text(
      data = labels |>
        semi_join(df |> distinct(Node), by = "Node"),
      aes(label = label,
          x = min_val + (max_val - min_val) * 0.05,
          y = max_val * 0.95),
      hjust = 0, vjust = 1, size = 2.5, color = "black", inherit.aes = FALSE
    ) +
    facet_wrap(~ Node, ncol = 4) +
    scale_x_continuous(limits = c(min_val, max_val)) +
    scale_y_continuous(limits = c(min_val, max_val)) +
    labs(
      title    = paste0(sensor_col, " vs ", ref_col, " — Hourly (ppm)"),
      subtitle = "Red: OLS regression  |  Dashed: 1:1 line",
      x        = paste(ref_col, "(ppm) — CRDS reference"),
      y        = paste(sensor_col, "(ppm) — OTICE sensor")
    ) +
    theme_bw(base_size = 10)
}

# -----------------------------------------------------------------------------
# 5. Scatter plots
# -----------------------------------------------------------------------------
p_scatter_raw  <- make_scatter(paired, "NH3_OTICE",      "NH3_CRDS_in", stats_table)
p_scatter_corr <- make_scatter(paired, "NH3_OTICE_corr", "NH3_CRDS_in", stats_table)

ggsave(file.path(output_dir, "scatter_NH3_OTICE_vs_CRDS.png"),
       p_scatter_raw,  width = 14, height = 10, dpi = 150)
ggsave(file.path(output_dir, "scatter_NH3_OTICE_corr_vs_CRDS.png"),
       p_scatter_corr, width = 14, height = 10, dpi = 150)
message("Scatter plots saved.")

# -----------------------------------------------------------------------------
# 6. Time series — CRDS + NH3_OTICE + NH3_OTICE_corr overlaid
# -----------------------------------------------------------------------------
ts_data <- paired |>
  select(datetime_hour, Node, NH3_OTICE, NH3_OTICE_corr) |>
  pivot_longer(c(NH3_OTICE, NH3_OTICE_corr),
               names_to = "source", values_to = "NH3_ppm")

p_ts <- ggplot() +
  geom_line(
    data  = crds |> filter(!is.na(NH3_CRDS_in)),
    aes(x = datetime_hour, y = NH3_CRDS_in, color = "CRDS"),
    linewidth = 0.9
  ) +
  geom_line(
    data  = ts_data |> filter(source == "NH3_OTICE", !is.na(NH3_ppm)),
    aes(x = datetime_hour, y = NH3_ppm, group = Node, color = "OTICE (raw)"),
    alpha = 0.45, linewidth = 0.4
  ) +
  geom_line(
    data  = ts_data |> filter(source == "NH3_OTICE_corr", !is.na(NH3_ppm)),
    aes(x = datetime_hour, y = NH3_ppm, group = Node, color = "OTICE (barn-corr)"),
    alpha = 0.45, linewidth = 0.4
  ) +
  scale_color_manual(values = c(
    "CRDS"              = "black",
    "OTICE (raw)"       = "steelblue",
    "OTICE (barn-corr)" = "tomato"
  )) +
  labs(
    title = "NH3 Time Series — CRDS vs OTICE (raw & barn-corrected)",
    x     = NULL, y = "NH3 (ppm)", color = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top")

ggsave(file.path(output_dir, "timeseries_NH3.png"),
       p_ts, width = 14, height = 5, dpi = 150)
message("Time series saved.")

# -----------------------------------------------------------------------------
# 7. Helper: Bland-Altman plot for one comparison pair
# -----------------------------------------------------------------------------
make_ba_plot <- function(df, sensor_col, ref_col) {
  ba <- df |>
    filter(!is.na(.data[[sensor_col]]), !is.na(.data[[ref_col]])) |>
    mutate(
      mean_val = (.data[[ref_col]] + .data[[sensor_col]]) / 2,
      diff_val = .data[[sensor_col]] - .data[[ref_col]]
    ) |>
    group_by(Node) |>
    mutate(
      mean_diff = mean(diff_val, na.rm = TRUE),
      sd_diff   = sd(diff_val,   na.rm = TRUE),
      loa_upper = mean_diff + 1.96 * sd_diff,
      loa_lower = mean_diff - 1.96 * sd_diff
    ) |>
    ungroup()

  ggplot(ba, aes(x = mean_val, y = diff_val)) +
    geom_point(alpha = 0.3, size = 0.7, color = "steelblue") +
    geom_hline(aes(yintercept = mean_diff),  color = "red",    linetype = "solid") +
    geom_hline(aes(yintercept = loa_upper),  color = "orange", linetype = "dashed") +
    geom_hline(aes(yintercept = loa_lower),  color = "orange", linetype = "dashed") +
    geom_hline(yintercept = 0, color = "black", linetype = "dotted") +
    facet_wrap(~ Node, ncol = 4) +
    labs(
      title    = paste0("Bland-Altman — ", sensor_col, " minus ", ref_col, " (ppm)"),
      subtitle = "Red: mean bias  |  Orange dashed: \u00b11.96 SD limits of agreement",
      x        = paste0("Mean of ", sensor_col, " and ", ref_col, " (ppm)"),
      y        = paste0(sensor_col, " \u2212 ", ref_col, " (ppm)")
    ) +
    theme_bw(base_size = 10)
}

# -----------------------------------------------------------------------------
# 8. Bland-Altman plots
# -----------------------------------------------------------------------------
p_ba_raw  <- make_ba_plot(paired, "NH3_OTICE",      "NH3_CRDS_in")
p_ba_corr <- make_ba_plot(paired, "NH3_OTICE_corr", "NH3_CRDS_in")

ggsave(file.path(output_dir, "bland_altman_NH3_OTICE.png"),
       p_ba_raw,  width = 14, height = 10, dpi = 150)
ggsave(file.path(output_dir, "bland_altman_NH3_OTICE_corr.png"),
       p_ba_corr, width = 14, height = 10, dpi = 150)
message("Bland-Altman plots saved.")

# -----------------------------------------------------------------------------
# 9. Summary ranking
# -----------------------------------------------------------------------------
cat("\n===== SENSOR RANKING — NH3_OTICE vs CRDS (best R\u00b2 first) =====\n")
stats_table |>
  filter(comparison == "NH3_OTICE vs NH3_CRDS_in") |>
  mutate(performance = case_when(
    r2 >= 0.80 & abs(Bias_ppm) < 0.1 ~ "Good",
    r2 >= 0.60                         ~ "Moderate",
    TRUE                               ~ "Poor"
  )) |>
  arrange(desc(r2)) |>
  print(n = Inf)

cat("\n===== SENSOR RANKING — NH3_OTICE_corr vs CRDS (best R\u00b2 first) =====\n")
stats_table |>
  filter(comparison == "NH3_OTICE_corr vs NH3_CRDS_in") |>
  mutate(performance = case_when(
    r2 >= 0.80 & abs(Bias_ppm) < 0.1 ~ "Good",
    r2 >= 0.60                         ~ "Moderate",
    TRUE                               ~ "Poor"
  )) |>
  arrange(desc(r2)) |>
  print(n = Inf)

message("\nAll outputs saved to: ", output_dir)
message("Done.")
