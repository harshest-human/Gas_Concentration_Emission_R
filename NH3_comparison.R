# =============================================================================
# NH3 & CO2 Comparison: CRDS (three location groups) vs OTICE (node average)
# Period: 2025-09-29 to 2025-12-31
# =============================================================================
#
# CRDS location groups (location column normalised to lowercase):
#   "in"                 -> NH3_CRDS_in,    CO2_CRDS_in
#   "s"                  -> NH3_CRDS_s,     CO2_CRDS_s
#   all other locations  -> NH3_CRDS_hires, CO2_CRDS_hires
#
# OTICE nodes excluded: O04, O09, O13
# OTICE column mapping:
#   NH3_ppm      -> NH3_OTICE        (raw calibrated)
#   NH3_ppm_barn -> NH3_OTICE_corr   (barn-corrected)
#   CO2.RAW      -> CO2_OTICE
#   CO2.AVG_barn -> CO2_OTICE_corr
#
# Primary comparisons (OTICE node-avg vs CRDS_in):
#   NH3_OTICE / NH3_OTICE_corr vs NH3_CRDS_in
#   CO2_OTICE / CO2_OTICE_corr vs CO2_CRDS_in
# =============================================================================

library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(patchwork)

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
crds_dir  <- "processed_data/CRDS8_processed"
otice_dir <- "processed_data/OTICE_processed"
out_dir   <- "output"
dir.create(out_dir, showWarnings = FALSE)

study_start    <- as.POSIXct("2025-09-29 00:00:00", tz = "Europe/Berlin")
study_end      <- as.POSIXct("2025-12-31 23:59:59", tz = "Europe/Berlin")
OTICE_EXCLUDE  <- c("O04", "O09", "O13")

# =============================================================================
# 1. LOAD CRDS — all locations
#    location normalised to lowercase: "in", "s", numeric strings
# =============================================================================
crds_files <- list.files(crds_dir, pattern = "\\.csv$", full.names = TRUE)
message("CRDS files found: ", length(crds_files))

crds_raw <- lapply(crds_files, function(f) {
  df <- read.csv(f, stringsAsFactors = FALSE)
  df$location <- tolower(trimws(as.character(df$location)))
  df
}) |>
  bind_rows() |>
  mutate(
    DATE.TIME = as.POSIXct(DATE.TIME, format = "%Y-%m-%d %H:%M:%S", tz = "Europe/Berlin"),
    NH3       = as.numeric(NH3),
    CO2       = as.numeric(CO2)
  ) |>
  filter(!is.na(DATE.TIME), !is.na(NH3)) |>
  filter(DATE.TIME >= study_start, DATE.TIME <= study_end) |>
  distinct(analyzer, DATE.TIME, location, .keep_all = TRUE)

message("CRDS steps after filtering: ", nrow(crds_raw))
message("  Devices:   ", paste(sort(unique(crds_raw$analyzer)), collapse = ", "))
message("  Locations: ", paste(sort(unique(crds_raw$location)), collapse = ", "))

# Helper: hourly mean for a location subset; averages across CRDS devices
# when two instruments overlap in the same hour
make_crds_hourly <- function(df) {
  df |>
    mutate(datetime_hour = floor_date(DATE.TIME, "hour")) |>
    group_by(datetime_hour, analyzer) |>
    summarise(NH3 = mean(NH3, na.rm = TRUE),
              CO2 = mean(CO2, na.rm = TRUE),
              .groups = "drop") |>
    group_by(datetime_hour) |>
    summarise(NH3 = mean(NH3), CO2 = mean(CO2), .groups = "drop")
}

# "in" group also tracks CRDS_device for statistics splitting
crds_hourly_in <- crds_raw |>
  filter(location == "in") |>
  mutate(datetime_hour = floor_date(DATE.TIME, "hour")) |>
  group_by(datetime_hour, analyzer) |>
  summarise(NH3 = mean(NH3, na.rm = TRUE),
            CO2 = mean(CO2, na.rm = TRUE),
            n   = n(),
            .groups = "drop") |>
  group_by(datetime_hour) |>
  summarise(NH3_CRDS_in  = mean(NH3),
            CO2_CRDS_in  = mean(CO2),
            CRDS_device  = paste(sort(unique(analyzer)), collapse = "+"),
            .groups = "drop")

crds_hourly_s <- crds_raw |>
  filter(location == "s") |>
  make_crds_hourly() |>
  rename(NH3_CRDS_s = NH3, CO2_CRDS_s = CO2)

crds_hourly_hires <- crds_raw |>
  filter(!location %in% c("in", "s")) |>
  make_crds_hourly() |>
  rename(NH3_CRDS_hires = NH3, CO2_CRDS_hires = CO2)

message("CRDS 'in'    hourly rows: ", nrow(crds_hourly_in))
message("CRDS 's'     hourly rows: ", nrow(crds_hourly_s))
message("CRDS 'hires' hourly rows: ", nrow(crds_hourly_hires))

# =============================================================================
# 2. LOAD OTICE — exclude nodes O04, O09, O13
# =============================================================================
otice_files <- list.files(otice_dir, pattern = "^min_calibrated.*\\.csv$", full.names = TRUE)
message("\nOTICE files found: ", length(otice_files))

otice_raw <- lapply(otice_files, read.csv, stringsAsFactors = FALSE) |>
  bind_rows() |>
  filter(!Node %in% OTICE_EXCLUDE) |>
  mutate(
    Datetime_Berlin = as.POSIXct(Datetime_Berlin, format = "%Y-%m-%d %H:%M:%S",
                                 tz = "Europe/Berlin"),
    NH3_OTICE       = suppressWarnings(as.numeric(NH3_ppm)),
    NH3_OTICE_corr  = suppressWarnings(as.numeric(NH3_ppm_barn)),
    CO2_OTICE       = suppressWarnings(as.numeric(CO2.RAW)),
    CO2_OTICE_corr  = suppressWarnings(as.numeric(CO2.AVG_barn))
  ) |>
  filter(!is.na(Datetime_Berlin)) |>
  filter(Datetime_Berlin >= study_start, Datetime_Berlin <= study_end) |>
  distinct(Node, Datetime_Berlin, .keep_all = TRUE)

message("OTICE rows after filtering: ", nrow(otice_raw))
message("  Nodes retained: ", paste(sort(unique(otice_raw$Node)), collapse = ", "))

# Nodes active per month
node_summary <- otice_raw |>
  filter(!is.na(NH3_OTICE_corr)) |>
  mutate(month = format(Datetime_Berlin, "%Y-%m")) |>
  group_by(month) |>
  summarise(nodes = paste(sort(unique(Node)), collapse = ", "), .groups = "drop")

cat("\nOTICE nodes active by month (NH3_OTICE_corr available):\n")
print(node_summary, n = Inf)

safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

# Hourly mean per node
otice_node_hourly <- otice_raw |>
  mutate(datetime_hour = floor_date(Datetime_Berlin, "hour")) |>
  group_by(datetime_hour, Node) |>
  summarise(
    NH3_OTICE      = safe_mean(NH3_OTICE),
    NH3_OTICE_corr = safe_mean(NH3_OTICE_corr),
    CO2_OTICE      = safe_mean(CO2_OTICE),
    CO2_OTICE_corr = safe_mean(CO2_OTICE_corr),
    n_minutes      = n(),
    .groups        = "drop"
  )

# Mean across available nodes per hour
otice_hourly <- otice_node_hourly |>
  group_by(datetime_hour) |>
  summarise(
    NH3_OTICE      = safe_mean(NH3_OTICE),
    NH3_OTICE_corr = safe_mean(NH3_OTICE_corr),
    CO2_OTICE      = safe_mean(CO2_OTICE),
    CO2_OTICE_corr = safe_mean(CO2_OTICE_corr),
    .groups        = "drop"
  )

message("OTICE hourly rows: ", nrow(otice_hourly))

# =============================================================================
# 3. JOIN
#    Anchor on CRDS_in hours; left-join S and hires (may be NA);
#    inner-join OTICE (keep only hours where both instruments are present)
# =============================================================================
comparison <- crds_hourly_in |>
  left_join(crds_hourly_s,     by = "datetime_hour") |>
  left_join(crds_hourly_hires, by = "datetime_hour") |>
  inner_join(otice_hourly,     by = "datetime_hour") |>
  arrange(datetime_hour) |>
  select(datetime_hour,
         NH3_CRDS_in, NH3_CRDS_s, NH3_CRDS_hires,
         CO2_CRDS_in, CO2_CRDS_s, CO2_CRDS_hires,
         NH3_OTICE, NH3_OTICE_corr,
         CO2_OTICE, CO2_OTICE_corr,
         CRDS_device)   # kept internally for stats; excluded from CSV

message("\nPaired hourly rows: ", nrow(comparison))
message("  Range: ", min(comparison$datetime_hour), " to ", max(comparison$datetime_hour))

# =============================================================================
# 4. STATISTICS (primary reference: CRDS_in)
# =============================================================================
calc_stats <- function(df, sensor_col, ref_col, label) {
  y  <- df[[sensor_col]]
  x  <- df[[ref_col]]
  ok <- !is.na(y) & !is.na(x)
  if (sum(ok) < 5) {
    message("  Skipping '", label, "' (", sensor_col, " vs ", ref_col,
            ") — fewer than 5 rows")
    return(NULL)
  }
  fit <- lm(y[ok] ~ x[ok])
  res <- y[ok] - x[ok]
  tibble(
    period      = label,
    comparison  = paste0(sensor_col, " vs ", ref_col),
    n_hours     = sum(ok),
    R2          = round(summary(fit)$r.squared, 4),
    slope       = round(unname(coef(fit)[2]),   4),
    intercept   = round(unname(coef(fit)[1]),   4),
    RMSE        = round(sqrt(mean(res^2)),       4),
    MAE         = round(mean(abs(res)),          4),
    Bias        = round(mean(res),               4),
    mean_CRDS   = round(mean(x[ok]), 3),
    mean_sensor = round(mean(y[ok]), 3)
  )
}

pairs <- list(
  list(sensor = "NH3_OTICE",      ref = "NH3_CRDS_in"),
  list(sensor = "NH3_OTICE_corr", ref = "NH3_CRDS_in"),
  list(sensor = "CO2_OTICE",      ref = "CO2_CRDS_in"),
  list(sensor = "CO2_OTICE_corr", ref = "CO2_CRDS_in")
)

stats_overall <- lapply(pairs, function(p) {
  calc_stats(comparison, p$sensor, p$ref, "Overall 2025-09-29 to 2025-12-31")
}) |> bind_rows()

stats_by_device <- lapply(pairs, function(p) {
  comparison |>
    group_by(CRDS_device) |>
    group_modify(~ calc_stats(.x, p$sensor, p$ref, label = .y$CRDS_device)) |>
    ungroup() |>
    select(-CRDS_device)
}) |> bind_rows()

stats_table <- bind_rows(stats_overall, stats_by_device) |>
  arrange(comparison, period)

cat("\n=============================================================================\n")
cat("OTICE node-avg vs CRDS_in — NH3 & CO2 (with/without barn correction)\n")
cat("=============================================================================\n\n")
print(stats_table, n = Inf, width = 140)

# =============================================================================
# 5. COVERAGE SUMMARY BY MONTH
# =============================================================================
coverage <- comparison |>
  mutate(month = format(datetime_hour, "%Y-%m")) |>
  group_by(month) |>
  summarise(
    paired_hours        = n(),
    CRDS_device         = paste(sort(unique(CRDS_device)), collapse = "/"),
    mean_NH3_CRDS_in    = round(mean(NH3_CRDS_in,    na.rm = TRUE), 3),
    mean_NH3_CRDS_s     = round(mean(NH3_CRDS_s,     na.rm = TRUE), 3),
    mean_NH3_CRDS_hires = round(mean(NH3_CRDS_hires, na.rm = TRUE), 3),
    mean_NH3_OTICE      = round(mean(NH3_OTICE,      na.rm = TRUE), 3),
    mean_NH3_OTICE_corr = round(mean(NH3_OTICE_corr, na.rm = TRUE), 3),
    .groups = "drop"
  )

cat("\nMonthly coverage summary:\n")
print(coverage, n = Inf, width = 140)

# =============================================================================
# 6. PLOTS — node-average OTICE vs CRDS_in
# =============================================================================

scatter_panel <- function(df, sensor_col, ref_col, title) {
  st  <- stats_table |>
    filter(comparison == paste0(sensor_col, " vs ", ref_col),
           period == "Overall 2025-09-29 to 2025-12-31")
  lbl <- paste0("n=", st$n_hours, "\nR\u00b2=", st$R2,
                "\nSlope=", st$slope, "\nRMSE=", st$RMSE, "\nBias=", st$Bias)
  d    <- df |> filter(!is.na(.data[[ref_col]]), !is.na(.data[[sensor_col]]))
  lims <- range(c(d[[ref_col]], d[[sensor_col]]), na.rm = TRUE)
  pad  <- diff(lims) * 0.05
  lims <- c(lims[1] - pad, lims[2] + pad)
  ggplot(d, aes(x = .data[[ref_col]], y = .data[[sensor_col]])) +
    geom_point(alpha = 0.4, size = 1, color = "steelblue") +
    geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "black", linewidth = 0.5) +
    annotate("text",
             x = lims[1] + diff(lims) * 0.03,
             y = lims[2] - diff(lims) * 0.03,
             label = lbl, hjust = 0, vjust = 1, size = 3.2) +
    scale_x_continuous(limits = lims) +
    scale_y_continuous(limits = lims) +
    labs(title = title,
         x = paste(ref_col, "(ppm)"), y = paste(sensor_col, "(ppm)")) +
    theme_bw(base_size = 11)
}

p_nh3_scatter <-
  scatter_panel(comparison, "NH3_OTICE",      "NH3_CRDS_in", "NH3: raw OTICE vs CRDS_in") +
  scatter_panel(comparison, "NH3_OTICE_corr", "NH3_CRDS_in", "NH3: barn-corr OTICE vs CRDS_in") +
  plot_annotation(
    title    = "NH3 — OTICE node average vs CRDS_in (hourly, ppm)",
    subtitle = "Red: OLS regression  |  Dashed: 1:1 line",
    theme    = theme(plot.title = element_text(size = 13, face = "bold"))
  )
ggsave(file.path(out_dir, "scatter_NH3_nodeavg_vs_CRDS.png"),
       p_nh3_scatter, width = 12, height = 6, dpi = 150)
message("NH3 scatter saved.")

p_co2_scatter <-
  scatter_panel(comparison, "CO2_OTICE",      "CO2_CRDS_in", "CO2: raw OTICE vs CRDS_in") +
  scatter_panel(comparison, "CO2_OTICE_corr", "CO2_CRDS_in", "CO2: barn-corr OTICE vs CRDS_in") +
  plot_annotation(
    title    = "CO2 — OTICE node average vs CRDS_in (hourly, ppm)",
    subtitle = "Red: OLS regression  |  Dashed: 1:1 line",
    theme    = theme(plot.title = element_text(size = 13, face = "bold"))
  )
ggsave(file.path(out_dir, "scatter_CO2_nodeavg_vs_CRDS.png"),
       p_co2_scatter, width = 12, height = 6, dpi = 150)
message("CO2 scatter saved.")

ts_long <- comparison |>
  select(datetime_hour, NH3_CRDS_in, NH3_OTICE, NH3_OTICE_corr) |>
  pivot_longer(-datetime_hour, names_to = "source", values_to = "NH3_ppm") |>
  mutate(source = factor(source,
    levels = c("NH3_CRDS_in", "NH3_OTICE", "NH3_OTICE_corr"),
    labels = c("CRDS (in)", "OTICE raw (node avg)", "OTICE barn-corr (node avg)")
  ))
p_ts <- ggplot(ts_long |> filter(!is.na(NH3_ppm)),
               aes(x = datetime_hour, y = NH3_ppm, color = source)) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = c("CRDS (in)"                  = "black",
                                "OTICE raw (node avg)"       = "steelblue",
                                "OTICE barn-corr (node avg)" = "tomato")) +
  labs(title = "NH3 Time Series — CRDS_in vs OTICE node average",
       x = NULL, y = "NH3 (ppm)", color = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "top")
ggsave(file.path(out_dir, "timeseries_NH3_nodeavg.png"),
       p_ts, width = 14, height = 5, dpi = 150)
message("NH3 time series saved.")

ba_panel <- function(df, sensor_col, ref_col, title) {
  d    <- df |> filter(!is.na(.data[[sensor_col]]), !is.na(.data[[ref_col]]))
  diff <- d[[sensor_col]] - d[[ref_col]]
  avg  <- (d[[sensor_col]] + d[[ref_col]]) / 2
  m    <- mean(diff);  s <- sd(diff)
  ggplot(data.frame(avg = avg, diff = diff), aes(x = avg, y = diff)) +
    geom_point(alpha = 0.35, size = 0.9, color = "steelblue") +
    geom_hline(yintercept = m,             color = "red",    linewidth = 0.7) +
    geom_hline(yintercept = m + 1.96 * s,  color = "orange", linetype = "dashed") +
    geom_hline(yintercept = m - 1.96 * s,  color = "orange", linetype = "dashed") +
    geom_hline(yintercept = 0,             color = "black",  linetype = "dotted") +
    annotate("text", x = Inf, y = m,
             label = paste0("  bias = ", round(m, 3)),
             hjust = 1, vjust = -0.5, size = 3, color = "red") +
    annotate("text", x = Inf, y = m + 1.96 * s,
             label = paste0("  +1.96 SD = ", round(m + 1.96 * s, 3)),
             hjust = 1, vjust = -0.5, size = 3, color = "darkorange") +
    annotate("text", x = Inf, y = m - 1.96 * s,
             label = paste0("  \u22121.96 SD = ", round(m - 1.96 * s, 3)),
             hjust = 1, vjust = 1.5,  size = 3, color = "darkorange") +
    labs(title = title,
         x = paste0("Mean of ", sensor_col, " and ", ref_col, " (ppm)"),
         y = paste0(sensor_col, " \u2212 ", ref_col, " (ppm)")) +
    theme_bw(base_size = 11)
}
p_ba <-
  ba_panel(comparison, "NH3_OTICE",      "NH3_CRDS_in", "Bland-Altman: NH3 raw") +
  ba_panel(comparison, "NH3_OTICE_corr", "NH3_CRDS_in", "Bland-Altman: NH3 barn-corr") +
  plot_annotation(
    title    = "Bland-Altman — OTICE node average vs CRDS_in (NH3)",
    subtitle = "Red: mean bias  |  Orange dashed: \u00b11.96 SD limits of agreement",
    theme    = theme(plot.title = element_text(size = 13, face = "bold"))
  )
ggsave(file.path(out_dir, "bland_altman_NH3_nodeavg.png"),
       p_ba, width = 14, height = 6, dpi = 150)
message("Bland-Altman saved.")

# =============================================================================
# 7. TIME SERIES — CRDS_hires vs OTICE_corr (daily x-axis scale)
# =============================================================================

# Shared x-axis: one major label per week, minor gridline every 24 hours (1 day)
x_scale_daily <- scale_x_datetime(
  date_breaks       = "1 week",
  date_minor_breaks = "1 day",
  date_labels       = "%d %b",
  expand            = expansion(add = c(0, 0))
)

ts_theme <- theme_bw(base_size = 11) +
  theme(
    legend.position   = "top",
    legend.key.width  = unit(1.5, "cm"),
    panel.grid.major  = element_line(colour = "grey80", linewidth = 0.4),
    panel.grid.minor  = element_line(colour = "grey92", linewidth = 0.25),
    axis.text.x       = element_text(angle = 45, hjust = 1),
    plot.title        = element_text(face = "bold", size = 12)
  )

# --- NH3 panel ---------------------------------------------------------------
nh3_hires <- comparison |>
  select(datetime_hour, NH3_CRDS_hires, NH3_OTICE_corr) |>
  pivot_longer(-datetime_hour, names_to = "source", values_to = "ppm") |>
  mutate(source = factor(source,
    levels = c("NH3_CRDS_hires", "NH3_OTICE_corr"),
    labels = c("CRDS hires",     "OTICE barn-corr (node avg)")
  ))

p_nh3_hires <- ggplot(nh3_hires |> filter(!is.na(ppm)),
                      aes(x = datetime_hour, y = ppm,
                          color = source, linewidth = source)) +
  geom_line() +
  x_scale_daily +
  scale_color_manual(
    values = c("CRDS hires"               = "#1a1a2e",
               "OTICE barn-corr (node avg)"= "#e63946")
  ) +
  scale_linewidth_manual(
    values = c("CRDS hires"               = 0.7,
               "OTICE barn-corr (node avg)"= 0.55),
    guide  = "none"
  ) +
  labs(
    title = "NH3 — CRDS hires vs OTICE barn-corrected (node average)",
    x     = NULL,
    y     = "NH\u2083 (ppm)",
    color = NULL
  ) +
  ts_theme

# --- CO2 panel ---------------------------------------------------------------
co2_hires <- comparison |>
  select(datetime_hour, CO2_CRDS_hires, CO2_OTICE_corr) |>
  pivot_longer(-datetime_hour, names_to = "source", values_to = "ppm") |>
  mutate(source = factor(source,
    levels = c("CO2_CRDS_hires", "CO2_OTICE_corr"),
    labels = c("CRDS hires",     "OTICE barn-corr (node avg)")
  ))

p_co2_hires <- ggplot(co2_hires |> filter(!is.na(ppm)),
                      aes(x = datetime_hour, y = ppm,
                          color = source, linewidth = source)) +
  geom_line() +
  x_scale_daily +
  scale_color_manual(
    values = c("CRDS hires"               = "#1a1a2e",
               "OTICE barn-corr (node avg)"= "#457b9d")
  ) +
  scale_linewidth_manual(
    values = c("CRDS hires"               = 0.7,
               "OTICE barn-corr (node avg)"= 0.55),
    guide  = "none"
  ) +
  labs(
    title = "CO\u2082 — CRDS hires vs OTICE barn-corrected (node average)",
    x     = NULL,
    y     = "CO\u2082 (ppm)",
    color = NULL
  ) +
  ts_theme

# --- Combined (stacked) -------------------------------------------------------
p_hires_ts <- p_nh3_hires / p_co2_hires +
  plot_annotation(
    caption = "X-axis: major labels every week  |  minor gridlines every 24 h",
    theme   = theme(plot.caption = element_text(size = 8, colour = "grey50"))
  )

ggsave(file.path(out_dir, "timeseries_hires_vs_OTICEcorr.png"),
       p_hires_ts, width = 20, height = 10, dpi = 150)
message("CRDS hires vs OTICE_corr time series saved.")

# =============================================================================
# 8. DAILY AVERAGES — CRDS_hires vs OTICE_corr
#    Aggregate hourly paired data to 24-hour means; days with < 12 paired hours
#    are excluded to avoid misleading daily averages from sparse data.
# =============================================================================
comparison_daily <- comparison |>
  mutate(date = as.Date(datetime_hour, tz = "Europe/Berlin")) |>
  group_by(date) |>
  summarise(
    NH3_CRDS_hires = safe_mean(NH3_CRDS_hires),
    NH3_OTICE_corr = safe_mean(NH3_OTICE_corr),
    CO2_CRDS_hires = safe_mean(CO2_CRDS_hires),
    CO2_OTICE_corr = safe_mean(CO2_OTICE_corr),
    n_hours        = n(),
    .groups        = "drop"
  ) |>
  filter(n_hours >= 12)

message("\nDaily paired rows (>= 12 hrs): ", nrow(comparison_daily))

# --- Daily statistics ---------------------------------------------------------
calc_stats_daily <- function(sensor_col, ref_col) {
  y  <- comparison_daily[[sensor_col]]
  x  <- comparison_daily[[ref_col]]
  ok <- !is.na(y) & !is.na(x)
  if (sum(ok) < 5) return(NULL)
  fit <- lm(y[ok] ~ x[ok])
  res <- y[ok] - x[ok]
  tibble(
    comparison  = paste0(sensor_col, " vs ", ref_col),
    n_days      = sum(ok),
    R2          = round(summary(fit)$r.squared, 4),
    slope       = round(unname(coef(fit)[2]),   4),
    intercept   = round(unname(coef(fit)[1]),   4),
    RMSE        = round(sqrt(mean(res^2)),       4),
    MAE         = round(mean(abs(res)),          4),
    Bias        = round(mean(res),               4)
  )
}

daily_stats <- bind_rows(
  calc_stats_daily("NH3_OTICE_corr", "NH3_CRDS_hires"),
  calc_stats_daily("CO2_OTICE_corr", "CO2_CRDS_hires")
)

cat("\n===== Daily statistics (CRDS_hires vs OTICE_corr) =====\n")
print(daily_stats, width = 120)

write.csv(comparison_daily, file.path(out_dir, "comparison_daily.csv"), row.names = FALSE)
write.csv(daily_stats,      file.path(out_dir, "daily_stats.csv"),      row.names = FALSE)
message("Daily CSVs saved.")

# --- Daily scatter plots (side-by-side) ---------------------------------------
daily_scatter <- function(sensor_col, ref_col, gas_label, st) {
  d    <- comparison_daily |> filter(!is.na(.data[[ref_col]]), !is.na(.data[[sensor_col]]))
  lims <- range(c(d[[ref_col]], d[[sensor_col]]), na.rm = TRUE)
  pad  <- diff(lims) * 0.05
  lims <- c(lims[1] - pad, lims[2] + pad)
  lbl  <- paste0("n=", st$n_days, "\nR\u00b2=", st$R2,
                 "\nSlope=", st$slope, "\nRMSE=", st$RMSE, "\nBias=", st$Bias)
  ggplot(d, aes(x = .data[[ref_col]], y = .data[[sensor_col]])) +
    geom_point(alpha = 0.7, size = 2, color = "steelblue") +
    geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "black", linewidth = 0.5) +
    annotate("text",
             x = lims[1] + diff(lims) * 0.03,
             y = lims[2] - diff(lims) * 0.03,
             label = lbl, hjust = 0, vjust = 1, size = 3.2) +
    scale_x_continuous(limits = lims) +
    scale_y_continuous(limits = lims) +
    labs(
      title = paste0(gas_label, ": daily mean OTICE_corr vs CRDS_hires"),
      x     = paste0(ref_col,    " (ppm)"),
      y     = paste0(sensor_col, " (ppm)")
    ) +
    theme_bw(base_size = 11)
}

p_daily_scatter <-
  daily_scatter("NH3_OTICE_corr", "NH3_CRDS_hires", "NH3",
                daily_stats |> filter(comparison == "NH3_OTICE_corr vs NH3_CRDS_hires")) +
  daily_scatter("CO2_OTICE_corr", "CO2_CRDS_hires", "CO2",
                daily_stats |> filter(comparison == "CO2_OTICE_corr vs CO2_CRDS_hires")) +
  plot_annotation(
    title    = "Daily mean — OTICE barn-corr (node avg) vs CRDS hires",
    subtitle = "Each point = one calendar day (days with \u2265 12 paired hours)  |  Red: OLS  |  Dashed: 1:1",
    theme    = theme(plot.title = element_text(size = 13, face = "bold"))
  )

ggsave(file.path(out_dir, "scatter_daily_CRDS_hires_vs_OTICEcorr.png"),
       p_daily_scatter, width = 12, height = 6, dpi = 150)
message("Daily scatter saved.")

# --- Daily time series (stacked) ----------------------------------------------
x_scale_weekly <- scale_x_date(
  date_breaks       = "2 weeks",
  date_minor_breaks = "1 week",
  date_labels       = "%d %b",
  expand            = expansion(add = c(1, 1))
)

daily_ts_theme <- theme_bw(base_size = 11) +
  theme(
    legend.position  = "top",
    legend.key.width = unit(1.2, "cm"),
    panel.grid.major = element_line(colour = "grey80", linewidth = 0.4),
    panel.grid.minor = element_line(colour = "grey92", linewidth = 0.25),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    plot.title       = element_text(face = "bold", size = 12)
  )

p_daily_nh3 <- ggplot(comparison_daily |> filter(!is.na(NH3_CRDS_hires) | !is.na(NH3_OTICE_corr))) +
  geom_line(aes(x = date, y = NH3_CRDS_hires,  color = "CRDS hires"),        linewidth = 0.8) +
  geom_point(aes(x = date, y = NH3_CRDS_hires, color = "CRDS hires"),        size = 1.8) +
  geom_line(aes(x = date, y = NH3_OTICE_corr,  color = "OTICE barn-corr"),   linewidth = 0.7) +
  geom_point(aes(x = date, y = NH3_OTICE_corr, color = "OTICE barn-corr"),   size = 1.8) +
  scale_color_manual(values = c("CRDS hires" = "#1a1a2e", "OTICE barn-corr" = "#e63946")) +
  x_scale_weekly +
  labs(
    title = "NH\u2083 — Daily mean: CRDS hires vs OTICE barn-corrected (node avg)",
    x     = NULL, y = "NH\u2083 (ppm)", color = NULL
  ) +
  daily_ts_theme

p_daily_co2 <- ggplot(comparison_daily |> filter(!is.na(CO2_CRDS_hires) | !is.na(CO2_OTICE_corr))) +
  geom_line(aes(x = date, y = CO2_CRDS_hires,  color = "CRDS hires"),       linewidth = 0.8) +
  geom_point(aes(x = date, y = CO2_CRDS_hires, color = "CRDS hires"),       size = 1.8) +
  geom_line(aes(x = date, y = CO2_OTICE_corr,  color = "OTICE barn-corr"),  linewidth = 0.7) +
  geom_point(aes(x = date, y = CO2_OTICE_corr, color = "OTICE barn-corr"),  size = 1.8) +
  scale_color_manual(values = c("CRDS hires" = "#1a1a2e", "OTICE barn-corr" = "#457b9d")) +
  x_scale_weekly +
  labs(
    title = "CO\u2082 — Daily mean: CRDS hires vs OTICE barn-corrected (node avg)",
    x     = NULL, y = "CO\u2082 (ppm)", color = NULL
  ) +
  daily_ts_theme

p_daily_ts <- p_daily_nh3 / p_daily_co2 +
  plot_annotation(
    caption = "X-axis: major labels every 2 weeks  |  minor gridlines every week",
    theme   = theme(plot.caption = element_text(size = 8, colour = "grey50"))
  )

ggsave(file.path(out_dir, "timeseries_daily_CRDS_hires_vs_OTICEcorr.png"),
       p_daily_ts, width = 16, height = 10, dpi = 150)
message("Daily time series saved.")

# =============================================================================
# 9. SAVE
#    comparison_hourly.csv: only datetime + gas columns (no metadata)
# =============================================================================
comparison_out <- comparison |>
  select(datetime_hour,
         NH3_CRDS_in, NH3_CRDS_s, NH3_CRDS_hires,
         CO2_CRDS_in, CO2_CRDS_s, CO2_CRDS_hires,
         NH3_OTICE, NH3_OTICE_corr,
         CO2_OTICE, CO2_OTICE_corr)

write.csv(comparison_out,    file.path(out_dir, "comparison_hourly.csv"),  row.names = FALSE)
write.csv(stats_table,       file.path(out_dir, "comparison_stats.csv"),   row.names = FALSE)
write.csv(otice_node_hourly, file.path(out_dir, "OTICE_node_hourly.csv"),  row.names = FALSE)

cat("\nFiles saved to:", out_dir, "\n")
cat("  comparison_hourly.csv                      <- 10 columns: datetime + 3 CRDS groups + 4 OTICE\n")
cat("  comparison_stats.csv                       <- R2, RMSE, bias for all pairs vs CRDS_in\n")
cat("  OTICE_node_hourly.csv                      <- per-node hourly (nodes O04/O09/O13 excluded)\n")
cat("  comparison_daily.csv                       <- daily means (days >= 12 paired hours)\n")
cat("  daily_stats.csv                            <- R2, RMSE, bias at daily resolution\n")
cat("  timeseries_hires_vs_OTICEcorr.png         <- NH3 & CO2 CRDS_hires vs OTICE_corr (hourly)\n")
cat("  scatter_daily_CRDS_hires_vs_OTICEcorr.png <- daily scatter: NH3 & CO2\n")
cat("  timeseries_daily_CRDS_hires_vs_OTICEcorr.png <- daily time series: NH3 & CO2\n")
cat("\nColumn naming:\n")
cat("  NH3/CO2_CRDS_in     = CRDS location 'in'\n")
cat("  NH3/CO2_CRDS_s      = CRDS location 's' (south)\n")
cat("  NH3/CO2_CRDS_hires  = CRDS all other locations (high-res)\n")
cat("  NH3/CO2_OTICE       = OTICE raw calibrated\n")
cat("  NH3/CO2_OTICE_corr  = OTICE barn-corrected\n")
