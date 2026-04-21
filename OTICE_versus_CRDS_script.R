# =============================================================================
# OTICE versus CRDS
# Minute-level campaign comparison of OTICE node averages versus CRDS reference
# averages, based on the campaign metadata document.
# =============================================================================

library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(patchwork)

# -----------------------------------------------------------------------------
# 1. Paths and basic settings
# -----------------------------------------------------------------------------
timezone_local <- "Europe/Berlin"

otice_dir <- "processed_data/OTICE_processed"
crds_dir  <- "processed_data/CRDS8_processed"
out_dir   <- file.path("output", "OTICE_versus_CRDS")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "campaign_data"), showWarnings = FALSE)
dir.create(file.path(out_dir, "campaign_plots"), showWarnings = FALSE)

unlink(list.files(file.path(out_dir, "campaign_data"), full.names = TRUE), force = TRUE)
unlink(list.files(file.path(out_dir, "campaign_plots"), full.names = TRUE), force = TRUE)

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
                                     analyzer, node_locations) {
  bind_rows(lapply(names(node_locations), function(node_id) {
    tibble(
      period_id      = period_id,
      period_label   = period_label,
      start_time     = as.POSIXct(start_time, tz = timezone_local),
      end_time       = as.POSIXct(end_time,   tz = timezone_local),
      OTICE_location = node_id,
      analyzer       = analyzer,
      CRDS_location  = as.character(node_locations[[node_id]])
    )
  }))
}

calc_campaign_stats <- function(df, sensor_col, ref_col, gas_label) {
  ok <- !is.na(df[[sensor_col]]) & !is.na(df[[ref_col]])
  if (sum(ok) < 5) {
    return(NULL)
  }

  fit <- lm(df[[sensor_col]][ok] ~ df[[ref_col]][ok])
  residuals <- df[[sensor_col]][ok] - df[[ref_col]][ok]

  tibble(
    period_id    = df$period_id[1],
    period_label = df$period_label[1],
    gas          = gas_label,
    n_overlap    = sum(ok),
    R2           = round(summary(fit)$r.squared, 4),
    slope        = round(unname(coef(fit)[2]), 4),
    intercept    = round(unname(coef(fit)[1]), 4),
    RMSE         = round(sqrt(mean(residuals^2)), 4),
    MAE          = round(mean(abs(residuals)), 4),
    Bias         = round(mean(residuals), 4)
  )
}

make_campaign_plot <- function(df, stats_tbl) {
  nh3_stats <- stats_tbl |>
    filter(gas == "NH3")

  co2_stats <- stats_tbl |>
    filter(gas == "CO2")

  nh3_label <- if (nrow(nh3_stats) == 1) {
    paste0("n=", nh3_stats$n_overlap,
           "\nR2=", nh3_stats$R2,
           "\nSlope=", nh3_stats$slope,
           "\nBias=", nh3_stats$Bias)
  } else {
    "Insufficient overlap"
  }

  co2_label <- if (nrow(co2_stats) == 1) {
    paste0("n=", co2_stats$n_overlap,
           "\nR2=", co2_stats$R2,
           "\nSlope=", co2_stats$slope,
           "\nBias=", co2_stats$Bias)
  } else {
    "Insufficient overlap"
  }

  p_nh3 <- ggplot(df, aes(x = datetime_hour)) +
    geom_line(aes(y = CRDS_NH3, color = "CRDS average"), linewidth = 0.7, na.rm = TRUE) +
    geom_line(aes(y = OTICE_NH3, color = "OTICE average"), linewidth = 0.6, na.rm = TRUE) +
    geom_point(aes(y = CRDS_NH3, color = "CRDS average"), size = 0.6, alpha = 0.5, na.rm = TRUE) +
    geom_point(aes(y = OTICE_NH3, color = "OTICE average"), size = 0.6, alpha = 0.5, na.rm = TRUE) +
    annotate("text",
             x = min(df$datetime_hour, na.rm = TRUE),
             y = max(c(df$CRDS_NH3, df$OTICE_NH3), na.rm = TRUE),
             label = nh3_label,
             hjust = 0, vjust = 1, size = 3.2) +
    scale_color_manual(values = c("CRDS average" = "#1a1a2e",
                                  "OTICE average" = "#e76f51")) +
    scale_x_datetime(date_breaks = "1 day",
                     date_labels = "%d %b",
                     expand = expansion(mult = c(0.01, 0.02))) +
    labs(title = "NH3: hourly average OTICE nodes versus hourly average CRDS references",
         x = NULL,
         y = "NH3 (ppm)",
         color = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold"))

  p_co2 <- ggplot(df, aes(x = datetime_hour)) +
    geom_line(aes(y = CRDS_CO2, color = "CRDS average"), linewidth = 0.7, na.rm = TRUE) +
    geom_line(aes(y = OTICE_CO2, color = "OTICE average"), linewidth = 0.6, na.rm = TRUE) +
    geom_point(aes(y = CRDS_CO2, color = "CRDS average"), size = 0.6, alpha = 0.5, na.rm = TRUE) +
    geom_point(aes(y = OTICE_CO2, color = "OTICE average"), size = 0.6, alpha = 0.5, na.rm = TRUE) +
    annotate("text",
             x = min(df$datetime_hour, na.rm = TRUE),
             y = max(c(df$CRDS_CO2, df$OTICE_CO2), na.rm = TRUE),
             label = co2_label,
             hjust = 0, vjust = 1, size = 3.2) +
    scale_color_manual(values = c("CRDS average" = "#1a1a2e",
                                  "OTICE average" = "#2a9d8f")) +
    scale_x_datetime(date_breaks = "1 day",
                     date_labels = "%d %b",
                     expand = expansion(mult = c(0.01, 0.02))) +
    labs(title = "CO2: hourly average OTICE nodes versus hourly average CRDS references",
         x = NULL,
         y = "CO2 (ppm)",
         color = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold"))

  p_nh3 / p_co2
}

# -----------------------------------------------------------------------------
# 3. Campaign metadata from OTICE_CRDS_calibration_routine_Silvia.docx
# -----------------------------------------------------------------------------
campaign_reference_map <- bind_rows(
  build_campaign_reference(
    period_id = "C1",
    period_label = "2025-09-29 to 2025-10-06",
    start_time = "2025-09-29 00:00:00",
    end_time   = "2025-10-06 00:00:00",
    analyzer = "CRDS9",
    node_locations = list(
      "2" = c("24", "27"), "6" = c("24", "27"), "7" = c("24", "27"),
      "11" = c("24", "27"), "17" = c("24", "27"), "18" = c("24", "27"),
      "14" = c("24", "27")
    )
  ),
  build_campaign_reference(
    period_id = "C2",
    period_label = "2025-10-06 to 2025-10-21",
    start_time = "2025-10-06 00:00:00",
    end_time   = "2025-10-21 00:00:00",
    analyzer = "CRDS8",
    node_locations = list(
      "2" = "30", "6" = "30", "7" = "30", "11" = "30",
      "17" = "30", "18" = "30", "14" = "30", "5" = "30"
    )
  ),
  build_campaign_reference(
    period_id = "C3",
    period_label = "2025-10-21 to 2025-11-10",
    start_time = "2025-10-21 00:00:00",
    end_time   = "2025-11-10 00:00:00",
    analyzer = "CRDS8",
    node_locations = list(
      "2" = "30", "6" = "30", "7" = "30", "11" = "30",
      "17" = "30", "18" = "30", "14" = "30", "5" = "30"
    )
  ),
  build_campaign_reference(
    period_id = "C4",
    period_label = "2025-11-10 to 2025-11-12",
    start_time = "2025-11-10 00:00:00",
    end_time   = "2025-11-12 00:00:00",
    analyzer = "CRDS8",
    node_locations = list(
      "2" = "30", "6" = "30", "7" = "30", "11" = "30",
      "17" = "30", "18" = "30", "14" = "30", "3" = "30"
    )
  ),
  build_campaign_reference(
    period_id = "C5",
    period_label = "2025-11-12 to 2025-11-19",
    start_time = "2025-11-12 00:00:00",
    end_time   = "2025-11-19 00:00:00",
    analyzer = "CRDS8",
    node_locations = list(
      "2" = "30", "6" = "30", "7" = "30", "11" = "30",
      "17" = "30", "18" = "30", "14" = "30", "3" = "30"
    )
  ),
  build_campaign_reference(
    period_id = "C6",
    period_label = "2025-11-26 to 2025-12-01",
    start_time = "2025-11-26 00:00:00",
    end_time   = "2025-12-01 00:00:00",
    analyzer = "CRDS8",
    node_locations = list(
      "2" = "36", "6" = "30", "7" = "15", "11" = "45",
      "17" = "42", "18" = "12", "3" = "18"
    )
  ),
  build_campaign_reference(
    period_id = "C7",
    period_label = "2025-12-01 to 2025-12-09",
    start_time = "2025-12-01 00:00:00",
    end_time   = "2025-12-09 00:00:00",
    analyzer = "CRDS8",
    node_locations = list(
      "2" = "36", "6" = "30", "7" = "15", "11" = "45",
      "17" = "42", "18" = "12", "3" = "18"
    )
  ),
  build_campaign_reference(
    period_id = "C8",
    period_label = "2025-12-09 to 2025-12-31",
    start_time = "2025-12-09 00:00:00",
    end_time   = "2026-01-01 00:00:00",
    analyzer = "CRDS8",
    node_locations = list(
      "2" = "5", "6" = "4", "7" = "3", "11" = "7",
      "17" = "6", "18" = "1", "3" = "2"
    )
  )
)

campaign_periods <- campaign_reference_map |>
  distinct(period_id, period_label, start_time, end_time) |>
  arrange(start_time)

write.csv(campaign_reference_map,
          file.path(out_dir, "campaign_reference_map.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 4. Read and standardize OTICE minute data
# -----------------------------------------------------------------------------
otice_files <- list.files(otice_dir, pattern = "^min_calibrated.*\\.csv$", full.names = TRUE)
otice_files <- otice_files[file.info(otice_files)$size > 0]

OTICE_dataset <- lapply(otice_files, read.csv, stringsAsFactors = FALSE) |>
  bind_rows() |>
  transmute(
    DATE.TIME = as.POSIXct(Datetime_Berlin, format = "%Y-%m-%d %H:%M:%S", tz = timezone_local),
    location  = clean_otice_location(Node),
    analyzer  = as.character(Type),
    CO2       = suppressWarnings(as.numeric(CO2.RAW)),
    NH3       = suppressWarnings(as.numeric(NH3_ppm)),
    NH3_cal   = suppressWarnings(as.numeric(NH3_ppm_barn)),
    CO2_cal   = suppressWarnings(as.numeric(CO2.AVG_barn))
  ) |>
  filter(!is.na(DATE.TIME)) |>
  distinct(analyzer, location, DATE.TIME, .keep_all = TRUE) |>
  arrange(DATE.TIME, location)

message("OTICE rows in standardised dataset: ", nrow(OTICE_dataset))
message("OTICE locations: ", paste(sort(unique(OTICE_dataset$location)), collapse = ", "))

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
    location  = clean_crds_location(location),
    analyzer  = toupper(trimws(as.character(analyzer))),
    CO2       = suppressWarnings(as.numeric(CO2)),
    NH3       = suppressWarnings(as.numeric(NH3)),
    NH3_cal   = NA_real_,
    CO2_cal   = NA_real_
  ) |>
  filter(!is.na(DATE.TIME)) |>
  distinct(analyzer, location, DATE.TIME, .keep_all = TRUE) |>
  arrange(DATE.TIME, analyzer, location)

message("CRDS rows in standardised dataset: ", nrow(CRDS_dataset))
message("CRDS analyzers: ", paste(sort(unique(CRDS_dataset$analyzer)), collapse = ", "))

# Combined standardised dataset only to keep both sources in one comparable table.
combined_standardised_dataset <- bind_rows(
  OTICE_dataset |> mutate(dataset = "OTICE"),
  CRDS_dataset  |> mutate(dataset = "CRDS")
)

write.csv(combined_standardised_dataset,
          file.path(out_dir, "combined_standardised_dataset.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 6. Keep only campaign periods where OTICE and CRDS are defined in parallel
# -----------------------------------------------------------------------------
otice_campaign_data <- OTICE_dataset |>
  mutate(datetime_hour = floor_date(DATE.TIME, "hour")) |>
  inner_join(
    campaign_reference_map |>
      distinct(period_id, period_label, start_time, end_time, OTICE_location),
    by = c("location" = "OTICE_location"),
    relationship = "many-to-many"
  ) |>
  filter(DATE.TIME >= start_time, DATE.TIME < end_time)

crds_campaign_data <- CRDS_dataset |>
  mutate(datetime_hour = floor_date(DATE.TIME, "hour")) |>
  inner_join(
    campaign_reference_map |>
      distinct(period_id, period_label, start_time, end_time, analyzer, CRDS_location),
    by = c("analyzer", "location" = "CRDS_location"),
    relationship = "many-to-many"
  ) |>
  filter(DATE.TIME >= start_time, DATE.TIME < end_time)

message("OTICE rows inside campaign periods: ", nrow(otice_campaign_data))
message("CRDS rows inside campaign periods: ", nrow(crds_campaign_data))

# -----------------------------------------------------------------------------
# 7. Average available OTICE and CRDS locations for each campaign and hour
# -----------------------------------------------------------------------------
otice_campaign_average <- otice_campaign_data |>
  group_by(period_id, period_label, start_time, end_time, datetime_hour) |>
  summarise(
    OTICE_NH3        = safe_mean(NH3),
    OTICE_CO2        = safe_mean(CO2),
    OTICE_NH3_cal    = safe_mean(NH3_cal),
    OTICE_CO2_cal    = safe_mean(CO2_cal),
    n_OTICE_locations = n_distinct(location[!is.na(NH3) | !is.na(CO2)]),
    OTICE_locations   = paste(sort(unique(location[!is.na(NH3) | !is.na(CO2)])), collapse = ", "),
    .groups = "drop"
  )

crds_campaign_average <- crds_campaign_data |>
  group_by(period_id, period_label, start_time, end_time, datetime_hour) |>
  summarise(
    CRDS_NH3        = safe_mean(NH3),
    CRDS_CO2        = safe_mean(CO2),
    n_CRDS_locations = n_distinct(location[!is.na(NH3) | !is.na(CO2)]),
    CRDS_locations   = paste(sort(unique(location[!is.na(NH3) | !is.na(CO2)])), collapse = ", "),
    .groups = "drop"
  )

campaign_average_comparison <- full_join(
  otice_campaign_average,
  crds_campaign_average,
  by = c("period_id", "period_label", "start_time", "end_time", "datetime_hour")
) |>
  filter(
    !is.na(OTICE_NH3),
    !is.na(CRDS_NH3),
    !is.na(OTICE_CO2),
    !is.na(CRDS_CO2)
  ) |>
  arrange(start_time, datetime_hour)

write.csv(campaign_average_comparison,
          file.path(out_dir, "campaign_average_comparison.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 8. Keep only campaigns with real OTICE-versus-CRDS overlap
# -----------------------------------------------------------------------------
campaign_summary_all <- campaign_average_comparison |>
  group_by(period_id, period_label) |>
  summarise(
    total_hours = n(),
    overlap_hours = n(),
    .groups = "drop"
  )

valid_campaigns <- campaign_summary_all |>
  filter(overlap_hours > 0)

campaign_average_comparison_valid <- campaign_average_comparison |>
  semi_join(valid_campaigns, by = c("period_id", "period_label"))

write.csv(campaign_summary_all,
          file.path(out_dir, "campaign_summary_all_periods.csv"),
          row.names = FALSE)

write.csv(campaign_average_comparison_valid,
          file.path(out_dir, "campaign_average_comparison_valid_only.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 9. Campaign statistics and campaign-specific datasets
# -----------------------------------------------------------------------------
campaign_split <- split(campaign_average_comparison_valid, campaign_average_comparison_valid$period_id)

campaign_stats <- bind_rows(lapply(campaign_split, function(df) {
  bind_rows(
    calc_campaign_stats(df, "OTICE_NH3", "CRDS_NH3", "NH3"),
    calc_campaign_stats(df, "OTICE_CO2", "CRDS_CO2", "CO2")
  )
}))

write.csv(campaign_stats,
          file.path(out_dir, "campaign_stats.csv"),
          row.names = FALSE)

for (campaign_id in names(campaign_split)) {
  campaign_df <- campaign_split[[campaign_id]]
  campaign_label <- unique(campaign_df$period_label)

  write.csv(
    campaign_df,
    file.path(out_dir, "campaign_data",
              paste0(campaign_id, "_",
                     gsub("[^0-9A-Za-z]+", "_", campaign_label),
                     "_average_comparison.csv")),
    row.names = FALSE
  )
}

# -----------------------------------------------------------------------------
# 10. Plots for each campaign
# -----------------------------------------------------------------------------
for (campaign_id in names(campaign_split)) {
  campaign_df <- campaign_split[[campaign_id]]
  campaign_label <- unique(campaign_df$period_label)
  campaign_stats_this <- campaign_stats |>
    filter(period_id == campaign_id)

  if (nrow(campaign_df) == 0) {
    next
  }

  campaign_plot <- make_campaign_plot(campaign_df, campaign_stats_this) +
    plot_annotation(
      title = paste0(campaign_id, ": ", campaign_label),
      subtitle = "Hourly comparison of average available OTICE nodes and average available CRDS references",
      theme = theme(plot.title = element_text(size = 13, face = "bold"))
    )

  ggsave(
    filename = file.path(out_dir, "campaign_plots",
                         paste0(campaign_id, "_",
                                gsub("[^0-9A-Za-z]+", "_", campaign_label),
                                ".png")),
    plot = campaign_plot,
    width = 16,
    height = 10,
    dpi = 150
  )
}

# -----------------------------------------------------------------------------
# 11. Combined overview plots
# -----------------------------------------------------------------------------
overview_nh3 <- campaign_average_comparison_valid |>
  select(period_label, datetime_hour, CRDS_NH3, OTICE_NH3) |>
  pivot_longer(cols = c(CRDS_NH3, OTICE_NH3),
               names_to = "source",
               values_to = "NH3") |>
  mutate(source = recode(source,
                         CRDS_NH3 = "CRDS average",
                         OTICE_NH3 = "OTICE average"))

overview_co2 <- campaign_average_comparison_valid |>
  select(period_label, datetime_hour, CRDS_CO2, OTICE_CO2) |>
  pivot_longer(cols = c(CRDS_CO2, OTICE_CO2),
               names_to = "source",
               values_to = "CO2") |>
  mutate(source = recode(source,
                         CRDS_CO2 = "CRDS average",
                         OTICE_CO2 = "OTICE average"))

p_overview_nh3 <- ggplot(overview_nh3, aes(x = datetime_hour, y = NH3, color = source)) +
  geom_line(linewidth = 0.6, na.rm = TRUE) +
  facet_wrap(~ period_label, scales = "free_x", ncol = 2) +
  scale_color_manual(values = c("CRDS average" = "#1a1a2e",
                                "OTICE average" = "#e76f51")) +
  labs(title = "NH3 overview by campaign",
       x = NULL,
       y = "NH3 (ppm)",
       color = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold"))

p_overview_co2 <- ggplot(overview_co2, aes(x = datetime_hour, y = CO2, color = source)) +
  geom_line(linewidth = 0.6, na.rm = TRUE) +
  facet_wrap(~ period_label, scales = "free_x", ncol = 2) +
  scale_color_manual(values = c("CRDS average" = "#1a1a2e",
                                "OTICE average" = "#2a9d8f")) +
  labs(title = "CO2 overview by campaign",
       x = NULL,
       y = "CO2 (ppm)",
       color = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold"))

ggsave(
  filename = file.path(out_dir, "all_campaigns_average_comparison.png"),
  plot = p_overview_nh3 / p_overview_co2,
  width = 18,
  height = 12,
  dpi = 150
)

# -----------------------------------------------------------------------------
# 12. Console summary
# -----------------------------------------------------------------------------
cat("\n============================================================\n")
cat("OTICE versus CRDS campaign summary\n")
cat("============================================================\n")
print(campaign_summary_all, n = Inf)

cat("\nFiles written to:", out_dir, "\n")
cat("  campaign_reference_map.csv\n")
cat("  combined_standardised_dataset.csv\n")
cat("  campaign_average_comparison.csv\n")
cat("  campaign_average_comparison_valid_only.csv\n")
cat("  campaign_summary_all_periods.csv\n")
cat("  campaign_stats.csv\n")
cat("  campaign_data/*.csv\n")
cat("  campaign_plots/*.png\n")
cat("  all_campaigns_average_comparison.png\n")
