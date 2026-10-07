suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(lubridate)
  library(readr)
  library(tidyr)
})

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args_all[grepl("--file=", args_all)])
script_dir <- dirname(normalizePath(script_path[1], winslash = "/", mustWork = TRUE))
source(file.path(script_dir, "device_comparison_helpers.R"))

workflow_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
input_path <- file.path(workflow_dir, "data", "integrated", "device_comparison_2026081909_2026100610.csv")
output_dir <- file.path(workflow_dir, "result_data", "plots", "cubic_pronova_crds_continuation_20260819")
table_dir <- file.path(workflow_dir, "result_data", "tables")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

time_zone <- "Europe/Berlin"
minimum_hours <- 18L
device_levels <- c("CRDS8", "PRONOVA", "CUBIC")
device_colors <- c(CRDS8 = "#4B5563", PRONOVA = "#7570B3", CUBIC = "#1B9E77")

data <- read_csv(input_path, show_col_types = FALSE, col_types = cols(date_time = col_character())) |>
  mutate(
    date_time = with_tz(ymd_hms(date_time, quiet = TRUE), time_zone),
    analyzer = recode(device, crds = "CRDS8", cubic = "CUBIC", pronova = "PRONOVA"),
    date_day = as.Date(date_time, tz = time_zone)
  ) |>
  filter(analyzer %in% device_levels)

daily_stats <- data |>
  select(date_day, analyzer, delta_co2, delta_ch4, delta_nh3) |>
  pivot_longer(starts_with("delta_"), names_to = "gas", values_to = "value") |>
  group_by(date_day, analyzer, gas) |>
  summarise(
    n_hours = sum(is.finite(value)),
    mean = if_else(n_hours > 0, mean(value, na.rm = TRUE), NA_real_),
    sd = if_else(n_hours > 1, sd(value, na.rm = TRUE), NA_real_),
    .groups = "drop"
  ) |>
  mutate(
    sufficient_coverage = n_hours >= minimum_hours,
    mean = if_else(sufficient_coverage, mean, NA_real_),
    sd = if_else(sufficient_coverage, sd, NA_real_),
    lower = mean - sd,
    upper = mean + sd,
    analyzer = factor(analyzer, levels = device_levels),
    gas = factor(
      gas,
      levels = c("delta_co2", "delta_ch4", "delta_nh3"),
      labels = c("Delta CO2 (ppm)", "Delta CH4 (ppm)", "Delta NH3 (ppm)")
    )
  ) |>
  arrange(gas, analyzer, date_day) |>
  group_by(gas, analyzer) |>
  mutate(
    new_segment = is.na(mean) | is.na(lag(mean)) | is.na(lag(date_day)) | as.integer(date_day - lag(date_day)) != 1L,
    segment = cumsum(replace_na(new_segment, TRUE))
  ) |>
  ungroup()

plot_data <- filter(daily_stats, sufficient_coverage, is.finite(mean), is.finite(sd))

p <- ggplot(
  plot_data,
  aes(
    x = as.POSIXct(date_day, tz = time_zone),
    y = mean,
    colour = analyzer,
    fill = analyzer,
    group = interaction(analyzer, segment)
  )
) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.14, colour = NA) +
  geom_line(linewidth = 0.9) +
  facet_grid(gas ~ ., scales = "free_y", switch = "y") +
  scale_colour_manual(values = device_colors, drop = FALSE) +
  scale_fill_manual(values = device_colors, drop = FALSE) +
  scale_x_datetime(
    date_breaks = "1 week",
    date_minor_breaks = "1 day",
    date_labels = "%d-%m-%Y",
    expand = expansion(mult = c(0.01, 0.01)),
    timezone = time_zone
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  labs(
    title = "CUBIC, PRONOVA and CRDS8 daily delta concentrations",
    subtitle = "Daily line = mean; shaded band = mean +/- within-day SD; at least 18 valid hourly values required",
    x = NULL,
    y = NULL,
    colour = NULL,
    fill = NULL
  ) +
  device_plot_theme_classic() +
  theme(
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.position = "bottom"
  ) +
  guides(fill = "none", colour = guide_legend(nrow = 1))

tag <- "2026081909_2026100609"
png_path <- file.path(output_dir, paste0("device_delta_daily_mean_sd_", tag, ".png"))
pdf_path <- file.path(output_dir, paste0("device_delta_daily_mean_sd_", tag, ".pdf"))
stats_path <- file.path(table_dir, paste0("device_delta_daily_mean_sd_", tag, ".csv"))

ggsave(png_path, p, width = 12, height = 7.5, dpi = 200, bg = "white")
ggsave(pdf_path, p, width = 12, height = 7.5)
write_csv(daily_stats, stats_path, na = "NA")

cat("PNG:", normalizePath(png_path, winslash = "/"), "\n")
cat("PDF:", normalizePath(pdf_path, winslash = "/"), "\n")
cat("Daily statistics:", normalizePath(stats_path, winslash = "/"), "\n")
print(
  daily_stats |>
    group_by(analyzer, gas) |>
    summarise(
      first_day = min(date_day[sufficient_coverage], na.rm = TRUE),
      last_day = max(date_day[sufficient_coverage], na.rm = TRUE),
      plotted_days = sum(sufficient_coverage),
      excluded_days = sum(!sufficient_coverage),
      .groups = "drop"
    )
)
