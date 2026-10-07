library(dplyr)
library(ggplot2)
library(lubridate)
library(readr)
library(tidyr)

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args_all[grepl("--file=", args_all)])
script_dir <- dirname(normalizePath(script_path[1], winslash = "/", mustWork = TRUE))
source(file.path(script_dir, "device_comparison_helpers.R"))

workflow_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
input_path <- file.path(workflow_dir, "data", "integrated", "device_comparison_2026070600_2026082000.csv")
output_dir <- file.path(workflow_dir, "result_data", "plots", "cubic_pronova_crds8_since_20260706")
table_dir <- file.path(workflow_dir, "result_data", "tables")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

start_time <- ymd_hms("2026-07-06 00:00:00", tz = "Europe/Berlin")

data <- read_csv(input_path, show_col_types = FALSE) |>
  mutate(
    date_time = ymd_hms(date_time, tz = "Europe/Berlin"),
    analyzer = recode(device, crds = "CRDS8", cubic = "CUBIC", pronova = "PRONOVA")
  ) |>
  filter(date_time >= start_time, analyzer %in% c("CUBIC", "PRONOVA", "CRDS8"))

coverage <- data |>
  pivot_longer(c(delta_co2, delta_ch4, delta_nh3), names_to = "gas", values_to = "value") |>
  filter(!is.na(value)) |>
  group_by(analyzer, gas) |>
  summarise(first_hour = min(date_time), last_hour = max(date_time), n_hours = n(), .groups = "drop")

latest_time <- max(coverage$last_hour)
plot_data <- data |>
  pivot_longer(c(delta_co2, delta_ch4, delta_nh3), names_to = "gas", values_to = "value") |>
  filter(!is.na(value)) |>
  mutate(
    analyzer = factor(analyzer, levels = c("CRDS8", "PRONOVA", "CUBIC")),
    gas = factor(
      gas,
      levels = c("delta_co2", "delta_ch4", "delta_nh3"),
      labels = c("Delta CO2 (ppm)", "Delta CH4 (ppm)", "Delta NH3 (ppm)")
    )
  )

colors <- c(CRDS8 = "#4B5563", PRONOVA = "#7570B3", CUBIC = "#1B9E77")

p <- ggplot(plot_data, aes(date_time, value, color = analyzer, group = analyzer)) +
  geom_line(linewidth = 0.55, alpha = 0.95, na.rm = TRUE) +
  facet_grid(gas ~ ., scales = "free_y", switch = "y") +
  scale_color_manual(values = colors, drop = FALSE) +
  scale_x_datetime(
    limits = c(start_time, latest_time),
    date_breaks = "1 week",
    date_minor_breaks = "1 day",
    date_labels = "%d-%m-%Y",
    expand = expansion(mult = c(0, 0.01)),
    timezone = "Europe/Berlin"
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  labs(
    title = "CUBIC, PRONOVA and CRDS8 delta concentration trends",
    subtitle = paste0("06-07-2026 to latest available data (", format(latest_time, "%d-%m-%Y %H:%M"), ")"),
    x = NULL, y = NULL, color = NULL
  ) +
  device_plot_theme_classic() +
  theme(plot.subtitle = element_text(hjust = 0.5, size = 10))

tag <- paste0("20260706_", format(latest_time, "%Y%m%d%H"))
png_path <- file.path(output_dir, paste0("device_delta_trends_cubic_pronova_crds8_", tag, ".png"))
pdf_path <- file.path(output_dir, paste0("device_delta_trends_cubic_pronova_crds8_", tag, ".pdf"))
coverage_path <- file.path(table_dir, paste0("device_delta_coverage_cubic_pronova_crds8_", tag, ".csv"))

ggsave(png_path, p, width = 12, height = 7.2, dpi = 200)
ggsave(pdf_path, p, width = 12, height = 7.2)
write_csv(coverage, coverage_path)

cat("PNG:", normalizePath(png_path, winslash = "/"), "\n")
cat("PDF:", normalizePath(pdf_path, winslash = "/"), "\n")
cat("Coverage:", normalizePath(coverage_path, winslash = "/"), "\n")
print(coverage)
