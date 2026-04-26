# =============================================================================
# Device trend comparison
#
# Period: 2025-12-09 to 2025-12-22
# Devices: crds reference, logas_ndir, logas_tdlas, otice
# =============================================================================

library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(readr)
library(readxl)
library(purrr)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args_all[grepl(file_arg, args_all)])
project_dir <- if (length(script_path) > 0) {
  normalizePath(file.path(dirname(script_path[1]), "..", ".."), winslash = "/", mustWork = FALSE)
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

source(file.path(project_dir, "scripts", "utils", "remove_outliers_function.R"))
source(file.path(project_dir, "scripts", "utils", "indirect.CO2.balance function.R"))

period_start <- ymd_hms("2025-12-09 12:00:00")
period_end <- ymd_hms("2025-12-22 23:00:00")

out_dir <- file.path(project_dir, "workflows", "device_comparison", "result_data")
plot_dir <- file.path(out_dir, "plots")
table_dir <- file.path(out_dir, "tables")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

in_period <- function(x) {
  x >= period_start & x <= period_end
}

parse_hour <- function(x) {
  parsed <- suppressWarnings(ymd_hms(as.character(x), tz = "Europe/Berlin"))
  missing <- is.na(parsed)
  parsed[missing] <- suppressWarnings(parse_date_time(
    as.character(x),
    orders = c("ymd HMS", "ymd HM", "Ymd HMS", "Ymd HM"),
    tz = "Europe/Berlin"
  ))[missing]
  parsed
}

device_colors <- c(
  "crds" = "#6B7280",
  "logas_tdlas" = "#1B9E77",
  "logas_ndir" = "#7570B3",
  "otice" = "#C45A11"
)

rename_device_levels <- function(df) {
  df |>
    mutate(
      analyzer = recode(analyzer, "CRDS" = "crds", .default = analyzer),
      analyzer = factor(analyzer, levels = names(device_colors))
    )
}

device_trend_plot <- function(data, y, title_text) {
  facet_labels <- c(
    "delta_CO2" = "Delta CO2 (ppm)",
    "delta_CH4" = "Delta CH4 (ppm)",
    "delta_NH3" = "Delta NH3 (ppm)",
    "Q_vent" = "Q (m3 h-1 LU-1)",
    "e_CH4_ghLU" = "eCH4 (g h-1 LU-1)",
    "e_NH3_ghLU" = "eNH3 (g h-1 LU-1)"
  )

  plot_data <- data |>
    filter(var %in% y) |>
    group_by(DATE.TIME, analyzer, var) |>
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop") |>
    filter(!is.na(value)) |>
    mutate(
      var = factor(var, levels = y, labels = facet_labels[y])
    )

  ggplot(plot_data, aes(x = DATE.TIME, y = value, color = analyzer, group = analyzer)) +
    geom_line(linewidth = 0.8, alpha = 0.95, na.rm = TRUE) +
    facet_grid(var ~ ., scales = "free_y", switch = "y") +
    scale_color_manual(values = device_colors, drop = FALSE) +
    scale_x_datetime(
      limits = c(period_start, period_end),
      breaks = seq(period_start, period_end, by = "1 day"),
      date_labels = "%d-%m-%Y",
      expand = expansion(mult = c(0, 0.01))
    ) +
    scale_y_continuous(
      breaks = scales::pretty_breaks(n = 7),
      labels = scales::label_number(accuracy = 0.1, big.mark = "")
    ) +
    labs(title = title_text, x = NULL, y = NULL, color = NULL) +
    theme_classic(base_size = 15) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", hjust = 0.5, size = 17),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
      axis.text.y = element_text(size = 13),
      strip.text.y.left = element_text(size = 14),
      panel.border = element_rect(color = "black", fill = NA)
    )
}

# -----------------------------------------------------------------------------
# 1. crds reference
# -----------------------------------------------------------------------------
crds_files <- list.files(
  file.path(project_dir, "workflows", "crds_routine_cleaning", "clean_data", "crds_clean", "hourly_in_out_avg"),
  pattern = "\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)

crds_files <- crds_files[
  grepl("CRDS", basename(crds_files)) &
    !grepl("FTIR|LUFA|UB|ANECO|MBBM", basename(crds_files))
]

crds_data <- map_dfr(
  crds_files,
  ~ read_csv(.x, col_types = cols(.default = col_character()))
) |>
  mutate(
    DATE.HOUR = parse_hour(DATE.HOUR),
    across(any_of(c("CO2_in", "CO2_S", "CH4_in", "CH4_S", "NH3_in", "NH3_S")), as.numeric),
    delta_CO2 = CO2_in - CO2_S,
    delta_CH4 = CH4_in - CH4_S,
    delta_NH3 = NH3_in - NH3_S,
    analyzer = "crds"
  ) |>
  filter(in_period(DATE.HOUR)) |>
  select(DATE.HOUR, analyzer, delta_CO2, delta_CH4, delta_NH3) |>
  remove_outliers(exclude_cols = c("DATE.HOUR", "analyzer"), group_cols = c("DATE.HOUR"))

# -----------------------------------------------------------------------------
# 2. logas_ndir
# -----------------------------------------------------------------------------
logas_ndir_files <- list.files(
  path = file.path(project_dir, "workflows", "device_comparison", "raw_data", "logas_ndir_raw", "Messdaten"),
  pattern = "^differenzmessung_.*\\.txt$",
  full.names = TRUE
)

logas_ndir_read_log_file <- function(file) {
  tryCatch(
    read.table(
      file,
      header = TRUE,
      sep = "\t",
      dec = ",",
      check.names = FALSE,
      stringsAsFactors = FALSE,
      fill = TRUE,
      comment.char = "",
      colClasses = "character"
    ),
    error = function(e) NULL
  )
}

logas_ndir_data <- logas_ndir_files |>
  lapply(logas_ndir_read_log_file) |>
  bind_rows() |>
  mutate(
    DATE.TIME = dmy_hms(`Datum Uhrzeit`),
    DATE.HOUR = floor_date(DATE.TIME, "hour")
  ) |>
  group_by(DATE.HOUR) |>
  summarise(
    delta_CH4 = mean(as.numeric(gsub(",", ".", `CH4 in ppm`)), na.rm = TRUE),
    delta_CO2 = mean(as.numeric(gsub(",", ".", `CO2 in ppm`)), na.rm = TRUE),
    delta_NH3 = mean(as.numeric(gsub(",", ".", `NH3 in ppm`)), na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(analyzer = "logas_ndir") |>
  filter(in_period(DATE.HOUR)) |>
  select(DATE.HOUR, analyzer, delta_CO2, delta_CH4, delta_NH3) |>
  remove_outliers(exclude_cols = c("DATE.HOUR", "analyzer"), group_cols = c("DATE.HOUR"))

# -----------------------------------------------------------------------------
# 3. logas_tdlas
# -----------------------------------------------------------------------------
logas_tdlas_files <- list.files(
  path = file.path(project_dir, "workflows", "device_comparison", "raw_data", "logas_tdlas_raw"),
  pattern = "^[^~].*\\.xlsx$",
  full.names = TRUE
)

logas_tdlas_data <- logas_tdlas_files |>
  map_dfr(read_excel) |>
  mutate(
    Time = ymd_hms(Time),
    DATE.HOUR = floor_date(Time, "hour")
  ) |>
  pivot_longer(cols = any_of(c("CH4", "NH3", "CO2")), names_to = "gas", values_to = "value") |>
  mutate(gas_name = paste0(gas, ifelse(Type == 1, "_in", "_S"))) |>
  group_by(DATE.HOUR, gas_name) |>
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = gas_name, values_from = value) |>
  mutate(
    delta_CO2 = CO2_in - CO2_S,
    delta_CH4 = CH4_in - CH4_S,
    delta_NH3 = NH3_in - NH3_S,
    analyzer = "logas_tdlas"
  ) |>
  filter(in_period(DATE.HOUR)) |>
  select(DATE.HOUR, analyzer, delta_CO2, delta_CH4, delta_NH3) |>
  remove_outliers(exclude_cols = c("DATE.HOUR", "analyzer"), group_cols = c("DATE.HOUR"))

# -----------------------------------------------------------------------------
# 4. otice calibrated node average
# -----------------------------------------------------------------------------
otice_data <- read_csv(
  file.path(project_dir, "workflows", "device_comparison", "clean_data", "device_hourly", "otice_hourly_20251209_20251222.csv"),
  show_col_types = FALSE
) |>
  mutate(DATE.HOUR = parse_hour(DATE.HOUR)) |>
  filter(in_period(DATE.HOUR)) |>
  select(DATE.HOUR, analyzer, delta_CO2, delta_CH4, delta_NH3) |>
  remove_outliers(exclude_cols = c("DATE.HOUR", "analyzer"), group_cols = c("DATE.HOUR"))

# -----------------------------------------------------------------------------
# 5. Combine devices and calculate emissions
# -----------------------------------------------------------------------------
gas_data <- bind_rows(crds_data, logas_ndir_data, logas_tdlas_data, otice_data) |>
  rename_device_levels() |>
  arrange(DATE.HOUR)

animal_data <- read_excel(file.path(project_dir, "shared_data", "clean_data", "animal_clean", "animal_data_2025-01-10_2026-01-19.xlsx")) |>
  mutate(DATE.HOUR = parse_hour(DATE.HOUR)) |>
  filter(in_period(DATE.HOUR))

temp_data <- list.files(file.path(project_dir, "shared_data", "clean_data", "temp_rh_clean"), pattern = "\\.csv$", full.names = TRUE, recursive = TRUE) |>
  map_dfr(~ read_csv(.x, show_col_types = FALSE)) |>
  select(Date, T_inside) |>
  rename(temp_in = T_inside) |>
  mutate(
    DATE.TIME = mdy_hms(sub(" \\+0000$", "", Date)),
    DATE.HOUR = floor_date(DATE.TIME, "hour")
  ) |>
  group_by(DATE.HOUR) |>
  summarise(temp_in = mean(temp_in, na.rm = TRUE), .groups = "drop") |>
  filter(in_period(DATE.HOUR))

input_data <- gas_data |>
  left_join(animal_data, by = "DATE.HOUR", relationship = "many-to-many") |>
  left_join(temp_data, by = "DATE.HOUR", relationship = "many-to-many") |>
  rename(DATE.TIME = DATE.HOUR)

emission_data <- indirect.CO2.balance(input_data)
emission_reshaped <- reshaper(emission_data)

write_csv(gas_data, file.path(table_dir, "device_delta_hourly_20251209_20251222.csv"))
write_csv(input_data, file.path(table_dir, "emission_input_hourly_20251209_20251222.csv"))
write_csv(emission_data, file.path(table_dir, "emission_data_20251209_20251222.csv"))
write_csv(emission_reshaped, file.path(table_dir, "emission_reshaped_20251209_20251222.csv"))

daily_summary <- emission_data |>
  mutate(date_day = as.Date(DATE.TIME)) |>
  group_by(date_day, analyzer) |>
  summarise(
    across(c(delta_CO2, delta_CH4, delta_NH3, Q_vent, e_CH4_ghLU, e_NH3_ghLU), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

write_csv(daily_summary, file.path(table_dir, "daily_summary_20251209_20251222.csv"))

# -----------------------------------------------------------------------------
# 6. Plots
# -----------------------------------------------------------------------------
delta_trend_plot <- device_trend_plot(
  emission_reshaped,
  y = c("delta_CO2", "delta_CH4", "delta_NH3"),
  title_text = "Device delta concentration trends, 09-12-2025 12:00 to 22-12-2025"
)

emission_trend_plot <- device_trend_plot(
  emission_reshaped,
  y = c("Q_vent", "e_CH4_ghLU", "e_NH3_ghLU"),
  title_text = "Device emission trends, 09-12-2025 12:00 to 22-12-2025"
)

delta_errorbar_plot <- emierrorbarplot(emission_reshaped, y = c("delta_CO2", "delta_CH4", "delta_NH3")) +
  scale_color_manual(values = device_colors, drop = FALSE)
emission_errorbar_plot <- emierrorbarplot(emission_reshaped, y = c("Q_vent", "e_CH4_ghLU", "e_NH3_ghLU")) +
  scale_color_manual(values = device_colors, drop = FALSE)

ggsave(file.path(plot_dir, "device_delta_trends_20251209_20251222.png"), delta_trend_plot, width = 15, height = 9, dpi = 150)
ggsave(file.path(plot_dir, "device_emission_trends_20251209_20251222.png"), emission_trend_plot, width = 15, height = 9, dpi = 150)
ggsave(file.path(plot_dir, "device_delta_errorbars_20251209_20251222.png"), delta_errorbar_plot, width = 10, height = 8, dpi = 150)
ggsave(file.path(plot_dir, "device_emission_errorbars_20251209_20251222.png"), emission_errorbar_plot, width = 10, height = 8, dpi = 150)

cat("Wrote multi-device comparison outputs to:\n", normalizePath(out_dir), "\n")
cat("Rows by analyzer:\n")
print(count(gas_data, analyzer))
