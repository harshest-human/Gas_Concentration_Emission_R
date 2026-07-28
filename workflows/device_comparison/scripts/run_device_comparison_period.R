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
helper_dir <- if (length(script_path) > 0) dirname(normalizePath(script_path[1], winslash = "/", mustWork = FALSE)) else getwd()
source(file.path(helper_dir, "device_comparison_helpers.R"))

default_device_old <- getOption("device")
options(device = function(...) {
  pdf(file = tempfile(fileext = ".pdf"), ...)
})
on.exit(options(device = default_device_old), add = TRUE)

parse_cli_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    parts <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- parts[1]
    value <- if (length(parts) > 1) paste(parts[-1], collapse = "=") else TRUE
    out[[key]] <- value
  }
  out
}

cli_args <- parse_cli_args(commandArgs(trailingOnly = TRUE))

timezone_local <- "Europe/Berlin"
period_start <- ymd_hms(if (!is.null(cli_args$start)) cli_args$start else "2025-12-09 12:00:00", tz = timezone_local)
period_end <- ymd_hms(if (!is.null(cli_args$end)) cli_args$end else "2025-12-22 23:00:00", tz = timezone_local)
output_suffix <- if (!is.null(cli_args$suffix)) cli_args$suffix else "20251209_20251222"
plot_subdir <- if (!is.null(cli_args$plot_subdir)) cli_args$plot_subdir else "dec_2025"
label_text <- if (!is.null(cli_args$label)) cli_args$label else "09-12-2025 12:00 to 22-12-2025"
ndir_clean_name <- if (!is.null(cli_args$ndir_clean)) cli_args$ndir_clean else "logas_ndir_hourly_dec_2025.csv"
tdlas_clean_name <- if (!is.null(cli_args$tdlas_clean)) cli_args$tdlas_clean else "logas_tdlas_hourly_dec_2025.csv"
otice_clean_name <- if (!is.null(cli_args$otice_clean)) cli_args$otice_clean else "otice_hourly_for_comparison_dec_2025.csv"

in_period <- function(x) {
  x >= period_start & x <= period_end
}

finalize_device_data <- function(df) {
  df |>
    filter(in_period(DATE.HOUR)) |>
    select(DATE.HOUR, analyzer, delta_CO2, delta_CH4, delta_NH3) |>
    remove_outliers(exclude_cols = c("DATE.HOUR", "analyzer"), group_cols = c("DATE.HOUR"))
}

read_clean_device_file <- function(path) {
  if (!file.exists(path)) {
    stop("Missing clean input file: ", path)
  }

  read_csv(path, show_col_types = FALSE)
}

rename_device_levels <- function(df) {
  df |>
    mutate(
      analyzer = recode(analyzer, "CRDS" = "crds", .default = analyzer),
      analyzer = factor(analyzer, levels = names(device_colors))
    )
}

device_trend_plot <- function(data, y, title_text) {
  facet_labels <- c(
    delta_CO2 = "Delta CO2 (ppm)",
    delta_CH4 = "Delta CH4 (ppm)",
    delta_NH3 = "Delta NH3 (ppm)",
    Q_vent = "Q (m3 h-1 LU-1)",
    e_CH4_ghLU = "eCH4 (g h-1 LU-1)",
    e_NH3_ghLU = "eNH3 (g h-1 LU-1)"
  )

  plot_data <- data |>
    filter(var %in% y) |>
    group_by(DATE.TIME, analyzer, var) |>
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop") |>
    filter(!is.na(value)) |>
    mutate(
      var = factor(var, levels = y, labels = facet_labels[y]),
      analyzer = factor(analyzer, levels = names(device_colors), labels = unname(device_labels[names(device_colors)]))
    )

  ggplot(plot_data, aes(x = DATE.TIME, y = value, color = analyzer, group = analyzer)) +
    geom_line(linewidth = 0.8, alpha = 0.95, na.rm = TRUE) +
    facet_grid(var ~ ., scales = "free_y", switch = "y") +
    scale_color_manual(values = setNames(unname(device_colors), unname(device_labels[names(device_colors)])), drop = FALSE) +
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
    device_plot_theme_classic()
}

project_dir <- resolve_project_dir(require_crds = TRUE, require_utils = TRUE)
source(file.path(project_dir, "scripts", "utils", "remove_outliers_function.R"))
source(file.path(project_dir, "scripts", "utils", "indirect.CO2.balance function.R"))

result_dir <- file.path(project_dir, "workflows", "device_comparison", "result_data")
plot_dir <- file.path(result_dir, "plots")
table_dir <- file.path(result_dir, "tables")
plot_output_dir <- file.path(plot_dir, plot_subdir)
dir.create(plot_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

rplots_pdf_path <- file.path(project_dir, "Rplots.pdf")
if (file.exists(rplots_pdf_path)) {
  unlink(rplots_pdf_path, force = TRUE)
}

load_crds_hourly <- function(project_dir) {
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

  map_dfr(
    crds_files,
    ~ read_csv(.x, col_types = cols(.default = col_character()))
  ) |>
    mutate(
      DATE.HOUR = parse_datetime_local(DATE.HOUR),
      across(any_of(c("CO2_in", "CO2_S", "CH4_in", "CH4_S", "NH3_in", "NH3_S")), as.numeric)
    ) |>
    filter(in_period(DATE.HOUR))
}

crds_hourly <- load_crds_hourly(project_dir)

crds_data <- crds_hourly |>
  mutate(
    delta_CO2 = CO2_in - CO2_S,
    delta_CH4 = CH4_in - CH4_S,
    delta_NH3 = NH3_in - NH3_S,
    analyzer = "crds"
  ) |>
  finalize_device_data()

logas_ndir_data <- read_clean_device_file(
  file.path(project_dir, "workflows", "device_comparison", "clean_data", "logas_ndir_clean", ndir_clean_name)
) |>
  mutate(
    DATE.HOUR = parse_datetime_local(DATE.HOUR),
    analyzer = "logas_ndir",
    across(any_of(c("delta_CO2", "delta_CH4", "delta_NH3")), as.numeric)
  ) |>
  finalize_device_data()

logas_tdlas_data <- read_clean_device_file(
  file.path(project_dir, "workflows", "device_comparison", "clean_data", "logas_tdlas_clean", tdlas_clean_name)
) |>
  mutate(
    DATE.HOUR = parse_datetime_local(DATE.HOUR),
    analyzer = "logas_tdlas",
    across(any_of(c("delta_CO2", "delta_CH4", "delta_NH3")), as.numeric)
  ) |>
  finalize_device_data()

otice_data <- read_clean_device_file(
  file.path(project_dir, "workflows", "device_comparison", "clean_data", "otice_clean", otice_clean_name)
) |>
  mutate(
    DATE.HOUR = parse_datetime_local(DATE.HOUR),
    analyzer = "otice",
    across(any_of(c("delta_CO2", "delta_CH4", "delta_NH3")), as.numeric)
  ) |>
  finalize_device_data()

otice_raw_data <- read_clean_device_file(
  file.path(project_dir, "workflows", "device_comparison", "clean_data", "otice_clean", otice_clean_name)
) |>
  mutate(
    DATE.HOUR = parse_datetime_local(DATE.HOUR),
    analyzer = "otice_raw",
    delta_CO2 = as.numeric(delta_CO2_raw),
    delta_CH4 = NA_real_,
    delta_NH3 = as.numeric(delta_NH3_raw)
  ) |>
  finalize_device_data()

gas_data <- bind_rows(crds_data, logas_ndir_data, logas_tdlas_data, otice_data) |>
  rename_device_levels() |>
  arrange(DATE.HOUR)

animal_data <- read_excel(
  file.path(project_dir, "shared_data", "clean_data", "animal_clean", "animal_data_2025-01-10_2026-01-19.xlsx")
) |>
  mutate(DATE.HOUR = parse_datetime_local(DATE.HOUR)) |>
  filter(in_period(DATE.HOUR))

temp_data <- list.files(
  file.path(project_dir, "shared_data", "clean_data", "temp_rh_clean"),
  pattern = "\\.csv$",
  full.names = TRUE,
  recursive = TRUE
) |>
  map_dfr(~ read_csv(.x, show_col_types = FALSE)) |>
  select(Date, T_inside) |>
  rename(temp_in = T_inside) |>
  mutate(
    DATE.TIME = parse_datetime_local(sub(" \\+0000$", "", Date)),
    DATE.HOUR = floor_date(DATE.TIME, "hour")
  ) |>
  group_by(DATE.HOUR) |>
  summarise(temp_in = safe_mean(temp_in), .groups = "drop") |>
  filter(in_period(DATE.HOUR))

input_data <- gas_data |>
  left_join(animal_data, by = "DATE.HOUR", relationship = "many-to-many") |>
  left_join(temp_data, by = "DATE.HOUR", relationship = "many-to-many") |>
  rename(DATE.TIME = DATE.HOUR)

emission_data <- indirect.CO2.balance(input_data)
emission_reshaped <- reshaper(emission_data)

write_csv(gas_data, file.path(table_dir, paste0("device_delta_hourly_", output_suffix, ".csv")))
write_csv(input_data, file.path(table_dir, paste0("emission_input_hourly_", output_suffix, ".csv")))
write_csv(emission_data, file.path(table_dir, paste0("emission_data_", output_suffix, ".csv")))
write_csv(emission_reshaped, file.path(table_dir, paste0("emission_reshaped_", output_suffix, ".csv")))

daily_summary <- emission_data |>
  mutate(date_day = as.Date(DATE.TIME)) |>
  group_by(date_day, analyzer) |>
  summarise(
    across(c(delta_CO2, delta_CH4, delta_NH3, Q_vent, e_CH4_ghLU, e_NH3_ghLU), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

write_csv(daily_summary, file.path(table_dir, paste0("daily_summary_", output_suffix, ".csv")))

delta_summary_table <- bind_rows(crds_data, logas_ndir_data, logas_tdlas_data, otice_raw_data, otice_data) |>
  mutate(
    analyzer = as.character(analyzer),
    analyzer = recode(analyzer, !!!device_labels, .default = analyzer)
  ) |>
  pivot_longer(
    cols = c(delta_CO2, delta_CH4, delta_NH3),
    names_to = "variable",
    values_to = "value"
  ) |>
  filter(!is.na(value)) |>
  group_by(analyzer, variable) |>
  summarise(
    n_hours = n(),
    mean_value = mean(value, na.rm = TRUE),
    sd_value = sd(value, na.rm = TRUE),
    se_value = sd_value / sqrt(n_hours),
    summary_text = paste0(round(mean_value, 3), " +/- ", round(se_value, 3)),
    .groups = "drop"
  ) |>
  arrange(variable, analyzer)

precision_vs_crds_table <- bind_rows(crds_data, logas_ndir_data, logas_tdlas_data, otice_raw_data, otice_data) |>
  mutate(
    analyzer = as.character(analyzer),
    analyzer = recode(analyzer, !!!device_labels, .default = analyzer)
  ) |>
  pivot_longer(
    cols = c(delta_CO2, delta_CH4, delta_NH3),
    names_to = "variable",
    values_to = "value"
  ) |>
  select(DATE.HOUR, analyzer, variable, value) |>
  pivot_wider(names_from = analyzer, values_from = value) |>
  pivot_longer(
    cols = any_of(c("PRONOVA", "CUBIC", "OTICE raw", "OTICE")),
    names_to = "compare_analyzer",
    values_to = "compare_value"
  ) |>
  filter(!is.na(CRDS), !is.na(compare_value)) |>
  mutate(diff_vs_crds = compare_value - CRDS) |>
  group_by(compare_analyzer, variable) |>
  summarise(
    paired_hours = n(),
    mean_diff = mean(diff_vs_crds, na.rm = TRUE),
    sd_diff = sd(diff_vs_crds, na.rm = TRUE),
    se_diff = sd_diff / sqrt(paired_hours),
    mae = mean(abs(diff_vs_crds), na.rm = TRUE),
    rmse = sqrt(mean(diff_vs_crds^2, na.rm = TRUE)),
    summary_text = paste0(round(mean_diff, 3), " +/- ", round(se_diff, 3)),
    .groups = "drop"
  ) |>
  arrange(variable, compare_analyzer)

write_csv(delta_summary_table, file.path(table_dir, paste0("device_delta_summary_", output_suffix, ".csv")))
write_csv(precision_vs_crds_table, file.path(table_dir, paste0("device_delta_precision_vs_crds_", output_suffix, ".csv")))

delta_trend_plot <- device_trend_plot(
  emission_reshaped,
  y = c("delta_CO2", "delta_CH4", "delta_NH3"),
  title_text = paste("Device delta concentration trends,", label_text)
)

emission_trend_plot <- device_trend_plot(
  emission_reshaped,
  y = c("Q_vent", "e_CH4_ghLU", "e_NH3_ghLU"),
  title_text = paste("Device emission trends,", label_text)
)

tmp_plot_device <- tempfile(fileext = ".pdf")
pdf(tmp_plot_device)
delta_errorbar_plot <- emierrorbarplot(emission_reshaped, y = c("delta_CO2", "delta_CH4", "delta_NH3"))
emission_errorbar_plot <- emierrorbarplot(emission_reshaped, y = c("Q_vent", "e_CH4_ghLU", "e_NH3_ghLU"))
invisible(dev.off())
unlink(tmp_plot_device, force = TRUE)

ggsave(file.path(plot_output_dir, paste0("device_delta_trends_", output_suffix, ".png")), delta_trend_plot, width = 11, height = 6.5, dpi = 150)
ggsave(file.path(plot_output_dir, paste0("device_emission_trends_", output_suffix, ".png")), emission_trend_plot, width = 11, height = 6.5, dpi = 150)
ggsave(file.path(plot_output_dir, paste0("device_delta_errorbars_", output_suffix, ".png")), delta_errorbar_plot, width = 7.5, height = 6, dpi = 150)
ggsave(file.path(plot_output_dir, paste0("device_emission_errorbars_", output_suffix, ".png")), emission_errorbar_plot, width = 7.5, height = 6, dpi = 150)

if (file.exists(rplots_pdf_path)) {
  unlink(rplots_pdf_path, force = TRUE)
}

cat("Wrote device comparison outputs to:\n", normalizePath(result_dir, winslash = "/", mustWork = FALSE), "\n")
cat("Rows by analyzer:\n")
print(count(gas_data, analyzer))
