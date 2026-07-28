library(dplyr)
library(lubridate)
library(readr)
library(tidyr)
library(ggplot2)
library(purrr)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args_all[grepl(file_arg, args_all)])
helper_dir <- if (length(script_path) > 0) dirname(normalizePath(script_path[1], winslash = "/", mustWork = FALSE)) else getwd()
source(file.path(helper_dir, "device_comparison_helpers.R"))

campaign_start <- as.POSIXct("2025-09-29 00:00:00", tz = "Europe/Berlin")

load_crds_hourly <- function(project_dir) {
  crds_dir <- file.path(project_dir, "workflows", "crds_routine_cleaning", "clean_data", "crds_clean", "hourly_in_out_avg")
  crds_files <- list.files(crds_dir, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
  crds_files <- crds_files[
    grepl("CRDS", basename(crds_files)) &
      !grepl("FTIR|LUFA|UB|ANECO|MBBM", basename(crds_files))
  ]

  bind_rows(lapply(crds_files, function(path) {
    read_csv(path, col_types = cols(.default = col_character()), show_col_types = FALSE)
  })) %>%
    mutate(
      DATE.HOUR = parse_datetime_local(DATE.HOUR),
      CO2_in = as.numeric(CO2_in),
      CO2_S = as.numeric(CO2_S),
      CH4_in = as.numeric(CH4_in),
      CH4_S = as.numeric(CH4_S),
      NH3_in = as.numeric(NH3_in),
      NH3_S = as.numeric(NH3_S)
    ) %>%
    filter(!is.na(DATE.HOUR), DATE.HOUR >= campaign_start) %>%
    group_by(DATE.HOUR) %>%
    summarise(
      CO2_in_ppm = safe_mean(CO2_in),
      CH4_in_ppm = safe_mean(CH4_in),
      NH3_in_ppm = safe_mean(NH3_in),
      CO2_s_ppm = safe_mean(CO2_S),
      CH4_s_ppm = safe_mean(CH4_S),
      NH3_s_ppm = safe_mean(NH3_S),
      CO2_delta_ppm = safe_mean(CO2_in - CO2_S),
      CH4_delta_ppm = safe_mean(CH4_in - CH4_S),
      NH3_delta_ppm = safe_mean(NH3_in - NH3_S),
      .groups = "drop"
    ) %>%
    mutate(analyzer = "CRDS") %>%
    select(
      DATE.HOUR, analyzer,
      CO2_in_ppm, CH4_in_ppm, NH3_in_ppm,
      CO2_s_ppm, CH4_s_ppm, NH3_s_ppm,
      CO2_delta_ppm, CH4_delta_ppm, NH3_delta_ppm
    )
}

load_pronova_hourly <- function(project_dir) {
  read_csv(
    file.path(project_dir, "workflows", "device_comparison", "clean_data", "logas_ndir_campaign", "logas_ndir_campaign_hourly.csv"),
    show_col_types = FALSE
  ) %>%
    mutate(
      DATE.HOUR = parse_datetime_local(DATE.HOUR),
      CO2_in_ppm = NA_real_,
      CH4_in_ppm = NA_real_,
      NH3_in_ppm = NA_real_,
      CO2_s_ppm = NA_real_,
      CH4_s_ppm = NA_real_,
      NH3_s_ppm = NA_real_,
      CO2_delta_ppm = as.numeric(delta_CO2),
      CH4_delta_ppm = as.numeric(delta_CH4),
      NH3_delta_ppm = as.numeric(delta_NH3),
      analyzer = "PRONOVA"
    ) %>%
    filter(!is.na(DATE.HOUR), DATE.HOUR >= campaign_start) %>%
    select(
      DATE.HOUR, analyzer,
      CO2_in_ppm, CH4_in_ppm, NH3_in_ppm,
      CO2_s_ppm, CH4_s_ppm, NH3_s_ppm,
      CO2_delta_ppm, CH4_delta_ppm, NH3_delta_ppm
    )
}

load_cubic_hourly <- function(project_dir) {
  read_csv(
    file.path(project_dir, "workflows", "device_comparison", "clean_data", "logas_tdlas_campaign", "logas_tdlas_campaign_hourly.csv"),
    show_col_types = FALSE
  ) %>%
    mutate(
      DATE.HOUR = parse_datetime_local(DATE.HOUR),
      CO2_in_ppm = as.numeric(CO2_in),
      CH4_in_ppm = as.numeric(CH4_in),
      NH3_in_ppm = as.numeric(NH3_in),
      CO2_s_ppm = as.numeric(CO2_S),
      CH4_s_ppm = as.numeric(CH4_S),
      NH3_s_ppm = as.numeric(NH3_S),
      CO2_delta_ppm = as.numeric(delta_CO2),
      CH4_delta_ppm = as.numeric(delta_CH4),
      NH3_delta_ppm = as.numeric(delta_NH3),
      analyzer = "CUBIC"
    ) %>%
    filter(!is.na(DATE.HOUR), DATE.HOUR >= campaign_start) %>%
    select(
      DATE.HOUR, analyzer,
      CO2_in_ppm, CH4_in_ppm, NH3_in_ppm,
      CO2_s_ppm, CH4_s_ppm, NH3_s_ppm,
      CO2_delta_ppm, CH4_delta_ppm, NH3_delta_ppm
    )
}

load_otice_hourly <- function(project_dir) {
  read_csv(
    file.path(project_dir, "workflows", "device_comparison", "clean_data", "otice_campaign", "otice_campaign_hourly.csv"),
    show_col_types = FALSE
  ) %>%
    mutate(
      DATE.HOUR = parse_datetime_local(DATE.HOUR),
      delta_CO2_fitted = as.numeric(delta_CO2_fitted),
      delta_NH3_fitted = as.numeric(delta_NH3_fitted)
    ) %>%
    filter(!is.na(DATE.HOUR), DATE.HOUR >= campaign_start) %>%
    group_by(DATE.HOUR) %>%
    summarise(
      CO2_in_ppm = NA_real_,
      CH4_in_ppm = NA_real_,
      NH3_in_ppm = NA_real_,
      CO2_s_ppm = NA_real_,
      CH4_s_ppm = NA_real_,
      NH3_s_ppm = NA_real_,
      CO2_delta_ppm = safe_mean(delta_CO2_fitted),
      CH4_delta_ppm = NA_real_,
      NH3_delta_ppm = safe_mean(delta_NH3_fitted),
      .groups = "drop"
    ) %>%
    mutate(analyzer = "OTICE") %>%
    select(
      DATE.HOUR, analyzer,
      CO2_in_ppm, CH4_in_ppm, NH3_in_ppm,
      CO2_s_ppm, CH4_s_ppm, NH3_s_ppm,
      CO2_delta_ppm, CH4_delta_ppm, NH3_delta_ppm
    )
}

drop_empty_rows <- function(df) {
  df %>%
    filter(!(is.na(CO2_delta_ppm) & is.na(CH4_delta_ppm) & is.na(NH3_delta_ppm)))
}

prepare_for_export <- function(df) {
  df %>%
    mutate(DATE.HOUR = format(with_tz(DATE.HOUR, tzone = "Europe/Berlin"), "%Y-%m-%d %H:%M:%S"))
}

build_month_lookup <- function(df) {
  month_starts <- seq(
    from = floor_date(min(df$DATE.HOUR, na.rm = TRUE), "month"),
    to = floor_date(max(df$DATE.HOUR, na.rm = TRUE), "month"),
    by = "month"
  )

  tibble(
    month_id = format(month_starts, "%Y-%m"),
    month_label = format(month_starts, "%B %Y"),
    start = month_starts,
    end = (month_starts %m+% months(1)) - hours(1)
  )
}

build_pair_plot <- function(df, title_text) {
  analyzer_levels <- c("CRDS", "PRONOVA", "CUBIC", "OTICE")
  analyzer_labels <- c(
    "CRDS" = "CRDS",
    "PRONOVA" = "PRONOVA",
    "CUBIC" = "CUBIC",
    "OTICE" = "OTICE"
  )
  color_map <- c(
    "CRDS" = unname(device_colors["crds"]),
    "PRONOVA" = unname(device_colors["logas_ndir"]),
    "CUBIC" = unname(device_colors["logas_tdlas"]),
    "OTICE" = unname(device_colors["otice"])
  )

  plot_df <- df %>%
    pivot_longer(cols = c(CO2_delta_ppm, CH4_delta_ppm, NH3_delta_ppm), names_to = "gas", values_to = "value") %>%
    filter(!is.na(value)) %>%
    mutate(
      gas = factor(
        gas,
        levels = c("CO2_delta_ppm", "CH4_delta_ppm", "NH3_delta_ppm"),
        labels = c("Delta CO2 (ppm)", "Delta CH4 (ppm)", "Delta NH3 (ppm)")
      ),
      analyzer = factor(analyzer, levels = analyzer_levels)
    )

  ggplot(plot_df, aes(x = DATE.HOUR, y = value, color = analyzer)) +
    geom_line(linewidth = 0.6, na.rm = TRUE) +
    facet_grid(gas ~ ., scales = "free_y", switch = "y") +
    scale_color_manual(
      values = color_map,
      breaks = unique(as.character(na.omit(plot_df$analyzer))),
      labels = analyzer_labels[unique(as.character(na.omit(plot_df$analyzer)))],
      drop = TRUE
    ) +
    labs(title = title_text, x = NULL, y = NULL, color = NULL) +
    device_plot_theme_classic()
}

write_monthly_pair_files <- function(df, dataset_name, output_dir) {
  month_lookup <- build_month_lookup(df)
  dataset_dir <- file.path(output_dir, dataset_name)
  dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)

  for (i in seq_len(nrow(month_lookup))) {
    month_row <- month_lookup[i, ]
    month_df <- df %>%
      filter(DATE.HOUR >= month_row$start, DATE.HOUR <= month_row$end)

    if (nrow(month_df) == 0) {
      next
    }

    write_csv(
      prepare_for_export(month_df),
      file.path(dataset_dir, paste0(dataset_name, "_", gsub("-", "_", month_row$month_id), ".csv"))
    )
  }
}

write_monthly_pair_plots <- function(df, dataset_name, plot_dir) {
  month_lookup <- build_month_lookup(df)
  dataset_plot_dir <- file.path(plot_dir, dataset_name)
  dir.create(dataset_plot_dir, recursive = TRUE, showWarnings = FALSE)

  full_plot <- build_pair_plot(df, paste(dataset_name, "vs CRDS"))
  ggsave(
    file.path(dataset_plot_dir, paste0(dataset_name, "_full_campaign.png")),
    full_plot,
    width = 11,
    height = 7,
    dpi = 150
  )

  for (i in seq_len(nrow(month_lookup))) {
    month_row <- month_lookup[i, ]
    month_df <- df %>%
      filter(DATE.HOUR >= month_row$start, DATE.HOUR <= month_row$end)

    if (nrow(month_df) == 0) {
      next
    }

    month_plot <- build_pair_plot(df = month_df, title_text = paste(dataset_name, "vs CRDS,", month_row$month_label))
    ggsave(
      file.path(dataset_plot_dir, paste0(dataset_name, "_", gsub("-", "_", month_row$month_id), ".png")),
      month_plot,
      width = 11,
      height = 7,
      dpi = 150
    )
  }
}

project_dir <- resolve_project_dir(require_crds = TRUE, require_utils = TRUE)
result_table_dir <- file.path(project_dir, "workflows", "device_comparison", "result_data", "tables", "analyzer_comparison")
result_plot_dir <- file.path(project_dir, "workflows", "device_comparison", "result_data", "plots", "analyzer_comparison")
monthly_table_dir <- file.path(result_table_dir, "monthly")

dir.create(result_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(monthly_table_dir, recursive = TRUE, showWarnings = FALSE)

main_dataset <- bind_rows(
  load_crds_hourly(project_dir),
  load_pronova_hourly(project_dir),
  load_cubic_hourly(project_dir),
  load_otice_hourly(project_dir)
) %>%
  drop_empty_rows() %>%
  arrange(DATE.HOUR, analyzer)

write_csv(
  prepare_for_export(main_dataset),
  file.path(result_table_dir, "device_delta_hourly_main_2025_09_29_onward.csv")
)

pair_datasets <- list(
  pronova_vs_crds = main_dataset %>% filter(analyzer %in% c("CRDS", "PRONOVA")),
  cubic_vs_crds = main_dataset %>% filter(analyzer %in% c("CRDS", "CUBIC")),
  otice_vs_crds = main_dataset %>% filter(analyzer %in% c("CRDS", "OTICE"))
)

for (dataset_name in names(pair_datasets)) {
  dataset_df <- pair_datasets[[dataset_name]] %>% arrange(DATE.HOUR, analyzer)

  write_csv(
    prepare_for_export(dataset_df),
    file.path(result_table_dir, paste0(dataset_name, "_hourly.csv"))
  )

  write_monthly_pair_files(dataset_df, dataset_name, monthly_table_dir)
  write_monthly_pair_plots(dataset_df, dataset_name, result_plot_dir)
}

cat("Wrote main analyzer comparison dataset to:\n", normalizePath(file.path(result_table_dir, "device_delta_hourly_main_2025_09_29_onward.csv"), winslash = "/", mustWork = FALSE), "\n")
cat("Wrote analyzer-specific datasets to:\n", normalizePath(result_table_dir, winslash = "/", mustWork = FALSE), "\n")
cat("Wrote analyzer comparison plots to:\n", normalizePath(result_plot_dir, winslash = "/", mustWork = FALSE), "\n")
