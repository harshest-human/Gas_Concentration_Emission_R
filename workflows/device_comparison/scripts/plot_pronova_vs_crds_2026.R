library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(readr)
library(purrr)
library(tibble)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args_all[grepl(file_arg, args_all)])
helper_dir <- if (length(script_path) > 0) dirname(normalizePath(script_path[1], winslash = "/", mustWork = FALSE)) else getwd()
source(file.path(helper_dir, "device_comparison_helpers.R"))

load_crds_hourly <- function(project_dir, time_zone = "Europe/Berlin") {
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

  if (length(crds_files) == 0) {
    stop("No CRDS hourly files found under cleaned CRDS output.")
  }

  map_dfr(
    crds_files,
    ~ read_csv(.x, col_types = cols(.default = col_character()), show_col_types = FALSE)
  ) %>%
    mutate(
      DATE.HOUR = parse_datetime_local(DATE.HOUR, time_zone = time_zone),
      across(any_of(c("CO2_in", "CO2_S", "CH4_in", "CH4_S", "NH3_in", "NH3_S")), as.numeric)
    ) %>%
    filter(!is.na(DATE.HOUR)) %>%
    group_by(DATE.HOUR) %>%
    summarise(
      delta_CO2 = safe_mean(CO2_in - CO2_S),
      delta_CH4 = safe_mean(CH4_in - CH4_S),
      delta_NH3 = safe_mean(NH3_in - NH3_S),
      .groups = "drop"
    ) %>%
    arrange(DATE.HOUR)
}

load_pronova_hourly <- function(project_dir, time_zone = "Europe/Berlin") {
  ndir_path <- file.path(
    project_dir,
    "workflows", "device_comparison", "clean_data", "logas_ndir_campaign", "logas_ndir_campaign_hourly.csv"
  )

  if (!file.exists(ndir_path)) {
    stop("Missing cleaned Pronova campaign file: ", ndir_path)
  }

  read_csv(ndir_path, show_col_types = FALSE) %>%
    mutate(
      DATE.HOUR = parse_datetime_local(DATE.HOUR, time_zone = time_zone),
      across(any_of(c("delta_CO2", "delta_CH4", "delta_NH3")), as.numeric)
    ) %>%
    filter(!is.na(DATE.HOUR)) %>%
    select(DATE.HOUR, delta_CO2, delta_CH4, delta_NH3) %>%
    arrange(DATE.HOUR)
}

build_month_lookup <- function(min_dt, max_dt, time_zone = "Europe/Berlin") {
  if (is.na(min_dt) || is.na(max_dt)) {
    return(tibble(
      period_id   = character(),
      month_label = character(),
      start       = as.POSIXct(character(), tz = time_zone),
      end         = as.POSIXct(character(), tz = time_zone)
    ))
  }

  month_starts <- seq(
    from = floor_date(min_dt, "month"),
    to = floor_date(max_dt, "month"),
    by = "month"
  )

  tibble(
    period_id = tolower(format(month_starts, "%b_%Y")),
    month_label = format(month_starts, "%B %Y"),
    start = month_starts,
    end = (month_starts %m+% months(1)) - hours(1)
  )
}

filter_period <- function(df, start_time, end_time) {
  df %>%
    filter(DATE.HOUR >= start_time, DATE.HOUR <= end_time)
}

build_feasibility_table <- function(ndir_hourly, crds_hourly, month_lookup) {
  bind_rows(lapply(seq_len(nrow(month_lookup)), function(i) {
    period_row <- month_lookup[i, ]
    ndir_month <- filter_period(ndir_hourly, period_row$start, period_row$end)
    crds_month <- filter_period(crds_hourly, period_row$start, period_row$end)
    joined_month <- inner_join(
      ndir_month,
      crds_month,
      by = "DATE.HOUR",
      suffix = c("_ndir", "_crds")
    )

    tibble(
      period_id = period_row$period_id,
      month_label = period_row$month_label,
      ndir_hours = nrow(ndir_month),
      crds_hours = nrow(crds_month),
      overlap_hours = nrow(joined_month),
      feasible = overlap_hours > 0
    )
  }))
}

plot_month_comparison <- function(joined_month, month_label) {
  plot_df <- joined_month %>%
    select(
      DATE.HOUR,
      delta_CO2_ndir, delta_CO2_crds,
      delta_CH4_ndir, delta_CH4_crds,
      delta_NH3_ndir, delta_NH3_crds
    ) %>%
    pivot_longer(
      cols = -DATE.HOUR,
      names_to = c("gas", "device"),
      names_pattern = "(delta_[A-Z0-9]+)_(ndir|crds)",
      values_to = "value"
    ) %>%
    mutate(
      gas = recode(
        gas,
        delta_CO2 = "Delta CO2 (ppm)",
        delta_CH4 = "Delta CH4 (ppm)",
        delta_NH3 = "Delta NH3 (ppm)"
      ),
      device = recode(
        device,
        ndir = "PRONOVA",
        crds = "CRDS"
      )
    )

  ggplot(plot_df, aes(x = DATE.HOUR, y = value, color = device)) +
    geom_line(linewidth = 0.6, na.rm = TRUE) +
    facet_grid(gas ~ ., scales = "free_y", switch = "y") +
    scale_color_manual(values = c("PRONOVA" = "#7570B3", "CRDS" = "#4B5563")) +
    labs(
      title = paste("PRONOVA vs CRDS,", month_label),
      x = NULL,
      y = NULL,
      color = NULL
    ) +
    device_plot_theme_classic()
}

project_dir <- resolve_project_dir(require_crds = TRUE, require_utils = TRUE)
time_zone <- "Europe/Berlin"
workflow_dir <- file.path(project_dir, "workflows", "device_comparison")
clean_dir <- file.path(workflow_dir, "clean_data", "logas_ndir_clean")
table_dir <- file.path(workflow_dir, "result_data", "tables")
plot_dir <- file.path(workflow_dir, "result_data", "plots", "logas_ndir_vs_crds_2026")

dir.create(clean_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

ndir_hourly <- load_pronova_hourly(project_dir, time_zone = time_zone)
crds_hourly <- load_crds_hourly(project_dir, time_zone = time_zone)

analysis_year <- 2026
year_start <- as.POSIXct(sprintf("%d-01-01 00:00:00", analysis_year), tz = time_zone)
year_end <- as.POSIXct(sprintf("%d-12-31 23:00:00", analysis_year), tz = time_zone)

ndir_2026 <- ndir_hourly %>% filter(DATE.HOUR >= year_start, DATE.HOUR <= year_end)
crds_2026 <- crds_hourly %>% filter(DATE.HOUR >= year_start, DATE.HOUR <= year_end)

data_min <- suppressWarnings(min(c(ndir_2026$DATE.HOUR, crds_2026$DATE.HOUR), na.rm = TRUE))
data_max <- suppressWarnings(max(c(ndir_2026$DATE.HOUR, crds_2026$DATE.HOUR), na.rm = TRUE))

if (is.infinite(data_min) || is.infinite(data_max)) {
  stop("No PRONOVA or CRDS data found in ", analysis_year, " - nothing to compare.")
}

month_lookup <- build_month_lookup(data_min, data_max, time_zone = time_zone)

first_id <- month_lookup$period_id[1]
last_id <- month_lookup$period_id[nrow(month_lookup)]
first_month <- sub("_.*", "", first_id)
last_month <- sub("_.*", "", last_id)
year_str <- sub(".*_", "", last_id)
range_id <- if (first_id == last_id) first_id else paste0(first_month, "_to_", last_month, "_", year_str)
range_short <- if (first_id == last_id) first_month else paste0(first_month, "_to_", last_month)

ndir_clean_name <- paste0("logas_ndir_hourly_", range_id, ".csv")
feasibility_name <- paste0("logas_ndir_vs_crds_feasibility_", year_str, "_", range_short, ".csv")

write_csv(ndir_2026, file.path(clean_dir, ndir_clean_name))

feasibility_table <- build_feasibility_table(ndir_2026, crds_2026, month_lookup)
write_csv(feasibility_table, file.path(table_dir, feasibility_name))

for (i in seq_len(nrow(month_lookup))) {
  period_row <- month_lookup[i, ]
  ndir_month <- filter_period(ndir_2026, period_row$start, period_row$end)
  crds_month <- filter_period(crds_2026, period_row$start, period_row$end)
  joined_month <- inner_join(
    ndir_month,
    crds_month,
    by = "DATE.HOUR",
    suffix = c("_ndir", "_crds")
  ) %>%
    arrange(DATE.HOUR)

  write_csv(
    joined_month,
    file.path(table_dir, paste0("logas_ndir_vs_crds_", period_row$period_id, ".csv"))
  )

  if (nrow(joined_month) == 0) {
    next
  }

  month_plot <- plot_month_comparison(joined_month, period_row$month_label)

  ggsave(
    filename = file.path(plot_dir, paste0("logas_ndir_vs_crds_", period_row$period_id, ".png")),
    plot = month_plot,
    width = 11,
    height = 7,
    dpi = 150
  )
}

cat("Read PRONOVA campaign hourly data from:\n", normalizePath(file.path(project_dir, "workflows", "device_comparison", "clean_data", "logas_ndir_campaign", "logas_ndir_campaign_hourly.csv"), winslash = "/", mustWork = TRUE), "\n")
cat("Read CRDS hourly data from:\n", normalizePath(file.path(project_dir, "workflows", "crds_routine_cleaning", "clean_data", "crds_clean", "hourly_in_out_avg"), winslash = "/", mustWork = TRUE), "\n")
cat("Wrote feasibility table to:\n", normalizePath(file.path(table_dir, feasibility_name), winslash = "/", mustWork = FALSE), "\n")
print(feasibility_table)
