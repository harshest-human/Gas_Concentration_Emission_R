library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(readr)
library(tibble)

resolve_project_dir <- function() {
  project_marker <- function(path) {
    dir.exists(file.path(path, "workflows", "device_comparison")) &&
      dir.exists(file.path(path, "workflows", "crds_routine_cleaning")) &&
      dir.exists(file.path(path, "scripts", "utils"))
  }

  candidate_paths <- character()
  args_all <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- sub(file_arg, "", args_all[grepl(file_arg, args_all)])

  if (length(script_path) > 0) {
    candidate_paths <- c(
      candidate_paths,
      normalizePath(file.path(dirname(script_path[1]), "..", "..", ".."), winslash = "/", mustWork = FALSE)
    )
  }

  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    active_path <- tryCatch(
      rstudioapi::getActiveDocumentContext()$path,
      error = function(e) ""
    )

    if (nzchar(active_path)) {
      candidate_paths <- c(
        candidate_paths,
        normalizePath(file.path(dirname(active_path), "..", "..", ".."), winslash = "/", mustWork = FALSE)
      )
    }
  }

  candidate_paths <- c(
    candidate_paths,
    normalizePath(getwd(), winslash = "/", mustWork = FALSE),
    normalizePath(file.path(getwd(), "Gas_Concentration_Emission_R"), winslash = "/", mustWork = FALSE)
  )

  candidate_paths <- unique(candidate_paths[nzchar(candidate_paths)])

  for (path in candidate_paths) {
    if (project_marker(path)) {
      return(path)
    }
  }

  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

parse_crds_datetime <- function(x, time_zone = "Europe/Berlin") {
  x_chr <- as.character(x)
  parsed <- as.POSIXct(rep(NA_character_, length(x_chr)), tz = time_zone)

  formats <- c(
    "%Y/%m/%d %H:%M:%S",
    "%Y-%m-%d %H:%M:%S",
    "%Y-%m-%dT%H:%M:%SZ",
    "%Y-%m-%dT%H:%M:%S"
  )

  for (fmt in formats) {
    missing_idx <- which(is.na(parsed) & !is.na(x_chr))
    if (length(missing_idx) == 0) {
      break
    }

    parsed_try <- as.POSIXct(x_chr[missing_idx], format = fmt, tz = time_zone)
    parsed[missing_idx[!is.na(parsed_try)]] <- parsed_try[!is.na(parsed_try)]
  }

  parsed
}

parse_ndir_datetime <- function(x, time_zone = "Europe/Berlin") {
  parse_date_time(
    as.character(x),
    orders = c("dmy HMS", "dmy HM"),
    tz = time_zone,
    quiet = TRUE
  )
}

safe_mean <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }

  mean(x, na.rm = TRUE)
}

read_log_file <- function(file_path) {
  tryCatch(
    read.table(
      file_path,
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

read_crds_hourly <- function(crds_hourly_dir, time_zone = "Europe/Berlin") {
  crds_files <- list.files(
    path = crds_hourly_dir,
    recursive = TRUE,
    pattern = "\\.csv$",
    full.names = TRUE
  )

  if (length(crds_files) == 0) {
    stop("No CRDS hourly files found in: ", crds_hourly_dir)
  }

  bind_rows(lapply(crds_files, function(file_path) {
    read_csv(file_path, col_types = cols(.default = col_character()), show_col_types = FALSE)
  })) %>%
    mutate(
      DATE.HOUR = parse_crds_datetime(DATE.HOUR, time_zone = time_zone),
      CO2_in = as.numeric(CO2_in),
      CO2_S = as.numeric(CO2_S),
      CH4_in = as.numeric(CH4_in),
      CH4_S = as.numeric(CH4_S),
      NH3_in = as.numeric(NH3_in),
      NH3_S = as.numeric(NH3_S)
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

read_ndir_hourly <- function(ndir_raw_dir, time_zone = "Europe/Berlin") {
  ndir_files <- list.files(
    path = ndir_raw_dir,
    pattern = "^differenzmessung_.*\\.txt$",
    full.names = TRUE
  )

  if (length(ndir_files) == 0) {
    stop("No Pronova LoGAS raw files found in: ", ndir_raw_dir)
  }

  lapply(ndir_files, read_log_file) %>%
    bind_rows() %>%
    mutate(
      DATE.TIME = parse_ndir_datetime(`Datum Uhrzeit`, time_zone = time_zone),
      DATE.HOUR = floor_date(DATE.TIME, "hour")
    ) %>%
    filter(!is.na(DATE.HOUR)) %>%
    group_by(DATE.HOUR) %>%
    summarise(
      delta_CO2 = safe_mean(as.numeric(gsub(",", ".", `CO2 in ppm`))),
      delta_CH4 = safe_mean(as.numeric(gsub(",", ".", `CH4 in ppm`))),
      delta_NH3 = safe_mean(as.numeric(gsub(",", ".", `NH3 in ppm`))),
      .groups = "drop"
    ) %>%
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
    to   = floor_date(max_dt, "month"),
    by   = "month"
  )
  month_ends <- ceiling_date(month_starts, "month") - hours(1)

  tibble(
    period_id   = tolower(format(month_starts, "%b_%Y")),
    month_label = format(month_starts, "%B %Y"),
    start       = month_starts,
    end         = month_ends
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
  label_map <- c(
    delta_CO2_ndir = "Delta CO2 (ppm)",
    delta_CH4_ndir = "Delta CH4 (ppm)",
    delta_NH3_ndir = "Delta NH3 (ppm)"
  )

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
        ndir = "Pronova LoGAS",
        crds = "CRDS"
      )
    )

  ggplot(plot_df, aes(x = DATE.HOUR, y = value, color = device)) +
    geom_line(linewidth = 0.6, na.rm = TRUE) +
    facet_grid(gas ~ ., scales = "free_y", switch = "y") +
    scale_color_manual(values = c("Pronova LoGAS" = "#7570B3", "CRDS" = "#4B5563")) +
    labs(
      title = paste("Pronova LoGAS vs CRDS,", month_label),
      x = NULL,
      y = NULL,
      color = NULL
    ) +
    theme_classic(base_size = 15) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", hjust = 0.5, size = 17),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      strip.text.y.left = element_text(size = 13),
      panel.border = element_rect(color = "black", fill = NA)
    )
}

project_dir <- resolve_project_dir()
time_zone <- "Europe/Berlin"

ndir_raw_dir <- "D:/Data_Analysis_R/owncloud_sync_data/Pronova_LoGAS_raw/Messdaten"
crds_hourly_dir <- file.path(
  project_dir,
  "workflows", "crds_routine_cleaning", "clean_data", "crds_clean", "hourly_in_out_avg", "2026"
)

workflow_dir <- file.path(project_dir, "workflows", "device_comparison")
clean_dir <- file.path(workflow_dir, "clean_data", "logas_ndir_clean")
table_dir <- file.path(workflow_dir, "result_data", "tables")
plot_dir <- file.path(workflow_dir, "result_data", "plots", "logas_ndir_vs_crds_2026")

dir.create(clean_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

ndir_hourly <- read_ndir_hourly(ndir_raw_dir, time_zone = time_zone)
crds_hourly <- read_crds_hourly(crds_hourly_dir, time_zone = time_zone)

# Restrict to the analysis year (consistent with the script's title), then
# derive the month lookup from whatever data is actually present, so new
# months are picked up automatically as raw data accumulates.
analysis_year <- 2026
year_start <- as.POSIXct(sprintf("%d-01-01 00:00:00", analysis_year), tz = time_zone)
year_end   <- as.POSIXct(sprintf("%d-12-31 23:00:00", analysis_year), tz = time_zone)

ndir_2026 <- ndir_hourly %>% filter(DATE.HOUR >= year_start, DATE.HOUR <= year_end)
crds_2026 <- crds_hourly %>% filter(DATE.HOUR >= year_start, DATE.HOUR <= year_end)

data_min <- suppressWarnings(min(c(ndir_2026$DATE.HOUR, crds_2026$DATE.HOUR), na.rm = TRUE))
data_max <- suppressWarnings(max(c(ndir_2026$DATE.HOUR, crds_2026$DATE.HOUR), na.rm = TRUE))

if (is.infinite(data_min) || is.infinite(data_max)) {
  stop("No NDIR or CRDS data found in ", analysis_year, " - nothing to compare.")
}

month_lookup <- build_month_lookup(data_min, data_max, time_zone = time_zone)

# Dynamic output names: 'jan_to_may_2026' when the range spans multiple months,
# 'may_2026' when only a single month is present.
first_id    <- month_lookup$period_id[1]
last_id     <- month_lookup$period_id[nrow(month_lookup)]
first_month <- sub("_.*", "", first_id)
last_month  <- sub("_.*", "", last_id)
year_str    <- sub(".*_", "", last_id)
range_id    <- if (first_id == last_id) first_id else paste0(first_month, "_to_", last_month, "_", year_str)
range_short <- if (first_id == last_id) first_month else paste0(first_month, "_to_", last_month)

ndir_clean_name  <- paste0("logas_ndir_hourly_", range_id, ".csv")
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

cat("Read Pronova LoGAS raw data from:\n", normalizePath(ndir_raw_dir, winslash = "/", mustWork = TRUE), "\n")
cat("Read CRDS hourly data from:\n", normalizePath(crds_hourly_dir, winslash = "/", mustWork = TRUE), "\n")
cat("Wrote feasibility table to:\n", normalizePath(file.path(table_dir, feasibility_name), winslash = "/", mustWork = FALSE), "\n")
print(feasibility_table)
