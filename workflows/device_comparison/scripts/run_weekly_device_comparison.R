suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(lubridate)
  library(purrr)
  library(readr)
  library(readxl)
  library(tibble)
  library(tidyr)
})

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args_all[grepl("--file=", args_all)])
script_dir <- if (length(script_path)) {
  dirname(normalizePath(script_path[1], winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

workflow_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
project_dir <- normalizePath(file.path(workflow_dir, "..", ".."), winslash = "/", mustWork = TRUE)

time_zone <- "Europe/Berlin"
start_time <- ymd_hms("2025-12-09 00:00:00", tz = time_zone)
end_time <- ymd_hms("2026-10-08 00:00:00", tz = time_zone) # exclusive; includes 07 October

paths <- list(
  cubic = "D:/Data_Analysis_R/owncloud_sync_data/Cubic_raw",
  pronova = "D:/Data_Analysis_R/owncloud_sync_data/Pronova_LoGAS_raw/Messdaten",
  crds = file.path(project_dir, "workflows", "crds_routine_cleaning", "clean_data", "crds_clean", "hourly_in_out_avg"),
  hourly = file.path(workflow_dir, "result_data", "tables", "hourly"),
  daily = file.path(workflow_dir, "result_data", "tables", "daily"),
  weekly = file.path(workflow_dir, "result_data", "tables", "weekly"),
  qc = file.path(workflow_dir, "result_data", "tables", "qc"),
  plots = file.path(workflow_dir, "result_data", "plots", "weekly")
)

walk(c(paths$hourly, paths$daily, paths$weekly, paths$qc), dir.create, recursive = TRUE, showWarnings = FALSE)
walk(c("mean_sd", "mean_95ci", "median_iqr"), ~ dir.create(file.path(paths$plots, .x), recursive = TRUE, showWarnings = FALSE))

device_levels <- c("CRDS08", "PRONOVA", "CUBIC")
device_colors <- c(CRDS08 = "#4B5563", PRONOVA = "#7570B3", CUBIC = "#1B9E77")
gas_levels <- c("delta_co2", "delta_ch4", "delta_nh3")
gas_labels <- c(
  delta_co2 = "Delta CO2 (ppm)",
  delta_ch4 = "Delta CH4 (ppm)",
  delta_nh3 = "Delta NH3 (ppm)"
)

safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
safe_sd <- function(x) if (sum(is.finite(x)) < 2) NA_real_ else sd(x, na.rm = TRUE)

parse_local_datetime <- function(x) {
  parse_date_time(
    as.character(x),
    orders = c("ymd HMS", "ymd HM", "Ymd HMS", "Ymd HM", "Y/m/d HMS", "dmy HMS", "dmy HM"),
    tz = time_zone,
    quiet = TRUE
  )
}

collapse_sources <- function(x) paste(sort(unique(x[!is.na(x) & nzchar(x)])), collapse = ";")

# -----------------------------------------------------------------------------
# CUBIC: raw Excel -> hourly inside/outside -> hourly delta
# -----------------------------------------------------------------------------

read_cubic_hourly <- function() {
  files <- list.files(paths$cubic, pattern = "\\.xlsx$", full.names = TRUE)
  parsing_log <- list()

  raw <- map_dfr(files, function(path) {
    tryCatch(
      read_excel(path) |>
        mutate(source_file = basename(path)),
      error = function(e) {
        parsing_log[[length(parsing_log) + 1L]] <<- tibble(
          analyser = "CUBIC", source_file = basename(path), status = "failed", message = conditionMessage(e)
        )
        NULL
      }
    )
  })

  required <- c("Time", "Type", "CO2", "CH4", "NH3", "source_file")
  missing <- setdiff(required, names(raw))
  if (length(missing)) stop("CUBIC input is missing: ", paste(missing, collapse = ", "), call. = FALSE)

  hourly_long <- raw |>
    transmute(
      date_time = floor_date(parse_local_datetime(Time), "hour"),
      stream = case_when(
        suppressWarnings(as.numeric(Type)) == 1 ~ "in",
        !is.na(Type) ~ "out",
        TRUE ~ NA_character_
      ),
      CO2 = suppressWarnings(as.numeric(CO2)),
      CH4 = suppressWarnings(as.numeric(CH4)),
      NH3 = suppressWarnings(as.numeric(NH3)),
      source_file
    ) |>
    filter(date_time >= start_time, date_time < end_time, !is.na(stream)) |>
    pivot_longer(c(CO2, CH4, NH3), names_to = "gas", values_to = "value") |>
    group_by(date_time, stream, gas) |>
    summarise(
      value = safe_mean(value),
      n_raw = sum(is.finite(value)),
      source_files = collapse_sources(source_file),
      .groups = "drop"
    )

  wide <- hourly_long |>
    select(date_time, stream, gas, value) |>
    pivot_wider(names_from = c(gas, stream), values_from = value, names_glue = "{tolower(gas)}_{stream}")

  sources <- hourly_long |>
    group_by(date_time) |>
    summarise(source_files = collapse_sources(source_files), .groups = "drop")

  hourly <- wide |>
    left_join(sources, by = "date_time") |>
    mutate(
      analyser = "CUBIC",
      delta_co2 = co2_in - co2_out,
      delta_ch4 = ch4_in - ch4_out,
      delta_nh3 = nh3_in - nh3_out
    ) |>
    select(date_time, analyser, co2_in, co2_out, ch4_in, ch4_out, nh3_in, nh3_out,
           delta_co2, delta_ch4, delta_nh3, source_files) |>
    arrange(date_time)

  list(
    data = hourly,
    log = bind_rows(
      tibble(analyser = "CUBIC", source_file = basename(files), status = "read", message = NA_character_),
      bind_rows(parsing_log)
    )
  )
}

# -----------------------------------------------------------------------------
# PRONOVA: differential raw TXT -> hourly delta
# -----------------------------------------------------------------------------

read_pronova_file <- function(path) {
  read.table(
    path,
    header = TRUE,
    sep = "\t",
    dec = ",",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    fill = TRUE,
    comment.char = "",
    colClasses = "character"
  ) |>
    as_tibble() |>
    mutate(source_file = basename(path))
}

read_pronova_hourly <- function() {
  files <- list.files(paths$pronova, pattern = "^differenzmessung_.*\\.txt$", full.names = TRUE)
  parsing_log <- list()
  raw <- map_dfr(files, function(path) {
    tryCatch(
      read_pronova_file(path),
      error = function(e) {
        parsing_log[[length(parsing_log) + 1L]] <<- tibble(
          analyser = "PRONOVA", source_file = basename(path), status = "failed", message = conditionMessage(e)
        )
        NULL
      }
    )
  })

  required <- c("Datum Uhrzeit", "CO2 in ppm", "CH4 in ppm", "NH3 in ppm", "source_file")
  missing <- setdiff(required, names(raw))
  if (length(missing)) stop("PRONOVA input is missing: ", paste(missing, collapse = ", "), call. = FALSE)

  parse_comma <- function(x) parse_double(x, locale = locale(decimal_mark = ","))

  hourly <- raw |>
    transmute(
      date_time = floor_date(parse_date_time(`Datum Uhrzeit`, orders = c("dmy HMS", "dmy HM"), tz = time_zone, quiet = TRUE), "hour"),
      delta_co2 = parse_comma(`CO2 in ppm`),
      delta_ch4 = parse_comma(`CH4 in ppm`),
      delta_nh3 = parse_comma(`NH3 in ppm`),
      source_file
    ) |>
    filter(date_time >= start_time, date_time < end_time) |>
    group_by(date_time) |>
    summarise(
      delta_co2 = safe_mean(delta_co2),
      delta_ch4 = safe_mean(delta_ch4),
      delta_nh3 = safe_mean(delta_nh3),
      source_files = collapse_sources(source_file),
      .groups = "drop"
    ) |>
    mutate(
      analyser = "PRONOVA",
      co2_in = NA_real_, co2_out = NA_real_,
      ch4_in = NA_real_, ch4_out = NA_real_,
      nh3_in = NA_real_, nh3_out = NA_real_
    ) |>
    select(date_time, analyser, co2_in, co2_out, ch4_in, ch4_out, nh3_in, nh3_out,
           delta_co2, delta_ch4, delta_nh3, source_files) |>
    arrange(date_time)

  list(
    data = hourly,
    log = bind_rows(
      tibble(analyser = "PRONOVA", source_file = basename(files), status = "read", message = NA_character_),
      bind_rows(parsing_log)
    )
  )
}

# -----------------------------------------------------------------------------
# CRDS08: cleaned hourly routine output -> explicitly filtered CRDS08 reference
# -----------------------------------------------------------------------------

first_numeric_column <- function(data, candidates) {
  present <- intersect(candidates, names(data))
  if (!length(present)) return(rep(NA_real_, nrow(data)))
  result <- rep(NA_real_, nrow(data))
  for (column in present) {
    candidate <- suppressWarnings(as.numeric(data[[column]]))
    replace <- is.na(result) & !is.na(candidate)
    result[replace] <- candidate[replace]
  }
  result
}

read_crds08_hourly <- function() {
  files <- list.files(paths$crds, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
  files <- files[grepl("CRDS", basename(files)) & !grepl("FTIR|LUFA|UB|ANECO|MBBM", basename(files))]
  parsing_log <- list()

  standardized <- map_dfr(files, function(path) {
    tryCatch({
      raw <- read_csv(path, show_col_types = FALSE, col_types = cols(.default = col_character()))
      datetime_column <- intersect(c("DATE.HOUR", "DATE.TIME"), names(raw))[1]
      if (is.na(datetime_column)) stop("No DATE.HOUR or DATE.TIME column")

      date_time <- parse_local_datetime(raw[[datetime_column]])
      analyser_raw <- if ("analyzer" %in% names(raw)) raw$analyzer else NA_character_
      analyser_normalized <- toupper(gsub("[^A-Z0-9]", "", analyser_raw))
      filename_is_crds08 <- grepl("CRDS8", basename(path)) && !grepl("CRDS9_", basename(path))
      is_crds08 <- analyser_normalized == "CRDS8"
      is_crds08[is.na(is_crds08)] <- filename_is_crds08

      tibble(
        date_time = date_time,
        is_crds08 = is_crds08,
        co2_in = first_numeric_column(raw, c("CO2_ppm_in", "CO2_in")),
        co2_out = first_numeric_column(raw, c("CO2_ppm_S", "CO2_S")),
        ch4_in = first_numeric_column(raw, c("CH4_ppm_in", "CH4_in")),
        ch4_out = first_numeric_column(raw, c("CH4_ppm_S", "CH4_S")),
        nh3_in = first_numeric_column(raw, c("NH3_ppm_in", "NH3_in")),
        nh3_out = first_numeric_column(raw, c("NH3_ppm_S", "NH3_S")),
        source_file = basename(path)
      )
    }, error = function(e) {
      parsing_log[[length(parsing_log) + 1L]] <<- tibble(
        analyser = "CRDS08", source_file = basename(path), status = "failed", message = conditionMessage(e)
      )
      NULL
    })
  })

  hourly <- standardized |>
    filter(is_crds08, date_time >= start_time, date_time < end_time) |>
    group_by(date_time) |>
    summarise(
      across(c(co2_in, co2_out, ch4_in, ch4_out, nh3_in, nh3_out), safe_mean),
      source_files = collapse_sources(source_file),
      .groups = "drop"
    ) |>
    mutate(
      analyser = "CRDS08",
      delta_co2 = co2_in - co2_out,
      delta_ch4 = ch4_in - ch4_out,
      delta_nh3 = nh3_in - nh3_out
    ) |>
    select(date_time, analyser, co2_in, co2_out, ch4_in, ch4_out, nh3_in, nh3_out,
           delta_co2, delta_ch4, delta_nh3, source_files) |>
    arrange(date_time)

  selected_files <- unique(unlist(strsplit(hourly$source_files, ";", fixed = TRUE)))
  list(
    data = hourly,
    log = bind_rows(
      tibble(
        analyser = "CRDS08",
        source_file = basename(files),
        status = if_else(basename(files) %in% selected_files, "read_crds08", "excluded_non_crds08_or_outside_period"),
        message = NA_character_
      ),
      bind_rows(parsing_log)
    )
  )
}

# -----------------------------------------------------------------------------
# Combine the hourly analyzer rows and calculate daily statistics
# -----------------------------------------------------------------------------

message("Reading and cleaning CUBIC raw data...")
cubic <- read_cubic_hourly()
message("Reading and cleaning PRONOVA raw data...")
pronova <- read_pronova_hourly()
message("Reading cleaned CRDS08 reference data...")
crds08 <- read_crds08_hourly()

hourly <- bind_rows(crds08$data, pronova$data, cubic$data) |>
  mutate(
    analyser = factor(analyser, levels = device_levels),
    date = as.Date(date_time, tz = time_zone),
    iso_year = isoyear(date_time),
    iso_week = isoweek(date_time),
    week_start = as.Date(floor_date(date_time, unit = "week", week_start = 1), tz = time_zone),
    week_end = week_start + days(6)
  ) |>
  arrange(date_time, analyser)

duplicate_keys <- hourly |>
  count(date_time, analyser) |>
  filter(n > 1)
if (nrow(duplicate_keys)) stop("Duplicate date_time/analyser rows remain after cleaning.", call. = FALSE)

hourly_path <- file.path(paths$hourly, "device_hourly_20251209_20261007.csv")
write_csv(hourly, hourly_path, na = "NA")

hourly_long <- hourly |>
  select(date_time, date, iso_year, iso_week, week_start, week_end, analyser, all_of(gas_levels)) |>
  pivot_longer(all_of(gas_levels), names_to = "gas", values_to = "value")

daily <- hourly_long |>
  group_by(date, iso_year, iso_week, week_start, week_end, analyser, gas) |>
  summarise(
    n_hours = sum(is.finite(value)),
    coverage_percent = 100 * n_hours / 24,
    mean = if_else(n_hours > 0, mean(value, na.rm = TRUE), NA_real_),
    sd = safe_sd(value),
    se = if_else(n_hours > 1, sd / sqrt(n_hours), NA_real_),
    ci_multiplier = if (n_hours > 1) qt(0.975, df = n_hours - 1) else NA_real_,
    ci95_lower = mean - ci_multiplier * se,
    ci95_upper = mean + ci_multiplier * se,
    median = if_else(n_hours > 0, median(value, na.rm = TRUE), NA_real_),
    q1 = if_else(n_hours > 0, as.numeric(quantile(value, 0.25, na.rm = TRUE, names = FALSE)), NA_real_),
    q3 = if_else(n_hours > 0, as.numeric(quantile(value, 0.75, na.rm = TRUE, names = FALSE)), NA_real_),
    iqr = q3 - q1,
    minimum = if (n_hours > 0) min(value, na.rm = TRUE) else NA_real_,
    maximum = if (n_hours > 0) max(value, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  ) |>
  mutate(
    analyser = factor(analyser, levels = device_levels),
    gas = factor(gas, levels = gas_levels, labels = unname(gas_labels[gas_levels])),
    quality_flag = case_when(
      n_hours == 0 ~ "missing",
      n_hours < 18 ~ "partial_lt18h",
      n_hours < 24 ~ "partial_18_to_23h",
      TRUE ~ "complete_24h"
    )
  ) |>
  arrange(date, analyser, gas)

daily_path <- file.path(paths$daily, "device_daily_statistics_20251209_20261007.csv")
write_csv(daily, daily_path, na = "NA")

coverage <- daily |>
  group_by(iso_year, iso_week, week_start, week_end, analyser, gas) |>
  summarise(
    days_with_data = sum(n_hours > 0),
    complete_days = sum(n_hours == 24),
    partial_days = sum(n_hours > 0 & n_hours < 24),
    valid_hourly_values = sum(n_hours),
    .groups = "drop"
  )

coverage_path <- file.path(paths$weekly, "weekly_coverage_summary.csv")
write_csv(coverage, coverage_path, na = "NA")
write_csv(bind_rows(crds08$log, pronova$log, cubic$log), file.path(paths$qc, "source_file_log.csv"), na = "NA")
write_csv(daily |> select(date, analyser, gas, n_hours, coverage_percent, quality_flag), file.path(paths$qc, "daily_coverage.csv"), na = "NA")
write_csv(duplicate_keys, file.path(paths$qc, "duplicate_hourly_keys.csv"), na = "NA")

# -----------------------------------------------------------------------------
# Weekly calendar plots
# -----------------------------------------------------------------------------

add_segments <- function(plot_data, value_column) {
  plot_data |>
    arrange(gas, analyser, date) |>
    group_by(gas, analyser) |>
    mutate(
      plot_value = .data[[value_column]],
      gap = is.na(plot_value) | is.na(lag(plot_value)) | is.na(lag(date)) | as.integer(date - lag(date)) != 1L,
      segment = cumsum(replace_na(gap, TRUE))
    ) |>
    ungroup()
}

weekly_plot <- function(week_data, statistic, week_start_value, week_end_value) {
  specifications <- list(
    mean_sd = list(value = "mean", lower = "mean_minus_sd", upper = "mean_plus_sd",
                   subtitle = "Daily line = arithmetic mean; shaded band = mean +/- within-day SD"),
    mean_95ci = list(value = "mean", lower = "ci95_lower", upper = "ci95_upper",
                     subtitle = "Daily line = arithmetic mean; shaded band = 95% confidence interval of the mean")
  )
  spec <- specifications[[statistic]]

  display_start <- max(as.Date(start_time, tz = time_zone), week_start_value)
  display_end <- min(as.Date(end_time - seconds(1), tz = time_zone), week_end_value)

  complete_grid <- crossing(
    date = seq(display_start, display_end, by = "day"),
    analyser = factor(device_levels, levels = device_levels),
    gas = factor(unname(gas_labels[gas_levels]), levels = unname(gas_labels[gas_levels]))
  ) |>
    left_join(week_data, by = c("date", "analyser", "gas")) |>
    mutate(
      mean_minus_sd = mean - sd,
      mean_plus_sd = mean + sd,
      lower = pmax(0, .data[[spec$lower]]),
      upper = .data[[spec$upper]]
    ) |>
    add_segments(spec$value)

  ggplot(
    complete_grid,
    aes(
      x = date,
      y = .data[[spec$value]],
      colour = analyser,
      fill = analyser,
      group = interaction(analyser, segment)
    )
  ) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, colour = NA, na.rm = TRUE) +
    geom_line(linewidth = 0.9, na.rm = TRUE) +
    geom_point(size = 1.8, na.rm = TRUE) +
    facet_grid(gas ~ ., scales = "free_y", switch = "y") +
    scale_colour_manual(values = device_colors, drop = FALSE) +
    scale_fill_manual(values = device_colors, drop = FALSE) +
    scale_x_date(
      limits = c(display_start, display_end),
      breaks = seq(display_start, display_end, by = "day"),
      date_labels = "%a\n%d-%m",
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(
      limits = c(0, NA),
      breaks = scales::pretty_breaks(n = 5),
      expand = expansion(mult = c(0, 0.06))
    ) +
    labs(
      title = paste0(
        "Device comparison - ISO week ", isoyear(week_start_value), "-W",
        sprintf("%02d", isoweek(week_start_value))
      ),
      subtitle = paste0(
        format(display_start, "%d-%m-%Y"), " to ", format(display_end, "%d-%m-%Y"),
        " | ", spec$subtitle
      ),
      caption = "All available analyzers are retained; missing analyzer-days are left blank. Daily n and coverage are stored in the result tables.",
      x = NULL,
      y = NULL,
      colour = NULL,
      fill = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      plot.caption = element_text(hjust = 0, size = 8, colour = "grey35"),
      axis.text.x = element_text(size = 9),
      strip.placement = "outside",
      strip.background = element_rect(fill = "white", colour = "black"),
      strip.text.y.left = element_text(angle = 90, size = 10),
      panel.border = element_rect(colour = "black", fill = NA),
      legend.position = "bottom"
    ) +
    guides(fill = "none", colour = guide_legend(nrow = 1))
}

weekly_boxplot <- function(hourly_week_data, week_start_value, week_end_value) {
  display_start <- max(as.Date(start_time, tz = time_zone), week_start_value)
  display_end <- min(as.Date(end_time - seconds(1), tz = time_zone), week_end_value)

  plot_data <- hourly_week_data |>
    filter(date >= display_start, date <= display_end, is.finite(value)) |>
    mutate(
      analyser = factor(analyser, levels = device_levels),
      gas = factor(gas, levels = gas_levels, labels = unname(gas_labels[gas_levels])),
      date_label = factor(
        format(date, "%a\n%d-%m"),
        levels = format(seq(display_start, display_end, by = "day"), "%a\n%d-%m")
      )
    )

  ggplot(plot_data, aes(x = date_label, y = value, colour = analyser, fill = analyser)) +
    geom_boxplot(
      aes(group = interaction(date_label, analyser)),
      position = position_dodge2(width = 0.82, preserve = "single", padding = 0.12),
      width = 0.72,
      linewidth = 0.55,
      alpha = 0.18,
      outlier.alpha = 0.45,
      outlier.size = 1.1,
      na.rm = TRUE
    ) +
    facet_grid(gas ~ ., scales = "free_y", switch = "y") +
    scale_colour_manual(values = device_colors, drop = FALSE) +
    scale_fill_manual(values = device_colors, drop = FALSE) +
    scale_y_continuous(
      limits = c(0, NA),
      breaks = scales::pretty_breaks(n = 5),
      expand = expansion(mult = c(0, 0.06))
    ) +
    labs(
      title = paste0(
        "Device comparison - ISO week ", isoyear(week_start_value), "-W",
        sprintf("%02d", isoweek(week_start_value))
      ),
      subtitle = paste0(
        format(display_start, "%d-%m-%Y"), " to ", format(display_end, "%d-%m-%Y"),
        " | Daily boxplots: median, IQR, 1.5 x IQR whiskers and hourly outliers"
      ),
      caption = "Boxes use all available hourly values. Missing analyzer-days are left blank; daily n and coverage are stored in the result tables.",
      x = NULL,
      y = NULL,
      colour = NULL,
      fill = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      plot.caption = element_text(hjust = 0, size = 8, colour = "grey35"),
      axis.text.x = element_text(size = 9),
      strip.placement = "outside",
      strip.background = element_rect(fill = "white", colour = "black"),
      strip.text.y.left = element_text(angle = 90, size = 10),
      panel.border = element_rect(colour = "black", fill = NA),
      legend.position = "bottom"
    ) +
    guides(colour = guide_legend(nrow = 1), fill = "none")
}

weeks <- daily |>
  distinct(iso_year, iso_week, week_start, week_end) |>
  arrange(week_start)

plot_manifest <- list()
manifest_index <- 1L

for (i in seq_len(nrow(weeks))) {
  week <- weeks[i, ]
  week_data <- daily |>
    filter(iso_year == week$iso_year, iso_week == week$iso_week)
  hourly_week_data <- hourly_long |>
    filter(iso_year == week$iso_year, iso_week == week$iso_week)
  display_start <- max(as.Date(start_time, tz = time_zone), week$week_start)
  display_end <- min(as.Date(end_time - seconds(1), tz = time_zone), week$week_end)
  base_tag <- paste0(
    week$iso_year, "-W", sprintf("%02d", week$iso_week), "_",
    format(display_start, "%Y-%m-%d"), "_", format(display_end, "%Y-%m-%d")
  )

  for (statistic in c("mean_sd", "mean_95ci", "median_iqr")) {
    p <- if (statistic == "median_iqr") {
      weekly_boxplot(hourly_week_data, week$week_start, week$week_end)
    } else {
      weekly_plot(week_data, statistic, week$week_start, week$week_end)
    }
    filename <- paste0(base_tag, "_", statistic, ".png")
    output_path <- file.path(paths$plots, statistic, filename)
    ggsave(output_path, p, width = 11, height = 7.2, dpi = 180, bg = "white")
    plot_manifest[[manifest_index]] <- tibble(
      iso_year = week$iso_year,
      iso_week = week$iso_week,
      week_start = week$week_start,
      week_end = week$week_end,
      display_start = display_start,
      display_end = display_end,
      statistic = statistic,
      file = normalizePath(output_path, winslash = "/", mustWork = TRUE)
    )
    manifest_index <- manifest_index + 1L
  }
}

manifest <- bind_rows(plot_manifest)
manifest_path <- file.path(paths$weekly, "weekly_plot_manifest.csv")
write_csv(manifest, manifest_path, na = "NA")

cat("\nCompleted weekly device comparison.\n")
cat("Hourly data:", normalizePath(hourly_path, winslash = "/"), "\n")
cat("Daily statistics:", normalizePath(daily_path, winslash = "/"), "\n")
cat("Coverage:", normalizePath(coverage_path, winslash = "/"), "\n")
cat("Plot manifest:", normalizePath(manifest_path, winslash = "/"), "\n")
cat("Weekly plots:", nrow(manifest), "\n\n")

print(
  hourly |>
    group_by(analyser) |>
    summarise(
      first_hour = min(date_time),
      last_hour = max(date_time),
      hourly_rows = n(),
      .groups = "drop"
    )
)
