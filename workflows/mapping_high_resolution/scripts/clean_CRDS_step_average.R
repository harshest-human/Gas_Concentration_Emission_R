# CRDS8/CRDS9 step-average cleaning for the high-resolution mapping workflow.
#
# For each continuous MPVPosition step:
#   - use a 240-second measurement interval;
#   - discard the first 60 seconds as flushing time;
#   - average gas measurements from elapsed second 60 up to second 240.
#
# This script writes only one step-average CSV. It does not create hourly,
# KTBL, intermediate, or diagnostic CSV files.

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package 'data.table' is required.")
}

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
analyser_name <- if (length(args) >= 1L) toupper(args[1]) else "CRDS8"
start_text <- if (length(args) >= 2L) {
  args[2]
} else {
  "2024-05-16 17:27:00"
}
end_text <- if (length(args) >= 3L) {
  args[3]
} else {
  "2024-10-29 13:29:12"
}
location_map_name <- if (length(args) >= 4L) {
  tolower(args[4])
} else {
  "original"
}
make_plot <- if (length(args) >= 5L) {
  tolower(args[5]) %in% c("true", "yes", "1")
} else {
  TRUE
}

analyser_paths <- c(
  CRDS8 = "D:/Data_Analysis_R/CRDS_raw_recovered/Picarro_G2508/CRDS08_raw/2024",
  CRDS9 = "D:/Data_Analysis_R/CRDS_raw_recovered/Picarro_G2509/CRDS09_raw/2024"
)

if (!analyser_name %in% names(analyser_paths)) {
  stop("Analyser must be CRDS8 or CRDS9.")
}

input_dir <- unname(analyser_paths[analyser_name])
start_time <- as.POSIXct(
  start_text,
  format = "%Y-%m-%d %H:%M:%S",
  tz = "Europe/Berlin"
)
end_time <- as.POSIXct(
  end_text,
  format = "%Y-%m-%d %H:%M:%S",
  tz = "Europe/Berlin"
)

if (is.na(start_time) || is.na(end_time) || start_time > end_time) {
  stop("Invalid start or end time.")
}

output_dir <- paste0(
  "D:/Data_Analysis_R/Gas_Concentration_Emission_R/",
  "workflows/mapping_high_resolution/clean_data"
)
output_file <- file.path(
  output_dir,
  paste0(
    format(start_time, "%Y%m%d.%H%M%S"),
    "_",
    format(end_time, "%Y%m%d.%H%M%S"),
    "_",
    analyser_name,
    "_240s_step_average.csv"
  )
)

flush_seconds <- 60
interval_seconds <- 240
gas_columns <- c("CO2", "CH4", "NH3", "H2O")
required_columns <- c("DATE", "TIME", "MPVPosition", gas_columns)

extract_file_time <- function(path) {
  file_name <- basename(path)
  matched <- regexec("-(\\d{8})-(\\d{6})-", file_name)
  pieces <- regmatches(file_name, matched)[[1]]
  if (length(pieces) != 3L) {
    return(as.POSIXct(NA, tz = "Europe/Berlin"))
  }
  as.POSIXct(
    paste0(pieces[2], pieces[3]),
    format = "%Y%m%d%H%M%S",
    tz = "Europe/Berlin"
  )
}

all_files <- list.files(
  input_dir,
  pattern = "\\.dat$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(all_files) == 0L) {
  stop("No CRDS .dat files found in: ", input_dir)
}

file_times <- as.POSIXct(
  vapply(all_files, function(path) as.numeric(extract_file_time(path)), numeric(1)),
  origin = "1970-01-01",
  tz = "Europe/Berlin"
)

# Include two hours before the requested start so a file that straddles the
# boundary is not missed. Row-level filtering below enforces the exact window.
files_to_read <- all_files[
  !is.na(file_times) &
    file_times >= start_time - 2 * 60 * 60 &
    file_times <= end_time
]
files_to_read <- files_to_read[order(file_times[
  !is.na(file_times) &
    file_times >= start_time - 2 * 60 * 60 &
    file_times <= end_time
])]

if (length(files_to_read) == 0L) {
  stop("No CRDS files overlap the requested time range.")
}

message("Reading ", length(files_to_read), " ", analyser_name, " raw files...")

read_one_file <- function(path, index, total) {
  column_names <- names(fread(path, nrows = 0, showProgress = FALSE))
  missing_columns <- setdiff(required_columns, column_names)
  if (length(missing_columns) > 0L) {
    warning(
      "Skipping ", basename(path), "; missing: ",
      paste(missing_columns, collapse = ", ")
    )
    return(NULL)
  }

  result <- tryCatch(
    fread(
      path,
      select = required_columns,
      showProgress = FALSE
    ),
    error = function(error) {
      warning("Skipping unreadable file ", basename(path), ": ", error$message)
      NULL
    }
  )

  if (index %% 250L == 0L || index == total) {
    message("Read ", index, " of ", total, " files")
  }
  result
}

raw_parts <- lapply(
  seq_along(files_to_read),
  function(index) read_one_file(
    files_to_read[index],
    index,
    length(files_to_read)
  )
)
raw_parts <- Filter(Negate(is.null), raw_parts)

if (length(raw_parts) == 0L) {
  stop("None of the selected CRDS files could be read.")
}

raw <- rbindlist(raw_parts, use.names = TRUE, fill = TRUE)
rm(raw_parts)
gc()

raw[, DATE.TIME := as.POSIXct(
  paste(DATE, TIME),
  format = "%Y-%m-%d %H:%M:%OS",
  tz = "Europe/Berlin"
)]
raw <- raw[
  !is.na(DATE.TIME) &
    DATE.TIME >= start_time &
    DATE.TIME <= end_time
]
raw[, c("DATE", "TIME") := NULL]

# Remove recovered-file overlaps using the timestamp and measurement fields.
setorder(raw, DATE.TIME)
raw <- unique(
  raw,
  by = c("DATE.TIME", "MPVPosition", gas_columns)
)

# Keep valid, non-zero integer valve positions.
raw <- raw[
  !is.na(MPVPosition) &
    MPVPosition != 0 &
    MPVPosition == floor(MPVPosition)
]

# Start a new step when the valve position changes or when a data gap exceeds
# 10 seconds. The gap rule prevents separate runs at the same position from
# being combined.
setorder(raw, DATE.TIME)
raw[, time_gap := c(Inf, diff(as.numeric(DATE.TIME)))]
raw[, step_id := cumsum(
  seq_len(.N) == 1L |
    MPVPosition != shift(MPVPosition, fill = first(MPVPosition)) |
    time_gap > 10
)]
raw[, step_start := min(DATE.TIME), by = step_id]
raw[, elapsed_seconds := as.numeric(
  difftime(DATE.TIME, step_start, units = "secs")
)]

# Subdivide long continuous runs at the same MPVPosition into consecutive
# 240-second intervals. This ensures that a position held for hours or days
# still produces one average per complete 240-second interval.
raw[, interval_id := floor(elapsed_seconds / interval_seconds)]
raw[, interval_elapsed := elapsed_seconds - interval_id * interval_seconds]

# Retain nominally complete 240-second intervals and average the post-flush
# window. Valve transitions and sub-second sampling commonly leave an observed
# first-to-last span of 235-239 seconds, so 230 seconds is used as a conservative
# completeness threshold. Gas values are still restricted to elapsed seconds
# 60 through 240.
complete_intervals <- raw[
  ,
  .(
    interval_end = max(DATE.TIME),
    observed_span = max(interval_elapsed)
  ),
  by = .(step_id, interval_id, MPVPosition)
][observed_span >= interval_seconds - 10]

average_window <- raw[
  complete_intervals,
  on = .(step_id, interval_id, MPVPosition),
  nomatch = 0L
][
  interval_elapsed >= flush_seconds &
    interval_elapsed < interval_seconds
]

step_average <- average_window[
  ,
  .(
    DATE.TIME = min(interval_end),
    analyser = analyser_name,
    location = as.character(first(MPVPosition)),
    CO2 = mean(CO2, na.rm = TRUE),
    CH4 = mean(CH4, na.rm = TRUE),
    NH3 = mean(NH3, na.rm = TRUE) / 1000,
    H2O = mean(H2O, na.rm = TRUE)
  ),
  by = .(step_id, interval_id)
]

step_average <- step_average[
  DATE.TIME >= start_time &
    DATE.TIME <= end_time
]

if (location_map_name == "campaign2") {
  campaign2_locations <- list(
    CRDS9 = c(1, 3, 4, 6, 7, 9, 10, 12, 13, 15, 16, 18, 21, 22, 24, 25),
    CRDS8 = c(27, 28, 30, 31, 33, 34, 36, 37, 39, 42, 43, 45, 46, 48, 49, 51)
  )
  source_position <- as.integer(step_average[["location"]])
  valid_position <- !is.na(source_position) &
    source_position %in% seq_along(campaign2_locations[[analyser_name]])

  if (any(!valid_position)) {
    warning(
      "Discarding ", sum(!valid_position),
      " interval(s) outside MPVPositions 1:16"
    )
  }

  step_average <- step_average[valid_position]
  step_average[, location := as.character(
    campaign2_locations[[analyser_name]][source_position[valid_position]]
  )]
} else if (location_map_name != "original") {
  stop("Location map must be 'original' or 'campaign2'.")
}

setorder(step_average, DATE.TIME)
step_average[, c("step_id", "interval_id") := NULL]
setcolorder(
  step_average,
  c("DATE.TIME", "analyser", "location", "CO2", "CH4", "NH3", "H2O")
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(
  step_average,
  output_file,
  quote = FALSE,
  dateTimeAs = "write.csv"
)

message(
  "Wrote ", nrow(step_average), " ", analyser_name,
  " step averages to: ", output_file
)

# Plot the CRDS measurement location through time when requested.
if (isTRUE(make_plot)) {
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required to create the CRDS location plot.")
}

plot_dir <- paste0(
  "D:/Data_Analysis_R/Gas_Concentration_Emission_R/",
  "workflows/mapping_high_resolution/plots"
)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

location_breaks <- sort(unique(as.numeric(step_average[["location"]])))
location_plot <- ggplot2::ggplot(
  step_average,
  ggplot2::aes(
    x = DATE.TIME,
    y = as.numeric(location),
    colour = analyser
  )
) +
  ggplot2::geom_point(alpha = 0.65, size = 0.75) +
  ggplot2::scale_colour_manual(
    values = stats::setNames("#0072B2", analyser_name)
  ) +
  ggplot2::scale_y_continuous(
    breaks = location_breaks,
    labels = location_breaks,
    minor_breaks = NULL
  ) +
  ggplot2::labs(
    title = paste(analyser_name, "measurement locations over time"),
    subtitle = "240-second intervals; first 60 seconds flushed",
    x = "DATE.TIME",
    y = "Location (MPVPosition)",
    colour = "Analyser"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    legend.position = "top",
    panel.grid.minor = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)
  )

plot_file <- file.path(
  plot_dir,
  paste0(analyser_name, "_measurement_locations_over_time.png")
)
ggplot2::ggsave(
  plot_file,
  location_plot,
  width = 14,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)
message("Plot saved to: ", plot_file)
}
