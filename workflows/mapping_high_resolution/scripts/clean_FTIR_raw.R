# Clean raw FTIR exports for the high-resolution mapping workflow.
#
# Input files are detected from raw_data/FTIR1_raw and the analyser name is
# taken from FTIR1 or FTIR2 in each filename. The script keeps one row per raw
# measurement and selects only the fields needed by the mapping workflow.

script_args <- commandArgs(trailingOnly = FALSE)
script_flag <- grep("^--file=", script_args, value = TRUE)

if (length(script_flag) == 1L) {
  script_path <- normalizePath(sub("^--file=", "", script_flag))
  workflow_dir <- dirname(dirname(script_path))
} else {
  workflow_dir <- normalizePath(getwd())
}

input_dir <- file.path(workflow_dir, "raw_data", "FTIR1_raw")
output_dir <- file.path(workflow_dir, "clean_data")
plot_dir <- file.path(workflow_dir, "plots")
date_time_cutoff <- as.POSIXct(
  "2024-12-31 23:59:59",
  format = "%Y-%m-%d %H:%M:%S",
  tz = "Europe/Berlin"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

input_files <- list.files(
  input_dir,
  pattern = "\\.TXT$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(input_files) == 0L) {
  stop("No .TXT files found in: ", input_dir)
}

clean_one_ftir <- function(path) {
  analyser_match <- regmatches(
    basename(path),
    regexpr("FTIR[12]", basename(path), ignore.case = TRUE)
  )

  if (length(analyser_match) == 0L || analyser_match == "") {
    stop("Could not determine analyser from filename: ", basename(path))
  }

  analyser_name <- toupper(analyser_match)

  raw <- read.delim(
    path,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    fileEncoding = "windows-1252"
  )

  required <- c("Messstelle", "Datum", "Zeit", "CO2", "CH4", "NH3", "H2O")
  missing_columns <- setdiff(required, names(raw))
  if (length(missing_columns) > 0L) {
    stop(
      "Missing required column(s) in ", basename(path), ": ",
      paste(missing_columns, collapse = ", ")
    )
  }

  to_number <- function(x) {
    as.numeric(gsub(",", ".", trimws(x), fixed = TRUE))
  }

  date_time <- as.POSIXct(
    paste(trimws(raw[["Datum"]]), trimws(raw[["Zeit"]])),
    format = "%Y-%m-%d %H:%M:%S",
    tz = "Europe/Berlin"
  )

  valid_measurement <- !is.na(date_time)
  excluded_rows <- sum(!valid_measurement)
  if (excluded_rows > 0L) {
    message(
      basename(path), ": excluded ", excluded_rows,
      " embedded header/non-measurement row(s)"
    )
  }

  raw <- raw[valid_measurement, , drop = FALSE]
  date_time <- date_time[valid_measurement]

  # Keep measurements through 31 December 2024, inclusive.
  cutoff <- as.POSIXct("2025-01-01 00:00:00", tz = "Europe/Berlin")
  within_date_range <- date_time < cutoff
  message(
    basename(path), ": excluded ", sum(!within_date_range),
    " measurement row(s) after 2024-12-31"
  )
  raw <- raw[within_date_range, , drop = FALSE]
  date_time <- date_time[within_date_range]

  within_date_range <- date_time <= date_time_cutoff
  excluded_after_cutoff <- sum(!within_date_range)
  if (excluded_after_cutoff > 0L) {
    message(
      basename(path), ": excluded ", excluded_after_cutoff,
      " measurement row(s) after 2024-12-31 23:59:59"
    )
  }

  raw <- raw[within_date_range, , drop = FALSE]
  date_time <- date_time[within_date_range]

  cleaned <- data.frame(
    DATE.TIME = format(date_time, "%Y-%m-%d %H:%M:%S"),
    analyser = analyser_name,
    location = trimws(as.character(raw[["Messstelle"]])),
    CO2 = to_number(raw[["CO2"]]),
    CH4 = to_number(raw[["CH4"]]),
    NH3 = to_number(raw[["NH3"]]),
    H2O = to_number(raw[["H2O"]]),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  cleaned
}

cleaned_files <- lapply(input_files, clean_one_ftir)
names(cleaned_files) <- vapply(
  input_files,
  function(path) toupper(regmatches(
    basename(path),
    regexpr("FTIR[12]", basename(path), ignore.case = TRUE)
  )),
  character(1)
)

for (analyser_name in names(cleaned_files)) {
  output_path <- file.path(
    output_dir,
    paste0("2024_05_17_", analyser_name, "_clean.csv")
  )
  write.csv(cleaned_files[[analyser_name]], output_path, row.names = FALSE, na = "")
  message(analyser_name, ": wrote ", nrow(cleaned_files[[analyser_name]]),
          " rows to ", output_path)
}

combined <- do.call(rbind, cleaned_files)
row.names(combined) <- NULL
combined <- combined[order(combined[["DATE.TIME"]], combined[["analyser"]]), ]

combined_path <- file.path(output_dir, "2024_05_17_FTIR_clean_combined.csv")
write.csv(combined, combined_path, row.names = FALSE, na = "")
message("Combined: wrote ", nrow(combined), " rows to ", combined_path)

# Plot the location measured by each FTIR through time.
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required to create the FTIR location plot.")
}

plot_data <- combined
plot_data[["DATE.TIME"]] <- as.POSIXct(
  plot_data[["DATE.TIME"]],
  format = "%Y-%m-%d %H:%M:%S",
  tz = "Europe/Berlin"
)
plot_data[["location_numeric"]] <- as.numeric(plot_data[["location"]])
plot_data[["location_dodged"]] <- plot_data[["location_numeric"]] +
  ifelse(plot_data[["analyser"]] == "FTIR1", -0.16, 0.16)
location_breaks <- sort(unique(plot_data[["location_numeric"]]))

location_plot <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(x = DATE.TIME, y = location_dodged, colour = analyser)
) +
  ggplot2::geom_point(alpha = 0.55, size = 0.65) +
  ggplot2::scale_colour_manual(
    values = c("FTIR1" = "#0072B2", "FTIR2" = "#D55E00")
  ) +
  ggplot2::scale_y_continuous(
    breaks = location_breaks,
    labels = location_breaks,
    minor_breaks = NULL
  ) +
  ggplot2::labs(
    title = "FTIR measurement locations over time (vertically dodged)",
    subtitle = "FTIR1 is offset below and FTIR2 above each measurement location",
    x = "DATE.TIME",
    y = "Location (Messstelle)",
    colour = "Analyser"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    legend.position = "top",
    panel.grid.minor = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)
  )

plot_path <- file.path(plot_dir, "FTIR_measurement_locations_over_time.png")
ggplot2::ggsave(
  filename = plot_path,
  plot = location_plot,
  width = 14,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)
message("Plot saved to ", plot_path)
