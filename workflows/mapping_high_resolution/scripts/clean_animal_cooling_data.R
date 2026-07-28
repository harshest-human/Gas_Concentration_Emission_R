###############################################################################
# Clean animal-counting and cooling-system records for campaigns 3 and 4.
#
# Important rules:
#   - Total herd size is 58 cows.
#   - When Out_Cow is supplied: inside_cow = 58 - out_cow.
#   - The first two files supply Inside_cow directly; cow_e and cow_w are ignored.
#   - Yes/No cooling-system values become 1/0.
#   - Known entry errors "yess" and "dss" are both treated as Yes (1).
#   - Local Excel clock times are preserved and assigned Europe/Berlin.
#   - No hourly aggregation or gas-data merge is performed in this script.
###############################################################################

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package 'data.table' is required.")
}
if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Package 'readxl' is required.")
}

library(data.table)

workflow_dir <- paste0(
  "D:/Data_Analysis_R/Gas_Concentration_Emission_R/",
  "workflows/mapping_high_resolution"
)
input_dir <- file.path(
  workflow_dir,
  "raw_data",
  "Animal_Counting_Cooling_System_Data"
)
campaign3_dir <- file.path(workflow_dir, "clean_data", "3_campaign")
campaign4_dir <- file.path(workflow_dir, "clean_data", "4_campaign")
dir.create(campaign3_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(campaign4_dir, recursive = TRUE, showWarnings = FALSE)

herd_size <- 58

campaign3_start <- as.POSIXct(
  "2025-08-19 13:00:12",
  format = "%Y-%m-%d %H:%M:%S",
  tz = "Europe/Berlin"
)
campaign3_end <- as.POSIXct(
  "2025-11-19 09:06:42",
  format = "%Y-%m-%d %H:%M:%S",
  tz = "Europe/Berlin"
)
campaign4_start <- as.POSIXct(
  "2025-12-09 12:03:59",
  format = "%Y-%m-%d %H:%M:%S",
  tz = "Europe/Berlin"
)
campaign4_end <- as.POSIXct(
  "2026-07-21 02:38:20",
  format = "%Y-%m-%d %H:%M:%S",
  tz = "Europe/Berlin"
)

output_columns <- c(
  "DATE.TIME",
  "inside_cow",
  "out_cow",
  "east_vent",
  "mid_vent",
  "west_vent",
  "water_spray"
)


###############################################################################
##### REUSABLE CLEANING FUNCTIONS
###############################################################################

# Make the changing workbook column names consistent.
standardize_column_names <- function(x) {
  names_clean <- gsub("[^a-z0-9]+", "_", tolower(names(x)))
  names_clean[names_clean %in% c("mid_went", "mid_vent")] <- "mid_vent"
  names_clean[names_clean %in% c("west_went", "west_vent")] <- "west_vent"
  names_clean[names_clean %in% c("water_spray", "water_spry")] <- "water_spray"
  setnames(x, names_clean)
  x
}


# Preserve the displayed Excel wall-clock time and assign Europe/Berlin.
# readxl stores Excel datetimes as UTC instants; formatting in UTC recovers the
# values displayed in the workbook without introducing a one- or two-hour shift.
excel_datetime_to_local <- function(x) {
  text_value <- format(
    as.POSIXct(x, origin = "1970-01-01", tz = "UTC"),
    "%Y-%m-%d %H:%M:%S",
    tz = "UTC"
  )
  as.POSIXct(
    text_value,
    format = "%Y-%m-%d %H:%M:%S",
    tz = "Europe/Berlin"
  )
}


# Build DATE.TIME for each workbook schema.
extract_date_time <- function(data) {
  if ("date_time" %in% names(data)) {
    return(excel_datetime_to_local(data[["date_time"]]))
  }

  if (!"date" %in% names(data)) {
    stop("No Date.Time or Date column found.")
  }

  # Most files store a complete Excel datetime in the Date column.
  if (inherits(data[["date"]], "POSIXt")) {
    return(excel_datetime_to_local(data[["date"]]))
  }

  # The first workbook stores Date as text and Time in a separate Excel column.
  if (!"time" %in% names(data)) {
    stop("Text Date column found without a Time column.")
  }

  date_text <- trimws(as.character(data[["date"]]))

  # Correct the documented obvious date-entry errors in the first workbook.
  date_text[date_text == "17/02/2025"] <- "17/07/2025"
  date_text <- sub("/2024$", "/2025", date_text)

  time_text <- format(
    as.POSIXct(data[["time"]], origin = "1970-01-01", tz = "UTC"),
    "%H:%M:%S",
    tz = "UTC"
  )

  as.POSIXct(
    paste(date_text, time_text),
    format = "%d/%m/%Y %H:%M:%S",
    tz = "Europe/Berlin"
  )
}


# Convert cooling-system states to binary values.
clean_binary_state <- function(x, column_name, file_name) {
  value <- tolower(trimws(as.character(x)))
  result <- rep(NA_integer_, length(value))
  result[value == "no"] <- 0L
  result[value %in% c("yes", "yess", "dss")] <- 1L

  unknown <- sort(unique(value[!is.na(value) & !value %in% c(
    "no", "yes", "yess", "dss"
  )]))
  if (length(unknown) > 0L) {
    warning(
      basename(file_name), ": unknown ", column_name, " value(s): ",
      paste(unknown, collapse = ", ")
    )
  }
  result
}


# Clean one source workbook.
clean_one_workbook <- function(path) {
  data <- as.data.table(readxl::read_excel(path))
  data <- standardize_column_names(data)

  required_conditions <- c(
    "east_vent", "mid_vent", "west_vent", "water_spray"
  )
  missing_conditions <- setdiff(required_conditions, names(data))
  if (length(missing_conditions) > 0L) {
    stop(
      basename(path), " is missing: ",
      paste(missing_conditions, collapse = ", ")
    )
  }

  date_time <- extract_date_time(data)

  if ("inside_cow" %in% names(data)) {
    inside_cow <- as.numeric(data[["inside_cow"]])
    invalid_inside <- !is.na(inside_cow) & (
      inside_cow < 0 | inside_cow > herd_size
    )
    if (any(invalid_inside)) {
      warning(
        basename(path), ": set ", sum(invalid_inside),
        " impossible Inside_cow value(s) outside 0:58 to NA"
      )
      inside_cow[invalid_inside] <- NA_real_
    }
    out_cow <- herd_size - inside_cow
  } else if ("out_cow" %in% names(data)) {
    out_cow <- as.numeric(data[["out_cow"]])
    invalid_out <- !is.na(out_cow) & (
      out_cow < 0 | out_cow > herd_size
    )
    if (any(invalid_out)) {
      warning(
        basename(path), ": set ", sum(invalid_out),
        " impossible Out_Cow value(s) outside 0:58 to NA"
      )
      out_cow[invalid_out] <- NA_real_
    }
    inside_cow <- herd_size - out_cow
  } else {
    stop(basename(path), " has neither Inside_cow nor Out_Cow.")
  }

  cleaned <- data.table(
    DATE.TIME = date_time,
    inside_cow = inside_cow,
    out_cow = out_cow,
    east_vent = clean_binary_state(
      data[["east_vent"]], "east_vent", path
    ),
    mid_vent = clean_binary_state(
      data[["mid_vent"]], "mid_vent", path
    ),
    west_vent = clean_binary_state(
      data[["west_vent"]], "west_vent", path
    ),
    water_spray = clean_binary_state(
      data[["water_spray"]], "water_spray", path
    )
  )

  if (anyNA(cleaned[["DATE.TIME"]])) {
    warning(
      basename(path), ": discarded ",
      sum(is.na(cleaned[["DATE.TIME"]])),
      " row(s) with an invalid timestamp"
    )
    cleaned <- cleaned[!is.na(DATE.TIME)]
  }

  cleaned
}


###############################################################################
##### IMPORT AND CLEAN ALL AVAILABLE FILES
###############################################################################

input_files <- list.files(
  input_dir,
  pattern = "\\.xlsx$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)
input_files <- input_files[!startsWith(basename(input_files), "~$")]

if (length(input_files) == 0L) {
  stop("No .xlsx source files found in: ", input_dir)
}

animal_cooling <- rbindlist(
  lapply(input_files, clean_one_workbook),
  use.names = TRUE
)
setorder(animal_cooling, DATE.TIME)

# Keep one record per 15-minute timestamp if adjacent files overlap.
animal_cooling <- unique(animal_cooling, by = "DATE.TIME")
setcolorder(animal_cooling, output_columns)


###############################################################################
##### CAMPAIGN 3 OUTPUT
###############################################################################

campaign3 <- animal_cooling[
  DATE.TIME >= campaign3_start & DATE.TIME <= campaign3_end
]
campaign3_path <- file.path(
  campaign3_dir,
  "animal_cooling_campaign3.csv"
)
fwrite(
  campaign3,
  campaign3_path,
  quote = FALSE,
  dateTimeAs = "write.csv"
)


###############################################################################
##### CAMPAIGN 4 OUTPUT
###############################################################################

campaign4 <- animal_cooling[
  DATE.TIME >= campaign4_start & DATE.TIME <= campaign4_end
]
campaign4_path <- file.path(
  campaign4_dir,
  "animal_cooling_campaign4.csv"
)
fwrite(
  campaign4,
  campaign4_path,
  quote = FALSE,
  dateTimeAs = "write.csv"
)


###############################################################################
##### COMPLETION SUMMARY
###############################################################################

message("Read ", length(input_files), " animal/cooling workbooks")
message(
  "Campaign 3: wrote ", nrow(campaign3), " rows to: ", campaign3_path
)
message(
  "Campaign 4: wrote ", nrow(campaign4), " rows to: ", campaign4_path
)
message("No hourly aggregation or gas-data merge was performed.")
