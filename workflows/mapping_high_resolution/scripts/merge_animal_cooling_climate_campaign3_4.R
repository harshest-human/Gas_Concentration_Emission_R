###############################################################################
# Add five-minute climate observations to the Campaign 3 and Campaign 4
# animal/cooling datasets.
#
# Temporal handling:
#   - A regular five-minute sequence spanning each animal/cooling campaign
#     defines the output rows.
#   - Temperature/humidity timestamps containing +0000 are parsed as UTC and
#     converted to Europe/Berlin.
#   - USA 16 wind date/time fields are parsed as UTC and converted to
#     Europe/Berlin. The complete 02:00--02:55 sequence recorded during the
#     spring DST transition confirms that the mast logger clock is UTC.
#   - Wind observations are matched only at identical five-minute timestamps.
#   - Each 15-minute animal/cooling state is carried forward for at most
#     14 minutes 59 seconds, so it applies to its corresponding three
#     five-minute climate rows. No temporal averaging or interpolation is used.
#   - Missing climate or wind measurements remain NA. In particular, no wind
#     source file is available for April 2026.
#
# Source-preservation:
#   - Before the named campaign outputs are replaced, their original 15-minute
#     versions are copied once to files ending in "_15min_source.csv".
###############################################################################

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package 'data.table' is required.")
}

library(data.table)

workflow_dir <- paste0(
  "D:/Data_Analysis_R/Gas_Concentration_Emission_R/",
  "workflows/mapping_high_resolution"
)
climate_dir <- file.path(workflow_dir, "raw_data", "CLIMATE_DATA")

campaign_files <- c(
  campaign3 = file.path(
    workflow_dir,
    "clean_data",
    "3_campaign",
    "animal_cooling_campaign3.csv"
  ),
  campaign4 = file.path(
    workflow_dir,
    "clean_data",
    "4_campaign",
    "animal_cooling_campaign4.csv"
  )
)

animal_columns <- c(
  "inside_cow",
  "out_cow",
  "east_vent",
  "mid_vent",
  "west_vent",
  "water_spray"
)

climate_columns <- c(
  "T_Outside",
  "RH_Outside",
  "T_inside",
  "RH_inside"
)

output_columns <- c(
  "DATE.TIME",
  "DATE.TIME.UTC",
  "UTC_offset",
  animal_columns,
  climate_columns,
  "wind_speed",
  "wind_direction",
  "station_id"
)


###############################################################################
##### REUSABLE IMPORT AND MERGE FUNCTIONS
###############################################################################

# Display a UTC instant in Europe/Berlin without changing or reparsing the
# underlying instant. This preserves both occurrences of the autumn clock hour.
utc_to_berlin <- function(x) {
  x <- as.POSIXct(x, origin = "1970-01-01", tz = "UTC")
  attr(x, "tzone") <- "Europe/Berlin"
  x
}


# Read one temperature/humidity CSV.
read_trh_file <- function(path) {
  x <- fread(path, na.strings = c("", "NA", "NaN"))

  names_lower <- tolower(names(x))
  setnames(x, names(x), names_lower)
  setnames(
    x,
    old = intersect(c("t_outside", "t_out"), names(x)),
    new = rep("T_Outside", length(intersect(c("t_outside", "t_out"), names(x))))
  )
  setnames(
    x,
    old = intersect(c("rh_outside", "rh_out"), names(x)),
    new = rep("RH_Outside", length(intersect(c("rh_outside", "rh_out"), names(x))))
  )
  setnames(x, "t_inside", "T_inside", skip_absent = TRUE)
  setnames(x, "rh_inside", "RH_inside", skip_absent = TRUE)

  required <- c("date", climate_columns)
  missing <- setdiff(required, names(x))
  if (length(missing) > 0L) {
    stop(basename(path), " is missing: ", paste(missing, collapse = ", "))
  }

  # Files use month/day/year. Newer exports include the UTC suffix in one
  # Date field; older exports split the same UTC timestamp across date/time.
  if ("time" %in% names(x)) {
    timestamp_text <- paste(x[["date"]], x[["time"]])
    timestamp_format <- "%m/%d/%y %H:%M:%S"
  } else {
    timestamp_text <- x[["date"]]
    timestamp_format <- "%m/%d/%y %H:%M:%S %z"
  }

  timestamp_utc <- as.POSIXct(
    timestamp_text,
    format = timestamp_format,
    tz = "UTC"
  )

  if (anyNA(timestamp_utc)) {
    stop(
      basename(path), " contains ",
      sum(is.na(timestamp_utc)), " unparseable timestamp(s)."
    )
  }

  data.table(
    DATE.TIME = utc_to_berlin(timestamp_utc),
    T_Outside = as.numeric(x[["T_Outside"]]),
    RH_Outside = as.numeric(x[["RH_Outside"]]),
    T_inside = as.numeric(x[["T_inside"]]),
    RH_inside = as.numeric(x[["RH_inside"]])
  )
}


# Read one USA 16 wind CSV.
read_wind_file <- function(path) {
  x <- fread(path, na.strings = c("", "NA", "NaN"))

  required <- c(
    "station_id", "date", "time", "wind_speed", "wind_direction"
  )
  missing <- setdiff(required, names(x))
  if (length(missing) > 0L) {
    stop(basename(path), " is missing: ", paste(missing, collapse = ", "))
  }

  timestamp_utc <- as.POSIXct(
    paste(x[["date"]], x[["time"]]),
    format = "%d.%m.%Y %H:%M:%S",
    tz = "UTC"
  )

  if (anyNA(timestamp_utc)) {
    stop(
      basename(path), " contains ",
      sum(is.na(timestamp_utc)), " unparseable timestamp(s)."
    )
  }

  data.table(
    DATE.TIME = utc_to_berlin(timestamp_utc),
    wind_speed = as.numeric(x[["wind_speed"]]),
    wind_direction = as.numeric(x[["wind_direction"]]),
    station_id = as.character(x[["station_id"]])
  )
}


# Retain one row per timestamp and warn when overlapping source files disagree.
deduplicate_time <- function(x, source_name) {
  duplicate_times <- x[duplicated(DATE.TIME) | duplicated(DATE.TIME, fromLast = TRUE)]

  if (nrow(duplicate_times) > 0L) {
    warning(
      source_name, ": found ", uniqueN(duplicate_times[["DATE.TIME"]]),
      " overlapping timestamp(s); retaining the first record."
    )
  }

  setorder(x, DATE.TIME)
  unique(x, by = "DATE.TIME")
}


# Add animal/cooling states and climate observations to a five-minute timeline.
merge_one_campaign <- function(animal_path, climate, campaign_name) {
  backup_path <- sub(
    "\\.csv$",
    "_15min_source.csv",
    animal_path,
    ignore.case = TRUE
  )
  animal_source_path <- if (file.exists(backup_path)) backup_path else animal_path

  animal <- fread(
    animal_source_path,
    na.strings = c("", "NA", "NaN"),
    colClasses = c("DATE.TIME" = "character")
  )

  required <- c("DATE.TIME", animal_columns)
  missing <- setdiff(required, names(animal))
  if (length(missing) > 0L) {
    stop(basename(animal_path), " is missing: ", paste(missing, collapse = ", "))
  }

  animal[, DATE.TIME := as.POSIXct(
    DATE.TIME,
    format = "%Y-%m-%d %H:%M:%S",
    tz = "Europe/Berlin"
  )]
  setorder(animal, DATE.TIME)

  campaign_start <- min(animal[["DATE.TIME"]], na.rm = TRUE)
  campaign_end <- max(animal[["DATE.TIME"]], na.rm = TRUE)

  output <- data.table(
    DATE.TIME = seq(
      from = campaign_start,
      to = campaign_end,
      by = 5 * 60
    )
  )
  output <- merge(
    output,
    climate,
    by = "DATE.TIME",
    all.x = TRUE,
    sort = TRUE
  )

  # Retain a unique UTC timestamp and explicit offset because the autumn clock
  # change repeats 02:00--02:55 in Europe/Berlin.
  output[, `:=`(
    DATE.TIME.UTC = format(
      DATE.TIME,
      tz = "UTC",
      format = "%Y-%m-%d %H:%M:%S"
    ),
    UTC_offset = format(
      DATE.TIME,
      tz = "Europe/Berlin",
      format = "%z"
    )
  )]

  # Prepare local wall-clock representations for the autumn fallback only.
  output_wall_time <- as.POSIXct(
    format(
      output[["DATE.TIME"]],
      tz = "Europe/Berlin",
      format = "%Y-%m-%d %H:%M:%S"
    ),
    format = "%Y-%m-%d %H:%M:%S",
    tz = "UTC"
  )
  animal_wall_time <- as.POSIXct(
    format(
      animal[["DATE.TIME"]],
      tz = "Europe/Berlin",
      format = "%Y-%m-%d %H:%M:%S"
    ),
    format = "%Y-%m-%d %H:%M:%S",
    tz = "UTC"
  )

  # First match by true elapsed time. This handles the spring clock jump, where
  # 01:59:59 and 03:00:00 are only one second apart.
  animal_index <- findInterval(
    as.numeric(output[["DATE.TIME"]]),
    as.numeric(animal[["DATE.TIME"]])
  )
  valid_index <- animal_index > 0L

  elapsed_seconds <- rep(Inf, nrow(output))
  elapsed_seconds[valid_index] <- (
    as.numeric(output[["DATE.TIME"]][valid_index]) -
      as.numeric(animal[["DATE.TIME"]][animal_index[valid_index]])
  )
  valid_index <- valid_index & elapsed_seconds >= 0 & elapsed_seconds < 15 * 60

  # For still-unmatched rows, use local wall-clock matching. This assigns the
  # same animal state to both occurrences of the repeated autumn 02:00 hour.
  fallback_rows <- which(!valid_index)
  if (length(fallback_rows) > 0L) {
    fallback_index <- findInterval(
      as.numeric(output_wall_time[fallback_rows]),
      as.numeric(animal_wall_time)
    )
    fallback_valid <- fallback_index > 0L
    fallback_elapsed <- rep(Inf, length(fallback_rows))
    fallback_elapsed[fallback_valid] <- (
      as.numeric(output_wall_time[fallback_rows[fallback_valid]]) -
        as.numeric(animal_wall_time[fallback_index[fallback_valid]])
    )
    fallback_valid <- (
      fallback_valid &
        fallback_elapsed >= 0 &
        fallback_elapsed < 15 * 60
    )

    accepted_rows <- fallback_rows[fallback_valid]
    animal_index[accepted_rows] <- fallback_index[fallback_valid]
    valid_index[accepted_rows] <- TRUE
  }

  for (column_name in animal_columns) {
    set(output, j = column_name, value = rep(NA_real_, nrow(output)))
    set(
      output,
      i = which(valid_index),
      j = column_name,
      value = animal[[column_name]][animal_index[valid_index]]
    )
  }

  setcolorder(output, output_columns)
  setorder(output, DATE.TIME)

  if (!file.exists(backup_path)) {
    if (!file.copy(animal_path, backup_path, overwrite = FALSE)) {
      stop("Could not preserve source file: ", backup_path)
    }
  }

  fwrite(
    output,
    animal_path,
    quote = FALSE,
    dateTimeAs = "write.csv",
    na = "NA"
  )

  message(
    campaign_name, ": wrote ", nrow(output),
    " five-minute rows from ",
    format(min(output[["DATE.TIME"]]), "%Y-%m-%d %H:%M:%S"),
    " to ",
    format(max(output[["DATE.TIME"]]), "%Y-%m-%d %H:%M:%S")
  )
  message(
    campaign_name, ": rows without matched animal state = ",
    sum(!valid_index)
  )
  message(
    campaign_name, ": rows without wind speed = ",
    sum(is.na(output[["wind_speed"]]))
  )
  message(
    campaign_name, ": rows without temperature/humidity = ",
    sum(is.na(output[["T_Outside"]]))
  )

  invisible(output)
}


###############################################################################
##### BUILD THE FIVE-MINUTE CLIMATE MASTER DATASET
###############################################################################

trh_files <- list.files(
  climate_dir,
  pattern = "\\.csv$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)
trh_files <- trh_files[grepl("TRH", trh_files, ignore.case = TRUE)]

wind_files <- list.files(
  climate_dir,
  pattern = "\\.csv$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)
wind_files <- wind_files[grepl("WS_WD", wind_files, ignore.case = TRUE)]

if (length(trh_files) == 0L) {
  stop("No temperature/humidity CSV files found in: ", climate_dir)
}
if (length(wind_files) == 0L) {
  stop("No wind CSV files found in: ", climate_dir)
}

temperature_humidity <- rbindlist(
  lapply(trh_files, read_trh_file),
  use.names = TRUE
)
temperature_humidity <- deduplicate_time(
  temperature_humidity,
  "Temperature/humidity data"
)

wind <- rbindlist(
  lapply(wind_files, read_wind_file),
  use.names = TRUE
)
wind <- deduplicate_time(wind, "Wind data")

# Temperature/humidity and wind are matched at identical five-minute
# timestamps. The complete campaign timeline is created later from the
# animal/cooling campaign bounds, so missing climate intervals remain explicit.
climate <- merge(
  temperature_humidity,
  wind,
  by = "DATE.TIME",
  all.x = TRUE,
  sort = TRUE
)


###############################################################################
##### WRITE CAMPAIGN 3 AND CAMPAIGN 4 OUTPUTS
###############################################################################

campaign3 <- merge_one_campaign(
  campaign_files[["campaign3"]],
  climate,
  "Campaign 3"
)

campaign4 <- merge_one_campaign(
  campaign_files[["campaign4"]],
  climate,
  "Campaign 4"
)

message("No temporal averaging or wind interpolation was performed.")
