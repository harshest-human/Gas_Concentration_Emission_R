###############################################################################
# High-resolution barn gas mapping: campaigns 1 to 4
#
# This script:
#   1. cleans the raw FTIR exports with one reusable function;
#   2. reads and combines cleaned CRDS files with one reusable function;
#   3. applies the location layout used in each campaign;
#   4. preserves the special reference locations "s" and "out";
#   5. writes one CSV per campaign;
#   6. row-binds all campaigns and orders the final CSV by DATE.TIME; and
#   7. writes a text README documenting the decisions.
###############################################################################

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package 'data.table' is required.")
}

library(data.table)

workflow_dir <- paste0(
  "D:/Data_Analysis_R/Gas_Concentration_Emission_R/",
  "workflows/mapping_high_resolution"
)
clean_dir <- file.path(workflow_dir, "clean_data")
output_dir <- file.path(clean_dir, "campaign_combined")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_columns <- c(
  "DATE.TIME", "analyser", "location", "CO2", "CH4", "NH3", "H2O"
)


###############################################################################
##### REUSABLE FUNCTIONS
###############################################################################

# Convert DATE.TIME text to a Europe/Berlin POSIXct value.
parse_date_time <- function(x) {
  as.POSIXct(
    x,
    format = "%Y-%m-%d %H:%M:%S",
    tz = "Europe/Berlin"
  )
}


# Standardize the two non-numeric reference locations used in the campaigns.
#
# Source "in" is the internal/reference sampling line and is labelled "s".
# Source "S" is the outside-south sampling line and is labelled "out".
# Existing lowercase "s" and "out" labels are left unchanged.
standardize_reference_locations <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("in", "IN", "In")] <- "s"
  x[x == "S"] <- "out"
  x
}


# Clean one or more raw FTIR tab-delimited exports.
#
# The analyser is read from FTIR1/FTIR2 in the filename. Repeated embedded
# headers are removed by retaining only rows with a valid date and time.
# Gas values remain in the instrument units: gases in ppm and H2O in vol-%.
clean_ftir_raw <- function(input_files) {
  clean_one_file <- function(path) {
    analyser_match <- regmatches(
      basename(path),
      regexpr("FTIR[12]", basename(path), ignore.case = TRUE)
    )
    if (length(analyser_match) == 0L || analyser_match == "") {
      stop("Cannot identify FTIR analyser from: ", basename(path))
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

    required <- c(
      "Messstelle", "Datum", "Zeit", "CO2", "CH4", "NH3", "H2O"
    )
    missing_columns <- setdiff(required, names(raw))
    if (length(missing_columns) > 0L) {
      stop(
        basename(path), " is missing: ",
        paste(missing_columns, collapse = ", ")
      )
    }

    date_time <- as.POSIXct(
      paste(trimws(raw[["Datum"]]), trimws(raw[["Zeit"]])),
      format = "%Y-%m-%d %H:%M:%S",
      tz = "Europe/Berlin"
    )
    valid_measurement <- !is.na(date_time)
    raw <- raw[valid_measurement, , drop = FALSE]
    date_time <- date_time[valid_measurement]

    to_number <- function(x) {
      as.numeric(gsub(",", ".", trimws(x), fixed = TRUE))
    }

    data.table(
      DATE.TIME = date_time,
      analyser = analyser_name,
      location = trimws(as.character(raw[["Messstelle"]])),
      CO2 = to_number(raw[["CO2"]]),
      CH4 = to_number(raw[["CH4"]]),
      NH3 = to_number(raw[["NH3"]]),
      H2O = to_number(raw[["H2O"]])
    )
  }

  rbindlist(lapply(input_files, clean_one_file), use.names = TRUE)
}


# Read and combine cleaned CRDS CSV files.
#
# location_map is optional and maps MPVPosition values to campaign-specific
# barn locations. Positions not in location_map retain their source reference
# label, which is then standardized to "s" or "out".
combine_crds_files <- function(input_files, location_map = NULL) {
  read_one_file <- function(path) {
    column_names <- names(fread(path, nrows = 0, showProgress = FALSE))
    analyser_column <- if ("analyser" %in% column_names) {
      "analyser"
    } else if ("analyzer" %in% column_names) {
      "analyzer"
    } else {
      stop("No analyser/analyzer column in: ", basename(path))
    }

    required <- c(
      "DATE.TIME", analyser_column, "location", "CO2", "CH4", "NH3", "H2O"
    )
    if (!is.null(location_map)) {
      required <- c(required, "MPVPosition")
    }
    missing_columns <- setdiff(required, column_names)
    if (length(missing_columns) > 0L) {
      stop(
        basename(path), " is missing: ",
        paste(missing_columns, collapse = ", ")
      )
    }

    data <- fread(
      path,
      select = required,
      colClasses = list(
        character = c("DATE.TIME", analyser_column, "location")
      ),
      showProgress = FALSE
    )
    if (analyser_column == "analyzer") {
      setnames(data, "analyzer", "analyser")
    }

    data[, DATE.TIME := parse_date_time(DATE.TIME)]
    data[, location := standardize_reference_locations(location)]

    if (!is.null(location_map)) {
      mapped <- unname(location_map[as.character(data[["MPVPosition"]])])
      use_mapping <- !is.na(mapped)
      data[use_mapping, location := as.character(mapped[use_mapping])]
      data[, MPVPosition := NULL]
    }

    data[, ..output_columns]
  }

  combined <- rbindlist(lapply(input_files, read_one_file), use.names = TRUE)
  setorder(combined, DATE.TIME, analyser, location)
  combined
}


# Remap existing location numbers with a named vector.
remap_locations <- function(data, analyser_name, location_map) {
  rows <- data[["analyser"]] == analyser_name
  source_location <- as.character(data[["location"]][rows])
  mapped_location <- unname(location_map[source_location])

  if (anyNA(mapped_location)) {
    stop(
      "Unmapped ", analyser_name, " location(s): ",
      paste(sort(unique(source_location[is.na(mapped_location)])), collapse = ", ")
    )
  }
  data[rows, location := mapped_location]
  data
}


# Filter, row-bind, standardize, and deduplicate campaign data.
combine_campaign <- function(..., start_time = NULL, end_time = NULL) {
  data_list <- list(...)
  combined <- rbindlist(data_list, use.names = TRUE)
  combined[, location := standardize_reference_locations(location)]

  if (!is.null(start_time)) {
    combined <- combined[DATE.TIME >= parse_date_time(start_time)]
  }
  if (!is.null(end_time)) {
    combined <- combined[DATE.TIME <= parse_date_time(end_time)]
  }

  setorder(combined, DATE.TIME, analyser, location)
  combined <- unique(
    combined,
    by = c("DATE.TIME", "analyser", "location")
  )
  setcolorder(combined, output_columns)
  combined
}


# Write a campaign CSV using a consistent filename.
write_campaign_csv <- function(data, campaign_number) {
  path <- file.path(
    output_dir,
    paste0(
      "campaign_", campaign_number, "_",
      format(min(data[["DATE.TIME"]]), "%Y%m%d.%H%M%S"),
      "_",
      format(max(data[["DATE.TIME"]]), "%Y%m%d.%H%M%S"),
      ".csv"
    )
  )
  fwrite(data, path, quote = FALSE, dateTimeAs = "write.csv")
  message("Campaign ", campaign_number, ": wrote ", nrow(data), " rows")
  path
}


###############################################################################
##### CAMPAIGN 1: FTIR1 + FTIR2 + CRDS8 + CRDS9
###############################################################################

# Common measurement overlap used previously for the four analysers.
campaign1_start <- "2024-10-01 00:02:32"
campaign1_end <- "2024-10-24 14:38:58"

# FTIR raw exports are cleaned directly here.
campaign1_ftir_files <- file.path(
  workflow_dir,
  "raw_data",
  "FTIR_raw",
  c("2024_05_17_FTIR1.TXT", "2024_05_17_FTIR2.TXT")
)
campaign1_ftir <- clean_ftir_raw(campaign1_ftir_files)
campaign1_ftir <- campaign1_ftir[
  DATE.TIME >= parse_date_time(campaign1_start) &
    DATE.TIME <= parse_date_time(campaign1_end)
]

# Apply the original 51-location layout. FTIR1 position 10 is the special
# reference location "s" and is deliberately retained.
campaign1_ftir <- remap_locations(
  campaign1_ftir,
  "FTIR2",
  setNames(as.character(17:26), as.character(1:10))
)
campaign1_ftir <- remap_locations(
  campaign1_ftir,
  "FTIR1",
  setNames(c(as.character(43:51), "s"), as.character(1:10))
)

# CRDS step averages were already generated with 240-second intervals:
# first 60 seconds flushed, remaining observations averaged.
campaign1_crds_files <- file.path(
  clean_dir,
  "1_campaign",
  c(
    "20240516.172700_20241029.132912_CRDS8_240s_step_average.csv",
    "20240516.172700_20241029.132912_CRDS9_240s_step_average.csv"
  )
)
campaign1_crds <- combine_crds_files(campaign1_crds_files)
campaign1_crds <- remap_locations(
  campaign1_crds,
  "CRDS9",
  setNames(as.character(1:16), as.character(1:16))
)
campaign1_crds <- remap_locations(
  campaign1_crds,
  "CRDS8",
  setNames(as.character(27:42), as.character(1:16))
)

campaign1 <- combine_campaign(
  campaign1_ftir,
  campaign1_crds,
  start_time = campaign1_start,
  end_time = campaign1_end
)
campaign1_path <- write_campaign_csv(campaign1, 1)


###############################################################################
##### CAMPAIGN 2: CRDS8 + CRDS9, 32 LOCATIONS
###############################################################################

# The mid-height tubes and locations 19 and 40 were removed. These cleaned
# files already contain the confirmed serial CRDS9/CRDS8 location assignment.
campaign2_start <- "2024-11-16 00:06:47"
campaign2_end <- "2024-12-31 23:57:10"
campaign2_crds_files <- file.path(
  clean_dir,
  "2_campaign",
  c(
    "20241116.000000_20241231.235959_CRDS8_240s_step_average.csv",
    "20241116.000000_20241231.235959_CRDS9_240s_step_average.csv"
  )
)

campaign2_crds <- combine_crds_files(campaign2_crds_files)
campaign2 <- combine_campaign(
  campaign2_crds,
  start_time = campaign2_start,
  end_time = campaign2_end
)
campaign2_path <- write_campaign_csv(campaign2, 2)


###############################################################################
##### CAMPAIGN 3: SEGREGATED CLEAN CRDS8 + CRDS9 FILES
###############################################################################

# Use only the files in clean_data/3_campaign. Existing generated combined
# files are excluded so the script can be rerun safely.
campaign3_dir <- file.path(clean_dir, "3_campaign")
campaign3_crds_files <- list.files(
  campaign3_dir,
  pattern = "\\.csv$",
  full.names = TRUE,
  ignore.case = TRUE
)
campaign3_crds_files <- campaign3_crds_files[
  !grepl("_combined\\.csv$", basename(campaign3_crds_files))
]

campaign3_crds <- combine_crds_files(campaign3_crds_files)
campaign3 <- combine_campaign(campaign3_crds)
campaign3_path <- write_campaign_csv(campaign3, 3)


###############################################################################
##### CAMPAIGN 4: CRDS8, SEVEN BARN LOCATIONS + REFERENCES
###############################################################################

# Fourth-campaign MPVPosition mapping:
#   1 -> 12, 2 -> 18, 3 -> 15, 4 -> 30,
#   5 -> 36, 6 -> 42, 7 -> 45.
#
# MPVPosition 8 ("in") becomes "s" and MPVPosition 9 ("S") becomes "out".
campaign4_location_map <- c(
  "1" = "12",
  "2" = "18",
  "3" = "15",
  "4" = "30",
  "5" = "36",
  "6" = "42",
  "7" = "45"
)

campaign4_dir <- file.path(clean_dir, "4_campaign")
campaign4_crds_files <- list.files(
  campaign4_dir,
  pattern = "\\.csv$",
  full.names = TRUE,
  ignore.case = TRUE
)
campaign4_crds_files <- campaign4_crds_files[
  !grepl("_combined\\.csv$", basename(campaign4_crds_files))
]

campaign4_crds <- combine_crds_files(
  campaign4_crds_files,
  location_map = campaign4_location_map
)
campaign4 <- combine_campaign(campaign4_crds)
campaign4_path <- write_campaign_csv(campaign4, 4)


###############################################################################
##### FINAL: ROW-BIND CAMPAIGNS 1, 2, 3, AND 4
###############################################################################

all_campaigns <- rbindlist(
  list(campaign1, campaign2, campaign3, campaign4),
  use.names = TRUE
)
setorder(all_campaigns, DATE.TIME, analyser, location)
setcolorder(all_campaigns, output_columns)

final_csv_path <- file.path(
  output_dir,
  "all_campaigns_combined.csv"
)
fwrite(
  all_campaigns,
  final_csv_path,
  quote = FALSE,
  dateTimeAs = "write.csv"
)


###############################################################################
##### CAMPAIGN README
###############################################################################

reference_summary <- function(data) {
  references <- sort(unique(
    data[["location"]][data[["location"]] %in% c("s", "out")]
  ))
  if (length(references) == 0L) "none measured" else paste(references, collapse = ", ")
}

readme_lines <- c(
  "HIGH-RESOLUTION MAPPING: CAMPAIGN CLEANING README",
  "================================================",
  "",
  "Final columns",
  "-------------",
  "DATE.TIME, analyser, location, CO2, CH4, NH3, H2O",
  "",
  "General processing",
  "------------------",
  "- All timestamps use the Europe/Berlin timezone.",
  "- Final rows are ordered by DATE.TIME, analyser, and location.",
  "- Duplicate DATE.TIME + analyser + location records are removed within each campaign.",
  "- CO2, CH4, and NH3 are in ppm; H2O is in vol-%.",
  "- CRDS step-average files use a 240-second measurement interval.",
  "- For CRDS raw cleaning, the first 60 seconds are flushing time and are discarded;",
  "  the remaining observations in the interval are averaged.",
  "- Pre-cleaned CRDS files are not averaged or flushed a second time.",
  "- FTIR exports are cleaned by removing embedded headers and selecting the required gases.",
  "- FTIR measurements are retained at their instrument timestamps; no extra flush is applied.",
  "",
  "Reference locations",
  "-------------------",
  "- Source label 'in' is standardized to 's'.",
  "- Source label 'S' is standardized to 'out'.",
  "- Existing lowercase 's' and 'out' values are retained.",
  "- Reference rows are preserved and are not silently discarded.",
  "",
  "Campaign 1",
  "----------",
  paste0("- Period: ", campaign1_start, " to ", campaign1_end),
  "- Analysers: FTIR1, FTIR2, CRDS8, CRDS9.",
  "- CRDS9 positions 1:16 -> locations 1:16.",
  "- FTIR2 positions 1:10 -> locations 17:26.",
  "- CRDS8 positions 1:16 -> locations 27:42.",
  "- FTIR1 positions 1:9 -> locations 43:51; position 10 -> s.",
  paste0("- Reference locations present: ", reference_summary(campaign1), "."),
  paste0("- Output rows: ", nrow(campaign1), "."),
  "",
  "Campaign 2",
  "----------",
  paste0("- Period: ", campaign2_start, " to ", campaign2_end),
  "- Analysers: CRDS8 and CRDS9.",
  "- Mid-height locations plus locations 19 and 40 were removed from the layout.",
  "- CRDS9 MPVPositions 1:16 -> 1,3,4,6,7,9,10,12,13,15,16,18,21,22,24,25.",
  "- CRDS8 MPVPositions 1:16 -> 27,28,30,31,33,34,36,37,39,42,43,45,46,48,49,51.",
  paste0("- Reference locations present: ", reference_summary(campaign2), "."),
  paste0("- Output rows: ", nrow(campaign2), "."),
  "",
  "Campaign 3",
  "----------",
  "- Source: only CSV files segregated in clean_data/3_campaign.",
  "- Numeric saved locations are preserved.",
  "- Source 'in' and 'S' locations are retained as s and out.",
  "- Most files are 240-second averages; the supplied 450avg file is used as supplied.",
  paste0(
    "- Period: ",
    format(min(campaign3[["DATE.TIME"]]), "%Y-%m-%d %H:%M:%S"),
    " to ",
    format(max(campaign3[["DATE.TIME"]]), "%Y-%m-%d %H:%M:%S"),
    "."
  ),
  paste0("- Reference locations present: ", reference_summary(campaign3), "."),
  paste0("- Output rows: ", nrow(campaign3), "."),
  "",
  "Campaign 4",
  "----------",
  "- Analyser: CRDS8.",
  "- MPVPosition mapping: 1->12, 2->18, 3->15, 4->30, 5->36, 6->42, 7->45.",
  "- MPVPosition 8 ('in') -> s; MPVPosition 9 ('S') -> out.",
  "- Overlapping source exports are deduplicated.",
  paste0(
    "- Period: ",
    format(min(campaign4[["DATE.TIME"]]), "%Y-%m-%d %H:%M:%S"),
    " to ",
    format(max(campaign4[["DATE.TIME"]]), "%Y-%m-%d %H:%M:%S"),
    "."
  ),
  paste0("- Reference locations present: ", reference_summary(campaign4), "."),
  paste0("- Output rows: ", nrow(campaign4), "."),
  "",
  "Final combined dataset",
  "----------------------",
  paste0("- Total rows: ", nrow(all_campaigns), "."),
  paste0(
    "- Overall period: ",
    format(min(all_campaigns[["DATE.TIME"]]), "%Y-%m-%d %H:%M:%S"),
    " to ",
    format(max(all_campaigns[["DATE.TIME"]]), "%Y-%m-%d %H:%M:%S"),
    "."
  ),
  paste0("- Final CSV: ", final_csv_path),
  "",
  "Campaign CSV files",
  "------------------",
  paste0("- ", campaign1_path),
  paste0("- ", campaign2_path),
  paste0("- ", campaign3_path),
  paste0("- ", campaign4_path)
)

readme_path <- file.path(output_dir, "campaign_readme.txt")
writeLines(readme_lines, readme_path, useBytes = TRUE)

message("Final combined rows: ", nrow(all_campaigns))
message("Final CSV: ", final_csv_path)
message("README: ", readme_path)
