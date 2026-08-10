###############################################################################
# High-resolution gas mapping: unified preparation for Campaigns 1--4
#
# This is the single preparation script for analyser-wise cleaning, sampling-
# point mapping, campaign-wise aggregation, and final row-binding.
#
# Terminology:
#   sampling.point = physical sampling point in or around the barn
#   MPVPosition    = CRDS valve-selector position
#   Messstelle     = FTIR selector position
#
# Existing accepted data are never overwritten. All outputs are written to
# clean_data/prepared_campaigns_v01 for validation before adoption.
###############################################################################

suppressPackageStartupMessages(library(data.table))

workflow <- paste0(
  "D:/Data_Analysis_R/Gas_Concentration_Emission_R/",
  "workflows/mapping_high_resolution"
)
metadata_dir <- file.path(workflow, "meta_data")
output_root <- file.path(workflow, "clean_data", "prepared_campaigns_v01")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

tz_local <- "Europe/Berlin"
gas_columns <- c("CO2", "CH4", "NH3", "H2O")
ftir_divisors <- c(CO2 = 1.06, CH4 = 1.06, NH3 = 1.09, H2O = 1.00)


###############################################################################
##### GENERAL HELPERS
###############################################################################

to_number <- function(x) {
  suppressWarnings(as.numeric(gsub(",", ".", trimws(as.character(x)),
                                fixed = TRUE)))
}

parse_local_time <- function(x) {
  x <- trimws(as.character(x))
  out <- as.POSIXct(x, format = "%Y-%m-%d %H:%M:%S", tz = tz_local)
  missing <- is.na(out)
  if (any(missing)) {
    out[missing] <- as.POSIXct(
      x[missing], format = "%d.%m.%Y %H:%M:%S", tz = tz_local
    )
  }
  out
}

normalise_analyser <- function(x) {
  toupper(gsub("[^A-Za-z0-9]", "", trimws(as.character(x))))
}

assert_columns <- function(x, required, label) {
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop(label, " is missing: ", paste(missing, collapse = ", "))
  }
}

read_mapping <- function(campaign_number) {
  path <- file.path(
    metadata_dir,
    paste0("campaign", campaign_number,
           "_sampling_point_analyser_mapping.csv")
  )
  map <- fread(path, colClasses = "character")
  assert_columns(
    map,
    c("sampling.point", "analyser", "Messstelle_or_MPVPosition"),
    basename(path)
  )
  map[, analyser := normalise_analyser(analyser)]
  setnames(map, "Messstelle_or_MPVPosition", "source.position")
  unique(map)
}

add_instrument_scale_columns <- function(x) {
  for (g in gas_columns) {
    raw_name <- paste0(g, "_raw")
    corrected_name <- paste0(g, "_corr")
    x[, (raw_name) := get(g)]
    x[, (corrected_name) := fifelse(
      grepl("^FTIR", analyser),
      get(g) / unname(ftir_divisors[g]),
      get(g)
    )]
    x[, (g) := get(corrected_name)]
  }
  x[, correction_applied := as.integer(grepl("^FTIR", analyser))]
  x
}

filter_period <- function(x, start, end) {
  start <- as.POSIXct(start, tz = tz_local)
  end <- as.POSIXct(end, tz = tz_local)
  x[DATE.TIME >= start & DATE.TIME <= end]
}


###############################################################################
##### FTIR RAW READER
###############################################################################

read_ftir_raw <- function(files) {
  read_one <- function(path) {
    analyser <- regmatches(
      basename(path), regexpr("FTIR[12]", basename(path), ignore.case = TRUE)
    )
    if (!length(analyser) || !nzchar(analyser)) {
      stop("Cannot identify FTIR1/FTIR2 from filename: ", basename(path))
    }
    analyser <- toupper(analyser)

    raw <- fread(
      path, sep = "\t", encoding = "Latin-1", fill = TRUE,
      check.names = TRUE, showProgress = FALSE
    )
    assert_columns(
      raw, c("Messstelle", "Datum", "Zeit", gas_columns), basename(path)
    )

    out <- data.table(
      DATE.TIME = parse_local_time(paste(raw$Datum, raw$Zeit)),
      analyser = analyser,
      source.variable = "Messstelle",
      source.position = trimws(as.character(raw$Messstelle))
    )
    for (g in gas_columns) out[, (g) := to_number(raw[[g]])]
    out[!is.na(DATE.TIME) & nzchar(source.position)]
  }

  rbindlist(lapply(files, read_one), use.names = TRUE, fill = TRUE)
}


###############################################################################
##### CRDS RAW READER AND STEP AVERAGER
###############################################################################

read_crds_raw <- function(files, analyser) {
  analyser <- normalise_analyser(analyser)

  read_one <- function(path) {
    header <- names(fread(path, nrows = 0, showProgress = FALSE))
    required_gases <- intersect(gas_columns, header)
    if (!all(gas_columns %in% required_gases) ||
        !"MPVPosition" %in% header) {
      stop("Required CRDS fields are missing from: ", basename(path))
    }

    if ("DATE.TIME" %in% header) {
      keep <- c("DATE.TIME", "MPVPosition", gas_columns)
      x <- fread(path, select = keep, showProgress = FALSE)
      x[, DATE.TIME := parse_local_time(DATE.TIME)]
    } else {
      assert_columns(as.list(setNames(header, header)), c("DATE", "TIME"),
                     basename(path))
      keep <- c("DATE", "TIME", "MPVPosition", gas_columns)
      x <- fread(path, select = keep, showProgress = FALSE)
      x[, DATE.TIME := parse_local_time(paste(DATE, TIME))]
      x[, c("DATE", "TIME") := NULL]
    }

    x[, analyser := analyser]
    x[, source.variable := "MPVPosition"]
    x[, source.position := trimws(as.character(MPVPosition))]
    x[, MPVPosition := NULL]
    for (g in gas_columns) set(x, j = g, value = to_number(x[[g]]))
    x[!is.na(DATE.TIME) & nzchar(source.position)]
  }

  setorder(rbindlist(lapply(files, read_one), use.names = TRUE), DATE.TIME)
}

average_crds_steps <- function(raw, flush_seconds = 60,
                               nominal_step_seconds = 240,
                               maximum_gap_seconds = 10) {
  setorder(raw, analyser, DATE.TIME)
  raw[, gap := as.numeric(difftime(DATE.TIME, shift(DATE.TIME), units = "secs")),
      by = analyser]
  raw[, new_step := is.na(gap) | gap > maximum_gap_seconds |
        source.position != shift(source.position), by = analyser]
  raw[, step.id := cumsum(new_step), by = analyser]
  raw[, elapsed.seconds := as.numeric(
    difftime(DATE.TIME, min(DATE.TIME), units = "secs")
  ), by = .(analyser, step.id)]

  averaged <- raw[
    elapsed.seconds >= flush_seconds & elapsed.seconds < nominal_step_seconds,
    c(
      list(
        DATE.TIME = max(DATE.TIME),
        source.variable = first(source.variable),
        source.position = first(source.position),
        flush.seconds = flush_seconds,
        averaging.seconds = nominal_step_seconds - flush_seconds,
        n.raw = .N
      ),
      lapply(.SD, mean, na.rm = TRUE)
    ),
    by = .(analyser, step.id),
    .SDcols = gas_columns
  ]
  averaged[, step.id := NULL]
  averaged
}


###############################################################################
##### PRE-CLEANED CRDS STEP-AVERAGE READER
###############################################################################

read_crds_step_average <- function(files) {
  read_one <- function(path) {
    header <- names(fread(path, nrows = 0, showProgress = FALSE))
    analyser_field <- if ("analyser" %in% header) "analyser" else "analyzer"
    selector_field <- if ("MPVPosition" %in% header) {
      "MPVPosition"
    } else if ("location" %in% header) {
      "location"
    } else {
      stop(basename(path), " has neither MPVPosition nor location.")
    }
    assert_columns(
      as.list(setNames(header, header)),
      c("DATE.TIME", analyser_field, selector_field, gas_columns),
      basename(path)
    )
    x <- fread(path, showProgress = FALSE)
    setnames(x, analyser_field, "analyser")
    x[, `:=`(
      DATE.TIME = parse_local_time(DATE.TIME),
      analyser = normalise_analyser(analyser),
      source.variable = "MPVPosition",
      source.position = trimws(as.character(get(selector_field)))
    )]
    for (g in gas_columns) set(x, j = g, value = to_number(x[[g]]))
    x <- x[!is.na(DATE.TIME)]
    x[, .SD, .SDcols = c(
      "DATE.TIME", "analyser", "source.variable", "source.position",
      gas_columns, intersect(c("location", "dwell_seconds"), names(x))
    )]
  }
  rbindlist(lapply(files, read_one), use.names = TRUE, fill = TRUE)
}


###############################################################################
##### SAMPLING-POINT MAPPING AND OUTPUT
###############################################################################

apply_sampling_map <- function(x, campaign_number,
                               preserve_existing_location = FALSE) {
  map <- read_mapping(campaign_number)

  if (preserve_existing_location && "location" %in% names(x)) {
    x[, sampling.point := trimws(as.character(location))]
    x[sampling.point == "", sampling.point := NA_character_]
  } else {
    x[, sampling.point := NA_character_]
  }

  unresolved <- is.na(x$sampling.point)
  if (any(unresolved)) {
    mapped <- map[x[unresolved],
      on = .(analyser, source.position),
      .(sampling.point = x.sampling.point),
      allow.cartesian = TRUE
    ]
    if (nrow(mapped) != sum(unresolved)) {
      stop("Ambiguous mapping in Campaign ", campaign_number,
           "; inspect repeated analyser/source-position assignments.")
    }
    x[unresolved, sampling.point := mapped$sampling.point]
  }

  if (anyNA(x$sampling.point)) {
    bad <- unique(x[is.na(sampling.point), .(analyser, source.position)])
    stop("Unmapped selector positions in Campaign ", campaign_number, ": ",
         paste(paste(bad$analyser, bad$source.position, sep = ":"),
               collapse = ", "))
  }

  x[, campaign := campaign_number]
  x[, location := NULL]
  add_instrument_scale_columns(x)
}

write_campaign_outputs <- function(x, campaign_number) {
  campaign_dir <- file.path(output_root, paste0("campaign_", campaign_number))
  analyser_dir <- file.path(campaign_dir, "analyser_wise")
  dir.create(analyser_dir, recursive = TRUE, showWarnings = FALSE)

  detailed_columns <- c(
    "DATE.TIME", "campaign", "analyser", "source.variable",
    "source.position", "sampling.point",
    paste0(gas_columns, "_raw"), paste0(gas_columns, "_corr"), gas_columns,
    "correction_applied"
  )
  x <- x[, ..detailed_columns]
  setorder(x, DATE.TIME, analyser, sampling.point)
  # Overlapping source exports can contain the same physical observation.
  # Retain the first record after deterministic chronological ordering.
  x <- unique(x, by = c("DATE.TIME", "analyser", "sampling.point"))

  for (a in sort(unique(x$analyser))) {
    fwrite(
      x[analyser == a],
      file.path(analyser_dir, paste0(a, "_clean_mapped.csv")),
      dateTimeAs = "write.csv"
    )
  }

  compact <- x[, .(
    DATE.TIME, campaign, analyser, sampling.point, CO2, CH4, NH3, H2O
  )]
  fwrite(
    compact,
    file.path(campaign_dir, paste0("campaign_", campaign_number,
                                  "_combined.csv")),
    dateTimeAs = "write.csv"
  )
  compact
}


###############################################################################
##### EXECUTIVE RUNS
###############################################################################

ftir_files <- list.files(
  file.path(workflow, "raw_data", "FTIR_raw"),
  pattern = "FTIR[12].*[.]TXT$", full.names = TRUE, ignore.case = TRUE
)
ftir_all <- read_ftir_raw(ftir_files)

# Campaign 1: summer and autumn Campaign 1 intervals.
c1_crds_files <- c(
  list.files(
    file.path(workflow, "clean_data", "1_campaign",
              "recovered_crds_june_august_v03"),
    pattern = "recovered_step_average_.*_v03[.]csv$",
    recursive = TRUE, full.names = TRUE
  ),
  file.path(workflow, "clean_data", "1_campaign", "recovered_crds_v01",
            "campaign1_CRDS8_CRDS9_recovered_step_average_v01.csv")
)
c1_crds_files <- c1_crds_files[file.exists(c1_crds_files)]
c1_crds <- read_crds_step_average(c1_crds_files)
c1_crds <- c1_crds[is.na(dwell_seconds) | dwell_seconds == 240]
c1_ftir <- ftir_all[
  (DATE.TIME >= as.POSIXct("2024-06-01", tz = tz_local) &
     DATE.TIME <= as.POSIXct("2024-08-31 23:59:59", tz = tz_local)) |
  (DATE.TIME >= as.POSIXct("2024-10-01", tz = tz_local) &
     DATE.TIME <= as.POSIXct("2024-10-24 23:59:59", tz = tz_local))
]
c1 <- rbindlist(list(c1_crds, c1_ftir), use.names = TRUE, fill = TRUE)
c1 <- apply_sampling_map(c1, 1)
campaign_1 <- write_campaign_outputs(c1, 1)

# Campaign 2: retain points 19 and 40; FTIR2 positions 1:3 map to 49, 51, 52.
c2_crds_files <- file.path(
  workflow, "clean_data", "2_campaign",
  c(
    "20241116.000000_20241231.235959_CRDS8_240s_step_average.csv",
    "20241116.000000_20241231.235959_CRDS9_240s_step_average.csv"
  )
)
c2_crds <- read_crds_step_average(c2_crds_files)
# These accepted step-average files pre-date this unified workflow and contain
# the former mapped sampling point instead of MPVPosition. Because the former
# mapping was one-to-one, invert it to recover MPVPosition before applying the
# revised Campaign 2 metadata mapping.
legacy_c2_points <- list(
  CRDS9 = c(1, 3, 4, 6, 7, 9, 10, 12, 13, 15, 16, 18, 21, 22, 24, 25),
  CRDS8 = c(27, 28, 30, 31, 33, 34, 36, 37, 39, 42, 43, 45, 46, 48, 49, 51)
)
for (a in names(legacy_c2_points)) {
  inverse <- setNames(as.character(1:16), as.character(legacy_c2_points[[a]]))
  rows <- c2_crds$analyser == a
  c2_crds[rows, source.position := unname(inverse[source.position])]
}
if (anyNA(c2_crds$source.position)) {
  stop("Campaign 2 legacy sampling points could not all be mapped to MPVPosition.")
}
c2_ftir <- filter_period(
  ftir_all[analyser == "FTIR2"],
  "2024-11-16 00:00:00", "2024-12-31 23:59:59"
)
c2 <- rbindlist(list(c2_crds, c2_ftir), use.names = TRUE, fill = TRUE)
c2 <- apply_sampling_map(c2, 2)
campaign_2 <- write_campaign_outputs(c2, 2)

# Campaign 3: supplied step-average files already retain physical locations.
c3_files <- list.files(
  file.path(workflow, "clean_data", "3_campaign"),
  pattern = "ATB_.*CRDS[89][.]csv$", full.names = TRUE
)
c3 <- read_crds_step_average(c3_files)
c3 <- apply_sampling_map(c3, 3, preserve_existing_location = TRUE)
campaign_3 <- write_campaign_outputs(c3, 3)

# Campaign 4: MPVPosition is mapped through the campaign metadata table.
c4_files <- list.files(
  file.path(workflow, "clean_data", "4_campaign"),
  pattern = "ATB_.*CRDS8[.]csv$", full.names = TRUE
)
c4 <- read_crds_step_average(c4_files)
c4 <- apply_sampling_map(c4, 4)
campaign_4 <- write_campaign_outputs(c4, 4)

all_campaigns <- rbindlist(
  list(campaign_1, campaign_2, campaign_3, campaign_4),
  use.names = TRUE, fill = TRUE
)
setorder(all_campaigns, DATE.TIME, campaign, analyser, sampling.point)
fwrite(
  all_campaigns,
  file.path(output_root, "all_campaigns_combined.csv"),
  dateTimeAs = "write.csv"
)

validation <- all_campaigns[, .(
  first = min(DATE.TIME, na.rm = TRUE), last = max(DATE.TIME, na.rm = TRUE),
  observations = .N,
  sampling.points = uniqueN(sampling.point)
), by = .(campaign, analyser)]
fwrite(validation, file.path(output_root, "preparation_validation_summary.csv"),
       dateTimeAs = "write.csv")

message("Prepared Campaigns 1--4 in: ", output_root)
