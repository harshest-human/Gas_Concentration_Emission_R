###############################################################################
# Sensitivity preparation: 180 s CRDS flush + final 60 s average
#
# Purpose
#   1. Re-read CRDS8 and CRDS9 per-second raw data for Campaigns 1--4.
#   2. Identify ordered multiplexer cycles and individual valve dwells.
#   3. Flush elapsed seconds 0--<180 and average seconds 180--<240 only.
#   4. Timestamp each retained CRDS value at the final whole minute of the
#      four-minute window.
#   5. Harmonise CRDS and FTIR records into one campaign-wise data frame.
#
# This is a sensitivity workflow. It does not overwrite accepted clean data.
# Output data are written to:
#   clean_data/crds_180s_flush_60s_sensitivity_v01
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(lubridate)
})

workflow <- paste0(
  "D:/Data_Analysis_R/Gas_Concentration_Emission_R/",
  "workflows/mapping_high_resolution"
)
metadata_dir <- file.path(workflow, "meta_data")
output_root <- file.path(
  workflow, "clean_data", "crds_180s_flush_60s_sensitivity_v01"
)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

tz_local <- "Europe/Berlin"
gas_columns <- c("CO2", "CH4", "NH3", "H2O")
flush_seconds <- 180L
average_end_seconds <- 240L
minimum_average_rows <- 30L
maximum_gap_seconds <- 10

crds_direct_roots <- c(
  CRDS8 = paste0(
    "D:/Data_Analysis_R/CRDS_raw_recovered/CRDS8_raw_recovered/",
    "Picarro_G2508/CRDS08_raw"
  ),
  CRDS9 = paste0(
    "D:/Data_Analysis_R/CRDS_raw_recovered/CRDS9_raw_recovered/",
    "Picarro_G2509/CRDS09_raw"
  )
)
crds_zip_roots <- c(
  CRDS8 = paste0(
    "D:/Data_Analysis_R/CRDS_raw_recovered/",
    "Picarro_G2508_backup_all"
  ),
  CRDS9 = paste0(
    "D:/Data_Analysis_R/CRDS_raw_recovered/",
    "Picarro_G2509_backup_all"
  )
)

campaign_periods <- list(
  `1` = list(
    c("2024-06-01 00:00:00", "2024-08-31 23:59:59"),
    c("2024-10-01 00:00:00", "2024-10-24 23:59:59")
  ),
  `2` = list(c("2024-11-16 00:00:00", "2024-12-31 23:59:59")),
  `3` = list(c("2025-08-19 12:00:00", "2025-11-19 12:00:00")),
  `4` = list(c("2025-12-09 12:00:00", "2026-07-21 03:00:00"))
)

expected_positions <- list(
  `1` = list(CRDS8 = 1:16, CRDS9 = 1:16),
  `2` = list(CRDS8 = 1:16, CRDS9 = 1:16),
  `3` = list(CRDS8 = 1:16, CRDS9 = 1:16),
  `4` = list(CRDS8 = 1:9)
)

###############################################################################
##### General helpers
###############################################################################

parse_local_time <- function(x) {
  as.POSIXct(x, format = "%Y-%m-%d %H:%M:%OS", tz = tz_local)
}

parse_filename_time <- function(path) {
  value <- sub(
    ".*-(20[0-9]{6})-([0-9]{6})-DataLog.*", "\\1\\2", basename(path)
  )
  bad <- !grepl("^[0-9]{14}$", value)
  value[bad] <- NA_character_
  as.POSIXct(value, format = "%Y%m%d%H%M%S", tz = tz_local)
}

parse_backup_time <- function(path) {
  value <- sub(
    ".*Backup_(20[0-9]{6})_([0-9]{6})[.]zip$", "\\1\\2", basename(path)
  )
  bad <- !grepl("^[0-9]{14}$", value)
  value[bad] <- NA_character_
  as.POSIXct(value, format = "%Y%m%d%H%M%S", tz = tz_local)
}

in_periods <- function(time, periods, buffer_seconds = 0) {
  keep <- rep(FALSE, length(time))
  for (period in periods) {
    start <- as.POSIXct(period[1], tz = tz_local) - buffer_seconds
    end <- as.POSIXct(period[2], tz = tz_local) + buffer_seconds
    keep <- keep | (!is.na(time) & time >= start & time <= end)
  }
  keep
}

filter_exact_periods <- function(x, periods) {
  keep <- rep(FALSE, nrow(x))
  for (period in periods) {
    start <- as.POSIXct(period[1], tz = tz_local)
    end <- as.POSIXct(period[2], tz = tz_local)
    keep <- keep | (x$DATE.TIME >= start & x$DATE.TIME <= end)
  }
  x[keep]
}

read_mapping <- function(campaign) {
  path <- file.path(
    metadata_dir,
    paste0("campaign", campaign, "_sampling_point_analyser_mapping.csv")
  )
  map <- fread(path, colClasses = "character")
  setnames(map, "Messstelle_or_MPVPosition", "source.position")
  map[, analyser := toupper(analyser)]
  map
}

normalise_sampling_point <- function(x) {
  x <- trimws(as.character(x))
  x[tolower(x) == "in"] <- "ring_in"
  x[tolower(x) == "s"] <- "s"
  x[x == "52"] <- "s"
  x
}

map_sampling_points <- function(x, campaign) {
  map <- read_mapping(campaign)

  if (campaign == 3L) {
    # CRDS9 positions 15 and 16 changed assignment on 28 August 2025.
    map <- map[!(analyser == "CRDS9" & source.position %in% c("15", "16"))]
  }

  x[, sampling.point := NA_character_]
  x[map, on = .(analyser, source.position),
    sampling.point := i.sampling.point]

  if (campaign == 3L) {
    switch_time <- as.POSIXct("2025-08-28 12:00:00", tz = tz_local)
    x[analyser == "CRDS9" & source.position == "15" & DATE.TIME < switch_time,
      sampling.point := "25"]
    x[analyser == "CRDS9" & source.position == "16" & DATE.TIME < switch_time,
      sampling.point := "27"]
    x[analyser == "CRDS9" & source.position == "15" & DATE.TIME >= switch_time,
      sampling.point := "ring_in"]
    x[analyser == "CRDS9" & source.position == "16" & DATE.TIME >= switch_time,
      sampling.point := "s"]
  }

  x[, sampling.point := normalise_sampling_point(sampling.point)]
  if (anyNA(x$sampling.point)) {
    bad <- unique(x[is.na(sampling.point), .(analyser, source.position)])
    stop(
      "Unmapped positions in Campaign ", campaign, ": ",
      paste(paste(bad$analyser, bad$source.position, sep = ":"),
            collapse = ", ")
    )
  }
  x
}

###############################################################################
##### Raw CRDS readers
###############################################################################

read_crds_dat <- function(paths, analyser) {
  if (!length(paths)) return(data.table())
  keep <- c("DATE", "TIME", "MPVPosition", gas_columns)

  parts <- lapply(paths, function(path) {
    tryCatch(
      fread(path, select = keep, showProgress = FALSE),
      error = function(e) {
        warning("Skipped unreadable raw file: ", path, " [", conditionMessage(e), "]")
        NULL
      }
    )
  })
  parts <- Filter(Negate(is.null), parts)
  if (!length(parts)) return(data.table())

  x <- rbindlist(parts, use.names = TRUE, fill = TRUE)
  x[, DATE.TIME := parse_local_time(paste(DATE, TIME))]
  x[, c("DATE", "TIME") := NULL]
  x <- x[!is.na(DATE.TIME) & !is.na(MPVPosition) & is.finite(MPVPosition)]

  # The analogue selector signal briefly ramps between positions. Retain only
  # settled values within 0.02 of an integer selector position.
  x[, position.integer := as.integer(round(MPVPosition))]
  x <- x[
    position.integer >= 1L & position.integer <= 16L &
      abs(MPVPosition - position.integer) <= 0.02
  ]
  x[, `:=`(
    MPVPosition = NULL,
    source.position = as.character(position.integer),
    position.integer = NULL,
    analyser = analyser
  )]
  for (g in gas_columns) set(x, j = g, value = as.numeric(x[[g]]))
  setorder(x, DATE.TIME)
  unique(x, by = c("DATE.TIME", "source.position", gas_columns))
}

summarise_complete_cycles <- function(raw, analyser, positions) {
  if (!nrow(raw)) return(list(averages = data.table(), audit = data.table()))
  setorder(raw, DATE.TIME)

  raw[, gap.seconds := c(Inf, diff(as.numeric(DATE.TIME)))]
  raw[, new.dwell := .I == 1L |
        source.position != shift(source.position, fill = first(source.position)) |
        gap.seconds > maximum_gap_seconds]
  raw[, dwell.id := cumsum(new.dwell)]
  raw[, elapsed.seconds := as.numeric(
    difftime(DATE.TIME, min(DATE.TIME), units = "secs")
  ), by = dwell.id]

  dwell <- raw[, .(
    dwell.start = min(DATE.TIME),
    dwell.end = max(DATE.TIME),
    observed.seconds = as.numeric(
      difftime(max(DATE.TIME), min(DATE.TIME), units = "secs")
    ),
    source.position = first(source.position),
    n.raw.total = .N
  ), by = dwell.id]
  setorder(dwell, dwell.start)
  dwell[, position.integer := as.integer(source.position)]
  dwell[, cycle.id := cumsum(
    position.integer == min(positions) &
      shift(position.integer, fill = 0L) != min(positions)
  )]
  dwell[, transition.ok := position.integer == fifelse(
    shift(position.integer, fill = max(positions)) == max(positions),
    min(positions), shift(position.integer, fill = 0L) + 1L
  )]

  cycle <- dwell[cycle.id > 0, .(
    cycle.start = min(dwell.start),
    cycle.end = max(dwell.end),
    n.steps = .N,
    n.positions = uniqueN(position.integer),
    positions.complete = setequal(unique(position.integer), positions),
    order.score = mean(transition.ok),
    complete.window.fraction = mean(observed.seconds >= 230)
  ), by = cycle.id]
  cycle[, cyclical.experiment :=
          positions.complete & order.score >= 0.75 &
          complete.window.fraction >= 0.75]

  dwell <- cycle[dwell, on = "cycle.id"]
  raw <- dwell[raw, on = "dwell.id"]
  selected <- raw[
    cyclical.experiment == TRUE &
      observed.seconds >= 230 &
      elapsed.seconds >= flush_seconds &
      elapsed.seconds < average_end_seconds
  ]

  averages <- selected[, c(list(
    DATE.TIME = as.POSIXct(
      floor(as.numeric(max(DATE.TIME)) / 60) * 60,
      origin = "1970-01-01", tz = tz_local
    ),
    analyser = analyser,
    source.position = first(source.position),
    dwell.seconds.observed = first(observed.seconds),
    flush.seconds = flush_seconds,
    averaging.seconds = average_end_seconds - flush_seconds,
    n.raw = .N,
    cycle.id = first(cycle.id)
  ), lapply(.SD, mean, na.rm = TRUE)),
  by = dwell.id, .SDcols = gas_columns]
  averages <- averages[n.raw >= minimum_average_rows]
  averages[, dwell.id := NULL]

  # Picarro reports NH3 in ppb. The harmonised data use ppm.
  averages[, NH3 := NH3 / 1000]

  audit <- dwell[, .(
    analyser = analyser,
    cycles = uniqueN(cycle.id[cycle.id > 0]),
    valid.cycles = uniqueN(cycle.id[cyclical.experiment == TRUE]),
    dwells = .N,
    complete.240s.windows = sum(observed.seconds >= 230),
    retained.averages = nrow(averages)
  )]
  list(averages = averages, audit = audit)
}

process_chunk_stream <- function(chunks, analyser, positions, periods,
                                 label) {
  carry <- data.table()
  all.averages <- list()
  chunk.audit <- list()

  for (i in seq_along(chunks)) {
    message(label, ": chunk ", i, "/", length(chunks))
    current <- chunks[[i]]()
    current <- filter_exact_periods(current, periods)
    combined <- rbindlist(list(carry, current), use.names = TRUE, fill = TRUE)
    if (!nrow(combined)) next
    setorder(combined, DATE.TIME)
    combined <- unique(
      combined, by = c("DATE.TIME", "source.position", gas_columns)
    )

    position.integer <- as.integer(combined$source.position)
    cycle.starts <- which(
      position.integer == min(positions) &
        shift(position.integer, fill = 0L) != min(positions)
    )

    is.final <- i == length(chunks)
    if (!is.final && length(cycle.starts)) {
      last.start <- cycle.starts[length(cycle.starts)]
      process.raw <- if (last.start > 1L) combined[seq_len(last.start - 1L)] else data.table()
      carry <- combined[last.start:nrow(combined)]
    } else if (!is.final) {
      process.raw <- data.table()
      carry <- combined
    } else {
      process.raw <- combined
      carry <- data.table()
    }

    result <- summarise_complete_cycles(process.raw, analyser, positions)
    if (nrow(result$averages)) all.averages[[length(all.averages) + 1L]] <- result$averages
    if (nrow(result$audit)) {
      result$audit[, chunk := i]
      chunk.audit[[length(chunk.audit) + 1L]] <- result$audit
    }
    rm(current, combined, process.raw, result)
    invisible(gc())
  }

  list(
    averages = unique(
      rbindlist(all.averages, use.names = TRUE, fill = TRUE),
      by = c("DATE.TIME", "analyser", "source.position")
    ),
    audit = rbindlist(chunk.audit, use.names = TRUE, fill = TRUE)
  )
}

make_direct_chunks <- function(root, analyser, periods, chunk_size = 24L) {
  files <- list.files(root, pattern = "[.]dat$", recursive = TRUE,
                      full.names = TRUE, ignore.case = TRUE)
  times <- parse_filename_time(files)
  files <- files[in_periods(times, periods, buffer_seconds = 7200)]
  times <- parse_filename_time(files)
  files <- files[order(times, basename(files))]
  files <- files[!duplicated(basename(files))]
  groups <- split(files, ceiling(seq_along(files) / chunk_size))
  lapply(groups, function(paths) {
    force(paths)
    function() read_crds_dat(paths, analyser)
  })
}

make_zip_chunks <- function(root, analyser, periods) {
  archives <- list.files(root, pattern = "[.]zip$", recursive = TRUE,
                         full.names = TRUE, ignore.case = TRUE)
  times <- parse_backup_time(archives)
  archives <- archives[in_periods(times, periods, buffer_seconds = 24 * 3600)]
  times <- parse_backup_time(archives)
  archives <- archives[order(times, basename(archives))]

  seen.members <- new.env(parent = emptyenv())
  lapply(archives, function(archive) {
    force(archive)
    function() {
      listing <- tryCatch(unzip(archive, list = TRUE), error = function(e) NULL)
      if (is.null(listing) || !nrow(listing)) return(data.table())
      members <- listing$Name[grepl("[.]dat$", listing$Name, ignore.case = TRUE)]
      if (!length(members)) return(data.table())
      base <- basename(members)
      unseen <- !vapply(base, exists, logical(1), envir = seen.members,
                        inherits = FALSE)
      members <- members[unseen]
      base <- base[unseen]
      if (!length(members)) return(data.table())
      for (name in base) assign(name, TRUE, envir = seen.members)

      member.times <- parse_filename_time(base)
      members <- members[in_periods(member.times, periods, buffer_seconds = 7200)]
      if (!length(members)) return(data.table())

      temp.root <- tempfile(paste0("crds_", analyser, "_"))
      dir.create(temp.root, recursive = TRUE)
      on.exit(unlink(temp.root, recursive = TRUE, force = TRUE), add = TRUE)
      tryCatch(
        unzip(archive, files = members, exdir = temp.root),
        error = function(e) character()
      )
      extracted <- list.files(
        temp.root, pattern = "[.]dat$", recursive = TRUE, full.names = TRUE,
        ignore.case = TRUE
      )
      read_crds_dat(extracted, analyser)
    }
  })
}

process_direct <- function(campaign, analyser, periods) {
  chunks <- make_direct_chunks(crds_direct_roots[analyser], analyser, periods)
  process_chunk_stream(
    chunks, analyser, expected_positions[[as.character(campaign)]][[analyser]],
    periods, paste0("Campaign ", campaign, " ", analyser, " direct")
  )
}

process_zip <- function(campaign, analyser, periods) {
  chunks <- make_zip_chunks(crds_zip_roots[analyser], analyser, periods)
  process_chunk_stream(
    chunks, analyser, expected_positions[[as.character(campaign)]][[analyser]],
    periods, paste0("Campaign ", campaign, " ", analyser, " zip")
  )
}

###############################################################################
##### FTIR reader and instrument harmonisation
###############################################################################

read_ftir_files <- function() {
  files <- list.files(
    file.path(workflow, "raw_data", "FTIR_raw"),
    pattern = "FTIR[12].*[.]TXT$", full.names = TRUE, ignore.case = TRUE
  )
  parts <- lapply(files, function(path) {
    analyser <- toupper(regmatches(
      basename(path), regexpr("FTIR[12]", basename(path), ignore.case = TRUE)
    ))
    raw <- fread(path, sep = "\t", encoding = "Latin-1", fill = TRUE,
                 check.names = TRUE, showProgress = FALSE)
    date.text <- paste(raw$Datum, raw$Zeit)
    date.time <- as.POSIXct(
      date.text, format = "%Y-%m-%d %H:%M:%OS", tz = tz_local
    )
    missing.time <- is.na(date.time)
    date.time[missing.time] <- as.POSIXct(
      date.text[missing.time], format = "%d.%m.%Y %H:%M:%OS",
      tz = tz_local
    )
    out <- data.table(
      DATE.TIME = date.time,
      analyser = analyser,
      source.position = trimws(as.character(raw$Messstelle)),
      flush.seconds = NA_integer_,
      averaging.seconds = NA_integer_,
      n.raw = NA_integer_,
      cycle.id = NA_integer_,
      dwell.seconds.observed = NA_real_
    )
    for (g in gas_columns) {
      out[, (g) := suppressWarnings(as.numeric(
        gsub(",", ".", trimws(as.character(raw[[g]])), fixed = TRUE)
      ))]
    }
    out[!is.na(DATE.TIME) & nzchar(source.position)]
  })
  rbindlist(parts, use.names = TRUE, fill = TRUE)
}

add_scale_columns <- function(x) {
  ftir_divisors <- c(CO2 = 1.06, CH4 = 1.06, NH3 = 1.09, H2O = 1.00)
  x[, correction.applied := analyser %chin% c("FTIR1", "FTIR2") & campaign == 1L]
  for (g in gas_columns) {
    raw.name <- paste0(g, "_raw")
    corr.name <- paste0(g, "_corr")
    x[, (raw.name) := get(g)]
    x[, (corr.name) := fifelse(
      correction.applied,
      get(g) / unname(ftir_divisors[g]),
      get(g)
    )]
    x[, (g) := get(corr.name)]
  }
  x
}

###############################################################################
##### Executive processing
###############################################################################

results <- list()
audits <- list()

# Campaigns 1 and 2 are fully available as recovered direct DAT files.
for (campaign in 1:2) {
  periods <- campaign_periods[[as.character(campaign)]]
  for (analyser in c("CRDS8", "CRDS9")) {
    result <- process_direct(campaign, analyser, periods)
    result$averages[, campaign := campaign]
    result$audit[, campaign := campaign]
    results[[paste(campaign, analyser, sep = "_")]] <- result$averages
    audits[[paste(campaign, analyser, sep = "_")]] <- result$audit
  }
}

# Campaign 3 uses direct recovered files until 29 September 2025 and backup
# archives thereafter. The periods do not overlap, preventing duplicate cycles.
c3.direct.period <- list(c("2025-08-19 12:00:00", "2025-09-29 23:59:59"))
c3.zip.period <- list(c("2025-09-30 00:00:00", "2025-11-19 12:00:00"))
for (analyser in c("CRDS8", "CRDS9")) {
  direct <- process_direct(3L, analyser, c3.direct.period)
  zipped <- process_zip(3L, analyser, c3.zip.period)
  combined <- unique(
    rbindlist(list(direct$averages, zipped$averages), use.names = TRUE,
              fill = TRUE),
    by = c("DATE.TIME", "analyser", "source.position")
  )
  combined[, campaign := 3L]
  audit <- rbindlist(list(direct$audit, zipped$audit), use.names = TRUE,
                     fill = TRUE)
  audit[, campaign := 3L]
  results[[paste(3, analyser, sep = "_")]] <- combined
  audits[[paste(3, analyser, sep = "_")]] <- audit
}

# Campaign 4 uses CRDS8 only and is recovered from backup archives.
c4 <- process_zip(4L, "CRDS8", campaign_periods[["4"]])
c4$averages[, campaign := 4L]
c4$audit[, campaign := 4L]
results[["4_CRDS8"]] <- c4$averages
audits[["4_CRDS8"]] <- c4$audit

crds <- rbindlist(results, use.names = TRUE, fill = TRUE)
setorder(crds, campaign, analyser, DATE.TIME)

# Add FTIR observations to Campaigns 1 and 2 without altering their native
# timestamps. FTIR records are already selector-level results.
ftir <- read_ftir_files()
ftir.campaigns <- list()
for (campaign in 1:2) {
  z <- filter_exact_periods(
    copy(ftir), campaign_periods[[as.character(campaign)]]
  )
  if (campaign == 2L) z <- z[analyser == "FTIR2"]
  z[, campaign := campaign]
  ftir.campaigns[[as.character(campaign)]] <- z
}
ftir <- rbindlist(ftir.campaigns, use.names = TRUE, fill = TRUE)

combined <- rbindlist(list(crds, ftir), use.names = TRUE, fill = TRUE)
combined <- combined[!is.na(DATE.TIME)]
combined[, source.variable := fifelse(
  grepl("^CRDS", analyser), "MPVPosition", "Messstelle"
)]

mapped <- list()
for (campaign.id in sort(unique(combined$campaign))) {
  mapped[[as.character(campaign.id)]] <- map_sampling_points(
    copy(combined[campaign == campaign.id]), campaign.id
  )
}
combined <- rbindlist(mapped, use.names = TRUE, fill = TRUE)
combined <- add_scale_columns(combined)
setorder(combined, campaign, DATE.TIME, analyser, sampling.point)
combined <- unique(
  combined, by = c("campaign", "DATE.TIME", "analyser", "sampling.point")
)

column.order <- c(
  "DATE.TIME", "campaign", "analyser", "source.variable",
  "source.position", "sampling.point",
  paste0(gas_columns, "_raw"), paste0(gas_columns, "_corr"), gas_columns,
  "correction.applied", "flush.seconds", "averaging.seconds", "n.raw",
  "dwell.seconds.observed", "cycle.id"
)
setcolorder(combined, column.order)

for (campaign.id in 1:4) {
  campaign.dir <- file.path(output_root, paste0("campaign_", campaign.id))
  dir.create(campaign.dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(
    combined[campaign == campaign.id],
    file.path(campaign.dir, paste0("campaign_", campaign.id, "_combined.csv")),
    dateTimeAs = "write.csv"
  )
}
fwrite(
  combined, file.path(output_root, "all_campaigns_combined.csv"),
  dateTimeAs = "write.csv"
)

audit <- rbindlist(audits, use.names = TRUE, fill = TRUE)
fwrite(audit, file.path(output_root, "crds_chunk_cycle_audit.csv"),
       dateTimeAs = "write.csv")

summary <- combined[, .(
  first = min(DATE.TIME),
  last = max(DATE.TIME),
  observations = .N,
  sampling.points = uniqueN(sampling.point),
  median.raw.per.average = as.numeric(median(n.raw, na.rm = TRUE))
), by = .(campaign, analyser)]
summary[!is.finite(median.raw.per.average), median.raw.per.average := NA_real_]
fwrite(summary, file.path(output_root, "preparation_summary.csv"),
       dateTimeAs = "write.csv")

message("Sensitivity data written to: ", output_root)
