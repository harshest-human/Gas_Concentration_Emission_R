# Recover and clean Campaign 1 CRDS8/CRDS9 measurements.
#
# The script identifies sequential 1:16 MPV cycles, distinguishes nominal
# 180-s and 240-s valve dwells, removes the first 60 s of each dwell, and
# averages the remaining 120 or 180 s. Single-position and irregular periods
# remain in the audit table but are excluded from the campaign data.

library(data.table)

tz_local <- "Europe/Berlin"
args <- commandArgs(trailingOnly = TRUE)
month_tag <- if (length(args)) args[1] else "2024-06"
month_start_date <- as.Date(paste0(month_tag, "-01"))
campaign_start <- as.POSIXct(month_start_date, tz = tz_local)
campaign_end <- as.POSIXct(seq(month_start_date, by = "1 month", length.out = 2)[2],
                           tz = tz_local) - 1
workflow <- "D:/Data_Analysis_R/Gas_Concentration_Emission_R/workflows/mapping_high_resolution"
out_root <- file.path(workflow, "clean_data", "1_campaign",
                      "recovered_crds_june_august_v03", month_tag)
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

sources <- c(
  CRDS8 = "D:/Data_Analysis_R/CRDS_raw_recovered/CRDS8_raw_recovered",
  CRDS9 = "D:/Data_Analysis_R/CRDS_raw_recovered/CRDS9_raw_recovered"
)
maps <- list(
  CRDS8 = setNames(as.character(27:42), as.character(1:16)),
  CRDS9 = setNames(as.character(1:16), as.character(1:16))
)
gas <- c("CO2", "CH4", "NH3", "H2O")

file_time <- function(x) {
  z <- sub(".*-(\\d{8})-(\\d{6})-.*", "\\1\\2", basename(x))
  as.POSIXct(z, "%Y%m%d%H%M%S", tz = tz_local)
}

read_campaign_files <- function(root) {
  f <- list.files(root, "\\.dat$", recursive = TRUE, full.names = TRUE,
                  ignore.case = TRUE)
  ft <- file_time(f)
  # Picarro files may span an hour; retain a two-hour leading buffer.
  f <- f[!is.na(ft) & ft >= campaign_start - 7200 & ft <= campaign_end]
  message("Reading ", length(f), " files from ", root)
  parts <- lapply(f, function(path) {
    tryCatch(fread(path, select = c("DATE", "TIME", "MPVPosition", gas),
                   showProgress = FALSE), error = function(e) NULL)
  })
  rbindlist(Filter(Negate(is.null), parts), fill = TRUE)
}

clean_one <- function(analyser) {
  x <- read_campaign_files(sources[[analyser]])
  x[, DATE.TIME := as.POSIXct(paste(DATE, TIME),
                              "%Y-%m-%d %H:%M:%OS", tz = tz_local)]
  x <- x[DATE.TIME >= campaign_start & DATE.TIME <= campaign_end]
  # The multiplexer analogue signal ramps between integer positions for a few
  # seconds. Keep only settled readings close to an integer valve position.
  x[, position_integer := as.integer(round(MPVPosition))]
  x <- x[position_integer %between% c(1L, 16L) &
           abs(MPVPosition - position_integer) <= 0.02]
  x[, MPVPosition := position_integer]
  x[, position_integer := NULL]
  x[, c("DATE", "TIME") := NULL]
  setorder(x, DATE.TIME)
  x <- unique(x, by = c("DATE.TIME", "MPVPosition", gas))

  # A dwell ends at a valve transition or a gap longer than 10 s.
  x[, gap_s := c(Inf, diff(as.numeric(DATE.TIME)))]
  x[, dwell_id := cumsum(.I == 1L | MPVPosition != shift(MPVPosition,
      fill = first(MPVPosition)) | gap_s > 10)]
  x[, elapsed_s := as.numeric(DATE.TIME - min(DATE.TIME)), by = dwell_id]

  dwell <- x[, .(
    dwell_start = min(DATE.TIME), dwell_end = max(DATE.TIME),
    observed_s = as.numeric(max(DATE.TIME) - min(DATE.TIME)),
    n_raw = .N, MPVPosition = first(MPVPosition)
  ), by = dwell_id]
  setorder(dwell, dwell_start)

  # Choose the closest supported protocol when the observed span is plausible.
  dwell[, nominal_s := fifelse(observed_s >= 150 & observed_s < 210, 180L,
                        fifelse(observed_s >= 210 & observed_s <= 270, 240L, NA_integer_))]
  dwell[, transition_ok := MPVPosition ==
          fifelse(shift(MPVPosition, fill = 16L) == 16L, 1L,
                  shift(MPVPosition, fill = 0L) + 1L)]
  dwell[, cycle_id := cumsum(MPVPosition == 1L & shift(MPVPosition,
      fill = 0L) != 1L)]
  cycle <- dwell[cycle_id > 0, .(
    n_steps = .N, n_positions = uniqueN(MPVPosition),
    order_score = mean(transition_ok),
    valid_dwell_fraction = mean(!is.na(nominal_s)),
    cycle_start = min(dwell_start), cycle_end = max(dwell_end)
  ), by = cycle_id]
  cycle[, cyclical_experiment := n_positions >= 12 & n_steps >= 12 &
          order_score >= 0.75 & valid_dwell_fraction >= 0.75]
  dwell <- cycle[dwell, on = "cycle_id"]
  dwell[is.na(cyclical_experiment), cyclical_experiment := FALSE]
  dwell[, exclusion_reason := fcase(
    !cyclical_experiment, "not_a_complete_ordered_1_to_16_cycle",
    is.na(nominal_s), "unsupported_or_incomplete_dwell",
    default = "included"
  )]

  x <- dwell[x, on = "dwell_id"]
  retained <- x[cyclical_experiment & !is.na(nominal_s) &
                  elapsed_s >= 60 & elapsed_s < nominal_s]
  avg <- retained[, c(list(
    DATE.TIME = as.POSIXct(floor(as.numeric(max(DATE.TIME)) / 60) * 60,
                           origin = "1970-01-01", tz = tz_local),
    analyser = analyser, MPVPosition = first(MPVPosition),
    location = unname(maps[[analyser]][as.character(first(MPVPosition))]),
    dwell_seconds = first(nominal_s), flush_seconds = 60L,
    averaging_seconds = first(nominal_s) - 60L, cycle_id = first(cycle_id),
    n_raw = .N),
    lapply(.SD, mean, na.rm = TRUE)), by = dwell_id, .SDcols = gas]
  # Picarro NH3 is recorded in ppb; manuscript concentrations use ppm.
  avg[, NH3 := NH3 / 1000]
  setcolorder(avg, c("DATE.TIME", "analyser", "MPVPosition", "location",
                     gas, "dwell_seconds", "flush_seconds",
                     "averaging_seconds", "cycle_id", "n_raw", "dwell_id"))
  fwrite(dwell, file.path(out_root, paste0(analyser, "_dwell_cycle_audit.csv")))
  avg
}

clean <- rbindlist(lapply(names(sources), clean_one), fill = TRUE)
setorder(clean, DATE.TIME, analyser, MPVPosition)
fwrite(clean, file.path(out_root,
  paste0("campaign1_CRDS8_CRDS9_recovered_step_average_", month_tag, "_v03.csv")),
       dateTimeAs = "write.csv")

summary <- clean[, .(
  first_time = min(DATE.TIME), last_time = max(DATE.TIME), steps = .N,
  cycles = uniqueN(cycle_id), dwell_180 = sum(dwell_seconds == 180),
  dwell_240 = sum(dwell_seconds == 240)
), by = analyser]
fwrite(summary, file.path(out_root,
  paste0("campaign1_recovery_summary_", month_tag, "_v03.csv")))
print(summary)

# Compare CRDS timing with the two FTIR schedules and create the corrected,
# four-analyser Campaign 1 table used by Manuscript 3.
ftir_files <- file.path(workflow, "clean_data", "1_campaign", c(
  "2024_05_17_FTIR1_clean.csv", "2024_05_17_FTIR2_clean.csv"))
ftir <- rbindlist(lapply(ftir_files, fread), fill = TRUE)
ftir[, DATE.TIME := as.POSIXct(DATE.TIME, "%Y-%m-%d %H:%M:%S", tz = tz_local)]
ftir <- ftir[DATE.TIME >= campaign_start & DATE.TIME <= campaign_end]
ftir[, MPVPosition := as.integer(location)]
ftir[, location := as.character(location)]
ftir[analyser == "FTIR2", location := as.character(16L + MPVPosition)]
ftir[analyser == "FTIR1" & MPVPosition <= 9L,
     location := as.character(42L + MPVPosition)]
ftir[analyser == "FTIR1" & MPVPosition == 10L, location := "s"]

timing_audit <- rbindlist(lapply(split(ftir, by = "analyser"), function(z) {
  setorder(z, DATE.TIME)
  dwell_gaps <- as.numeric(diff(z$DATE.TIME), units = "secs")
  dwell_gaps <- dwell_gaps[dwell_gaps < 600]
  data.table(analyser = first(z$analyser), rows = nrow(z),
             first_time = min(z$DATE.TIME), last_time = max(z$DATE.TIME),
             median_dwell_s = median(dwell_gaps, na.rm = TRUE),
             pct_dwells_near_240_s = 100 * mean(dwell_gaps >= 230 &
                                                 dwell_gaps <= 250, na.rm=TRUE))
}))
timing_audit <- rbind(timing_audit, summary[, .(
  analyser, rows = steps, first_time, last_time,
  median_dwell_s = 240,
  pct_dwells_near_240_s = 100 * dwell_240 / steps
)], fill = TRUE)
fwrite(timing_audit, file.path(out_root,
  paste0("campaign1_FTIR_CRDS_timing_audit_", month_tag, "_v03.csv")),
       dateTimeAs = "write.csv")

crds <- copy(clean)
crds <- crds[, .(DATE.TIME, analyser, MPVPosition, location, CO2, CH4, NH3, H2O)]
all <- rbindlist(list(ftir[, .(DATE.TIME, analyser, MPVPosition, location,
                              CO2, CH4, NH3, H2O)], crds), fill = TRUE)
setorder(all, DATE.TIME, analyser, location)
all[, `:=`(
  CO2_raw = CO2, CH4_raw = CH4, NH3_raw = NH3, H2O_raw = H2O,
  CO2_corr = fifelse(grepl("FTIR", analyser), CO2 / 1.06, CO2),
  CH4_corr = fifelse(grepl("FTIR", analyser), CH4 / 1.06, CH4),
  NH3_corr = fifelse(grepl("FTIR", analyser), NH3 / 1.09, NH3),
  H2O_corr = H2O,
  correction_applied = as.integer(grepl("FTIR", analyser)),
  correction_basis = fifelse(grepl("FTIR", analyser),
    "FTIR harmonised to CRDS: CO2/1.06; CH4/1.06; NH3/1.09",
    "CRDS retained without correction")
)]
all[, c("CO2", "CH4", "NH3", "H2O") := NULL]
setcolorder(all, c("DATE.TIME", "analyser", "MPVPosition", "location",
  "CO2_raw", "CH4_raw", "NH3_raw", "H2O_raw", "CO2_corr", "CH4_corr",
  "NH3_corr", "H2O_corr", "correction_applied", "correction_basis"))
fwrite(all, file.path(out_root,
  paste0("campaign1_recovered_four_analyser_CRDS_scale_corrected_",
         month_tag, "_v03.csv")),
  dateTimeAs = "write.csv")
