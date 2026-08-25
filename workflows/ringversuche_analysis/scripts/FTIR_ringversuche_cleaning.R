# =============================================================================
# FTIR_ringversuche_cleaning.R
#
# Cleans Gasmet CX4000 FTIR raw outputs for the Ringversuche
# 2025-04-08..2025-04-15 window and writes V5-format clean data that
# Ringversuche_analysis_script.R can consume.
#
# Outputs per analyzer:
#   clean_data/20250408-15_concentrations_7.5minute_data/
#     20250408-15_7.5_avg_<lab>_<analyzer>.csv  - cycle-level (7.5 min)
#   clean_data/Version_9/
#   20250408-15_long_<lab>_<analyzer>.csv     - hourly, one row per (hour, location)
#   20250408-15_wide_<lab>_<analyzer>.csv     - hourly, pivoted by location
#
# Labs handled (different raw layouts):
#   ATB  FTIR.1  - English headers (Line/Date/Time), ppm + vol-%; the source
#                  TXT contains multiple Gasmet exports concatenated, so we
#                  read with fill=Inf and drop the repeated header rows.
#                  Line: 1=N, 2=in, 3=S
#   LUFA FTIR.2  - German headers (Messstelle/Datum/Zeit), ppm + vol-%; uses
#                  Messstelle as the line identifier (1=in, 2=S, 3=N).
#                  No unit conversion required.
#   MBBM FTIR.3  - English headers, comma-as-decimal; uses time-cycle
#                  inference (no reliable Line column).
#   ANECO FTIR.4 - 8 daily German RESULTS_DDMMYY.TXT files in mixed units
#                  (g/m3 + vol-% + mg/m3); uses time-cycle inference.
#
# Time-cycle labs (MBBM, ANECO): location_cycle = c("in","N","in","S"),
# 450 s per step, 180 s flush, 270 s averaged.
# =============================================================================

# ---- libraries --------------------------------------------------------------
suppressPackageStartupMessages({
        library(tidyverse)
        library(lubridate)
        library(data.table)
        library(readr)
})

# ---- paths & helpers --------------------------------------------------------
proj_root  <- "D:/Data_Analysis_R/Gas_Concentration_Emission_R/workflows/ringversuche_analysis"
utils_dir  <- "D:/Data_Analysis_R/Gas_Concentration_Emission_R/workflows/utils"

source(file.path(utils_dir, "round to interval function.R"))     # round_to_interval()
source(file.path(utils_dir, "remove_outliers_function.R"))        # remove_outliers()

# ---- config -----------------------------------------------------------------
start_time <- ymd_hms("2025-04-08 12:00:00", tz = "UTC")
end_time   <- ymd_hms("2025-04-14 13:00:00", tz = "UTC")

flush_sec      <- 180
interval_sec   <- 450
location_cycle <- c("in", "N", "in", "S")    # 4 * 450 s = 30 min super-cycle

out_version  <- "Version_9"
out_dir      <- file.path(proj_root, "clean_data", out_version)
interval_out_dir <- file.path(
        proj_root, "clean_data", "20250408-15_concentrations_7.5minute_data"
)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(interval_out_dir, showWarnings = FALSE, recursive = TRUE)
file_stub    <- "20250408-15"

# ---- raw readers ------------------------------------------------------------

#' Read an English-header Gasmet RESULTS.TXT (Line/Date/Time/H2O/CO2/NH3/CH4/...).
#' Handles the case where multiple exports have been concatenated into one file
#' (the second header row gets re-treated as data; we filter it out below).
read_ftir_results_english <- function(path, gases = c("CO2","CH4","NH3","H2O","N2O")) {
        d <- fread(path, header = TRUE, fill = Inf, na.strings = c("","NA"))
        # Pick only the columns we care about (by exact name).
        keep <- intersect(c("Line","Date","Time", gases), names(d))
        d <- d[, ..keep]
        # Drop header-repeat rows where Line/Date are the literal header strings.
        d <- d[grepl("^\\d+$", as.character(Line)) & grepl("^\\d{4}-\\d{2}-\\d{2}$", as.character(Date))]
        # Coerce types
        d[, Line := as.character(Line)]
        d[, DATE.TIME := as.POSIXct(paste(Date, Time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")]
        d[, c("Date","Time") := NULL]
        for (g in intersect(gases, names(d))) d[, (g) := as.numeric(get(g))]
        as_tibble(d)
}

#' Read a German-header LUFA-style TXT (Messstelle/Datum/Zeit, ppm + vol-%).
#' The file has duplicate "Messstelle" headers (one is the line identifier,
#' another is a redundant readback) - we only use the first (column 1).
read_ftir_results_german_lufa <- function(path, gases = c("CO2","CH4","NH3","H2O","N2O")) {
        d <- fread(path, header = TRUE, fill = Inf, encoding = "Latin-1",
                   na.strings = c("","NA"), check.names = FALSE)
        # The file has Messstelle in col 1 *and* somewhere later; we want col 1.
        # Rename the FIRST Messstelle to "Line" by position.
        ms_idx <- which(names(d) == "Messstelle")[1]
        if (length(ms_idx) == 0) stop("No 'Messstelle' column in ", basename(path))
        setnames(d, old = ms_idx, new = "Line")
        # Drop any later Messstelle columns to avoid duplicate-name confusion
        dupes <- which(names(d) == "Messstelle")
        if (length(dupes) > 0) d[, (names(d)[dupes]) := NULL]
        # Now rename Datum/Zeit
        setnames(d, c("Datum","Zeit"), c("Date","Time"), skip_absent = TRUE)
        # Pick gas columns by exact name
        keep <- intersect(c("Line","Date","Time", gases), names(d))
        d <- d[, ..keep]
        d <- d[grepl("^\\d+$", as.character(Line))]
        d[, Line := as.character(Line)]
        d[, DATE.TIME := as.POSIXct(paste(Date, Time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")]
        d[, c("Date","Time") := NULL]
        for (g in intersect(gases, names(d))) d[, (g) := as.numeric(get(g))]
        as_tibble(d)
}

#' Read an MBBM-style TXT (English headers, comma-as-decimal numbers).
read_ftir_results_mbbm <- function(path, gases = c("CO2","CH4","NH3","H2O","N2O")) {
        d <- read.delim(path, stringsAsFactors = FALSE, fill = TRUE, check.names = FALSE)
        keep <- intersect(c("Date","Time", gases), names(d))
        d <- d[, keep, drop = FALSE]
        d <- d[grepl("^\\d{4}-\\d{2}-\\d{2}$", d$Date), , drop = FALSE]
        d$DATE.TIME <- as.POSIXct(paste(d$Date, d$Time),
                                  format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
        d$Date <- NULL; d$Time <- NULL
        for (g in intersect(gases, names(d))) {
                d[[g]] <- as.numeric(gsub(",", ".", d[[g]]))
        }
        as_tibble(d) %>% select(DATE.TIME, all_of(intersect(gases, names(d))))
}

#' Read all 8 ANECO daily RESULTS_DDMMYY.TXT files; return a long tibble
#' with DATE.TIME plus mixed-unit raw columns.
read_aneco_daily_results <- function(dir) {
        files <- list.files(dir, pattern = "^RESULTS_\\d{6}\\.TXT$", full.names = TRUE)
        if (length(files) == 0) stop("No ANECO RESULTS files in ", dir)
        purrr::map_dfr(files, function(p) {
                dt <- fread(p, sep = "\t", header = TRUE, fill = Inf,
                            encoding = "Latin-1", na.strings = c("","NA"))
                need <- c("Datum","Zeit","Water vapor H2O",
                          "Carbon dioxide CO2","NH3","CH4")
                if (any(!need %in% names(dt))) {
                        stop(sprintf("ANECO %s missing: %s", basename(p),
                                     paste(setdiff(need, names(dt)), collapse=", ")))
                }
                # Drop header-repeat rows
                dt <- dt[grepl("^\\d{4}-\\d{2}-\\d{2}$", as.character(Datum))]
                tibble(
                        DATE.TIME = as.POSIXct(paste(dt$Datum, dt$Zeit),
                                               format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
                        H2O_gm3   = as.numeric(dt[["Water vapor H2O"]]),
                        CO2_vol   = as.numeric(dt[["Carbon dioxide CO2"]]),
                        NH3_mgm3  = as.numeric(dt[["NH3"]]),
                        CH4_mgm3  = as.numeric(dt[["CH4"]])
                )
        })
}

# ---- shared aggregators -----------------------------------------------------

#' Time-cycle averaging: anchor a fixed location_cycle at start_time, drop
#' the first flush_sec of each step, average the rest.
cycle_average_by_time <- function(df_sec, start_time, end_time,
                                  flush_sec = 180, interval_sec = 450,
                                  location_cycle = c("in","N","in","S")) {
        time_seq <- tibble(DATE.TIME = seq(start_time, end_time, by = "1 sec"))
        df_full  <- time_seq %>% left_join(df_sec, by = "DATE.TIME") %>%
                mutate(
                        step_index        = floor(as.numeric(difftime(DATE.TIME, start_time, units = "secs")) / interval_sec),
                        interval_start    = start_time + step_index * interval_sec,
                        seconds_into_step = as.numeric(difftime(DATE.TIME, interval_start, units = "secs")),
                        location          = location_cycle[(step_index %% length(location_cycle)) + 1]
                )
        gas_cols <- intersect(c("CO2","CH4","NH3","H2O","N2O"), names(df_full))
        df_full %>%
                filter(seconds_into_step >= flush_sec & seconds_into_step < interval_sec) %>%
                group_by(interval_start, location) %>%
                summarise(
                        DATE.TIME = max(interval_start) + interval_sec,
                        across(all_of(gas_cols), ~ mean(.x, na.rm = TRUE)),
                        .groups = "drop"
                ) %>%
                select(-interval_start) %>%
                arrange(DATE.TIME)
}

#' Line-based averaging for labs that have a reliable Line/Messstelle column.
cycle_average_by_line <- function(df_sec, line_to_loc) {
        gas_cols <- intersect(c("CO2","CH4","NH3","H2O","N2O"), names(df_sec))
        df_sec %>%
                filter(.data$Line %in% names(line_to_loc)) %>%
                mutate(location = unname(line_to_loc[as.character(.data$Line)])) %>%
                group_by(DATE.TIME, location) %>%
                summarise(across(all_of(gas_cols), ~ mean(.x, na.rm = TRUE)),
                          .groups = "drop")
}

# ---- V5 output formatters ---------------------------------------------------

make_hourly_long_v5 <- function(df_cycle, lab_name, analyzer_name) {
        gas_cols <- intersect(c("CO2","CH4","NH3","H2O"), names(df_cycle))
        if (nrow(df_cycle) == 0) {
                return(tibble(DATE.TIME = character(), location = character(),
                              lab = character(), analyzer = character(),
                              !!!setNames(rep(list(numeric()), length(gas_cols)), gas_cols)))
        }
        df_cycle %>%
                mutate(DATE.TIME = floor_date(DATE.TIME, unit = "hour")) %>%
                group_by(DATE.TIME, location) %>%
                summarise(across(all_of(gas_cols), ~ mean(.x, na.rm = TRUE)),
                          .groups = "drop") %>%
                mutate(
                        lab       = lab_name,
                        analyzer  = analyzer_name,
                        DATE.TIME = format(DATE.TIME, "%d/%m/%Y %H:%M")
                ) %>%
                select(DATE.TIME, location, lab, analyzer, all_of(gas_cols)) %>%
                arrange(DATE.TIME, location)
}

make_hourly_wide_v5 <- function(df_long, lab_name) {
        if (nrow(df_long) == 0) {
                return(tibble(DATE.TIME = character(), analyzer = character()))
        }
        gas_cols <- intersect(c("CO2","CH4","NH3","H2O"), names(df_long))
        df_long %>%
                select(-lab) %>%
                pivot_wider(
                        names_from  = location,
                        values_from = all_of(gas_cols),
                        names_sep   = "_"
                ) %>%
                rename_with(
                        .cols = -any_of(c("DATE.TIME", "analyzer")),
                        .fn   = function(x) paste0(x, "_", lab_name)
                ) %>%
                arrange(DATE.TIME)
}

write_outputs <- function(df_cycle, lab_name, analyzer_name, out_dir,
                          interval_out_dir, file_stub) {
        gas_cols <- intersect(c("CO2","CH4","NH3","H2O","N2O"), names(df_cycle))
        # Snap to 7.5-min grid.
        cycle_df <- df_cycle %>%
                mutate(DATE.TIME = round_to_interval(DATE.TIME, interval_sec = 450))

        # Raw descriptive stats per (location, gas) BEFORE any outlier removal
        # (handover §2.3, "pre-Stage-1 7.5-min data"). Goes into the M&M raw-
        # stats paragraph and feeds Annex Table A1.
        raw_stats <- cycle_df %>%
                pivot_longer(any_of(gas_cols), names_to = "gas", values_to = "value") %>%
                filter(!is.na(value), is.finite(value)) %>%
                group_by(location, gas) %>%
                summarise(n      = n(),
                          mean   = mean(value),
                          median = median(value),
                          min    = min(value),
                          max    = max(value),
                          sd     = sd(value),
                          .groups = "drop") %>%
                mutate(lab = lab_name, analyzer = analyzer_name, .before = 1)
        raw_stats_path <- file.path(out_dir,
                                    sprintf("raw_stats_%s_%s.csv",
                                            lab_name, analyzer_name))
        write_excel_csv(raw_stats, raw_stats_path)
        cat(sprintf("  -> %s  (raw 7.5-min stats, %d rows)\n",
                    basename(raw_stats_path), nrow(raw_stats)))

        # Stage 1: 1.5 x IQR per location on the 7.5-min data; logged to CSV.
        stage1_path <- file.path(out_dir,
                                 sprintf("stage1_dropouts_%s_%s.csv",
                                         lab_name, analyzer_name))
        cat(sprintf("  outlier removal (per location):\n"))
        cycle_df <- remove_outliers(cycle_df,
                                    group_cols     = c("location"),
                                    summary_path   = stage1_path,
                                    analyzer_label = analyzer_name,
                                    lab_label      = lab_name)
        df_cycle <- cycle_df

        cycle_out <- df_cycle %>%
                mutate(lab = lab_name, analyzer = analyzer_name) %>%
                select(DATE.TIME, location, lab, analyzer, all_of(gas_cols))

        cycle_path <- file.path(interval_out_dir,
                                sprintf("%s_7.5_avg_%s_%s.csv",
                                        file_stub, lab_name, analyzer_name))
        write_excel_csv(cycle_out, cycle_path)
        cat(sprintf("  -> %s  (%d rows)\n", basename(cycle_path), nrow(cycle_out)))

        long_df <- make_hourly_long_v5(df_cycle, lab_name, analyzer_name)
        long_path <- file.path(out_dir, sprintf("%s_long_%s_%s.csv",
                                                file_stub, lab_name, analyzer_name))
        write_excel_csv(long_df, long_path)
        cat(sprintf("  -> %s  (%d rows)\n", basename(long_path), nrow(long_df)))

        wide_df <- make_hourly_wide_v5(long_df, lab_name)
        wide_path <- file.path(out_dir, sprintf("%s_wide_%s_%s.csv",
                                                file_stub, lab_name, analyzer_name))
        write_excel_csv(wide_df, wide_path)
        cat(sprintf("  -> %s  (%d rows)\n", basename(wide_path), nrow(wide_df)))
}

filter_window <- function(df, start_time, end_time) {
        df %>% filter(DATE.TIME >= start_time & DATE.TIME <= end_time)
}

# =============================================================================
# 1)  ATB FTIR.1
# =============================================================================
cat("\n============================================================\n")
cat("Processing ATB FTIR.1\n")
cat("============================================================\n")

atb_sec <- read_ftir_results_english(
        file.path(proj_root, "raw_data/ATB/FTIR.1/FTIR.1_raw_recovered/FTIR2_RESULTS.TXT")
) %>% filter_window(start_time, end_time)
cat(sprintf("  raw rows in window: %d\n", nrow(atb_sec)))
cat(sprintf("  Line value counts: %s\n",
            paste(names(table(atb_sec$Line)), table(atb_sec$Line), sep="=", collapse=", ")))

atb_cycle <- cycle_average_by_line(atb_sec,
                                   line_to_loc = c("1"="N", "2"="in", "3"="S"))
write_outputs(atb_cycle, "ATB", "FTIR.1", out_dir, interval_out_dir, file_stub)

# =============================================================================
# 2)  LUFA FTIR.2
# =============================================================================
cat("\n============================================================\n")
cat("Processing LUFA FTIR.2\n")
cat("============================================================\n")

lufa_files <- list.files(file.path(proj_root, "raw_data/LUFA/FTIR.2"),
                         pattern = "Ringversuch.*LUFA\\.TXT$",
                         full.names = TRUE)
if (length(lufa_files) == 0) stop("No LUFA Ringversuch TXT found.")
cat(sprintf("  using: %s\n", basename(lufa_files[1])))

lufa_sec <- read_ftir_results_german_lufa(lufa_files[1]) %>%
        filter_window(start_time, end_time)
cat(sprintf("  raw rows in window: %d\n", nrow(lufa_sec)))
cat(sprintf("  Line value counts: %s\n",
            paste(names(table(lufa_sec$Line)), table(lufa_sec$Line), sep="=", collapse=", ")))

lufa_cycle <- cycle_average_by_line(lufa_sec,
                                    line_to_loc = c("1"="in", "2"="S", "3"="N"))
write_outputs(lufa_cycle, "LUFA", "FTIR.2", out_dir, interval_out_dir, file_stub)

# =============================================================================
# 3)  MBBM FTIR.3
# =============================================================================
cat("\n============================================================\n")
cat("Processing MBBM FTIR.3\n")
cat("============================================================\n")

mbbm_path <- file.path(proj_root,
                       "raw_data/MBBM/FTIR.3/2025-06-02_FTIR_Ringversuche_RESULTS_MBBM.TXT")
mbbm_sec <- read_ftir_results_mbbm(mbbm_path) %>%
        filter_window(start_time, end_time)
cat(sprintf("  raw rows in window: %d\n", nrow(mbbm_sec)))

mbbm_cycle <- cycle_average_by_time(mbbm_sec, start_time, end_time,
                                    flush_sec, interval_sec, location_cycle)
write_outputs(mbbm_cycle, "MBBM", "FTIR.3", out_dir, interval_out_dir, file_stub)

# =============================================================================
# 4)  ANECO FTIR.4 legacy export (intentionally excluded: noFTIR4old)
# =============================================================================
cat("\n============================================================\n")
cat("Processing ANECO FTIR.4\n")
cat("============================================================\n")

aneco_dir <- file.path(proj_root, "raw_data/ANECO/ANECO_roh/ANECO_calcamet_results")
aneco_files <- list.files(aneco_dir, pattern = "^RESULTS_\\d{6}\\.TXT$", full.names = TRUE)

include_ftir4_old <- FALSE
if (include_ftir4_old && length(aneco_files) > 0) {
  aneco_raw <- read_aneco_daily_results(aneco_dir)
  cat(sprintf("  raw rows total: %d\n", nrow(aneco_raw)))

  # Unit conversions
  R_const <- 8.314472
  T_const <- 273.15
  P_const <- 101325
  aneco_sec <- aneco_raw %>%
          mutate(
                  H2O = (H2O_gm3 * 100 * R_const * T_const) / (18.015 * P_const),
                  CO2 = CO2_vol * 10000,
                  NH3 = (NH3_mgm3 * 1000 * R_const * T_const) / (17.031 * P_const),
                  CH4 = (CH4_mgm3 * 1000 * R_const * T_const) / (16.04  * P_const)
          ) %>%
          select(DATE.TIME, CO2, CH4, NH3, H2O) %>%
          filter_window(start_time, end_time)
  cat(sprintf("  raw rows in window: %d\n", nrow(aneco_sec)))

  aneco_cycle <- cycle_average_by_time(aneco_sec, start_time, end_time,
                                       flush_sec, interval_sec, location_cycle)
  write_outputs(aneco_cycle, "ANECO", "FTIR.4_old", out_dir,
                interval_out_dir, file_stub)
} else if (include_ftir4_old) {
  cat("  WARN: No daily ANECO RESULTS files (RESULTS_DDMMYY.TXT) found. Skipping FTIR.4.\n")
} else {
  cat("  skipped by configuration (noFTIR4old).\n")
}

# =============================================================================
# 4b) ANECO FTIR.4 - V2 (recovered from spectrum)
# =============================================================================
cat("\n============================================================\n")
cat("Processing ANECO FTIR.4 - V2 (recovered from spectrum)\n")
cat("============================================================\n")

aneco_v2_path <- file.path(proj_root,
                           "raw_data/ANECO/ANECO_roh/2025-04-08_2025-04-15_RESULTS_ANECO_FTIR.4.TXT")
if (file.exists(aneco_v2_path)) {
  aneco_v2_raw <- fread(aneco_v2_path, sep = "\t", header = TRUE, fill = Inf,
                        encoding = "Latin-1", na.strings = c("","NA"))

  # Extract Date, Time and gas columns (already in ppm)
  aneco_v2_sec <- tibble(
    DATE.TIME = as.POSIXct(paste(aneco_v2_raw$Date, aneco_v2_raw$Time),
                           format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    CO2 = as.numeric(aneco_v2_raw$CO2),
    CH4 = as.numeric(aneco_v2_raw$CH4),
    NH3 = as.numeric(aneco_v2_raw$NH3),
    H2O = as.numeric(aneco_v2_raw$H2O)
  ) %>%
    filter_window(start_time, end_time)

  cat(sprintf("  raw rows in window: %d\n", nrow(aneco_v2_sec)))

  aneco_v2_cycle <- cycle_average_by_time(aneco_v2_sec, start_time, end_time,
                                          flush_sec, interval_sec, location_cycle)
  write_outputs(aneco_v2_cycle, "ANECO", "FTIR.4", out_dir,
                interval_out_dir, file_stub)
} else {
  cat(sprintf("  WARN: file not found: %s\n", basename(aneco_v2_path)))
}

# ---- copy LUFA CRDS.2 V5 -> V6 (no raw .dat) --------------------------------
v5_dir <- file.path(proj_root, "clean_data", "Version_5")
lufa_crds_files <- c(
        sprintf("%s_long_LUFA_CRDS.2.csv", file_stub),
        sprintf("%s_wide_LUFA_CRDS.2.csv", file_stub)
)
cat("\n--- copying LUFA CRDS.2 V5 -> V6 (no raw .dat for LUFA CRDS) ---\n")
for (f in lufa_crds_files) {
        src <- file.path(v5_dir, f)
        dst <- file.path(out_dir, f)
        if (file.exists(src) && !file.exists(dst)) {
                file.copy(src, dst, overwrite = FALSE)
                cat(sprintf("  copied %s\n", f))
        } else if (file.exists(dst)) {
                cat(sprintf("  already present: %s\n", f))
        } else {
                cat(sprintf("  WARN: %s not found in V5\n", f))
        }
}

cat("\nDone. Run analysis with data_version <- \"Version_6\".\n")
