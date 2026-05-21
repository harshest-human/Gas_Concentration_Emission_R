# =============================================================================
# CRDS_ringversuche_cleaning.R
#
# Cleans Picarro CRDS .dat files for the Ringversuche 2025-04-08..2025-04-15
# window and writes V5-format clean data that the Ringversuche_analysis_script.R
# can consume.
#
# Outputs per analyzer (in clean_data/Version_6/):
#   20250408-15_7.5_avg_<lab>_<analyzer>.csv  - cycle-level (7.5 min) intermediate
#   20250408-15_long_<lab>_<analyzer>.csv     - hourly, one row per (hour, location)
#   20250408-15_wide_<lab>_<analyzer>.csv     - hourly, pivoted by location
#
# Long format (matches existing Version_5 schema):
#   DATE.TIME (DD/MM/YYYY HH:MM), MPVPosition, location, lab, analyzer,
#   CO2, CH4, NH3, H2O   (gas in ppm except H2O in vol-%)
#
# Wide format:
#   DATE.TIME, analyzer, CO2_in_<lab>, CO2_S_<lab>, CO2_N_<lab>,
#                        CH4_in_<lab>, CH4_S_<lab>, CH4_N_<lab>,
#                        NH3_in_<lab>, NH3_S_<lab>, NH3_N_<lab>,
#                        H2O_in_<lab>, H2O_S_<lab>, H2O_N_<lab>
#
# LUFA CRDS.2 is skipped here (no raw .dat files available) - its existing
# Version_5 long + wide files are auto-copied into Version_6 so the folder is
# complete. LUFA FTIR.2 is NOT skipped - it is handled in the FTIR cleaning
# script.
#
# Cleaning logic (from Picarro_CRDS_data_cleaning_script.R::piclean):
#   - Read every .dat, keep DATE/TIME/MPVPosition + requested gases.
#   - Treat each contiguous run of identical MPVPosition as one "step".
#   - Drop the first <flush> seconds of each step (line still rinsing),
#     then average the next (interval - flush) seconds.
#   - Map MPVPosition -> human-readable location using a per-lab table.
#   - Divide NH3 by 1000 (analyzer reports ppb instead of ppm).
#
# Per-lab cycle parameters confirmed from Picarro_CRDS_data_analysis_script.R:
#   ATB CRDS.1  -  positions 1,2,3  ->  N, in, S
#   UB  CRDS.3  -  positions 8,1,9  ->  N, in, S    [analyzer label CORRECTED;
#                                                    old script had CRDS.2]
#   LUFA CRDS.2 -  positions 3,1,2  ->  N, in, S    [SKIPPED, see above]
#   flush = 180 s, interval = 450 s   (i.e. 7.5-minute cycles)
# =============================================================================

# ---- libraries --------------------------------------------------------------
suppressPackageStartupMessages({
        library(tidyverse)
        library(lubridate)
        library(data.table)
        library(readr)
})

# ---- helpers ----------------------------------------------------------------
proj_root  <- "D:/Data_Analysis_R/Gas_Concentration_Emission_R/workflows/ringversuche_analysis"
helpers_dir_crds <- file.path(proj_root, "Picarro-G2508_CRDS_gas_measurement")
helpers_dir_ftir <- file.path(proj_root, "GasmetCX4000_FTIR_Gas_Measurement")

source(file.path(helpers_dir_crds, "Picarro_CRDS_data_cleaning_script.R"))  # piclean()
source(file.path(helpers_dir_crds, "remove_outliers_function.R"))           # remove_outliers()
source(file.path(helpers_dir_ftir, "round to interval function.R"))         # round_to_interval()

# ---- config -----------------------------------------------------------------
start_time <- "2025-04-08 12:00:00"
end_time   <- "2025-04-14 13:00:00"   # 1 h overrun keeps the last hour intact

flush_sec    <- 180
interval_sec <- 450
gases        <- c("CO2", "CH4", "NH3", "H2O", "N2O")

out_version  <- "Version_9"
out_dir      <- file.path(proj_root, "clean_data", out_version)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

file_stub    <- "20250408-15"   # matches existing V5 naming

# Per-lab config. Edit input_path / MPVPosition.levels here if a future
# deployment re-routes ports.
lab_configs <- list(
        list(
                lab            = "ATB",
                analyzer       = "CRDS.1",
                input_path     = "D:/Data_Analysis_R/CRDS_raw_recovered/Picarro_G2508/CRDS08_raw/2025/04",
                mpv_levels     = c("1", "2", "3"),
                location_levels = c("N", "in", "S"),
                time_offset_sec = 0
        ),
        list(
                lab            = "UB",
                analyzer       = "CRDS.3",   # corrected from old "CRDS.2"
                input_path     = file.path(proj_root, "raw_data/UniBonn/CRDS.3"),
                mpv_levels     = c("8", "1", "9"),
                location_levels = c("N", "in", "S"),
                time_offset_sec = 0
        )
        # LUFA is skipped - see auto-copy block at the bottom
        # ,list(
        #         lab            = "LUFA",
        #         analyzer       = "CRDS.2",   # corrected from old "CRDS.3"
        #         input_path     = "<path-to-LUFA-CRDS-raw>",
        #         mpv_levels     = c("3", "1", "2"),
        #         location_levels = c("N", "in", "S"),
        #         time_offset_sec = -20        # clock-drift correction
        # )
)

# ---- output formatters (V5 schema) ------------------------------------------

#' Convert a piclean result (cycle-level) into the V5 long-format hourly file.
#' One row per (hour, location). Keeps a representative MPVPosition per row.
make_hourly_long_v5 <- function(df_cycle, lab_name, analyzer_name) {
        df_cycle %>%
                mutate(DATE.TIME = floor_date(DATE.TIME, unit = "hour")) %>%
                group_by(DATE.TIME, location) %>%
                summarise(
                        MPVPosition = as.character(first(MPVPosition)),
                        CO2 = mean(CO2, na.rm = TRUE),
                        CH4 = mean(CH4, na.rm = TRUE),
                        NH3 = mean(NH3, na.rm = TRUE),
                        H2O = if ("H2O" %in% names(.)) mean(H2O, na.rm = TRUE) else NA_real_,
                        .groups = "drop"
                ) %>%
                mutate(
                        lab      = lab_name,
                        analyzer = analyzer_name,
                        DATE.TIME = format(DATE.TIME, "%d/%m/%Y %H:%M")
                ) %>%
                select(DATE.TIME, MPVPosition, location, lab, analyzer,
                       CO2, CH4, NH3, H2O) %>%
                arrange(DATE.TIME, location)
}

#' Convert long V5 -> wide V5 format for one lab.
make_hourly_wide_v5 <- function(df_long, lab_name) {
        df_long %>%
                select(-MPVPosition, -lab) %>%
                pivot_wider(
                        names_from  = location,
                        values_from = c(CO2, CH4, NH3, H2O),
                        names_sep   = "_"
                ) %>%
                rename_with(
                        .cols = !c(DATE.TIME, analyzer),
                        .fn   = ~ paste0(.x, "_", lab_name)
                ) %>%
                arrange(DATE.TIME)
}

# ---- main loop --------------------------------------------------------------
for (cfg in lab_configs) {

        cat("\n============================================================\n")
        cat(sprintf("Processing %s %s\n", cfg$lab, cfg$analyzer))
        cat(sprintf("  input:  %s\n", cfg$input_path))
        cat("============================================================\n")

        # 1) cycle-level cleaning via piclean (writes its own working CSV in cwd)
        owd <- setwd(out_dir); on.exit(setwd(owd), add = TRUE)
        cycle_df <- piclean(
                input_path         = cfg$input_path,
                gas                = gases,
                start_time         = start_time,
                end_time           = end_time,
                flush              = flush_sec,
                interval           = interval_sec,
                MPVPosition.levels = cfg$mpv_levels,
                location.levels    = cfg$location_levels,
                lab                = cfg$lab,
                analyzer           = cfg$analyzer
        )
        setwd(owd)

        # 2) Force UTC tz + optional clock-drift correction.
        # piclean's as.POSIXct() doesn't set a timezone, so DATE.TIME comes
        # back interpreted in the system's local TZ (CEST/CET). Picarro
        # analyzers log in UTC, so use force_tz to relabel without shifting
        # the wall-clock value (otherwise output is off by the local UTC
        # offset, e.g. -2 h on a CEST machine).
        cycle_df <- cycle_df %>%
                mutate(DATE.TIME = lubridate::force_tz(DATE.TIME, "UTC") +
                                   cfg$time_offset_sec)

        # 3) snap timestamps onto the 7.5-min grid
        cycle_df <- cycle_df %>%
                mutate(DATE.TIME = round_to_interval(DATE.TIME, interval_sec = 450))

        # 3b) outlier removal on 7.5-min cycles, per location (Tukey 1.5*IQR).
        # Run before any hourly aggregation so the hourly mean isn't dragged
        # by a single off-cycle reading.
        cat(sprintf("  outlier removal (per location):\n"))
        cycle_df <- remove_outliers(
                cycle_df,
                exclude_cols = c("step_id", "measuring.time"),
                group_cols   = c("location")
        )

        # 4) write the 7.5-min intermediate
        cycle_out <- cycle_df %>%
                select(DATE.TIME, MPVPosition, location, lab, analyzer,
                       any_of(c("CO2", "CH4", "NH3", "H2O", "N2O")))
        cycle_path <- file.path(out_dir,
                                sprintf("%s_7.5_avg_%s_%s.csv",
                                        file_stub, cfg$lab, cfg$analyzer))
        write_excel_csv(cycle_out, cycle_path)
        cat(sprintf("  -> %s  (%d rows)\n", basename(cycle_path), nrow(cycle_out)))

        # 5) hourly long (V5)
        long_df <- make_hourly_long_v5(cycle_df, cfg$lab, cfg$analyzer)
        long_path <- file.path(out_dir,
                               sprintf("%s_long_%s_%s.csv",
                                       file_stub, cfg$lab, cfg$analyzer))
        write_excel_csv(long_df, long_path)
        cat(sprintf("  -> %s  (%d rows)\n", basename(long_path), nrow(long_df)))

        # 6) hourly wide (V5)
        wide_df <- make_hourly_wide_v5(long_df, cfg$lab)
        wide_path <- file.path(out_dir,
                               sprintf("%s_wide_%s_%s.csv",
                                       file_stub, cfg$lab, cfg$analyzer))
        write_excel_csv(wide_df, wide_path)
        cat(sprintf("  -> %s  (%d rows)\n", basename(wide_path), nrow(wide_df)))
}

# ---- LUFA: copy V5 long + wide into V6 so the folder is complete ------------
v5_dir <- file.path(proj_root, "clean_data", "Version_5")
lufa_files <- c(
        sprintf("%s_long_LUFA_CRDS.2.csv", file_stub),
        sprintf("%s_wide_LUFA_CRDS.2.csv", file_stub)
)
cat("\n--- copying LUFA CRDS.2 V5 -> V6 (cleaning skipped for LUFA) ---\n")
for (f in lufa_files) {
        src <- file.path(v5_dir, f)
        dst <- file.path(out_dir, f)
        if (file.exists(src)) {
                file.copy(src, dst, overwrite = TRUE)
                cat(sprintf("  copied %s\n", f))
        } else {
                cat(sprintf("  WARN: %s not found in Version_5 - skipped\n", f))
        }
}

cat("\nDone. Run the analysis script with data_version <- \"Version_6\".\n")
