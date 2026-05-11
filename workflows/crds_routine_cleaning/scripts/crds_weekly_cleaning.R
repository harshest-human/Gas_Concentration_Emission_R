##### Function Development ####

piclean <- function(input_path, gas, start_time, end_time, flush, interval,
                    MPVPosition.levels = NULL, location.levels = NULL,
                    lab = NULL, analyzer = NULL,
                    output_dir = getwd()) {
  library(data.table)
  library(lubridate)

  required_cols    <- c("DATE", "TIME", "MPVPosition")
  all_needed_cols  <- unique(c(required_cols, gas))

  # Step 1: Load only .dat files inside the requested date window.
  # Picarro layout: <input_path>/datalog_user/YYYY/MM/DD/*.dat. We construct
  # the day-folder paths from start_time/end_time and only descend into those,
  # so the whole history isn't scanned. If the layout isn't found, we fall back
  # to a full recursive scan (compatible with the old raw_data folder).
  appendData <- function() {
    start_date <- as.Date(as.POSIXct(start_time))
    end_date   <- as.Date(as.POSIXct(end_time))
    day_seq    <- seq(start_date, end_date, by = "day")

    base <- if (dir.exists(file.path(input_path, "datalog_user"))) {
      file.path(input_path, "datalog_user")
    } else {
      input_path
    }

    day_dirs <- file.path(
      base,
      format(day_seq, "%Y"),
      format(day_seq, "%m"),
      format(day_seq, "%d")
    )
    day_dirs <- day_dirs[dir.exists(day_dirs)]

    if (length(day_dirs) > 0) {
      cat(sprintf("Scanning %d day folder(s) in %s\n", length(day_dirs), base))
      dat_files <- unlist(lapply(day_dirs, function(d) {
        list.files(d, pattern = "\\.dat$", full.names = TRUE)
      }), use.names = FALSE)
    } else {
      cat("No date-structured day folders found; falling back to recursive scan.\n")
      dat_files <- list.files(input_path, recursive = TRUE,
                              pattern = "\\.dat$", full.names = TRUE)
    }

    total_files <- length(dat_files)
    cat(sprintf("Found %d .dat file(s) for date range\n", total_files))
    data_list <- vector("list", total_files)

    for (i in seq_along(dat_files)) {
      current_data <- tryCatch(fread(dat_files[i]), error = function(e) NULL)
      if (is.null(current_data)) {
        cat(sprintf("Warning: Failed to read file: %s\n", dat_files[i]))
        next
      }
      existing_cols <- intersect(all_needed_cols, names(current_data))
      if (length(existing_cols) == 0) {
        cat(sprintf("Warning: No matching columns in file: %s\n", dat_files[i]))
        next
      }
      data_list[[i]] <- current_data[, ..existing_cols]
      cat(sprintf("Processed file %d of %d\n", i, total_files))
    }
    data_list[!vapply(data_list, is.null, logical(1))]
  }

  list_of_data <- appendData()

  if (length(list_of_data) == 0) {
    stop("No data frames created. Check input files and column names.")
  }

  merged_data <- rbindlist(list_of_data, fill = TRUE)

  # Step 2: Merge date and time & filter early
  cat("Merging DATE and TIME column...\n")
  merged_data[, DATE.TIME := as.POSIXct(paste(DATE, TIME),
                                        format = "%Y-%m-%d %H:%M:%S")]
  merged_data[, c("DATE", "TIME") := NULL]

  if (!is.null(start_time) && !is.null(end_time)) {
    cat("Filtering data between:\n",
        "Start:", start_time, "\n",
        "End:  ", end_time, "\n")
    merged_data <- merged_data[
      DATE.TIME >= as.POSIXct(start_time) &
      DATE.TIME <= as.POSIXct(end_time)
    ]
  }

  # Step 3: Remove invalid MPVPositions, then build step_id on clean integers
  cat("Removing non-integer and zero MPVPosition values...\n")
  merged_data <- merged_data[
    !is.na(MPVPosition) & MPVPosition %% 1 == 0 & MPVPosition != 0
  ]
  setorder(merged_data, DATE.TIME)
  merged_data[, step_id := cumsum(c(1L, as.integer(diff(MPVPosition) != 0)))]

  # Step 4: Apply MPVPosition & location levels if provided
  if (!is.null(MPVPosition.levels)) {
    merged_data <- merged_data[MPVPosition %in% MPVPosition.levels]
    merged_data[, MPVPosition := factor(MPVPosition, levels = MPVPosition.levels)]
  }

  if (!is.null(location.levels) && !is.null(MPVPosition.levels)) {
    merged_data[, location := factor(
      location.levels[as.integer(MPVPosition)],
      levels = location.levels
    )]
  }

  # Step 5: Flush — mark the first `flush` rows of each step as NA on gas cols
  setorder(merged_data, step_id, MPVPosition, DATE.TIME)
  merged_data[, `:=`(
    timestamp_for_step = last(DATE.TIME),
    time_rank          = seq_len(.N)
  ), by = .(step_id, MPVPosition)]

  gas_present <- intersect(gas, names(merged_data))
  for (g in gas_present) {
    merged_data[time_rank <= flush, (g) := NA_real_]
  }

  # Step 6: Summarise per step (only first `interval` rows feed the mean)
  merged_data <- merged_data[time_rank <= interval]

  summarized <- merged_data[, c(
    list(
      DATE.TIME      = last(timestamp_for_step),
      measuring.time = sum(time_rank > flush)
    ),
    lapply(.SD, mean, na.rm = TRUE)
  ),
  by = .(step_id, MPVPosition, location),
  .SDcols = gas_present]

  # NH3 unit conversion (kept bug-for-bug from the dplyr version: if NH3 isn't
  # present, the column is added as NA — that's the original behavior)
  if ("NH3" %in% names(summarized)) {
    summarized[, NH3 := NH3 / 1000]
  } else {
    summarized[, NH3 := NA_real_]
  }

  # Step 8: Add lab & analyzer info
  cat("Adding lab and analyzer info...\n")
  if (!is.null(lab))      summarized[, lab := ..lab]
  if (!is.null(analyzer)) summarized[, analyzer := ..analyzer]

  # Step 9: Save output
  cat("Saving data to CSV...\n")
  start_str <- format(as.POSIXct(start_time), "%Y%m%d.%H%M")
  end_str   <- format(as.POSIXct(end_time),   "%Y%m%d.%H%M")
  file_name <- paste0(start_str, "_", end_str, "_",
                      lab, "_", interval, "avg_", analyzer, ".csv")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  full_output_path <- file.path(output_dir, file_name)

  fwrite(summarized, full_output_path, quote = FALSE, dateTimeAs = "write.csv")

  cat("Data successfully processed and saved to:", full_output_path, "\n")
  cat("Dataframe contains", ncol(summarized), "variables\n")

  return(summarized)
}


reshape_crds <- function(df_raw) {
  library(data.table)
  library(lubridate)

  df_raw          <- as.data.table(df_raw)
  valid_locations <- c("in", "S", "N")
  gases           <- c("CO2", "CH4", "NH3", "H2O", "N2O")
  gases_present   <- intersect(gases, names(df_raw))

  df_hourly_long <- df_raw[
    location %in% valid_locations,
    lapply(.SD, mean, na.rm = TRUE),
    by = .(DATE.HOUR = floor_date(DATE.TIME, "hour"),
           location, lab, analyzer),
    .SDcols = gases_present
  ]

  df_dummy <- CJ(
    DATE.HOUR = seq(min(df_hourly_long$DATE.HOUR),
                    max(df_hourly_long$DATE.HOUR), by = "hour"),
    analyzer  = unique(df_hourly_long$analyzer),
    location  = unique(df_hourly_long$location),
    lab       = unique(df_hourly_long$lab)
  )

  merged <- merge(df_dummy, df_hourly_long,
                  by = c("DATE.HOUR", "location", "analyzer", "lab"),
                  all.x = TRUE)

  df_hourly_wide <- dcast(
    merged,
    DATE.HOUR + analyzer + lab ~ location,
    value.var = gases_present
  )
  setorder(df_hourly_wide, DATE.HOUR)

  make_final <- function(df_wide, suffix, site_code_label) {
    df_wide <- as.data.table(df_wide)
    cols    <- paste0(gases, "_", suffix)
    cols    <- intersect(cols, names(df_wide))

    df_site <- df_wide[, c("DATE.HOUR", "analyzer", "lab", cols),
                       with = FALSE]
    if (length(cols) > 0) {
      setnames(df_site, cols, gsub(paste0("_", suffix), "", cols))
    }

    for (g in gases) {
      if (!g %in% names(df_site)) df_site[, (g) := NA_real_]
    }

    df_out <- df_site[, .(
      LOCATION                 = "Gross Kreutz",
      SITE_CODE                = site_code_label,
      MEASUREMENT_PERIOD       = "",
      DATE_TIME                = DATE.HOUR,
      MEASUREMENT_POINT_CODE   = lab,
      GAS_MEASUREMENT_METHOD   = paste0("LINE_", analyzer),
      CO2_DRY_PPM              = CO2,
      CH4_DRY_PPM              = CH4,
      NH3_DRY_PPM              = NH3,
      H2O_VOL_PCT              = H2O,
      N2O_DRY_PPM              = N2O,
      NO2_DRY_PPM              = "",
      NO_DRY_PPM               = "",
      CO_DRY_PPM               = "",
      COMMENTS                 = "",
      AVG_TOTAL_NUMBER_ANIMALS = "",
      AVG_WEIGHT_KG            = "",
      AVG_DAYS_IN_PREGNANCY    = "",
      AVG_MILK_YIELD_CORR_KGAD = "",
      IM_AIR_TEMPERATURE_C     = ""
    )]
    setorder(df_out, DATE_TIME)

    setattr(df_out, "start_dt", min(df_out$DATE_TIME, na.rm = TRUE))
    setattr(df_out, "end_dt",   max(df_out$DATE_TIME, na.rm = TRUE))

    df_out
  }

  list(
    hourly_wide     = df_hourly_wide,
    final_IN        = make_final(df_hourly_wide, "in", "IN"),
    final_OUT_South = make_final(df_hourly_wide, "S",  "OUT_South"),
    final_OUT_North = make_final(df_hourly_wide, "N",  "OUT_North")
  )
}


##### Orchestrator ####

autocrds <- function(input_path, output_path,
                     gas, start_time, end_time, flush, interval,
                     MPVPosition.levels = NULL, location.levels = NULL,
                     lab = NULL, analyzer = NULL,
                     sites = c("IN", "S", "N"),
                     plot = TRUE) {

  library(data.table)
  library(writexl)
  library(ggplot2)

  cat("\n--- AUTOCRDS START ---\n")

  # Output subdirectories
  dir_step   <- file.path(output_path, "step_avg")
  dir_hourly <- file.path(output_path, "hourly_in_out_avg")
  dir_ktbl   <- file.path(output_path, "KTBL_database")
  dir_plots  <- file.path(output_path, "plots")
  for (d in c(dir_step, dir_hourly, dir_ktbl, dir_plots)) {
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
  }

  # 1. Per-step summary (piclean writes its CSV into dir_step)
  step_summary <- piclean(
    input_path         = input_path,
    gas                = gas,
    start_time         = start_time,
    end_time           = end_time,
    flush              = flush,
    interval           = interval,
    MPVPosition.levels = MPVPosition.levels,
    location.levels    = location.levels,
    lab                = lab,
    analyzer           = analyzer,
    output_dir         = dir_step
  )

  # 2. Reshape to hourly wide + per-site KTBL tables
  reshape_out <- reshape_crds(step_summary)

  # 3. Hourly in/out CSV
  start_str <- format(as.POSIXct(start_time), "%Y%m%d")
  end_str   <- format(as.POSIXct(end_time),   "%Y%m%d")
  hourly_path <- file.path(
    dir_hourly,
    paste0("H_", analyzer, "+9_", start_str, "_", end_str, ".csv")
  )
  fwrite(reshape_out$hourly_wide, hourly_path, quote = FALSE,
         dateTimeAs = "write.csv")
  cat("Wrote hourly in/out:", hourly_path, "\n")

  # 4. KTBL Excel per requested site
  site_to_info <- list(
    IN = list(key = "final_IN",        label = "IN"),
    S  = list(key = "final_OUT_South", label = "OUT_South"),
    N  = list(key = "final_OUT_North", label = "OUT_North")
  )

  for (site in sites) {
    info <- site_to_info[[site]]
    if (is.null(info)) {
      warning("Unknown site '", site, "' — expected one of IN, S, N. Skipping.")
      next
    }
    tbl <- reshape_out[[ info$key ]]
    if (is.null(tbl) || nrow(tbl) == 0) {
      warning("No rows for site '", site, "'. Skipping.")
      next
    }
    tbl_out <- copy(tbl)
    tbl_out[, DATE_TIME := format(DATE_TIME, "%Y-%m-%d %H:%M:%S")]

    ktbl_path <- file.path(
      dir_ktbl,
      paste0("GrossKreutz_", info$label, "_", start_str, "_", end_str, ".xlsx")
    )
    write_xlsx(as.data.frame(tbl_out), ktbl_path)
    cat("Wrote KTBL", info$label, ":", ktbl_path, "\n")
  }

  # 5. Plots — one PNG per gas, faceted by 24h segment
  if (isTRUE(plot)) {
    cat("Generating plots...\n")

    plot_dt <- copy(step_summary)
    plot_dt[, location := as.character(location)]
    plot_dt <- plot_dt[location %in% c("in", "S", "N")]
    setorder(plot_dt, DATE.TIME)
    plot_dt[, time_diff := c(0, diff(as.numeric(DATE.TIME)))]
    plot_dt[, segment := cumsum(time_diff > 86400)]

    gas_present     <- intersect(gas, names(plot_dt))
    location_colors <- c("in" = "black", "S" = "grey50", "N" = "steelblue")

    for (g in gas_present) {
      p <- ggplot(plot_dt,
                  aes(x = DATE.TIME, y = .data[[g]], color = location)) +
        geom_line() +
        scale_color_manual(values = location_colors) +
        facet_wrap(~segment, scales = "free_x") +
        labs(y = paste0(g, " (ppm)"), x = "Time") +
        theme_minimal()

      ggsave(
        file.path(dir_plots,
                  paste0(g, "_", start_str, "_", end_str, "_timeseries.png")),
        p, width = 12, height = 6
      )
    }
    cat("Plots saved\n")
  }

  cat("--- AUTOCRDS COMPLETE ---\n")

  invisible(list(
    step_summary = step_summary,
    reshape      = reshape_out
  ))
}


##### 20260309_20260509 #####
CRDS8_20260309_20260509 <- autocrds(
  input_path  = "D:/Data_Analysis_R/owncloud_sync_data/CRDS08_raw",
  output_path = "D:/Data_Analysis_R/Gas_Concentration_Emission_R/workflows/crds_routine_cleaning/clean_data/crds_clean",
  gas         = c("CO2", "CH4", "NH3", "H2O", "N2O"),
  start_time  = "2026-03-09 08:35:00",
  end_time    = "2026-05-09 01:54:17",
  flush       = 60,
  interval    = 240,
  MPVPosition.levels = as.character(1:9),
  location.levels    = c("1","2","3","4","5","6","7","in","S"),
  lab         = "ATB",
  analyzer    = "CRDS8",
  sites       = c("IN", "S")
)
