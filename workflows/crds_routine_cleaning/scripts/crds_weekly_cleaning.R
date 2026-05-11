##### Function Develipment ####
crds_routine_process <- function(
    input_path,
    gas = c("CO2", "CH4", "NH3", "H2O", "N2O"),
    start_time,
    end_time,
    flush = 60,
    interval = 240,
    MPVPosition.levels = as.character(1:9),
    location_map = c(
      "1"="1","2"="2","3"="3",
      "4"="4","5"="5","6"="6",
      "7"="7","8"="in","9"="S"
    ),
    lab = "ATB",
    analyzer = "CRDS8",
    output_path
) {
  
  library(data.table)
  library(lubridate)
  library(ggplot2)
  library(writexl)
  
  cat("\n--- CRDS ROUTINE START ---\n")
  
  # ---------------- OUTPUT PATHS ----------------
  year_str  <- format(as.POSIXct(start_time), "%Y")
  start_str <- format(as.POSIXct(start_time), "%Y%m%d")
  end_str   <- format(as.POSIXct(end_time), "%Y%m%d")
  
  dir_hourly <- file.path(output_path, "hourly_in_out_avg", year_str)
  dir_step   <- file.path(output_path, "step_avg")
  dir_ktbl   <- file.path(output_path, "KTBL_database")
  dir_plots  <- file.path(output_path, "plots")
  
  dir.create(dir_hourly, TRUE, FALSE)
  dir.create(dir_step, TRUE, FALSE)
  dir.create(dir_ktbl, TRUE, FALSE)
  dir.create(dir_plots, TRUE, FALSE)
  
  cat("Directories ready\n")
  
  # ---------------- LOAD DATA ----------------
  files <- list.files(input_path, full.names = TRUE, recursive = TRUE, pattern = "\\.dat$")
  files <- files[file.exists(files)]
  
  cat("Files found:", length(files), "\n")
  
  dt <- rbindlist(lapply(seq_along(files), function(i) {
    cat("Reading file", i, "of", length(files), "\n")
    tryCatch(fread(files[i]), error = function(e) NULL)
  }), fill = TRUE)
  
  setDT(dt)
  
  cat("Rows loaded:", nrow(dt), "\n")
  
  # ---------------- TIMESTAMP ----------------
  if (!"timestamp" %in% names(dt)) {
    if (all(c("DATE", "TIME") %in% names(dt))) {
      dt[, timestamp := as.POSIXct(
        paste(DATE, TIME),
        format = "%Y-%m-%d %H:%M:%S",
        tz = "UTC"
      )]
      dt[, c("DATE", "TIME") := NULL]
    } else {
      stop("No timestamp or DATE/TIME columns found")
    }
  }
  
  # ---------------- FILTER TIME ----------------
  dt <- dt[
    timestamp >= as.POSIXct(start_time) &
      timestamp <= as.POSIXct(end_time)
  ]
  
  setorder(dt, timestamp)
  
  cat("After filtering:", nrow(dt), "\n")
  
  # ---------------- METADATA FIX (CRITICAL) ----------------
  dt[, lab := lab]
  dt[, analyzer := analyzer]
  
  # ---------------- LOCATION MAP ----------------
  dt[, MPVPosition := as.character(MPVPosition)]
  dt[, location := location_map[MPVPosition]]
  
  # ---------------- CYCLE + FLUSH ----------------
  dt[, cycle_id := rleid(MPVPosition)]
  dt[, t_in_cycle := as.numeric(timestamp - first(timestamp)), by = cycle_id]
  
  dt_clean <- dt[t_in_cycle > flush]
  
  cat("After flush:", nrow(dt_clean), "\n")
  
  # =====================================================
  # FIXED AGGREGATION (NO dplyr, NO BUGS)
  # =====================================================
  
  cat("Aggregating...\n")
  
  dt_clean[, time_bin := as.numeric(timestamp) %/% interval * interval]
  dt_clean[, timestamp_bin := as.POSIXct(time_bin, origin = "1970-01-01", tz = "UTC")]
  
  agg <- dt_clean[
    ,
    c(
      list(
        timestamp = timestamp_bin,
        location = location,
        lab = lab,
        analyzer = analyzer
      ),
      lapply(.SD, mean, na.rm = TRUE)
    ),
    by = .(timestamp_bin, location, lab, analyzer),
    .SDcols = gas
  ]
  
  setnames(agg, "timestamp_bin", "timestamp")
  
  cat("Aggregated rows:", nrow(agg), "\n")
  
  # ---------------- NH3 CONVERSION ----------------
  if ("NH3" %in% names(agg)) {
    cat("Converting NH3 ppb -> ppm\n")
    agg[, NH3 := NH3 / 1000]  # convert ppb → ppm
  }
  
  # ---------------- 24h SEGMENTS ----------------
  setorder(agg, timestamp)
  agg[, time_diff := c(0, diff(as.numeric(timestamp)))]
  agg[, segment := cumsum(time_diff > 86400)]
  
  # =====================================================
  # RESHAPE (PURE data.table VERSION OF reshape_crds)
  # =====================================================
  
  cat("Reshaping hourly structure...\n")
  
  valid_locations <- c("in", "S", "N")
  
  long <- melt(
    agg[location %in% valid_locations],
    id.vars = c("timestamp", "location", "lab", "analyzer", "segment"),
    measure.vars = gas,
    variable.name = "gas",
    value.name = "value"
  )
  
  long_hour <- long[
    ,
    .(value = mean(value, na.rm = TRUE)),
    by = .(
      DATE.HOUR = floor_date(timestamp, "hour"),
      location,
      lab,
      analyzer,
      gas
    )
  ]
  
  dummy <- CJ(
    DATE.HOUR = seq(min(long_hour$DATE.HOUR),
                    max(long_hour$DATE.HOUR),
                    by = "hour"),
    location = unique(long_hour$location),
    lab = unique(long_hour$lab),
    analyzer = unique(long_hour$analyzer),
    gas = gas
  )
  
  wide <- dcast(
    merge(dummy, long_hour, all.x = TRUE),
    DATE.HOUR ~ location + gas,
    value.var = "value"
  )
  
  cat("Reshape complete\n")
  
  # =====================================================
  # EXPORT FILES
  # =====================================================
  
  fwrite(
    wide,
    file.path(dir_hourly,
              paste0("H_", analyzer, "_", start_str, "_", end_str, ".csv"))
  )
  
  write_xlsx(
    agg,
    file.path(dir_step,
              paste0(start_str, "_", end_str, "_step.xlsx"))
  )
  
  write_xlsx(
    agg,
    file.path(dir_ktbl,
              paste0("GrossKreutz_", start_str, "_to_", end_str, ".xlsx"))
  )
  
  cat("Files exported\n")
  
  # =====================================================
  # PLOTS (with styling rules)
  # =====================================================
  
  cat("Generating plots...\n")
  
  for (g in gas) {
    
    p <- ggplot(agg, aes(x = timestamp, y = get(g), color = location)) +
      geom_line() +
      scale_color_manual(values = c(
        in = "black",
        S = "grey50",
        N = "steelblue"
      )) +
      facet_wrap(~segment, scales = "free_x") +
      labs(y = paste0(g, " (ppm)"), x = "Time") +
      theme_minimal()
    
    ggsave(
      file.path(dir_plots, paste0(g, "_timeseries.png")),
      p,
      width = 12,
      height = 6
    )
  }
  
  cat("Plots saved\n")
  cat("--- CRDS COMPLETE ---\n")
  
  return(list(
    raw = dt,
    cleaned = dt_clean,
    aggregated = agg,
    hourly_wide = wide
  ))
}


##### 20260501_20260509 #####
CRDS8_20260501_20260509 <- crds_routine_process(
  
  input_path = "D:/Data_Analysis_R/Gas_Concentration_Emission_R/workflows/crds_routine_cleaning/raw_data",
  
  start_time = "2026-05-01 02:34:56",
  end_time   = "2026-05-09 01:54:17",
  
  output_path = "D:/Data_Analysis_R/Gas_Concentration_Emission_R/workflows/crds_routine_cleaning/clean_data/crds_clean"
)

