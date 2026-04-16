#### Develpoment of function ####
harshcleaner <- function(
    input_path,
    gas = c("CO2", "CH4", "NH3", "H2O", "N2O"),
    start_time = NULL,
    end_time = NULL,
    flush = 0,
    interval = Inf,
    MPVPosition.levels = NULL,
    location.levels = NULL,
    valid_locations = c("in", "S", "N"),
    lab = NULL,
    analyzer = NULL,
    tz = "Europe/Berlin",
    save_csv = TRUE,
    output_dir = getwd(),
    prefix = NULL) 
{
  
  # ---------------------------
  # 1. Check packages
  # ---------------------------
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.")
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package 'tidyr' is required.")
  }
  if (!requireNamespace("lubridate", quietly = TRUE)) {
    stop("Package 'lubridate' is required.")
  }
  
  # ---------------------------
  # 2. Basic checks
  # ---------------------------
  if (!dir.exists(input_path)) {
    stop("input_path does not exist: ", input_path)
  }
  
  if (!is.numeric(flush) || length(flush) != 1 || flush < 0) {
    stop("flush must be one non-negative number.")
  }
  
  if (!is.numeric(interval) || length(interval) != 1 || interval <= 0) {
    stop("interval must be one positive number.")
  }
  
  if (!is.null(MPVPosition.levels) && !is.null(location.levels)) {
    if (length(MPVPosition.levels) != length(location.levels)) {
      stop("MPVPosition.levels and location.levels must have the same length.")
    }
  }
  
  if (!is.null(start_time)) {
    start_time <- as.POSIXct(start_time, tz = tz)
  }
  
  if (!is.null(end_time)) {
    end_time <- as.POSIXct(end_time, tz = tz)
  }
  
  # ---------------------------
  # 3. Columns needed
  # ---------------------------
  required_cols <- c("DATE", "TIME", "MPVPosition")
  all_needed_cols <- unique(c(required_cols, gas))
  
  # ---------------------------
  # 4. Find .dat files
  # ---------------------------
  dat_files <- list.files(
    path = input_path,
    pattern = "\\.dat$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  if (length(dat_files) == 0) {
    stop("No .dat files found in input_path.")
  }
  
  message("Found ", length(dat_files), " .dat files.")
  
  # ---------------------------
  # 5. Read all files
  # ---------------------------
  data_list <- list()
  
  for (i in seq_along(dat_files)) {
    this_file <- dat_files[i]
    
    tmp <- tryCatch(
      read.table(this_file, header = TRUE, stringsAsFactors = FALSE),
      error = function(e) {
        warning("Could not read file: ", this_file)
        return(NULL)
      }
    )
    
    if (is.null(tmp)) {
      next
    }
    
    existing_cols <- all_needed_cols[all_needed_cols %in% names(tmp)]
    
    if (length(existing_cols) == 0) {
      warning("No required columns found in file: ", basename(this_file))
      next
    }
    
    tmp <- tmp[, existing_cols, drop = FALSE]
    
    missing_cols <- setdiff(all_needed_cols, names(tmp))
    if (length(missing_cols) > 0) {
      for (col_name in missing_cols) {
        tmp[[col_name]] <- NA
      }
    }
    
    tmp <- tmp[, all_needed_cols, drop = FALSE]
    data_list[[length(data_list) + 1]] <- tmp
    
    message("Processed file ", i, " of ", length(dat_files))
  }
  
  if (length(data_list) == 0) {
    stop("No usable files were read.")
  }
  
  raw_data <- dplyr::bind_rows(data_list)
  
  # ---------------------------
  # 6. Create DATE.TIME
  # ---------------------------
  raw_data <- raw_data |>
    dplyr::mutate(
      DATE.TIME = as.POSIXct(
        paste(DATE, TIME),
        format = "%Y-%m-%d %H:%M:%S",
        tz = tz
      )
    ) |>
    dplyr::select(-DATE, -TIME)
  
  if (all(is.na(raw_data$DATE.TIME))) {
    stop("DATE.TIME could not be parsed. Please check DATE and TIME format.")
  }
  
  # ---------------------------
  # 7. Filter by time
  # ---------------------------
  if (!is.null(start_time)) {
    raw_data <- raw_data |>
      dplyr::filter(DATE.TIME >= start_time)
  }
  
  if (!is.null(end_time)) {
    raw_data <- raw_data |>
      dplyr::filter(DATE.TIME <= end_time)
  }
  
  if (nrow(raw_data) == 0) {
    stop("No rows left after time filtering.")
  }
  
  # ---------------------------
  # 8. Clean MPVPosition
  # ---------------------------
  raw_data <- raw_data |>
    dplyr::mutate(MPVPosition = suppressWarnings(as.numeric(MPVPosition))) |>
    dplyr::filter(!is.na(MPVPosition)) |>
    dplyr::filter(MPVPosition %% 1 == 0) |>
    dplyr::filter(MPVPosition != 0) |>
    dplyr::arrange(DATE.TIME)
  
  if (nrow(raw_data) == 0) {
    stop("No rows left after MPVPosition cleaning.")
  }
  
  # ---------------------------
  # 9. Optional MPV filter
  # ---------------------------
  if (!is.null(MPVPosition.levels)) {
    raw_data <- raw_data |>
      dplyr::filter(MPVPosition %in% MPVPosition.levels)
    
    if (nrow(raw_data) == 0) {
      stop("No rows left after filtering MPVPosition.levels.")
    }
  }
  
  # ---------------------------
  # 10. Add location mapping
  # ---------------------------
  if (!is.null(MPVPosition.levels) && !is.null(location.levels)) {
    map_df <- data.frame(
      MPVPosition = MPVPosition.levels,
      location = location.levels,
      stringsAsFactors = FALSE
    )
    
    raw_data <- raw_data |>
      dplyr::left_join(map_df, by = "MPVPosition")
  } else {
    raw_data$location <- NA_character_
  }
  
  # ---------------------------
  # 11. Convert gas columns
  # ---------------------------
  for (g in gas) {
    if (g %in% names(raw_data)) {
      raw_data[[g]] <- suppressWarnings(as.numeric(raw_data[[g]]))
    }
  }
  
  # ---------------------------
  # 12. Create step_id
  # New step when MPVPosition changes
  # ---------------------------
  raw_data <- raw_data |>
    dplyr::arrange(DATE.TIME) |>
    dplyr::mutate(
      step_change = c(TRUE, diff(MPVPosition) != 0),
      step_id = cumsum(step_change)
    ) |>
    dplyr::select(-step_change)
  
  # ---------------------------
  # 13. Apply flush and summarise
  # ---------------------------
  step_data <- raw_data |>
    dplyr::group_by(step_id, MPVPosition, location) |>
    dplyr::arrange(DATE.TIME, .by_group = TRUE) |>
    dplyr::mutate(
      time_rank = dplyr::row_number(),
      timestamp_for_step = dplyr::last(DATE.TIME)
    ) |>
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(gas),
        ~ dplyr::if_else(time_rank <= flush, NA_real_, .x)
      )
    ) |>
    dplyr::filter(time_rank <= interval) |>
    dplyr::summarise(
      DATE.TIME = dplyr::last(timestamp_for_step),
      measuring.time = sum(time_rank > flush & time_rank <= interval),
      dplyr::across(
        dplyr::all_of(gas),
        ~ if (all(is.na(.x))) NA_real_ else mean(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    )
  
  # ---------------------------
  # 14. Optional NH3 conversion
  # Use only if NH3 raw unit is ppb
  # ---------------------------
  if ("NH3" %in% names(step_data)) {
    step_data <- step_data |>
      dplyr::mutate(NH3 = NH3 / 1000)
  }
  
  # ---------------------------
  # 15. Add metadata
  # ---------------------------
  if (!is.null(lab)) {
    step_data$lab <- lab
  }
  
  if (!is.null(analyzer)) {
    step_data$analyzer <- analyzer
  }
  
  # ---------------------------
  # 16. Hourly long table
  # ---------------------------
  gases_present <- intersect(c("CO2", "CH4", "NH3", "H2O", "N2O"), names(step_data))
  
  hourly_long <- step_data |>
    dplyr::filter(location %in% valid_locations) |>
    dplyr::mutate(DATE.HOUR = lubridate::floor_date(DATE.TIME, unit = "hour")) |>
    dplyr::group_by(DATE.HOUR, location, lab, analyzer) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(gases_present),
        ~ if (all(is.na(.x))) NA_real_ else mean(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    )
  
  # ---------------------------
  # 17. Hourly wide table
  # ---------------------------
  if (nrow(hourly_long) > 0) {
    dummy_grid <- tidyr::expand_grid(
      DATE.HOUR = seq(
        from = min(hourly_long$DATE.HOUR, na.rm = TRUE),
        to = max(hourly_long$DATE.HOUR, na.rm = TRUE),
        by = "hour"
      ),
      analyzer = unique(hourly_long$analyzer),
      location = unique(hourly_long$location),
      lab = unique(hourly_long$lab)
    )
    
    hourly_wide <- dummy_grid |>
      dplyr::left_join(
        hourly_long,
        by = c("DATE.HOUR", "location", "analyzer", "lab")
      ) |>
      tidyr::pivot_wider(
        names_from = location,
        values_from = dplyr::all_of(gases_present),
        names_sep = "_"
      ) |>
      dplyr::arrange(DATE.HOUR)
  } else {
    hourly_wide <- data.frame()
  }
  
  # ---------------------------
  # 18. Helper to create final site table
  # ---------------------------
  make_final_site <- function(df_wide, suffix, site_code_label, lab_value, analyzer_value) {
    gases_std <- c("CO2", "CH4", "NH3", "H2O", "N2O")
    wanted_cols <- paste0(gases_std, "_", suffix)
    present_cols <- wanted_cols[wanted_cols %in% names(df_wide)]
    
    df_site <- df_wide |>
      dplyr::select(DATE.HOUR, analyzer, lab, dplyr::all_of(present_cols)) |>
      dplyr::rename_with(
        ~ gsub(paste0("_", suffix), "", .x),
        dplyr::all_of(present_cols)
      )
    
    for (g in gases_std) {
      if (!g %in% names(df_site)) {
        df_site[[g]] <- NA_real_
      }
    }
    
    df_out <- df_site |>
      dplyr::mutate(
        LOCATION = "Gross Kreutz",
        SITE_CODE = site_code_label,
        MEASUREMENT_PERIOD = "",
        DATE_TIME = DATE.HOUR,
        MEASUREMENT_POINT_CODE = lab_value,
        GAS_MEASUREMENT_METHOD = paste0("LINE_", analyzer_value),
        CO2_DRY_PPM = CO2,
        CH4_DRY_PPM = CH4,
        NH3_DRY_PPM = NH3,
        H2O_VOL_PCT = H2O,
        N2O_DRY_PPM = N2O,
        NO2_DRY_PPM = "",
        NO_DRY_PPM = "",
        CO_DRY_PPM = "",
        COMMENTS = "",
        AVG_TOTAL_NUMBER_ANIMALS = "",
        AVG_WEIGHT_KG = "",
        AVG_DAYS_IN_PREGNANCY = "",
        AVG_MILK_YIELD_CORR_KGAD = "",
        IM_AIR_TEMPERATURE_C = ""
      ) |>
      dplyr::select(
        LOCATION, SITE_CODE, MEASUREMENT_PERIOD, DATE_TIME,
        MEASUREMENT_POINT_CODE, GAS_MEASUREMENT_METHOD,
        CO2_DRY_PPM, CH4_DRY_PPM, NH3_DRY_PPM, H2O_VOL_PCT,
        N2O_DRY_PPM, NO2_DRY_PPM, NO_DRY_PPM, CO_DRY_PPM,
        COMMENTS, AVG_TOTAL_NUMBER_ANIMALS, AVG_WEIGHT_KG,
        AVG_DAYS_IN_PREGNANCY, AVG_MILK_YIELD_CORR_KGAD,
        IM_AIR_TEMPERATURE_C
      ) |>
      dplyr::arrange(DATE_TIME)
    
    return(df_out)
  }
  
  # ---------------------------
  # 19. Final outputs
  # ---------------------------
  if (nrow(hourly_wide) > 0) {
    final_IN <- make_final_site(hourly_wide, "in", "IN", lab, analyzer)
    final_OUT_South <- make_final_site(hourly_wide, "S", "OUT_South", lab, analyzer)
    final_OUT_North <- make_final_site(hourly_wide, "N", "OUT_North", lab, analyzer)
  } else {
    final_IN <- data.frame()
    final_OUT_South <- data.frame()
    final_OUT_North <- data.frame()
  }
  
  # ---------------------------
  # 20. Save csv files
  # ---------------------------
  if (isTRUE(save_csv)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    start_str <- if (!is.null(start_time)) format(start_time, "%Y%m%d.%H%M") else "startNA"
    end_str <- if (!is.null(end_time)) format(end_time, "%Y%m%d.%H%M") else "endNA"
    
    name_parts <- c(start_str, end_str, lab, analyzer, paste0(interval, "avg"), prefix)
    name_parts <- name_parts[!is.na(name_parts)]
    name_parts <- name_parts[name_parts != ""]
    base_name <- paste(name_parts, collapse = "_")
    
    utils::write.csv(
      step_data,
      file.path(output_dir, paste0(base_name, "_step.csv")),
      row.names = FALSE,
      quote = FALSE
    )
    
    utils::write.csv(
      hourly_long,
      file.path(output_dir, paste0(base_name, "_hourly_long.csv")),
      row.names = FALSE,
      quote = FALSE
    )
    
    utils::write.csv(
      hourly_wide,
      file.path(output_dir, paste0(base_name, "_hourly_wide.csv")),
      row.names = FALSE,
      quote = FALSE
    )
    
    utils::write.csv(
      final_IN,
      file.path(output_dir, paste0(base_name, "_final_IN.csv")),
      row.names = FALSE,
      quote = FALSE
    )
    
    utils::write.csv(
      final_OUT_South,
      file.path(output_dir, paste0(base_name, "_final_OUT_South.csv")),
      row.names = FALSE,
      quote = FALSE
    )
    
    utils::write.csv(
      final_OUT_North,
      file.path(output_dir, paste0(base_name, "_final_OUT_North.csv")),
      row.names = FALSE,
      quote = FALSE
    )
    
    message("Files saved to: ", output_dir)
  }
  
  # ---------------------------
  # 21. Return results
  # ---------------------------
  return(
    list(
      raw_data = raw_data,
      step_data = step_data,
      hourly_long = hourly_long,
      hourly_wide = hourly_wide,
      final_IN = final_IN,
      final_OUT_South = final_OUT_South,
      final_OUT_North = final_OUT_North
    )
  )
}
