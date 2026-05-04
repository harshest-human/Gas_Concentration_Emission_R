###### libraries ######
library(dplyr)
library(lubridate)
library(tidyr)
library(readr)
library(writexl)

resolve_project_dir <- function() {
  project_marker <- function(path) {
    dir.exists(file.path(path, "workflows", "crds_routine_cleaning")) &&
      dir.exists(file.path(path, "scripts", "cleaning"))
  }

  candidate_paths <- character()
  args_all <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- sub(file_arg, "", args_all[grepl(file_arg, args_all)])

  if (length(script_path) > 0) {
    candidate_paths <- c(
      candidate_paths,
      normalizePath(file.path(dirname(script_path[1]), "..", "..", ".."), winslash = "/", mustWork = FALSE)
    )
  }

  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    active_path <- tryCatch(
      rstudioapi::getActiveDocumentContext()$path,
      error = function(e) ""
    )

    if (nzchar(active_path)) {
      candidate_paths <- c(
        candidate_paths,
        normalizePath(file.path(dirname(active_path), "..", "..", ".."), winslash = "/", mustWork = FALSE)
      )
    }
  }

  candidate_paths <- c(
    candidate_paths,
    normalizePath(getwd(), winslash = "/", mustWork = FALSE),
    normalizePath(file.path(getwd(), "Gas_Concentration_Emission_R"), winslash = "/", mustWork = FALSE)
  )

  candidate_paths <- unique(candidate_paths[nzchar(candidate_paths)])

  for (path in candidate_paths) {
    if (project_marker(path)) {
      return(path)
    }
  }

  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

format_datetime_text <- function(x, time_zone = "Europe/Berlin") {
  if (length(x) == 0) {
    return(character())
  }

  format(
    as.POSIXct(x, tz = time_zone),
    tz = time_zone,
    usetz = FALSE,
    format = "%Y-%m-%d %H:%M:%S"
  )
}

format_datetime_columns <- function(data, columns, time_zone = "Europe/Berlin") {
  for (column_name in columns) {
    if (column_name %in% names(data)) {
      data[[column_name]] <- format_datetime_text(data[[column_name]], time_zone = time_zone)
    }
  }

  data
}

extract_script_end_times <- function(script_paths, time_zone = "Europe/Berlin") {
  end_times <- as.POSIXct(character(), tz = time_zone)

  for (script_path in script_paths[file.exists(script_paths)]) {
    lines <- readLines(script_path, warn = FALSE)
    matches <- regmatches(
      lines,
      gregexpr('end_time\\s*=\\s*"[^"]+"', lines, perl = TRUE)
    )

    matched_values <- unlist(matches, use.names = FALSE)
    if (length(matched_values) == 0) {
      next
    }

    datetime_strings <- sub('.*"([^"]+)".*', "\\1", matched_values)
    parsed <- as.POSIXct(datetime_strings, tz = time_zone, format = "%Y-%m-%d %H:%M:%S")
    end_times <- c(end_times, parsed[!is.na(parsed)])
  }

  end_times
}

extract_filename_end_times <- function(step_dir, time_zone = "Europe/Berlin") {
  files <- list.files(step_dir, pattern = "\\.csv$", full.names = FALSE)

  if (length(files) == 0) {
    return(as.POSIXct(character(), tz = time_zone))
  }

  end_strings <- sub("^\\d{8}\\.\\d{4}-(\\d{8}\\.\\d{4}).*$", "\\1", files, perl = TRUE)
  valid_strings <- end_strings[grepl("^\\d{8}\\.\\d{4}$", end_strings)]

  if (length(valid_strings) == 0) {
    return(as.POSIXct(character(), tz = time_zone))
  }

  as.POSIXct(valid_strings, tz = time_zone, format = "%Y%m%d.%H%M")
}

find_last_cleaned_end_time <- function(log_path, step_dir, reference_scripts, time_zone = "Europe/Berlin") {
  if (file.exists(log_path)) {
    run_log <- suppressMessages(read_csv(log_path, show_col_types = FALSE))

    if ("end_time" %in% names(run_log)) {
      parsed <- as.POSIXct(run_log$end_time, tz = time_zone, format = "%Y-%m-%d %H:%M:%S")
      parsed <- parsed[!is.na(parsed)]
      if (length(parsed) > 0) {
        return(max(parsed))
      }
    }
  }

  script_end_times <- extract_script_end_times(reference_scripts, time_zone = time_zone)
  if (length(script_end_times) > 0) {
    return(max(script_end_times))
  }

  filename_end_times <- extract_filename_end_times(step_dir, time_zone = time_zone)
  if (length(filename_end_times) > 0) {
    return(max(filename_end_times))
  }

  as.POSIXct(NA, tz = time_zone)
}

prompt_datetime <- function(prompt_text, default_value = NULL, time_zone = "Europe/Berlin") {
  repeat {
    default_text <- if (!is.null(default_value) && !is.na(default_value)) {
      format_datetime_text(default_value, time_zone = time_zone)
    } else {
      NULL
    }

    message_text <- prompt_text
    if (!is.null(default_text) && nzchar(default_text)) {
      message_text <- paste0(message_text, " [", default_text, "]")
    }
    message_text <- paste0(message_text, ": ")

    user_input <- readline(message_text)
    if (!nzchar(user_input) && !is.null(default_text)) {
      user_input <- default_text
    }

    parsed <- as.POSIXct(user_input, tz = time_zone, format = "%Y-%m-%d %H:%M:%S")
    if (!is.na(parsed)) {
      return(parsed)
    }

    message("Please enter the timestamp as YYYY-MM-DD HH:MM:SS")
  }
}

select_input_folder <- function(default_root) {
  chosen_dir <- NULL

  if (.Platform$OS.type == "windows" && interactive()) {
    chosen_dir <- tryCatch(
      utils::choose.dir(default = default_root, caption = "Select the CRDS raw-data folder for this run"),
      error = function(e) NULL
    )
  }

  if (is.null(chosen_dir) || is.na(chosen_dir) || !nzchar(chosen_dir)) {
    chosen_dir <- readline(
      paste0(
        "Enter raw-data folder path",
        if (nzchar(default_root)) paste0(" [", default_root, "]") else "",
        ": "
      )
    )

    if (!nzchar(chosen_dir)) {
      chosen_dir <- default_root
    }
  }

  normalizePath(chosen_dir, winslash = "/", mustWork = TRUE)
}

append_run_log <- function(log_path, run_record) {
  if (file.exists(log_path)) {
    existing_log <- suppressMessages(read_csv(log_path, show_col_types = FALSE))
  } else {
    existing_log <- tibble()
  }

  existing_log <- mutate(existing_log, across(everything(), as.character))
  run_record <- mutate(run_record, across(everything(), as.character))

  updated_log <- bind_rows(existing_log, run_record)
  write_csv(updated_log, log_path)
}

parse_dat_file_timestamp <- function(file_paths, time_zone = "Europe/Berlin") {
  file_names <- basename(file_paths)
  matched_text <- sub(".*-(\\d{8})-(\\d{6})-DataLog_User\\.dat$", "\\1 \\2", file_names, perl = TRUE)
  is_valid <- grepl("^\\d{8} \\d{6}$", matched_text)

  parsed <- as.POSIXct(rep(NA_character_, length(file_paths)), tz = time_zone)
  parsed[is_valid] <- as.POSIXct(
    matched_text[is_valid],
    tz = time_zone,
    format = "%Y%m%d %H%M%S"
  )

  parsed
}

find_relevant_dat_files <- function(input_path, start_time, end_time, time_zone = "Europe/Berlin") {
  dat_files <- list.files(
    path = input_path,
    recursive = TRUE,
    pattern = "\\.dat$",
    full.names = TRUE
  )

  if (length(dat_files) == 0) {
    return(character())
  }

  file_timestamps <- parse_dat_file_timestamp(dat_files, time_zone = time_zone)
  ordering <- order(file_timestamps, na.last = TRUE)
  dat_files <- dat_files[ordering]
  file_timestamps <- file_timestamps[ordering]

  valid_idx <- which(!is.na(file_timestamps))

  if (length(valid_idx) == 0) {
    return(dat_files)
  }

  matching_idx <- which(
    !is.na(file_timestamps) &
      file_timestamps >= (start_time - hours(2)) &
      file_timestamps <= (end_time + hours(2))
  )

  if (length(matching_idx) == 0) {
    matching_idx <- which(
      !is.na(file_timestamps) &
        file_timestamps >= floor_date(start_time, "day") &
        file_timestamps <= ceiling_date(end_time, "day")
    )
  }

  if (length(matching_idx) == 0) {
    return(dat_files)
  }

  candidate_idx <- unique(c(
    max(1, min(matching_idx) - 1),
    matching_idx,
    min(length(dat_files), max(matching_idx) + 1)
  ))

  dat_files[candidate_idx]
}

###### CRDS_cleaner function ######
crds_cleaner <- function(
    input_path,
    dat_files = NULL,
    gas = c("CO2", "CH4", "NH3", "H2O", "N2O"),
    start_time,
    end_time,
    flush,
    interval,
    MPVPosition.levels = NULL,
    location.levels = NULL,
    valid_locations = c("in", "S", "N"),
    lab = NULL,
    analyzer = NULL,
    time_zone = "Europe/Berlin",
    nh3_divisor = 1000,
    site_name = "Gross Kreutz") {

  required_cols <- c("DATE", "TIME", "MPVPosition")
  all_needed_cols <- unique(c(required_cols, gas))

  if (is.null(dat_files)) {
    dat_files <- list.files(
      path = input_path,
      recursive = TRUE,
      pattern = "\\.dat$",
      full.names = TRUE
    )
  }

  total_files <- length(dat_files)

  if (total_files == 0) {
    stop("No .dat files found in: ", input_path)
  }

  start_time <- as.POSIXct(start_time, tz = time_zone)
  end_time <- as.POSIXct(end_time, tz = time_zone)

  append_data <- function() {
    data_frames <- list()

    for (i in seq_along(dat_files)) {
      current_data <- tryCatch(
        read.table(dat_files[i], header = TRUE, stringsAsFactors = FALSE),
        error = function(e) NULL
      )

      if (is.null(current_data)) {
        warning("Could not read file: ", dat_files[i])
        next
      }

      existing_cols <- all_needed_cols[all_needed_cols %in% colnames(current_data)]

      if (length(existing_cols) == 0) {
        warning("No matching columns in file: ", dat_files[i])
        next
      }

      current_data <- current_data[, existing_cols, drop = FALSE]

      missing_cols <- setdiff(all_needed_cols, names(current_data))
      for (missing_col in missing_cols) {
        current_data[[missing_col]] <- NA
      }

      current_data <- current_data[, all_needed_cols, drop = FALSE]
      data_frames[[length(data_frames) + 1]] <- current_data

      cat(sprintf("Processed file %d of %d\n", i, total_files))
    }

    data_frames
  }

  list_of_data_frames <- append_data()

  if (length(list_of_data_frames) == 0) {
    stop("No usable data frames were created. Check input files and column names.")
  }

  raw_data <- bind_rows(list_of_data_frames)

  cat("Merging DATE and TIME column...\n")
  raw_data <- raw_data %>%
    mutate(
      DATE.TIME = as.POSIXct(
        paste(DATE, TIME),
        format = "%Y-%m-%d %H:%M:%S",
        tz = time_zone
      ),
      MPVPosition = suppressWarnings(as.numeric(MPVPosition))
    ) %>%
    select(-DATE, -TIME)

  if (all(is.na(raw_data$DATE.TIME))) {
    stop("DATE.TIME could not be parsed. Please check DATE and TIME format.")
  }

  for (gas_name in gas) {
    if (gas_name %in% names(raw_data)) {
      raw_data[[gas_name]] <- suppressWarnings(as.numeric(raw_data[[gas_name]]))
    }
  }

  cat("Filtering data between:\n",
      "Start:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n",
      "End:  ", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")

  cleaned_data <- raw_data %>%
    filter(
      DATE.TIME > start_time,
      DATE.TIME <= end_time,
      !is.na(MPVPosition),
      MPVPosition %% 1 == 0,
      MPVPosition != 0
    ) %>%
    arrange(DATE.TIME)

  if (!is.null(MPVPosition.levels)) {
    mpv_numeric <- suppressWarnings(as.numeric(MPVPosition.levels))

    cleaned_data <- cleaned_data %>%
      filter(MPVPosition %in% mpv_numeric)

    if (!is.null(location.levels)) {
      location_map <- data.frame(
        MPVPosition = mpv_numeric,
        location = location.levels,
        stringsAsFactors = FALSE
      )

      cleaned_data <- cleaned_data %>%
        left_join(location_map, by = "MPVPosition")
    }
  }

  if (!"location" %in% names(cleaned_data)) {
    cleaned_data$location <- NA_character_
  }

  if (nrow(cleaned_data) == 0) {
    stop("No rows left after filtering.")
  }

  cat("Cleaning MPVPosition and building step_id...\n")
  step_data <- cleaned_data %>%
    mutate(step_id = cumsum(c(TRUE, diff(MPVPosition) != 0))) %>%
    group_by(step_id, MPVPosition, location) %>%
    arrange(DATE.TIME, .by_group = TRUE) %>%
    mutate(
      timestamp_for_step = last(DATE.TIME),
      time_rank = row_number()
    ) %>%
    mutate(across(all_of(gas), ~ if_else(time_rank <= flush, NA_real_, .x))) %>%
    filter(time_rank <= interval) %>%
    summarise(
      DATE.TIME = last(timestamp_for_step),
      measuring.time = sum(time_rank > flush & time_rank <= interval),
      across(
        all_of(gas),
        ~ if (all(is.na(.x))) NA_real_ else mean(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    )

  if ("NH3" %in% names(step_data) && !is.null(nh3_divisor) && nh3_divisor != 1) {
    step_data <- step_data %>%
      mutate(NH3 = NH3 / nh3_divisor)
  }

  cat("Adding lab and analyzer info...\n")
  if (!is.null(lab)) step_data$lab <- lab
  if (!is.null(analyzer)) step_data$analyzer <- analyzer

  gases_present <- intersect(gas, names(step_data))

  hourly_long <- step_data %>%
    filter(location %in% valid_locations) %>%
    mutate(DATE.HOUR = floor_date(DATE.TIME, "hour")) %>%
    group_by(DATE.HOUR, location, lab, analyzer) %>%
    summarise(
      across(
        all_of(gases_present),
        ~ if (all(is.na(.x))) NA_real_ else mean(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    )

  if (nrow(hourly_long) > 0) {
    hourly_dummy <- expand_grid(
      DATE.HOUR = seq(
        min(hourly_long$DATE.HOUR),
        max(hourly_long$DATE.HOUR),
        by = "hour"
      ),
      analyzer = unique(hourly_long$analyzer),
      location = unique(hourly_long$location),
      lab = unique(hourly_long$lab)
    )

    hourly_wide <- hourly_dummy %>%
      left_join(
        hourly_long,
        by = c("DATE.HOUR", "location", "analyzer", "lab")
      ) %>%
      pivot_wider(
        names_from = location,
        values_from = all_of(gases_present),
        names_sep = "_"
      ) %>%
      arrange(DATE.HOUR)
  } else {
    hourly_wide <- data.frame()
  }

  make_final <- function(df_wide, suffix, site_code_label) {
    gases_export <- c("CO2", "CH4", "NH3", "H2O", "N2O")
    cols <- paste0(gases_export, "_", suffix)
    cols <- cols[cols %in% names(df_wide)]

    if (nrow(df_wide) == 0) {
      return(data.frame())
    }

    df_site <- df_wide %>%
      select(DATE.HOUR, analyzer, lab, all_of(cols)) %>%
      rename_with(~ gsub(paste0("_", suffix), "", .x), all_of(cols))

    for (gas_name in gases_export) {
      if (!gas_name %in% names(df_site)) df_site[[gas_name]] <- NA_real_
    }

    df_site %>%
      mutate(
        LOCATION = site_name,
        SITE_CODE = site_code_label,
        MEASUREMENT_PERIOD = "",
        DATE_TIME = DATE.HOUR,
        MEASUREMENT_POINT_CODE = lab,
        GAS_MEASUREMENT_METHOD = paste0("LINE_", analyzer),
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
      ) %>%
      select(
        LOCATION, SITE_CODE, MEASUREMENT_PERIOD, DATE_TIME,
        MEASUREMENT_POINT_CODE, GAS_MEASUREMENT_METHOD,
        CO2_DRY_PPM, CH4_DRY_PPM, NH3_DRY_PPM, H2O_VOL_PCT,
        N2O_DRY_PPM, NO2_DRY_PPM, NO_DRY_PPM, CO_DRY_PPM,
        COMMENTS, AVG_TOTAL_NUMBER_ANIMALS, AVG_WEIGHT_KG,
        AVG_DAYS_IN_PREGNANCY, AVG_MILK_YIELD_CORR_KGAD,
        IM_AIR_TEMPERATURE_C
      ) %>%
      arrange(DATE_TIME)
  }

  final_IN <- make_final(hourly_wide, "in", "IN")
  final_OUT_South <- make_final(hourly_wide, "S", "OUT_South")
  final_OUT_North <- make_final(hourly_wide, "N", "OUT_North")

  return(list(
    raw_data = raw_data,
    cleaned_data = cleaned_data,
    step_data = step_data,
    hourly_long = hourly_long,
    hourly_wide = hourly_wide,
    final_IN = final_IN,
    final_OUT_South = final_OUT_South,
    final_OUT_North = final_OUT_North
  ))
}

run_crds_period <- function(
    input_path,
    start_time,
    end_time,
    step_dir,
    hourly_dir,
    export_dir,
    time_zone = "Europe/Berlin") {
  start_stamp <- format(as.POSIXct(start_time, tz = time_zone), "%Y%m%d.%H%M")
  end_stamp <- format(as.POSIXct(end_time, tz = time_zone), "%Y%m%d.%H%M")
  start_date <- format(as.POSIXct(start_time, tz = time_zone), "%Y-%m-%d")
  end_date <- format(as.POSIXct(end_time, tz = time_zone), "%Y-%m-%d")

  step_file <- paste0(start_stamp, "-", end_stamp, "_ATB_240avg_CRDS8.csv")
  hourly_file <- paste0("H_CRDS8+9_", format(as.POSIXct(start_time, tz = time_zone), "%Y%m%d"), "_", format(as.POSIXct(end_time, tz = time_zone), "%Y%m%d"), ".csv")
  in_file <- paste0("GrossKreutz_IN_", start_date, "_to_", end_date, ".xlsx")
  out_south_file <- paste0("GrossKreutz_OUT_South_", start_date, "_to_", end_date, ".xlsx")
  out_north_file <- paste0("GrossKreutz_OUT_North_", start_date, "_to_", end_date, ".xlsx")
  selected_dat_files <- find_relevant_dat_files(
    input_path = input_path,
    start_time = as.POSIXct(start_time, tz = time_zone),
    end_time = as.POSIXct(end_time, tz = time_zone),
    time_zone = time_zone
  )

  if (length(selected_dat_files) == 0) {
    stop("No relevant .dat files found for the requested time range in: ", input_path)
  }

  message("Reading ", length(selected_dat_files), " relevant .dat files for the requested time range.")

  result <- crds_cleaner(
    input_path = input_path,
    dat_files = selected_dat_files,
    gas = c("CO2", "CH4", "NH3", "H2O", "N2O"),
    start_time = start_time,
    end_time = end_time,
    flush = 60,
    interval = 240,
    MPVPosition.levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9"),
    location.levels = c("1", "2", "3", "4", "5", "6", "7", "in", "S"),
    valid_locations = c("in", "S", "N"),
    lab = "ATB",
    analyzer = "CRDS8",
    time_zone = time_zone
  )

  step_export <- format_datetime_columns(result$step_data, "DATE.TIME", time_zone = time_zone)
  hourly_export <- format_datetime_columns(result$hourly_wide, "DATE.HOUR", time_zone = time_zone)
  final_in_export <- format_datetime_columns(result$final_IN, "DATE_TIME", time_zone = time_zone)
  final_out_south_export <- format_datetime_columns(result$final_OUT_South, "DATE_TIME", time_zone = time_zone)
  final_out_north_export <- format_datetime_columns(result$final_OUT_North, "DATE_TIME", time_zone = time_zone)

  write.csv(
    step_export,
    file.path(step_dir, step_file),
    row.names = FALSE,
    quote = FALSE
  )

  write_csv(
    hourly_export,
    file.path(hourly_dir, hourly_file)
  )

  write_xlsx(
    final_in_export,
    file.path(export_dir, in_file)
  )

  write_xlsx(
    final_out_south_export,
    file.path(export_dir, out_south_file)
  )

  write_xlsx(
    final_out_north_export,
    file.path(export_dir, out_north_file)
  )

  list(
    result = result,
    files = list(
      step_file = step_file,
      hourly_file = hourly_file,
      in_file = in_file,
      out_south_file = out_south_file,
      out_north_file = out_north_file
    )
  )
}

###### Weekly routine ######
project_dir <- resolve_project_dir()
workflow_dir <- file.path(project_dir, "workflows", "crds_routine_cleaning")
clean_dir <- file.path(workflow_dir, "clean_data", "crds_clean")
step_dir <- file.path(clean_dir, "step_avg")
hourly_dir <- file.path(clean_dir, "hourly_in_out_avg", "2026")
export_dir <- file.path(workflow_dir, "result_data", "reference_exports")
meta_dir <- file.path(workflow_dir, "meta_data")
log_path <- file.path(meta_dir, "cleaning_runs.csv")
default_raw_root <- "D:/Data_Analysis_R/owncloud_sync_data/CRDS08_raw"
reference_scripts <- c(
  file.path(workflow_dir, "scripts", "crds_weekly_cleaning.R"),
  file.path(project_dir, "scripts", "cleaning", "crds_cleaning_script.R")
)

dir.create(step_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(hourly_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(meta_dir, recursive = TRUE, showWarnings = FALSE)

last_end_time <- find_last_cleaned_end_time(
  log_path = log_path,
  step_dir = step_dir,
  reference_scripts = reference_scripts
)

if (!is.na(last_end_time)) {
  recommended_start_time <- last_end_time + 1
  message("Last cleaned end time found: ", format_datetime_text(last_end_time))
  message("Recommended next start_time to avoid overlap: ", format_datetime_text(recommended_start_time))
} else {
  recommended_start_time <- NULL
  message("No previous cleaning end time could be detected automatically.")
}

input_path <- select_input_folder(default_raw_root)
dat_count <- length(list.files(input_path, pattern = "\\.dat$", recursive = TRUE, full.names = TRUE))

if (dat_count == 0) {
  stop("No .dat files found in selected folder: ", input_path)
}

message("Selected raw-data folder: ", input_path)
message("Found ", dat_count, " .dat files.")
message("Tip: selecting a parent folder such as .../datalog_user/2026 is more robust than selecting a single month folder.")

start_time <- prompt_datetime(
  "Type start_time as YYYY-MM-DD HH:MM:SS",
  default_value = recommended_start_time
)
end_time <- prompt_datetime(
  "Type end_time as YYYY-MM-DD HH:MM:SS"
)

if (end_time <= start_time) {
  stop("end_time must be later than start_time.")
}

run_output <- run_crds_period(
  input_path = input_path,
  start_time = start_time,
  end_time = end_time,
  step_dir = step_dir,
  hourly_dir = hourly_dir,
  export_dir = export_dir
)

run_record <- tibble(
  run_timestamp = format_datetime_text(Sys.time()),
  selected_raw_folder = input_path,
  start_time = format_datetime_text(start_time),
  end_time = format_datetime_text(end_time),
  step_file = run_output$files$step_file,
  hourly_file = run_output$files$hourly_file,
  in_file = run_output$files$in_file,
  out_south_file = run_output$files$out_south_file,
  out_north_file = run_output$files$out_north_file
)

append_run_log(log_path, run_record)

message("Cleaning completed.")
message("Step output: ", file.path(step_dir, run_output$files$step_file))
message("Hourly output: ", file.path(hourly_dir, run_output$files$hourly_file))
message("Run logged in: ", log_path)
