resolve_project_dir <- function(require_crds = FALSE, require_utils = FALSE) {
  project_marker <- function(path) {
    has_device_comparison <- dir.exists(file.path(path, "workflows", "device_comparison"))
    has_crds <- !require_crds || dir.exists(file.path(path, "workflows", "crds_routine_cleaning"))
    has_utils <- !require_utils || dir.exists(file.path(path, "scripts", "utils"))
    has_device_comparison && has_crds && has_utils
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
    active_path <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
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

campaign_start_default <- as.POSIXct("2025-09-29 00:00:00", tz = "Europe/Berlin")

device_colors <- c(
  crds = "#4B5563",
  logas_ndir = "#7570B3",
  logas_tdlas = "#1B9E77",
  otice = "#C45A11"
)

device_labels <- c(
  crds = "CRDS",
  logas_ndir = "PRONOVA",
  logas_tdlas = "CUBIC",
  otice = "OTICE",
  otice_raw = "OTICE raw"
)

parse_mixed_datetime <- function(x, time_zone = "Europe/Berlin") {
  x_chr <- as.character(x)
  parsed <- as.POSIXct(rep(NA_character_, length(x_chr)), tz = time_zone)
  formats <- c(
    "%Y-%m-%dT%H:%M:%SZ",
    "%Y-%m-%d %H:%M:%S",
    "%Y/%m/%d %H:%M:%S",
    "%Y-%m-%dT%H:%M:%S"
  )

  for (fmt in formats) {
    idx <- which(is.na(parsed) & !is.na(x_chr))
    if (length(idx) == 0) break
    trial <- as.POSIXct(x_chr[idx], format = fmt, tz = time_zone)
    parsed[idx[!is.na(trial)]] <- trial[!is.na(trial)]
  }

  parsed
}

parse_datetime_local <- function(x, time_zone = "Europe/Berlin") {
  lubridate::parse_date_time(
    as.character(x),
    orders = c(
      "ymd HMS", "ymd HM",
      "Ymd HMS", "Ymd HM",
      "Y/m/d HMS", "Y/m/d HM",
      "dmy HMS", "dmy HM",
      "mdy HMS", "mdy HM"
    ),
    tz = time_zone,
    quiet = TRUE
  )
}

safe_mean <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_rmse <- function(obs, pred) {
  sqrt(mean((pred - obs)^2, na.rm = TRUE))
}

safe_mae <- function(obs, pred) {
  mean(abs(pred - obs), na.rm = TRUE)
}

safe_r2 <- function(obs, pred) {
  if (length(obs) < 2 || sd(obs, na.rm = TRUE) == 0 || sd(pred, na.rm = TRUE) == 0) {
    return(NA_real_)
  }
  cor(obs, pred, use = "complete.obs")^2
}

label_analyzer <- function(x) {
  dplyr::recode(x, !!!device_labels, .default = x)
}

device_plot_theme_classic <- function(
  base_size = 13,
  title_size = 15,
  x_text_size = 10,
  y_text_size = 10,
  strip_size = 11,
  legend_text_size = 10
) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = legend_text_size),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = title_size),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = x_text_size),
      axis.text.y = ggplot2::element_text(size = y_text_size),
      strip.text = ggplot2::element_text(size = strip_size),
      strip.text.y.left = ggplot2::element_text(size = strip_size),
      panel.border = ggplot2::element_rect(color = "black", fill = NA)
    )
}

device_plot_theme_minimal <- function(
  base_size = 11,
  title_size = 14,
  subtitle_size = 10,
  strip_size = 10,
  legend_text_size = 10
) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = title_size),
      plot.subtitle = ggplot2::element_text(size = subtitle_size),
      strip.text = ggplot2::element_text(face = "bold", size = strip_size),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = legend_text_size),
      panel.grid.minor = ggplot2::element_blank()
    )
}
