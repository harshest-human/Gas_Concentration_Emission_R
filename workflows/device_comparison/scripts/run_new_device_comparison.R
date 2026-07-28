suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(lubridate); library(purrr)
  library(readr); library(readxl); library(tidyr); library(tibble)
})

script_dir <- {
  arg <- commandArgs(trailingOnly = FALSE)
  path <- sub("--file=", "", arg[grepl("--file=", arg)])
  if (length(path)) dirname(normalizePath(path[1], winslash = "/")) else normalizePath(getwd(), winslash = "/")
}
r_dir <- normalizePath(file.path(script_dir, "..", "R"), winslash = "/", mustWork = TRUE)
invisible(lapply(c("config.R", "utils.R", "process_devices.R", "load_crds.R", "integrate.R", "plot.R"), function(x) source(file.path(r_dir, x))))

parse_args <- function(args) {
  values <- list(start_time = NULL, end_time = NULL, devices = c("crds", "cubic", "pronova", "otice"), emissions = TRUE)
  for (arg in args) {
    bits <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    if (length(bits) == 2 && bits[1] %in% names(values)) values[[bits[1]]] <- bits[2]
  }
  values$devices <- strsplit(values$devices, ",", fixed = TRUE)[[1]]
  values$emissions <- tolower(as.character(values$emissions)) %in% c("true", "1", "yes")
  values
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(args$start_time) || is.null(args$end_time)) {
  stop("Usage: Rscript run_new_device_comparison.R --start_time='2026-07-01 00:00:00' --end_time='2026-07-10 00:00:00' [--devices=crds,cubic,pronova,otice] [--emissions=true]", call. = FALSE)
}

project_dir <- normalizePath(file.path(script_dir, "..", "..", ".."), winslash = "/", mustWork = TRUE)
config <- device_comparison_config(project_dir)
window <- parse_window(args$start_time, args$end_time, config$timezone)
processors <- list(crds = load_crds_hourly, cubic = process_cubic_hourly, pronova = process_pronova_hourly, otice = process_otice_hourly)
unknown <- setdiff(args$devices, names(processors))
if (length(unknown)) stop("Unknown devices: ", paste(unknown, collapse = ", "), call. = FALSE)

hourly <- purrr::map_dfr(args$devices, function(device) {
  message("Processing ", device, "...")
  data <- processors[[device]](window$start, window$end, config)
  path <- write_hourly(data, device, window$start, window$end, config)
  message("  wrote ", path)
  data
})

main_df <- build_main_df(hourly) |> add_context_data(config)
final_df <- if (args$emissions) calculate_emissions(main_df, config$project_dir) else main_df
dir.create(config$integrated, recursive = TRUE, showWarnings = FALSE)
tag <- paste0(format(window$start, "%Y%m%d%H"), "_", format(window$end, "%Y%m%d%H"))
main_path <- file.path(config$integrated, paste0("device_comparison_", tag, ".csv"))
write_csv(final_df, main_path, na = "")
message("Wrote integrated dataset: ", main_path)
