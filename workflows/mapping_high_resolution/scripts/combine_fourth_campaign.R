# Combine fourth-campaign CRDS8 files, remap MPVPositions 1:7 to barn
# locations, remove overlapping exports, and create a mean +/- SE plot.

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package 'data.table' is required.")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required.")
}

library(data.table)

workflow_dir <- paste0(
  "D:/Data_Analysis_R/Gas_Concentration_Emission_R/",
  "workflows/mapping_high_resolution"
)
input_dir <- file.path(workflow_dir, "clean_data", "4_campaign")
plot_dir <- file.path(workflow_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

required_input <- c(
  "DATE.TIME", "analyzer", "MPVPosition", "CO2", "CH4", "NH3", "H2O"
)
output_columns <- c(
  "DATE.TIME", "analyser", "location", "CO2", "CH4", "NH3", "H2O"
)

input_files <- list.files(
  input_dir,
  pattern = "\\.csv$",
  full.names = TRUE,
  ignore.case = TRUE
)
input_files <- input_files[
  !grepl("_fourth_campaign_combined\\.csv$", basename(input_files))
]

if (length(input_files) == 0L) {
  stop("No fourth-campaign CSV files found in: ", input_dir)
}

location_map <- c(
  "1" = 12,
  "2" = 18,
  "3" = 15,
  "4" = 30,
  "5" = 36,
  "6" = 42,
  "7" = 45
)

read_one <- function(path) {
  available <- names(fread(path, nrows = 0, showProgress = FALSE))
  missing_columns <- setdiff(required_input, available)
  if (length(missing_columns) > 0L) {
    stop(
      basename(path), " is missing: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  data <- fread(
    path,
    select = required_input,
    colClasses = list(character = c("DATE.TIME", "analyzer")),
    showProgress = FALSE
  )
  data <- data[MPVPosition %in% 1:7]
  data[, location := as.character(
    unname(location_map[as.character(MPVPosition)])
  )]
  data[, DATE.TIME := as.POSIXct(
    DATE.TIME,
    format = "%Y-%m-%d %H:%M:%S",
    tz = "Europe/Berlin"
  )]
  setnames(data, "analyzer", "analyser")
  data[, MPVPosition := NULL]
  setcolorder(data, output_columns)
  data
}

combined_raw <- rbindlist(lapply(input_files, read_one), use.names = TRUE)
raw_rows <- nrow(combined_raw)
setorder(combined_raw, DATE.TIME, analyser, location)

# Overlapping source exports contain the same measurement periods. Keep one
# record per timestamp, analyser, and mapped location.
combined <- unique(
  combined_raw,
  by = c("DATE.TIME", "analyser", "location")
)
duplicate_rows_removed <- raw_rows - nrow(combined)
setcolorder(combined, output_columns)

output_file <- file.path(
  input_dir,
  paste0(
    format(min(combined[["DATE.TIME"]]), "%Y%m%d.%H%M%S"),
    "_",
    format(max(combined[["DATE.TIME"]]), "%Y%m%d.%H%M%S"),
    "_fourth_campaign_combined.csv"
  )
)
fwrite(
  combined,
  output_file,
  quote = FALSE,
  dateTimeAs = "write.csv"
)

plot_long <- melt(
  combined,
  id.vars = c("DATE.TIME", "analyser", "location"),
  measure.vars = c("CO2", "CH4", "NH3"),
  variable.name = "gas",
  value.name = "concentration"
)
location_levels <- sort(unique(as.integer(plot_long[["location"]])))
plot_long[, location := factor(
  as.integer(location),
  levels = location_levels
)]

plot_summary <- plot_long[
  !is.na(concentration),
  .(
    mean_concentration = mean(concentration),
    se_concentration = if (.N > 1L) {
      stats::sd(concentration) / sqrt(.N)
    } else {
      NA_real_
    }
  ),
  by = .(gas, location)
]

concentration_plot <- ggplot2::ggplot(
  plot_summary,
  ggplot2::aes(x = location, y = mean_concentration)
) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = mean_concentration - se_concentration,
      ymax = mean_concentration + se_concentration
    ),
    width = 0.3,
    linewidth = 0.6,
    colour = "#009E73"
  ) +
  ggplot2::geom_point(size = 2.2, colour = "#009E73") +
  ggplot2::facet_grid(
    rows = ggplot2::vars(gas),
    scales = "free_y",
    switch = "y"
  ) +
  ggplot2::labs(
    title = "Fourth campaign gas concentrations by barn location",
    subtitle = paste0(
      format(min(combined[["DATE.TIME"]]), "%Y-%m-%d %H:%M:%S"),
      " to ",
      format(max(combined[["DATE.TIME"]]), "%Y-%m-%d %H:%M:%S")
    ),
    x = "Location",
    y = "Concentration"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(
      colour = "grey35",
      fill = NA,
      linewidth = 0.7
    ),
    panel.spacing.y = grid::unit(1.2, "lines"),
    strip.placement = "outside",
    strip.background = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5)
  )

plot_file <- file.path(
  plot_dir,
  "fourth_campaign_CO2_CH4_NH3_mean_SE_by_location.png"
)
ggplot2::ggsave(
  filename = plot_file,
  plot = concentration_plot,
  width = 11,
  height = 11,
  units = "in",
  dpi = 300,
  bg = "white"
)

message("Read ", length(input_files), " fourth-campaign files")
message("Removed ", duplicate_rows_removed, " overlapping row(s)")
message("Wrote ", nrow(combined), " combined rows to: ", output_file)
message("Plot saved to: ", plot_file)
