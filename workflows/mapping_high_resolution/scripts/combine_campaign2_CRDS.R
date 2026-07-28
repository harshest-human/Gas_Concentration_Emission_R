# Combine the mapped CRDS8 and CRDS9 step averages for campaign 2 and create
# a mean +/- SE concentration plot by barn location.

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
clean_dir <- file.path(workflow_dir, "clean_data")
plot_dir <- file.path(workflow_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

input_files <- c(
  CRDS8 = "20241116.000000_20241231.235959_CRDS8_240s_step_average.csv",
  CRDS9 = "20241116.000000_20241231.235959_CRDS9_240s_step_average.csv"
)
required_columns <- c(
  "DATE.TIME", "analyser", "location", "CO2", "CH4", "NH3", "H2O"
)

read_one <- function(analyser_name, file_name) {
  data <- fread(
    file.path(clean_dir, file_name),
    colClasses = list(character = c("DATE.TIME", "analyser", "location"))
  )
  if (!identical(names(data), required_columns)) {
    stop(file_name, " does not have the expected columns.")
  }
  data[, DATE.TIME := as.POSIXct(
    DATE.TIME,
    format = "%Y-%m-%d %H:%M:%S",
    tz = "Europe/Berlin"
  )]
  data[, analyser := analyser_name]
  data
}

data_list <- Map(read_one, names(input_files), unname(input_files))
names(data_list) <- names(input_files)

# Use the strict intersection of both observed campaign ranges.
common_start <- as.POSIXct(
  max(vapply(
    data_list,
    function(data) as.numeric(min(data[["DATE.TIME"]])),
    numeric(1)
  )),
  origin = "1970-01-01",
  tz = "Europe/Berlin"
)
common_end <- as.POSIXct(
  min(vapply(
    data_list,
    function(data) as.numeric(max(data[["DATE.TIME"]])),
    numeric(1)
  )),
  origin = "1970-01-01",
  tz = "Europe/Berlin"
)

data_list <- lapply(
  data_list,
  function(data) data[
    DATE.TIME >= common_start & DATE.TIME <= common_end
  ]
)

combined <- rbindlist(data_list, use.names = TRUE)
combined[, location := as.character(location)]
setorder(combined, DATE.TIME, analyser, location)
setcolorder(combined, required_columns)

output_file <- file.path(
  clean_dir,
  paste0(
    format(common_start, "%Y%m%d.%H%M%S"),
    "_",
    format(common_end, "%Y%m%d.%H%M%S"),
    "_CRDS8_CRDS9_campaign2_combined.csv"
  )
)
fwrite(
  combined,
  output_file,
  quote = FALSE,
  dateTimeAs = "write.csv"
)

vertical_groups <- list(
  top = c(1, 4, 7, 10, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 46, 49),
  bottom = c(3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51)
)

plot_data <- copy(combined)
plot_data[, location_number := as.integer(location)]
plot_data[, vertical_group := fcase(
  location_number %in% vertical_groups$top, "top",
  location_number %in% vertical_groups$bottom, "bottom",
  default = NA_character_
)]

if (anyNA(plot_data[["vertical_group"]])) {
  stop("One or more locations were not assigned to top or bottom.")
}

plot_long <- melt(
  plot_data,
  id.vars = c("DATE.TIME", "analyser", "location_number", "vertical_group"),
  measure.vars = c("CO2", "CH4", "NH3"),
  variable.name = "gas",
  value.name = "concentration"
)

location_levels <- sort(unique(plot_long[["location_number"]]))
plot_long[, location := factor(location_number, levels = location_levels)]
plot_long[, vertical_group := factor(
  vertical_group,
  levels = c("top", "bottom")
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
  by = .(gas, location, vertical_group)
]

concentration_plot <- ggplot2::ggplot(
  plot_summary,
  ggplot2::aes(
    x = location,
    y = mean_concentration,
    colour = vertical_group
  )
) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = mean_concentration - se_concentration,
      ymax = mean_concentration + se_concentration
    ),
    width = 0.35,
    linewidth = 0.55
  ) +
  ggplot2::geom_point(size = 2) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(gas),
    scales = "free_y",
    switch = "y"
  ) +
  ggplot2::scale_colour_manual(
    values = c("top" = "#0072B2", "bottom" = "#009E73")
  ) +
  ggplot2::labs(
    title = "Campaign 2 CRDS gas concentrations by barn location",
    subtitle = paste0(
      format(common_start, "%Y-%m-%d %H:%M:%S"),
      " to ",
      format(common_end, "%Y-%m-%d %H:%M:%S")
    ),
    x = "Location",
    y = "Concentration",
    colour = "Vertical group"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    legend.position = "top",
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
  "campaign2_CRDS8_CRDS9_CO2_CH4_NH3_mean_SE_by_location.png"
)
ggplot2::ggsave(
  filename = plot_file,
  plot = concentration_plot,
  width = 14,
  height = 11,
  units = "in",
  dpi = 300,
  bg = "white"
)

message(
  "Common overlap: ",
  format(common_start, "%Y-%m-%d %H:%M:%S"),
  " to ",
  format(common_end, "%Y-%m-%d %H:%M:%S")
)
message("Wrote ", nrow(combined), " combined rows to: ", output_file)
message("Plot saved to: ", plot_file)
