# Combine CRDS8, CRDS9, FTIR1, and FTIR2 over their common October 2024
# measurement period and remap analyser-specific locations.

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package 'data.table' is required.")
}

library(data.table)

workflow_dir <- paste0(
  "D:/Data_Analysis_R/Gas_Concentration_Emission_R/",
  "workflows/mapping_high_resolution"
)
clean_dir <- file.path(workflow_dir, "clean_data")

input_files <- c(
  CRDS8 = "20240516.172700_20241029.132912_CRDS8_240s_step_average.csv",
  CRDS9 = "20240516.172700_20241029.132912_CRDS9_240s_step_average.csv",
  FTIR1 = "2024_05_17_FTIR1_clean.csv",
  FTIR2 = "2024_05_17_FTIR2_clean.csv"
)

required_columns <- c(
  "DATE.TIME", "analyser", "location", "CO2", "CH4", "NH3", "H2O"
)

read_clean_data <- function(analyser_name, file_name) {
  path <- file.path(clean_dir, file_name)
  data <- fread(
    path,
    colClasses = list(character = c("DATE.TIME", "analyser", "location"))
  )

  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0L) {
    stop(
      basename(path), " is missing: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  data <- data[, ..required_columns]
  data[, DATE.TIME := as.POSIXct(
    DATE.TIME,
    format = "%Y-%m-%d %H:%M:%S",
    tz = "Europe/Berlin"
  )]
  data[, analyser := analyser_name]
  data
}

data_list <- Map(read_clean_data, names(input_files), unname(input_files))
names(data_list) <- names(input_files)

# First restrict every analyser to October 2024.
october_start <- as.POSIXct("2024-10-01 00:00:00", tz = "Europe/Berlin")
november_start <- as.POSIXct("2024-11-01 00:00:00", tz = "Europe/Berlin")
data_list <- lapply(
  data_list,
  function(data) data[
    DATE.TIME >= october_start & DATE.TIME < november_start
  ]
)

if (any(vapply(data_list, nrow, integer(1)) == 0L)) {
  empty_analysers <- names(data_list)[
    vapply(data_list, nrow, integer(1)) == 0L
  ]
  stop(
    "No October 2024 data for: ",
    paste(empty_analysers, collapse = ", ")
  )
}

# Use the strict intersection of the four observed October time ranges.
common_start <- max(vapply(
  data_list,
  function(data) as.numeric(min(data[["DATE.TIME"]])),
  numeric(1)
))
common_end <- min(vapply(
  data_list,
  function(data) as.numeric(max(data[["DATE.TIME"]])),
  numeric(1)
))
common_start <- as.POSIXct(
  common_start,
  origin = "1970-01-01",
  tz = "Europe/Berlin"
)
common_end <- as.POSIXct(
  common_end,
  origin = "1970-01-01",
  tz = "Europe/Berlin"
)

data_list <- lapply(
  data_list,
  function(data) data[
    DATE.TIME >= common_start & DATE.TIME <= common_end
  ]
)

location_maps <- list(
  CRDS9 = setNames(as.character(1:16), as.character(1:16)),
  FTIR2 = setNames(as.character(17:26), as.character(1:10)),
  CRDS8 = setNames(as.character(27:42), as.character(1:16)),
  FTIR1 = setNames(as.character(c(43:51, "s")), as.character(1:10))
)

for (analyser_name in names(data_list)) {
  data <- data_list[[analyser_name]]
  source_location <- as.character(data[["location"]])
  mapped_location <- unname(location_maps[[analyser_name]][source_location])

  if (anyNA(mapped_location)) {
    unmapped <- sort(unique(source_location[is.na(mapped_location)]))
    stop(
      "Unmapped ", analyser_name, " location(s): ",
      paste(unmapped, collapse = ", ")
    )
  }

  data[, location := mapped_location]
  data_list[[analyser_name]] <- data
}

combined <- rbindlist(data_list, use.names = TRUE)
combined[, location := factor(
  location,
  levels = c(as.character(1:51), "s")
)]

# Mid-level tubes were removed from all four devices. Exclude those numbered
# barn locations while retaining the separate "s" sampling location.
discarded_mid_locations <- c(
  2, 5, 8, 11, 14, 17, 20, 23, 26,
  29, 32, 35, 38, 41, 44, 47, 50
)
remaining_numbered_locations <- setdiff(1:51, discarded_mid_locations)
combined <- combined[
  !as.character(location) %in% as.character(discarded_mid_locations)
]

setorder(combined, DATE.TIME, location, analyser)
setcolorder(combined, required_columns)

output_file <- file.path(
  clean_dir,
  paste0(
    format(common_start, "%Y%m%d.%H%M%S"),
    "_",
    format(common_end, "%Y%m%d.%H%M%S"),
    "_FTIR1_FTIR2_CRDS8_CRDS9_combined.csv"
  )
)

fwrite(
  combined,
  output_file,
  quote = FALSE,
  dateTimeAs = "write.csv"
)

message(
  "Common overlap: ",
  format(common_start, "%Y-%m-%d %H:%M:%S"),
  " to ",
  format(common_end, "%Y-%m-%d %H:%M:%S")
)
message("Wrote ", nrow(combined), " combined rows to: ", output_file)

# Plot CO2, CH4, and NH3 by remapped barn location. Location "s" is excluded.
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required to create the concentration plot.")
}

vertical_groups <- list(
  top = c(1, 4, 7, 10, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 46, 49),
  mid = c(2, 5, 8, 11, 14, 17, 20, 23, 26, 29, 32, 35, 38, 41, 44, 47, 50),
  bottom = c(3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51)
)

plot_data <- copy(combined[as.character(location) != "s"])
plot_data[, location_number := as.integer(as.character(location))]
plot_data[, vertical_group := fcase(
  location_number %in% vertical_groups$top, "top",
  location_number %in% vertical_groups$mid, "mid",
  location_number %in% vertical_groups$bottom, "bottom",
  default = NA_character_
)]

if (anyNA(plot_data[["vertical_group"]])) {
  stop("One or more plotted locations were not assigned to a vertical group.")
}

plot_long <- melt(
  plot_data,
  id.vars = c("DATE.TIME", "analyser", "location_number", "vertical_group"),
  measure.vars = c("CO2", "CH4", "NH3"),
  variable.name = "gas",
  value.name = "concentration"
)
plot_long[, location := factor(location_number, levels = 1:51)]
plot_long[, vertical_group := factor(
  vertical_group,
  levels = c("top", "mid", "bottom")
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
  ggplot2::geom_point(
    size = 1.9
  ) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(gas),
    scales = "free_y",
    switch = "y"
  ) +
  ggplot2::scale_colour_manual(
    values = c(
      "top" = "#0072B2",
      "mid" = "#E69F00",
      "bottom" = "#009E73"
    )
  ) +
  ggplot2::labs(
    title = "Gas concentrations by barn location",
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
    axis.text.x = ggplot2::element_text(
      angle = 0,
      hjust = 0.5,
      size = 7
    )
  )

plot_dir <- file.path(workflow_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
plot_file <- file.path(
  plot_dir,
  "CO2_CH4_NH3_by_location_vertical_groups.png"
)
ggplot2::ggsave(
  filename = plot_file,
  plot = concentration_plot,
  width = 15,
  height = 11,
  units = "in",
  dpi = 300,
  bg = "white"
)
message("Plot saved to: ", plot_file)
