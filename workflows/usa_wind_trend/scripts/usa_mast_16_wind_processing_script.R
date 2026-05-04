library(dplyr)
library(ggplot2)
library(grid)
library(lubridate)
library(readr)
library(tidyr)

# ------------------------------------------------------------
# Small helper functions
# ------------------------------------------------------------
# Calculate horizontal wind speed from the u and v components.
# We use only the horizontal wind here because wind roses are based on horizontal flow, while w is kept separately in the output data.
uv_to_ws <- function(u, v) {
  sqrt((u ^ 2) + (v ^ 2))
}

# Convert u and v into meteorological wind direction.
# The result is the direction the wind comes from, in degrees clockwise from north.
uv_to_wd <- function(u, v) {
  wd <- (270 - atan2(v, u) * 180 / pi) %% 360
  wd[uv_to_ws(u, v) == 0] <- NA_real_
  wd
}

# Group wind direction into 8 compass sectors for plotting.
deg_to_compass8 <- function(direction_deg) {
  labels <- c("N", "NE", "E", "SE", "S", "SW", "W", "NW")
  width <- 360 / length(labels)
  idx <- floor(((direction_deg %% 360) + width / 2) / width) %% length(labels) + 1
  labels[idx]
}

format_date_range <- function(start_date, end_date) {
  paste0(
    format(as.Date(start_date), "%d.%m.%Y"),
    " to ",
    format(as.Date(end_date), "%d.%m.%Y")
  )
}

# Build one windrose figure with several panels.
# The same frequency-ring scale is used in all panels so the figures are easy to compare across seasons or days.
build_faceted_windrose_plot <- function(
  data,
  facet_var,
  facet_nrow,
  facet_ncol,
  fill_values,
  radial_breaks
) {
  direction_levels <- c("N", "NE", "E", "SE", "S", "SW", "W", "NW")
  direction_axis_labels <- c(
    "N\n0°", "NE\n45°", "E\n90°", "SE\n135°",
    "S\n180°", "SW\n225°", "W\n270°", "NW\n315°"
  )
  speed_bin_levels <- c("0.0-0.5", "0.5-1.0", "1.0-2.0", "2.0-3.0", "3.0-4.0", ">=4.0")

  windrose_data <- data |>
    count({{ facet_var }}, direction_8, speed_bin) |>
    complete(
      {{ facet_var }},
      direction_8 = factor(direction_levels, levels = direction_levels),
      speed_bin = factor(speed_bin_levels, levels = speed_bin_levels),
      fill = list(n = 0)
    ) |>
    group_by({{ facet_var }}) |>
    mutate(percentage = if (sum(n) > 0) n / sum(n) * 100 else 0) |>
    ungroup() |>
    mutate(
      direction_8 = factor(as.character(direction_8), levels = direction_levels),
      speed_bin = factor(as.character(speed_bin), levels = speed_bin_levels),
      facet_label = as.character({{ facet_var }})
    )

  plot_obj <- ggplot(windrose_data, aes(x = direction_8, y = percentage, fill = speed_bin)) +
    geom_col(width = 1, color = "white", linewidth = 0.3) +
    coord_polar(start = -pi / 8) +
    facet_wrap(vars(facet_label), nrow = facet_nrow, ncol = facet_ncol) +
    scale_fill_manual(values = fill_values, drop = FALSE) +
    scale_y_continuous(
      breaks = radial_breaks,
      labels = function(x) paste0(x, "%"),
      limits = c(0, max(radial_breaks)),
      expand = expansion(mult = c(0, 0.02))
    ) +
    scale_x_discrete(labels = direction_axis_labels, drop = FALSE) +
    labs(
      x = NULL,
      y = "Frequency (%)",
      fill = expression(paste("Wind speed (m ", s^{-1}, ")"))
    ) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme_minimal(base_size = 18) +
    theme(
      strip.text = element_text(size = 16, face = "bold", lineheight = 1.25, margin = margin(b = 6)),
      axis.text.x = element_text(size = 12.5, face = "bold", color = "#333333"),
      axis.text.y = element_text(size = 11, color = "#5F5F5F"),
      axis.title.y = element_text(size = 15, face = "bold"),
      panel.background = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "#B3B3B3", linewidth = 0.35),
      panel.grid.major.y = element_line(color = "#B3B3B3", linewidth = 0.35),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 13),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.direction = "horizontal",
      legend.key.width = unit(1.2, "cm"),
      panel.spacing.x = unit(0.5, "lines"),
      panel.spacing.y = unit(2.0, "lines")
    )

  print(plot_obj)
  plot_obj
}

# ------------------------------------------------------------
# Paths and settings
# ------------------------------------------------------------
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args_all[grepl(file_arg, args_all)])

find_workflow_dir <- function() {
  candidate_starts <- unique(c(
    if (length(script_path) > 0) dirname(script_path[1]) else character(0),
    getwd(),
    "D:/Data_Analysis_R/Gas_Concentration_Emission_R/workflows/usa_wind_trend/scripts"
  ))

  for (start_dir in candidate_starts) {
    current_dir <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)

    for (i in seq_len(6)) {
      workflow_candidate <- if (basename(current_dir) == "scripts") {
        dirname(current_dir)
      } else if (basename(current_dir) == "usa_wind_trend") {
        current_dir
      } else if (basename(dirname(current_dir)) == "usa_wind_trend") {
        dirname(current_dir)
      } else {
        file.path(current_dir, "workflows", "usa_wind_trend")
      }

      raw_candidate <- file.path(
        workflow_candidate,
        "raw_data",
        "cloud.gr_2024_2025",
        "20240101_20250825_USA_mast_16_cloudgr_raw.csv"
      )

      if (file.exists(raw_candidate)) {
        return(normalizePath(workflow_candidate, winslash = "/", mustWork = FALSE))
      }

      parent_dir <- dirname(current_dir)
      if (identical(parent_dir, current_dir)) {
        break
      }
      current_dir <- parent_dir
    }
  }

  stop(
    "Could not locate the usa_wind_trend workflow folder. ",
    "Please run the script from within the project or keep the standard workflow structure."
  )
}

workflow_dir <- find_workflow_dir()
script_dir <- file.path(workflow_dir, "scripts")
project_dir <- normalizePath(file.path(workflow_dir, "..", ".."), winslash = "/", mustWork = FALSE)

timezone_local <- "Europe/Berlin"
node_label <- "USA_mast_16"

raw_file_path <- file.path(
  workflow_dir,
  "raw_data",
  "cloud.gr_2024_2025",
  "20240101_20250825_USA_mast_16_cloudgr_raw.csv"
)

clean_output_dir <- file.path(workflow_dir, "clean_data")
result_root_dir <- file.path(workflow_dir, "result_data")
result_dir <- file.path(result_root_dir, "plots")

dir.create(clean_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

rplots_pdf_path <- file.path(project_dir, "Rplots.pdf")
if (file.exists(rplots_pdf_path)) {
  unlink(rplots_pdf_path, force = TRUE)
}

if (!file.exists(raw_file_path)) {
  stop("Raw USA wind file not found: ", raw_file_path)
}

# ------------------------------------------------------------
# Data processing
# ------------------------------------------------------------
raw_data <- read_csv(
  raw_file_path,
  show_col_types = FALSE,
  col_types = cols(
    Date.Time = col_character(),
    node = col_character(),
    u_west_east = col_double(),
    v_south_north = col_double(),
    w_sky_soil = col_double()
  )
) |>
  mutate(
    datetime_utc = ymd_hms(Date.Time, tz = "UTC"),
    datetime_local = with_tz(datetime_utc, tzone = timezone_local),
    node = if_else(is.na(node) | node == "", node_label, node)
  ) |>
  filter(!is.na(datetime_local)) |>
  arrange(datetime_local)

# Aggregate component data to hourly means before deriving wind speed and direction.
# This keeps the output aligned with hourly summaries derived from the vector components.
hourly_data <- raw_data |>
  mutate(datetime_hour = floor_date(datetime_local, unit = "hour")) |>
  group_by(datetime_hour) |>
  summarise(
    node = dplyr::first(node),
    u = mean(u_west_east, na.rm = TRUE),
    v = mean(v_south_north, na.rm = TRUE),
    w = mean(w_sky_soil, na.rm = TRUE),
    observations_in_hour = n(),
    .groups = "drop"
  ) |>
  mutate(
    ws = uv_to_ws(u, v),
    wd = uv_to_wd(u, v),
    direction_8 = factor(deg_to_compass8(wd), levels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW")),
    speed_bin = cut(
      ws,
      breaks = c(0, 0.5, 1, 2, 3, 4, Inf),
      labels = c("0.0-0.5", "0.5-1.0", "1.0-2.0", "2.0-3.0", "3.0-4.0", ">=4.0"),
      include.lowest = TRUE,
      right = FALSE
    ),
    date_only = as.Date(datetime_hour, tz = timezone_local)
  )

hourly_export_path <- file.path(
  clean_output_dir,
  "20240101_20250825_USA_mast_16_hourly_uvw_wd_ws.csv"
)

write_csv(
  hourly_data |>
    transmute(
      datetime_hour = format(datetime_hour, "%Y-%m-%d %H:%M:%S"),
      node,
      u,
      v,
      w,
      wd,
      ws
    ),
  hourly_export_path
)

# ------------------------------------------------------------
# Plot subsets
# ------------------------------------------------------------
historical_periods <- tibble::tribble(
  ~season_name, ~start_date, ~end_date,
  "Summer 2024", as.Date("2024-06-01"), as.Date("2024-08-31"),
  "Autumn 2024", as.Date("2024-09-01"), as.Date("2024-11-30"),
  "Winter 2024/25", as.Date("2024-12-01"), as.Date("2025-02-28"),
  "Spring 2025", as.Date("2025-03-01"), as.Date("2025-04-07")
)

measurement_dates <- seq.Date(as.Date("2025-04-08"), as.Date("2025-04-15"), by = "day")

# Palette categorical contrast
windrose_colors <- c(
  "0.0-0.5" = "#9E0142",
  "0.5-1.0" = "#D53E4F",
  "1.0-2.0" = "#FDAE61",
  "2.0-3.0" = "#E6F598",
  "3.0-4.0" = "#66C2A5",
  ">=4.0" = "#3288BD"
)

historical_data <- dplyr::bind_rows(lapply(seq_len(nrow(historical_periods)), function(i) {
  period_row <- historical_periods[i, ]
  season_name <- as.character(period_row$season_name[[1]])
  start_date <- as.Date(period_row$start_date[[1]], origin = "1970-01-01")
  end_date <- as.Date(period_row$end_date[[1]], origin = "1970-01-01")

  hourly_data |>
    filter(date_only >= start_date, date_only <= end_date) |>
    mutate(season_label = paste0(season_name, "\n", format_date_range(start_date, end_date)))
}))

seasonal_plot <- build_faceted_windrose_plot(
  data = historical_data,
  facet_var = season_label,
  facet_nrow = 1,
  facet_ncol = 4,
  fill_values = windrose_colors,
  radial_breaks = c(10, 20, 30, 40)
)

ggsave(
  filename = file.path(result_dir, "historical_seasonal_windroses_summer2024_to_spring2025.png"),
  plot = seasonal_plot,
  width = 18.5,
  height = 6.4,
  units = "in",
  dpi = 500,
  bg = "white"
)

daily_panel_data <- dplyr::bind_rows(lapply(measurement_dates, function(measurement_date) {
  measurement_date <- as.Date(measurement_date, origin = "1970-01-01")
  hourly_data |>
    filter(date_only == measurement_date) |>
    mutate(day_label = format(measurement_date, "%d.%m.%Y"))
}))

daily_plot <- build_faceted_windrose_plot(
  data = daily_panel_data,
  facet_var = day_label,
  facet_nrow = 2,
  facet_ncol = 4,
  fill_values = windrose_colors,
  radial_breaks = c(20, 40, 60, 80)
)

ggsave(
  filename = file.path(result_dir, "ringversuche_daily_windroses_20250408_20250415.png"),
  plot = daily_plot,
  width = 15.8,
  height = 8.8,
  units = "in",
  dpi = 500,
  bg = "white"
)
