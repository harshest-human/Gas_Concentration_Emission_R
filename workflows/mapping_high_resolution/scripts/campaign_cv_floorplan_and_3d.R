###############################################################################
# Campaign CV maps and interactive 3D supplement
#
# Outputs:
#   - CV summary and reusable spatial-coordinate CSV files
#   - one six-panel floor-plan CV map for each campaign
#   - one interactive 3D HTML file for all campaigns
#   - one interactive 3D HTML file for each campaign
#
# The cleaned campaign files are read but never modified.
###############################################################################

required_packages <- c(
  "data.table", "ggplot2", "png", "plotly", "htmlwidgets"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(plotly)
})


###############################################################################
##### PATHS
###############################################################################

workflow_dir <- paste0(
  "D:/Data_Analysis_R/Gas_Concentration_Emission_R/",
  "workflows/mapping_high_resolution"
)
campaign_dir <- file.path(workflow_dir, "clean_data", "campaign_combined")
floorplan_file <- file.path(
  workflow_dir,
  "Manuscript_3_Mapping_high_resolution_concentration_Latex_draft",
  "figures",
  "Fig1a_campaign1_sampling_layout_no_ringline.png"
)
output_data_dir <- file.path(
  workflow_dir,
  "clean_data",
  "campaign_cv_spatial"
)
output_plot_dir <- file.path(
  workflow_dir,
  "plots",
  "campaign_cv_spatial"
)
output_html_dir <- file.path(output_plot_dir, "interactive_3d")

dir.create(output_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_html_dir, recursive = TRUE, showWarnings = FALSE)

# htmlwidgets uses Pandoc to create portable self-contained HTML files. Use the
# Pandoc bundled with RStudio when it is not already available on PATH.
rstudio_pandoc <- paste0(
  "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
)
if (
  !nzchar(Sys.which("pandoc")) &&
    file.exists(file.path(rstudio_pandoc, "pandoc.exe"))
) {
  Sys.setenv(RSTUDIO_PANDOC = rstudio_pandoc)
}

campaign_files <- c(
  "Campaign 1" = "campaign_1_20241001.000232_20241024.143858.csv",
  "Campaign 2" = "campaign_2_20241116.000647_20241231.235710.csv",
  "Campaign 3" = "campaign_3_20250819.130012_20251119.090642.csv",
  "Campaign 4" = "campaign_4_20251209.120359_20260721.023820.csv"
)

if (!file.exists(floorplan_file)) {
  stop("Floor-plan image not found: ", floorplan_file)
}


###############################################################################
##### REUSABLE LOCATION COORDINATES
###############################################################################

# Pixel centres were registered against the 1739 x 904 manuscript floor plan.
# Each horizontal position contains top, middle, and bottom sampling locations.
# The slightly separated map x-coordinate keeps the three height markers visible
# on the two-dimensional floor plan. The 3D x/y coordinates use the common
# centre of each vertical sampling group.
group_coordinates <- data.table(
  horizontal_position = 1:17,
  group_x_pixel = c(
    1043, 1308, 1575, 1575, 1308, 1043, 775, 775, 1043,
    1308, 1575, 1575, 1308, 1043, 775, 510, 510
  ),
  group_y_pixel = c(
    663, 663, 663, 498, 498, 498, 498, 339, 339,
    339, 339, 171, 171, 171, 171, 339, 171
  )
)

coordinates <- data.table(location_number = 1:51)
coordinates[, `:=`(
  horizontal_position = ((location_number - 1L) %/% 3L) + 1L,
  height = c("top", "middle", "bottom")[
    ((location_number - 1L) %% 3L) + 1L
  ]
)]
coordinates <- merge(
  coordinates,
  group_coordinates,
  by = "horizontal_position",
  all.x = TRUE,
  sort = FALSE
)
coordinates[, map_x_pixel := group_x_pixel + fifelse(
  height == "top", -28,
  fifelse(height == "middle", 0, 28)
)]
coordinates[, map_y_pixel := group_y_pixel]

# Approximate metre coordinates are obtained by scaling the dimensioned plan:
# group x = 510...1575 pixels corresponds to 0...38.8 m;
# group y = 663...171 pixels corresponds to 0...17.6 m.
coordinates[, `:=`(
  x_m = (group_x_pixel - 510) / (1575 - 510) * 38.8,
  y_m = (663 - group_y_pixel) / (663 - 171) * 17.6,
  z_layer = c(top = 3, middle = 2, bottom = 1)[height],
  z_label = factor(
    height,
    levels = c("bottom", "middle", "top")
  ),
  coordinate_basis = "Scaled from the dimensioned manuscript floor plan"
)]
setorder(coordinates, location_number)
fwrite(
  coordinates,
  file.path(output_data_dir, "location_coordinate_lookup.csv")
)


###############################################################################
##### CALCULATE CAMPAIGN-LEVEL MEAN, SD, AND CV
###############################################################################

campaign_data <- rbindlist(
  lapply(names(campaign_files), function(campaign_name) {
    file_path <- file.path(campaign_dir, campaign_files[[campaign_name]])
    if (!file.exists(file_path)) {
      stop("Campaign file not found: ", file_path)
    }
    x <- fread(file_path, na.strings = c("", "NA", "NaN"))
    x[, campaign := campaign_name]
    x
  }),
  use.names = TRUE,
  fill = TRUE
)

campaign_data[, location := trimws(as.character(location))]
for (gas in c("CO2", "CH4", "NH3")) {
  campaign_data[, (gas) := suppressWarnings(as.numeric(get(gas)))]
  campaign_data[!is.finite(get(gas)) | get(gas) <= 0, (gas) := NA_real_]
}

campaign_data[, `:=`(
  CH4_CO2 = fifelse(CO2 > 0 & CH4 > 0, CH4 / CO2, NA_real_),
  NH3_CO2 = fifelse(CO2 > 0 & NH3 > 0, NH3 / CO2, NA_real_),
  NH3_CH4 = fifelse(CH4 > 0 & NH3 > 0, NH3 / CH4, NA_real_)
)]

variable_names <- c(
  "CO2", "CH4", "NH3", "CH4_CO2", "NH3_CO2", "NH3_CH4"
)
variable_labels <- c(
  CO2 = "CO2",
  CH4 = "CH4",
  NH3 = "NH3",
  CH4_CO2 = "CH4 / CO2",
  NH3_CO2 = "NH3 / CO2",
  NH3_CH4 = "NH3 / CH4"
)

long_data <- melt(
  campaign_data,
  id.vars = c("campaign", "analyser", "location"),
  measure.vars = variable_names,
  variable.name = "variable",
  value.name = "value"
)

cv_summary <- long_data[
  is.finite(value) & value > 0,
  .(
    n = .N,
    mean = mean(value),
    sd = sd(value),
    cv_pct = 100 * sd(value) / mean(value),
    analysers = paste(sort(unique(analyser)), collapse = ", ")
  ),
  by = .(campaign, location, variable)
]
cv_summary[, cv_class := cut(
  pmin(cv_pct, 100),
  breaks = c(-Inf, 20, 40, 60, 80, Inf),
  labels = c("0-20%", ">20-40%", ">40-60%", ">60-80%", ">80-100%"),
  right = TRUE
)]
cv_summary[, location_number := suppressWarnings(as.integer(location))]
cv_summary[, variable_label := unname(variable_labels[variable])]

internal_summary <- merge(
  cv_summary[!is.na(location_number)],
  coordinates,
  by = "location_number",
  all.x = TRUE,
  sort = FALSE
)
reference_summary <- cv_summary[is.na(location_number)]

setorder(internal_summary, campaign, variable, location_number)
setorder(reference_summary, campaign, location, variable)

fwrite(
  internal_summary,
  file.path(output_data_dir, "campaign_location_variable_cv_summary.csv")
)
fwrite(
  reference_summary,
  file.path(output_data_dir, "campaign_reference_location_cv_summary.csv")
)


###############################################################################
##### STATIC FLOOR-PLAN CV MAPS
###############################################################################

floorplan <- png::readPNG(floorplan_file)
image_width <- dim(floorplan)[2]
image_height <- dim(floorplan)[1]
# Remove the original height legend below the plan. The CV maps provide their
# own height-shape and CV-colour legends.
map_image_height <- 770L
floorplan_map <- floorplan[seq_len(map_image_height), , , drop = FALSE]
# De-emphasise the original coloured sampling symbols already embedded in the
# floor plan, so only the campaign-specific CV markers carry colour.
floorplan_luminance <- (
  0.299 * floorplan_map[, , 1] +
    0.587 * floorplan_map[, , 2] +
    0.114 * floorplan_map[, , 3]
)
floorplan_luminance <- 0.35 * floorplan_luminance + 0.65
for (channel in 1:3) {
  floorplan_map[, , channel] <- floorplan_luminance
}

cv_colours <- c(
  "0-20%" = "#2C7BB6",
  ">20-40%" = "#1A9641",
  ">40-60%" = "#F9D423",
  ">60-80%" = "#F28E2B",
  ">80-100%" = "#D7191C"
)
height_shapes <- c(top = 24, middle = 22, bottom = 21)

for (campaign_name in names(campaign_files)) {
  map_data <- copy(internal_summary[campaign == campaign_name])
  map_data[, variable := factor(
    variable,
    levels = variable_names,
    labels = unname(variable_labels[variable_names])
  )]

  map_plot <- ggplot() +
    annotation_raster(
      floorplan_map,
      xmin = 0,
      xmax = image_width,
      ymin = 0,
      ymax = map_image_height
    ) +
    geom_point(
      data = map_data,
      aes(
        x = map_x_pixel,
        y = map_image_height - map_y_pixel,
        fill = cv_class,
        shape = height
      ),
      size = 3.2,
      colour = "black",
      stroke = 0.45
    ) +
    facet_wrap(~ variable, ncol = 2) +
    scale_fill_manual(
      values = cv_colours,
      drop = FALSE,
      name = "CV"
    ) +
    scale_shape_manual(
      values = height_shapes,
      name = "Height"
    ) +
    coord_fixed(
      xlim = c(0, image_width),
      ylim = c(0, map_image_height),
      expand = FALSE,
      clip = "on"
    ) +
    labs(
      title = paste0(campaign_name, ": temporal coefficient of variation"),
      subtitle = "Colour represents CV; shape represents sampling height",
      x = NULL,
      y = NULL
    ) +
    theme_void(base_size = 10) +
    theme(
      strip.text = element_text(face = "bold"),
      strip.background = element_rect(
        fill = "grey95",
        colour = "grey40",
        linewidth = 0.3
      ),
      legend.position = "bottom",
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(colour = "grey30"),
      panel.spacing = grid::unit(0.4, "lines"),
      plot.margin = margin(5, 5, 15, 5)
    )

  file_stub <- gsub(" ", "_", tolower(campaign_name))
  ggsave(
    file.path(
      output_plot_dir,
      paste0(file_stub, "_floorplan_cv_gases_and_ratios.png")
    ),
    map_plot,
    width = 14,
    height = 11,
    dpi = 400,
    bg = "white"
  )
}


###############################################################################
##### INTERACTIVE ROTATABLE 3D POINT CLOUD
###############################################################################

# The z coordinate is deliberately categorical. The top sampling points were
# 60 cm below a sloping roof, so absolute top heights must not be fabricated.
plotly_colourscale <- list(
  c(0.00, "#2C7BB6"), c(0.20, "#2C7BB6"),
  c(0.20, "#1A9641"), c(0.40, "#1A9641"),
  c(0.40, "#F9D423"), c(0.60, "#F9D423"),
  c(0.60, "#F28E2B"), c(0.80, "#F28E2B"),
  c(0.80, "#D7191C"), c(1.00, "#D7191C")
)

make_3d_plot <- function(plot_data, plot_title) {
  plot_data <- copy(plot_data)
  plot_data[, selector := paste(campaign, variable_label, sep = " | ")]
  selector_levels <- unique(plot_data$selector)

  p <- plot_ly()
  for (selector_index in seq_along(selector_levels)) {
    selector_value <- selector_levels[selector_index]
    d <- plot_data[selector == selector_value]
    hover_text <- paste0(
      "<b>", d$campaign, "</b>",
      "<br>Variable: ", d$variable_label,
      "<br>Location: ", d$location,
      "<br>Height: ", d$height,
      "<br>Mean: ", signif(d$mean, 5),
      "<br>SD: ", signif(d$sd, 5),
      "<br>CV: ", sprintf("%.1f%%", d$cv_pct),
      "<br>n: ", d$n,
      "<br>Analyser(s): ", d$analysers
    )
    p <- add_trace(
      p,
      data = d,
      x = ~x_m,
      y = ~y_m,
      z = ~z_layer,
      type = "scatter3d",
      mode = "markers",
      text = hover_text,
      hoverinfo = "text",
      visible = selector_index == 1L,
      marker = list(
        size = 6,
        color = pmin(d$cv_pct, 100),
        colorscale = plotly_colourscale,
        cmin = 0,
        cmax = 100,
        showscale = selector_index == 1L,
        colorbar = list(title = "CV (%)"),
        symbol = unname(c(
          top = "diamond",
          middle = "square",
          bottom = "circle"
        )[d$height]),
        line = list(color = "#303030", width = 1)
      ),
      name = selector_value,
      showlegend = FALSE
    )
  }

  # A combined selector avoids conflicting visibility states from two separate
  # Plotly menus while still giving access to every campaign-variable pair.
  menu_buttons <- lapply(seq_along(selector_levels), function(i) {
    list(
      method = "update",
      args = list(
        list(
          visible = seq_along(selector_levels) == i,
          "marker.showscale" = seq_along(selector_levels) == i
        ),
        list(title = paste(plot_title, selector_levels[i], sep = "<br>"))
      ),
      label = selector_levels[i]
    )
  })

  layout(
    p,
    title = paste(plot_title, selector_levels[1], sep = "<br>"),
    scene = list(
      xaxis = list(title = "Barn length (m)", range = c(-1, 40)),
      yaxis = list(title = "Barn width (m)", range = c(-1, 19)),
      zaxis = list(
        title = "Sampling layer",
        tickmode = "array",
        tickvals = c(1, 2, 3),
        ticktext = c("Bottom", "Middle", "Top"),
        range = c(0.7, 3.3)
      ),
      aspectmode = "manual",
      aspectratio = list(x = 2.2, y = 1, z = 0.55),
      camera = list(eye = list(x = 1.5, y = 1.5, z = 1.0))
    ),
    updatemenus = list(list(
      type = "dropdown",
      x = 0,
      y = 1.12,
      xanchor = "left",
      yanchor = "top",
      buttons = menu_buttons
    )),
    margin = list(l = 20, r = 20, b = 20, t = 110)
  )
}

all_3d <- make_3d_plot(
  internal_summary,
  "Temporal CV in the barn"
)
htmlwidgets::saveWidget(
  all_3d,
  file.path(output_html_dir, "interactive_CV_3D_all_campaigns.html"),
  selfcontained = TRUE,
  title = "Interactive 3D campaign CV"
)

for (campaign_name in names(campaign_files)) {
  campaign_3d <- make_3d_plot(
    internal_summary[campaign == campaign_name],
    paste0(campaign_name, " temporal CV")
  )
  file_stub <- gsub(" ", "_", tolower(campaign_name))
  htmlwidgets::saveWidget(
    campaign_3d,
    file.path(
      output_html_dir,
      paste0(file_stub, "_interactive_CV_3D.html")
    ),
    selfcontained = TRUE,
    title = paste0(campaign_name, " interactive 3D CV")
  )
}

# The HTML files are self-contained. htmlwidgets may leave temporary dependency
# directories after Pandoc has embedded them; remove only those generated
# sibling directories after confirming that every HTML file exists.
expected_html <- c(
  file.path(output_html_dir, "interactive_CV_3D_all_campaigns.html"),
  file.path(
    output_html_dir,
    paste0(
      gsub(" ", "_", tolower(names(campaign_files))),
      "_interactive_CV_3D.html"
    )
  )
)
if (all(file.exists(expected_html))) {
  generated_dependency_dirs <- list.dirs(
    output_html_dir,
    full.names = TRUE,
    recursive = FALSE
  )
  generated_dependency_dirs <- generated_dependency_dirs[
    grepl("_files$", generated_dependency_dirs)
  ]
  if (length(generated_dependency_dirs) > 0L) {
    unlink(generated_dependency_dirs, recursive = TRUE, force = TRUE)
  }
}


###############################################################################
##### ANALYSIS RECORD
###############################################################################

record <- c(
  "Campaign CV floor-plan and 3D spatial supplement",
  paste0("Run time: ", Sys.time()),
  paste0("Input directory: ", campaign_dir),
  paste0("Floor plan: ", floorplan_file),
  "CV definition: 100 * SD / arithmetic mean.",
  "Ratios were calculated at each measurement timestamp before summarising.",
  "CV display scale: 0-100%, in equal 20-percentage-point classes.",
  "Values above 100% are retained numerically and capped at 100% for colour.",
  "The same CV scale is used for all campaigns, gases, and ratios.",
  "Only positive finite concentrations and ratios were included.",
  paste0(
    "Internal numeric sampling locations mapped: ",
    uniqueN(internal_summary$location_number)
  ),
  paste0(
    "Reference labels retained in a separate summary: ",
    paste(sort(unique(reference_summary$location)), collapse = ", ")
  ),
  paste0(
    "Reference labels are not plotted spatially because their meanings differ ",
    "among campaigns and they are not ordinary internal point locations."
  ),
  paste0(
    "3D z values are categorical layers, not fabricated absolute top heights."
  ),
  paste0("Output data directory: ", output_data_dir),
  paste0("Output plot directory: ", output_plot_dir),
  "",
  capture.output(sessionInfo())
)
writeLines(record, file.path(output_data_dir, "campaign_cv_spatial_readme.txt"))

message("Campaign CV maps and interactive 3D outputs completed.")
