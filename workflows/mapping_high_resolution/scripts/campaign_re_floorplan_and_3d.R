###############################################################################
# Campaign relative-error floor-plan maps and interactive 3D supplement
#
# Campaigns 1 and 2:
#   baseline = equally weighted mean of all internal numeric location means.
#
# Campaigns 3 and 4:
#   baseline = mean at location "s", which represents the "in" ring line in
#   these campaigns.
#
# Signed RE is retained in the output table. Absolute RE is mapped because the
# requested common display scale is 0-100%.
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
coordinate_file <- file.path(
  workflow_dir,
  "clean_data",
  "campaign_cv_spatial",
  "location_coordinate_lookup.csv"
)
floorplan_file <- file.path(
  workflow_dir,
  "Manuscript_3_Mapping_high_resolution_concentration_Latex_draft",
  "figures",
  "Fig1a_campaign1_sampling_layout_no_ringline.png"
)
output_data_dir <- file.path(
  workflow_dir,
  "clean_data",
  "campaign_re_spatial"
)
output_plot_dir <- file.path(
  workflow_dir,
  "plots",
  "campaign_re_spatial"
)
output_html_dir <- file.path(output_plot_dir, "interactive_3d")

dir.create(output_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_html_dir, recursive = TRUE, showWarnings = FALSE)

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

if (!file.exists(coordinate_file)) {
  stop(
    "Coordinate lookup not found. Run campaign_cv_floorplan_and_3d.R first."
  )
}
if (!file.exists(floorplan_file)) {
  stop("Floor-plan image not found: ", floorplan_file)
}


###############################################################################
##### READ DATA AND CALCULATE TIME-RESOLVED RATIOS
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
# Harmonise Campaign 1 FTIR values to the CRDS scale using the 24-hour
# co-located field comparison before calculating cross-instrument means and RE.
campaign_data[
  campaign == "Campaign 1" & grepl("^FTIR", analyser),
  `:=`(
    CO2 = CO2 / 1.06,
    CH4 = CH4 / 1.06,
    NH3 = NH3 / 1.09
  )
]
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

location_means <- long_data[
  is.finite(value) & value > 0,
  .(
    n = .N,
    location_mean = mean(value),
    location_sd = sd(value),
    analysers = paste(sort(unique(analyser)), collapse = ", ")
  ),
  by = .(campaign, location, variable)
]
location_means[
  ,
  location_number := suppressWarnings(as.integer(location))
]


###############################################################################
##### DEFINE CAMPAIGN-SPECIFIC BASELINES
###############################################################################

# Equal location weighting prevents locations with more observations from
# dominating the Campaign 1 and Campaign 2 barn-wide baseline.
internal_baselines <- location_means[
  campaign %in% c("Campaign 1", "Campaign 2") &
    !is.na(location_number),
  .(
    baseline_mean = mean(location_mean),
    baseline_locations = .N,
    baseline_location = "Mean of all internal location means",
    baseline_type = "internal_location_mean"
  ),
  by = .(campaign, variable)
]

# In Campaigns 3 and 4, the cleaned label "s" is the source "in" ring line.
ringline_baselines <- location_means[
  campaign %in% c("Campaign 3", "Campaign 4") &
    location == "s",
  .(
    baseline_mean = location_mean,
    baseline_locations = 1L,
    baseline_location = "s (source label 'in'; ring line)",
    baseline_type = "ringline_in"
  ),
  by = .(campaign, variable)
]

baseline_table <- rbindlist(
  list(internal_baselines, ringline_baselines),
  use.names = TRUE
)

expected_baselines <- CJ(
  campaign = names(campaign_files),
  variable = variable_names,
  unique = TRUE
)
missing_baselines <- baseline_table[
  expected_baselines,
  on = .(campaign, variable)
][is.na(baseline_mean)]
if (nrow(missing_baselines) > 0L) {
  stop(
    "Missing baseline(s): ",
    paste(
      paste(missing_baselines$campaign, missing_baselines$variable),
      collapse = ", "
    )
  )
}

internal_means <- location_means[!is.na(location_number)]
re_summary <- merge(
  internal_means,
  baseline_table,
  by = c("campaign", "variable"),
  all.x = TRUE,
  sort = FALSE
)
re_summary[, `:=`(
  signed_re_pct = 100 * (location_mean - baseline_mean) / baseline_mean,
  absolute_re_pct = abs(
    100 * (location_mean - baseline_mean) / baseline_mean
  )
)]
re_summary[, display_re_pct := pmin(absolute_re_pct, 100)]
re_summary[, re_class := cut(
  display_re_pct,
  breaks = c(-Inf, 20, 40, 60, 80, Inf),
  labels = c("0-20%", ">20-40%", ">40-60%", ">60-80%", ">80-100%"),
  right = TRUE
)]
re_summary[, variable_label := unname(variable_labels[variable])]

coordinates <- fread(coordinate_file)
re_summary <- merge(
  re_summary,
  coordinates,
  by = "location_number",
  all.x = TRUE,
  sort = FALSE
)
setorder(re_summary, campaign, variable, location_number)
setorder(baseline_table, campaign, variable)

fwrite(
  re_summary,
  file.path(output_data_dir, "campaign_location_relative_error_summary.csv")
)
fwrite(
  baseline_table,
  file.path(output_data_dir, "campaign_relative_error_baselines.csv")
)


###############################################################################
##### STATIC FLOOR-PLAN ABSOLUTE-RE MAPS
###############################################################################

floorplan <- png::readPNG(floorplan_file)
image_width <- dim(floorplan)[2]
map_image_height <- 770L
floorplan_map <- floorplan[seq_len(map_image_height), , , drop = FALSE]
floorplan_luminance <- (
  0.299 * floorplan_map[, , 1] +
    0.587 * floorplan_map[, , 2] +
    0.114 * floorplan_map[, , 3]
)
floorplan_luminance <- 0.35 * floorplan_luminance + 0.65
for (channel in 1:3) {
  floorplan_map[, , channel] <- floorplan_luminance
}

re_colours <- c(
  "0-20%" = "#2C7BB6",
  ">20-40%" = "#1A9641",
  ">40-60%" = "#F9D423",
  ">60-80%" = "#F28E2B",
  ">80-100%" = "#D7191C"
)
height_shapes <- c(top = 24, middle = 22, bottom = 21)

for (campaign_name in names(campaign_files)) {
  map_data <- copy(re_summary[campaign == campaign_name])
  map_data[, variable := factor(
    variable,
    levels = variable_names,
    labels = unname(variable_labels[variable_names])
  )]
  baseline_note <- if (campaign_name %in% c("Campaign 1", "Campaign 2")) {
    "Baseline: equally weighted mean of internal location means"
  } else {
    "Baseline: 'in' ring line (cleaned location label s)"
  }

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
        fill = re_class,
        shape = height
      ),
      size = 3.2,
      colour = "black",
      stroke = 0.45
    ) +
    facet_wrap(~ variable, ncol = 2) +
    scale_fill_manual(
      values = re_colours,
      drop = FALSE,
      name = "|RE|"
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
      title = paste0(
        campaign_name,
        ": absolute relative percentage error"
      ),
      subtitle = paste0(
        baseline_note,
        "; colour scale capped at 100%"
      ),
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
      paste0(file_stub, "_floorplan_RE_gases_and_ratios.png")
    ),
    map_plot,
    width = 14,
    height = 11,
    dpi = 400,
    bg = "white"
  )
}


###############################################################################
##### INTERACTIVE ROTATABLE 3D ABSOLUTE-RE POINT CLOUD
###############################################################################

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
      "<br>Location mean: ", signif(d$location_mean, 5),
      "<br>Baseline mean: ", signif(d$baseline_mean, 5),
      "<br>Signed RE: ", sprintf("%+.1f%%", d$signed_re_pct),
      "<br>Absolute RE: ", sprintf("%.1f%%", d$absolute_re_pct),
      "<br>Baseline: ", d$baseline_location,
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
        color = d$display_re_pct,
        colorscale = plotly_colourscale,
        cmin = 0,
        cmax = 100,
        showscale = selector_index == 1L,
        colorbar = list(title = "|RE| (%)"),
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
  re_summary,
  "Absolute relative error in the barn"
)
htmlwidgets::saveWidget(
  all_3d,
  file.path(output_html_dir, "interactive_RE_3D_all_campaigns.html"),
  selfcontained = TRUE,
  title = "Interactive 3D campaign relative error"
)

for (campaign_name in names(campaign_files)) {
  campaign_3d <- make_3d_plot(
    re_summary[campaign == campaign_name],
    paste0(campaign_name, " absolute relative error")
  )
  file_stub <- gsub(" ", "_", tolower(campaign_name))
  htmlwidgets::saveWidget(
    campaign_3d,
    file.path(
      output_html_dir,
      paste0(file_stub, "_interactive_RE_3D.html")
    ),
    selfcontained = TRUE,
    title = paste0(campaign_name, " interactive 3D RE")
  )
}

expected_html <- c(
  file.path(output_html_dir, "interactive_RE_3D_all_campaigns.html"),
  file.path(
    output_html_dir,
    paste0(
      gsub(" ", "_", tolower(names(campaign_files))),
      "_interactive_RE_3D.html"
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
  "Campaign relative-error floor-plan and 3D spatial supplement",
  paste0("Run time: ", Sys.time()),
  paste0("Input directory: ", campaign_dir),
  "Signed RE definition: 100 * (location mean - baseline) / baseline.",
  "Absolute RE is mapped; signed and absolute values are retained in the CSV.",
  paste0(
    "Campaigns 1-2 baseline: equally weighted mean of internal location means."
  ),
  "Campaigns 3-4 baseline: location s, originating from source label 'in'.",
  "Display scale: 0-100%, in equal 20-percentage-point classes.",
  "Values above 100% are retained numerically and capped at 100% for colour.",
  "Ratios were calculated at each measurement timestamp before averaging.",
  "Only positive finite concentrations and ratios were included.",
  "3D z values are categorical layers, not fabricated absolute top heights.",
  paste0("Output data directory: ", output_data_dir),
  paste0("Output plot directory: ", output_plot_dir),
  "",
  capture.output(sessionInfo())
)
writeLines(record, file.path(output_data_dir, "campaign_re_spatial_readme.txt"))

message("Campaign RE maps and interactive 3D outputs completed.")
