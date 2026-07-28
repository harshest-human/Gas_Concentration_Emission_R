###############################################################################
##### MANUSCRIPT 1: CAMPAIGNS 1 AND 2
##### HIGH-RESOLUTION GAS MAPPING AND SAMPLING-DENSITY REDUCTION
###############################################################################
#
# Outputs are written to:
#   clean_data/manuscript1/
#   plots/manuscript1/
#
# Campaign 1 uses the FTIR concentrations harmonised to the CRDS scale.
# Location "s" (south background; conceptual location 52) is retained in the
# clean analytical dataset but excluded from internal spatial comparisons.
#
# Vertical-group colours required for all manuscript figures:
#   top    = orange
#   middle = green3
#   bottom = steelblue1
###############################################################################

required_packages <- c("data.table", "ggplot2", "patchwork")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Install required package(s): ", paste(missing_packages, collapse = ", "))
}

library(data.table)
library(ggplot2)
library(patchwork)

workflow_dir <- paste0(
  "D:/Data_Analysis_R/Gas_Concentration_Emission_R/",
  "workflows/mapping_high_resolution"
)
clean_output_dir <- file.path(workflow_dir, "clean_data", "manuscript1")
plot_output_dir <- file.path(workflow_dir, "plots", "manuscript1")
dir.create(clean_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_output_dir, recursive = TRUE, showWarnings = FALSE)

campaign1_path <- file.path(
  workflow_dir,
  "clean_data",
  "1_campaign",
  paste0(
    "20241001.000232_20241024.143858_",
    "campaign1_raw_and_CRDS_scale_corrected.csv"
  )
)
campaign2_path <- file.path(
  workflow_dir,
  "clean_data",
  "2_campaign",
  "20241116.000647_20241231.235710_CRDS8_CRDS9_campaign2_combined.csv"
)

top_locations <- seq(1L, 49L, by = 3L)
mid_locations <- seq(2L, 50L, by = 3L)
bottom_locations <- seq(3L, 51L, by = 3L)
internal_locations <- 1:51
reduced_locations <- c(top_locations, bottom_locations)

vertical_colours <- c(
  "top" = "orange",
  "mid" = "green3",
  "bottom" = "steelblue1"
)

parse_time <- function(x) {
  as.POSIXct(x, format = "%Y-%m-%d %H:%M:%S", tz = "Europe/Berlin")
}

assign_vertical_group <- function(location) {
  location_number <- suppressWarnings(as.integer(as.character(location)))
  fcase(
    location_number %in% top_locations, "top",
    location_number %in% mid_locations, "mid",
    location_number %in% bottom_locations, "bottom",
    as.character(location) == "s", "out",
    default = NA_character_
  )
}

assign_horizontal_position <- function(location) {
  location_number <- suppressWarnings(as.integer(as.character(location)))
  fifelse(
    is.na(location_number),
    NA_integer_,
    as.integer(ceiling(location_number / 3))
  )
}

floor_to_seconds <- function(x, seconds) {
  as.POSIXct(
    floor(as.numeric(x) / seconds) * seconds,
    origin = "1970-01-01",
    tz = "Europe/Berlin"
  )
}

mean_ci <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  mean_value <- if (n > 0L) mean(x) else NA_real_
  se <- if (n > 1L) sd(x) / sqrt(n) else NA_real_
  multiplier <- if (n > 1L) qt(0.975, df = n - 1L) else NA_real_
  list(
    mean = mean_value,
    lower = mean_value - multiplier * se,
    upper = mean_value + multiplier * se,
    n_blocks = n
  )
}

lin_ccc <- function(observed, predicted) {
  keep <- is.finite(observed) & is.finite(predicted)
  observed <- observed[keep]
  predicted <- predicted[keep]
  covariance <- cov(observed, predicted)
  2 * covariance / (
    var(observed) + var(predicted) +
      (mean(observed) - mean(predicted))^2
  )
}

gas_long <- function(data, columns, value_name = "concentration") {
  melt(
    data,
    id.vars = setdiff(names(data), columns),
    measure.vars = columns,
    variable.name = "gas",
    value.name = value_name
  )
}


###############################################################################
##### PREPARE CAMPAIGN 1 AND CAMPAIGN 2 ANALYTICAL DATA
###############################################################################

campaign1 <- fread(
  campaign1_path,
  colClasses = list(
    character = c("DATE.TIME", "analyser", "location")
  ),
  showProgress = FALSE
)
campaign1[, DATE.TIME := parse_time(DATE.TIME)]
campaign1[, campaign := "Campaign 1"]
campaign1[, vgroup := assign_vertical_group(location)]
campaign1[, horizontal_position := assign_horizontal_position(location)]

# Non-positive concentration entries cannot represent valid mole fractions and
# would make the gas ratios undefined. They are retained in the source CSV but
# represented as missing values in the manuscript analytical dataset.
# CO2 below 300 ppm is physically implausible for barn/ambient air and denotes
# an invalid transition or analyser record rather than a valid concentration.
campaign1[
  CO2_raw < 300 | !is.finite(CO2_raw),
  c("CO2_raw", "CO2_corr") := list(NA_real_, NA_real_)
]
campaign1[
  CH4_raw <= 0 | !is.finite(CH4_raw),
  c("CH4_raw", "CH4_corr") := list(NA_real_, NA_real_)
]
campaign1[
  NH3_raw <= 0 | !is.finite(NH3_raw),
  c("NH3_raw", "NH3_corr") := list(NA_real_, NA_real_)
]
campaign1[
  H2O_raw <= 0 | !is.finite(H2O_raw),
  c("H2O_raw", "H2O_corr") := list(NA_real_, NA_real_)
]
campaign1[, `:=`(
  CH4_CO2_pct = fifelse(
    CO2_corr > 0 & CH4_corr > 0,
    100 * CH4_corr / CO2_corr,
    NA_real_
  ),
  NH3_CO2_pct = fifelse(
    CO2_corr > 0 & NH3_corr > 0,
    100 * NH3_corr / CO2_corr,
    NA_real_
  )
)]

campaign2 <- fread(
  campaign2_path,
  colClasses = list(
    character = c("DATE.TIME", "analyser", "location")
  ),
  showProgress = FALSE
)
campaign2[, DATE.TIME := parse_time(DATE.TIME)]
campaign2[, campaign := "Campaign 2"]
campaign2[, vgroup := assign_vertical_group(location)]
campaign2[, horizontal_position := assign_horizontal_position(location)]
campaign2[CO2 < 300 | !is.finite(CO2), CO2 := NA_real_]
campaign2[CH4 <= 0 | !is.finite(CH4), CH4 := NA_real_]
campaign2[NH3 <= 0 | !is.finite(NH3), NH3 := NA_real_]
campaign2[H2O <= 0 | !is.finite(H2O), H2O := NA_real_]
campaign2[, `:=`(
  CO2_raw = CO2,
  CH4_raw = CH4,
  NH3_raw = NH3,
  H2O_raw = H2O,
  CO2_corr = CO2,
  CH4_corr = CH4,
  NH3_corr = NH3,
  H2O_corr = H2O,
  correction_applied = 0L,
  correction_basis = "CRDS retained without correction",
  CH4_CO2_pct = fifelse(
    CO2 > 0 & CH4 > 0,
    100 * CH4 / CO2,
    NA_real_
  ),
  NH3_CO2_pct = fifelse(
    CO2 > 0 & NH3 > 0,
    100 * NH3 / CO2,
    NA_real_
  )
)]
campaign2[, c("CO2", "CH4", "NH3", "H2O") := NULL]

analytical_columns <- c(
  "DATE.TIME", "campaign", "analyser", "location",
  "horizontal_position", "vgroup",
  "CO2_raw", "CH4_raw", "NH3_raw", "H2O_raw",
  "CO2_corr", "CH4_corr", "NH3_corr", "H2O_corr",
  "CH4_CO2_pct", "NH3_CO2_pct",
  "correction_applied", "correction_basis"
)
setcolorder(campaign1, analytical_columns)
setcolorder(campaign2, analytical_columns)

manuscript_data <- rbindlist(
  list(campaign1, campaign2),
  use.names = TRUE,
  fill = TRUE
)
setorder(manuscript_data, DATE.TIME, analyser, location)

fwrite(
  manuscript_data,
  file.path(clean_output_dir, "manuscript1_campaign1_2_analytical_data.csv"),
  quote = TRUE,
  dateTimeAs = "write.csv"
)


###############################################################################
##### FIGURE 2: CAMPAIGN 1 CONCENTRATIONS AT 51 INTERNAL LOCATIONS
###############################################################################

# Daily location means are the independent blocks used for descriptive
# confidence intervals. The south-background location is excluded.
c1_internal <- campaign1[location != "s"]
c1_internal[, date := as.IDate(DATE.TIME, tz = "Europe/Berlin")]

c1_daily <- c1_internal[
  ,
  .(
    CO2 = mean(CO2_corr, na.rm = TRUE),
    CH4 = mean(CH4_corr, na.rm = TRUE),
    NH3 = mean(NH3_corr, na.rm = TRUE)
  ),
  by = .(date, location, horizontal_position, vgroup)
]

c1_daily_long <- melt(
  c1_daily,
  id.vars = c("date", "location", "horizontal_position", "vgroup"),
  measure.vars = c("CO2", "CH4", "NH3"),
  variable.name = "gas",
  value.name = "concentration"
)

c1_location_summary <- c1_daily_long[
  ,
  mean_ci(concentration),
  by = .(gas, location, horizontal_position, vgroup)
]
c1_location_summary[, location_number := as.integer(location)]
c1_location_summary[, gas := factor(gas, levels = c("CO2", "CH4", "NH3"))]
c1_location_summary[, vgroup := factor(
  vgroup,
  levels = c("top", "mid", "bottom")
)]

gas_labels <- c(
  "CO2" = expression(CO[2]~"(ppm)"),
  "CH4" = expression(CH[4]~"(ppm)"),
  "NH3" = expression(NH[3]~"(ppm)")
)

figure2 <- ggplot(
  c1_location_summary,
  aes(x = location_number, y = mean, colour = vgroup)
) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0,
    linewidth = 0.35,
    alpha = 0.75
  ) +
  geom_point(size = 1.8) +
  facet_wrap(
    ~gas,
    ncol = 1,
    scales = "free_y",
    labeller = as_labeller(gas_labels)
  ) +
  scale_colour_manual(values = vertical_colours, guide = "none") +
  scale_x_continuous(
    breaks = seq(1, 51, by = 2),
    minor_breaks = NULL
  ) +
  labs(x = "Internal sampling location", y = NULL) +
  theme_classic(base_size = 9) +
  theme(
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold"),
    panel.spacing = grid::unit(0.8, "lines"),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

ggsave(
  file.path(
    plot_output_dir,
    "Fig2_campaign1_corrected_concentrations_51_locations.png"
  ),
  figure2,
  width = 180,
  height = 190,
  units = "mm",
  dpi = 600,
  bg = "white"
)


###############################################################################
##### FIGURE 3: SAMPLING-CONFIGURATION DEVIATION FROM FULL DESIGN
###############################################################################

# Two-hour blocks provide at least one observation from almost every sampling
# location despite the sequential analyser cycles. Only complete 51-location
# blocks are used as the full-density reference.
c1_internal[, block_2h := floor_to_seconds(DATE.TIME, 2L * 60L * 60L)]
c1_block_location <- c1_internal[
  ,
  .(
    CO2 = mean(CO2_corr, na.rm = TRUE),
    CH4 = mean(CH4_corr, na.rm = TRUE),
    NH3 = mean(NH3_corr, na.rm = TRUE)
  ),
  by = .(block_2h, location, vgroup)
]
c1_block_location[, location_number := as.integer(location)]

complete_blocks <- c1_block_location[
  ,
  .(n_locations = uniqueN(location_number)),
  by = block_2h
][n_locations == 51L, block_2h]

c1_complete <- c1_block_location[block_2h %in% complete_blocks]

configuration_map <- list(
  "Top only (17)" = top_locations,
  "Middle only (17)" = mid_locations,
  "Bottom only (17)" = bottom_locations,
  "Top + middle (34)" = c(top_locations, mid_locations),
  "Middle + bottom (34)" = c(mid_locations, bottom_locations),
  "Top + bottom (34)" = reduced_locations,
  "Full design (51)" = internal_locations
)

c1_complete_long <- melt(
  c1_complete,
  id.vars = c("block_2h", "location", "location_number", "vgroup"),
  measure.vars = c("CO2", "CH4", "NH3"),
  variable.name = "gas",
  value.name = "concentration"
)

configuration_results <- rbindlist(
  lapply(names(configuration_map), function(configuration_name) {
    retained <- configuration_map[[configuration_name]]
    c1_complete_long[
      location_number %in% retained,
      .(
        estimate = mean(concentration, na.rm = TRUE),
        heterogeneity = sd(log(concentration), na.rm = TRUE)
      ),
      by = .(block_2h, gas)
    ][, configuration := configuration_name]
  })
)

reference_results <- configuration_results[
  configuration == "Full design (51)",
  .(
    block_2h,
    gas,
    reference_mean = estimate,
    reference_heterogeneity = heterogeneity
  )
]
configuration_results <- merge(
  configuration_results,
  reference_results,
  by = c("block_2h", "gas"),
  all.x = TRUE
)
configuration_results[, relative_deviation_pct := 100 * (
  estimate - reference_mean
) / reference_mean]
configuration_results[, heterogeneity_deviation_pct := 100 * (
  heterogeneity - reference_heterogeneity
) / reference_heterogeneity]

configuration_levels <- names(configuration_map)
configuration_results[, configuration := factor(
  configuration,
  levels = configuration_levels
)]
configuration_results[, gas := factor(gas, levels = c("CO2", "CH4", "NH3"))]

figure3 <- ggplot(
  configuration_results[configuration != "Full design (51)"],
  aes(x = configuration, y = relative_deviation_pct)
) +
  annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = -10,
    ymax = 10,
    fill = "#F6E8B1",
    alpha = 0.35
  ) +
  annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = -5,
    ymax = 5,
    fill = "#D9EAD3",
    alpha = 0.55
  ) +
  geom_hline(yintercept = 0, linewidth = 0.45) +
  geom_boxplot(
    width = 0.62,
    outlier.shape = NA,
    linewidth = 0.4,
    fill = "white"
  ) +
  coord_cartesian(ylim = quantile(
    configuration_results[
      configuration != "Full design (51)"
    ][["relative_deviation_pct"]],
    probs = c(0.005, 0.995),
    na.rm = TRUE
  )) +
  facet_wrap(
    ~gas,
    ncol = 1,
    scales = "free_y",
    labeller = as_labeller(gas_labels)
  ) +
  scale_x_discrete(
    labels = c(
      "Top only (17)" = "Top",
      "Middle only (17)" = "Middle",
      "Bottom only (17)" = "Bottom",
      "Top + middle (34)" = "Top + middle",
      "Middle + bottom (34)" = "Middle + bottom",
      "Top + bottom (34)" = "Top + bottom"
    )
  ) +
  labs(
    x = "Sampling configuration",
    y = "Deviation from full-density mean (%)"
  ) +
  theme_classic(base_size = 9) +
  theme(
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold"),
    panel.spacing = grid::unit(0.8, "lines"),
    axis.text.x = element_text(angle = 15, hjust = 1)
  )

ggsave(
  file.path(
    plot_output_dir,
    "Fig3_campaign1_sampling_configuration_deviation.png"
  ),
  figure3,
  width = 180,
  height = 190,
  units = "mm",
  dpi = 600,
  bg = "white"
)


###############################################################################
##### FIGURE 4: FULL 51-LOCATION VERSUS REDUCED 34-LOCATION DESIGN
###############################################################################

agreement <- configuration_results[
  configuration == "Top + bottom (34)",
  .(
    block_2h,
    gas,
    reduced_mean = estimate,
    full_mean = reference_mean,
    reduced_heterogeneity = heterogeneity,
    full_heterogeneity = reference_heterogeneity
  )
]

agreement_metrics <- agreement[
  ,
  .(
    n_blocks = .N,
    mean_bias = mean(reduced_mean - full_mean, na.rm = TRUE),
    relative_bias_pct = 100 * mean(
      (reduced_mean - full_mean) / full_mean,
      na.rm = TRUE
    ),
    rmse = sqrt(mean((reduced_mean - full_mean)^2, na.rm = TRUE)),
    R2 = cor(reduced_mean, full_mean, use = "complete.obs")^2,
    CCC = lin_ccc(full_mean, reduced_mean),
    heterogeneity_relative_bias_pct = 100 * mean(
      (reduced_heterogeneity - full_heterogeneity) /
        full_heterogeneity,
      na.rm = TRUE
    ),
    heterogeneity_CCC = lin_ccc(
      full_heterogeneity,
      reduced_heterogeneity
    )
  ),
  by = gas
]

agreement_plot_data <- merge(
  agreement,
  agreement_metrics[, .(gas, CCC, relative_bias_pct)],
  by = "gas"
)
agreement_plot_data[, annotation := sprintf(
  "CCC = %.3f; bias = %.2f%%",
  CCC,
  relative_bias_pct
)]
agreement_plot_data[, gas := factor(
  gas,
  levels = c("CO2", "CH4", "NH3")
)]

figure4 <- ggplot(
  agreement_plot_data,
  aes(x = full_mean, y = reduced_mean)
) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.5) +
  geom_point(size = 1.3, alpha = 0.55, colour = "black") +
  facet_wrap(
    ~gas,
    scales = "free",
    labeller = as_labeller(gas_labels)
  ) +
  geom_text(
    data = unique(
      agreement_plot_data[, .(gas, annotation)]
    ),
    aes(x = -Inf, y = Inf, label = annotation),
    inherit.aes = FALSE,
    hjust = -0.05,
    vjust = 1.25,
    size = 2.7
  ) +
  labs(
    x = "Full 51-location mean (ppm)",
    y = "Reduced 34-location mean (ppm)"
  ) +
  theme_classic(base_size = 9) +
  theme(
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold"),
    panel.spacing = grid::unit(0.8, "lines")
  )

ggsave(
  file.path(
    plot_output_dir,
    "Fig4_campaign1_full_vs_reduced_agreement.png"
  ),
  figure4,
  width = 180,
  height = 75,
  units = "mm",
  dpi = 600,
  bg = "white"
)


###############################################################################
##### FIGURE 5: CAMPAIGN 2 LONGER-TERM REDUCED-DESIGN CONCENTRATIONS
###############################################################################

c2_internal <- campaign2[location != "s" & vgroup %in% c("top", "bottom")]
c2_internal[, date := as.IDate(DATE.TIME, tz = "Europe/Berlin")]

c2_daily <- c2_internal[
  ,
  .(
    CO2 = mean(CO2_corr, na.rm = TRUE),
    CH4 = mean(CH4_corr, na.rm = TRUE),
    NH3 = mean(NH3_corr, na.rm = TRUE)
  ),
  by = .(date, location, horizontal_position, vgroup)
]
c2_daily_long <- melt(
  c2_daily,
  id.vars = c("date", "location", "horizontal_position", "vgroup"),
  measure.vars = c("CO2", "CH4", "NH3"),
  variable.name = "gas",
  value.name = "concentration"
)
c2_location_summary <- c2_daily_long[
  ,
  mean_ci(concentration),
  by = .(gas, location, horizontal_position, vgroup)
]
c2_location_summary[, location_number := as.integer(location)]
c2_location_summary[, gas := factor(gas, levels = c("CO2", "CH4", "NH3"))]
c2_location_summary[, vgroup := factor(
  vgroup,
  levels = c("top", "bottom")
)]

figure5 <- ggplot(
  c2_location_summary,
  aes(x = location_number, y = mean, colour = vgroup)
) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0,
    linewidth = 0.35,
    alpha = 0.75
  ) +
  geom_point(size = 1.8) +
  facet_wrap(
    ~gas,
    ncol = 1,
    scales = "free_y",
    labeller = as_labeller(gas_labels)
  ) +
  scale_colour_manual(values = vertical_colours, guide = "none") +
  scale_x_continuous(
    breaks = sort(unique(c2_location_summary[["location_number"]])),
    minor_breaks = NULL
  ) +
  labs(x = "Internal sampling location", y = NULL) +
  theme_classic(base_size = 9) +
  theme(
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold"),
    panel.spacing = grid::unit(0.8, "lines"),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

ggsave(
  file.path(
    plot_output_dir,
    "Fig5_campaign2_concentrations_reduced_design.png"
  ),
  figure5,
  width = 180,
  height = 190,
  units = "mm",
  dpi = 600,
  bg = "white"
)


###############################################################################
##### MANUSCRIPT TABLES AND REPRODUCIBILITY OUTPUTS
###############################################################################

campaign_summary <- manuscript_data[
  ,
  .(
    start = min(DATE.TIME, na.rm = TRUE),
    end = max(DATE.TIME, na.rm = TRUE),
    rows = .N,
    analysers = paste(sort(unique(analyser)), collapse = ", "),
    internal_locations = uniqueN(location[location != "s"]),
    south_background_rows = sum(location == "s"),
    top_locations = uniqueN(location[vgroup == "top"]),
    mid_locations = uniqueN(location[vgroup == "mid"]),
    bottom_locations = uniqueN(location[vgroup == "bottom"])
  ),
  by = campaign
]

configuration_summary <- configuration_results[
  configuration != "Full design (51)",
  .(
    n_blocks = .N,
    mean_relative_deviation_pct = mean(
      relative_deviation_pct,
      na.rm = TRUE
    ),
    median_absolute_deviation_pct = median(
      abs(relative_deviation_pct),
      na.rm = TRUE
    ),
    p95_absolute_deviation_pct = quantile(
      abs(relative_deviation_pct),
      0.95,
      na.rm = TRUE
    )
  ),
  by = .(gas, configuration)
]

gas_summary <- manuscript_data[
  location != "s",
  .(
    mean_CO2 = mean(CO2_corr, na.rm = TRUE),
    sd_CO2 = sd(CO2_corr, na.rm = TRUE),
    mean_CH4 = mean(CH4_corr, na.rm = TRUE),
    sd_CH4 = sd(CH4_corr, na.rm = TRUE),
    mean_NH3 = mean(NH3_corr, na.rm = TRUE),
    sd_NH3 = sd(NH3_corr, na.rm = TRUE),
    mean_CH4_CO2_pct = mean(CH4_CO2_pct, na.rm = TRUE),
    mean_NH3_CO2_pct = mean(NH3_CO2_pct, na.rm = TRUE)
  ),
  by = .(campaign, vgroup)
]

fwrite(
  campaign_summary,
  file.path(clean_output_dir, "Table_campaign_summary.csv"),
  quote = TRUE,
  dateTimeAs = "write.csv"
)
fwrite(
  gas_summary,
  file.path(clean_output_dir, "Table_gas_summary_by_campaign_height.csv"),
  quote = TRUE
)
fwrite(
  configuration_summary,
  file.path(clean_output_dir, "Table_sampling_configuration_performance.csv"),
  quote = TRUE
)
fwrite(
  agreement_metrics,
  file.path(clean_output_dir, "Table_full_vs_reduced_agreement.csv"),
  quote = TRUE
)
fwrite(
  configuration_results,
  file.path(clean_output_dir, "campaign1_configuration_block_results.csv"),
  quote = TRUE,
  dateTimeAs = "write.csv"
)

writeLines(
  c(
    "MANUSCRIPT 1 ANALYSIS SUMMARY",
    "=============================",
    paste0("Campaign 1 complete two-hour blocks: ", length(complete_blocks)),
    paste0("Campaign 1 rows: ", nrow(campaign1)),
    paste0("Campaign 2 rows: ", nrow(campaign2)),
    paste0("South-background label: s (conceptual location 52)"),
    "Spatial comparisons exclude s.",
    "FTIR correction: CO2/1.06, CH4/1.06, NH3/1.09.",
    "CRDS and H2O values are unchanged.",
    "Error bars in location plots use daily location means as blocks.",
    "Sampling configurations are compared in complete two-hour blocks."
  ),
  file.path(clean_output_dir, "analysis_readme.txt")
)

message("Manuscript 1 analytical rows: ", nrow(manuscript_data))
message("Complete Campaign 1 two-hour blocks: ", length(complete_blocks))
message("Clean outputs: ", clean_output_dir)
message("Plot outputs: ", plot_output_dir)
