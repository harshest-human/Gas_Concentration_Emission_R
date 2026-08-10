###############################################################################
# Manuscript 3: plotting and statistical tables only
#
# Scientific emphasis for this version:
#   - median-based representativeness is primary;
#   - ratios remain dimensionless (not percentages);
#   - Campaigns 1 and 2 use common y scales within each response;
#   - both mean +/- SD and median +/- SD are shown;
#   - Campaign 2 entropy is calculated for all gases and ratios;
#   - CV and relative error are deliberately excluded.
#
# This script does not read, edit, or write the LaTeX manuscript.
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(lubridate)
})

workflow <- paste0(
  "D:/Data_Analysis_R/Gas_Concentration_Emission_R/",
  "workflows/mapping_high_resolution"
)
input_root <- file.path(workflow, "clean_data", "prepared_campaigns_v01")
output_root <- file.path(workflow, "clean_data", "manuscript3_analysis_v01")
table_dir <- file.path(output_root, "tables")
audit_dir <- file.path(output_root, "audit")
plot_dir <- file.path(workflow, "plots", "manuscript3_analysis_v01")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

responses <- c("CO2", "NH3_CO2", "CH4", "CH4_CO2", "NH3", "NH3_CH4")
response_labels <- c(
  CO2 = "CO2", NH3_CO2 = "NH3 / CO2",
  CH4 = "CH4", CH4_CO2 = "CH4 / CO2",
  NH3 = "NH3", NH3_CH4 = "NH3 / CH4"
)
height_colours <- c(
  top = "orange", middle = "green3", bottom = "steelblue1",
  reference = "grey45"
)


###############################################################################
##### READ AND HARMONISE CAMPAIGNS 1 AND 2
###############################################################################

read_campaign <- function(number) {
  path <- file.path(
    input_root, paste0("campaign_", number),
    paste0("campaign_", number, "_combined.csv")
  )
  x <- fread(path, colClasses = c(sampling.point = "character"))
  x[, DATE.TIME := ymd_hms(DATE.TIME, tz = "Europe/Berlin", quiet = TRUE)]
  x[, campaign := NULL]
  x[, campaign := paste("Campaign", number)]

  # Canonical reference spelling. Numeric sampling point 52 is retained as 52
  # because it is the explicitly assigned Campaign 2 sampling point.
  x[tolower(sampling.point) == "in", sampling.point := "ring_in"]
  x[tolower(sampling.point) == "s", sampling.point := "s"]
  x
}

x <- rbindlist(list(read_campaign(1), read_campaign(2)), use.names = TRUE)
for (g in c("CO2", "CH4", "NH3", "H2O")) {
  x[!is.finite(get(g)) | get(g) <= 0, (g) := NA_real_]
}
x[, `:=`(
  CH4_CO2 = fifelse(CO2 > 0, CH4 / CO2, NA_real_),
  NH3_CO2 = fifelse(CO2 > 0, NH3 / CO2, NA_real_),
  NH3_CH4 = fifelse(CH4 > 0, NH3 / CH4, NA_real_)
)]

x[, sampling.point.numeric := suppressWarnings(as.integer(sampling.point))]
x[, internal := !is.na(sampling.point.numeric) &
                  sampling.point.numeric >= 1L & sampling.point.numeric <= 51L]
x[, height := fifelse(
  !internal, "reference",
  c("top", "middle", "bottom")[(sampling.point.numeric - 1L) %% 3L + 1L]
)]
x[, height := factor(height, levels = names(height_colours))]
x[, campaign := factor(campaign, levels = c("Campaign 1", "Campaign 2"))]

long <- melt(
  x,
  id.vars = c("DATE.TIME", "campaign", "analyser", "sampling.point",
              "sampling.point.numeric", "internal", "height"),
  measure.vars = responses,
  variable.name = "response", value.name = "value"
)[is.finite(value) & value > 0]
long[, response := factor(response, levels = responses)]

fwrite(
  x,
  file.path(output_root, "manuscript3_campaign1_2_analysis_data.csv"),
  dateTimeAs = "write.csv"
)


###############################################################################
##### DESCRIPTIVE TABLES
###############################################################################

describe <- function(v) {
  v <- v[is.finite(v)]
  data.table(
    observations = length(v), mean = mean(v), median = median(v),
    sd = sd(v), mad = mad(v), minimum = min(v), q1 = quantile(v, 0.25),
    q3 = quantile(v, 0.75), maximum = max(v)
  )
}

location_descriptives <- long[, describe(value),
  by = .(campaign, response, sampling.point, sampling.point.numeric,
         internal, height)]
setorder(location_descriptives, campaign, response, sampling.point.numeric)
fwrite(location_descriptives,
       file.path(table_dir, "Table_01_sampling_point_descriptives.csv"))

campaign_descriptives <- long[internal == TRUE, describe(value),
  by = .(campaign, response)]
fwrite(campaign_descriptives,
       file.path(table_dir, "Table_02_campaign_descriptives.csv"))


###############################################################################
##### MEDIAN-BASED REPRESENTATIVENESS
###############################################################################

# Location medians are first calculated within hourly blocks. The equally
# weighted median across the available internal sampling points is the hourly
# network reference. This prevents locations with more raw rows from receiving
# greater weight.
hourly_location <- long[internal == TRUE,
  .(location.median = median(value), location.mean = mean(value), raw.rows = .N),
  by = .(
    campaign, response,
    hour = floor_date(DATE.TIME, "hour"),
    sampling.point, sampling.point.numeric, height
  )]

hourly_location[, `:=`(
  network.median = median(location.median),
  network.mean = mean(location.mean),
  locations.in.block = uniqueN(sampling.point)
), by = .(campaign, response, hour)]
hourly_location[, `:=`(
  median.difference = location.median - network.median,
  mean.difference = location.mean - network.mean
)]

representativeness <- hourly_location[, .(
  hourly.blocks = .N,
  median.bias = median(median.difference),
  median.absolute.difference = median(abs(median.difference)),
  median.mean.absolute.difference = mean(abs(median.difference)),
  median.root.mean.square.difference = sqrt(mean(median.difference^2)),
  q90.absolute.difference = quantile(abs(median.difference), 0.90),
  mean.bias = mean(mean.difference),
  mean.absolute.difference = median(abs(mean.difference)),
  mean.root.mean.square.difference = sqrt(mean(mean.difference^2)),
  spearman.rho = suppressWarnings(cor(
    location.median, network.median, method = "spearman",
    use = "pairwise.complete.obs"
  ))
), by = .(campaign, response, sampling.point, sampling.point.numeric, height)]
representativeness[, median.rank := frank(
  median.absolute.difference, ties.method = "min"
), by = .(campaign, response)]
representativeness[, mean.rank := frank(
  mean.absolute.difference, ties.method = "min"
), by = .(campaign, response)]
representativeness[, rank.change.mean.minus.median := mean.rank - median.rank]
setorder(representativeness, campaign, response, median.rank,
         sampling.point.numeric)
fwrite(representativeness,
       file.path(table_dir, "Table_03_median_representativeness_ranking.csv"))

sp25_blocks <- hourly_location[sampling.point.numeric == 25L]
sp25_summary <- sp25_blocks[, {
  difference <- location.median - network.median
  test <- tryCatch(
    wilcox.test(location.median, network.median, paired = TRUE, exact = FALSE),
    error = function(e) NULL
  )
  .(
    hourly.blocks = .N,
    network.median = median(network.median),
    SP25.median = median(location.median),
    median.difference = median(difference),
    median.absolute.difference = median(abs(difference)),
    mean.estimator.absolute.difference = median(abs(
      location.mean - network.mean
    )),
    q1.difference = quantile(difference, 0.25),
    q3.difference = quantile(difference, 0.75),
    spearman.rho = suppressWarnings(cor(
      location.median, network.median, method = "spearman",
      use = "pairwise.complete.obs"
    )),
    paired.wilcoxon.p = if (is.null(test)) NA_real_ else test$p.value
  )
}, by = .(campaign, response)]
fwrite(sp25_summary,
       file.path(table_dir, "Table_04_SP25_median_vs_network_median.csv"))


###############################################################################
##### SHANNON ENTROPY: CONCENTRATIONS AND RATIOS
###############################################################################

normalised_entropy <- function(v, breaks) {
  v <- v[is.finite(v)]
  if (length(v) < 8L || length(unique(breaks)) < 3L) return(NA_real_)
  groups <- cut(v, breaks = breaks, include.lowest = TRUE)
  p <- as.numeric(table(groups)) / length(v)
  p <- p[p > 0]
  max(0, -sum(p * log2(p)) / log2(length(breaks) - 1L))
}

# Four fixed-width bins are defined from the pooled range of each
# campaign-response combination, then applied identically to every SP.
entropy <- hourly_location[, {
  pooled_breaks <- seq(
    min(location.median), max(location.median), length.out = 5L
  )
  .SD[, .(
    hourly.blocks = .N,
    entropy.normalised = normalised_entropy(location.median, pooled_breaks)
  ), by = .(sampling.point, sampling.point.numeric, height)]
}, by = .(campaign, response)]
setorder(entropy, campaign, response, sampling.point.numeric)
fwrite(entropy, file.path(table_dir,
  "Table_05_Shannon_entropy_concentrations_and_ratios.csv"))
fwrite(entropy[campaign == "Campaign 2"], file.path(table_dir,
  "Table_05b_Campaign2_Shannon_entropy_concentrations_and_ratios.csv"))


###############################################################################
##### PLOTTING HELPERS
###############################################################################

facet_labels <- labeller(
  response = as_labeller(response_labels),
  campaign = label_value
)

common_theme <- theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = "black"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 0, size = 6)
  )

save_plot <- function(name, plot, width = 15, height = 15) {
  ggsave(file.path(plot_dir, paste0(name, ".png")), plot,
         width = width, height = height, dpi = 350, bg = "white")
  ggsave(file.path(plot_dir, paste0(name, ".pdf")), plot,
         width = width, height = height, bg = "white")
}

# One facet grid is used so each response row has the same y scale in both
# campaign columns. A common pooled 0.5--99.5% display range is applied within
# each response to prevent a few extreme values from obscuring the boxes.
plot_limits <- long[internal == TRUE, .(
  display.low = quantile(value, 0.005),
  display.high = quantile(value, 0.995)
), by = response]
fwrite(plot_limits, file.path(table_dir, "Table_06_plot_display_limits.csv"))

long_plot <- plot_limits[long[internal == TRUE], on = "response"]
long_plot <- long_plot[value >= display.low & value <= display.high]

p_box <- ggplot(long_plot,
  aes(factor(sampling.point.numeric), value, fill = height)) +
  geom_boxplot(outlier.shape = NA, width = 0.72, linewidth = 0.25) +
  facet_grid(response ~ campaign, scales = "free_y", labeller = facet_labels) +
  scale_fill_manual(values = height_colours, name = "Height") +
  labs(
    x = "Sampling point", y = "Concentration or dimensionless ratio",
    title = "Sampling-point distributions and medians",
    subtitle = "Common 0.5-99.5% display range per response across both campaigns"
  ) + common_theme
save_plot("Fig_01_boxplots_median_campaign1_2_common_ratio_scales", p_box)

summary_points <- long[internal == TRUE, .(
  mean = mean(value), median = median(value), sd = sd(value), observations = .N
), by = .(campaign, response, sampling.point.numeric, height)]
summary_points <- plot_limits[summary_points, on = "response"]
summary_points[, `:=`(
  mean.low.display = pmax(mean - sd, display.low),
  mean.high.display = pmin(mean + sd, display.high),
  median.low.display = pmax(median - sd, display.low),
  median.high.display = pmin(median + sd, display.high)
)]

p_mean_sd <- ggplot(summary_points,
  aes(sampling.point.numeric, mean, colour = height)) +
  geom_errorbar(aes(ymin = mean.low.display, ymax = mean.high.display),
                width = 0.15, linewidth = 0.3) +
  geom_point(size = 1.3) +
  facet_grid(response ~ campaign, scales = "free_y", labeller = facet_labels) +
  scale_colour_manual(values = height_colours, name = "Height") +
  scale_x_continuous(breaks = seq(1, 51, 2)) +
  labs(x = "Sampling point", y = "Mean +/- SD",
       title = "Mean and standard deviation by sampling point",
       subtitle = "SD bars clipped only for display to the common pooled response range") + common_theme
save_plot("Fig_02_mean_plus_minus_SD_campaign1_2", p_mean_sd)

p_median_sd <- ggplot(summary_points,
  aes(sampling.point.numeric, median, colour = height)) +
  geom_errorbar(aes(ymin = median.low.display, ymax = median.high.display),
                width = 0.15, linewidth = 0.3) +
  geom_point(size = 1.3) +
  facet_grid(response ~ campaign, scales = "free_y", labeller = facet_labels) +
  scale_colour_manual(values = height_colours, name = "Height") +
  scale_x_continuous(breaks = seq(1, 51, 2)) +
  labs(x = "Sampling point", y = "Median +/- SD",
       title = "Median and standard deviation by sampling point",
       subtitle = "SD bars clipped only for display to the common pooled response range") + common_theme
save_plot("Fig_03_median_plus_minus_SD_campaign1_2", p_median_sd)

representative_limits <- hourly_location[, .(
  display.low = quantile(median.difference, 0.005),
  display.high = quantile(median.difference, 0.995)
), by = response]
representative_plot <- representative_limits[hourly_location, on = "response"]
representative_plot <- representative_plot[
  median.difference >= display.low & median.difference <= display.high
]

p_representative <- ggplot(representative_plot,
  aes(factor(sampling.point.numeric), median.difference, fill = height)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.35) +
  geom_boxplot(outlier.shape = NA, width = 0.72, linewidth = 0.25) +
  facet_grid(response ~ campaign, scales = "free_y", labeller = facet_labels) +
  scale_fill_manual(values = height_colours, name = "Height") +
  labs(
    x = "Sampling point",
    y = "Hourly sampling-point median - hourly network median",
    title = "Median-based sampling-point representativeness",
    subtitle = "Common 0.5-99.5% display range per response across both campaigns"
  ) + common_theme
save_plot("Fig_04_median_representativeness_boxplots", p_representative)

entropy_c2 <- entropy[campaign == "Campaign 2"]
p_entropy_c2 <- ggplot(entropy_c2,
  aes(sampling.point.numeric, response, fill = entropy.normalised)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.2f", entropy.normalised)), size = 2.1) +
  scale_x_continuous(breaks = sort(unique(entropy_c2$sampling.point.numeric))) +
  scale_y_discrete(labels = response_labels) +
  scale_fill_viridis_c(limits = c(0, 1), name = "Normalised\nShannon entropy") +
  coord_fixed(ratio = 1.25) +
  labs(
    x = "Sampling point", y = NULL,
    title = "Campaign 2 Shannon entropy of concentrations and ratios"
  ) + common_theme +
  theme(axis.text.x = element_text(size = 7))
save_plot("Fig_05_Campaign2_Shannon_entropy_gases_and_ratios",
          p_entropy_c2, width = 13.5, height = 4.2)


###############################################################################
##### AUDIT
###############################################################################

audit <- data.table(
  decision = c(
    "physical_location_name", "ring_line_label", "south_label",
    "primary_representativeness", "ratio_units", "ratio_scale",
    "variability_display", "plot_display_range", "entropy_inputs",
    "excluded_metrics"
  ),
  definition = c(
    "sampling.point", "in is standardised to ring_in",
    "s and S are standardised to lowercase s",
    "hourly sampling-point median compared with equally weighted hourly network median",
    "dimensionless; not multiplied by 100",
    "shared within each response across Campaigns 1 and 2",
    "mean +/- SD and median +/- SD are both reported",
    "pooled 0.5-99.5% range per response; full values retained in tables",
    "hourly sampling-point medians for CO2, CH4, NH3 and three ratios",
    "CV and relative error are not calculated or plotted in this version"
  )
)
fwrite(audit, file.path(audit_dir, "analysis_decisions.csv"))
writeLines(capture.output(sessionInfo()), file.path(audit_dir, "sessionInfo.txt"))

message("Manuscript 3 plots and tables written to: ", output_root,
        " and ", plot_dir)
