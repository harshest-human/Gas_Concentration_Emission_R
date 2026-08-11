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
  library(lme4)
  library(emmeans)
})
emm_options(lmer.df = "asymptotic")

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
    sd = sd(v), cv.percent = 100 * sd(v) / mean(v), mad = mad(v),
    minimum = min(v), q1 = quantile(v, 0.25),
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
##### HYPOTHESIS 1: VERTICAL CONCENTRATION GRADIENT
###############################################################################

# The instruments visited sampling points sequentially. Clock-aligned two-hour
# blocks approximate a complete network cycle and are therefore the inferential
# unit; raw four-minute rows are not treated as simultaneous observations.
h1 <- long[internal == TRUE & response %chin% c("CO2", "CH4", "NH3")]
h1[, `:=`(
  block.2h = floor_date(DATE.TIME, "2 hours"),
  horizontal.position = (sampling.point.numeric - 1L) %/% 3L + 1L
)]
h1_blocks <- h1[, .(
  block.median = median(value), block.mean = mean(value), raw.rows = .N
), by = .(campaign, response, block.2h, sampling.point.numeric,
          horizontal.position, height)]

# Only blocks containing every expected height at a horizontal position enter
# a paired vertical comparison (3 heights in C1; top and bottom in C2).
expected_heights <- data.table(
  campaign = factor(c("Campaign 1", "Campaign 2"), levels = levels(x$campaign)),
  expected = c(3L, 2L)
)
h1_blocks <- expected_heights[h1_blocks, on = "campaign"]
h1_complete <- h1_blocks[, if (uniqueN(height) == first(expected)) .SD,
  by = .(campaign, response, block.2h, horizontal.position)]

height_descriptives <- h1_complete[, describe(block.median),
  by = .(campaign, response, height)]
fwrite(height_descriptives,
       file.path(table_dir, "Table_07_H1_height_descriptives_2h_blocks.csv"))

fit_height_model <- function(d) {
  d <- copy(d[is.finite(block.median) & block.median > 0])
  d[, height := droplevels(height)]
  if (nrow(d) < 30L || uniqueN(d$height) < 2L) return(NULL)
  model <- lmer(
    log(block.median) ~ height + factor(horizontal.position) +
      (1 | block.2h), data = d, REML = TRUE
  )
  means <- as.data.table(confint(emmeans(model, ~ height), type = "response"))
  setnames(means, "response", "estimated.concentration")
  contrasts <- as.data.table(summary(
    pairs(emmeans(model, ~ height), adjust = "holm"),
    infer = c(TRUE, TRUE), type = "response"
  ))
  list(means = means, contrasts = contrasts,
       singular = isSingular(model, tol = 1e-4), observations = nrow(d),
       blocks = uniqueN(d$block.2h))
}

model_results <- list()
for (cc in levels(x$campaign)) for (rr in c("CO2", "CH4", "NH3")) {
  result <- fit_height_model(h1_complete[campaign == cc & response == rr])
  if (is.null(result)) next
  result$means[, `:=`(campaign = cc, gas = rr,
                      model.observations = result$observations,
                      model.blocks = result$blocks,
                      singular.fit = result$singular)]
  result$contrasts[, `:=`(campaign = cc, response = rr,
                          model.observations = result$observations,
                          model.blocks = result$blocks,
                          singular.fit = result$singular)]
  model_results[[paste(cc, rr, sep = "_")]] <- result
}
h1_model_means <- rbindlist(lapply(model_results, `[[`, "means"), fill = TRUE)
h1_model_contrasts <- rbindlist(lapply(model_results, `[[`, "contrasts"), fill = TRUE)
fwrite(h1_model_means,
       file.path(table_dir, "Table_08_H1_mixed_model_height_estimates.csv"))
fwrite(h1_model_contrasts,
       file.path(table_dir, "Table_09_H1_mixed_model_height_contrasts.csv"))


###############################################################################
##### HYPOTHESIS 1: PRECISION ACROSS TEMPORAL RESOLUTIONS
###############################################################################

aggregate_precision <- function(d, unit) {
  z <- copy(d)
  if (unit == "2-hour") z[, period := floor_date(DATE.TIME, "2 hours")]
  if (unit == "day") z[, period := floor_date(DATE.TIME, "day")]
  if (unit == "week") z[, period := floor_date(DATE.TIME, "week", week_start = 1)]
  if (unit == "campaign") z[, period := min(DATE.TIME)]
  location_period <- z[, .(estimate = median(value), raw.rows = .N),
    by = .(campaign, response, period, sampling.point.numeric)]
  location_period[, .(
    sampling.points = .N,
    spatial.median = median(estimate),
    spatial.mad = mad(estimate),
    robust.cv.percent = 100 * mad(estimate) / median(estimate),
    spatial.minimum = min(estimate), spatial.maximum = max(estimate),
    spatial.range.percent = 100 * (max(estimate) - min(estimate)) / median(estimate)
  ), by = .(campaign, response, period)][, temporal.resolution := unit]
}

precision <- rbindlist(lapply(
  c("2-hour", "day", "week", "campaign"),
  function(unit) aggregate_precision(h1, unit)
), fill = TRUE)
fwrite(precision,
       file.path(table_dir, "Table_10_H1_precision_by_temporal_resolution.csv"))
precision_summary <- precision[, .(
  periods = .N, median.sampling.points = as.numeric(median(sampling.points)),
  median.robust.cv.percent = median(robust.cv.percent, na.rm = TRUE),
  q1.robust.cv.percent = quantile(robust.cv.percent, 0.25, na.rm = TRUE),
  q3.robust.cv.percent = quantile(robust.cv.percent, 0.75, na.rm = TRUE),
  median.spatial.range.percent = median(spatial.range.percent, na.rm = TRUE)
), by = .(campaign, response, temporal.resolution)]
fwrite(precision_summary,
       file.path(table_dir, "Table_11_H1_precision_resolution_summary.csv"))


###############################################################################
##### HYPOTHESIS 1: SPATIAL ENTROPY AND PCA DIAGNOSTICS
###############################################################################

spatial_entropy <- h1_complete[, {
  v <- block.median[is.finite(block.median) & block.median > 0]
  p <- v / sum(v)
  .(sampling.points = length(v), entropy = -sum(p * log(p)),
    mixing.index = if (length(v) > 1L) -sum(p * log(p)) / log(length(v)) else NA_real_)
}, by = .(campaign, response, block.2h)]
fwrite(spatial_entropy,
       file.path(table_dir, "Table_12_H1_spatial_entropy_2h_blocks.csv"))

run_pca <- function(d, campaign_name, gas_name) {
  w <- dcast(d[campaign == campaign_name & response == gas_name],
             block.2h ~ sampling.point.numeric, value.var = "block.median")
  if (ncol(w) < 4L) return(NULL)
  matrix_data <- as.matrix(w[, -"block.2h"])
  keep_columns <- colSums(is.finite(matrix_data)) >= 10L
  matrix_data <- matrix_data[, keep_columns, drop = FALSE]
  complete <- complete.cases(matrix_data)
  if (sum(complete) < 10L || ncol(matrix_data) < 3L) return(NULL)
  fit <- prcomp(matrix_data[complete, , drop = FALSE], center = TRUE, scale. = TRUE)
  loadings <- data.table(
    sampling.point.numeric = as.integer(rownames(fit$rotation)),
    PC1 = fit$rotation[, 1], PC2 = fit$rotation[, 2]
  )
  loadings[, height := c("top", "middle", "bottom")[(sampling.point.numeric - 1L) %% 3L + 1L]]
  loadings[, `:=`(campaign = campaign_name, response = gas_name,
                  complete.blocks = sum(complete),
                  PC1.variance = summary(fit)$importance[2, 1],
                  PC2.variance = summary(fit)$importance[2, 2])]
  loadings
}
# PCA is a three-height structural diagnostic and is consequently restricted
# to Campaign 1. Campaign 2 deliberately omitted the middle height.
pca_loadings <- rbindlist(lapply(c("CO2", "CH4", "NH3"), function(rr) {
  run_pca(h1_complete, "Campaign 1", rr)
}), fill = TRUE)
fwrite(pca_loadings, file.path(table_dir, "Table_13_H1_PCA_loadings.csv"))


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

# Four pooled quantile bins are defined for each campaign-response combination
# and then applied identically to every SP. This is robust to extreme gas-ratio
# values while preserving a common comparison basis among locations.
entropy <- hourly_location[, {
  pooled_breaks <- unique(as.numeric(quantile(
    location.median, probs = seq(0, 1, 0.25), na.rm = TRUE
  )))
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

# Spatial entropy describes the distribution across sampling points within an
# aligned two-hour network block. Normalisation by log(N) permits comparison
# between campaigns with different numbers of available sampling points.
entropy_2h_locations <- long[internal == TRUE, .(
  location.median = median(value), raw.rows = .N
), by = .(
  campaign, response, block.2h = floor_date(DATE.TIME, "2 hours"),
  sampling.point.numeric
)]
spatial_entropy_all <- entropy_2h_locations[, {
  v <- location.median[is.finite(location.median) & location.median > 0]
  p <- v / sum(v)
  h <- -sum(p * log(p))
  .(
    sampling.points = length(v), entropy = h,
    entropy.maximum = log(length(v)),
    entropy.normalised = if (length(v) > 1L) h / log(length(v)) else NA_real_,
    information.deficit.percent = if (length(v) > 1L) 100 * (1 - h / log(length(v))) else NA_real_
  )
}, by = .(campaign, response, block.2h)]
spatial_entropy_all[, expected.sampling.points := max(sampling.points),
                    by = .(campaign, response)]
spatial_entropy_all[, coverage.percent := 100 * sampling.points /
                      expected.sampling.points]
spatial_entropy_all[, valid.coverage := coverage.percent >= 80]
fwrite(spatial_entropy_all, file.path(table_dir,
  "Table_14_spatial_Shannon_entropy_gases_and_ratios_2h.csv"))

entropy_summary <- spatial_entropy_all[valid.coverage == TRUE, .(
  blocks = .N,
  median.sampling.points = as.numeric(median(sampling.points)),
  median.entropy.normalised = median(entropy.normalised, na.rm = TRUE),
  q1.entropy.normalised = quantile(entropy.normalised, 0.25, na.rm = TRUE),
  q3.entropy.normalised = quantile(entropy.normalised, 0.75, na.rm = TRUE),
  median.information.deficit.percent = median(information.deficit.percent, na.rm = TRUE)
), by = .(campaign, response)]
fwrite(entropy_summary, file.path(table_dir,
  "Table_15_spatial_Shannon_entropy_summary.csv"))


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

p_entropy_temporal <- ggplot(entropy,
  aes(sampling.point.numeric, response, fill = entropy.normalised)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.2f", entropy.normalised)), size = 2.1) +
  facet_grid(campaign ~ ., scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = sort(unique(entropy$sampling.point.numeric))) +
  scale_y_discrete(labels = response_labels) +
  scale_fill_viridis_c(limits = c(0, 1), name = "Normalised\nShannon entropy") +
  labs(
    x = "Sampling point", y = NULL,
    title = "Temporal Shannon entropy of concentrations and ratios",
    subtitle = "Common pooled bins within each campaign-response combination"
  ) + common_theme +
  theme(axis.text.x = element_text(size = 7))
save_plot("Fig_05_temporal_Shannon_entropy_campaign1_2",
          p_entropy_temporal, width = 15, height = 8)

h1_model_means[, gas := factor(gas, levels = c("CO2", "CH4", "NH3"))]
p_h1_height <- ggplot(h1_model_means,
  aes(height, estimated.concentration, colour = height)) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                width = 0.14, linewidth = 0.45) +
  geom_point(size = 2) +
  facet_grid(gas ~ campaign, scales = "free_y") +
  scale_colour_manual(values = height_colours, guide = "none") +
  labs(
    x = "Sampling height", y = "Model-estimated concentration (95% CI)",
    title = "Hypothesis 1: vertical concentration gradient",
    subtitle = paste(
      "Campaign 1 tests top, middle and bottom; Campaign 2 validates",
      "the retained top-bottom contrast"
    )
  ) + common_theme
save_plot("Fig_06_H1_mixed_model_vertical_gradient", p_h1_height,
          width = 9, height = 9)

precision_summary[, temporal.resolution := factor(
  temporal.resolution, levels = c("2-hour", "day", "week", "campaign")
)]
p_precision <- ggplot(precision_summary,
  aes(temporal.resolution, median.robust.cv.percent,
      group = response, colour = response)) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 1.8) +
  facet_wrap(~ campaign) +
  labs(
    x = "Temporal aggregation", y = "Median spatial robust CV (%)",
    colour = "Gas", title = "Spatial precision across temporal resolutions",
    subtitle = "Robust CV = 100 x MAD of sampling-point medians / network median"
  ) + common_theme
save_plot("Fig_07_H1_precision_temporal_resolution", p_precision,
          width = 9, height = 5.5)

p_pca <- ggplot(pca_loadings,
  aes(PC1, PC2, colour = height, label = sampling.point.numeric)) +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_vline(xintercept = 0, colour = "grey80") +
  geom_point(size = 1.8) +
  geom_text(nudge_y = 0.004, size = 2.2, check_overlap = TRUE) +
  facet_wrap(~ response, scales = "free") +
  scale_colour_manual(values = height_colours, name = "Height") +
  labs(x = "PC1 loading", y = "PC2 loading",
       title = "Campaign 1 PCA loadings by sampling height") + common_theme
save_plot("Fig_08_H1_Campaign1_PCA_height_diagnostic", p_pca,
          width = 10, height = 4.5)

p_spatial_entropy <- ggplot(spatial_entropy_all[valid.coverage == TRUE],
  aes(block.2h, entropy.normalised, colour = response)) +
  geom_line(linewidth = 0.25, alpha = 0.65) +
  facet_grid(response ~ campaign, scales = "free_x", labeller = facet_labels) +
  scale_y_continuous(breaks = seq(0.94, 1, 0.02)) +
  coord_cartesian(ylim = c(0.94, 1)) +
  labs(
    x = "Two-hour sampling-cycle block", y = "Normalised spatial entropy",
    title = "Spatial Shannon entropy over time",
    subtitle = "Blocks with >=80% network coverage; y-axis zoomed to show differences near one"
  ) + common_theme + theme(legend.position = "none")
save_plot("Fig_09_spatial_Shannon_entropy_time_series_campaign1_2",
          p_spatial_entropy, width = 15, height = 12)

p_spatial_entropy_distribution <- ggplot(spatial_entropy_all[valid.coverage == TRUE],
  aes(response, entropy.normalised, fill = response)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.35) +
  facet_wrap(~ campaign) +
  scale_x_discrete(labels = response_labels) +
  scale_y_continuous(breaks = seq(0.94, 1, 0.02)) +
  coord_cartesian(ylim = c(0.94, 1)) +
  labs(
    x = NULL, y = "Normalised spatial entropy",
    title = "Distribution of spatial Shannon entropy",
    subtitle = "Aligned two-hour blocks with >=80% network coverage; y-axis zoomed"
  ) + common_theme +
  theme(legend.position = "none", axis.text.x = element_text(size = 7))
save_plot("Fig_10_spatial_Shannon_entropy_distribution_campaign1_2",
          p_spatial_entropy_distribution, width = 11, height = 5.5)


###############################################################################
##### AUDIT
###############################################################################

audit <- data.table(
  decision = c(
    "physical_location_name", "ring_line_label", "south_label",
    "primary_representativeness", "ratio_units", "ratio_scale",
    "variability_display", "plot_display_range", "entropy_inputs",
    "spatial_entropy", "temporal_entropy",
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
    "two-hour sampling-point medians for CO2, CH4, NH3 and three ratios",
    "normalised across sampling points within aligned two-hour blocks using log(N)",
    "normalised across time using common pooled campaign-response bins",
    "CV and relative error are not calculated or plotted in this version"
  )
)
fwrite(audit, file.path(audit_dir, "analysis_decisions.csv"))
writeLines(capture.output(sessionInfo()), file.path(audit_dir, "sessionInfo.txt"))

message("Manuscript 3 plots and tables written to: ", output_root,
        " and ", plot_dir)
