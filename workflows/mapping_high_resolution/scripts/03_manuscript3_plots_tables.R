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
  library(patchwork)
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
# Preserve the complete prepared gas table for the CO2-balance calculation.
# The plotting/statistical copy below converts non-positive values to NA, but
# the emissions audit must retain them so validity rules can be applied later.
x_emission_source <- copy(x)
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
# Campaign 2 FTIR2 operated only during three short sessions. Its internal
# points 49 and 51 are retained for transparent descriptive auditing but are
# ineligible for Campaign 2 inferential statistics. The primary Campaign 2
# analysis therefore uses the 32 continuously cycled CRDS sampling points.
x[, analysis.eligible := internal &
    !(campaign == "Campaign 2" & !grepl("^CRDS", analyser))]
x[, height := fifelse(
  !internal, "reference",
  c("top", "middle", "bottom")[(sampling.point.numeric - 1L) %% 3L + 1L]
)]
x[, height := factor(height, levels = names(height_colours))]
x[, campaign := factor(campaign, levels = c("Campaign 1", "Campaign 2"))]

long <- melt(
  x,
  id.vars = c("DATE.TIME", "campaign", "analyser", "sampling.point",
              "sampling.point.numeric", "internal", "analysis.eligible",
              "height"),
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
         internal, analysis.eligible, height)]
setorder(location_descriptives, campaign, response, sampling.point.numeric)
fwrite(location_descriptives,
       file.path(table_dir, "Table_01_sampling_point_descriptives.csv"))

campaign_descriptives <- long[analysis.eligible == TRUE, describe(value),
  by = .(campaign, response)]
fwrite(campaign_descriptives,
       file.path(table_dir, "Table_02_campaign_descriptives.csv"))

eligibility_audit <- unique(x[, .(
  campaign, analyser, sampling.point, sampling.point.numeric,
  internal, analysis.eligible
)])
setorder(eligibility_audit, campaign, analysis.eligible,
         sampling.point.numeric, analyser)
fwrite(eligibility_audit,
       file.path(table_dir, "Table_00_sampling_point_analysis_eligibility.csv"))


###############################################################################
##### HYPOTHESIS 1: VERTICAL CONCENTRATION GRADIENT
###############################################################################

# The instruments visited sampling points sequentially. Clock-aligned two-hour
# blocks approximate a complete network cycle and are therefore the inferential
# unit; raw four-minute rows are not treated as simultaneous observations.
h1 <- long[analysis.eligible == TRUE & response %chin% c("CO2", "CH4", "NH3")]
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
hourly_location <- long[analysis.eligible == TRUE,
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
entropy_2h_locations <- long[analysis.eligible == TRUE, .(
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
##### STANDARD CO2-BALANCE VENTILATION AND EMISSION ESTIMATES
###############################################################################

# Every aligned row is retained. Invalid or unavailable calculations remain in
# the CSV with the function's validity flags and NA estimates; no subsequent
# filtering is performed here.
co2_balance_function <- file.path(
  dirname(workflow), "utils", "indirect.CO2.balance function.R"
)
if (!file.exists(co2_balance_function)) {
  stop("CO2-balance function was not found: ", co2_balance_function)
}
source(co2_balance_function, local = environment())
if (!exists("indirect.CO2.balance.st", mode = "function")) {
  stop("indirect.CO2.balance.st() was not defined by: ", co2_balance_function)
}

median_keep_all <- function(v) {
  # Missing values cannot contribute to a median, but zero and negative values
  # are deliberately retained. Return NA only when the entire group is missing.
  if (all(is.na(v))) return(NA_real_)
  median(v, na.rm = TRUE)
}

emission_2h <- x_emission_source[, .(
  CO2 = median_keep_all(CO2),
  CH4 = median_keep_all(CH4),
  NH3 = median_keep_all(NH3),
  gas.raw.rows = .N,
  gas.nonmissing.CO2 = sum(!is.na(CO2)),
  gas.nonmissing.CH4 = sum(!is.na(CH4)),
  gas.nonmissing.NH3 = sum(!is.na(NH3))
), by = .(
  campaign, analyser,
  block.2h = floor_date(DATE.TIME, "2 hours"),
  sampling.point
)]
emission_2h[, sampling.point.numeric := suppressWarnings(as.integer(sampling.point))]

# Campaign 1 used `s` as the outside line. Campaign 2 FTIR2 selector position 3
# maps to sampling point 52 and is used here as the campaign outside reference.
emission_2h[, reference.role := fcase(
  campaign == "Campaign 1" & tolower(sampling.point) == "s", "outside",
  campaign == "Campaign 2" & sampling.point == "52", "outside",
  !is.na(sampling.point.numeric) & sampling.point.numeric %between% c(1L, 51L),
  "inside_sampling_point",
  default = "unclassified"
)]
emission_2h[, analysis.eligible :=
  reference.role == "inside_sampling_point" &
    !(campaign == "Campaign 2" & !grepl("^CRDS", analyser))]

outside_2h <- emission_2h[reference.role == "outside", .(
  CO2_ppm_out = median_keep_all(CO2),
  CH4_ppm_out = median_keep_all(CH4),
  NH3_ppm_out = median_keep_all(NH3),
  outside.raw.rows = sum(gas.raw.rows),
  outside.source = paste(sort(unique(sampling.point)), collapse = "+")
), by = .(campaign, block.2h)]

inside_points <- emission_2h[reference.role == "inside_sampling_point", .(
  campaign, block.2h, estimate.level = "sampling_point", sampling.point,
  analyser, analysis.eligible,
  sampling.points.in.estimate = 1L,
  CO2_ppm_in = CO2, CH4_ppm_in = CH4, NH3_ppm_in = NH3,
  inside.raw.rows = gas.raw.rows,
  inside.nonmissing.CO2 = gas.nonmissing.CO2,
  inside.nonmissing.CH4 = gas.nonmissing.CH4,
  inside.nonmissing.NH3 = gas.nonmissing.NH3
)]

inside_network <- emission_2h[analysis.eligible == TRUE, .(
  estimate.level = "network_median",
  sampling.point = "network_median",
  analyser = "eligible_network",
  analysis.eligible = TRUE,
  sampling.points.in.estimate = uniqueN(sampling.point),
  CO2_ppm_in = median_keep_all(CO2),
  CH4_ppm_in = median_keep_all(CH4),
  NH3_ppm_in = median_keep_all(NH3),
  inside.raw.rows = sum(gas.raw.rows),
  inside.nonmissing.CO2 = sum(!is.na(CO2)),
  inside.nonmissing.CH4 = sum(!is.na(CH4)),
  inside.nonmissing.NH3 = sum(!is.na(NH3))
), by = .(campaign, block.2h)]

emission_input <- rbindlist(
  list(inside_points, inside_network), use.names = TRUE, fill = TRUE
)
emission_input <- outside_2h[emission_input,
  on = .(campaign, block.2h)]

weather_path <- file.path(
  workflow, "clean_data", "manuscript3_weather",
  "dwd_potsdam_03987_hourly_2024_campaign_periods.csv"
)
if (!file.exists(weather_path)) {
  stop("DWD weather input was not found: ", weather_path)
}
weather_emission <- fread(weather_path)
weather_emission[, DATE.TIME := with_tz(
  ymd_hms(DATE.TIME, tz = "UTC", quiet = TRUE), "Europe/Berlin"
)]
weather_2h <- weather_emission[, .(
  temp_dwd = if (all(is.na(temperature_C))) NA_real_ else
    mean(temperature_C, na.rm = TRUE),
  pressure_dwd_Pa = NA_real_,
  dwd.hourly.rows = .N,
  dwd.temperature.rows = sum(!is.na(temperature_C))
), by = .(block.2h = floor_date(DATE.TIME, "2 hours"))]

emission_input <- weather_2h[emission_input, on = "block.2h"]
setorder(emission_input, campaign, block.2h, estimate.level, sampling.point)

emission_result <- indirect.CO2.balance.st(
  as.data.frame(emission_input),
  dwd_temp_col = "temp_dwd",
  pressure_col = NULL,
  n_dairy_cows = 58,
  cow_weight_kg = 500,
  milk_kg_cow_d = 31.21864407,
  pregnancy_day = 119.5510204,
  min_delta_co2_ppm = 0,
  low_delta_co2_ppm = 200,
  annualisation_hours = 8760
)
emission_result <- as.data.table(emission_result)
emission_result[, campaign.result.role := fcase(
  campaign == "Campaign 1", "primary_campaign_result",
  campaign == "Campaign 2" & estimate.level == "sampling_point" &
    analysis.eligible == FALSE, "retained_audit_excluded_from_inference",
  campaign == "Campaign 2", "exploratory_only_incomplete_outside_reference",
  default = "unclassified"
)]
emission_result[, calculation.row := .I]
setcolorder(emission_result, c(
  "calculation.row", "campaign", "block.2h", "estimate.level",
  "sampling.point", "sampling.points.in.estimate"
))

fwrite(
  emission_result,
  file.path(output_root,
            "H2_standard_CO2_balance_ventilation_emissions_all_rows.csv"),
  dateTimeAs = "write.csv"
)
fwrite(
  emission_result[estimate.level == "network_median"],
  file.path(output_root,
            "H2_standard_CO2_balance_network_median_all_rows.csv"),
  dateTimeAs = "write.csv"
)

emission_audit <- emission_result[, .(
  rows = .N,
  rows.with.outside = sum(is.finite(CO2_ppm_out)),
  rows.with.temperature = sum(is.finite(temp_dwd)),
  rows.positive.delta.CO2 = sum(delta_co2_positive, na.rm = TRUE),
  rows.low.delta.CO2 = sum(delta_co2_lt_200ppm, na.rm = TRUE),
  rows.finite.Q = sum(is.finite(Q_vent_m3_h_barn)),
  rows.negative.CH4.enhancement = sum(delta_ch4_negative, na.rm = TRUE),
  rows.negative.NH3.enhancement = sum(delta_nh3_negative, na.rm = TRUE),
  rows.finite.CH4.emission = sum(is.finite(e_CH4_kg_year_LU)),
  rows.finite.NH3.emission = sum(is.finite(e_NH3_kg_year_LU))
), by = .(campaign, estimate.level)]
fwrite(emission_audit,
       file.path(table_dir, "Table_16_H2_CO2_balance_unfiltered_audit.csv"))

campaign2_exclusion_audit <- unique(emission_result[
  campaign == "Campaign 2" & estimate.level == "sampling_point",
  .(sampling.point, analyser, analysis.eligible, campaign.result.role)
])
setorder(campaign2_exclusion_audit, analysis.eligible, sampling.point)
fwrite(campaign2_exclusion_audit,
       file.path(table_dir, "Table_17_Campaign2_sampling_point_eligibility.csv"))


###############################################################################
##### CAMPAIGN 1 SAMPLING CONFIGURATION AND DENSITY EFFECTS ON EMISSIONS
###############################################################################

add_emission_validity <- function(z) {
  z <- as.data.table(z)
  z[, delta.CO2.validity := fcase(
    !is.finite(delta_CO2_ppm), "missing",
    delta_CO2_ppm <= 0, "non-positive",
    delta_CO2_ppm < 70, "below-70-ppm",
    delta_CO2_ppm < 200, "70-to-199-ppm-sensitivity",
    default = "at-least-200-ppm-primary"
  )]
  z[, primary.valid := delta.CO2.validity == "at-least-200-ppm-primary" &
      is.finite(Q_vent_m3_h_LU)]
  z[, sensitivity.valid := delta_CO2_ppm >= 70 &
      is.finite(Q_vent_m3_h_LU)]
  z
}

emission_result <- add_emission_validity(emission_result)
# Rewrite the two unfiltered files with the validity classification included.
fwrite(emission_result, file.path(
  output_root, "H2_standard_CO2_balance_ventilation_emissions_all_rows.csv"
), dateTimeAs = "write.csv")
fwrite(emission_result[estimate.level == "network_median"], file.path(
  output_root, "H2_standard_CO2_balance_network_median_all_rows.csv"
), dateTimeAs = "write.csv")

c1_config_source <- copy(emission_2h[
  campaign == "Campaign 1" & reference.role == "inside_sampling_point"
])
c1_config_source[, `:=`(
  height = c("top", "middle", "bottom")[(sampling.point.numeric - 1L) %% 3L + 1L],
  horizontal.position = (sampling.point.numeric - 1L) %/% 3L + 1L
)]

configuration_membership <- rbindlist(list(
  c1_config_source[, .(block.2h, sampling.point, CO2, CH4, NH3,
                       configuration = "all_51")],
  c1_config_source[height == "top", .(block.2h, sampling.point, CO2, CH4, NH3,
                                      configuration = "top_17")],
  c1_config_source[height == "middle", .(block.2h, sampling.point, CO2, CH4, NH3,
                                         configuration = "middle_17")],
  c1_config_source[height == "bottom", .(block.2h, sampling.point, CO2, CH4, NH3,
                                         configuration = "bottom_17")],
  c1_config_source[height != "middle", .(block.2h, sampling.point, CO2, CH4, NH3,
                                         configuration = "top_bottom_34")],
  c1_config_source[height != "middle" & !sampling.point.numeric %in% c(19L, 40L),
    .(block.2h, sampling.point, CO2, CH4, NH3,
      configuration = "campaign2_like_32")]
), use.names = TRUE)

configuration_input <- configuration_membership[, .(
  sampling.points.in.estimate = uniqueN(sampling.point),
  CO2_ppm_in = median_keep_all(CO2),
  CH4_ppm_in = median_keep_all(CH4),
  NH3_ppm_in = median_keep_all(NH3)
), by = .(block.2h, configuration)]
configuration_input[, campaign := "Campaign 1"]
configuration_input <- outside_2h[campaign == "Campaign 1"][
  configuration_input, on = .(campaign, block.2h)
]
configuration_input <- weather_2h[configuration_input, on = "block.2h"]

configuration_result <- indirect.CO2.balance.st(
  as.data.frame(configuration_input), dwd_temp_col = "temp_dwd",
  n_dairy_cows = 58, cow_weight_kg = 500,
  milk_kg_cow_d = 31.21864407, pregnancy_day = 119.5510204,
  min_delta_co2_ppm = 0, low_delta_co2_ppm = 200,
  annualisation_hours = 8760
)
configuration_result <- add_emission_validity(configuration_result)
configuration_result <- as.data.table(configuration_result)

baseline <- configuration_result[configuration == "all_51", .(
  block.2h,
  Q.reference = Q_vent_m3_h_LU,
  CH4.reference = e_CH4_ghLU,
  NH3.reference = e_NH3_ghLU,
  reference.primary.valid = primary.valid
)]
configuration_result <- baseline[configuration_result, on = "block.2h"]
configuration_result[, `:=`(
  Q.relative.error.percent = 100 * (Q_vent_m3_h_LU - Q.reference) / Q.reference,
  CH4.relative.error.percent = 100 * (e_CH4_ghLU - CH4.reference) / CH4.reference,
  NH3.relative.error.percent = 100 * (e_NH3_ghLU - NH3.reference) / NH3.reference,
  CH4_CO2 = CH4_ppm_in / CO2_ppm_in,
  NH3_CO2 = NH3_ppm_in / CO2_ppm_in,
  NH3_CH4 = NH3_ppm_in / CH4_ppm_in
)]
fwrite(configuration_result, file.path(
  output_root, "H2_Campaign1_sampling_configuration_emissions_all_rows.csv"
), dateTimeAs = "write.csv")

configuration_summary <- configuration_result[
  primary.valid == TRUE & reference.primary.valid == TRUE,
  .(
    valid.blocks = .N,
    median.sampling.points = as.numeric(median(sampling.points.in.estimate)),
    Q.median = median(Q_vent_m3_h_LU),
    Q.abs.RE.median = median(abs(Q.relative.error.percent), na.rm = TRUE),
    Q.abs.RE.q95 = quantile(abs(Q.relative.error.percent), 0.95, na.rm = TRUE),
    CH4.median = median(e_CH4_ghLU, na.rm = TRUE),
    CH4.abs.RE.median = median(abs(CH4.relative.error.percent), na.rm = TRUE),
    CH4.abs.RE.q95 = quantile(abs(CH4.relative.error.percent), 0.95, na.rm = TRUE),
    NH3.median = median(e_NH3_ghLU, na.rm = TRUE),
    NH3.abs.RE.median = median(abs(NH3.relative.error.percent), na.rm = TRUE),
    NH3.abs.RE.q95 = quantile(abs(NH3.relative.error.percent), 0.95, na.rm = TRUE)
  ), by = configuration
]
fwrite(configuration_summary, file.path(
  table_dir, "Table_18_Campaign1_sampling_configuration_emission_effects.csv"
))

# Balanced-density simulation: randomly select horizontal positions and retain
# all three heights at each selected position. Each subset is fixed across time,
# representing a physically installed reduced network rather than resampling
# points independently at every block.
set.seed(20260812)
density_levels <- c(1L, 2L, 4L, 8L, 12L, 17L)
density_design <- rbindlist(lapply(density_levels, function(k) {
  reps <- if (k == 17L) 1L else 100L
  rbindlist(lapply(seq_len(reps), function(rep_id) data.table(
    horizontal.position = sample(1:17, k),
    horizontal.positions = k,
    sampling.points = 3L * k,
    replicate = rep_id
  )))
}))

# The explicit loop is easier to audit than a cartesian expansion and keeps
# memory bounded for the 100 fixed-network replicates per density.
density_networks <- unique(
  density_design[, .(horizontal.positions, sampling.points, replicate)]
)
density_input_list <- vector("list", nrow(density_networks))
for (i in seq_len(nrow(density_networks))) {
  design_i <- density_networks[i]
  selected <- density_design[
    horizontal.positions == design_i$horizontal.positions &
      replicate == design_i$replicate,
    horizontal.position
  ]
  density_input_list[[i]] <- c1_config_source[
    horizontal.position %in% selected,
    .(
      CO2_ppm_in = median_keep_all(CO2),
      CH4_ppm_in = median_keep_all(CH4),
      NH3_ppm_in = median_keep_all(NH3)
    ), by = block.2h
  ][, `:=`(
    horizontal.positions = design_i$horizontal.positions,
    sampling.points.in.estimate = design_i$sampling.points,
    replicate = design_i$replicate
  )]
}
density_input <- rbindlist(density_input_list, use.names = TRUE)
density_input[, campaign := "Campaign 1"]
density_input <- outside_2h[campaign == "Campaign 1"][
  density_input, on = .(campaign, block.2h)
]
density_input <- weather_2h[density_input, on = "block.2h"]
density_result <- indirect.CO2.balance.st(
  as.data.frame(density_input), dwd_temp_col = "temp_dwd",
  n_dairy_cows = 58, cow_weight_kg = 500,
  milk_kg_cow_d = 31.21864407, pregnancy_day = 119.5510204,
  min_delta_co2_ppm = 0, low_delta_co2_ppm = 200,
  annualisation_hours = 8760
)
density_result <- add_emission_validity(density_result)
density_result <- baseline[as.data.table(density_result), on = "block.2h"]
density_result[, `:=`(
  Q.abs.RE = abs(100 * (Q_vent_m3_h_LU - Q.reference) / Q.reference),
  CH4.abs.RE = abs(100 * (e_CH4_ghLU - CH4.reference) / CH4.reference),
  NH3.abs.RE = abs(100 * (e_NH3_ghLU - NH3.reference) / NH3.reference)
)]
fwrite(density_result, file.path(
  output_root, "H2_Campaign1_sampling_density_simulation_all_rows.csv"
), dateTimeAs = "write.csv")

density_summary <- melt(
  density_result[primary.valid == TRUE & reference.primary.valid == TRUE],
  id.vars = c("horizontal.positions", "sampling.points.in.estimate", "replicate"),
  measure.vars = c("Q.abs.RE", "CH4.abs.RE", "NH3.abs.RE"),
  variable.name = "outcome", value.name = "absolute.relative.error.percent"
)[is.finite(absolute.relative.error.percent), .(
  median.absolute.RE = median(absolute.relative.error.percent),
  q1.absolute.RE = quantile(absolute.relative.error.percent, 0.25),
  q3.absolute.RE = quantile(absolute.relative.error.percent, 0.75),
  q95.absolute.RE = quantile(absolute.relative.error.percent, 0.95)
), by = .(horizontal.positions, sampling.points.in.estimate, replicate, outcome)]
fwrite(density_summary, file.path(
  table_dir, "Table_19_Campaign1_sampling_density_simulation.csv"
))

ratio_outcome_correlations <- melt(
  configuration_result[
    configuration == "all_51" & primary.valid == TRUE,
    .(block.2h, CH4_CO2, NH3_CO2, NH3_CH4,
      Q_vent_m3_h_LU, e_CH4_ghLU, e_NH3_ghLU)
  ],
  id.vars = c("block.2h", "Q_vent_m3_h_LU", "e_CH4_ghLU", "e_NH3_ghLU"),
  measure.vars = c("CH4_CO2", "NH3_CO2", "NH3_CH4"),
  variable.name = "ratio", value.name = "ratio.value"
)
ratio_outcome_correlations <- melt(
  ratio_outcome_correlations,
  id.vars = c("block.2h", "ratio", "ratio.value"),
  measure.vars = c("Q_vent_m3_h_LU", "e_CH4_ghLU", "e_NH3_ghLU"),
  variable.name = "outcome", value.name = "outcome.value"
)[is.finite(ratio.value) & ratio.value > 0 &
    is.finite(outcome.value) & outcome.value > 0, .(
  observations = .N,
  spearman.rho = cor(ratio.value, outcome.value, method = "spearman"),
  pearson.r = cor(log(ratio.value), log(outcome.value), method = "pearson")
), by = .(ratio, outcome)]
fwrite(ratio_outcome_correlations, file.path(
  table_dir, "Table_20_Campaign1_ratio_ventilation_emission_correlations.csv"
))


###############################################################################
##### REPRESENTATIVE HEIGHT/SP AND GAS-COUPLING DIAGNOSTICS
###############################################################################

block_location <- long[analysis.eligible == TRUE, .(
  estimate = median(value)
), by = .(
  campaign, response, block.2h = floor_date(DATE.TIME, "2 hours"),
  sampling.point.numeric, height
)]

# Leave-one-SP-out reference prevents a candidate point contributing to its own
# benchmark. It is evaluated in Campaign 1, where all 51 points were installed.
c1_loo <- block_location[campaign == "Campaign 1", {
  v <- estimate
  ref <- vapply(seq_along(v), function(i) {
    other <- v[-i]
    if (all(!is.finite(other))) NA_real_ else median(other, na.rm = TRUE)
  }, numeric(1))
  .(sampling.point.numeric, height, estimate,
    reference.leave.one.out = ref)
}, by = .(response, block.2h)]
c1_loo[, `:=`(
  difference = estimate - reference.leave.one.out,
  absolute.percent.error = 100 * abs(estimate - reference.leave.one.out) /
    reference.leave.one.out
)]

sp_representativeness <- c1_loo[
  is.finite(estimate) & is.finite(reference.leave.one.out), {
    fit <- lm(estimate ~ reference.leave.one.out)
    .(
      blocks = .N,
      bias = median(difference),
      normalised.MAE.percent = median(absolute.percent.error),
      RMSE = sqrt(mean(difference^2)),
      q90.absolute.percent.error = quantile(absolute.percent.error, 0.90),
      spearman.rho = cor(estimate, reference.leave.one.out, method = "spearman"),
      intercept = unname(coef(fit)[1]), slope = unname(coef(fit)[2])
    )
  }, by = .(response, sampling.point.numeric, height)
]
sp_representativeness[, `:=`(
  rank.MAE = frank(normalised.MAE.percent, ties.method = "average") / .N,
  rank.bias = frank(abs(bias), ties.method = "average") / .N,
  rank.slope = frank(abs(1 - slope), ties.method = "average") / .N,
  rank.correlation = frank(1 - spearman.rho, ties.method = "average") / .N
), by = response]
sp_representativeness[, response.score := rowMeans(.SD),
  .SDcols = c("rank.MAE", "rank.bias", "rank.slope", "rank.correlation")]
sp_composite <- sp_representativeness[, .(
  composite.score = mean(response.score),
  median.normalised.MAE.percent = median(normalised.MAE.percent),
  median.absolute.bias = median(abs(bias)),
  minimum.covered.blocks = min(blocks)
), by = .(sampling.point.numeric, height)]
setorder(sp_composite, composite.score)
sp_composite[, overall.rank := .I]
fwrite(sp_representativeness, file.path(
  table_dir, "Table_21_Campaign1_SP_representativeness_by_response.csv"
))
fwrite(sp_composite, file.path(
  table_dir, "Table_22_Campaign1_SP_composite_representativeness.csv"
))

# Each height is compared with the median of the other two heights, avoiding a
# self-including whole-network reference.
c1_height_block <- block_location[campaign == "Campaign 1", .(
  height.estimate = median(estimate),
  height.points = uniqueN(sampling.point.numeric)
), by = .(response, block.2h, height)]
c1_height_loo <- c1_height_block[, {
  v <- height.estimate
  ref <- vapply(seq_along(v), function(i) median(v[-i], na.rm = TRUE), numeric(1))
  .(height, height.estimate, reference.other.heights = ref)
}, by = .(response, block.2h)]
c1_height_loo[, `:=`(
  difference = height.estimate - reference.other.heights,
  absolute.percent.error = 100 * abs(height.estimate - reference.other.heights) /
    reference.other.heights
)]
height_representativeness <- c1_height_loo[, .(
  blocks = .N,
  median.bias = median(difference),
  normalised.MAE.percent = median(absolute.percent.error),
  RMSE = sqrt(mean(difference^2)),
  q90.absolute.percent.error = quantile(absolute.percent.error, 0.90),
  spearman.rho = cor(height.estimate, reference.other.heights,
                     method = "spearman")
), by = .(response, height)]
fwrite(height_representativeness, file.path(
  table_dir, "Table_23_Campaign1_height_representativeness.csv"
))

# Gas coupling uses observed two-hour SP medians, not model-derived
# concentrations. Slopes and R2 are descriptive diagnostics of co-variation.
gas_block_wide <- dcast(
  block_location[response %chin% c("CO2", "CH4", "NH3")],
  campaign + block.2h + sampling.point.numeric + height ~ response,
  value.var = "estimate"
)
gas_pairs <- rbindlist(list(
  gas_block_wide[, .(campaign, block.2h, sampling.point.numeric, height,
                     pair = "CH4_vs_CO2", x = CO2, y = CH4)],
  gas_block_wide[, .(campaign, block.2h, sampling.point.numeric, height,
                     pair = "NH3_vs_CO2", x = CO2, y = NH3)],
  gas_block_wide[, .(campaign, block.2h, sampling.point.numeric, height,
                     pair = "NH3_vs_CH4", x = CH4, y = NH3)]
))[is.finite(x) & x > 0 & is.finite(y) & y > 0]
gas_coupling <- gas_pairs[, {
  fit <- lm(log(y) ~ log(x))
  ci <- confint(fit)[2, ]
  .(
    observations = .N, slope = unname(coef(fit)[2]),
    slope.lower95 = ci[1], slope.upper95 = ci[2],
    R2 = summary(fit)$r.squared,
    spearman.rho = cor(x, y, method = "spearman")
  )
}, by = .(campaign, pair, height)]
fwrite(gas_coupling, file.path(
  table_dir, "Table_24_gas_coupling_by_campaign_height.csv"
))


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
plot_limits <- long[analysis.eligible == TRUE, .(
  display.low = quantile(value, 0.005),
  display.high = quantile(value, 0.995)
), by = response]
fwrite(plot_limits, file.path(table_dir, "Table_06_plot_display_limits.csv"))

long_plot <- plot_limits[long[analysis.eligible == TRUE], on = "response"]
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

summary_points <- long[analysis.eligible == TRUE, .(
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

# Main observed-data figure: campaigns are stacked as rows and each campaign
# uses the agreed 3 x 2 response arrangement. These are empirical boxplots, not
# model-derived estimates.
make_campaign_boxplot <- function(campaign_name) {
  ggplot(long_plot[campaign == campaign_name],
    aes(factor(sampling.point.numeric), value, fill = height)) +
    geom_boxplot(outlier.shape = NA, width = 0.72, linewidth = 0.22) +
    facet_wrap(~ response, ncol = 2, scales = "free_y",
               labeller = as_labeller(response_labels)) +
    scale_fill_manual(values = height_colours, name = "Height") +
    labs(x = "Sampling point", y = "Observed concentration or ratio",
         title = campaign_name) + common_theme +
    theme(axis.text.x = element_text(size = 5), legend.position = "bottom")
}
p_observed_boxes <- make_campaign_boxplot("Campaign 1") /
  make_campaign_boxplot("Campaign 2") +
  plot_annotation(
    title = "Observed gas concentrations and ratios at each sampling point",
    subtitle = paste(
      "Campaign 2 contains the 32 statistically eligible CRDS points;",
      "FTIR-connected SP49 and SP51 are excluded"
    )
  )
save_plot("Fig_11_observed_SP_boxplots_campaign1_2_main",
          p_observed_boxes, width = 15, height = 18)

pair_labels <- c(
  CH4_vs_CO2 = "CH4 versus CO2",
  NH3_vs_CO2 = "NH3 versus CO2",
  NH3_vs_CH4 = "NH3 versus CH4"
)
set.seed(20260812)
gas_pairs_plot <- gas_pairs[, if (.N > 12000L) .SD[sample(.N, 12000L)] else .SD,
                            by = .(campaign, pair)]
p_gas_coupling <- ggplot(gas_pairs_plot,
  aes(x, y, colour = height)) +
  geom_point(alpha = 0.08, size = 0.35) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.7) +
  facet_grid(pair ~ campaign, scales = "free", labeller = labeller(
    pair = as_labeller(pair_labels), campaign = label_value
  )) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = height_colours, name = "Height") +
  labs(
    x = "Observed two-hour median of predictor gas (log scale)",
    y = "Observed two-hour median of response gas (log scale)",
    title = "Gas coupling by campaign and sampling height",
    subtitle = "Lines describe co-variation and are not proof of complete physical mixing"
  ) + common_theme
save_plot("Fig_12_gas_coupling_observed_two_hour_medians",
          p_gas_coupling, width = 12, height = 10)

p_height_rep <- ggplot(c1_height_loo,
  aes(height, absolute.percent.error, fill = height)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.3) +
  facet_wrap(~ response, ncol = 2, scales = "free_y",
             labeller = as_labeller(response_labels)) +
  scale_fill_manual(values = height_colours, guide = "none") +
  coord_cartesian(ylim = c(0, quantile(
    c1_height_loo$absolute.percent.error, 0.98, na.rm = TRUE
  ))) +
  labs(x = "Height", y = "Absolute difference from other-height median (%)",
       title = "A  Representative height") + common_theme

p_sp_rep <- ggplot(sp_composite,
  aes(sampling.point.numeric, composite.score, colour = height)) +
  geom_line(colour = "grey75", linewidth = 0.35) +
  geom_point(size = 1.7) +
  geom_point(data = sp_composite[overall.rank <= 5L |
                                  sampling.point.numeric == 25L],
             shape = 21, fill = "white", stroke = 0.8, size = 3) +
  geom_text(data = sp_composite[overall.rank <= 5L |
                                 sampling.point.numeric == 25L],
            aes(label = sampling.point.numeric), nudge_y = 0.02,
            size = 3, show.legend = FALSE) +
  scale_colour_manual(values = height_colours, name = "Height") +
  scale_x_continuous(breaks = 1:51) +
  labs(x = "Sampling point", y = "Composite representativeness score",
       title = "B  Representative sampling point",
       subtitle = "Lower scores indicate closer leave-one-out network agreement") +
  common_theme + theme(axis.text.x = element_text(size = 5))

p_representative_main <- p_height_rep / p_sp_rep +
  plot_annotation(title = "Campaign 1 representative height and sampling point")
save_plot("Fig_13_representative_height_and_SP",
          p_representative_main, width = 12, height = 13)

configuration_order <- c(
  "all_51", "top_bottom_34", "campaign2_like_32",
  "top_17", "middle_17", "bottom_17"
)
configuration_plot <- melt(
  configuration_result[primary.valid == TRUE],
  id.vars = c("block.2h", "configuration"),
  measure.vars = c("Q_vent_m3_h_LU", "e_CH4_ghLU", "e_NH3_ghLU"),
  variable.name = "outcome", value.name = "value"
)[is.finite(value)]
configuration_plot[, configuration := factor(configuration,
                                               levels = configuration_order)]
outcome_labels <- c(
  Q_vent_m3_h_LU = "Ventilation rate Q (m3 h-1 LU-1)",
  e_CH4_ghLU = "CH4 emission (g h-1 LU-1)",
  e_NH3_ghLU = "NH3 emission (g h-1 LU-1)"
)
p_configuration_emission <- ggplot(configuration_plot,
  aes(configuration, value, fill = configuration)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.3) +
  facet_wrap(~ outcome, ncol = 1, scales = "free_y",
             labeller = as_labeller(outcome_labels)) +
  scale_fill_viridis_d(option = "C", guide = "none") +
  labs(
    x = "Sampling configuration", y = NULL,
    title = "Campaign 1 ventilation and emissions by sampling configuration",
    subtitle = "Primary estimates require indoor-outdoor delta CO2 >= 200 ppm"
  ) + common_theme +
  theme(axis.text.x = element_text(angle = 0, size = 8))
save_plot("Fig_14_Campaign1_Q_emissions_sampling_configuration",
          p_configuration_emission, width = 10, height = 10)

density_plot <- copy(density_summary)
density_plot[, outcome := factor(
  outcome, levels = c("Q.abs.RE", "CH4.abs.RE", "NH3.abs.RE"),
  labels = c("Ventilation rate Q", "CH4 emission", "NH3 emission")
)]
p_density <- ggplot(density_plot,
  aes(factor(sampling.points.in.estimate), median.absolute.RE,
      fill = outcome)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.3) +
  facet_wrap(~ outcome, ncol = 1, scales = "free_y") +
  scale_fill_viridis_d(option = "D", guide = "none") +
  labs(
    x = "Number of sampling points (balanced three-height network)",
    y = "Median absolute relative error versus full 51-SP network (%)",
    title = "Effect of sampling density on ventilation and emission estimates",
    subtitle = "One hundred fixed random networks per density; delta CO2 >= 200 ppm"
  ) + common_theme
save_plot("Fig_15_Campaign1_sampling_density_Q_emission_error",
          p_density, width = 10, height = 9)

p_ratio_outcome <- ggplot(ratio_outcome_correlations,
  aes(outcome, ratio, fill = spearman.rho)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", spearman.rho)), size = 4) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, limits = c(-1, 1), name = "Spearman rho") +
  scale_x_discrete(labels = c(
    Q_vent_m3_h_LU = "Q", e_CH4_ghLU = "CH4 emission",
    e_NH3_ghLU = "NH3 emission"
  )) +
  scale_y_discrete(labels = c(
    CH4_CO2 = "CH4/CO2", NH3_CO2 = "NH3/CO2", NH3_CH4 = "NH3/CH4"
  )) +
  labs(x = NULL, y = NULL,
       title = "Relationship of gas ratios with ventilation and emissions",
       subtitle = "Campaign 1 full 51-SP network; algebraic coupling requires cautious interpretation") +
  common_theme
save_plot("Fig_16_Campaign1_ratio_Q_emission_correlations",
          p_ratio_outcome, width = 8, height = 5.5)


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
