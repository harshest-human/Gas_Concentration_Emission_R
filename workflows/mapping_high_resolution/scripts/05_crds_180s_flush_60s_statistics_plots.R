###############################################################################
# Statistical sensitivity analysis for the 180 s flush + 60 s CRDS average
#
# Inputs
#   clean_data/crds_180s_flush_60s_sensitivity_v01/all_campaigns_combined.csv
#
# Outputs
#   descriptive statistics, leave-one-out relative errors, coefficients of
#   variation, crossed mixed-model variance components, height likelihood-
#   ratio tests, and manuscript-ready PNG figures.
#
# Rules
#   - Internal numeric sampling points are used for spatial statistics.
#   - Missing, zero and negative values are excluded response-wise.
#   - No IQR or other statistical outlier removal is applied.
#   - Ratios are dimensionless and are not multiplied by 100.
#   - Two-hour clock-aligned blocks synchronise the sequential SP network.
#   - Shannon entropy is deliberately not calculated in this workflow.
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(lubridate)
  library(ggplot2)
  library(patchwork)
  library(lme4)
})

workflow <- paste0(
  "D:/Data_Analysis_R/Gas_Concentration_Emission_R/",
  "workflows/mapping_high_resolution"
)
input_root <- file.path(
  workflow, "clean_data", "crds_180s_flush_60s_sensitivity_v01"
)
statistics_dir <- file.path(input_root, "statistics")
plot_dir <- file.path(
  workflow, "plots", "crds_180s_flush_60s_sensitivity_v01"
)
dir.create(statistics_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

tz_local <- "Europe/Berlin"
responses <- c(
  "CO2", "NH3_CO2", "CH4", "CH4_CO2", "NH3", "NH3_CH4"
)
response_labels <- c(
  CO2 = "CO[2]~'(ppm)'",
  NH3_CO2 = "NH[3]/CO[2]",
  CH4 = "CH[4]~'(ppm)'",
  CH4_CO2 = "CH[4]/CO[2]",
  NH3 = "NH[3]~'(ppm)'",
  NH3_CH4 = "NH[3]/CH[4]"
)
vertical_colours <- c(
  top = "orange", middle = "green3", bottom = "steelblue1"
)

###############################################################################
##### Read and derive responses
###############################################################################

gas <- fread(
  file.path(input_root, "all_campaigns_combined.csv"),
  colClasses = c("DATE.TIME" = "character")
)
gas[, DATE.TIME := as.POSIXct(
  DATE.TIME, format = "%Y-%m-%d %H:%M:%S", tz = tz_local
)]
stopifnot(!anyNA(gas$DATE.TIME))

gas[, sampling.point.numeric := suppressWarnings(as.integer(sampling.point))]
gas[, internal.point := !is.na(sampling.point.numeric)]
gas[, horizontal.position := (sampling.point.numeric + 2L) %/% 3L]
gas[, vertical.group := fcase(
  sampling.point.numeric %% 3L == 1L, "top",
  sampling.point.numeric %% 3L == 2L, "middle",
  sampling.point.numeric %% 3L == 0L, "bottom",
  default = NA_character_
)]
gas[, vertical.group := factor(
  vertical.group, levels = c("top", "middle", "bottom")
)]
gas[, block.2h := floor_date(DATE.TIME, "2 hours")]

# Ratios are formed only from positive constituent concentrations.
gas[, NH3_CO2 := fifelse(NH3 > 0 & CO2 > 0, NH3 / CO2, NA_real_)]
gas[, CH4_CO2 := fifelse(CH4 > 0 & CO2 > 0, CH4 / CO2, NA_real_)]
gas[, NH3_CH4 := fifelse(NH3 > 0 & CH4 > 0, NH3 / CH4, NA_real_)]

long <- melt(
  gas[internal.point == TRUE],
  id.vars = c(
    "DATE.TIME", "campaign", "analyser", "sampling.point",
    "sampling.point.numeric", "horizontal.position", "vertical.group",
    "block.2h"
  ),
  measure.vars = responses,
  variable.name = "response", value.name = "value"
)
long[, response := factor(response, levels = responses)]
long <- long[is.finite(value) & value > 0]

fwrite(
  long[, .(
    campaign, DATE.TIME, block.2h, analyser, sampling.point,
    sampling.point.numeric, horizontal.position, vertical.group,
    response, value
  )],
  file.path(statistics_dir, "analysis_long_positive_internal.csv"),
  dateTimeAs = "write.csv"
)

###############################################################################
##### Descriptive statistics
###############################################################################

sp_summary <- long[, .(
  n = .N,
  mean = mean(value),
  median = median(value),
  minimum = min(value),
  maximum = max(value),
  sd = sd(value),
  cv_percent = 100 * sd(value) / abs(mean(value)),
  mad = mad(value, constant = 1)
), by = .(
  campaign, response, sampling.point, sampling.point.numeric,
  horizontal.position, vertical.group
)]

network_summary <- long[, .(
  n = .N,
  sampling.points = uniqueN(sampling.point),
  mean = mean(value),
  median = median(value),
  minimum = min(value),
  maximum = max(value),
  sd = sd(value),
  cv_percent = 100 * sd(value) / abs(mean(value)),
  mad = mad(value, constant = 1)
), by = .(campaign, response)]

vertical_summary <- long[, .(
  n = .N,
  sampling.points = uniqueN(sampling.point),
  mean = mean(value),
  median = median(value),
  minimum = min(value),
  maximum = max(value),
  sd = sd(value),
  cv_percent = 100 * sd(value) / abs(mean(value)),
  mad = mad(value, constant = 1)
), by = .(campaign, response, vertical.group)]

fwrite(sp_summary, file.path(statistics_dir, "sampling_point_descriptive.csv"))
fwrite(network_summary, file.path(statistics_dir, "network_descriptive.csv"))
fwrite(vertical_summary, file.path(statistics_dir, "vertical_group_descriptive.csv"))

###############################################################################
##### Block-wise leave-one-out relative error
###############################################################################

# Median aggregation within SP and two-hour block prevents an analyser or SP
# with repeated observations from receiving greater weight.
block_sp <- long[, .(
  value = median(value),
  observations = .N
), by = .(
  campaign, response, block.2h, sampling.point, sampling.point.numeric,
  horizontal.position, vertical.group
)]

block_sp[, network.points := .N, by = .(campaign, response, block.2h)]
block_sp[, loo.reference := vapply(
  seq_len(.N),
  function(i) {
    other <- value[-i]
    if (length(other) < 2L) return(NA_real_)
    median(other, na.rm = TRUE)
  },
  numeric(1)
), by = .(campaign, response, block.2h)]
block_sp[, loo.re_percent := 100 * (value - loo.reference) / loo.reference]
block_sp[, loo.absolute.re_percent := abs(loo.re_percent)]
block_sp <- block_sp[
  network.points >= 3L & is.finite(loo.absolute.re_percent)
]

loo_summary <- block_sp[, .(
  blocks = .N,
  mean_signed_re_percent = mean(loo.re_percent),
  median_signed_re_percent = median(loo.re_percent),
  mean_absolute_re_percent = mean(loo.absolute.re_percent),
  median_absolute_re_percent = median(loo.absolute.re_percent),
  p95_absolute_re_percent = as.numeric(
    quantile(loo.absolute.re_percent, 0.95, names = FALSE)
  )
), by = .(
  campaign, response, sampling.point, sampling.point.numeric,
  horizontal.position, vertical.group
)]

fwrite(block_sp, file.path(statistics_dir, "blockwise_leave_one_out_re.csv"),
       dateTimeAs = "write.csv")
fwrite(loo_summary, file.path(statistics_dir, "sampling_point_leave_one_out_re.csv"))

sp_complete <- merge(
  sp_summary, loo_summary,
  by = c(
    "campaign", "response", "sampling.point", "sampling.point.numeric",
    "horizontal.position", "vertical.group"
  ),
  all.x = TRUE
)
fwrite(sp_complete, file.path(
  statistics_dir, "sampling_point_all_statistics.csv"
))

###############################################################################
##### Crossed mixed models: spatial versus temporal variance
###############################################################################

fit_variance_model <- function(z) {
  z <- copy(z[is.finite(value) & value > 0])
  if (nrow(z) < 100L || uniqueN(z$sampling.point) < 3L ||
      uniqueN(z$block.2h) < 3L) return(NULL)
  z[, `:=`(
    sampling.point.factor = factor(sampling.point),
    block.factor = factor(block.2h),
    log.value = log(value)
  )]
  tryCatch(
    lmer(
      log.value ~ 1 + (1 | sampling.point.factor) + (1 | block.factor),
      data = z, REML = TRUE,
      control = lmerControl(
        optimizer = "bobyqa", calc.derivs = FALSE,
        check.conv.singular = "ignore"
      )
    ),
    error = function(e) NULL
  )
}

variance_results <- list()
model_diagnostics <- list()
groups <- unique(block_sp[, .(campaign, response)])
for (i in seq_len(nrow(groups))) {
  cid <- groups$campaign[i]
  response_id <- as.character(groups$response[i])
  z <- block_sp[campaign == cid & response == response_id]
  model <- fit_variance_model(z)
  if (is.null(model)) next

  vc <- as.data.table(VarCorr(model))[, .(grp, vcov)]
  spatial <- vc[grp == "sampling.point.factor", sum(vcov)]
  temporal <- vc[grp == "block.factor", sum(vcov)]
  residual <- vc[grp == "Residual", sum(vcov)]
  total <- spatial + temporal + residual
  variance_results[[length(variance_results) + 1L]] <- data.table(
    campaign = cid,
    response = response_id,
    spatial_variance = spatial,
    temporal_variance = temporal,
    residual_variance = residual,
    spatial_percent = 100 * spatial / total,
    temporal_percent = 100 * temporal / total,
    residual_percent = 100 * residual / total,
    observations = nobs(model),
    sampling.points = uniqueN(z$sampling.point),
    time.blocks = uniqueN(z$block.2h),
    singular = isSingular(model, tol = 1e-5)
  )
  model_diagnostics[[length(model_diagnostics) + 1L]] <- data.table(
    campaign = cid, response = response_id,
    AIC = AIC(model), BIC = BIC(model), logLik = as.numeric(logLik(model))
  )
}
variance_partition <- rbindlist(variance_results, fill = TRUE)
model_diagnostics <- rbindlist(model_diagnostics, fill = TRUE)
fwrite(variance_partition, file.path(
  statistics_dir, "mixed_model_spatial_temporal_variance.csv"
))
fwrite(model_diagnostics, file.path(
  statistics_dir, "mixed_model_diagnostics.csv"
))

###############################################################################
##### Vertical-group likelihood-ratio tests
###############################################################################

height_tests <- list()
for (i in seq_len(nrow(groups))) {
  cid <- groups$campaign[i]
  response_id <- as.character(groups$response[i])
  z <- copy(block_sp[campaign == cid & response == response_id])
  z <- z[is.finite(value) & value > 0 & !is.na(vertical.group)]
  if (nrow(z) < 100L || uniqueN(z$vertical.group) < 2L ||
      uniqueN(z$horizontal.position) < 3L) next
  z[, `:=`(
    log.value = log(value),
    vertical.group = droplevels(factor(vertical.group)),
    horizontal.factor = factor(horizontal.position),
    block.factor = factor(block.2h)
  )]
  m0 <- tryCatch(lmer(
    log.value ~ 1 + (1 | horizontal.factor) + (1 | block.factor),
    data = z, REML = FALSE,
    control = lmerControl(optimizer = "bobyqa", calc.derivs = FALSE)
  ), error = function(e) NULL)
  m1 <- tryCatch(lmer(
    log.value ~ vertical.group + (1 | horizontal.factor) + (1 | block.factor),
    data = z, REML = FALSE,
    control = lmerControl(optimizer = "bobyqa", calc.derivs = FALSE)
  ), error = function(e) NULL)
  if (is.null(m0) || is.null(m1)) next
  comparison <- anova(m0, m1)
  height_tests[[length(height_tests) + 1L]] <- data.table(
    campaign = cid, response = response_id,
    likelihood_ratio_chisq = comparison$Chisq[2],
    df = comparison$Df[2],
    p_value = comparison$`Pr(>Chisq)`[2],
    AIC_without_height = AIC(m0),
    AIC_with_height = AIC(m1),
    observations = nrow(z)
  )
}
height_tests <- rbindlist(height_tests, fill = TRUE)
fwrite(height_tests, file.path(
  statistics_dir, "vertical_group_likelihood_ratio_tests.csv"
))

###############################################################################
##### Comparison with the former 60 s flush + 180 s average
###############################################################################

old_path <- file.path(
  workflow, "clean_data", "prepared_campaigns_v01",
  "all_campaigns_combined.csv"
)
if (file.exists(old_path)) {
  old <- fread(old_path, colClasses = c("DATE.TIME" = "character"))
  old <- old[grepl("^CRDS", analyser)]
  old[, sampling.point.numeric := suppressWarnings(as.integer(sampling.point))]
  old <- old[!is.na(sampling.point.numeric)]
  old[, `:=`(
    NH3_CO2 = fifelse(NH3 > 0 & CO2 > 0, NH3 / CO2, NA_real_),
    CH4_CO2 = fifelse(CH4 > 0 & CO2 > 0, CH4 / CO2, NA_real_),
    NH3_CH4 = fifelse(NH3 > 0 & CH4 > 0, NH3 / CH4, NA_real_)
  )]
  old_long <- melt(
    old, id.vars = c("campaign", "analyser"), measure.vars = responses,
    variable.name = "response", value.name = "value"
  )[is.finite(value) & value > 0]
  old_summary <- old_long[, .(
    old_mean = mean(value), old_median = median(value), old_n = .N
  ), by = .(campaign, analyser, response)]
  new_summary <- long[grepl("^CRDS", analyser), .(
    new_mean = mean(value), new_median = median(value), new_n = .N
  ), by = .(campaign, analyser, response)]
  flush_comparison <- merge(
    old_summary, new_summary,
    by = c("campaign", "analyser", "response"), all = TRUE
  )
  flush_comparison[, `:=`(
    mean_change_percent = 100 * (new_mean - old_mean) / old_mean,
    median_change_percent = 100 * (new_median - old_median) / old_median
  )]
  fwrite(flush_comparison, file.path(
    statistics_dir, "old_60s_vs_new_180s_flush_comparison.csv"
  ))
}

###############################################################################
##### Compact summaries used by the standalone sensitivity report
###############################################################################

report_network_gases <- network_summary[
  response %in% c("CO2", "CH4", "NH3")
][order(campaign, match(response, c("CO2", "CH4", "NH3")))]
fwrite(report_network_gases, file.path(
  statistics_dir, "report_network_gas_summary.csv"
))

report_variance_gases <- variance_partition[
  response %in% c("CO2", "CH4", "NH3")
][order(campaign, match(response, c("CO2", "CH4", "NH3")))]
fwrite(report_variance_gases, file.path(
  statistics_dir, "report_mixed_model_gas_variance.csv"
))

report_sp_metrics <- sp_complete[, .(
  sampling_points = uniqueN(sampling.point),
  median_sp_cv_percent = median(cv_percent, na.rm = TRUE),
  maximum_sp_cv_percent = max(cv_percent, na.rm = TRUE),
  median_sp_loo_absolute_re_percent = median(
    median_absolute_re_percent, na.rm = TRUE
  ),
  maximum_sp_loo_absolute_re_percent = max(
    median_absolute_re_percent, na.rm = TRUE
  )
), by = .(campaign, response)]
fwrite(report_sp_metrics, file.path(
  statistics_dir, "report_sampling_point_metric_summary.csv"
))

if (exists("flush_comparison")) {
  report_flush_gases <- flush_comparison[
    response %in% c("CO2", "CH4", "NH3")
  ][, .(
    analysers = uniqueN(analyser),
    median_mean_change_percent = median(mean_change_percent, na.rm = TRUE),
    minimum_mean_change_percent = min(mean_change_percent, na.rm = TRUE),
    maximum_mean_change_percent = max(mean_change_percent, na.rm = TRUE),
    median_median_change_percent = median(
      median_change_percent, na.rm = TRUE
    )
  ), by = .(campaign, response)]
  fwrite(report_flush_gases, file.path(
    statistics_dir, "report_old_vs_new_flush_gas_summary.csv"
  ))
}

###############################################################################
##### Figures
###############################################################################

theme_manuscript <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey94"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

sp_summary[, response := factor(response, levels = responses)]
for (cid in 1:4) {
  z <- sp_summary[campaign == cid]
  p <- ggplot(
    z,
    aes(
      x = sampling.point.numeric, y = mean, colour = vertical.group,
      group = vertical.group
    )
  ) +
    geom_errorbar(
      aes(ymin = pmax(0, mean - sd), ymax = mean + sd),
      width = 0.25, linewidth = 0.35
    ) +
    geom_point(size = 1.6) +
    facet_wrap(
      ~response, ncol = 2, scales = "free_y",
      labeller = as_labeller(response_labels, default = label_parsed)
    ) +
    scale_colour_manual(values = vertical_colours, drop = FALSE) +
    scale_x_continuous(
      breaks = sort(unique(z$sampling.point.numeric)),
      expand = expansion(mult = c(0.015, 0.015))
    ) +
    labs(
      x = "Sampling point", y = "Mean +/- SD",
      colour = "Vertical group",
      title = paste0("Campaign ", cid, ": mean +/- SD by sampling point"),
      subtitle = "CRDS values use a 180 s flush and the following 60 s average"
    ) +
    theme_manuscript +
    theme(axis.text.x = element_text(size = 6))
  ggsave(
    file.path(plot_dir, paste0("Fig_mean_SD_campaign", cid, ".png")),
    p, width = 13, height = 12, dpi = 300, bg = "white"
  )
}

# CV and leave-one-out RE heatmaps use a common 0--100% display range. Values
# above 100% remain available in the CSV and are capped only for colour display.
cv_plot <- copy(sp_summary)
cv_plot[, display_value := pmin(cv_percent, 100)]
p_cv <- ggplot(
  cv_plot,
  aes(x = sampling.point.numeric, y = response, fill = display_value)
) +
  geom_tile(colour = "white", linewidth = 0.15) +
  facet_wrap(~campaign, ncol = 1, scales = "free_x") +
  scale_fill_gradientn(
    colours = c("#2166ac", "#67a9cf", "#1a9850", "#fee08b", "#f46d43", "#b2182b"),
    limits = c(0, 100), name = "CV (%)\n(capped at 100)"
  ) +
  scale_x_continuous(breaks = 1:51) +
  scale_y_discrete(labels = c(
    CO2 = "CO2", NH3_CO2 = "NH3/CO2", CH4 = "CH4",
    CH4_CO2 = "CH4/CO2", NH3 = "NH3", NH3_CH4 = "NH3/CH4"
  )) +
  labs(x = "Sampling point", y = NULL, title = "Sampling-point coefficient of variation") +
  theme_manuscript +
  theme(axis.text.x = element_text(size = 6), panel.grid = element_blank())
ggsave(file.path(plot_dir, "Fig_CV_heatmaps_campaign1_4.png"), p_cv,
       width = 13, height = 12, dpi = 300, bg = "white")

loo_plot <- copy(loo_summary)
loo_plot[, response := factor(response, levels = responses)]
loo_plot[, display_value := pmin(median_absolute_re_percent, 100)]
p_loo <- ggplot(
  loo_plot,
  aes(x = sampling.point.numeric, y = response, fill = display_value)
) +
  geom_tile(colour = "white", linewidth = 0.15) +
  facet_wrap(~campaign, ncol = 1, scales = "free_x") +
  scale_fill_gradientn(
    colours = c("#2166ac", "#67a9cf", "#1a9850", "#fee08b", "#f46d43", "#b2182b"),
    limits = c(0, 100), name = "Median absolute\nLOO RE (%)"
  ) +
  scale_x_continuous(breaks = 1:51) +
  scale_y_discrete(labels = c(
    CO2 = "CO2", NH3_CO2 = "NH3/CO2", CH4 = "CH4",
    CH4_CO2 = "CH4/CO2", NH3 = "NH3", NH3_CH4 = "NH3/CH4"
  )) +
  labs(
    x = "Sampling point", y = NULL,
    title = "Block-wise leave-one-out relative error",
    subtitle = "Reference: median of all other internal SPs in each two-hour block"
  ) +
  theme_manuscript +
  theme(axis.text.x = element_text(size = 6), panel.grid = element_blank())
ggsave(file.path(plot_dir, "Fig_LOO_RE_heatmaps_campaign1_4.png"), p_loo,
       width = 13, height = 12, dpi = 300, bg = "white")

variance_long <- melt(
  variance_partition,
  id.vars = c("campaign", "response"),
  measure.vars = c("spatial_percent", "temporal_percent", "residual_percent"),
  variable.name = "component", value.name = "percent"
)
variance_long[, component := factor(
  component,
  levels = c("residual_percent", "spatial_percent", "temporal_percent"),
  labels = c("Residual", "Spatial SP", "Temporal two-hour block")
)]
variance_long[, response := factor(response, levels = responses)]
p_variance <- ggplot(
  variance_long,
  aes(x = response, y = percent, fill = component)
) +
  geom_col(width = 0.75) +
  facet_wrap(~campaign, ncol = 2) +
  scale_fill_manual(values = c(
    "Residual" = "grey70", "Spatial SP" = "#d95f02",
    "Temporal two-hour block" = "#1b9e77"
  )) +
  scale_x_discrete(labels = c(
    CO2 = "CO2", NH3_CO2 = "NH3/CO2", CH4 = "CH4",
    CH4_CO2 = "CH4/CO2", NH3 = "NH3", NH3_CH4 = "NH3/CH4"
  )) +
  labs(
    x = NULL, y = "Variance component (%)", fill = NULL,
    title = "Crossed mixed-model variance partition",
    subtitle = "Models fitted to log-transformed two-hour sampling-point medians"
  ) +
  theme_manuscript +
  theme(axis.text.x = element_text(angle = 0, size = 8))
ggsave(file.path(plot_dir, "Fig_mixed_model_variance_partition.png"),
       p_variance, width = 11, height = 7.5, dpi = 300, bg = "white")

message("Statistics written to: ", statistics_dir)
message("Figures written to: ", plot_dir)
