##### Manuscript 3 v12: unified analysis and figure pipeline ####

# This script rebuilds the manuscript analysis using internal locations only.
# Location "s"/52 is retained in immutable source files but is excluded here.

library(data.table)
library(ggplot2)
library(lubridate)
library(lme4)
library(lmerTest)
library(emmeans)
library(patchwork)

set.seed(20260807)

workflow <- "D:/Data_Analysis_R/Gas_Concentration_Emission_R/workflows/mapping_high_resolution"
out <- file.path(workflow, "clean_data", "manuscript3_v12", "statistical_tables")
fig <- file.path(workflow, "plots", "manuscript3_v12", "main")
man_fig <- file.path(workflow,
  "Manuscript_3_Mapping_high_resolution_concentration_Latex_draft", "figures")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(fig, recursive = TRUE, showWarnings = FALSE)
dir.create(man_fig, recursive = TRUE, showWarnings = FALSE)

height_colours <- c(top = "orange", middle = "green3", bottom = "steelblue1")
response_order <- c("CO2", "NH3_CO2", "CH4", "CH4_CO2", "NH3", "NH3_CH4")
response_labels <- c(
  CO2 = expression(CO[2]), CH4 = expression(CH[4]), NH3 = expression(NH[3]),
  NH3_CO2 = expression(NH[3]/CO[2]), CH4_CO2 = expression(CH[4]/CO[2]),
  NH3_CH4 = expression(NH[3]/CH[4])
)

##### Read and harmonise campaign data ####

c1_summer_files <- list.files(
  file.path(workflow, "clean_data", "1_campaign", "recovered_crds_june_august_v03"),
  pattern = "four_analyser_CRDS_scale_corrected_.*_v03[.]csv$",
  recursive = TRUE, full.names = TRUE
)
c1_autumn_file <- file.path(
  workflow, "clean_data", "1_campaign",
  "20241001.000232_20241024.143858_campaign1_raw_and_CRDS_scale_corrected.csv"
)
c2_file <- file.path(
  workflow, "clean_data", "2_campaign",
  "20241116.000647_20241231.235710_CRDS8_CRDS9_campaign2_combined.csv"
)
stopifnot(length(c1_summer_files) == 3L, file.exists(c1_autumn_file), file.exists(c2_file))

# Keep only the intended 240 s CRDS switching steps.  The recovered step-average
# files retain the original dwell duration, which allows shorter QA/configuration
# records to be removed before any summary statistic or model is calculated.
c1_step_files <- list.files(
  file.path(workflow, "clean_data", "1_campaign", "recovered_crds_june_august_v03"),
  pattern = "campaign1_CRDS8_CRDS9_recovered_step_average_.*_v03[.]csv$",
  recursive = TRUE, full.names = TRUE
)
stopifnot(length(c1_step_files) == 3L)
short_step_keys <- rbindlist(lapply(c1_step_files, function(file) {
  z <- fread(file, colClasses = c(DATE.TIME = "character", location = "character"))
  z[, DATE.TIME := ymd_hms(DATE.TIME, tz = "Europe/Berlin", quiet = TRUE)]
  z[dwell_seconds < 240, .(DATE.TIME, analyser, location)]
}), use.names = TRUE, fill = TRUE)
short_step_keys <- unique(short_step_keys)

read_c1 <- function(file, period) {
  z <- fread(file, colClasses = c(DATE.TIME = "character", location = "character"))
  z[, DATE.TIME := ymd_hms(DATE.TIME, tz = "Europe/Berlin", quiet = TRUE)]
  z[, `:=`(
    period = period, campaign = "Campaign 1",
    CO2 = CO2_corr, CH4 = CH4_corr, NH3 = NH3_corr
  )]
  z[, .(DATE.TIME, campaign, period, analyser, location, CO2, CH4, NH3)]
}

read_c2 <- function(file) {
  z <- fread(file, colClasses = c(DATE.TIME = "character", location = "character"))
  z[, DATE.TIME := ymd_hms(DATE.TIME, tz = "Europe/Berlin", quiet = TRUE)]
  z[, period := fifelse(month(DATE.TIME) == 11L,
    "C2 late autumn (Nov)", "C2 winter (Dec)")]
  z[, campaign := "Campaign 2"]
  z[, .(DATE.TIME, campaign, period, analyser, location, CO2, CH4, NH3)]
}

x <- rbindlist(list(
  rbindlist(lapply(c1_summer_files, read_c1, period = "C1 summer (Jun-Aug)")),
  read_c1(c1_autumn_file, "C1 autumn (Oct)"),
  read_c2(c2_file)
), use.names = TRUE, fill = TRUE)

if (nrow(short_step_keys)) {
  x <- x[!short_step_keys, on = .(DATE.TIME, analyser, location)]
}

x[, location_n := suppressWarnings(as.integer(location))]
x <- x[!is.na(DATE.TIME) & !is.na(location_n) & location_n %between% c(1L, 51L)]
x[CO2 <= 0, CO2 := NA_real_]
x[CH4 <= 0, CH4 := NA_real_]
x[NH3 <= 0, NH3 := NA_real_]
x[, `:=`(
  horizontal_position = ceiling(location_n / 3),
  height = factor(c("top", "middle", "bottom")[(location_n - 1L) %% 3L + 1L],
                  levels = c("bottom", "middle", "top")),
  height_m = c(NA_real_, 3.6, 2.6)[1L]
)]
x[height == "middle", height_m := 3.6]
x[height == "bottom", height_m := 2.6]
# Top points were defined relative to the roof: 0.60 m below the local roof.
x[, `:=`(
  season = factor(fcase(
    period == "C1 summer (Jun-Aug)", "summer",
    period %chin% c("C1 autumn (Oct)", "C2 late autumn (Nov)"), "autumn",
    default = "winter"), levels = c("summer", "autumn", "winter")),
  date = as.Date(DATE.TIME, tz = "Europe/Berlin"),
  hour_utc = floor_date(with_tz(DATE.TIME, "UTC"), "hour"),
  CH4_CO2 = CH4 / CO2,
  NH3_CO2 = NH3 / CO2,
  NH3_CH4 = NH3 / CH4
)]

period_levels <- c("C1 summer (Jun-Aug)", "C1 autumn (Oct)",
                   "C2 late autumn (Nov)", "C2 winter (Dec)")
x[, period := factor(period, levels = period_levels)]

responses <- c("CO2", "CH4", "NH3", "NH3_CO2", "CH4_CO2", "NH3_CH4")
long <- melt(x,
  id.vars = c("DATE.TIME", "hour_utc", "date", "campaign", "period", "season",
              "analyser", "location_n", "horizontal_position", "height"),
  measure.vars = responses, variable.name = "response", value.name = "value"
)[is.finite(value) & value > 0]
long[, response := factor(response, levels = response_order)]

fwrite(x, file.path(out, "manuscript3_internal_only_analysis_dataset.csv"))
fwrite(x[, .(rows = .N, locations = uniqueN(location_n), dates = uniqueN(date)),
         by = .(campaign, period, season)],
       file.path(out, "Table_01_analysis_period_audit.csv"))

##### Hourly location data and campaign-1 spatial summaries ####

hourly <- long[, .(value = mean(value, na.rm = TRUE), raw_n = .N),
  by = .(period, season, date, hour_utc, response, location_n,
         horizontal_position, height)]

c1_hourly <- hourly[period %in% c("C1 summer (Jun-Aug)", "C1 autumn (Oct)")]
fig2_summary <- c1_hourly[, .(
  n_hours = .N,
  mean = mean(value), sd = sd(value), median = median(value),
  cv_pct = 100 * sd(value) / abs(mean(value))
), by = .(response, location_n, horizontal_position, height)]
fwrite(fig2_summary, file.path(out, "Table_02_campaign1_location_mean_SD_CV.csv"))

facet_labeller <- as_labeller(response_labels, default = label_parsed)

# Confidence intervals are based on daily location means, rather than treating
# autocorrelated hourly observations as independent replicates.
fig2_ci_summary <- c1_hourly[, .(daily_mean = mean(value)),
  by = .(period, date, response, location_n, horizontal_position, height)
][, .(
  n_days = .N,
  mean = mean(daily_mean),
  sd_day = sd(daily_mean),
  se = sd(daily_mean) / sqrt(.N)
), by = .(response, location_n, horizontal_position, height)]
fig2_ci_summary[, `:=`(
  ci_low = mean - qt(0.975, pmax(n_days - 1L, 1L)) * se,
  ci_high = mean + qt(0.975, pmax(n_days - 1L, 1L)) * se
)]
fwrite(fig2_ci_summary, file.path(out, "Table_02b_campaign1_location_mean_95CI.csv"))

p_mean <- ggplot(fig2_ci_summary,
  aes(horizontal_position, mean, colour = height, group = height)) +
  geom_line(linewidth = 0.55) +
  geom_errorbar(aes(ymin = pmax(0, ci_low), ymax = ci_high),
                width = 0.14, linewidth = 0.35, alpha = 0.75) +
  geom_point(size = 1.8) +
  facet_wrap(~response, scales = "free_y", ncol = 3, labeller = facet_labeller) +
  scale_x_continuous(breaks = 1:17, limits = c(0.6, 17.4)) +
  scale_colour_manual(values = height_colours, name = "Height") +
  labs(x = "Horizontal position", y = "Campaign 1 mean (95% CI)") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

cv_plot_data <- copy(fig2_summary)
cv_plot_data[, `:=`(cv_display = pmin(cv_pct, 100), exceeds_100 = cv_pct > 100)]
p_cv <- ggplot(cv_plot_data,
  aes(horizontal_position, cv_display, fill = height, group = height)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.72) +
  geom_point(data = cv_plot_data[exceeds_100 == TRUE], aes(y = 99),
             position = position_dodge(width = 0.78), shape = 17, size = 2) +
  facet_wrap(~response, ncol = 3, labeller = facet_labeller) +
  scale_x_continuous(breaks = 1:17, limits = c(0.45, 17.55)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = c(0, 0)) +
  scale_fill_manual(values = height_colours, name = "Height") +
  labs(x = "Horizontal position", y = "Coefficient of variation (%)") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

##### Three single-height configurations relative to the dense mean ####

c1 <- x[period %in% c("C1 summer (Jun-Aug)", "C1 autumn (Oct)")]
c1[, block_2h := as.POSIXct(floor(as.numeric(DATE.TIME) / 7200) * 7200,
                            origin = "1970-01-01", tz = "Europe/Berlin")]
block <- c1[, lapply(.SD, mean, na.rm = TRUE),
  by = .(period, block_2h, location_n, horizontal_position, height),
  .SDcols = responses]
complete_blocks <- block[, .(locations = uniqueN(location_n)),
                         by = .(period, block_2h)][locations == 51L]
block <- block[complete_blocks, on = .(period, block_2h)]

block_long <- melt(block,
  id.vars = c("period", "block_2h", "location_n", "horizontal_position", "height"),
  measure.vars = responses, variable.name = "response", value.name = "value"
)[is.finite(value)]
block_long[, response := factor(response, levels = response_order)]
block_long[, dense_mean := mean(value), by = .(period, block_2h, response)]
height_cfg <- block_long[, .(estimate = mean(value), dense_mean = first(dense_mean)),
  by = .(period, block_2h, response, height)]
height_cfg[, RE_pct := 100 * (estimate - dense_mean) / dense_mean]
fwrite(height_cfg, file.path(out, "Table_03_single_height_relative_errors_blocks.csv"))
fwrite(height_cfg[, .(
  blocks = .N, mean_RE_pct = mean(RE_pct), median_RE_pct = median(RE_pct),
  mean_absolute_RE_pct = mean(abs(RE_pct)),
  p95_absolute_RE_pct = quantile(abs(RE_pct), 0.95)
), by = .(period, response, height)],
file.path(out, "Table_03b_single_height_relative_error_summary.csv"))

p_height <- ggplot(height_cfg, aes(height, RE_pct, fill = height)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
  geom_boxplot(outlier.alpha = 0.08, width = 0.68) +
  facet_wrap(~response, scales = "free_y", ncol = 3, labeller = facet_labeller) +
  scale_fill_manual(values = height_colours, guide = "none") +
  labs(x = "Single-height configuration",
       y = "Relative deviation from the 51-location mean (%)") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank())

##### Shannon entropy with square spatial tiles ####

entropy4 <- function(v) {
  v <- v[is.finite(v)]
  if (length(v) < 8L || diff(range(v)) == 0) return(NA_real_)
  bins <- cut(v, breaks = seq(min(v), max(v), length.out = 5), include.lowest = TRUE)
  p <- as.numeric(table(bins)) / length(v)
  p <- p[p > 0]
  -sum(p * log2(p)) / log2(4)
}

entropy <- block_long[, .(
  n_blocks = .N, entropy_normalised = entropy4(value)
), by = .(response, location_n, horizontal_position, height)]
entropy[, response := factor(as.character(response), levels = c(
  "CO2", "NH3_CO2",
  "CH4", "CH4_CO2",
  "NH3", "NH3_CH4"
))]
entropy[, height_y := c(bottom = 1, middle = 2, top = 3)[as.character(height)]]
fwrite(entropy, file.path(out, "Table_04_four_bin_Shannon_entropy.csv"))

p_entropy <- ggplot(entropy,
  aes(horizontal_position, height_y, fill = entropy_normalised)) +
  geom_tile(width = 0.94, height = 0.94, colour = "white", linewidth = 0.25) +
  geom_text(aes(label = location_n), size = 2.3) +
  facet_wrap(~response, ncol = 2, labeller = facet_labeller) +
  scale_x_continuous(breaks = NULL, labels = NULL,
                     limits = c(0.5, 17.5), expand = c(0, 0)) +
  scale_y_continuous(breaks = 1:3, labels = c("bottom", "middle", "top"),
                     limits = c(0.5, 3.5), expand = c(0, 0)) +
  scale_fill_viridis_c(limits = c(0, 1), name = "Normalised\nShannon entropy") +
  coord_fixed(ratio = 1) +
  labs(x = NULL, y = "Height") +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(), legend.position = "bottom")

##### Spatiotemporal precision across hour, day, week and period ####

scale_blocks <- list(
  hour = function(t) floor_date(t, "hour"),
  day = function(t) floor_date(t, "day"),
  week = function(t) floor_date(t, "week", week_start = 1)
)

precision_blocks <- rbindlist(lapply(names(scale_blocks), function(scale_name) {
  z <- copy(long)
  z[, time_block := scale_blocks[[scale_name]](DATE.TIME)]
  loc <- z[, .(value = mean(value)),
    by = .(period, response, time_block, location_n, horizontal_position, height)]
  loc[, .(
    locations = uniqueN(location_n), spatial_mean = mean(value),
    spatial_sd = sd(value), spatial_cv_pct = 100 * sd(value) / abs(mean(value))
  ), by = .(period, response, time_block)][, scale := scale_name][]
}))

period_location <- long[, .(value = mean(value)),
  by = .(period, response, location_n, horizontal_position, height)]
period_precision <- period_location[, .(
  locations = uniqueN(location_n), spatial_mean = mean(value),
  spatial_sd = sd(value), spatial_cv_pct = 100 * sd(value) / abs(mean(value)),
  time_block = as.POSIXct(NA)
), by = .(period, response)][, scale := "period"]
precision_blocks <- rbind(precision_blocks, period_precision, fill = TRUE)
precision_blocks[, scale := factor(scale, levels = c("hour", "day", "week", "period"))]
fwrite(precision_blocks, file.path(out, "Table_05_spatial_precision_by_temporal_scale.csv"))
precision_summary <- precision_blocks[, .(
  blocks = .N, median_spatial_CV_pct = median(spatial_cv_pct, na.rm = TRUE),
  q25 = quantile(spatial_cv_pct, 0.25, na.rm = TRUE),
  q75 = quantile(spatial_cv_pct, 0.75, na.rm = TRUE)
), by = .(period, response, scale)]
fwrite(precision_summary, file.path(out, "Table_05b_spatial_precision_summary.csv"))

# Repeated location-level CV supports a formal response-by-scale test.
location_precision <- rbindlist(lapply(names(scale_blocks), function(scale_name) {
  z <- copy(long)
  z[, time_block := scale_blocks[[scale_name]](DATE.TIME)]
  z <- z[, .(value = mean(value)),
    by = .(period, response, time_block, location_n)]
  z[, .(n_blocks = .N, cv_pct = 100 * sd(value) / abs(mean(value))),
    by = .(period, response, location_n)][, scale := scale_name][]
}))
location_precision <- location_precision[is.finite(cv_pct) & cv_pct > 0]
location_precision[, `:=`(
  scale = factor(scale, levels = c("hour", "day", "week")),
  location_period = interaction(period, location_n, drop = TRUE)
)]
precision_model <- lmer(log(cv_pct) ~ response * scale + (1 | location_period),
                        data = location_precision, REML = FALSE)
precision_anova <- as.data.table(anova(precision_model), keep.rownames = "term")
fwrite(precision_anova, file.path(out, "Table_06_precision_response_scale_mixed_model.csv"))

##### Vertical-gradient mixed models ####

model_one_response <- function(d, formula, model_name) {
  d <- copy(d[is.finite(value) & value > 0])
  d[, `:=`(horizontal_f = factor(horizontal_position), date_f = factor(date))]
  fit <- lmer(formula, data = d, REML = FALSE)
  a <- as.data.table(anova(fit), keep.rownames = "term")
  a[, `:=`(response = as.character(first(d$response)), model = model_name)]
  emm <- emmeans(fit, ~ height | period, type = "response")
  em <- as.data.table(summary(emm, infer = c(TRUE, TRUE)))
  em[, `:=`(response = as.character(first(d$response)), model = model_name)]
  con <- as.data.table(summary(pairs(emm, adjust = "tukey"), infer = c(TRUE, TRUE)))
  con[, `:=`(response = as.character(first(d$response)), model = model_name)]
  list(anova = a, means = em, contrasts = con)
}

dense_models <- list()
for (resp in response_order) {
  d <- hourly[response == resp & period %in% c("C1 summer (Jun-Aug)", "C1 autumn (Oct)")]
  d[, period := droplevels(period)]
  dense_models[[resp]] <- model_one_response(
    d, log(value) ~ height * period + (1 | horizontal_f) + (1 | date_f),
    "Campaign 1 dense: height x period")
}

common_tb <- Reduce(intersect, lapply(split(
  hourly[height %in% c("top", "bottom")]$location_n,
  hourly[height %in% c("top", "bottom")]$period), unique))
tb_models <- list()
for (resp in response_order) {
  d <- hourly[response == resp & height %in% c("top", "bottom") &
                location_n %in% common_tb]
  d[, height := droplevels(height)]
  tb_models[[resp]] <- model_one_response(
    d, log(value) ~ height * period + (1 | horizontal_f) + (1 | date_f),
    "Shared top-bottom: height x four periods")
}

write_model_component <- function(models, component, file) {
  fwrite(rbindlist(lapply(models, `[[`, component), fill = TRUE), file.path(out, file))
}
write_model_component(dense_models, "anova", "Table_07_dense_vertical_model_ANOVA.csv")
write_model_component(dense_models, "means", "Table_07b_dense_vertical_emmeans.csv")
write_model_component(dense_models, "contrasts", "Table_07c_dense_vertical_contrasts.csv")
write_model_component(tb_models, "anova", "Table_08_top_bottom_period_model_ANOVA.csv")
write_model_component(tb_models, "means", "Table_08b_top_bottom_period_emmeans.csv")
write_model_component(tb_models, "contrasts", "Table_08c_top_bottom_period_contrasts.csv")

gradient_position <- c1_hourly[, .(mean = mean(value)),
  by = .(period, response, horizontal_position, height)]
gradient_wide <- dcast(gradient_position,
  period + response + horizontal_position ~ height, value.var = "mean")
gradient_wide[, monotonic_increase := top > middle & middle > bottom]
gradient_consistency <- gradient_wide[, .(
  horizontal_positions = .N,
  monotonic_positions = sum(monotonic_increase, na.rm = TRUE),
  monotonic_pct = 100 * mean(monotonic_increase, na.rm = TRUE),
  median_top_bottom_pct = median(100 * (top - bottom) / bottom, na.rm = TRUE),
  median_middle_bottom_pct = median(100 * (middle - bottom) / bottom, na.rm = TRUE)
), by = .(period, response)]
fwrite(gradient_wide, file.path(out, "Table_09_horizontal_gradient_checks.csv"))
fwrite(gradient_consistency, file.path(out, "Table_09b_gradient_consistency_summary.csv"))

##### Height-by-period descriptive summaries ####

height_period_summary <- hourly[, .(
  n_hour_location = .N, mean = mean(value), median = median(value),
  sd = sd(value), cv_pct = 100 * sd(value) / abs(mean(value))
), by = .(period, season, response, height)]
fwrite(height_period_summary, file.path(out, "Table_10_height_period_descriptives.csv"))

##### Exploratory westerly-wind analysis without external concentration points ####

weather_file <- file.path(workflow, "clean_data", "manuscript3_weather",
  "dwd_potsdam_03987_hourly_2024_campaign_periods.csv")
weather <- fread(weather_file)
weather[, hour_utc := ymd_hms(DATE.TIME, tz = "UTC", quiet = TRUE)]
weather[, westerly := wind_direction_deg >= 225 & wind_direction_deg < 315]

internal_hour <- hourly[, .(value = mean(value)), by = .(period, response, hour_utc)]
wind_data <- merge(internal_hour,
  weather[, .(hour_utc, wind_speed_ms, wind_direction_deg, temperature_C,
              relative_humidity_pct, precipitation_mm, westerly)],
  by = "hour_utc", all = FALSE)
wind_data <- wind_data[response %in% c("CO2", "CH4", "NH3") & value > 0]

wind_summary <- wind_data[, .(
  n_hours = .N, mean = mean(value), median = median(value), sd = sd(value)
), by = .(period, response, westerly)]
wind_tests <- wind_data[, {
  if (uniqueN(westerly) < 2L) list(statistic = NA_real_, p_value = NA_real_)
  else {
    wt <- wilcox.test(value ~ westerly, exact = FALSE)
    list(statistic = unname(wt$statistic), p_value = wt$p.value)
  }
}, by = .(period, response)]

wind_models <- wind_data[, {
  fit <- lm(log(value) ~ wind_speed_ms + westerly + temperature_C + relative_humidity_pct)
  cf <- as.data.table(summary(fit)$coefficients, keep.rownames = "term")
  setnames(cf, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"),
           c("estimate", "std_error", "t_value", "p_value"))
  cf
}, by = .(period, response)]
fwrite(wind_summary, file.path(out, "Table_11_westerly_wind_descriptives.csv"))
fwrite(wind_tests, file.path(out, "Table_11b_westerly_wind_tests.csv"))
fwrite(wind_models, file.path(out, "Table_11c_westerly_adjusted_models.csv"))

##### Manuscript v12: network-versus-SP25 statistics and uniform figures ####

v12_out <- file.path(workflow, "clean_data", "manuscript3_v12")
v12_fig <- file.path(workflow, "plots", "manuscript3_v12")
dir.create(file.path(v12_out, "statistical_tables"), recursive = TRUE,
           showWarnings = FALSE)
dir.create(file.path(v12_out, "audit"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(v12_fig, "main"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(v12_fig, "supplementary"), recursive = TRUE,
           showWarnings = FALSE)

daily_location <- long[, .(daily_mean = mean(value)),
  by = .(campaign, date, response, location_n, horizontal_position, height)]
daily_network <- daily_location[, .(daily_mean = mean(daily_mean)),
  by = .(campaign, date, response)]
daily_sp25 <- daily_location[location_n == 25L,
  .(campaign, date, response, daily_mean)]

describe_vector <- function(v) {
  v <- v[is.finite(v)]
  n <- length(v)
  se <- if (n > 1L) sd(v) / sqrt(n) else NA_real_
  critical <- if (n > 1L) qt(0.975, n - 1L) else NA_real_
  list(
    mean = mean(v), sd = sd(v), cv_pct = 100 * sd(v) / abs(mean(v)),
    median = median(v), mad = mad(v, constant = 1),
    minimum = min(v), q1 = quantile(v, 0.25, names = FALSE),
    q3 = quantile(v, 0.75, names = FALSE), maximum = max(v),
    ci95_low = mean(v) - critical * se,
    ci95_high = mean(v) + critical * se,
    observations = n, measurement_days = uniqueN(v)
  )
}

network_stats <- daily_network[, describe_vector(daily_mean),
  by = .(campaign, response)][, scope := "Whole network"]
sp25_stats <- daily_sp25[, describe_vector(daily_mean),
  by = .(campaign, response)][, scope := "SP 25"]
summary_wide_source <- rbindlist(list(network_stats, sp25_stats), fill = TRUE)

summary_long <- melt(summary_wide_source,
  id.vars = c("campaign", "response", "scope"),
  variable.name = "statistic", value.name = "value")
summary_long[, response := factor(as.character(response), levels = response_order)]
summary_long[, statistic := factor(statistic, levels = c(
  "mean", "sd", "cv_pct", "median", "mad", "minimum", "q1", "q3",
  "maximum", "ci95_low", "ci95_high", "observations", "measurement_days"
))]
setorder(summary_long, response, statistic, campaign, scope)
summary_table <- dcast(summary_long, response + statistic ~ campaign + scope,
                       value.var = "value")
fwrite(summary_long, file.path(v12_out, "statistical_tables",
  "Table_v12_descriptive_statistics_long.csv"))
fwrite(summary_table, file.path(v12_out, "statistical_tables",
  "Table_v12_descriptive_statistics_full_page.csv"))

sp25_comparison <- merge(daily_network, daily_sp25,
  by = c("campaign", "date", "response"), suffixes = c("_network", "_sp25"))
sp25_comparison[, RE_pct := 100 * (daily_mean_sp25 - daily_mean_network) /
  daily_mean_network]
sp25_performance <- sp25_comparison[, .(
  paired_days = .N,
  mean_signed_RE_pct = mean(RE_pct),
  RE_se_pct = sd(RE_pct) / sqrt(.N),
  equivalence_90CI_low_pct = mean(RE_pct) - qt(0.95, .N - 1L) * sd(RE_pct) / sqrt(.N),
  equivalence_90CI_high_pct = mean(RE_pct) + qt(0.95, .N - 1L) * sd(RE_pct) / sqrt(.N),
  mean_absolute_RE_pct = mean(abs(RE_pct)),
  RMSE_pct = sqrt(mean(RE_pct^2)),
  median_RE_pct = median(RE_pct),
  RE_sd_pct = sd(RE_pct),
  days_within_5pct = 100 * mean(abs(RE_pct) <= 5),
  days_within_10pct = 100 * mean(abs(RE_pct) <= 10),
  spearman_rho = cor(daily_mean_sp25, daily_mean_network, method = "spearman")
), by = .(campaign, response)]
sp25_performance[, equivalent_within_5pct :=
  equivalence_90CI_low_pct > -5 & equivalence_90CI_high_pct < 5]
fwrite(sp25_performance, file.path(v12_out, "statistical_tables",
  "Table_v12_SP25_network_performance.csv"))

re_table_rows <- rbindlist(list(
  sp25_performance[, .(campaign, response, scope = "Whole network",
                       statistic = "SP25_RE_pct", value = 0)],
  sp25_performance[, .(campaign, response, scope = "SP 25",
                       statistic = "SP25_RE_pct", value = mean_signed_RE_pct)]
))
summary_long <- rbindlist(list(
  as.data.table(summary_long)[, statistic := as.character(statistic)],
  re_table_rows
), use.names = TRUE, fill = TRUE)
summary_long[, response := factor(as.character(response), levels = response_order)]
summary_long[, statistic := factor(as.character(statistic), levels = c(
  "mean", "sd", "cv_pct", "median", "mad", "minimum", "q1", "q3",
  "maximum", "ci95_low", "ci95_high", "observations", "measurement_days",
  "SP25_RE_pct"
))]
setorder(summary_long, response, statistic, campaign, scope)
summary_table <- dcast(summary_long, response + statistic ~ campaign + scope,
                       value.var = "value")
fwrite(summary_long, file.path(v12_out, "statistical_tables",
  "Table_v12_descriptive_statistics_long.csv"))
fwrite(summary_table, file.path(v12_out, "statistical_tables",
  "Table_v12_descriptive_statistics_full_page.csv"))

man_table_dir <- file.path(workflow,
  "Manuscript_3_Mapping_high_resolution_concentration_Latex_draft", "tables")
dir.create(man_table_dir, recursive = TRUE, showWarnings = FALSE)

metric_labels <- c(
  mean = "Mean", sd = "SD", cv_pct = "CV (\\%)", median = "Median",
  mad = "MAD", minimum = "Minimum", q1 = "Q1", q3 = "Q3",
  maximum = "Maximum", ci95_low = "95\\% CI lower",
  ci95_high = "95\\% CI upper", observations = "Daily estimates",
  measurement_days = "Measurement days", SP25_RE_pct = "SP25 RE (\\%)"
)
response_tex <- c(
  CO2 = "$\\mathrm{CO_2}$", NH3_CO2 = "$\\mathrm{NH_3}/\\mathrm{CO_2}$",
  CH4 = "$\\mathrm{CH_4}$", CH4_CO2 = "$\\mathrm{CH_4}/\\mathrm{CO_2}$",
  NH3 = "$\\mathrm{NH_3}$", NH3_CH4 = "$\\mathrm{NH_3}/\\mathrm{CH_4}$"
)
format_table_value <- function(value, response, statistic) {
  if (!is.finite(value)) return("--")
  if (statistic %chin% c("observations", "measurement_days"))
    return(format(round(value), big.mark = ",", scientific = FALSE))
  digits <- if (response %chin% c("CO2", "CH4", "NH3")) 2L else 5L
  if (statistic %chin% c("cv_pct", "SP25_RE_pct")) digits <- 2L
  formatC(value, digits = digits, format = "f")
}

table_columns <- c("Campaign 1_Whole network", "Campaign 1_SP 25",
                   "Campaign 2_Whole network", "Campaign 2_SP 25")
latex_lines <- c(
  "\\begin{sidewaystable}[p]", "\\centering", "\\scriptsize",
  "\\caption{Descriptive statistics for daily whole-network estimates and sampling point (SP) 25. Gas ratios are dimensionless and are not multiplied by 100. Whole-network daily estimates weight each available SP equally. Relative error (RE) is signed for SP 25 against the contemporaneous network mean.}",
  "\\label{tab:full_descriptive_statistics}",
  "\\setlength{\\tabcolsep}{3.2pt}"
)
for (idx in seq_along(response_order)) {
  response_name <- response_order[idx]
  d <- summary_table[as.character(response) == response_name]
  block <- c(
    paste0("\\begin{minipage}[t]{0.48\\textwidth}\\centering\\textbf{",
           response_tex[[response_name]], "}\\par\\vspace{1mm}"),
    "\\begin{tabular}{lrrrr}", "\\toprule",
    "Statistic & C1 network & C1 SP25 & C2 network & C2 SP25\\\\", "\\midrule"
  )
  for (metric in levels(summary_long$statistic)) {
    row <- d[as.character(statistic) == metric]
    values <- if (nrow(row)) vapply(table_columns, function(column_name)
      format_table_value(row[[column_name]], response_name, metric), character(1))
      else rep("--", 4)
    block <- c(block, paste(metric_labels[[metric]], paste(values, collapse = " & "),
                            sep = " & "), "\\\\")
  }
  block <- c(block, "\\bottomrule", "\\end{tabular}", "\\end{minipage}")
  latex_lines <- c(latex_lines, block,
    if (idx %% 2L == 1L) "\\hfill" else if (idx < length(response_order)) "\\par\\vspace{2mm}" else "")
}
latex_lines <- c(latex_lines, "\\end{sidewaystable}")
writeLines(latex_lines, file.path(man_table_dir,
  "Table_v12_full_descriptive_statistics.tex"))

format_p <- function(p) {
  if (!is.finite(p)) return("--")
  if (p < 0.001) return("$<0.001$")
  formatC(p, digits = 3, format = "f")
}
extract_p <- function(model_list, response_name, term_name) {
  a <- model_list[[response_name]]$anova
  value <- a[term == term_name][["Pr(>F)"]]
  if (length(value)) value[1] else NA_real_
}
model_table_lines <- c(
  "\\begin{table}[htbp]", "\\centering", "\\small",
  "\\caption{Omnibus tests from the linear mixed-effects models. The dense model used Campaign 1 and all three heights; the shared model used top and bottom SPs available across the four seasonal periods. Models were fitted to log-transformed responses with horizontal position and date as random intercepts.}",
  "\\label{tab:mixed_model_tests}",
  "\\begin{tabular}{lrrrrrr}", "\\toprule",
  "& \\multicolumn{3}{c}{Campaign 1 dense} & \\multicolumn{3}{c}{Shared top--bottom}\\\\",
  "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}",
  "Response & Height & Period & Height $\\times$ period & Height & Period & Height $\\times$ period\\\\",
  "\\midrule"
)
for (response_name in response_order) {
  values <- c(
    extract_p(dense_models, response_name, "height"),
    extract_p(dense_models, response_name, "period"),
    extract_p(dense_models, response_name, "height:period"),
    extract_p(tb_models, response_name, "height"),
    extract_p(tb_models, response_name, "period"),
    extract_p(tb_models, response_name, "height:period")
  )
  model_table_lines <- c(model_table_lines,
    paste(response_tex[[response_name]], paste(vapply(values, format_p, character(1)),
      collapse = " & "), sep = " & "), "\\\\")
}
model_table_lines <- c(model_table_lines, "\\bottomrule", "\\end{tabular}",
                       "\\end{table}")
writeLines(model_table_lines, file.path(man_table_dir,
  "Table_v12_mixed_model_omnibus.tex"))

sp25_table_lines <- c(
  "\\begin{table}[htbp]", "\\centering", "\\scriptsize",
  "\\caption{Performance and equivalence of SP25 against the equally weighted daily whole-network mean. Equivalence required the 90\\% confidence interval (CI) for signed relative error (RE) to lie entirely within $-5$ to $+5$\\%.}",
  "\\label{tab:sp25_equivalence}",
  "\\begin{tabular}{llrrrrrr}", "\\toprule",
  "Campaign & Response & Signed RE & 90\\% CI & MAE & RMSE & $\\rho$ & Equivalent\\\\",
  "& & (\\%) & (\\%) & (\\%) & (\\%) & & $\\pm5$\\%\\\\", "\\midrule"
)
for (campaign_name in c("Campaign 1", "Campaign 2")) {
  for (response_name in response_order) {
    row <- sp25_performance[campaign == campaign_name &
                              as.character(response) == response_name]
    if (!nrow(row)) next
    ci_text <- sprintf("%.2f to %.2f", row$equivalence_90CI_low_pct,
                       row$equivalence_90CI_high_pct)
    sp25_table_lines <- c(sp25_table_lines, paste(
      campaign_name, response_tex[[response_name]],
      sprintf("%.2f", row$mean_signed_RE_pct), ci_text,
      sprintf("%.2f", row$mean_absolute_RE_pct), sprintf("%.2f", row$RMSE_pct),
      sprintf("%.3f", row$spearman_rho),
      ifelse(row$equivalent_within_5pct, "Yes", "No"), sep = " & "), "\\\\")
  }
  if (campaign_name == "Campaign 1") sp25_table_lines <- c(sp25_table_lines, "\\midrule")
}
sp25_table_lines <- c(sp25_table_lines, "\\bottomrule", "\\end{tabular}",
                      "\\end{table}")
writeLines(sp25_table_lines, file.path(man_table_dir,
  "Table_v12_SP25_equivalence.tex"))

location_stats <- daily_location[, {
  z <- describe_vector(daily_mean)
  .(days = z$observations, mean = z$mean, sd = z$sd, cv_pct = z$cv_pct,
    median = z$median, mad = z$mad, ci95_low = z$ci95_low,
    ci95_high = z$ci95_high)
}, by = .(campaign, response, location_n, horizontal_position, height)]
reference_means <- location_stats[location_n == 25L,
  .(campaign, response, reference_mean = mean)]
location_stats <- reference_means[location_stats,
  on = .(campaign, response)]
location_stats[, `:=`(
  RE_pct = 100 * (mean - reference_mean) / reference_mean,
  absolute_RE_pct = abs(100 * (mean - reference_mean) / reference_mean),
  relative_MAD_pct = 100 * mad / abs(median),
  relative_CI_halfwidth_pct = 100 * ((ci95_high - ci95_low) / 2) / abs(mean)
)]
fwrite(location_stats, file.path(v12_out, "statistical_tables",
  "Table_v12_all_locations_statistics_and_SP25_RE.csv"))

response_facet_order <- factor(response_order, levels = response_order)
location_stats[, response := factor(as.character(response),
                                     levels = levels(response_facet_order))]

p_ci_campaign <- function(campaign_name) {
  ggplot(location_stats[campaign == campaign_name],
    aes(location_n, mean, colour = height)) +
    geom_errorbar(aes(ymin = ci95_low, ymax = ci95_high),
                  width = 0.16, linewidth = 0.35) +
    geom_point(size = 1.65) +
    facet_wrap(~response, scales = "free_y", ncol = 2,
               labeller = facet_labeller) +
    scale_x_continuous(breaks = seq(1, 51, by = 2), limits = c(0.5, 51.5)) +
    scale_colour_manual(values = height_colours, name = "Height") +
    labs(x = "Sampling point", y = "Mean (95% CI)", title = campaign_name) +
    theme_bw(base_size = 9) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
}

metric_colours <- c("#2166ac", "#1fa187", "#7ad151", "#fdae61", "#b2182b")
metric_breaks <- c(0, 25, 50, 75, 100)
p_metric_heatmap <- function(campaign_name, metric, legend_title) {
  d <- copy(location_stats[campaign == campaign_name])
  d[, display_value := pmin(pmax(get(metric), 0), 100)]
  d[, height_y := c(bottom = 1, middle = 2, top = 3)[as.character(height)]]
  ggplot(d, aes(horizontal_position, height_y, fill = display_value)) +
    geom_tile(width = 0.94, height = 0.94, colour = "white", linewidth = 0.25) +
    geom_text(aes(label = location_n), size = 2.25) +
    facet_wrap(~response, ncol = 2, labeller = facet_labeller) +
    scale_x_continuous(breaks = NULL, labels = NULL,
                       limits = c(0.5, 17.5), expand = c(0, 0)) +
    scale_y_continuous(breaks = 1:3, labels = c("bottom", "middle", "top"),
                       limits = c(0.5, 3.5), expand = c(0, 0)) +
    scale_fill_gradientn(colours = metric_colours, limits = c(0, 100),
                         breaks = metric_breaks, name = legend_title) +
    coord_fixed(ratio = 1) +
    labs(x = NULL, y = "Height", title = campaign_name) +
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank(), legend.position = "bottom")
}

v12_plots <- list(
  C1_mean_CI = p_ci_campaign("Campaign 1"),
  C2_mean_CI = p_ci_campaign("Campaign 2"),
  C1_RE_SP25 = p_metric_heatmap("Campaign 1", "absolute_RE_pct", "Absolute RE (%)"),
  C2_RE_SP25 = p_metric_heatmap("Campaign 2", "absolute_RE_pct", "Absolute RE (%)"),
  C1_CV = p_metric_heatmap("Campaign 1", "cv_pct", "CV (%)"),
  C2_CV = p_metric_heatmap("Campaign 2", "cv_pct", "CV (%)"),
  C1_MAD = p_metric_heatmap("Campaign 1", "relative_MAD_pct", "Relative MAD (%)"),
  C2_MAD = p_metric_heatmap("Campaign 2", "relative_MAD_pct", "Relative MAD (%)"),
  C1_CI = p_metric_heatmap("Campaign 1", "relative_CI_halfwidth_pct", "Relative CI (%)"),
  C2_CI = p_metric_heatmap("Campaign 2", "relative_CI_halfwidth_pct", "Relative CI (%)")
)

for (nm in names(v12_plots)) {
  dimensions <- if (grepl("mean_CI", nm)) c(11.2, 10.2) else c(10.5, 8.4)
  target <- file.path(v12_fig, "main", paste0("Fig_v12_", nm, ".png"))
  ggsave(target, v12_plots[[nm]], width = dimensions[1], height = dimensions[2],
         dpi = 350, bg = "white")
  file.copy(target, file.path(man_fig, basename(target)), overwrite = TRUE)
}

fwrite(data.table(
  item = c("ratios", "reference", "heatmap_scale", "CI_basis"),
  definition = c(
    "Dimensionless; not multiplied by 100",
    "Sampling point 25; central pragmatic candidate",
    "0-100%; values above 100% clipped only for display",
    "95% t interval calculated from daily sampling-point means"
  )), file.path(v12_out, "audit", "analysis_decisions_v12.csv"))

##### Save figures ####

save_both <- function(filename, plot, width, height) {
  for (d in c(fig, man_fig)) {
    ggsave(file.path(d, paste0(filename, ".png")), plot,
           width = width, height = height, dpi = 350, bg = "white")
  }
}

save_both("Fig_v12_dense_means_95CI", p_mean, 12.2, 7.4)
save_both("Fig_v12_dense_CV", p_cv, 12.2, 7.4)
save_both("Fig_v12_single_height_relative_error", p_height, 11.5, 7.1)
save_both("Fig_v12_Shannon_entropy", p_entropy, 10.5, 8.4)

writeLines(capture.output(sessionInfo()), file.path(v12_out, "audit", "sessionInfo_v12.txt"))
writeLines(c(
  "Manuscript 3 v12 unified analysis",
  "Internal locations only: numeric locations 1-51; source location s/52 excluded.",
  "Campaign 1 summer: June-August 2024; Campaign 1 autumn: October 2024.",
  "Campaign 2 was split into late autumn (November) and winter (December).",
  "Middle height: 3.6 m above floor; bottom height: 2.6 m above floor.",
  "Top height: 0.60 m below the local roof, not represented by one absolute elevation.",
  "Gas ratios are dimensionless and are not multiplied by 100.",
  "Westerly DWD sector: 225 <= direction < 315 degrees; regional context only.",
  "DWD wind is not barn-local velocity or turbulence intensity."
), file.path(v12_out, "audit", "analysis_readme_v12.txt"))

cat("Internal rows:", nrow(x), "\n")
cat("Complete dense two-hour blocks:", nrow(complete_blocks), "\n")
cat("Common top/bottom locations across four periods:", length(common_tb), "\n")
print(gradient_consistency)
print(wind_summary)
