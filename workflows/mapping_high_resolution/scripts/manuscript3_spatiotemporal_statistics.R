###############################################################################
# Manuscript 3: inferential analysis of spatial and temporal gas heterogeneity
#
# Primary question:
#   How strongly and in what spatial/temporal structure do CO2, CH4, and NH3
#   concentrations vary within the naturally ventilated dairy barn?
#
# Secondary question:
#   Can sampling heights be reduced after the heterogeneity structure has been
#   quantified?
#
# No source data are modified. All tables, figures, model summaries, and the
# session record are written to dedicated manuscript statistical-output folders.
###############################################################################

required_packages <- c(
  "data.table", "ggplot2", "lme4", "lmerTest", "emmeans", "mgcv",
  "broom.mixed", "boot", "DescTools", "lubridate"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing R package(s): ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(mgcv)
  library(lubridate)
})

set.seed(20260729)
options(contrasts = c("contr.sum", "contr.poly"))


###############################################################################
##### PATHS AND CONSTANTS
###############################################################################

workflow_dir <- paste0(
  "D:/Data_Analysis_R/Gas_Concentration_Emission_R/",
  "workflows/mapping_high_resolution"
)
input_dir <- file.path(workflow_dir, "clean_data", "manuscript1")
input_file <- file.path(
  input_dir,
  "manuscript1_campaign1_2_analytical_data.csv"
)
configuration_file <- file.path(
  input_dir,
  "campaign1_configuration_block_results.csv"
)

table_dir <- file.path(input_dir, "statistical_analysis")
figure_dir <- file.path(
  workflow_dir,
  "plots",
  "manuscript1",
  "statistical_analysis"
)
diagnostic_dir <- file.path(figure_dir, "diagnostics")
model_dir <- file.path(table_dir, "models")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(diagnostic_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

gas_names <- c("CO2", "CH4", "NH3")
ratio_names <- c("CH4_CO2", "NH3_CO2", "NH3_CH4")
response_names <- c(gas_names, ratio_names)

height_levels <- c("bottom", "mid", "top")
height_colours <- c(
  top = "orange",
  mid = "green3",
  bottom = "steelblue1"
)

top_locations <- seq(1L, 49L, by = 3L)
mid_locations <- seq(2L, 50L, by = 3L)
bottom_locations <- seq(3L, 51L, by = 3L)

barn_timezone <- "Europe/Berlin"
day_start_hour <- 6L
night_start_hour <- 18L

bootstrap_replicates <- 2000L
equivalence_margin_primary <- 5
equivalence_margin_sensitivity <- 10


###############################################################################
##### HELPERS
###############################################################################

safe_numeric <- function(x) suppressWarnings(as.numeric(x))


height_from_location <- function(location) {
  location_number <- suppressWarnings(as.integer(as.character(location)))
  fifelse(
    location_number %in% top_locations, "top",
    fifelse(
      location_number %in% mid_locations, "mid",
      fifelse(location_number %in% bottom_locations, "bottom", NA_character_)
    )
  )
}


percent_from_log <- function(x) 100 * (exp(x) - 1)


bootstrap_mean_ci <- function(x, replicates = bootstrap_replicates) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) {
    return(list(mean = NA_real_, lower = NA_real_, upper = NA_real_))
  }
  boot_object <- boot::boot(
    data = x,
    statistic = function(values, index) mean(values[index]),
    R = replicates
  )
  interval <- tryCatch(
    boot::boot.ci(boot_object, type = "perc")$percent[4:5],
    error = function(e) quantile(
      boot_object$t,
      probs = c(0.025, 0.975),
      na.rm = TRUE
    )
  )
  list(
    mean = mean(x),
    lower = unname(interval[1]),
    upper = unname(interval[2])
  )
}


model_diagnostics <- function(model, response, model_family, model_data = NULL) {
  residual_values <- residuals(model)
  fitted_values <- fitted(model)
  lag1_correlation <- NA_real_
  time_column <- if (
    !is.null(model_data) && "DATE.TIME" %in% names(model_data)
  ) {
    "DATE.TIME"
  } else if (
    !is.null(model_data) && "block_2h" %in% names(model_data)
  ) {
    "block_2h"
  } else {
    NA_character_
  }
  if (
    !is.null(model_data) &&
      "location" %in% names(model_data) &&
      !is.na(time_column) &&
      nrow(model_data) == length(residual_values)
  ) {
    residual_data <- data.table(
      location = model_data$location,
      observation_time = model_data[[time_column]],
      residual = residual_values
    )
    setorder(residual_data, location, observation_time)
    lag1_by_location <- residual_data[
      ,
      .(
        lag1 = if (.N >= 3L) {
          cor(residual[-1L], residual[-.N], use = "complete.obs")
        } else {
          NA_real_
        }
      ),
      by = location
    ]
    lag1_correlation <- median(lag1_by_location$lag1, na.rm = TRUE)
  }
  qq_correlation <- cor(
    sort(residual_values),
    qnorm(ppoints(length(residual_values))),
    use = "complete.obs"
  )
  data.table(
    response = response,
    model_family = model_family,
    n = length(residual_values),
    residual_mean = mean(residual_values, na.rm = TRUE),
    residual_sd = sd(residual_values, na.rm = TRUE),
    residual_fitted_correlation = cor(
      residual_values,
      fitted_values,
      use = "complete.obs"
    ),
    median_location_lag1_residual_correlation = lag1_correlation,
    qq_correlation = qq_correlation,
    singular = if (inherits(model, "merMod")) {
      lme4::isSingular(model, tol = 1e-4)
    } else {
      NA
    }
  )
}


save_diagnostic_plot <- function(model, response, model_family) {
  residual_values <- residuals(model)
  fitted_values <- fitted(model)
  output_path <- file.path(
    diagnostic_dir,
    paste0("diagnostic_", model_family, "_", response, ".png")
  )
  png(output_path, width = 1800, height = 900, res = 180)
  par(mfrow = c(1, 2))
  plot(
    fitted_values,
    residual_values,
    pch = 16,
    cex = 0.45,
    col = adjustcolor("black", alpha.f = 0.25),
    xlab = "Fitted value",
    ylab = "Residual",
    main = paste(response, "residuals versus fitted")
  )
  abline(h = 0, lty = 2, col = "red")
  qqnorm(
    residual_values,
    pch = 16,
    cex = 0.45,
    col = adjustcolor("black", alpha.f = 0.25),
    main = paste(response, "normal Q-Q")
  )
  qqline(residual_values, col = "red")
  dev.off()
  invisible(output_path)
}


extract_lmer_tests <- function(model, response, model_family) {
  test_table <- as.data.table(
    as.data.frame(anova(model, type = 3)),
    keep.rownames = "term"
  )
  setnames(
    test_table,
    old = intersect(
      c("NumDF", "DenDF", "F value", "Pr(>F)"),
      names(test_table)
    ),
    new = c("num_df", "den_df", "F_value", "p_value")[
      match(
        intersect(c("NumDF", "DenDF", "F value", "Pr(>F)"), names(test_table)),
        c("NumDF", "DenDF", "F value", "Pr(>F)")
      )
    ]
  )
  test_table[, `:=`(
    response = response,
    model_family = model_family
  )]
  setcolorder(
    test_table,
    c(
      "response", "model_family", "term",
      setdiff(names(test_table), c("response", "model_family", "term"))
    )
  )
  test_table
}


extract_height_contrasts <- function(model, response, model_family) {
  marginal <- emmeans::emmeans(model, ~ height, weights = "equal")
  pairwise_result <- as.data.table(
    summary(
      pairs(marginal, adjust = "tukey"),
      infer = c(TRUE, TRUE)
    )
  )
  lower_column <- intersect(c("lower.CL", "asymp.LCL"), names(pairwise_result))[1]
  upper_column <- intersect(c("upper.CL", "asymp.UCL"), names(pairwise_result))[1]
  pairwise_result[, `:=`(
    response = response,
    model_family = model_family,
    percent_difference = percent_from_log(estimate),
    percent_lower = percent_from_log(get(lower_column)),
    percent_upper = percent_from_log(get(upper_column)),
    adjustment = "Tukey"
  )]
  pairwise_result
}


fit_spatial_model <- function(model_data, response, model_family) {
  subset_data <- model_data[
    is.finite(get(response))
  ]
  subset_data[, response_value := get(response)]

  model <- lmerTest::lmer(
    response_value ~ height * horizontal_position + (1 | block_2h),
    data = subset_data,
    REML = FALSE
  )

  saveRDS(
    model,
    file.path(model_dir, paste0("model_", model_family, "_", response, ".rds"))
  )
  save_diagnostic_plot(model, response, model_family)

  list(
    model = model,
    tests = extract_lmer_tests(model, response, model_family),
    contrasts = extract_height_contrasts(model, response, model_family),
    diagnostics = model_diagnostics(
      model,
      response,
      model_family,
      subset_data
    )
  )
}


###############################################################################
##### IMPORT AND VALIDATE
###############################################################################

if (!file.exists(input_file)) {
  stop("Analytical input file not found: ", input_file)
}
if (!file.exists(configuration_file)) {
  stop("Configuration result file not found: ", configuration_file)
}

data <- fread(
  input_file,
  na.strings = c("", "NA", "NaN"),
  colClasses = c("DATE.TIME" = "character")
)

required_columns <- c(
  "DATE.TIME", "campaign", "analyser", "location", "horizontal_position",
  "vgroup", "CO2_corr", "CH4_corr", "NH3_corr"
)
missing_columns <- setdiff(required_columns, names(data))
if (length(missing_columns) > 0L) {
  stop("Analytical dataset is missing: ", paste(missing_columns, collapse = ", "))
}

data[, DATE.TIME := as.POSIXct(
  DATE.TIME,
  format = "%Y-%m-%d %H:%M:%S",
  tz = barn_timezone
)]
data[, location := as.character(location)]
data[, location_number := suppressWarnings(as.integer(location))]
data[, height := factor(
  fifelse(is.na(vgroup), height_from_location(location), as.character(vgroup)),
  levels = height_levels
)]
data[, horizontal_position := factor(
  as.integer(as.character(horizontal_position)),
  levels = 1:17
)]
data[, date := as.Date(DATE.TIME, tz = barn_timezone)]
data[, block_2h := floor_date(DATE.TIME, unit = "2 hours")]
data[, hour := (
  hour(DATE.TIME) +
    minute(DATE.TIME) / 60 +
    second(DATE.TIME) / 3600
)]

for (gas in gas_names) {
  source_column <- paste0(gas, "_corr")
  data[, (gas) := fifelse(
    is.finite(get(source_column)) & get(source_column) > 0,
    safe_numeric(get(source_column)),
    NA_real_
  )]
}

data[, `:=`(
  CH4_CO2 = fifelse(CH4 > 0 & CO2 > 0, CH4 / CO2, NA_real_),
  NH3_CO2 = fifelse(NH3 > 0 & CO2 > 0, NH3 / CO2, NA_real_),
  NH3_CH4 = fifelse(NH3 > 0 & CH4 > 0, NH3 / CH4, NA_real_)
)]

internal <- data[location != "s" & !is.na(location_number)]

validation_summary <- internal[
  ,
  .(
    rows = .N,
    start = min(DATE.TIME, na.rm = TRUE),
    end = max(DATE.TIME, na.rm = TRUE),
    analysers = paste(sort(unique(analyser)), collapse = ", "),
    locations = uniqueN(location),
    missing_CO2 = sum(is.na(CO2)),
    missing_CH4 = sum(is.na(CH4)),
    missing_NH3 = sum(is.na(NH3))
  ),
  by = campaign
]
fwrite(
  validation_summary,
  file.path(table_dir, "Table_data_validation_summary.csv")
)


###############################################################################
##### CAMPAIGN 1 COMPLETE TWO-HOUR LOCATION BLOCKS
###############################################################################

campaign1 <- internal[campaign == "Campaign 1"]

block_location <- campaign1[
  ,
  .(
    CO2 = mean(CO2, na.rm = TRUE),
    CH4 = mean(CH4, na.rm = TRUE),
    NH3 = mean(NH3, na.rm = TRUE)
  ),
  by = .(
    block_2h,
    date,
    location,
    location_number,
    horizontal_position,
    height
  )
]

for (gas in gas_names) {
  set(
    block_location,
    i = which(!is.finite(block_location[[gas]])),
    j = gas,
    value = NA_real_
  )
}

complete_blocks <- block_location[
  complete.cases(CO2, CH4, NH3),
  .(n_locations = uniqueN(location)),
  by = block_2h
][n_locations == 51L, block_2h]

complete <- block_location[block_2h %in% complete_blocks]
complete[, `:=`(
  log_CO2 = log(CO2),
  log_CH4 = log(CH4),
  log_NH3 = log(NH3),
  log_CH4_CO2 = log(CH4 / CO2),
  log_NH3_CO2 = log(NH3 / CO2),
  log_NH3_CH4 = log(NH3 / CH4)
)]

model_response_map <- c(
  CO2 = "log_CO2",
  CH4 = "log_CH4",
  NH3 = "log_NH3",
  CH4_CO2 = "log_CH4_CO2",
  NH3_CO2 = "log_NH3_CO2",
  NH3_CH4 = "log_NH3_CH4"
)


###############################################################################
##### PRIMARY SPATIAL MODELS
###############################################################################

spatial_results <- lapply(names(model_response_map), function(response_label) {
  fit_spatial_model(
    complete,
    model_response_map[[response_label]],
    if (response_label %in% gas_names) "spatial_concentration" else "spatial_ratio"
  )
})
names(spatial_results) <- names(model_response_map)

spatial_tests <- rbindlist(lapply(spatial_results, `[[`, "tests"), fill = TRUE)
spatial_contrasts <- rbindlist(
  lapply(spatial_results, `[[`, "contrasts"),
  fill = TRUE
)
spatial_diagnostics <- rbindlist(
  lapply(spatial_results, `[[`, "diagnostics"),
  fill = TRUE
)

spatial_tests[, response := names(model_response_map)[
  match(response, unname(model_response_map))
]]
spatial_contrasts[, response := names(model_response_map)[
  match(response, unname(model_response_map))
]]
spatial_diagnostics[, response := names(model_response_map)[
  match(response, unname(model_response_map))
]]

fwrite(spatial_tests, file.path(table_dir, "Table_spatial_model_tests.csv"))
fwrite(
  spatial_contrasts,
  file.path(table_dir, "Table_spatial_height_contrasts.csv")
)


###############################################################################
##### DESCRIPTIVE HEIGHT AND RATIO SUMMARIES
###############################################################################

height_summary <- complete[
  ,
  .(
    mean_CO2 = mean(CO2, na.rm = TRUE),
    sd_CO2 = sd(CO2, na.rm = TRUE),
    mean_CH4 = mean(CH4, na.rm = TRUE),
    sd_CH4 = sd(CH4, na.rm = TRUE),
    mean_NH3 = mean(NH3, na.rm = TRUE),
    sd_NH3 = sd(NH3, na.rm = TRUE),
    CH4_CO2_pct = 100 * mean(CH4 / CO2, na.rm = TRUE),
    NH3_CO2_pct = 100 * mean(NH3 / CO2, na.rm = TRUE),
    NH3_CH4_pct = 100 * mean(NH3 / CH4, na.rm = TRUE)
  ),
  by = height
]
fwrite(
  height_summary,
  file.path(table_dir, "Table_campaign1_height_and_ratio_summary.csv")
)


###############################################################################
##### HETEROGENEITY AND CONFIGURATION EQUIVALENCE
###############################################################################

configuration <- fread(
  configuration_file,
  na.strings = c("", "NA", "NaN"),
  colClasses = c("block_2h" = "character")
)
configuration[, block_2h := as.POSIXct(
  block_2h,
  format = "%Y-%m-%d %H:%M:%S",
  tz = barn_timezone
)]

heterogeneity_summary <- configuration[
  ,
  {
    estimate_ci <- bootstrap_mean_ci(heterogeneity)
    deviation_ci <- bootstrap_mean_ci(heterogeneity_deviation_pct)
    list(
      n_blocks = .N,
      mean_sd_log = estimate_ci$mean,
      mean_sd_log_ci_lower = estimate_ci$lower,
      mean_sd_log_ci_upper = estimate_ci$upper,
      mean_heterogeneity_deviation_pct = deviation_ci$mean,
      heterogeneity_deviation_ci_lower = deviation_ci$lower,
      heterogeneity_deviation_ci_upper = deviation_ci$upper
    )
  },
  by = .(gas, configuration)
]
fwrite(
  heterogeneity_summary,
  file.path(table_dir, "Table_heterogeneity_by_configuration.csv")
)

equivalence_results <- configuration[
  configuration == "Top + bottom (34)",
  {
    test90 <- t.test(relative_deviation_pct, conf.level = 0.90)
    list(
      n_blocks = .N,
      mean_relative_difference_pct = mean(relative_deviation_pct),
      ci90_lower = test90$conf.int[1],
      ci90_upper = test90$conf.int[2],
      equivalent_margin_5pct = (
        test90$conf.int[1] > -equivalence_margin_primary &&
          test90$conf.int[2] < equivalence_margin_primary
      ),
      equivalent_margin_10pct = (
        test90$conf.int[1] > -equivalence_margin_sensitivity &&
          test90$conf.int[2] < equivalence_margin_sensitivity
      ),
      proportion_within_5pct = mean(
        abs(relative_deviation_pct) <= equivalence_margin_primary
      ),
      proportion_within_10pct = mean(
        abs(relative_deviation_pct) <= equivalence_margin_sensitivity
      )
    )
  },
  by = gas
]
fwrite(
  equivalence_results,
  file.path(table_dir, "Table_top_bottom_equivalence.csv")
)


###############################################################################
##### DAY/NIGHT CLASSIFICATION AND MIXED MODELS
###############################################################################

# Fixed clock-time classification requested by the author. "Day" represents
# 06:00-17:59 and "night" represents 18:00-05:59 in Europe/Berlin local time.
# These labels describe operational clock-time periods, not solar daylight and
# darkness.
complete[, day_night := factor(
  fifelse(
    hour(block_2h) >= day_start_hour &
      hour(block_2h) < night_start_hour,
    "day",
    "night"
  ),
  levels = c("night", "day")
)]

day_night_tests <- list()
day_night_contrasts <- list()
day_night_diagnostics <- list()

for (gas in gas_names) {
  response_column <- paste0("log_", gas)
  model_data <- complete[is.finite(get(response_column))]
  model_data[, response_value := get(response_column)]
  model_data[, date_factor := factor(date)]

  model <- lmerTest::lmer(
    response_value ~ day_night * height + horizontal_position +
      (1 | date_factor) + (1 | block_2h),
    data = model_data,
    REML = FALSE
  )
  saveRDS(model, file.path(model_dir, paste0("model_day_night_", gas, ".rds")))
  save_diagnostic_plot(model, gas, "day_night")

  day_night_tests[[gas]] <- extract_lmer_tests(model, gas, "day_night")

  marginal <- emmeans::emmeans(model, ~ day_night | height)
  contrast_table <- as.data.table(
    summary(
      pairs(marginal, reverse = TRUE, adjust = "holm"),
      infer = c(TRUE, TRUE)
    )
  )
  lower_column <- intersect(c("lower.CL", "asymp.LCL"), names(contrast_table))[1]
  upper_column <- intersect(c("upper.CL", "asymp.UCL"), names(contrast_table))[1]
  contrast_table[, `:=`(
    response = gas,
    percent_difference = percent_from_log(estimate),
    percent_lower = percent_from_log(get(lower_column)),
    percent_upper = percent_from_log(get(upper_column)),
    adjustment = "Holm"
  )]
  day_night_contrasts[[gas]] <- contrast_table
  day_night_diagnostics[[gas]] <- model_diagnostics(
    model,
    gas,
    "day_night",
    model_data
  )
}

day_night_tests <- rbindlist(day_night_tests, fill = TRUE)
day_night_contrasts <- rbindlist(day_night_contrasts, fill = TRUE)
day_night_diagnostics <- rbindlist(day_night_diagnostics, fill = TRUE)

fwrite(
  day_night_tests,
  file.path(table_dir, "Table_day_night_model_tests.csv")
)
fwrite(
  day_night_contrasts,
  file.path(table_dir, "Table_day_night_height_contrasts.csv")
)


###############################################################################
##### CYCLIC HOUR-OF-DAY GAMS
###############################################################################

complete[, `:=`(
  hour = (
    hour(block_2h) +
      minute(block_2h) / 60
  ),
  date_factor = factor(date)
)]

hour_tests <- list()
hour_predictions <- list()
hour_diagnostics <- list()

for (gas in gas_names) {
  response_column <- paste0("log_", gas)
  model_data <- complete[is.finite(get(response_column))]
  model_data[, response_value := get(response_column)]

  model <- mgcv::gam(
    response_value ~
      height +
      s(hour, by = height, bs = "cc", k = 8) +
      s(horizontal_position, bs = "re") +
      s(date_factor, bs = "re"),
    data = model_data,
    method = "REML",
    knots = list(hour = c(0, 24))
  )
  saveRDS(model, file.path(model_dir, paste0("model_hour_gam_", gas, ".rds")))
  save_diagnostic_plot(model, gas, "hour_gam")

  summary_model <- summary(model)
  smooth_table <- as.data.table(
    as.data.frame(summary_model$s.table),
    keep.rownames = "term"
  )
  setnames(
    smooth_table,
    old = intersect(c("edf", "Ref.df", "F", "p-value"), names(smooth_table)),
    new = c("edf", "reference_df", "F_value", "p_value")[
      match(
        intersect(c("edf", "Ref.df", "F", "p-value"), names(smooth_table)),
        c("edf", "Ref.df", "F", "p-value")
      )
    ]
  )
  smooth_table[, response := gas]
  hour_tests[[gas]] <- smooth_table

  prediction_grid <- CJ(
    hour = seq(0, 23.75, by = 0.25),
    height = factor(height_levels, levels = height_levels)
  )
  prediction_grid[, `:=`(
    horizontal_position = factor("1", levels = levels(model_data$horizontal_position)),
    date_factor = factor(
      levels(model_data$date_factor)[1],
      levels = levels(model_data$date_factor)
    )
  )]

  prediction <- predict(
    model,
    newdata = prediction_grid,
    type = "link",
    se.fit = TRUE,
    exclude = c("s(horizontal_position)", "s(date_factor)")
  )
  prediction_grid[, `:=`(
    response = gas,
    fit_log = as.numeric(prediction$fit),
    se_log = as.numeric(prediction$se.fit)
  )]
  prediction_grid[, `:=`(
    fit = exp(fit_log),
    lower = exp(fit_log - 1.96 * se_log),
    upper = exp(fit_log + 1.96 * se_log)
  )]
  hour_predictions[[gas]] <- prediction_grid
  hour_diagnostics[[gas]] <- model_diagnostics(
    model,
    gas,
    "hour_gam",
    model_data
  )
}

hour_tests <- rbindlist(hour_tests, fill = TRUE)
hour_predictions <- rbindlist(hour_predictions, fill = TRUE)
hour_predictions[, response := factor(response, levels = gas_names)]
hour_diagnostics <- rbindlist(hour_diagnostics, fill = TRUE)

fwrite(hour_tests, file.path(table_dir, "Table_hour_gam_tests.csv"))
fwrite(
  hour_predictions,
  file.path(table_dir, "Table_hour_gam_predictions.csv")
)


###############################################################################
##### MANUSCRIPT FIGURES
###############################################################################

# Vertical profiles across horizontal positions.
spatial_plot_data <- complete[
  ,
  .(
    mean = mean(c(CO2, CH4, NH3), na.rm = TRUE)
  )
]

concentration_long <- melt(
  complete,
  id.vars = c("block_2h", "horizontal_position", "height"),
  measure.vars = gas_names,
  variable.name = "gas",
  value.name = "concentration"
)
concentration_summary <- concentration_long[
  is.finite(concentration),
  .(
    mean = mean(concentration),
    se = sd(concentration) / sqrt(.N)
  ),
  by = .(gas, horizontal_position, height)
]
concentration_summary[, `:=`(
  lower = mean - 1.96 * se,
  upper = mean + 1.96 * se
)]
concentration_summary[, gas := factor(gas, levels = gas_names)]

p_vertical <- ggplot(
  concentration_summary,
  aes(
    x = as.integer(as.character(horizontal_position)),
    y = mean,
    colour = height,
    group = height
  )
) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper, fill = height),
    alpha = 0.12,
    colour = NA
  ) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.8) +
  facet_wrap(~ gas, scales = "free_y", ncol = 1) +
  scale_colour_manual(values = height_colours) +
  scale_fill_manual(values = height_colours) +
  scale_x_continuous(breaks = 1:17) +
  labs(
    x = "Horizontal position",
    y = "Concentration (ppm)",
    colour = "Height",
    fill = "Height"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.spacing = grid::unit(0.8, "lines"),
    legend.position = "top"
  )

ggsave(
  file.path(figure_dir, "Fig_spatial_vertical_profiles.png"),
  p_vertical,
  width = 9,
  height = 10,
  dpi = 400
)

# Spatial ratio profiles.
ratio_long <- complete[
  ,
  .(
    block_2h,
    horizontal_position,
    height,
    CH4_CO2 = 100 * CH4 / CO2,
    NH3_CO2 = 100 * NH3 / CO2,
    NH3_CH4 = 100 * NH3 / CH4
  )
]
ratio_long <- melt(
  ratio_long,
  id.vars = c("block_2h", "horizontal_position", "height"),
  measure.vars = ratio_names,
  variable.name = "ratio",
  value.name = "ratio_pct"
)
ratio_summary <- ratio_long[
  is.finite(ratio_pct),
  .(
    mean = mean(ratio_pct),
    se = sd(ratio_pct) / sqrt(.N)
  ),
  by = .(ratio, horizontal_position, height)
]
ratio_summary[, `:=`(
  lower = mean - 1.96 * se,
  upper = mean + 1.96 * se
)]
ratio_summary[, ratio := factor(ratio, levels = ratio_names)]

p_ratio <- ggplot(
  ratio_summary,
  aes(
    x = as.integer(as.character(horizontal_position)),
    y = mean,
    colour = height,
    group = height
  )
) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper, fill = height),
    alpha = 0.12,
    colour = NA
  ) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.8) +
  facet_wrap(
    ~ ratio,
    scales = "free_y",
    ncol = 1,
    labeller = as_labeller(c(
      CH4_CO2 = "CH4 / CO2",
      NH3_CO2 = "NH3 / CO2",
      NH3_CH4 = "NH3 / CH4"
    ))
  ) +
  scale_colour_manual(values = height_colours) +
  scale_fill_manual(values = height_colours) +
  scale_x_continuous(breaks = 1:17) +
  labs(
    x = "Horizontal position",
    y = "Ratio (%)",
    colour = "Height",
    fill = "Height"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.spacing = grid::unit(0.8, "lines"),
    legend.position = "top"
  )

ggsave(
  file.path(figure_dir, "Fig_spatial_gas_ratios.png"),
  p_ratio,
  width = 9,
  height = 10,
  dpi = 400
)

# Cyclic hour predictions.
p_hour <- ggplot(
  hour_predictions,
  aes(x = hour, y = fit, colour = height, fill = height)
) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.14, colour = NA) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ response, scales = "free_y", ncol = 1) +
  scale_colour_manual(values = height_colours) +
  scale_fill_manual(values = height_colours) +
  scale_x_continuous(breaks = seq(0, 24, by = 3), limits = c(0, 24)) +
  labs(
    x = "Hour of day",
    y = "Predicted concentration (ppm)",
    colour = "Height",
    fill = "Height"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.spacing = grid::unit(0.8, "lines"),
    legend.position = "top"
  )

ggsave(
  file.path(figure_dir, "Fig_hour_of_day_profiles.png"),
  p_hour,
  width = 9,
  height = 10,
  dpi = 400
)

# Heterogeneity across configurations.
configuration[, gas := factor(gas, levels = gas_names)]
p_heterogeneity <- ggplot(
  configuration,
  aes(x = configuration, y = heterogeneity, fill = configuration)
) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.4) +
  facet_wrap(~ gas, scales = "free_y", ncol = 1) +
  labs(
    x = "Sampling configuration",
    y = "SD of log-transformed location concentrations"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none",
    panel.spacing = grid::unit(0.8, "lines")
  )

ggsave(
  file.path(figure_dir, "Fig_heterogeneity_by_configuration.png"),
  p_heterogeneity,
  width = 11,
  height = 10,
  dpi = 400
)


###############################################################################
##### FINAL DIAGNOSTICS AND ANALYSIS RECORD
###############################################################################

all_diagnostics <- rbindlist(
  list(
    spatial_diagnostics,
    day_night_diagnostics,
    hour_diagnostics
  ),
  fill = TRUE
)
fwrite(
  all_diagnostics,
  file.path(table_dir, "Table_model_diagnostics.csv")
)

analysis_record <- c(
  "Manuscript 3 spatiotemporal statistical analysis",
  paste0("Run time: ", Sys.time()),
  paste0("Input: ", input_file),
  paste0("Campaign 1 complete two-hour blocks: ", length(complete_blocks)),
  paste0("Campaign 1 model rows: ", nrow(complete)),
  paste0(
    "Fixed clock-time classification: day = ",
    sprintf("%02d:00", day_start_hour),
    "-",
    sprintf("%02d:00", night_start_hour),
    "; night = ",
    sprintf("%02d:00", night_start_hour),
    "-",
    sprintf("%02d:00", day_start_hour),
    " (Europe/Berlin)."
  ),
  "Day/night labels represent operational time periods, not solar light status.",
  paste0("Block bootstrap replicates: ", bootstrap_replicates),
  paste0("Primary equivalence margin: +/-", equivalence_margin_primary, "%"),
  paste0(
    "Sensitivity equivalence margin: +/-",
    equivalence_margin_sensitivity, "%"
  ),
  "Concentration responses: log(CO2), log(CH4), log(NH3)",
  "Ratio responses: log(CH4/CO2), log(NH3/CO2), log(NH3/CH4)",
  "No background-corrected enhancement-ratio regression was performed.",
  "No gas concentration was interpolated.",
  "",
  capture.output(sessionInfo())
)
writeLines(
  analysis_record,
  file.path(table_dir, "statistical_analysis_readme.txt")
)

message("Statistical analysis completed.")
message("Tables: ", table_dir)
message("Figures: ", figure_dir)
