##### External wind, south-background and internal concentration analysis ####

# Purpose: test exploratory, time-aligned associations between regional DWD
# weather and Campaign 1 gas concentrations. Hourly DWD wind speed is not a
# measurement of turbulence intensity or barn-level inlet velocity.

library(data.table)
library(lubridate)
library(ggplot2)

set.seed(20260803)

workflow_dir <- "D:/Data_Analysis_R/Gas_Concentration_Emission_R/workflows/mapping_high_resolution"
campaign_dir <- file.path(workflow_dir, "clean_data", "1_campaign")
weather_file <- file.path(workflow_dir, "clean_data", "manuscript3_weather",
                          "dwd_potsdam_03987_hourly_2024_campaign_periods.csv")
output_dir <- file.path(workflow_dir, "clean_data", "manuscript3_weather")
plot_dir <- file.path(workflow_dir, "plots", "manuscript3_weather")
manuscript_figure_dir <- file.path(
  workflow_dir, "Manuscript_3_Mapping_high_resolution_concentration_Latex_draft", "figures"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manuscript_figure_dir, recursive = TRUE, showWarnings = FALSE)

summer_files <- list.files(
  file.path(campaign_dir, "recovered_crds_june_august_v03"),
  pattern = "four_analyser_CRDS_scale_corrected_.*_v03[.]csv$",
  recursive = TRUE, full.names = TRUE
)
autumn_file <- file.path(
  campaign_dir,
  "20241001.000232_20241024.143858_campaign1_raw_and_CRDS_scale_corrected.csv"
)
stopifnot(length(summer_files) == 3L, file.exists(autumn_file), file.exists(weather_file))

read_campaign <- function(file, period) {
  x <- fread(file)
  required <- c("DATE.TIME", "location", "CO2_corr", "CH4_corr", "NH3_corr")
  stopifnot(all(required %in% names(x)))
  x[, DATE.TIME := ymd_hms(DATE.TIME, tz = "Europe/Berlin", quiet = TRUE)]
  x <- x[!is.na(DATE.TIME)]
  x[, `:=`(period = period, location = as.character(location))]
  melt(
    x,
    id.vars = c("DATE.TIME", "period", "location"),
    measure.vars = c("CO2_corr", "CH4_corr", "NH3_corr"),
    variable.name = "gas", value.name = "concentration"
  )[, gas := sub("_corr$", "", gas)][]
}

gas <- rbindlist(list(
  rbindlist(lapply(summer_files, read_campaign, period = "Campaign 1: summer")),
  read_campaign(autumn_file, "Campaign 1: autumn")
), use.names = TRUE)
gas <- gas[is.finite(concentration)]

# DWD timestamps are UTC. Gas timestamps are parsed as local civil time and
# converted to UTC before hourly alignment, including daylight-saving changes.
gas[, hour_utc := floor_date(with_tz(DATE.TIME, "UTC"), "hour")]

background <- gas[location == "s", .(
  background_mean = mean(concentration, na.rm = TRUE),
  background_median = median(concentration, na.rm = TRUE),
  n_background = .N
), by = .(period, gas, hour_utc)]

inside <- gas[location != "s", .(
  inside_mean = mean(concentration, na.rm = TRUE),
  inside_median = median(concentration, na.rm = TRUE),
  n_inside = .N,
  n_locations = uniqueN(location)
), by = .(period, gas, hour_utc)]

hourly <- merge(inside, background, by = c("period", "gas", "hour_utc"), all = FALSE)
hourly[, `:=`(
  enhancement_mean = inside_mean - background_mean,
  enhancement_median = inside_median - background_median
)]

weather <- fread(weather_file)
weather[, hour_utc := ymd_hms(DATE.TIME, tz = "UTC", quiet = TRUE)]
sector_labels <- c("N", "NE", "E", "SE", "S", "SW", "W", "NW")
weather[, wind_sector := factor(
  sector_labels[((wind_direction_deg + 22.5) %/% 45) %% 8 + 1],
  levels = sector_labels
)]
hourly <- merge(
  hourly,
  weather[, .(hour_utc, wind_speed_ms, wind_direction_deg, wind_sector,
              temperature_C, relative_humidity_pct, precipitation_mm)],
  by = "hour_utc", all = FALSE
)
hourly[, date_utc := as.Date(hour_utc, tz = "UTC")]

block_boot_spearman <- function(x, y, block, reps = 300L) {
  ok <- is.finite(x) & is.finite(y) & !is.na(block)
  x <- x[ok]; y <- y[ok]; block <- block[ok]
  estimate <- suppressWarnings(cor(x, y, method = "spearman"))
  blocks <- unique(block)
  if (length(blocks) < 3L || length(x) < 8L) {
    return(c(rho = estimate, lower = NA_real_, upper = NA_real_, n = length(x), blocks = length(blocks)))
  }
  boot <- replicate(reps, {
    sampled <- sample(blocks, length(blocks), replace = TRUE)
    idx <- unlist(lapply(sampled, function(b) which(block == b)), use.names = FALSE)
    suppressWarnings(cor(x[idx], y[idx], method = "spearman"))
  })
  ci <- quantile(boot[is.finite(boot)], c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  c(rho = estimate, lower = ci[1], upper = ci[2], n = length(x), blocks = length(blocks))
}

outcomes <- c("background_mean", "inside_mean", "enhancement_mean")
outcome_labels <- c(
  background_mean = "South background (s)",
  inside_mean = "Internal spatial mean",
  enhancement_mean = "Inside - background"
)
correlations <- rbindlist(lapply(unique(hourly$period), function(p) {
  rbindlist(lapply(unique(hourly$gas), function(g) {
    d <- hourly[period == p & gas == g]
    rbindlist(lapply(outcomes, function(outcome) {
      z <- block_boot_spearman(d$wind_speed_ms, d[[outcome]], d$date_utc)
      data.table(period = p, gas = g, outcome = outcome,
                 rho = z["rho"], ci_lower = z["lower"], ci_upper = z["upper"],
                 n_hours = as.integer(z["n"]), n_dates = as.integer(z["blocks"]))
    }))
  }))
}))

# Lag definition: positive lag means that wind precedes the concentration
# outcome. Lags are exploratory and are not interpreted causally.
lag_correlations <- rbindlist(lapply(unique(hourly$period), function(p) {
  rbindlist(lapply(unique(hourly$gas), function(g) {
    d <- copy(hourly[period == p & gas == g])
    rbindlist(lapply(outcomes, function(outcome) {
      rbindlist(lapply(-6:6, function(lag_h) {
        shifted <- d[, .(hour_utc = hour_utc + hours(lag_h), wind_lagged = wind_speed_ms)]
        m <- merge(d[, c("hour_utc", outcome), with = FALSE], shifted, by = "hour_utc")
        data.table(period = p, gas = g, outcome = outcome, lag_hours = lag_h,
                   rho = suppressWarnings(cor(m$wind_lagged, m[[outcome]], method = "spearman")),
                   n_hours = nrow(m))
      }))
    }))
  }))
}))
lag_correlations[, outcome_label := factor(outcome_labels[outcome], levels = outcome_labels)]

sector_summary <- hourly[, .(
  n_hours = .N,
  background_mean = mean(background_mean, na.rm = TRUE),
  inside_mean = mean(inside_mean, na.rm = TRUE),
  enhancement_mean = mean(enhancement_mean, na.rm = TRUE)
), by = .(period, gas, wind_sector)][order(period, gas, wind_sector)]

fwrite(hourly, file.path(output_dir, "campaign1_hourly_gas_background_dwd_aligned.csv"))
fwrite(correlations, file.path(output_dir, "campaign1_dwd_wind_spearman_block_bootstrap.csv"))
fwrite(lag_correlations, file.path(output_dir, "campaign1_dwd_wind_lag_correlations.csv"))
fwrite(sector_summary, file.path(output_dir, "campaign1_dwd_wind_sector_summary.csv"))

correlations[, outcome_label := factor(outcome_labels[outcome], levels = outcome_labels)]

p1 <- ggplot(correlations, aes(rho, gas, colour = outcome_label)) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "grey55") +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), width = 0.18,
                orientation = "y", position = position_dodge(width = 0.55),
                linewidth = 0.55) +
  geom_point(position = position_dodge(width = 0.55), size = 2.4) +
  facet_wrap(~period, ncol = 1) +
  scale_colour_brewer(palette = "Dark2", name = NULL) +
  coord_cartesian(xlim = c(-1, 1)) +
  labs(x = "Spearman correlation with regional hourly wind speed (95% date-block bootstrap CI)",
       y = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

best_lags <- lag_correlations[, .SD[which.max(abs(rho))], by = .(period, gas, outcome)]
p2 <- ggplot(lag_correlations, aes(lag_hours, rho, colour = outcome_label)) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey55") +
  geom_line(linewidth = 0.65) +
  geom_point(data = best_lags, size = 1.8) +
  facet_grid(gas ~ period) +
  scale_colour_brewer(palette = "Dark2", name = NULL) +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(x = "Lag (h; positive = wind precedes concentration)",
       y = "Spearman correlation") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

for (dir in c(plot_dir, manuscript_figure_dir)) {
  ggsave(file.path(dir, "campaign1_dwd_wind_background_correlations.pdf"), p1,
         width = 7.2, height = 5.2, device = cairo_pdf)
  ggsave(file.path(dir, "campaign1_dwd_wind_lag_correlations.pdf"), p2,
         width = 7.4, height = 6.0, device = cairo_pdf)
}

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo_external_wind.txt"))

cat("Aligned hourly rows:", nrow(hourly), "\n")
print(correlations)
cat("\nLargest absolute exploratory lag correlations:\n")
print(best_lags[order(period, gas, outcome)])
