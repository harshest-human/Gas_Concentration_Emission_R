##### DWD regional meteorology for Manuscript 3 ##############################

# Purpose:
#   1. Read hourly DWD historical observations from Potsdam station 03987.
#   2. Retain the 2024 periods relevant to Campaigns 1 and 2.
#   3. Produce study-period wind roses.
#   4. Produce an integrated temperature, relative-humidity, and rainfall plot.
#
# Scientific boundary:
# DWD Potsdam describes the regional meteorological setting approximately
# 18.5 km from the barn. It is not a substitute for wind measured at the
# on-site 10 m mast and must not be described as barn-inlet airflow.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(lubridate)
  library(patchwork)
})

workflow <- "D:/Data_Analysis_R/Gas_Concentration_Emission_R/workflows/mapping_high_resolution"
dwd_root <- file.path(workflow, "raw_data", "DWD_Potsdam_03987")
clean_dir <- file.path(workflow, "clean_data", "manuscript3_weather")
plot_dir <- file.path(workflow, "plots", "manuscript3_weather")
figure_dir <- file.path(
  workflow, "Manuscript_3_Mapping_high_resolution_concentration_Latex_draft",
  "figures"
)
dir.create(clean_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

find_product <- function(pattern) {
  files <- list.files(dwd_root, pattern = pattern, recursive = TRUE,
                      full.names = TRUE)
  if (length(files) != 1L) {
    stop("Expected one DWD product file for pattern ", pattern,
         "; found ", length(files))
  }
  files
}

read_dwd <- function(path) {
  x <- fread(path, sep = ";", strip.white = TRUE, na.strings = "-999")
  setnames(x, trimws(names(x)))
  if ("eor" %in% names(x)) x[, eor := NULL]
  x
}

wind <- read_dwd(find_product("^produkt_ff_stunde_.*\\.txt$"))
thermo <- read_dwd(find_product("^produkt_tf_stunde_.*\\.txt$"))
rain <- read_dwd(find_product("^produkt_rr_stunde_.*\\.txt$"))

parse_dwd_time <- function(x) {
  # DWD hourly climate observations use UTC timestamps.
  with_tz(ymd_h(as.character(x), tz = "UTC"), "Europe/Berlin")
}

wind[, DATE.TIME := parse_dwd_time(MESS_DATUM)]
thermo[, DATE.TIME := parse_dwd_time(MESS_DATUM)]
rain[, DATE.TIME := parse_dwd_time(MESS_DATUM)]

weather <- merge(
  wind[, .(DATE.TIME, wind_speed_ms = F, wind_direction_deg = D)],
  thermo[, .(DATE.TIME, temperature_C = TT_STD,
             relative_humidity_pct = RF_STD)],
  by = "DATE.TIME", all = TRUE
)
weather <- merge(
  weather,
  rain[, .(DATE.TIME, precipitation_mm = R1, rain_indicator = RS_IND)],
  by = "DATE.TIME", all = TRUE
)

weather <- weather[
  DATE.TIME >= as.POSIXct("2024-06-01 00:00:00", tz = "Europe/Berlin") &
    DATE.TIME <= as.POSIXct("2024-12-31 23:59:59", tz = "Europe/Berlin")
]

weather[, study_period := fcase(
  DATE.TIME >= as.POSIXct("2024-06-01 00:00:00", tz = "Europe/Berlin") &
    DATE.TIME <= as.POSIXct("2024-08-31 23:59:59", tz = "Europe/Berlin"),
  "Campaign 1: summer",
  DATE.TIME >= as.POSIXct("2024-10-01 00:00:00", tz = "Europe/Berlin") &
    DATE.TIME <= as.POSIXct("2024-10-24 23:59:59", tz = "Europe/Berlin"),
  "Campaign 1: autumn",
  DATE.TIME >= as.POSIXct("2024-11-16 00:00:00", tz = "Europe/Berlin") &
    DATE.TIME <= as.POSIXct("2024-12-31 23:59:59", tz = "Europe/Berlin"),
  "Campaign 2: late autumn/winter",
  default = NA_character_
)]

period_levels <- c(
  "Campaign 1: summer",
  "Campaign 1: autumn",
  "Campaign 2: late autumn/winter"
)
weather[, study_period := factor(study_period, levels = period_levels)]

fwrite(weather, file.path(clean_dir,
                          "dwd_potsdam_03987_hourly_2024_campaign_periods.csv"))

##### Figure 1: seasonal study-period wind roses #############################

wind_plot <- weather[
  !is.na(study_period) & !is.na(wind_speed_ms) &
    !is.na(wind_direction_deg) & wind_direction_deg >= 0
]

# Sixteen 22.5-degree directional sectors centred on the compass directions.
wind_plot[, direction_sector :=
            (floor(((wind_direction_deg + 11.25) %% 360) / 22.5) * 22.5)]
wind_plot[, speed_class := cut(
  wind_speed_ms,
  breaks = c(-Inf, 1, 2, 3, 5, 7, Inf),
  labels = c("<1", "1--<2", "2--<3", "3--<5", "5--<7", "≥7"),
  right = FALSE
)]

rose <- wind_plot[, .N, by = .(study_period, direction_sector, speed_class)]
rose[, frequency_pct := 100 * N / sum(N), by = study_period]

speed_colours <- c(
  "<1" = "#deebf7", "1--<2" = "#9ecae1", "2--<3" = "#6baed6",
  "3--<5" = "#3182bd", "5--<7" = "#08519c", "≥7" = "#08306b"
)

p_rose <- ggplot(rose, aes(x = direction_sector, y = frequency_pct,
                            fill = speed_class)) +
  geom_col(width = 22.5, colour = "grey30", linewidth = 0.12) +
  coord_polar(start = -pi / 16, direction = 1) +
  facet_wrap(~study_period, nrow = 1) +
  scale_x_continuous(
    breaks = seq(0, 315, 45),
    labels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW"),
    limits = c(-11.25, 348.75)
  ) +
  scale_fill_manual(values = speed_colours, drop = FALSE,
                    name = expression("Wind speed (m "*s^-1*")")) +
  labs(y = "Frequency (%)", x = NULL,
       caption = paste0(
         "DWD Potsdam station 03987 (regional context; approximately 18.5 km ",
         "from the barn). Direction denotes the meteorological direction from ",
         "which the wind originated."
       )) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey94"),
    legend.position = "bottom",
    plot.caption = element_text(hjust = 0, size = 8)
  )

rose_file <- file.path(plot_dir, "Fig_DWD_seasonal_wind_roses.png")
ggsave(rose_file, p_rose, width = 12.5, height = 5.1, dpi = 250,
       bg = "white")
file.copy(rose_file, file.path(figure_dir,
                              "Fig_DWD_seasonal_wind_roses.png"),
          overwrite = TRUE)

##### Figure 2: temperature, humidity, and precipitation #####################

daily <- weather[!is.na(study_period), .(
  temperature_mean_C = mean(temperature_C, na.rm = TRUE),
  temperature_min_C = min(temperature_C, na.rm = TRUE),
  temperature_max_C = max(temperature_C, na.rm = TRUE),
  humidity_mean_pct = mean(relative_humidity_pct, na.rm = TRUE),
  precipitation_sum_mm = sum(precipitation_mm, na.rm = TRUE)
), by = .(date = as.Date(DATE.TIME), study_period)]

period_fill <- c(
  "Campaign 1: summer" = "#f4a261",
  "Campaign 1: autumn" = "#bc6c25",
  "Campaign 2: late autumn/winter" = "#457b9d"
)

p_temp <- ggplot(daily, aes(date, temperature_mean_C,
                            colour = study_period, fill = study_period)) +
  geom_ribbon(aes(ymin = temperature_min_C, ymax = temperature_max_C),
              alpha = 0.16, colour = NA) +
  geom_line(linewidth = 0.55) +
  scale_colour_manual(values = period_fill, guide = "none") +
  scale_fill_manual(values = period_fill, guide = "none") +
  labs(x = NULL, y = "Air temperature (°C)") +
  theme_bw(base_size = 10) + theme(legend.position = "none")

p_rh <- ggplot(daily, aes(date, humidity_mean_pct, colour = study_period)) +
  geom_line(linewidth = 0.55) +
  scale_colour_manual(values = period_fill, guide = "none") +
  coord_cartesian(ylim = c(35, 100)) +
  labs(x = NULL, y = "Relative humidity (%)") +
  theme_bw(base_size = 10) + theme(legend.position = "none")

p_rain <- ggplot(daily, aes(date, precipitation_sum_mm, fill = study_period)) +
  geom_col(width = 0.85) +
  scale_fill_manual(values = period_fill, name = "Study period") +
  labs(x = "Date in 2024", y = "Precipitation (mm d⁻¹)",
       caption = paste0(
         "Daily summaries from DWD Potsdam station 03987. Shaded temperature ",
         "bands show the daily minimum--maximum range. Missing intervals between ",
         "study periods are intentionally retained."
       )) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        plot.caption = element_text(hjust = 0, size = 8))

p_weather <- p_temp / p_rh / p_rain +
  plot_layout(heights = c(1, 1, 0.9))

weather_file <- file.path(
  plot_dir, "Fig_DWD_temperature_humidity_precipitation.png"
)
ggsave(weather_file, p_weather, width = 11.2, height = 8.0, dpi = 250,
       bg = "white")
file.copy(weather_file, file.path(
  figure_dir, "Fig_DWD_temperature_humidity_precipitation.png"
), overwrite = TRUE)

##### Reproducibility summary ################################################

summary_table <- weather[!is.na(study_period), .(
  first_time = min(DATE.TIME),
  last_time = max(DATE.TIME),
  hourly_rows = .N,
  valid_wind = sum(!is.na(wind_speed_ms) & !is.na(wind_direction_deg)),
  valid_temperature = sum(!is.na(temperature_C)),
  valid_humidity = sum(!is.na(relative_humidity_pct)),
  precipitation_total_mm = sum(precipitation_mm, na.rm = TRUE)
), by = study_period]
fwrite(summary_table, file.path(clean_dir,
                                "dwd_potsdam_03987_period_summary.csv"))

writeLines(capture.output(sessionInfo()),
           file.path(clean_dir, "sessionInfo_dwd_weather.txt"))

message("Saved DWD weather data and manuscript figures.")
