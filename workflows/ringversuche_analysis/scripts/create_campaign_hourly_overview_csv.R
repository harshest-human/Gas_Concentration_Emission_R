library(dplyr)
library(tidyr)
library(readr)
library(stringr)

base_dir <- "D:/Data_Analysis_R/Gas_Concentration_Emission_R/workflows/ringversuche_analysis"
tables_dir <- file.path(base_dir, "result_data/tables/Version_11")

input_path <- file.path(tables_dir, "20250408-15_input_combined.csv")
loss_path <- file.path(tables_dir, "pre_hourly_outlier_loss.csv")
out_path <- file.path(tables_dir, "campaign_hourly_overview_table1.csv")

input_combined <- read_csv(input_path, show_col_types = FALSE)
pre_hourly_loss <- read_csv(loss_path, show_col_types = FALSE)

location_map <- c(
  "in" = "Indoor",
  "N" = "OutdoorNE",
  "S" = "OutdoorSW"
)

gas_hourly <- input_combined %>%
  pivot_longer(
    cols = matches("^(CO2|CH4|NH3)_ppm_(N|S|in)$"),
    names_to = c("variable", "location_code"),
    names_pattern = "^(CO2|CH4|NH3)_ppm_(N|S|in)$",
    values_to = "value"
  ) %>%
  mutate(
    source = "gas_concentration",
    location = recode(location_code, !!!location_map)
  ) %>%
  group_by(source, lab, analyser, location, variable) %>%
  summarise(
    expected_hourly_timestamps = n(),
    hourly_observations = sum(is.finite(value)),
    hourly_missing = sum(!is.finite(value)),
    hourly_mean = mean(value, na.rm = TRUE),
    hourly_min = min(value, na.rm = TRUE),
    hourly_max = max(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    pre_hourly_loss %>%
      filter(variable %in% c("CO2", "CH4", "NH3")) %>%
      group_by(lab, analyser, variable) %>%
      summarise(
        raw_7p5min_observations = max(n_before, na.rm = TRUE),
        outliers_removed_before_hourly = max(n_removed, na.rm = TRUE),
        .groups = "drop"
      ),
    by = c("lab", "analyser", "variable")
  ) %>%
  mutate(
    raw_7p5min_observations = coalesce(raw_7p5min_observations, 0),
    outliers_removed_before_hourly = coalesce(outliers_removed_before_hourly, 0),
    raw_resolution_min = 7.5,
    hourly_resolution_min = 60
  ) %>%
  select(
    source, lab, analyser, location, variable,
    raw_7p5min_observations, raw_resolution_min,
    outliers_removed_before_hourly,
    expected_hourly_timestamps, hourly_observations, hourly_missing,
    hourly_resolution_min, hourly_mean, hourly_min, hourly_max
  )

meta_hourly <- bind_rows(
  input_combined %>%
    distinct(DATE.TIME, n_dairycows_in) %>%
    transmute(
      source = "animal",
      lab = "meta",
      analyser = "meta",
      location = "Indoor",
      variable = "n_dairycows_in",
      value = n_dairycows_in
  ),
  input_combined %>%
    distinct(DATE.TIME, m_weight) %>%
    transmute(
      source = "animal",
      lab = "meta",
      analyser = "meta",
      location = "meta",
      variable = "m_weight",
      value = m_weight
  ),
  input_combined %>%
    distinct(DATE.TIME, p_pregnancy_day) %>%
    transmute(
      source = "animal",
      lab = "meta",
      analyser = "meta",
      location = "meta",
      variable = "p_pregnancy_day",
      value = p_pregnancy_day
  ),
  input_combined %>%
    distinct(DATE.TIME, Y1_milk_prod) %>%
    transmute(
      source = "animal",
      lab = "meta",
      analyser = "meta",
      location = "meta",
      variable = "Y1_milk_prod",
      value = Y1_milk_prod
  ),
  input_combined %>%
    distinct(DATE.TIME, temp_N) %>%
    transmute(
      source = "climate",
      lab = "meta",
      analyser = "meta",
      location = "OutdoorNE",
      variable = "temp_N",
      value = temp_N
  ),
  input_combined %>%
    distinct(DATE.TIME, RH_N) %>%
    transmute(
      source = "climate",
      lab = "meta",
      analyser = "meta",
      location = "OutdoorNE",
      variable = "RH_N",
      value = RH_N
  ),
  input_combined %>%
    distinct(DATE.TIME, temp_in) %>%
    transmute(
      source = "climate",
      lab = "meta",
      analyser = "meta",
      location = "Indoor",
      variable = "temp_in",
      value = temp_in
  ),
  input_combined %>%
    distinct(DATE.TIME, RH_in) %>%
    transmute(
      source = "climate",
      lab = "meta",
      analyser = "meta",
      location = "Indoor",
      variable = "RH_in",
      value = RH_in
  ),
  input_combined %>%
    distinct(DATE.TIME, wd_mst) %>%
    transmute(
      source = "wind",
      lab = "meta",
      analyser = "Ultrasonic_Anemometer",
      location = "Mast",
      variable = "wd_mst",
      value = wd_mst
  ),
  input_combined %>%
    distinct(DATE.TIME, ws_mst) %>%
    transmute(
      source = "wind",
      lab = "meta",
      analyser = "Ultrasonic_Anemometer",
      location = "Mast",
      variable = "ws_mst",
      value = ws_mst
    )
) %>%
  group_by(source, lab, analyser, location, variable) %>%
  summarise(
    raw_7p5min_observations = n(),
    raw_resolution_min = 60,
    outliers_removed_before_hourly = 0,
    expected_hourly_timestamps = n(),
    hourly_observations = sum(is.finite(value)),
    hourly_missing = sum(!is.finite(value)),
    hourly_resolution_min = 60,
    hourly_mean = mean(value, na.rm = TRUE),
    hourly_min = min(value, na.rm = TRUE),
    hourly_max = max(value, na.rm = TRUE),
    .groups = "drop"
  )

overview_table <- bind_rows(gas_hourly, meta_hourly) %>%
  mutate(
    across(c(hourly_mean, hourly_min, hourly_max), ~ round(.x, 3))
  ) %>%
  arrange(source, lab, analyser, location, variable)

write_csv(overview_table, out_path, na = "")
message("Wrote: ", out_path)
