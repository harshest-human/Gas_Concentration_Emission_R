library(dplyr)
library(tidyr)
library(readr)
library(lubridate)
library(ggplot2)
library(purrr)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args_all[grepl(file_arg, args_all)])
helper_dir <- if (length(script_path) > 0) dirname(normalizePath(script_path[1], winslash = "/", mustWork = FALSE)) else getwd()
source(file.path(helper_dir, "device_comparison_helpers.R"))

calc_metrics <- function(df, ref_col, cmp_col) {
  dat <- df %>%
    select(all_of(c(ref_col, cmp_col))) %>%
    rename(ref = all_of(ref_col), cmp = all_of(cmp_col)) %>%
    filter(!is.na(ref), !is.na(cmp))

  if (nrow(dat) < 2) {
    return(tibble(
      n_pairs = nrow(dat),
      mean_ref = NA_real_,
      mean_cmp = NA_real_,
      rd_percent = NA_real_,
      bias = NA_real_,
      mae = NA_real_,
      rmse = NA_real_,
      slope = NA_real_,
      intercept = NA_real_,
      r2 = NA_real_
    ))
  }

  fit <- lm(cmp ~ ref, data = dat)
  coef_fit <- coef(fit)

  tibble(
    n_pairs = nrow(dat),
    mean_ref = mean(dat$ref, na.rm = TRUE),
    mean_cmp = mean(dat$cmp, na.rm = TRUE),
    rd_percent = 100 * (mean(dat$cmp, na.rm = TRUE) - mean(dat$ref, na.rm = TRUE)) / mean(dat$ref, na.rm = TRUE),
    bias = mean(dat$cmp - dat$ref, na.rm = TRUE),
    mae = safe_mae(dat$ref, dat$cmp),
    rmse = safe_rmse(dat$ref, dat$cmp),
    slope = unname(coef_fit[["ref"]]),
    intercept = unname(coef_fit[["(Intercept)"]]),
    r2 = safe_r2(dat$ref, dat$cmp)
  )
}

label_gas <- function(x) {
  recode(
    x,
    delta_CO2 = "Delta CO2",
    delta_CH4 = "Delta CH4",
    delta_NH3 = "Delta NH3"
  )
}

build_joined <- function(device_delta, analyzer_id, start_time = NULL, end_time = NULL) {
  ref_df <- device_delta %>%
    filter(analyzer == "crds") %>%
    select(DATE.HOUR, delta_CO2, delta_CH4, delta_NH3)

  cmp_df <- device_delta %>%
    filter(analyzer == analyzer_id) %>%
    select(DATE.HOUR, delta_CO2, delta_CH4, delta_NH3)

  joined <- inner_join(ref_df, cmp_df, by = "DATE.HOUR", suffix = c("_ref", "_cmp"))

  if (!is.null(start_time)) {
    joined <- joined %>% filter(DATE.HOUR >= start_time)
  }
  if (!is.null(end_time)) {
    joined <- joined %>% filter(DATE.HOUR <= end_time)
  }

  joined
}

build_metric_table <- function(device_delta, analyzers, start_time = NULL, end_time = NULL, window_label = "campaign") {
  gases <- c("delta_CO2", "delta_CH4", "delta_NH3")

  bind_rows(lapply(analyzers, function(analyzer_id) {
    joined <- build_joined(device_delta, analyzer_id, start_time, end_time)

    bind_rows(lapply(gases, function(gas_id) {
      ref_col <- paste0(gas_id, "_ref")
      cmp_col <- paste0(gas_id, "_cmp")

      if (!cmp_col %in% names(joined)) {
        return(NULL)
      }

      calc_metrics(joined, ref_col, cmp_col) %>%
        mutate(
          window = window_label,
          analyzer = label_analyzer(analyzer_id),
          gas = label_gas(gas_id)
        )
    }))
  })) %>%
    select(window, analyzer, gas, n_pairs, mean_ref, mean_cmp, rd_percent, bias, mae, rmse, slope, intercept, r2)
}

build_overlap_summary <- function(device_delta, analyzers, start_time = NULL, end_time = NULL, window_label = "campaign") {
  gases <- c("delta_CO2", "delta_CH4", "delta_NH3")

  bind_rows(lapply(analyzers, function(analyzer_id) {
    joined <- build_joined(device_delta, analyzer_id, start_time, end_time)

    bind_rows(lapply(gases, function(gas_id) {
      ref_col <- paste0(gas_id, "_ref")
      cmp_col <- paste0(gas_id, "_cmp")
      if (!all(c(ref_col, cmp_col) %in% names(joined))) {
        return(NULL)
      }
      tibble(
        window = window_label,
        analyzer = label_analyzer(analyzer_id),
        gas = label_gas(gas_id),
        n_pairs = sum(!is.na(joined[[ref_col]]) & !is.na(joined[[cmp_col]])),
        first_hour = suppressWarnings(min(joined$DATE.HOUR[!is.na(joined[[ref_col]]) & !is.na(joined[[cmp_col]])], na.rm = TRUE)),
        last_hour = suppressWarnings(max(joined$DATE.HOUR[!is.na(joined[[ref_col]]) & !is.na(joined[[cmp_col]])], na.rm = TRUE))
      ) %>%
        mutate(
          first_hour = ifelse(is.infinite(first_hour), NA, as.character(first_hour)),
          last_hour = ifelse(is.infinite(last_hour), NA, as.character(last_hour))
        )
    }))
  }))
}

build_scatter_plot <- function(device_delta, analyzers, start_time, end_time) {
  plot_df <- bind_rows(lapply(analyzers, function(analyzer_id) {
    joined <- build_joined(device_delta, analyzer_id, start_time, end_time)
    joined %>%
      transmute(
        DATE.HOUR,
        analyzer = label_analyzer(analyzer_id),
        `Delta CO2_ref` = delta_CO2_ref,
        `Delta CO2_cmp` = delta_CO2_cmp,
        `Delta CH4_ref` = delta_CH4_ref,
        `Delta CH4_cmp` = delta_CH4_cmp,
        `Delta NH3_ref` = delta_NH3_ref,
        `Delta NH3_cmp` = delta_NH3_cmp
      ) %>%
      pivot_longer(
        cols = -c(DATE.HOUR, analyzer),
        names_to = c("gas", "series"),
        names_sep = "_",
        values_to = "value"
      ) %>%
      pivot_wider(names_from = series, values_from = value) %>%
      filter(!is.na(ref), !is.na(cmp))
  }))

  ggplot(plot_df, aes(x = ref, y = cmp)) +
    geom_abline(slope = 1, intercept = 0, color = "#9CA3AF", linetype = "dashed") +
    geom_point(color = "#005CA9", alpha = 0.35, size = 1.4) +
    geom_smooth(method = "lm", se = FALSE, color = "#99C13C", linewidth = 0.8) +
    facet_grid(gas ~ analyzer, scales = "free") +
    labs(
      title = "December 2025 overlap: analyzer vs CRDS reference",
      subtitle = "Hourly paired delta concentrations during the four-device overlap window",
      x = "CRDS",
      y = "Compared analyzer"
    ) +
    device_plot_theme_minimal()
}

build_rd_plot <- function(metric_df) {
  ggplot(metric_df, aes(x = gas, y = rd_percent, fill = analyzer)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65) +
    geom_hline(yintercept = 0, color = "#374151", linewidth = 0.5) +
    scale_fill_manual(values = c(
      "PRONOVA" = "#6A5ACD",
      "CUBIC" = "#169B62",
      "OTICE" = "#C45A11"
    )) +
    labs(
      title = "Relative difference vs CRDS",
      subtitle = "December 2025 paired hourly overlap",
      x = NULL,
      y = "RD (%)",
      fill = NULL
    ) +
    device_plot_theme_minimal()
}

project_dir <- resolve_project_dir()
workflow_dir <- file.path(project_dir, "workflows", "device_comparison")
table_dir <- file.path(workflow_dir, "result_data", "tables")
plot_dir <- file.path(workflow_dir, "result_data", "plots", "cigr_conference")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

device_delta <- read_csv(file.path(table_dir, "device_delta_campaign_all.csv"), show_col_types = FALSE) %>%
  mutate(
    DATE.HOUR = parse_mixed_datetime(DATE.HOUR),
    month_id = format(DATE.HOUR, "%Y-%m")
  ) %>%
  filter(!is.na(DATE.HOUR))

coverage_summary <- device_delta %>%
  mutate(analyzer = label_analyzer(analyzer)) %>%
  group_by(analyzer) %>%
  summarise(
    first_hour = min(DATE.HOUR, na.rm = TRUE),
    last_hour = max(DATE.HOUR, na.rm = TRUE),
    total_rows = n(),
    .groups = "drop"
  )

dec_start <- as.POSIXct("2025-12-09 12:00:00", tz = "Europe/Berlin")
dec_end <- as.POSIXct("2025-12-22 23:00:00", tz = "Europe/Berlin")

analyzers <- c("logas_ndir", "logas_tdlas", "otice")

metrics_campaign <- build_metric_table(device_delta, analyzers, window_label = "campaign")
metrics_dec <- build_metric_table(device_delta, analyzers, dec_start, dec_end, "dec_2025_overlap")

overlap_campaign <- build_overlap_summary(device_delta, analyzers, window_label = "campaign")
overlap_dec <- build_overlap_summary(device_delta, analyzers, dec_start, dec_end, "dec_2025_overlap")

write_csv(coverage_summary, file.path(table_dir, "cigr_coverage_summary.csv"))
write_csv(metrics_campaign, file.path(table_dir, "cigr_metrics_campaign.csv"))
write_csv(metrics_dec, file.path(table_dir, "cigr_metrics_dec_2025.csv"))
write_csv(bind_rows(overlap_campaign, overlap_dec), file.path(table_dir, "cigr_overlap_summary.csv"))

ggsave(
  file.path(plot_dir, "cigr_scatter_dec_2025.png"),
  build_scatter_plot(device_delta, analyzers, dec_start, dec_end),
  width = 13,
  height = 8,
  dpi = 180
)

ggsave(
  file.path(plot_dir, "cigr_rd_dec_2025.png"),
  build_rd_plot(metrics_dec %>% filter(!is.na(rd_percent), analyzer != "CRDS")),
  width = 10,
  height = 5.8,
  dpi = 180
)

cat("Wrote conference coverage and metric tables to:\n", normalizePath(table_dir, winslash = "/", mustWork = FALSE), "\n")
cat("Wrote conference plots to:\n", normalizePath(plot_dir, winslash = "/", mustWork = FALSE), "\n")
