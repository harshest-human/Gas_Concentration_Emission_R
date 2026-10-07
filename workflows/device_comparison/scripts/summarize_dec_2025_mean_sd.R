suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
})

script_dir <- {
  args_all <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- sub(file_arg, "", args_all[grepl(file_arg, args_all)])
  if (length(script_path) > 0) {
    dirname(normalizePath(script_path[1], winslash = "/", mustWork = FALSE))
  } else {
    normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  }
}

workflow_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
table_dir <- file.path(workflow_dir, "result_data", "tables")
plot_dir <- file.path(workflow_dir, "result_data", "plots", "dec_2025")
input_path <- file.path(table_dir, "emission_data_20251209_20251222.csv")

if (!file.exists(input_path)) stop("Missing December emission data: ", input_path, call. = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

device_order <- c("CRDS", "CUBIC", "PRONOVA", "OTICE")
device_colors <- c(CRDS = "#59636F", CUBIC = "#1FA77D", PRONOVA = "#7774BA", OTICE = "#C65D0E")
device_shapes <- c(CRDS = 2, CUBIC = 0, PRONOVA = 1, OTICE = 15)

parameter_info <- tibble::tribble(
  ~variable,    ~parameter,       ~group,      ~unit,
  "delta_CO2", "Delta CO2",      "delta",     "ppm",
  "delta_CH4", "Delta CH4",      "delta",     "ppm",
  "delta_NH3", "Delta NH3",      "delta",     "ppm",
  "Q_vent",    "Q",              "emission",  "m3 h-1 LU-1",
  "e_CH4_ghLU", "eCH4",          "emission",  "g h-1 LU-1",
  "e_NH3_ghLU", "eNH3",          "emission",  "g h-1 LU-1"
)

errorbar_zero_plot <- function(data, variables, labels) {
  plot_data <- data |>
    select(analyzer, any_of(variables)) |>
    pivot_longer(-analyzer, names_to = "variable", values_to = "value") |>
    filter(is.finite(value)) |>
    group_by(analyzer, variable) |>
    summarise(
      mean = mean(value),
      sd = sd(value),
      n = n(),
      se = sd / sqrt(n),
      .groups = "drop"
    ) |>
    mutate(
      analyzer = factor(analyzer, levels = device_order),
      variable = factor(variable, levels = variables, labels = unname(labels[variables]))
    )

  ggplot(plot_data, aes(x = analyzer, y = mean, colour = analyzer, shape = analyzer)) +
    geom_point(size = 3, stroke = 1) +
    geom_errorbar(aes(ymin = pmax(0, mean - se), ymax = mean + se), width = 0.2, linewidth = 0.65) +
    facet_grid(variable ~ ., scales = "free_y", switch = "y") +
    scale_colour_manual(values = device_colors, drop = FALSE) +
    scale_shape_manual(values = device_shapes, drop = FALSE) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.06))) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      strip.placement = "outside",
      strip.background = element_rect(fill = "white", colour = "black"),
      strip.text.y.left = element_text(angle = 90, size = 11),
      panel.border = element_rect(colour = "black", fill = NA),
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.margin = margin(8, 15, 8, 8)
    ) +
    guides(colour = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))
}

data <- read_csv(input_path, show_col_types = FALSE) |>
  mutate(
    analyzer = recode(
      as.character(analyzer),
      crds = "CRDS", logas_tdlas = "CUBIC", logas_ndir = "PRONOVA",
      otice = "OTICE", `CRDS` = "CRDS", `CUBIC` = "CUBIC",
      `PRONOVA` = "PRONOVA", `OTICE` = "OTICE",
      .default = as.character(analyzer)
    ),
    analyzer = factor(analyzer, levels = device_order)
  )

hourly_long <- data |>
  select(DATE.TIME, analyzer, any_of(parameter_info$variable)) |>
  pivot_longer(any_of(parameter_info$variable), names_to = "variable", values_to = "value")

daily_long <- hourly_long |>
  mutate(date_day = as.Date(DATE.TIME)) |>
  group_by(date_day, analyzer, variable) |>
  summarise(value = if (all(is.na(value))) NA_real_ else mean(value, na.rm = TRUE), .groups = "drop")

summary_long <- daily_long |>
  left_join(parameter_info, by = "variable") |>
  group_by(group, variable, parameter, unit, analyzer) |>
  summarise(
    n = sum(!is.na(value)),
    mean = if_else(n > 0, mean(value, na.rm = TRUE), NA_real_),
    sd = if_else(n > 1, sd(value, na.rm = TRUE), NA_real_),
    .groups = "drop"
  ) |>
  group_by(variable) |>
  mutate(
    crds_mean = mean[analyzer == "CRDS"][1],
    percent_difference_vs_crds = if_else(
      analyzer == "CRDS" | is.na(mean) | is.na(crds_mean) | crds_mean == 0,
      if_else(analyzer == "CRDS", 0, NA_real_),
      100 * (mean - crds_mean) / crds_mean
    )
  ) |>
  ungroup() |>
  arrange(match(variable, parameter_info$variable), analyzer)

write_csv(summary_long, file.path(table_dir, "dec_2025_device_mean_sd_long.csv"), na = "NA")

format_number <- function(x, digits) {
  mapply(function(value, places) {
    if (is.na(value)) "NA" else formatC(value, format = "f", digits = places, big.mark = "")
  }, x, digits, USE.NAMES = FALSE)
}

presentation <- summary_long |>
  mutate(
    digits = case_when(variable == "delta_CO2" ~ 1L, variable == "Q_vent" ~ 1L, TRUE ~ 2L),
    value_text = if_else(
      is.na(mean),
      "NA",
      paste0(
        format_number(mean, digits), " +/- ", format_number(sd, digits),
        if_else(analyzer == "CRDS", "", paste0("\n(", sprintf("%+.1f", percent_difference_vs_crds), "%)"))
      )
    )
  ) |>
  select(parameter, analyzer, value_text) |>
  pivot_wider(names_from = analyzer, values_from = value_text) |>
  arrange(match(parameter, parameter_info$parameter))

write_csv(presentation, file.path(table_dir, "dec_2025_device_mean_sd_presentation.csv"), na = "NA")

moving_summary <- function(df, hours = 24, min_observations = 6) {
  df <- arrange(df, DATE.TIME)
  half_window <- hours * 3600 / 2
  timestamp <- as.numeric(df$DATE.TIME)
  stats <- lapply(seq_len(nrow(df)), function(i) {
    values <- df$value[abs(timestamp - timestamp[i]) <= half_window]
    values <- values[is.finite(values)]
    if (length(values) < min_observations) return(c(mean = NA_real_, sd = NA_real_, n = length(values)))
    c(mean = mean(values), sd = if (length(values) > 1) sd(values) else NA_real_, n = length(values))
  })
  bind_cols(df, as_tibble(do.call(rbind, stats)))
}

smooth_time_plot <- function(hourly_df, group_name, title_text) {
  plot_df <- hourly_df |>
    left_join(parameter_info, by = "variable") |>
    filter(group == group_name) |>
    group_by(analyzer, variable) |>
    group_modify(~ moving_summary(.x, hours = 24, min_observations = 6)) |>
    ungroup() |>
    mutate(
      parameter = factor(parameter, levels = parameter_info$parameter[parameter_info$group == group_name]),
      analyzer = factor(analyzer, levels = device_order),
      lower = mean - sd,
      upper = mean + sd
    )

  ggplot(plot_df, aes(x = DATE.TIME, y = mean, colour = analyzer, fill = analyzer)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, colour = NA, na.rm = TRUE) +
    geom_line(linewidth = 0.8, na.rm = TRUE) +
    facet_wrap(~parameter, scales = "free_y", ncol = 1, strip.position = "left") +
    scale_colour_manual(values = device_colors, drop = FALSE) +
    scale_fill_manual(values = device_colors, drop = FALSE) +
    scale_x_datetime(date_breaks = "2 days", date_labels = "%d-%m-%Y") +
    labs(title = title_text, subtitle = "Line: centred 24-hour moving mean; shaded band: moving mean +/- one SD", x = NULL, y = NULL) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
      plot.subtitle = element_text(hjust = 0.5, colour = "grey35"),
      strip.placement = "outside",
      strip.background = element_rect(fill = "white", colour = "black"),
      strip.text.y.left = element_text(angle = 90, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      panel.border = element_rect(colour = "black", fill = NA),
      plot.margin = margin(10, 28, 10, 10)
    ) +
    guides(colour = guide_legend(title = NULL), fill = "none")
}

render_summary_table <- function(table_df, parameters, title_text, output_path) {
  display <- table_df |>
    filter(parameter %in% parameters) |>
    mutate(parameter = factor(parameter, levels = parameters)) |>
    arrange(parameter) |>
    mutate(parameter = as.character(parameter))
  display <- display[, c("parameter", device_order)]
  names(display)[1] <- "Parameter"

  long <- display |>
    mutate(row = row_number()) |>
    pivot_longer(-row, names_to = "column", values_to = "label") |>
    mutate(
      col = match(column, names(display)),
      x = col - 0.5,
      y = nrow(display) - row + 0.5
    )
  headers <- tibble(column = names(display), col = seq_along(display), x = col - 0.5, y = nrow(display) + 0.5)

  p <- ggplot() +
    geom_rect(aes(xmin = 0, xmax = ncol(display), ymin = nrow(display), ymax = nrow(display) + 1), fill = "#F4F4F4", colour = NA) +
    geom_vline(xintercept = 0:ncol(display), colour = "grey35", linewidth = 0.5) +
    geom_hline(yintercept = 0:(nrow(display) + 1), colour = "grey35", linewidth = 0.5) +
    geom_text(data = headers, aes(x = x, y = y, label = column), fontface = "bold", size = 4.2) +
    geom_text(data = filter(long, column == "Parameter"), aes(x = x, y = y, label = label), fontface = "bold", lineheight = 1.15, size = 3.8) +
    geom_text(data = filter(long, column != "Parameter"), aes(x = x, y = y, label = label), lineheight = 1.15, size = 3.8) +
    coord_cartesian(xlim = c(0, ncol(display)), ylim = c(0, nrow(display) + 1), expand = FALSE) +
    labs(title = title_text, subtitle = "Mean +/- SD; parentheses show difference from CRDS mean") +
    theme_void() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15), plot.subtitle = element_text(hjust = 0.5, colour = "grey35", margin = margin(b = 8)))

  ggsave(output_path, p, width = 9, height = 1.15 + 1.25 * nrow(display), dpi = 200, bg = "white")
}

delta_plot <- smooth_time_plot(hourly_long, "delta", "Smoothed device delta concentrations, 09-12-2025 to 22-12-2025")
emission_plot <- smooth_time_plot(hourly_long, "emission", "Smoothed ventilation and emissions, 09-12-2025 to 22-12-2025")

delta_zero_plot <- errorbar_zero_plot(
  data,
  c("delta_CO2", "delta_CH4", "delta_NH3"),
  c(delta_CO2 = "Delta CO2 (ppm)", delta_CH4 = "Delta CH4 (ppm)", delta_NH3 = "Delta NH3 (ppm)")
)
emission_zero_plot <- errorbar_zero_plot(
  data,
  c("Q_vent", "e_CH4_ghLU", "e_NH3_ghLU"),
  c(Q_vent = "Q (m3 h-1 LU-1)", e_CH4_ghLU = "eCH4 (g h-1 LU-1)", e_NH3_ghLU = "eNH3 (g h-1 LU-1)")
)

ggsave(file.path(plot_dir, "device_delta_smoothed_mean_sd_20251209_20251222.png"), delta_plot, width = 11, height = 7.5, dpi = 200, bg = "white")
ggsave(file.path(plot_dir, "device_emission_smoothed_mean_sd_20251209_20251222.png"), emission_plot, width = 11, height = 7.5, dpi = 200, bg = "white")
ggsave(file.path(plot_dir, "device_delta_errorbars_zero_y_20251209_20251222.png"), delta_zero_plot, width = 7.5, height = 6, dpi = 200, bg = "white")
ggsave(file.path(plot_dir, "device_emission_errorbars_zero_y_20251209_20251222.png"), emission_zero_plot, width = 7.5, height = 6, dpi = 200, bg = "white")

render_summary_table(
  presentation,
  c("Delta CO2", "Delta CH4", "Delta NH3"),
  "Delta concentration summary, 09-12-2025 to 22-12-2025",
  file.path(plot_dir, "device_delta_mean_sd_table_20251209_20251222.png")
)
render_summary_table(
  presentation,
  c("Q", "eCH4", "eNH3"),
  "Ventilation and emission summary, 09-12-2025 to 22-12-2025",
  file.path(plot_dir, "device_emission_mean_sd_table_20251209_20251222.png")
)

cat("Wrote December mean/SD summaries to:\n", normalizePath(plot_dir, winslash = "/"), "\n")
print(summary_long)
