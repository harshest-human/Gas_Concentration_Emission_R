# OTICE v5 node comparison: combine and visualise raw sensor data

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(lubridate)
  library(purrr)
  library(scales)
})

tz_local <- "Europe/Berlin"
selected_cols <- c(
  "bridge_time", "sensor_time", "nodeId", "nh3", "nh3T", "nh3RH",
  "co2T", "co2RH", "co2Raw"
)

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
workflow_dir <- if (length(script_arg)) {
  normalizePath(file.path(dirname(script_arg[1]), ".."), winslash = "/", mustWork = TRUE)
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

raw_dir <- file.path(
  workflow_dir, "raw_data", "otice_raw", "2026",
  "OTICE_v5_IOT_bridge_dashboard"
)
clean_dir <- file.path(workflow_dir, "clean_data", "otice_v5_iot_bridge_dashboard")
plot_dir <- file.path(workflow_dir, "plots", "otice_v5_iot_bridge_dashboard")
report_dir <- file.path(workflow_dir, "result_data", "reports")
walk(c(clean_dir, plot_dir, report_dir), dir.create, recursive = TRUE, showWarnings = FALSE)

files <- list.files(raw_dir, pattern = "\\.csv$", full.names = TRUE, ignore.case = TRUE)
if (!length(files)) stop("No CSV files found in: ", raw_dir)

read_node <- function(path) {
  x <- read_csv(path, show_col_types = FALSE, na = c("", "NA", "NaN"))
  missing_cols <- setdiff(selected_cols, names(x))
  if (length(missing_cols)) stop(basename(path), " is missing: ", paste(missing_cols, collapse = ", "))

  x <- x |>
    select(all_of(selected_cols)) |>
    mutate(
      bridge_time = ymd_hms(bridge_time, tz = tz_local, quiet = TRUE),
      sensor_time = ymd_hms(sensor_time, tz = tz_local, quiet = TRUE),
      nodeId = as.integer(nodeId),
      across(c(nh3, nh3T, nh3RH, co2T, co2RH, co2Raw), as.numeric)
    )

  # The filename is authoritative for a physically absent CO2 sensor.
  if (grepl("no_CO2_Sensor", basename(path), ignore.case = TRUE)) {
    x <- x |> mutate(across(c(co2T, co2RH, co2Raw), ~ NA_real_))
  }
  x
}

all_nodes <- map_dfr(files, read_node) |>
  filter(!is.na(bridge_time), !is.na(nodeId)) |>
  arrange(bridge_time, nodeId)

write_csv(all_nodes, file.path(clean_dir, "otice_all_nodes_selected.csv"), na = "")
saveRDS(all_nodes, file.path(clean_dir, "otice_all_nodes_selected.rds"))

long_raw <- all_nodes |>
  select(bridge_time, nodeId, NH3 = nh3, CO2 = co2Raw) |>
  pivot_longer(c(NH3, CO2), names_to = "gas", values_to = "value") |>
  mutate(gas = factor(gas, levels = c("NH3", "CO2")))

node_palette <- hue_pal()(n_distinct(all_nodes$nodeId))
names(node_palette) <- sort(unique(as.character(all_nodes$nodeId)))

raw_plot <- ggplot(long_raw, aes(bridge_time, value, colour = factor(nodeId), group = nodeId)) +
  geom_line(linewidth = 0.7, alpha = 0.85, na.rm = TRUE) +
  facet_grid(. ~ gas, scales = "free_y", labeller = as_labeller(c(NH3 = "NH3 (raw)", CO2 = "CO2 raw"))) +
  scale_colour_manual(values = node_palette) +
  scale_x_datetime(date_labels = "%d %b\n%H:%M", date_breaks = "6 hours", expand = expansion(mult = c(0.01, 0.02))) +
  labs(x = "Bridge time (Europe/Berlin)", y = "Raw sensor value", colour = "Node ID",
       title = "Raw OTICE sensor values over time") +
  theme_bw(base_size = 16) +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 17),
    strip.text = element_text(size = 16),
    plot.title = element_text(size = 20),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 15),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(plot_dir, "otice_raw_nh3_co2.png"), raw_plot, width = 13, height = 5.5, dpi = 300)
ggsave(file.path(plot_dir, "otice_raw_nh3_co2.pdf"), raw_plot, width = 13, height = 5.5)

fmt_value <- function(x) ifelse(is.na(x), "--", formatC(x, digits = 2, format = "f"))

sample_display <- all_nodes |>
  group_by(nodeId) |>
  slice_min(bridge_time, n = 1, with_ties = FALSE) |>
  ungroup() |>
  arrange(nodeId) |>
  transmute(
    sensor_time = format(sensor_time, "%Y-%m-%d %H:%M"),
    nodeId = as.character(nodeId),
    nh3 = fmt_value(nh3),
    co2Raw = fmt_value(co2Raw)
  )

sample_tex_rows <- apply(sample_display, 1, paste, collapse = " & ") |>
  paste0(" \\\\")

report_lines <- c(
  "\\documentclass[11pt]{article}",
  "\\usepackage[a4paper,margin=2.2cm]{geometry}",
  "\\usepackage{graphicx}",
  "\\usepackage{float}",
  "\\usepackage{url}",
  "\\title{Short Test Report: OTICE v5 Sensor-Node Comparison}",
  "\\author{Automated R workflow}",
  "\\date{\\today}",
  "\\begin{document}",
  "\\maketitle",
  "\\section{Objective and method}",
  paste0("Measurements from ", length(files), " OTICE node files were combined row-wise. The retained fields were bridge and sensor time, node ID, NH3, NH3 temperature and relative humidity, and CO2 temperature, relative humidity and raw signal. Bridge time was interpreted in Europe/Berlin time. Raw NH3 and CO2 values were plotted using independent y-axis scales."),
  "\\subsection{Laboratory span-test setup}",
  "The accompanying laboratory note identifies the instruments as EE894 CO$_2$ sensors and TB600 NH$_3$ sensors, tested using a HovaCAL\\textsuperscript{\\textregistered} N 422-SP gas-mixture system. Seven new LoRa-enabled OTICE v5 nodes were installed in an acrylic chamber measuring $50\\times50\\times50$ cm (approximately 125 L). One legacy OTICE v4 node without LoRa was also placed in the chamber. An additional EE894 CO$_2$ sensor, of the same model used with the OTICE v5 nodes, was connected to this OTICE v4 node. No additional NH$_3$ sensor was installed on the OTICE v4 node. The HovaCAL evaporator was set to 150$^{\\circ}$C, the heated line to 50$^{\\circ}$C, and the total flow set point to 5000 mL min$^{-1}$. Source-gas concentrations were 500 ppm NH$_3$ and 10\\% CO$_2$ (100,000 ppm); chamber set points were 12 ppm NH$_3$ and 2,000 ppm CO$_2$, with 1.2\\% water during the wet phase.",
  "\\begin{figure}[H]\\centering",
  "\\includegraphics[width=0.78\\textwidth]{../../meta_data/photos/hovacal_setup.jpeg}",
  "\\caption{HovaCAL N 422-SP gas-mixture system beside the acrylic chamber containing the LoRa-enabled OTICE v5 nodes and the legacy OTICE v4 node without LoRa.}",
  "\\end{figure}",
  "\\subsection{Test sequence}",
  "The chamber was first flushed with dry synthetic air for 10 minutes, from 12:50 to 13:00. During flushing, water and both analyte-gas valves were closed and the synthetic-air valve was fully open. For chamber filling, water was set to 1.2\\%, the NH$_3$ and CO$_2$ valves were opened, and the total flow remained 5000 mL min$^{-1}$. The chamber outlet was left open from 13:00 to 13:05, then closed while the chamber filled from 13:05 to 13:35 (approximately 25--30 minutes).",
  "\\begin{figure}[H]\\centering",
  "\\includegraphics[width=0.88\\textwidth]{../../meta_data/photos/screen_image.jpeg}",
  "\\caption{HovaCAL operating screen during the wet span test, showing the total-flow, water, NH$_3$, CO$_2$, evaporator, and heated-line settings.}",
  "\\end{figure}",
  "\\begin{figure}[H]\\centering",
  "\\includegraphics[width=0.44\\textwidth]{../../meta_data/photos/chamber_1.jpeg}\\hfill",
  "\\includegraphics[width=0.44\\textwidth]{../../meta_data/photos/chamber_2.jpeg}",
  "\\caption{Front and oblique views of the acrylic chamber containing the OTICE sensor-node test configuration.}",
  "\\end{figure}",
  "\\section{Results}",
  paste0("The combined OTICE v5 dataset contains ", nrow(all_nodes), " observations from ", n_distinct(all_nodes$nodeId), " LoRa-enabled nodes, spanning ", format(min(all_nodes$bridge_time), "%Y-%m-%d %H:%M"), " to ", format(max(all_nodes$bridge_time), "%Y-%m-%d %H:%M"), ". Node ID 60 does not show a CO$_2$ reading because its CO$_2$ sensor was connected to the legacy OTICE v4 node. The additional OTICE v4 CO$_2$ graph is therefore shown separately in Figure~5 and is not included in this combined dataframe."),
  "\\begin{table}[H]",
  "\\centering",
  "\\caption{Example of the final combined dataframe, showing the earliest available record from each node ID.}",
  "\\begin{tabular}{lrrr}",
  "\\hline",
  "sensor\\_time & nodeId & nh3 & co2Raw \\\\",
  "\\hline",
  sample_tex_rows,
  "\\hline",
  "\\end{tabular}",
  "\\end{table}",
  "\\begin{figure}[H]\\centering",
  "\\includegraphics[width=\\textwidth]{../../plots/otice_v5_iot_bridge_dashboard/otice_raw_nh3_co2.pdf}",
  "\\caption{Raw NH3 and CO2 sensor values.}",
  "\\end{figure}",
  "\\begin{figure}[H]\\centering",
  "\\includegraphics[width=\\textwidth]{../../meta_data/photos/otice_v4_ee894_co2_dashboard.png}",
  "\\caption{Dashboard record from the additional EE894 CO$_2$ sensor connected to the legacy OTICE v4 node. The panels show the averaged and raw CO$_2$ channels. No additional NH$_3$ sensor was connected to OTICE v4.}",
  "\\end{figure}",
  "\\section{Inference and recommended checks}",
  "\\paragraph{CO$_2$ channel.}",
  "The chamber CO$_2$ set point was 2,000 ppm, whereas the OTICE v5 \\texttt{co2Raw} maxima were approximately 95.6 (node 50), 435.0 (node 51), 199.3 (node 52), and 337.2 (node 56); nodes 58 and 61 reported zero, and node 60 had no CO$_2$ values in the combined OTICE v5 dataset. Thus, none of the available OTICE v5 CO$_2$ channels approached the chamber set point. By contrast, the additional EE894 sensor connected to the legacy OTICE v4 node reported values on the expected 2,000 ppm scale (Figure~5; displayed range approximately 1,960--2,100 ppm and latest value 1,963 ppm). Because the same CO$_2$ sensor model was used, the comparison makes a sensor-only explanation less likely and points instead to a possible scaling, configuration, firmware, or backend-conversion issue in the OTICE v5 signal chain. This remains a working hypothesis rather than a confirmed root cause. It should be tested by comparing the direct EE894 output with the value before and after each OTICE v5 firmware/backend conversion step and by checking the configured units and correction factors.",
  "\\paragraph{NH$_3$ channel.}",
  "The exported NH$_3$ values are labelled as raw readings and span 0--1,551. The present files do not contain a documented transfer function, so the data alone cannot determine whether division by 100, division by 1,000, or another calibration equation is required to obtain ppm. \\textbf{Action for Ash:} please confirm the TB600 output definition, unit conversion, and any firmware scaling applied to the OTICE v5 \\texttt{nh3} field. Although the first available records occur at different times after the chamber-filling period, the exports contain only 100 records per node and are reportedly affected by the new 24-hour cloud-upload window. Consequently, the apparent multi-hour delay cannot yet be assigned uniquely to NH$_3$ sensor wake-up or warm-up; sensor stabilization, delayed upload, and export truncation are all plausible contributors. A follow-up test should log the local sensor output continuously from power-on while recording the corresponding cloud-arrival timestamps. The monitoring interface is available at \\url{https://iot.adaptiveagrotech.com/iot_monitoring.php}.",
  "\\end{document}"
)

writeLines(report_lines, file.path(report_dir, "otice_v5_node_comparison_report.tex"), useBytes = TRUE)
message("Completed OTICE comparison. Outputs:\n  ", clean_dir, "\n  ", plot_dir, "\n  ", report_dir)
