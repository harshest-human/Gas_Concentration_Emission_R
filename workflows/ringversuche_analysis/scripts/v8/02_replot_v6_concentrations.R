######################################################################
# 02_replot_v6_concentrations.R   (v8)
# ---------------------------------------------------------------------
# Re-renders the v6 concentration / delta / emission plots with NE/SW
# labels, optional y >= 0 clip on absolute & delta concentration, and
# facet row spacing for multi-row free-y composites. Mirrors the v6
# filenames so the manuscript build can swap V6 for V8 paths.
#
# Strategy:
#   1) Source v6 function definitions only (lines 1..1122 of
#      Ringversuche_analysis_script.R - everything before the analysis
#      section). This avoids triggering a fresh v6 run.
#   2) Read the v6-already-reshaped long tables from
#      result_data/tables/Version_6/ as inputs (concentration_reshaped,
#      emission_reshaped). These contain the canonical schema:
#      var, value, location, analyzer, DATE.TIME, day.
#   3) Relabel `location` strings to NE/SW.
#   4) Call each v6 plot function with `locations = c("NE background",
#      "SW background")` (or with "Barn inside" prepended).
#   5) Post-mod each ggplot with apply_v8_theme(clip_y0=..., spacing).
#   6) Save with the original filename to result_data/plots/Version_8/.
######################################################################

if (!exists("relabel_ne_sw")) source("01_relabel_NE_SW.R")

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(lubridate)
  library(tidyr); library(ggplot2); library(patchwork)
  library(scales); library(viridis); library(ggcorrplot); library(ggridges)
  library(rstatix); library(multcompView)
  library(kableExtra); library(knitr); library(reshape2)
})

v8_base       <- "D:/Data_Analysis_R/Gas_Concentration_Emission_R/workflows/ringversuche_analysis"
v8_v6_script  <- file.path(v8_base, "scripts", "Ringversuche_analysis_script.R")
v8_v6_tables  <- file.path(v8_base, "result_data/tables/Version_6")
v8_plots_out  <- file.path(v8_base, "result_data/plots/Version_8")
dir.create(v8_plots_out, showWarnings = FALSE, recursive = TRUE)

# ---- 1) Source v6 function definitions only ------------------------
v8_source_v6_funcs <- function() {
  lines <- readLines(v8_v6_script, warn = FALSE)
  # Split point: first occurrence of "######## Import Data" line
  import_idx <- grep("^######## Import Data", lines)[1]
  if (is.na(import_idx)) stop("could not find ## Import Data marker in v6 script")
  src <- paste(lines[seq_len(import_idx - 1L)], collapse = "\n")
  eval(parse(text = src), envir = globalenv())
  message("v8: sourced v6 function defs (1..", import_idx - 1L, ")")
}

# ---- 2) Load v6 reshaped data and relabel --------------------------
v8_load_inputs <- function() {
  cr <- read_csv(file.path(v8_v6_tables, "20250408-15_ringversuche_concentration_reshaped.csv"),
                 show_col_types = FALSE, guess_max = 50000)
  er <- read_csv(file.path(v8_v6_tables, "20250408-15_ringversuche_emission_reshaped.csv"),
                 show_col_types = FALSE, guess_max = 50000)
  # ensure DATE.TIME parses
  if (!inherits(cr$DATE.TIME, "POSIXct")) cr <- cr %>% mutate(DATE.TIME = as.POSIXct(DATE.TIME))
  if (!inherits(er$DATE.TIME, "POSIXct")) er <- er %>% mutate(DATE.TIME = as.POSIXct(DATE.TIME))
  # add a day column if not present (some v6 plots expect it)
  if (!"day" %in% names(cr)) cr$day <- factor(as.Date(cr$DATE.TIME))
  if (!"day" %in% names(er)) er$day <- factor(as.Date(er$DATE.TIME))
  # relabel location column to NE/SW
  cr <- relabel_ne_sw_strings(cr)
  er <- relabel_ne_sw_strings(er)
  list(concentration_reshaped = cr, emission_reshaped = er)
}

# ---- 3) Save helper that applies v8 cosmetics ---------------------
v8_ggsave <- function(plot, filename, width, height, clip_y0 = FALSE,
                      panel_spacing_pt = 20) {
  p2 <- tryCatch(apply_v8_theme(plot, clip_y0 = clip_y0,
                                panel_spacing_pt = panel_spacing_pt),
                 error = function(e) plot)
  out <- file.path(v8_plots_out, filename)
  suppressMessages(suppressWarnings(
    ggsave(out, p2, width = width, height = height, dpi = 300)
  ))
  message("v8: wrote plot ", filename)
}

# ---- 4) Run the v6 plot pipeline against relabeled inputs ----------
v8_run_v6_plots <- function() {
  v8_source_v6_funcs()
  inp <- v8_load_inputs()
  concentration_reshaped <- inp$concentration_reshaped
  emission_reshaped      <- inp$emission_reshaped

  ne_sw <- c("NE background", "SW background")
  in_ne_sw <- c("Barn inside", "NE background", "SW background")

  # ---------- HSD / RE matrix plots (multi-gas pairwise_diff) -------
  abs_cd  <- pairwise_diff(concentration_reshaped,
                           vars = c("CO2_mgm3","CH4_mgm3","NH3_mgm3"),
                           locations = in_ne_sw)
  delta_cd <- pairwise_diff(concentration_reshaped,
                            vars = c("delta_CO2","delta_CH4","delta_NH3"),
                            locations = ne_sw)
  vent_ed  <- pairwise_diff(emission_reshaped,
                            vars = c("Q_vent","e_CH4_ghLU","e_NH3_ghLU"),
                            locations = ne_sw)
  v8_ggsave(abs_cd,   "abs_concentration_diff.png",   12, 6, clip_y0 = TRUE)
  v8_ggsave(delta_cd, "delta_concentration_diff.png", 9,  6, clip_y0 = TRUE)
  v8_ggsave(vent_ed,  "vent_emission_diff.png",       9,  6, clip_y0 = FALSE)

  # ---------- absolute concentration trend / errorbar --------------
  c_trend  <- emitrendplot(concentration_reshaped,
                           y = c("CO2_mgm3","CH4_mgm3","NH3_mgm3"),
                           plot_err = FALSE)
  c_err    <- emierrorbarplot(concentration_reshaped,
                              y = c("CO2_mgm3","CH4_mgm3","NH3_mgm3"),
                              plot_err = FALSE)
  v8_ggsave(c_trend, "c_trend_plot.png",   12, 8, clip_y0 = TRUE)
  v8_ggsave(c_err,   "c_errorbarplot.png", 12, 8, clip_y0 = TRUE)

  # ---------- delta concentration trend / errorbar -----------------
  d_trend <- emitrendplot(concentration_reshaped,
                          y = c("delta_CO2","delta_CH4","delta_NH3"),
                          plot_err = FALSE)
  d_err   <- emierrorbarplot(concentration_reshaped,
                             y = c("delta_CO2","delta_CH4","delta_NH3"),
                             plot_err = FALSE)
  v8_ggsave(d_trend, "d_trend_plot.png",   12, 8, clip_y0 = TRUE)
  v8_ggsave(d_err,   "d_errorbarplot.png", 12, 8, clip_y0 = TRUE)

  # ---------- vent + emission trend / errorbar (NOT clipped) -------
  q_trend <- emitrendplot(emission_reshaped,
                          y = c("Q_vent","e_CH4_ghLU","e_NH3_ghLU"),
                          plot_err = FALSE)
  q_err   <- emierrorbarplot(emission_reshaped,
                             y = c("Q_vent","e_CH4_ghLU","e_NH3_ghLU"),
                             plot_err = FALSE)
  v8_ggsave(q_trend, "q_e_trend_plot.png",   12, 8, clip_y0 = FALSE)
  v8_ggsave(q_err,   "q_e_errorbarplot.png", 12, 8, clip_y0 = FALSE)

  # ---------- Bland-Altman composites (ATB pair as Analyzer A) -----
  v8_ba_save <- function(filename, var_pair, ana_pair, locs, w, h) {
    plots <- lapply(seq_along(var_pair), function(i) {
      v <- var_pair[[i]]; loc <- locs[[i]]
      bland_altman_plot(emission_reshaped, var_filter = v,
                        analyzer_pair = ana_pair, location_filter = loc)
    })
    plots <- lapply(plots, function(p) p + theme(plot.margin = margin(10,10,10,10)))
    p <- patchwork::wrap_plots(plots, ncol = length(plots), nrow = 1)
    v8_ggsave(p, filename, w, h, clip_y0 = FALSE)
  }
  # e_BlandAltman_AnalyzerA: 4 panels  CH4_NE, NH3_NE, CH4_SW, NH3_SW
  v8_ba_save("e_BlandAltman_AnalyzerA.png",
             list("e_CH4_ghLU","e_NH3_ghLU","e_CH4_ghLU","e_NH3_ghLU"),
             c("FTIR.1","CRDS.1"),
             list("NE background","NE background","SW background","SW background"),
             16, 4)
  # q_BlandAltman_AnalyzerA: 2 panels  Q_NE, Q_SW
  v8_ba_save("q_BlandAltman_AnalyzerA.png",
             list("Q_vent","Q_vent"),
             c("FTIR.1","CRDS.1"),
             list("NE background","SW background"),
             10, 4)
  # Analyzer B (ATB FTIR vs LUFA CRDS pair, mimicking v6 pattern)
  v8_ba_save("e_BlandAltman_AnalyzerB.png",
             list("e_CH4_ghLU","e_NH3_ghLU","e_CH4_ghLU","e_NH3_ghLU"),
             c("FTIR.2","CRDS.2"),
             list("NE background","NE background","SW background","SW background"),
             16, 4)
  v8_ba_save("q_BlandAltman_AnalyzerB.png",
             list("Q_vent","Q_vent"),
             c("FTIR.2","CRDS.2"),
             list("NE background","SW background"),
             10, 4)

  # ---------- corrgrams ---------------------------------------------
  d_cg  <- tryCatch(emicorrgram(concentration_reshaped,
                                target_variables = c("delta_CO2","delta_CH4","delta_NH3"),
                                locations = ne_sw),
                    error = function(e) { message("emicorrgram delta: ", e$message); NULL })
  if (!is.null(d_cg))  v8_ggsave(d_cg,  "d_corrgram.png",   10, 6, clip_y0 = FALSE)
  qe_cg <- tryCatch(emicorrgram(emission_reshaped,
                                target_variables = c("e_CH4_ghLU","e_NH3_ghLU","Q_vent"),
                                locations = ne_sw),
                    error = function(e) { message("emicorrgram emission: ", e$message); NULL })
  if (!is.null(qe_cg)) v8_ggsave(qe_cg, "q_e_corrgram.png", 10, 6, clip_y0 = FALSE)

  # ---------- box+jitter per analyzer per location (NEW in v8) -----
  # Replaces the v6 errorbar mean+SE plots; full boxplots show the
  # per-cycle distribution rather than just the summary point.
  cb <- v8_box_plot(concentration_reshaped,
                    vars = c("CO2_mgm3","CH4_mgm3","NH3_mgm3"),
                    locations = in_ne_sw)
  db <- v8_box_plot(concentration_reshaped,
                    vars = c("delta_CO2","delta_CH4","delta_NH3"),
                    locations = ne_sw)
  qb <- v8_box_plot(emission_reshaped,
                    vars = c("Q_vent","e_CH4_ghLU","e_NH3_ghLU"),
                    locations = ne_sw)
  v8_ggsave(cb, "c_boxplot.png",   12, 9, clip_y0 = TRUE)
  v8_ggsave(db, "d_boxplot.png",   10, 9, clip_y0 = TRUE)
  v8_ggsave(qb, "q_e_boxplot.png", 10, 9, clip_y0 = FALSE)

  # CV heatmaps removed in v8 (per reviewer feedback - not needed).

  # ---------- CV scatter (emipointplot) ----------------------------
  v8_pp_save <- function(filename, df, vars, w = 10, h = 7) {
    p <- tryCatch(emipointplot(df, vars = vars, locations = ne_sw),
                  error = function(e) { message("emipointplot: ", e$message); NULL })
    if (!is.null(p)) v8_ggsave(p, filename, w, h, clip_y0 = FALSE)
  }
  v8_pp_save("d_cvplot.png",    concentration_reshaped,
             c("delta_CO2","delta_CH4","delta_NH3"))
  v8_pp_save("q_e_cvplot.png",  emission_reshaped,
             c("Q_vent","e_CH4_ghLU","e_NH3_ghLU"))

  # ---------- weather trend (no N/S labels, just copy through) -----
  src <- file.path(v8_base, "result_data/plots/Version_6/weather_trendplot.png")
  if (file.exists(src)) {
    file.copy(src, file.path(v8_plots_out, "weather_trendplot.png"), overwrite = TRUE)
    message("v8: copied weather_trendplot.png unchanged (no N/S labels).")
  }

  invisible(NULL)
}

if (sys.nframe() == 0L) v8_run_v6_plots()
