######################################################################
# 03_replot_v7_outputs.R   (v8)
# ---------------------------------------------------------------------
# Re-renders the v7 figures with NE/SW labels and y >= 0 clipping
# where appropriate. Strategy:
#   1) Source v7 setup (sets v7_to_wide, helpers, paths).
#   2) Override v7_to_wide / v7_save_plot to use NE/SW location names
#      and Version_8 plot dir respectively.
#   3) Load relabeled concentration_reshaped & emission_reshaped.
#   4) Source v7 modules that PRODUCE PLOTS (02, 03, 04, 09, 10) one
#      by one. Sourcing rebuilds the plots from the relabeled data.
#   5) Move/relabel the resulting tables (handled in 04_relabel...).
#
# Modules without plots (05, 06, 07, 08, 11) are skipped here; their
# tables are picked up by 04_relabel_v6_v7_tables.R.
######################################################################

if (!exists("relabel_ne_sw")) source("01_relabel_NE_SW.R")

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
  library(patchwork); library(stringr)
})

v8_base       <- "D:/Data_Analysis_R/Gas_Concentration_Emission_R/workflows/ringversuche_analysis"
v8_v7_dir     <- file.path(v8_base, "scripts", "v7")
v8_v6_tables  <- file.path(v8_base, "result_data/tables/Version_6")
v8_plots_out  <- file.path(v8_base, "result_data/plots/Version_8")
dir.create(v8_plots_out, showWarnings = FALSE, recursive = TRUE)

v8_run_v7_plots <- function() {
  # ---- 1) source v7 setup ----
  old_wd <- getwd(); on.exit(setwd(old_wd), add = TRUE)
  setwd(v8_v7_dir)
  source("01_setup_v7.R", local = FALSE)

  # ---- 2) override v7 helpers in-memory (do NOT touch source files) ----
  ne_sw_locs <- c("NE background", "SW background")

  v7_to_wide_orig <- v7_to_wide
  assign("v7_to_wide", function(df_long, var_name, locs = ne_sw_locs) {
    v7_to_wide_orig(df_long, var_name, locs)
  }, envir = globalenv())

  v7_save_plot_orig <- v7_save_plot
  assign("v7_save_plot", function(plot, filename, width = 8, height = 6, dpi = 300) {
    p2 <- tryCatch(apply_v8_theme(plot,
                                  clip_y0 = grepl("regression_(delta|CO2|CH4|NH3)|concentration",
                                                  filename, ignore.case = TRUE),
                                  panel_spacing_pt = 20),
                   error = function(e) plot)
    out <- file.path(v8_plots_out, filename)
    suppressMessages(suppressWarnings(
      ggplot2::ggsave(out, p2, width = width, height = height, dpi = dpi)
    ))
    message("v8: wrote v7-style plot ", filename)
    invisible(out)
  }, envir = globalenv())

  v7_write_csv_orig <- v7_write_csv
  # redirect tables to a v8 temp dir to avoid clobbering v7 tables
  v8_v7_tmp_tables <- file.path(v8_base, "result_data/tables/Version_8/_tmp_v7_rerun")
  dir.create(v8_v7_tmp_tables, showWarnings = FALSE, recursive = TRUE)
  assign("v7_write_csv", function(df, filename) {
    path <- file.path(v8_v7_tmp_tables, filename)
    readr::write_excel_csv(df, path)
    message("v8: (v7 rerun temp) wrote ", filename)
    invisible(path)
  }, envir = globalenv())

  # ---- 3) load v6 reshaped data + relabel locations ----
  cr <- read_csv(file.path(v8_v6_tables, "20250408-15_ringversuche_concentration_reshaped.csv"),
                 show_col_types = FALSE, guess_max = 50000)
  er <- read_csv(file.path(v8_v6_tables, "20250408-15_ringversuche_emission_reshaped.csv"),
                 show_col_types = FALSE, guess_max = 50000)
  ew <- read_csv(file.path(v8_v6_tables, "20250408-15_ringversuche_emission_result.csv"),
                 show_col_types = FALSE, guess_max = 50000)
  if (!inherits(cr$DATE.TIME, "POSIXct")) cr <- cr %>% mutate(DATE.TIME = as.POSIXct(DATE.TIME))
  if (!inherits(er$DATE.TIME, "POSIXct")) er <- er %>% mutate(DATE.TIME = as.POSIXct(DATE.TIME))
  if (!inherits(ew$DATE.TIME, "POSIXct")) ew <- ew %>% mutate(DATE.TIME = as.POSIXct(DATE.TIME))
  cr <- relabel_ne_sw_strings(cr)
  er <- relabel_ne_sw_strings(er)
  ew <- relabel_ne_sw_cols(ew)   # wide format: rename _N->_NE, _S->_SW columns
  assign("concentration_reshaped", cr, envir = globalenv())
  assign("emission_reshaped",      er, envir = globalenv())
  assign("emission_result",        ew, envir = globalenv())

  # ---- 4) source v7 plot-producing modules ----
  modules <- c(
    "02_regression_pairs.R",
    "03_lins_ccc.R",
    "04_bland_altman_relative.R",
    # 09 and 10 also produce plots; 09 needs input_combined which has wd column
    # so we attempt them but trap errors to keep the pipeline running.
    "09_wind_conditional_ba.R",
    "10_variance_decomposition.R"
  )
  for (m in modules) {
    res <- tryCatch({
      message("v8: re-running ", m, " ...")
      source(m, local = FALSE)
      "ok"
    }, error = function(e) {
      message("v8: WARN ", m, " failed: ", e$message)
      conditionMessage(e)
    })
  }

  # ---- 5) restore originals (defensive) ----
  assign("v7_to_wide",   v7_to_wide_orig,   envir = globalenv())
  assign("v7_save_plot", v7_save_plot_orig, envir = globalenv())
  assign("v7_write_csv", v7_write_csv_orig, envir = globalenv())

  invisible(NULL)
}

if (sys.nframe() == 0L) v8_run_v7_plots()
