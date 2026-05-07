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

  # ---- 4a) v8 regression plot builder (post-source rerender) ----
  # We source v7 module 02 below to get the regression TABLE, but its
  # for-loop overwrites our function override; so we re-render the 6
  # regression plots AFTER the module finishes, with math titles.
  v8_render_regression <- function() {
    var_set <- c("delta_CO2","delta_CH4","delta_NH3",
                 "Q_vent","e_CH4_ghLU","e_NH3_ghLU")
    for (vn in var_set) {
      df_use <- if (vn %in% c("Q_vent","e_CH4_ghLU","e_NH3_ghLU")) emission_reshaped else concentration_reshaped
      var_expr <- v8_var_expr_short[[vn]]
      plots <- lapply(v7_pairs_cross, function(pr) {
        wide <- v7_to_wide(df_use, vn)
        v <- v7_pair_vectors(wide, pr[1], pr[2])
        if (is.null(v)) return(ggplot2::ggplot() + ggplot2::theme_void())
        d <- tibble::tibble(x = v$x, y = v$y)
        ttl <- parse(text = sprintf("%s ~ ': ' ~ '%s vs %s'",
                                    var_expr, pr[1], pr[2]))[[1]]
        ggplot2::ggplot(d, ggplot2::aes(x = x, y = y)) +
          ggplot2::geom_point(alpha = 0.4, size = 0.7) +
          ggplot2::geom_abline(slope = 1, intercept = 0,
                               linetype = "dashed", colour = "grey40") +
          ggplot2::geom_smooth(method = "lm", formula = y ~ x,
                               se = TRUE, colour = "red", linewidth = 0.6) +
          ggplot2::labs(x = pr[1], y = pr[2], title = ttl) +
          ggplot2::theme_classic(base_size = 9) +
          ggplot2::theme(plot.title = ggplot2::element_text(size = 9))
      })
      v7_save_plot(patchwork::wrap_plots(plots, ncol = 3),
                   paste0("02_regression_", vn, ".png"),
                   width = 11, height = 11)
    }
  }

  # ---- 4b) source v7 plot-producing modules ----
  # Lin's CCC needs a y-axis re-order and parsed math labels; instead
  # of touching the v7 source we override the ggplot it builds AFTER
  # module 03 finishes (the code uses v7_save_plot which we already
  # wrapped, so we intercept the plot at v7_save_plot time when the
  # filename is the CCC lattice).
  v7_save_plot_outer <- v7_save_plot
  ccc_var_order <- c("Q_vent","e_CH4_ghLU","e_NH3_ghLU",
                     "delta_CO2","delta_CH4","delta_NH3",
                     "CO2_mgm3","CH4_mgm3","NH3_mgm3")
  ccc_var_labs  <- v8_var_expr_short[ccc_var_order]

  assign("v7_save_plot", function(plot, filename, width = 8, height = 6, dpi = 300) {
    if (filename == "03_lins_ccc_lattice.png" && inherits(plot, "ggplot")) {
      plot <- plot +
        ggplot2::scale_y_discrete(
          limits = ccc_var_order,
          labels = function(x) parse(text = unname(ccc_var_labs[x]))
        ) +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 11))
    }
    v7_save_plot_outer(plot, filename, width = width, height = height, dpi = dpi)
  }, envir = globalenv())

  modules <- c(
    "02_regression_pairs.R",
    "03_lins_ccc.R",
    "04_bland_altman_relative.R",
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

  # ---- 4c) re-render regression scatters with v8 math-expression titles
  message("v8: re-rendering regression scatter with math titles...")
  tryCatch(v8_render_regression(),
           error = function(e) message("v8: regression rerender WARN: ", e$message))

  # ---- 5) restore originals (defensive) ----
  assign("v7_to_wide",   v7_to_wide_orig,   envir = globalenv())
  assign("v7_save_plot", v7_save_plot_orig, envir = globalenv())
  assign("v7_write_csv", v7_write_csv_orig, envir = globalenv())

  invisible(NULL)
}

if (sys.nframe() == 0L) v8_run_v7_plots()
