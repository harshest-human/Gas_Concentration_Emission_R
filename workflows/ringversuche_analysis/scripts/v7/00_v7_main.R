######################################################################
# 00_v7_main.R
# ---------------------------------------------------------------------
# Top-level driver for Version_7 (reviewer-revision extension to v6).
#
# How it works:
#   1. Source 01_setup_v7.R  -> defines paths, libs, helpers.
#   2. Source the v6 main analysis script up through the data-import
#      and emission-computation block. We do this by sourcing the
#      whole v6 file but stopping just before the heavy plotting
#      section so we don't re-render every v6 figure. We control
#      this with the env var `RUN_V6_PLOTS` (set FALSE here).
#   3. Source each v7 module in numerical order.
#
# The v6 driver writes its outputs to result_data/plots/Version_6/
# and result_data/tables/Version_6/. v7 only ever writes to
# Version_7 sub-directories. v6 outputs are never touched.
######################################################################

# ---- Phase 1: setup --------------------------------------------------
# Discover this script's directory robustly across source(), Rscript,
# and RStudio's source-on-save.
v7_dir <- tryCatch({
        f <- sys.frame(1)$ofile
        if (is.null(f) || !nzchar(f)) stop("no ofile")
        dirname(normalizePath(f))
}, error = function(e) {
        # Fallback to absolute path. Adjust if you move the project.
        file.path("D:/Data_Analysis_R/Gas_Concentration_Emission_R",
                  "workflows/ringversuche_analysis/scripts/v7")
})
source(file.path(v7_dir, "01_setup_v7.R"))

# ---- Phase 2: load v6 datasets --------------------------------------
# We source the v6 main script. It will run end-to-end. To avoid
# wasting time and disk space on regenerating Version_6 plots in this
# overnight run, set the option below to TRUE only when you want a
# full v6 + v7 re-run.
#
# This driver is a SCRIPT not a function, so we rely on the v6 main
# script to leave `concentration_reshaped`, `emission_reshaped`, and
# `emission_result` in the global environment. Those three frames
# are what every v7 module consumes.

v6_main <- file.path(v7_dir, "..", "Ringversuche_analysis_script.R")
v6_main <- normalizePath(v6_main, mustWork = TRUE)
message("v7 driver: sourcing v6 main script at ", v6_main)
source(v6_main, echo = FALSE)

stopifnot(exists("concentration_reshaped"),
          exists("emission_reshaped"),
          exists("emission_result"))

# ---- Phase 3: run v7 modules in order -------------------------------
v7_modules <- c(
        "02_regression_pairs.R",
        "03_lins_ccc.R",
        "04_bland_altman_relative.R",
        "05_range_dependence.R",
        "06_extended_pcc.R",
        "07_ftir2_excluded.R",
        "08_meteo_regression.R",
        "09_wind_conditional_ba.R",
        "10_variance_decomposition.R",
        "11_cross_correlation.R"
)

for (mod in v7_modules) {
        message("\n===== v7 module: ", mod, " =====")
        source(file.path(v7_dir, mod), echo = FALSE)
}

message("\nAll v7 modules complete.")
message("Tables in: ", v7_tables_dir)
message("Plots  in: ", v7_plots_dir)
