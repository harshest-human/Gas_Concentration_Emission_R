######################################################################
# 01_setup_v7.R
# ---------------------------------------------------------------------
# Common paths, library loads and helpers used by all v7 modules.
# Sourcing this file is idempotent: it only initialises objects that
# are missing in the current R session.
######################################################################

# ---- Paths -----------------------------------------------------------
v7_base_dir   <- "D:/Data_Analysis_R/Gas_Concentration_Emission_R/workflows/ringversuche_analysis"
v7_data_dir   <- file.path(v7_base_dir, "clean_data",   "Version_6")  # use v6 cleaned data as input
v7_meta_dir   <- file.path(v7_base_dir, "meta_data")
v7_plots_dir  <- file.path(v7_base_dir, "result_data", "plots",  "Version_7")
v7_tables_dir <- file.path(v7_base_dir, "result_data", "tables", "Version_7")

dir.create(v7_plots_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(v7_tables_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Period (matches v6 driver) -------------------------------------
v7_start_time <- as.POSIXct("2025-04-08 12:00:00", tz = "UTC")
v7_end_time   <- as.POSIXct("2025-04-14 12:00:00", tz = "UTC")

# ---- Library loads (silent, but report any missing) -----------------
v7_libs <- c(
        "tidyverse",  "lubridate",  "dplyr",      "tidyr",       "readr",
        "ggplot2",    "patchwork",  "scales",     "stringr",
        "DescTools",  "rstatix",    "broom",      "broom.mixed",
        "lme4",       "lmerTest",   "performance",
        "Hmisc",      "psych",      "epiR",       "boot"
)
.v7_missing <- v7_libs[!vapply(v7_libs, requireNamespace, logical(1L), quietly = TRUE)]
if (length(.v7_missing) > 0L) {
        message("v7 setup: the following packages are missing and should be installed:\n  ",
                paste(.v7_missing, collapse = ", "))
}
suppressPackageStartupMessages({
        invisible(lapply(setdiff(v7_libs, .v7_missing), library, character.only = TRUE))
})

# ---- Analyzer aesthetics (mirror v6) --------------------------------
v7_analyzer_levels <- c("FTIR.1","FTIR.2","FTIR.3","FTIR.4",
                        "CRDS.1","CRDS.2","CRDS.3","baseline")
v7_analyzer_colors <- c(
        "FTIR.1"   = "#1b9e77", "FTIR.2"   = "#d95f02", "FTIR.3"   = "#7570b3",
        "FTIR.4"   = "#e7298a", "CRDS.1"   = "#66a61e", "CRDS.2"   = "#e6ab02",
        "CRDS.3"   = "#a6761d", "baseline" = "black"
)
v7_analyzer_shapes <- c(
        "FTIR.1" = 0,  "FTIR.2" = 1,  "FTIR.3" = 2,  "FTIR.4" = 5,
        "CRDS.1" = 15, "CRDS.2" = 19, "CRDS.3" = 17, "baseline" = 4
)

# ---- Reference analyzer pairs (lab-internal pairs) ------------------
# A = ATB, B = LUFA, C = UB. ANECO has FTIR only (FTIR.4); MBBM has FTIR only (FTIR.3).
v7_pairs_lab <- list(
        "ATB"  = c("FTIR.1", "CRDS.1"),
        "LUFA" = c("FTIR.2", "CRDS.2"),
        "UB"   = c("FTIR.3", "CRDS.3")  # NB: UB pair uses MBBM's FTIR.3 as proxy for cross-lab CRDS.3
)
# Cross-platform pairs - FTIR vs CRDS, all combinations
v7_pairs_cross <- list(
        c("FTIR.1","CRDS.1"), c("FTIR.1","CRDS.2"), c("FTIR.1","CRDS.3"),
        c("FTIR.2","CRDS.1"), c("FTIR.2","CRDS.2"), c("FTIR.2","CRDS.3"),
        c("FTIR.3","CRDS.1"), c("FTIR.3","CRDS.2"), c("FTIR.3","CRDS.3"),
        c("FTIR.4","CRDS.1"), c("FTIR.4","CRDS.2"), c("FTIR.4","CRDS.3")
)
# FTIR-FTIR
v7_pairs_ftir <- list(
        c("FTIR.1","FTIR.2"), c("FTIR.1","FTIR.3"), c("FTIR.1","FTIR.4"),
        c("FTIR.2","FTIR.3"), c("FTIR.2","FTIR.4"), c("FTIR.3","FTIR.4")
)
# CRDS-CRDS
v7_pairs_crds <- list(
        c("CRDS.1","CRDS.2"), c("CRDS.1","CRDS.3"), c("CRDS.2","CRDS.3")
)
v7_pairs_all <- c(v7_pairs_cross, v7_pairs_ftir, v7_pairs_crds)

# ---- Helper: pivot a v6 emission_reshaped df to wide on analyzer ----
# Returns a wide tibble keyed by DATE.TIME + location, with one numeric
# column per analyzer for the requested variable. Used by every module
# that compares analyzer pairs.
v7_to_wide <- function(df_long, var_name, locs = c("North background","South background")) {
        df_long %>%
                dplyr::filter(.data$var == var_name,
                              .data$location %in% locs,
                              .data$analyzer != "baseline") %>%
                dplyr::select(DATE.TIME, location, analyzer, value) %>%
                tidyr::pivot_wider(names_from = analyzer, values_from = value,
                                   values_fn = ~ mean(.x, na.rm = TRUE))
}

# ---- Helper: for an analyzer pair, return aligned vectors -----------
v7_pair_vectors <- function(df_wide, a1, a2) {
        if (!all(c(a1, a2) %in% names(df_wide))) return(NULL)
        out <- df_wide %>%
                dplyr::select(DATE.TIME, location, dplyr::all_of(c(a1, a2))) %>%
                dplyr::filter(complete.cases(.))
        if (nrow(out) < 5L) return(NULL)
        list(x = out[[a1]], y = out[[a2]],
             location = out$location, time = out$DATE.TIME, n = nrow(out))
}

# ---- Helper: nice rounding helper for tables ------------------------
v7_round <- function(x, digits = 3L) {
        if (is.numeric(x)) signif(x, digits) else x
}

# ---- Helper: write a CSV and a console message ----------------------
v7_write_csv <- function(df, filename) {
        path <- file.path(v7_tables_dir, filename)
        readr::write_excel_csv(df, path)
        message("v7: wrote table -> ", filename, "  (", nrow(df), " rows)")
        invisible(path)
}

# ---- Helper: ggsave wrapper -----------------------------------------
v7_save_plot <- function(plot, filename, width = 8, height = 6, dpi = 300) {
        path <- file.path(v7_plots_dir, filename)
        suppressMessages(ggplot2::ggsave(path, plot, width = width, height = height, dpi = dpi))
        message("v7: wrote plot  -> ", filename)
        invisible(path)
}

# ---- Wind-direction sector helper -----------------------------------
# 0/360 = N, 90 = E, 180 = S, 270 = W. Returns 4-letter sector.
v7_wd_sector <- function(wd_deg) {
        if (any(is.na(wd_deg))) {
                ifelse(is.na(wd_deg), NA_character_,
                       v7_wd_sector_inner(wd_deg))
        } else {
                v7_wd_sector_inner(wd_deg)
        }
}
v7_wd_sector_inner <- function(wd) {
        wd <- (wd %% 360)
        cuts <- c(45, 135, 225, 315)
        labels <- c("N", "E", "S", "W")
        out <- rep("N", length(wd))
        out[wd >= 45  & wd < 135] <- "E"
        out[wd >= 135 & wd < 225] <- "S"
        out[wd >= 225 & wd < 315] <- "W"
        out
}

message("v7 setup complete: outputs to ", v7_plots_dir, " and ", v7_tables_dir)
