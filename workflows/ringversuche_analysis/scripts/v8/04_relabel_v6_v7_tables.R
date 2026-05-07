######################################################################
# 04_relabel_v6_v7_tables.R   (v8)
# ---------------------------------------------------------------------
# Reads every CSV under result_data/tables/Version_6 and Version_7,
# applies relabel_ne_sw (column-suffix + string-value swap),
# writes mirrored filenames into result_data/tables/Version_8/.
#
# Source CSVs are NEVER overwritten.
######################################################################

if (!exists("relabel_ne_sw")) source("01_relabel_NE_SW.R")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

v8_base    <- "D:/Data_Analysis_R/Gas_Concentration_Emission_R/workflows/ringversuche_analysis"
v8_v6_in   <- file.path(v8_base, "result_data/tables/Version_6")
v8_v7_in   <- file.path(v8_base, "result_data/tables/Version_7")
v8_tab_out <- file.path(v8_base, "result_data/tables/Version_8")
dir.create(v8_tab_out, showWarnings = FALSE, recursive = TRUE)

v8_relabel_one_csv <- function(path_in, path_out) {
  if (!file.exists(path_in)) return(invisible(FALSE))
  df <- suppressMessages(read_csv(path_in, show_col_types = FALSE,
                                  guess_max = 50000))
  df2 <- tryCatch(relabel_ne_sw(df), error = function(e) {
    message("v8: relabel failed on ", basename(path_in), ": ", e$message,
            " - copying as-is.")
    df
  })
  write_excel_csv(df2, path_out)
  message("v8: wrote ", basename(path_out),
          "  (", nrow(df2), " rows x ", ncol(df2), " cols)")
  invisible(TRUE)
}

v8_relabel_all_tables <- function() {
  csvs <- c(
    list.files(v8_v6_in, pattern = "\\.csv$", full.names = TRUE),
    list.files(v8_v7_in, pattern = "\\.csv$", full.names = TRUE)
  )
  for (f in csvs) {
    v8_relabel_one_csv(f, file.path(v8_tab_out, basename(f)))
  }
  invisible(NULL)
}

# Run when sourced standalone
if (sys.nframe() == 0L) v8_relabel_all_tables()
