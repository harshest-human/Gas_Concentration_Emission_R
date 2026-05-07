######################################################################
# 05_pipeline.R   (v8 driver)
# ---------------------------------------------------------------------
# Top-level driver for v8. Sources 01..04, runs the negative-value
# investigation, writes README_v8.md.
######################################################################

# Always source from this directory so relative paths in 02..04 work.
v8_script_dir <- "D:/Data_Analysis_R/Gas_Concentration_Emission_R/workflows/ringversuche_analysis/scripts/v8"
old_wd <- getwd(); setwd(v8_script_dir); on.exit(setwd(old_wd), add = TRUE)

source("01_relabel_NE_SW.R", local = FALSE)
message("================================================================")
message("v8 pipeline: phase 1 - relabel V6+V7 tables -> Version_8/tables")
message("================================================================")
source("04_relabel_v6_v7_tables.R", local = FALSE)
v8_relabel_all_tables()

message("================================================================")
message("v8 pipeline: phase 2 - re-render v6 plots -> Version_8/plots")
message("================================================================")
source("02_replot_v6_concentrations.R", local = FALSE)
v8_run_v6_plots()

message("================================================================")
message("v8 pipeline: phase 3 - re-render v7 plots")
message("================================================================")
source("03_replot_v7_outputs.R", local = FALSE)
v8_run_v7_plots()

message("================================================================")
message("v8 pipeline: phase 4 - negative-value inventory")
message("================================================================")
v8_neg_value_inventory <- function() {
  suppressPackageStartupMessages({ library(readr); library(dplyr) })
  out_path <- file.path(v8_base, "result_data/tables/Version_8/negative_value_inventory.csv")
  cr <- read_csv(file.path(v8_v6_tables, "20250408-15_ringversuche_concentration_reshaped.csv"),
                 show_col_types = FALSE, guess_max = 50000)
  cr <- relabel_ne_sw_strings(cr)
  neg <- cr %>%
    filter(is.finite(value), value < 0) %>%
    mutate(gas = case_when(
      grepl("CO2",  var) ~ "CO2",
      grepl("CH4",  var) ~ "CH4",
      grepl("NH3",  var) ~ "NH3",
      TRUE               ~ var
    )) %>%
    group_by(var, gas, analyzer, location) %>%
    summarise(n_negative = n(),
              min_value  = min(value, na.rm = TRUE),
              mean_value = round(mean(value, na.rm = TRUE), 4),
              .groups = "drop")
  total <- cr %>% group_by(var, analyzer, location) %>%
    summarise(n_total = n(), .groups = "drop")
  out <- neg %>%
    left_join(total, by = c("var", "analyzer", "location")) %>%
    mutate(pct_negative = round(100 * n_negative / n_total, 2)) %>%
    select(var, gas, analyzer, location,
           n_negative, n_total, pct_negative,
           min_value, mean_value) %>%
    arrange(desc(n_negative))
  readr::write_excel_csv(out, out_path)
  message("v8: wrote negative_value_inventory.csv  (", nrow(out), " rows)")
  invisible(out)
}
v8_neg_value_inventory()

message("================================================================")
message("v8 pipeline: phase 5 - README_v8.md")
message("================================================================")
v8_write_readme <- function() {
  rd_path <- file.path(v8_script_dir, "README_v8.md")

  v6_plots <- list.files(file.path(v8_base, "result_data/plots/Version_6"),
                         pattern = "\\.png$") |> sort()
  v7_plots <- list.files(file.path(v8_base, "result_data/plots/Version_7"),
                         pattern = "\\.png$") |> sort()
  v8_plots <- list.files(file.path(v8_base, "result_data/plots/Version_8"),
                         pattern = "\\.png$") |> sort()
  v6_tabs  <- list.files(file.path(v8_base, "result_data/tables/Version_6"),
                         pattern = "\\.csv$") |> sort()
  v7_tabs  <- list.files(file.path(v8_base, "result_data/tables/Version_7"),
                         pattern = "\\.csv$") |> sort()
  v8_tabs  <- list.files(file.path(v8_base, "result_data/tables/Version_8"),
                         pattern = "\\.csv$") |> sort()

  fmt_li <- function(x) if (length(x) == 0) "(none)\n" else paste0("- `", x, "`\n", collapse = "")

  txt <- paste0(
"# Version_8 - NE/SW relabel + concentration y-clip + facet spacing\n",
"\n",
"Manuscript: **YBENG-D-26-00430** (Biosystems Engineering, first revision).\n",
"Reviewer comments addressed: **R2.18** (NE/SW relabel) and **R1.12**\n",
"(physically valid concentration y-axis).\n",
"\n",
"## Transformations applied\n",
"\n",
"1. `_N -> _NE` and `_S -> _SW` in column-name suffixes; `north background\n",
"   -> NE background` and `south background -> SW background` in string\n",
"   values, axis labels, legends, captions, and CSV row/column headers.\n",
"   Source data in `clean_data/Version_6/` is **untouched**.\n",
"2. `coord_cartesian(ylim = c(0, NA))` applied to every plot showing a\n",
"   gas concentration (absolute or delta). Sub-zero readings (analyser\n",
"   noise / quantification artefacts) are suppressed from the figure but\n",
"   retained in the underlying CSV. See `negative_value_inventory.csv`\n",
"   in `result_data/tables/Version_8/` for the count by gas x analyser x\n",
"   location.\n",
"3. `theme(panel.spacing.y = unit(20, 'pt'))` applied to every multi-row\n",
"   free-y faceted plot so different gas units don't visually merge.\n",
"\n",
"## Suggested manuscript caption note (concentration figures)\n",
"\n",
"> Note: y-axis clipped at 0; sub-zero readings (analyser noise near\n",
"> the limit of detection or numerical artefacts at near-zero deltas)\n",
"> are suppressed from the plot; values are retained in the underlying\n",
"> data and tabulated in Table S-NV (`negative_value_inventory.csv`).\n",
"\n",
"## Inventory\n",
"\n",
"### Version_6 plots (input, ", length(v6_plots), " files)\n", fmt_li(v6_plots),
"\n### Version_7 plots (input, ", length(v7_plots), " files)\n", fmt_li(v7_plots),
"\n### Version_8 plots (output, ", length(v8_plots), " files)\n", fmt_li(v8_plots),
"\n### Version_6 tables (input, ", length(v6_tabs), " files)\n", fmt_li(v6_tabs),
"\n### Version_7 tables (input, ", length(v7_tabs), " files)\n", fmt_li(v7_tabs),
"\n### Version_8 tables (output, ", length(v8_tabs), " files)\n", fmt_li(v8_tabs),
"\n## Pipeline files\n",
"\n",
"- `01_relabel_NE_SW.R` - relabel functions, theme helper, write_csv shim\n",
"- `02_replot_v6_concentrations.R` - re-renders v6 plots into Version_8\n",
"- `03_replot_v7_outputs.R` - re-renders v7 plots into Version_8 by\n",
"   sourcing v7 modules with `v7_to_wide` defaults and `v7_save_plot`\n",
"   destination wrapped in-memory.\n",
"- `04_relabel_v6_v7_tables.R` - mirrors V6 + V7 CSVs into Version_8 with\n",
"   relabel_ne_sw applied to columns and string values.\n",
"- `05_pipeline.R` - top-level driver. Source this file to rebuild v8.\n",
"\n",
"## Stopping-rule items / manual handling needed\n",
"\n",
"- `weather_trendplot.png`: contains no NE/SW labels (temperature, wind\n",
"   axes); copied through from V6 unchanged.\n",
"- `06_extended_pcc`, `07_ftir2_excluded`, `08_meteo_regression`,\n",
"   `11_cross_correlation`: tables only, relabelled by 04. No plots in v7.\n",
"- The wind-direction sector axis in `09_wind_conditional_BA_summary.png`\n",
"   uses N/E/S/W as compass labels (per the brief, those stay - they are\n",
"   not the background-line N/S labels).\n",
"\n",
"## Provenance\n",
"\n",
"Author: Harsh Sahu (ATB Potsdam).\n",
"Pipeline runs over `clean_data/Version_6/`, source PNGs and CSVs in\n",
"V6/V7 unchanged; outputs in `result_data/{plots,tables}/Version_8/`.\n"
  )
  writeLines(txt, rd_path)
  message("v8: wrote README_v8.md")
}
v8_write_readme()

message("================================================================")
message("v8 pipeline: DONE")
message("================================================================")
