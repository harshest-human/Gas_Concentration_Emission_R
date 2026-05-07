# Version_8 - NE/SW relabel + concentration y-clip + facet spacing

Manuscript: **YBENG-D-26-00430** (Biosystems Engineering, first revision).
Reviewer comments addressed: **R2.18** (NE/SW relabel) and **R1.12**
(physically valid concentration y-axis).

## Transformations applied

1. `_N -> _NE` and `_S -> _SW` in column-name suffixes; `north background
   -> NE background` and `south background -> SW background` in string
   values, axis labels, legends, captions, and CSV row/column headers.
   Source data in `clean_data/Version_6/` is **untouched**.
2. `coord_cartesian(ylim = c(0, NA))` applied to every plot showing a
   gas concentration (absolute or delta). Sub-zero readings (analyser
   noise / quantification artefacts) are suppressed from the figure but
   retained in the underlying CSV. See `negative_value_inventory.csv`
   in `result_data/tables/Version_8/` for the count by gas x analyser x
   location.
3. `theme(panel.spacing.y = unit(20, 'pt'))` applied to every multi-row
   free-y faceted plot so different gas units don't visually merge.

## Suggested manuscript caption note (concentration figures)

> Note: y-axis clipped at 0; sub-zero readings (analyser noise near
> the limit of detection or numerical artefacts at near-zero deltas)
> are suppressed from the plot; values are retained in the underlying
> data and tabulated in Table S-NV (`negative_value_inventory.csv`).

## Inventory

### Version_6 plots (input, 24 files)
- `abs_concentration_diff.png`
- `c_errorbarplot.png`
- `c_trend_plot.png`
- `d_CH4_heatmap.png`
- `d_CO2_heatmap.png`
- `d_corrgram.png`
- `d_cvplot.png`
- `d_errorbarplot.png`
- `d_NH3_heatmap.png`
- `d_trend_plot.png`
- `delta_concentration_diff.png`
- `e_BlandAltman_AnalyzerA.png`
- `e_BlandAltman_AnalyzerB.png`
- `e_CH4_heatmap.png`
- `e_NH3_heatmap.png`
- `q_BlandAltman_AnalyzerA.png`
- `q_BlandAltman_AnalyzerB.png`
- `q_e_corrgram.png`
- `q_e_cvplot.png`
- `q_e_errorbarplot.png`
- `q_e_trend_plot.png`
- `q_heatmap.png`
- `vent_emission_diff.png`
- `weather_trendplot.png`

### Version_7 plots (input, 12 files)
- `02_regression_delta_CH4.png`
- `02_regression_delta_CO2.png`
- `02_regression_delta_NH3.png`
- `02_regression_e_CH4_ghLU.png`
- `02_regression_e_NH3_ghLU.png`
- `02_regression_Q_vent.png`
- `03_lins_ccc_lattice.png`
- `04_BA_relative_lab_ATB.png`
- `04_BA_relative_lab_LUFA.png`
- `04_BA_relative_lab_UB.png`
- `09_wind_conditional_BA_summary.png`
- `10_variance_decomposition.png`

### Version_8 plots (output, 33 files)
- `02_regression_delta_CH4.png`
- `02_regression_delta_CO2.png`
- `02_regression_delta_NH3.png`
- `02_regression_e_CH4_ghLU.png`
- `02_regression_e_NH3_ghLU.png`
- `02_regression_Q_vent.png`
- `03_lins_ccc_lattice.png`
- `04_BA_relative_lab_ATB.png`
- `04_BA_relative_lab_LUFA.png`
- `04_BA_relative_lab_UB.png`
- `09_wind_conditional_BA_summary.png`
- `10_variance_decomposition.png`
- `abs_concentration_diff.png`
- `c_boxplot.png`
- `c_errorbarplot.png`
- `c_trend_plot.png`
- `d_boxplot.png`
- `d_corrgram.png`
- `d_cvplot.png`
- `d_errorbarplot.png`
- `d_trend_plot.png`
- `delta_concentration_diff.png`
- `e_BlandAltman_AnalyzerA.png`
- `e_BlandAltman_AnalyzerB.png`
- `q_BlandAltman_AnalyzerA.png`
- `q_BlandAltman_AnalyzerB.png`
- `q_e_boxplot.png`
- `q_e_corrgram.png`
- `q_e_cvplot.png`
- `q_e_errorbarplot.png`
- `q_e_trend_plot.png`
- `vent_emission_diff.png`
- `weather_trendplot.png`

### Version_6 tables (input, 6 files)
- `20250408-15_ringversuche_concentration_reshaped.csv`
- `20250408-15_ringversuche_emission_reshaped.csv`
- `20250408-15_ringversuche_emission_result.csv`
- `20250408-15_ringversuche_input_combined_data.csv`
- `concentration_stat.csv`
- `emission_stat.csv`

### Version_7 tables (input, 12 files)
- `02_regression_pairs.csv`
- `03_lins_ccc.csv`
- `04_bland_altman_relative.csv`
- `05_range_dependence.csv`
- `06_extended_pcc.csv`
- `07_ftir2_sensitivity_ccc.csv`
- `07_ftir2_sensitivity_regression.csv`
- `07_ftir2_sensitivity_summary.csv`
- `08_meteo_regression.csv`
- `09_wind_conditional_BA.csv`
- `10_variance_decomposition.csv`
- `11_cross_correlation.csv`

### Version_8 tables (output, 19 files)
- `02_regression_pairs.csv`
- `03_lins_ccc.csv`
- `04_bland_altman_relative.csv`
- `05_range_dependence.csv`
- `06_extended_pcc.csv`
- `07_ftir2_sensitivity_ccc.csv`
- `07_ftir2_sensitivity_regression.csv`
- `07_ftir2_sensitivity_summary.csv`
- `08_meteo_regression.csv`
- `09_wind_conditional_BA.csv`
- `10_variance_decomposition.csv`
- `11_cross_correlation.csv`
- `20250408-15_ringversuche_concentration_reshaped.csv`
- `20250408-15_ringversuche_emission_reshaped.csv`
- `20250408-15_ringversuche_emission_result.csv`
- `20250408-15_ringversuche_input_combined_data.csv`
- `concentration_stat.csv`
- `emission_stat.csv`
- `negative_value_inventory.csv`

## Pipeline files

- `01_relabel_NE_SW.R` - relabel functions, theme helper, write_csv shim
- `02_replot_v6_concentrations.R` - re-renders v6 plots into Version_8
- `03_replot_v7_outputs.R` - re-renders v7 plots into Version_8 by
   sourcing v7 modules with `v7_to_wide` defaults and `v7_save_plot`
   destination wrapped in-memory.
- `04_relabel_v6_v7_tables.R` - mirrors V6 + V7 CSVs into Version_8 with
   relabel_ne_sw applied to columns and string values.
- `05_pipeline.R` - top-level driver. Source this file to rebuild v8.

## Stopping-rule items / manual handling needed

- `weather_trendplot.png`: contains no NE/SW labels (temperature, wind
   axes); copied through from V6 unchanged.
- `06_extended_pcc`, `07_ftir2_excluded`, `08_meteo_regression`,
   `11_cross_correlation`: tables only, relabelled by 04. No plots in v7.
- The wind-direction sector axis in `09_wind_conditional_BA_summary.png`
   uses N/E/S/W as compass labels (per the brief, those stay - they are
   not the background-line N/S labels).

## Provenance

Author: Harsh Sahu (ATB Potsdam).
Pipeline runs over `clean_data/Version_6/`, source PNGs and CSVs in
V6/V7 unchanged; outputs in `result_data/{plots,tables}/Version_8/`.

