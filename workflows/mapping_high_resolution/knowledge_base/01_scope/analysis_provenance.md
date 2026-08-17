# Manuscript 3 analysis provenance

## Pipeline

```text
raw analyser files
  -> scripts/01_prepare_all_campaign_gas_data.R
  -> clean_data/prepared_campaigns_v01/
  -> scripts/03_manuscript3_plots_tables.R
  -> clean_data/manuscript3_analysis_v01/tables/
  -> plots/manuscript3_analysis_v01/
  -> manuscript figures and reported values
```

Regional meteorological context is prepared separately:

```text
DWD source data
  -> scripts/02_manuscript3_dwd_weather.R
  -> plots/manuscript3_weather/
  -> manuscript wind-rose figure
```

## Principal analysis outputs

| Output | Purpose |
|---|---|
| `Table_00_sampling_point_analysis_eligibility.csv` | SP eligibility and coverage audit |
| `Table_02_campaign_descriptives.csv` | campaign-level descriptive statistics |
| `Table_07_H1_height_descriptives_2h_blocks.csv` | vertical-group summaries from two-hour blocks |
| `Table_08_H1_mixed_model_height_estimates.csv` | model-derived height estimates |
| `Table_09_H1_mixed_model_height_contrasts.csv` | height contrasts and multiplicity adjustment |
| `Table_10_H1_precision_by_temporal_resolution.csv` | precision across temporal aggregation scales |
| `Table_17_SP_representativeness_CV_robustCV.csv` | individual-SP temporal variability |
| `Table_28_standard_CO2_balance_hourly_annualised_summary.csv` | standard-input ventilation and emission summaries |
| `Table_29_full_network_campaign_descriptives.csv` | full-network campaign statistics |
| `Table_30_vertical_group_effects_and_accuracy.csv` | vertical effect and network-accuracy metrics |
| `Table_31_spatial_temporal_variance_partition.csv` | temporal and persistent spatial variance components |
| `Table_32_diurnal_each_SP_descriptives.csv` | hourly summaries for every SP |
| `Table_33_diurnal_functional_zone_summary.csv` | hourly functional-zone summaries |

## Reproducibility rule

Numerical manuscript claims must be traceable to an R-generated table. Console output and manually calculated values are not accepted as sources. If analysis changes, rerun the active script and update this map before editing manuscript numbers.
