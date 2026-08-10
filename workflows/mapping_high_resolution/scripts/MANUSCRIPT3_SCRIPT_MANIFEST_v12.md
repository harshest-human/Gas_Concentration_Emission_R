# Manuscript 3 analysis-script manifest (v12)

## Canonical script

- `manuscript3_pipeline_v12.R`: the single executable pipeline for Campaigns 1 and 2, ratios, descriptive statistics, SP25 comparisons, multiscale precision, mixed-effects models, Shannon entropy, meteorological associations, tables, and figures.
- `dwd_weather_figures_manuscript3.R`: retained as the second canonical script because it creates the DWD weather input and the two meteorological figures used by the manuscript.

## Supporting upstream scripts retained

These create accepted clean inputs and remain outside the manuscript-analysis consolidation:

- `clean_FTIR_raw.R`
- `clean_CRDS_step_average.R`
- `process_campaign1_instrument_correction.R`
- `recover_campaign1_crds.R`
- `recover_campaign1_crds_june_august.R`
- `combine_campaign2_CRDS.R`

## Redundant manuscript-analysis candidates for archival

The following overlap with the canonical v12 pipeline. They must not be deleted until the v12 numerical audit is approved:

- `manuscript1_campaign1_2_analysis.R`
- `manuscript3_campaign1_june_august_analysis_v03.R`
- `manuscript3_campaign1_recovered_analysis.R`
- `manuscript3_external_wind_background_analysis.R`
- `manuscript3_multiscale_reference_analysis_v04.R`
- `manuscript3_spatiotemporal_statistics.R`

## Archival rule

After approval, redundant scripts should be moved together to `scripts/archive/manuscript3_pre_v12/`. Git history, rather than additional filename suffixes, will track subsequent corrections. A new versioned output directory should be created only when a scientific analysis is frozen.

## Cleanup completed

On 2026-08-10, the six redundant analysis scripts listed above, their four superseded clean-output directories, four superseded plot directories, and 17 obsolete manuscript figure copies were moved to the Windows Recycle Bin. They remain recoverable until the Recycle Bin is emptied. Current v12 files and DWD dependencies were retained.
