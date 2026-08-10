# Manuscript 3 analysis-script manifest (v12)

## Canonical scripts

- `manuscript3_pipeline_v12.R`: the single executable pipeline for Campaigns 1 and 2, ratios, descriptive statistics, SP25 comparisons, multiscale precision, mixed-effects models, Shannon entropy, meteorological associations, tables, and figures.
- `dwd_weather_figures_manuscript3.R`: retained as the second canonical script because it creates the DWD weather input and the two meteorological figures used by the manuscript.

All manuscript concentration outputs---including time series, relative error
(RE), coefficient of variation (CV), confidence intervals, median absolute
deviation, entropy, statistical models, summary tables, and manuscript
figures---must be produced through `manuscript3_pipeline_v12.R`. New analyses
should be added as functions or clearly labelled sections in this script, not
as additional manuscript-analysis scripts.

## Supporting upstream scripts retained

These create accepted clean inputs and remain outside the manuscript-analysis consolidation:

- `clean_FTIR_raw.R`
- `clean_CRDS_step_average.R`
- `process_campaign1_instrument_correction.R`
- `recover_campaign1_crds.R`
- `recover_campaign1_crds_june_august.R`
- `combine_campaign2_CRDS.R`

## Version-control policy

- `main.tex` is the single canonical manuscript source. Scientific and
  editorial changes are tracked with Git commits instead of copied TeX files.
- The two canonical R scripts above are edited in place and tracked with Git.
- Generated tables, clean derivatives, and figures are reproducible outputs;
  they must not be edited manually.
- A Git tag may mark a submitted or otherwise frozen scientific version. A
  copied TeX or R version is created only when an external submission system
  explicitly requires a frozen standalone package.

## Cleanup completed

On 2026-08-10, six redundant Manuscript 3 analysis scripts, four superseded
clean-output directories, four superseded plot directories, and 17 obsolete
manuscript figure copies were moved to the Windows Recycle Bin. They remain
recoverable until the Recycle Bin is emptied. The two canonical scripts and
all upstream cleaning dependencies were retained.
