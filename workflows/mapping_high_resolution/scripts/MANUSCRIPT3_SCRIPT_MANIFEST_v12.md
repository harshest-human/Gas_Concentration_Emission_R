# Manuscript 3 analysis-script manifest (v12)

## Canonical scripts

- `01_prepare_all_campaign_gas_data.R`: the unified upstream gas-preparation
  script. It reads FTIR raw exports and CRDS step-average inputs, provides a
  raw CRDS reader and 60-s-flush averaging function, applies the campaign
  metadata mappings, writes analyser-wise traceable files, writes one combined
  CSV per campaign, and row-binds all four campaigns. Its validation outputs
  are isolated in `clean_data/prepared_campaigns_v01`.
- `03_manuscript3_plots_tables.R`: the current Manuscript 3 plotting and table
  pipeline for Campaigns 1 and 2. Median-based representativeness is primary;
  mean-based performance is retained for comparison. It produces shared-scale
  gas and ratio boxplots, mean +/- SD and median +/- SD figures, median-based
  sampling-point rankings, and Campaign 2 gas-and-ratio Shannon entropy. It
  deliberately excludes CV and relative error and does not modify LaTeX.
- `dwd_weather_figures_manuscript3.R`: retained as the second canonical script because it creates the DWD weather input and the two meteorological figures used by the manuscript.

New Manuscript 3 concentration analyses, tables, and plots must be added as
functions or clearly labelled sections in `03_manuscript3_plots_tables.R`, not
as additional analysis scripts. CV and relative error remain deferred until
they are explicitly reinstated in the study plan.

`manuscript3_pipeline_v12.R` is retained temporarily only to reproduce the
currently compiled pre-revision manuscript. It is superseded for new analysis
by `03_manuscript3_plots_tables.R` and can be recycled once the new figures and
tables are approved for manuscript integration.

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
- The canonical R scripts above are edited in place and tracked with Git.
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

A second cleanup on 2026-08-10 moved two broken legacy R scripts, the copied
TeX-version directory, superseded `manuscript1` outputs, LaTeX build
intermediates, and unused manuscript-figure copies to the Windows Recycle Bin.
The canonical scripts, accepted campaign inputs, current v12 outputs,
`main.tex`, `main.pdf`, and all figures referenced by `main.tex` were verified
afterwards.
