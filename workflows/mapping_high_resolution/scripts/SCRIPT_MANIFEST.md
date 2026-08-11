# High-resolution mapping script manifest

## Active gas workflow

- `01_prepare_all_campaign_gas_data.R`: reads FTIR and CRDS inputs, applies
  flushing and instrument correction, maps analyser selectors to
  `sampling.point`, saves analyser-wise and campaign-wise CSVs, and row-binds
  Campaigns 1--4.
- `03_manuscript3_plots_tables.R`: current Campaign 1--2 Manuscript 3 analysis.
  It produces median-based representativeness tables, mean and median with SD,
  shared-scale concentration and ratio plots, and Shannon entropy. It does not
  modify LaTeX and currently excludes CV and relative error.
- `dwd_weather_figures_manuscript3.R`: prepares DWD regional meteorology and
  produces the Manuscript 3 weather figures.

## Active Campaign 3--4 supporting workflow

- `clean_animal_cooling_data.R`: cleans animal counts and cooling states.
- `merge_animal_cooling_climate_campaign3_4.R`: merges animal/cooling records
  with five-minute barn climate and mast weather.
- `campaign_cv_floorplan_and_3d.R`: retained for future Campaign 3--4 spatial
  visualisation.
- `campaign_re_floorplan_and_3d.R`: retained for future Campaign 3--4 spatial
  visualisation.
- `campaign_cv_re_floorplan_3d_node.js`: retained as the existing HTML/3D
  renderer until the Manuscript 4 workflow is consolidated.
- `build_manuscript4_research_synopsis_html.js`: retained for the existing
  Manuscript 4 synopsis HTML.

## Archived workflow

Superseded R scripts and old plot outputs were moved on 2026-08-11 to:

```text
archive/pre_unified_2026-08-11/scripts
archive/pre_unified_2026-08-11/plots
```

The archived scripts remain version-controlled and recoverable. Archived plots
remain local and are not pushed to GitHub. Current outputs remain in
`plots/manuscript3_analysis_v01` and `plots/manuscript3_weather`.

## Version-control policy

- Edit the active scripts in place; use Git commits and tags for versions.
- Do not create copied R or TeX files for routine revisions.
- Generated data, tables and plots are reproducible outputs and are not edited
  manually.
