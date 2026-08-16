# Paper 3: Mapping gas concentrations

This folder contains the LaTeX draft for Manuscript 3 on Campaigns 1 and 2.

## Build

Compile `main.tex` with a LaTeX installation that provides `natbib`,
`siunitx`, `subcaption`, `lineno`, and `booktabs`:

```text
pdflatex main
bibtex main
pdflatex main
pdflatex main
```

## Reproducible analysis

The unified concentration-analysis, statistics, table, and plotting script is:

```text
D:/Data_Analysis_R/Gas_Concentration_Emission_R/workflows/
mapping_high_resolution/scripts/manuscript3_pipeline_v12.R
```

The separate DWD preparation and meteorological-figure script is:

```text
D:/Data_Analysis_R/Gas_Concentration_Emission_R/workflows/
mapping_high_resolution/scripts/02_manuscript3_dwd_weather.R
```

Clean manuscript data and tables are written to:

```text
D:/Data_Analysis_R/Gas_Concentration_Emission_R/workflows/
mapping_high_resolution/clean_data/manuscript3_v12
```

Plots are written to:

```text
D:/Data_Analysis_R/Gas_Concentration_Emission_R/workflows/
mapping_high_resolution/plots/manuscript3_v12
```

`main.tex` is the single canonical manuscript source. Git commits and tags,
rather than copied TeX files, provide version history.

## Items requiring author confirmation

- Final author list and contribution roles.
- Funding project names and grant numbers.
- Public data-repository location and DOI.
- Exact proprietary analyser model, manufacturer, city, and country details.
- Whether the corresponding author should remain David Janke.
- Whether Campaign 2 locations 19 and 40 should be described as already
  excluded because of fan proximity and sampling-line leakage.

The journal currently requests Microsoft Word source files by default.
The LaTeX draft can be converted after the scientific content is stable.

## Citation, reference, and prose policy

- Cite every scholarly source with a `natbib` command such as `\citet{}` or
  `\citep{}`. Do not type author--year citations manually.
- Load `natbib` before `hyperref`. Every in-text citation must link to its
  corresponding bibliography entry in the compiled PDF.
- Every citation key must resolve in `references.bib`. Each reference should
  include a DOI when one exists, or a stable source URL otherwise.
- Compile with BibTeX and confirm that no undefined citation or reference
  warnings remain before committing manuscript changes.
- Use full stops rather than semicolons in manuscript prose.
