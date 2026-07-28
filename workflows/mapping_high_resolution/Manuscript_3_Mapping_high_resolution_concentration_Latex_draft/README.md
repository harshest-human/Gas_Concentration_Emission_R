# Paper 3: Mapping gas concentrations

This folder contains the LaTeX draft for Manuscript 1 on Campaigns 1 and 2.

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

The data-processing and plotting script is:

```text
D:/Data_Analysis_R/Gas_Concentration_Emission_R/workflows/
mapping_high_resolution/scripts/manuscript1_campaign1_2_analysis.R
```

Clean manuscript data and tables are written to:

```text
D:/Data_Analysis_R/Gas_Concentration_Emission_R/workflows/
mapping_high_resolution/clean_data/manuscript1
```

Plots are written to:

```text
D:/Data_Analysis_R/Gas_Concentration_Emission_R/workflows/
mapping_high_resolution/plots/manuscript1
```

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
