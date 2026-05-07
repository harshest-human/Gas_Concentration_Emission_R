# Version_7 — reviewer-driven extension of v6

This folder adds the statistical analyses requested by the Biosystems
Engineering reviewers for manuscript **YBENG-D-26-00430**
("Exploring the analytical limits of the CO2-based tracer balancing
indirect method: a round-robin study"). Resubmission deadline: 20 Jun 2026.

v7 does **not** change v6. It re-uses v6 inputs (`clean_data/Version_6/`)
and v6 derived datasets (`indirect.CO2.balance()`, `reshaper()`, etc.)
and writes new outputs to:
- `result_data/plots/Version_7/`
- `result_data/tables/Version_7/`

## Directory contents

| File | Tier | Reviewer comment(s) addressed | What it does |
|------|------|------------------------------|---------------|
| `00_v7_main.R`              | -     | -                          | Top-level driver. Sources v6 functions, then runs every module below. |
| `01_setup_v7.R`             | -     | -                          | Paths, library loads, helper utilities common to all modules. |
| `02_regression_pairs.R`     | A.1   | R2-L249, R2-L353, R1-L301  | OLS regression per analyzer pair (slope, intercept, R^2, SEE, p). |
| `03_lins_ccc.R`             | A.2   | R2-L353, R2-L391           | Lin's concordance correlation coefficient with bootstrap 95 % CI. |
| `04_bland_altman_relative.R`| A.3   | R2-L353, R1-L327           | Bland-Altman plotted as RELATIVE differences (% of pair mean). |
| `05_range_dependence.R`     | A.4   | R1-L327, R2-L391           | Regression of |Diff| vs Mean to test heteroscedasticity. |
| `06_extended_pcc.R`         | A.5   | R1-L301, R2-L391           | PCC matrix per variable, per location, per time-resolution. |
| `07_ftir2_excluded.R`       | A.6   | Editor note, R2-L391       | Re-runs all comparisons with FTIR.2 NH3 column excluded. |
| `08_meteo_regression.R`     | B.1   | R1-L335 (wind influence)   | Multivariate regression: emission rates ~ wd + ws + T + RH. |
| `09_wind_conditional_ba.R`  | B.2   | R1-L335                    | Bland-Altman conditioned on wind direction sectors (N/E/S/W). |
| `10_variance_decomposition.R`| B.3  | R1-L301                    | lme4 random-effects: variance share by analyzer / day / residual. |
| `11_cross_correlation.R`    | B.4   | R2-L353                    | Cross-correlation function (CCF) of analyzer time series, lag detection. |

Tier labels follow `Analyses_brief.md`
(`D:\Doctorate\Manuscripts\Paper_2_Achievable_precision_by_Indirect_CO2_method_VERA\Analyses_brief.md`).

## How to run

```r
source("00_v7_main.R")
```

The driver script sources `01_setup_v7.R` first, then sources the v6
main analysis script up to the function definitions (using
`v6_load_functions = TRUE` flag). It then sources each numbered module
in order. Modules are independent and can also be run individually
provided `01_setup_v7.R` has been sourced.

## Reviewer-comment cross-reference

| Module(s) | Manuscript line / context |
|-----------|---------------------------|
| 02, 03, 04, 05 | R2-L249/353/391 — "Bland-Altman alone is not enough; please add regression and Lin's CCC". |
| 06 | R1-L301 — "Show the correlation more rigorously per variable, not just lumped". |
| 07 | Editor & R2 — "Be transparent about FTIR.2 NH3 outlier; show the analysis with and without". |
| 08, 09 | R1-L335 — "Discuss wind direction / wind speed effects on the background gradients". |
| 10 | R1-L301 — "How much of the spread is between analyzers vs. between days?". |
| 11 | R2-L353 — "Are the analyzers in phase, or is one running ahead/behind?". |

## Outputs

Plots are PNG (300 dpi) for figures, PDF for multi-panel composites.
Tables are tab-aware CSVs (`write_excel_csv`) so they open cleanly in
both Excel and LibreOffice.

## Provenance

Author: Harsh Sahu (ATB Potsdam).
Created: 7 May 2026, in response to first-revision review.
v6 baseline at this commit was branch `main` HEAD.
