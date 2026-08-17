# Manuscript 3 project state

Last updated: 2026-08-17

## Active manuscript

`Manuscript_3_Mapping_high_resolution_concentration_Latex_draft/LaTex_drafts/20260813_Manuscript3_draft.tex`

The bibliography is `references.bib` in the same directory. The draft compiles successfully with manual `pdflatex`, `bibtex`, `pdflatex`, `pdflatex` passes. `latexmk` is unavailable locally because MiKTeX has no Perl engine.

## Scientific framing

The research question is how horizontal position, vertical level and temporal aggregation affect measured CO2, CH4 and NH3 concentrations and their pairwise ratios in an operational naturally ventilated dairy barn.

The current null hypothesis is that the long-term distributions of the measured gases do not differ among horizontal SPs or vertical levels. Effects on CO2-balance ventilation and emission estimates are engineering consequences rather than a second hypothesis.

Campaign 1 is the primary dense field experiment. Campaign 2 is a constrained seasonal extension.

## Campaign definitions

- Campaign 1: 1 June to 31 August 2024 and 1 to 24 October 2024.
- Campaign 1 network: 51 internal SPs across 17 horizontal positions and top, middle and bottom levels.
- Campaign 1 analysers: CRDS8, CRDS9, FTIR1 and FTIR2.
- Campaign 2: 16 November to 31 December 2024.
- Campaign 2 network: 32 top-and-bottom SPs measured by CRDS8 and CRDS9 after FTIR unavailability.
- SP19 and SP40 remain included.

## Processing decisions

- Campaign 1 FTIR harmonisation: CO2 and CH4 divided by 1.06. NH3 divided by 1.09.
- CRDS valve dwell: 240 s. The first 60 s are excluded and the remaining 180 s are averaged.
- Three-minute sampling intervals are excluded from the current manuscript analysis.
- No 300 ppm lower threshold is applied to CO2.
- Missing, zero and negative concentrations are excluded response-wise from concentration statistics.
- No IQR or other statistical outlier filter is applied.
- Pairwise gas ratios are dimensionless.
- Clock-aligned two-hour blocks are used for contemporaneous network comparisons.
- The outside line is used where required for Campaign 1 balance calculations. It is not treated as an internal SP.
- Campaign 2 lacks sufficient contemporaneous outside data for comparative emission inference.

## Active reproducible scripts

1. `scripts/01_prepare_all_campaign_gas_data.R`
2. `scripts/02_manuscript3_dwd_weather.R`
3. `scripts/03_manuscript3_plots_tables.R`
4. `scripts/literature_ingest.ps1`

The first three scripts remain the only active Manuscript 3 data and statistical workflow. Older scripts remain in the archive.

## Literature status

- The Introduction and Materials and methods now cite Apostolico et al. (2026), Declerck et al. (2025) and Janke et al. (2022).
- Declerck (2025) and Janke et al. (2022) have primary-source-verified evidence notes.
- Other central papers are registered but still require individual evidence-note verification.
- The Declerck report must not be described as proving a sampling-height bias. It compared complete measurement methods.
- The direct method in Janke et al. (2022) was treated as a conditional reference under selected stable cross-flow conditions, not as a universal reference.

## Current editing state

The Introduction and relevant Materials and methods passages were revised on 2026-08-17. Results and Discussion were not revised during that operation. Existing user changes in the manuscript and bibliography remain uncommitted.

The literature workflow was reorganised on 2026-08-17. The source registry contains the central Manuscript 3 literature. Local searchable text has been generated for Declerck2025 and Janke2022. The extracted text and its machine-specific manifest are ignored by Git.

Last repository commit observed when this file was created:

```text
a9b08be309400516a6762a9d65a7ec5d2397cd4c
```

## Next task

Continue supervisor-led review of the revised Introduction and Materials and methods. Do not rewrite Results or Discussion unless explicitly requested. When another paper is introduced, register it, extract it once, create a verified evidence note, update the appropriate synthesis map and only then amend the manuscript.
