# NotebookLM research protocol for Manuscript 3

Notebook: `Uncertainty and Errors in Livestock Gas Emission Measurement`

Date: 2026-08-03

## Pass 1: evidence audit

NotebookLM was instructed to audit the imported Gemini draft rather than to
write manuscript prose. The requested ledger distinguished concentration,
ventilation-rate and emission errors; direct and indirect applicability;
numerical results; limitations; and confidence. It was explicitly told that
sampling was sequential, the dense-network mean was a reference estimate,
Campaign 2 was not an independent dense-network validation, no CFD was
conducted, and concentration representativeness was not emission
representativeness.

The audit covered Calvet (2013), Van Buggenhout (2009), Janke (2020), Mendes
(2015), D'Urso (2022, 2024), De Vogeleer (2017), Saha (2013), Doumbia (2021),
Rom (2010), Chen (2026), and analyser-performance sources.

## Pass 2: evidence-conservative synthesis

NotebookLM was provided with the fixed experimental design and verified R
results. It was asked to draft Introduction, Methods-improvement and Discussion
text in UK English and passive voice, using only notebook-supported citations.
The prompt prohibited claims of simultaneous mapping, absolute ground truth,
universal sensor optimality, completed CFD validation, high predictive
accuracy, and equivalence between concentration and emission sampling.

The numerical results supplied to NotebookLM included the 743 complete blocks,
configuration CCC/bias/RMSE, multiscale errors, variance fractions,
grouped-date test R-squared values, and Campaign-specific biases at location
25. This prevented the earlier imported draft's unsupported claim of model
R-squared above 0.90.

## Human/Codex curation rule

NotebookLM wording was not copied automatically. Every proposed statement was
checked against the current R outputs, the study design, publisher metadata,
and the Biosystems Engineering revision checklist. Unsupported causal wording
was recast as a plausible mechanism or removed.

