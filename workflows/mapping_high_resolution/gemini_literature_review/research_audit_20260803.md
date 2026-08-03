# Scientific audit of Gemini/NotebookLM material

## Corrections applied

- Sequential valve switching was retained; the dataset was not described as a
  simultaneous three-dimensional snapshot.
- The 51-location mean was called a dense-network reference estimate rather
  than ground truth.
- The Campaign 1 FTIR factors were described as empirical cross-platform
  alignment to the CRDS operational scale, not absolute calibration.
- Campaign 2 was described as an operational application of the reduced design,
  not an independent 51-versus-32 validation.
- Location 25 was restricted to a site-specific concentration approximation;
  no claim of universal or emission representativeness was retained.
- Chen (2026) was used as conceptual support for information-based sampling.
  The present entropy calculation was identified as an empirical adaptation,
  not a replication of the CFD-informed optimisation algorithm.
- Grouped-date predictive performance was reported as moderate. The imported
  Gemini claim of R-squared greater than 0.90 was rejected.
- Mechanisms involving wind, fans, openings, manure surfaces and animal sources
  were separated from effects directly demonstrated by the present data.

## Important cross-check against a primary source

NotebookLM's first audit stated that Van Buggenhout et al. (2009) supported a
maximum error of 33% and that the earlier 86% value should be removed. The
publisher abstract instead states that ventilation-rate errors for internal
sampling positions could rise to 86% of the actual ventilation rate (DOI:
10.1016/j.biosystemseng.2009.04.018). The manuscript therefore retains 86%,
but explicitly scopes it to a tracer-gas ventilation-rate estimate in a
mechanically ventilated test installation. This discrepancy demonstrates why
NotebookLM citations were treated as leads rather than final verification.

## Interpretive boundaries retained

- Gas ratios are proportionality diagnostics, not source apportionment or
  emission factors.
- Weekly aggregation reduces transient location error but cannot demonstrate
  spatial homogeneity.
- Whole-campaign CV includes the complete temporal range and is not a fourth
  stage in a monotonic averaging sequence.
- Agreement in the barn mean can occur through compensating vertical errors;
  preservation of mean and preservation of spatial spread are separate tests.
- Any emission application additionally requires representative airflow at a
  defined control surface.

