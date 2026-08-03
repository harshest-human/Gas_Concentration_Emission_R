# External wind, turbulence and background concentration evidence audit

Date: 3 August 2026

## Scope

This audit was produced from a focused query to the user's NotebookLM notebook,
*Uncertainty and Errors in Livestock Gas Emission Measurement* (126 uploaded
sources), followed by an empirical analysis of Campaign 1 gas data and hourly
DWD Potsdam observations. NotebookLM was instructed to use only uploaded
sources and to distinguish concentration, concentration enhancement,
ventilation rate and emission rate.

## Core source evidence retained

- Saha et al. (2013) measured external wind, temperature, humidity, and gas
  concentrations in a naturally ventilated dairy building. Wind speed and
  direction affected air change, and wind conditions were associated with
  spatial concentrations. The study did not supply a simple correlation
  coefficient that can be transferred to the present data.
- König et al. (2018) showed that the indoor/outdoor sampling-point choice and
  the resulting air-exchange estimate depended on approach direction. This
  supports wind-sector-aware representativeness, not a universal reference
  point.
- Janke et al. (2020, Biosystems Engineering) compared concentration sampling
  strategies across wind directions and showed that lateral and cross-wise
  flow produced different spatial concentration behaviour and emission
  estimates.
- Janke et al. (2020, Sensors) measured velocity and tracer concentration in a
  wind-tunnel model. It supports vertically resolved sampling, but its stable
  controlled inflow is not equivalent to the occupied barn's temporal weather.
- Schmithausen et al. (2018) used wind velocity and cross-flow criteria when
  retaining periods for emission calculation. The criterion concerns the
  validity of an emission method, not a direct wind-concentration correlation.
- Fiedler et al. (2014) reported a negative CO2-wind relationship, but a source
  phrase rendered as `p = -0.7` is statistically ambiguous and was not repeated
  numerically in the manuscript.

## Empirical result in the present study

Campaign 1 summer and October data were hourly aligned with DWD Potsdam wind.
The south-background point `s`, internal spatial mean, and inside-minus-`s`
enhancement were analysed separately. Zero-lag Spearman correlations between
regional wind speed and internal means were weak in summer (CO2 -0.14, CH4
-0.17, NH3 -0.20) and stronger in October (-0.43, -0.41 and -0.53). October
enhancement correlations were -0.40, -0.39 and -0.56; background correlations
were only -0.13, -0.27 and -0.10. Date-block bootstrap intervals excluded zero
for all October internal and enhancement outcomes.

These are exploratory associations. They do not prove a turbulent dilution
mechanism. The DWD station is approximately 18.5 km from the barn; hourly DWD
wind is neither local opening velocity nor turbulence intensity. No 2024
high-frequency on-site mast file was available in the workflow. Sequential gas
sampling also prevents an instantaneous spatial map.

## Manuscript decisions

- Report the measured negative regional-wind associations, including the
  period dependence and weaker association at `s`.
- Describe dilution/air exchange only as a compatible mechanism.
- Do not call DWD hourly wind "turbulence" or infer turbulent flux.
- Do not transfer published emission or ventilation effects directly to gas
  concentration.
- Preserve the separation between `s`, internal mean and enhancement.
- Treat selected lag maxima as descriptive because adjacent lags are dependent
  and multiple lags were scanned.

## Reproducibility

The analysis is implemented in
`scripts/manuscript3_external_wind_background_analysis.R`. Generated aligned
data and result tables are written to `clean_data/manuscript3_weather`; figures
are written to `plots/manuscript3_weather` and the manuscript figure directory.
Generated data and figures are intentionally excluded from version control.
