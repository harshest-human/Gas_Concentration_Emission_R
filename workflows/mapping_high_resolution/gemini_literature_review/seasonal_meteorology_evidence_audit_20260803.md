# Seasonal and meteorological framing audit

## Scope

NotebookLM notebook: *Uncertainty and Errors in Livestock Gas Emission Measurement*  
Sources selected in the notebook: 126  
Audit date: 3 August 2026

The audit was requested to support a conservative revision of Manuscript 3
after recovery of an additional dense-network period from 1--24 October 2024.
The literature was queried for evidence concerning wind, temperature,
humidity, precipitation, seasonal repetition, and spatial representativeness.

## Evidence retained for the manuscript

- Saha et al. (2013) supports the conclusion that external wind speed and
  approach direction can change sampling-point concentrations and air exchange
  in a naturally ventilated dairy building. This justifies treating wind as a
  potential modifier of spatial patterns, but not transferring numerical
  effects directly to the Gross Kreutz barn.
- Ngwabie et al. (2009) supports measurement across contrasting periods because
  gas concentrations and emissions varied temporally and seasonally in a
  naturally ventilated dairy barn.
- Ngwabie et al. (2011) supports animal activity and indoor air temperature as
  explanatory variables for methane and ammonia dynamics. Its results concern
  emission rates and cannot be assumed to be identical to spatial
  concentration effects.
- Fiedler et al. (2014) supports wind-sector-specific internal airflow and
  concentration structure. It also shows why regional station wind cannot be
  treated as barn-local airflow.
- Janke et al. (2020) supports evaluating sampling strategies under contrasting
  seasonal and approach-flow regimes, while preserving the distinction between
  concentration representativeness and emission representativeness.
- Qu et al. (2021), not 2015 as initially labelled in the NotebookLM response,
  supports temperature and relative humidity as environmental covariates in a
  meta-analysis of dairy-barn emission rates. It does not establish their
  effects on the present barn's concentration field.

## Modelling implications

- Represent wind direction circularly with sine and cosine components, or with
  pre-defined sectors tied to barn orientation.
- Permit a non-linear wind-speed response.
- Include temperature and relative humidity as physically justified
  covariates.
- Treat rain initially as an event indicator and through pre-specified lags;
  the reviewed literature did not justify rain as an independent instantaneous
  driver of indoor spatial concentration.
- Prefer the on-site 10 m mast for wind. DWD observations provide regional
  meteorological context, gap support, and sensitivity analysis.
- Keep summer (June--August) and autumn (October) as distinct periods within
  Campaign 1. Campaign 2 represents late autumn/winter but also changes the
  analyser set and sampling network.

## Claims explicitly rejected

- Seasonal causation cannot be inferred from Campaign 1 versus Campaign 2.
- DWD station measurements cannot be described as the airflow entering the
  barn.
- The dense-network mean is not an absolute or volume-weighted truth.
- A representative concentration location is not automatically representative
  of emission flux.
- Sequential multiplexed sampling is not an instantaneous spatial snapshot.
- Meteorological models must not be described as completed until their results
  have been generated and validated.

