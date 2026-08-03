# Campaigns 1–2 multiscale and reference-location analysis protocol (v04)

## Primary question

Do sampling locations retain systematic concentration signatures through time,
or do their differences converge towards the contemporaneous barn mean as the
averaging duration increases?

## Secondary reference-location objective

Identify a physically defensible single reference location that best represents
the spatial network across gases, time blocks and both campaigns. Physical
eligibility is applied before statistical ranking.

## Pre-specified physical screen

- Location 40: excluded because of confirmed sampling-line leakage.
- Location 19: excluded because of confirmed local fan influence.
- Conservative fan-corridor sensitivity screen: horizontal triplets 19–21 and
  43–45, identified from the experimental floor plan.
- Only locations measured in both campaigns can be cross-campaign candidates.
- Final engineering acceptance requires author confirmation against the exact
  fan axes and intended emission-measurement plane.

## Analysis units and leakage prevention

- Concentrations are summarised by location within contemporaneous clock-time
  blocks. The barn reference is recalculated inside every block.
- Cross-validation is grouped by complete date; rows from one date cannot occur
  in both training and testing sets.
- No randomly split row-level validation is reported.
- Absolute minima and maxima are descriptive only; 5th and 95th percentiles are
  used as robust extremes.

## Statistical analyses

1. Hourly, daily, weekly and whole-campaign mean, median, minimum, maximum, SD,
   MAD, CV, and 5th and 95th percentiles. Whole-campaign statistics are
   calculated directly from all valid observations, not averaged from weeks.
   Scale-dependent SD, MAD and CV are tested using repeated-location Friedman
   tests.
2. Mixed-effects variance decomposition of log concentrations and ratios, with
   random intercepts for location and date.
3. Temporal convergence of absolute relative error at 1, 2, 4, 8, 12, 24, 72
   and 168 h.
4. Daily rank stability using Spearman correlations and Kendall's W.
5. Cyclic-hour GAMs with location and date effects.
6. PCA and correlation-distance hierarchical clustering of location profiles.
7. Four-bin Shannon entropy and discretised mutual information with the
   contemporaneous barn mean.
8. Change-point screening of daily barn mean and spatial CV using recursive
   regression-tree segmentation. This is exploratory, not confirmatory.
9. Date-grouped cross-validation comparing a linear model, GAM and regression
   tree for predicting location relative error.

## Reference-location score

For every common, physically eligible location and gas, calculate absolute
bias, MAE, RMSE, Spearman correlation with the block barn mean, temporal SD of
relative error, mutual information and coverage. Metrics are converted to
within-gas percentile ranks and averaged equally. No gas receives a larger
weight unless a later engineering objective justifies it.

The primary selection uses the conservative fan-corridor exclusion. A
sensitivity result ranks all otherwise valid common locations.

## Interpretation boundary

The selected location is a pragmatic representative of concentration, not an
absolute truth and not automatically an optimal emission outlet location. Its
suitability for emission calculations must also consider airflow measurement,
outlet geometry and the intended mass-balance boundary.
