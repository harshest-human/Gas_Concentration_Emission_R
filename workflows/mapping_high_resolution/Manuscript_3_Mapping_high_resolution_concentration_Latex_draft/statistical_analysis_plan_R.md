# Statistical analysis plan for Manuscript 3

## 1. Scientific focus

The primary objective is to quantify and interpret the vertical, horizontal,
and temporal heterogeneity of CO2, CH4, and NH3 concentrations within the
naturally ventilated dairy barn.

The evaluation of reduced sampling configurations is secondary. It will be
presented as a practical consequence of the observed spatial structure rather
than as the central scientific hypothesis.

Campaign 1 provides the full 51-location, three-height reference dataset.
Campaign 2 provides a longer-term application of the reduced top-and-bottom
design. Campaign 2 is not an independent replication of the full 51-location
experiment.


## 2. Input data

Primary analytical dataset:

```text
D:/Data_Analysis_R/Gas_Concentration_Emission_R/workflows/
mapping_high_resolution/clean_data/manuscript1/
manuscript1_campaign1_2_analytical_data.csv
```

Existing reproducible preparation script:

```text
D:/Data_Analysis_R/Gas_Concentration_Emission_R/workflows/
mapping_high_resolution/scripts/manuscript1_campaign1_2_analysis.R
```

Required variables:

- `DATE.TIME`
- `campaign`
- `analyser`
- `location`
- `horizontal_position`
- `vgroup`
- `CO2`
- `CH4`
- `NH3`
- Campaign 1 corrected concentrations
- `CH4_CO2_pct`
- `NH3_CO2_pct`
- `NH3_CH4_pct`, to be calculated when the denominator is positive

The south-background location `s` will be retained in the analytical data but
excluded from internal spatial models.


## 3. Software and R packages

Analyses will be conducted in R 4.5.3 or later.

Core packages:

```r
library(data.table)
library(ggplot2)
library(lme4)
library(lmerTest)
library(emmeans)
library(mgcv)
library(performance)
library(DHARMa)
library(broom.mixed)
library(parameters)
library(effectsize)
library(boot)
library(DescTools)
```

Optional packages:

```r
library(glmmTMB)   # Alternative model with temporal correlation structure
library(nlme)      # AR(1) correlation where required
library(equivalence)
```

Package versions will be recorded with `sessionInfo()` or `renv`.


## 4. General preprocessing

1. Parse `DATE.TIME` in `Europe/Berlin`.
2. Exclude the south-background location `s` from internal heterogeneity
   models.
3. Retain only positive, physically valid concentrations.
4. Use the CRDS-scale corrected Campaign 1 values:

   ```r
   CO2_corr = CO2_raw / 1.06
   CH4_corr = CH4_raw / 1.06
   NH3_corr = NH3_raw / 1.09
   ```

5. CRDS concentrations and H2O remain unchanged.
6. Define:

   ```r
   height = factor(vgroup, levels = c("bottom", "mid", "top"))
   horizontal_position = factor(horizontal_position, levels = 1:17)
   location = factor(location)
   date = as.Date(DATE.TIME)
   hour = hour(DATE.TIME) + minute(DATE.TIME) / 60
   block_2h = floor_date(DATE.TIME, "2 hours")
   ```

7. Use complete two-hour Campaign 1 blocks containing all 51 internal
   locations for direct configuration comparisons.
8. Do not replace missing measurements or interpolate gas concentrations.


## 5. Transformations

The primary models will use natural-log-transformed concentrations:

```r
log_CO2 = log(CO2)
log_CH4 = log(CH4)
log_NH3 = log(NH3)
```

Ratio responses:

```r
log_CH4_CO2 = log(CH4 / CO2)
log_NH3_CO2 = log(NH3 / CO2)
log_NH3_CH4 = log(NH3 / CH4)
```

The NH3/CH4 ratio will be used to compare the relative spatial behaviour of
the two principal animal- and manure-related gases without CO2 as the
normalising gas.

Model estimates will be back-transformed and reported as geometric-mean ratios
or percentage differences:

```r
percent_difference = 100 * (exp(estimate) - 1)
```


## 6. Descriptive spatial analysis

For each campaign, gas, location, and height:

- Number of observations
- Arithmetic mean
- Median
- Standard deviation
- Interquartile range
- Minimum and maximum
- Coefficient of variation
- 95% confidence interval based on independent daily location means

Campaign 1 maps will display all 51 internal locations.

Colours:

```r
height_colours <- c(
  top = "orange",
  mid = "green3",
  bottom = "steelblue1"
)
```

Primary visualisations:

1. Location means and 95% confidence intervals for each gas.
2. Vertical profiles at each horizontal position.
3. Heat map or barn-layout map of location-level deviations from the barn
   mean.
4. Block-level distributions of spatial heterogeneity.


## 7. Primary mixed-effects models: spatial heterogeneity

### 7.1 Analysis unit

Campaign 1 observations will first be averaged by:

```r
block_2h, location, horizontal_position, height
```

This reduces imbalance caused by sequential analyser switching and provides
comparable spatial records within common two-hour periods.

### 7.2 Model

A separate model will be fitted for each gas:

```r
lmer(
  log_concentration ~ height * horizontal_position +
    (1 | block_2h),
  data = campaign1_complete
)
```

Primary tests:

- Overall height effect
- Overall horizontal-position effect
- Height × horizontal-position interaction

The interaction is the main test of whether the vertical gradient changes
across the barn.

The interaction will be retained when scientifically meaningful even if its
global p-value is above 0.05, because it represents the spatial structure under
investigation.

### 7.3 Inference

- Type III tests with Satterthwaite degrees of freedom
- Estimated marginal means using `emmeans`
- Tukey-adjusted top–middle, top–bottom, and middle–bottom contrasts
- Position-specific vertical contrasts where supported
- Back-transformed percentage differences and 95% confidence intervals

Position-specific comparisons will be treated as secondary and corrected using
the Holm method within each gas.


## 8. Quantification of heterogeneity

Within each complete two-hour block, calculate across locations:

```r
mean_concentration
sd_concentration
cv_concentration
iqr_concentration
range_concentration
sd_log_concentration
```

`sd_log_concentration` will be the primary heterogeneity metric because it is
less dominated by high absolute concentrations and is compatible with
log-normal concentration distributions.

Heterogeneity will be summarised separately for:

- Full 51-location design
- Top only
- Middle only
- Bottom only
- Top + middle
- Middle + bottom
- Top + bottom

Block bootstrap confidence intervals will be generated by resampling complete
two-hour blocks rather than individual rows.


## 9. Ratio models

Separate Campaign 1 mixed models will be fitted for:

```r
log(CH4 / CO2)
log(NH3 / CO2)
log(NH3 / CH4)
```

Model:

```r
lmer(
  log_ratio ~ height * horizontal_position +
    (1 | block_2h),
  data = campaign1_complete
)
```

Interpretation:

- A spatial effect remaining after normalisation to CO2 indicates gas-specific
  behaviour rather than a uniform ventilation-driven dilution effect.
- Spatial variation in NH3/CH4 indicates that NH3 and CH4 do not respond
  proportionally to the same source and transport conditions.
- Results will be reported as percentage differences in the ratio.
- Ratios will not be interpreted as emission ratios because no
  background-corrected enhancement-ratio regression is planned.


## 10. Hour-of-day analysis

Hour-of-day effects will be evaluated using cyclic generalised additive mixed
models rather than treating hour as a linear predictor.

Separate model for each gas:

```r
gam(
  log_concentration ~
    height +
    s(hour, bs = "cc", k = 8) +
    s(hour, by = height, bs = "cc", k = 8) +
    s(horizontal_position, bs = "re") +
    s(date, bs = "re"),
  method = "REML",
  data = campaign_data,
  knots = list(hour = c(0, 24))
)
```

The same model structure may be used for all three log-ratio responses.

Tests and outputs:

- Overall cyclic hour effect
- Height-specific cyclic effects
- Difference curves between heights with simultaneous 95% confidence bands
- Predicted 24-hour profiles for each gas and height

If residual autocorrelation remains material, an AR(1)-capable model in
`gamm()`, `nlme`, or `glmmTMB` will be used.


## 11. Fixed clock-time analysis

The author-defined operational periods will be:

```r
day   = 06:00-17:59
night = 18:00-05:59
```

Times will be interpreted in the `Europe/Berlin` time zone. These categories
will be reported as fixed clock-time periods and will not be described as
solar daylight and darkness.

Separate mixed model for each gas:

```r
lmer(
  log_concentration ~ day_night * height +
    horizontal_position +
    (1 | date),
  data = campaign_data
)
```

Primary quantities:

- Day-versus-night percentage difference
- Day/night × height interaction
- Estimated marginal means and 95% confidence intervals

This analysis is associative. Campaigns 1 and 2 do not contain the complete
animal, wind, and cooling-system covariates required for causal attribution.


## 12. Campaign comparison

Campaign 1 and Campaign 2 differ in calendar period, analyser deployment, and
sampling design. Campaign effects are therefore not interpreted as isolated
seasonal effects.

The persistence of the top–bottom contrast will be evaluated with:

```r
lmer(
  log_concentration ~ campaign * height +
    horizontal_position +
    (1 | date),
  data = top_bottom_common
)
```

Only locations or horizontal positions represented comparably in both
campaigns will be included.

Results will be described as a comparison between measurement periods, not as
a definitive seasonal effect.


## 13. Sampling-configuration analysis

For every complete Campaign 1 two-hour block, compare each reduced
configuration with the full 51-location reference.

Metrics:

- Mean bias
- Mean relative bias
- Median absolute deviation
- RMSE
- 95th percentile of absolute relative deviation
- Coefficient of determination
- Lin's concordance correlation coefficient
- Agreement in `sd_log_concentration`

Uncertainty:

- Block bootstrap 95% confidence intervals
- Resampling unit: complete two-hour block

Equivalence analysis:

- Primary equivalence margin: ±5% relative difference
- Sensitivity margin: ±10%
- Two one-sided tests or a 90% confidence interval for the paired block-level
  relative difference

The top-and-bottom configuration will be considered practically equivalent
only if its entire equivalence interval falls within the prespecified margin.


## 14. Model diagnostics

For every fitted model:

1. Residual-versus-fitted plot
2. Normal Q–Q plot
3. Random-effect Q–Q plot
4. Residual temporal autocorrelation
5. Heteroscedasticity by height and horizontal position
6. Influential block and location diagnostics
7. Singular-fit check
8. Comparison of observed and model-predicted values

Corrective actions:

- Retain log transformation as the default.
- Use variance structures or robust/block-bootstrap inference if
  heteroscedasticity remains.
- Use AR(1) correlation if temporal autocorrelation is material.
- Do not delete observations solely to improve model assumptions.
- Any excluded observation must have a documented physical or measurement
  justification.


## 15. Multiplicity control

- Tukey adjustment for the three pairwise height comparisons
- Holm adjustment for position-specific contrasts within each gas
- Separate inferential families for CO2, CH4, NH3, CH4/CO2, NH3/CO2, and
  NH3/CH4
- Exact p-values will be reported; emphasis will remain on effect sizes and
  confidence intervals

The nominal significance level is:

```r
alpha <- 0.05
```


## 16. Sensitivity analyses

1. Arithmetic concentration versus log concentration.
2. One-hour, two-hour, and four-hour block definitions.
3. Analyses with and without locations 19 and 40.
4. Complete-block analysis versus all available blocks.
5. Campaign 1 corrected FTIR values versus raw values, shown only as a
   measurement-harmonisation sensitivity analysis.
6. Fixed horizontal-position effects versus random horizontal-position
   effects.
7. Day/night classification with and without twilight observations.


## 17. Planned manuscript tables

1. Campaign design and retained analytical records.
2. Descriptive concentrations and ratios by gas and height.
3. Mixed-model tests for height, horizontal position, and their interaction.
4. Back-transformed height contrasts with 95% confidence intervals.
5. Sampling-configuration bias, error, concordance, and equivalence results.
6. Sensitivity-analysis summary.


## 18. Planned manuscript figures

1. Barn plan and cross-sectional sampling layout.
2. Location-level concentration means and 95% confidence intervals.
3. Vertical gradients across the 17 horizontal positions.
4. Spatially resolved CH4/CO2, NH3/CO2, and NH3/CH4 ratios.
5. Predicted cyclic hour-of-day profiles by gas and height.
6. Full-versus-reduced configuration agreement.
7. Distribution of block-level relative error for reduced configurations.


## 19. Reporting principles

- Report estimates, uncertainty, and practical magnitude before p-values.
- Use UK English.
- Distinguish statistical association from mechanism.
- Interpret spatial patterns in relation to barn layout, animal-occupied zones,
  manure handling, openings, and fans only where the measurement design
  supports the interpretation.
- Do not describe Campaign 1–2 differences as purely seasonal.
- Present sampling reduction as an applied outcome of the heterogeneity
  analysis.


## 20. Decisions required before model execution

The following choices must be confirmed before final inferential models are
run:

3. Whether the ±5% equivalence margin is acceptable for all three gases.
4. Whether horizontal positions should be treated primarily as fixed effects
   or as a random sample of possible barn positions.
5. Whether Campaign 2 should be included in hour-of-day modelling or retained
   solely as longer-term confirmation of the top–bottom pattern.
