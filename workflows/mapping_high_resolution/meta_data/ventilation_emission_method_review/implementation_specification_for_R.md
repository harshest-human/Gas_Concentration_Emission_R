# Implementation specification for R (no analysis executed)

## Input tables and columns

### Gas observations

`timestamp` (POSIXct Europe/Berlin plus preserved UTC), `campaign`, `analyser`, `sampling_point`, `point_role`, `CO2_ppm`, `CH4_ppm`, `NH3_ppm`, `H2O_volpct`, `gas_basis`, `point_x`, `point_y`, `point_z_m`, `distance_opening_m`, `line_recovery`, `calibration_id`.

### Animals/production

`timestamp_start`, `timestamp_end`, `animal_category`, `n_inside`, `mean_live_mass_kg`, `milk_kg_cow_d`, `pregnancy_day`, `weight_gain_kg_d`, `feed_ME_MJ_kgDM`, `value_source`, `value_uncertainty`.

### Environment/operations

`timestamp`, `T_inside_C`, `RH_inside_pct`, `T_outside_C`, `RH_outside_pct`, `pressure_Pa`, `wind_speed_m_s`, `wind_direction_deg`, `fan_on`, `sprinkler_on`, `manure_scraper_on`, `slurry_event`, `milking_or_cows_out`.

### Metadata

Floor/manure class, control-volume membership, molar masses, analyser accuracy/precision/drift, tube delay, eligible point flags and source citations.

## Processing order

1. Parse source timestamps, convert to UTC for joins, preserve Europe/Berlin clock time and DST fold indicator.
2. Apply analyser-specific flush/delay rules already established upstream; do not re-clean source data here.
3. Join animal and operations states by bounded interval; calculate within-hour time-weighted `n_inside` and animal inputs. Flag stale carries.
4. Determine wet/dry status from metadata. If wet, `x_dry=x_wet/(1-H2O_volpct/100)` independently for every gas/point/time. No correction if instrument already reports dry.
5. Aggregate to hourly point values using duration weighting; require a pre-specified minimum temporal coverage.
6. Assign wind sector and upwind incoming point(s). Calculate spatial primary and sensitivity estimators from eligible points only. Do not silently substitute a missing background.
7. Calculate `LU=sum(n_inside*mean_live_mass_kg)/500`.
8. Calculate category heat and CO2 production with the exact VERA Annex H equation appropriate to category; do not apply the milking-cow equation to heifers. Sum categories and add/encode manure contribution exactly once.
9. `delta_CO2_ppm=CO2_inside_dry_ppm-CO2_background_dry_ppm`; `Q_m3_h=P_CO2_m3_h*1e6/delta_CO2_ppm` for positive deltas. Preserve raw negative/zero deltas, flag and leave Q missing for the physical primary estimate.
10. Convert gas mole fractions to mass at hourly measured pressure and temperature: `C_mg_m3=ppm*1e-6*pressure_Pa*M_g_mol/(R*T_K)*1000`; equivalently with M in g/mol, multiply by `1e-3`.
11. `E_g_h=Q_m3_h*(C_inside_mg_m3-C_bg_mg_m3)/1000`; preserve signed gas enhancement.
12. Normalise: `Q_cow=Q/n_inside`, `Q_LU=Q/LU`, `E_g_h_cow=E/n_inside`, `E_g_h_LU=E/LU`; set missing if denominator <=0.
13. `rate_equivalent_kg_y_LU=E_g_h_LU*8760/1000`, labelled explicitly. For a defensible annual estimate integrate hourly barn emissions and divide by time-integrated LU exposure using pre-specified seasonal/operation weights.

## Required output columns

Keys and provenance; spatial strategy; all input summaries; `LU`, `heat_W`, `temperature_factor`, `activity_factor`, `P_CO2_animals_m3_h`, `P_CO2_manure_m3_h`, `P_CO2_total_m3_h`, dry/wet concentrations, enhancements, `Q_m3_h`, `Q_m3_h_cow`, `Q_m3_h_LU`, CH4/NH3 mg m-3, barn/cow/LU emissions, rate-equivalents, Monte Carlo quantiles, structural strategy range, and every validity flag.

## Filters and validation

- Never use capacity 58 in place of observed hourly occupancy without an explicit imputation flag.
- Never apply an automatic 200 ppm VERA cutoff; flag `<200 ppm` and conduct sensitivity analysis.
- Require same basis and compatible timestamps for inside/background and gases.
- Check `0<=H2O<100`, `pressure>0`, `T_K>0`, `0<=n_inside<=58`, non-negative live mass/production.
- Detect impossible concentration units, duplicate keys, analyser discontinuities, background exceeding all indoor points, Q non-finite/extreme, and unrepresented wind sectors.
- Cross-check tracer-ratio gas emission against `Q*deltaC`; results must agree numerically.
- Compare arithmetic mean, median, wind-paired and candidate-point results; report ratios and rank reversals.

## Monte Carlo

At least 10,000 draws per hour or a convergence-verified smaller number. Draw correlated calibration/drift errors by analyser-day, residual precision by observation, spatial points by hierarchical bootstrap, occupancy as discrete uncertainty, live mass/milk/pregnancy and manure factor from documented distributions only. Recompute basis correction, PCO2, Q and emissions each draw. Save seed, distribution table/version, median, SD, 2.5/97.5 percentiles and invalid-draw share. Do not invent distribution widths; unresolved widths remain required decisions.
