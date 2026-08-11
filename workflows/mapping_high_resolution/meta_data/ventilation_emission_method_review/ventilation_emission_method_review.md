# Methodological review: ventilation and gaseous emissions from the Groß Kreutz dairy barn

## Scope and hierarchy of evidence

This review defines a reproducible method; it does not calculate study results. Evidence is classified as: **requirement** (VERA/standard), **recommended method** (peer-reviewed or protocol guidance), **default** (only if farm data are unavailable), or **study decision** (must be fixed before analysis). VERA v3:2018 is the principal measurement protocol. CIGR (2002) and Pedersen et al. (2008) define animal heat/CO2 production. Janke et al. (2020) provides the most directly relevant spatial-sampling evidence for a naturally ventilated dairy barn. IPCC and EMEP/EEA are inventory comparators, not substitutes for the measured hourly method.

## A. Livestock unit (LU/GV)

Use the German/KTBL definition **1 GV = 500 kg live mass**. KTBL's EmiDaT report states that barn live mass is the sum-product of animal number and category live mass, divided by 500 kg (KTBL, 2021, p. 35). Janke et al. define 1 LU identically (2020, p. 19, Eq. 8). Thus, at hour \(t\):

\[
LU_t=\frac{\sum_k N_{k,t}\bar m_{k,t}}{500\ \mathrm{kg}}.
\]

For a homogeneous milking-cow group with no individual weights, use \(LU_t=N_t\bar m/500\), where \(\bar m\) is a documented herd/category mean. VERA's **650 kg** milking-cow value is a protocol spreadsheet default, not a measured fact (VERA, Annex H, pp. 51-52). Report results using measured/category live mass if available and a sensitivity analysis around any assumed mass. Do not equate one cow with one LU.

Changing animal numbers must enter both total CO2 production and the hourly denominator. Use time-weighted occupancy within each hour when movements occur; do not use the licensed capacity of 58 unless it equals observed occupancy. Set animal-derived Q and per-LU output to missing when no cows are inside.

## B. Ventilation rate by CO2 mass balance

For a well-mixed steady control volume:

\[
Q_t=\frac{P_{CO_2,t}}{x_{CO_2,in,t}-x_{CO_2,bg,t}},
\]

where \(Q\) is m3 air h-1, \(P_{CO_2}\) is m3 CO2 h-1 on the same gas-volume reference basis, and \(x\) is dimensionless volume fraction (m3 m-3). With ppm, \(x=ppm\times10^{-6}\), hence

\[
Q_t=\frac{P_{CO_2,t}\,10^6}{CO_{2,in,t}^{ppm}-CO_{2,bg,t}^{ppm}}.
\]

The equivalent VERA tracer-ratio equation is \(E_g=P_{CO_2}\Delta C_g/\Delta C_{CO_2}\) (VERA, s. 7.4.2.3, p. 27). If both production and concentration differences are volume quantities on a consistent basis, conversion to CO2 mass concentration is unnecessary. If mass concentration is used, convert the production term to mass on exactly the same temperature, pressure and moisture basis.

Concentrations must be consistent. Preferred implementation: convert analyser-reported wet mole fractions to dry mole fractions,

\[
x_{dry}=\frac{x_{wet}}{1-x_{H_2O,wet}},
\]

using simultaneous H2O at the same sampling point; if the analyser already reports dry concentrations, do not correct again. Record this from analyser documentation. The ideal-gas conversion is

\[
C_g[\mathrm{mg\,m^{-3}}]=x_g[\mathrm{ppm}]\frac{pM_g}{RT}\times10^{-3},
\]

with \(p\) Pa, \(T\) K, \(M\) g mol-1 and \(R=8.314462618\) J mol-1 K-1. This uses actual local T and p; a standard-condition conversion is permissible only if both flow and concentrations are referred to that same condition.

### Spatial representation

VERA recommends at least one inside point per 10 m barn length, controlled flow/dust filters, points not affected by obstacles, and open-path length averages where available (Annex F, pp. 46-47). For cross-flow, it recommends the middle of the barn at >=3 m; points near outlets should be >=2 m from openings. Incoming air should normally be sampled >=5 m from the barn; nearby sources require additional points. With rapidly varying wind, simultaneous measurement on both sides is recommended and the lower side concentration identifies inlet air (p. 47).

Janke et al. (2020) found large Q and NH3 differences among plausible strategies and showed that including upwind low indoor concentrations can inflate Q. Therefore one point must not be declared representative without validation. **Proposed study method:** calculate a primary robust spatial estimator (pre-specified median of eligible internal points) and sensitivity estimates using (i) arithmetic mean, (ii) wind-sector inlet/outlet pairing, and (iii) candidate representative points. Because a median is not specifically prescribed by VERA, label it a study decision. Use synchronous or tightly time-aligned concentration summaries; never combine measurements from substantially different air-flow states.

### Background, thresholds and intervals

Use measured incoming air, not a generic atmospheric constant. Select upwind/background locations using wind direction and exclude locally contaminated outside points only by a pre-specified rule. VERA explicitly states that a minimum \(\Delta CO_2=200\) ppm **is not required**, because removing low differences systematically underestimates Q and emissions (Annex F, p. 47). Nevertheless, uncertainty explodes as \(\Delta CO_2\to0\). Retain all positive enhancements meeting instrument-specific signal-to-noise criteria, flag values below 200 ppm, and report sensitivity with and without such observations; do not call 200 ppm a VERA validity requirement.

VERA Annex H uses daily-mean heat production and 24 h mean concentrations for its annual test calculation. Janke et al. calculated hourly values. For this high-resolution study, the recommended output interval is hourly because animal count, wind and operation vary within day; derive it from time-weighted, aligned observations and additionally aggregate valid hourly results to daily/seasonal summaries.

Total tracer production must include animals plus manure. VERA incorporates a 10% manure-pit contribution through \(P_C=0.20\) m3 CO2 h-1 HPU-1 for partly slatted floors versus 0.18 for closed floors, but requires measurement or artificial tracer where deep litter, indoor storage, flushing or reduced slurry surfaces make the default unsuitable (Annex H, pp. 51-52). Combustion, humans, calves and connected spaces are separate CO2 sources and must be quantified or affected hours excluded. This barn's manure pit and calf area therefore prevent blind use of an animal-only balance.

## C. Animal CO2 production

VERA Annex H gives for a milking cow:

\[
\Phi_{20}=5.6m^{0.75}+22Y_1+1.6\times10^{-5}p^3\quad[\mathrm W],
\]
\[
f_T=\frac{1000+4(20-t_i)}{1000},\qquad
P_{CO_2,cow}=P_C\frac{\Phi_{20}f_T}{1000}\quad[\mathrm{m^3\,h^{-1}}].
\]

Here \(m\) is kg live mass, \(Y_1\) kg milk cow-1 d-1, \(p\) days pregnant, and \(t_i\) degrees C. \(P_C=0.18\) (closed floor; animal/house contribution stated by VERA) or 0.20 (partly slatted; includes 10% manure-pit contribution) m3 CO2 h-1 HPU-1. VERA's spreadsheet equations use 0.20; floor/manure-system choice is unresolved for Groß Kreutz. The 650 kg and 160 d pregnancy values are defaults only. Sum category-specific production multiplied by actual hourly category counts.

Janke et al. (2020, Eqs. 2-6) instead used 0.185 m3 CO2 h-1 HPU-1 and an hourly activity factor \(A\); its heat equation is the same maintenance + milk + pregnancy form. VERA says activity is not included in CIGR rules and an activity correction is optional (Table 9, p. 28). The activity amplitude and minimum-time parameters must not be invented: obtain them from the applicable Pedersen/CIGR table or omit activity in the primary VERA implementation and treat the Janke activity model as sensitivity analysis.

Total heat production (W) is metabolic heat; sensible heat is only the convective/radiative fraction and is not used to infer respiratory CO2. HPU is 1000 W total heat at 20 degrees C. The CO2/heat factor embeds an assumed respiratory quotient and empirical animal-house relation; do not multiply it by a second respiratory quotient.

Required measured inputs are hourly occupancy/category, indoor temperature, live mass/category mean, milk yield and pregnancy day. Feed energy, dry-matter intake and growth are required for heifer equations or alternative energy-balance models, but not the VERA milking-cow equation. Farm-specific values take precedence. Defaults must be propagated as uncertain inputs and reported separately.

## D. CH4 and NH3 emissions

After Q is estimated:

\[
E_{g,t}=Q_t(C_{g,in,t}-C_{g,bg,t}),
\]

where Q is m3 h-1 and concentrations are g m-3, yielding g h-1 (Janke et al., 2020, Eq. 7). Molar masses: CH4 = 16.04246 g mol-1, NH3 = 17.03052 g mol-1 and CO2 = 44.0095 g mol-1. Use the ideal-gas formula above; ppb is identical with \(10^{-9}\) rather than \(10^{-6}\).

The same spatial estimator/pairing must be used consistently for tracer and pollutant, because VERA requires tracer and pollutant concentrations at the same sampling points (s. 7.4.2.3, p. 27). Prefer wind-dependent inlet/outlet pairing when simultaneous perimeter data support it; otherwise use the validated spatial estimator and quantify between-strategy uncertainty. Negative enhancements are physically possible as noise, advection or contaminated background. Preserve signed values and a flag; do not truncate to zero in the primary data. Exclude only on a documented instrument/background failure rule and report the negative proportion.

Report whole-barn net CH4 and NH3. It includes all sources within the control volume (enteric animals plus manure for CH4; principally manure/urine surfaces for NH3) and must not be labelled enteric-only or manure-only without source partitioning.

## E. Scaling and reporting

\[
E_{g,cow}=E_{g,barn}/N_t,
\]
\[
E_{g,LU}=E_{g,barn}/LU_t=E_{g,barn}\frac{500}{\sum N_tm_t},
\]
\[
E_{g,annual}=E_{g,LU}[\mathrm{g\,h^{-1}\,LU^{-1}}]\frac{8760}{1000}
\quad[\mathrm{kg\,year^{-1}\,LU^{-1}}].
\]

The last equation is a rate-equivalent only. It is scientifically defensible as an annual emission factor only for a probability/sample-weighted annual mean covering seasons, diurnal states, management, occupancy and missingness. Simple multiplication of a short campaign mean by 8760 is not an annual inventory estimate. Prefer integration of hourly barn emissions over observed time plus pre-specified seasonal/operational weights. When cows are outside, remove their CO2 production and live mass from the barn balance; do not calculate a barn animal-normalised rate at zero occupancy. Report barn occupancy hours and distinguish barn-only from whole-farm annual emissions.

## F. IPCC comparison

IPCC Volume 4 Chapter 10 addresses inventory categories **3.A.1 Enteric Fermentation (cattle)** and **3.A.2 Manure Management (CH4)** using annual animal populations, gross-energy/volatile-solids and management-system factors. These factors benchmark annual source-specific CH4, not short-term whole-barn flux. Compare only after matching year, animal class, mass/production and boundary; never insert an IPCC factor into the CO2 balance.

NH3 is an air pollutant rather than an IPCC greenhouse-gas inventory species. Use EMEP/EEA Guidebook 2023, category 3.B Manure Management, or the German inventory method for annual NH3 benchmarks. Inventory housing factors may exclude/include storage and land spreading differently from the measured barn control volume; state boundaries side-by-side.

## G. QA/QC and uncertainty

VERA requires prevention of adsorption, condensation, leaks and blockage; delay/rise/drying times and cross-sensitivities; operation away from detection limits; calibration/maintenance consistent with ISO/IEC 17025; and pre-test sampling-line recovery checks with certified gas (ss. 7.4.2.1 and Annex G, pp. 26, 50-51).

Required QA/QC: certified zero/span and drift records for every analyser; inter-analyser collocation; line recovery/leak tests; flushing/response-time verification; simultaneous H2O basis; pressure/T records; clock/time-zone audit including DST; temporal coverage; background plausibility and upwind status; animal-count freshness; fan/sprinkler/manure/milking flags; sampling-point eligibility; and removal reasons retained without deleting raw values.

For independent inputs to \(Q=P/\Delta x\):

\[
\left(\frac{u_Q}{Q}\right)^2=\left(\frac{u_P}{P}\right)^2+
\left(\frac{u_{\Delta x}}{\Delta x}\right)^2,
\quad
u_{\Delta x}^2=u_{in}^2+u_{bg}^2-2\,cov(in,bg).
\]

For \(E=Q\Delta C_g\), include covariance because Q and emissions share measurements:

\[
u_E^2=(\Delta C_g)^2u_Q^2+Q^2u_{\Delta C_g}^2+2Q\Delta C_g\,cov(Q,\Delta C_g).
\]

Recommended Monte Carlo: for each hour draw correlated analyser errors/drift, spatial locations (hierarchical bootstrap), animal count, live mass, milk/pregnancy/default parameters, manure factor and T/p/H2O; recompute the entire chain; retain signed draws; summarise median, 2.5/97.5 percentiles and invalid-draw proportion. Add between-sampling-strategy variability as a structural uncertainty, not merely instrument error.

Retain flags: `valid_timestamp`, `coverage_ok`, `calibration_ok`, `line_recovery_ok`, `analyser_collocation_ok`, `basis_known`, `background_available`, `background_upwind`, `spatial_coverage_ok`, `animal_count_fresh`, `delta_co2_positive`, `delta_co2_lt_200ppm`, `delta_gas_negative`, `nonanimal_co2_event`, `fan_on`, `sprinkler_on`, `manure_event`, `milking_or_cows_out`, `steady_state_proxy_ok`, `q_plausible`, and `primary_valid`.

## Strongest recommended study method

Use VERA's CIGR/Pedersen CO2-production chain with hourly observed occupancy and measured production inputs; explicitly model the barn manure contribution; calculate hourly Q from dry-basis, simultaneous, spatially representative indoor and true incoming-air CO2; calculate signed whole-barn net CH4/NH3 by mass balance; normalise by hourly live mass/500 kg; and quantify structural uncertainty across spatial strategies. Annualise only using seasonally and operationally representative weighted hours.

## Principal references

- VERA (2018). *VERA Test Protocol for Livestock Housing and Management Systems, Version 3:2018-09*. https://www.vera-verification.eu/app/uploads/sites/9/2019/05/VERA_Testprotocol_Housing_v3_2018.pdf
- Pedersen, S. & Sällvik, K. (eds) (2002). *Heat and moisture production at animal and house levels*. CIGR Section II, 4th report. https://www.cigr.org/sites/default/files/documets/CIGR_4TH_WORK_GR.pdf
- Pedersen, S. et al. (2008). Carbon dioxide production in animal houses: a literature review. *Agricultural Engineering International*, X, Manuscript BC 08 008.
- Janke, D. et al. (2020). Calculation of ventilation rates and ammonia emissions: Comparison of sampling strategies for a naturally ventilated dairy barn. *Biosystems Engineering*, 198, 15-30. https://doi.org/10.1016/j.biosystemseng.2020.07.004
- Schmithausen, A.J. et al. (2018). Quantification of methane and ammonia emissions in a naturally ventilated barn. *Animals*, 8, 75. https://doi.org/10.3390/ani8050075
- Saha, C.K. et al. (2013). The effect of external wind speed and direction... *Biosystems Engineering*, 114, 267-278. https://doi.org/10.1016/j.biosystemseng.2012.12.002
- KTBL (2021). *Aktuelle rechtliche Rahmenbedingungen für die Tierhaltung 2021*, p. 35. https://www.ktbl.de/fileadmin/user_upload/Allgemeines/Download/Tagungen_2021/ARR/ARR_2021.pdf
- IPCC (2019). *2019 Refinement*, Vol. 4, Ch. 10. https://www.ipcc-nggip.iges.or.jp/public/2019rf/pdf/4_Volume4/19R_V4_Ch10_Livestock.pdf
- EMEP/EEA (2023). *Air pollutant emission inventory guidebook 2023*, 3.B Manure Management. https://www.eea.europa.eu/en/analysis/publications/emep-eea-guidebook-2023
