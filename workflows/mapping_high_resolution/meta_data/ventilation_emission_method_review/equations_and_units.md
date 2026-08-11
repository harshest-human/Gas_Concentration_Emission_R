# Equation chain and dimensional checks

## 1. Occupancy and livestock units

\[LU_t=\sum_kN_{k,t}m_k/(500\ \mathrm{kg})\]

Units: animals x kg animal-1 / kg = LU (dimensionless reference unit).

## 2. Milking-cow total heat and CO2 production

\[\Phi_{20}=5.6m^{0.75}+22Y+1.6\times10^{-5}p^3\quad[W\ cow^{-1}]\]
\[f_T=(1000+4(20-t_i))/1000\]
\[P_{CO2,t}=\sum_kN_{k,t}P_{C,k}\Phi_{20,k}f_T/1000\]

Units: HPU = W/(1000 W); (m3 CO2 h-1 HPU-1) x HPU = m3 CO2 h-1. Source: VERA (2018), Annex H, pp. 51-52. Do not use the milking-cow equation for heifers without the category-specific equation.

## 3. Moisture basis

\[x_{dry}=x_{wet}/(1-x_{H2O,wet})\]

Units: mol mol-1 / mol mol-1 = mol mol-1. Apply only if input is wet basis.

## 4. Ventilation

\[Q_t=P_{CO2,t}/(x_{CO2,in,t}-x_{CO2,bg,t})\]
\[Q_t=P_{CO2,t}10^6/(CO2_{in}^{ppm}-CO2_{bg}^{ppm})\]

Units: m3 CO2 h-1 / (m3 CO2 m-3 air) = m3 air h-1.

## 5. Mole fraction to mass concentration

\[C_g[mg\ m^{-3}]=ppm_g\times10^{-6}\times pM_g/(RT)\times10^3\]

Equivalently, when \(M\) is entered in g mol-1, the last combined multiplier is \(10^{-3}\). Units: mol mol-1 x Pa x g mol-1 /(Pa m3 mol-1 K-1 x K) x 1000 mg g-1 = mg m-3. For ppb replace 10-6 by 10-9. Use \(M_{CH4}=16.04246\), \(M_{NH3}=17.03052\), \(M_{CO2}=44.0095\) g mol-1.

## 6. Gas emission

\[E_{g,barn}=Q(C_{g,in}-C_{g,bg})\]

If C is g m-3: m3 h-1 x g m-3 = g h-1. Source: Janke et al. (2020), Eq. 7, p. 19.

## 7. Normalisation and annual rate-equivalent

\[E_{cow}=E_{barn}/N\]
\[E_{LU}=E_{barn}/LU\]
\[E_{annual,LU}=E_{LU}\times8760/1000\]

Units respectively: g h-1 cow-1; g h-1 LU-1; kg year-1 LU-1. The final calculation is not an annual emission factor unless the mean hourly rate is annually representative.

## 8. Uncertainty

\[(u_Q/Q)^2=(u_P/P)^2+(u_{\Delta x}/\Delta x)^2\]
\[u_E^2=(\Delta C)^2u_Q^2+Q^2u_{\Delta C}^2+2Q\Delta C\,cov(Q,\Delta C)\]

Use Monte Carlo for non-linearity, correlations and small concentration enhancements.

## Fictional worked example (not study data)

Assume 50 cows, a deliberately fictional category-mean live mass 625 kg, total tracer production 8.00 m3 CO2 h-1, dry CO2 of 900 ppm inside and 450 ppm incoming, CH4 of 30 and 2 ppm, NH3 of 5.0 and 0.2 ppm, T=293.15 K and p=101325 Pa.

1. \(LU=50\times625/500=62.5\ LU\).
2. \(Q=8.00\times10^6/(900-450)=17,777.8\ m^3 h^{-1}\); \(Q/N=355.6\ m^3 h^{-1}cow^{-1}\); \(Q/LU=284.4\ m^3 h^{-1}LU^{-1}\).
3. Fictional CH4 enhancement is 28 ppm. At stated T,p, \(\Delta C_{CH4}=28\times10^{-6}pM/(RT)=0.01867\ g\ m^{-3}\). Thus \(E_{CH4}=331.9\ g\ h^{-1}\), \(5.31\ g\ h^{-1}LU^{-1}\), and the rate-equivalent is \(46.5\ kg\ year^{-1}LU^{-1}\).
4. Fictional NH3 enhancement is 4.8 ppm, or 0.003399 g m-3; \(E_{NH3}=60.4\ g\ h^{-1}\), \(0.966\ g\ h^{-1}LU^{-1}\), rate-equivalent \(8.46\ kg\ year^{-1}LU^{-1}\).

These values illustrate units only; the annual multipliers are not defensible study estimates.
