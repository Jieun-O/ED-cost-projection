# Heat-attributable emergency department admission costs in South Korea

Analysis code for the estimation and projection of heat-attributable emergency
department (ED) admission costs across 250 si-gun-gu districts in South Korea
(2010-2023), and projections to 2039 under SSP2-4.5 and SSP5-8.5.

The workflow covers four steps: (1) two-stage distributed lag non-linear
modelling of the temperature-admission association, (2) projection of
attributable admissions and costs under future climate and demographic
scenarios, (3) Shapley decomposition of the projected change into four drivers,
and (4) calculation of attributable costs and attributable fractions.

---

## Scripts

| Script | Purpose | Main output |
|---|---|---|
| `code/01_twostage.R` | District-specific DLNM (time-stratified case-crossover, conditional quasi-Poisson) and multivariate meta-regression | `meta.xlsx`, `result_gnm.RData` |
| `code/02_projection.R` | Attributable admissions and costs by district, period, cause, age group and GCM, with 1,000 Monte Carlo simulations | `proj_{ssp}_mod{a}_{b}_dec{a}_{b}.RData` |
| `code/03_decomposition.R` | GCM ensemble aggregation and Shapley decomposition into temperature, population ageing, population size and cost per admission | `decomp_*.csv`, decomposition figures |
| `code/04_attributable.R` | Attributable costs (AC) and attributable fractions (AF), by cause and overall | `af_*.csv`, main projection figures |

---

## Data

National Health Insurance Service (NHIS) claims data are held
under a data use agreement and cannot be redistributed, so no example data are
provided. All scripts assume the derived datasets below in `data/`, and write
to `output/`. The structures are documented so that the code can be adapted to
equivalent data from other settings.

### 1. `d_ts_cause_all_2010_2023.RData` (object `d_ts`)

Daily district-level time series, 250 districts x 5,113 days (2010-01-01 to
2023-12-31). Used in `01_twostage.R`.

| Column | Type | Description |
|---|---|---|
| `sgg_h_v2019` | integer | District (si-gun-gu) code, 2019 boundaries |
| `date` | Date | Calendar date |
| `tmean` | numeric | Daily mean temperature, ERA5-Land (degrees C) |
| `shum` | numeric | Daily mean specific humidity, NCEP CFSv2 (kg/kg) |
| `rhum`, `ahum` | numeric | Relative and absolute humidity (sensitivity analyses) |
| `holiday` | integer | 1 for public holidays, 0 otherwise |
| `n_icd_{cause}_{age}` | integer | Daily ED admission counts by cause and age group |

`cause` is one of `i` (circulatory), `j` (respiratory), `n` (genitourinary),
`f` (mental and behavioural), `AB` (infectious), `E` (endocrine and metabolic),
`G` (nervous system); ICD-10 primary diagnosis only.
`age` is one of `u65` (<65 years), `6574` (65-74), `o75` (>=75).

### 2. `era5_noaa_kor_1023.RData` (object `weather_list_sgg`)

Observed daily temperature by district, 2010-2023, as a list of data frames
with columns `region` (district code), `date`, `tmean`. Used in
`02_projection.R` to derive district-specific temperature percentiles (knots,
centring value and basis boundaries).

### 3. `tmean_scn_{ssp}_cal_{gcm}.RData` (object `tmean_scn_mod`)

Bias-corrected daily mean temperature by district, 2010-2039, one file per SSP
scenario and per GCM (19 CMIP6 models). Columns: `sgg_h` (district code),
`date`, `tmean`. `{ssp}` is `245`, `370` or `585`. The GCM list is read from
the file names in `data/`.

### 4. `doy_list.RData` (object `doy_list`)

Observed ED admissions and costs by district, cause, month and age group,
averaged over 2010-2023 and expanded to day of year. A list with one element
per cause; each element has columns `sgg_h_v2019`, `icd`, `month`, `doy`,
`n_{age}` (daily admissions), `cost_{age}` (daily cost, KRW),
`cost_pc_{age}` (cost per admission, KRW).

### 5. `doy_mth_kor.xlsx`

Nationwide monthly admissions and costs by cause and age group, used to
replace district-level cost per admission with the national monthly average.
Columns: `month`, `icd`, `sub` (age group), `n`, `cost`.

### 6. `cost_est.RData` (objects `cpcdf`, `cpc_growth_sum_df`)

Projected cost per admission relative to the baseline period, by cause, age
group and period. `cpcdf` holds the ratio columns
`cpc_ratio_{age}_{scenario}`, where `{scenario}` is `lin1023` (log-linear
trend over 2010-2023, main analysis), `linex21` (excluding 2020-2021) or
`damp` (damped trend, phi = 0.95). Costs are deflated to 2023 using the
general consumer price index.

### 7. `pop_est.RData` (object `popestlist`)

District-level population by year (2010-2039) and age group under the low,
medium (`med`, main analysis) and high fertility variants of the Statistics
Korea projections. Each element has columns `sgg_h`, `year`, `scn_dec`
(`pop`, `age`, `pop+age`) and `pop_{age}`. `scn_dec` fixes either the age
structure or the total population size at baseline values, and is used for the
decomposition scenarios.

### 8. `KOSIS_pop_proj.RData` (object `poplist`)

National population by year and age group, observed (`obs`) and projected
(`med`), used as the person-year denominator.

### 9. `sgg_info_bh_v2019.xlsx`

District code key linking climate projection grids to 2019 administrative
boundaries.

---

## Analytical settings

Warm season (May-October) only. Exposure-response modelled with a natural
cubic spline of temperature (internal knots at the 50th and 90th percentiles
of the district-specific distribution) and a lag of 5 days with one internal
knot at equally spaced values on the log scale. Effects are centred at the
75th percentile of the district-specific warm-season distribution, and only
days above that value contribute to the attributable burden; negative
attributable values are truncated to zero.

Periods are 2010-2014, 2015-2019, 2020-2023 (baseline, 14 years pooled) and
2024-2029, 2030-2034, 2035-2039 (projection). Uncertainty is propagated with
1,000 Monte Carlo samples from the multivariate normal distribution of the
pooled coefficients; reported intervals are 2.5th and 97.5th percentiles of the
resulting empirical distribution. Costs are converted from KRW to USD at
1,305.6625 KRW/USD (2023 average).

The four decomposition factors are temperature change (`heat`), population
ageing (`age`), population size (`pop`) and cost per admission (`cpc`). All
2^4 = 16 factor states are projected, and Shapley values are computed as the
weighted average of the marginal contribution of each factor across all
orderings.
