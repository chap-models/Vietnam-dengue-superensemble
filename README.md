# Vietnam Dengue Superensemble

A CHAP integration of the Vietnam dengue superensemble model originally developed by
[Colón-González et al. (2021)](https://github.com/FelipeJColon/paper_dengue_superensemble).

> Colón-González FJ, et al. *Probabilistic seasonal dengue forecasting in Vietnam:
> A modelling study using superensembles.* PLOS Medicine, 2021.
> <https://doi.org/10.1371/journal.pmed.1003542>

## Model description

The model is a Bayesian Model Averaging (BMA) superensemble of five INLA candidate
models. All five share the same spatiotemporal random-effect structure:

- **BYM** (Besag-York-Mollié) spatial random effect per province
- **IID** year-by-province random effect
- **AR1** month-by-province random effect
They differ in which environmental fixed effects are included:

| Model | Fixed covariates (beyond loglag and land cover) |
|---|---|
| f1 | shum02, wind_speed, dtr02, nino3403 |
| f2 | shum02, dtr02 |
| f3 | tmin02, dtr02, nino3403 |
| f4 | tmin02, tmax02 |
| f5 | wind_speed |

All fixed effects beyond land cover and the log-lag term are computed as rolling means
over the combined historic + forecast series before model fitting:

| Covariate | Source variable | Operation | Window |
|---|---|---|---|
| `tmin02` | `minimum_temperature` | mean | 3 months |
| `tmax02` | `maximum_temperature` | mean | 3 months |
| `dtr02` | `maximum_temperature − minimum_temperature` | mean | 3 months |
| `shum02` | `specific_surface_humidity` | mean | 3 months |
| `nino3403` | `nino34_anomaly` | mean | 4 months (captures ENSO lead time) |

`precipitation_amount_per_day` was also smoothed with a 3-month rolling mean (`pre02`) in
the original code, but it does not appear in any of the five model formulas and was
therefore excluded from this integration.

All models include a **seasonal index shift of −6 months**: `tsdatetime` is shifted back
by 6 months before deriving `ID.year` and `ID.month`, so that the full dengue
epidemiological season falls within a single model "year". Without the shift, a season
spanning two calendar years would be split across two `ID.year` levels, undermining the
IID year-by-province random effect. The 6-month value was chosen for Vietnam's seasonal
pattern. **This offset will vary by country** and should be revisited if the model is
applied elsewhere.

BMA weights are computed as a 50/50 average of softmax-normalised marginal likelihoods
and DIC scores.

## Differences from the original code

The original scripts (`04_Fit_models.R` etc.) were a batch analysis pipeline; this adaptation is a CHAP-callable service. The model formulas (f1–f5), all hyperparameters, and all preprocessing transformations are taken directly from the original. The following are intentional departures:

- **INLA strategy** — original uses `strategy="gaussian", int.strategy="eb"`; this adaptation uses `strategy="simplified.laplace"` (with a `laplace` fallback) because the Gaussian strategy segfaults on the INLA devel build in Docker. `simplified.laplace` is a strictly more accurate approximation and does not change the model.
- **f0 and multi-step lag propagation** [up for discussion] — the original fits a baseline model (f0, no climate covariates) iteratively across six leads to propagate `loglag` forward: predicted counts from lead *t* become the lag input for lead *t+1*. This is necessary when the full six-month forecast is produced in one batch. Here, CHAP drives the loop itself by calling `predict_chap` once per month and appending the previous prediction to the historic data before the next call. The lag is therefore always available as the last observed (or last predicted) row, and f0 is not needed.
- **Rolling-mean warm-up NAs** [this should be checked] — the original fills forecast covariates from a pre-computed ensemble-member table. Here, the first 2–3 NA rows per province are back-filled with the first available value (`fill(.direction="up")`), since CHAP provides a single forecast member with no overlap.
- **BMA DIC weighting** [this must also be checked in detail] — the original applies `reweight(dics)` (softmax on raw DIC values, which would give higher weight to *worse* models). This adaptation uses `reweight(-dics)` so that lower DIC correctly receives higher weight.
- **Output format** — the original produces weighted summary statistics and epidemic-threshold metrics. This adaptation draws 1000 posterior predictive samples as required by CHAP.

## Climate forecast data

The original paper used seasonal climate forecasts from **GloSea5** (UK Met Office
Global Seasonal forecasting system, version 5), with 42 ensemble members at
0.83° × 0.56° resolution. For the forecast period, rolling-mean covariates
(`tmin02`, `tmax02`, `shum02`, `dtr02`, `nino3403`) are computed from these forecast
values rather than from observed climate data. The original code loads one file per
ensemble member and processes them in parallel, which is why the preprocessing
has an `ensmember` dimension.

In this CHAP integration, CHAP supplies the future climate data directly as a single
`future_data.csv` — the ensemble averaging or member selection is handled upstream by
CHAP before the model is called.

## How the original code handles multi-month forecasts (f0)

The original pipeline produces all six forecast months in a single batch run. Because
`loglag` (log of the previous month's case count) is a fixed covariate in all models,
it needs a plausible value for every forecast month before f1–f5 can be fitted. The
original solves this with a baseline model **f0** — identical random-effect structure to
f1–f5 but no climate covariates — fitted iteratively six times:

1. Fit f0 on data up to lead 1; extract the posterior mean for that forecast row as `lag1` for lead 2.
2. Fit f0 on data up to lead 2; extract the posterior mean as `lag1` for lead 3.
3. …repeat through lead 6.

Once all six `loglag` values have been chained together this way, f1–f5 are each fitted
once on the full six-month window with those propagated lags already in place.

In this CHAP integration f0 is not needed because CHAP calls `predict_chap` one month
at a time and appends the previous prediction to the historic data before each call, so
the lag is always the last row of the observed series.

## Forecast horizon

This implementation is designed to forecast **1 month ahead**. CHAP loops it to produce
multi-month forecasts by calling `predict_chap` repeatedly, each time with updated
historic data that includes the previous prediction as the most recent observation.

## Required covariates

| CHAP column | Description |
|---|---|
| `population` | Province population (used as offset) |
| `minimum_temperature` | Monthly minimum temperature (°C) |
| `maximum_temperature` | Monthly maximum temperature (°C) |
| `nino34_anomaly` | Niño 3.4 SST anomaly |
| `specific_surface_humidity` | Specific surface humidity (kg/kg) |
| `wind_speed` | Mean wind speed (m/s) |
| `periurban_landcover` | Fraction of peri-urban land cover |
| `urban_landcover` | Fraction of urban land cover |

A GeoJSON file containing province boundaries is also required and is passed as the
`geojson` argument to both `train` and `predict` entry points.

## Known warnings

INLA may emit warnings of the form:

```
*** WARNING *** GMRFLib_2order_approx: rescue NAN/INF values in logl
```

These occur when the log-likelihood returns NaN or Inf for extreme hyperparameter
configurations during INLA's Laplace approximation. The immediate cause is NAs in the
rolling-mean covariates (`tmin02`, `shum02`, etc.) that arise during the warm-up period
at the start of each province's time series — any NA covariate makes the linear
predictor NA, which propagates to the likelihood. The preprocessing step back-fills these warm-up NAs with the first available rolling
mean per province before passing data to INLA, so the warnings should not appear in
normal use. An alternative is to drop the warm-up rows entirely (filter out any row
where a rolling-mean covariate is NA after computing the windows); this loses the first
3–4 months of training data per province but avoids any imputation. If they do appear on real data, they are non-fatal: INLA substitutes a
large negative value and continues, and the final predictions remain valid.

Similarly, `vb.correction` divergence warnings are non-fatal. They indicate that INLA's
variational Bayes correction step was skipped and the standard Laplace approximation
was used instead.

## Local development

```bash
# From the repo root
Rscript isolated_run.R
```

This runs against the synthetic five-province example data in `example_data/` and
writes predictions to `example_data/predictions.csv`.
