# Column mapping (CHAP name → internal name):
#   disease_cases → dengue_cases   (count outcome)
#   population    → population     (offset, unchanged)
#   location      → areaid         (spatial unit)
#   time_period   → tsdatetime     (parsed to Date)

options(warn = 1)

library(INLA)
library(dplyr)
library(tidyr)
library(lubridate)
library(zoo)
library(spdep)
library(sf)
source("lib.R")

predict_chap <- function(model_fn, hist_fn, future_fn, preds_fn, geojson_fn) {

  # ------------------------------------------------------------------
  # 1. Load data
  # ------------------------------------------------------------------
  historic_df <- read.csv(hist_fn, stringsAsFactors = FALSE)
  future_df   <- read.csv(future_fn, stringsAsFactors = FALSE)

  # Multi-month forecasts must be driven by repeated single-month calls from CHAP.
  loc_col       <- intersect(c("location", "areaid"), names(future_df))[1]
  future_counts <- table(future_df[[loc_col]])
  if (any(future_counts != 1L)) {
    bad <- names(future_counts[future_counts != 1L])
    stop("Future data must contain exactly 1 row per province. ",
         "Got multiple rows for: ", paste(bad, collapse = ", "))
  }

  future_df$disease_cases <- NA_integer_
  historic_df$is_future   <- FALSE
  future_df$is_future     <- TRUE

  df <- bind_rows(historic_df, future_df)

  # ------------------------------------------------------------------
  # 2. Rename CHAP columns to internal names
  # ------------------------------------------------------------------
  if ("location" %in% names(df))
    names(df)[names(df) == "location"]      <- "areaid"
  if ("disease_cases" %in% names(df))
    names(df)[names(df) == "disease_cases"] <- "dengue_cases"

  df$tsdatetime <- ymd(paste0(df$time_period, "-01"))

  # ------------------------------------------------------------------
  # 3. Sort (required for rolling means and lags)
  # ------------------------------------------------------------------
  df <- df %>% arrange(areaid, tsdatetime)

  # ------------------------------------------------------------------
  # 4. Derived variable: diurnal temperature range
  # ------------------------------------------------------------------
  df$dtr <- df$maximum_temperature - df$minimum_temperature

  # ------------------------------------------------------------------
  # 5. Rolling means within province (model-specific preprocessing)
  #    3-month window: temperature, humidity, dtr
  #    4-month window: Niño 3.4 (captures ENSO lead time)
  # ------------------------------------------------------------------
  df <- df %>%
    group_by(areaid) %>%
    mutate(
      tmin02   = rollapply(minimum_temperature,          3, mean, fill = NA, align = "right"),
      tmax02   = rollapply(maximum_temperature,          3, mean, fill = NA, align = "right"),
      shum02   = rollapply(specific_surface_humidity,    3, mean, fill = NA, align = "right"),
      dtr02    = rollapply(dtr,                          3, mean, fill = NA, align = "right"),
      nino3403 = rollapply(nino34_anomaly,               4, mean, fill = NA, align = "right")
    ) %>%
    # Back-fill warm-up NAs with the first available value per province.
    fill(tmin02, tmax02, shum02, dtr02, nino3403, .direction = "up") %>%
    ungroup()

  # ------------------------------------------------------------------
  # 6. Seasonal index shift (−6 months): aligns temporal random effects
  #    with Vietnam's epidemiological season.
  # ------------------------------------------------------------------
  df$date2    <- df$tsdatetime %m-% months(6)
  df$ID.year  <- year(df$date2)
  df$ID.month <- month(df$date2)

  # Normalise year index to start at 1
  df$ID.year <- df$ID.year - min(df$ID.year) + 1L

  # ------------------------------------------------------------------
  # 7. Lagged dengue cases
  # ------------------------------------------------------------------
  # LOCF is applied inline to dengue_cases only for lag computation;
  # disease_cases itself is not modified.
  df <- df %>%
    group_by(areaid) %>%
    mutate(dengueL1 = dplyr::lag(zoo::na.locf(dengue_cases, na.rm = FALSE), 1)) %>%
    ungroup()

  df$loglag <- log1p(df$dengueL1)

  # NA loglag means no historic data precedes the forecast row, so the model cannot run.
  if (any(is.na(df$loglag[df$is_future]))) {
    stop("loglag is NA for prediction rows: historic data must cover at least 1 month per province.")
  }

  # ------------------------------------------------------------------
  # 8. Numeric IDs for INLA random effects
  # ------------------------------------------------------------------
  df$areaid   <- factor(df$areaid)
  df$ID.area  <- as.integer(df$areaid)
  df$ID.area1 <- df$ID.area   # replicate index for AR1 group
  df$ID.area2 <- df$ID.area   # replicate index for IID year group
  df$ID.month1 <- as.integer(df$ID.month)
  df$ID.year   <- as.integer(df$ID.year)

  # ------------------------------------------------------------------
  # 9. Spatial adjacency graph from GeoJSON
  #    The GeoJSON must contain a "province" property matching areaid.
  # ------------------------------------------------------------------
  map <- st_read(geojson_fn, quiet = TRUE)

  # Identify province ID column (first non-geometry column)
  id_col <- setdiff(names(map), attr(map, "sf_column"))[1]

  # Re-order GeoJSON rows to match areaid factor levels so poly2nb
  # indices correspond to ID.area integers.
  prov_levels <- levels(df$areaid)
  row_order   <- match(prov_levels, map[[id_col]])
  if (any(is.na(row_order))) {
    stop("GeoJSON is missing provinces: ",
         paste(prov_levels[is.na(row_order)], collapse = ", "))
  }
  map <- map[row_order, ]

  nb      <- poly2nb(map, queen = FALSE)
  adj_tmp <- tempfile(fileext = ".graph")
  nb2INLA(adj_tmp, nb)
  vnm.adj <- adj_tmp

  # ------------------------------------------------------------------
  # 10. Model formulas (identical to Colón-González et al. 2021)
  #     Shared: BYM spatial + IID year-by-province + AR1 month-by-province.
  #     Differ in environmental fixed effects (see README).
  # ------------------------------------------------------------------
  pc3  <- list(prior = "pc.prec", param = c(3, 0.01))
  pc_r <- list(prior = "pc.cor1", param = c(0.5, 0.75))

  f1 <- dengue_cases ~ loglag +
    f(ID.area,  model = "bym", graph = vnm.adj,
      adjust.for.con.comp = FALSE, constr = TRUE, scale.model = TRUE,
      hyper = list(prec.unstruct = pc3, prec.spatial = pc3)) +
    f(ID.year,  model = "iid", hyper = list(prec = pc3),
      group = ID.area2,
      control.group = list(model = "iid", hyper = list(prec = pc3))) +
    f(ID.month1, model = "ar1", hyper = list(prec = pc3, rho = pc_r),
      group = ID.area1,
      control.group = list(model = "iid", hyper = list(prec = pc3))) +
    periurban_landcover + urban_landcover +
    shum02 + wind_speed + dtr02 + nino3403

  f2 <- dengue_cases ~ loglag +
    f(ID.area,  model = "bym", graph = vnm.adj,
      adjust.for.con.comp = FALSE, constr = TRUE, scale.model = TRUE,
      hyper = list(prec.unstruct = pc3, prec.spatial = pc3)) +
    f(ID.year,  model = "iid", hyper = list(prec = pc3),
      group = ID.area2,
      control.group = list(model = "iid", hyper = list(prec = pc3))) +
    f(ID.month1, model = "ar1", hyper = list(prec = pc3, rho = pc_r),
      group = ID.area1,
      control.group = list(model = "iid", hyper = list(prec = pc3))) +
    periurban_landcover + urban_landcover +
    shum02 + dtr02

  f3 <- dengue_cases ~ loglag +
    f(ID.area,  model = "bym", graph = vnm.adj,
      adjust.for.con.comp = FALSE, constr = TRUE, scale.model = TRUE,
      hyper = list(prec.unstruct = pc3, prec.spatial = pc3)) +
    f(ID.year,  model = "iid", hyper = list(prec = pc3),
      group = ID.area2,
      control.group = list(model = "iid", hyper = list(prec = pc3))) +
    f(ID.month1, model = "ar1", hyper = list(prec = pc3, rho = pc_r),
      group = ID.area1,
      control.group = list(model = "iid", hyper = list(prec = pc3))) +
    periurban_landcover + urban_landcover +
    tmin02 + dtr02 + nino3403

  f4 <- dengue_cases ~ loglag +
    f(ID.area,  model = "bym", graph = vnm.adj,
      adjust.for.con.comp = FALSE, constr = TRUE, scale.model = TRUE,
      hyper = list(prec.unstruct = pc3, prec.spatial = pc3)) +
    f(ID.year,  model = "iid", hyper = list(prec = pc3),
      group = ID.area2,
      control.group = list(model = "iid", hyper = list(prec = pc3))) +
    f(ID.month1, model = "ar1", hyper = list(prec = pc3, rho = pc_r),
      group = ID.area1,
      control.group = list(model = "iid", hyper = list(prec = pc3))) +
    periurban_landcover + urban_landcover +
    tmin02 + tmax02

  f5 <- dengue_cases ~ loglag +
    f(ID.area,  model = "bym", graph = vnm.adj,
      adjust.for.con.comp = FALSE, constr = TRUE, scale.model = TRUE,
      hyper = list(prec.unstruct = pc3, prec.spatial = pc3)) +
    f(ID.year,  model = "iid", hyper = list(prec = pc3),
      group = ID.area2,
      control.group = list(model = "iid", hyper = list(prec = pc3))) +
    f(ID.month1, model = "ar1", hyper = list(prec = pc3, rho = pc_r),
      group = ID.area1,
      control.group = list(model = "iid", hyper = list(prec = pc3))) +
    periurban_landcover + urban_landcover +
    wind_speed

  formulas <- list(f1, f2, f3, f4, f5)

  # ------------------------------------------------------------------
  # 11. Fit all candidate models
  # ------------------------------------------------------------------
  fit_model <- function(formula, idx) {
    cat(sprintf("Fitting model %d ...\n", idx))
    result <- tryCatch(
      inla(
        formula,
        family            = "nbinomial",
        data              = df,
        offset            = log(pmax(df$population, 1)),
        control.predictor = list(compute = TRUE, link = 1),
        control.compute   = list(dic = TRUE, config = TRUE,
                                 return.marginals = FALSE),
        control.inla      = list(strategy = "simplified.laplace"),
        verbose           = FALSE,
        safe              = FALSE
      ),
      error = function(e) {
        cat(sprintf("Model %d failed with simplified.laplace, retrying with laplace ...\n", idx))
        inla(
          formula,
          family            = "nbinomial",
          data              = df,
          offset            = log(pmax(df$population, 1)),
          control.predictor = list(compute = TRUE, link = 1),
          control.compute   = list(dic = TRUE, config = TRUE,
                                   return.marginals = FALSE),
          control.inla      = list(strategy = "laplace"),
          verbose           = FALSE,
          safe              = FALSE
        )
      }
    )
    cat(sprintf("Model %d done (DIC = %.1f)\n", idx, result$dic$dic))
    result
  }

  myModels <- mapply(fit_model, formulas, seq_along(formulas), SIMPLIFY = FALSE)

  # ------------------------------------------------------------------
  # 12. BMA weights: 50 % marginal likelihood + 50 % DIC
  # ------------------------------------------------------------------
  mliks   <- get.mliks(myModels)
  dics    <- get.dics(myModels)
  # reweight() is softmax; negate DIC so lower (better fit) gets higher weight.
  weights <- reweight(mliks) * 0.5 + reweight(-dics) * 0.5

  # ------------------------------------------------------------------
  # 13. Posterior predictive samples (1000 total from BMA mixture)
  # ------------------------------------------------------------------
  s        <- 1000L
  idx.pred <- which(df$is_future)
  mpred    <- length(idx.pred)

  n_per_model <- round(weights * s)
  # Adjust for integer rounding so samples sum to exactly s
  rounding_adj <- s - sum(n_per_model)
  n_per_model[which.max(weights)] <- n_per_model[which.max(weights)] + rounding_adj

  y.pred <- matrix(NA_integer_, mpred, 0)

  for (m in seq_along(myModels)) {
    nm <- n_per_model[m]
    if (nm < 1L) next

    xx    <- inla.posterior.sample(nm, myModels[[m]])
    xx.s  <- inla.posterior.sample.eval(
      function(idx.pred) c(theta[1], Predictor[idx.pred]),
      xx,
      idx.pred = idx.pred
    )

    y_m <- apply(xx.s, 2, function(col) {
      as.integer(rnbinom(mpred, mu = exp(col[-1]), size = exp(col[1])))
    })
    if (is.vector(y_m)) y_m <- matrix(y_m, nrow = mpred)
    y.pred <- cbind(y.pred, y_m)
  }

  # ------------------------------------------------------------------
  # 14. Write output
  # ------------------------------------------------------------------
  out <- data.frame(
    time_period = format(df$tsdatetime[idx.pred], "%Y-%m"),
    location    = as.character(df$areaid[idx.pred]),
    y.pred,
    stringsAsFactors = FALSE
  )
  colnames(out) <- c("time_period", "location", paste0("sample_", 0L:(s - 1L)))
  write.csv(out, preds_fn, row.names = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 4) {
  predict_chap(
    model_fn   = args[1],
    hist_fn    = args[2],
    future_fn  = args[3],
    preds_fn   = args[4],
    geojson_fn = if (length(args) >= 5) args[5] else ""
  )
}
