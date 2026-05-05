source("train.R")
source("predict.R")

train_chap(
  train_fn   = "input/historic_data.csv",
  model_fn   = "model",
  geojson_fn = "input/historic_data.geojson"
)

predict_chap(
  model_fn   = "model",
  hist_fn    = "input/historic_data.csv",
  future_fn  = "input/future_data.csv",
  preds_fn   = "input/predictions.csv",
  geojson_fn = "input/historic_data.geojson"
)
