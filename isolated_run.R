source("train.R")
source("predict.R")

train_chap(
  train_fn   = "example_data/historic_data.csv",
  model_fn   = "model",
  geojson_fn = "example_data/historic_data.geojson"
)

predict_chap(
  model_fn   = "model",
  hist_fn    = "example_data/historic_data.csv",
  future_fn  = "example_data/future_data.csv",
  preds_fn   = "example_data/predictions.csv",
  geojson_fn = "example_data/historic_data.geojson"
)
