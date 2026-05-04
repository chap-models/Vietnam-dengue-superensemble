options(warn = 1)

train_chap <- function(train_fn, model_fn, geojson_fn) {
  # The superensemble is fit fresh on each predict call; no training artifact is saved.
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2) {
  train_chap(
    train_fn   = args[1],
    model_fn   = args[2],
    geojson_fn = if (length(args) >= 3) args[3] else ""
  )
}
