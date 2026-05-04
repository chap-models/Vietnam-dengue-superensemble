# Helper functions for the Vietnam dengue superensemble.
# Derived from FelipeJColon/paper_dengue_superensemble (00_Functions.R).

get.mliks <- function(models) {
  sapply(models, function(m) m$mlik[1, 1])
}

get.dics <- function(models) {
  sapply(models, function(m) m$dic$dic)
}

# Softmax normalisation of log-scores (numerically stable).
reweight <- function(x) {
  x <- x - max(x, na.rm = TRUE)
  w <- exp(x)
  w / sum(w, na.rm = TRUE)
}
