commandArgs <- function(trailingOnly = TRUE) {
  if (trailingOnly) {
    return(c("71", "baseline"))
  }
}

source("scripts/R/analysis/15_save_snapshot.R")