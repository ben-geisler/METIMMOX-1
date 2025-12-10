commandArgs <- function(trailingOnly = TRUE) {
  if (trailingOnly) {
    return(c("90", "fixed"))
  }
}

source("scripts/R/analysis/15_save_snapshot.R")