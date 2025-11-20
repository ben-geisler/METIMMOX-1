commandArgs <- function(trailingOnly = TRUE) {
  if (trailingOnly) {
    return(c("67", "fixed"))
  }
}

source("scripts/R/analysis/15_save_snapshot.R")