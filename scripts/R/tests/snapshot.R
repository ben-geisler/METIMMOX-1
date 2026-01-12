commandArgs <- function(trailingOnly = TRUE) {
  if (trailingOnly) {
    return(c("18", "fixed"))
  }
}

source("scripts/R/analysis/15_save_snapshot.R")