commandArgs <- function(trailingOnly = TRUE) {
  if (trailingOnly) {
    return(c("118", "baseline"))
  }
}

source("analysis/13_save_snapshot.R")