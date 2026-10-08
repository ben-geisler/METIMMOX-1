#' Mask base::commandArgs() so that 13_save_snapshot.R receives the fixed
#' arguments issue 118 and "baseline" when sourced from this helper.
commandArgs <- function(trailingOnly = TRUE) {
  if (trailingOnly) {
    return(c("118", "baseline"))
  }
}

source("analysis/13_save_snapshot.R")