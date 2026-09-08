# ===============================================================================
# PFS ENDPOINT DERIVATION - SINGLE SOURCE OF TRUTH (issue #149)
# ===============================================================================
# The trial export records two separate things:
#   - "Progression exit" (renamed `Progression` in 02): exit from the study
#     because of radiological progression only, and
#   - "Days until progression" (`PFSwk` after conversion): time to that exit,
#     which for patients who died without a recorded progression is simply the
#     last progression-free assessment.
#
# METIMMOX (Ree et al. 2024, doi:10.1038/s41416-024-02696-6) defines PFS as the
# time to first progression on active therapy OR death from any cause, whichever
# came first. A partitioned survival model also requires this definition: the
# progression-free state means alive AND progression-free, so a PFS curve that
# censors deaths overstates that state and lets OS and PFS cross.
#
# derive_pfs_endpoint() therefore recodes, in place,
#   Progression <- 1 if progressed OR died, else 0
#   PFSwk       <- time to progression if progressed,
#                  time to death       if died without progression,
#                  follow-up time      otherwise (censored)
# and preserves the original trial variables as `ProgressionExit` and `TTPwk`
# (time to progression, deaths censored) for reference. The function is
# idempotent: re-applying it to an already-derived data frame is a no-op.
#
# Called by 02_setup_and_global_variables.R (the analysis pipeline) and by any
# report that reads data/tidy/METIMMOX.rds directly.
# ===============================================================================

#' Derive the composite PFS endpoint (progression or death)
#'
#' @param data Data frame with numeric columns `Progression` (progression-exit
#'   flag), `PFSwk` (time to progression in weeks), `Death` (0/1) and `OSwk`
#'   (time to death or last follow-up in weeks).
#' @return The same data frame with `Progression` and `PFSwk` recoded to the
#'   composite endpoint, plus `ProgressionExit` and `TTPwk` holding the
#'   original values.
#' @export
derive_pfs_endpoint <- function(data) {
  required <- c("Progression", "PFSwk", "Death", "OSwk")
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop("derive_pfs_endpoint(): missing column(s): ",
         paste(missing, collapse = ", "))
  }

  # Preserve the raw trial variables once; on repeat calls reuse them so the
  # derivation is idempotent.
  if (!"ProgressionExit" %in% names(data)) {
    data$ProgressionExit <- as.numeric(data$Progression)
  }
  if (!"TTPwk" %in% names(data)) {
    data$TTPwk <- as.numeric(data$PFSwk)
  }

  progressed <- data$ProgressionExit == 1
  died <- as.numeric(data$Death) == 1

  if (anyNA(progressed) || anyNA(died)) {
    stop("derive_pfs_endpoint(): Progression/Death flags contain NA.")
  }
  if (any(data$TTPwk > data$OSwk, na.rm = TRUE)) {
    stop("derive_pfs_endpoint(): time to progression exceeds OS time for ",
         sum(data$TTPwk > data$OSwk, na.rm = TRUE), " patient(s).")
  }

  data$Progression <- as.numeric(progressed | died)
  data$PFSwk <- ifelse(progressed, data$TTPwk,
                       ifelse(died, as.numeric(data$OSwk), data$TTPwk))

  data
}

message("PFS endpoint helper loaded from pfs_endpoint.R")
