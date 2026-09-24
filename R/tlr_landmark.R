# ===============================================================================
# TLR LANDMARK COHORTS - SINGLE SOURCE OF TRUTH (issue #155)
# ===============================================================================
# Target lesion reduction (TLR) is read at the first on-treatment CT, so any
# analysis that stratifies survival by TLR from randomisation conditions on
# having survived, progression-free, to that scan (guarantee-time bias). The
# landmark analyses in clinical_effectiveness.qmd, clin_effect_figure1.qmd and
# the DAG association reports restrict to patients still at risk at a landmark
# and measure time from the landmark onward.
#
# Three landmark definitions are provided:
#   week9  - fixed landmark at week 9 (primary, retained from earlier versions)
#   scan   - per-patient landmark at the patient's own first on-treatment CT
#            (CT1wk, derived in 02_setup_and_global_variables.R)
#   week12 - fixed landmark at week 12, after the latest first scan (week 12.1)
#
# A patient enters an endpoint-specific cohort only if the endpoint time is
# strictly greater than the landmark (event-free and uncensored AT the
# landmark). Under the per-patient landmark this excludes the patients whose
# progression is recorded at the first scan itself (PFS time equals the scan
# time). Under the fixed week-9 landmark some of those first-scan progressors
# survive the cut (their scan fell after week 9) and some patients have their
# TLR read after the landmark; build_tlr_landmark_cohorts() counts both so the
# reports can state them.
# ===============================================================================

TLR_LANDMARK_SPECS <- list(
  week9 = list(
    id = "week9", type = "fixed", week = 9,
    label = "Fixed week 9 (primary)",
    short_label = "Week 9"
  ),
  scan = list(
    id = "scan", type = "scan", week = NA_real_,
    label = "Per-patient first on-treatment CT date",
    short_label = "Scan date"
  ),
  week12 = list(
    id = "week12", type = "fixed", week = 12,
    label = "Fixed week 12",
    short_label = "Week 12"
  )
)

#' Landmark specification by id
#'
#' @param landmark One of "week9", "scan", "week12".
tlr_landmark_spec <- function(landmark = "week9") {
  landmark <- match.arg(landmark, names(TLR_LANDMARK_SPECS))
  TLR_LANDMARK_SPECS[[landmark]]
}

#' Landmark time (weeks from randomisation) for every row of a TLR cohort
#'
#' @param df TLR-complete data frame; needs `CT1wk` for the per-patient landmark.
#' @param spec A landmark specification from `tlr_landmark_spec()`.
#' @return Numeric vector of landmark weeks, one per row.
tlr_landmark_time <- function(df, spec) {
  if (spec$type == "fixed") {
    return(rep(spec$week, nrow(df)))
  }
  if (!"CT1wk" %in% names(df)) {
    stop("tlr_landmark_time(): the per-patient landmark needs the CT1wk column ",
         "(first on-treatment CT week, from 02_setup_and_global_variables.R).")
  }
  if (anyNA(df$CT1wk)) {
    stop("tlr_landmark_time(): CT1wk is NA for ", sum(is.na(df$CT1wk)),
         " row(s); TLR cohorts must be restricted to patients with a first ",
         "on-treatment scan before building the per-patient landmark.")
  }
  as.numeric(df$CT1wk)
}

#' Build endpoint-specific TLR landmark cohorts
#'
#' @param df TLR-complete data frame with `OSwk`, `Death`, `PFSwk`,
#'   `Progression` and a TLR indicator (`tlr_num` or `tlr`); `CT1wk` is
#'   required for the per-patient landmark and for the first-scan diagnostics.
#'   If `tlr_num` and `Rx_num` are present the ridge product term `tlr_x_rx`
#'   is added to both cohorts.
#' @param landmark One of "week9" (default), "scan", "week12".
#' @return A list with the specification (`spec`), the OS cohort (`os`, with
#'   `lm_wk` and `OSwk_lm`), the PFS cohort (`pfs`, with `lm_wk` and
#'   `PFSwk_lm`), an `attrition` table and a `diagnostics` list.
build_tlr_landmark_cohorts <- function(df, landmark = "week9") {
  spec <- tlr_landmark_spec(landmark)
  required <- c("OSwk", "Death", "PFSwk", "Progression", "ProgressionExit")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop("build_tlr_landmark_cohorts(): missing column(s): ",
         paste(missing, collapse = ", "))
  }
  tlr_col <- if ("tlr_num" %in% names(df)) "tlr_num" else "tlr"
  if (!tlr_col %in% names(df) || anyNA(df[[tlr_col]])) {
    stop("build_tlr_landmark_cohorts(): df must be the TLR-complete cohort ",
         "(no NA in ", tlr_col, ").")
  }

  lm_wk <- tlr_landmark_time(df, spec)
  at_risk_os <- df$OSwk > lm_wk
  at_risk_pfs <- df$PFSwk > lm_wk

  os <- df[at_risk_os, , drop = FALSE]
  os$lm_wk <- lm_wk[at_risk_os]
  os$OSwk_lm <- os$OSwk - os$lm_wk

  pfs <- df[at_risk_pfs, , drop = FALSE]
  pfs$lm_wk <- lm_wk[at_risk_pfs]
  pfs$PFSwk_lm <- pfs$PFSwk - pfs$lm_wk

  if (all(c("tlr_num", "Rx_num") %in% names(df))) {
    os$tlr_x_rx <- os$tlr_num * os$Rx_num
    pfs$tlr_x_rx <- pfs$tlr_num * pfs$Rx_num
  }

  # First-scan progressors: progression recorded at the first on-treatment CT
  # (PFS time equals the scan time). They are TLR-negative by construction
  # because progression and TLR are read from the same scan.
  has_scan <- "CT1wk" %in% names(df) && !anyNA(df$CT1wk)
  first_scan_prog <- if (has_scan) {
    df$Progression == 1 & abs(df$PFSwk - df$CT1wk) < 1e-8
  } else {
    rep(NA, nrow(df))
  }
  scan_after_lm <- if (has_scan) df$CT1wk > lm_wk else rep(NA, nrow(df))

  pfs_excl <- df[!at_risk_pfs, , drop = FALSE]
  os_excl <- df[!at_risk_os, , drop = FALSE]

  diagnostics <- list(
    landmark = spec$id,
    n_total = nrow(df),
    n_os = nrow(os),
    n_pfs = nrow(pfs),
    n_os_excluded = nrow(os_excl),
    n_pfs_excluded = nrow(pfs_excl),
    n_os_excl_deaths = sum(os_excl$Death == 1),
    n_pfs_excl_progression_exit = sum(pfs_excl$ProgressionExit == 1),
    n_pfs_excl_progression_exit_tlr_neg = sum(pfs_excl$ProgressionExit == 1 &
                                          pfs_excl[[tlr_col]] == 0),
    n_pfs_excl_death_pf = sum(pfs_excl$Progression == 1 & pfs_excl$ProgressionExit == 0),
    n_pfs_excl_censored = sum(pfs_excl$Progression == 0),
    n_first_scan_progressors = if (has_scan) sum(first_scan_prog) else NA_integer_,
    n_first_scan_progressors_kept = if (has_scan) {
      sum(first_scan_prog & at_risk_pfs)
    } else NA_integer_,
    n_scan_after_landmark_kept = if (has_scan) {
      sum(scan_after_lm & at_risk_pfs)
    } else NA_integer_,
    min_pfs_lm = if (nrow(pfs) > 0) min(pfs$PFSwk_lm) else NA_real_
  )

  attrition <- data.frame(
    Endpoint = c("Overall survival", "Progression-free survival"),
    N_total = c(nrow(df), nrow(df)),
    N_landmark = c(nrow(os), nrow(pfs)),
    N_excluded = c(nrow(os_excl), nrow(pfs_excl)),
    Events_landmark = c(sum(os$Death), sum(pfs$Progression)),
    stringsAsFactors = FALSE
  )

  list(spec = spec, landmark_wk = lm_wk, os = os, pfs = pfs,
       attrition = attrition, diagnostics = diagnostics)
}

#' Build all three landmark cohorts
#'
#' @inheritParams build_tlr_landmark_cohorts
#' @return Named list (week9, scan, week12) of `build_tlr_landmark_cohorts()` results.
build_all_tlr_landmark_cohorts <- function(df) {
  lapply(TLR_LANDMARK_SPECS, function(spec) build_tlr_landmark_cohorts(df, spec$id))
}

# ===============================================================================
# SHARED BASE COHORT FOR EVERY TLR LANDMARK ANALYSIS (issue #184)
# ===============================================================================
# TLR (target lesion reduction) is positive when the sum of target-lesion
# diameters at the first on-treatment CT is at least 10% below baseline
# (TLRcat: ratio <= 0.9 in 01_data_prep.R). The landmark analyses of
# clinical_effectiveness.qmd, clin_effect_figure1.qmd and the DAG association
# reports all start from the TLR-classified patients WITHIN the complete-case
# cohort of the clinical and economic survival models, so the cohort (not only
# the landmark rule) is shared. Patients with a TLR value but a missing
# TMB/BRAF result are therefore not in any landmark cohort. Non-landmark TLR
# analyses (for example the logistic models into TLR in the DAG tests) may use
# every patient with an observed TLR; they do not call this helper.

#' Complete-case variables of the clinical and economic survival models
TLR_BASE_COHORT_VARS <- c("Age", "sex", "Rx", "crp", "tmb_braf",
                          "OSwk", "Death", "PFSwk", "Progression")

#' TLR-classified patients within the complete-case cohort
#'
#' @param data Data frame produced by scripts 02 and 03 (or a subset of it).
#' @param tlr_col Name of the TLR indicator column (default `tlr`).
#' @return The rows of `data` that are complete on `TLR_BASE_COHORT_VARS` and
#'   have an observed TLR, in their original order; pass the result to
#'   `build_tlr_landmark_cohorts()` or `build_all_tlr_landmark_cohorts()`.
tlr_landmark_base_cohort <- function(data, tlr_col = "tlr") {
  needed <- c(TLR_BASE_COHORT_VARS, tlr_col)
  missing <- setdiff(needed, names(data))
  if (length(missing) > 0) {
    stop("tlr_landmark_base_cohort(): missing column(s): ",
         paste(missing, collapse = ", "))
  }
  complete <- stats::complete.cases(data[, TLR_BASE_COHORT_VARS, drop = FALSE])
  data[complete & !is.na(data[[tlr_col]]), , drop = FALSE]
}

message("TLR landmark helpers loaded from tlr_landmark.R")
