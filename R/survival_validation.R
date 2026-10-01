# ===============================================================================
# SURVIVAL EXTRAPOLATION VALIDATION (issue #162)
# ===============================================================================
# External benchmarks and the pre-specified tolerance for comparing the modelled
# standard-of-care curves with published first-line mCRC survival. The benchmark
# file and the rule below were committed before any comparison was computed;
# see data/external/README.md for provenance.
# ===============================================================================

SURVIVAL_BENCHMARK_PATH <- "data/external/survival_benchmarks.csv"

SURVIVAL_BENCHMARK_COLUMNS <- c(
  "id", "source", "doi", "population", "era", "msi_selection", "median_age",
  "endpoint", "measure", "time_years", "estimate", "lower", "upper", "n",
  "events", "ci_method", "role", "note"
)

#' Read and check the published survival benchmarks
#'
#' @param path Path to the benchmark CSV (relative to the repository root).
#' @return Data frame with one row per benchmark.
read_survival_benchmarks <- function(path = here::here(SURVIVAL_BENCHMARK_PATH)) {
  b <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  problems <- check_survival_benchmarks(b)
  if (length(problems) > 0) {
    stop("Invalid survival benchmark file ", path, ":\n  ",
         paste(problems, collapse = "\n  "))
  }
  b
}

#' Contract of the benchmark table; returns a character vector of problems
check_survival_benchmarks <- function(b) {
  problems <- character(0)
  missing_cols <- setdiff(SURVIVAL_BENCHMARK_COLUMNS, names(b))
  if (length(missing_cols) > 0) {
    return(paste("missing columns:", paste(missing_cols, collapse = ", ")))
  }
  if (anyDuplicated(b$id)) problems <- c(problems, "duplicate id")
  if (any(!nzchar(b$doi) | is.na(b$doi))) problems <- c(problems, "row without DOI")
  if (!all(b$endpoint %in% c("OS", "PFS"))) problems <- c(problems, "endpoint not OS/PFS")
  if (!all(b$measure %in% c("median_months", "survival"))) {
    problems <- c(problems, "measure not median_months/survival")
  }
  if (!all(b$role %in% c("comparable", "lower_bound"))) {
    problems <- c(problems, "role not comparable/lower_bound")
  }
  if (!all(b$ci_method %in% c("reported", "none"))) {
    problems <- c(problems, "ci_method not reported/none")
  }
  if (any(!is.finite(b$estimate))) problems <- c(problems, "non-finite estimate")
  surv <- b$measure == "survival"
  if (any(surv & (!is.finite(b$time_years) | b$time_years <= 0))) {
    problems <- c(problems, "survival row without a positive time_years")
  }
  if (any(surv & (b$estimate < 0 | b$estimate > 1))) {
    problems <- c(problems, "survival estimate outside [0, 1]")
  }
  has_ci <- is.finite(b$lower) & is.finite(b$upper)
  if (any(has_ci != (b$ci_method == "reported"))) {
    problems <- c(problems, "ci_method disagrees with the presence of lower/upper")
  }
  if (any(has_ci & !(b$lower <= b$estimate & b$estimate <= b$upper))) {
    problems <- c(problems, "estimate outside its own interval")
  }
  if (any(b$role == "comparable" & !has_ci)) {
    problems <- c(problems, "comparable row without a reported CI")
  }
  problems
}

#' The pre-specified tolerance, as printed in the report
survival_benchmark_tolerance <- function() {
  paste(
    "For a comparable benchmark the model value must lie within the benchmark's",
    "reported 95% confidence interval (inclusive); for a lower-bound benchmark the",
    "model value must be at least the benchmark estimate. Medians are compared in",
    "months and survival proportions at the stated time. A value the model cannot",
    "provide is not evaluable. A failed check is reported, not enforced."
  )
}

#' Apply the pre-specified tolerance to one benchmark
#'
#' @param model_value Model median (months) or survival proportion matching the row.
#' @param row One row of the benchmark table (list or one-row data frame).
#' @return One of "Within CI", "Outside CI", "Above lower bound",
#'   "Below lower bound", "Not evaluable".
evaluate_benchmark <- function(model_value, row) {
  if (length(model_value) != 1 || !is.finite(model_value)) return("Not evaluable")
  if (identical(row$role, "comparable")) {
    if (!is.finite(row$lower) || !is.finite(row$upper)) return("Not evaluable")
    if (model_value >= row$lower && model_value <= row$upper) "Within CI" else "Outside CI"
  } else if (identical(row$role, "lower_bound")) {
    if (model_value >= row$estimate) "Above lower bound" else "Below lower bound"
  } else {
    stop("Unknown benchmark role: ", row$role)
  }
}

#' TRUE for verdicts that meet the tolerance, FALSE for failures, NA otherwise
benchmark_passed <- function(verdict) {
  ifelse(verdict %in% c("Within CI", "Above lower bound"), TRUE,
         ifelse(verdict %in% c("Outside CI", "Below lower bound"), FALSE, NA))
}
