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

# ===============================================================================
# KAPLAN-MEIER AND MODEL LANDMARKS
# ===============================================================================

WEEKS_PER_YEAR <- 52
MONTHS_PER_WEEK <- 12 / 52
KM_FEW_AT_RISK <- 5

#' Kaplan-Meier survival at landmark years, with log-log 95% CIs
#'
#' A landmark after the last observed time is not estimable ("NE") unless the
#' curve has already reached zero, in which case it stays zero.
#'
#' @param time Event/censoring times in weeks.
#' @param event 1 = event, 0 = censored.
#' @param years Landmarks in years.
#' @return Data frame: years, surv, lower, upper, n_risk, status.
km_landmarks <- function(time, event, years = c(1, 2, 3, 5)) {
  fit <- survival::survfit(survival::Surv(time, event) ~ 1, conf.type = "log-log")
  weeks <- years * WEEKS_PER_YEAR
  s <- summary(fit, times = weeks, extend = TRUE)
  out <- data.frame(years = years, surv = s$surv, lower = s$lower,
                    upper = s$upper, n_risk = s$n.risk, status = "Estimated",
                    stringsAsFactors = FALSE)
  last_time <- max(time)
  last_surv <- utils::tail(fit$surv, 1)
  beyond <- weeks > last_time
  if (last_surv == 0) {
    out$status[beyond] <- "All events before landmark"
    out$surv[beyond] <- 0
    out$lower[beyond] <- NA
    out$upper[beyond] <- NA
  } else {
    out$status[beyond] <- "NE (beyond follow-up)"
    out[beyond, c("surv", "lower", "upper")] <- NA
  }
  # Before the first event the estimate is 1 and the CI degenerate (1 to 1)
  no_event <- !beyond & out$surv == 1
  out$status[no_event] <- "No events before landmark"
  out[no_event, c("lower", "upper")] <- NA
  few <- !beyond & !no_event & out$n_risk < KM_FEW_AT_RISK
  out$status[few] <- paste0("Fewer than ", KM_FEW_AT_RISK, " at risk")
  out
}

#' Kaplan-Meier median with 95% CI, in months
km_median_months <- function(time, event) {
  fit <- survival::survfit(survival::Surv(time, event) ~ 1, conf.type = "log-log")
  tab <- summary(fit)$table
  c(median = unname(tab["median"]), lower = unname(tab["0.95LCL"]),
    upper = unname(tab["0.95UCL"])) * MONTHS_PER_WEEK
}

#' Value of a weekly curve (index 1 = week 0) at landmark years; NA beyond it
curve_at_years <- function(curve, years) {
  idx <- round(years * WEEKS_PER_YEAR) + 1
  ifelse(idx <= length(curve), curve[pmin(idx, length(curve))], NA_real_)
}

#' Median of a weekly survival curve in months (linear interpolation between
#' weeks); NA when the curve does not reach 0.5 within its horizon
curve_median_months <- function(curve) {
  k <- which(curve <= 0.5)[1]
  if (is.na(k)) return(NA_real_)
  if (k == 1) return(0)
  w <- (k - 2) + (curve[k - 1] - 0.5) / (curve[k - 1] - curve[k])
  w * MONTHS_PER_WEEK
}

#' Population-averaged curves of one endpoint for coefficient draws
#'
#' Uses the closed-form gamma evaluation of `direct_gamma_pop_avg()`
#' (R/pfs_os_violation_diagnostics.R) when the fit is gamma, otherwise rebuilds
#' each draw with `sampled_survival_models()` and predicts with
#' `predict_pop_avg()`. Each patient keeps the `Rx` in `newdata`.
#'
#' @param component Joint sampling component (`get_joint_sampling_models()`).
#' @param endpoint "os" or "pfs".
#' @param newdata Patients to average over.
#' @param idx Draw indices.
#' @param times Weekly grid.
#' @return Matrix, one row per draw and one column per time.
draw_curves <- function(component, endpoint, newdata, idx, times) {
  fit <- component[[paste0("original_", endpoint)]]
  draws <- component$draws[[endpoint]]
  if (identical(fit$dlist$name, "gamma")) {
    rows <- lapply(idx, function(i) {
      direct_gamma_pop_avg(fit, draws[i, ], newdata, times)
    })
  } else {
    rows <- lapply(idx, function(i) {
      predict_pop_avg(sampled_survival_models(component, i)[[endpoint]]$model,
                      newdata, times)
    })
  }
  do.call(rbind, rows)
}

#' 95% percentile interval of landmark survival and medians over draws
#'
#' @param mat Output of `draw_curves()`.
#' @param years Landmark years.
#' @return Data frame with one row per landmark plus one "median" row
#'   (`years = NA`, values in months).
draw_landmark_intervals <- function(mat, years) {
  q <- function(x) stats::quantile(x, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  land <- t(vapply(years, function(y) q(apply(mat, 1, curve_at_years, years = y)),
                   numeric(2)))
  med <- q(apply(mat, 1, curve_median_months))
  data.frame(years = c(years, NA), lower = c(land[, 1], med[1]),
             upper = c(land[, 2], med[2]),
             measure = c(rep("survival", length(years)), "median_months"))
}
