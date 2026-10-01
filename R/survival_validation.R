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
  df_benchmarks <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  v_problems <- check_survival_benchmarks(df_benchmarks)
  if (length(v_problems) > 0) {
    stop("Invalid survival benchmark file ", path, ":\n  ",
         paste(v_problems, collapse = "\n  "))
  }
  df_benchmarks
}

#' Contract of the benchmark table; returns a character vector of problems
check_survival_benchmarks <- function(b) {
  v_problems <- character(0)
  v_missing_cols <- setdiff(SURVIVAL_BENCHMARK_COLUMNS, names(b))
  if (length(v_missing_cols) > 0) {
    return(paste("missing columns:", paste(v_missing_cols, collapse = ", ")))
  }
  if (anyDuplicated(b$id)) v_problems <- c(v_problems, "duplicate id")
  if (any(!nzchar(b$doi) | is.na(b$doi))) v_problems <- c(v_problems, "row without DOI")
  if (!all(b$endpoint %in% c("OS", "PFS"))) v_problems <- c(v_problems, "endpoint not OS/PFS")
  if (!all(b$measure %in% c("median_months", "survival"))) {
    v_problems <- c(v_problems, "measure not median_months/survival")
  }
  if (!all(b$role %in% c("comparable", "lower_bound"))) {
    v_problems <- c(v_problems, "role not comparable/lower_bound")
  }
  if (!all(b$ci_method %in% c("reported", "none"))) {
    v_problems <- c(v_problems, "ci_method not reported/none")
  }
  if (any(!is.finite(b$estimate))) v_problems <- c(v_problems, "non-finite estimate")
  v_is_survival <- b$measure == "survival"
  if (any(v_is_survival & (!is.finite(b$time_years) | b$time_years <= 0))) {
    v_problems <- c(v_problems, "survival row without a positive time_years")
  }
  if (any(v_is_survival & (b$estimate < 0 | b$estimate > 1))) {
    v_problems <- c(v_problems, "survival estimate outside [0, 1]")
  }
  v_has_ci <- is.finite(b$lower) & is.finite(b$upper)
  if (any(v_has_ci != (b$ci_method == "reported"))) {
    v_problems <- c(v_problems, "ci_method disagrees with the presence of lower/upper")
  }
  if (any(v_has_ci & !(b$lower <= b$estimate & b$estimate <= b$upper))) {
    v_problems <- c(v_problems, "estimate outside its own interval")
  }
  if (any(b$role == "comparable" & !v_has_ci)) {
    v_problems <- c(v_problems, "comparable row without a reported CI")
  }
  v_problems
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
  v_weeks <- years * WEEKS_PER_YEAR
  s <- summary(fit, times = v_weeks, extend = TRUE)
  df_out <- data.frame(years = years, surv = s$surv, lower = s$lower,
                    upper = s$upper, n_risk = s$n.risk, status = "Estimated",
                    stringsAsFactors = FALSE)
  last_time <- max(time)
  last_surv <- utils::tail(fit$surv, 1)
  v_beyond <- v_weeks > last_time
  if (last_surv == 0) {
    df_out$status[v_beyond] <- "All events before landmark"
    df_out$surv[v_beyond] <- 0
    df_out$lower[v_beyond] <- NA
    df_out$upper[v_beyond] <- NA
  } else {
    df_out$status[v_beyond] <- "NE (beyond follow-up)"
    df_out[v_beyond, c("surv", "lower", "upper")] <- NA
  }
  # Before the first event the estimate is 1 and the CI degenerate (1 to 1)
  v_no_event <- !v_beyond & df_out$surv == 1
  df_out$status[v_no_event] <- "No events before landmark"
  df_out[v_no_event, c("lower", "upper")] <- NA
  v_few <- !v_beyond & !v_no_event & df_out$n_risk < KM_FEW_AT_RISK
  df_out$status[v_few] <- paste0("Fewer than ", KM_FEW_AT_RISK, " at risk")
  df_out
}

#' Kaplan-Meier median with 95% CI, in months
km_median_months <- function(time, event) {
  fit <- survival::survfit(survival::Surv(time, event) ~ 1, conf.type = "log-log")
  v_km_table <- summary(fit)$table
  c(median = unname(v_km_table["median"]), lower = unname(v_km_table["0.95LCL"]),
    upper = unname(v_km_table["0.95UCL"])) * MONTHS_PER_WEEK
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
  m_draws <- component$draws[[endpoint]]
  if (identical(fit$dlist$name, "gamma")) {
    rows <- lapply(idx, function(i) {
      direct_gamma_pop_avg(fit, m_draws[i, ], newdata, times)
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
  m_landmarks <- t(vapply(years, function(y) q(apply(mat, 1, curve_at_years, years = y)),
                   numeric(2)))
  v_median <- q(apply(mat, 1, curve_median_months))
  data.frame(years = c(years, NA), lower = c(m_landmarks[, 1], v_median[1]),
             upper = c(m_landmarks[, 2], v_median[2]),
             measure = c(rep("survival", length(years)), "median_months"))
}

# ===============================================================================
# BENCHMARK COMPARISON UNDER UNCERTAINTY (added after review, 1 October 2026)
# ===============================================================================
# The pre-specified rule above judges the point estimate of the parametric
# survival model. The headline assessment adopted after review asks whether
# the published value is compatible with the model's 95% coefficient-draw
# interval: for a comparable benchmark the published estimate must lie inside
# that interval; for a lower bound the interval's upper limit must reach it.

#' Parametric-model value matching a benchmark row (median in months, or
#' survival at `time_years`)
benchmark_model_value <- function(curve, row) {
  if (row$measure == "median_months") curve_median_months(curve)
  else curve_at_years(curve, row$time_years)
}

#' 95% percentile interval over coefficient-draw curves for a benchmark row
#'
#' @param mat Output of `draw_curves()` (one row per draw).
benchmark_model_interval <- function(mat, row) {
  v <- if (row$measure == "median_months") apply(mat, 1, curve_median_months)
       else apply(mat, 1, curve_at_years, years = row$time_years)
  stats::quantile(v, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
}

#' Is a published benchmark consistent with the model's 95% interval?
#'
#' @return "Consistent", "Not consistent" or "Not evaluable".
benchmark_consistent <- function(model_lower, model_upper, row) {
  if (!is.finite(model_lower) || !is.finite(model_upper)) return("Not evaluable")
  ok <- if (identical(row$role, "comparable")) {
    row$estimate >= model_lower && row$estimate <= model_upper
  } else if (identical(row$role, "lower_bound")) {
    model_upper >= row$estimate
  } else {
    stop("Unknown benchmark role: ", row$role)
  }
  if (ok) "Consistent" else "Not consistent"
}

#' Short display label for a pre-specified point-estimate verdict
point_estimate_label <- function(verdict) {
  v_labels <- c("Within CI" = "Within CI", "Outside CI" = "Outside CI",
              "Above lower bound" = "At or above bound",
              "Below lower bound" = "Below bound", "Not evaluable" = "--")
  unname(v_labels[verdict])
}
