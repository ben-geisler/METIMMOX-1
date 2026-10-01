# Synthetic checks of the survival-validation helpers (issue #162).
# Requires no confidential data or production caches; run from the repo root:
#   Rscript tests/test_survival_validation.R
suppressPackageStartupMessages(library(survival))
source("R/survival_validation.R")

# ---- Benchmark file contract ---------------------------------------------------
b <- read_survival_benchmarks("data/external/survival_benchmarks.csv")
stopifnot(nrow(b) > 0, all(nzchar(b$doi)), !anyDuplicated(b$id),
          any(b$role == "comparable"), any(b$role == "lower_bound"))
bad <- b
bad$lower[1] <- bad$estimate[1] + 1
stopifnot(length(check_survival_benchmarks(bad)) > 0)
bad <- b
bad$role[1] <- "upper_bound"
stopifnot(length(check_survival_benchmarks(bad)) > 0)
bad <- b
bad$doi[2] <- ""
stopifnot(length(check_survival_benchmarks(bad)) > 0)
cat("PASS: benchmark file contract\n")

# ---- Tolerance -------------------------------------------------------------------
comp <- list(role = "comparable", estimate = 20, lower = 17, upper = 23)
low <- list(role = "lower_bound", estimate = 0.10, lower = NA, upper = NA)
stopifnot(
  evaluate_benchmark(17, comp) == "Within CI",      # inclusive bounds
  evaluate_benchmark(23, comp) == "Within CI",
  evaluate_benchmark(16.99, comp) == "Outside CI",
  evaluate_benchmark(23.01, comp) == "Outside CI",
  evaluate_benchmark(0.10, low) == "Above lower bound",
  evaluate_benchmark(0.099, low) == "Below lower bound",
  evaluate_benchmark(NA, comp) == "Not evaluable",
  evaluate_benchmark(NA, low) == "Not evaluable",
  identical(benchmark_passed(c("Within CI", "Outside CI", "Above lower bound",
                               "Below lower bound", "Not evaluable")),
            c(TRUE, FALSE, TRUE, FALSE, NA))
)
cat("PASS: pre-specified tolerance\n")

# ---- Kaplan-Meier landmarks ------------------------------------------------------
set.seed(162)
n <- 200
t_event <- rexp(n, rate = log(2) / 80)          # median 80 weeks
t_cens <- runif(n, 50, 200)                     # follow-up ends by week 200
time <- pmin(t_event, t_cens)
event <- as.numeric(t_event <= t_cens)
lm <- km_landmarks(time, event, years = c(1, 2, 3, 5))
ref <- summary(survfit(Surv(time, event) ~ 1, conf.type = "log-log"),
               times = c(52, 104))
stopifnot(
  max(abs(lm$surv[1:2] - ref$surv)) < 1e-12,
  max(abs(lm$lower[1:2] - ref$lower)) < 1e-12,
  max(abs(lm$upper[1:2] - ref$upper)) < 1e-12,
  identical(lm$n_risk[1:2], ref$n.risk),
  max(time) < 260,
  is.na(lm$surv[4]), lm$status[4] == "NE (beyond follow-up)"
)
# Few at risk is flagged, not suppressed
lm3 <- lm[lm$years == 3, ]
if (max(time) >= 156) {
  stopifnot(is.finite(lm3$surv),
            (lm3$n_risk < KM_FEW_AT_RISK) == grepl("Fewer", lm3$status))
}
# A curve that reaches zero stays at zero beyond the last time
lm0 <- km_landmarks(c(10, 20, 30), c(1, 1, 1), years = c(1, 2))
stopifnot(all(lm0$surv == 0), all(lm0$status == "All events before landmark"))
# No event before the landmark: estimate 1, no (degenerate) CI
lm1 <- km_landmarks(c(60, 70, 120, 150), c(0, 1, 1, 0), years = c(1, 2))
stopifnot(lm1$surv[1] == 1, is.na(lm1$lower[1]), is.na(lm1$upper[1]),
          lm1$status[1] == "No events before landmark",
          lm1$surv[2] < 1, is.finite(lm1$lower[2]))
cat("PASS: Kaplan-Meier landmarks\n")

med <- km_median_months(time, event)
km_tab <- summary(survfit(Surv(time, event) ~ 1, conf.type = "log-log"))$table
stopifnot(abs(med[["median"]] - km_tab[["median"]] * 12 / 52) < 1e-12,
          med[["lower"]] <= med[["median"]], med[["median"]] <= med[["upper"]])
cat("PASS: Kaplan-Meier median in months\n")

# ---- Model curve landmarks and median --------------------------------------------
lambda <- log(2) / 80
curve <- exp(-lambda * (0:520))
stopifnot(
  abs(curve_at_years(curve, 1) - exp(-lambda * 52)) < 1e-15,
  abs(curve_at_years(curve, 10) - exp(-lambda * 520)) < 1e-15,
  is.na(curve_at_years(curve, 11)),
  abs(curve_median_months(curve) - 80 * 12 / 52) < 1 * 12 / 52,
  is.na(curve_median_months(rep(1, 521))),
  curve_median_months(c(1, 0.4, 0.3)) < 1 * 12 / 52
)
cat("PASS: model curve landmarks and median\n")

# ---- Draw intervals ------------------------------------------------------------
rates <- log(2) / c(60, 70, 80, 90, 100)
mat <- t(sapply(rates, function(r) exp(-r * (0:520))))
iv <- draw_landmark_intervals(mat, years = c(1, 2))
stopifnot(nrow(iv) == 3, iv$measure[3] == "median_months",
          all(iv$lower <= iv$upper),
          abs(iv$lower[1] - quantile(exp(-rates * 52), 0.025, names = FALSE)) < 1e-12)
cat("PASS: draw percentile intervals\n")

# ---- Consistency with the model's 95% interval (added after review) ------------
stopifnot(
  benchmark_consistent(16, 24, comp) == "Consistent",          # 20 inside
  benchmark_consistent(21, 24, comp) == "Not consistent",      # 20 below interval
  benchmark_consistent(0.01, 0.10, low) == "Consistent",       # upper reaches bound
  benchmark_consistent(0.01, 0.09, low) == "Not consistent",
  benchmark_consistent(NA, 0.2, low) == "Not evaluable",
  identical(point_estimate_label(c("Within CI", "Below lower bound", "Above lower bound")),
            c("Within CI", "Below bound", "At or above bound"))
)
row_med <- list(measure = "median_months", time_years = NA)
row_s1 <- list(measure = "survival", time_years = 1)
stopifnot(
  abs(benchmark_model_value(curve, row_med) - curve_median_months(curve)) < 1e-12,
  abs(benchmark_model_value(curve, row_s1) - exp(-lambda * 52)) < 1e-15,
  abs(benchmark_model_interval(mat, row_s1)[2] -
        quantile(exp(-rates * 52), 0.975, names = FALSE)) < 1e-12
)
cat("PASS: consistency with the model interval\n")

cat("All survival-validation tests passed.\n")
