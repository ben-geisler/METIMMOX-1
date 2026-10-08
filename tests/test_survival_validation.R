# Synthetic checks of the survival-validation helpers (issue #162).
# Requires no confidential data or production caches; run from the repo root:
#   Rscript tests/test_survival_validation.R
suppressPackageStartupMessages(library(survival))
source("R/survival_validation.R")

# ---- Benchmark file contract ---------------------------------------------------
# The benchmark file lives in data/external/, which the public repository does
# not carry (README, "Data"); its contract is checked wherever the file exists.
benchmark_file <- "data/external/survival_benchmarks.csv"
if (file.exists(benchmark_file)) {
  df_benchmarks <- read_survival_benchmarks(benchmark_file)
  stopifnot(nrow(df_benchmarks) > 0, all(nzchar(df_benchmarks$doi)), !anyDuplicated(df_benchmarks$id),
            any(df_benchmarks$role == "comparable"), any(df_benchmarks$role == "lower_bound"))
  df_bad <- df_benchmarks
  df_bad$lower[1] <- df_bad$estimate[1] + 1
  stopifnot(length(check_survival_benchmarks(df_bad)) > 0)
  df_bad <- df_benchmarks
  df_bad$role[1] <- "upper_bound"
  stopifnot(length(check_survival_benchmarks(df_bad)) > 0)
  df_bad <- df_benchmarks
  df_bad$doi[2] <- ""
  stopifnot(length(check_survival_benchmarks(df_bad)) > 0)
  cat("PASS: benchmark file contract\n")
} else {
  cat("SKIP: benchmark file contract (", benchmark_file, " absent)\n", sep = "")
}

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
v_t_event <- rexp(n, rate = log(2) / 80)          # median 80 weeks
v_t_cens <- runif(n, 50, 200)                     # follow-up ends by week 200
v_time <- pmin(v_t_event, v_t_cens)
v_event <- as.numeric(v_t_event <= v_t_cens)
df_lm <- km_landmarks(v_time, v_event, years = c(1, 2, 3, 5))
ref <- summary(survfit(Surv(v_time, v_event) ~ 1, conf.type = "log-log"),
               times = c(52, 104))
stopifnot(
  max(abs(df_lm$surv[1:2] - ref$surv)) < 1e-12,
  max(abs(df_lm$lower[1:2] - ref$lower)) < 1e-12,
  max(abs(df_lm$upper[1:2] - ref$upper)) < 1e-12,
  identical(df_lm$n_risk[1:2], ref$n.risk),
  max(v_time) < 260,
  is.na(df_lm$surv[4]), df_lm$status[4] == "NE (beyond follow-up)"
)
# Few at risk is flagged, not suppressed
df_lm3 <- df_lm[df_lm$years == 3, ]
if (max(v_time) >= 156) {
  stopifnot(is.finite(df_lm3$surv),
            (df_lm3$n_risk < KM_FEW_AT_RISK) == grepl("Fewer", df_lm3$status))
}
# A curve that reaches zero stays at zero beyond the last time
df_lm0 <- km_landmarks(c(10, 20, 30), c(1, 1, 1), years = c(1, 2))
stopifnot(all(df_lm0$surv == 0), all(df_lm0$status == "All events before landmark"))
# No event before the landmark: estimate 1, no (degenerate) CI
df_lm1 <- km_landmarks(c(60, 70, 120, 150), c(0, 1, 1, 0), years = c(1, 2))
stopifnot(df_lm1$surv[1] == 1, is.na(df_lm1$lower[1]), is.na(df_lm1$upper[1]),
          df_lm1$status[1] == "No events before landmark",
          df_lm1$surv[2] < 1, is.finite(df_lm1$lower[2]))
cat("PASS: Kaplan-Meier landmarks\n")

v_med <- km_median_months(v_time, v_event)
v_km_tab <- summary(survfit(Surv(v_time, v_event) ~ 1, conf.type = "log-log"))$table
stopifnot(abs(v_med[["median"]] - v_km_tab[["median"]] * 12 / 52) < 1e-12,
          v_med[["lower"]] <= v_med[["median"]], v_med[["median"]] <= v_med[["upper"]])
cat("PASS: Kaplan-Meier median in months\n")

# ---- Model curve landmarks and median --------------------------------------------
lambda <- log(2) / 80
v_curve <- exp(-lambda * (0:520))
stopifnot(
  abs(curve_at_years(v_curve, 1) - exp(-lambda * 52)) < 1e-15,
  abs(curve_at_years(v_curve, 10) - exp(-lambda * 520)) < 1e-15,
  is.na(curve_at_years(v_curve, 11)),
  abs(curve_median_months(v_curve) - 80 * 12 / 52) < 1 * 12 / 52,
  is.na(curve_median_months(rep(1, 521))),
  curve_median_months(c(1, 0.4, 0.3)) < 1 * 12 / 52
)
cat("PASS: model curve landmarks and median\n")

# ---- Draw intervals ------------------------------------------------------------
v_rates <- log(2) / c(60, 70, 80, 90, 100)
m_curves <- t(sapply(v_rates, function(r) exp(-r * (0:520))))
df_iv <- draw_landmark_intervals(m_curves, years = c(1, 2))
stopifnot(nrow(df_iv) == 3, df_iv$measure[3] == "median_months",
          all(df_iv$lower <= df_iv$upper),
          abs(df_iv$lower[1] - quantile(exp(-v_rates * 52), 0.025, names = FALSE)) < 1e-12)
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
  abs(benchmark_model_value(v_curve, row_med) - curve_median_months(v_curve)) < 1e-12,
  abs(benchmark_model_value(v_curve, row_s1) - exp(-lambda * 52)) < 1e-15,
  abs(benchmark_model_interval(m_curves, row_s1)[2] -
        quantile(exp(-v_rates * 52), 0.975, names = FALSE)) < 1e-12
)
cat("PASS: consistency with the model interval\n")

cat("All survival-validation tests passed.\n")
