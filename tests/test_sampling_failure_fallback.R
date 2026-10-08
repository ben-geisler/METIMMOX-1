# Regression test for the survival coefficient sampler's failure handling
# (issues #156, #159, #186).
#
# Since #159 sample_survival_coefficients() fits OS and PFS, estimates their
# joint covariance with a paired patient bootstrap
# (estimate_joint_survival_covariance(), which calls flexsurv::flexsurvreg()
# directly) and draws joint multivariate-normal coefficients. Individual draws
# cannot fail; the failure modes are a fit without a usable covariance matrix and
# too many failed bootstrap pairs, and both must stop before any draw is made.
# The test fits real flexsurvreg() models on small simulated data; no trial data.

suppressMessages({
  library(survival)
  library(flexsurv)
})
source("R/joint_survival_sampling.R")

# Evaluate only the sampler so this focused test does not run the
# cache-generation side effects in 06_sampling.R.
sampling_expressions <- parse("analysis/06_sampling.R")
#' The single top-level `name <- ...` expression of analysis/06_sampling.R, so a
#' production helper can be evaluated without running the script.
extract_assignment <- function(name) {
  v_matches <- vapply(sampling_expressions, function(expr) {
    is.call(expr) &&
      identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name(name))
  }, logical(1))
  stopifnot(sum(v_matches) == 1L)
  sampling_expressions[[which(v_matches)]]
}

test_env <- new.env(parent = globalenv())
eval(extract_assignment("SAMPLING_METHOD"), envir = test_env)
eval(extract_assignment("sample_survival_coefficients"), envir = test_env)

# Simulated cohort with dependent endpoints: PFS is the earlier of progression
# and death, and both share a patient frailty, so OS and PFS coefficients
# estimated on the same resample covary.
set.seed(186)
n_patients <- 200L
df_sim <- data.frame(age = rnorm(n_patients))
v_frailty <- rnorm(n_patients)  # unobserved prognosis shared by both endpoints
v_death <- rweibull(n_patients, shape = 1.3, scale = exp(1.2 - 0.3 * df_sim$age + 0.8 * v_frailty))
v_progression <- rweibull(n_patients, shape = 1.3, scale = exp(0.4 - 0.3 * df_sim$age + 0.8 * v_frailty))
v_censor <- runif(n_patients, 2, 6)
df_sim$OS <- pmin(v_death, v_censor)
df_sim$Death <- as.integer(v_death <= v_censor)
df_sim$PFS <- pmin(v_death, v_progression, v_censor)
df_sim$Progression <- as.integer(pmin(v_death, v_progression) <= v_censor)
f_os <- Surv(OS, Death) ~ age
f_pfs <- Surv(PFS, Progression) ~ age

#' Run the sampler on the simulated cohort with a small bootstrap.
run_sampler <- function(n_samples = 2000L, n_bootstrap = 60L, data = df_sim,
                        formula_os = f_os, formula_pfs = f_pfs) {
  invisible(capture.output(result <- test_env$sample_survival_coefficients(
    formula_os = formula_os, formula_pfs = formula_pfs, data = data,
    dist_os = "weibull", dist_pfs = "weibull",
    n_samples = n_samples, seed = 11L, n_bootstrap = n_bootstrap
  )))
  result
}

# 1. Converged fits: every draw is finite, none fails, and the OS and PFS draws
#    are jointly (not independently) distributed.
result <- run_sampler()
v_names <- c("shape", "scale", "age")
v_cor <- vapply(v_names, function(p) cor(result$draws$os[, p], result$draws$pfs[, p]),
                numeric(1))
stopifnot(
  identical(result$method, test_env$SAMPLING_METHOD),
  result$n_failed == 0L,
  is.null(result$samples),
  identical(dim(result$draws$os), c(2000L, 3L)),
  identical(dim(result$draws$pfs), c(2000L, 3L)),
  identical(colnames(result$draws$os), v_names),
  identical(colnames(result$draws$pfs), v_names),
  all(is.finite(result$draws$os)), all(is.finite(result$draws$pfs)),
  result$joint_covariance$n_attempted == 60L,
  result$joint_covariance$n_failed == 0L,
  # Paired resamples carry the dependence between the endpoints.
  all(v_cor > 0.3),
  nrow(result$coefficient_summary) == 6L
)
cat("PASS: converged fits give finite, jointly correlated OS/PFS draws with n_failed = 0",
    sprintf("(draw correlations %s).\n", paste(round(v_cor, 2), collapse = ", ")))

# 2. A fit without a covariance matrix must stop before any draw is made.
draw_calls <- 0L
test_env$draw_joint_survival_coefficients <- function(...) {
  draw_calls <<- draw_calls + 1L
  draw_joint_survival_coefficients(...)
}
test_env$flexsurvreg <- function(formula, ...) {
  fit <- flexsurv::flexsurvreg(formula, ...)
  if (identical(formula, f_pfs)) fit$cov <- NULL
  fit
}
failure <- tryCatch(run_sampler(n_samples = 10L), error = identity)
stopifnot(
  inherits(failure, "error"),
  grepl("covariance is unavailable", conditionMessage(failure)),
  draw_calls == 0L
)
rm("flexsurvreg", envir = test_env)
cat("PASS: a fit without a covariance matrix stops before any draw.\n")

# 3. More than 10% failed bootstrap pairs must stop generation. A binary
#    covariate carried by a single patient is absent from about a third of the
#    resamples, so those refits fail the coefficient checks.
df_rare <- df_sim
df_rare$marker <- 0
df_rare$marker[which(df_rare$Death == 1L)[1]] <- 1
failure <- tryCatch(
  run_sampler(n_samples = 10L, data = df_rare,
              formula_os = Surv(OS, Death) ~ age + marker,
              formula_pfs = Surv(PFS, Progression) ~ age + marker),
  error = identity
)
stopifnot(
  inherits(failure, "error"),
  grepl("Paired covariance bootstrap failed", conditionMessage(failure)),
  grepl("limit 10%", conditionMessage(failure), fixed = TRUE),
  draw_calls == 0L
)
cat("PASS: more than 10% failed bootstrap pairs stops before any draw.\n")

cat("Sampling failure handling regression test passed.\n")
