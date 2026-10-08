# ===============================================================================
# TEST: multivariate-normal survival sampling and PSA draw/model alignment
# (issue #156)
# ===============================================================================
# 1. A flexsurvreg object rebuilt from a sampled coefficient vector predicts
#    exactly what the closed-form survival function gives for that vector, and
#    the accessor returns the legacy sample shape from the draw matrices.
# 2. get_joint_sampling_models() resolves the single joint component and the
#    legacy per-biomarker layout.
# 3. build_psa_obj() records, per retained row, the cached model that produced
#    it, and aligns additional (model-derived) columns by that index rather
#    than by the draw number.
# 4. The cache fingerprint changes when a formula, the data, a distribution,
#    n_samples, the seed or the method changes, and not otherwise.
#
# Uses simulated data only; no trial data or caches are needed.
# Run: "C:\Program Files\R\R-4.3.2\bin\x64\Rscript.exe" tests/test_mvn_sampling.R
# ===============================================================================

suppressPackageStartupMessages({
  library(survival)
  library(flexsurv)
  library(mvtnorm)
})
suppressMessages({
  source("R/model_configs.R")
  source("R/cache_paths.R")
  source("R/prediction_functions.R")
  source("R/parameter_distributions.R")
  source("R/psa_functions.R")
})

# ---------------------------------------------------------------------------
# 1. Rebuilt model == closed-form survival for the drawn coefficients
# ---------------------------------------------------------------------------
set.seed(42)
n <- 120
df_sim <- data.frame(
  Age = round(runif(n, 45, 80)),
  Rx = factor(sample(c("control", "experimental"), n, TRUE),
              levels = c("control", "experimental")),
  crp = factor(rbinom(n, 1, 0.4), levels = c(0, 1))
)
v_lp <- -4 + 0.01 * df_sim$Age + 0.3 * (df_sim$Rx == "experimental") - 0.4 * (df_sim$crp == "1")
df_sim$time <- rgamma(n, shape = 0.9, rate = exp(v_lp))
df_sim$event <- as.integer(df_sim$time < 250)
df_sim$time <- pmin(df_sim$time, 250)

fit <- flexsurvreg(Surv(time, event) ~ Age + Rx + crp * Rx, data = df_sim, dist = "gamma")

set.seed(7)
v_draw <- mvtnorm::rmvnorm(1, fit$opt$par, fit$cov)[1, ]
names(v_draw) <- rownames(fit$res)
sampled <- build_sampled_flexsurv_model(fit, v_draw)

df_newdata <- df_sim[1:3, ]
v_times <- c(26, 52, 104)
pred <- predict(sampled, newdata = df_newdata, type = "survival", times = v_times)
m_pred <- extract_all_survival_probabilities(pred)

X <- model.matrix(~ Age + Rx + crp * Rx, df_newdata)[, -1, drop = FALSE]
v_cov_names <- rownames(fit$res)[fit$covpars]
stopifnot(identical(colnames(X), v_cov_names))
m_manual <- sapply(seq_len(nrow(df_newdata)), function(i) {
  rate <- exp(v_draw["rate"] + sum(v_draw[v_cov_names] * X[i, ]))
  1 - pgamma(v_times, shape = exp(v_draw["shape"]), rate = rate)
})
stopifnot(isTRUE(all.equal(unname(m_pred), unname(m_manual), tolerance = 1e-10)))

# The fitted model itself is untouched and predicts differently from the draw.
m_orig <- extract_all_survival_probabilities(
  predict(fit, newdata = df_newdata, type = "survival", times = v_times)
)
stopifnot(!isTRUE(all.equal(unname(m_orig), unname(m_pred))))

# Accessor: legacy sample shape from the draw matrices, index validation.
component <- list(
  method = "mvn_v1",
  draws = list(os = rbind(draw = v_draw, v_draw + 0.01), pfs = rbind(v_draw - 0.02, draw = v_draw)),
  original_os = fit, original_pfs = fit,
  dist_os = "gamma", dist_pfs = "gamma", n_samples = 2L, seed = 123L
)
s2 <- sampled_survival_models(component, 2L)
stopifnot(
  setequal(names(s2), c("os", "pfs", "failed")),
  identical(s2$failed, FALSE),
  isTRUE(all.equal(unname(s2$os$coefficients), unname(v_draw + 0.01))),
  isTRUE(all.equal(unname(s2$pfs$coefficients), unname(v_draw))),
  inherits(s2$os$model, "flexsurvreg"),
  identical(s2$os$dist, "gamma"),
  inherits(tryCatch(sampled_survival_models(component, 3L), error = identity), "error"),
  inherits(tryCatch(sampled_survival_models(component, 0L), error = identity), "error")
)
legacy <- list(samples = list(list(os = list(coefficients = 1), failed = TRUE)))
stopifnot(isTRUE(sampled_survival_models(legacy, 1L)$failed))
cat("PASS: sampled flexsurvreg objects reproduce closed-form survival\n")

# ---------------------------------------------------------------------------
# 2. Joint component resolution
# ---------------------------------------------------------------------------
joint_layout <- list(joint = component, biomarkers = c("crp", "tmb_braf"))
legacy_layout <- setNames(list(component, component), get_biomarkers())
stopifnot(
  identical(get_joint_sampling_models(joint_layout), component),
  identical(get_joint_sampling_models(legacy_layout), component),
  inherits(tryCatch(get_joint_sampling_models(list(other = 1)), error = identity), "error")
)
cat("PASS: joint sampling component resolves for both cache layouts\n")

# ---------------------------------------------------------------------------
# 3. model_idx alignment in build_psa_obj()
# ---------------------------------------------------------------------------
# Draw 2 fails with its own cached model and succeeds with model 3; the
# model-derived column must then carry model 3's value in row 2.
#' Stub model: control cost 100 + sim_idx; errors when draw 2 is paired with its
#' own economic draw, so build_psa_obj() must re-run it with another model.
model_fun <- function(params, time_horizon, cl, determpsa, return_traces, sim_idx) {
  df_result <- data.frame(Strategy = "control", Cost = 100 + sim_idx, Effect = 1)
  attr(df_result, "fallback_used") <- FALSE
  if (identical(sim_idx, 2L) && identical(params$cost_draw, params$draws_for_2)) {
    stop("model 2 cannot be used with draw 2")
  }
  df_result
}
#' Identity stub for the legacy diagnostic-price sync; returns `params` unchanged.
sync_biomarker_test_costs <- function(params) params
distributions <- list(cost_draw = list(dist = "unif", min = 0, max = 1))
draws_for_2 <- generate_psa_samples(distributions, 4L, seed = 5L)$cost_draw[2]
base <- list(draws_for_2 = draws_for_2)
df_coef_by_model <- data.frame(b_model = c(10, 20, 30, 40))
v_built <- suppressWarnings(capture.output(
  psa_build <- build_psa_obj(
    l_params_base = base, param_distributions = distributions,
    strategies = "control", time_horizon = 1, cl = 1, n_sim = 4L, seed = 5L,
    additional_params = df_coef_by_model, fallback_threshold = 1
  )
))
stopifnot(
  identical(psa_build$psa_results$model_idx, c(1L, 3L, 3L, 4L)),
  identical(psa_build$psa_params$model_idx, c(1L, 3L, 3L, 4L)),
  identical(psa_build$psa_obj$model_idx, c(1L, 3L, 3L, 4L)),
  identical(psa_build$psa_params$b_model, c(10, 30, 30, 40)),
  identical(psa_build$psa_params$sim, 1:4),
  identical(as.numeric(psa_build$psa_obj$cost$control), c(101, 103, 103, 104))
)
cat("PASS: replaced PSA draws carry the index of the model that produced them\n")

# ---------------------------------------------------------------------------
# 4. Cache fingerprint sensitivity
# ---------------------------------------------------------------------------
base_args <- list(
  formulas = list(os = Surv(time, event) ~ Age + Rx, pfs = Surv(time, event) ~ Age + Rx),
  data = df_sim, distributions = list(os = "gamma", pfs = "gamma"),
  n_samples = 10L, seed = 123L, method = "mvn_v1"
)
#' Sampling-cache fingerprint string for a list of sampling_cache_fingerprint()
#' arguments.
fp <- function(args) do.call(sampling_cache_fingerprint, args)$fingerprint
reference <- fp(base_args)
#' Fingerprint of the base arguments with the named arguments in `...` replaced.
change <- function(...) {
  updates <- list(...)
  a <- base_args
  a[names(updates)] <- updates
  fp(a)
}
stopifnot(
  identical(fp(base_args), reference),
  change(n_samples = 11L) != reference,
  change(seed = 124L) != reference,
  change(method = "bootstrap_v1") != reference,
  change(distributions = list(os = "weibull", pfs = "gamma")) != reference,
  change(formulas = list(os = Surv(time, event) ~ Age, pfs = base_args$formulas$pfs)) != reference,
  change(data = df_sim[-1, ]) != reference,
  change(data = df_sim[, rev(names(df_sim))]) == reference   # column order is irrelevant
)
#' PSA-cache fingerprint string for the parameter list `params`, all other
#' inputs held fixed.
psa_fp <- function(params) psa_cache_fingerprint(
  reference, params, distributions, "control", 100L, 1L, 520, 1 / 52)$fingerprint
stopifnot(
  psa_fp(list(a = 1, b = 2)) == psa_fp(list(b = 2, a = 1)),
  psa_fp(list(a = 1, b = 2)) != psa_fp(list(a = 1, b = 3))
)
cat("PASS: cache fingerprints track formulas, data, distributions, n, seed and method\n")

cat("\nAll multivariate-normal sampling tests passed.\n")
