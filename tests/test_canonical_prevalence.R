# Issue #166: identical treatment has zero increments under prevalence changes.
# Synthetic data only; no production cache is accessed.
source("R/model_configs.R")
source("R/prediction_functions.R")
source("R/model_fun.R")
source("R/calculate_outcomes.R")
source("R/parameter_distributions.R")
source("R/psa_functions.R")
#' Mock predict() method: exponential survival whose per-patient rate rises with
#' age, CRP and TMB/BRAF; returns flexsurv-style nested `.pred` survival frames.
predict.mock_survival_model <- function(object, newdata, type, times, ...) {
  v_rate <- (0.01 + newdata$Age / 3000 + as.numeric(as.character(newdata$crp)) / 100 +
            as.numeric(as.character(newdata$tmb_braf)) / 50) * object$rate
  structure(list(.pred = lapply(v_rate, function(r)
    data.frame(.pred_survival = exp(-r * times)))), class = "data.frame", row.names = seq_along(v_rate))
}
data <- data.frame(ID = 1:9, Age = seq(45, 85, 5), sex = factor(rep(0:1, length.out = 9)),
  Rx = factor(rep(c("control", "experimental"), length.out = 9)),
  crp = c(0,0,0,0,1,1,1,1,1), tmb_braf = c(0,0,0,1,0,1,1,1,1),
  OSwk = c(rep(100, 8), NA), Death = 1, PFSwk = 50, Progression = 1)
data_complete <- economic_prediction_population(data)
stopifnot(nrow(data_complete) == 8L)
models <- list(os = structure(list(rate = 1), class = "mock_survival_model"),
               pfs = structure(list(rate = 2), class = "mock_survival_model"))
strategies_df <- data.frame(id = get_strategies(), prevalence = c(1, 0.5, 0.5))
preds <- generate_population_averaged_predictions(models, strategies_df, data_complete, 0:52, quiet = TRUE)
params <- list(time_horizon = 52, cl = 1 / 52,  # grid owned by params (#172)
  dr_costs = 0.04, dr_effects = 0.04, u_np = 0.73, u_p = 0.59, u_decrement = 0.14,
  c_drug_nivo = 0, c_drug_FLOX = 100, c_test_CT = 10, c_test_blood = 5,
  c_test_CRP = 0, c_test_NGS = 0, c_test_biomarker = list(crp = 0, tmb_braf = 0),
  c_other_visit = 10, c_other_baseline = 100, c_other_follow = 20, c_other_last = 1000)
for (name in c("l_nivo", "l_FLOX_exp", "l_FLOX_control", "l_CT", "l_blood", "l_visit"))
  params[[name]] <- rep(1, 53)
params <- set_population_predictions(params, preds)
joint <- joint_biomarker_population(data_complete)
params$joint_counts <- joint$counts
for (cell in names(joint$probabilities)) params[[paste0("p_joint_", cell)]] <- unname(joint$probabilities[cell])
# Load only the production PSA helper, without running script 06 or its cache.
for (expr in parse("analysis/06_sampling.R")) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name("generate_psa_population_averaged_predictions"))) eval(expr)
}
n_samples <- 1L
sampling_models <- list(joint = list(samples = list(list(os = list(model = models$os),
  pfs = list(model = models$pfs), failed = FALSE))))
#' Assert that all strategies of a model_fun() result have identical cost and
#' effect (to 1e-10) and that no PSA fallback was used.
assert_equal_strategies <- function(result) {
  stopifnot(max(abs(result$Cost - result$Cost[1])) < 1e-10,
            max(abs(result$Effect - result$Effect[1])) < 1e-10,
            !isTRUE(attr(result, "fallback_used")))
}
for (crp in c(0.4, 0.5, 0.6)) for (tmb in c(0.4, 0.5, 0.6)) {
  p <- params; p$p_crp <- crp; p$p_tmb_braf <- tmb
  v_weights <- model_population_weights(p, data_complete)
  stopifnot(abs(sum(v_weights) - 1) < 1e-12,
    abs(sum(v_weights * data_complete$crp) - crp) < 1e-10,
    abs(sum(v_weights * data_complete$tmb_braf) - tmb) < 1e-10)
  df_deterministic <- model_fun(p, time_horizon = 52)
  df_probabilistic <- model_fun(p, time_horizon = 52, determpsa = "psa", sim_idx = 1)
  assert_equal_strategies(df_deterministic)
  assert_equal_strategies(df_probabilistic)
  stopifnot(isTRUE(all.equal(df_deterministic, df_probabilistic, tolerance = 1e-10)))
}
# Endpoint-missing patients cannot leak into the PSA through the full data.
df_before <- model_fun(params, time_horizon = 52, determpsa = "psa", sim_idx = 1)
data$Age[9] <- 200
df_after <- model_fun(params, time_horizon = 52, determpsa = "psa", sim_idx = 1)
stopifnot(identical(df_before, df_after))
p <- params; p$p_tmb_braf <- 0.6
df_changed <- model_fun(p, time_horizon = 52)
stopifnot(abs(df_changed$Effect[1] - df_before$Effect[1]) > 1e-4)
config <- configure_parameter_distributions(params)
df_draws <- generate_psa_samples(config$distributions, 20000, seed = 166L)
v_joint_keys <- paste0("p_joint_", c("00", "01", "10", "11"))
stopifnot(max(abs(rowSums(df_draws[v_joint_keys]) - 1)) < 1e-12,
  all(as.matrix(df_draws[v_joint_keys]) >= 0),
  identical(df_draws$p_crp, df_draws$p_joint_10 + df_draws$p_joint_11),
  identical(df_draws$p_tmb_braf, df_draws$p_joint_01 + df_draws$p_joint_11),
  max(abs(colMeans(df_draws[v_joint_keys]) - joint$probabilities)) < 0.005,
  abs(cov(df_draws$p_crp, df_draws$p_tmb_braf) - (0.375 - 0.25) / 9) < 0.001,
  identical(config$groups$prevalence, v_joint_keys[1:3]))
for (i in 1:10) {
  p <- params
  for (name in names(config$distributions)) p[[name]] <- df_draws[[name]][i]
  assert_equal_strategies(model_fun(p, time_horizon = 52))
  assert_equal_strategies(model_fun(p, time_horizon = 52, determpsa = "psa", sim_idx = 1))
}
cat("PASS: one population, zero identical-treatment increments under marginal/joint changes, no endpoint-missingness leak, coherent joint PSA draws.\n")
