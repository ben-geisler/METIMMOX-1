# Regression tests for deterministic stochastic blocks and cache seed metadata.
# Run from the repository root with Rscript. No trial data or caches are needed.

source("R/parameter_distributions.R")  # apply_derived_psa_parameters()
source("R/psa_functions.R")
source("R/joint_survival_sampling.R")

param_distributions <- list(
  lognormal = list(dist = "lnorm", meanlog = 0, sdlog = 0.2),
  probability = list(dist = "beta", shape1 = 2, shape2 = 3),
  cost = list(dist = "gamma", shape = 4, rate = 2),
  prevalence = list(dist = "unif", min = 0.1, max = 0.9)
)

# PSA draws must be independent of the incoming RNG state, including the
# difference between cache-hit and cache-miss execution paths.
df_psa_reference <- generate_psa_samples(param_distributions, 25, seed = 123L)
set.seed(987L)
for (i in seq_len(20)) {
  sample.int(100, 100, replace = TRUE)
}
df_psa_after_rng_use <- generate_psa_samples(param_distributions, 25, seed = 123L)
df_psa_other_seed <- generate_psa_samples(param_distributions, 25, seed = 124L)

stopifnot(
  identical(df_psa_reference, df_psa_after_rng_use),
  !identical(unname(as.matrix(df_psa_reference)),
             unname(as.matrix(df_psa_other_seed))),
  identical(attr(df_psa_reference, "seed"), 123L),
  psa_samples_seed_matches(df_psa_reference, 123L),
  !psa_samples_seed_matches(df_psa_reference, 124L),
  !psa_samples_seed_matches(data.frame(sim = seq_len(25)), 123L)
)

# Evaluate only the sampling function and cache helper so the analysis script's
# cache-generation code is not executed.
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

sampling_env <- new.env(parent = globalenv())
eval(extract_assignment("SAMPLING_METHOD"), envir = sampling_env)
eval(extract_assignment("sample_survival_coefficients"), envir = sampling_env)
eval(extract_assignment("sampling_cache_seed_matches"), envir = sampling_env)

# A minimal stand-in for a converged flexsurvreg fit: two baseline parameters
# and one covariate effect, all optimised, with a proper covariance matrix.
#' Stand-in for flexsurvreg(): a converged three-parameter fit (shape, rate,
#' predictor) with a diagonal covariance matrix; ignores its arguments.
mock_fit <- function(...) {
  v_par_names <- c("shape", "rate", "predictor")
  v_est <- c(0.1, -3, 0.2)
  m_res <- matrix(v_est, ncol = 4, nrow = 3,
                dimnames = list(v_par_names, c("est", "L95%", "U95%", "se")))
  list(
    res = m_res, res.t = m_res, coefficients = setNames(v_est, v_par_names),
    opt = list(par = v_est), cov = diag(c(0.01, 0.04, 0.02)),
    optpars = 1:3, fixedpars = integer(0),
    dlist = list(pars = c("shape", "rate"),
                 inv.transforms = list(exp, exp))
  )
}
sampling_env$flexsurvreg <- mock_fit
# Isolate the normal draw RNG from covariance estimation in this unit test;
# the paired-refit estimator is exercised in test_joint_survival_sampling.R.
sampling_env$estimate_joint_survival_covariance <- function(...) {
  list(covariance = diag(rep(c(0.01, 0.04, 0.02), 2)))
}

sampling_args <- list(
  formula_os = response ~ predictor,
  formula_pfs = response ~ predictor,
  data = data.frame(predictor = seq_len(20)),
  n_samples = 5
)

sampling_reference <- do.call(
  sampling_env$sample_survival_coefficients,
  c(sampling_args, list(seed = 123L))
)
invisible(runif(100))
sampling_after_rng_use <- do.call(
  sampling_env$sample_survival_coefficients,
  c(sampling_args, list(seed = 123L))
)
sampling_other_seed <- do.call(
  sampling_env$sample_survival_coefficients,
  c(sampling_args, list(seed = 124L))
)

stopifnot(
  identical(sampling_reference$draws, sampling_after_rng_use$draws),
  !identical(sampling_reference$draws, sampling_other_seed$draws),
  identical(dim(sampling_reference$draws$os), c(5L, 3L)),
  identical(sampling_reference$seed, 123L),
  identical(sampling_reference$method, "mvn_joint_v2"),
  identical(sampling_reference$n_failed, 0L)
)

# The seed check resolves the joint component (issue #156 layout) and falls
# back to a legacy per-biomarker component through get_joint_sampling_models().
sampling_env$get_joint_sampling_models <- function(sampling_models) {
  if (!is.null(sampling_models$joint)) return(sampling_models$joint)
  sampling_models[["crp"]]
}
stopifnot(
  sampling_env$sampling_cache_seed_matches(list(joint = list(seed = 123L)), 123L),
  sampling_env$sampling_cache_seed_matches(list(crp = list(seed = 123L)), 123L),
  !sampling_env$sampling_cache_seed_matches(list(joint = list()), 123L),
  !sampling_env$sampling_cache_seed_matches(list(joint = list(seed = 124L)), 123L)
)

# The regression EVPPI point estimate is deterministic; its Monte Carlo
# standard error comes from random draws of the regression coefficients and
# must therefore follow the shared seed (issue #152).
evppi_env <- new.env(parent = globalenv())
sys.source("R/evppi_functions.R", envir = evppi_env)

n_evppi <- 400L
set.seed(1L)
v_evppi_theta <- rnorm(n_evppi)
psa_obj <- list(
  cost = cbind(strategy_one = rep(0, n_evppi),
               strategy_two = rep(0, n_evppi)),
  effect = cbind(strategy_one = rep(0, n_evppi),
                 strategy_two = v_evppi_theta + rnorm(n_evppi, sd = 0.5))
)
df_psa_params <- data.frame(
  parameter_one = v_evppi_theta,
  parameter_two = rnorm(n_evppi)
)

#' Two-parameter group EVPPI (B = 50) under `seed` after unrelated RNG use;
#' returns c(estimate, se).
run_grouped_evppi <- function(seed) {
  runif(50)
  result <- NULL
  capture.output(result <- evppi_env$calculate_evppi_regression(
    psa_obj = psa_obj,
    psa_params = df_psa_params,
    param_names = names(df_psa_params),
    wtp = 1,
    B = 50,
    seed = seed
  ))
  c(estimate = result$evppi, se = result$evppi_se)
}

v_evppi_reference <- run_grouped_evppi(123L)
v_evppi_repeat <- run_grouped_evppi(123L)
v_evppi_other_seed <- run_grouped_evppi(124L)

stopifnot(
  all(is.finite(v_evppi_reference)),
  identical(v_evppi_reference, v_evppi_repeat),
  identical(v_evppi_reference[["estimate"]], v_evppi_other_seed[["estimate"]]),
  !identical(v_evppi_reference[["se"]], v_evppi_other_seed[["se"]])
)

# Public wrappers must retain the documented default and forward it to their
# stochastic dependencies.
interface_env <- new.env(parent = globalenv())
sys.source("R/evppi_functions.R", envir = interface_env)
sys.source("R/scenario_analysis.R", envir = interface_env)

v_seeded_interfaces <- c(
  "calculate_evppi_regression",
  "run_evppi_analysis",
  "run_scenario_psa",
  "run_scenario_evppi",
  "run_scenario_analysis",
  "run_all_scenarios"
)
stopifnot(all(vapply(v_seeded_interfaces, function(name) {
  identical(formals(interface_env[[name]])$seed, 123L)
}, logical(1))))

v_forwarding_wrappers <- setdiff(v_seeded_interfaces, "calculate_evppi_regression")
stopifnot(all(vapply(v_forwarding_wrappers, function(name) {
  grepl("seed = seed", paste(deparse(body(interface_env[[name]])), collapse = " "),
        fixed = TRUE)
}, logical(1))))

# Report-level stochastic blocks must use the shared analysis seed. The
# permutation helper is shared by the DAG report and its Table S1 vignette, so
# seeding inside the helper also makes their independently rendered p-values
# identical.
assoc_env <- new.env(parent = globalenv())
sys.source("R/assoc_tests.R", envir = assoc_env)

stopifnot(
  identical(formals(assoc_env$run_permutation_test)$seed, 123L),
  grepl(
    "set.seed(seed)",
    paste(deparse(body(assoc_env$run_permutation_test)), collapse = " "),
    fixed = TRUE
  )
)

dag_report <- paste(
  readLines("reports/dag_associations.qmd", warn = FALSE),
  collapse = "\n"
)
dag_table <- paste(
  readLines("reports/vignettes/clin_effect_table_s3.qmd", warn = FALSE),
  collapse = "\n"
)
figure_one <- paste(
  readLines("reports/vignettes/figure1.qmd", warn = FALSE),
  collapse = "\n"
)

stopifnot(
  grepl("seed = analysis_seed", dag_report, fixed = TRUE),
  grepl("seed = analysis_seed", dag_table, fixed = TRUE),
  grepl(
    "set.seed(analysis_seed)\n  boot_indices <- sample",
    figure_one,
    fixed = TRUE
  )
)

cat("RNG reproducibility and cache seed validation tests passed.\n")
