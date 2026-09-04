# Regression tests for deterministic stochastic blocks and cache seed metadata.
# Run from the repository root with Rscript. No trial data or caches are needed.

source("scripts/R/functions/psa_functions.R")

param_distributions <- list(
  lognormal = list(dist = "lnorm", meanlog = 0, sdlog = 0.2),
  probability = list(dist = "beta", shape1 = 2, shape2 = 3),
  cost = list(dist = "gamma", shape = 4, rate = 2),
  prevalence = list(dist = "unif", min = 0.1, max = 0.9)
)

# PSA draws must be independent of the incoming RNG state, including the
# difference between cache-hit and cache-miss execution paths.
psa_reference <- generate_psa_samples(param_distributions, 25, seed = 123L)
set.seed(987L)
for (i in seq_len(20)) {
  sample.int(100, 100, replace = TRUE)
}
psa_after_rng_use <- generate_psa_samples(param_distributions, 25, seed = 123L)
psa_other_seed <- generate_psa_samples(param_distributions, 25, seed = 124L)

stopifnot(
  identical(psa_reference, psa_after_rng_use),
  !identical(unname(as.matrix(psa_reference)),
             unname(as.matrix(psa_other_seed))),
  identical(attr(psa_reference, "seed"), 123L),
  psa_samples_seed_matches(psa_reference, 123L),
  !psa_samples_seed_matches(psa_reference, 124L),
  !psa_samples_seed_matches(data.frame(sim = seq_len(25)), 123L)
)

# Evaluate only the sampling function and cache helper so the analysis script's
# cache-generation code is not executed.
sampling_expressions <- parse("scripts/R/analysis/06_sampling.R")
extract_assignment <- function(name) {
  matches <- vapply(sampling_expressions, function(expr) {
    is.call(expr) &&
      identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name(name))
  }, logical(1))
  stopifnot(sum(matches) == 1L)
  sampling_expressions[[which(matches)]]
}

sampling_env <- new.env(parent = baseenv())
eval(extract_assignment("sample_correlated_survival"), envir = sampling_env)
eval(extract_assignment("sampling_cache_seed_matches"), envir = sampling_env)

sampling_env$flexsurvreg <- function(...) {
  list(coefficients = c(intercept = 0))
}

sampling_args <- list(
  formula_os = response ~ predictor,
  formula_pfs = response ~ predictor,
  data = data.frame(predictor = seq_len(20)),
  n_samples = 5
)

sampling_reference <- do.call(
  sampling_env$sample_correlated_survival,
  c(sampling_args, list(seed = 123L))
)
invisible(runif(100))
sampling_after_rng_use <- do.call(
  sampling_env$sample_correlated_survival,
  c(sampling_args, list(seed = 123L))
)
sampling_other_seed <- do.call(
  sampling_env$sample_correlated_survival,
  c(sampling_args, list(seed = 124L))
)

resample_indices <- function(result) {
  lapply(result$samples, `[[`, "resample_idx")
}

stopifnot(
  identical(resample_indices(sampling_reference),
            resample_indices(sampling_after_rng_use)),
  !identical(resample_indices(sampling_reference),
             resample_indices(sampling_other_seed)),
  identical(sampling_reference$seed, 123L)
)

valid_sampling_cache <- list(
  control = list(seed = 123L),
  crp = list(seed = 123L)
)
legacy_sampling_cache <- list(
  control = list(),
  crp = list()
)
mismatched_sampling_cache <- list(
  control = list(seed = 123L),
  crp = list(seed = 124L)
)

stopifnot(
  sampling_env$sampling_cache_seed_matches(
    valid_sampling_cache, c("control", "crp"), 123L
  ),
  !sampling_env$sampling_cache_seed_matches(
    legacy_sampling_cache, c("control", "crp"), 123L
  ),
  !sampling_env$sampling_cache_seed_matches(
    mismatched_sampling_cache, c("control", "crp"), 123L
  )
)

# The regression EVPPI point estimate is deterministic; its Monte Carlo
# standard error comes from random draws of the regression coefficients and
# must therefore follow the shared seed (issue #152).
evppi_env <- new.env(parent = globalenv())
sys.source("scripts/R/functions/evppi_functions.R", envir = evppi_env)

n_evppi <- 400L
set.seed(1L)
evppi_theta <- rnorm(n_evppi)
psa_obj <- list(
  cost = cbind(strategy_one = rep(0, n_evppi),
               strategy_two = rep(0, n_evppi)),
  effect = cbind(strategy_one = rep(0, n_evppi),
                 strategy_two = evppi_theta + rnorm(n_evppi, sd = 0.5))
)
psa_params <- data.frame(
  parameter_one = evppi_theta,
  parameter_two = rnorm(n_evppi)
)

run_grouped_evppi <- function(seed) {
  runif(50)
  result <- NULL
  capture.output(result <- evppi_env$calculate_evppi_regression(
    psa_obj = psa_obj,
    psa_params = psa_params,
    param_names = names(psa_params),
    wtp = 1,
    B = 50,
    seed = seed
  ))
  c(estimate = result$evppi, se = result$evppi_se)
}

evppi_reference <- run_grouped_evppi(123L)
evppi_repeat <- run_grouped_evppi(123L)
evppi_other_seed <- run_grouped_evppi(124L)

stopifnot(
  all(is.finite(evppi_reference)),
  identical(evppi_reference, evppi_repeat),
  identical(evppi_reference[["estimate"]], evppi_other_seed[["estimate"]]),
  !identical(evppi_reference[["se"]], evppi_other_seed[["se"]])
)

# Public wrappers must retain the documented default and forward it to their
# stochastic dependencies.
interface_env <- new.env(parent = globalenv())
sys.source("scripts/R/functions/evppi_functions.R", envir = interface_env)
sys.source("scripts/R/functions/scenario_analysis.R", envir = interface_env)

seeded_interfaces <- c(
  "calculate_evppi_regression",
  "run_evppi_analysis",
  "run_scenario_psa",
  "run_scenario_evppi",
  "run_scenario_analysis",
  "run_all_scenarios"
)
stopifnot(all(vapply(seeded_interfaces, function(name) {
  identical(formals(interface_env[[name]])$seed, 123L)
}, logical(1))))

forwarding_wrappers <- setdiff(seeded_interfaces, "calculate_evppi_regression")
stopifnot(all(vapply(forwarding_wrappers, function(name) {
  grepl("seed = seed", paste(deparse(body(interface_env[[name]])), collapse = " "),
        fixed = TRUE)
}, logical(1))))

# Report-level stochastic blocks must use the shared analysis seed. The
# permutation helper is shared by the DAG report and its Table S1 vignette, so
# seeding inside the helper also makes their independently rendered p-values
# identical.
assoc_env <- new.env(parent = globalenv())
sys.source("scripts/R/functions/assoc_tests.R", envir = assoc_env)

stopifnot(
  identical(formals(assoc_env$run_permutation_test)$seed, 123L),
  grepl(
    "set.seed(seed)",
    paste(deparse(body(assoc_env$run_permutation_test)), collapse = " "),
    fixed = TRUE
  )
)

dag_report <- paste(
  readLines("scripts/QMD/report/dag_associations.qmd", warn = FALSE),
  collapse = "\n"
)
dag_table <- paste(
  readLines("scripts/QMD/vignettes/clin_effect_table_s2.qmd", warn = FALSE),
  collapse = "\n"
)
figure_one <- paste(
  readLines("scripts/QMD/vignettes/figure1.qmd", warn = FALSE),
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
