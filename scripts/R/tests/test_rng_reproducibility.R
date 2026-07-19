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
sampling_expressions <- parse("scripts/R/analysis/08_sampling.R")
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

# Capture the random row subset used by grouped EVPPI calculations.
evppi_expressions <- parse("scripts/R/functions/evppi_functions.R")
is_evppi_function <- vapply(evppi_expressions, function(expr) {
  is.call(expr) &&
    identical(expr[[1]], as.name("<-")) &&
    identical(expr[[2]], as.name("calculate_evppi_improved"))
}, logical(1))
stopifnot(sum(is_evppi_function) == 1L)

evppi_env <- new.env(parent = globalenv())
evppi_env$sampled_indices <- NULL
evppi_env$sample <- eval(quote(function(x, size, ...) {
  indices <- base::sample(x, size, ...)
  sampled_indices <<- indices
  indices
}), envir = evppi_env)
eval(evppi_expressions[[which(is_evppi_function)]], envir = evppi_env)

n_evppi <- 30L
effect_one <- seq(0, 1, length.out = n_evppi)
psa_obj <- list(
  cost = cbind(strategy_one = rep(0, n_evppi),
               strategy_two = rep(0, n_evppi)),
  effect = cbind(strategy_one = effect_one,
                 strategy_two = rev(effect_one))
)
psa_params <- data.frame(
  parameter_one = effect_one,
  parameter_two = cos(seq(0, 2 * pi, length.out = n_evppi))
)

run_grouped_evppi <- function(seed) {
  runif(50)
  capture.output(evppi_env$calculate_evppi_improved(
    psa_obj = psa_obj,
    psa_params = psa_params,
    param_names = names(psa_params),
    wtp = 1,
    n_inner = 5,
    n_grid = 8,
    seed = seed
  ))
  evppi_env$sampled_indices
}

evppi_reference <- run_grouped_evppi(123L)
evppi_repeat <- run_grouped_evppi(123L)
evppi_other_seed <- run_grouped_evppi(124L)

stopifnot(
  identical(evppi_reference, evppi_repeat),
  !identical(evppi_reference, evppi_other_seed)
)

# Public wrappers must retain the documented default and forward it to their
# stochastic dependencies.
interface_env <- new.env(parent = globalenv())
sys.source("scripts/R/functions/evppi_functions.R", envir = interface_env)
sys.source("scripts/R/functions/scenario_analysis.R", envir = interface_env)

seeded_interfaces <- c(
  "calculate_evppi_improved",
  "run_evppi_analysis",
  "run_scenario_psa",
  "run_scenario_evppi",
  "run_scenario_analysis",
  "run_all_scenarios"
)
stopifnot(all(vapply(seeded_interfaces, function(name) {
  identical(formals(interface_env[[name]])$seed, 123L)
}, logical(1))))

forwarding_wrappers <- setdiff(seeded_interfaces, "calculate_evppi_improved")
stopifnot(all(vapply(forwarding_wrappers, function(name) {
  grepl("seed = seed", paste(deparse(body(interface_env[[name]])), collapse = " "),
        fixed = TRUE)
}, logical(1))))

cat("RNG reproducibility and cache seed validation tests passed.\n")
