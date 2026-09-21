# Cache contracts shared by producers and read-only consumers (issues #167/169).
# Canonical patient inputs remain in the ignored local cache, never in logs.
canonical_sampling_data <- function(data) {
  lapply(data, function(x) {
    if (is.factor(x)) return(list(values = enc2utf8(as.character(x)),
      levels = enc2utf8(levels(x)), ordered = is.ordered(x),
      contrasts = attr(x, "contrasts")))
    if (is.character(x)) return(enc2utf8(x))
    if (is.integer(x) || is.double(x)) return(as.numeric(x))
    if (is.logical(x)) return(as.logical(x))
    stop("Unsupported sampling column class: ", paste(class(x), collapse = ", "))
  })
}

cache_runtime_info <- function() {
  list(R = R.version.string, libraries = .libPaths(),
       locale = Sys.getlocale(), session = utils::sessionInfo())
}

# Parse source rather than hash bytes: comments, line endings and srcrefs do not
# change calculations. Include loaded definitions to detect interactive edits.
calculation_identity <- function(stage = c("psa", "sampling", "evppi", "scenario")) {
  stage <- match.arg(stage)
  files <- switch(stage,
    sampling = c("analysis/06_sampling.R", "R/model_configs.R", "R/joint_survival_sampling.R"),
    psa = c("R/model_fun.R", "R/calculate_outcomes.R", "R/prediction_functions.R",
            "R/psa_functions.R", "R/parameter_distributions.R", "R/model_configs.R",
            "R/prediction_population.R",
            "analysis/06_sampling.R", "R/joint_survival_sampling.R"),
    evppi = "R/evppi_functions.R",
    scenario = c("R/scenario_analysis.R", "R/evppi_functions.R"))
  code <- lapply(files, function(file) {
    exprs <- parse(here::here(file), keep.source = FALSE)
    definitions <- list()
    for (expr in exprs) {
      if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
          is.symbol(expr[[2]]) && is.call(expr[[3]]) &&
          identical(expr[[3]][[1]], as.name("function"))) {
        name <- as.character(expr[[2]])
        fn <- get0(name, envir = .GlobalEnv, mode = "function", inherits = FALSE)
        if (is.null(fn)) fn <- eval(expr[[3]], envir = baseenv())
        definitions[[name]] <- list(formals = paste(deparse(formals(fn)), collapse = "\n"),
                                    body = paste(deparse(body(fn)), collapse = "\n"))
      }
    }
    list(source = paste(deparse(exprs), collapse = "\n"), definitions = definitions)
  })
  names(code) <- files
  packages <- switch(stage, sampling = c("flexsurv", "survival", "mvtnorm"),
    psa = c("flexsurv", "survival", "dampack", "mvtnorm"),
    evppi = c("voi", "mgcv", "dampack"), scenario = c("voi", "mgcv", "dampack"))
  list(code = vapply(code, cache_fingerprint, character(1)),
       packages = setNames(vapply(packages, function(p)
         as.character(utils::packageVersion(p)), character(1)), packages),
       R = paste(R.version$major, R.version$minor, sep = "."),
       contrasts = getOption("contrasts"), rng = RNGkind())
}

report_fingerprint_differences <- function(cached_inputs, expected_inputs) {
  for (key in union(names(expected_inputs), names(cached_inputs))) {
    if (!identical(cached_inputs[[key]], expected_inputs[[key]])) {
      cat("  - ", key, " differs\n", sep = "")
      if (key == "data_columns") {
        for (column in union(names(cached_inputs[[key]]), names(expected_inputs[[key]]))) {
          old <- cached_inputs[[key]][column]
          new <- expected_inputs[[key]][column]
          if (!identical(old, new)) cat("    ", column, ": ", old, " -> ", new, "\n", sep = "")
        }
      }
    }
  }
}

assert_sampling_regeneration_allowed <- function() {
  if (!isTRUE(getOption("metimmox.sampling_allow_regenerate", TRUE)) ||
      isTRUE(getOption("knitr.in.progress"))) {
    stop("Sampling cache missing or stale; regeneration is disabled for this reader. ",
         "Run analysis/06_sampling.R in the analysis pipeline first.", call. = FALSE)
  }
}

psa_result_hash <- function(psa_obj) cache_fingerprint(list(
  cost = psa_obj$cost, effect = psa_obj$effect, strategies = psa_obj$strategies,
  sim = psa_obj$sim, model_idx = psa_obj$model_idx))

bind_psa_pair <- function(psa_obj, psa_params) {
  psa_obj$sim <- psa_params$sim
  metadata <- list(schema = "psa_pair_v1", fingerprint = psa_obj$fingerprint,
    params_hash = cache_fingerprint(psa_params), outcomes_hash = psa_result_hash(psa_obj))
  # A content-addressed generation ID binds the two files even after a partial save.
  metadata$generation_id <- cache_fingerprint(metadata)
  psa_obj$pair_metadata <- metadata
  attr(psa_params, "pair_metadata") <- metadata
  validate_psa_pair(psa_obj, psa_params)
  list(psa_obj = psa_obj, psa_params = psa_params)
}

validate_psa_pair <- function(psa_obj, psa_params) {
  fail <- function(why) stop("Invalid PSA cache pair: ", why,
                             ". Rerun analysis/10_PSA.R.", call. = FALSE)
  meta <- psa_obj$pair_metadata
  if (is.null(meta) || !identical(meta, attr(psa_params, "pair_metadata")))
    fail("missing or mixed generation metadata")
  if (!identical(meta$fingerprint, psa_obj$fingerprint) ||
      !identical(meta$params_hash, cache_fingerprint(psa_params)) ||
      !identical(meta$outcomes_hash, psa_result_hash(psa_obj))) fail("content hash mismatch")
  valid_ids <- function(x) is.numeric(x) && length(x) == psa_obj$n_sim &&
    all(is.finite(x) & x >= 1 & x == floor(x))
  if (!valid_ids(psa_params$sim) || anyDuplicated(psa_params$sim) ||
      !valid_ids(psa_params$model_idx) ||
      !identical(psa_obj$sim, psa_params$sim) ||
      !identical(psa_obj$model_idx, psa_params$model_idx) ||
      nrow(psa_params) != psa_obj$n_sim || nrow(psa_obj$cost) != psa_obj$n_sim ||
      nrow(psa_obj$effect) != psa_obj$n_sim) fail("draw/model indices are not aligned")
  invisible(TRUE)
}

evppi_cache_fingerprint <- function(psa_obj, wtp, config, seed, population,
                                    settings = list(method = "gam", B = 1000L, se = TRUE)) {
  inputs <- list(schema = "evppi_v2", psa_fingerprint = psa_obj$fingerprint,
    psa_generation = psa_obj$pair_metadata$generation_id,
    psa_result = psa_result_hash(psa_obj), wtp = wtp, config = config,
    seed = seed, population = population, settings = settings,
    implementation = calculation_identity("evppi"))
  list(fingerprint = cache_fingerprint(inputs), inputs = inputs)
}

scenario_cache_fingerprint <- function(psa_fingerprint, scenarios, config, seed, population) {
  inputs <- list(schema = "scenario_v2", psa_fingerprint = psa_fingerprint,
    scenarios = scenarios, config = config, seed = seed, population = population,
    implementation = calculation_identity("scenario"))
  list(fingerprint = cache_fingerprint(inputs), inputs = inputs)
}

# One complete-case economic population (#166); params also fingerprints the
# explicit prediction population, curve basis and target weights. Unused
# clinical-only columns stay outside this identity.
prediction_population_identity <- function() {
  variables <- unique(c("ID", unlist(lapply(c("os", "pfs"), function(outcome)
    all.vars(get_strategy_formula(get_biomarkers()[1], outcome))))))
  setNames(lapply(c("data_complete"), function(name) {
    population <- get0(name, envir = .GlobalEnv, inherits = FALSE)
    if (!is.data.frame(population)) return(NULL)
    columns <- sort(intersect(variables, names(population)))
    cache_fingerprint(canonical_sampling_data(population[, columns, drop = FALSE]))
  }), c("data_complete"))
}

current_psa_fingerprint <- function() psa_cache_fingerprint(
  sampling_models$fingerprint, l_params_base, param_distributions, strategies,
  n_sim, analysis_seed, time_horizon, cl)

current_population_inputs <- function() list(annual_incidence = annual_incidence_norway,
  research_horizon = research_horizon_years, discount_rate = discount_rate_research)

current_evppi_config <- function() {
  config <- configure_parameter_distributions(l_params_base)
  configure_evppi_analysis(config$distributions, config$groups)
}

current_scenario_fingerprint <- function() scenario_cache_fingerprint(
  current_psa_fingerprint()$fingerprint,
  define_scenarios(WTP, l_params_base$c_drug_nivo), current_evppi_config(),
  analysis_seed, current_population_inputs())

load_current_psa_cache <- function() load_psa_cache(required = TRUE,
  expected_fingerprint = current_psa_fingerprint(), strategies = strategies,
  n_sim = n_sim, sampling_fingerprint = sampling_models$fingerprint)

load_evppi_cache <- function(psa_obj, expected, wtp, path = evppi_path()) {
  if (!file.exists(path)) stop("Missing EVPPI cache. Run analysis/11_EVPPIs.R.")
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  if (!identical(e$evppi_fingerprint, expected$fingerprint))
    stop("Stale EVPPI cache: input fingerprint mismatch. Run analysis/11_EVPPIs.R.")
  nmb <- as.matrix(psa_obj$effect) * wtp - as.matrix(psa_obj$cost)
  evpi <- mean(apply(nmb, 1, max)) - max(colMeans(nmb))
  if (!isTRUE(all.equal(e$evpi_manual, evpi, tolerance = 1e-8)))
    stop("Cached EVPI disagrees with the current PSA and WTP. Run analysis/11_EVPPIs.R.")
  as.list(e)
}

load_scenario_cache <- function(expected = current_scenario_fingerprint(),
                                path = scenario_evppi_path()) {
  if (!file.exists(path)) stop("Missing scenario cache. Run analysis/12_scenario_EVPPIs.R.")
  cache <- readRDS(path)
  if (!identical(cache$fingerprint, expected$fingerprint))
    stop("Stale scenario cache: input fingerprint mismatch. Run analysis/12_scenario_EVPPIs.R.")
  if (!identical(cache$scenarios, expected$inputs$scenarios) ||
      !setequal(names(cache$all_scenario_results), cache$scenarios$scenario_id))
    stop("Scenario cache definitions/results do not match.")
  for (result in cache$all_scenario_results) validate_psa_pair(result$psa_obj, result$psa_params)
  if (!identical(cache$result_hash, scenario_result_hash(cache)))
    stop("Scenario cache result content does not match its saved identity. Run analysis/12_scenario_EVPPIs.R.")
  cache
}

scenario_result_hash <- function(cache) cache_fingerprint(list(
  results = cache$all_scenario_results, compiled = cache$evppi_all_scenarios))
