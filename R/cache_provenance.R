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

#' Paper 3 calculation sources: the `voi` category of calculation_identity()
#' (issue #205)
VOI_CALCULATION_FILES <- c("R/voi_price_axis.R", "R/voi_survival_parameters.R", "R/voi_policies.R",
                           "R/voi_incidence.R", "R/design_information.R", "R/evidence_update.R",
                           "R/evsi_functions.R", "R/voi_cache.R", "R/voi_designs.R",
                           "R/voi_research_value.R", "R/voi_evidence_psa.R", "analysis/14_VOI.R")

#' Packages whose versions enter the `voi` identity: the EVPPI and EVSI
#' regressions (voi, mgcv; earth and dbarts for the step 1a cross-checks), the
#' gamma likelihood and survival fits (flexsurv, survival) and the
#' multivariate normal draws (mvtnorm)
VOI_CALCULATION_PACKAGES <- c("voi", "mgcv", "earth", "dbarts", "flexsurv", "survival", "mvtnorm")

#' Calculation sources and packages of one cache category
#'
#' The single list of the code and packages hashed into each cache's
#' fingerprint. The `psa`, `sampling`, `evppi` and `scenario` lists are those
#' of issues #156-#169; `voi` (#205) lists every paper 3 calculation source, so
#' that only the value-of-information caches depend on them and an edit of
#' value-of-information code never invalidates the PSA, EVPPI or scenario
#' caches.
#'
#' @param stage Cache category: "psa", "sampling", "evppi", "scenario" or "voi".
#' @return List with `files` (paths relative to the repository root) and
#'   `packages` (names of the packages whose versions are recorded).
calculation_identity_sources <- function(stage = c("psa", "sampling", "evppi", "scenario", "voi")) {
  stage <- match.arg(stage)
  v_files <- switch(stage,
    sampling = c("analysis/06_sampling.R", "R/model_configs.R", "R/joint_survival_sampling.R"),
    psa = c("R/model_fun.R", "R/calculate_outcomes.R", "R/prediction_functions.R",
            "R/psa_functions.R", "R/parameter_distributions.R", "R/model_configs.R",
            "R/prediction_population.R",
            "analysis/06_sampling.R", "R/joint_survival_sampling.R"),
    evppi = "R/evppi_functions.R",
    scenario = c("R/scenario_analysis.R", "R/evppi_functions.R"),
    voi = VOI_CALCULATION_FILES)
  v_packages <- switch(stage, sampling = c("flexsurv", "survival", "mvtnorm"),
    psa = c("flexsurv", "survival", "dampack", "mvtnorm"),
    evppi = c("voi", "mgcv", "dampack"), scenario = c("voi", "mgcv", "dampack"),
    voi = VOI_CALCULATION_PACKAGES)
  list(files = v_files, packages = v_packages)
}

#' Calculation identity of one cache category
#'
#' Parses each source rather than hashing its bytes, so comments, line endings
#' and srcrefs do not change the identity, and hashes the loaded definition of
#' every function the file defines, so an interactive edit invalidates the
#' cache. Also records the package versions, the R version, the contrasts and
#' the RNG kind.
#'
#' @param stage Cache category (see calculation_identity_sources()).
#' @param root Repository root the source paths are relative to; the option
#'   `metimmox.calculation_root` overrides the default (the sensitivity tests
#'   hash an edited copy of the sources through it).
#' @return List with `code` (one hash per file), `packages`, `R`, `contrasts`
#'   and `rng`.
calculation_identity <- function(stage = c("psa", "sampling", "evppi", "scenario", "voi"),
                                 root = getOption("metimmox.calculation_root", here::here())) {
  stage <- match.arg(stage)
  sources <- calculation_identity_sources(stage)
  v_files <- sources$files
  code <- lapply(v_files, function(file) {
    exprs <- parse(file.path(root, file), keep.source = FALSE)
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
  names(code) <- v_files
  v_packages <- sources$packages
  list(code = vapply(code, cache_fingerprint, character(1)),
       packages = setNames(vapply(v_packages, function(p)
         as.character(utils::packageVersion(p)), character(1)), v_packages),
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

#' Content hash of the outcomes of a PSA object
#'
#' Covers cost, QALYs, strategies, draw and model indices and, since issues
#' #197 and #200, the life-year matrices and the cell-level component.
#'
#' @param psa_obj PSA object.
#' @return Character fingerprint.
psa_result_hash <- function(psa_obj) cache_fingerprint(list(
  cost = psa_obj$cost, effect = psa_obj$effect, strategies = psa_obj$strategies,
  sim = psa_obj$sim, model_idx = psa_obj$model_idx,
  ly_disc = psa_obj$ly_disc, ly_undisc = psa_obj$ly_undisc, cells = psa_obj$cells))

# Schema of a bound PSA pair. psa_pair_v1 (#167) carried cost and QALYs only;
# psa_pair_v2 adds life-years (#197) and the cell-level outcomes (#200).
PSA_PAIR_SCHEMA <- "psa_pair_v2"

#' Describe why a PSA object lacks the extended outcomes of psa_pair_v2
#'
#' @param psa_obj PSA object.
#' @param n_sim Expected number of draws (rows).
#' @return NULL when life-years and cells are present and well formed,
#'   otherwise a phrase naming the problem.
psa_extended_component_problem <- function(psa_obj, n_sim = psa_obj$n_sim) {
  v_strategies <- psa_obj$strategies
  for (key in c("ly_disc", "ly_undisc")) {
    m_ly <- psa_obj[[key]]
    if (is.null(m_ly)) return(paste0("no ", key, " life-year matrix (#197): life-year draws missing"))
    # One finite column per strategy and one row per draw, aligned with the
    # cost table (issue #197).
    if (nrow(m_ly) != n_sim || !identical(colnames(m_ly), as.character(v_strategies)) ||
        (!is.null(psa_obj$cost) && !identical(colnames(m_ly), colnames(psa_obj$cost))) ||
        any(!is.finite(as.matrix(m_ly)))) {
      return(paste0(key, " life-year draws missing or not aligned with the cost draws ",
                    "(not a finite draw x strategy matrix)"))
    }
  }
  cells <- psa_obj$cells
  if (is.null(cells)) return("no cell-level component (#200)")
  if (!identical(cells$schema, "psa_cells_v1")) return("unknown cell-component schema")
  v_dim <- dim(cells$values)
  if (length(v_dim) != 4L || v_dim[1] != n_sim || v_dim[2] != length(cells$cells) ||
      !identical(dimnames(cells$values)[[2]], cells$cells) ||
      !identical(dimnames(cells$values)[[3]], cells$arms) ||
      !identical(dimnames(cells$values)[[4]], cells$outcomes))
    return("cell values are not a draw x cell x arm x outcome array")
  m_weights <- as.matrix(cells$weights)
  if (!identical(dim(m_weights), c(as.integer(n_sim), length(cells$cells))) ||
      any(!is.finite(m_weights)) || any(m_weights < 0) ||
      max(abs(rowSums(m_weights) - 1)) > 1e-10) return("cell weights are not draw-wise probabilities")
  v_observed <- colSums(m_weights > 0) > 0
  if (any(!is.finite(cells$values[, v_observed, , ]))) return("cell values are not finite")
  if (length(cells$nivolumab_price) != n_sim || any(!is.finite(cells$nivolumab_price)))
    return("cell nivolumab price is not one finite value per draw")
  NULL
}

#' Bind a PSA object and its parameter draws into one generation
#'
#' Refuses a PSA object without the life-years and cell-level outcomes of
#' psa_pair_v2: a legacy object is regenerated by analysis/10_PSA.R, never
#' stamped.
#'
#' @param psa_obj PSA object from build_psa_obj().
#' @param psa_params Retained PSA parameter rows.
#' @return List with the bound `psa_obj` and `psa_params`.
bind_psa_pair <- function(psa_obj, psa_params) {
  problem <- psa_extended_component_problem(psa_obj)
  if (!is.null(problem)) {
    stop("Refusing to bind a legacy PSA object: ", problem,
         ". Rerun analysis/10_PSA.R with the current model.", call. = FALSE)
  }
  psa_obj$sim <- psa_params$sim
  metadata <- list(schema = PSA_PAIR_SCHEMA, fingerprint = psa_obj$fingerprint,
    params_hash = cache_fingerprint(psa_params), outcomes_hash = psa_result_hash(psa_obj))
  # A content-addressed generation ID binds the two files even after a partial save.
  metadata$generation_id <- cache_fingerprint(metadata)
  psa_obj$pair_metadata <- metadata
  attr(psa_params, "pair_metadata") <- metadata
  validate_psa_pair(psa_obj, psa_params)
  list(psa_obj = psa_obj, psa_params = psa_params)
}

#' Validate a bound PSA pair
#'
#' Stops on missing or mixed generation metadata, a legacy pair (schema
#' psa_pair_v1, without life-years and cell-level outcomes), a content-hash
#' mismatch or misaligned draw/model indices.
#'
#' @param psa_obj PSA object.
#' @param psa_params PSA parameter rows.
#' @return TRUE invisibly.
validate_psa_pair <- function(psa_obj, psa_params) {
  fail <- function(why) stop("Invalid PSA cache pair: ", why,
                             ". Rerun analysis/10_PSA.R.", call. = FALSE)
  meta <- psa_obj$pair_metadata
  if (is.null(meta) || !identical(meta, attr(psa_params, "pair_metadata")))
    fail("missing or mixed generation metadata")
  if (!identical(meta$schema, PSA_PAIR_SCHEMA))
    fail(paste0("legacy pair (schema ", meta$schema, ") without the life-years (#197) ",
                "and cell-level outcomes (#200) of ", PSA_PAIR_SCHEMA))
  problem <- psa_extended_component_problem(psa_obj)
  if (!is.null(problem)) fail(problem)
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
  v_variables <- unique(c("ID", unlist(lapply(c("os", "pfs"), function(outcome)
    all.vars(get_strategy_formula(get_biomarkers()[1], outcome))))))
  setNames(lapply(c("data_complete"), function(name) {
    df_population <- get0(name, envir = .GlobalEnv, inherits = FALSE)
    if (!is.data.frame(df_population)) return(NULL)
    v_columns <- sort(intersect(v_variables, names(df_population)), method = "radix")
    cache_fingerprint(canonical_sampling_data(df_population[, v_columns, drop = FALSE]))
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
  m_nmb <- as.matrix(psa_obj$effect) * wtp - as.matrix(psa_obj$cost)
  evpi <- mean(apply(m_nmb, 1, max)) - max(colMeans(m_nmb))
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

# ===============================================================================
# UPSTREAM PROVENANCE OF THE VALUE-OF-INFORMATION CACHES (issue #205)
# ===============================================================================
# Every paper 3 cache (step 1a results, evidence states, EVSI and research-value
# results) records, under `upstream` in its fingerprint inputs, the fingerprints
# of the caches it was built from, and its implementation through the `voi`
# category of calculation_identity(). Loaders compare the saved identity with
# the one built from the current caches and name the differing inputs, like
# load_evppi_cache() and report_fingerprint_differences().

#' Fingerprints of the caches a value-of-information result is built from
#'
#' @param sampling_models Sampling cache (script 06), or NULL.
#' @param psa_obj Bound PSA object (psa_pair_v2), or NULL.
#' @param scenario_cache Scenario cache (script 12), or NULL.
#' @param evidence_states Named list of evidence states (each with a
#'   `fingerprint`) or named character vector of their fingerprints, or NULL.
#' @return Named list with the elements supplied: `sampling` (fingerprint),
#'   `psa` (`fingerprint`, `generation_id`, `schema`), `scenario`
#'   (`fingerprint`, `result_hash` and `psa_generations`, the generation ID of
#'   every embedded PSA pair by scenario id) and `evidence_states`
#'   (fingerprint by label). Stops when a supplied cache carries no
#'   fingerprint.
voi_upstream_fingerprints <- function(sampling_models = NULL, psa_obj = NULL, scenario_cache = NULL,
                                      evidence_states = NULL) {
  need <- function(value, what, instruction) {
    if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
      stop("Value-of-information input without provenance: ", what, " is missing. ", instruction,
           call. = FALSE)
    }
    value
  }
  upstream <- list()
  if (!is.null(sampling_models)) {
    upstream$sampling <- need(sampling_models$fingerprint, "the sampling-cache fingerprint",
                              "Rerun analysis/06_sampling.R.")
  }
  if (!is.null(psa_obj)) {
    instruction <- "Rerun analysis/10_PSA.R."
    upstream$psa <- list(
      fingerprint = need(psa_obj$fingerprint, "the PSA input fingerprint", instruction),
      generation_id = need(psa_obj$pair_metadata$generation_id, "the PSA pair generation ID", instruction),
      schema = need(psa_obj$pair_metadata$schema, "the PSA pair schema", instruction))
  }
  if (!is.null(scenario_cache)) {
    instruction <- "Rerun analysis/12_scenario_EVPPIs.R."
    upstream$scenario <- list(
      fingerprint = need(scenario_cache$fingerprint, "the scenario-cache fingerprint", instruction),
      result_hash = need(scenario_cache$result_hash, "the scenario-cache result hash", instruction),
      psa_generations = vapply(scenario_cache$all_scenario_results, function(result)
        need(result$psa_obj$pair_metadata$generation_id, "a scenario PSA pair generation ID", instruction),
        character(1)))
  }
  if (!is.null(evidence_states)) {
    if (is.null(names(evidence_states)) || any(!nzchar(names(evidence_states)))) {
      stop("evidence_states must be named by evidence-state label", call. = FALSE)
    }
    upstream$evidence_states <- vapply(evidence_states, function(state)
      need(if (is.list(state)) state$fingerprint else state, "an evidence-state fingerprint",
           "Rebuild the evidence state (build_evidence_state())."), character(1))
  }
  if (!length(upstream)) stop("voi_upstream_fingerprints() needs at least one upstream cache", call. = FALSE)
  upstream
}

#' Names of the inputs in which two cache identities differ
#'
#' Descends into named lists and named vectors (`upstream`, `implementation`,
#' the per-file code hashes), so a difference is named down to the cache or
#' file, e.g. `upstream$psa$generation_id` or
#' `implementation$voi$code$R/voi_policies.R`.
#'
#' @param saved Inputs recorded in the cache.
#' @param expected Inputs built from the current caches and code.
#' @param prefix Prefix of the names (used in the recursion).
#' @return Character vector of differing input names (empty when identical).
voi_identity_differences <- function(saved, expected, prefix = "") {
  is_nested <- function(x) !is.null(names(x)) && !is.data.frame(x) &&
    (is.list(x) || (is.atomic(x) && length(x) >= 1L))
  v_differences <- character(0)
  for (key in union(names(expected), names(saved))) {
    old <- saved[[key]]
    new <- expected[[key]]
    if (identical(old, new)) next
    label <- paste0(prefix, key)
    v_differences <- c(v_differences, if (is_nested(old) && is_nested(new)) {
      voi_identity_differences(as.list(old), as.list(new), paste0(label, "$"))
    } else label)
  }
  v_differences
}

#' Validate the saved identity of a value-of-information cache
#'
#' @param saved_fingerprint Fingerprint stored in the cache.
#' @param saved_inputs Fingerprint inputs stored in the cache (NULL when the
#'   cache recorded none).
#' @param expected Identity built from the current caches and code (list with
#'   `fingerprint` and `inputs`, as returned by the VoI fingerprint functions).
#' @param what Description of the cache for the message.
#' @param instruction What to do about a stale cache (appended to the message).
#' @param required TRUE stops on a stale cache; FALSE gives a message and
#'   returns FALSE (for producers that then recompute).
#' @return TRUE when current; FALSE (with a message) when stale and not
#'   required.
check_voi_identity <- function(saved_fingerprint, saved_inputs, expected, what, instruction,
                               required = TRUE) {
  why <- NULL
  if (!is.null(saved_inputs) && !identical(cache_fingerprint(saved_inputs), saved_fingerprint)) {
    why <- "its recorded inputs do not match its saved fingerprint"
  } else if (!identical(saved_fingerprint, expected$fingerprint)) {
    v_differences <- if (is.null(saved_inputs)) character(0) else
      voi_identity_differences(saved_inputs, expected$inputs)
    why <- if (length(v_differences)) paste0("inputs differ from the current ones: ",
                                             paste(v_differences, collapse = ", ")) else
      "fingerprint differs (no recorded inputs to compare)"
  }
  if (is.null(why)) return(TRUE)
  msg <- paste0("Stale ", what, ": ", why, ". ", instruction)
  if (required) stop(msg, call. = FALSE)
  message(msg)
  FALSE
}
