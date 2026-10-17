# ===============================================================================
# CACHE PATHS AND FINGERPRINTS OF THE PAPER 3 VOI ENGINES (issues #202, #205)
# ===============================================================================
# The implementation of every value-of-information cache is the `voi` category
# of calculation_identity() (R/cache_provenance.R: VOI_CALCULATION_FILES, every
# paper 3 calculation source including this file and analysis/14_VOI.R), so an
# edit of any of them invalidates every VoI result and none of the PSA, EVPPI
# or scenario caches. Each fingerprint also records, under `upstream`, the
# fingerprints of the caches it was built from (voi_upstream_fingerprints());
# the loaders name the differing inputs (check_voi_identity()).
#
# Every cache written here lives next to the base caches (cache_dir()) under
# its own name; save_evidence_state() refuses to write to the base sampling
# cache path, so an evidence state can never replace E0.
#
# Requires R/cache_paths.R (cache_dir(), cache_fingerprint(),
# sampling_cache_path(), and through it R/cache_provenance.R).
# ===============================================================================

#' Calculation identity of the VoI engines
#'
#' Thin wrapper kept for callers: the single definition is the `voi` category
#' of calculation_identity() (issue #205).
#'
#' @return List with `code` (hash per file), `packages`, `R`, `contrasts` and
#'   `rng`.
voi_engine_identity <- function() calculation_identity("voi")

#' Hash of a reference population used by the maps g_k and the designs
#'
#' @param reference_population Data frame of reference patients.
#' @param weights Fixed weights, or NULL.
#' @return Character hash (cache_fingerprint() of the columns and weights;
#'   values never leave the hash).
reference_population_hash <- function(reference_population, weights = NULL) {
  cache_fingerprint(list(population = as.data.frame(reference_population), weights = weights))
}

#' Fingerprint of an evidence state
#'
#' The upstream cache of an evidence state is the base (E0) sampling cache,
#' recorded as `base_sampling_fingerprint`; the implementation is the `voi`
#' category of calculation_identity() (#205).
#'
#' @param base_fingerprint Fingerprint of the base (E0) sampling cache.
#' @param label Evidence-state label.
#' @param evidence Evidence rows used (data frame or NULL).
#' @param evidence_md5 md5 of the evidence file.
#' @param specs List of trial_summary_spec() objects.
#' @param sceptical Sceptical-prior specification or NULL.
#' @param tau Heterogeneity settings.
#' @param reference_population,weights Reference patients and weights.
#' @param n_samples,seed Draws of the regenerated component.
#' @param tol,step Relinearisation tolerance and quadrature spacing.
#' @return List with `fingerprint` and its `inputs`.
evidence_state_fingerprint <- function(base_fingerprint, label, evidence, evidence_md5, specs, sceptical,
                                       tau, reference_population, weights, n_samples, seed, tol, step) {
  inputs <- list(schema = "voi_evidence_state_v2", base_sampling_fingerprint = base_fingerprint,
                 label = label, evidence = if (is.null(evidence)) NULL else as.data.frame(evidence),
                 evidence_md5 = evidence_md5, specs = lapply(specs, unclass), sceptical = sceptical,
                 tau = tau, reference_population = reference_population_hash(reference_population, weights),
                 n_samples = as.integer(n_samples), seed = as.integer(seed), tol = tol, step = step,
                 implementation = calculation_identity("voi"))
  list(fingerprint = cache_fingerprint(inputs), inputs = inputs)
}

#' Path of an evidence-state sampling cache
#'
#' @param label Evidence-state label (letters, digits, "_" and "-").
#' @param n Number of draws.
#' @param directory Cache directory.
#' @return File path "sampling_models_n<n>_evidence_<label>.rds".
evidence_state_cache_path <- function(label, n, directory = cache_dir()) {
  if (!grepl("^[A-Za-z0-9_-]+$", label)) stop("Evidence-state labels may contain letters, digits, '_' and '-' only")
  file.path(directory, paste0("sampling_models_n", n, "_evidence_", label, ".rds"))
}

#' Save an evidence state without touching the base sampling cache
#'
#' @param state `state` element of build_evidence_state().
#' @param path Destination (default evidence_state_cache_path()).
#' @return The path, invisibly.
save_evidence_state <- function(state, path = evidence_state_cache_path(state$evidence_state,
                                                                        state$joint$n_samples)) {
  v_base <- normalizePath(sampling_cache_path(state$joint$n_samples, dirname(path)), mustWork = FALSE)
  if (identical(normalizePath(path, mustWork = FALSE), v_base)) {
    stop("Refusing to overwrite the base sampling cache with an evidence state")
  }
  if (!identical(state$joint$method, "mvn_joint_v2") || is.null(state$fingerprint)) {
    stop("An evidence state needs method mvn_joint_v2 and a fingerprint")
  }
  if (!dir.exists(dirname(path))) dir.create(dirname(path), recursive = TRUE)
  saveRDS(state, path)
  invisible(path)
}

#' Load an evidence state and check its fingerprint
#'
#' @param path Cache path.
#' @param expected Expected identity: evidence_state_fingerprint() built from
#'   the current base sampling cache (a list with `fingerprint` and `inputs`,
#'   so that the differing inputs are named), or its fingerprint alone.
#' @param required TRUE stops on an absent or stale state, naming the
#'   differing inputs; FALSE (the default, for the producer that rebuilds it)
#'   returns NULL with a message.
#' @return The state, or NULL (with a message) when absent or stale and not
#'   required.
load_evidence_state <- function(path, expected, required = FALSE) {
  instruction <- "Rebuild it with build_evidence_state() and save_evidence_state()."
  if (!file.exists(path)) {
    if (required) stop("Missing evidence-state cache ", basename(path), ". ", instruction, call. = FALSE)
    message("No evidence-state cache at ", path)
    return(NULL)
  }
  state <- readRDS(path)
  if (!is.list(expected)) expected <- list(fingerprint = expected, inputs = NULL)
  saved_inputs <- if (is.null(expected$inputs)) NULL else state$fingerprint_inputs
  if (!check_voi_identity(state$fingerprint, saved_inputs, expected,
                          paste0("evidence-state cache ", basename(path)), instruction, required)) {
    return(NULL)
  }
  state
}

#' Path of an EVSI results cache
#'
#' Also the path of the research-value results of script 14
#' (`evsi_results_research_<label>.rds`).
#'
#' @param label Label (e.g. the evidence state and price).
#' @param directory Cache directory.
#' @return File path "evsi_results_<label>.rds".
evsi_cache_path <- function(label, directory = cache_dir()) {
  if (!grepl("^[A-Za-z0-9_-]+$", label)) stop("EVSI cache labels may contain letters, digits, '_' and '-' only")
  file.path(directory, paste0("evsi_results_", label, ".rds"))
}

#' Fingerprint of an EVSI run
#'
#' @param upstream Fingerprints of the caches the NMB matrix and the draws come
#'   from (voi_upstream_fingerprints(); must include `psa`, and
#'   `evidence_states` for an evidence-state run).
#' @param designs List of design specifications (evsi_design()).
#' @param n_grid Sample sizes.
#' @param wtp,price Threshold and nivolumab price of the NMB matrix.
#' @param seed,B Seed and regression-coefficient draws of the SEs.
#' @param methods EVSI methods.
#' @return List with `fingerprint` and its `inputs`.
evsi_fingerprint <- function(upstream, designs, n_grid, wtp, price, seed, B, methods) {
  if (!is.list(upstream) || is.null(upstream$psa)) {
    stop("evsi_fingerprint() needs upstream$psa (voi_upstream_fingerprints(psa_obj = ...))", call. = FALSE)
  }
  inputs <- list(schema = "voi_evsi_v2", upstream = upstream,
                 designs = lapply(designs, function(d) list(name = d$name, targets = d$targets,
                                                             covariances = d$covariances)),
                 n_grid = n_grid, wtp = wtp, price = price, seed = as.integer(seed), B = as.integer(B),
                 methods = methods, implementation = calculation_identity("voi"))
  list(fingerprint = cache_fingerprint(inputs), inputs = inputs)
}

#' Save value-of-information results with their identity
#'
#' Layout of the EVSI and research-value caches: `fingerprint`, `inputs`
#' (including `upstream` and `implementation`), `results` and
#' `creation_time`.
#'
#' @param results Results object.
#' @param identity Identity (list with `fingerprint` and `inputs`), e.g. from
#'   evsi_fingerprint() or voi_research_fingerprint().
#' @param path Cache path (evsi_cache_path()).
#' @return The path, invisibly.
save_voi_results <- function(results, identity, path) {
  if (is.null(identity$fingerprint) || is.null(identity$inputs$upstream)) {
    stop("Value-of-information results need an identity with upstream fingerprints (#205)", call. = FALSE)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(fingerprint = identity$fingerprint, inputs = identity$inputs, results = results,
               creation_time = Sys.time()), path)
  invisible(path)
}

#' Load value-of-information results and validate their identity
#'
#' Compares the saved identity with `expected`, built from the current caches
#' and code, and names the differing inputs (e.g. `upstream$psa$generation_id`).
#'
#' @param path Cache path.
#' @param expected Current identity (list with `fingerprint` and `inputs`).
#' @param what Description of the results for the message.
#' @param instruction What to do about a missing or stale cache.
#' @param required TRUE (readers) stops on an absent or stale cache; FALSE
#'   (the producer, which then recomputes) returns NULL, with a message when
#'   the cache is stale.
#' @return The saved `results`, or NULL when absent or stale and not required.
load_voi_results <- function(path, expected, what = "value-of-information results",
                             instruction = "Rerun analysis/14_VOI.R.", required = TRUE) {
  if (!file.exists(path)) {
    if (required) stop("Missing ", what, " (", basename(path), "). ", instruction, call. = FALSE)
    return(NULL)
  }
  saved <- readRDS(path)
  if (!check_voi_identity(saved$fingerprint, saved$inputs, expected,
                          paste0(what, " (", basename(path), ")"), instruction, required)) {
    return(NULL)
  }
  saved$results
}
