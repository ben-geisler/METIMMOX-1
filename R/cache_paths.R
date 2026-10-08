# ===============================================================================
# CACHE FILE PATHS - SINGLE SOURCE OF TRUTH
# ===============================================================================
# Centralizes the utility-source label and every cached-object file path so the
# analysis scripts and reports resolve identical locations. All managed caches
# live under data/tidy/. Sourced by 02_setup_and_global_variables.R (and, as a
# fallback, by cea_helpers.R).
#
# Cache-file naming:
#   sampling_models_n{n_samples}_full.rds     survival-sampling cache (08)
#   psa_obj_{label}.rds                        PSA dampack object       (12/15)
#   psa_params_{label}.rds                     PSA parameter draws      (12/15)
#   evppi_results_{label}.RData                EVPPI results            (13)
#   scenario_evppi_results_{label}.rds         scenario EVPPI           (14)
# where {label} is "ipd" or "correct".
# ===============================================================================

#' Directory holding tidy caches
#'
#' @return Absolute path to the tidy-cache directory.
#' @export
cache_dir <- function() getOption("metimmox.cache_dir", here::here("data", "tidy"))

#' Path to the snapshot directory or a file within it
#'
#' @param filename Optional snapshot filename. When omitted, returns the
#'   snapshot directory.
#' @param directory Directory containing snapshots.
#' @return Absolute path to the snapshot directory or requested file.
#' @export
snapshot_path <- function(filename = NULL,
                          directory = here::here("data", "snapshots")) {
  if (is.null(filename)) directory else file.path(directory, filename)
}

#' Resolve the utility-source label ("ipd" or "correct")
#'
#' @param label Explicit label, or NULL to use the global utility_source_label
#'   (falling back to "ipd" when unset).
#' @return Character scalar.
#' @export
resolve_util_label <- function(label = NULL) {
  if (!is.null(label)) return(label)
  if (exists("utility_source_label", inherits = TRUE)) {
    return(get("utility_source_label", inherits = TRUE))
  }
  "ipd"
}

#' Path to the correlated survival-sampling cache
#'
#' @param n Sample size; defaults to the global n_samples.
#' @param directory Directory containing the cache.
#' @export
sampling_cache_path <- function(n = NULL, directory = cache_dir()) {
  if (is.null(n)) {
    if (!exists("n_samples", inherits = TRUE)) {
      stop("n_samples is not set; pass n explicitly.")
    }
    n <- get("n_samples", inherits = TRUE)
  }
  file.path(directory, paste0("sampling_models_n", n, "_full.rds"))
}

#' Path to the PSA results cache (dampack PSA object)
#' @param label Utility-source label (see resolve_util_label).
#' @param directory Directory containing the cache.
#' @export
psa_obj_path <- function(label = NULL, directory = cache_dir()) {
  file.path(directory, paste0("psa_obj_", resolve_util_label(label), ".rds"))
}

#' Path to the PSA parameter-sample cache
#' @inheritParams psa_obj_path
#' @export
psa_params_path <- function(label = NULL, directory = cache_dir()) {
  file.path(directory, paste0("psa_params_", resolve_util_label(label), ".rds"))
}

#' Path to the EVPPI results cache
#' @inheritParams psa_obj_path
#' @export
evppi_path <- function(label = NULL, directory = cache_dir()) {
  file.path(directory, paste0("evppi_results_", resolve_util_label(label), ".RData"))
}

#' Path to the scenario-EVPPI results cache
#'
#' @param label Utility-source label (see resolve_util_label).
#' @param directory Directory containing the cache.
#' @export
scenario_evppi_path <- function(label = NULL, directory = cache_dir()) {
  file.path(
    directory,
    paste0("scenario_evppi_results_", resolve_util_label(label), ".rds")
  )
}

# ===============================================================================
# CACHE FINGERPRINTS AND PROVENANCE (issue #156)
# ===============================================================================
# The file name of a cache encodes only n_samples (sampling) or the utility
# source (PSA). Before issue #156 the load checks compared distribution names and
# the RNG seed, so a cache built from an earlier formula, data set or base
# parameter list was silently reused. Each cache now stores a fingerprint of the
# inputs that determine its content; the loading script recomputes it and
# regenerates on any mismatch.

#' Hash an arbitrary R object into a short, stable string
#'
#' @param x Object to hash (formulas are deparsed first so environments do
#'   not enter the hash; factors are hashed as character).
#' @return Character scalar.
#' @export
cache_fingerprint <- function(x) {
  canonical <- function(obj) {
    if (inherits(obj, "formula")) return(paste(deparse(obj), collapse = " "))
    # Radix (C-order) sorting keeps hashes independent of LC_COLLATE (#185).
    if (is.data.frame(obj)) {
      obj <- obj[, order(names(obj), method = "radix"), drop = FALSE]
      return(lapply(as.list(obj), function(col) {
        if (is.factor(col)) as.character(col) else col
      }))
    }
    if (is.list(obj)) {
      out <- lapply(obj, canonical)
      if (!is.null(names(obj))) out <- out[order(names(out), method = "radix")]
      return(out)
    }
    if (is.function(obj)) return(paste(deparse(obj), collapse = "\n"))
    obj
  }
  # Serialization v2 materialises ALTREP vectors, unlike hashing their internal
  # representation. Hash canonical values, independent of lazy materialisation.
  digest::digest(serialize(canonical(x), NULL, version = 2),
                 algo = "sha256", serialize = FALSE)
}

#' Fingerprint of the inputs that determine the survival-sampling cache
#'
#' @param formulas Named list of OS and PFS formulas.
#' @param data Data frame the models are fitted on.
#' @param distributions List with \code{os} and \code{pfs} distribution names.
#' @param n_samples Number of draws.
#' @param seed RNG seed.
#' @param method Sampling method label.
#' @return List with the \code{fingerprint} and its \code{inputs}.
#' @export
sampling_cache_fingerprint <- function(formulas, data, distributions,
                                       n_samples, seed, method) {
  # Radix (C-order) sorting: default collation depends on LC_COLLATE (#185).
  v_columns <- sort(method = "radix", unique(c(intersect("ID", names(data)),
                           unlist(lapply(formulas, all.vars)))))
  v_missing <- setdiff(v_columns, names(data))
  if (length(v_missing)) stop("Missing sampling columns: ", paste(v_missing, collapse = ", "))
  n_rows <- nrow(data)
  data <- canonical_sampling_data(data[, v_columns, drop = FALSE])
  inputs <- list(
    schema = "sampling_v2",
    formulas = lapply(formulas, function(f) paste(deparse(f), collapse = " ")),
    data_hash = cache_fingerprint(data),
    data_columns = vapply(data, cache_fingerprint, character(1)),
    n_rows = n_rows,
    distributions = distributions,
    n_samples = as.integer(n_samples),
    seed = as.integer(seed),
    method = method,
    implementation = calculation_identity("sampling")
  )
  list(fingerprint = cache_fingerprint(inputs), inputs = inputs,
       canonical_data = data, runtime = cache_runtime_info())
}

#' Fingerprint of the inputs that determine the PSA cache
#'
#' @param sampling_fingerprint Fingerprint of the sampling cache used.
#' @param l_params_base Base-case parameter list (every element enters the
#'   hash, including the survival curves and schedules).
#' @param param_distributions Named list of PSA distributions.
#' @param strategies Strategy identifiers.
#' @param n_sim Requested number of simulations.
#' @param seed PSA RNG seed.
#' @param time_horizon Model horizon in cycles.
#' @param cl Cycle length.
#' @return List with the \code{fingerprint} and its \code{inputs}.
#' @export
psa_cache_fingerprint <- function(sampling_fingerprint, l_params_base,
                                  param_distributions, strategies, n_sim,
                                  seed, time_horizon, cl) {
  inputs <- list(
    schema = "psa_v2",
    sampling_fingerprint = sampling_fingerprint,
    implementation = calculation_identity("psa"),
    prediction_population = prediction_population_identity(),
    params_hash = cache_fingerprint(l_params_base),
    distributions_hash = cache_fingerprint(param_distributions),
    strategies = strategies,
    n_sim = as.integer(n_sim),
    seed = as.integer(seed),
    time_horizon = time_horizon,
    cl = cl
  )
  list(fingerprint = cache_fingerprint(inputs), inputs = inputs)
}

#' Provenance (md5 and modification time) of a cache file
#'
#' @param path File path.
#' @return List with \code{path}, \code{exists}, \code{md5} and \code{mtime}
#'   (NA when the file is absent).
#' @export
cache_file_provenance <- function(path) {
  if (!file.exists(path)) {
    return(list(path = path, exists = FALSE, md5 = NA_character_,
                mtime = as.POSIXct(NA)))
  }
  list(
    path = path,
    exists = TRUE,
    md5 = unname(tools::md5sum(path)),
    mtime = file.info(path)$mtime
  )
}

source(here::here("R/cache_provenance.R"))
message("Cache-path helpers loaded from cache_paths.R")
