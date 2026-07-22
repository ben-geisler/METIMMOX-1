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
cache_dir <- function() here::here("data", "tidy")

#' Path to the snapshot directory or a file within it
#'
#' @param filename Optional snapshot filename. When omitted, returns the
#'   snapshot directory.
#' @param directory Directory containing snapshots.
#' @return Absolute path to the snapshot directory or requested file.
#' @export
snapshot_path <- function(filename = NULL,
                          directory = here::here("data", "output", "snapshots")) {
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

message("Cache-path helpers loaded from cache_paths.R")
