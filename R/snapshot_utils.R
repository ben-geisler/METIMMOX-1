# Snapshot Utility Functions
# Helper functions for managing single-model result snapshots for bug fix impact analysis.

if (!exists("snapshot_path", mode = "function")) {
  source(here::here("R/cache_paths.R"))
}

#' Get Git Commit Hash
#'
#' Retrieves the current git commit hash (short version).
#'
#' @return Character string with short commit hash, or "unknown" if unavailable.
get_git_commit <- function() {
  tryCatch({
    system("git rev-parse --short HEAD", intern = TRUE)
  }, error = function(e) {
    warning("Could not retrieve git commit hash")
    "unknown"
  })
}

#' Get Git Commit Time
#'
#' @return POSIXct commit time of HEAD, or NA if unavailable.
get_git_commit_time <- function() {
  tryCatch({
    epoch <- suppressWarnings(as.numeric(system("git log -1 --format=%ct HEAD", intern = TRUE)))
    if (length(epoch) != 1 || is.na(epoch)) return(as.POSIXct(NA))
    as.POSIXct(epoch, origin = "1970-01-01", tz = "")
  }, error = function(e) as.POSIXct(NA))
}

#' Collect cache provenance for a snapshot
#'
#' @param label Utility-source label.
#' @param n Sampling-cache sample size.
#' @return Named list of cache_file_provenance() results for the sampling, PSA
#'   object, PSA parameter, EVPPI and scenario-EVPPI caches.
collect_cache_provenance <- function(label = NULL, n = NULL) {
  list(
    sampling = cache_file_provenance(sampling_cache_path(n)),
    psa_obj = cache_file_provenance(psa_obj_path(label)),
    psa_params = cache_file_provenance(psa_params_path(label)),
    evppi = cache_file_provenance(evppi_path(label)),
    scenario_evppi = cache_file_provenance(scenario_evppi_path(label))
  )
}

#' Describe the PSA provenance of a before/after snapshot pair
#'
#' @param before,after Snapshot objects.
#' @return Character scalar suitable for printing in a comparison.
describe_psa_provenance <- function(before, after) {
  b <- before$metadata$caches
  a <- after$metadata$caches
  describe_one <- function(x, label) {
    if (is.null(x)) {
      return(sprintf("%s snapshot predates issue #156 and did not record the PSA cache md5",
                     label))
    }
    sprintf("%s PSA cache md5 %s (modified %s%s)", label, x$psa_obj$md5,
            format(x$psa_obj$mtime, "%Y-%m-%d %H:%M"),
            if (isTRUE(x$psa_regenerated)) ", regenerated inside the snapshot run" else "")
  }
  if (is.null(b) || is.null(a)) {
    return(paste0("PSA provenance: ", describe_one(b, "before"), "; ",
                  describe_one(a, "after"), ". Identity of the two PSA files ",
                  "cannot be established from metadata."))
  }
  same <- identical(b$psa_obj$md5, a$psa_obj$md5)
  regen <- isTRUE(a$psa_regenerated)
  msg <- sprintf(
    "PSA cache md5 before %s (modified %s), after %s (modified %s).",
    b$psa_obj$md5, format(b$psa_obj$mtime, "%Y-%m-%d %H:%M"),
    a$psa_obj$md5, format(a$psa_obj$mtime, "%Y-%m-%d %H:%M")
  )
  if (same) {
    msg <- paste(msg, "The two snapshots copied the same PSA cache, so any",
                 "PSA comparison shows no change by construction.")
  } else if (regen) {
    msg <- paste(msg, "The after PSA was regenerated inside the snapshot run",
                 "because the cache predated the fix commit.")
  }
  msg
}

#' Get Package Versions
#'
#' Retrieves version numbers for key R packages used in the analysis.
#'
#' @return Named list of package versions.
get_package_versions <- function() {
  key_packages <- c(
    "dampack", "flexsurv", "survival", "dplyr", "ggplot2",
    "here", "mvtnorm", "Matrix", "darthtools", "gems"
  )

  versions <- list()
  for (pkg in key_packages) {
    versions[[pkg]] <- if (requireNamespace(pkg, quietly = TRUE)) {
      as.character(packageVersion(pkg))
    } else {
      "not installed"
    }
  }

  versions
}

#' Get Current Economic Model Metadata
#'
#' @return List describing the single economic model, or NULL if unavailable.
get_economic_model_metadata <- function() {
  if (!exists("get_current_model_config", mode = "function", inherits = TRUE)) {
    return(NULL)
  }

  config <- get_current_model_config()
  list(
    label = config$label,
    description = config$description,
    strategies = config$strategies,
    biomarkers = config$biomarkers
  )
}

#' Collect Metadata
#'
#' Collects comprehensive metadata about the current single-model analysis run.
#'
#' @param issue_number GitHub issue number (character or numeric).
#' @return List containing metadata.
collect_metadata <- function(issue_number, caches = NULL) {
  list(
    issue_number = as.character(issue_number),
    git_commit = get_git_commit(),
    git_commit_time = get_git_commit_time(),
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    # Cache provenance (issue #156): md5 and modification time of every cache
    # the snapshot copied, plus whether the PSA/EVPPI were regenerated inside
    # the snapshot run. Without this, a "fixed" snapshot could carry a PSA file
    # copied unchanged from before the fix and report "no PSA change".
    caches = caches,
    r_version = R.version.string,
    packages = get_package_versions(),
    model = get_economic_model_metadata(),
    parameters = list(
      WTP = if (exists("WTP")) WTP else NA,
      time_horizon = if (exists("time_horizon")) time_horizon else NA,
      n_samples = if (exists("n_samples")) n_samples else NA,
      n_sim = if (exists("n_sim")) n_sim else NA,
      cl = if (exists("cl")) cl else NA,
      dr = if (exists("dr")) dr else NA,
      DSA_mult = if (exists("DSA_mult")) DSA_mult else NA,
      UTILITY_SOURCE = if (exists("UTILITY_SOURCE")) UTILITY_SOURCE else NA,
      utility_source_label = if (exists("utility_source_label")) utility_source_label else NA
    )
  )
}

#' Generate Snapshot Filename
#'
#' Creates standardized filenames used by 13_save_snapshot.R.
#'
#' @param type File type ("snapshot" or "psa").
#' @param issue_number GitHub issue number.
#' @param status Snapshot status ("baseline" or "fixed").
#' @param commit Git commit hash.
#' @return Character string with filename.
generate_snapshot_filename <- function(type, issue_number, status,
                                       commit = get_git_commit()) {
  if (!type %in% c("snapshot", "psa")) {
    stop("type must be 'snapshot' or 'psa'.")
  }
  if (!status %in% c("baseline", "fixed")) {
    stop("status must be 'baseline' or 'fixed'.")
  }

  paste0(type, "_", issue_number, "_", status, "_", commit, ".rds")
}

#' Parse Snapshot Filename
#'
#' Supports current status-based filenames and older timestamp-based filenames.
#'
#' @param filename Snapshot or PSA filename.
#' @return Data frame with parsed fields.
parse_snapshot_filename <- function(filename) {
  base <- basename(filename)
  no_ext <- sub("\\.rds$", "", base)
  parts <- strsplit(no_ext, "_")[[1]]

  empty <- data.frame(
    filename = base,
    type = NA_character_,
    issue = NA_character_,
    status = NA_character_,
    timestamp = NA_character_,
    commit = NA_character_,
    stringsAsFactors = FALSE
  )

  if (length(parts) < 4 || !parts[1] %in% c("snapshot", "psa")) {
    return(empty)
  }

  type <- parts[1]
  issue <- parts[2]
  token <- parts[3]

  if (token %in% c("baseline", "fixed")) {
    status <- token
    timestamp <- NA_character_
    commit <- paste(parts[4:length(parts)], collapse = "_")
  } else if (length(parts) >= 5) {
    status <- NA_character_
    timestamp <- paste(parts[3], parts[4], sep = "_")
    commit <- paste(parts[5:length(parts)], collapse = "_")
  } else {
    status <- NA_character_
    timestamp <- NA_character_
    commit <- paste(parts[4:length(parts)], collapse = "_")
  }

  data.frame(
    filename = base,
    type = type,
    issue = issue,
    status = status,
    timestamp = timestamp,
    commit = commit,
    stringsAsFactors = FALSE
  )
}

#' List Snapshots for Issue
#'
#' Lists all snapshot and PSA files for a given issue number.
#'
#' @param issue_number GitHub issue number.
#' @param snapshots_dir Directory containing snapshots.
#' @return Data frame with snapshot information.
list_snapshots <- function(issue_number,
                           snapshots_dir = snapshot_path()) {
  if (!dir.exists(snapshots_dir)) {
    message("Snapshot directory not found: ", snapshots_dir)
    return(data.frame(
      filename = character(0),
      type = character(0),
      issue = character(0),
      status = character(0),
      timestamp = character(0),
      commit = character(0),
      mtime = as.POSIXct(character(0))
    ))
  }

  pattern <- paste0("^(snapshot|psa)_", issue_number, "_.*\\.rds$")
  files <- list.files(snapshots_dir, pattern = pattern, full.names = FALSE)

  if (length(files) == 0) {
    message("No snapshots found for issue #", issue_number)
    return(data.frame(
      filename = character(0),
      type = character(0),
      issue = character(0),
      status = character(0),
      timestamp = character(0),
      commit = character(0),
      mtime = as.POSIXct(character(0))
    ))
  }

  snapshot_info <- do.call(rbind, lapply(files, parse_snapshot_filename))
  snapshot_info$mtime <- file.info(file.path(snapshots_dir, snapshot_info$filename))$mtime
  snapshot_info <- snapshot_info[order(snapshot_info$mtime), ]
  rownames(snapshot_info) <- NULL

  snapshot_info
}

#' Load Snapshot
#'
#' Loads a snapshot file and validates its single-model structure.
#'
#' @param filename Snapshot filename without path.
#' @param snapshots_dir Directory containing snapshots.
#' @return Snapshot object.
load_snapshot <- function(filename,
                          snapshots_dir = snapshot_path()) {
  filepath <- file.path(snapshots_dir, filename)

  if (!file.exists(filepath)) {
    stop("Snapshot file not found: ", filepath)
  }

  snapshot <- readRDS(filepath)

  if (grepl("^snapshot_", filename)) {
    required_components <- c("metadata", "base_results", "icer_obj",
                             "nmb_at_wtp", "psa_summary")
    missing <- setdiff(required_components, names(snapshot))
    if (length(missing) > 0) {
      warning("Snapshot missing components: ", paste(missing, collapse = ", "))
    }
    if ("models" %in% names(snapshot)) {
      warning("Snapshot contains obsolete multi-model 'models' component.")
    }
  }

  if (grepl("^psa_", filename) && !inherits(snapshot, "psa")) {
    warning("PSA file does not contain a dampack psa object.")
  }

  snapshot
}

#' Find Matching Snapshot Pair
#'
#' Finds the matching PSA file for a snapshot file.
#'
#' @param snapshot_file Snapshot filename.
#' @param snapshots_dir Directory containing snapshots.
#' @return List with snapshot and psa filenames, or NULL if no match.
find_snapshot_pair <- function(snapshot_file,
                               snapshots_dir = snapshot_path()) {
  info <- parse_snapshot_filename(snapshot_file)

  if (is.na(info$type) || info$type != "snapshot") {
    stop("Expected a snapshot filename, got: ", snapshot_file)
  }

  psa_file <- if (!is.na(info$status)) {
    paste0("psa_", info$issue, "_", info$status, "_", info$commit, ".rds")
  } else if (!is.na(info$timestamp)) {
    paste0("psa_", info$issue, "_", info$timestamp, "_", info$commit, ".rds")
  } else {
    NULL
  }

  if (!is.null(psa_file) && file.exists(file.path(snapshots_dir, psa_file))) {
    return(list(snapshot = snapshot_file, psa = psa_file))
  }

  warning("No matching PSA file found for ", snapshot_file)
  NULL
}

#' Select Snapshots for Comparison
#'
#' Selects baseline and fixed single-model snapshot pairs for an issue.
#'
#' @param issue_number GitHub issue number.
#' @param snapshots_dir Directory containing snapshots.
#' @return List with before and after snapshot pairs.
select_snapshots_for_comparison <- function(issue_number,
                                            snapshots_dir = snapshot_path()) {
  snapshot_info <- list_snapshots(issue_number, snapshots_dir)
  snapshot_files <- snapshot_info[snapshot_info$type == "snapshot", ]

  if (nrow(snapshot_files) == 0) {
    stop("No snapshots found for issue #", issue_number)
  }

  has_status_pair <- any(snapshot_files$status == "baseline", na.rm = TRUE) &&
    any(snapshot_files$status == "fixed", na.rm = TRUE)

  if (has_status_pair) {
    baseline_files <- snapshot_files[snapshot_files$status == "baseline", ]
    fixed_files <- snapshot_files[snapshot_files$status == "fixed", ]

    before_file <- baseline_files$filename[which.max(baseline_files$mtime)]
    after_file <- fixed_files$filename[which.max(fixed_files$mtime)]
    message("Selected latest baseline and fixed snapshots for issue #", issue_number)
  } else {
    if (nrow(snapshot_files) == 1) {
      stop("Only one snapshot found for issue #", issue_number,
           ". Need at least two snapshots to compare.")
    }

    message("No baseline/fixed status pair found; defaulting to earliest vs latest snapshot.")
    snapshot_files <- snapshot_files[order(snapshot_files$mtime), ]
    before_file <- snapshot_files$filename[1]
    after_file <- snapshot_files$filename[nrow(snapshot_files)]
  }

  before_pair <- find_snapshot_pair(before_file, snapshots_dir)
  after_pair <- find_snapshot_pair(after_file, snapshots_dir)

  if (is.null(before_pair) || is.null(after_pair)) {
    stop("Could not find matching snapshot/PSA pairs")
  }

  list(before = before_pair, after = after_pair)
}
