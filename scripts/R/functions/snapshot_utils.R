# Snapshot Utility Functions
# Helper functions for managing result snapshots for bug fix impact analysis

#' Get Git Commit Hash
#'
#' Retrieves the current git commit hash (short version)
#'
#' @return Character string with short commit hash, or "unknown" if not in git repo
get_git_commit <- function() {
  tryCatch({
    commit <- system("git rev-parse --short HEAD", intern = TRUE)
    return(commit)
  }, error = function(e) {
    warning("Could not retrieve git commit hash")
    return("unknown")
  })
}

#' Get Package Versions
#'
#' Retrieves version numbers for key R packages used in the analysis
#'
#' @return Named list of package versions
get_package_versions <- function() {
  key_packages <- c(
    "dampack", "flexsurv", "survival", "dplyr", "ggplot2",
    "here", "mvtnorm", "Matrix", "darthtools", "gems"
  )

  versions <- list()
  for (pkg in key_packages) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      versions[[pkg]] <- as.character(packageVersion(pkg))
    } else {
      versions[[pkg]] <- "not installed"
    }
  }

  return(versions)
}

#' Collect Metadata
#'
#' Collects comprehensive metadata about the current analysis run
#'
#' @param issue_number GitHub issue number (character or numeric)
#' @return List containing all metadata
collect_metadata <- function(issue_number) {
  metadata <- list(
    issue_number = as.character(issue_number),
    git_commit = get_git_commit(),
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    r_version = R.version.string,
    packages = get_package_versions(),
    parameters = list(
      WTP = if(exists("WTP")) WTP else NA,
      time_horizon = if(exists("time_horizon")) time_horizon else NA,
      n_samples = if(exists("n_samples")) n_samples else NA,
      n_sim = if(exists("n_sim")) n_sim else NA,
      cl = if(exists("cl")) cl else NA,
      dr = if(exists("dr")) dr else NA,
      DSA_mult = if(exists("DSA_mult")) DSA_mult else NA,
      USE_BOTH_MODELS = if(exists("USE_BOTH_MODELS")) USE_BOTH_MODELS else NA,
      MODEL_STRUCTURE = if(exists("MODEL_STRUCTURE")) MODEL_STRUCTURE else NA
    )
  )

  return(metadata)
}

#' Generate Snapshot Filename
#'
#' Creates standardized filename for snapshot files
#'
#' @param type File type ("snapshot" or "psa")
#' @param issue_number GitHub issue number
#' @param timestamp Timestamp string (default: current time)
#' @param commit Git commit hash (default: current commit)
#' @return Character string with filename
generate_snapshot_filename <- function(type, issue_number,
                                      timestamp = format(Sys.time(), "%Y%m%d_%H%M%S"),
                                      commit = get_git_commit()) {
  filename <- paste0(type, "_", issue_number, "_", timestamp, "_", commit, ".rds")
  return(filename)
}

#' List Snapshots for Issue
#'
#' Lists all snapshot files for a given issue number
#'
#' @param issue_number GitHub issue number
#' @param snapshots_dir Directory containing snapshots (default: data/output/snapshots/)
#' @return Data frame with snapshot information (filename, type, timestamp, commit)
list_snapshots <- function(issue_number, snapshots_dir = here::here("data", "output", "snapshots")) {
  pattern <- paste0("^(snapshot|psa)_", issue_number, "_")
  files <- list.files(snapshots_dir, pattern = pattern, full.names = FALSE)

  if (length(files) == 0) {
    message("No snapshots found for issue #", issue_number)
    return(data.frame(filename = character(0), type = character(0),
                     timestamp = character(0), commit = character(0)))
  }

  # Parse filenames to extract information
  snapshot_info <- do.call(rbind, lapply(files, function(f) {
    parts <- strsplit(f, "_")[[1]]
    type <- parts[1]
    # timestamp is parts 3 and 4
    timestamp <- paste(parts[3], parts[4], sep = "_")
    # commit is part 5 without .rds
    commit <- gsub("\\.rds$", "", parts[5])

    data.frame(
      filename = f,
      type = type,
      timestamp = timestamp,
      commit = commit,
      stringsAsFactors = FALSE
    )
  }))

  # Sort by timestamp
  snapshot_info <- snapshot_info[order(snapshot_info$timestamp), ]

  return(snapshot_info)
}

#' Load Snapshot
#'
#' Loads a snapshot file and validates its structure
#'
#' @param filename Snapshot filename (without path)
#' @param snapshots_dir Directory containing snapshots
#' @return Snapshot object (list)
load_snapshot <- function(filename, snapshots_dir = here::here("data", "output", "snapshots")) {
  filepath <- file.path(snapshots_dir, filename)

  if (!file.exists(filepath)) {
    stop("Snapshot file not found: ", filepath)
  }

  snapshot <- readRDS(filepath)

  # Validate structure for snapshot files
  if (grepl("^snapshot_", filename)) {
    required_components <- c("metadata", "base_results", "icer_obj", "nmb_at_wtp", "psa_summary")
    missing <- setdiff(required_components, names(snapshot))
    if (length(missing) > 0) {
      warning("Snapshot missing components: ", paste(missing, collapse = ", "))
    }
  }

  # Validate structure for PSA files (accept single psa object or named list)
  if (grepl("^psa_", filename)) {
    if (!("psa" %in% class(snapshot)) && !is.list(snapshot)) {
      warning("PSA file does not contain a valid dampack psa object or model list")
    }
  }

  return(snapshot)
}

#' Check if a snapshot contains multi-model results
#'
#' @param snapshot A loaded snapshot object
#' @return Logical TRUE if snapshot contains models list
is_multimodel_snapshot <- function(snapshot) {
  "models" %in% names(snapshot)
}

#' Get model names from a snapshot (multi or single)
#'
#' @param snapshot A loaded snapshot object
#' @return Character vector of model names
get_snapshot_model_names <- function(snapshot) {
  if (is_multimodel_snapshot(snapshot)) {
    return(names(snapshot$models))
  }
  ms <- snapshot$metadata$parameters$MODEL_STRUCTURE
  if (is.null(ms) || is.na(ms)) return("joint")
  c("joint", "focused", "separate")[ms + 1]
}

#' Find Matching Snapshot Pair
#'
#' Finds matching snapshot and PSA files by timestamp and commit
#'
#' @param snapshot_file Snapshot filename
#' @param snapshots_dir Directory containing snapshots
#' @return List with snapshot and psa filenames, or NULL if no match
find_snapshot_pair <- function(snapshot_file, snapshots_dir = here::here("data", "output", "snapshots")) {
  # Extract timestamp and commit from snapshot filename
  parts <- strsplit(snapshot_file, "_")[[1]]
  issue <- parts[2]
  timestamp <- parts[3]
  time <- parts[4]
  commit <- gsub("\\.rds$", "", parts[5])

  # Construct expected PSA filename
  psa_file <- paste0("psa_", issue, "_", timestamp, "_", time, "_", commit, ".rds")

  # Check if PSA file exists
  if (file.exists(file.path(snapshots_dir, psa_file))) {
    return(list(snapshot = snapshot_file, psa = psa_file))
  } else {
    warning("No matching PSA file found for ", snapshot_file)
    return(NULL)
  }
}

#' Select Snapshots for Comparison
#'
#' Interactive function to select two snapshots for comparison
#'
#' @param issue_number GitHub issue number
#' @param snapshots_dir Directory containing snapshots
#' @return List with two snapshot pairs (before and after)
select_snapshots_for_comparison <- function(issue_number,
                                           snapshots_dir = here::here("data", "output", "snapshots")) {
  snapshot_info <- list_snapshots(issue_number, snapshots_dir)

  # Filter to only snapshot files (not psa)
  snapshot_files <- snapshot_info[snapshot_info$type == "snapshot", ]

  if (nrow(snapshot_files) == 0) {
    stop("No snapshots found for issue #", issue_number)
  }

  if (nrow(snapshot_files) == 1) {
    stop("Only one snapshot found for issue #", issue_number,
         ". Need at least two snapshots to compare.")
  }

  if (nrow(snapshot_files) == 2) {
    # Exactly two snapshots - use them
    message("Found 2 snapshots for issue #", issue_number)
    before_file <- snapshot_files$filename[1]
    after_file <- snapshot_files$filename[2]
  } else {
    # More than two snapshots - let user choose or default to first and last
    message("Found ", nrow(snapshot_files), " snapshots for issue #", issue_number)
    message("\nAvailable snapshots (chronological order):")
    for (i in 1:nrow(snapshot_files)) {
      message(sprintf("  %d. %s (commit: %s)", i, snapshot_files$timestamp[i],
                     snapshot_files$commit[i]))
    }

    message("\nDefaulting to earliest vs. latest (1 vs. ", nrow(snapshot_files), ")")
    message("To compare different snapshots, modify the compare_snapshots.R script.")

    before_file <- snapshot_files$filename[1]
    after_file <- snapshot_files$filename[nrow(snapshot_files)]
  }

  # Find matching PSA files
  before_pair <- find_snapshot_pair(before_file, snapshots_dir)
  after_pair <- find_snapshot_pair(after_file, snapshots_dir)

  if (is.null(before_pair) || is.null(after_pair)) {
    stop("Could not find matching snapshot pairs")
  }

  return(list(before = before_pair, after = after_pair))
}
