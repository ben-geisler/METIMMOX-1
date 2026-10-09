# ===============================================================================
# RUN ALL REGRESSION TESTS (issue #165)
# ===============================================================================
# Runs every tests/test_*.R in its own Rscript process, records exit status and
# run time, and prints a summary table. A test passes when its process exits
# with status 0.
#
# Caches: every test runs against a temporary COPY of data/tidy/ (option
# metimmox.cache_dir set through a temporary R_PROFILE_USER), so no test can
# write to or delete the live caches. The md5 of every live cache file is
# compared before and after the run; any change fails the runner.
#
# Trial data: tests listed in REQUIRES_TRIAL_DATA need the confidential
# data/tidy/METIMMOX.rds and are reported as SKIPPED when it is absent. Other
# tests that read trial data or caches skip those sections themselves.
#
# Known failures: KNOWN_FAILURES lists tests documented as failing in the test
# protocol (docs/validation/README.md of the private repository). They are
# reported as XFAIL and do not fail the run; a known failure that passes is
# reported as XPASS and does fail the run, so the list cannot go stale. The
# list is empty: test_psa_basecase_alignment.R left it
# when #180 retired its five-MCSE criterion, test_sampling_failure_fallback.R
# when #186 rewrote it for the #159 sampler.
#
# Usage (from the repository root):
#   "C:\Program Files\R\R-4.3.2\bin\x64\Rscript.exe" tests/run_all.R
#   ... tests/run_all.R --only test_pfs_endpoint,test_psa_summary
#   ... tests/run_all.R --no-data      # skip REQUIRES_TRIAL_DATA tests
#
# Exit status 0 when every test passes, skips or fails as expected; 1 otherwise.
# The summary is also written to data/output/test_run_all.csv (ignored by git).
# ===============================================================================

KNOWN_FAILURES <- character(0)

# Tests that stop (rather than skip) without data/tidy/METIMMOX.rds.
REQUIRES_TRIAL_DATA <- c(
  "test_survival_ordering.R",
  "test_enriched_survival_ordering.R"
)

#' Parse the runner's command-line options
#'
#' @param args Character vector of trailing command-line arguments.
#' @return A list with `only` (test file names or NULL) and `no_data` (logical).
parse_runner_args <- function(args) {
  only <- NULL
  only_pos <- match("--only", args)
  if (!is.na(only_pos)) {
    if (only_pos == length(args)) stop("--only needs a comma-separated list of tests.")
    only <- trimws(strsplit(args[only_pos + 1], ",", fixed = TRUE)[[1]])
    only <- ifelse(grepl("\\.R$", only), only, paste0(only, ".R"))
  }
  list(only = only, no_data = "--no-data" %in% args)
}

#' Copy the cache directory to a temporary location
#'
#' @param v_files Paths of the live cache files to copy.
#' @param directory Destination directory (created if needed).
#' @return The destination directory, invisibly.
copy_caches <- function(v_files, directory) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(v_files, directory, overwrite = TRUE, copy.date = TRUE)
  if (!all(ok)) stop("Could not copy cache files: ", paste(basename(v_files[!ok]), collapse = ", "))
  invisible(directory)
}

#' Run one test file in a fresh Rscript process
#'
#' @param test_file Path of the test script, relative to the repository root.
#' @param rscript Path of the Rscript executable.
#' @param profile Path of the temporary R profile that sets the cache directory.
#' @param log_file Path that receives the combined stdout and stderr.
#' @return A list with the exit `status` and the elapsed `seconds`.
run_one_test <- function(test_file, rscript, profile, log_file) {
  old_profile <- Sys.getenv("R_PROFILE_USER", unset = NA)
  Sys.setenv(R_PROFILE_USER = profile)
  on.exit(if (is.na(old_profile)) Sys.unsetenv("R_PROFILE_USER")
          else Sys.setenv(R_PROFILE_USER = old_profile), add = TRUE)
  t0 <- Sys.time()
  status <- suppressWarnings(system2(rscript, shQuote(test_file),
                                     stdout = log_file, stderr = log_file))
  list(status = as.integer(status), seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

#' Classify a test result
#'
#' @param status Integer exit status, or NA when the test was not run.
#' @param test_name File name of the test.
#' @return One of "PASS", "FAIL", "XFAIL", "XPASS", "SKIPPED".
classify_result <- function(status, test_name) {
  if (is.na(status)) return("SKIPPED")
  known <- test_name %in% KNOWN_FAILURES
  if (status == 0L) if (known) "XPASS" else "PASS" else if (known) "XFAIL" else "FAIL"
}

repo_root <- here::here()
setwd(repo_root)
opts <- parse_runner_args(commandArgs(trailingOnly = TRUE))
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")

v_tests <- sort(list.files("tests", pattern = "^test_.*\\.R$", full.names = TRUE))
if (!is.null(opts$only)) {
  unknown <- setdiff(opts$only, basename(v_tests))
  if (length(unknown) > 0) stop("Unknown test(s): ", paste(unknown, collapse = ", "))
  v_tests <- v_tests[basename(v_tests) %in% opts$only]
}
stale_known <- setdiff(KNOWN_FAILURES, basename(v_tests))
if (is.null(opts$only) && length(stale_known) > 0) {
  stop("KNOWN_FAILURES names tests that do not exist: ", paste(stale_known, collapse = ", "))
}

has_trial_data <- file.exists(file.path("data", "tidy", "METIMMOX.rds")) && !opts$no_data

# Live caches: everything in data/tidy except the trial data itself.
live_dir <- file.path(repo_root, "data", "tidy")
v_live_caches <- setdiff(list.files(live_dir, full.names = TRUE, pattern = "\\.(rds|RData)$"),
                         file.path(live_dir, "METIMMOX.rds"))
v_md5_before <- tools::md5sum(v_live_caches)

work_dir <- file.path(tempdir(), "metimmox_run_all")
cache_copy <- file.path(work_dir, "cache")
log_dir <- file.path(work_dir, "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
if (length(v_live_caches) > 0) copy_caches(v_live_caches, cache_copy) else dir.create(cache_copy, showWarnings = FALSE)

profile <- file.path(work_dir, "Rprofile_tests")
user_profile <- Sys.getenv("R_PROFILE_USER", unset = file.path(Sys.getenv("HOME"), ".Rprofile"))
writeLines(c(
  if (file.exists(user_profile)) sprintf("source(%s)", deparse(normalizePath(user_profile, winslash = "/"))),
  sprintf("options(metimmox.cache_dir = %s)", deparse(normalizePath(cache_copy, winslash = "/")))
), profile)

cat(sprintf("Running %d test(s); caches copied to %s\n", length(v_tests), cache_copy))
if (!has_trial_data) cat("Trial data not used: tests that require it are skipped.\n")

df_results <- data.frame(test = basename(v_tests), status = NA_integer_, result = NA_character_,
                         seconds = NA_real_, log = NA_character_, stringsAsFactors = FALSE)
for (i in seq_along(v_tests)) {
  test_name <- df_results$test[i]
  if (!has_trial_data && test_name %in% REQUIRES_TRIAL_DATA) {
    df_results$result[i] <- "SKIPPED"
    cat(sprintf("  %-46s SKIPPED (needs trial data)\n", test_name))
    next
  }
  log_file <- file.path(log_dir, sub("\\.R$", ".log", test_name))
  run <- run_one_test(v_tests[i], rscript, profile, log_file)
  df_results$status[i] <- run$status
  df_results$seconds[i] <- round(run$seconds, 1)
  df_results$log[i] <- log_file
  df_results$result[i] <- classify_result(run$status, test_name)
  cat(sprintf("  %-46s %-7s %6.1f s\n", test_name, df_results$result[i], run$seconds))
}

v_md5_after <- tools::md5sum(v_live_caches)
v_changed <- names(v_md5_before)[v_md5_before != v_md5_after | is.na(v_md5_after)]

cat("\nSummary:\n")
print(table(factor(df_results$result, levels = c("PASS", "FAIL", "XFAIL", "XPASS", "SKIPPED"))))
v_bad <- df_results$test[df_results$result %in% c("FAIL", "XPASS")]
for (test_name in v_bad) {
  log_file <- df_results$log[df_results$test == test_name]
  cat(sprintf("\n--- %s (last 20 lines of %s) ---\n", test_name, log_file))
  cat(tail(readLines(log_file, warn = FALSE), 20), sep = "\n")
}
if (length(v_changed) > 0) {
  cat("\nLIVE CACHE MODIFIED during the run: ", paste(basename(v_changed), collapse = ", "), "\n")
}

dir.create(file.path("data", "output"), recursive = TRUE, showWarnings = FALSE)
write.csv(df_results[, c("test", "status", "result", "seconds")],
          file.path("data", "output", "test_run_all.csv"), row.names = FALSE)

quit(save = "no", status = if (length(v_bad) > 0 || length(v_changed) > 0) 1L else 0L)
