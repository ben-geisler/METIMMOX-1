# Aggregate-only durable record; run after finish.ps1 completes.
source("R/cache_paths.R")
source("R/snapshot_utils.R")
utility_source_label <- "correct"
n_samples <- 5000
pair <- select_snapshots_for_comparison("181")
before <- load_snapshot(pair$before$snapshot)
after <- load_snapshot(pair$after$snapshot)
before_psa <- load_snapshot(pair$before$psa)
after_psa <- load_snapshot(pair$after$psa)
stopifnot(!identical(before$metadata$caches$psa_obj$md5, after$metadata$caches$psa_obj$md5),
          !identical(before_psa$fingerprint, after_psa$fingerprint),
          !isTRUE(all.equal(before$base_results, after$base_results)))
sampling <- readRDS(sampling_cache_path())
evppi <- new.env()
load(evppi_path(), envir = evppi)
report <- capture.output({
  cat("Issue #181 regenerated results\n")
  cat("Baseline:", pair$before$snapshot, "\nFixed:", pair$after$snapshot, "\n")
  cat(describe_psa_provenance(before, after), "\n")
  cat("Sampling fingerprint:", sampling$fingerprint, "\n")
  cat("Baseline PSA fingerprint:", before_psa$fingerprint, "\n")
  cat("Fixed PSA fingerprint:", after_psa$fingerprint, "\n")
  cat("Deterministic results:\n")
  print(after$base_results)
  cat("PSA summary:\n")
  print(after$psa_summary)
  cat("EVPI per patient at EUR 51,000:", evppi$evpi_manual, "\n")
  cat("EVPPI rows:", nrow(evppi$evppi_results), "\n")
  cat("Paired bootstrap attempted/successful/failed:\n")
  print(sampling$joint$joint_covariance[c("n_attempted", "n_successful", "n_failed")])
})
writeLines(report, "validation/issue181_2026-09-24/results.txt")
cat(paste(report, collapse = "\n"), "\n")
