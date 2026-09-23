# No economic-model change is intended in the publication-artifact fix.
source(here::here("R/snapshot_utils.R"))
source(here::here("R/cache_paths.R"))
pairs <- select_snapshots_for_comparison(168)
before <- load_snapshot(pairs$before$snapshot)
after <- load_snapshot(pairs$after$snapshot)
before_psa <- load_snapshot(pairs$before$psa)
after_psa <- load_snapshot(pairs$after$psa)
checks <- list(
  baseline_file = pairs$before$snapshot,
  fixed_file = pairs$after$snapshot,
  baseline_time = as.character(before$metadata$timestamp),
  fixed_time = as.character(after$metadata$timestamp),
  baseline_commit = before$metadata$git_commit,
  fixed_commit = after$metadata$git_commit,
  n_sim = after_psa$n_sim,
  deterministic_max_cost_difference = max(abs(after$base_results$Cost - before$base_results$Cost)),
  deterministic_max_qaly_difference = max(abs(after$base_results$Effect - before$base_results$Effect)),
  psa_max_cost_difference = max(abs(as.matrix(after_psa$cost) - as.matrix(before_psa$cost))),
  psa_max_qaly_difference = max(abs(as.matrix(after_psa$effect) - as.matrix(before_psa$effect))),
  same_psa_input_fingerprint = identical(before_psa$fingerprint, after_psa$fingerprint),
  same_model_indices = identical(before_psa$model_idx, after_psa$model_idx),
  before_evpi = before$evppi$evpi,
  after_evpi = after$evppi$evpi,
  baseline_psa_regenerated = before$metadata$caches$psa_regenerated,
  fixed_psa_regenerated = after$metadata$caches$psa_regenerated)
stopifnot(checks$n_sim == 5000L, checks$same_psa_input_fingerprint,
  checks$same_model_indices, checks$deterministic_max_cost_difference <= 1e-10,
  checks$deterministic_max_qaly_difference <= 1e-10,
  checks$psa_max_cost_difference <= 1e-7, checks$psa_max_qaly_difference <= 1e-11,
  checks$before_evpi == 0, checks$after_evpi == 0,
  isTRUE(checks$fixed_psa_regenerated))
jsonlite::write_json(checks, here::here("validation/issue168_2026-09-22/snapshot_comparison.json"),
  pretty = TRUE, auto_unbox = TRUE, digits = 16)
print(checks)
cat("PASS: deterministic results unchanged; all 5,000 freshly regenerated PSA rows agree within numerical tolerance; EVPI remains zero.\n")
