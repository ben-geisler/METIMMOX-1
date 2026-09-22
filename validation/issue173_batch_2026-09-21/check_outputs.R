# Run from the repository root after recompute.R and the fixed snapshot.
invisible(parse("validation/issue173_batch_2026-09-21/recompute.R"))
source("R/model_configs.R")
diag <- readRDS("data/tidy/pfs_os_violations_n5000.rds")
psa <- readRDS("data/tidy/psa_obj_correct.rds")
params <- readRDS("data/tidy/psa_params_correct.rds")
stopifnot(diag$n_draws == psa$n_sim, psa$n_sim == nrow(params),
          identical(diag$fingerprint_inputs$sampling_fingerprint,
                    psa$fingerprint_inputs$sampling_fingerprint))
pc <- diag$per_curve
curve_qalys <- function(curve) {
  rows <- pc[pc$curve == curve, ]
  stopifnot(nrow(rows) == psa$n_sim, !anyDuplicated(rows$draw),
            setequal(rows$draw, params$sim))
  rows$qaly_clamped[match(params$sim, rows$draw)]
}
reconstructed <- matrix(NA_real_, nrow(params), ncol(psa$effectiveness),
  dimnames = list(NULL, colnames(psa$effectiveness)))
control <- get_control_strategy()
reconstructed[, control] <- curve_qalys(control)
for (bm in get_biomarkers()) {
  prevalence <- params[[paste0("p_", bm)]]
  reconstructed[, bm] <- prevalence * curve_qalys(paste0(bm, "_pos")) +
    (1 - prevalence) * curve_qalys(paste0(bm, "_neg"))
}
error <- max(abs(reconstructed - as.matrix(psa$effectiveness)))
stopifnot(is.finite(error), error < 1e-11)
cat("PASS: all", psa$n_sim, "draws' diagnostic QALYs match the PSA; max difference", error, "\n")

before <- readRDS("data/output/snapshots/snapshot_173_baseline_1b50e8d.rds")
after <- readRDS("data/output/snapshots/snapshot_173_fixed_1b50e8d.rds")
sampling_hash <- before$metadata$caches$sampling$md5
psa_before <- before$metadata$caches$psa_obj$md5
psa_after <- after$metadata$caches$psa_obj$md5
stopifnot(is.character(sampling_hash), length(sampling_hash) == 1L,
          nzchar(sampling_hash), identical(sampling_hash, after$metadata$caches$sampling$md5),
          identical(sampling_hash, unname(tools::md5sum("data/tidy/sampling_models_n5000_full.rds"))),
          is.character(psa_before), length(psa_before) == 1L, nzchar(psa_before),
          !identical(psa_before, psa_after),
          identical(psa_after, unname(tools::md5sum("data/tidy/psa_obj_correct.rds"))))
cat("PASS: snapshots have changed economic results and identical survival sampling.\n")
