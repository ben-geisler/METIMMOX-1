# ===============================================================================
# TEST: interaction coefficients extract from real flexsurv coefficient names
# ===============================================================================
# extract_interaction_coefficients() previously matched "^crp.*:Rx" (and before
# that "crp.*:Rx"). Neither ever matched, because the economic formula
#   ~ Age + sex + Rx + crp * Rx + tmb_braf * Rx
# puts Rx ahead of the biomarker in the term expansion, so R labels the
# interaction "Rx:crp". With flexsurvreg appending factor levels, the actual
# coefficient is "RxExperimental arm:crp1".
#
# The mismatch was silent: every b_*_rx_* row was dropped from the EVPPI and
# scenario-EVPPI tables behind a warning, so the value of information on the
# treatment-by-biomarker interaction was absent from published results.
#
# Run: "C:\Program Files\R\R-4.3.2\bin\x64\Rscript.exe" \
#        scripts/R/tests/test_interaction_coefficient_extraction.R
# ===============================================================================

source(here::here("scripts", "R", "functions", "model_configs.R"))
source(here::here("scripts", "R", "functions", "evppi_functions.R"))

# Coefficient names exactly as flexsurvreg emits them for this model. Verified
# against data/tidy/sampling_models_n5000_full.rds.
real_coef_names <- c(
  "shape", "rate", "Age", "sex1", "RxExperimental arm", "crp1", "tmb_braf1",
  "RxExperimental arm:crp1", "RxExperimental arm:tmb_braf1"
)

make_coefs <- function(interaction_values) {
  coefs <- c(0.5, 0.02, -0.01, 0.1, 0.3, 0.2, 0.15,
             interaction_values[["crp"]], interaction_values[["tmb_braf"]])
  names(coefs) <- real_coef_names
  coefs
}

n_sim <- 4L
one_sample <- function(k) {
  coefs <- make_coefs(c(crp = 0.100 * k, tmb_braf = -0.200 * k))
  list(os = list(coefficients = coefs), pfs = list(coefficients = coefs))
}

sampling_models <- list(
  control  = list(samples = lapply(seq_len(n_sim), one_sample)),
  crp      = list(samples = lapply(seq_len(n_sim), one_sample)),
  tmb_braf = list(samples = lapply(seq_len(n_sim), one_sample))
)

extracted <- suppressWarnings(
  extract_interaction_coefficients(sampling_models, n_sim)
)

# ---------------------------------------------------------------------------
# 1. All four interaction columns are produced and fully populated
# ---------------------------------------------------------------------------
expected_cols <- c("b_crp_rx_os", "b_crp_rx_pfs",
                   "b_tmb_braf_rx_os", "b_tmb_braf_rx_pfs")
stopifnot(
  !is.null(extracted),
  all(expected_cols %in% names(extracted)),
  nrow(extracted) == n_sim
)
for (col in expected_cols) {
  stopifnot(!any(is.na(extracted[[col]])))
}

# ---------------------------------------------------------------------------
# 2. Each biomarker picks up ITS OWN coefficient, not the other's
# ---------------------------------------------------------------------------
stopifnot(
  isTRUE(all.equal(extracted$b_crp_rx_os, 0.100 * seq_len(n_sim))),
  isTRUE(all.equal(extracted$b_tmb_braf_rx_os, -0.200 * seq_len(n_sim)))
)

# ---------------------------------------------------------------------------
# 3. Term ORDER must not matter: the reversed labelling extracts identically
# ---------------------------------------------------------------------------
reversed <- real_coef_names
reversed[8] <- "crp1:RxExperimental arm"
reversed[9] <- "tmb_braf1:RxExperimental arm"
rev_sample <- function(k) {
  coefs <- c(0.5, 0.02, -0.01, 0.1, 0.3, 0.2, 0.15, 0.100 * k, -0.200 * k)
  names(coefs) <- reversed
  list(os = list(coefficients = coefs), pfs = list(coefficients = coefs))
}
rev_models <- list(
  control  = list(samples = lapply(seq_len(n_sim), rev_sample)),
  crp      = list(samples = lapply(seq_len(n_sim), rev_sample)),
  tmb_braf = list(samples = lapply(seq_len(n_sim), rev_sample))
)
rev_extracted <- suppressWarnings(
  extract_interaction_coefficients(rev_models, n_sim)
)
stopifnot(
  isTRUE(all.equal(rev_extracted$b_crp_rx_os, extracted$b_crp_rx_os)),
  isTRUE(all.equal(rev_extracted$b_tmb_braf_rx_os, extracted$b_tmb_braf_rx_os))
)

# ---------------------------------------------------------------------------
# 4. A main effect alone must NOT be mistaken for an interaction
# ---------------------------------------------------------------------------
main_only <- c(0.5, 0.02, -0.01, 0.1, 0.3, 0.2, 0.15)
names(main_only) <- real_coef_names[1:7]
no_int <- list(
  control  = list(samples = list(list(os = list(coefficients = main_only),
                                      pfs = list(coefficients = main_only)))),
  crp      = list(samples = list(list(os = list(coefficients = main_only),
                                      pfs = list(coefficients = main_only)))),
  tmb_braf = list(samples = list(list(os = list(coefficients = main_only),
                                      pfs = list(coefficients = main_only))))
)
no_int_result <- suppressWarnings(extract_interaction_coefficients(no_int, 1L))
stopifnot(is.null(no_int_result) || all(is.na(no_int_result$b_crp_rx_os)))

# ---------------------------------------------------------------------------
# 5. Issue #156 layout: one joint component with coefficient-draw matrices
#    gives the same columns as the legacy per-biomarker sample lists
# ---------------------------------------------------------------------------
draw_matrix <- do.call(rbind, lapply(seq_len(n_sim), function(k) {
  make_coefs(c(crp = 0.100 * k, tmb_braf = -0.200 * k))
}))
joint_models <- list(
  joint = list(method = "mvn_v1", draws = list(os = draw_matrix, pfs = draw_matrix)),
  biomarkers = c("crp", "tmb_braf")
)
joint_extracted <- suppressWarnings(extract_interaction_coefficients(joint_models, n_sim))
stopifnot(
  setequal(names(joint_extracted), expected_cols),
  isTRUE(all.equal(joint_extracted$b_crp_rx_os, extracted$b_crp_rx_os)),
  isTRUE(all.equal(joint_extracted$b_tmb_braf_rx_pfs, extracted$b_tmb_braf_rx_pfs))
)

cat("PASS: treatment-by-biomarker interaction coefficients extract correctly.\n")
