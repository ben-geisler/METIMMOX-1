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
#        tests/test_interaction_coefficient_extraction.R
# ===============================================================================

source(here::here("R", "model_configs.R"))
source(here::here("R", "evppi_functions.R"))

# Coefficient names exactly as flexsurvreg emits them for this model. Verified
# against data/tidy/sampling_models_n5000_full.rds.
v_real_coef_names <- c(
  "shape", "rate", "Age", "sex1", "RxExperimental arm", "crp1", "tmb_braf1",
  "RxExperimental arm:crp1", "RxExperimental arm:tmb_braf1"
)

#' Named coefficient vector in flexsurvreg order with the given CRP and TMB/BRAF
#' interaction values (named list or vector with elements crp, tmb_braf).
make_coefs <- function(interaction_values) {
  v_coefs <- c(0.5, 0.02, -0.01, 0.1, 0.3, 0.2, 0.15,
             interaction_values[["crp"]], interaction_values[["tmb_braf"]])
  names(v_coefs) <- v_real_coef_names
  v_coefs
}

n_sim <- 4L
#' Legacy bootstrap sample k: OS and PFS fits whose interactions scale with k.
one_sample <- function(k) {
  v_coefs <- make_coefs(c(crp = 0.100 * k, tmb_braf = -0.200 * k))
  list(os = list(coefficients = v_coefs), pfs = list(coefficients = v_coefs))
}

sampling_models <- list(
  control  = list(samples = lapply(seq_len(n_sim), one_sample)),
  crp      = list(samples = lapply(seq_len(n_sim), one_sample)),
  tmb_braf = list(samples = lapply(seq_len(n_sim), one_sample))
)

df_extracted <- suppressWarnings(
  extract_interaction_coefficients(sampling_models, n_sim)
)

# ---------------------------------------------------------------------------
# 1. All four interaction columns are produced and fully populated
# ---------------------------------------------------------------------------
v_expected_cols <- c("b_crp_rx_os", "b_crp_rx_pfs",
                   "b_tmb_braf_rx_os", "b_tmb_braf_rx_pfs")
stopifnot(
  !is.null(df_extracted),
  all(v_expected_cols %in% names(df_extracted)),
  nrow(df_extracted) == n_sim
)
for (col in v_expected_cols) {
  stopifnot(!any(is.na(df_extracted[[col]])))
}

# ---------------------------------------------------------------------------
# 2. Each biomarker picks up ITS OWN coefficient, not the other's
# ---------------------------------------------------------------------------
stopifnot(
  isTRUE(all.equal(df_extracted$b_crp_rx_os, 0.100 * seq_len(n_sim))),
  isTRUE(all.equal(df_extracted$b_tmb_braf_rx_os, -0.200 * seq_len(n_sim)))
)

# ---------------------------------------------------------------------------
# 3. Term ORDER must not matter: the reversed labelling extracts identically
# ---------------------------------------------------------------------------
v_reversed <- v_real_coef_names
v_reversed[8] <- "crp1:RxExperimental arm"
v_reversed[9] <- "tmb_braf1:RxExperimental arm"
#' As one_sample(), but with the interaction terms labelled biomarker-first
#' ("crp1:RxExperimental arm").
rev_sample <- function(k) {
  v_coefs <- c(0.5, 0.02, -0.01, 0.1, 0.3, 0.2, 0.15, 0.100 * k, -0.200 * k)
  names(v_coefs) <- v_reversed
  list(os = list(coefficients = v_coefs), pfs = list(coefficients = v_coefs))
}
rev_models <- list(
  control  = list(samples = lapply(seq_len(n_sim), rev_sample)),
  crp      = list(samples = lapply(seq_len(n_sim), rev_sample)),
  tmb_braf = list(samples = lapply(seq_len(n_sim), rev_sample))
)
df_rev_extracted <- suppressWarnings(
  extract_interaction_coefficients(rev_models, n_sim)
)
stopifnot(
  isTRUE(all.equal(df_rev_extracted$b_crp_rx_os, df_extracted$b_crp_rx_os)),
  isTRUE(all.equal(df_rev_extracted$b_tmb_braf_rx_os, df_extracted$b_tmb_braf_rx_os))
)

# ---------------------------------------------------------------------------
# 4. A main effect alone must NOT be mistaken for an interaction
# ---------------------------------------------------------------------------
v_main_only <- c(0.5, 0.02, -0.01, 0.1, 0.3, 0.2, 0.15)
names(v_main_only) <- v_real_coef_names[1:7]
no_int <- list(
  control  = list(samples = list(list(os = list(coefficients = v_main_only),
                                      pfs = list(coefficients = v_main_only)))),
  crp      = list(samples = list(list(os = list(coefficients = v_main_only),
                                      pfs = list(coefficients = v_main_only)))),
  tmb_braf = list(samples = list(list(os = list(coefficients = v_main_only),
                                      pfs = list(coefficients = v_main_only))))
)
df_no_int_result <- suppressWarnings(extract_interaction_coefficients(no_int, 1L))
stopifnot(is.null(df_no_int_result) || all(is.na(df_no_int_result$b_crp_rx_os)))

# ---------------------------------------------------------------------------
# 5. Issue #156 layout: one joint component with coefficient-draw matrices
#    gives the same columns as the legacy per-biomarker sample lists
# ---------------------------------------------------------------------------
m_draw <- do.call(rbind, lapply(seq_len(n_sim), function(k) {
  make_coefs(c(crp = 0.100 * k, tmb_braf = -0.200 * k))
}))
joint_models <- list(
  joint = list(method = "mvn_v1", draws = list(os = m_draw, pfs = m_draw)),
  biomarkers = c("crp", "tmb_braf")
)
df_joint_extracted <- suppressWarnings(extract_interaction_coefficients(joint_models, n_sim))
stopifnot(
  setequal(names(df_joint_extracted), v_expected_cols),
  isTRUE(all.equal(df_joint_extracted$b_crp_rx_os, df_extracted$b_crp_rx_os)),
  isTRUE(all.equal(df_joint_extracted$b_tmb_braf_rx_pfs, df_extracted$b_tmb_braf_rx_pfs))
)

cat("PASS: treatment-by-biomarker interaction coefficients extract correctly.\n")
