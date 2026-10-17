# ===============================================================================
# TEST: survival coefficients and cell effects as EVPPI inputs (issue #199)
# ===============================================================================
# extract_survival_coefficients() returns every coefficient of the joint OS and
# PFS models for the cached draw behind each PSA row (`model_idx`), and
# add_cell_effects() derives the treatment effect in each joint biomarker cell.
# A synthetic sampling component with the flexsurvreg coefficient names of the
# economic model (no trial data or caches) checks that:
#
#   1. rows follow model_idx, not the PSA row number, for both endpoints;
#   2. coefficients are classified by structure (baseline, Rx main effect,
#      prognostic, interaction) whatever the order of the interaction terms;
#   3. the cell effects equal Rx plus the interactions that apply;
#   4. the interactions agree with extract_interaction_coefficients();
#   5. build_survival_evppi_inputs() appends the columns and stops when the
#      PSA's own interaction columns are not aligned with model_idx;
#   6. the EVPPI groups have the documented members;
#   7. invalid indices and unrecognised coefficients stop.
#
# Run: "C:\Program Files\R\R-4.3.2\bin\x64\Rscript.exe" \
#        tests/test_survival_coefficient_extraction.R
# ===============================================================================

source(here::here("R", "model_configs.R"))
source(here::here("R", "evppi_functions.R"))
source(here::here("R", "prediction_functions.R"))
source(here::here("R", "voi_survival_parameters.R"))

#' Stop unless `expr` raises an error whose message matches `pattern`.
expect_error_matching <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  if (!inherits(error, "error") || !grepl(pattern, conditionMessage(error))) {
    stop("Expected an error matching '", pattern, "', got: ",
         if (inherits(error, "error")) conditionMessage(error) else "no error")
  }
  invisible(TRUE)
}

v_coef_names <- c("shape", "rate", "Age", "sex1", "RxExperimental arm", "crp1", "tmb_braf1",
                  "RxExperimental arm:crp1", "RxExperimental arm:tmb_braf1")
n_cached <- 20L
#' Draw matrix whose entry (k, j) encodes draw k, coefficient j and the endpoint.
make_draws <- function(offset, v_names = v_coef_names) {
  m <- outer(seq_len(n_cached), seq_along(v_names), function(k, j) offset + k + j / 100)
  colnames(m) <- v_names
  m
}
m_os <- make_draws(0)
m_pfs <- make_draws(1000)
sampling_models <- list(joint = list(method = "mvn_joint_v2", draws = list(os = m_os, pfs = m_pfs)),
                        biomarkers = c("crp", "tmb_braf"), fingerprint = "synthetic")

# ---------------------------------------------------------------------------
# 1-2. Alignment by model_idx and structural classification
# ---------------------------------------------------------------------------
v_model_idx <- c(5L, 5L, 17L, 1L, 20L, 3L)
df_coefs <- extract_survival_coefficients(sampling_models, n_cached, v_model_idx)
df_map <- attr(df_coefs, "coefficient_map")
stopifnot(nrow(df_coefs) == length(v_model_idx), ncol(df_coefs) == 18,
          identical(attr(df_coefs, "model_idx"), v_model_idx),
          all(df_coefs$coef_os_shape == m_os[v_model_idx, "shape"]),
          all(df_coefs$coef_pfs_rx_x_tmb_braf == m_pfs[v_model_idx, "RxExperimental arm:tmb_braf1"]),
          all(df_coefs$coef_os_Age == m_os[v_model_idx, "Age"]),
          identical(df_map$role[df_map$outcome == "os"],
                    c("baseline", "baseline", "baseline", "baseline", "rx_main",
                      "prognostic", "prognostic", "interaction", "interaction")),
          identical(df_map$biomarker[df_map$role == "interaction"],
                    c("crp", "tmb_braf", "crp", "tmb_braf")))

# Reversed interaction labels ("crp1:RxExperimental arm") classify identically.
v_reversed <- v_coef_names
v_reversed[8:9] <- c("crp1:RxExperimental arm", "tmb_braf1:RxExperimental arm")
reversed_models <- list(joint = list(draws = list(os = make_draws(0, v_reversed),
                                                  pfs = make_draws(1000, v_reversed))))
df_rev <- extract_survival_coefficients(reversed_models, n_cached, v_model_idx)
stopifnot(identical(names(df_rev), names(df_coefs)),
          isTRUE(all.equal(df_rev, df_coefs, check.attributes = FALSE)))

# The legacy per-sample layout gives the same columns.
legacy_models <- list(crp = list(samples = lapply(seq_len(n_cached), function(k) {
  list(os = list(coefficients = m_os[k, ]), pfs = list(coefficients = m_pfs[k, ]))
})))
df_legacy <- extract_survival_coefficients(legacy_models, n_cached, v_model_idx)
stopifnot(isTRUE(all.equal(df_legacy, df_coefs, check.attributes = FALSE)))
cat("PASS: coefficients follow model_idx and are classified by structure\n")

# ---------------------------------------------------------------------------
# 3. Cell effects are the sums of Rx and the applicable interactions
# ---------------------------------------------------------------------------
df_cells <- add_cell_effects(df_coefs)
for (outcome in c("os", "pfs")) {
  v_rx <- df_coefs[[paste0("coef_", outcome, "_rx")]]
  v_crp <- df_coefs[[paste0("coef_", outcome, "_rx_x_crp")]]
  v_tmb <- df_coefs[[paste0("coef_", outcome, "_rx_x_tmb_braf")]]
  stopifnot(
    isTRUE(all.equal(df_cells[[paste0("cell_effect_00_", outcome)]], v_rx)),
    isTRUE(all.equal(df_cells[[paste0("cell_effect_10_", outcome)]], v_rx + v_crp)),
    isTRUE(all.equal(df_cells[[paste0("cell_effect_01_", outcome)]], v_rx + v_tmb)),
    isTRUE(all.equal(df_cells[[paste0("cell_effect_11_", outcome)]], v_rx + v_crp + v_tmb)))
}
cat("PASS: cell effects equal Rx plus the interactions of the positive biomarkers\n")

# ---------------------------------------------------------------------------
# 4. Agreement with extract_interaction_coefficients() (draw order)
# ---------------------------------------------------------------------------
v_extraction_log <- suppressWarnings(utils::capture.output(
  df_int <- extract_interaction_coefficients(sampling_models, n_cached)))
df_all_draws <- extract_survival_coefficients(sampling_models, n_cached)
stopifnot(all(df_int$b_crp_rx_os == df_all_draws$coef_os_rx_x_crp),
          all(df_int$b_tmb_braf_rx_pfs == df_all_draws$coef_pfs_rx_x_tmb_braf))
cat("PASS: interactions agree with extract_interaction_coefficients()\n")

# ---------------------------------------------------------------------------
# 5. PSA inputs: appended by model_idx, misalignment detected
# ---------------------------------------------------------------------------
df_psa_params <- data.frame(sim = seq_along(v_model_idx), u_np = 0.7, model_idx = v_model_idx)
df_psa_params$b_crp_rx_os <- m_os[v_model_idx, "RxExperimental arm:crp1"]
df_inputs <- build_survival_evppi_inputs(df_psa_params, sampling_models)
stopifnot(nrow(df_inputs) == nrow(df_psa_params), ncol(df_inputs) == ncol(df_psa_params) + 26,
          all(df_inputs$coef_os_rx == m_os[v_model_idx, "RxExperimental arm"]))
df_misaligned <- df_psa_params
df_misaligned$b_crp_rx_os <- m_os[seq_along(v_model_idx), "RxExperimental arm:crp1"]
expect_error_matching(build_survival_evppi_inputs(df_misaligned, sampling_models),
                      "not aligned by model_idx")
expect_error_matching(build_survival_evppi_inputs(df_psa_params[, c("sim", "u_np")], sampling_models),
                      "no model_idx")
cat("PASS: PSA inputs are aligned by model_idx and misalignment stops\n")

# ---------------------------------------------------------------------------
# 6. Groups
# ---------------------------------------------------------------------------
groups <- attr(df_inputs, "survival_groups")
stopifnot(
  identical(names(groups), c("survival_all", "cell_effect_00", "cell_effect_10", "cell_effect_01",
                             "cell_effect_11", "rx_main", "prognostic", "baseline_survival")),
  length(groups$survival_all) == 18, !any(grepl("^cell_effect", groups$survival_all)),
  identical(groups$cell_effect_10, c("cell_effect_10_os", "cell_effect_10_pfs")),
  identical(groups$rx_main, c("coef_os_rx", "coef_pfs_rx")),
  setequal(groups$prognostic, c("coef_os_crp", "coef_os_tmb_braf", "coef_pfs_crp", "coef_pfs_tmb_braf")),
  setequal(groups$baseline_survival, paste0("coef_", rep(c("os", "pfs"), each = 4), "_",
                                             c("shape", "rate", "Age", "sex1"))),
  identical(evppi_method_label("gam", 2), "gam_tensor"),
  identical(evppi_method_label("gam", 18), "gam_additive"),
  identical(evppi_method_label("bart", 18), "bart"))
cat("PASS: survival EVPPI groups have the documented members\n")

# ---------------------------------------------------------------------------
# 7. Invalid input stops
# ---------------------------------------------------------------------------
expect_error_matching(extract_survival_coefficients(sampling_models, n_cached, c(1L, 21L)),
                      "model_idx must be integers")
expect_error_matching(extract_survival_coefficients(sampling_models, n_cached + 5L),
                      "hold 20 rows")
v_odd <- v_coef_names
v_odd[8] <- "Age:sex1"
expect_error_matching(extract_survival_coefficients(
  list(joint = list(draws = list(os = make_draws(0, v_odd), pfs = m_pfs))), n_cached),
  "Unrecognised survival coefficient")
cat("PASS: invalid indices and unrecognised coefficients stop\n")

cat("PASS: test_survival_coefficient_extraction.R\n")
