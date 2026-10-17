# ===============================================================================
# TEST: PSA on the draws of an evidence state (paper 3, E2; issue #203)
# ===============================================================================
# Synthetic (gamma fits on simulated patients, an mvn_joint_v2 component of 8
# draws, 12 prediction patients in all four joint biomarker cells; no trial
# data, no live cache):
#   1. A state without evidence (posterior = prior, same seed) regenerates the
#      E0 draws exactly, and its PSA reproduces the E0 PSA exactly: cost,
#      QALYs, life-years, cells and model_idx; every economic draw is shared.
#   2. A state with evidence (shifted mean, smaller covariance) keeps the E0
#      economic draws but changes the survival outcomes; the PSA is a bound
#      psa_pair_v2 stamped with the state's PSA fingerprint and label.
#   3. with_sampling_models() swaps the survival draws for the evaluation only
#      and restores them afterwards, also on error, and removes a binding it
#      created; check_common_economic_draws() stops on a differing draw.
#   4. Cache: load_or_run_evidence_psa() runs and saves, then reads the same
#      PSA back; a changed evidence state makes the loader name
#      upstream$evidence_states and reruns; the E0 PSA path is refused.
#   5. evidence_state_comparison(): the resample check equals the state's PSA
#      when the state is E0 itself.
#
# Run: "C:\Program Files\R\R-4.3.2\bin\x64\Rscript.exe" tests/test_voi_evidence_psa.R
# ===============================================================================

suppressPackageStartupMessages({
  library(dampack)
  library(survival)
  library(flexsurv)
})
for (f in c("model_configs", "prediction_functions", "model_fun", "calculate_outcomes", "parameter_distributions",
            "psa_functions", "cache_paths", "cache_provenance", "evppi_functions", "joint_survival_sampling",
            "pfs_os_violation_diagnostics", "design_information", "voi_price_axis", "voi_policies", "voi_cache",
            "evidence_update", "voi_evidence_psa")) {
  source(here::here("R", paste0(f, ".R")))
}

#' Stop unless `expr` raises an error whose message matches `pattern`.
expect_error_matching <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  if (!inherits(error, "error") || !grepl(pattern, conditionMessage(error))) {
    stop("Expected an error matching '", pattern, "', got: ",
         if (inherits(error, "error")) conditionMessage(error) else "no error")
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Synthetic fits, joint component, population and parameters
# ---------------------------------------------------------------------------
set.seed(31)
n_fit <- 300
df_fit <- data.frame(Age = stats::rnorm(n_fit, 62, 9), sex = factor(stats::rbinom(n_fit, 1, 0.45)),
                     Rx = factor(sample(c("control", "experimental"), n_fit, TRUE),
                                 levels = c("control", "experimental")),
                     crp = stats::rbinom(n_fit, 1, 0.35), tmb_braf = stats::rbinom(n_fit, 1, 0.45))
v_eta <- log(0.02) + 0.01 * (df_fit$Age - 62) - 0.2 * (df_fit$Rx == "experimental") -
  0.3 * df_fit$crp * (df_fit$Rx == "experimental")
v_t <- stats::rgamma(n_fit, 1.6, exp(v_eta)); v_c <- stats::runif(n_fit, 30, 150)
df_fit$OSwk <- pmin(v_t, v_c); df_fit$Death <- as.numeric(v_t <= v_c)
v_t2 <- pmin(v_t, stats::rgamma(n_fit, 1.2, exp(v_eta + 0.8)))
df_fit$PFSwk <- pmin(v_t2, v_c); df_fit$Progression <- as.numeric(v_t2 <= v_c)
rhs <- "Age + sex + Rx + crp * Rx + tmb_braf * Rx"
fit_os <- flexsurvreg(stats::as.formula(paste("Surv(OSwk, Death) ~", rhs)), data = df_fit, dist = "gamma")
fit_pfs <- flexsurvreg(stats::as.formula(paste("Surv(PFSwk, Progression) ~", rhs)), data = df_fit, dist = "gamma")
m_ro <- covariance_root(fit_os$cov); m_rp <- covariance_root(fit_pfs$cov)
m_joint <- rbind(cbind(fit_os$cov, 0.4 * m_ro %*% m_rp), cbind(0.4 * m_rp %*% m_ro, fit_pfs$cov))
n_draw <- 8L
seed_draws <- 5L
component <- list(method = "mvn_joint_v2", joint_covariance = list(covariance = m_joint),
                  draws = draw_joint_survival_coefficients(fit_os, fit_pfs, m_joint, n_draw, seed_draws),
                  original_os = fit_os, original_pfs = fit_pfs, n_samples = n_draw, n_failed = 0L,
                  dist_os = "gamma", dist_pfs = "gamma", seed = seed_draws)
dimnames(component$joint_covariance$covariance) <- list(joint_coefficient_names(component),
                                                        joint_coefficient_names(component))

horizon <- 104
time_points <- seq(0, horizon)
df_population <- data.frame(
  ID = 1:12, Age = c(45, 52, 60, 71, 48, 66, 74, 55, 69, 50, 63, 77),
  sex = factor(c(0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0)),
  Rx = factor(rep(c("control", "experimental"), 6), levels = c("control", "experimental")),
  crp = c(0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1), tmb_braf = c(0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1),
  OSwk = 100, Death = 1, PFSwk = 50, Progression = 1)
df_population <- economic_prediction_population(df_population)
joint_cells <- joint_biomarker_population(df_population)
params <- list(time_horizon = horizon, cl = 1 / 52, dr_costs = 0.04, dr_effects = 0.04,
               u_np = 0.73, u_p = 0.59, u_decrement = 0.14, c_drug_nivo = 13923, c_drug_FLOX = 100,
               c_test_CT = 150, c_test_blood = 20, c_test_CRP = 16, c_test_NGS = 2518,
               c_other_visit = 213.46, c_other_baseline = 213.46, c_other_follow = 213.46,
               c_other_pp = 0, c_other_last = 5000)
params <- c(params, build_treatment_schedules(horizon + 1))
params$joint_counts <- joint_cells$counts
v_mass <- joint_cells$counts / sum(joint_cells$counts)
for (k in seq_along(v_mass)) params[[paste0("p_joint_", rownames(joint_cells$cells)[k])]] <- v_mass[[k]]
params$p_crp <- sum(v_mass[c(3, 4)])
params$p_tmb_braf <- sum(v_mass[c(2, 4)])
preds <- generate_population_averaged_predictions(list(os = fit_os, pfs = fit_pfs), NULL, df_population,
                                                  time_points, quiet = TRUE,
                                                  weights = model_population_weights(params, df_population))
params <- set_population_predictions(params, preds)
# Later treatment sequences (issue #208).
source(here::here("R", "later_sequences.R"))
later <- later_sequence_inputs(df_population, data.frame(
  ID = df_population$ID, Seq3 = rep(c(1, 1, 0), 4), Seq3wk = rep(c(40, 60, NA), 4),
  Seq4 = rep(c(0, 1, 0), 4), Seq4wk = rep(c(NA, 80, NA), 4)))
params$later_sequences <- later$later_sequences
params[names(later$shares)] <- later$shares
config <- configure_parameter_distributions(params)
v_strategies <- get_strategies()
psa_seed <- 7L

# E0: the base sampling cache and its PSA.
n_samples <- n_draw
sampling_models <- list(joint = component, fingerprint = "synthetic-sampling-e0")
sampling_e0 <- sampling_models
#' Build and bind the PSA on the current `sampling_models` (quietly).
build_e0 <- function() {
  build <- NULL
  capture.output(build <- suppressWarnings(build_psa_obj(params, config$distributions, v_strategies, horizon,
                                                         1 / 52, n_draw, seed = psa_seed, extended_outcomes = TRUE)))
  build$psa_obj$fingerprint <- "synthetic-psa-e0"
  bind_psa_pair(build$psa_obj, build$psa_params)
}
pair_e0 <- build_e0()
#' Run the evidence-state PSA quietly with the synthetic inputs.
run_state_psa <- function(state, base_params = pair_e0$psa_params) {
  out <- NULL
  capture.output(out <- suppressWarnings(run_evidence_psa(state, params, config$distributions, v_strategies,
                                                          horizon, 1 / 52, n_draw, psa_seed,
                                                          base_psa_params = base_params)))
  out
}
#' Largest absolute difference between the outcomes of two PSA objects.
psa_outcome_gap <- function(a, b) {
  max(abs(as.matrix(a$cost) - as.matrix(b$cost)), abs(as.matrix(a$effectiveness) - as.matrix(b$effectiveness)),
      abs(as.matrix(a$ly_disc) - as.matrix(b$ly_disc)), abs(a$cells$values - b$cells$values))
}

# ---------------------------------------------------------------------------
# 1. A state without evidence reproduces E0 exactly
# ---------------------------------------------------------------------------
prior <- joint_prior_moments(component)
joint_none <- regenerate_evidence_component(component, prior$mean, prior$cov, n_draw, seed_draws)
state_none <- list(joint = joint_none, fingerprint = "synthetic-state-none", evidence_state = "E2none",
                   base_sampling_fingerprint = sampling_e0$fingerprint)
stopifnot(max(abs(joint_draw_matrix(joint_none) - joint_draw_matrix(component))) < 1e-12)
psa_none <- run_state_psa(state_none)
stopifnot(psa_outcome_gap(psa_none$psa_obj, pair_e0$psa_obj) < 1e-8,
          identical(psa_none$psa_obj$model_idx, pair_e0$psa_obj$model_idx),
          psa_none$diagnostics$draws_shared_with_e0 == n_draw,
          identical(psa_none$psa_obj$evidence_state, "E2none"),
          identical(sampling_models, sampling_e0))
cat(sprintf("1. no evidence: state draws equal E0's; its PSA reproduces E0 (max |diff| %.1e), all %d economic draws shared\n",
            psa_outcome_gap(psa_none$psa_obj, pair_e0$psa_obj), n_draw))

# ---------------------------------------------------------------------------
# 2. A state with evidence: same economic draws, other survival outcomes
# ---------------------------------------------------------------------------
v_shift <- prior$mean
v_int <- "pfs::Rxexperimental:crp"
stopifnot(v_int %in% names(v_shift))
v_shift[v_int] <- v_shift[v_int] - 0.5
joint_shift <- regenerate_evidence_component(component, v_shift, prior$cov * 0.5, n_draw, seed_draws)
state_shift <- list(joint = joint_shift, fingerprint = "synthetic-state-shift", evidence_state = "E2shift",
                    base_sampling_fingerprint = sampling_e0$fingerprint)
psa_shift <- run_state_psa(state_shift)
v_econ <- names(config$distributions)
stopifnot(psa_shift$diagnostics$draws_shared_with_e0 == psa_shift$psa_obj$n_sim,
          identical(psa_shift$psa_params[v_econ], pair_e0$psa_params[match(psa_shift$psa_params$sim,
                                                                           pair_e0$psa_params$sim), v_econ]),
          psa_outcome_gap(psa_shift$psa_obj, pair_e0$psa_obj) > 1e-6,
          identical(psa_shift$psa_obj$fingerprint,
                    psa_cache_fingerprint("synthetic-state-shift", params, config$distributions, v_strategies,
                                          n_draw, psa_seed, horizon, 1 / 52)$fingerprint),
          identical(psa_shift$psa_obj$pair_metadata$schema, PSA_PAIR_SCHEMA),
          identical(sampling_models, sampling_e0))
validate_psa_pair(psa_shift$psa_obj, psa_shift$psa_params)
stopifnot(mean(joint_draw_matrix(joint_shift)[, v_int]) < mean(joint_draw_matrix(component)[, v_int]))
cat("2. evidence: E0 economic draws kept, survival outcomes changed, bound pair with the state's fingerprint\n")

# ---------------------------------------------------------------------------
# 3. with_sampling_models() and the shared-draw check
# ---------------------------------------------------------------------------
seen <- with_sampling_models(state_shift, get("sampling_models", envir = globalenv())$evidence_state,
                             envir = globalenv())
stopifnot(identical(seen, "E2shift"), identical(sampling_models, sampling_e0))
expect_error_matching(with_sampling_models(state_shift, stop("boom"), envir = globalenv()), "boom")
stopifnot(identical(sampling_models, sampling_e0))
env_empty <- new.env()
stopifnot(identical(with_sampling_models(state_shift, get("sampling_models", envir = env_empty)$evidence_state,
                                         envir = env_empty), "E2shift"),
          !exists("sampling_models", envir = env_empty, inherits = FALSE))
expect_error_matching(with_sampling_models(list(fingerprint = "x"), 1), "no joint sampling component")
df_tampered <- pair_e0$psa_params
df_tampered[[v_econ[1]]][2] <- df_tampered[[v_econ[1]]][2] * 1.01
expect_error_matching(check_common_economic_draws(df_tampered, psa_shift$psa_params, v_econ),
                      "does not reuse the E0 economic draws")
cat("3. the survival draws are swapped for the evaluation only (restored on error); differing economic draws stop\n")

# ---------------------------------------------------------------------------
# 4. Cache: run and save, read back, stale on a changed state
# ---------------------------------------------------------------------------
cache_tmp <- file.path(tempdir(), "voi_evidence_psa_cache")
dir.create(cache_tmp, showWarnings = FALSE)
path_shift <- evidence_psa_cache_path("E2shift", "synthetic", cache_tmp)
unlink(path_shift)
#' load_or_run_evidence_psa() quietly with the synthetic inputs.
cached_psa <- function(state, path) {
  out <- NULL
  capture.output(out <- suppressMessages(suppressWarnings(load_or_run_evidence_psa(
    state, params, config$distributions, v_strategies, horizon, 1 / 52, n_draw, psa_seed,
    base_sampling_models = sampling_e0, base_psa_obj = pair_e0$psa_obj, base_psa_params = pair_e0$psa_params,
    path = path))))
  out
}
first <- cached_psa(state_shift, path_shift)
second <- cached_psa(state_shift, path_shift)
stopifnot(!first$from_cache, file.exists(path_shift), second$from_cache,
          identical(second$psa_obj$pair_metadata$generation_id, first$psa_obj$pair_metadata$generation_id),
          psa_outcome_gap(second$psa_obj, psa_shift$psa_obj) == 0,
          identical(basename(path_shift), "psa_evidence_E2shift_synthetic.rds"))
state_changed <- state_shift
state_changed$fingerprint <- "synthetic-state-shift-v2"
identity_changed <- evidence_psa_identity(
  state_changed, psa_cache_fingerprint(state_changed$fingerprint, params, config$distributions, v_strategies,
                                       n_draw, psa_seed, horizon, 1 / 52)$inputs, sampling_e0, pair_e0$psa_obj)
expect_error_matching(load_evidence_psa(path_shift, identity_changed, required = TRUE),
                      "upstream\\$evidence_states\\$E2shift")
third <- cached_psa(state_changed, path_shift)
stopifnot(!third$from_cache)
expect_error_matching(cached_psa(state_shift, psa_obj_path(directory = cache_tmp)), "Refusing to use the E0 PSA path")
expect_error_matching(evidence_psa_cache_path("E2 bad"), "letters, digits")
unlink(cache_tmp, recursive = TRUE)
cat("4. cache: saved and read back by identity; a changed evidence state is named and rerun; E0 path refused\n")

# ---------------------------------------------------------------------------
# 5. The resample check against the state's own PSA
# ---------------------------------------------------------------------------
pairing_e0 <- policy_price_pairing(policy_space_from_psa(pair_e0$psa_obj, params))
pairing_none <- policy_price_pairing(policy_space_from_psa(psa_none$psa_obj, params))
df_check <- evidence_state_comparison(pairing_e0, seq_len(pair_e0$psa_obj$n_sim), pairing_none, c(0, 13923), 51000)
pairing_shift <- policy_price_pairing(policy_space_from_psa(psa_shift$psa_obj, params))
df_check_shift <- evidence_state_comparison(pairing_e0, seq_len(pair_e0$psa_obj$n_sim), pairing_shift, 0, 51000)
stopifnot(nrow(df_check) == 2, max(abs(df_check$evpi_resampled - df_check$evpi_psa)) < 1e-6,
          max(df_check$max_abs_inmb_difference) < 1e-6, df_check_shift$max_abs_inmb_difference > 1e-6)
cat("5. resample check: equal to the state's PSA when the state is E0, different under evidence\n")

cat("\nAll evidence-state PSA tests passed.\n")
