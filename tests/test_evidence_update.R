# ===============================================================================
# TEST: conjugate evidence updating of the coefficient draws (issue #202, step 3)
# ===============================================================================
# Synthetic data only (simulated patients without identifiers, a synthetic
# joint gamma component in the sampling-cache layout, synthetic evidence
# files in a temporary directory):
#
#   1. conjugate_update(): an infinite-variance likelihood returns the prior
#      exactly; a 2-parameter normal-normal example equals the precision-form
#      closed form (1e-12); two observations at once equal two sequential
#      updates.
#   2. Sceptical prior: SD log(2) / qnorm(0.95) gives P(rate ratio < 0.5) = 5%;
#      the calibrated variant stays centred at 0, reproduces the target
#      shrinkage (1e-8) and, for a coefficient uncorrelated with the others,
#      the closed form sd^2 = s V / (1 - s).
#   3. Cox limit g_k: for a homogeneous exponential population it equals the
#      log rate ratio (1e-4) at 1:1 and 2:1 allocation, with and without
#      accrual and a maintenance landmark; it does not depend on the order of
#      the reference patients.
#   4. Shared heterogeneity: with many identical trials the posterior variance
#      of the observed combination approaches the floor set by tau_shared, not
#      zero.
#   5. Evidence states on the synthetic component: an external all-comers PFS
#      and OS result (nonlinear g_k, relinearised to convergence) and a
#      sceptical prior on the four interactions pass the importance-weight
#      check; the ESS is reported; a wrong conjugate mean stops the check.
#   6. Regenerated draws: method "mvn_joint_v2", posterior-centred fits,
#      sampled_survival_models() works, draw means within 4 SE of the
#      posterior mean; the fingerprint changes with tau; save_evidence_state()
#      refuses the base sampling-cache path and loads only on a matching
#      fingerprint.
#   7. External-evidence reader: a missing required column stops; an 80% CI
#      is converted with z = 1.2816; NA estimates and rows without a CI are
#      skipped with a message; allocation ratios are parsed as
#      experimental:control. The real table, when present, passes the
#      contract, and its rows with a lower limit printed as 0 (infinite SE)
#      are skipped by the likelihood and excluded, by name, by the selection.
#   8. SLOW (METIMMOX_SLOW_TESTS=1 only): g_k against coxph() fitted to large
#      trials simulated from three draws (and one maintenance design).
#   9. select_evidence_rows() (ruling (2) of 11 October 2026) on a synthetic
#      table that exercises every rule: single-arm and non-randomised designs
#      excluded first (#214), pMMR/MSS over ITT, the ITT row's dMMR
#      share recorded, latest data cut (undated cuts last), 95% over 80% for
#      the same analysis, exclusion of interaction P values, non-HR and
#      adjusted-HR rows, zero lower limits, secondary duplicates, rows without
#      follow-up or allocation, other populations and biomarker subgroups;
#      subgroups returned apart with only a TMB >= 9 mut/Mb-or-BRAF definition
#      matching; the real table gives at most one row per trial, population
#      and endpoint and no matching subgroup.
#  10. Ruling (6): an external result four prior SDs of g_k away with an SE
#      a tenth of that SD stops the importance-weight check with its
#      discrepancy message; one SD away with SE equal to that SD passes.
#
# Run: "C:\Program Files\R\R-4.3.2\bin\x64\Rscript.exe" tests/test_evidence_update.R
# ===============================================================================

suppressPackageStartupMessages({
  library(survival)
  library(flexsurv)
})
source(here::here("R", "pfs_os_violation_diagnostics.R"))
source(here::here("R", "joint_survival_sampling.R"))
source(here::here("R", "prediction_functions.R"))
source(here::here("R", "voi_survival_parameters.R"))
source(here::here("R", "cache_paths.R"))
source(here::here("R", "design_information.R"))
source(here::here("R", "voi_cache.R"))
source(here::here("R", "evidence_update.R"))

#' Stop unless `expr` raises an error whose message matches `pattern`.
expect_error_matching <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  if (!inherits(error, "error") || !grepl(pattern, conditionMessage(error))) {
    stop("Expected an error matching '", pattern, "', got: ",
         if (inherits(error, "error")) conditionMessage(error) else "no error")
  }
  invisible(TRUE)
}

#' Simulated patients with the covariates of the economic model (no identifiers).
simulate_patients <- function(n, seed) {
  set.seed(seed)
  data.frame(Age = stats::rnorm(n, 62, 9), sex = factor(stats::rbinom(n, 1, 0.45)),
             Rx = factor(stats::rbinom(n, 1, 0.5), labels = c("Control", "Experimental arm")),
             crp = factor(stats::rbinom(n, 1, 0.35)), tmb_braf = factor(stats::rbinom(n, 1, 0.45)))
}

#' Synthetic joint gamma component in the sampling-cache layout.
make_component <- function(n = 300, seed = 11, n_draws = 2000, rho = 0.4) {
  df <- simulate_patients(n, seed)
  v_eta <- log(0.02) + 0.01 * (df$Age - 62) - 0.2 * (df$Rx != "Control") + 0.3 * (df$crp == "1") -
    0.3 * (df$crp == "1") * (df$Rx != "Control")
  v_t <- stats::rgamma(n, 1.6, exp(v_eta))
  v_c <- stats::runif(n, 30, 150)
  df$OSwk <- pmin(v_t, v_c); df$Death <- as.numeric(v_t <= v_c)
  v_t2 <- pmin(v_t, stats::rgamma(n, 1.2, exp(v_eta + 0.8)))
  df$PFSwk <- pmin(v_t2, v_c); df$Prog <- as.numeric(v_t2 <= v_c)
  rhs <- "Age + sex + Rx + crp * Rx + tmb_braf * Rx"
  fit_os <- flexsurvreg(stats::as.formula(paste("Surv(OSwk, Death) ~", rhs)), data = df, dist = "gamma")
  fit_pfs <- flexsurvreg(stats::as.formula(paste("Surv(PFSwk, Prog) ~", rhs)), data = df, dist = "gamma")
  m_ro <- covariance_root(fit_os$cov); m_rp <- covariance_root(fit_pfs$cov)
  m_joint <- rbind(cbind(fit_os$cov, rho * m_ro %*% m_rp), cbind(rho * m_rp %*% m_ro, fit_pfs$cov))
  v_names <- c(paste0("os::", rownames(fit_os$res)), paste0("pfs::", rownames(fit_pfs$res)))
  dimnames(m_joint) <- list(v_names, v_names)
  list(component = list(method = "mvn_joint_v2", joint_covariance = list(covariance = m_joint),
                        draws = draw_joint_survival_coefficients(fit_os, fit_pfs, m_joint, n_draws, seed),
                        original_os = fit_os, original_pfs = fit_pfs, n_samples = n_draws, n_failed = 0L,
                        dist_os = "gamma", dist_pfs = "gamma", seed = seed),
       data = df)
}

# ---------------------------------------------------------------------------
# 1. conjugate_update()
# ---------------------------------------------------------------------------
v_m <- c(a = 0.3, b = -0.2)
m_s <- matrix(c(0.5, 0.2, 0.2, 0.3), 2, dimnames = list(names(v_m), names(v_m)))
prior_only <- conjugate_update(v_m, m_s, c(1, 0), 1.5, Inf)
stopifnot(identical(prior_only$mean, v_m), identical(prior_only$cov, m_s))
m_a <- rbind(c(1, 0), c(1, 1)); v_y <- c(0.1, 0.4); m_v <- diag(c(0.2, 0.4))
post <- conjugate_update(v_m, m_s, m_a, v_y, m_v)
m_prec <- solve(m_s) + t(m_a) %*% solve(m_v) %*% m_a
v_closed <- solve(m_prec, solve(m_s) %*% v_m + t(m_a) %*% solve(m_v) %*% v_y)
stopifnot(max(abs(post$mean - v_closed)) < 1e-12, max(abs(post$cov - solve(m_prec))) < 1e-12)
step1 <- conjugate_update(v_m, m_s, m_a[1, ], v_y[1], m_v[1, 1])
step2 <- conjugate_update(step1$mean, step1$cov, m_a[2, ], v_y[2], m_v[2, 2])
stopifnot(max(abs(step2$mean - post$mean)) < 1e-12, max(abs(step2$cov - post$cov)) < 1e-12)
expect_error_matching(conjugate_update(v_m, m_s, m_a, v_y, matrix(c(Inf, 0.1, 0.1, 1), 2)), "uncorrelated")
cat("1. conjugate update: infinite variance returns the prior; normal-normal closed form (1e-12)\n")

# ---------------------------------------------------------------------------
# 2. Sceptical prior
# ---------------------------------------------------------------------------
sd_sceptical <- sceptical_prior_sd()
stopifnot(abs(sd_sceptical - log(2) / qnorm(0.95)) < 1e-15,
          abs(stats::pnorm(log(0.5), 0, sd_sceptical) - 0.05) < 1e-12)
v_m3 <- c("os::x" = -1.2, "os::y" = 0.4, "os::z" = 0.1)
m_s3 <- diag(c(0.6, 0.3, 0.2)); dimnames(m_s3) <- list(names(v_m3), names(v_m3))
m_s3["os::y", "os::z"] <- m_s3["os::z", "os::y"] <- 0.05
cal <- calibrate_sceptical_sd(v_m3, m_s3, "os::x", "os::x", target_shrinkage = 0.3)
stopifnot(abs(cal$shrinkage - 0.3) < 1e-8, abs(cal$sd^2 - 0.3 * 0.6 / 0.7) < 1e-8,
          abs(cal$update$mean[["os::x"]] - 0.3 * v_m3[["os::x"]]) < 1e-8)
cal2 <- calibrate_sceptical_sd(v_m3, m_s3, c("os::x", "os::y"), "os::x", target_shrinkage = 0.5)
stopifnot(abs(cal2$shrinkage - 0.5) < 1e-8, cal2$probability_below_half > 0, cal2$probability_below_half < 0.5)
terms <- sceptical_prior_terms(c("os::x", "os::y"), names(v_m3), 0.4)
stopifnot(all(terms$y == 0))
expect_error_matching(calibrate_sceptical_sd(v_m3, m_s3, "os::x", "os::x", 1.2), "target_shrinkage")
cat(sprintf("2. sceptical prior SD %.5f (P(rate ratio < 0.5) = 5%%); calibrated SD reproduces the shrinkage\n",
            sd_sceptical))

# ---------------------------------------------------------------------------
# 3. Cox limit for a homogeneous exponential population
# ---------------------------------------------------------------------------
synthetic <- make_component()
component <- synthetic$component
df_reference_all <- synthetic$data[, c("Age", "sex", "Rx", "crp", "tmb_braf")]
set.seed(7)
df_reference <- df_reference_all[sample.int(nrow(df_reference_all), 60), ]
prior <- joint_prior_moments(component)
theta_exp <- prior$mean
theta_exp["pfs::shape"] <- 0
theta_exp["pfs::rate"] <- log(0.02)
df_one <- df_reference[1, , drop = FALSE]
df_one$crp <- factor("0", levels = c("0", "1")); df_one$tmb_braf <- factor("1", levels = c("0", "1"))
true_log_rr <- theta_exp[["pfs::RxExperimental arm"]] + theta_exp[["pfs::RxExperimental arm:tmb_braf1"]]
for (cfg in list(list(p = 0.5, F = 60, A = 0, L = 0), list(p = 2 / 3, F = 40, A = 100, L = 0),
                 list(p = 2 / 3, F = 40, A = 100, L = 16))) {
  spec <- trial_summary_spec("exp", "pfs", cfg$p, cfg$F, cfg$A, cfg$L)
  g <- make_trial_summary_map(spec, component, df_one)
  if (abs(g(theta_exp) - true_log_rr) > 1e-4) stop("Cox limit ", g(theta_exp), " differs from the log rate ratio ", true_log_rr)
}
spec_all <- trial_summary_spec("all-comers PFS", "pfs", 2 / 3, 60, 80)
g_all <- make_trial_summary_map(spec_all, component, df_reference)
g_perm <- make_trial_summary_map(spec_all, component, df_reference[rev(seq_len(nrow(df_reference))), ])
stopifnot(abs(g_all(prior$mean) - g_perm(prior$mean)) < 1e-12)
cat("3. Cox limit equals the log rate ratio under proportional hazards (1e-4) and ignores patient order\n")

# ---------------------------------------------------------------------------
# 4. Shared heterogeneity is a floor
# ---------------------------------------------------------------------------
n_trials <- 200
df_obs <- data.frame(label = paste0("trial", seq_len(n_trials)), trial = paste0("trial", seq_len(n_trials)),
                     endpoint = "pfs", se = 0.1, type = "external", stringsAsFactors = FALSE)
tau_s <- 0.15
m_v_many <- evidence_covariance(df_obs, tau_shared = tau_s, tau_trial = 0.05)
v_a <- c(1, 0); prior_var <- m_s[1, 1]
post_many <- conjugate_update(v_m, m_s, matrix(v_a, n_trials, 2, byrow = TRUE), rep(0.1, n_trials), m_v_many)
floor_var <- 1 / (1 / prior_var + 1 / tau_s^2)
stopifnot(post_many$cov[1, 1] >= floor_var - 1e-12, post_many$cov[1, 1] < 1.01 * floor_var + 1e-4)
cat(sprintf("4. 200 identical trials: posterior variance %.5f against the tau_shared floor %.5f\n",
            post_many$cov[1, 1], floor_var))

# ---------------------------------------------------------------------------
# 5. Evidence states and the importance-weight check
# ---------------------------------------------------------------------------
df_evidence <- data.frame(
  label = c("Trial A | ITT | PFS", "Trial A | ITT | OS"), trial = "Trial A", population = "ITT",
  endpoint = c("pfs", "os"), log_hr = log(c(0.8, 0.95)), se = c(0.12, 0.15), z = 1.96,
  p_experimental = 2 / 3, follow_up_weeks = 90, maintenance = FALSE, source_row = 1:2,
  accrual_weeks = 70, stringsAsFactors = FALSE)
state_e2 <- build_evidence_state(component, df_reference, evidence = df_evidence, tau_shared = 0.05,
                                 tau_trial = 0.05, within_trial_correlation = 0.3, label = "E2_synthetic",
                                 base_fingerprint = "synthetic", step = 2)
stopifnot(all(state_e2$check$table$agree), state_e2$update$iterations >= 2,
          state_e2$check$ess > 0, state_e2$check$ess <= nrow(component$draws$os))
v_moved <- abs(state_e2$observations$predicted_posterior - state_e2$observations$observed) <
  abs(state_e2$observations$predicted_prior - state_e2$observations$observed)
stopifnot(all(v_moved))
cat(sprintf("5. E2 (two nonlinear g_k): %d linearisations, ESS %.0f of %d, max |z| %.2f (threshold %.2f)\n",
            state_e2$update$iterations, state_e2$check$ess, nrow(component$draws$os),
            max(state_e2$check$table$z), state_e2$check$threshold))
v_interactions <- joint_coefficient_names(component)[grepl(":", sub("^[a-z]+::", "", joint_coefficient_names(component)))]
stopifnot(length(v_interactions) == 4)
state_e1 <- build_evidence_state(component, df_reference, sceptical = list(coefs = v_interactions, sd = sd_sceptical),
                                 label = "E1_synthetic", base_fingerprint = "synthetic")
direct <- apply_sceptical_prior(prior$mean, prior$cov, v_interactions, sd_sceptical)
stopifnot(max(abs(state_e1$update$mean - direct$mean)) < 1e-10, state_e1$update$iterations <= 2,
          all(state_e1$check$table$agree))
cat(sprintf("5. E1 (sceptical prior on %d interactions): ESS %.0f, max |z| %.2f\n", length(v_interactions),
            state_e1$check$ess, max(state_e1$check$table$z)))
obs_e1 <- lapply(v_interactions, function(nm) list(label = nm, y = 0, a = stats::setNames(1, nm)))
m_v_e1 <- diag(sd_sceptical^2, length(v_interactions))
expect_error_matching(importance_weight_check(joint_draw_matrix(component), obs_e1, m_v_e1,
                                              conjugate_mean = prior$mean, prior_mean = prior$mean),
                      "disagree")
cat("5. a wrong conjugate mean stops the importance-weight check\n")

# ---------------------------------------------------------------------------
# 6. Regenerated draws and the evidence-state cache
# ---------------------------------------------------------------------------
joint_e2 <- state_e2$state$joint
stopifnot(identical(joint_e2$method, "mvn_joint_v2"),
          max(abs(c(joint_e2$original_os$opt$par, joint_e2$original_pfs$opt$par) - state_e2$update$mean)) < 1e-12)
model_draw <- sampled_survival_models(joint_e2, 5)
stopifnot(inherits(model_draw$os$model, "flexsurvreg"), isFALSE(model_draw$failed))
m_post_draws <- joint_draw_matrix(joint_e2)
v_se <- sqrt(diag(state_e2$update$cov) / nrow(m_post_draws))
stopifnot(all(abs(colMeans(m_post_draws) - state_e2$update$mean) < 4 * v_se))
state_e2b <- build_evidence_state(component, df_reference, evidence = df_evidence, tau_shared = 0.1,
                                  tau_trial = 0.05, label = "E2_synthetic", base_fingerprint = "synthetic", step = 2)
stopifnot(!identical(state_e2$state$fingerprint, state_e2b$state$fingerprint))
cache_tmp <- file.path(tempdir(), "voi_evidence_cache")
dir.create(cache_tmp, showWarnings = FALSE)
expect_error_matching(save_evidence_state(state_e2$state, sampling_cache_path(joint_e2$n_samples, cache_tmp)),
                      "Refusing to overwrite")
path_e2 <- save_evidence_state(state_e2$state, evidence_state_cache_path("E2_synthetic", joint_e2$n_samples, cache_tmp))
stopifnot(!is.null(load_evidence_state(path_e2, state_e2$state$fingerprint)),
          is.null(suppressMessages(load_evidence_state(path_e2, state_e2b$state$fingerprint))))
cat("6. regenerated mvn_joint_v2 draws centred at the posterior; cache refuses the base path\n")

# ---------------------------------------------------------------------------
# 7. External-evidence reader
# ---------------------------------------------------------------------------
df_csv <- data.frame(
  trial = c("Trial A", "Trial A", "Trial B", "Trial C", "Trial D"),
  population = c("ITT", "ITT", "ITT", "ITT", "ITT"),
  endpoint = c("PFS (BICR)", "OS", "PFS", "PFS", "ORR"),
  effect_measure = c("HR", "HR", "HR", "HR", "proportion_experimental_percent"),
  estimate = c(0.8, NA, 0.7, 0.9, 50), ci_lower = c(0.6, NA, 0.55, NA, 40), ci_upper = c(1.07, NA, 0.89, NA, 60),
  ci_level = c(80, NA, 95, NA, 95),
  allocation_ratio = c("2:1 (experimental:control)", "2:1 (experimental:control)", "1:2 (control:experimental)",
                       "1:1", "1:1"),
  randomisation_landmark = c(NA, NA, NA, "after induction", NA),
  median_follow_up_months = c("12 (IQR 8-16)", "12", "24", "10.5", "10"), stringsAsFactors = FALSE)
path_csv <- file.path(tempdir(), "synthetic_external_evidence.csv")
utils::write.csv(df_csv, path_csv, row.names = FALSE)
df_read <- read_external_evidence(path_csv)
stopifnot(nchar(attr(df_read, "md5")) == 32)
v_messages <- character(0)
df_like <- withCallingHandlers(external_evidence_likelihood(df_read),
                               message = function(m) { v_messages <<- c(v_messages, conditionMessage(m)); invokeRestart("muffleMessage") })
stopifnot(nrow(df_like) == 2, identical(df_like$trial, c("Trial A", "Trial B")),
          abs(df_like$se[1] - (log(1.07) - log(0.6)) / (2 * 1.2816)) < 1e-4,
          abs(df_like$z[1] - 1.2816) < 1e-4, abs(df_like$se[2] - (log(0.89) - log(0.55)) / (2 * qnorm(0.975))) < 1e-12,
          abs(df_like$p_experimental[1] - 2 / 3) < 1e-12, abs(df_like$p_experimental[2] - 2 / 3) < 1e-12,
          abs(df_like$follow_up_weeks[1] - 12 * 52 / 12) < 1e-12,
          any(grepl("no estimate", v_messages)), any(grepl("no confidence interval", v_messages)),
          any(grepl("effect measure", v_messages)))
stopifnot(identical(parse_allocation_ratio(c("1:1", "3:2 (experimental:control)", NA)), c(0.5, 0.6, NA)))
utils::write.csv(df_csv[, setdiff(names(df_csv), "ci_level")], path_csv, row.names = FALSE)
expect_error_matching(read_external_evidence(path_csv), "lacks required column\\(s\\): ci_level")
df_bad <- df_csv; df_bad$ci_upper[1] <- 0.7
utils::write.csv(df_bad, path_csv, row.names = FALSE)
expect_error_matching(read_external_evidence(path_csv), "estimate outside its CI")
path_real <- here::here(EXTERNAL_EVIDENCE_PATH)
if (file.exists(path_real)) {
  raw_real <- read_external_evidence(path_real)
  df_real <- suppressMessages(external_evidence_likelihood(raw_real))
  # The likelihood rows that enter an update come from the selection rule; rows
  # whose lower limit prints as 0 (infinite SE) are excluded by it and named.
  sel7 <- select_evidence_rows(raw_real)
  v_zero <- which(!is.na(raw_real$ci_lower) & raw_real$ci_lower <= 0)
  stopifnot(nrow(df_real) > 0, !any(df_real$source_row %in% v_zero),
            nrow(sel7$selected) > 0, all(is.finite(sel7$selected$se)), all(sel7$selected$se > 0),
            all(sel7$selected$p_experimental > 0 & sel7$selected$p_experimental < 1),
            all(v_zero %in% sel7$excluded$source_row),
            all(grepl("lower limit prints as 0",
                      sel7$excluded$exclusion_reason[match(v_zero, sel7$excluded$source_row)])))
  cat(sprintf(paste0("7. external-evidence contract: synthetic checks pass; the real table gives %d usable HR rows; ",
                     "the selection excludes the %d rows with a zero lower limit (rows %s)\n"),
              nrow(df_real), length(v_zero), paste(v_zero, collapse = ", ")))
} else {
  cat("7. external-evidence contract: synthetic checks pass (real table absent)\n")
}

# ---------------------------------------------------------------------------
# 8. SLOW: g_k against Cox fits to simulated trials
# ---------------------------------------------------------------------------
if (identical(Sys.getenv("METIMMOX_SLOW_TESTS"), "1")) {
  #' Simulate one large trial from a draw for a trial_summary_spec().
  simulate_cox_trial <- function(spec, theta, n, seed) {
    set.seed(seed)
    fit <- component[[paste0("original_", spec$endpoint)]]
    par <- split_joint_draw(theta, component)
    v_w <- rep(1, nrow(df_reference))
    if (spec$landmark_weeks > 0) {
      v_eta_sel <- gamma_linear_predictor(component$original_pfs, par$pfs, df_reference,
                                          model_rx_levels(component$original_pfs)[1])
      v_w <- stats::pgamma(spec$landmark_weeks, exp(par$pfs[["shape"]]), exp(v_eta_sel), lower.tail = FALSE)
    }
    df <- df_reference[sample.int(nrow(df_reference), n, replace = TRUE, prob = v_w), ]
    v_arm <- stats::rbinom(n, 1, spec$p_experimental)
    df$Rx <- factor(model_rx_levels(fit)[v_arm + 1], levels = model_rx_levels(fit))
    v_eta <- gamma_linear_predictor(fit, par[[spec$endpoint]], df)
    k <- exp(par[[spec$endpoint]][["shape"]])
    s_land <- stats::pgamma(spec$landmark_weeks, k, exp(v_eta), lower.tail = FALSE)
    v_time <- stats::qgamma(stats::runif(n) * s_land, k, exp(v_eta), lower.tail = FALSE) - spec$landmark_weeks
    v_cens <- spec$min_follow_up_weeks + spec$accrual_weeks - stats::runif(n, 0, spec$accrual_weeks)
    data.frame(time = pmin(v_time, v_cens), status = as.numeric(v_time <= v_cens), arm = v_arm)
  }
  df_slow <- data.frame()
  v_specs <- list(spec_all, trial_summary_spec("all-comers OS", "os", 0.5, 100, 60),
                  trial_summary_spec("maintenance PFS", "pfs", 2 / 3, 50, 80, landmark_weeks = 16))
  for (spec in v_specs) for (idx in c(1, 2, 3)) {
    theta <- joint_draw_matrix(component)[idx, ]
    g <- make_trial_summary_map(spec, component, df_reference)
    df_sim <- simulate_cox_trial(spec, theta, 40000, 1000 + idx)
    cox <- survival::coxph(survival::Surv(time, status) ~ arm, data = df_sim)
    df_slow <- rbind(df_slow, data.frame(spec = spec$label, draw = idx, g = g(theta), cox = unname(coef(cox)),
                                         se = sqrt(vcov(cox)[1, 1])))
  }
  df_slow$z <- (df_slow$cox - df_slow$g) / df_slow$se
  print(df_slow, row.names = FALSE)
  if (any(abs(df_slow$z) > 3.5)) stop("g_k disagrees with the Cox fits to simulated trials")
  cat("8. SLOW: g_k agrees with coxph() on simulated trials within 3.5 SE\n")
} else {
  cat("8. SLOW g_k simulation check skipped (set METIMMOX_SLOW_TESTS=1)\n")
}


# ---------------------------------------------------------------------------
# 9. Evidence row selection (ruling (2) of 11 October 2026)
# ---------------------------------------------------------------------------
df_sel <- data.frame(
  trial = c(rep("Trial A", 14), rep("Trial B", 4), "Trial C"),
  population = c("ITT", "pMMR subgroup", "pMMR subgroup", "pMMR subgroup", "ITT", "ITT", "ITT", "pMMR subgroup",
                 "pMMR subgroup", "pMMR subgroup", "pMMR subgroup", "pMMR cohort, TMB subgroups", "pMMR cohort",
                 "pMMR subgroup", "ITT", "ITT", "ITT", "ITT", "Per-protocol"),
  subgroup_definition = c(NA, "proficient mismatch repair", "proficient mismatch repair",
                          "proficient mismatch repair", NA, NA, NA, "proficient mismatch repair",
                          "pMMR, TMB-high (>= 10 mutations/Mb, FoundationOne CDx)", "pMMR, TMB-high (>= 9 mut/Mb)",
                          "pMMR, TMB >= 9 mut/Mb or any BRAF mutation", "TMB-high (>= 10 mut/Mb)",
                          "pMMR, Immunoscore IC high", "proficient mismatch repair", NA, NA, NA, NA, NA),
  endpoint = c("PFS", "PFS", "PFS", "PFS (updated)", "OS", "OS", "ORR", "OS", "PFS", "PFS", "PFS", "OS", "PFS",
               "PFS", "PFS", "PFS", "PFS", "OS", "PFS"),
  effect_measure = c("HR", "HR", "HR", "HR", "HR", "HR (adjusted)", "proportion_experimental_percent", "HR", "HR",
                     "HR", "HR", "treatment x TMB interaction p-value", "HR", "HR", "HR", "HR", "HR", "HR", "HR"),
  estimate = c(0.8, 0.78, 0.77, 0.75, 0.9, 0.88, 50, NA, 0.05, 0.5, 0.45, 0.04, 0.5, 0.78, 0.9, 0.92, 0.95, 0.85, 0.7),
  ci_lower = c(0.6, 0.62, 0.55, 0.5, 0.7, 0.68, 40, NA, 0, 0.2, 0.2, NA, 0.25, 0.62, 0.7, 0.72, 0.75, 0.6, 0.5),
  ci_upper = c(1.07, 0.97, 1.08, 1.12, 1.15, 1.13, 60, NA, 0.4, 1.2, 1.0, NA, 1.0, 0.97, 1.15, 1.17, 1.2, 1.2, 0.98),
  ci_level = c(95, 80, 95, 95, 95, 95, 95, NA, 95, 95, 95, NA, 95, 80, 95, 95, 95, 95, 95),
  allocation_ratio = c(rep("2:1 (experimental:control)", 14), rep("1:1", 3), NA, "1:1"),
  randomisation_landmark = c(rep(NA, 14), rep("after induction", 4), NA),
  median_follow_up_months = c(rep("20", 14), "18", "20.3", NA, "20.3", "12"),
  data_cut = c("2022-06-30", "2022-06-30", "2022-06-30", "2020-01-01", "2022-06-30", "2022-06-30", "2022-06-30",
               "2022-06-30", "2022-06-30", "2022-06-30", "2022-06-30", "2022-06-30", "2022-06-30", "2022-06-30",
               "not stated (after 2019-05-31)", "2019-05-31", "2019-05-31", "2019-05-31", "2021-01-01"),
  read_from = c(rep("full text", 13), "secondary full text", "abstract", "full text", "full text", "registry",
                "full text"),
  note = c(rep(NA, 13), "duplicate of the primary row, do not use", rep(NA, 5)),
  dmmr_share_in_population = c("10/200 dMMR (4 control, 6 experimental)", NA, NA, NA,
                               "10/200 dMMR (4 control, 6 experimental)", rep(NA, 9), "7 MSI of 378 tested",
                               "7 MSI of 378 tested (5/252, 2/126)", NA, NA, NA),
  stringsAsFactors = FALSE)
path_sel <- file.path(tempdir(), "synthetic_evidence_selection.csv")
utils::write.csv(df_sel, path_sel, row.names = FALSE)
raw_sel <- read_external_evidence(path_sel)
sel <- select_evidence_rows(raw_sel, subgroups = TRUE)
reason_of <- function(row) sel$excluded$exclusion_reason[sel$excluded$source_row == row]
stopifnot(identical(sel$selected$source_row, c(5L, 3L, 16L)),
          identical(sel$selected$endpoint, c("os", "pfs", "pfs")),
          identical(sel$selected$population_class, c("itt", "pmmr", "itt")),
          !anyDuplicated(paste(sel$selected$trial, sel$selected$endpoint)),
          setequal(c(sel$selected$source_row, sel$excluded$source_row), seq_len(nrow(df_sel))),
          # pMMR preferred, ITT then not used
          grepl("ITT row not used: the pMMR/MSS row 3", reason_of(1)),
          # 95% preferred over 80% for the same analysis and cut; earlier cut dropped
          grepl("^80% interval; row 3", reason_of(2)), grepl("earlier data cut", reason_of(4)),
          # ITT used where the pMMR row has no estimate; its dMMR share recorded
          grepl("no estimate", reason_of(8)), abs(sel$selected$dmmr_share[1] - 0.05) < 1e-12,
          identical(sel$diagnostics$source_row, c(5L, 16L)), abs(sel$diagnostics$dmmr_share[2] - 7 / 378) < 1e-12,
          grepl("secondary analysis", reason_of(6)), grepl("not a hazard ratio", reason_of(7)),
          grepl("lower limit prints as 0", reason_of(9)), grepl("interaction P value", reason_of(12)),
          grepl("secondary duplicate", reason_of(14)), grepl("data cut not stated", reason_of(15)),
          grepl("no median follow-up", reason_of(17)), grepl("no allocation", reason_of(18)),
          grepl("neither ITT nor pMMR", reason_of(19)), grepl("biomarker subgroup \\(tmb\\)", reason_of(10)),
          grepl("biomarker subgroup \\(immunoscore\\)", reason_of(13)),
          # maintenance landmark
          sel$selected$landmark_weeks[3] == EVIDENCE_MAINTENANCE_LANDMARK_WEEKS, all(sel$selected$landmark_weeks[1:2] == 0),
          grepl("pMMR/MSS population; data cut 2022-06-30; 95% interval", sel$selected$selection_reason[2]))
# Subgroups: returned separately; only the combined TMB-or-BRAF definition at 9 mut/Mb matches.
sg <- sel$subgroups
stopifnot(identical(sort(sg$source_row), c(9L, 10L, 11L, 12L, 13L)),
          identical(sg$source_row[sg$matches_metimmox], 11L),
          grepl("10 mut/Mb \\(FoundationOne CDx\\)", sg$match_reason[sg$source_row == 9]),
          grepl("subset", sg$match_reason[sg$source_row == 10]),
          grepl("not a METIMMOX biomarker", sg$match_reason[sg$source_row == 13]),
          !sg$usable[sg$source_row == 9], grepl("lower limit", sg$unusable_reason[sg$source_row == 9]),
          is.null(select_evidence_rows(raw_sel)$subgroups))
# Maintenance rows need a landmark; the likelihood format gives the same selection.
sel_nm <- select_evidence_rows(raw_sel, maintenance_landmark_weeks = NA)
stopifnot(identical(sel_nm$selected$source_row, c(5L, 3L)),
          grepl("maintenance design", sel_nm$excluded$exclusion_reason[sel_nm$excluded$source_row == 16]))
sel_lik <- select_evidence_rows(suppressMessages(external_evidence_likelihood(raw_sel)))
# (Without the subgroup, cut and source columns the likelihood format cannot
# tell subgroups or cuts apart, so only the one-row-per-trial-and-endpoint
# guarantee is checked.)
stopifnot(nrow(sel_lik$selected) == 3, 5L %in% sel_lik$selected$source_row,
          !anyDuplicated(paste(sel_lik$selected$trial, sel_lik$selected$endpoint)))
# A single-arm or non-randomised design (#214) is excluded first, even when the
# row is otherwise a usable hazard ratio; "randomised" designs are unaffected.
df_sa <- df_sel[c(5, 5, 5), ]
df_sa$trial <- c("Trial D", "Trial E", "Trial F")
df_sa$design <- c("randomised open-label phase 2", "single-arm phase 1b/2 (no control arm)",
                  "HR from a pooled non-randomised cross-trial comparison")
sel_sa <- select_evidence_rows(df_sa)
stopifnot(identical(sel_sa$selected$trial, "Trial D"),
          identical(sel_sa$excluded$trial, c("Trial E", "Trial F")),
          all(grepl("^single-arm or non-randomised design", sel_sa$excluded$exclusion_reason)),
          identical(suppressMessages(external_evidence_likelihood(df_sa))$trial, "Trial D"),
          !grepl(EVIDENCE_NONRANDOMISED_PATTERN, "randomised double-blind placebo-controlled phase 2", ignore.case = TRUE))
# build_evidence_state() selects by default and keeps a selection as given.
stopifnot(!is.null(state_e2$selection), nrow(state_e2$selection$selected) == 2)
cat(paste0("9. evidence selection: pMMR over ITT, latest cut, 95% interval, exclusions with reasons; subgroups apart; ",
           "single-arm and non-randomised designs excluded first\n"))
if (file.exists(path_real)) {
  raw_real9 <- read_external_evidence(path_real)
  sel_real <- select_evidence_rows(raw_real9, subgroups = TRUE)
  tab <- table(sel_real$selected$endpoint)
  # Single-arm and non-randomised rows (POCHI, MEDITREME; #214) are documented
  # as excluded with the design reason and never selected.
  v_nr <- which(!is.na(raw_real9$design) & grepl(EVIDENCE_NONRANDOMISED_PATTERN, raw_real9$design, ignore.case = TRUE))
  stopifnot(!any(v_nr %in% sel_real$selected$source_row), all(v_nr %in% sel_real$excluded$source_row),
            all(grepl("^single-arm or non-randomised design",
                      sel_real$excluded$exclusion_reason[match(v_nr, sel_real$excluded$source_row)])))
  stopifnot(!anyDuplicated(paste(sel_real$selected$trial, sel_real$selected$population_class,
                                 sel_real$selected$endpoint)),
            !anyDuplicated(paste(sel_real$selected$trial, sel_real$selected$endpoint)),
            all(is.finite(sel_real$selected$se)), all(sel_real$selected$se > 0),
            nrow(sel_real$selected) + nrow(sel_real$excluded) == nrow(read_external_evidence(path_real)),
            !any(sel_real$subgroups$matches_metimmox))
  # Every trial with a usable pMMR/MSS row for an endpoint uses it.
  v_pmmr_avail <- unique(paste(sel_real$selected$trial[sel_real$selected$population_class == "pmmr"],
                               sel_real$selected$endpoint[sel_real$selected$population_class == "pmmr"]))
  stopifnot(all(!grepl("^ITT row not used", sel_real$excluded$exclusion_reason) |
                  paste(sel_real$excluded$trial, ifelse(grepl("^OS", sel_real$excluded$endpoint), "os", "pfs")) %in%
                    v_pmmr_avail))
  print(sel_real$selected[, c("trial", "endpoint", "source_row", "population_class", "selection_reason")],
        row.names = FALSE)
  cat(sprintf(paste0("9. real table: %d selected rows (%s), %d excluded (%d single-arm or non-randomised), ",
                     "%d subgroup rows (none matches METIMMOX)\n"),
              nrow(sel_real$selected), paste(names(tab), tab, sep = " ", collapse = ", "),
              nrow(sel_real$excluded), length(v_nr), nrow(sel_real$subgroups)))
}

# ---------------------------------------------------------------------------
# 10. A genuinely poor linearisation must stop the step (ruling (6) of
#     11 October 2026: the multiplicity-adjusted importance-weight check is
#     accepted as the linearisation check, with a test that it fails when the
#     linearisation is poor)
# ---------------------------------------------------------------------------
spec_far <- trial_summary_spec("Trial F | ITT | PFS", "pfs", 2 / 3, 60, 70, step = 2)
g_far <- make_trial_summary_map(spec_far, component, df_reference)
v_jac <- numerical_jacobian(function(theta) g_far(theta), prior$mean)
sd_g <- sqrt(as.numeric(v_jac %*% prior$cov %*% t(v_jac)))
g0 <- g_far(prior$mean)
evidence_at <- function(shift_sd, se) {
  data.frame(label = "Trial F | ITT | PFS", trial = "Trial F", population = "ITT", endpoint = "pfs",
             log_hr = g0 + shift_sd * sd_g, se = se, z = 1.96, p_experimental = 2 / 3, follow_up_weeks = 60,
             maintenance = FALSE, source_row = 1L, accrual_weeks = 70, stringsAsFactors = FALSE)
}
# Far: four prior SDs of g away with an SE a tenth of that SD; the external
# likelihood sits where the prior draws have almost no mass, the weights
# collapse and the weighted mean cannot follow the conjugate mean.
expect_error_matching(build_evidence_state(component, df_reference, evidence = evidence_at(-4, 0.1 * sd_g),
                                           label = "far", base_fingerprint = "synthetic", step = 2),
                      "Conjugate and importance-weighted posterior means disagree")
# Moderate: one prior SD away (towards a smaller benefit) with an SE equal to
# that SD passes. (Shifts of 0.5-1 SD towards a larger benefit give max z of
# 3.1-3.7 on pfs::shape here, just above the 3.02 threshold, for differences
# of about 0.05 posterior SD: the check also reacts to the curvature of g_k.)
state_moderate <- build_evidence_state(component, df_reference, evidence = evidence_at(1, sd_g),
                                       label = "moderate", base_fingerprint = "synthetic", step = 2)
stopifnot(all(state_moderate$check$table$agree), state_moderate$check$ess > 100)
cat(sprintf("10. linearisation check (ruling 6): a result 4 prior SDs away (SE %.3f) stops; 1 SD away (SE = SD) passes (ESS %.0f)\n",
            0.1 * sd_g, state_moderate$check$ess))

cat("\nAll evidence-update tests passed.\n")
