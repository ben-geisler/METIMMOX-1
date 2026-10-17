# ===============================================================================
# TEST: design information of the joint gamma models (issue #202, paper 3 step 2)
# ===============================================================================
# Synthetic data only (simulated patients without identifiers; no trial data,
# no caches):
#
#   1. Derivatives: the analytic event derivatives and the Richardson
#      finite-difference censored derivatives of the gamma log-likelihood in
#      (log shape, eta) equal numDeriv derivatives of dgamma/pgamma (1e-6).
#   2. Observed information by arm x cell: the strata sum to the numerical
#      Hessian of flexsurv's own log-likelihood at the estimates (scaled
#      difference < 1e-6) and to the inverse of the fit's covariance (< 5e-3:
#      flexsurv's covariance comes from optim's numerical Hessian), for a
#      right-censored and a left-truncated fit; stratum counts add up.
#   3. Expected information per patient equals the Monte Carlo mean of the
#      observed information over simulated event and censoring times (shape
#      below and above one, accrual on and off, landmark 0 and 4) within four
#      Monte Carlo SEs plus 0.3%.
#   4. Cox heuristic: under an exponential draw (shape 1, proportional hazards)
#      and the saturated Rx x CRP model, the interaction variance with the
#      shape held fixed equals sum(1 / expected events) over the four arm x
#      CRP strata within 1e-3; with the shape estimated it is at least that.
#   5. Enriched designs: TMB/BRAF-positive only identifies {Rx:crp,
#      cell_effect_01}, CRP-positive only {Rx:tmb_braf, cell_effect_10}; the
#      information has rank 7 of 9 and its null space aliases rate with
#      tmb_braf (or crp) and Rx with the interaction; a target that is not
#      identifiable stops.
#   6. Scaling: with fixed accrual the information is linear in n; with a
#      fixed accrual rate the information per patient grows with n.
#   7. Landmark randomisation reweights the cell mix by P(PFS > L).
#   8. Cross-endpoint covariance: "recolour" keeps the whitened OS/PFS
#      correlation of the prior targets; "independent" sets it to zero.
#   9. SLOW (METIMMOX_SLOW_TESTS=1 only): V_d(n) against the sampling
#      covariance of flexsurvreg refits over 500 simulated trials per n
#      (n = 200, 400, 800; the all-comers design randomised at week 4).
#
# Run: "C:\Program Files\R\R-4.3.2\bin\x64\Rscript.exe" tests/test_design_information.R
#      (slow part: set METIMMOX_SLOW_TESTS=1 first; about 15 minutes)
# ===============================================================================

suppressPackageStartupMessages({
  library(survival)
  library(flexsurv)
})
source(here::here("R", "pfs_os_violation_diagnostics.R"))
source(here::here("R", "joint_survival_sampling.R"))
source(here::here("R", "voi_survival_parameters.R"))
source(here::here("R", "design_information.R"))

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
make_component <- function(n = 300, seed = 11, n_draws = 500, rho = 0.4,
                           formula_rhs = "Age + sex + Rx + crp * Rx + tmb_braf * Rx") {
  df <- simulate_patients(n, seed)
  v_eta <- log(0.02) + 0.01 * (df$Age - 62) - 0.2 * (df$Rx != "Control") + 0.3 * (df$crp == "1") -
    0.3 * (df$crp == "1") * (df$Rx != "Control")
  v_t <- stats::rgamma(n, 1.6, exp(v_eta))
  v_c <- stats::runif(n, 30, 150)
  df$OSwk <- pmin(v_t, v_c); df$Death <- as.numeric(v_t <= v_c)
  v_t2 <- pmin(v_t, stats::rgamma(n, 1.2, exp(v_eta + 0.8)))
  df$PFSwk <- pmin(v_t2, v_c); df$Prog <- as.numeric(v_t2 <= v_c)
  fit_os <- flexsurvreg(stats::as.formula(paste("Surv(OSwk, Death) ~", formula_rhs)), data = df, dist = "gamma")
  fit_pfs <- flexsurvreg(stats::as.formula(paste("Surv(PFSwk, Prog) ~", formula_rhs)), data = df, dist = "gamma")
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

synthetic <- make_component()
component <- synthetic$component
df_reference <- synthetic$data[, c("Age", "sex", "Rx", "crp", "tmb_braf")]

# ---------------------------------------------------------------------------
# 1. Derivatives against numDeriv
# ---------------------------------------------------------------------------
v_a <- log(c(0.7, 1.6)); v_eta <- c(log(0.02), log(0.05)); v_t <- c(3, 40)
for (a in v_a) for (eta in v_eta) for (t in v_t) {
  for (ev in 0:1) {
    f_l <- function(p) if (ev == 1) stats::dgamma(t, exp(p[1]), exp(p[2]), log = TRUE) else
      stats::pgamma(t, exp(p[1]), exp(p[2]), lower.tail = FALSE, log.p = TRUE)
    m_h <- numDeriv::hessian(f_l, c(a, eta))
    v_g <- numDeriv::grad(f_l, c(a, eta))
    m_d <- gamma_loglik_derivatives(a, eta, t, ev)
    v_ours <- c(m_d[, "a"], m_d[, "eta"], m_d[, "aa"], m_d[, "aeta"], m_d[, "etaeta"])
    v_ref <- c(v_g, m_h[1, 1], m_h[1, 2], m_h[2, 2])
    if (max(abs(v_ours - v_ref) / pmax(1, abs(v_ref))) > 1e-6) {
      stop("Derivative mismatch (event = ", ev, ", a = ", a, ", eta = ", eta, ", t = ", t, "): ",
           paste(signif(v_ours - v_ref, 3), collapse = ", "))
    }
  }
}
cat("1. gamma log-likelihood derivatives agree with numDeriv (1e-6)\n")

# ---------------------------------------------------------------------------
# 2. Observed information by stratum
# ---------------------------------------------------------------------------
check_fit_information <- function(fit, df, formula, label) {
  s <- stratum_observed_information(fit)
  nll <- function(p) {
    refit <- flexsurvreg(formula, data = df, dist = "gamma", inits = c(exp(p[1:2]), p[-(1:2)]), fixedpars = TRUE)
    -refit$loglik
  }
  m_h <- numDeriv::hessian(nll, fit$opt$par)
  dimnames(m_h) <- dimnames(s$total)
  d_hess <- scaled_matrix_difference(s$total, m_h)
  stopifnot(d_hess < 1e-6, s$scaled_difference < 5e-3,
            sum(s$counts$n) == nrow(df), length(s$strata) == 8,
            isTRUE(all.equal(Reduce(`+`, s$strata), s$total)))
  cat(sprintf("2. %s: strata sum to the numerical Hessian (scaled difference %.1e) and to solve(vcov) (%.1e)\n",
              label, d_hess, s$scaled_difference))
  s
}
formula_os <- Surv(OSwk, Death) ~ Age + sex + Rx + crp * Rx + tmb_braf * Rx
s_os <- check_fit_information(component$original_os, synthetic$data, formula_os, "right-censored OS fit")
stopifnot(sum(s_os$counts$events) == sum(synthetic$data$Death))

df_trunc <- synthetic$data
df_trunc$entry <- pmin(df_trunc$OSwk / 2, 4)
formula_trunc <- Surv(entry, OSwk, Death) ~ Age + sex + Rx + crp * Rx + tmb_braf * Rx
fit_trunc <- flexsurvreg(formula_trunc, data = df_trunc, dist = "gamma")
invisible(check_fit_information(fit_trunc, df_trunc, formula_trunc, "left-truncated fit"))

# ---------------------------------------------------------------------------
# 3. Expected information against Monte Carlo
# ---------------------------------------------------------------------------
set.seed(3)
for (cfg in list(c(a = log(0.7), L = 0, F = 60, A = 100), c(a = log(1.8), L = 4, F = 50, A = 80),
                 c(a = log(1.2), L = 0, F = 40, A = 0))) {
  a <- cfg[["a"]]; landmark <- cfg[["L"]]; follow <- cfg[["F"]]; accrual <- cfg[["A"]]
  for (eta in c(log(0.015), log(0.03))) {
    m_e <- expected_gamma_information(a, eta, follow, accrual, landmark, step = 0.25)
    n_mc <- 200000
    k <- exp(a)
    s_land <- stats::pgamma(landmark, k, exp(eta), lower.tail = FALSE)
    v_time <- stats::qgamma(stats::runif(n_mc) * s_land, k, exp(eta), lower.tail = FALSE)
    v_cens <- follow + accrual - stats::runif(n_mc, 0, accrual)
    v_obs <- pmin(v_time, landmark + v_cens); v_ev <- as.numeric(v_time <= landmark + v_cens)
    m_d <- gamma_loglik_derivatives(a, eta, v_obs, v_ev)
    if (landmark > 0) {
      m_l <- gamma_censored_derivatives(a, eta, landmark)
      m_d[, c("aa", "aeta", "etaeta")] <- sweep(m_d[, c("aa", "aeta", "etaeta")], 2, m_l[1, c("aa", "aeta", "etaeta")])
    }
    for (j in c("aa", "aeta", "etaeta")) {
      mc <- -mean(m_d[, j]); se <- stats::sd(m_d[, j]) / sqrt(n_mc)
      if (abs(m_e[1, j] - mc) > 4 * se + 0.003 * abs(mc)) {
        stop("Expected information ", j, " = ", m_e[1, j], " differs from Monte Carlo ", mc, " (SE ", se, ")")
      }
    }
    if (abs(m_e[1, "events"] - mean(v_ev)) > 4 * sqrt(mean(v_ev) * (1 - mean(v_ev)) / n_mc) + 0.002) {
      stop("Expected events differ from Monte Carlo")
    }
  }
}
cat("3. expected information equals the Monte Carlo mean of the observed information (4 SE + 0.3%)\n")

# ---------------------------------------------------------------------------
# 4. Cox heuristic for an interaction under an exponential draw
# ---------------------------------------------------------------------------
exp_synth <- make_component(n = 300, seed = 21, formula_rhs = "Rx * crp")
exp_component <- exp_synth$component
draw_exp <- component_mean_draw(exp_component)
draw_exp$os[["shape"]] <- 0
draw_exp$pfs[["shape"]] <- 0
design_crp <- trial_design("crp_interaction", enrolled_cells = c("0", "1"), min_follow_up_weeks = 52,
                           accrual_weeks = 104, biomarkers = "crp")
info_exp <- design_information(design_crp, 400, draw_exp, exp_component,
                               exp_synth$data[, c("Rx", "crp")], cross_endpoint = "independent")
for (outcome in c("os", "pfs")) {
  m_info <- info_exp$endpoints[[outcome]]$information
  v_par <- rownames(m_info)
  int <- v_par[grepl(":", v_par)]
  var_fixed_shape <- solve(m_info[-1, -1])[int, int]
  var_free_shape <- solve(m_info)[int, int]
  heuristic <- cox_interaction_heuristic(info_exp, outcome, "crp")
  if (abs(var_fixed_shape / heuristic - 1) > 1e-3) {
    stop(outcome, ": shape-fixed interaction variance ", var_fixed_shape, " vs sum(1/d) ", heuristic)
  }
  stopifnot(var_free_shape >= var_fixed_shape * (1 - 1e-9))
  cat(sprintf("4. %s: interaction variance with shape fixed / sum(1/d) = %.6f; with shape estimated %.4f\n",
              outcome, var_fixed_shape / heuristic, var_free_shape / heuristic))
}

# ---------------------------------------------------------------------------
# 5. Enriched designs identify only part of the coefficients
# ---------------------------------------------------------------------------
draw_mean <- component_mean_draw(component)
design_tmb <- trial_design("tmb_braf_enriched", enrolled_cells = c("01", "11"), min_follow_up_weeks = 52,
                           accrual_weeks = 104)
info_tmb <- design_information(design_tmb, 80, draw_mean, component, df_reference)
stopifnot(identical(rownames(info_tmb$targets),
                    c("os::rx_x_crp", "os::cell_effect_01", "pfs::rx_x_crp", "pfs::cell_effect_01")),
          info_tmb$endpoints$os$rank == 7, info_tmb$endpoints$pfs$rank == 7,
          all(eigen(info_tmb$covariance, only.values = TRUE)$values > 0))
m_null <- info_tmb$endpoints$os$null_space
# The null space is spanned by (rate - tmb_braf1) and (Rx - Rx:tmb_braf1).
m_alias <- cbind(c(rate = 1, tmb_braf1 = -1), c(1, -1))
v_alias_rows <- list(c("rate", "tmb_braf1"), c("RxExperimental arm", "RxExperimental arm:tmb_braf1"))
m_basis <- matrix(0, nrow(m_null), 2, dimnames = list(rownames(m_null), NULL))
m_basis[v_alias_rows[[1]], 1] <- c(1, -1); m_basis[v_alias_rows[[2]], 2] <- c(1, -1)
stopifnot(max(abs(qr.resid(qr(m_basis), m_null))) < 1e-6)
design_crp_pos <- trial_design("crp_enriched", enrolled_cells = c("10", "11"), min_follow_up_weeks = 52,
                               accrual_weeks = 104)
info_crp_pos <- design_information(design_crp_pos, 80, draw_mean, component, df_reference)
stopifnot(identical(rownames(info_crp_pos$targets),
                    c("os::rx_x_tmb_braf", "os::cell_effect_10", "pfs::rx_x_tmb_braf", "pfs::cell_effect_10")))
m_bad <- matrix(0, 1, 9, dimnames = list("os::rx", paste0("os::", rownames(component$original_os$res))))
m_bad[1, "os::RxExperimental arm"] <- 1
expect_error_matching(design_information(design_tmb, 80, draw_mean, component, df_reference, targets = m_bad),
                      "does not identify: os::rx")
cat("5. enriched designs: rank 7 of 9, aliased pairs in the null space, identifiable targets only\n")

# ---------------------------------------------------------------------------
# 6. Scaling in n
# ---------------------------------------------------------------------------
design_all <- trial_design("all_comers", min_follow_up_weeks = 52, accrual_weeks = 104)
info_200 <- design_information(design_all, 200, draw_mean, component, df_reference)
info_400 <- design_information(design_all, 400, draw_mean, component, df_reference)
stopifnot(max(abs(info_400$endpoints$os$information - 2 * info_200$endpoints$os$information)) <
            1e-9 * max(abs(info_400$endpoints$os$information)),
          max(abs(info_200$covariance - 2 * info_400$covariance)) < 1e-9)
design_rate <- trial_design("all_comers_rate", min_follow_up_weeks = 52, accrual_rate = 2)
info_r100 <- design_information(design_rate, 100, draw_mean, component, df_reference)
info_r400 <- design_information(design_rate, 400, draw_mean, component, df_reference)
stopifnot(info_r100$accrual_weeks == 50, info_r400$accrual_weeks == 200,
          all(diag(info_r400$endpoints$os$information) / 400 > diag(info_r100$endpoints$os$information) / 100))
cat("6. information linear in n at fixed accrual; per-patient information grows with a fixed accrual rate\n")

# ---------------------------------------------------------------------------
# 7. Landmark randomisation reweights the cell mix
# ---------------------------------------------------------------------------
design_week4 <- trial_design("all_comers_week4", min_follow_up_weeks = 52, accrual_weeks = 104,
                             randomisation_week = 4)
df_mix0 <- design_patient_mix(design_all, df_reference, NULL, component, draw_mean)
df_mix4 <- design_patient_mix(design_week4, df_reference, NULL, component, draw_mean)
stopifnot(abs(sum(df_mix4$weight) - 1) < 1e-12, all(df_mix4$selection < 1), all(df_mix0$selection == 1))
v_ratio <- (df_mix4$weight / df_mix0$weight) / df_mix4$selection
stopifnot(max(v_ratio) - min(v_ratio) < 1e-10)
cat("7. randomisation at week 4 reweights the mix by P(PFS > 4) under control\n")

# ---------------------------------------------------------------------------
# 8. Cross-endpoint covariance
# ---------------------------------------------------------------------------
m_t <- info_400$targets
m_prior <- m_t[, colnames(component$joint_covariance$covariance)] %*% component$joint_covariance$covariance %*%
  t(m_t[, colnames(component$joint_covariance$covariance)])
i_os <- 1:3; i_pfs <- 4:6
whiten <- function(m) covariance_root(m[i_os, i_os], inverse = TRUE) %*% m[i_os, i_pfs] %*%
  covariance_root(m[i_pfs, i_pfs], inverse = TRUE)
stopifnot(max(abs(whiten(info_400$covariance) - whiten(m_prior))) < 1e-8)
info_ind <- design_information(design_all, 400, draw_mean, component, df_reference, cross_endpoint = "independent")
stopifnot(all(info_ind$covariance[i_os, i_pfs] == 0),
          max(abs(info_ind$covariance[i_os, i_os] - info_400$covariance[i_os, i_os])) < 1e-12)
cat("8. recoloured cross-endpoint covariance keeps the whitened prior correlation\n")

# ---------------------------------------------------------------------------
# 9. SLOW: V_d(n) against full simulation
# ---------------------------------------------------------------------------
if (identical(Sys.getenv("METIMMOX_SLOW_TESTS"), "1")) {
  df_slow <- data.frame()
  for (n in c(200, 400, 800)) {
    info_n <- design_information(design_week4, n, draw_mean, component, df_reference)
    for (outcome in c("os", "pfs")) {
      if (outcome == "os" && n != 400) next
      sim <- simulate_design_covariance(design_week4, n, draw_mean, component, df_reference,
                                        outcome = outcome, n_trials = 500, seed = 202L + n)
      v_targets <- rownames(sim$covariance)
      v_ratio <- diag(sim$covariance) / diag(info_n$covariance)[v_targets]
      tol <- 4 * sqrt(2 / sim$n_ok)
      df_slow <- rbind(df_slow, data.frame(n = n, outcome = outcome, target = v_targets, ratio = v_ratio,
                                           tolerance = tol, n_ok = sim$n_ok, n_failed = sim$n_failed))
    }
  }
  print(df_slow, row.names = FALSE)
  if (any(abs(df_slow$ratio - 1) > df_slow$tolerance)) stop("V_d(n) disagrees with the simulated covariance")
  cat("9. SLOW: simulated variances / V_d(n) within 4 SEs of 1 at n = 200, 400, 800\n")
} else {
  cat("9. SLOW full-simulation check skipped (set METIMMOX_SLOW_TESTS=1)\n")
}

cat("\nAll design-information tests passed.\n")
