# ===============================================================================
# TEST: EVSI by design (issue #202, paper 3 step 4)
# ===============================================================================
# Synthetic PSA outputs only (no trial data, no caches):
#
#   1. Normal-normal, linear NMB, one target: theta ~ N(mu, s0^2),
#      T ~ N(theta, sigma^2 / n); EVSI = E[max(0, b E(theta | T))] - max(0, b mu)
#      with E(theta | T) ~ N(mu, s0^4 / (s0^2 + sigma^2 / n)), i.e. the
#      expected value of the preposterior maximum. Both the GAM regression and
#      the Gaussian approximation agree within 3 MCSE + 3%.
#   2. Two targets, three policies, linear NMB: both methods against the
#      exact preposterior value (Monte Carlo over the normal preposterior
#      means, 2e6 draws), within 3 MCSE + 4%.
#   3. EVSI(0) = 0 ("zero_information").
#   4. EVSI is nondecreasing in n within two combined SEs and tends to the
#      EVPPI of the observed target (voi::evppi) for large n.
#   5. Zero EVPI gives EVSI exactly 0 with method "zero_evpi_bound" and SE 0
#      (#168); a small positive EVPI gives NA ("not_estimated").
#   6. Seeded reproducibility: the same seed reproduces the results exactly;
#      another seed changes them.
#   7. A per-draw covariance array with identical slices equals the common
#      matrix; a missing sample size stops.
#   8. Integration: evsi_design_from_information() on a synthetic joint gamma
#      component (all-comers design, six treatment targets) runs both methods
#      with finite results not above EVPI beyond Monte Carlo error.
#
# Run: "C:\Program Files\R\R-4.3.2\bin\x64\Rscript.exe" tests/test_evsi.R
# ===============================================================================

suppressPackageStartupMessages({
  library(survival)
  library(flexsurv)
})
source(here::here("R", "evppi_functions.R"))
source(here::here("R", "joint_survival_sampling.R"))
source(here::here("R", "pfs_os_violation_diagnostics.R"))
source(here::here("R", "voi_survival_parameters.R"))
source(here::here("R", "design_information.R"))
source(here::here("R", "evsi_functions.R"))

#' Stop unless `expr` raises an error whose message matches `pattern`.
expect_error_matching <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  if (!inherits(error, "error") || !grepl(pattern, conditionMessage(error))) {
    stop("Expected an error matching '", pattern, "', got: ",
         if (inherits(error, "error")) conditionMessage(error) else "no error")
  }
  invisible(TRUE)
}

#' Stop unless an EVSI row is within k MCSE plus a relative tolerance of `exact`.
expect_close <- function(df_rows, exact, k = 3, rel = 0.03, label = "") {
  for (i in seq_len(nrow(df_rows))) {
    tol <- k * df_rows$mcse[i] + rel * max(abs(exact[i]), 1)
    if (!is.finite(df_rows$evsi[i]) || abs(df_rows$evsi[i] - exact[i]) > tol) {
      stop(label, ": ", df_rows$method[i], " at n = ", df_rows$n[i], " gives ", signif(df_rows$evsi[i], 5),
           " against ", signif(exact[i], 5), " (tolerance ", signif(tol, 3), ")")
    }
  }
  invisible(TRUE)
}

#' One-target design with V = sigma^2 / n.
scalar_design <- function(n_grid, sigma = 1) {
  evsi_design("scalar", matrix(1, 1, 1, dimnames = list("os::x", "os::x")),
              stats::setNames(lapply(n_grid, function(n) matrix(sigma^2 / n)), n_grid))
}

#' E[max(0, X)] for X ~ N(mu, s^2).
expected_positive_part <- function(mu, s) mu * stats::pnorm(mu / s) + s * stats::dnorm(mu / s)

# ---------------------------------------------------------------------------
# 1. Normal-normal closed form, one target
# ---------------------------------------------------------------------------
set.seed(5)
n_sim <- 10000
v_theta <- stats::rnorm(n_sim, -0.1, 0.3)
m_nmb <- cbind(soc = 0, new = 1000 * (v_theta + stats::rnorm(n_sim, 0, 0.2)))
m_theta <- cbind("os::x" = v_theta)
v_n <- c(0, 2, 10, 50)
df_1 <- evsi_by_design(m_nmb, m_theta, scalar_design(v_n[-1]), v_n, B = 500, seed = 123L)
mu_hat <- mean(v_theta); s0_hat <- stats::sd(v_theta)
exact_1 <- vapply(v_n, function(n) if (n == 0) 0 else
  1000 * expected_positive_part(mu_hat, sqrt(s0_hat^4 / (s0_hat^2 + 1 / n))) - max(0, 1000 * mu_hat), numeric(1))
for (method in EVSI_METHODS) {
  v_rows <- df_1$method == method
  expect_close(df_1[v_rows, ], exact_1[match(df_1$n[v_rows], v_n)], label = "closed form")
}
print(cbind(df_1[, c("n", "method", "evsi", "mcse")], exact = exact_1[match(df_1$n, v_n)]), row.names = FALSE)
cat("1. one-target normal-normal closed form reproduced by both methods (3 MCSE + 3%)\n")

# ---------------------------------------------------------------------------
# 2. Two targets, three policies, linear NMB
# ---------------------------------------------------------------------------
set.seed(9)
m_sigma <- matrix(c(0.09, 0.03, 0.03, 0.06), 2)
m_th2 <- mvtnorm::rmvnorm(5000, c(-0.05, 0.02), m_sigma)
colnames(m_th2) <- c("os::a", "pfs::b")
m_b <- rbind(c(0, 0), c(800, 200), c(-300, 900))
v_c <- c(0, 10, -15)
m_nmb2 <- sweep(m_th2 %*% t(m_b), 2, v_c, "+") + matrix(stats::rnorm(5000 * 3, 0, 50), 5000)
m_v_unit <- diag(c(0.5, 0.8))
m_id2 <- diag(2); dimnames(m_id2) <- list(colnames(m_th2), colnames(m_th2))
design_2 <- evsi_design("two_targets", m_id2,
                        list("5" = m_v_unit / 5, "40" = m_v_unit / 40))
df_2 <- evsi_by_design(m_nmb2, m_th2, design_2, c(5, 40), B = 500)
v_mu2 <- colMeans(m_th2); m_s2 <- stats::cov(m_th2)
exact_2 <- vapply(c(5, 40), function(n) {
  m_r <- m_s2 %*% solve(m_s2 + m_v_unit / n) %*% m_s2
  set.seed(1)
  m_post <- mvtnorm::rmvnorm(2e6, v_mu2, m_r)
  m_val <- sweep(m_post %*% t(m_b), 2, v_c, "+")
  mean(apply(m_val, 1, max)) - max(v_c + as.vector(m_b %*% v_mu2))
}, numeric(1))
for (method in EVSI_METHODS) {
  v_rows <- df_2$method == method
  expect_close(df_2[v_rows, ], exact_2, rel = 0.04, label = "two targets")
}
print(cbind(df_2[, c("n", "method", "evsi", "mcse")], exact = rep(exact_2, each = 2)), row.names = FALSE)
cat("2. two-target, three-policy linear example reproduced by both methods (3 MCSE + 4%)\n")

# ---------------------------------------------------------------------------
# 3. EVSI(0) = 0
# ---------------------------------------------------------------------------
v_zero <- df_1[df_1$n == 0, ]
stopifnot(nrow(v_zero) == 1, v_zero$evsi == 0, v_zero$mcse == 0, v_zero$method == "zero_information")
cat("3. EVSI(0) = 0\n")

# ---------------------------------------------------------------------------
# 4. Monotone in n and tending to the EVPPI of the observed target
# ---------------------------------------------------------------------------
v_n4 <- c(1, 4, 16, 64, 256, 1e5)
df_4 <- evsi_by_design(m_nmb[1:5000, ], m_theta[1:5000, , drop = FALSE], scalar_design(v_n4), v_n4,
                       methods = "gam_regression", B = 500)
v_drop <- df_4$evsi[-1] - df_4$evsi[-nrow(df_4)]
v_tol <- 2 * sqrt(df_4$mcse[-1]^2 + df_4$mcse[-nrow(df_4)]^2)
stopifnot(all(v_drop >= -v_tol))
set.seed(123L)
evppi <- voi::evppi(m_nmb[1:5000, ], data.frame(x = v_theta[1:5000]), pars = "x", se = TRUE, B = 500)
last <- df_4[nrow(df_4), ]
stopifnot(abs(last$evsi - evppi$evppi) <= 3 * sqrt(last$mcse^2 + evppi$se^2))
cat(sprintf("4. EVSI nondecreasing in n; EVSI(1e5) = %.2f against EVPPI %.2f\n", last$evsi, evppi$evppi))

# ---------------------------------------------------------------------------
# 5. Zero and small EVPI
# ---------------------------------------------------------------------------
m_nmb_dom <- cbind(soc = 1000 + abs(m_nmb[, 2]), new = m_nmb[, 2])
df_5 <- evsi_by_design(m_nmb_dom, m_theta, scalar_design(c(2, 10)), c(0, 2, 10), B = 100)
stopifnot(all(df_5$evsi == 0), all(df_5$mcse == 0), all(df_5$method == "zero_evpi_bound"), all(df_5$evpi == 0))
m_nmb_small <- cbind(soc = rep(0, 100), new = c(0.5, rep(-1e3, 99)))
df_5b <- evsi_by_design(m_nmb_small, m_theta[1:100, , drop = FALSE], scalar_design(2), 2, B = 100)
stopifnot(is.na(df_5b$evsi), df_5b$method == "not_estimated", nzchar(df_5b$error))
cat("5. zero EVPI gives EVSI 0 (zero_evpi_bound); EVPI <= 0.01 is not estimated (NA)\n")

# ---------------------------------------------------------------------------
# 6. Seeded reproducibility
# ---------------------------------------------------------------------------
df_6a <- evsi_by_design(m_nmb[1:3000, ], m_theta[1:3000, , drop = FALSE], scalar_design(10), 10, B = 200, seed = 42L)
df_6b <- evsi_by_design(m_nmb[1:3000, ], m_theta[1:3000, , drop = FALSE], scalar_design(10), 10, B = 200, seed = 42L)
df_6c <- evsi_by_design(m_nmb[1:3000, ], m_theta[1:3000, , drop = FALSE], scalar_design(10), 10, B = 200, seed = 43L)
stopifnot(identical(df_6a, df_6b),
          df_6a$evsi[df_6a$method == "gam_regression"] != df_6c$evsi[df_6c$method == "gam_regression"])
cat("6. identical results under one seed; another seed changes them\n")

# ---------------------------------------------------------------------------
# 7. Per-draw covariance arrays; missing sample sizes
# ---------------------------------------------------------------------------
a_cov <- array(0.1, dim = c(3000, 1, 1))
design_array <- evsi_design("array", matrix(1, 1, 1, dimnames = list("os::x", "os::x")), list("10" = a_cov))
df_7 <- evsi_by_design(m_nmb[1:3000, ], m_theta[1:3000, , drop = FALSE], design_array, 10, B = 200, seed = 42L)
stopifnot(max(abs(df_7$evsi - df_6a$evsi)) < 1e-8, max(abs(df_7$mcse - df_6a$mcse)) < 1e-8)
expect_error_matching(evsi_by_design(m_nmb[1:3000, ], m_theta[1:3000, , drop = FALSE], scalar_design(10), 20, B = 50),
                      "no covariance for n = 20")
cat("7. per-draw covariance array with identical slices equals the common matrix\n")

# ---------------------------------------------------------------------------
# 8. Integration with design_information()
# ---------------------------------------------------------------------------
set.seed(11)
n_pat <- 300
df_pat <- data.frame(Age = stats::rnorm(n_pat, 62, 9), sex = factor(stats::rbinom(n_pat, 1, 0.45)),
                     Rx = factor(stats::rbinom(n_pat, 1, 0.5), labels = c("Control", "Experimental arm")),
                     crp = factor(stats::rbinom(n_pat, 1, 0.35)), tmb_braf = factor(stats::rbinom(n_pat, 1, 0.45)))
v_eta <- log(0.02) + 0.01 * (df_pat$Age - 62) - 0.2 * (df_pat$Rx != "Control") - 0.3 * (df_pat$crp == "1") *
  (df_pat$Rx != "Control")
v_t <- stats::rgamma(n_pat, 1.6, exp(v_eta)); v_cens <- stats::runif(n_pat, 30, 150)
df_pat$OSwk <- pmin(v_t, v_cens); df_pat$Death <- as.numeric(v_t <= v_cens)
v_t2 <- pmin(v_t, stats::rgamma(n_pat, 1.2, exp(v_eta + 0.8)))
df_pat$PFSwk <- pmin(v_t2, v_cens); df_pat$Prog <- as.numeric(v_t2 <= v_cens)
fit_os <- flexsurvreg(Surv(OSwk, Death) ~ Age + sex + Rx + crp * Rx + tmb_braf * Rx, data = df_pat, dist = "gamma")
fit_pfs <- flexsurvreg(Surv(PFSwk, Prog) ~ Age + sex + Rx + crp * Rx + tmb_braf * Rx, data = df_pat, dist = "gamma")
m_joint <- as.matrix(Matrix::bdiag(fit_os$cov, fit_pfs$cov))
v_joint_names <- c(paste0("os::", rownames(fit_os$res)), paste0("pfs::", rownames(fit_pfs$res)))
dimnames(m_joint) <- list(v_joint_names, v_joint_names)
component <- list(method = "mvn_joint_v2", joint_covariance = list(covariance = m_joint),
                  draws = draw_joint_survival_coefficients(fit_os, fit_pfs, m_joint, 3000, 11L),
                  original_os = fit_os, original_pfs = fit_pfs, n_samples = 3000L, dist_os = "gamma",
                  dist_pfs = "gamma", seed = 11L)
design_s1 <- trial_design("S1", min_follow_up_weeks = 52, accrual_weeks = 104, randomisation_week = 4)
design_evsi <- evsi_design_from_information(design_s1, c(100, 400), component_mean_draw(component), component,
                                            df_pat[, c("Age", "sex", "Rx", "crp", "tmb_braf")])
stopifnot(nrow(design_evsi$targets) == 6, length(design_evsi$covariances) == 2)
m_draws <- cbind(component$draws$os, component$draws$pfs)
colnames(m_draws) <- v_joint_names
v_cell_crp <- m_draws[, "pfs::RxExperimental arm"] + m_draws[, "pfs::RxExperimental arm:crp1"]
m_nmb8 <- cbind(soc = 0, crp_guided = -4000 * v_cell_crp - 600 + stats::rnorm(3000, 0, 100))
df_8 <- evsi_by_design(m_nmb8, m_draws, design_evsi, c(0, 100, 400), B = 200)
v_pos <- df_8$n > 0
stopifnot(all(is.finite(df_8$evsi)), all(df_8$evsi[v_pos] <= df_8$evpi[v_pos] + 3 * df_8$mcse[v_pos]),
          all(df_8$evsi[v_pos] >= -3 * df_8$mcse[v_pos]))
print(df_8[, c("design", "n", "method", "evsi", "mcse", "evpi")], row.names = FALSE)
cat("8. evsi_design_from_information() and evsi_by_design() run on a synthetic joint gamma component\n")

cat("\nAll EVSI tests passed.\n")
