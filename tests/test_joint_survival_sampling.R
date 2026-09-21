# Synthetic covariance and paired-refit contracts for issue #159; no trial data.
library(survival)
source("R/joint_survival_sampling.R")
set.seed(812)
n <- 4000L
x <- cbind(a = rnorm(n), b = rnorm(n))
y <- cbind(c = 0.8 * x[, 1] + 0.6 * rnorm(n), d = -0.5 * x[, 2] + rnorm(n))
v1 <- matrix(c(1, 0.2, 0.2, 2), 2)
v2 <- matrix(c(3, -0.4, -0.4, 1), 2)
estimate <- calibrate_joint_covariance(x, y, v1, v2)$covariance
stopifnot(max(abs(estimate[1:2, 1:2] - v1)) < 1e-12,
          max(abs(estimate[3:4, 3:4] - v2)) < 1e-12,
          min(eigen(estimate, symmetric = TRUE)$values) > 0,
          estimate[1, 3] > 1, estimate[2, 4] < -0.4)
independent <- calibrate_joint_covariance(x, y[sample.int(n), ], v1, v2)$covariance
stopifnot(max(abs(independent[1:2, 3:4])) < 0.1)
cat("PASS: calibrated cross-covariance preserves fitted marginal blocks and dependence.\n")

# Identical endpoints must use identical patient resamples. Separate endpoint
# resampling would destroy the equality of bootstrap estimates in this probe.
set.seed(81)
d <- data.frame(age = rnorm(150))
d$time <- rexp(150, exp(-1 + 0.25 * d$age))
d$event <- as.integer(d$time <= 5)
d$time <- pmin(d$time, 5)
f <- Surv(time, event) ~ age
fit <- flexsurv::flexsurvreg(f, data = d, dist = "exp")
get_cov <- function(seed) estimate_joint_survival_covariance(f, f, d, fit, fit,
  "exp", "exp", n_bootstrap = 40L, seed = seed)
a <- get_cov(45L)
invisible(runif(30))
b <- get_cov(45L)
c <- get_cov(46L)
stopifnot(identical(a$bootstrap_draws, b$bootstrap_draws),
          !identical(a$bootstrap_draws, c$bootstrap_draws),
          identical(a$bootstrap_draws$os, a$bootstrap_draws$pfs),
          a$n_successful + a$n_failed == 40L,
          max(abs(a$covariance[1:2, 3:4] - fit$cov)) < 1e-10)
draws <- draw_joint_survival_coefficients(fit, fit, a$covariance, 10000L, 123L)
stopifnot(max(abs(draws$os - draws$pfs)) < 1e-6,
          max(abs(cov(draws$os) - fit$cov)) < 0.002,
          identical(draws, draw_joint_survival_coefficients(fit, fit, a$covariance, 10000L, 123L)))
bad <- try(calibrate_joint_covariance(x[1:3, ], y[1:3, ], v1, v2), silent = TRUE)
stopifnot(inherits(bad, "try-error"))
bad <- try(covariance_root(matrix(c(1, 2, 2, 1), 2)), silent = TRUE)
stopifnot(inherits(bad, "try-error"))
cat("PASS: paired refits, joint draws, RNG reproducibility and invalid covariance guards.\n")
