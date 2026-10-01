# Joint OS/PFS normal approximation (issue #159).
# Paired patient bootstraps estimate dependence only. Block whitening and
# recolouring retain each fitted model's observed-information covariance.

covariance_root <- function(x, inverse = FALSE) {
  if (!is.matrix(x) || nrow(x) != ncol(x) || any(!is.finite(x)))
    stop("A finite square covariance matrix is required")
  eig <- eigen((x + t(x)) / 2, symmetric = TRUE)
  if (min(eig$values) <= max(eig$values) * 1e-12)
    stop("Covariance matrix is not positive definite")
  power <- if (inverse) -0.5 else 0.5
  eig$vectors %*% (eig$values^power * t(eig$vectors))
}

calibrate_joint_covariance <- function(bootstrap_os, bootstrap_pfs, cov_os, cov_pfs) {
  if (nrow(bootstrap_os) != nrow(bootstrap_pfs) ||
      nrow(bootstrap_os) <= ncol(bootstrap_os) + ncol(bootstrap_pfs))
    stop("Too few paired bootstrap fits for joint covariance estimation")
  if (ncol(bootstrap_os) != nrow(cov_os) || ncol(bootstrap_pfs) != nrow(cov_pfs))
    stop("Bootstrap and fitted covariance dimensions differ")
  m_bootstrap_cov <- stats::cov(cbind(bootstrap_os, bootstrap_pfs))
  ios <- seq_len(ncol(bootstrap_os))
  ipfs <- length(ios) + seq_len(ncol(bootstrap_pfs))
  # If B is the paired bootstrap covariance, the cross-block is
  # V_os^(1/2) B_os^(-1/2) B_os,pfs B_pfs^(-1/2) V_pfs^(1/2).
  # Congruence with a block-diagonal transform preserves positive
  # semidefiniteness; no clipping of correlations or marginal variances.
  m_cross <- covariance_root(cov_os) %*%
    covariance_root(m_bootstrap_cov[ios, ios, drop = FALSE], inverse = TRUE) %*%
    m_bootstrap_cov[ios, ipfs, drop = FALSE] %*%
    covariance_root(m_bootstrap_cov[ipfs, ipfs, drop = FALSE], inverse = TRUE) %*%
    covariance_root(cov_pfs)
  m_joint <- rbind(cbind(cov_os, m_cross), cbind(t(m_cross), cov_pfs))
  nm <- c(paste0("os::", colnames(bootstrap_os)),
          paste0("pfs::", colnames(bootstrap_pfs)))
  dimnames(m_joint) <- list(nm, nm)
  if (min(eigen(m_joint, symmetric = TRUE, only.values = TRUE)$values) < -1e-10)
    stop("Estimated joint covariance is not positive semidefinite")
  list(covariance = m_joint, bootstrap_covariance = m_bootstrap_cov)
}

estimate_joint_survival_covariance <- function(formula_os, formula_pfs, data,
                                                original_os, original_pfs,
                                                dist_os, dist_pfs,
                                                n_bootstrap = 1000L, seed = 123L,
                                                max_failure_rate = 0.1) {
  fits <- list(os = original_os, pfs = original_pfs)
  par_names <- lapply(fits, function(fit) rownames(fit$res)[fit$optpars])
  for (fit in fits) {
    if (is.null(fit$cov)) stop("Fitted model covariance is unavailable")
    covariance_root(fit$cov)
  }
  if (length(n_bootstrap) != 1L || !is.finite(n_bootstrap) ||
      n_bootstrap != as.integer(n_bootstrap) || n_bootstrap <= sum(lengths(par_names)))
    stop("n_bootstrap must exceed the total number of estimated coefficients")
  boot <- lapply(par_names, function(nm) matrix(NA_real_, n_bootstrap, length(nm),
                                              dimnames = list(NULL, nm)))
  errors <- rep(NA_character_, n_bootstrap)
  set.seed(seed)
  for (b in seq_len(n_bootstrap)) {
    # One patient index vector serves BOTH endpoints; endpoint resampling
    # separately would erase the cross-covariance this procedure estimates.
    v_rows <- sample.int(nrow(data), nrow(data), replace = TRUE)
    pair <- tryCatch({
      os <- suppressWarnings(flexsurv::flexsurvreg(formula_os, data = data[v_rows, ], dist = dist_os))
      pfs <- suppressWarnings(flexsurv::flexsurvreg(formula_pfs, data = data[v_rows, ], dist = dist_pfs))
      candidate <- list(os = os, pfs = pfs)
      for (outcome in names(candidate)) {
        fit <- candidate[[outcome]]
        if (!identical(fit$opt$convergence, 0L) || any(!is.finite(fit$opt$par)) ||
            !identical(rownames(fit$res)[fit$optpars], par_names[[outcome]]) ||
            is.null(fit$cov)) stop(outcome, " fit failed convergence/parameter checks")
        covariance_root(fit$cov)
      }
      lapply(candidate, function(fit) fit$opt$par)
    }, error = function(e) {
      errors[b] <<- conditionMessage(e)
      NULL
    })
    if (!is.null(pair)) for (outcome in names(pair)) boot[[outcome]][b, ] <- pair[[outcome]]
    if (b %% 100L == 0L)
      cat(sprintf("  Paired covariance bootstrap: %d/%d attempted; %d failed pairs\n",
                  b, n_bootstrap, sum(!is.na(errors[seq_len(b)]))))
  }
  v_valid <- is.na(errors)
  if (mean(!v_valid) > max_failure_rate)
    stop(sprintf("Paired covariance bootstrap failed for %d/%d pairs (limit %.0f%%)",
                 sum(!v_valid), n_bootstrap, 100 * max_failure_rate))
  boot <- lapply(boot, function(x) x[v_valid, , drop = FALSE])
  estimate <- calibrate_joint_covariance(boot$os, boot$pfs, original_os$cov, original_pfs$cov)
  c(estimate, list(bootstrap_draws = boot, n_attempted = n_bootstrap,
                  n_successful = sum(v_valid), n_failed = sum(!v_valid),
                  failure_reasons = table(errors, useNA = "no"),
                  seed = seed, resampling = "paired patients, unstratified",
                  calibration = "block whitening; fitted marginal covariances"))
}

draw_joint_survival_coefficients <- function(original_os, original_pfs, covariance,
                                              n_samples, seed) {
  fits <- list(os = original_os, pfs = original_pfs)
  v_mean <- unlist(lapply(fits, function(fit) fit$opt$par), use.names = FALSE)
  if (!identical(dim(covariance), c(length(v_mean), length(v_mean))) ||
      any(!is.finite(covariance))) stop("Invalid joint coefficient covariance")
  # Separate reproducible draw stream: bootstrap RNG consumption never alters
  # the normal variates. A zero covariance returns the fitted coefficients.
  set.seed(seed)
  m_joint <- mvtnorm::rmvnorm(n_samples, v_mean, covariance)
  offset <- 0L
  lapply(fits, function(fit) {
    nm <- rownames(fit$res)
    m_draws <- matrix(rep(fit$res.t[, "est"], each = n_samples), nrow = n_samples,
                    dimnames = list(NULL, nm))
    idx <- offset + seq_along(fit$optpars)
    m_draws[, fit$optpars] <- m_joint[, idx, drop = FALSE]
    offset <<- offset + length(fit$optpars)
    m_draws
  })
}
