# Read-only diagnostics for issue #179; no production cache writer is called.

population_reference_survival <- function(age, sex, life_table, years) {
  stopifnot(length(age) == length(sex), all(is.finite(age)),
            all(sex %in% c("0", "1")), all(years >= 0),
            all(is.finite(life_table$qx)),
            all(life_table$qx >= 0 & life_table$qx <= 1))
  vapply(seq_along(age), function(i) {
    tab <- life_table[life_table$sex == if (sex[i] == "0") "female" else "male", ]
    tab <- tab[order(tab$age), ]
    stopifnot(!anyDuplicated(tab$age), all(diff(tab$age) == 1), min(age[i]) >= min(tab$age),
              max(age[i] + years) < max(tab$age) + 1)
    # Integrate only intervals actually reached (avoids 0 * Inf at qx = 1).
    vapply(years, function(t) {
      duration <- pmax(0, pmin(age[i] + t, tab$age + 1) - pmax(age[i], tab$age))
      used <- duration > 0
      exp(sum(log1p(-tab$qx[used]) * duration[used]))
    }, numeric(1))
  }, numeric(length(years)))
}

extrapolation_tail_summary <- function(rmst, qaly, fractions = c(.01, .05, .10)) {
  stopifnot(length(rmst) == length(qaly), length(qaly) > 1,
            all(is.finite(rmst)), all(is.finite(qaly)), sum(qaly) > 0,
            all(fractions > 0 & fractions < 1))
  order_id <- order(rmst, decreasing = TRUE, method = "radix")
  do.call(rbind, lapply(fractions, function(f) {
    k <- ceiling(length(qaly) * f)
    stopifnot(k < length(qaly))
    top <- order_id[seq_len(k)]
    data.frame(fraction = f, removed = k, cutoff = min(rmst[top]),
      contribution = sum(qaly[top]) / sum(qaly), mean = mean(qaly),
      trimmed = mean(qaly[-top]), shift = mean(qaly) - mean(qaly[-top]))
  }))
}

read_extrapolation_ordering <- function(path, sampling, population, params, base) {
  if (!file.exists(path)) stop("Ordering diagnostic cache is required; rendering cannot rebuild it.")
  x <- readRDS(path)
  expected <- list(sampling_fingerprint = sampling$fingerprint,
    n_draws = as.integer(nrow(params)), time_horizon = base$time_horizon, cl = base$cl,
    u_np = base$u_np, u_p = base$u_p, dr_effects = base$dr_effects,
    prediction_population = cache_fingerprint(population),
    psa_parameters = cache_fingerprint(params), version = "pfs_os_violations_v3_trapezoid")
  if (!identical(x$fingerprint, cache_fingerprint(expected)))
    stop("Ordering diagnostic cache is stale; no cache was changed.")
  x
}

reconstruct_extrapolation_curves <- function(sampling, psa, draws, base, ordering) {
  population <- base$prediction_population
  component <- get_joint_sampling_models(sampling)
  times <- seq(0, base$time_horizon)
  strategies <- psa$strategies
  n <- psa$n_sim
  curves <- setNames(lapply(strategies, function(s) list(
    os = matrix(NA_real_, n, length(times)),
    pfs = matrix(NA_real_, n, length(times)))), strategies)
  crossing <- matrix(FALSE, n, length(strategies), dimnames = list(NULL, strategies))
  # Include sim as well as model_idx: replacement models may recur under
  # different population weights. A model-index-only join would multiply rows.
  pc <- ordering$per_curve
  key <- paste(pc$draw, pc$model_idx, pc$curve)
  stopifnot(!anyDuplicated(key))
  max_qaly_error <- 0
  for (i in seq_len(n)) {
    p <- base
    for (nm in intersect(names(draws), names(p))) p[[nm]] <- draws[[nm]][i]
    weights <- model_population_weights(p, population)
    sg <- psa_violation_subgroups(population, weights = weights)
    raw <- psa_draw_subgroup_curves(component, psa$model_idx[i], sg, times)
    get_cross <- function(cn) {
      j <- match(paste(psa$sim[i], psa$model_idx[i], cn), key)
      if (is.na(j)) stop("Missing sim/model_idx/curve ordering diagnostic.")
      flag <- any(raw[[cn]]$pfs > raw[[cn]]$os)
      if (!identical(flag, pc$violated[j])) stop("Ordering flag disagrees with reconstructed curves.")
      flag
    }
    for (s in strategies) {
      cn <- if (s == "control") s else paste0(s, c("_pos", "_neg"))
      mix <- if (s == "control") 1 else {
        pos <- sum(weights[as.character(population[[s]]) == "1"])
        c(pos, 1 - pos)
      }
      os <- Reduce(`+`, Map(function(k, w) raw[[k]]$os * w, cn, mix))
      pfs <- Reduce(`+`, Map(function(k, w) pmin(raw[[k]]$pfs, raw[[k]]$os) * w, cn, mix))
      curves[[s]]$os[i, ] <- os
      curves[[s]]$pfs[i, ] <- pfs
      crossing[i, s] <- any(vapply(cn, get_cross, logical(1)))
      q <- restricted_mean_survival((pfs * p$u_np + (os - pfs) * p$u_p) /
        (1 + p$dr_effects)^(times * p$cl), p$cl)
      max_qaly_error <- max(max_qaly_error, abs(q - psa$effectiveness[i, s]))
    }
    if (i %% 500 == 0) message("Reconstructed ", i, "/", n, " PSA rows")
  }
  if (!is.finite(max_qaly_error) || max_qaly_error > 1e-10)
    stop("Reconstructed curves do not reproduce cached PSA QALYs: ", max_qaly_error)
  list(curves = curves, crossing = crossing, max_qaly_error = max_qaly_error)
}
