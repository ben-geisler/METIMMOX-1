# Read-only diagnostics for issue #179; no production cache writer is called.

population_reference_survival <- function(age, sex, life_table, years) {
  stopifnot(length(age) == length(sex), all(is.finite(age)),
            all(sex %in% c("0", "1")), all(years >= 0),
            all(is.finite(life_table$qx)),
            all(life_table$qx >= 0 & life_table$qx <= 1))
  vapply(seq_along(age), function(i) {
    df_tab <- life_table[life_table$sex == if (sex[i] == "0") "female" else "male", ]
    df_tab <- df_tab[order(df_tab$age), ]
    stopifnot(!anyDuplicated(df_tab$age), all(diff(df_tab$age) == 1), min(age[i]) >= min(df_tab$age),
              max(age[i] + years) < max(df_tab$age) + 1)
    # Integrate only intervals actually reached (avoids 0 * Inf at qx = 1).
    vapply(years, function(t) {
      v_duration <- pmax(0, pmin(age[i] + t, df_tab$age + 1) - pmax(age[i], df_tab$age))
      v_used <- v_duration > 0
      exp(sum(log1p(-df_tab$qx[v_used]) * v_duration[v_used]))
    }, numeric(1))
  }, numeric(length(years)))
}

extrapolation_tail_summary <- function(rmst, qaly, fractions = c(.01, .05, .10)) {
  stopifnot(length(rmst) == length(qaly), length(qaly) > 1,
            all(is.finite(rmst)), all(is.finite(qaly)), sum(qaly) > 0,
            all(fractions > 0 & fractions < 1))
  v_order_id <- order(rmst, decreasing = TRUE, method = "radix")
  do.call(rbind, lapply(fractions, function(f) {
    k <- ceiling(length(qaly) * f)
    stopifnot(k < length(qaly))
    v_top <- v_order_id[seq_len(k)]
    data.frame(fraction = f, removed = k, cutoff = min(rmst[v_top]),
      contribution = sum(qaly[v_top]) / sum(qaly), mean = mean(qaly),
      trimmed = mean(qaly[-v_top]), shift = mean(qaly) - mean(qaly[-v_top]))
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
  v_times <- seq(0, base$time_horizon)
  v_strategies <- psa$strategies
  n <- psa$n_sim
  curves <- setNames(lapply(v_strategies, function(s) list(
    os = matrix(NA_real_, n, length(v_times)),
    pfs = matrix(NA_real_, n, length(v_times)))), v_strategies)
  m_crossing <- matrix(FALSE, n, length(v_strategies), dimnames = list(NULL, v_strategies))
  # Include sim as well as model_idx: replacement models may recur under
  # different population weights. A model-index-only join would multiply rows.
  df_pc <- ordering$per_curve
  v_key <- paste(df_pc$draw, df_pc$model_idx, df_pc$curve)
  stopifnot(!anyDuplicated(v_key))
  max_qaly_error <- 0
  for (i in seq_len(n)) {
    p <- base
    for (nm in intersect(names(draws), names(p))) p[[nm]] <- draws[[nm]][i]
    v_weights <- model_population_weights(p, population)
    sg <- psa_violation_subgroups(population, weights = v_weights)
    raw <- psa_draw_subgroup_curves(component, psa$model_idx[i], sg, v_times)
    get_cross <- function(cn) {
      j <- match(paste(psa$sim[i], psa$model_idx[i], cn), v_key)
      if (is.na(j)) stop("Missing sim/model_idx/curve ordering diagnostic.")
      flag <- any(raw[[cn]]$pfs > raw[[cn]]$os)
      if (!identical(flag, df_pc$violated[j])) stop("Ordering flag disagrees with reconstructed curves.")
      flag
    }
    for (s in v_strategies) {
      cn <- if (s == "control") s else paste0(s, c("_pos", "_neg"))
      v_mix <- if (s == "control") 1 else {
        pos <- sum(v_weights[as.character(population[[s]]) == "1"])
        c(pos, 1 - pos)
      }
      v_os <- Reduce(`+`, Map(function(k, w) raw[[k]]$os * w, cn, v_mix))
      v_pfs <- Reduce(`+`, Map(function(k, w) pmin(raw[[k]]$pfs, raw[[k]]$os) * w, cn, v_mix))
      curves[[s]]$os[i, ] <- v_os
      curves[[s]]$pfs[i, ] <- v_pfs
      m_crossing[i, s] <- any(vapply(cn, get_cross, logical(1)))
      q <- restricted_mean_survival((v_pfs * p$u_np + (v_os - v_pfs) * p$u_p) /
        (1 + p$dr_effects)^(v_times * p$cl), p$cl)
      max_qaly_error <- max(max_qaly_error, abs(q - psa$effectiveness[i, s]))
    }
    if (i %% 500 == 0) message("Reconstructed ", i, "/", n, " PSA rows")
  }
  if (!is.finite(max_qaly_error) || max_qaly_error > 1e-10)
    stop("Reconstructed curves do not reproduce cached PSA QALYs: ", max_qaly_error)
  list(curves = curves, crossing = m_crossing, max_qaly_error = max_qaly_error)
}
