# ===============================================================================
# PFS > OS ordering violations in the PSA: per-draw diagnostics
#
# The PSA draws the OS and PFS coefficient vectors independently (issue #156),
# so a sampled pair can put PFS above OS on part of the weekly grid. model_fun()
# then clamps PFS to OS (pmin) for that curve and emits a
# `survival_ordering_warning`; 10_PSA.R only counts those warnings, it does not
# keep the curves. The functions here re-predict the five subgroup curves of
# every draw from the sampling cache and quantify the violation before and
# after the clamp:
#
#   * n_violations, first/last violating week, max excess (PFS - OS);
#   * excess area  = integral of max(PFS - OS, 0) over the horizon, in
#     patient-years (trapezoidal rule on the weekly grid). After the clamp it
#     is zero by construction, so "after" is reported as the difference the
#     clamp makes to the progressed-state area;
#   * progressed-state area = integral of (OS - PFS) raw (can be negative) and
#     of max(OS - PFS, 0) clamped, in patient-years;
#   * QALY impact = discounted QALYs the raw curves would have produced minus
#     the clamped QALYs. Within the model's cycle sum this is exactly
#     sum_t (PFS_t - OS_t)^+ (u_np - u_p) cl w_t, because the raw curves give
#     p_pf = PFS and p_p = OS - PFS < 0 on the violating cycles.
#
# The curves are the same population-averaged predictions as
# generate_psa_population_averaged_predictions() (script 06): control = every
# complete-case patient with Rx = control; biomarker-positive = the positives
# with Rx = experimental; biomarker-negative = the negatives with Rx = control.
# They are computed in closed form from the gamma AFT parameters instead of
# through flexsurv::predict(), which is far slower; the closed form is
# checked against the pipeline helper by validate_direct_curves().
# ===============================================================================

#' Design matrix for a fitted flexsurvreg model on new data (no intercept)
#'
#' @param model A flexsurvreg object.
#' @param newdata Data frame with the model covariates.
#' @return Numeric matrix whose columns are in the model's covariate order
#'   (`rownames(model$res)[model$covpars]`).
flexsurv_design_matrix <- function(model, newdata) {
  tt <- delete.response(terms(model$all.formulae[[model$dlist$location]]))
  mf <- model.frame(tt, newdata, xlev = model$xlev, na.action = na.pass)
  X <- model.matrix(tt, mf, contrasts.arg = model$contrasts)
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  cov_names <- rownames(model$res)[model$covpars]
  if (!setequal(colnames(X), cov_names)) {
    stop("Design matrix columns (", paste(colnames(X), collapse = ", "),
         ") do not match model covariates (", paste(cov_names, collapse = ", "), ")")
  }
  X[, cov_names, drop = FALSE]
}

#' Population-averaged gamma survival curve for one coefficient draw
#'
#' @param model Original flexsurvreg gamma fit (supplies formula and levels).
#' @param draw Named coefficient vector on flexsurvreg's optimisation scale
#'   (log shape, log rate, covariate effects on log rate).
#' @param newdata Patients to predict for (their `Rx` already set).
#' @param times Time points (weeks).
#' @return Numeric vector of length `length(times)`: mean survival over patients.
direct_gamma_pop_avg <- function(model, draw, newdata, times, weights = rep(1, nrow(newdata))) {
  if (model$dlist$name != "gamma") {
    stop("direct_gamma_pop_avg() supports the gamma distribution only; got ",
         model$dlist$name)
  }
  X <- flexsurv_design_matrix(model, newdata)
  keep <- stats::complete.cases(X)
  if (!all(keep)) stop("The explicit prediction population contains incomplete predictors.")
  shape <- exp(draw[["shape"]])
  eta <- as.vector(draw[["rate"]] + X %*% draw[colnames(X)])
  rate <- exp(eta)
  # S_i(t) = 1 - F_gamma(t; shape, rate_i); average over patients.
  surv <- outer(rate, times, function(r, t) {
    stats::pgamma(t, shape = shape, rate = r, lower.tail = FALSE)
  })
  weighted_survival_average(t(surv), weights)
}

#' Subgroup definitions mirroring model_fun() / script 06
#'
#' @param data_complete Complete-case cohort.
#' @param biomarkers Economic biomarker names.
#' @return Named list of `list(data, rx, label, strategy, subgroup)`.
psa_violation_subgroups <- function(data_complete, biomarkers = get_biomarkers(),
                                    weights = rep(1 / nrow(data_complete), nrow(data_complete))) {
  rx_levels <- levels(data_complete$Rx)
  control_strategy <- get_control_strategy()
  out <- list()
  out[[control_strategy]] <- list(
    data = data_complete, weights = weights, rx = rx_levels[1],
    label = "Standard of care (control)", strategy = control_strategy,
    subgroup = "control"
  )
  for (bm in biomarkers) {
    status <- as.numeric(as.character(data_complete[[bm]]))
    disp <- unname(strategy_display_name(bm))
    out[[paste0(bm, "_pos")]] <- list(
      data = data_complete[status == 1, , drop = FALSE], weights = weights[status == 1], rx = rx_levels[2],
      label = paste0(disp, " positive (experimental)"), strategy = bm,
      subgroup = "positive"
    )
    out[[paste0(bm, "_neg")]] <- list(
      data = data_complete[status == 0, , drop = FALSE], weights = weights[status == 0], rx = rx_levels[1],
      label = paste0(disp, " negative (control)"), strategy = bm,
      subgroup = "negative"
    )
  }
  out
}

#' Predict the OS and PFS curves of every subgroup for one draw
#'
#' @param component Joint sampling component (`get_joint_sampling_models()`).
#' @param idx Draw index.
#' @param subgroups Output of `psa_violation_subgroups()`.
#' @param times Weekly grid.
#' @return List (per subgroup) of `list(os, pfs)`.
psa_draw_subgroup_curves <- function(component, idx, subgroups, times) {
  os_draw <- component$draws$os[idx, ]
  pfs_draw <- component$draws$pfs[idx, ]
  lapply(subgroups, function(sg) {
    nd <- sg$data
    nd$Rx <- factor(sg$rx, levels = levels(sg$data$Rx))
    list(
      os = direct_gamma_pop_avg(component$original_os, os_draw, nd, times, sg$weights),
      pfs = direct_gamma_pop_avg(component$original_pfs, pfs_draw, nd, times, sg$weights)
    )
  })
}

#' Check the closed-form curves against the pipeline's predict()-based helper
#'
#' @param component Joint sampling component.
#' @param subgroups Output of `psa_violation_subgroups()`.
#' @param data_original The `data` frame model_fun() passes for biomarker
#'   subgroups (may contain incomplete patients, which predict() drops).
#' @param data_complete Complete-case cohort used for the control curve.
#' @param times Weekly grid.
#' @param draws Draw indices to check.
#' @param tol Maximum absolute difference allowed.
#' @return Data frame of max absolute differences per draw, curve and outcome.
validate_direct_curves <- function(component, subgroups, data_original,
                                   data_complete, times, draws = 1:3,
                                   tol = 1e-8, psa_params = NULL, params = NULL) {
  control_strategy <- get_control_strategy()
  rows <- list()
  add_row <- function(idx, curve, outcome, ref, direct) {
    rows[[length(rows) + 1]] <<- data.frame(
      draw = idx, curve = curve, outcome = outcome,
      max_abs_diff = max(abs(ref - direct)), stringsAsFactors = FALSE
    )
  }
  for (idx in draws) {
    weights <- rep(1 / nrow(data_complete), nrow(data_complete))
    model_idx <- idx
    if (!is.null(psa_params)) {
      p <- params
      for (nm in intersect(names(psa_params), names(p))) p[[nm]] <- psa_params[[nm]][idx]
      weights <- model_population_weights(p, data_complete)
      model_idx <- psa_params$model_idx[idx]
    }
    draw_subgroups <- psa_violation_subgroups(data_complete, weights = weights)
    direct <- psa_draw_subgroup_curves(component, model_idx, draw_subgroups, times)
    for (outcome in c("os", "pfs")) {
      ref_control <- generate_psa_population_averaged_predictions(
        component, biomarker_name = NULL, outcome = outcome,
        sample_idx = model_idx, data_original = data_complete, time_points = times, weights = weights
      )
      add_row(idx, control_strategy, outcome, ref_control,
              direct[[control_strategy]][[outcome]])
      for (bm in get_biomarkers()) {
        ref <- generate_psa_population_averaged_predictions(
          component, biomarker_name = bm, outcome = outcome,
          sample_idx = model_idx, data_original = data_complete, time_points = times, weights = weights
        )
        add_row(idx, paste0(bm, "_pos"), outcome, ref$positive,
                direct[[paste0(bm, "_pos")]][[outcome]])
        add_row(idx, paste0(bm, "_neg"), outcome, ref$negative,
                direct[[paste0(bm, "_neg")]][[outcome]])
      }
    }
  }
  out <- do.call(rbind, rows)
  if (any(!is.finite(out$max_abs_diff)) || max(out$max_abs_diff) > tol) {
    stop("Closed-form gamma curves differ from the pipeline predictions ",
         "(max |diff| = ", signif(max(out$max_abs_diff), 3), ", tol = ", tol, ")")
  }
  out
}

#' Trapezoidal area of a weekly curve in years
#'
#' @param y Values on the weekly grid.
#' @param cl Cycle length in years.
trapezoid_area_years <- function(y, cl) {
  n <- length(y)
  cl * sum((y[-1] + y[-n]) / 2)
}

#' Violation metrics for one OS/PFS pair
#'
#' @param os,pfs Survival curves on the weekly grid.
#' @param times Weekly grid (weeks).
#' @param cl Cycle length in years.
#' @param v_dw_e Effect discount weights (length of the grid).
#' @param u_np,u_p State utilities.
#' @return One-row data frame.
pfs_os_violation_metrics <- function(os, pfs, times, cl, v_dw_e, u_np, u_p) {
  excess <- pmax(pfs - os, 0)
  viol <- excess > 0
  n_viol <- sum(viol)
  gap_raw <- os - pfs
  gap_clamped <- pmax(gap_raw, 0)
  data.frame(
    violated = n_viol > 0,
    n_violations = n_viol,
    first_week = if (n_viol > 0) min(times[viol]) else NA_real_,
    last_week = if (n_viol > 0) max(times[viol]) else NA_real_,
    max_excess = if (n_viol > 0) max(excess) else 0,
    week_max_excess = if (n_viol > 0) times[which.max(excess)] else NA_real_,
    excess_area_years = trapezoid_area_years(excess, cl),
    pp_area_raw_years = trapezoid_area_years(gap_raw, cl),
    pp_area_clamped_years = trapezoid_area_years(gap_clamped, cl),
    pf_area_raw_years = trapezoid_area_years(pfs, cl),
    pf_area_clamped_years = trapezoid_area_years(pmin(pfs, os), cl),
    os_area_years = trapezoid_area_years(os, cl),
    qaly_raw_minus_clamped = sum(excess * (u_np - u_p) * cl * v_dw_e),
    qaly_clamped = sum((pmin(pfs, os) * u_np + gap_clamped * u_p) * cl * v_dw_e)
  )
}

#' Cache file for the violation diagnostics
#'
#' @param n Number of draws (default `n_samples`).
#' @param directory Cache directory.
pfs_os_violation_cache_path <- function(n = NULL, directory = cache_dir()) {
  if (is.null(n)) n <- get("n_samples", envir = globalenv())
  file.path(directory, paste0("pfs_os_violations_n", n, ".rds"))
}

#' Run (or load) the per-draw violation diagnostics for the whole sampling cache
#'
#' @param sampling_models Sampling cache list (script 06).
#' @param data_complete Complete-case cohort.
#' @param data_original The `data` frame (for the validation step only).
#' @param params Base-case parameter list (utilities, discount rate, cl).
#' @param time_horizon Horizon in weeks.
#' @param cache_path Where to store the result; regenerated when the stored
#'   fingerprint differs from the current sampling cache and utilities.
#' @param n_draws Draws to evaluate (default: all).
#' @param verbose Print progress.
#' @return List: `per_curve` (draw x curve data frame), `curves` (subgroup
#'   metadata), `timing` (share of draws violating at each week and mean excess
#'   by week, per curve), `examples` (raw curves of the worst draw per curve),
#'   `validation`, `fingerprint`, `n_draws`, `elapsed`.
run_pfs_os_violation_diagnostics <- function(sampling_models, data_complete,
                                             data_original, params,
                                             time_horizon,
                                             cache_path = pfs_os_violation_cache_path(),
                                             n_draws = NULL, verbose = TRUE,
                                             psa_params = NULL) {
  component <- get_joint_sampling_models(sampling_models)
  if (is.null(component$draws)) {
    stop("The sampling cache is not a multivariate-normal draw cache (issue #156)")
  }
  n_total <- if (is.null(psa_params)) nrow(component$draws$os) else nrow(psa_params)
  if (is.null(n_draws)) n_draws <- n_total
  n_draws <- as.integer(n_draws)
  if (n_draws < 1 || n_draws > n_total) stop("n_draws must be in 1..", n_total)

  times <- seq(0, time_horizon)
  cl <- params$cl
  v_dw_e <- 1 / (1 + params$dr_effects)^(times * cl)
  fp_inputs <- list(
    sampling_fingerprint = sampling_models$fingerprint,
    n_draws = n_draws, time_horizon = time_horizon, cl = cl,
    u_np = params$u_np, u_p = params$u_p, dr_effects = params$dr_effects,
    prediction_population = cache_fingerprint(data_complete),
    psa_parameters = cache_fingerprint(psa_params),
    version = "pfs_os_violations_v2"
  )
  fingerprint <- cache_fingerprint(fp_inputs)

  if (file.exists(cache_path)) {
    cached <- readRDS(cache_path)
    if (identical(cached$fingerprint, fingerprint)) {
      if (verbose) cat("Loaded PFS/OS violation diagnostics from", cache_path, "\n")
      return(cached)
    }
    if (verbose) cat("Violation diagnostics cache is stale; regenerating\n")
  }

  subgroups <- psa_violation_subgroups(data_complete)
  validation <- validate_direct_curves(component, subgroups, data_complete,
    data_complete, times, draws = seq_len(min(3L, n_draws)), psa_params = psa_params, params = params)
  curve_names <- names(subgroups)
  n_t <- length(times)
  viol_count <- matrix(0, nrow = length(curve_names), ncol = n_t,
                       dimnames = list(curve_names, NULL))
  excess_sum <- viol_count
  worst <- setNames(vector("list", length(curve_names)), curve_names)
  worst_area <- setNames(rep(-Inf, length(curve_names)), curve_names)
  rows <- vector("list", n_draws * length(curve_names))
  k <- 0L
  t0 <- Sys.time()
  for (idx in seq_len(n_draws)) {
    draw_params <- params
    model_idx <- draw_id <- idx
    draw_subgroups <- subgroups
    if (!is.null(psa_params)) {
      for (nm in intersect(names(psa_params), names(draw_params))) draw_params[[nm]] <- psa_params[[nm]][idx]
      draw_subgroups <- psa_violation_subgroups(data_complete,
        weights = model_population_weights(draw_params, data_complete))
      model_idx <- psa_params$model_idx[idx]
      draw_id <- psa_params$sim[idx]
    }
    curves <- psa_draw_subgroup_curves(component, model_idx, draw_subgroups, times)
    for (cn in curve_names) {
      os <- curves[[cn]]$os
      pfs <- curves[[cn]]$pfs
      m <- pfs_os_violation_metrics(os, pfs, times, cl, v_dw_e,
                                    draw_params$u_np, draw_params$u_p)
      k <- k + 1L
      rows[[k]] <- cbind(data.frame(draw = draw_id, model_idx = model_idx, curve = cn,
                                    stringsAsFactors = FALSE), m)
      if (m$violated) {
        excess <- pmax(pfs - os, 0)
        viol_count[cn, ] <- viol_count[cn, ] + (excess > 0)
        excess_sum[cn, ] <- excess_sum[cn, ] + excess
        if (m$excess_area_years > worst_area[[cn]]) {
          worst_area[[cn]] <- m$excess_area_years
          worst[[cn]] <- data.frame(draw = draw_id, curve = cn, week = times,
                                    os = os, pfs = pfs,
                                    stringsAsFactors = FALSE)
        }
      }
    }
    if (verbose && idx %% 500 == 0) {
      cat(sprintf("  draw %d of %d (%.0f s)\n", idx, n_draws,
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    }
  }
  per_curve <- do.call(rbind, rows)
  rownames(per_curve) <- NULL
  timing <- do.call(rbind, lapply(curve_names, function(cn) {
    data.frame(curve = cn, week = times,
               share_violating = viol_count[cn, ] / n_draws,
               mean_excess = excess_sum[cn, ] / n_draws,
               stringsAsFactors = FALSE)
  }))
  curves_meta <- data.frame(
    curve = curve_names,
    label = vapply(subgroups, `[[`, "", "label"),
    strategy = vapply(subgroups, `[[`, "", "strategy"),
    subgroup = vapply(subgroups, `[[`, "", "subgroup"),
    n_patients = vapply(subgroups, function(s) nrow(s$data), 1L),
    stringsAsFactors = FALSE
  )
  rownames(curves_meta) <- NULL
  result <- list(
    per_curve = per_curve, curves = curves_meta, timing = timing,
    examples = do.call(rbind, worst[!vapply(worst, is.null, TRUE)]),
    validation = validation, fingerprint = fingerprint,
    fingerprint_inputs = fp_inputs, n_draws = n_draws,
    elapsed = as.numeric(difftime(Sys.time(), t0, units = "secs")),
    creation_time = Sys.time()
  )
  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(result, cache_path)
  if (verbose) cat(sprintf("Saved diagnostics for %d draws to %s (%.0f s)\n",
                           n_draws, cache_path, result$elapsed))
  result
}
