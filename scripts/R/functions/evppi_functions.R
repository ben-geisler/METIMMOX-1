# EVPPI calculation functions for value of information analysis
#
# ===============================================================================
# ESTIMATOR (issue #152)
# ===============================================================================
# The expected value of partial perfect information (EVPPI) is estimated by
# nonparametric regression of each strategy's incremental net monetary benefit
# (NMB) on the parameter(s) of interest (Strong, Oakley & Brennan, Med Decis
# Making 2014;34:311-326), as implemented by voi::evppi():
#
#   EVPPI(theta) = E_theta[ max_d E(NMB_d | theta) ] - max_d E(NMB_d)
#
# E(NMB_d | theta) is the fitted value of a generalized additive model (GAM) of
# NMB_d - NMB_ref on theta. The Monte Carlo standard error is the standard
# deviation of the EVPPI recomputed for B draws of the regression coefficients
# from their asymptotic posterior distribution (voi's `se = TRUE`).
#
# Method choice (evppi_gam_formula()):
#   - up to 4 parameters: tensor-product cubic regression spline of all
#     parameters, voi's default (`te(x1, ..., xd, bs = "cr")`, with k = 4 per
#     margin for 4 parameters);
#   - more than 4 parameters: additive cubic regression splines,
#     `s(x1, bs = "cr") + ... + s(xd, bs = "cr")`. Tensor products are
#     infeasible beyond 4 dimensions and voi's default there (a Gaussian
#     process) does not return a standard error. Every current group with more
#     than 4 members is a cost group, and NMB is exactly linear and additive in
#     the unit costs, so the additive form is well specified for them.
#
# Why the previous estimator was replaced: it averaged NMB over the 1000
# nearest of 5000 PSA draws at 500 quantile grid points and took the maximum
# over strategies. Because the neighbourhoods included draws unevenly across
# the parameter range, the mean of the dominant strategy's NMB shifted between
# grid points, and that shift, not a change in the optimal decision, was
# reported as EVPPI. A pure-noise parameter returned a non-zero value, groups
# came out below their members, and values were then capped at EVPI. The
# regression estimator gives approximately zero for a parameter that never
# changes the optimal decision, and run_evppi_analysis() asserts that every
# group EVPPI is at least the EVPPI of its largest member (within Monte Carlo
# error).
# ===============================================================================

#' Cost and effect matrices of a PSA object
#'
#' dampack PSA objects store effects under `effectiveness`; hand-built test
#' objects use `effect`. Accept both without relying on partial `$` matching.
#'
#' @param psa_obj PSA object (dampack `psa` or a list with cost/effect matrices)
#' @return List with numeric matrices `cost` and `effect`
psa_outcome_matrices <- function(psa_obj) {
  effect <- if (!is.null(psa_obj[["effectiveness"]])) {
    psa_obj[["effectiveness"]]
  } else {
    psa_obj[["effect"]]
  }
  if (is.null(psa_obj[["cost"]]) || is.null(effect)) {
    stop("psa_obj must contain 'cost' and 'effectiveness' (or 'effect') matrices")
  }
  list(cost = as.matrix(psa_obj[["cost"]]), effect = as.matrix(effect))
}

#' Expected value of perfect information from an NMB matrix
#'
#' @param nmb_matrix Numeric matrix, one row per PSA draw, one column per strategy
#' @return EVPI (same units as NMB)
calculate_evpi_from_nmb <- function(nmb_matrix) {
  mean(apply(nmb_matrix, 1, max)) - max(colMeans(nmb_matrix))
}

#' GAM formula (right-hand side) used for an EVPPI regression
#'
#' @param param_names Character vector of parameter names
#' @return NULL to use voi's default tensor-product formula (<= 4 parameters),
#'   otherwise a character string with additive cubic regression splines
evppi_gam_formula <- function(param_names) {
  if (length(param_names) <= 4) {
    return(NULL)
  }
  paste(sprintf("s(%s, bs = 'cr')", param_names), collapse = " + ")
}

#' Calculate EVPPI for one parameter or parameter group by regression
#'
#' Wraps voi::evppi() (Strong, Oakley & Brennan 2014). See the file header for
#' the estimator and the method choice.
#'
#' @param psa_obj PSA object created by dampack::make_psa_obj (or a list with
#'   `cost` and `effect` matrices)
#' @param psa_params Data frame with PSA parameter samples, row-aligned with
#'   the PSA outcomes
#' @param param_names Character vector of parameter names (one entry for a
#'   single-parameter EVPPI, several for a joint EVPPI)
#' @param wtp Willingness-to-pay threshold
#' @param method voi::evppi() method (default "gam")
#' @param se Whether to compute the Monte Carlo standard error
#' @param B Number of regression-coefficient draws used for the standard error
#' @param seed RNG seed set before the coefficient draws so that standard
#'   errors are reproducible
#' @param gam_formula Optional right-hand side of the GAM formula; defaults to
#'   evppi_gam_formula(param_names)
#' @param ... Further arguments passed to voi::evppi()
#' @return List with `evppi`, `evppi_se`, `evpi`, `param_names`, `method`,
#'   `gam_formula`, `n_sim` (draws used) and, on failure, `error`. A failed
#'   estimate is NA, never zero.
calculate_evppi_regression <- function(psa_obj, psa_params, param_names, wtp,
                                       method = "gam", se = TRUE, B = 1000,
                                       seed = 123L, gam_formula = NULL, ...) {

  cat("Calculating EVPPI for parameter(s):", paste(param_names, collapse = ", "), "\n")

  outcomes <- psa_outcome_matrices(psa_obj)
  cost_matrix <- outcomes$cost
  effect_matrix <- outcomes$effect
  n_sim <- nrow(cost_matrix)

  result <- list(
    evppi = NA_real_,
    evppi_se = NA_real_,
    evpi = NA_real_,
    param_names = param_names,
    method = method,
    gam_formula = NA_character_,
    n_sim = n_sim,
    error = NULL
  )
  fail <- function(msg) {
    cat("  EVPPI not estimated:", msg, "\n")
    result$error <- msg
    result
  }

  if (nrow(psa_params) != n_sim) {
    return(fail(sprintf("psa_params has %d rows but the PSA has %d draws",
                        nrow(psa_params), n_sim)))
  }

  missing_params <- setdiff(param_names, colnames(psa_params))
  if (length(missing_params) > 0) {
    return(fail(paste("Parameter(s) not found in psa_params:",
                      paste(missing_params, collapse = ", "))))
  }

  inputs <- psa_params[, param_names, drop = FALSE]
  for (param_name in param_names) {
    if (!is.numeric(inputs[[param_name]])) {
      return(fail(paste("Parameter", param_name, "is not numeric")))
    }
  }

  # Keep draws with complete outcomes and complete inputs so the regression
  # rows stay aligned with the PSA draws.
  complete_rows <- complete.cases(cost_matrix) & complete.cases(effect_matrix) &
    complete.cases(inputs)
  if (sum(complete_rows) < n_sim) {
    cat("  Note: using", sum(complete_rows), "of", n_sim,
        "draws with complete outcomes and inputs\n")
    cost_matrix <- cost_matrix[complete_rows, , drop = FALSE]
    effect_matrix <- effect_matrix[complete_rows, , drop = FALSE]
    inputs <- inputs[complete_rows, , drop = FALSE]
  }
  result$n_sim <- nrow(cost_matrix)
  if (result$n_sim < 50) {
    return(fail(paste("Too few complete draws for regression:", result$n_sim)))
  }

  for (param_name in param_names) {
    if (length(unique(inputs[[param_name]])) < 5) {
      return(fail(paste("Parameter", param_name, "has fewer than 5 distinct values")))
    }
  }

  nmb_matrix <- effect_matrix * wtp - cost_matrix
  result$evpi <- calculate_evpi_from_nmb(nmb_matrix)

  if (identical(method, "gam") && is.null(gam_formula)) {
    gam_formula <- evppi_gam_formula(param_names)
  }
  result$gam_formula <- if (is.null(gam_formula)) {
    if (identical(method, "gam")) "voi default (tensor product)" else NA_character_
  } else {
    gam_formula
  }

  voi_args <- list(
    outputs = nmb_matrix,
    inputs = inputs,
    pars = param_names,
    method = method,
    se = se,
    B = B,
    ...
  )
  if (!is.null(gam_formula)) {
    voi_args$gam_formula <- gam_formula
  }

  fit <- tryCatch({
    set.seed(seed)  # coefficient draws for the standard error
    do.call(voi::evppi, voi_args)
  }, error = function(e) e)

  if (inherits(fit, "error")) {
    return(fail(conditionMessage(fit)))
  }

  result$evppi <- as.numeric(fit$evppi[1])
  result$evppi_se <- if (se && "se" %in% names(fit)) as.numeric(fit$se[1]) else NA_real_

  cat(sprintf("  EVPPI: %.4f  (SE %s; %.1f%% of EVPI %.4f; method %s, n = %d)\n",
              result$evppi,
              if (is.finite(result$evppi_se)) sprintf("%.4f", result$evppi_se) else "NA",
              if (result$evpi > 0) 100 * result$evppi / result$evpi else NA_real_,
              result$evpi, method, result$n_sim))

  result
}

#' Check that every group EVPPI is at least the EVPPI of its largest member
#'
#' Perfect information on a set of parameters is worth at least as much as
#' perfect information on any one of them. Regression estimates carry Monte
#' Carlo error, so a group is flagged only when it falls short of its largest
#' member by more than `tol_se` combined standard errors, with a floor of
#' `tol_frac_evpi` of EVPI for estimates whose standard errors are tiny or
#' unavailable.
#'
#' @param evppi_results Data frame returned by run_evppi_analysis()
#' @param param_groups Named list of parameter groups
#' @param tol_se Number of combined standard errors tolerated
#' @param tol_frac_evpi Absolute tolerance floor as a fraction of EVPI
#' @return Data frame with one row per group present in `evppi_results`:
#'   `group`, `group_evppi`, `group_se`, `max_member`, `member_evppi`,
#'   `member_se`, `shortfall` (member minus group), `tolerance`, `violation`
check_evppi_group_consistency <- function(evppi_results, param_groups,
                                          tol_se = 2, tol_frac_evpi = 0.01) {
  empty <- data.frame(
    group = character(), group_evppi = numeric(), group_se = numeric(),
    max_member = character(), member_evppi = numeric(), member_se = numeric(),
    shortfall = numeric(), tolerance = numeric(), violation = logical(),
    stringsAsFactors = FALSE
  )
  if (is.null(param_groups) || length(param_groups) == 0 ||
      is.null(evppi_results) || nrow(evppi_results) == 0) {
    return(empty)
  }

  rows <- lapply(names(param_groups), function(group_name) {
    group_row <- evppi_results[
      evppi_results$parameter == paste0("[GROUP] ", group_name), , drop = FALSE]
    members <- intersect(param_groups[[group_name]], evppi_results$parameter)
    if (nrow(group_row) == 0 || length(members) == 0) {
      return(NULL)
    }
    member_rows <- evppi_results[match(members, evppi_results$parameter), , drop = FALSE]
    member_evppi <- member_rows$evppi
    if (all(is.na(member_evppi)) || is.na(group_row$evppi[1])) {
      return(NULL)
    }
    idx <- which.max(replace(member_evppi, is.na(member_evppi), -Inf))

    group_se <- group_row$evppi_se[1]
    member_se <- member_rows$evppi_se[idx]
    se_tol <- if (is.finite(group_se) && is.finite(member_se)) {
      tol_se * sqrt(group_se^2 + member_se^2)
    } else {
      0
    }
    evpi <- group_row$evpi[1]
    tolerance <- max(se_tol, tol_frac_evpi * if (is.finite(evpi)) evpi else 0)
    shortfall <- member_evppi[idx] - group_row$evppi[1]

    data.frame(
      group = group_name,
      group_evppi = group_row$evppi[1],
      group_se = group_se,
      max_member = members[idx],
      member_evppi = member_evppi[idx],
      member_se = member_se,
      shortfall = shortfall,
      tolerance = tolerance,
      violation = shortfall > tolerance,
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) {
    return(empty)
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Stop if any group EVPPI falls below its largest member beyond tolerance
#'
#' @inheritParams check_evppi_group_consistency
#' @return The consistency table, invisibly
assert_evppi_group_consistency <- function(evppi_results, param_groups,
                                           tol_se = 2, tol_frac_evpi = 0.01) {
  consistency <- check_evppi_group_consistency(
    evppi_results, param_groups, tol_se = tol_se, tol_frac_evpi = tol_frac_evpi
  )
  bad <- consistency[consistency$violation, , drop = FALSE]
  if (nrow(bad) > 0) {
    stop(
      "Group EVPPI below its largest member beyond Monte Carlo tolerance:\n",
      paste(sprintf(
        "  [GROUP] %s = %.3f (SE %.3f) < %s = %.3f (SE %.3f); shortfall %.3f > tolerance %.3f",
        bad$group, bad$group_evppi, bad$group_se, bad$max_member,
        bad$member_evppi, bad$member_se, bad$shortfall, bad$tolerance
      ), collapse = "\n")
    )
  }
  invisible(consistency)
}

#' Run complete EVPPI analysis for all parameters and parameter groups
#'
#' @param psa_obj PSA object from dampack
#' @param psa_params Data frame with PSA parameter samples
#' @param wtp Willingness-to-pay threshold
#' @param evppi_params Vector of parameter names to analyze individually
#' @param param_groups Optional named list of parameter groups for joint EVPPI
#' @param seed RNG seed used for the standard-error coefficient draws
#' @param B Number of coefficient draws for the standard errors
#' @param assert_groups Stop when a group EVPPI falls below its largest member
#'   beyond Monte Carlo tolerance (see assert_evppi_group_consistency()); when
#'   FALSE the violation is reported as a warning instead
#' @return Data frame with one row per parameter and per group (prefixed
#'   "[GROUP] "): `parameter`, `evppi`, `evppi_se`, `evpi`,
#'   `evppi_percent_of_evpi`, `method`, `n_params`, `n_sim`, `error`. The
#'   group-consistency table is attached as attribute "group_consistency".
run_evppi_analysis <- function(psa_obj, psa_params, wtp, evppi_params,
                               param_groups = NULL, seed = 123L, B = 1000,
                               assert_groups = TRUE) {

  outcomes <- psa_outcome_matrices(psa_obj)
  nmb_matrix <- outcomes$effect * wtp - outcomes$cost
  evpi_manual <- calculate_evpi_from_nmb(nmb_matrix)

  cat("\n=== EVPPI Analysis at WTP =", wtp, "===\n")
  cat("Total EVPI:", round(evpi_manual, 4), "\n")
  cat("Estimator: voi::evppi() nonparametric regression (GAM);",
      "standard errors from B =", B, "coefficient draws\n\n")

  empty_results <- data.frame(
    parameter = character(),
    evppi = numeric(),
    evppi_se = numeric(),
    evpi = numeric(),
    evppi_percent_of_evpi = numeric(),
    method = character(),
    n_params = integer(),
    n_sim = integer(),
    error = character(),
    stringsAsFactors = FALSE
  )

  # Only proceed if we have meaningful EVPI
  if (evpi_manual <= 0.01) {
    cat("EVPI is too small (", round(evpi_manual, 4), ") - skipping EVPPI analysis\n")
    return(empty_results)
  }

  # Filter to only include parameters that exist and have variation
  available_params <- intersect(evppi_params, colnames(psa_params))

  params_with_variation <- character()
  for (param in available_params) {
    param_values <- psa_params[[param]]
    if (is.numeric(param_values) && !all(is.na(param_values))) {
      param_values <- param_values[!is.na(param_values)]
      if (length(unique(param_values)) >= 5 && sd(param_values) > 0) {
        params_with_variation <- c(params_with_variation, param)
      }
    }
  }

  cat("Parameters with sufficient variation:",
      paste(params_with_variation, collapse = ", "), "\n\n")

  result_row <- function(label, result) {
    evppi_percent <- if (is.finite(result$evpi) && result$evpi > 0) {
      100 * result$evppi / result$evpi
    } else {
      NA_real_
    }
    data.frame(
      parameter = label,
      evppi = result$evppi,
      evppi_se = result$evppi_se,
      evpi = result$evpi,
      evppi_percent_of_evpi = evppi_percent,
      method = result$method,
      n_params = length(result$param_names),
      n_sim = result$n_sim,
      error = if (is.null(result$error)) "" else result$error,
      stringsAsFactors = FALSE
    )
  }

  evppi_results <- empty_results

  # Calculate EVPPI for each parameter individually
  for (param in params_with_variation) {
    result <- calculate_evppi_regression(psa_obj, psa_params, param,
                                         wtp = wtp, se = TRUE, B = B, seed = seed)
    evppi_results <- rbind(evppi_results, result_row(param, result))
  }

  # Sort individual parameters by EVPPI value (descending, NA last)
  if (nrow(evppi_results) > 0) {
    evppi_results <- evppi_results[order(-evppi_results$evppi, na.last = TRUE), ]
  }

  # Calculate EVPPI for parameter groups (if provided)
  if (!is.null(param_groups) && length(param_groups) > 0) {
    cat("\n=== Grouped Parameter EVPPI ===\n")

    for (group_name in names(param_groups)) {
      group_params <- intersect(param_groups[[group_name]], params_with_variation)

      if (length(group_params) >= 2) {
        cat("\nCalculating EVPPI for group:", group_name,
            "(", paste(group_params, collapse = ", "), ")\n")
        result <- calculate_evppi_regression(psa_obj, psa_params, group_params,
                                             wtp = wtp, se = TRUE, B = B, seed = seed)
        evppi_results <- rbind(evppi_results,
                               result_row(paste0("[GROUP] ", group_name), result))
      }
    }
  }
  rownames(evppi_results) <- NULL

  failed <- evppi_results$parameter[nzchar(evppi_results$error)]
  if (length(failed) > 0) {
    warning("EVPPI could not be estimated for: ", paste(failed, collapse = ", "),
            " (see the 'error' column; these rows are NA, not zero)")
  }

  # Perfect information on a group is worth at least as much as on any member.
  consistency <- check_evppi_group_consistency(evppi_results, param_groups)
  attr(evppi_results, "group_consistency") <- consistency
  if (any(consistency$violation)) {
    if (assert_groups) {
      assert_evppi_group_consistency(evppi_results, param_groups)
    } else {
      warning("Group EVPPI below its largest member for: ",
              paste(consistency$group[consistency$violation], collapse = ", "))
    }
  } else if (nrow(consistency) > 0) {
    cat("\nGroup consistency check passed for", nrow(consistency),
        "group(s): every group EVPPI >= largest member within tolerance\n")
  }

  evppi_results
}

#' Calculate Population-Level EVPPI
#'
#' Scales per-patient EVPPI values to the eligible population level
#' using annual incidence, research horizon, and discounting.
#'
#' @param evppi_per_patient Numeric. Per-patient EVPPI value in EUR.
#' @param annual_incidence Numeric. Annual eligible patients (default: 1500).
#' @param research_horizon Numeric. Years of research value (default: 10).
#' @param discount_rate Numeric. Annual discount rate for research benefits (default: 0.035).
#' @return Numeric. Population-level EVPPI in millions of EUR.
#'
#' @details
#' The calculation accounts for the time value of research by discounting
#' future benefits. The formula is:
#'   Population EVPPI = Per-patient EVPPI × Annual incidence ×
#'                      Sum of discount factors over research horizon
#'
#' Example: For EUR 1,000 per-patient EVPPI with 1,500 patients/year over 10 years:
#'   Sum of discount factors ≈ 8.317 (using 3.5% discount rate)
#'   Population EVPPI = EUR 1,000 × 1,500 × 8.317 = EUR 12.48 million
#'
#' The scaling is linear, so the same function converts a per-patient
#' standard error to a population-level standard error.
#'
#' @export
calculate_population_evppi <- function(evppi_per_patient,
                                       annual_incidence = 1500,
                                       research_horizon = 10,
                                       discount_rate = 0.035) {
  # Calculate sum of discount factors over research horizon
  # (Present value of an annuity of 1 unit per year for n years)
  years <- 1:research_horizon
  discount_factors <- 1 / (1 + discount_rate)^years
  sum_discount_factors <- sum(discount_factors)

  # Scale per-patient value to population
  population_evppi_eur <- evppi_per_patient * annual_incidence * sum_discount_factors

  # Convert to millions
  population_evppi_millions <- population_evppi_eur / 1e6

  return(population_evppi_millions)
}


#' Extract treatment-biomarker interaction coefficients from sampling models
#'
#' Extracts the biomarker:Rx interaction coefficients from resampled survival
#' models and returns them as additional columns aligned with PSA iterations.
#' These can be appended to psa_params for EVPPI analysis of treatment effect
#' modification parameters.
#'
#' @param sampling_models List of resampled models (global variable from 06_sampling.R)
#' @param n_sim Number of PSA iterations (must match length of sampling_models$*$samples)
#' @return Data frame with n_sim rows and columns for each extracted interaction
#'   coefficient (b_<biomarker>_rx_<outcome>). Returns NULL if extraction fails.
extract_interaction_coefficients <- function(sampling_models, n_sim) {
  if (is.null(sampling_models)) {
    warning("sampling_models is NULL - cannot extract interaction coefficients")
    return(NULL)
  }

  biomarkers <- get_biomarkers()
  outcomes <- c("os", "pfs")

  result <- data.frame(row.names = seq_len(n_sim))

  cat("\n=== Extracting interaction coefficients from sampling models ===\n")

  # Since issue #156 the cache holds one joint component (`joint`) shared by
  # every biomarker strategy, with the coefficient draws stored as matrices.
  # Legacy caches and test fixtures hold one component per biomarker key with
  # per-sample coefficient vectors; both layouts are read here.
  joint_component <- sampling_models$joint
  component_for <- function(biomarker) {
    if (!is.null(joint_component)) return(joint_component)
    sampling_models[[biomarker]]
  }

  # Coefficient vector of draw i for one outcome, or NULL when unavailable.
  draw_coefficients <- function(component, outcome, i) {
    if (!is.null(component$draws)) {
      m <- component$draws[[outcome]]
      if (is.null(m) || i > nrow(m)) return(NULL)
      return(m[i, ])
    }
    sample_i <- component$samples[[i]]
    if (is.null(sample_i)) return(NULL)
    sample_i[[outcome]]$coefficients
  }

  for (biomarker in biomarkers) {
    strategy_models <- component_for(biomarker)
    if (is.null(strategy_models)) {
      cat("  Warning:", biomarker, "not found in sampling_models - skipping\n")
      next
    }

    for (outcome in outcomes) {
      col_name <- paste0("b_", biomarker, "_rx_", outcome)
      values <- numeric(n_sim)
      n_extracted <- 0
      n_failed <- 0
      n_missing <- 0

      for (i in seq_len(n_sim)) {
        coefs <- draw_coefficients(strategy_models, outcome, i)

        if (is.null(coefs)) {
          values[i] <- NA
          n_failed <- n_failed + 1
          next
        }

        # Find the treatment-by-biomarker interaction coefficient.
        #
        # Term ORDER is not stable and must not be assumed. The economic
        # formula is "... + Rx + crp * Rx + tmb_braf * Rx", so Rx appears
        # before the biomarker in the term expansion and R labels the
        # interaction "Rx:crp", not "crp:Rx". Combined with flexsurvreg
        # appending factor levels, the real name is
        # "RxExperimental arm:crp1". Anchored patterns such as
        # "^crp.*:Rx" (or the earlier "crp.*:Rx") never match it, which
        # silently dropped every b_*_rx_* row from the EVPPI tables.
        #
        # Match structurally instead: an interaction term (contains ":")
        # mentioning both Rx and the biomarker, each at a term boundary so
        # one biomarker ID cannot match another's as a substring.
        coef_names <- names(coefs)
        matching_names <- coef_names[
          grepl(":", coef_names, fixed = TRUE) &
            grepl(paste0("(^|:)", biomarker), coef_names) &
            grepl("(^|:)Rx", coef_names)
        ]

        if (length(matching_names) > 0) {
          values[i] <- coefs[matching_names[1]]
          n_extracted <- n_extracted + 1
        } else {
          values[i] <- NA
          n_missing <- n_missing + 1
        }
      }

      result[[col_name]] <- values

      cat(sprintf("  %s: extracted=%d, failed=%d, missing=%d\n",
                  col_name, n_extracted, n_failed, n_missing))

      if (n_extracted / n_sim < 0.9) {
        warning("Only ", round(n_extracted / n_sim * 100, 1),
                "% of iterations have valid ", col_name, " values")
      }
    }
  }

  # Remove columns that are entirely NA (coefficient not in model formula)
  all_na_cols <- vapply(result, function(x) all(is.na(x)), logical(1))
  if (any(all_na_cols)) {
    cat("  Removing all-NA columns:",
        paste(names(result)[all_na_cols], collapse = ", "), "\n")
    result <- result[, !all_na_cols, drop = FALSE]
  }

  # Report summary statistics
  cat("\nExtracted coefficient summary:\n")
  for (col in names(result)) {
    vals <- result[[col]]
    valid_vals <- vals[!is.na(vals)]
    if (length(valid_vals) > 0) {
      cat(sprintf("  %s: mean=%.4f, sd=%.4f, range=[%.4f, %.4f], n_valid=%d\n",
                  col, mean(valid_vals), sd(valid_vals),
                  min(valid_vals), max(valid_vals), length(valid_vals)))
    }
  }

  if (ncol(result) == 0) {
    warning("No interaction coefficients could be extracted")
    return(NULL)
  }

  return(result)
}


#' Get interaction EVPPI parameter names for the economic model
#'
#' @return Character vector of parameter names (e.g., "b_crp_rx_os")
get_interaction_evppi_params <- function() {
  biomarkers <- get_biomarkers()
  outcomes <- c("os", "pfs")

  params <- character(0)
  for (biomarker in biomarkers) {
    for (outcome in outcomes) {
      params <- c(params, paste0("b_", biomarker, "_rx_", outcome))
    }
  }

  return(params)
}


#' Add interaction parameter groups for EVPPI analysis
#'
#' Extends existing param_groups with one group per biomarker (its OS and PFS
#' treatment-by-biomarker interaction coefficients) and, when more than one
#' biomarker is configured, a joint group `interaction_all` holding every
#' interaction coefficient. The joint group is what Figure 5 and Table S7
#' report as "Biomarker-treatment interaction"; group EVPPIs are never summed.
#'
#' @param existing_groups Named list of existing parameter groups
#' @return Updated named list with additional interaction groups
add_interaction_param_groups <- function(existing_groups) {
  interaction_params <- get_interaction_evppi_params()

  if (length(interaction_params) == 0) {
    return(existing_groups)
  }

  biomarkers <- get_biomarkers()
  n_biomarker_groups <- 0L

  for (biomarker in biomarkers) {
    bio_params <- grep(paste0("^b_", biomarker, "_rx_"),
                       interaction_params, value = TRUE)
    if (length(bio_params) >= 2) {
      existing_groups[[paste0("interaction_", biomarker)]] <- bio_params
      n_biomarker_groups <- n_biomarker_groups + 1L
    }
  }

  if (n_biomarker_groups > 1L && length(interaction_params) >= 2) {
    existing_groups[["interaction_all"]] <- interaction_params
  }

  return(existing_groups)
}


#' Assemble the parameter vector and groups used by EVPPI analyses
#'
#' @param param_distributions Named non-survival parameter distributions
#' @param param_groups Named non-survival parameter groups
#' @param include_interactions Whether to add treatment-biomarker interactions
#' @param available_params Optional names available in a realized PSA data frame
#' @return List with `params` and `groups`
configure_evppi_analysis <- function(param_distributions, param_groups,
                                     include_interactions = TRUE,
                                     available_params = NULL) {
  params <- names(param_distributions)
  groups <- param_groups

  if (include_interactions) {
    interaction_params <- get_interaction_evppi_params()
    if (!is.null(available_params)) {
      interaction_params <- intersect(interaction_params, available_params)
    }
    params <- unique(c(params, interaction_params))
    groups <- add_interaction_param_groups(groups)
  }

  list(params = params, groups = groups)
}
