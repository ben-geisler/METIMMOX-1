# ===============================================================================
# SINGLE-MODEL CEA HELPER FUNCTIONS
# ===============================================================================
# Functions for running and summarizing the single joint economic model.
#
# Prerequisites before calling run_basecase():
#   1. Global variables are set (source 02_setup_and_global_variables.R)
#   2. Biomarker strategies are defined (source 03_biomarker_strategies.R)
#   3. Survival predictions and l_params_base are available (source 04/05)
#   4. model_fun.R and calculate_outcomes.R are loaded
# ===============================================================================

if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(here, dplyr, ggplot2, scales, dampack)

# Source model configuration if needed. The configuration file defines the single
# joint economic model and the economic strategies included in CEA.
if (!exists("get_current_model_config") || !exists("get_strategies")) {
  model_configs_path <- here::here("R/model_configs.R")
  if (file.exists(model_configs_path)) {
    source(model_configs_path)
  } else {
    stop("model_configs.R not found. Please ensure R/model_configs.R exists.")
  }
}

# Source cache-path helpers if needed (utility label + cache file locations).
if (!exists("psa_obj_path") || !exists("resolve_util_label")) {
  cache_paths_path <- here::here("R/cache_paths.R")
  if (file.exists(cache_paths_path)) {
    source(cache_paths_path)
  } else {
    stop("cache_paths.R not found. Please ensure R/cache_paths.R exists.")
  }
}

#' Run base case analysis for the single joint economic model
#'
#' @param params Parameter list. Defaults to global l_params_base.
#' @param verbose Logical, print progress messages.
#' @return List containing base_results, icer_obj, and model_config.
run_basecase <- function(params = NULL, verbose = TRUE) {

  if (is.null(params)) {
    if (!exists("l_params_base", inherits = TRUE)) {
      stop("l_params_base not found. Source 04/05 before calling run_basecase().")
    }
    params <- get("l_params_base", inherits = TRUE)
  }

  if (!exists("model_fun", mode = "function", inherits = TRUE)) {
    stop("model_fun() not found. Source R/model_fun.R first.")
  }

  if (verbose) {
    config <- get_current_model_config()
    cat("Running base case for:", config$label, "\n")
    cat("Strategies:", paste(get_strategies(), collapse = ", "), "\n")
  }

  base_results <- model_fun(params)

  icer_obj <- dampack::calculate_icers(
    cost = base_results$Cost,
    effect = base_results$Effect,
    strategies = base_results$Strategy
  )

  if (verbose) {
    cat("Base case complete. Strategies:",
        paste(base_results$Strategy, collapse = ", "), "\n")
  }

  list(
    base_results = base_results,
    icer_obj = icer_obj,
    model_config = get_current_model_config()
  )
}

#' Calculate pairwise ICERs against the control strategy
#'
#' Convention (issue #153): every ICER in the returned table is a PAIRWISE
#' comparison of one strategy against standard of care, not a frontier ICER.
#' A strategy can be "Dominated" on the dampack efficiency frontier (more costly
#' and less effective than another biomarker strategy) while still having a
#' finite pairwise ICER versus standard of care. The frontier status is
#' therefore attached alongside the pairwise status as `Frontier_Status`
#' (dampack codes: ND = on the frontier, D = dominated, ED = extendedly
#' dominated) so that reports can state both without contradiction.
#'
#' @param results Data frame with Strategy, Cost, and Effect columns
#' @return Results with incremental outcomes versus control, the pairwise ICER,
#'   the pairwise `Status` ("Reference", "Pairwise ICER vs SoC", "Dominated by
#'   SoC", "Cost-saving vs SoC"), and the dampack `Frontier_Status`.
calculate_pairwise_icers <- function(results) {
  control <- results[results$Strategy == get_control_strategy(), ]
  if (nrow(control) != 1) stop("Control strategy must appear exactly once.")
  out <- transform(
    results,
    Inc_Cost = Cost - control$Cost,
    Inc_Effect = Effect - control$Effect
  )
  out$ICER <- out$Inc_Cost / out$Inc_Effect
  out$ICER[out$Strategy == get_control_strategy()] <- NA_real_
  out$Status <- "Pairwise ICER vs SoC"
  out$Status[out$Inc_Cost > 0 & out$Inc_Effect <= 0] <- "Dominated by SoC"
  out$Status[out$Inc_Cost < 0 & out$Inc_Effect > 0] <- "Cost-saving vs SoC"
  out$Status[out$Strategy == get_control_strategy()] <- "Reference"
  out$Frontier_Status <- frontier_status(out)
  out
}

#' dampack efficiency-frontier status for a set of strategies
#'
#' @param results Data frame with Strategy, Cost, and Effect columns
#' @return Character vector aligned with `results$Strategy`: "ND", "D", or "ED"
frontier_status <- function(results) {
  icers <- dampack::calculate_icers(
    cost = results$Cost,
    effect = results$Effect,
    strategies = results$Strategy
  )
  as.character(icers$Status)[match(results$Strategy, icers$Strategy)]
}

#' Human-readable frontier label for tables
frontier_label <- function(status) {
  labels <- c(ND = "On frontier", D = "Dominated", ED = "Extendedly dominated")
  out <- unname(labels[as.character(status)])
  out[is.na(out)] <- "--"
  out
}

#' Validate a cached PSA object against what the current session expects
#'
#' Mirrors the checks that 10_PSA.R performs before it accepts an existing
#' cache (issue #156), so that a report cannot render a stale PSA (issue #157).
#' Every expectation is optional; only the supplied ones are checked, but the
#' structural checks (failed-draw policy when `psa_failed_draw_policy()` is
#' loaded, presence of `model_idx`) always run because a cache failing them
#' would be regenerated by 10_PSA.R regardless of the inputs.
#'
#' @param psa_obj dampack PSA object as saved by 10_PSA.R.
#' @param expected_fingerprint Either the list returned by
#'   `psa_cache_fingerprint()` (fingerprint plus the inputs, so the differing
#'   input can be named in the error) or a bare fingerprint string.
#' @param strategies Expected strategy vector.
#' @param n_sim Expected number of requested PSA draws.
#' @param sampling_fingerprint Fingerprint of the sampling cache the PSA must
#'   have been built from (`sampling_models$fingerprint`).
#' @param cache_file Path used in messages.
#' @return Invisibly TRUE; stops with an actionable message otherwise.
validate_psa_cache <- function(psa_obj,
                               expected_fingerprint = NULL,
                               strategies = NULL,
                               n_sim = NULL,
                               sampling_fingerprint = NULL,
                               cache_file = "PSA cache") {
  rerun <- paste("Rerun analysis/10_PSA.R, then 11_EVPPIs.R and",
                 "12_scenario_EVPPIs.R, before rendering.")
  fail <- function(...) {
    stop("Stale PSA cache ", cache_file, ": ", ..., " ", rerun, call. = FALSE)
  }

  if (!is.null(strategies) &&
      !identical(as.character(psa_obj$strategies), as.character(strategies))) {
    fail("strategy set mismatch (cached: ",
         paste(psa_obj$strategies, collapse = ", "), "; expected: ",
         paste(strategies, collapse = ", "), ").")
  }
  if (!is.null(n_sim)) {
    cached_requested <- if (!is.null(psa_obj$requested_n_sim)) {
      psa_obj$requested_n_sim
    } else {
      psa_obj$n_sim
    }
    if (!isTRUE(cached_requested == n_sim)) {
      fail("requested n_sim mismatch (cached: ", cached_requested,
           "; expected: ", n_sim, ").")
    }
  }
  if (exists("psa_failed_draw_policy", mode = "function", inherits = TRUE) &&
      !identical(psa_obj$failed_draw_policy, psa_failed_draw_policy())) {
    fail("legacy failed-draw handling policy.")
  }
  if (is.null(psa_obj$model_idx)) {
    fail("no per-row model_idx (cache predates issue #156).")
  }
  if (!is.null(sampling_fingerprint) &&
      !identical(psa_obj$fingerprint_inputs$sampling_fingerprint,
                 sampling_fingerprint)) {
    fail("built from a different sampling cache (sampling fingerprint mismatch).")
  }
  if (!is.null(expected_fingerprint)) {
    expected_string <- if (is.list(expected_fingerprint)) {
      expected_fingerprint$fingerprint
    } else {
      expected_fingerprint
    }
    if (!identical(psa_obj$fingerprint, expected_string)) {
      differing <- character(0)
      if (is.list(expected_fingerprint) && !is.null(psa_obj$fingerprint_inputs)) {
        for (key in names(expected_fingerprint$inputs)) {
          if (!identical(psa_obj$fingerprint_inputs[[key]],
                         expected_fingerprint$inputs[[key]])) {
            differing <- c(differing, key)
          }
        }
      }
      detail <- if (length(differing) > 0) {
        paste0(" (differing inputs: ", paste(differing, collapse = ", "), ")")
      } else {
        ""
      }
      fail("input fingerprint mismatch", detail, ".")
    }
  }
  invisible(TRUE)
}

#' Load the single-model PSA cache
#'
#' Without expectations the loader behaves as before (a missing file gives a
#' warning and NULL) apart from the structural checks of
#' `validate_psa_cache()` and mandatory outcome/parameter pair validation. Both
#' files must exist and belong to the same generation. Reports pass
#' `required = TRUE` and the expected
#' fingerprint, strategies, n_sim and sampling fingerprint so that they stop on
#' an absent or stale cache instead of rendering placeholders (issue #157).
#'
#' @param util_label Utility-source label. Defaults to global utility_source_label.
#' @param directory Directory containing PSA cache files.
#' @param verbose Logical, print progress messages.
#' @param required Logical; stop (rather than warn and return NULL) when the
#'   cache file is absent.
#' @param expected_fingerprint,strategies,n_sim,sampling_fingerprint Optional
#'   expectations forwarded to `validate_psa_cache()`.
#' @return dampack PSA object, or NULL if the cache is not available and
#'   `required = FALSE`.
load_psa_cache <- function(util_label = NULL,
                           directory = cache_dir(),
                           verbose = TRUE,
                           required = FALSE,
                           expected_fingerprint = NULL,
                           strategies = NULL,
                           n_sim = NULL,
                           sampling_fingerprint = NULL) {

  cache_file <- psa_obj_path(util_label, directory)

  if (!file.exists(cache_file)) {
    msg <- paste0("PSA cache not found: ", cache_file,
                  ". Run analysis/10_PSA.R first.")
    if (isTRUE(required)) stop(msg, call. = FALSE)
    warning(msg)
    return(NULL)
  }

  psa_obj <- readRDS(cache_file)
  paired_path <- psa_params_path(util_label, directory)
  if (!file.exists(paired_path)) {
    stop("PSA parameter cache missing: ", paired_path,
         ". Rerun analysis/10_PSA.R.", call. = FALSE)
  }
  paired_params <- readRDS(paired_path)
  validate_psa_pair(psa_obj, paired_params)
  validate_psa_cache(
    psa_obj,
    expected_fingerprint = expected_fingerprint,
    strategies = strategies,
    n_sim = n_sim,
    sampling_fingerprint = sampling_fingerprint,
    cache_file = cache_file
  )

  if (verbose) {
    cat("Loaded PSA cache:", cache_file, "\n")
    if (!is.null(psa_obj$strategies)) {
      cat("Strategies:", paste(psa_obj$strategies, collapse = ", "), "\n")
    }
  }

  psa_obj
}

#' Load the single-model PSA parameter cache
#'
#' @param util_label Utility-source label. Defaults to global utility_source_label.
#' @param directory Directory containing PSA cache files.
#' @param verbose Logical, print progress messages.
#' @param required Logical; stop (rather than warn and return NULL) when the
#'   cache file is absent.
#' @param seed Expected RNG seed recorded on the draws (the "seed" attribute);
#'   checked with `psa_samples_seed_matches()` when that function is loaded,
#'   otherwise by identity.
#' @param n_sim Expected number of rows (the effective `psa_obj$n_sim`).
#' @return Data frame of PSA parameters, or NULL if the cache is not available
#'   and `required = FALSE`.
load_psa_params_cache <- function(util_label = NULL,
                                  directory = cache_dir(),
                                  verbose = TRUE,
                                  required = FALSE,
                                  seed = NULL,
                                  n_sim = NULL) {

  cache_file <- psa_params_path(util_label, directory)

  if (!file.exists(cache_file)) {
    msg <- paste0("PSA parameter cache not found: ", cache_file,
                  ". Run analysis/10_PSA.R first.")
    if (isTRUE(required)) stop(msg, call. = FALSE)
    warning(msg)
    return(NULL)
  }

  psa_params <- readRDS(cache_file)
  paired_path <- psa_obj_path(util_label, directory)
  if (!file.exists(paired_path)) {
    stop("PSA outcome cache missing: ", paired_path,
         ". Rerun analysis/10_PSA.R.", call. = FALSE)
  }
  paired_obj <- readRDS(paired_path)
  validate_psa_pair(paired_obj, psa_params)

  rerun <- "Rerun analysis/10_PSA.R before rendering."
  if (!is.null(seed)) {
    seed_ok <- if (exists("psa_samples_seed_matches", mode = "function",
                          inherits = TRUE)) {
      psa_samples_seed_matches(psa_params, seed)
    } else {
      identical(attr(psa_params, "seed"), seed)
    }
    if (!isTRUE(seed_ok)) {
      stop("Stale PSA parameter cache ", cache_file,
           ": RNG seed missing or mismatched (expected ", seed, "). ",
           rerun, call. = FALSE)
    }
  }
  if (!is.null(n_sim) && nrow(psa_params) != n_sim) {
    stop("Stale PSA parameter cache ", cache_file, ": ", nrow(psa_params),
         " rows but the PSA object has ", n_sim, " draws. ", rerun,
         call. = FALSE)
  }

  if (verbose) {
    cat("Loaded PSA parameter cache:", cache_file, "\n")
  }

  psa_params
}

#' Create a single-model CEAC plot
#'
#' @param psa_obj dampack PSA object.
#' @param wtp_range Vector of willingness-to-pay thresholds.
#' @param wtp_line Reference willingness-to-pay threshold.
#' @return ggplot object.
create_ceac_plot <- function(psa_obj,
                             wtp_range = seq(0, 100000, by = 2000),
                             wtp_line = NULL) {

  if (is.null(wtp_line)) {
    wtp_line <- if (exists("WTP", inherits = TRUE)) get("WTP", inherits = TRUE) else 51000
  }

  if (is.null(psa_obj)) {
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5,
                      label = "No CEAC data available", size = 5) +
             theme_void())
  }

  ceac_data <- dampack::ceac(wtp = wtp_range, psa = psa_obj)
  names(ceac_data)[names(ceac_data) == "Proportion"] <- "Probability"

  ggplot(ceac_data, aes(x = WTP, y = Probability, color = Strategy)) +
    geom_line(linewidth = 0.8) +
    geom_vline(xintercept = wtp_line, linetype = "dotted",
               color = "gray40", linewidth = 0.5) +
    scale_x_continuous(labels = scales::dollar_format(prefix = "EUR", scale = 0.001, suffix = "K")) +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
    labs(
      x = "Willingness-to-Pay Threshold (EUR/QALY)",
      y = "Probability of Being Cost-Effective",
      color = "Strategy"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
}

#' Create a PSA summary table for the single model
#'
#' @param psa_obj dampack PSA object.
#' @param wtp WTP threshold for probability of cost-effectiveness.
#' Intervals are the 2.5th/97.5th percentiles (type 7) of retained draws,
#' with increments calculated within each draw versus standard of care.
#' ICER is the ratio of mean increments, never the mean of draw-level ratios.
#' @return Numeric means, percentile limits, pairwise ICER/status and Prob_CE.
create_psa_summary_table <- function(psa_obj, wtp = NULL) {

  if (is.null(wtp)) {
    wtp <- if (exists("WTP", inherits = TRUE)) get("WTP", inherits = TRUE) else 51000
  }

  if (is.null(psa_obj)) {
    return(data.frame(
      Strategy = character(),
      Mean_Cost = numeric(),
      Mean_QALY = numeric(),
      Cost_Lower = numeric(), Cost_Upper = numeric(),
      QALY_Lower = numeric(), QALY_Upper = numeric(),
      Mean_Inc_Cost = numeric(), Mean_Inc_QALY = numeric(),
      Inc_Cost_Lower = numeric(), Inc_Cost_Upper = numeric(),
      Inc_QALY_Lower = numeric(), Inc_QALY_Upper = numeric(),
      ICER = numeric(), Status = character(),
      Prob_CE = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  cost_matrix <- as.matrix(psa_obj$cost)
  effect_matrix <- if (!is.null(psa_obj$effect)) {
    as.matrix(psa_obj$effect)
  } else {
    as.matrix(psa_obj$effectiveness)
  }
  control <- which(psa_obj$strategies == get_control_strategy())
  if (length(control) != 1L) stop("Control strategy must appear exactly once.")
  if (!identical(dim(cost_matrix), dim(effect_matrix)) ||
      ncol(cost_matrix) != length(psa_obj$strategies) || nrow(cost_matrix) < 1L ||
      anyDuplicated(psa_obj$strategies) ||
      any(!is.finite(cost_matrix)) || any(!is.finite(effect_matrix))) {
    stop("PSA summary requires aligned, finite cost/effect draws for every strategy.")
  }
  ceac_obj <- dampack::ceac(wtp = c(wtp - 1, wtp, wtp + 1), psa = psa_obj)

  summary_list <- lapply(seq_along(psa_obj$strategies), function(i) {
    strategy <- psa_obj$strategies[i]
    ceac_row <- ceac_obj[ceac_obj$WTP == wtp & ceac_obj$Strategy == strategy, ]
    prob_ce <- if (nrow(ceac_row) > 0) ceac_row$Proportion[1] else NA_real_
    inc_cost <- cost_matrix[, i] - cost_matrix[, control]
    inc_qaly <- effect_matrix[, i] - effect_matrix[, control]
    dc <- mean(inc_cost)
    de <- mean(inc_qaly)
    interval <- function(x) unname(quantile(x, c(0.025, 0.975), type = 7))
    cost_ci <- interval(cost_matrix[, i])
    qaly_ci <- interval(effect_matrix[, i])
    dc_ci <- interval(inc_cost)
    de_ci <- interval(inc_qaly)
    status <- if (i == control) "Reference" else if (dc >= 0 && de <= 0 && (dc > 0 || de < 0)) {
      "Dominated by SoC"
    } else if (dc <= 0 && de >= 0 && (dc < 0 || de > 0)) {
      "Cost-saving vs SoC"
    } else if (dc == 0 && de == 0) "Equivalent to SoC" else "Pairwise ICER vs SoC"

    data.frame(
      Strategy = strategy,
      Mean_Cost = mean(cost_matrix[, i]),
      Mean_QALY = mean(effect_matrix[, i]),
      Cost_Lower = cost_ci[1], Cost_Upper = cost_ci[2],
      QALY_Lower = qaly_ci[1], QALY_Upper = qaly_ci[2],
      Mean_Inc_Cost = dc, Mean_Inc_QALY = de,
      Inc_Cost_Lower = dc_ci[1], Inc_Cost_Upper = dc_ci[2],
      Inc_QALY_Lower = de_ci[1], Inc_QALY_Upper = de_ci[2],
      ICER = if (i == control || de == 0) NA_real_ else dc / de,
      Status = status,
      Prob_CE = prob_ce,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, summary_list)
}

#' Load scenario EVPPI results
#'
#' @param results_file Path to the .rds file containing scenario EVPPI results.
#' @return Scenario EVPPI result object, or NULL if unavailable.
load_scenario_evppi_results <- function(results_file = NULL) {

  if (is.null(results_file)) {
    results_file <- scenario_evppi_path()
  }

  if (!file.exists(results_file)) {
    warning("Scenario EVPPI results file not found: ", results_file)
    return(NULL)
  }

  readRDS(results_file)
}

#' Compare PSA strategy means with base-case values
#'
#' The PSA and the base case use the same joint survival model (issue #151),
#' but E[f(theta)] need not equal f(E[theta]) in a nonlinear model. This methods
#' diagnostic has no mean-alignment pass/fail threshold (issue #180). It
#' reports, per strategy and outcome, the base-case value, the PSA mean, its
#' Monte Carlo standard error, and the difference expressed in standard errors.
#'
#' @param psa_obj dampack PSA object.
#' @param base_results Data frame with Strategy, Cost, and Effect columns
#'   (the base-case output of \code{run_basecase()}).
#' @return Data frame with columns Strategy, Outcome, Base_Case, PSA_Mean,
#'   PSA_SE, Diff, Diff_SE, ordered by strategy then outcome (Cost, QALYs).
create_psa_basecase_comparison <- function(psa_obj, base_results) {
  if (is.null(psa_obj)) {
    return(data.frame(
      Strategy = character(), Outcome = character(),
      Base_Case = numeric(), PSA_Mean = numeric(), PSA_SE = numeric(),
      Diff = numeric(), Diff_SE = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  cost_matrix <- as.matrix(psa_obj$cost)
  effect_matrix <- if (!is.null(psa_obj$effect)) {
    as.matrix(psa_obj$effect)
  } else {
    as.matrix(psa_obj$effectiveness)
  }
  psa_strategies <- psa_obj$strategies
  missing <- setdiff(psa_strategies, base_results$Strategy)
  if (length(missing) > 0L) {
    stop("Base-case results lack strategies present in the PSA: ",
         paste(missing, collapse = ", "))
  }

  rows <- lapply(seq_along(psa_strategies), function(i) {
    strategy <- psa_strategies[i]
    base_row <- base_results[base_results$Strategy == strategy, , drop = FALSE]
    outcome_values <- list(
      Cost = list(base = base_row$Cost[1], draws = cost_matrix[, i]),
      QALYs = list(base = base_row$Effect[1], draws = effect_matrix[, i])
    )
    do.call(rbind, lapply(names(outcome_values), function(outcome) {
      draws <- outcome_values[[outcome]]$draws
      draws <- draws[is.finite(draws)]
      psa_mean <- mean(draws)
      psa_se <- stats::sd(draws) / sqrt(length(draws))
      diff <- psa_mean - outcome_values[[outcome]]$base
      data.frame(
        Strategy = strategy,
        Outcome = outcome,
        Base_Case = outcome_values[[outcome]]$base,
        PSA_Mean = psa_mean,
        PSA_SE = psa_se,
        Diff = diff,
        Diff_SE = diff / psa_se,
        stringsAsFactors = FALSE
      )
    }))
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Compare PSA incremental means (versus control) with base-case increments
#'
#' Companion to \code{create_psa_basecase_comparison()}. Because every PSA
#' draw predicts all strategies from one joint coefficient draw, a shift that is
#' common to every strategy cancels in the
#' increments; a control-specific shift, such as the separate age/sex control
#' model removed in issue #151, does not. The standard error is that of the
#' paired incremental draws, so it reflects the within-draw correlation.
#'
#' @param psa_obj dampack PSA object.
#' @param base_results Data frame with Strategy, Cost, and Effect columns.
#' @return Data frame with columns Comparator, Outcome, Base_Case, PSA_Mean,
#'   PSA_SE, Diff, Diff_SE, one row per biomarker strategy and outcome.
create_psa_incremental_comparison <- function(psa_obj, base_results) {
  control <- get_control_strategy()
  if (is.null(psa_obj)) {
    return(data.frame(
      Comparator = character(), Outcome = character(),
      Base_Case = numeric(), PSA_Mean = numeric(), PSA_SE = numeric(),
      Diff = numeric(), Diff_SE = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  cost_matrix <- as.matrix(psa_obj$cost)
  effect_matrix <- if (!is.null(psa_obj$effect)) {
    as.matrix(psa_obj$effect)
  } else {
    as.matrix(psa_obj$effectiveness)
  }
  colnames(cost_matrix) <- colnames(effect_matrix) <- psa_obj$strategies
  if (!control %in% psa_obj$strategies) {
    stop("PSA object has no control strategy '", control, "'")
  }
  base_control <- base_results[base_results$Strategy == control, , drop = FALSE]

  comparators <- setdiff(psa_obj$strategies, control)
  rows <- lapply(comparators, function(strategy) {
    base_row <- base_results[base_results$Strategy == strategy, , drop = FALSE]
    outcome_values <- list(
      Cost = list(
        base = base_row$Cost[1] - base_control$Cost[1],
        draws = cost_matrix[, strategy] - cost_matrix[, control]
      ),
      QALYs = list(
        base = base_row$Effect[1] - base_control$Effect[1],
        draws = effect_matrix[, strategy] - effect_matrix[, control]
      )
    )
    do.call(rbind, lapply(names(outcome_values), function(outcome) {
      draws <- outcome_values[[outcome]]$draws
      draws <- draws[is.finite(draws)]
      psa_mean <- mean(draws)
      psa_se <- stats::sd(draws) / sqrt(length(draws))
      diff <- psa_mean - outcome_values[[outcome]]$base
      data.frame(
        Comparator = strategy,
        Outcome = outcome,
        Base_Case = outcome_values[[outcome]]$base,
        PSA_Mean = psa_mean,
        PSA_SE = psa_se,
        Diff = diff,
        Diff_SE = diff / psa_se,
        stringsAsFactors = FALSE
      )
    }))
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

message("Single-model CEA helper functions loaded successfully.")
