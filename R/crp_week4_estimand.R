# ===============================================================================
# WEEK-4 CRP ESTIMAND DIAGNOSTICS (issue #175)
# ===============================================================================
# Deterministic diagnostics for the estimand of the CRP-guided strategy. The
# economic model keeps its complete-case cohort and its randomisation-time
# survival contrast (option (c) of issue #175); these functions quantify both
# choices without changing any model code or cache:
#
#   1. Missingness sensitivity. The patients without a week-4 CRP are restored
#      to the complete-case cohort under each combination of assigned CRP
#      status and assigned TMB/BRAF status (the latter only where TMB/BRAF is
#      unknown). Each scenario refits the joint OS/PFS models with the pipeline
#      formulas and the base-case selected distributions, rebuilds the common
#      prediction population, prevalences and joint biomarker cells, checks
#      OS >= PFS, and runs the deterministic model.
#   2. Pre-decision accrual. Costs and QALYs accrued before the first
#      nivolumab dose, from the deterministic traces, using the trapezoidal
#      interval weights of calculate_outcomes().
#
# Requires (sourced by the report setup): scripts 02-05 (data, l_params_base,
# models), R/model_fun.R, R/calculate_outcomes.R, R/prediction_functions.R,
# R/para_model_fit.R and R/cea_helpers.R. Nothing here writes to the global
# environment or to any cache.
# ===============================================================================

#' Assignment scenarios for patients without a week-4 CRP
#'
#' @return Data frame with one row per scenario: `id`, `label`, the CRP value
#'   assigned to every restored patient and the TMB/BRAF value assigned to the
#'   restored patients whose TMB/BRAF status is unknown.
crp_missingness_scenarios <- function() {
  df_grid <- expand.grid(tmb_braf = 0:1, crp = 0:1)[, c("crp", "tmb_braf")]
  sign <- function(x) ifelse(x == 1, "+", "-")
  data.frame(
    id = paste0("crp", df_grid$crp, "_tmb", df_grid$tmb_braf),
    label = paste0("CRP", sign(df_grid$crp), ", TMB/BRAF", sign(df_grid$tmb_braf)),
    crp = df_grid$crp,
    tmb_braf = df_grid$tmb_braf,
    stringsAsFactors = FALSE
  )
}

#' Rows of patients without a week-4 CRP who have every other model variable
#'
#' Only patients whose covariates (other than the biomarkers) and both
#' endpoints are complete can be restored; the function stops otherwise, so a
#' restored patient is never silently dropped again.
missing_crp_rows <- function(data) {
  v_missing <- which(is.na(data$crp))
  v_other <- c("Age", "sex", "Rx", "OSwk", "Death", "PFSwk", "Progression")
  if (!all(complete.cases(data[v_missing, v_other, drop = FALSE]))) {
    stop("A patient without a week-4 CRP lacks another model variable.")
  }
  v_missing
}

#' Assign biomarker values to the patients without a week-4 CRP
#'
#' @param data Prepared trial data (after scripts 02-04).
#' @param crp_value CRP status (0/1) assigned to every restored patient.
#' @param tmb_braf_value TMB/BRAF status (0/1) assigned to restored patients
#'   whose TMB/BRAF status is unknown; known values are kept.
#' @return `data` with the assigned values, preserving the column classes and
#'   factor levels used by the fitted models.
assign_missing_crp <- function(data, crp_value, tmb_braf_value) {
  v_rows <- missing_crp_rows(data)
  assign_value <- function(x, target_rows, value) {
    if (is.factor(x)) {
      v_out <- as.character(x)
      v_out[target_rows] <- as.character(value)
      factor(v_out, levels = levels(x))
    } else {
      x[target_rows] <- value
      x
    }
  }
  v_unknown_tmb <- v_rows[is.na(data$tmb_braf[v_rows])]
  data$crp <- assign_value(data$crp, v_rows, crp_value)
  data$tmb_braf <- assign_value(data$tmb_braf, v_unknown_tmb, tmb_braf_value)
  data
}

#' Refit the joint OS/PFS models on a population with fixed distributions
#'
#' Uses fit_all_direct() (the script 04 fitting path) for the pipeline formulas
#' and the distributions already selected for the base case; distributions are
#' not re-selected.
fit_fixed_distribution_models <- function(population, os_distribution,
                                          pfs_distribution,
                                          formulas = get_model_formulas()$full) {
  fit_one <- function(formula, distribution) {
    fit_all_direct(fit_data = population, fit_formula = formula,
                   fit_dists = distribution)$fitted_models[[distribution]]
  }
  fits <- list(os = fit_one(formulas$os, os_distribution),
               pfs = fit_one(formulas$pfs, pfs_distribution))
  v_converged <- vapply(fits, function(fit) {
    inherits(fit, "flexsurvreg") && identical(as.integer(fit$opt$convergence), 0L) &&
      !is.null(fit$cov) && !anyNA(fit$cov) && all(is.finite(fit$res[, "est"]))
  }, logical(1))
  list(models = fits, converged = v_converged)
}

#' Population, fits, ordering check and deterministic results for one cohort
#'
#' @param population Complete-case prediction population to fit and predict on.
#' @param base_params Base-case parameter list (l_params_base); a modified copy
#'   is returned, the argument is not changed.
#' @param os_distribution,pfs_distribution Selected base-case distributions.
#' @param time_points Weekly prediction grid (0, ..., time_horizon).
#' @return List with the population summary, fitted models, convergence flags,
#'   ordering check, parameter list and deterministic results (NULL when the
#'   fit failed or OS >= PFS is violated, in which case `status` says why).
run_cohort_model <- function(population, base_params, os_distribution,
                             pfs_distribution, time_points) {
  v_biomarkers <- get_biomarkers()
  status_of <- function(b) as.numeric(as.character(population[[b]]))
  v_prevalences <- vapply(v_biomarkers, function(b) mean(status_of(b)), numeric(1))
  v_arm <- as.integer(population$Rx) - 1L
  df_summary <- data.frame(
    n = nrow(population),
    n_control = sum(v_arm == 0L),
    n_experimental = sum(v_arm == 1L),
    deaths = sum(population$Death),
    pfs_events = sum(population$Progression),
    t(v_prevalences),
    stringsAsFactors = FALSE
  )
  fit <- fit_fixed_distribution_models(population, os_distribution, pfs_distribution)
  out <- list(summary = df_summary, models = fit$models, converged = fit$converged,
              ordering = NULL, params = NULL, results = NULL, pairwise = NULL,
              status = "ok")
  if (!all(fit$converged)) {
    out$status <- paste("fit failed:", paste(names(fit$converged)[!fit$converged],
                                             collapse = ", "))
    return(out)
  }
  df_strategies_s <- transform(get_strategy_metadata(get_strategies()),
    prevalence = c(1, v_prevalences)[match(get_strategies(),
                                          c(get_control_strategy(), v_biomarkers))],
    n_patients = nrow(population))
  predictions <- generate_population_averaged_predictions(
    models = fit$models, strategies_df = df_strategies_s,
    data_complete = population, time_points = time_points,
    prevalences = setNames(df_strategies_s$prevalence, df_strategies_s$id),
    quiet = TRUE)
  out$ordering <- check_population_survival_ordering(predictions)
  if (!out$ordering$ordered) {
    out$status <- sprintf("OS >= PFS violated at %s curve-time points",
                          out$ordering$n_violations)
    return(out)
  }
  # Parameter list as script 05 builds it: schedules, prices and utilities are
  # unchanged; curves, prevalences, population and joint cells are replaced.
  params <- set_population_predictions(base_params, predictions)
  joint <- joint_biomarker_population(population)
  params$joint_counts <- joint$counts
  for (cell in names(joint$probabilities)) {
    params[[paste0("p_joint_", cell)]] <- unname(joint$probabilities[cell])
  }
  out$params <- params
  out$results <- run_basecase(params, verbose = FALSE)$base_results
  out$pairwise <- calculate_pairwise_icers(out$results)
  out
}

#' Missingness sensitivity analysis for the patients without a week-4 CRP
#'
#' @param data Prepared trial data after scripts 02-04 (biomarkers as factors).
#' @param base_params Base-case parameter list (l_params_base).
#' @param models Script 04 `models` object; the selected distributions are
#'   read from `models$best_fit`.
#' @param time_points Weekly prediction grid.
#' @return List with `scenarios` (definitions), `runs` (one run_cohort_model()
#'   result per scenario), `reproduction` (the same code applied to the
#'   unmodified complete-case cohort) and `restored_rows`.
run_crp_missingness_sensitivity <- function(data, base_params, models, time_points) {
  os_distribution <- models$best_fit$os_distribution
  pfs_distribution <- models$best_fit$pfs_distribution
  df_base_population <- economic_prediction_population(data)
  v_restored <- missing_crp_rows(data)
  df_scenarios <- crp_missingness_scenarios()
  runs <- lapply(seq_len(nrow(df_scenarios)), function(i) {
    df_assigned <- assign_missing_crp(data, df_scenarios$crp[i], df_scenarios$tmb_braf[i])
    df_population <- economic_prediction_population(df_assigned)
    v_expected_ids <- union(df_base_population$ID, data$ID[v_restored])
    if (!setequal(df_population$ID, v_expected_ids) ||
        nrow(df_population) != nrow(df_base_population) + length(v_restored)) {
      stop("Scenario ", df_scenarios$id[i], " did not add exactly the patients ",
           "without a week-4 CRP to the complete-case cohort.")
    }
    run_cohort_model(df_population, base_params, os_distribution,
                     pfs_distribution, time_points)
  })
  names(runs) <- df_scenarios$id
  reproduction <- run_cohort_model(df_base_population, base_params, os_distribution,
                                   pfs_distribution, time_points)
  list(scenarios = df_scenarios, runs = runs, reproduction = reproduction,
       restored_rows = v_restored, os_distribution = os_distribution,
       pfs_distribution = pfs_distribution)
}

#' Range of the scenario results for each guided strategy
#'
#' @param sensitivity Result of run_crp_missingness_sensitivity().
#' @param base_pairwise calculate_pairwise_icers() of the deterministic base case.
#' @param wtp Willingness-to-pay threshold for the highest-NMB strategy.
#' @return Data frame with one row per guided strategy: the number of evaluated
#'   scenarios, the smallest and largest finite pairwise ICER and their
#'   scenario IDs, the number of scenarios per pairwise status and on the
#'   frontier, and the base-case ICER and status. Attribute `optimal` holds the
#'   highest-NMB strategy of the base case and of every evaluated scenario.
summarise_crp_missingness <- function(sensitivity, base_pairwise, wtp) {
  evaluated <- Filter(function(run) !is.null(run$pairwise), sensitivity$runs)
  v_guided <- setdiff(get_strategies(), get_control_strategy())
  df_out <- do.call(rbind, lapply(v_guided, function(strategy) {
    df_rows <- do.call(rbind, lapply(names(evaluated), function(id) {
      df_pw <- evaluated[[id]]$pairwise
      cbind(scenario = id, df_pw[df_pw$Strategy == strategy, ], stringsAsFactors = FALSE)
    }))
    df_finite <- df_rows[df_rows$Status == "Pairwise ICER vs SoC", , drop = FALSE]
    df_base_row <- base_pairwise[base_pairwise$Strategy == strategy, ]
    pick <- function(f) if (nrow(df_finite)) df_finite[f(df_finite$ICER), , drop = FALSE] else NULL
    df_lo <- pick(which.min)
    df_hi <- pick(which.max)
    data.frame(
      Strategy = strategy,
      n_evaluated = nrow(df_rows),
      n_finite = nrow(df_finite),
      icer_min = if (is.null(df_lo)) NA_real_ else df_lo$ICER,
      icer_min_scenario = if (is.null(df_lo)) NA_character_ else df_lo$scenario,
      icer_max = if (is.null(df_hi)) NA_real_ else df_hi$ICER,
      icer_max_scenario = if (is.null(df_hi)) NA_character_ else df_hi$scenario,
      n_dominated = sum(df_rows$Status == "Dominated by SoC"),
      n_cost_saving = sum(df_rows$Status == "Cost-saving vs SoC"),
      n_frontier = sum(df_rows$Frontier_Status == "ND"),
      base_icer = df_base_row$ICER,
      base_status = df_base_row$Status,
      stringsAsFactors = FALSE
    )
  }))
  optimal_of <- function(results) results$Strategy[which.max(results$Effect * wtp - results$Cost)]
  df_base_results <- base_pairwise[, c("Strategy", "Cost", "Effect")]
  attr(df_out, "optimal") <- c(base = optimal_of(df_base_results),
                            vapply(evaluated, function(run) optimal_of(run$results), character(1)))
  df_out
}

#' Costs and QALYs accrued before a decision week, from deterministic traces
#'
#' Recomputes every subgroup with calculate_outcomes() from the traces of
#' model_fun(return_traces = TRUE) and splits each component at `decision_week`:
#'   - QALYs and ongoing progressed-state costs (trapezoidal interval flows):
#'     the intervals [0, 1], ..., [decision_week - 1, decision_week] weeks, i.e.
#'     the full trapezoidal weight of grid points 0, ..., decision_week - 1 and
#'     half the weight of grid point decision_week;
#'   - scheduled and one-time charges (drugs, tests, visits, baseline and
#'     diagnostic test): grid points 0, ..., decision_week - 1; the charges at
#'     grid point decision_week (the first nivolumab dose) are excluded;
#'   - end-of-life event costs: deaths in the intervals ending at weeks
#'     1, ..., decision_week (charged at those grid points).
#'
#' @return List with `strategy` (per-strategy pre-decision and full-horizon
#'   QALYs and costs, with the diagnostic-test cost separated), `check` (the
#'   largest absolute difference between the recomputed full-horizon totals and
#'   model_fun()) and `schedule` (checks that the schedules are identical before
#'   `decision_week` and that nivolumab starts at `decision_week`).
pre_decision_outcomes <- function(params, decision_week = 4,
                                  time_horizon = params$time_horizon,
                                  cl = params$cl) {
  n <- time_horizon + 1
  if (decision_week < 1 || decision_week >= time_horizon) {
    stop("decision_week must lie strictly inside the time horizon.")
  }
  v_before <- seq_len(decision_week)                  # grid points 0..decision_week-1
  at_decision <- decision_week + 1                  # grid point decision_week
  v_flow_share <- rep(0, n)
  v_flow_share[v_before] <- 1
  v_flow_share[at_decision] <- 0.5
  v_events <- seq(2, decision_week + 1)               # deaths in (0, decision_week]

  schedule <- list(
    nivo_before = sum(params$l_nivo[v_before]),
    flox_identical = identical(params$l_FLOX_exp[v_before], params$l_FLOX_control[v_before]),
    first_nivo_week = which(params$l_nivo == 1)[1] - 1
  )

  run <- model_fun(params, time_horizon = time_horizon, cl = cl,
                   determpsa = "det", return_traces = TRUE)
  weights <- discount_weights(params, n, cl)

  split_outcome <- function(states, treatment_type, biomarker) {
    o <- calculate_outcomes(params, states$p_pf, states$p_p, states$p_d,
                            treatment_type, biomarker, weights$cost,
                            weights$effect, cl)
    diagnostic <- if (is.null(biomarker)) 0 else params[[biomarker_cost_key(biomarker)[[1]]]]
    v_scheduled <- (o$drug_costs + o$test_costs + o$visit_costs) * weights$cost
    v_flows <- (o$follow_up_costs + o$progressed_imaging_costs +
               o$post_progression_costs) * weights$cost
    v_events_cost <- o$end_life_costs * weights$cost
    c(qalys_pre = sum(o$qalys_discounted * v_flow_share),
      cost_pre = sum(v_scheduled[v_before]) + sum(v_flows * v_flow_share) + sum(v_events_cost[v_events]),
      diagnostic_pre = diagnostic,
      qalys_full = o$qalys_total,
      cost_full = o$costs_total)
  }

  rows <- lapply(get_strategies(), function(strategy) {
    trace <- run$traces[[strategy]]
    if (strategy == get_control_strategy()) {
      v_values <- split_outcome(trace, "standard", NULL)
    } else {
      prevalence <- trace$biomarker_prevalence
      v_values <- prevalence * split_outcome(trace$positive, "experimental", strategy) +
        (1 - prevalence) * split_outcome(trace$negative, "standard", strategy)
    }
    data.frame(Strategy = strategy, t(v_values), stringsAsFactors = FALSE)
  })
  df_out <- do.call(rbind, rows)
  df_model_results <- run$results[match(df_out$Strategy, run$results$Strategy), ]
  check <- max(abs(c(df_out$qalys_full - df_model_results$Effect,
                     df_out$cost_full - df_model_results$Cost)))
  list(strategy = df_out, check = check, schedule = schedule,
       decision_week = decision_week)
}

#' Increments of pre-decision and full-horizon outcomes versus standard of care
#'
#' @param pre Result of pre_decision_outcomes().
#' @return Data frame for each guided strategy: pre-decision and full-horizon
#'   increments, the diagnostic-test and survival-dependent parts of the
#'   pre-decision cost increment, and each pre-decision increment as a share of
#'   the full-horizon increment.
pre_decision_increments <- function(pre) {
  x <- pre$strategy
  df_control <- x[x$Strategy == get_control_strategy(), ]
  df_guided <- x[x$Strategy != get_control_strategy(), ]
  df_out <- data.frame(
    Strategy = df_guided$Strategy,
    inc_qalys_pre = df_guided$qalys_pre - df_control$qalys_pre,
    inc_qalys_full = df_guided$qalys_full - df_control$qalys_full,
    inc_cost_pre = df_guided$cost_pre - df_control$cost_pre,
    inc_diagnostic_pre = df_guided$diagnostic_pre - df_control$diagnostic_pre,
    inc_cost_full = df_guided$cost_full - df_control$cost_full,
    stringsAsFactors = FALSE
  )
  df_out$inc_cost_pre_survival <- df_out$inc_cost_pre - df_out$inc_diagnostic_pre
  df_out$share_qalys <- df_out$inc_qalys_pre / df_out$inc_qalys_full
  df_out$share_cost <- df_out$inc_cost_pre / df_out$inc_cost_full
  df_out$share_cost_survival <- df_out$inc_cost_pre_survival / df_out$inc_cost_full
  df_out
}
