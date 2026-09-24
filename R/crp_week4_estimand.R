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
  grid <- expand.grid(tmb_braf = 0:1, crp = 0:1)[, c("crp", "tmb_braf")]
  sign <- function(x) ifelse(x == 1, "+", "-")
  data.frame(
    id = paste0("crp", grid$crp, "_tmb", grid$tmb_braf),
    label = paste0("CRP", sign(grid$crp), ", TMB/BRAF", sign(grid$tmb_braf)),
    crp = grid$crp,
    tmb_braf = grid$tmb_braf,
    stringsAsFactors = FALSE
  )
}

#' Rows of patients without a week-4 CRP who have every other model variable
#'
#' Only patients whose covariates (other than the biomarkers) and both
#' endpoints are complete can be restored; the function stops otherwise, so a
#' restored patient is never silently dropped again.
missing_crp_rows <- function(data) {
  missing <- which(is.na(data$crp))
  other <- c("Age", "sex", "Rx", "OSwk", "Death", "PFSwk", "Progression")
  if (!all(complete.cases(data[missing, other, drop = FALSE]))) {
    stop("A patient without a week-4 CRP lacks another model variable.")
  }
  missing
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
  rows <- missing_crp_rows(data)
  assign_value <- function(x, target_rows, value) {
    if (is.factor(x)) {
      out <- as.character(x)
      out[target_rows] <- as.character(value)
      factor(out, levels = levels(x))
    } else {
      x[target_rows] <- value
      x
    }
  }
  unknown_tmb <- rows[is.na(data$tmb_braf[rows])]
  data$crp <- assign_value(data$crp, rows, crp_value)
  data$tmb_braf <- assign_value(data$tmb_braf, unknown_tmb, tmb_braf_value)
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
  converged <- vapply(fits, function(fit) {
    inherits(fit, "flexsurvreg") && identical(as.integer(fit$opt$convergence), 0L) &&
      !is.null(fit$cov) && !anyNA(fit$cov) && all(is.finite(fit$res[, "est"]))
  }, logical(1))
  list(models = fits, converged = converged)
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
  biomarkers <- get_biomarkers()
  status_of <- function(b) as.numeric(as.character(population[[b]]))
  prevalences <- vapply(biomarkers, function(b) mean(status_of(b)), numeric(1))
  arm <- as.integer(population$Rx) - 1L
  summary <- data.frame(
    n = nrow(population),
    n_control = sum(arm == 0L),
    n_experimental = sum(arm == 1L),
    deaths = sum(population$Death),
    pfs_events = sum(population$Progression),
    t(prevalences),
    stringsAsFactors = FALSE
  )
  fit <- fit_fixed_distribution_models(population, os_distribution, pfs_distribution)
  out <- list(summary = summary, models = fit$models, converged = fit$converged,
              ordering = NULL, params = NULL, results = NULL, pairwise = NULL,
              status = "ok")
  if (!all(fit$converged)) {
    out$status <- paste("fit failed:", paste(names(fit$converged)[!fit$converged],
                                             collapse = ", "))
    return(out)
  }
  strategies_df_s <- transform(get_strategy_metadata(get_strategies()),
    prevalence = c(1, prevalences)[match(get_strategies(),
                                          c(get_control_strategy(), biomarkers))],
    n_patients = nrow(population))
  predictions <- generate_population_averaged_predictions(
    models = fit$models, strategies_df = strategies_df_s,
    data_complete = population, time_points = time_points,
    prevalences = setNames(strategies_df_s$prevalence, strategies_df_s$id),
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
  base_population <- economic_prediction_population(data)
  restored <- missing_crp_rows(data)
  scenarios <- crp_missingness_scenarios()
  runs <- lapply(seq_len(nrow(scenarios)), function(i) {
    assigned <- assign_missing_crp(data, scenarios$crp[i], scenarios$tmb_braf[i])
    population <- economic_prediction_population(assigned)
    expected_ids <- union(base_population$ID, data$ID[restored])
    if (!setequal(population$ID, expected_ids) ||
        nrow(population) != nrow(base_population) + length(restored)) {
      stop("Scenario ", scenarios$id[i], " did not add exactly the patients ",
           "without a week-4 CRP to the complete-case cohort.")
    }
    run_cohort_model(population, base_params, os_distribution,
                     pfs_distribution, time_points)
  })
  names(runs) <- scenarios$id
  reproduction <- run_cohort_model(base_population, base_params, os_distribution,
                                   pfs_distribution, time_points)
  list(scenarios = scenarios, runs = runs, reproduction = reproduction,
       restored_rows = restored, os_distribution = os_distribution,
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
  guided <- setdiff(get_strategies(), get_control_strategy())
  out <- do.call(rbind, lapply(guided, function(strategy) {
    rows <- do.call(rbind, lapply(names(evaluated), function(id) {
      pw <- evaluated[[id]]$pairwise
      cbind(scenario = id, pw[pw$Strategy == strategy, ], stringsAsFactors = FALSE)
    }))
    finite <- rows[rows$Status == "Pairwise ICER vs SoC", , drop = FALSE]
    base_row <- base_pairwise[base_pairwise$Strategy == strategy, ]
    pick <- function(f) if (nrow(finite)) finite[f(finite$ICER), , drop = FALSE] else NULL
    lo <- pick(which.min)
    hi <- pick(which.max)
    data.frame(
      Strategy = strategy,
      n_evaluated = nrow(rows),
      n_finite = nrow(finite),
      icer_min = if (is.null(lo)) NA_real_ else lo$ICER,
      icer_min_scenario = if (is.null(lo)) NA_character_ else lo$scenario,
      icer_max = if (is.null(hi)) NA_real_ else hi$ICER,
      icer_max_scenario = if (is.null(hi)) NA_character_ else hi$scenario,
      n_dominated = sum(rows$Status == "Dominated by SoC"),
      n_cost_saving = sum(rows$Status == "Cost-saving vs SoC"),
      n_frontier = sum(rows$Frontier_Status == "ND"),
      base_icer = base_row$ICER,
      base_status = base_row$Status,
      stringsAsFactors = FALSE
    )
  }))
  optimal_of <- function(results) results$Strategy[which.max(results$Effect * wtp - results$Cost)]
  base_results <- base_pairwise[, c("Strategy", "Cost", "Effect")]
  attr(out, "optimal") <- c(base = optimal_of(base_results),
                            vapply(evaluated, function(run) optimal_of(run$results), character(1)))
  out
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
  before <- seq_len(decision_week)                  # grid points 0..decision_week-1
  at_decision <- decision_week + 1                  # grid point decision_week
  flow_share <- rep(0, n)
  flow_share[before] <- 1
  flow_share[at_decision] <- 0.5
  events <- seq(2, decision_week + 1)               # deaths in (0, decision_week]

  schedule <- list(
    nivo_before = sum(params$l_nivo[before]),
    flox_identical = identical(params$l_FLOX_exp[before], params$l_FLOX_control[before]),
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
    scheduled <- (o$drug_costs + o$test_costs + o$visit_costs) * weights$cost
    flows <- (o$follow_up_costs + o$post_progression_costs) * weights$cost
    events_cost <- o$end_life_costs * weights$cost
    c(qalys_pre = sum(o$qalys_discounted * flow_share),
      cost_pre = sum(scheduled[before]) + sum(flows * flow_share) + sum(events_cost[events]),
      diagnostic_pre = diagnostic,
      qalys_full = o$qalys_total,
      cost_full = o$costs_total)
  }

  rows <- lapply(get_strategies(), function(strategy) {
    trace <- run$traces[[strategy]]
    if (strategy == get_control_strategy()) {
      values <- split_outcome(trace, "standard", NULL)
    } else {
      prevalence <- trace$biomarker_prevalence
      values <- prevalence * split_outcome(trace$positive, "experimental", strategy) +
        (1 - prevalence) * split_outcome(trace$negative, "standard", strategy)
    }
    data.frame(Strategy = strategy, t(values), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  model_results <- run$results[match(out$Strategy, run$results$Strategy), ]
  check <- max(abs(c(out$qalys_full - model_results$Effect,
                     out$cost_full - model_results$Cost)))
  list(strategy = out, check = check, schedule = schedule,
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
  control <- x[x$Strategy == get_control_strategy(), ]
  guided <- x[x$Strategy != get_control_strategy(), ]
  out <- data.frame(
    Strategy = guided$Strategy,
    inc_qalys_pre = guided$qalys_pre - control$qalys_pre,
    inc_qalys_full = guided$qalys_full - control$qalys_full,
    inc_cost_pre = guided$cost_pre - control$cost_pre,
    inc_diagnostic_pre = guided$diagnostic_pre - control$diagnostic_pre,
    inc_cost_full = guided$cost_full - control$cost_full,
    stringsAsFactors = FALSE
  )
  out$inc_cost_pre_survival <- out$inc_cost_pre - out$inc_diagnostic_pre
  out$share_qalys <- out$inc_qalys_pre / out$inc_qalys_full
  out$share_cost <- out$inc_cost_pre / out$inc_cost_full
  out$share_cost_survival <- out$inc_cost_pre_survival / out$inc_cost_full
  out
}
