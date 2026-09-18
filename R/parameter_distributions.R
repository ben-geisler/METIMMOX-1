# ===============================================================================
# NON-SURVIVAL UNCERTAINTY PARAMETERS
# ===============================================================================
# The table below is the single source of truth for PSA, DSA, and EVPPI
# parameter membership. Add a parameter once here; downstream analyses derive
# their parameter vectors from the specification.
#
# Columns (issue #154):
#   psa      TRUE if the parameter is sampled in (or derived within) the PSA.
#   dsa      TRUE if the parameter is varied in the one-way DSA. Published
#            tariffs and list prices are psa = FALSE, dsa = TRUE: their value is
#            known, so sampling them adds decision uncertainty that does not
#            exist and puts "value of research on the FLOX tariff" at the top of
#            the EVPPI ranking, but their level still matters and is tested
#            deterministically.
#   derived  TRUE if the PSA column is a deterministic function of other sampled
#            parameters (apply_derived_psa_parameters()) rather than a draw.
#            Derived parameters are model inputs, so they take part in the DSA,
#            but they are excluded from EVPPI groups because they are exactly
#            collinear with the parameters they are derived from.

parameter_distribution_spec <- function() {
  biomarker_cost_parameters <- unique(unname(biomarker_cost_key()))
  prevalence_parameters <- unname(biomarker_prevalence_key())

  # Published list prices and laboratory/procedure tariffs: fixed in the PSA,
  # varied in the one-way DSA only.
  fixed_costs <- data.frame(
    parameter = c("c_drug_nivo", "c_drug_FLOX", "c_test_CT", "c_test_blood",
                  biomarker_cost_parameters),
    distribution = NA_character_,
    cv = NA_real_,
    group = c("drug_costs", "drug_costs", "test_costs", "test_costs",
              rep("test_costs", length(biomarker_cost_parameters))),
    psa = FALSE,
    order = c(10:11, 20:21, 22 + seq_along(biomarker_cost_parameters)),
    stringsAsFactors = FALSE
  )
  # Visit, baseline, follow-up and end-of-life costs bundle genuine resource-use
  # uncertainty (contact frequency, terminal-care intensity), so they remain
  # probabilistic.
  other_costs <- data.frame(
    parameter = c("c_other_visit", "c_other_baseline", "c_other_follow",
                  "c_other_last"),
    distribution = "gamma",
    cv = 0.20,
    group = "other_costs",
    psa = TRUE,
    order = 30:33,
    stringsAsFactors = FALSE
  )
  # u_p is not sampled directly. Independent beta draws for u_np and u_p put
  # 16% of draws (48% under UTILITY_SOURCE = 0) in the implausible region
  # u_p > u_np. Sampling a non-negative decrement instead makes u_p <= u_np hold
  # by construction, and E[u_decrement] = u_np - u_p leaves the marginal mean of
  # u_p unchanged.
  utilities <- data.frame(
    parameter = c("u_np", "u_decrement", "u_p"),
    distribution = c("beta", "gamma", "derived"),
    cv = c(0.15, 0.15, NA_real_),
    group = "utilities",
    psa = TRUE,
    order = 40:42,
    stringsAsFactors = FALSE
  )
  prevalence <- data.frame(
    parameter = c(paste0("p_joint_", c("00", "01", "10", "11")), prevalence_parameters),
    distribution = c(rep("dirichlet", 3), rep("derived", 1 + length(prevalence_parameters))),
    cv = NA_real_,
    group = "prevalence",
    psa = TRUE,
    order = 50 + seq_len(4 + length(prevalence_parameters)),
    stringsAsFactors = FALSE
  )

  spec <- rbind(fixed_costs, other_costs, utilities, prevalence)
  spec$derived <- spec$distribution %in% "derived"
  # DSA varies marginal prevalences, not individual joint masses. The sampled
  # utility decrement is also excluded because DSA varies state utilities.
  spec$dsa <- spec$parameter != "u_decrement" & !grepl("^p_joint_", spec$parameter)
  spec <- spec[order(spec$order), , drop = FALSE]
  spec$order <- NULL
  rownames(spec) <- NULL
  if (anyDuplicated(spec$parameter)) {
    stop("Parameter distribution specification contains duplicate names.")
  }
  spec
}

#' Look up the reporting group of any specified parameter
#'
#' Covers fixed as well as sampled parameters, so deterministic-analysis output
#' can be grouped even for parameters that have no PSA distribution.
#'
#' @param spec Parameter specification.
#' @return Named character vector mapping parameter to group.
parameter_group_lookup <- function(spec = parameter_distribution_spec()) {
  stats::setNames(spec$group, spec$parameter)
}

distribution_parameters <- function(mean_value, distribution, cv) {
  if (identical(distribution, "derived")) {
    return(list(dist = "derived"))
  }
  if (!is.finite(mean_value) || mean_value < 0) {
    stop("Distribution means must be finite and non-negative.")
  }
  if (identical(distribution, "gamma")) {
    shape <- 1 / cv^2
    return(list(dist = "gamma", shape = shape, rate = shape / mean_value))
  }
  if (identical(distribution, "beta")) {
    if (mean_value <= 0 || mean_value >= 1) {
      stop("Beta means must lie strictly between zero and one.")
    }
    variance <- (mean_value * cv)^2
    common <- mean_value * (1 - mean_value) / variance - 1
    if (common <= 0) {
      stop("Beta means and CVs must define positive shape parameters.")
    }
    return(list(
      dist = "beta",
      shape1 = mean_value * common,
      shape2 = (1 - mean_value) * common
    ))
  }
  stop("Unsupported parameter distribution: ", distribution)
}

create_parameter_distributions <- function(params,
                                           spec = parameter_distribution_spec()) {
  spec <- spec[spec$psa, , drop = FALSE]
  missing_parameters <- setdiff(spec$parameter, names(params))
  if (length(missing_parameters) > 0) {
    stop("Missing base-case parameters: ",
         paste(missing_parameters, collapse = ", "))
  }

  distributions <- lapply(seq_len(nrow(spec)), function(i) {
    if (spec$distribution[i] == "dirichlet") {
      counts <- params$joint_counts
      if (is.null(counts) || !identical(names(counts), c("00", "01", "10", "11")) ||
          any(!is.finite(counts)) || any(counts < 0) || sum(counts) <= 0)
        stop("Joint prevalence uncertainty requires complete-case joint cell counts.")
      return(list(dist = "dirichlet", alpha = counts,
                  component = sub("^p_joint_", "", spec$parameter[i])))
    }
    distribution_parameters(
      params[[spec$parameter[i]]], spec$distribution[i], spec$cv[i]
    )
  })
  stats::setNames(distributions, spec$parameter)
}

create_parameter_groups <- function(spec = parameter_distribution_spec()) {
  # EVPPI groups describe what the PSA actually samples. Fixed prices carry no
  # decision uncertainty, and derived columns are collinear with their inputs.
  spec <- spec[spec$psa & !spec$derived, , drop = FALSE]
  grouped <- split(spec$parameter, factor(
    spec$group,
    levels = c("drug_costs", "test_costs", "other_costs", "utilities", "prevalence")
  ))
  grouped <- grouped[lengths(grouped) > 0]

  cost_group_names <- c("drug_costs", "test_costs", "other_costs")
  sampled_cost_groups <- intersect(cost_group_names, names(grouped))
  # "all_costs" only earns a row when it spans more than one cost group;
  # otherwise it would duplicate that group in every EVPPI table.
  if (length(sampled_cost_groups) > 1) {
    grouped <- append(
      grouped,
      list(all_costs = spec$parameter[spec$group %in% cost_group_names]),
      after = match(utils::tail(sampled_cost_groups, 1), names(grouped))
    )
  }
  grouped
}

#' Fill in PSA columns that are deterministic functions of sampled parameters
#'
#' Applied by generate_psa_samples() after the sampling loop, so every consumer
#' of the PSA parameter table sees the same derived values.
#'
#' @param samples Data frame of sampled PSA parameters.
#' @return The data frame with derived columns added.
apply_derived_psa_parameters <- function(samples) {
  joint_keys <- paste0("p_joint_", c("00", "01", "10"))
  if (all(joint_keys %in% names(samples))) {
    samples$p_joint_11 <- pmax(0, 1 - rowSums(samples[joint_keys]))
    samples$p_crp <- samples$p_joint_10 + samples$p_joint_11
    samples$p_tmb_braf <- samples$p_joint_01 + samples$p_joint_11
  }
  if (all(c("u_np", "u_decrement") %in% names(samples))) {
    # Non-negative decrement, floored at zero: u_p <= u_np in every draw.
    samples$u_p <- pmax(samples$u_np - samples$u_decrement, 0)
  }
  samples
}

configure_parameter_distributions <- function(params) {
  spec <- parameter_distribution_spec()
  list(
    distributions = create_parameter_distributions(params, spec),
    groups = create_parameter_groups(spec),
    spec = spec
  )
}

# ===============================================================================
# DETERMINISTIC SENSITIVITY-ANALYSIS RANGES AND STRUCTURAL SCENARIOS
# ===============================================================================
# Single source of truth for the one-way DSA bounds (09_DSA.R, figure2.qmd) and
# for the input-parameter table (table_1.qmd), so the published ranges always
# correspond to an analysis that is actually run (issue #153).

#' One-way DSA bounds: base value +/- mult, capped to the parameter's support
#'
#' The DSA set is spec$dsa, not the PSA parameter set: unit prices that are
#' fixed in the PSA are still varied deterministically (issue #154).
#'
#' @param params Base-case parameter list (l_params_base).
#' @param spec Parameter specification from parameter_distribution_spec().
#' @param mult Proportional variation (default DSA_mult, +/-20%).
#' @return Data frame with columns pars, min, max.
build_dsa_ranges <- function(params, spec = parameter_distribution_spec(),
                             mult = DSA_mult) {
  spec <- spec[spec$dsa, , drop = FALSE]
  pars <- spec$parameter
  missing_parameters <- setdiff(pars, names(params))
  if (length(missing_parameters) > 0) {
    stop("Missing base-case parameters for the DSA: ",
         paste(missing_parameters, collapse = ", "))
  }
  ranges <- data.frame(
    pars = pars,
    min = unlist(params[pars]) * (1 - mult),
    max = unlist(params[pars]) * (1 + mult),
    stringsAsFactors = FALSE
  )
  cost_groups <- c("drug_costs", "test_costs", "other_costs")
  bounded_zero <- spec$group %in% c(cost_groups, "prevalence")
  bounded_one  <- spec$group %in% c("utilities", "prevalence")
  ranges$min[bounded_zero] <- pmax(ranges$min[bounded_zero], 0)
  ranges$max[bounded_one]  <- pmin(ranges$max[bounded_one], 1)
  if (any(ranges$min > ranges$max)) {
    stop("Invalid DSA ranges (min > max) for: ",
         paste(ranges$pars[ranges$min > ranges$max], collapse = ", "))
  }
  rownames(ranges) <- NULL
  ranges
}

#' Deterministic structural scenarios run alongside the one-way DSA
#'
#' Discount rate, model time horizon, post-progression treatment cost and the
#' second treatment sequence. The values here are what 09_DSA.R runs and what
#' table_1.qmd prints as ranges. `sides` names the endpoints that are actually
#' run; a one-sided scenario is one whose other endpoint is the base case.
#'
#' @param base_dr Base-case annual discount rate.
#' @param base_horizon_years Base-case time horizon in years.
#' @param base_pp_cost Base-case post-progression cost per quarter.
#' @return Named list; each element has label, unit, base, min, max, sides.
dsa_structural_scenarios <- function(base_dr, base_horizon_years,
                                     base_pp_cost = 0) {
  list(
    Discount_rate = list(
      label = "Discount rate (costs and QALYs)", unit = "rate",
      base = base_dr, min = 0, max = 0.08, sides = c("min", "max")
    ),
    Time_horizon = list(
      label = "Time horizon", unit = "years",
      base = base_horizon_years, min = 5, max = 20, sides = c("min", "max")
    ),
    # The base case costs nothing after progression beyond the quarterly
    # follow-up contact and the end-of-life cost (issue #154). EUR 5,000 per
    # quarter is an illustrative upper bound spanning a second-line systemic
    # therapy course plus imaging and visits; it is a scenario value, not a
    # costed Norwegian second-line pathway.
    Post_progression_cost = list(
      label = "Post-progression treatment cost (per quarter, progressed state)",
      unit = "EUR", base = base_pp_cost, min = base_pp_cost, max = 5000,
      sides = "max"
    ),
    # The base case gives every progression-free patient a second eight-cycle
    # sequence at weeks 24-38. The scenario removes it, holding survival at the
    # trial estimate, so it bounds the cost side of that assumption only.
    Second_sequence = list(
      label = "Second treatment sequence at weeks 24-38 (1 = given, 0 = omitted)",
      unit = "indicator", base = 1, min = 0, max = 1, sides = "min"
    )
  )
}

#' Endpoints of a structural scenario that 09_DSA.R runs
#'
#' @param scenario One element of dsa_structural_scenarios().
#' @return Character vector of endpoint names ("min", "max").
structural_scenario_sides <- function(scenario) {
  if (is.null(scenario$sides)) c("min", "max") else scenario$sides
}

# ---------------------------------------------------------------------------
# Applying a structural scenario to a base-case parameter list (issue #157)
# ---------------------------------------------------------------------------
# Moved here from 09_DSA.R so that the one-way DSA (09_DSA.R, OWSA.qmd,
# Table 1) and the deterministic scenario table (table_s8.qmd) run exactly the
# same scenario code. The discount rate is applied to costs and QALYs together.
# A different time horizon re-predicts the base-case survival curves on the
# new weekly grid from the selected joint models and rebuilds the schedule
# vectors with the same rules as script 05 (drug positions are fixed calendar
# weeks; CT every 12 weeks and blood tests every 4 weeks recur to the end of
# the horizon). The post-progression cost and second-sequence scenarios test
# the two scope assumptions raised in issue #154.

#' Rebuild the base-case parameter list for a different model horizon
#'
#' @param base_params Base-case parameter list (l_params_base).
#' @param horizon_weeks New horizon in weeks.
#' @param models list(os = <flexsurvreg>, pfs = <flexsurvreg>): the selected
#'   ordering-constrained best-fit models (models$best_fit$os / $pfs).
#' @param strategies_df Strategy table with `id` and `prevalence`.
#' @param data_complete Complete-case data used for population averaging.
#' @return Parameter list with curves and schedules on the new grid.
build_horizon_params <- function(base_params, horizon_weeks, models,
                                 strategies_df, data_complete) {
  if (is.null(models) || is.null(models$os) || is.null(models$pfs)) {
    stop("Fitted best-fit survival models are required for a time-horizon scenario.")
  }
  tp <- seq(0, horizon_weeks)
  n <- length(tp)
  preds <- generate_population_averaged_predictions(
    models = list(os = models$os, pfs = models$pfs),
    strategies_df = strategies_df,
    data_complete = data_complete,
    time_points = tp,
    prevalences = setNames(strategies_df$prevalence, strategies_df$id)
  )

  p <- set_population_predictions(base_params, preds)
  ctrl <- get_control_strategy()
  p$p_os  <- setNames(list(preds[[ctrl]]$os),  paste0(ctrl, "_OS"))
  p$p_pfs <- setNames(list(preds[[ctrl]]$pfs), paste0(ctrl, "_PFS"))
  for (bm in get_biomarkers()) {
    p$p_os[[paste0(bm, "_pos_OS")]]       <- preds[[bm]]$biomarker_positive$os
    p$p_os[[paste0(bm, "_neg_OS")]]       <- preds[[bm]]$biomarker_negative$os
    p$p_os[[paste0(bm, "_weighted_OS")]]  <- preds[[bm]]$os
    p$p_pfs[[paste0(bm, "_pos_PFS")]]     <- preds[[bm]]$biomarker_positive$pfs
    p$p_pfs[[paste0(bm, "_neg_PFS")]]     <- preds[[bm]]$biomarker_negative$pfs
    p$p_pfs[[paste0(bm, "_weighted_PFS")]] <- preds[[bm]]$pfs
  }

  schedule <- function(positions) {
    v <- rep(0, n)
    v[positions[positions <= n]] <- 1
    v
  }
  p$l_nivo         <- schedule(which(base_params$l_nivo == 1))
  p$l_FLOX_exp     <- schedule(which(base_params$l_FLOX_exp == 1))
  p$l_FLOX_control <- schedule(which(base_params$l_FLOX_control == 1))
  p$l_CT    <- schedule(c(1, seq(13, n, by = 12)))
  p$l_blood <- schedule(c(1, seq(5, n, by = 4)))
  p$l_visit <- as.numeric(p$l_nivo == 1 | p$l_FLOX_exp == 1 | p$l_FLOX_control == 1)
  p$l_visit[1] <- 1
  p$time_horizon <- horizon_weeks
  p
}

#' Remove the second treatment sequence
#'
#' Zero every administration from the given schedule position onwards and
#' rebuild the visit schedule from the remaining administrations. Monitoring
#' (CT, blood tests) is unchanged, and so is survival, so the scenario bounds
#' the cost of the assumption only.
#'
#' @param p Parameter list.
#' @param first_position First schedule position of the second sequence
#'   (position 25 = modeled week 24).
drop_second_sequence <- function(p, first_position = 25L) {
  zap <- function(v) {
    if (length(v) >= first_position) v[seq.int(first_position, length(v))] <- 0
    v
  }
  p$l_nivo         <- zap(p$l_nivo)
  p$l_FLOX_exp     <- zap(p$l_FLOX_exp)
  p$l_FLOX_control <- zap(p$l_FLOX_control)
  p$l_visit <- as.numeric(p$l_nivo == 1 | p$l_FLOX_exp == 1 | p$l_FLOX_control == 1)
  p$l_visit[1] <- 1
  p
}

#' Apply one endpoint of a structural scenario to a parameter list
#'
#' @param scenario_name Name of an entry of dsa_structural_scenarios().
#' @param value The endpoint value (scenario$min or scenario$max).
#' @param base_params Base-case parameter list.
#' @param base_horizon Base-case horizon in weeks (passed back unchanged for
#'   every scenario except Time_horizon).
#' @param models,strategies_df,data_complete Needed for Time_horizon only; see
#'   build_horizon_params().
#' @return list(params, time_horizon) ready for
#'   model_fun(params, time_horizon = time_horizon).
apply_structural_scenario <- function(scenario_name, value, base_params,
                                      base_horizon, models = NULL,
                                      strategies_df = NULL,
                                      data_complete = NULL) {
  p <- base_params
  horizon <- base_horizon
  if (scenario_name == "Discount_rate") {
    p$dr_costs <- value
    p$dr_effects <- value
  } else if (scenario_name == "Time_horizon") {
    horizon <- as.integer(round(value * 52))
    p <- build_horizon_params(base_params, horizon, models, strategies_df,
                              data_complete)
  } else if (scenario_name == "Post_progression_cost") {
    p$c_other_pp <- value
  } else if (scenario_name == "Second_sequence") {
    if (value == 0) p <- drop_second_sequence(p)
  } else {
    stop("Unknown structural scenario: ", scenario_name)
  }
  list(params = p, time_horizon = horizon)
}

#' Human-readable label for one structural scenario endpoint
#'
#' @param scenario_name Name of an entry of dsa_structural_scenarios().
#' @param value The endpoint value.
#' @return Character label used in Table S8 and the OWSA report.
structural_scenario_label <- function(scenario_name, value) {
  switch(scenario_name,
    Discount_rate = sprintf("Discount rate %g%% (costs and QALYs)", 100 * value),
    Time_horizon = sprintf("Time horizon %g years", value),
    Post_progression_cost = sprintf("Post-progression cost EUR %s per quarter",
                                    format(value, big.mark = ",", scientific = FALSE)),
    Second_sequence = if (value == 0) "Second treatment sequence omitted" else
      "Second treatment sequence given",
    paste(scenario_name, value)
  )
}
