# Load required packages
if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(here, dampack)

# Script 06 normally creates these objects. Recreate them cheaply when DSA is
# run on its own or from a report that intentionally skips survival resampling.
if (!exists("param_distributions") || !exists("param_groups")) {
  parameter_config <- configure_parameter_distributions(l_params_base)
  param_distributions <- parameter_config$distributions
  param_groups <- parameter_config$groups
  rm(parameter_config)
}

# The DSA varies more parameters than the PSA samples: unit prices are fixed in
# the PSA and tested deterministically here (issue #154), so the DSA reads the
# shared specification directly rather than the PSA distribution list.
param_spec <- parameter_distribution_spec()

# Load functions
source(here::here("scripts/R/functions/model_fun.R"))
source(here::here("scripts/R/functions/calculate_outcomes.R"))
source(here::here("scripts/R/functions/create_tornado_plot.R"))
source(here::here("scripts/R/functions/prediction_functions.R"))

# Use base case values as starting point
dsa_basecase <- l_params_base

# Calculate min/max ranges using DSA_mult (±20% variation), capped to the
# parameter support (costs and prevalences >= 0, utilities and prevalences <= 1).
# build_dsa_ranges() (parameter_distributions.R) is shared with table_1.qmd so
# the published ranges are the ranges actually run (issue #153). It validates
# min <= max (issue #48); min = max is allowed for zero-valued parameters.
dsa_ranges <- build_dsa_ranges(l_params_base, param_spec, mult = DSA_mult)
dsa_pars <- dsa_ranges$pars
prevalence_params <- param_groups$prevalence

# Deterministic structural scenarios (discount rate, time horizon,
# post-progression cost, second treatment sequence) run after the one-way
# parameter loop and reported in the OWSA report and Table 1.
dsa_scenarios <- dsa_structural_scenarios(
  base_dr = dr,
  base_horizon_years = time_horizon / 52,
  base_pp_cost = l_params_base$c_other_pp
)

# Print the ranges for verification
cat("Parameter ranges for DSA:\n")
print(dsa_ranges)

# Create a function to format model results consistently
format_results <- function(model_output, parameter = "base_case", value_type = "base", param_value = NA) {
  # Calculate NMB
  model_output$NMB <- model_output$Effect * WTP - model_output$Cost
  
  # Add parameter information
  model_output$Parameter <- parameter
  model_output$Value <- value_type
  model_output$ParamValue <- param_value
  
  return(model_output)
}

# Initialize an empty list to store all results
all_results <- list()
result_counter <- 1

# First run with base case
cat("Running base case...\n")
base_result <- model_fun(dsa_basecase)
all_results[[result_counter]] <- format_results(base_result)
result_counter <- result_counter + 1

# Identify optimal strategy in base case
base_nmb <- base_result$Effect * WTP - base_result$Cost
base_optimal <- base_result$Strategy[which.max(base_nmb)]
cat("Optimal strategy in base case:", base_optimal, "\n\n")

# For each parameter, run the model at min and max values
for (i in seq_len(nrow(dsa_ranges))) {
  param_name <- dsa_ranges$pars[i]
  cat("Testing parameter:", param_name, "- Min:", dsa_ranges$min[i],
      "Max:", dsa_ranges$max[i], "\n")

  for (side in c("min", "max")) {
    param_value <- dsa_ranges[[side]][i]
    params_side <- dsa_basecase
    params_side[[param_name]] <- param_value
    params_side <- sync_biomarker_test_costs(params_side)
    result_side <- model_fun(params_side)
    all_results[[result_counter]] <- format_results(
      result_side, param_name, side, param_value
    )
    result_counter <- result_counter + 1
  }
}

# ===============================================================================
# STRUCTURAL SCENARIOS (issues #153, #154)
# ===============================================================================
# Deterministic scenarios that Table 1 reports as ranges. The discount rate is
# applied to costs and QALYs together. A different time horizon re-predicts the
# base-case survival curves on the new weekly grid from the selected joint
# models and rebuilds the schedule vectors with the same rules as script 05
# (drug positions are fixed calendar weeks; CT every 12 weeks and blood tests
# every 4 weeks recur to the end of the horizon). The post-progression cost and
# second-sequence scenarios test the two scope assumptions raised in issue #154.

build_horizon_params <- function(base_params, horizon_weeks) {
  if (!exists("models") || is.null(models$best_fit$os) || is.null(models$best_fit$pfs)) {
    stop("Fitted best-fit survival models are required for a time-horizon scenario.")
  }
  tp <- seq(0, horizon_weeks)
  n <- length(tp)
  preds <- generate_population_averaged_predictions(
    models = list(os = models$best_fit$os, pfs = models$best_fit$pfs),
    strategies_df = strategies_df,
    data_complete = data_complete,
    time_points = tp,
    prevalences = setNames(strategies_df$prevalence, strategies_df$id)
  )

  p <- base_params
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

# Remove the second treatment sequence: zero every administration from the given
# schedule position onwards and rebuild the visit schedule from the remaining
# administrations. Monitoring (CT, blood tests) is unchanged, and so is
# survival, so the scenario bounds the cost of the assumption only.
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

run_structural_scenario <- function(scenario_name, value) {
  if (scenario_name == "Discount_rate") {
    p <- dsa_basecase
    p$dr_costs <- value
    p$dr_effects <- value
    return(model_fun(p))
  }
  if (scenario_name == "Time_horizon") {
    horizon_weeks <- as.integer(round(value * 52))
    p <- build_horizon_params(dsa_basecase, horizon_weeks)
    return(model_fun(p, time_horizon = horizon_weeks))
  }
  if (scenario_name == "Post_progression_cost") {
    p <- dsa_basecase
    p$c_other_pp <- value
    return(model_fun(p))
  }
  if (scenario_name == "Second_sequence") {
    p <- dsa_basecase
    if (value == 0) p <- drop_second_sequence(p)
    return(model_fun(p))
  }
  stop("Unknown structural scenario: ", scenario_name)
}

cat("\nRunning structural scenarios...\n")
dsa_scenario_status <- list()
for (scenario_name in names(dsa_scenarios)) {
  scenario <- dsa_scenarios[[scenario_name]]
  # One-sided scenarios declare only the endpoint that differs from base case.
  for (side in structural_scenario_sides(scenario)) {
    value <- scenario[[side]]
    cat("  ", scenario$label, "-", side, "=", value, scenario$unit, "\n")
    result_side <- tryCatch(
      run_structural_scenario(scenario_name, value),
      error = function(e) {
        cat("    Scenario failed:", conditionMessage(e), "\n")
        NULL
      }
    )
    dsa_scenario_status[[paste(scenario_name, side)]] <- !is.null(result_side)
    if (is.null(result_side)) next
    all_results[[result_counter]] <- format_results(
      result_side, scenario_name, side, value
    )
    result_counter <- result_counter + 1
  }
}

# Combine all results into a single data frame
cat("Combining results...\n")
dsa_results <- do.call(rbind, all_results)

# Add parameter group information. The lookup comes from the shared parameter
# specification rather than the PSA groups, so fixed unit prices are still
# grouped as drug/test costs in the tornado plots (issue #154).
group_lookup <- parameter_group_lookup(param_spec)
dsa_results$group <- sapply(dsa_results$Parameter, function(p) {
  if (p == "base_case") return("base_case")
  if (p %in% names(dsa_scenarios)) return("structural")
  if (p %in% names(group_lookup)) return(unname(group_lookup[p]))
  return("other")
})

# Calculate NMB differences from base case
cat("Calculating NMB differences...\n")
for (strat in strategies) {
  # Find base case NMB for this strategy
  base_nmb <- dsa_results$NMB[dsa_results$Parameter == "base_case" &
                                dsa_results$Strategy == strat]

  # Calculate differences
  dsa_results$NMB_diff[dsa_results$Strategy == strat] <-
    dsa_results$NMB[dsa_results$Strategy == strat] - base_nmb
}

# ===============================================================================
# SURVIVAL MODEL STRUCTURAL SENSITIVITY ANALYSIS
# ===============================================================================
# Tests all fitted parametric distributions (same distribution for OS and PFS)
# This addresses the gap in DSA where survival model choice was not varied
# See GitHub Issue #81
# ===============================================================================

cat("\nRunning survival model structural sensitivity analysis...\n")

# Determine if we have the necessary models
has_models <- FALSE
if (exists("models") && !is.null(models$full$os)) {
  has_models <- TRUE
}

if (has_models) {

  # Get list of distributions to test from the shared economic model
  distributions_tested <- names(models$full$os)
  distributions_tested <- distributions_tested[
    !vapply(models$full$os[distributions_tested], is.null, logical(1))
  ]

  cat("Distributions to test:", paste(distributions_tested, collapse = ", "), "\n")

  # Store results for each distribution
  model_sensitivity_results <- list()

  for (dist in distributions_tested) {
    cat("  Testing distribution:", dist, "\n")

    tryCatch({
      test_models <- list(
        os = models$full$os[[dist]],
        pfs = models$full$pfs[[dist]]
      )

      if (is.null(test_models$os) || is.null(test_models$pfs)) {
        cat("    Skipping - model not available for both OS and PFS\n")
        next
      }

      test_predictions <- generate_population_averaged_predictions(
        models = test_models,
        strategies_df = strategies_df,
        data_complete = data_complete,
        time_points = time_points,
        prevalences = setNames(strategies_df$prevalence, strategies_df$id)
      )

      # Build parameter list with new predictions
      test_params <- l_params_base

      # Update control survival curves
      control_strategy <- get_control_strategy()
      test_params$p_os[[paste0(control_strategy, "_OS")]] <-
        test_predictions[[control_strategy]]$os
      test_params$p_pfs[[paste0(control_strategy, "_PFS")]] <-
        test_predictions[[control_strategy]]$pfs

      # Update biomarker strategy survival curves
      biomarkers_to_loop <- get_biomarkers()
      for (biomarker in biomarkers_to_loop) {
        test_params$p_os[[paste0(biomarker, "_pos_OS")]] <- test_predictions[[biomarker]]$biomarker_positive$os
        test_params$p_os[[paste0(biomarker, "_neg_OS")]] <- test_predictions[[biomarker]]$biomarker_negative$os
        test_params$p_os[[paste0(biomarker, "_weighted_OS")]] <- test_predictions[[biomarker]]$os
        test_params$p_pfs[[paste0(biomarker, "_pos_PFS")]] <- test_predictions[[biomarker]]$biomarker_positive$pfs
        test_params$p_pfs[[paste0(biomarker, "_neg_PFS")]] <- test_predictions[[biomarker]]$biomarker_negative$pfs
        test_params$p_pfs[[paste0(biomarker, "_weighted_PFS")]] <- test_predictions[[biomarker]]$pfs
      }

      # Run model with this distribution's predictions
      dist_result <- model_fun(test_params, determpsa = "det")
      dist_result$NMB <- dist_result$Effect * WTP - dist_result$Cost
      dist_result$Distribution <- dist

      model_sensitivity_results[[dist]] <- dist_result

    }, error = function(e) {
      cat("    Error:", conditionMessage(e), "\n")
    })
  }

  # Combine results into data frame
  if (length(model_sensitivity_results) > 0) {
    model_sensitivity_df <- do.call(rbind, model_sensitivity_results)
    rownames(model_sensitivity_df) <- NULL

    # Identify best-fit distribution (used in base case)
    best_dist <- models$best_fit$os_distribution
    cat("\nBase case distribution:", best_dist, "\n")

    # Convert each distribution set to the same paired endpoint shape as DSA.
    model_endpoints <- do.call(rbind, lapply(strategies, function(strat) {
      strat_results <- model_sensitivity_df[model_sensitivity_df$Strategy == strat, ]
      if (nrow(strat_results) == 0) return(NULL)
      base_nmb_strat <- strat_results$NMB[strat_results$Distribution == best_dist]
      if (length(base_nmb_strat) == 0) {
        base_nmb_strat <- dsa_results$NMB[
          dsa_results$Parameter == "base_case" & dsa_results$Strategy == strat]
      }
      do.call(rbind, lapply(c("min", "max"), function(side) {
        endpoint <- if (side == "min") which.min(strat_results$NMB) else which.max(strat_results$NMB)
        data.frame(
          Strategy = strat, Cost = NA, Effect = NA,
          NMB = strat_results$NMB[endpoint], Parameter = "Survival_Model",
          Value = side, ParamValue = NA, group = "survival_model",
          NMB_diff = strat_results$NMB[endpoint] - base_nmb_strat,
          Distribution = strat_results$Distribution[endpoint],
          Base_NMB = base_nmb_strat, stringsAsFactors = FALSE
        )
      }))
    }))
    model_ranges <- summarise_param_ranges(model_endpoints)
    min_endpoints <- model_endpoints[model_endpoints$Value == "min", ]
    max_endpoints <- model_endpoints[model_endpoints$Value == "max", ]
    min_endpoints <- min_endpoints[match(model_ranges$Strategy, min_endpoints$Strategy), ]
    max_endpoints <- max_endpoints[match(model_ranges$Strategy, max_endpoints$Strategy), ]
    model_impact_summary <- transform(
      model_ranges,
      Base_NMB = min_endpoints$Base_NMB,
      Min_NMB = min_endpoints$NMB, Max_NMB = max_endpoints$NMB,
      Min_Dist = min_endpoints$Distribution,
      Max_Dist = max_endpoints$Distribution,
      group = "survival_model"
    )

    cat("\nSurvival model sensitivity analysis complete.\n")
    cat("NMB range by strategy:\n")
    print(model_impact_summary[, c("Strategy", "Min_Dist", "Min_NMB", "Max_Dist", "Max_NMB", "Range")])

    # Add the paired structural endpoints to the common DSA/tornado data.
    dsa_results <- rbind(dsa_results, model_endpoints[, names(dsa_results)])

  } else {
    cat("Warning: No survival model sensitivity results generated.\n")
    model_sensitivity_df <- NULL
    model_impact_summary <- NULL
  }

} else {
  cat("Warning: Required fitted models not found for the economic model.\n")
  cat("Run 04_parametric_survival_analysis.R first.\n")
  cat("Skipping survival model structural sensitivity analysis.\n")
  model_sensitivity_df <- NULL
  model_impact_summary <- NULL
}

# ===============================================================================
# END SURVIVAL MODEL STRUCTURAL SENSITIVITY ANALYSIS
# ===============================================================================

# Create tornado plots for each strategy
# First for the optimal strategy
par(mfrow = c(2, 2))  # Set up a 2x2 plot grid
optimal_tornado <- create_tornado_plot(base_optimal, dsa_results)

# Then for each biomarker strategy
biomarker_strategies <- get_biomarkers()
tornado_results <- list()

for (strat in biomarker_strategies) {
  # Skip if this is already the optimal strategy (to avoid duplication)
  if (strat == base_optimal) {
    cat("Skipping", strat, "as it is the optimal strategy and already plotted\n")
    next
  }
  
  # Create tornado plot for this strategy
  tornado_results[[strat]] <- create_tornado_plot(strat, dsa_results)
}

# Reset plot layout
par(mfrow = c(1, 1))

# Check if optimal strategy changes for any parameter values
cat("\nChecking if optimal strategy changes with parameter variations:\n")
changes_found <- FALSE

for (param_name in c(dsa_pars, names(dsa_scenarios))) {
  side_optimal <- setNames(character(2), c("min", "max"))
  for (side in names(side_optimal)) {
    side_results <- dsa_results[
      dsa_results$Parameter == param_name & dsa_results$Value == side, ]
    if (nrow(side_results) == 0) break
    side_optimal[side] <- side_results$Strategy[which.max(side_results$NMB)]
  }
  if (any(side_optimal == "")) next

  # Check if optimal strategy changes
  if (any(side_optimal != base_optimal)) {
    changes_found <- TRUE
    cat("Parameter:", param_name, "\n")
    cat("  Base optimal strategy:", base_optimal, "\n")
    cat("  Optimal at min value:", side_optimal["min"], "\n")
    cat("  Optimal at max value:", side_optimal["max"], "\n\n")
  }
}

if (!changes_found) {
  cat("The optimal strategy (", base_optimal, ") is robust to all parameter variations.\n")
}

# Create a summary table of parameter impact for all strategies
cat("\nSummary of parameter impact on NMB for all strategies:\n")
param_ranges <- summarise_param_ranges(dsa_results)
impact_summary <- do.call(rbind, lapply(strategies, function(strat) {
  top <- param_ranges[param_ranges$Strategy == strat, ]
  top <- head(top[order(-top$Range), ], 3)
  if (nrow(top) == 0) return(NULL)
  data.frame(Strategy = strat, Rank = seq_len(nrow(top)),
             Parameter = top$Parameter, Impact = top$Range)
}))

# Print impact summary
print(impact_summary)
