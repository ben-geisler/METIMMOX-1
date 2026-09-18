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
source(here::here("R/model_fun.R"))
source(here::here("R/calculate_outcomes.R"))
source(here::here("R/create_tornado_plot.R"))
source(here::here("R/prediction_functions.R"))

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
# Deterministic scenarios that Table 1 reports as ranges. The scenario code
# (discount rate applied to costs and QALYs together; time horizon re-predicted
# on the new weekly grid with the script-05 schedule rules; post-progression
# cost; second sequence omitted) lives in parameter_distributions.R as
# apply_structural_scenario(), shared with table_s8.qmd (issue #157), so the
# scenario table and the DSA run exactly the same scenarios.

horizon_models <- if (exists("models")) {
  list(os = models$best_fit$os, pfs = models$best_fit$pfs)
} else {
  NULL
}

run_structural_scenario <- function(scenario_name, value) {
  scenario <- apply_structural_scenario(
    scenario_name, value, dsa_basecase,
    base_horizon = time_horizon,
    models = horizon_models,
    strategies_df = strategies_df,
    data_complete = data_complete
  )
  model_fun(scenario$params, time_horizon = scenario$time_horizon)
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
# INCREMENTAL NMB VERSUS THE CONTROL STRATEGY (issue #156)
# ===============================================================================
# NMB_diff ranks parameters by how far they move each strategy's own NMB, so
# parameters that shift every strategy equally (utilities, end-of-life cost)
# top every tornado without affecting the decision. INMB_diff is the change in
# the strategy's incremental NMB versus the control strategy; it measures
# decision sensitivity and is the quantity Figure 2 plots.
control_strategy <- get_control_strategy()
add_incremental_nmb <- function(df) {
  key <- paste(df$Parameter, df$Value, df$ParamValue)
  control_rows <- df[df$Strategy == control_strategy, ]
  control_key <- paste(control_rows$Parameter, control_rows$Value, control_rows$ParamValue)
  df$INMB <- df$NMB - control_rows$NMB[match(key, control_key)]
  base_inmb <- df$INMB[df$Parameter == "base_case"]
  names(base_inmb) <- df$Strategy[df$Parameter == "base_case"]
  df$INMB_diff <- df$INMB - base_inmb[df$Strategy]
  df
}
dsa_results <- add_incremental_nmb(dsa_results)

# ===============================================================================
# SURVIVAL MODEL STRUCTURAL SENSITIVITY ANALYSIS
# ===============================================================================
# Tests every fitted parametric family with the same family for OS and PFS
# (issue #81). Since issue #156 a family whose OS/PFS pair violates OS >= PFS
# is not run through the model (model_fun() would stop on the ordering check)
# but is listed in dsa_distribution_status with the violation count, so the
# OWSA report can say which families were evaluated and why the others were
# not. The reference row is the ordering-constrained selected pair from
# script 04 (models$best_fit), which is the base case by construction; when
# the selected OS and PFS families differ, that pair is added explicitly
# because the same-family loop never produces it.
# ===============================================================================

cat("\nRunning survival model structural sensitivity analysis...\n")

has_models <- exists("models") && !is.null(models$full$os)

model_sensitivity_df <- NULL
model_impact_summary <- NULL
dsa_distribution_status <- NULL
dsa_results_inmb <- dsa_results

if (has_models) {

  distributions_tested <- names(models$full$os)
  distributions_tested <- distributions_tested[
    !vapply(models$full$os[distributions_tested], is.null, logical(1))
  ]
  cat("Distributions to test:", paste(distributions_tested, collapse = ", "), "\n")

  selected_pair <- list(
    os = models$best_fit$os_distribution,
    pfs = models$best_fit$pfs_distribution
  )
  selected_label <- if (identical(selected_pair$os, selected_pair$pfs)) {
    selected_pair$os
  } else {
    paste0(selected_pair$os, "/", selected_pair$pfs)
  }
  cat("Base case (selected) pair:", selected_label, "\n")

  pair_audit <- models$ordered_selection$pairs
  audit_row <- function(os_dist, pfs_dist) {
    if (is.null(pair_audit)) return(NULL)
    hit <- pair_audit[pair_audit$os_distribution == os_dist &
                        pair_audit$pfs_distribution == pfs_dist, , drop = FALSE]
    if (nrow(hit) == 0) NULL else hit[1, ]
  }

  predict_and_run <- function(test_models) {
    test_predictions <- generate_population_averaged_predictions(
      models = test_models,
      strategies_df = strategies_df,
      data_complete = data_complete,
      time_points = time_points,
      prevalences = setNames(strategies_df$prevalence, strategies_df$id),
      quiet = TRUE
    )
    test_params <- set_population_predictions(l_params_base, test_predictions)
    test_params$p_os[[paste0(control_strategy, "_OS")]] <-
      test_predictions[[control_strategy]]$os
    test_params$p_pfs[[paste0(control_strategy, "_PFS")]] <-
      test_predictions[[control_strategy]]$pfs
    for (biomarker in get_biomarkers()) {
      test_params$p_os[[paste0(biomarker, "_pos_OS")]] <- test_predictions[[biomarker]]$biomarker_positive$os
      test_params$p_os[[paste0(biomarker, "_neg_OS")]] <- test_predictions[[biomarker]]$biomarker_negative$os
      test_params$p_os[[paste0(biomarker, "_weighted_OS")]] <- test_predictions[[biomarker]]$os
      test_params$p_pfs[[paste0(biomarker, "_pos_PFS")]] <- test_predictions[[biomarker]]$biomarker_positive$pfs
      test_params$p_pfs[[paste0(biomarker, "_neg_PFS")]] <- test_predictions[[biomarker]]$biomarker_negative$pfs
      test_params$p_pfs[[paste0(biomarker, "_weighted_PFS")]] <- test_predictions[[biomarker]]$pfs
    }
    result <- model_fun(test_params, determpsa = "det")
    result$NMB <- result$Effect * WTP - result$Cost
    result
  }

  model_sensitivity_results <- list()
  status_rows <- list()

  for (dist in distributions_tested) {
    cat("  Testing distribution:", dist, "\n")
    test_models <- list(os = models$full$os[[dist]], pfs = models$full$pfs[[dist]])
    status <- data.frame(
      Distribution = dist, OS_Distribution = dist, PFS_Distribution = dist,
      Selected = identical(dist, selected_label),
      Evaluated = FALSE, N_Violations = NA_integer_, Reason = NA_character_,
      stringsAsFactors = FALSE
    )

    if (is.null(test_models$os) || is.null(test_models$pfs)) {
      status$Reason <- "model not fitted for both OS and PFS"
      cat("    Skipped:", status$Reason, "\n")
      status_rows[[dist]] <- status
      next
    }

    audit <- audit_row(dist, dist)
    if (!is.null(audit)) status$N_Violations <- audit$n_violations
    if (!is.null(audit) && !isTRUE(audit$ordered)) {
      status$Reason <- sprintf(
        "OS < PFS at %s curve-time points across the control and biomarker subgroup curves (%d weekly time points per curve); ordering constraint",
        if (is.na(audit$n_violations)) "an unknown number of" else audit$n_violations,
        length(time_points)
      )
      cat("    Skipped:", status$Reason, "\n")
      status_rows[[dist]] <- status
      next
    }

    dist_result <- tryCatch(predict_and_run(test_models), error = function(e) {
      status$Reason <<- conditionMessage(e)
      NULL
    })
    if (is.null(dist_result)) {
      cat("    Skipped:", status$Reason, "\n")
      status_rows[[dist]] <- status
      next
    }
    dist_result$Distribution <- dist
    model_sensitivity_results[[dist]] <- dist_result
    status$Evaluated <- TRUE
    status_rows[[dist]] <- status
  }

  # The selected pair is the base case; add it explicitly when it is not a
  # same-family pair so the reference row is always present.
  if (!selected_label %in% names(model_sensitivity_results)) {
    base_rows <- dsa_results[dsa_results$Parameter == "base_case",
                             c("Strategy", "Cost", "Effect", "NMB")]
    base_rows$Distribution <- selected_label
    model_sensitivity_results[[selected_label]] <- base_rows
    audit <- audit_row(selected_pair$os, selected_pair$pfs)
    status_rows[[selected_label]] <- data.frame(
      Distribution = selected_label, OS_Distribution = selected_pair$os,
      PFS_Distribution = selected_pair$pfs, Selected = TRUE, Evaluated = TRUE,
      N_Violations = if (is.null(audit)) NA_integer_ else audit$n_violations,
      Reason = "selected ordering-constrained pair (base case)",
      stringsAsFactors = FALSE
    )
  }

  dsa_distribution_status <- do.call(rbind, status_rows)
  rownames(dsa_distribution_status) <- NULL
  dsa_distribution_status <- dsa_distribution_status[
    order(!dsa_distribution_status$Selected, !dsa_distribution_status$Evaluated,
          dsa_distribution_status$Distribution), ]
  n_skipped <- sum(!dsa_distribution_status$Evaluated)
  cat(sprintf("\nDistributions evaluated: %d of %d; skipped: %d\n",
              sum(dsa_distribution_status$Evaluated), nrow(dsa_distribution_status), n_skipped))
  if (n_skipped > 0) {
    skipped <- dsa_distribution_status[!dsa_distribution_status$Evaluated, ]
    for (r in seq_len(nrow(skipped))) {
      cat("  -", skipped$Distribution[r], ":", skipped$Reason[r], "\n")
    }
  }

  if (length(model_sensitivity_results) > 0) {
    model_sensitivity_df <- do.call(rbind, model_sensitivity_results)
    rownames(model_sensitivity_df) <- NULL
    model_sensitivity_df$NMB <- model_sensitivity_df$Effect * WTP - model_sensitivity_df$Cost
    model_sensitivity_df$Selected <- model_sensitivity_df$Distribution == selected_label

    # Incremental NMB of each strategy versus control under the same family.
    ctrl <- model_sensitivity_df[model_sensitivity_df$Strategy == control_strategy, ]
    model_sensitivity_df$INMB <- model_sensitivity_df$NMB -
      ctrl$NMB[match(model_sensitivity_df$Distribution, ctrl$Distribution)]

    # Base-case anchors: the selected pair, taken from the base-case rows so the
    # anchor is the pair actually used, not the OS family alone.
    base_nmb_by_strategy <- setNames(
      dsa_results$NMB[dsa_results$Parameter == "base_case"],
      dsa_results$Strategy[dsa_results$Parameter == "base_case"]
    )
    base_inmb_by_strategy <- setNames(
      dsa_results$INMB[dsa_results$Parameter == "base_case"],
      dsa_results$Strategy[dsa_results$Parameter == "base_case"]
    )
    model_sensitivity_df$NMB_diff <- model_sensitivity_df$NMB -
      base_nmb_by_strategy[model_sensitivity_df$Strategy]
    model_sensitivity_df$INMB_diff <- model_sensitivity_df$INMB -
      base_inmb_by_strategy[model_sensitivity_df$Strategy]

    # Paired endpoints (lowest / highest value of the measure across the
    # evaluated families) in the same shape as the one-way DSA rows.
    survival_model_endpoints <- function(measure) {
      do.call(rbind, lapply(strategies, function(strat) {
        strat_results <- model_sensitivity_df[model_sensitivity_df$Strategy == strat, ]
        if (nrow(strat_results) == 0) return(NULL)
        do.call(rbind, lapply(c("min", "max"), function(side) {
          endpoint <- if (side == "min") which.min(strat_results[[measure]]) else
            which.max(strat_results[[measure]])
          data.frame(
            Strategy = strat, Cost = NA, Effect = NA,
            NMB = strat_results$NMB[endpoint], Parameter = "Survival_Model",
            Value = side, ParamValue = NA, group = "survival_model",
            NMB_diff = strat_results$NMB_diff[endpoint],
            INMB = strat_results$INMB[endpoint],
            INMB_diff = strat_results$INMB_diff[endpoint],
            Distribution = strat_results$Distribution[endpoint],
            Base_NMB = unname(base_nmb_by_strategy[strat]),
            stringsAsFactors = FALSE
          )
        }))
      }))
    }
    model_endpoints <- survival_model_endpoints("NMB")
    model_endpoints_inmb <- survival_model_endpoints("INMB")

    model_ranges <- summarise_param_ranges(model_endpoints)
    min_endpoints <- model_endpoints[model_endpoints$Value == "min", ]
    max_endpoints <- model_endpoints[model_endpoints$Value == "max", ]
    min_endpoints <- min_endpoints[match(model_ranges$Strategy, min_endpoints$Strategy), ]
    max_endpoints <- max_endpoints[match(model_ranges$Strategy, max_endpoints$Strategy), ]
    model_impact_summary <- transform(
      model_ranges,
      Base_NMB = min_endpoints$Base_NMB,
      Base_Dist = selected_label,
      Min_NMB = min_endpoints$NMB, Max_NMB = max_endpoints$NMB,
      Min_Dist = min_endpoints$Distribution,
      Max_Dist = max_endpoints$Distribution,
      N_Evaluated = sum(dsa_distribution_status$Evaluated),
      N_Candidates = nrow(dsa_distribution_status),
      group = "survival_model"
    )

    cat("\nSurvival model sensitivity analysis complete.\n")
    cat("NMB range by strategy:\n")
    print(model_impact_summary[, c("Strategy", "Base_Dist", "Min_Dist", "Min_NMB",
                                   "Max_Dist", "Max_NMB", "Range")])

    # Add the paired structural endpoints to the common DSA/tornado data. The
    # NMB tornado data take the NMB-chosen endpoints; the incremental tornado
    # data take the endpoints chosen on incremental NMB.
    dsa_results <- rbind(dsa_results, model_endpoints[, names(dsa_results)])
    dsa_results_inmb <- rbind(dsa_results_inmb,
                              model_endpoints_inmb[, names(dsa_results_inmb)])
  } else {
    cat("Warning: No survival model sensitivity results generated.\n")
  }

} else {
  cat("Warning: Required fitted models not found for the economic model.\n")
  cat("Run 04_parametric_survival_analysis.R first.\n")
  cat("Skipping survival model structural sensitivity analysis.\n")
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

# Incremental-NMB tornado plots for the guided strategies (issue #156): rank by
# decision sensitivity, i.e. the change in NMB versus the control strategy.
incremental_tornado_results <- list()
for (strat in biomarker_strategies) {
  incremental_tornado_results[[strat]] <-
    create_tornado_plot(strat, dsa_results_inmb, measure = "INMB_diff")
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
print(impact_summary)

# Same ranking on incremental NMB versus control (issue #156), guided strategies only.
cat("\nSummary of parameter impact on incremental NMB versus control:\n")
param_ranges_inmb <- summarise_param_ranges(dsa_results_inmb, measure = "INMB_diff")
incremental_impact_summary <- do.call(rbind, lapply(biomarker_strategies, function(strat) {
  top <- param_ranges_inmb[param_ranges_inmb$Strategy == strat, ]
  top <- head(top[order(-top$Range), ], 3)
  if (nrow(top) == 0) return(NULL)
  data.frame(Strategy = strat, Rank = seq_len(nrow(top)),
             Parameter = top$Parameter, Impact = top$Range)
}))
print(incremental_impact_summary)
