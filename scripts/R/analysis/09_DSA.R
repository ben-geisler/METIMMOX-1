# Load required packages
if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(here, dampack)

# Load functions
source(here::here("scripts/R/functions/model_fun.R"))
source(here::here("scripts/R/functions/calculate_outcomes.R"))
source(here::here("scripts/R/functions/create_tornado_plot.R"))
source(here::here("scripts/R/functions/prediction_functions.R"))

# Ensure consistent time indexing
if (!exists("time_points_length")) {
  time_points_length <- length(time_points)
}
# Validate time_points consistency
if (length(time_points) != time_points_length) {
  stop("time_points length inconsistency detected")
}

# Define parameters to vary in sensitivity analysis (matches EVPPI parameters)
dsa_pars <- c("c_drug_nivo", "c_drug_FLOX", "c_test_CRP", "c_test_NGS", "c_test_CT",
              "c_test_blood", "c_other_visit", "c_other_baseline",
              "c_other_follow", "c_other_last", "u_np", "u_p",
              "p_crp", "p_tmb_braf")

# Use base case values as starting point
dsa_basecase <- l_params_base

# Calculate min/max ranges using DSA_mult (±20% variation)
# But cap utility values at 1.0
dsa_ranges <- data.frame(
  pars = dsa_pars,
  min = unlist(l_params_base[dsa_pars]) * (1 - DSA_mult),
  max = unlist(l_params_base[dsa_pars]) * (1 + DSA_mult),
  stringsAsFactors = FALSE
)

# Cap utility values at 1.0
utility_params <- c("u_np", "u_p")
for (param in utility_params) {
  idx <- which(dsa_ranges$pars == param)
  if (length(idx) > 0) {
    dsa_ranges$max[idx] <- min(dsa_ranges$max[idx], 1.0)
  }
}

# Cap prevalence values between 0 and 1
prevalence_params <- c("p_crp", "p_tmb_braf")
for (param in prevalence_params) {
  idx <- which(dsa_ranges$pars == param)
  if (length(idx) > 0) {
    dsa_ranges$min[idx] <- max(dsa_ranges$min[idx], 0.0)
    dsa_ranges$max[idx] <- min(dsa_ranges$max[idx], 1.0)
  }
}

# Ensure cost parameters are non-negative (Issue #48)
cost_params <- c("c_drug_nivo", "c_drug_FLOX", "c_test_CRP", "c_test_NGS", "c_test_CT",
                 "c_test_blood", "c_other_visit", "c_other_baseline",
                 "c_other_follow", "c_other_last")
for (param in cost_params) {
  idx <- which(dsa_ranges$pars == param)
  if (length(idx) > 0) {
    dsa_ranges$min[idx] <- max(dsa_ranges$min[idx], 0)
  }
}

# Validate all ranges have min < max (Issue #48)
# Note: min = max is allowed for zero-valued parameters (e.g., c_test_blood = 0)
if (any(dsa_ranges$min > dsa_ranges$max)) {
  invalid <- dsa_ranges$pars[dsa_ranges$min > dsa_ranges$max]
  stop("Invalid DSA ranges (min > max) for: ", paste(invalid, collapse = ", "))
}

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
  param_min <- dsa_ranges$min[i]
  param_max <- dsa_ranges$max[i]
  
  cat("Testing parameter:", param_name, "- Min:", param_min, "Max:", param_max, "\n")
  
  # Create parameter sets for min and max values
  params_min <- dsa_basecase
  params_min[[param_name]] <- param_min
  if (param_name == "c_test_CRP") params_min$c_test_biomarker$crp <- param_min
  if (param_name == "c_test_NGS") params_min$c_test_biomarker$tmb_braf <- param_min
  
  params_max <- dsa_basecase
  params_max[[param_name]] <- param_max
  if (param_name == "c_test_CRP") params_max$c_test_biomarker$crp <- param_max
  if (param_name == "c_test_NGS") params_max$c_test_biomarker$tmb_braf <- param_max
  
  # Run model with min value
  result_min <- model_fun(params_min)
  all_results[[result_counter]] <- format_results(result_min, param_name, "min", param_min)
  result_counter <- result_counter + 1
  
  # Run model with max value
  result_max <- model_fun(params_max)
  all_results[[result_counter]] <- format_results(result_max, param_name, "max", param_max)
  result_counter <- result_counter + 1
}

# Combine all results into a single data frame
cat("Combining results...\n")
dsa_results <- do.call(rbind, all_results)

# Add parameter group information (param_groups from 06_sampling.R)
dsa_results$group <- sapply(dsa_results$Parameter, function(p) {
  if (p == "base_case") return("base_case")
  for (g in names(param_groups)) {
    if (p %in% param_groups[[g]] && g != "all_costs") return(g)
  }
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
      test_params$p_os$control_OS <- test_predictions$control$os
      test_params$p_pfs$control_PFS <- test_predictions$control$pfs

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

    # For each strategy, find min and max NMB across distributions
    model_impact_summary <- data.frame()
    for (strat in strategies) {
      strat_results <- model_sensitivity_df[model_sensitivity_df$Strategy == strat, ]

      if (nrow(strat_results) > 0) {
        base_nmb_strat <- strat_results$NMB[strat_results$Distribution == best_dist]

        # Handle case where base distribution might not be in results
        if (length(base_nmb_strat) == 0) {
          base_nmb_strat <- dsa_results$NMB[dsa_results$Parameter == "base_case" &
                                              dsa_results$Strategy == strat]
        }

        min_nmb <- min(strat_results$NMB)
        max_nmb <- max(strat_results$NMB)
        min_dist <- strat_results$Distribution[which.min(strat_results$NMB)]
        max_dist <- strat_results$Distribution[which.max(strat_results$NMB)]

        model_impact_summary <- rbind(model_impact_summary, data.frame(
          Strategy = strat,
          Parameter = "Survival_Model",
          Base_NMB = base_nmb_strat,
          Min_NMB = min_nmb,
          Max_NMB = max_nmb,
          Min_diff = min_nmb - base_nmb_strat,
          Max_diff = max_nmb - base_nmb_strat,
          Range = max_nmb - min_nmb,
          Min_Dist = min_dist,
          Max_Dist = max_dist,
          group = "survival_model",
          stringsAsFactors = FALSE
        ))
      }
    }

    cat("\nSurvival model sensitivity analysis complete.\n")
    cat("NMB range by strategy:\n")
    print(model_impact_summary[, c("Strategy", "Min_Dist", "Min_NMB", "Max_Dist", "Max_NMB", "Range")])

    # Add survival model results to dsa_results for tornado diagram
    for (strat in strategies) {
      strat_summary <- model_impact_summary[model_impact_summary$Strategy == strat, ]

      if (nrow(strat_summary) > 0) {
        # Add "min" entry
        dsa_results <- rbind(dsa_results, data.frame(
          Strategy = strat,
          Cost = NA,  # Not used for tornado
          Effect = NA,
          NMB = strat_summary$Min_NMB,
          Parameter = "Survival_Model",
          Value = "min",
          ParamValue = NA,
          group = "survival_model",
          NMB_diff = strat_summary$Min_diff,
          stringsAsFactors = FALSE
        ))

        # Add "max" entry
        dsa_results <- rbind(dsa_results, data.frame(
          Strategy = strat,
          Cost = NA,
          Effect = NA,
          NMB = strat_summary$Max_NMB,
          Parameter = "Survival_Model",
          Value = "max",
          ParamValue = NA,
          group = "survival_model",
          NMB_diff = strat_summary$Max_diff,
          stringsAsFactors = FALSE
        ))
      }
    }

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

for (param_name in dsa_pars) {
  # Get results for min value
  min_results <- dsa_results[dsa_results$Parameter == param_name & 
                               dsa_results$Value == "min", ]
  if (nrow(min_results) > 0) {
    min_optimal <- min_results$Strategy[which.max(min_results$NMB)]
  } else {
    next
  }
  
  # Get results for max value
  max_results <- dsa_results[dsa_results$Parameter == param_name & 
                               dsa_results$Value == "max", ]
  if (nrow(max_results) > 0) {
    max_optimal <- max_results$Strategy[which.max(max_results$NMB)]
  } else {
    next
  }
  
  # Check if optimal strategy changes
  if (min_optimal != base_optimal || max_optimal != base_optimal) {
    changes_found <- TRUE
    cat("Parameter:", param_name, "\n")
    cat("  Base optimal strategy:", base_optimal, "\n")
    cat("  Optimal at min value:", min_optimal, "\n")
    cat("  Optimal at max value:", max_optimal, "\n\n")
  }
}

if (!changes_found) {
  cat("The optimal strategy (", base_optimal, ") is robust to all parameter variations.\n")
}

# Create a summary table of parameter impact for all strategies
cat("\nSummary of parameter impact on NMB for all strategies:\n")
impact_summary <- data.frame()

# For each strategy, get the top 3 most impactful parameters
for (strat in strategies) {
  # Extract results for this strategy
  strat_results <- dsa_results[dsa_results$Strategy == strat & 
                                 dsa_results$Parameter != "base_case", ]
  
  # Calculate impact for each parameter
  param_impact <- data.frame()
  for (param in unique(strat_results$Parameter)) {
    min_row <- strat_results[strat_results$Parameter == param & 
                               strat_results$Value == "min", ]
    max_row <- strat_results[strat_results$Parameter == param & 
                               strat_results$Value == "max", ]
    
    # Skip if we don't have both min and max
    if (nrow(min_row) == 0 || nrow(max_row) == 0) next
    
    min_diff <- min_row$NMB_diff
    max_diff <- max_row$NMB_diff
    range <- abs(max_diff - min_diff)
    
    param_impact <- rbind(param_impact, data.frame(
      Parameter = param,
      Range = range
    ))
  }
  
  # Sort by impact
  if (nrow(param_impact) > 0) {
    param_impact <- param_impact[order(-param_impact$Range), ]
    
    # Get top 3 parameters (or fewer if there are fewer parameters)
    top_n <- min(3, nrow(param_impact))
    top_params <- param_impact[1:top_n, ]
    
    # Add to summary
    for (i in seq_len(top_n)) {
      impact_summary <- rbind(impact_summary, data.frame(
        Strategy = strat,
        Rank = i,
        Parameter = top_params$Parameter[i],
        Impact = top_params$Range[i]
      ))
    }
  }
}

# Print impact summary
print(impact_summary)
