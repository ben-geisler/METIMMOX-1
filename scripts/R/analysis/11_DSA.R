# Load required packages
if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(here, dampack)

# Load functions
source(here::here("scripts/R/functions/model_fun.R"))
source(here::here("scripts/R/functions/calculate_outcomes.R"))
source(here::here("scripts/R/functions/create_tornado_plot.R"))

# Ensure consistent time indexing
if (!exists("time_points_length")) {
  time_points_length <- length(time_points)
}
# Validate time_points consistency
if (length(time_points) != time_points_length) {
  stop("time_points length inconsistency detected")
}

# Define parameters to vary in sensitivity analysis
dsa_pars <- c("c_drug_nivo", "c_drug_FLOX", "c_test_NGS",
              "c_test_CT",  "u_np", "u_p", "c_other_last")

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
for (i in 1:nrow(dsa_ranges)) {
  param_name <- dsa_ranges$pars[i]
  param_min <- dsa_ranges$min[i]
  param_max <- dsa_ranges$max[i]
  
  cat("Testing parameter:", param_name, "- Min:", param_min, "Max:", param_max, "\n")
  
  # Create parameter sets for min and max values
  params_min <- dsa_basecase
  params_min[[param_name]] <- param_min
  
  params_max <- dsa_basecase
  params_max[[param_name]] <- param_max
  
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

# Create tornado plots for each strategy
# First for the optimal strategy
par(mfrow = c(2, 2))  # Set up a 2x2 plot grid
optimal_tornado <- create_tornado_plot(base_optimal, dsa_results)

# Then for each biomarker strategy
biomarker_strategies <- c("crp", "tlr", "tmb_braf")
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
    for (i in 1:top_n) {
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