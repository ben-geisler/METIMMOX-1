# Load required packages
if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(here, dampack)

# Load functions
source(here::here("R/model_fun.R"))
source(here::here("R/calculate_outcomes.R"))

# Run the model to get base case results
base_results <- model_fun(l_params_base)

# Display raw model results
print(base_results)

# Use dampack to calculate ICERs
icer_obj <- dampack::calculate_icers(
  cost = base_results$Cost,
  effect = base_results$Effect,
  strategies = base_results$Strategy
)

# Display cost-effectiveness table with ICERs
print(icer_obj)

# Use plot.icers method to create cost-effectiveness plane plot
# This works because icer_obj has class "icers"
ce_plot <- plot(
  icer_obj, 
  txtsize = 12,
  currency = "€",
  effect_units = "QALYs",
  label = "all",         # Label all strategies
  alpha = 0.8,           # Slightly transparent points
  n_x_ticks = 5,         # Number of x-axis ticks
  n_y_ticks = 6,         # Number of y-axis ticks
  xexpand = expansion(0.15),  # More padding around x-axis
  yexpand = expansion(0.15)   # More padding around y-axis
)

# Print the plot
print(ce_plot)

# Calculate Net Monetary Benefit at WTP threshold
wtp_threshold <- WTP  # Using the WTP value defined in section 1

# Calculate NMB manually
nmb_at_wtp <- data.frame(
  Strategy = base_results$Strategy,
  Cost = base_results$Cost,
  Effect = base_results$Effect,
  NMB = base_results$Effect * wtp_threshold - base_results$Cost
)
nmb_at_wtp <- nmb_at_wtp[order(-nmb_at_wtp$NMB), ]
print(nmb_at_wtp)

# Identify optimal strategy
cat("Optimal strategy at WTP threshold €", format(wtp_threshold, big.mark=","), ": ", nmb_at_wtp$Strategy[1], "\n")
optimal_strategy <- nmb_at_wtp$Strategy[1]

# If we need to see the efficiency frontier in a table format
frontier_strategies <- icer_obj$Strategy[!is.na(icer_obj$ICER) | icer_obj$Status == "Reference"]
cat("Strategies on the efficiency frontier: ", paste(frontier_strategies, collapse=", "), "\n")
#write.csv(icer_obj, "icer_obj.csv")