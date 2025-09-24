# Function to create tornado plot for a specific strategy
create_tornado_plot <- function(strategy_name, results_df) {
  cat("Creating tornado plot for strategy:", strategy_name, "\n")
  strat_results <- results_df[results_df$Strategy == strategy_name & 
                                results_df$Parameter != "base_case", ]
  
  # Reshape data for easier plotting
  tornado_data <- data.frame()
  for (param in unique(strat_results$Parameter)) {
    min_row <- strat_results[strat_results$Parameter == param & 
                               strat_results$Value == "min", ]
    max_row <- strat_results[strat_results$Parameter == param & 
                               strat_results$Value == "max", ]
    
    # Skip if we don't have both min and max
    if (nrow(min_row) == 0 || nrow(max_row) == 0) next
    
    min_diff <- min_row$NMB_diff
    max_diff <- max_row$NMB_diff
    
    tornado_data <- rbind(tornado_data, data.frame(
      Parameter = param,
      Min_diff = min_diff,
      Max_diff = max_diff,
      Range = abs(max_diff - min_diff)
    ))
  }
  
  # Sort by range for tornado plot (descending order)
  tornado_data <- tornado_data[order(-tornado_data$Range), ]
  
  # Make sure we have data to plot
  if (nrow(tornado_data) > 0) {
    # Create matrix for barplot - REVERSED order to get highest impact at the top
    barplot_data <- t(as.matrix(tornado_data[nrow(tornado_data):1, c("Min_diff", "Max_diff")]))
    colnames(barplot_data) <- tornado_data$Parameter[nrow(tornado_data):1]
    
    # Create the barplot
    barplot(
      barplot_data,
      beside = TRUE,
      horiz = TRUE,
      las = 1,  # Horizontal labels
      col = c("blue", "red"),
      main = paste("Tornado Plot for", strategy_name),
      xlab = "Change in Net Monetary Benefit (€)"
    )
    abline(v = 0, lty = 2)  # Add line at zero
    legend("bottomright", 
           legend = c("Min Value", "Max Value"), 
           fill = c("blue", "red"),
           cex = 0.8)
    
    # Return the tornado data for further analysis
    return(tornado_data)
  } else {
    cat("No data available for tornado plot for", strategy_name, "\n")
    return(NULL)
  }
}