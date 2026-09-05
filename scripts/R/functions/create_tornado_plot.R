#' Summarize paired DSA endpoints into parameter ranges
#'
#' @param results_df DSA results containing min/max rows and NMB differences
#' @param strategy_name Optional strategy to select
#' @return Data frame with one row per strategy and parameter
summarise_param_ranges <- function(results_df, strategy_name = NULL) {
  x <- results_df[results_df$Parameter != "base_case", ]
  if (!is.null(strategy_name)) x <- x[x$Strategy == strategy_name, ]
  groups <- split(x, interaction(x$Strategy, x$Parameter, drop = TRUE))
  rows <- lapply(groups, function(group) {
    endpoints <- match(c("min", "max"), group$Value)
    if (all(is.na(endpoints))) return(NULL)
    # A one-sided structural scenario runs only the endpoint that differs from
    # the base case (issue #154). The endpoint it omits IS the base case, whose
    # NMB difference is zero by definition, so the bar still spans base to
    # scenario rather than dropping out of the tornado altogether.
    diffs <- ifelse(is.na(endpoints), 0, group$NMB_diff[endpoints])
    data.frame(
      Strategy = group$Strategy[1], Parameter = group$Parameter[1],
      Min_diff = diffs[1], Max_diff = diffs[2],
      Range = abs(diff(diffs)), stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) {
    return(data.frame(Strategy = character(), Parameter = character(),
                      Min_diff = numeric(), Max_diff = numeric(), Range = numeric()))
  }
  do.call(rbind, rows)
}

# Function to create tornado plot for a specific strategy
create_tornado_plot <- function(strategy_name, results_df) {
  cat("Creating tornado plot for strategy:", strategy_name, "\n")
  tornado_data <- summarise_param_ranges(results_df, strategy_name)
  tornado_data <- tornado_data[, c("Parameter", "Min_diff", "Max_diff", "Range")]
  
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
