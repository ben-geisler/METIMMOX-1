#' Summarize paired DSA endpoints into parameter ranges
#'
#' @param results_df DSA results containing min/max rows and NMB differences
#' @param strategy_name Optional strategy to select
#' @param measure Column holding the change from base case: "NMB_diff" (per-
#'   strategy net monetary benefit) or "INMB_diff" (incremental NMB versus the
#'   control strategy; issue #156)
#' @return Data frame with one row per strategy and parameter
summarise_param_ranges <- function(results_df, strategy_name = NULL,
                                   measure = c("NMB_diff", "INMB_diff")) {
  measure <- match.arg(measure)
  if (!measure %in% names(results_df)) {
    stop("results_df has no '", measure, "' column")
  }
  x <- results_df[results_df$Parameter != "base_case", ]
  if (!is.null(strategy_name)) x <- x[x$Strategy == strategy_name, ]
  groups <- split(x, interaction(x$Strategy, x$Parameter, drop = TRUE))
  rows <- lapply(groups, function(group) {
    endpoints <- match(c("min", "max"), group$Value)
    if (all(is.na(endpoints))) return(NULL)
    # A one-sided structural scenario runs only the endpoint that differs from
    # the base case (issue #154). The endpoint it omits IS the base case, whose
    # difference is zero by definition, so the bar still spans base to
    # scenario rather than dropping out of the tornado altogether.
    diffs <- ifelse(is.na(endpoints), 0, group[[measure]][endpoints])
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

#' Tornado plot for one strategy
#'
#' @param strategy_name Strategy identifier
#' @param results_df DSA results (see summarise_param_ranges)
#' @param measure "NMB_diff" ranks parameters by how far they move the
#'   strategy's own NMB, which places parameters that move every strategy
#'   equally (utilities, end-of-life cost) at the top. "INMB_diff" ranks by the
#'   change in incremental NMB versus the control strategy, i.e. by decision
#'   sensitivity (issue #156).
#' @param title Optional plot title
#' @return The ranked tornado data (invisibly NULL when nothing to plot)
create_tornado_plot <- function(strategy_name, results_df,
                                measure = c("NMB_diff", "INMB_diff"),
                                title = NULL) {
  measure <- match.arg(measure)
  cat("Creating tornado plot for strategy:", strategy_name, "(", measure, ")\n")
  tornado_data <- summarise_param_ranges(results_df, strategy_name, measure)
  tornado_data <- tornado_data[, c("Parameter", "Min_diff", "Max_diff", "Range")]

  # Sort by range for tornado plot (descending order)
  tornado_data <- tornado_data[order(-tornado_data$Range), ]

  if (nrow(tornado_data) == 0) {
    cat("No data available for tornado plot for", strategy_name, "\n")
    return(NULL)
  }

  measure_label <- if (measure == "NMB_diff") {
    "Change in Net Monetary Benefit (€)"
  } else {
    "Change in incremental NMB versus standard of care (€)"
  }
  if (is.null(title)) {
    title <- if (measure == "NMB_diff") {
      paste("Tornado Plot for", strategy_name)
    } else {
      paste("Incremental NMB Tornado Plot for", strategy_name, "vs control")
    }
  }

  # Create matrix for barplot - REVERSED order to get highest impact at the top
  barplot_data <- t(as.matrix(tornado_data[nrow(tornado_data):1, c("Min_diff", "Max_diff")]))
  colnames(barplot_data) <- tornado_data$Parameter[nrow(tornado_data):1]

  barplot(
    barplot_data,
    beside = TRUE,
    horiz = TRUE,
    las = 1,  # Horizontal labels
    col = c("blue", "red"),
    main = title,
    xlab = measure_label
  )
  abline(v = 0, lty = 2)  # Add line at zero
  legend("bottomright",
         legend = c("Min Value", "Max Value"),
         fill = c("blue", "red"),
         cex = 0.8)

  tornado_data
}
