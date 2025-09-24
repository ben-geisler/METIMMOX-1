# ===============================================================================
# BIOMARKER-SPECIFIC SURVIVAL PLOTS  
# ===============================================================================
# This file creates separate OS and PFS plots for each biomarker showing
# biomarker-positive vs biomarker-negative survival curves
# ===============================================================================

if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(ggplot2, dplyr, scales)

# Modified plotting function for biomarker-specific plots
plot_biomarker_survival <- function(predictions, biomarker_name, outcome_type = "OS", time_points) {
  
  # Extract predictions for this biomarker
  biomarker_pred <- predictions[[biomarker_name]]
  
  if (is.null(biomarker_pred) || is.null(biomarker_pred$biomarker_positive)) {
    cat("No predictions available for", biomarker_name, "\n")
    return(NULL)
  }
  
  # Create plotting data
  if (outcome_type == "OS") {
    plot_data <- data.frame(
      time = rep(time_points, 3),
      survival = c(
        biomarker_pred$biomarker_positive$os,    # Biomarker+ (experimental treatment)
        biomarker_pred$biomarker_negative$os,    # Biomarker- (control treatment) 
        biomarker_pred$os                        # Population average
      ),
      group = rep(c(
        paste0(toupper(biomarker_name), "+"),
        paste0(toupper(biomarker_name), "-"), 
        "Population Average"
      ), each = length(time_points))
    )
    y_label <- "Overall Survival Probability"
    title_outcome <- "Overall Survival"
  } else {
    plot_data <- data.frame(
      time = rep(time_points, 3),
      survival = c(
        biomarker_pred$biomarker_positive$pfs,   # Biomarker+ (experimental treatment)
        biomarker_pred$biomarker_negative$pfs,   # Biomarker- (control treatment)
        biomarker_pred$pfs                       # Population average
      ),
      group = rep(c(
        paste0(toupper(biomarker_name), "+"),
        paste0(toupper(biomarker_name), "-"),
        "Population Average"
      ), each = length(time_points))
    )
    y_label <- "Progression-Free Survival Probability"  
    title_outcome <- "Progression-Free Survival"
  }
  
  # Define colors for consistency - using proper named vector syntax
  pos_name <- paste0(toupper(biomarker_name), "+")
  neg_name <- paste0(toupper(biomarker_name), "-")
  
  colors <- c("#E31A1C", "#1F78B4", "#33A02C")  # Red, Blue, Green
  names(colors) <- c(pos_name, neg_name, "Population Average")
  
  # Define line types
  linetypes <- c("solid", "solid", "dashed")
  names(linetypes) <- c(pos_name, neg_name, "Population Average")
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = time, y = survival, color = group, linetype = group)) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(values = colors) +
    scale_linetype_manual(values = linetypes) +
    theme_minimal() +
    labs(
      title = paste(title_outcome, "by", toupper(biomarker_name), "Status"),
      subtitle = paste("Prevalence:", round(biomarker_pred$prevalence * 100, 1), "%"),
      x = "Time (weeks)",
      y = y_label,
      color = "Group",
      linetype = "Group"
    ) +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    scale_x_continuous(breaks = seq(0, max(time_points), by = 52),
                       labels = seq(0, max(time_points), by = 52) / 52)
  
  return(p)
}

# Generate plots if predictions are available
if (exists("predictions") && length(predictions) > 0) {
  
  cat("\n=== GENERATING BIOMARKER-SPECIFIC SURVIVAL PLOTS ===\n")
  
  # Create storage list for plots
  survival_plots <- list()
  
  # Get biomarker names (excluding control)
  biomarkers <- strategies_df$id[strategies_df$id != "control"]
  
  # Generate plots for each biomarker and outcome
  for (biomarker in biomarkers) {
    cat("Creating plots for", toupper(biomarker), "...\n")
    
    # OS plot
    os_plot <- plot_biomarker_survival(
      predictions = predictions,
      biomarker_name = biomarker,
      outcome_type = "OS",
      time_points = time_points
    )
    
    if (!is.null(os_plot)) {
      survival_plots[[paste0(biomarker, "_os")]] <- os_plot
      cat("  - OS plot created\n")
    }
    
    # PFS plot
    pfs_plot <- plot_biomarker_survival(
      predictions = predictions, 
      biomarker_name = biomarker,
      outcome_type = "PFS",
      time_points = time_points
    )
    
    if (!is.null(pfs_plot)) {
      survival_plots[[paste0(biomarker, "_pfs")]] <- pfs_plot
      cat("  - PFS plot created\n")
    }
  }
  
  # Display plots
  cat("\nDisplaying survival plots:\n")
  for (plot_name in names(survival_plots)) {
    cat("Showing:", plot_name, "\n")
    print(survival_plots[[plot_name]])
    cat("\n")
  }
  
  # Save plots if directory exists or can be created
  plot_dir <- "survival_plots"
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }
  
  if (dir.exists(plot_dir)) {
    for (plot_name in names(survival_plots)) {
      filename <- paste0(plot_dir, "/", plot_name, "_biomarker_survival.png")
      ggsave(
        filename = filename,
        plot = survival_plots[[plot_name]],
        width = 10, height = 6, dpi = 300
      )
    }
    cat("All plots saved to", plot_dir, "directory\n")
  }
  
  cat("Total plots generated:", length(survival_plots), "\n")
  cat("=== BIOMARKER PLOTS COMPLETE ===\n\n")
  
} else {
  cat("Survival plots not available. Run main analysis first.\n")
}