# ===============================================================================
# BIOMARKER-SPECIFIC SURVIVAL PLOTS  
# ===============================================================================
# This file contains plotting functions for biomarker-specific survival curves.
# Plot generation calls have been moved to para_models.Rmd
# ===============================================================================

if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(ggplot2, dplyr, scales)

# Function for biomarker-specific survival plots
plot_biomarker_survival <- function(predictions, biomarker_name, outcome_type = "OS", time_points) {
  
  # Extract predictions for this biomarker
  biomarker_pred <- predictions[[biomarker_name]]
  
  if (is.null(biomarker_pred) || is.null(biomarker_pred$biomarker_positive)) {
    return(NULL)
  }
  
  # Create plotting data
  if (outcome_type == "OS") {
    df_plot <- data.frame(
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
    df_plot <- data.frame(
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
  
  v_colors <- c("#E31A1C", "#1F78B4", "#33A02C")  # Red, Blue, Green
  names(v_colors) <- c(pos_name, neg_name, "Population Average")
  
  # Define line types
  v_linetypes <- c("solid", "solid", "dashed")
  names(v_linetypes) <- c(pos_name, neg_name, "Population Average")
  
  # Create the plot
  p <- ggplot(df_plot, aes(x = time, y = survival, color = group, linetype = group)) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(values = v_colors) +
    scale_linetype_manual(values = v_linetypes) +
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

# Function to generate all biomarker plots
generate_biomarker_plots <- function(predictions, strategies_df, time_points) {
  
  if (!exists("predictions") || length(predictions) == 0) {
    return(NULL)
  }
  
  # Create storage list for plots
  survival_plots <- list()
  
  # Get economic biomarker names from the central configuration.
  v_biomarkers <- get_biomarkers()
  
  # Generate plots for each biomarker and outcome
  for (biomarker in v_biomarkers) {
    
    # OS plot
    os_plot <- plot_biomarker_survival(
      predictions = predictions,
      biomarker_name = biomarker,
      outcome_type = "OS",
      time_points = time_points
    )
    
    if (!is.null(os_plot)) {
      survival_plots[[paste0(biomarker, "_os")]] <- os_plot
      print(os_plot)
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
      print(pfs_plot)
    }
  }
  
  return(survival_plots)
}

# Auto-execute if predictions exist
if (exists("predictions") && exists("strategies_df") && exists("time_points")) {
  survival_plots <- generate_biomarker_plots(predictions, strategies_df, time_points)
}
