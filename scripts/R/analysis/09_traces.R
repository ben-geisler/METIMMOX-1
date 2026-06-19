# Ensure consistent time indexing
if (!exists("time_points_length")) {
  time_points_length <- length(time_points)
}
# Validate time_points consistency
if (length(time_points) != time_points_length) {
  stop("time_points length inconsistency detected")
}

# Generate traces by running the model with return_traces = TRUE
cat("Generating state occupancy traces for all strategies...\n")

# Generate traces by running the model once with return_traces = TRUE.
# model_fun() returns list(results, traces); `traces` is keyed by strategy
# (control, crp, tmb_braf) with weighted p_pf/p_p/p_d, and each biomarker entry
# additionally carries $positive/$negative subgroup traces.
model_output <- model_fun(l_params_base, time_horizon = time_horizon, cl = cl,
                          determpsa = "det", return_traces = TRUE, sim_idx = NULL)
traces <- model_output$traces

cat("Traces generated successfully.\n")

# Create faceted plot for all four strategies
cat("Creating faceted plot for all strategies...\n")

# Prepare data for all strategies
all_strategies_df <- data.frame()

for(strategy in strategies) {
  # Get the relevant traces
  p_pf <- traces[[strategy]]$p_pf
  p_p <- traces[[strategy]]$p_p
  p_d <- traces[[strategy]]$p_d
  
  # Extract data for plotting
  max_cycles <- min(10 * 52 + 1, time_points_length)
  cycles <- 0:(time_points_length-1)
  
  # Create data frame
  df <- data.frame(
    Year = cycles[1:max_cycles] / 52,
    PF = p_pf[1:max_cycles],
    Progressed = p_p[1:max_cycles],
    Dead = p_d[1:max_cycles],
    Strategy = strategy
  )
  
  # Convert to long format
  df_long <- reshape2::melt(
    df,
    id.vars = c("Year", "Strategy"),
    measure.vars = c("PF", "Progressed", "Dead"),
    variable.name = "State",
    value.name = "Proportion"
  )
  
  # Append to all strategies dataframe
  all_strategies_df <- rbind(all_strategies_df, df_long)
}

# Create strategy labels for faceting
# Labels are dynamic based on number of strategies (Model B excludes TLR)
strategy_labels <- if (length(strategies) == 3) {
  c("Standard of Care", "CRP Strategy", "TMB/BRAF Strategy")
} else {
  c("Standard of Care", "CRP Strategy", "TLR Strategy", "TMB/BRAF Strategy")
}
all_strategies_df$Strategy <- factor(all_strategies_df$Strategy,
                                     levels = strategies,
                                     labels = strategy_labels)

# Create faceted plot
all_strategies_plot <- ggplot(all_strategies_df, aes(x = Year, y = Proportion, fill = State)) +
  geom_area() +
  facet_wrap(~ Strategy, ncol = 2) +
  scale_fill_manual(values = c("PF" = "green", "Progressed" = "orange", "Dead" = "red")) +
  labs(title = "State Occupancy by Strategy",
       x = "Years",
       y = "Proportion of Cohort") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 9, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 11),
    strip.text = element_text(size = 10),
    strip.background = element_rect(fill = "lightgray"),
    panel.grid.major = element_line(color = "gray80"),
    panel.grid.minor = element_line(color = "gray90"),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  ) +
  scale_x_continuous(
    breaks = seq(0, 10, by = 2),
    minor_breaks = seq(0, 10, by = 1),
    limits = c(0, 10),
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    minor_breaks = seq(0, 1, by = 0.1),
    expand = c(0.01, 0)
  )

print(all_strategies_plot)

# Save the all strategies plot
cat("Saving all strategies plot as PNG file...\n")
ggsave("all_strategies_state_occupancy.png", plot = all_strategies_plot, 
       width = 10, height = 8, dpi = 300)

# Create faceted plots for biomarker positive and negative status for each biomarker
for(biomarker in biomarkers) {
  cat(paste0("Creating faceted plot for ", biomarker, " biomarker subgroups...\n"))
  
  # Prepare data for positive and negative subgroups
  pos_neg_df <- data.frame()
  
  # Add positive subgroup data
  pos_pf <- traces[[biomarker]]$positive$p_pf
  pos_p <- traces[[biomarker]]$positive$p_p
  pos_d <- traces[[biomarker]]$positive$p_d
  
  max_cycles <- min(10 * 52 + 1, time_points_length)
  cycles <- 0:(time_points_length-1)
  
  pos_df <- data.frame(
    Year = cycles[1:max_cycles] / 52,
    PF = pos_pf[1:max_cycles],
    Progressed = pos_p[1:max_cycles],
    Dead = pos_d[1:max_cycles],
    Status = "Positive"
  )
  
  pos_df_long <- reshape2::melt(
    pos_df,
    id.vars = c("Year", "Status"),
    measure.vars = c("PF", "Progressed", "Dead"),
    variable.name = "State",
    value.name = "Proportion"
  )
  
  pos_neg_df <- rbind(pos_neg_df, pos_df_long)
  
  # Add negative subgroup data
  neg_pf <- traces[[biomarker]]$negative$p_pf
  neg_p <- traces[[biomarker]]$negative$p_p
  neg_d <- traces[[biomarker]]$negative$p_d
  
  max_cycles <- min(10 * 52 + 1, time_points_length)
  cycles <- 0:(time_points_length-1)
  
  neg_df <- data.frame(
    Year = cycles[1:max_cycles] / 52,
    PF = neg_pf[1:max_cycles],
    Progressed = neg_p[1:max_cycles],
    Dead = neg_d[1:max_cycles],
    Status = "Negative"
  )
  
  neg_df_long <- reshape2::melt(
    neg_df,
    id.vars = c("Year", "Status"),
    measure.vars = c("PF", "Progressed", "Dead"),
    variable.name = "State",
    value.name = "Proportion"
  )
  
  pos_neg_df <- rbind(pos_neg_df, neg_df_long)
  
  # Create faceted plot for this biomarker
  biomarker_title <- switch(biomarker,
                            "crp" = "CRP",
                            "tlr" = "TLR",
                            "tmb_braf" = "TMB/BRAF")
  
  pos_neg_plot <- ggplot(pos_neg_df, aes(x = Year, y = Proportion, fill = State)) +
    geom_area() +
    facet_wrap(~ Status, ncol = 2) +
    scale_fill_manual(values = c("PF" = "green", "Progressed" = "orange", "Dead" = "red")) +
    labs(title = paste0(biomarker_title, " Biomarker - State Occupancy by Status"),
         x = "Years",
         y = "Proportion of Cohort") +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 9, color = "black"),
      axis.text.y = element_text(size = 9, color = "black"),
      axis.title = element_text(size = 11),
      strip.text = element_text(size = 10),
      strip.background = element_rect(fill = "lightgray"),
      panel.grid.major = element_line(color = "gray80"),
      panel.grid.minor = element_line(color = "gray90"),
      plot.title = element_text(size = 14, hjust = 0.5),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10)
    ) +
    scale_x_continuous(
      breaks = seq(0, 10, by = 2),
      minor_breaks = seq(0, 10, by = 1),
      limits = c(0, 10),
      expand = c(0.01, 0)
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      minor_breaks = seq(0, 1, by = 0.1),
      expand = c(0.01, 0)
    )
  
  print(pos_neg_plot)
  
  # Save the biomarker plot
  filename <- paste0(biomarker, "_biomarker_state_occupancy.png")
  cat("Saving", biomarker, "biomarker plot as", filename, "...\n")
  ggsave(filename, plot = pos_neg_plot, width = 10, height = 6, dpi = 300)
}

# Create a directory for plots if it doesn't exist
if (!dir.exists("plots")) {
  cat("Creating 'plots' directory...\n")
  dir.create("plots")
  
  # Move files to the plots directory
  file.copy("all_strategies_state_occupancy.png", "plots/")
  file.remove("all_strategies_state_occupancy.png")
  
  for (biomarker in biomarkers) {
    filename <- paste0(biomarker, "_biomarker_state_occupancy.png")
    file.copy(filename, "plots/")
    file.remove(filename)
  }
  
  cat("All plot files moved to 'plots' directory\n")
}