# Load required packages
if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(here, ggplot2, reshape2)

# Load functions
source(here::here("R/model_fun.R"))
source(here::here("R/calculate_outcomes.R"))

# Generate traces by running the model with return_traces = TRUE
cat("Generating state occupancy traces for all strategies...\n")

# Save generated trace figures in a predictable location
figs_dir <- here::here("outputs", "figs")
if (!dir.exists(figs_dir)) {
  cat("Creating 'figs' directory...\n")
  dir.create(figs_dir, recursive = TRUE)
}

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
df_all_strategies <- data.frame()

for(strategy in strategies) {
  # Get the relevant traces
  p_pf <- traces[[strategy]]$p_pf
  p_p <- traces[[strategy]]$p_p
  p_d <- traces[[strategy]]$p_d
  
  # Extract data for plotting
  max_cycles <- min(10 * 52 + 1, time_points_length)
  v_cycles <- 0:(time_points_length-1)
  
  # Create data frame
  df_strategy_trace <- data.frame(
    Year = v_cycles[1:max_cycles] / 52,
    PF = p_pf[1:max_cycles],
    Progressed = p_p[1:max_cycles],
    Dead = p_d[1:max_cycles],
    Strategy = strategy
  )
  
  # Convert to long format
  df_long <- reshape2::melt(
    df_strategy_trace,
    id.vars = c("Year", "Strategy"),
    measure.vars = c("PF", "Progressed", "Dead"),
    variable.name = "State",
    value.name = "Proportion"
  )
  
  # Append to all strategies dataframe
  df_all_strategies <- rbind(df_all_strategies, df_long)
}

# Create strategy labels for faceting from keyed configuration metadata. Kept
# local: the global strategy_labels (R/report_format.R) maps id -> report label
# and must not be overwritten by this script (issue #187).
v_trace_labels <- get_strategy_metadata(strategies)$trace_label
df_all_strategies$Strategy <- factor(df_all_strategies$Strategy,
                                     levels = strategies,
                                     labels = v_trace_labels)

# Create faceted plot
all_strategies_plot <- ggplot(df_all_strategies, aes(x = Year, y = Proportion, fill = State)) +
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
all_strategies_filename <- file.path(figs_dir, "traces_all_strategies_state_occupancy.png")
cat("Saving all strategies plot as", all_strategies_filename, "...\n")
ggsave(all_strategies_filename, plot = all_strategies_plot,
       width = 10, height = 8, dpi = 300)

# Create faceted plots for biomarker positive and negative status for each biomarker
for(biomarker in biomarkers) {
  cat(paste0("Creating faceted plot for ", biomarker, " biomarker subgroups...\n"))
  
  # Prepare data for positive and negative subgroups
  df_pos_neg <- data.frame()
  
  # Add positive subgroup data
  v_pos_pf <- traces[[biomarker]]$positive$p_pf
  v_pos_p <- traces[[biomarker]]$positive$p_p
  v_pos_d <- traces[[biomarker]]$positive$p_d
  
  max_cycles <- min(10 * 52 + 1, time_points_length)
  v_cycles <- 0:(time_points_length-1)
  
  df_pos <- data.frame(
    Year = v_cycles[1:max_cycles] / 52,
    PF = v_pos_pf[1:max_cycles],
    Progressed = v_pos_p[1:max_cycles],
    Dead = v_pos_d[1:max_cycles],
    Status = "Positive"
  )
  
  df_pos_long <- reshape2::melt(
    df_pos,
    id.vars = c("Year", "Status"),
    measure.vars = c("PF", "Progressed", "Dead"),
    variable.name = "State",
    value.name = "Proportion"
  )
  
  df_pos_neg <- rbind(df_pos_neg, df_pos_long)
  
  # Add negative subgroup data
  v_neg_pf <- traces[[biomarker]]$negative$p_pf
  v_neg_p <- traces[[biomarker]]$negative$p_p
  v_neg_d <- traces[[biomarker]]$negative$p_d
  
  max_cycles <- min(10 * 52 + 1, time_points_length)
  v_cycles <- 0:(time_points_length-1)
  
  df_neg <- data.frame(
    Year = v_cycles[1:max_cycles] / 52,
    PF = v_neg_pf[1:max_cycles],
    Progressed = v_neg_p[1:max_cycles],
    Dead = v_neg_d[1:max_cycles],
    Status = "Negative"
  )
  
  df_neg_long <- reshape2::melt(
    df_neg,
    id.vars = c("Year", "Status"),
    measure.vars = c("PF", "Progressed", "Dead"),
    variable.name = "State",
    value.name = "Proportion"
  )
  
  df_pos_neg <- rbind(df_pos_neg, df_neg_long)
  
  # Create faceted plot for this biomarker
  biomarker_title <- switch(biomarker,
                            "crp" = "CRP",
                            "tlr" = "TLR",
                            "tmb_braf" = "TMB/BRAF")
  
  pos_neg_plot <- ggplot(df_pos_neg, aes(x = Year, y = Proportion, fill = State)) +
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
  filename <- file.path(figs_dir, paste0("traces_", biomarker, "_biomarker_state_occupancy.png"))
  cat("Saving", biomarker, "biomarker plot as", filename, "...\n")
  ggsave(filename, plot = pos_neg_plot, width = 10, height = 6, dpi = 300)
}
