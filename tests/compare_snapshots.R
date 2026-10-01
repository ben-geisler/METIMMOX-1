# Compare Analysis Snapshots
# This script compares two single-model snapshots to assess the impact of bug fixes.

if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(here, dplyr, ggplot2, gridExtra)

source(here::here("R/snapshot_utils.R"))

cat("\n=== Compare Analysis Snapshots ===\n\n")
cat("This script compares two snapshots to assess the impact of bug fixes.\n\n")

v_args <- commandArgs(trailingOnly = TRUE)
issue_number <- if (length(v_args)) v_args[1] else readline(prompt = "Enter GitHub issue number: ")

if (issue_number == "" || is.na(as.numeric(issue_number))) {
  stop("Invalid issue number. Please provide a numeric issue number.")
}

cat("\nSearching for snapshots for issue #", issue_number, "...\n", sep = "")
snapshot_pairs <- select_snapshots_for_comparison(issue_number)

cat("\nLoading snapshots...\n")

before_snapshot <- load_snapshot(snapshot_pairs$before$snapshot)
after_snapshot <- load_snapshot(snapshot_pairs$after$snapshot)
before_psa <- load_snapshot(snapshot_pairs$before$psa)
after_psa <- load_snapshot(snapshot_pairs$after$psa)

if (!inherits(before_psa, "psa") || !inherits(after_psa, "psa")) {
  stop("Snapshot comparison expects single-model dampack PSA objects.")
}

cat("  Before: ", snapshot_pairs$before$snapshot, "\n", sep = "")
cat("          Commit: ", before_snapshot$metadata$git_commit,
    ", Timestamp: ", before_snapshot$metadata$timestamp, "\n", sep = "")
cat("  After:  ", snapshot_pairs$after$snapshot, "\n", sep = "")
cat("          Commit: ", after_snapshot$metadata$git_commit,
    ", Timestamp: ", after_snapshot$metadata$timestamp, "\n", sep = "")

cat("\n=== Comparison Analysis ===\n\n")

# 1. Compare base case results
cat("1. Base Case Results Comparison\n")
cat("--------------------------------\n")

df_icer_before <- before_snapshot$icer_obj
df_icer_after <- after_snapshot$icer_obj

df_comparison_base <- merge(
  df_icer_before[, c("Strategy", "Cost", "Effect", "Inc_Cost", "Inc_Effect", "ICER")],
  df_icer_after[, c("Strategy", "Cost", "Effect", "Inc_Cost", "Inc_Effect", "ICER")],
  by = "Strategy",
  suffixes = c("_before", "_after")
)

df_comparison_base$Cost_diff <- df_comparison_base$Cost_after - df_comparison_base$Cost_before
df_comparison_base$Effect_diff <- df_comparison_base$Effect_after - df_comparison_base$Effect_before
df_comparison_base$Inc_Cost_diff <- df_comparison_base$Inc_Cost_after - df_comparison_base$Inc_Cost_before
df_comparison_base$Inc_Effect_diff <- df_comparison_base$Inc_Effect_after - df_comparison_base$Inc_Effect_before
df_comparison_base$Cost_pct_change <- (df_comparison_base$Cost_diff / df_comparison_base$Cost_before) * 100
df_comparison_base$Effect_pct_change <- (df_comparison_base$Effect_diff / df_comparison_base$Effect_before) * 100
df_comparison_base$ICER_diff <- ifelse(
  is.na(df_comparison_base$ICER_before) | is.na(df_comparison_base$ICER_after),
  NA,
  df_comparison_base$ICER_after - df_comparison_base$ICER_before
)

print(df_comparison_base)

# 2. Compare NMB and optimal strategy
cat("\n2. Net Monetary Benefit Comparison\n")
cat("-----------------------------------\n")

df_nmb_before <- before_snapshot$nmb_at_wtp
df_nmb_after <- after_snapshot$nmb_at_wtp

df_comparison_nmb <- merge(
  df_nmb_before[, c("Strategy", "NMB")],
  df_nmb_after[, c("Strategy", "NMB")],
  by = "Strategy",
  suffixes = c("_before", "_after")
)

df_comparison_nmb$NMB_diff <- df_comparison_nmb$NMB_after - df_comparison_nmb$NMB_before
df_comparison_nmb <- df_comparison_nmb[order(-df_comparison_nmb$NMB_after), ]

print(df_comparison_nmb)

optimal_before <- df_nmb_before$Strategy[which.max(df_nmb_before$NMB)]
optimal_after <- df_nmb_after$Strategy[which.max(df_nmb_after$NMB)]

cat("\nOptimal strategy (before):", optimal_before, "\n")
cat("Optimal strategy (after): ", optimal_after, "\n")
if (optimal_before != optimal_after) {
  cat("*** OPTIMAL STRATEGY CHANGED ***\n")
}

# 3. Compare PSA summary statistics
cat("\n3. PSA Summary Statistics Comparison\n")
cat("-------------------------------------\n")

df_psa_sum_before <- before_snapshot$psa_summary
df_psa_sum_after <- after_snapshot$psa_summary

# Provenance (issue #156): say whether the two PSA files are the same cache.
cat(describe_psa_provenance(before_snapshot, after_snapshot), "\n\n")

if (is.null(df_psa_sum_before) || is.null(df_psa_sum_after)) {
  cat("PSA summary unavailable in one or both snapshots.\n")
  df_comparison_psa <- data.frame()
} else {
  df_comparison_psa <- merge(
    df_psa_sum_before[, c("Strategy", "meanCost", "sdCost", "meanEffect", "sdEffect",
                       "Cost_2.5_percent", "Cost_97.5_percent",
                       "Effect_2.5_percent", "Effect_97.5_percent")],
    df_psa_sum_after[, c("Strategy", "meanCost", "sdCost", "meanEffect", "sdEffect",
                      "Cost_2.5_percent", "Cost_97.5_percent",
                      "Effect_2.5_percent", "Effect_97.5_percent")],
    by = "Strategy",
    suffixes = c("_before", "_after")
  )

  df_comparison_psa$meanCost_diff <- df_comparison_psa$meanCost_after - df_comparison_psa$meanCost_before
  df_comparison_psa$meanEffect_diff <- df_comparison_psa$meanEffect_after - df_comparison_psa$meanEffect_before

  print(df_comparison_psa[, c("Strategy", "meanCost_before", "meanCost_after", "meanCost_diff",
                           "meanEffect_before", "meanEffect_after", "meanEffect_diff")])
}

# 4. Compare metadata
cat("\n4. Metadata Comparison\n")
cat("----------------------\n")

cat("Git commits:\n")
cat("  Before:", before_snapshot$metadata$git_commit, "\n")
cat("  After: ", after_snapshot$metadata$git_commit, "\n")

cat("\nTimestamps:\n")
cat("  Before:", before_snapshot$metadata$timestamp, "\n")
cat("  After: ", after_snapshot$metadata$timestamp, "\n")

cat("\nR versions:\n")
cat("  Before:", before_snapshot$metadata$r_version, "\n")
cat("  After: ", after_snapshot$metadata$r_version, "\n")

before_params <- before_snapshot$metadata$parameters
after_params <- after_snapshot$metadata$parameters
v_param_names <- union(names(before_params), names(after_params))
v_param_changes <- c()

for (param in v_param_names) {
  before_val <- before_params[[param]]
  after_val <- after_params[[param]]
  if (!identical(before_val, after_val)) {
    v_param_changes <- c(v_param_changes, param)
  }
}

if (length(v_param_changes) > 0) {
  cat("\nParameters that changed:\n")
  for (param in v_param_changes) {
    cat("  ", param, ": ",
        before_params[[param]], " -> ",
        after_params[[param]], "\n", sep = "")
  }
} else {
  cat("\nNo parameter changes detected.\n")
}

# 5. Create PSA scatter plots
cat("\n5. Creating PSA Scatter Plots\n")
cat("------------------------------\n")

#' Effect matrix of a dampack PSA object, reading `$effect` or the legacy
#' `$effectiveness` element. Returns a numeric matrix (draws x strategies).
get_psa_effect_matrix <- function(psa_obj) {
  if (!is.null(psa_obj$effect)) {
    as.matrix(psa_obj$effect)
  } else {
    as.matrix(psa_obj$effectiveness)
  }
}

m_cost_before <- as.matrix(before_psa$cost)
m_effect_before <- get_psa_effect_matrix(before_psa)
m_cost_after <- as.matrix(after_psa$cost)
m_effect_after <- get_psa_effect_matrix(after_psa)

df_psa_before <- data.frame(
  Strategy = rep(before_psa$strategies, each = nrow(m_cost_before)),
  Cost = as.vector(m_cost_before),
  Effect = as.vector(m_effect_before),
  Snapshot = "Before"
)

df_psa_after <- data.frame(
  Strategy = rep(after_psa$strategies, each = nrow(m_cost_after)),
  Cost = as.vector(m_cost_after),
  Effect = as.vector(m_effect_after),
  Snapshot = "After"
)

df_psa_combined <- rbind(df_psa_before, df_psa_after)

plot_before <- ggplot(df_psa_before, aes(x = Effect, y = Cost, color = Strategy)) +
  geom_point(alpha = 0.3, size = 1) +
  labs(title = "Before", x = "Effect (QALYs)", y = "Cost (EUR)") +
  theme_minimal() +
  theme(legend.position = "bottom")

plot_after <- ggplot(df_psa_after, aes(x = Effect, y = Cost, color = Strategy)) +
  geom_point(alpha = 0.3, size = 1) +
  labs(title = "After", x = "Effect (QALYs)", y = "Cost (EUR)") +
  theme_minimal() +
  theme(legend.position = "bottom")

plot_overlaid <- ggplot(df_psa_combined,
                        aes(x = Effect, y = Cost, color = Strategy, shape = Snapshot)) +
  geom_point(alpha = 0.3, size = 1.5) +
  scale_shape_manual(values = c("Before" = 1, "After" = 16)) +
  labs(title = "Before vs After Comparison",
       x = "Effect (QALYs)",
       y = "Cost (EUR)") +
  theme_minimal() +
  theme(legend.position = "bottom")

output_dir <- here::here("outputs", "temp")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cat("Saving plots to outputs/temp/...\n")

ggsave(
  filename = file.path(output_dir, paste0("psa_comparison_sidebyside_issue", issue_number, ".png")),
  plot = gridExtra::grid.arrange(plot_before, plot_after, ncol = 2),
  width = 12, height = 5, dpi = 300
)

ggsave(
  filename = file.path(output_dir, paste0("psa_comparison_overlaid_issue", issue_number, ".png")),
  plot = plot_overlaid,
  width = 8, height = 6, dpi = 300
)

cat("  Saved: psa_comparison_sidebyside_issue", issue_number, ".png\n", sep = "")
cat("  Saved: psa_comparison_overlaid_issue", issue_number, ".png\n", sep = "")

# 6. Save comparison results
cat("\n6. Saving Comparison Results\n")
cat("----------------------------\n")

comparison_results <- list(
  issue_number = issue_number,
  before_metadata = before_snapshot$metadata,
  after_metadata = after_snapshot$metadata,
  base_case_comparison = df_comparison_base,
  nmb_comparison = df_comparison_nmb,
  psa_comparison = df_comparison_psa,
  optimal_strategy_before = optimal_before,
  optimal_strategy_after = optimal_after,
  parameter_changes = v_param_changes
)

comparison_filename <- paste0("comparison_issue", issue_number, "_",
                              before_snapshot$metadata$git_commit, "_vs_",
                              after_snapshot$metadata$git_commit, ".rds")
comparison_path <- file.path(output_dir, comparison_filename)

saveRDS(comparison_results, comparison_path)
cat("Saved comparison results to:", comparison_path, "\n")

cat("\n=== Comparison Complete ===\n")
cat("\nTo generate a formatted PDF report, run:\n")
cat("  quarto render reports/technical/bug_fix_impact.qmd\n")
cat("\n")
