# Compare Analysis Snapshots
# This script compares two single-model snapshots to assess the impact of bug fixes.

if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(here, dplyr, ggplot2, gridExtra)

source(here::here("scripts/R/functions/snapshot_utils.R"))

cat("\n=== Compare Analysis Snapshots ===\n\n")
cat("This script compares two snapshots to assess the impact of bug fixes.\n\n")

issue_number <- readline(prompt = "Enter GitHub issue number: ")

if (issue_number == "" || is.na(as.numeric(issue_number))) {
  stop("Invalid issue number. Please provide a numeric issue number.")
}

cat("\nSearching for snapshots for issue #", issue_number, "...\n", sep = "")
snapshot_pairs <- select_snapshots_for_comparison(issue_number)

cat("\nLoading snapshots...\n")

snapshots_dir <- here::here("data", "output", "snapshots")
before_snapshot <- load_snapshot(snapshot_pairs$before$snapshot, snapshots_dir)
after_snapshot <- load_snapshot(snapshot_pairs$after$snapshot, snapshots_dir)
before_psa <- load_snapshot(snapshot_pairs$before$psa, snapshots_dir)
after_psa <- load_snapshot(snapshot_pairs$after$psa, snapshots_dir)

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

before_icer <- before_snapshot$icer_obj
after_icer <- after_snapshot$icer_obj

comparison_base <- merge(
  before_icer[, c("Strategy", "Cost", "Effect", "Inc_Cost", "Inc_Effect", "ICER")],
  after_icer[, c("Strategy", "Cost", "Effect", "Inc_Cost", "Inc_Effect", "ICER")],
  by = "Strategy",
  suffixes = c("_before", "_after")
)

comparison_base$Cost_diff <- comparison_base$Cost_after - comparison_base$Cost_before
comparison_base$Effect_diff <- comparison_base$Effect_after - comparison_base$Effect_before
comparison_base$Inc_Cost_diff <- comparison_base$Inc_Cost_after - comparison_base$Inc_Cost_before
comparison_base$Inc_Effect_diff <- comparison_base$Inc_Effect_after - comparison_base$Inc_Effect_before
comparison_base$Cost_pct_change <- (comparison_base$Cost_diff / comparison_base$Cost_before) * 100
comparison_base$Effect_pct_change <- (comparison_base$Effect_diff / comparison_base$Effect_before) * 100
comparison_base$ICER_diff <- ifelse(
  is.na(comparison_base$ICER_before) | is.na(comparison_base$ICER_after),
  NA,
  comparison_base$ICER_after - comparison_base$ICER_before
)

print(comparison_base)

# 2. Compare NMB and optimal strategy
cat("\n2. Net Monetary Benefit Comparison\n")
cat("-----------------------------------\n")

before_nmb <- before_snapshot$nmb_at_wtp
after_nmb <- after_snapshot$nmb_at_wtp

comparison_nmb <- merge(
  before_nmb[, c("Strategy", "NMB")],
  after_nmb[, c("Strategy", "NMB")],
  by = "Strategy",
  suffixes = c("_before", "_after")
)

comparison_nmb$NMB_diff <- comparison_nmb$NMB_after - comparison_nmb$NMB_before
comparison_nmb <- comparison_nmb[order(-comparison_nmb$NMB_after), ]

print(comparison_nmb)

optimal_before <- before_nmb$Strategy[which.max(before_nmb$NMB)]
optimal_after <- after_nmb$Strategy[which.max(after_nmb$NMB)]

cat("\nOptimal strategy (before):", optimal_before, "\n")
cat("Optimal strategy (after): ", optimal_after, "\n")
if (optimal_before != optimal_after) {
  cat("*** OPTIMAL STRATEGY CHANGED ***\n")
}

# 3. Compare PSA summary statistics
cat("\n3. PSA Summary Statistics Comparison\n")
cat("-------------------------------------\n")

before_psa_sum <- before_snapshot$psa_summary
after_psa_sum <- after_snapshot$psa_summary

if (is.null(before_psa_sum) || is.null(after_psa_sum)) {
  cat("PSA summary unavailable in one or both snapshots.\n")
  comparison_psa <- data.frame()
} else {
  comparison_psa <- merge(
    before_psa_sum[, c("Strategy", "meanCost", "sdCost", "meanEffect", "sdEffect",
                       "Cost_2.5_percent", "Cost_97.5_percent",
                       "Effect_2.5_percent", "Effect_97.5_percent")],
    after_psa_sum[, c("Strategy", "meanCost", "sdCost", "meanEffect", "sdEffect",
                      "Cost_2.5_percent", "Cost_97.5_percent",
                      "Effect_2.5_percent", "Effect_97.5_percent")],
    by = "Strategy",
    suffixes = c("_before", "_after")
  )

  comparison_psa$meanCost_diff <- comparison_psa$meanCost_after - comparison_psa$meanCost_before
  comparison_psa$meanEffect_diff <- comparison_psa$meanEffect_after - comparison_psa$meanEffect_before

  print(comparison_psa[, c("Strategy", "meanCost_before", "meanCost_after", "meanCost_diff",
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
param_names <- union(names(before_params), names(after_params))
param_changes <- c()

for (param in param_names) {
  before_val <- before_params[[param]]
  after_val <- after_params[[param]]
  if (!identical(before_val, after_val)) {
    param_changes <- c(param_changes, param)
  }
}

if (length(param_changes) > 0) {
  cat("\nParameters that changed:\n")
  for (param in param_changes) {
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

get_psa_effect_matrix <- function(psa_obj) {
  if (!is.null(psa_obj$effect)) {
    as.matrix(psa_obj$effect)
  } else {
    as.matrix(psa_obj$effectiveness)
  }
}

before_cost_matrix <- as.matrix(before_psa$cost)
before_effect_matrix <- get_psa_effect_matrix(before_psa)
after_cost_matrix <- as.matrix(after_psa$cost)
after_effect_matrix <- get_psa_effect_matrix(after_psa)

before_psa_data <- data.frame(
  Strategy = rep(before_psa$strategies, each = nrow(before_cost_matrix)),
  Cost = as.vector(before_cost_matrix),
  Effect = as.vector(before_effect_matrix),
  Snapshot = "Before"
)

after_psa_data <- data.frame(
  Strategy = rep(after_psa$strategies, each = nrow(after_cost_matrix)),
  Cost = as.vector(after_cost_matrix),
  Effect = as.vector(after_effect_matrix),
  Snapshot = "After"
)

combined_psa_data <- rbind(before_psa_data, after_psa_data)

plot_before <- ggplot(before_psa_data, aes(x = Effect, y = Cost, color = Strategy)) +
  geom_point(alpha = 0.3, size = 1) +
  labs(title = "Before", x = "Effect (QALYs)", y = "Cost (EUR)") +
  theme_minimal() +
  theme(legend.position = "bottom")

plot_after <- ggplot(after_psa_data, aes(x = Effect, y = Cost, color = Strategy)) +
  geom_point(alpha = 0.3, size = 1) +
  labs(title = "After", x = "Effect (QALYs)", y = "Cost (EUR)") +
  theme_minimal() +
  theme(legend.position = "bottom")

plot_overlaid <- ggplot(combined_psa_data,
                        aes(x = Effect, y = Cost, color = Strategy, shape = Snapshot)) +
  geom_point(alpha = 0.3, size = 1.5) +
  scale_shape_manual(values = c("Before" = 1, "After" = 16)) +
  labs(title = "Before vs After Comparison",
       x = "Effect (QALYs)",
       y = "Cost (EUR)") +
  theme_minimal() +
  theme(legend.position = "bottom")

output_dir <- here::here("output", "temp")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cat("Saving plots to output/temp/...\n")

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
  base_case_comparison = comparison_base,
  nmb_comparison = comparison_nmb,
  psa_comparison = comparison_psa,
  optimal_strategy_before = optimal_before,
  optimal_strategy_after = optimal_after,
  parameter_changes = param_changes
)

comparison_filename <- paste0("comparison_issue", issue_number, "_",
                              before_snapshot$metadata$git_commit, "_vs_",
                              after_snapshot$metadata$git_commit, ".rds")
comparison_path <- file.path(output_dir, comparison_filename)

saveRDS(comparison_results, comparison_path)
cat("Saved comparison results to:", comparison_path, "\n")

cat("\n=== Comparison Complete ===\n")
cat("\nTo generate a formatted PDF report, run:\n")
cat("  quarto render scripts/QMD/technical_docs/bug_fix_impact.qmd\n")
cat("\n")
