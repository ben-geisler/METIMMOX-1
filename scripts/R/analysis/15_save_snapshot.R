# Save Snapshot for Bug Fix Impact Analysis
#
# PURPOSE:
#   Track the impact of bug fixes on cost-effectiveness analysis results by
#   creating snapshots before (baseline) and after (fixed) the bug fix.
#
# USAGE:
#   Rscript 15_save_snapshot.R <issue_number> <baseline|fixed>
#
# WORKFLOW:
#   1. Before fixing bug: Rscript 15_save_snapshot.R 64 baseline
#      Alternative:
#      Set the command line arguments
#commandArgs <- function(trailingOnly = TRUE) {
#  if (trailingOnly) {
#    return(c("62", "baseline"))
#  }
#}
#      Then source the script
#source("scripts/R/analysis/15_save_snapshot.R")
#   2. Fix the bug and commit changes
#   3. After fixing bug:  Rscript 15_save_snapshot.R 64 fixed
#   4. Generate report:   quarto render scripts/QMD/technical_docs/bug_fix_impact.qmd
#
# OUTPUT FILES:
#   - snapshot_<issue>_<status>_<commit>.rds  (base case + metadata)
#   - psa_<issue>_<status>_<commit>.rds       (PSA object with 5000 sims)
#
# NOTES:
#   - Runs scripts 02, 03, 06, 07, 08, 10, 12 (takes 20-60 min for PSA)
#   - Uses sampling cache if available (much faster on subsequent runs)
#   - Snapshots saved to data/tidy/snapshots/

# Load required packages
if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(here)

# Load snapshot utilities
source(here::here("scripts/R/functions/snapshot_utils.R"))

# Get GitHub issue number and status from command line or prompt
cat("\n=== Save Analysis Snapshot ===\n\n")
cat("This script will run the complete analysis and save a snapshot.\n")
cat("The snapshot will be used to compare results before and after bug fixes.\n\n")
cat("Usage: Rscript 15_save_snapshot.R <issue_number> <baseline|fixed>\n\n")

# Try to get issue number and status from command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2) {
  issue_number <- args[1]
  snapshot_status <- args[2]
  cat("Issue number:", issue_number, "\n")
  cat("Snapshot status:", snapshot_status, "\n\n")
} else if (length(args) == 1) {
  issue_number <- args[1]
  snapshot_status <- readline(prompt = "Enter snapshot status (baseline/fixed): ")
  cat("\n")
} else {
  issue_number <- readline(prompt = "Enter GitHub issue number: ")
  snapshot_status <- readline(prompt = "Enter snapshot status (baseline/fixed): ")
  cat("\n")
}

if (issue_number == "" || is.na(suppressWarnings(as.numeric(issue_number)))) {
  stop("Invalid issue number. Please provide a numeric issue number.")
}

if (!snapshot_status %in% c("baseline", "fixed")) {
  stop("Invalid snapshot status. Must be 'baseline' or 'fixed'.")
}

cat("=== Running Analysis Scripts ===\n")

# Save parameters as environment variables (survives rm(list = ls()))
Sys.setenv(SNAPSHOT_ISSUE_NUMBER = issue_number)
Sys.setenv(SNAPSHOT_STATUS = snapshot_status)

# Run analysis scripts in sequence
cat("\n[1/7] Running 02_setup_and_global_variables.R...\n")
source(here::here("scripts/R/analysis/02_setup_and_global_variables.R"))

cat("[2/7] Running 03_biomarker_strategies.R...\n")
source(here::here("scripts/R/analysis/03_biomarker_strategies.R"))

cat("[3/7] Running 06_parametric_survival analysis.R...\n")
source(here::here("scripts/R/analysis/06_parametric_survival analysis.R"))

cat("[4/7] Running 07_basecase_input_parameters.R...\n")
source(here::here("scripts/R/analysis/07_basecase_input_parameters.R"))

cat("[5/7] Running 08_sampling.R (may take time if cache doesn't exist)...\n")
source(here::here("scripts/R/analysis/08_sampling.R"))

cat("[6/7] Running 10_basecase_analysis.R...\n")
source(here::here("scripts/R/analysis/10_basecase_analysis.R"))

cat("[7/7] Running 12_PSA.R...\n")
source(here::here("scripts/R/analysis/12_PSA.R"))

cat("\n=== Analysis Complete ===\n\n")

# Reload snapshot utilities and parameters (script 02 clears workspace)
source(here::here("scripts/R/functions/snapshot_utils.R"))
issue_number <- Sys.getenv("SNAPSHOT_ISSUE_NUMBER")
snapshot_status <- Sys.getenv("SNAPSHOT_STATUS")

# Check that required objects exist
cat("Validating required objects...\n")
required_objects <- c("base_results", "icer_obj", "nmb_at_wtp", "psa_obj")
missing_objects <- setdiff(required_objects, ls())

if (length(missing_objects) > 0) {
  stop("Missing required objects: ", paste(missing_objects, collapse = ", "),
       "\nAnalysis may have failed. Check for errors above.")
}

# Collect metadata
cat("Collecting metadata...\n")
metadata <- collect_metadata(issue_number)

# Get PSA summary with confidence intervals
cat("Calculating PSA summary statistics...\n")
psa_summary_full <- summary(psa_obj, calc_sds = TRUE, prob = c(0.025, 0.975))

# Calculate 95% credible intervals from PSA object directly
cost_quantiles <- apply(psa_obj$cost, 2, quantile, probs = c(0.025, 0.975))
effect_quantiles <- apply(psa_obj$effectiveness, 2, quantile, probs = c(0.025, 0.975))

# Create summary with credible intervals
psa_summary_with_ci <- data.frame(
  Strategy = psa_summary_full$Strategy,
  meanCost = psa_summary_full$meanCost,
  sdCost = psa_summary_full$sdCost,
  meanEffect = psa_summary_full$meanEffect,
  sdEffect = psa_summary_full$sdEffect,
  Cost_2.5_percent = cost_quantiles[1, ],
  Cost_97.5_percent = cost_quantiles[2, ],
  Effect_2.5_percent = effect_quantiles[1, ],
  Effect_97.5_percent = effect_quantiles[2, ],
  stringsAsFactors = FALSE
)

# Create snapshot object
cat("Creating snapshot object...\n")
snapshot <- list(
  metadata = metadata,
  base_results = base_results,
  icer_obj = as.data.frame(icer_obj),  # Ensure it's a data frame
  nmb_at_wtp = nmb_at_wtp,
  psa_summary = psa_summary_with_ci
)

# Generate filenames with baseline/fixed status
commit <- get_git_commit()

snapshot_filename <- paste0("snapshot_", issue_number, "_", snapshot_status, "_", commit, ".rds")
psa_filename <- paste0("psa_", issue_number, "_", snapshot_status, "_", commit, ".rds")

snapshots_dir <- here::here("data", "tidy", "snapshots")

# Save snapshot file
cat("\nSaving snapshot file...\n")
snapshot_path <- file.path(snapshots_dir, snapshot_filename)
saveRDS(snapshot, snapshot_path)
cat("  Saved:", snapshot_path, "\n")

# Save PSA object
cat("Saving PSA object...\n")
psa_path <- file.path(snapshots_dir, psa_filename)
saveRDS(psa_obj, psa_path)
cat("  Saved:", psa_path, "\n")

# Print summary
cat("\n=== Snapshot Saved Successfully ===\n")
cat("Issue number:", issue_number, "\n")
cat("Git commit:", commit, "\n")
cat("Timestamp:", metadata$timestamp, "\n")
cat("Snapshot file:", snapshot_filename, "\n")
cat("PSA file:", psa_filename, "\n")

cat("\nSnapshot contains:\n")
cat("  - Base case results (", nrow(base_results), " strategies)\n", sep = "")
cat("  - ICER analysis\n")
cat("  - Net Monetary Benefit at WTP = €", format(WTP, big.mark = ","), "\n", sep = "")
cat("  - PSA summary (", psa_obj$n_sim, " simulations)\n", sep = "")
cat("  - Metadata (git commit, R version, package versions, parameters)\n")

if (snapshot_status == "baseline") {
  cat("\nNext steps:\n")
  cat("  1. Fix the bug and commit your changes\n")
  cat("  2. Run: Rscript scripts/R/analysis/15_save_snapshot.R", issue_number, "fixed\n")
  cat("  3. Generate report: quarto render scripts/QMD/technical_docs/bug_fix_impact.qmd -P issue:", issue_number, "\n")
} else {
  cat("\nTo generate the bug fix impact report, run:\n")
  cat("  quarto render scripts/QMD/technical_docs/bug_fix_impact.qmd -P issue:", issue_number, "\n")
}

cat("\n=== Done ===\n\n")
