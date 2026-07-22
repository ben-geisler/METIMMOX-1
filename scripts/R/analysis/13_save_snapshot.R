# Save Snapshot for Bug Fix Impact Analysis
#
# PURPOSE:
#   Track the impact of bug fixes on cost-effectiveness analysis results by
#   creating snapshots before (baseline) and after (fixed) the bug fix.
#   Saves results for the single joint economic survival model.
#
# USAGE:
#   Rscript 13_save_snapshot.R <issue_number> <baseline|fixed>
#
# WORKFLOW:
#   1. Before fixing bug: Rscript 13_save_snapshot.R 64 baseline
#      Alternative:
#      Set the command line arguments
#commandArgs <- function(trailingOnly = TRUE) {
#  if (trailingOnly) {
#    return(c("62", "baseline"))
#  }
#}
#      Then source the script
#source("scripts/R/analysis/13_save_snapshot.R")
#   2. Fix the bug and commit changes
#   3. After fixing bug:  Rscript 13_save_snapshot.R 64 fixed
#   4. Generate report:   quarto render scripts/QMD/technical_docs/bug_fix_impact.qmd
#
# OUTPUT FILES:
#   - snapshot_<issue>_<status>_<commit>.rds  (base case + metadata)
#   - psa_<issue>_<status>_<commit>.rds       (PSA object, if available)
#
# NOTES:
#   - Runs base case for the single joint economic model
#   - Loads pre-computed PSA cache (psa_obj_{ipd|correct}.rds)
#   - Uses sampling cache if available (much faster on subsequent runs)
#   - Snapshots saved to data/output/snapshots/

# Load required packages
if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(here, dampack)

# Load snapshot utilities
source(here::here("scripts/R/functions/snapshot_utils.R"))

# Get GitHub issue number and status from command line or prompt
cat("\n=== Save Analysis Snapshot ===\n\n")
cat("This script will run the complete analysis and save a snapshot.\n")
cat("The snapshot will be used to compare results before and after bug fixes.\n\n")
cat("Usage: Rscript 13_save_snapshot.R <issue_number> <baseline|fixed>\n\n")

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

# ============================================================================
# Step 1: Run model-independent setup
# ============================================================================
cat("=== Running Setup Scripts ===\n")

cat("\n[1/5] Running 02_setup_and_global_variables.R...\n")
source(here::here("scripts/R/analysis/02_setup_and_global_variables.R"))

cat("[2/5] Running 03_biomarker_strategies.R...\n")
source(here::here("scripts/R/analysis/03_biomarker_strategies.R"))

cat("[3/5] Running 04_parametric_survival_analysis.R...\n")
source(here::here("scripts/R/analysis/04_parametric_survival_analysis.R"))

cat("[4/5] Running 05_basecase_input_parameters.R...\n")
source(here::here("scripts/R/analysis/05_basecase_input_parameters.R"))

cat("[5/5] Running 06_sampling.R (may take time if cache doesn't exist)...\n")
source(here::here("scripts/R/analysis/06_sampling.R"))

# ============================================================================
# Step 2: Source helper functions
# ============================================================================
cat("\n=== Loading Helper Functions ===\n")
source(here::here("scripts/R/functions/model_fun.R"))
source(here::here("scripts/R/functions/calculate_outcomes.R"))
source(here::here("scripts/R/functions/cea_helpers.R"))
source(here::here("scripts/R/functions/psa_functions.R"))
source(here::here("scripts/R/functions/snapshot_utils.R"))

# ============================================================================
# Step 3: Run analysis for the single economic model
# ============================================================================
cat("\n=== Running Single-Model Analysis ===\n")

# Helper: compute PSA summary with CIs from a psa object
compute_psa_summary <- function(psa_obj) {
  psa_summary_full <- summary(psa_obj, calc_sds = TRUE,
                              prob = c(0.025, 0.975))
  cost_matrix <- as.matrix(psa_obj$cost)
  effect_matrix <- if (!is.null(psa_obj$effect)) {
    as.matrix(psa_obj$effect)
  } else {
    as.matrix(psa_obj$effectiveness)
  }
  cost_q <- apply(cost_matrix, 2, quantile, probs = c(0.025, 0.975))
  effect_q <- apply(effect_matrix, 2, quantile, probs = c(0.025, 0.975))
  data.frame(
    Strategy = psa_summary_full$Strategy,
    meanCost = psa_summary_full$meanCost,
    sdCost = psa_summary_full$sdCost,
    meanEffect = psa_summary_full$meanEffect,
    sdEffect = psa_summary_full$sdEffect,
    Cost_2.5_percent = cost_q[1, ],
    Cost_97.5_percent = cost_q[2, ],
    Effect_2.5_percent = effect_q[1, ],
    Effect_97.5_percent = effect_q[2, ],
    stringsAsFactors = FALSE
  )
}

# Run base case
base_result <- run_basecase(verbose = TRUE)
base_results <- base_result$base_results
icer_obj <- base_result$icer_obj

# Calculate NMB
nmb_at_wtp <- data.frame(
  Strategy = base_results$Strategy,
  Cost = base_results$Cost,
  Effect = base_results$Effect,
  NMB = base_results$Effect * WTP - base_results$Cost
)
nmb_at_wtp <- nmb_at_wtp[order(-nmb_at_wtp$NMB), ]

# Load PSA cache, or generate if not available
cache_file_obj <- psa_obj_path(utility_source_label)
cache_file_params <- psa_params_path(utility_source_label)

psa_obj <- NULL
psa_summary <- NULL

psa_obj <- load_psa_cache(util_label = utility_source_label, verbose = TRUE)

if (is.null(psa_obj)) {
  if (exists("sampling_models") && !is.null(sampling_models)) {
    cat("PSA cache not found - generating PSA...\n")

    psa_build <- build_psa_obj(
      l_params_base = l_params_base,
      param_distributions = param_distributions,
      strategies = strategies,
      time_horizon = time_horizon,
      cl = cl,
      n_sim = n_sim,
      seed = analysis_seed,
      currency = "EUR"
    )
    psa_obj <- psa_build$psa_obj
    psa_params_gen <- psa_build$psa_params
    rm(psa_build)

    tryCatch({
      saveRDS(psa_obj, cache_file_obj)
      saveRDS(psa_params_gen, cache_file_params)
      cat("PSA cache saved to:", cache_file_obj, "\n")
    }, error = function(e) {
      cat("Warning: Failed to save PSA cache:", e$message, "\n")
    })
  } else {
    cat("WARNING: PSA cache not found and sampling_models is unavailable.\n")
    cat("Snapshot will be saved without PSA summary.\n")
  }
}

if (!is.null(psa_obj)) {
  psa_summary <- compute_psa_summary(psa_obj)
}

cat("\n=== Single-Model Analysis Complete ===\n")

# ============================================================================
# Step 4: Collect metadata and build snapshot
# ============================================================================
cat("\nCollecting metadata...\n")
metadata <- collect_metadata(issue_number)

cat("Creating snapshot object...\n")
snapshot <- list(
  metadata = metadata,
  base_results = base_results,
  icer_obj = as.data.frame(icer_obj),
  nmb_at_wtp = nmb_at_wtp,
  psa_summary = psa_summary
)

# ============================================================================
# Step 5: Save files
# ============================================================================
commit <- get_git_commit()
snapshot_filename <- paste0("snapshot_", issue_number, "_",
                            snapshot_status, "_", commit, ".rds")
psa_filename <- paste0("psa_", issue_number, "_",
                        snapshot_status, "_", commit, ".rds")
cat("\nSaving snapshot file...\n")
snapshot_file <- snapshot_path(snapshot_filename)
saveRDS(snapshot, snapshot_file)
cat("  Saved:", snapshot_file, "\n")

cat("Saving PSA object...\n")
psa_path <- snapshot_path(psa_filename)
if (!is.null(psa_obj)) {
  saveRDS(psa_obj, psa_path)
  cat("  Saved:", psa_path, "\n")
} else {
  cat("  WARNING: No PSA object available - PSA file not saved\n")
}

# ============================================================================
# Print summary
# ============================================================================
cat("\n=== Snapshot Saved Successfully ===\n")
cat("Issue number:", issue_number, "\n")
cat("Git commit:", commit, "\n")
cat("Timestamp:", metadata$timestamp, "\n")
cat("Snapshot file:", snapshot_filename, "\n")
cat("PSA file:", psa_filename, "\n")

cat("\nSnapshot contains:\n")
cat("  - Base case results: ", nrow(base_results), " strategies",
    if (!is.null(psa_summary)) " + PSA" else "",
    "\n", sep = "")
cat("  - Metadata (git commit, R version, package versions, parameters)\n")

if (snapshot_status == "baseline") {
  cat("\nNext steps:\n")
  cat("  1. Fix the bug and commit your changes\n")
  cat("  2. Run: Rscript scripts/R/analysis/13_save_snapshot.R",
      issue_number, "fixed\n")
  cat("  3. Generate report: quarto render",
      "scripts/QMD/technical_docs/bug_fix_impact.qmd\n")
} else {
  cat("\nTo generate the bug fix impact report, run:\n")
  cat("  quarto render scripts/QMD/technical_docs/bug_fix_impact.qmd\n")
}

cat("\n=== Done ===\n\n")
