# Save Snapshot for Bug Fix Impact Analysis
#
# PURPOSE:
#   Track the impact of bug fixes on cost-effectiveness analysis results by
#   creating snapshots before (baseline) and after (fixed) the bug fix.
#   Saves results for ALL three model structures (A=joint, B=focused, C=separate).
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
#   - snapshot_<issue>_<status>_<commit>.rds  (base case + metadata for all 3 models)
#   - psa_<issue>_<status>_<commit>.rds       (PSA objects for all 3 models)
#
# NOTES:
#   - Runs base case for all 3 model structures (joint, focused, separate)
#   - Loads pre-computed PSA caches (psa_obj_joint.rds, etc.)
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

# Save parameters as environment variables (survives rm(list = ls()))
Sys.setenv(SNAPSHOT_ISSUE_NUMBER = issue_number)
Sys.setenv(SNAPSHOT_STATUS = snapshot_status)

# ============================================================================
# Step 1: Run model-independent setup
# ============================================================================
cat("=== Running Setup Scripts ===\n")

cat("\n[1/3] Running 02_setup_and_global_variables.R...\n")
source(here::here("scripts/R/analysis/02_setup_and_global_variables.R"))

cat("[2/3] Running 03_biomarker_strategies.R...\n")
source(here::here("scripts/R/analysis/03_biomarker_strategies.R"))

cat("[3/3] Running 08_sampling.R (may take time if cache doesn't exist)...\n")
source(here::here("scripts/R/analysis/08_sampling.R"))

# Save WTP for NMB calculation (survives re-sourcing)
Sys.setenv(SNAPSHOT_WTP = as.character(WTP))

# ============================================================================
# Step 2: Source helper functions
# ============================================================================
cat("\n=== Loading Helper Functions ===\n")
source(here::here("scripts/R/functions/model_fun.R"))
source(here::here("scripts/R/functions/calculate_outcomes.R"))
source(here::here("scripts/R/functions/multi_model_cea.R"))
source(here::here("scripts/R/functions/snapshot_utils.R"))

# Recover environment variables
issue_number <- Sys.getenv("SNAPSHOT_ISSUE_NUMBER")
snapshot_status <- Sys.getenv("SNAPSHOT_STATUS")
wtp_val <- as.numeric(Sys.getenv("SNAPSHOT_WTP"))

# ============================================================================
# Step 3: Run analysis for all 3 model structures
# ============================================================================
cat("\n=== Running Multi-Model Analysis ===\n")

model_labels <- c("joint", "focused", "separate")
models_results <- list()
psa_objects <- list()

# Helper: compute PSA summary with CIs from a psa object
compute_psa_summary <- function(psa_obj) {
  psa_summary_full <- summary(psa_obj, calc_sds = TRUE,
                              prob = c(0.025, 0.975))
  cost_q <- apply(psa_obj$cost, 2, quantile, probs = c(0.025, 0.975))
  effect_q <- apply(psa_obj$effectiveness, 2, quantile,
                    probs = c(0.025, 0.975))
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

for (i in seq_along(model_labels)) {
  ms <- i - 1  # MODEL_STRUCTURE: 0, 1, 2
  label <- model_labels[i]
  cat("\n--- [Model ", i, "/3: ", label, "] ---\n", sep = "")

  # Run base case (re-sources 06 + 07 internally)
  tryCatch({
    bc <- run_basecase_for_structure(ms, verbose = TRUE)
    br <- bc$base_results
    io <- bc$icer_obj

    # Calculate NMB
    nmb <- data.frame(
      Strategy = br$Strategy,
      Cost = br$Cost,
      Effect = br$Effect,
      NMB = br$Effect * wtp_val - br$Cost
    )
    nmb <- nmb[order(-nmb$NMB), ]

    # Load PSA cache
    psa_obj <- load_psa_cache_for_structure(ms)
    psa_summary <- NULL
    if (!is.null(psa_obj)) {
      psa_summary <- compute_psa_summary(psa_obj)
      psa_objects[[label]] <- psa_obj
    } else {
      cat("  WARNING: PSA cache not found for", label, "- skipping PSA\n")
    }

    models_results[[label]] <- list(
      base_results = br,
      icer_obj = as.data.frame(io),
      nmb_at_wtp = nmb,
      psa_summary = psa_summary
    )
    cat("  Done:", label, "\n")

  }, error = function(e) {
    cat("  ERROR running", label, ":", conditionMessage(e), "\n")
    cat("  Skipping this model structure.\n")
  })
}

cat("\n=== Multi-Model Analysis Complete ===\n")
cat("Models completed:", paste(names(models_results), collapse = ", "), "\n")

if (length(models_results) == 0) {
  stop("No models completed successfully. Cannot save snapshot.")
}

# ============================================================================
# Step 4: Collect metadata and build snapshot
# ============================================================================
cat("\nCollecting metadata...\n")
metadata <- collect_metadata(issue_number)

# Use first available model as the primary (legacy top-level fields)
primary_label <- names(models_results)[1]
primary <- models_results[[primary_label]]

cat("Creating snapshot object...\n")
snapshot <- list(
  metadata = metadata,
  # Legacy top-level fields (backward compatibility)
  base_results = primary$base_results,
  icer_obj = primary$icer_obj,
  nmb_at_wtp = primary$nmb_at_wtp,
  psa_summary = primary$psa_summary,
  # Multi-model results
  models = models_results
)

# ============================================================================
# Step 5: Save files
# ============================================================================
commit <- get_git_commit()
snapshot_filename <- paste0("snapshot_", issue_number, "_",
                            snapshot_status, "_", commit, ".rds")
psa_filename <- paste0("psa_", issue_number, "_",
                        snapshot_status, "_", commit, ".rds")
snapshots_dir <- here::here("data", "output", "snapshots")

cat("\nSaving snapshot file...\n")
snapshot_path <- file.path(snapshots_dir, snapshot_filename)
saveRDS(snapshot, snapshot_path)
cat("  Saved:", snapshot_path, "\n")

cat("Saving PSA objects...\n")
psa_path <- file.path(snapshots_dir, psa_filename)
if (length(psa_objects) > 0) {
  saveRDS(psa_objects, psa_path)
  cat("  Saved:", psa_path, "(", length(psa_objects), "models)\n")
} else {
  cat("  WARNING: No PSA objects available - PSA file not saved\n")
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
for (label in names(models_results)) {
  mr <- models_results[[label]]
  cat("  - ", label, ": ", nrow(mr$base_results), " strategies",
      if (!is.null(mr$psa_summary)) " + PSA" else "",
      "\n", sep = "")
}
cat("  - Metadata (git commit, R version, package versions, parameters)\n")

if (snapshot_status == "baseline") {
  cat("\nNext steps:\n")
  cat("  1. Fix the bug and commit your changes\n")
  cat("  2. Run: Rscript scripts/R/analysis/15_save_snapshot.R",
      issue_number, "fixed\n")
  cat("  3. Generate report: quarto render",
      "scripts/QMD/technical_docs/bug_fix_impact.qmd\n")
} else {
  cat("\nTo generate the bug fix impact report, run:\n")
  cat("  quarto render scripts/QMD/technical_docs/bug_fix_impact.qmd\n")
}

cat("\n=== Done ===\n\n")
