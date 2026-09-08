# ===============================================================================
# REFERENCE VALUES FOR POPULATION-REPRESENTATIVE PREDICTIONS
# ===============================================================================
# This file outputs the reference values used for estimated marginal means
# in population-level survival predictions
# ===============================================================================

if (exists("ref_age") && exists("predictions")) {
  
  cat("\n=== REFERENCE VALUES FOR POPULATION PREDICTIONS ===\n")
  cat("These values represent the 'typical patient' used for predictions:\n\n")
  
  # Continuous variables (means)
  cat("Continuous Variables:\n")
  cat("- Age (mean):", round(ref_age, 2), "years\n")
  
  # Categorical variables (modal categories)
  cat("\nCategorical Variables (modal categories):\n")
  cat("- Sex:", ref_sex, "\n")
  cat("- CRP:", ref_crp, "\n") 
  cat("- TLR:", ref_tlr, "\n")
  cat("- TMB_BRAF:", ref_tmb_braf, "\n")
  
  # Treatment assignments by strategy
  cat("\nTreatment Assignments:\n")
  cat("- Control treatment:", ctrl_rx, "\n")
  cat("- Experimental treatment:", exp_rx, "\n")
  
  # Data completeness
  cat("\nData Completeness:\n")
  cat("- Complete cases used:", nrow(data_complete), "observations\n")
  cat("- Total original observations:", nrow(data), "\n")
  cat("- Proportion complete:", round(nrow(data_complete)/nrow(data) * 100, 1), "%\n")
  
  # Biomarker prevalences (from strategies_df)
  cat("\nBiomarker Prevalences (for population weighting):\n")
  for (strategy in strategies_df$id[strategies_df$id != "control"]) {
    prevalence <- strategies_df$prevalence[strategies_df$id == strategy]
    cat("- ", toupper(strategy), ":", round(prevalence * 100, 1), "%\n")
  }
  
  # Model selection results
  cat("\nSelected Best-Fitting Models:\n")
  if (exists("models") && !is.null(models$best_fit$os)) {
    cat("- OS model:", models$best_fit$os_structure, "-", models$best_fit$os_distribution, 
        "(AIC =", round(models$best_fit$os_aic, 2), ")\n")
  }
  if (exists("models") && !is.null(models$best_fit$pfs)) {
    cat("- PFS model:", models$best_fit$pfs_structure, "-", models$best_fit$pfs_distribution,
        "(AIC =", round(models$best_fit$pfs_aic, 2), ")\n")
  }
  
  cat("\n=== END REFERENCE VALUES ===\n\n")
  
} else {
  cat("Reference values not available. Run main analysis first.\n")
}