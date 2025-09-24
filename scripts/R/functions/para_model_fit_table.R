# ===============================================================================
# PARAMETRIC MODEL FIT COMPARISON TABLE
# ===============================================================================
# This file creates comprehensive tables comparing all fitted models
# across structures and distributions using AIC and BIC
# ===============================================================================

if (exists("models") && !is.null(models$full)) {
  
  cat("\n=== PARAMETRIC MODEL FIT COMPARISON TABLES ===\n")
  
  # Create combined table for OS models
  cat("\nOVERALL SURVIVAL (OS) MODELS:\n")
  cat(paste(rep("=", 50), collapse = ""), "\n")
  
  os_comparison <- rbind(
    data.frame(Structure = "Full", models$full$os_ic),
    data.frame(Structure = "No Age/Sex", models$no_age_sex$os_ic)
  )
  
  # Sort by AIC for easy comparison
  os_comparison <- os_comparison[order(os_comparison$AIC, na.last = TRUE), ]
  
  # Add ranking
  valid_aic <- !is.na(os_comparison$AIC)
  os_comparison$AIC_Rank <- NA
  os_comparison$AIC_Rank[valid_aic] <- rank(os_comparison$AIC[valid_aic])
  
  # Format for display
  os_comparison$AIC <- round(os_comparison$AIC, 2)
  os_comparison$BIC <- round(os_comparison$BIC, 2)
  
  # Mark best model
  if (exists("models") && !is.null(models$best_fit$os_aic)) {
    best_idx <- which.min(abs(os_comparison$AIC - models$best_fit$os_aic))
    if (length(best_idx) > 0) {
      os_comparison$Best <- ""
      os_comparison$Best[best_idx] <- "*** SELECTED ***"
    }
  }
  
  print(os_comparison)
  
  # Create combined table for PFS models  
  cat("\nPROGRESSION-FREE SURVIVAL (PFS) MODELS:\n")
  cat(paste(rep("=", 50), collapse = ""), "\n")
  
  pfs_comparison <- rbind(
    data.frame(Structure = "Full", models$full$pfs_ic),
    data.frame(Structure = "No Age/Sex", models$no_age_sex$pfs_ic)
  )
  
  # Sort by AIC for easy comparison
  pfs_comparison <- pfs_comparison[order(pfs_comparison$AIC, na.last = TRUE), ]
  
  # Add ranking
  valid_aic <- !is.na(pfs_comparison$AIC)
  pfs_comparison$AIC_Rank <- NA
  pfs_comparison$AIC_Rank[valid_aic] <- rank(pfs_comparison$AIC[valid_aic])
  
  # Format for display
  pfs_comparison$AIC <- round(pfs_comparison$AIC, 2)
  pfs_comparison$BIC <- round(pfs_comparison$BIC, 2)
  
  # Mark best model
  if (exists("models") && !is.null(models$best_fit$pfs_aic)) {
    best_idx <- which.min(abs(pfs_comparison$AIC - models$best_fit$pfs_aic))
    if (length(best_idx) > 0) {
      pfs_comparison$Best <- ""
      pfs_comparison$Best[best_idx] <- "*** SELECTED ***"
    }
  }
  
  print(pfs_comparison)
  
  # Summary statistics
  cat("\nPARAMETRIC MODEL FITTING SUMMARY:\n")
  cat(paste(rep("=", 35), collapse = ""), "\n")
  
  # Count successful fits
  os_successful <- sum(sapply(models$full$os, function(x) inherits(x, "flexsurvreg"))) +
    sum(sapply(models$no_age_sex$os, function(x) inherits(x, "flexsurvreg")))
  
  pfs_successful <- sum(sapply(models$full$pfs, function(x) inherits(x, "flexsurvreg"))) +
    sum(sapply(models$no_age_sex$pfs, function(x) inherits(x, "flexsurvreg")))
  
  total_attempted <- length(distributions_to_test) * 2 * 2  # distributions × structures × outcomes
  
  cat("- Total models attempted:", total_attempted, "\n")
  cat("- OS models successfully fitted:", os_successful, "\n")
  cat("- PFS models successfully fitted:", pfs_successful, "\n")
  cat("- Overall success rate:", round((os_successful + pfs_successful) / total_attempted * 100, 1), "%\n")
  
  # Best model summary
  if (exists("models") && !is.null(models$best_fit$os)) {
    cat("\nSELECTED MODELS FOR PREDICTIONS:\n")
    cat("- OS: ", models$best_fit$os_structure, " structure, ", models$best_fit$os_distribution, 
        " distribution (AIC = ", round(models$best_fit$os_aic, 2), ")\n", sep = "")
    cat("- PFS: ", models$best_fit$pfs_structure, " structure, ", models$best_fit$pfs_distribution,
        " distribution (AIC = ", round(models$best_fit$pfs_aic, 2), ")\n", sep = "")
  }
  
  cat("\n=== END PARAMETRIC MODEL FIT TABLES ===\n\n")
  
  # Store tables for further use
  para_model_fit_tables <- list(
    os_comparison = os_comparison,
    pfs_comparison = pfs_comparison
  )
  
} else {
  cat("Parametric model fit tables not available. Run main analysis first.\n")
}