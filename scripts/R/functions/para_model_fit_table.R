# ===============================================================================
# PARAMETRIC MODEL FIT COMPARISON FUNCTIONS
# ===============================================================================
# This file contains functions for creating model comparison tables.
# Display logic has been moved to para_models.Rmd
# ===============================================================================

# Function to create model comparison tables
create_model_comparison_tables <- function(models) {
  
  if (is.null(models) || is.null(models$full)) {
    return(NULL)
  }
  
  # Create combined table for OS models
  os_comparison <- rbind(
    data.frame(Structure = "Full", models$full$os_ic),
    data.frame(Structure = "No Age/Sex", models$no_age_sex$os_ic)
  )
  
  # Create combined table for PFS models  
  pfs_comparison <- rbind(
    data.frame(Structure = "Full", models$full$pfs_ic),
    data.frame(Structure = "No Age/Sex", models$no_age_sex$pfs_ic)
  )
  
  return(list(
    os_comparison = os_comparison,
    pfs_comparison = pfs_comparison
  ))
}

# Function to calculate model fitting summary statistics
calculate_model_summary <- function(models, distributions_to_test) {
  
  if (is.null(models)) {
    return(NULL)
  }
  
  # Count successful fits
  os_successful <- sum(sapply(models$full$os, function(x) inherits(x, "flexsurvreg"))) +
    sum(sapply(models$no_age_sex$os, function(x) inherits(x, "flexsurvreg")))
  
  pfs_successful <- sum(sapply(models$full$pfs, function(x) inherits(x, "flexsurvreg"))) +
    sum(sapply(models$no_age_sex$pfs, function(x) inherits(x, "flexsurvreg")))
  
  total_attempted <- length(distributions_to_test) * 2 * 2  # distributions × structures × outcomes
  
  return(list(
    total_attempted = total_attempted,
    os_successful = os_successful,
    pfs_successful = pfs_successful,
    success_rate = (os_successful + pfs_successful) / total_attempted * 100
  ))
}