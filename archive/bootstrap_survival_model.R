# Function to bootstrap a survival model using normboot.flexsurvreg
bootstrap_survival_model <- function(formula, data, dist = "weibull", n_boot = n_samples) {
  # Fit the original model
  original_model <- flexsurvreg(formula, data = data, dist = dist)
  
  # Use normboot.flexsurvreg to generate parameter samples
  # Include the original data as newdata to provide covariate information
  tryCatch({
    boot_results <- normboot.flexsurvreg(original_model, B = n_boot, 
                                         transform = TRUE, newdata = data)
  }, error = function(e) {
    cat("Bootstrap failed for model:", conditionMessage(e), "\n")
    # Fallback to parametric uncertainty
    return(NULL)
  })
  
  # Extract the matrix of bootstrapped coefficients
  coefficients_matrix <- t(boot_results)
  
  return(list(
    original = original_model,
    bootstrap = boot_results,
    coefficients_matrix = coefficients_matrix
  ))
}