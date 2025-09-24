# ===============================================================================
# PARAMETRIC MODEL FITTING FUNCTIONS FOR SURVIVAL ANALYSIS
# ===============================================================================
# This file contains helper functions for fitting and evaluating 
# parametric survival models across different model structures and distributions
# ===============================================================================

# Load required packages
if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(survival, flexsurv, dplyr)

# fit_all_direct function by Frederick
# Fits multiple parametric distributions to survival data
fit_all_direct <- function(fit_data,
                           fit_formula, # the survival formula
                           fit_dists = c("exponential",
                                         "gamma",
                                         "genf",
                                         "gengamma",
                                         "gompertz",
                                         "weibull",
                                         "weibullph",
                                         "llogis",
                                         "lognormal"),
                           strip_bkgr_dat = TRUE # to delete data from fitted objects
){
  
  # Load packages
  require(survival)
  require(flexsurv)
  
  # Define names to automatically create named list of fits
  names(fit_dists) <- fit_dists
  
  # Extract data object name
  data_name <- deparse(substitute(fit_data))
  
  # Fit KM ignoring any subgroups
  km_all <- eval(as.call(list(quote(survival::survfit),
                              formula = update(fit_formula, ~ 1),
                              data = quote(fit_data))))
  
  # Fit KM with strata - but remove interaction terms first
  # Extract main effects only for survfit
  formula_terms <- terms(fit_formula)
  main_effects <- attr(formula_terms, "term.labels")
  # Remove interaction terms (those containing ":")
  main_effects <- main_effects[!grepl(":", main_effects)]
  
  if(length(main_effects) > 0) {
    # Create formula with main effects only
    km_formula <- as.formula(paste(deparse(fit_formula[[2]]), "~", paste(main_effects, collapse = " + ")))
    km_std <- eval(as.call(list(quote(survival::survfit),
                                formula = km_formula,
                                data = quote(fit_data))))
  } else {
    # If no main effects, use the intercept-only model
    km_std <- km_all
  }
  
  # Fit all standard models (these can handle interactions)
  fits <- lapply(fit_dists,
                 function(x){
                   tryCatch(expr = eval(as.call(list(quote(flexsurvreg),
                                                     formula = fit_formula,
                                                     data = quote(fit_data),
                                                     dist = x))),
                            error = function(e) paste("Error message:", e$message))
                 })
  
  # Strip background data from fitted models
  if(strip_bkgr_dat) {
    for (i in 1:length(fit_dists)) {
      tryCatch(fits[[i]][["data"]] <- NULL,
               error = function(e) NULL)
    }
  }
  
  # Create result list
  res <- list(data_name = data_name,
              km_all = km_all,
              km_std = km_std,
              fitted_models = fits)
}

# Extract AIC and BIC from fitted models
# Returns a data frame with model fit statistics
extract_ic_single <- function(fit_results) {
  ic_data <- data.frame(
    Distribution = character(),
    AIC = numeric(),
    BIC = numeric(),
    stringsAsFactors = FALSE
  )
  
  fitted_models <- fit_results$fitted_models
  
  for (dist_name in names(fitted_models)) {
    fit_obj <- fitted_models[[dist_name]]
    
    aic_val <- NA
    bic_val <- NA
    
    if (inherits(fit_obj, "flexsurvreg")) {
      # Try AIC 
      aic_result <- try(AIC(fit_obj), silent = TRUE)
      if (!inherits(aic_result, "try-error")) {
        aic_val <- aic_result
      }
      
      # Calculate BIC manually: BIC = -2*loglik + k*log(n)
      tryCatch({
        loglik <- fit_obj$loglik
        k <- length(fit_obj$coefficients)
        n <- fit_obj$N
        
        if (!is.null(loglik) && !is.null(k) && !is.null(n)) {
          bic_val <- -2 * loglik + k * log(n)
        }
      }, error = function(e) NULL)
    }
    
    ic_data <- rbind(ic_data, data.frame(
      Distribution = dist_name,
      AIC = aic_val,
      BIC = bic_val,
      stringsAsFactors = FALSE
    ))
  }
  
  return(ic_data)
}

# Find best fitting model based on information criterion
# Returns list with best distribution and criterion value
find_best_model <- function(ic_data, criterion = "AIC") {
  valid_data <- ic_data[!is.na(ic_data[[criterion]]), ]
  if (nrow(valid_data) == 0) {
    return(NULL)
  }
  
  best_idx <- which.min(valid_data[[criterion]])
  return(list(
    distribution = valid_data$Distribution[best_idx],
    criterion_value = valid_data[[criterion]][best_idx],
    criterion = criterion
  ))
}

# Get modal (most frequent) category for categorical variables
get_modal_category <- function(x) {
  tbl <- table(x)
  names(tbl)[which.max(tbl)]
}