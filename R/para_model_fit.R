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

# fit_all_direct function by Frederick / fits multiple parametric distributions to survival data
# NB: I had to turn off strip_bkgr_dat for compatibility with later functions
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
                           strip_bkgr_dat = FALSE # to delete data from fitted objects
){
  
  # Load packages
  #require(survival)
  #require(flexsurv)
  
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
  v_main_effects <- attr(formula_terms, "term.labels")
  # Remove interaction terms (those containing ":")
  v_main_effects <- v_main_effects[!grepl(":", v_main_effects)]
  
  if(length(v_main_effects) > 0) {
    # Create formula with main effects only
    km_formula <- as.formula(paste(deparse(fit_formula[[2]]), "~", paste(v_main_effects, collapse = " + ")))
    km_std <- eval(as.call(list(quote(survival::survfit),
                                formula = km_formula,
                                data = quote(fit_data))))
  } else {
    # If no main effects, use the intercept-only model
    km_std <- km_all
  }
  
  # Fit all standard models (these can handle interactions). Warnings raised by
  # a candidate fit (e.g. a non-positive-definite Hessian) are recorded per
  # distribution in `fit_warnings` instead of propagating; the caller decides
  # whether they matter (script 04 re-raises them for the selected families).
  fit_warnings <- setNames(vector("list", length(fit_dists)), fit_dists)
  fits <- lapply(fit_dists,
                 function(x){
                   withCallingHandlers(
                     tryCatch(expr = eval(as.call(list(quote(flexsurvreg),
                                                       formula = fit_formula,
                                                       data = quote(fit_data),
                                                       dist = x))),
                              error = function(e) paste("Error message:", e$message)),
                     warning = function(w) {
                       fit_warnings[[x]] <<- c(fit_warnings[[x]], conditionMessage(w))
                       invokeRestart("muffleWarning")
                     })
                 })
  
  # Strip background data from fitted models
  if(strip_bkgr_dat) {
    for (i in seq_along(fit_dists)) {
      tryCatch(fits[[i]][["data"]] <- NULL,
               error = function(e) NULL)
    }
  }
  
  # Create result list
  res <- list(data_name = data_name,
              km_all = km_all,
              km_std = km_std,
              fitted_models = fits,
              fit_warnings = fit_warnings)
}

# Extract AIC and BIC from fitted models
# Returns a data frame with model fit statistics
extract_ic_single <- function(fit_results) {
  df_ic <- data.frame(
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
    
    df_ic <- rbind(df_ic, data.frame(
      Distribution = dist_name,
      AIC = aic_val,
      BIC = bic_val,
      stringsAsFactors = FALSE
    ))
  }
  
  return(df_ic)
}

# Find best fitting model based on information criterion / returns list with best distribution and IC
find_best_model <- function(ic_data, criterion = "AIC") {
  df_valid <- ic_data[!is.na(ic_data[[criterion]]), ]
  if (nrow(df_valid) == 0) {
    return(NULL)
  }
  
  best_idx <- which.min(df_valid[[criterion]])
  return(list(
    distribution = df_valid$Distribution[best_idx],
    criterion_value = df_valid[[criterion]][best_idx],
    criterion = criterion
  ))
}

# Get modal (most frequent) category for categorical variables
get_modal_category <- function(x) {
  tbl <- table(x)
  names(tbl)[which.max(tbl)]
}

# Select OS and PFS distributions jointly subject to OS >= PFS.
#
# `ordering_check` is deliberately injected by the analysis script so this
# fitting helper remains independent of the population-prediction implementation.
# It must accept two endpoint prediction objects (OS, PFS) and return a list
# containing `ordered`, `n_violations`, and `max_pfs_minus_os`.
find_best_ordered_model_pair <- function(os_candidates, pfs_candidates,
                                         ordering_check,
                                         criterion = "AIC") {
  v_required_columns <- c("formula_set", "distribution", criterion,
                        "ordering_data")
  if (!all(v_required_columns %in% names(os_candidates)) ||
      !all(v_required_columns %in% names(pfs_candidates))) {
    stop("Candidate tables must contain: ",
         paste(v_required_columns, collapse = ", "))
  }

  rows <- vector("list", nrow(os_candidates) * nrow(pfs_candidates))
  row_idx <- 0L

  for (os_idx in seq_len(nrow(os_candidates))) {
    for (pfs_idx in seq_len(nrow(pfs_candidates))) {
      row_idx <- row_idx + 1L
      check <- tryCatch(
        ordering_check(
          os_candidates$ordering_data[[os_idx]],
          pfs_candidates$ordering_data[[pfs_idx]]
        ),
        error = function(e) list(
          ordered = FALSE,
          n_violations = NA_integer_,
          max_pfs_minus_os = Inf,
          error = conditionMessage(e)
        )
      )

      rows[[row_idx]] <- data.frame(
        os_formula_set = os_candidates$formula_set[os_idx],
        os_distribution = os_candidates$distribution[os_idx],
        os_ic = os_candidates[[criterion]][os_idx],
        pfs_formula_set = pfs_candidates$formula_set[pfs_idx],
        pfs_distribution = pfs_candidates$distribution[pfs_idx],
        pfs_ic = pfs_candidates[[criterion]][pfs_idx],
        combined_ic = os_candidates[[criterion]][os_idx] +
          pfs_candidates[[criterion]][pfs_idx],
        ordered = isTRUE(check$ordered),
        n_violations = check$n_violations,
        max_pfs_minus_os = check$max_pfs_minus_os,
        check_error = if (is.null(check$error)) NA_character_ else check$error,
        stringsAsFactors = FALSE
      )
    }
  }

  df_pairs <- do.call(rbind, rows)
  df_pairs <- df_pairs[order(df_pairs$combined_ic), ]
  rownames(df_pairs) <- NULL
  df_feasible <- df_pairs[df_pairs$ordered, ]

  if (nrow(df_feasible) == 0) {
    return(list(
      selected = NULL,
      pairs = df_pairs,
      unconstrained_combined_ic = min(df_pairs$combined_ic),
      ic_penalty = Inf,
      criterion = criterion
    ))
  }

  df_selected <- df_feasible[which.min(df_feasible$combined_ic), , drop = FALSE]
  list(
    selected = df_selected,
    pairs = df_pairs,
    unconstrained_combined_ic = min(df_pairs$combined_ic),
    ic_penalty = df_selected$combined_ic[[1]] - min(df_pairs$combined_ic),
    criterion = criterion
  )
}
