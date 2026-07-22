# ===============================================================================
# SHARED COX MODEL FITTING AND EXTRACTION HELPERS
# ===============================================================================

#' Fit a Firth-corrected Cox model with profile-likelihood confidence intervals
fit_firth_cox <- function(formula, data, label = "model") {
  tryCatch(
    coxphf::coxphf(
      formula = formula,
      data = data,
      firth = TRUE,
      pl = TRUE,
      maxit = 500,
      maxstep = 0.1
    ),
    error = function(e) {
      warning(paste("coxphf failed for", label, ":", e$message))
      NULL
    }
  )
}

#' Extract a full result table from a coxphf model
#'
#' `ci.lower` and `ci.upper` from coxphf are already on the hazard-ratio scale.
extract_coxphf_results <- function(model, reorder_interactions = TRUE) {
  if (is.null(model)) return(NULL)

  results <- data.frame(
    Term = names(model$coefficients),
    HR = exp(model$coefficients),
    CI_Lower = model$ci.lower,
    CI_Upper = model$ci.upper,
    P_Value = model$prob,
    stringsAsFactors = FALSE
  )

  if (isTRUE(reorder_interactions)) {
    is_interaction <- grepl(":", results$Term)
    results <- rbind(results[is_interaction, ], results[!is_interaction, ])
  }

  rownames(results) <- NULL
  results
}

#' Extract the first matching interaction HR and profile-likelihood interval
get_interaction_hr <- function(results_df, pattern) {
  empty <- list(HR = NA_real_, CI_Lower = NA_real_, CI_Upper = NA_real_)
  if (is.null(results_df)) return(empty)

  idx <- grep(pattern, results_df$Term)
  if (length(idx) == 0) return(empty)

  list(
    HR = results_df$HR[idx[1]],
    CI_Lower = results_df$CI_Lower[idx[1]],
    CI_Upper = results_df$CI_Upper[idx[1]]
  )
}

#' Extract the first matching interaction HR and standard error from ridge Cox
extract_ridge_interaction <- function(model, term_pattern) {
  if (is.null(model)) return(list(HR = NA_real_, SE = NA_real_))

  coef_table <- summary(model)$coefficients
  idx <- grep(term_pattern, rownames(coef_table))
  if (length(idx) == 0) return(list(HR = NA_real_, SE = NA_real_))

  list(
    HR = exp(coef_table[idx[1], "coef"]),
    SE = coef_table[idx[1], "se(coef)"]
  )
}
