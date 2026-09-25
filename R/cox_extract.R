# ===============================================================================
# SHARED COX MODEL FITTING AND EXTRACTION HELPERS
# ===============================================================================

#' Iteration limit of fit_firth_cox()
FIRTH_COX_MAXIT <- 500

#' Fit a Firth-corrected Cox model with profile-likelihood confidence intervals
#'
#' Warnings that coxphf raises during the fit (numerical problems, iteration
#' limit reached) still propagate, and are also kept in the "fit_warnings"
#' attribute of the returned fit for coxphf_convergence() (issue #176). The
#' iteration limit is kept in the "maxit" attribute.
fit_firth_cox <- function(formula, data, label = "model") {
  fit_warnings <- character()
  fit <- tryCatch(
    withCallingHandlers(
      coxphf::coxphf(
        formula = formula,
        data = data,
        firth = TRUE,
        pl = TRUE,
        maxit = FIRTH_COX_MAXIT,
        maxstep = 0.1
      ),
      warning = function(w) fit_warnings <<- c(fit_warnings, conditionMessage(w))
    ),
    error = function(e) {
      warning(paste("coxphf failed for", label, ":", e$message))
      NULL
    }
  )
  if (!is.null(fit)) {
    attr(fit, "fit_warnings") <- fit_warnings
    attr(fit, "maxit") <- FIRTH_COX_MAXIT
  }
  fit
}

#' Convergence diagnostics of a coxphf fit (issue #176)
#'
#' coxphf returns the Newton-Raphson iterations of the penalised-likelihood
#' estimation (`iter`) and, with profile-likelihood inference, the iterations
#' of each confidence limit and penalised likelihood ratio test (`iter.ci`,
#' one row per coefficient; -9 marks a numerical error). A limit or p-value
#' whose iteration count reaches `maxit` is set to NA by coxphf.
#'
#' @param model A coxphf fit, or NULL for a failed fit.
#' @param maxit Iteration limit; defaults to the "maxit" attribute set by
#'   fit_firth_cox(), then to the `maxit` argument of the fit's call, then to
#'   the coxphf default.
#' @return One-row data frame: N, events, estimation iterations, largest
#'   profile-likelihood iteration count, iteration limit, number of warnings
#'   recorded during the fit and a status ("Converged", "Reached maxit",
#'   "Numerical warning" or "Fit failed").
coxphf_convergence <- function(model, maxit = NULL) {
  if (is.null(model)) {
    return(data.frame(N = NA_integer_, Events = NA_integer_, Iterations = NA_real_,
                      Max_PL_iterations = NA_real_, Maxit = NA_real_,
                      Warnings = NA_integer_, Status = "Fit failed",
                      stringsAsFactors = FALSE))
  }
  if (is.null(maxit)) maxit <- attr(model, "maxit")
  if (is.null(maxit)) maxit <- model$call$maxit
  if (is.null(maxit)) maxit <- formals(coxphf::coxphf)$maxit
  if (!is.numeric(maxit) || length(maxit) != 1L) {
    stop("coxphf_convergence(): cannot determine the iteration limit; pass maxit.")
  }
  iter_ci <- model$iter.ci
  reached <- model$iter >= maxit || (!is.null(iter_ci) && any(iter_ci >= maxit))
  pl_error <- !is.null(iter_ci) && any(iter_ci == -9)
  warns <- attr(model, "fit_warnings")
  if (is.null(warns)) warns <- character()
  status <- if (reached) "Reached maxit" else
    if (pl_error || length(warns) > 0) "Numerical warning" else "Converged"
  y <- model$y
  data.frame(N = nrow(y), Events = as.integer(sum(y[, "status"])),
             Iterations = model$iter,
             Max_PL_iterations = if (is.null(iter_ci)) NA_real_ else max(iter_ci),
             Maxit = maxit, Warnings = length(warns), Status = status,
             stringsAsFactors = FALSE)
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
