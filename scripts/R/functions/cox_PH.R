# Cox Proportional Hazards Functions for METIMMOX-1 Analysis

# Function to calculate median from survival probabilities
calculate_median_survival <- function(surv_probs, time_points) {
  # Find first time point where survival drops below 0.5
  median_idx <- which(surv_probs < 0.5)[1]
  if (is.na(median_idx)) {
    # If survival never drops below 0.5, return the last time point
    return(time_points[length(time_points)])
  } else {
    return(time_points[median_idx])
  }
}

# Convert weeks to months (divide by 4.345)
weeks_to_months <- function(weeks) {
  return(weeks / 4.345)
}

# Function to extract HR and CI for interaction terms
extract_hr_ci <- function(cox_summary, biomarker_name) {
  coef_table <- cox_summary$coefficients
  conf_int <- cox_summary$conf.int
  
  # Look for the interaction term with the biomarker
  # Pattern should match: biomarkerLevel:RxLevel or Rx:biomarkerLevel
  pattern1 <- paste0("^", biomarker_name, ".*:.*Rx")
  pattern2 <- paste0("^Rx.*:.*", biomarker_name)
  
  matching_rows1 <- grep(pattern1, rownames(coef_table), value = FALSE)
  matching_rows2 <- grep(pattern2, rownames(coef_table), value = FALSE)
  matching_rows <- c(matching_rows1, matching_rows2)
  
  if (length(matching_rows) > 0) {
    # Take the first matching interaction term
    idx <- matching_rows[1]
    hr <- conf_int[idx, "exp(coef)"]
    ci_lower <- conf_int[idx, "lower .95"]
    ci_upper <- conf_int[idx, "upper .95"]
    p_value <- coef_table[idx, "Pr(>|z|)"]
    
    return(list(hr = hr, ci_lower = ci_lower, ci_upper = ci_upper, p_value = p_value))
  } else {
    return(list(hr = NA, ci_lower = NA, ci_upper = NA, p_value = NA))
  }
}

# Function to extract HR with fallback logic
extract_hr_with_fallback <- function(cox_model, cox_summary, biomarker_name) {
  # First try the interaction term approach
  result <- extract_hr_ci(cox_summary, biomarker_name)
  
  # If no interaction found, try to find any term with the biomarker name
  if (is.na(result$hr)) {
    biomarker_terms <- grep(biomarker_name, names(coef(cox_model)), ignore.case = TRUE, value = TRUE)
    if (length(biomarker_terms) > 0) {
      cat("Found", biomarker_name, "terms:", biomarker_terms, "\n")
      # Use the last matching term (often the interaction)
      idx <- which(names(coef(cox_model)) == biomarker_terms[length(biomarker_terms)])
      if (length(idx) > 0) {
        result <- list(
          hr = exp(coef(cox_model)[idx]),
          ci_lower = exp(confint(cox_model)[idx, 1]),
          ci_upper = exp(confint(cox_model)[idx, 2]),
          p_value = cox_summary$coefficients[idx, "Pr(>|z|)"]
        )
      }
    }
  }
  
  return(result)
}

# Function to format HR with CI for display
format_hr_ci <- function(hr_object) {
  if (is.na(hr_object$hr)) {
    return("NA")
  } else {
    return(sprintf("%.2f (%.2f-%.2f)", hr_object$hr, hr_object$ci_lower, hr_object$ci_upper))
  }
}