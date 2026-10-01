# ===============================================================================
# SHARED DAG ASSOCIATION-TEST HELPERS
# ===============================================================================

complete_data <- function(df, vars) {
  df[stats::complete.cases(df[, vars, drop = FALSE]), , drop = FALSE]
}

fmt_num <- function(x, digits = 2) {
  ifelse(is.na(x), "--", formatC(x, digits = digits, format = "f"))
}

fmt_p <- function(x) {
  ifelse(
    is.na(x), "--",
    ifelse(x < 0.001, "<0.001", formatC(x, digits = 3, format = "f"))
  )
}

fmt_hr <- function(hr, lower, upper) {
  if (any(is.na(c(hr, lower, upper)))) return("--")
  sprintf("%.2f (95%% CI %.2f to %.2f)", hr, lower, upper)
}

fmt_or <- function(or, lower, upper) {
  if (any(is.na(c(or, lower, upper)))) return("--")
  sprintf("%.2f (95%% CI %.2f to %.2f)", or, lower, upper)
}

run_wilcox_test <- function(df, group_var, value_var, item, test_label,
                            group_labels) {
  df_complete <- complete_data(df, c(group_var, value_var))
  df_complete[[group_var]] <- factor(df_complete[[group_var]])
  wt <- stats::wilcox.test(df_complete[[value_var]] ~ df_complete[[group_var]], exact = FALSE)
  v_medians <- tapply(df_complete[[value_var]], df_complete[[group_var]], median, na.rm = TRUE)

  tibble::tibble(
    Item = item,
    Test = test_label,
    N = nrow(df_complete),
    Effect = paste0(
      "Median ", value_var, ": ",
      group_labels[1], "=", fmt_num(unname(v_medians[1]), 1),
      ", ", group_labels[2], "=", fmt_num(unname(v_medians[2]), 1)
    ),
    p_value = unname(wt$p.value)
  )
}

run_fisher_test <- function(df, x_var, y_var, item, test_label) {
  df_complete <- complete_data(df, c(x_var, y_var))
  ft <- stats::fisher.test(table(df_complete[[x_var]], df_complete[[y_var]]))

  tibble::tibble(
    Item = item,
    Test = test_label,
    N = nrow(df_complete),
    Effect = fmt_or(unname(ft$estimate), ft$conf.int[1], ft$conf.int[2]),
    p_value = unname(ft$p.value)
  )
}

run_mh_test <- function(df, x_var, y_var, strata_var, item, test_label) {
  df_complete <- complete_data(df, c(x_var, y_var, strata_var))
  mh <- stats::mantelhaen.test(
    table(df_complete[[x_var]], df_complete[[y_var]], df_complete[[strata_var]])
  )

  tibble::tibble(
    Item = item,
    Test = test_label,
    N = nrow(df_complete),
    Effect = fmt_or(unname(mh$estimate), mh$conf.int[1], mh$conf.int[2]),
    p_value = unname(mh$p.value)
  )
}

run_firth_cox_test <- function(df, formula, term, item, test_label) {
  df_complete <- complete_data(df, all.vars(formula))
  fit <- coxphf::coxphf(
    formula = formula, data = df_complete, firth = TRUE, pl = TRUE,
    maxit = 100, maxstep = 0.1
  )

  tibble::tibble(
    Item = item,
    Test = test_label,
    N = nrow(df_complete),
    Effect = fmt_hr(
      exp(unname(fit$coefficients[term])),
      unname(fit$ci.lower[term]),
      unname(fit$ci.upper[term])
    ),
    p_value = unname(fit$prob[term])
  )
}

run_firth_logistic_test <- function(df, formula, term, item, test_label) {
  df_complete <- complete_data(df, all.vars(formula))
  fit <- logistf::logistf(
    formula = formula, data = df_complete, firth = TRUE, pl = TRUE
  )

  tibble::tibble(
    Item = item,
    Test = test_label,
    N = nrow(df_complete),
    Effect = fmt_or(
      exp(unname(fit$coefficients[term])),
      exp(unname(fit$ci.lower[term])),
      exp(unname(fit$ci.upper[term]))
    ),
    p_value = unname(fit$prob[term])
  )
}

run_permutation_test <- function(df, formula, item, test_label, seed = 123L) {
  df_complete <- complete_data(df, all.vars(formula))
  set.seed(seed)
  perm <- coin::independence_test(
    formula,
    data = df_complete,
    distribution = coin::approximate(nresample = 9999)
  )

  tibble::tibble(
    Item = item,
    Test = test_label,
    N = nrow(df_complete),
    Effect = paste0(
      "Standardized Z = ",
      fmt_num(as.numeric(coin::statistic(perm, type = "standardized")), 2)
    ),
    p_value = as.numeric(coin::pvalue(perm))
  )
}

edge_support_label <- function(p_value) {
  ifelse(
    is.na(p_value), "Not estimable",
    ifelse(
      p_value < 0.05,
      "Association detected (p < 0.05)",
      "No association detected (p >= 0.05)"
    )
  )
}

ci_support_label <- function(p_value, randomization_check) {
  ifelse(
    is.na(p_value), "Not estimable",
    ifelse(
      p_value < 0.05,
      ifelse(
        randomization_check,
        "Arm imbalance in realised sample",
        "CI contradicted by data"
      ),
      ifelse(
        randomization_check,
        "Compatible with arm balance",
        "Compatible with DAG-implied CI"
      )
    )
  )
}
