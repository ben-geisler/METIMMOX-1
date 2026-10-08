# ===============================================================================
# SUPPLEMENTARY TABLE S3 (clinical effectiveness): CONDITIONAL INDEPENDENCE
# ===============================================================================
# The testable conditional independencies implied by the canonical DAG, then
# those that hold by construction, then the omitted interaction edges into TLR
# (effect modification of T -> TLR). The statistics come from run_dag_ci_tests()
# and run_dag_effect_modification_tests() (dag_association_tests.R); this file
# only orders, labels and annotates them (issue #192). Needs assoc_tests.R,
# supp_table_labels.R and dag_association_tests.R.

SUPP_S3_OMITTED_SECTION <- "Omitted edges: effect modification of T -> TLR"

SUPP_S3_TEST_FAMILIES <- c(
  "Wilcoxon rank-sum" = "Wilcoxon rank-sum test",
  "Fisher's exact" = "Fisher's exact test",
  "Mantel-Haenszel" = "Mantel-Haenszel test",
  "Conditional permutation" = "conditional permutation test (coin, 9,999 resamples)",
  "Firth logistic" = "Firth logistic regression"
)

#' Build Supplementary Table S3 (conditional independencies) and its footnotes
#'
#' @param ci_results Output of `run_dag_ci_tests()`.
#' @param em_results Output of `run_dag_effect_modification_tests()`.
#' @param dd Output of `prepare_dag_data()`.
#' @return List with `raw` (ASCII item labels, unformatted estimates, intervals,
#'   p-values and the check flags), `table` (Section, Item, Test, N,
#'   Effect (95% CI), P, Conclusion) and `footnotes` (character vector).
build_supp_table_s3 <- function(ci_results, em_results, dd) {
  cov <- attr(ci_results, "coverage")
  df_ci <- rbind(ci_results[!ci_results$Definitional, ],
                 ci_results[ci_results$Definitional, ])
  df_ci$Section <- "Conditional independence"
  df_em <- em_results
  df_em$Section <- SUPP_S3_OMITTED_SECTION
  df_em$Randomization_Check <- FALSE
  df_em$Stratum_Balance <- FALSE
  df_em$Definitional <- FALSE

  v_cols <- c("Section", "Item", "Test", "N", "Estimate", "CI_Lower", "CI_Upper",
              "p_value", "Effect", "Conclusion", "Randomization_Check",
              "Stratum_Balance", "Definitional")
  df_raw <- as.data.frame(rbind(df_ci[, v_cols], df_em[, v_cols]),
                          stringsAsFactors = FALSE)
  names(df_raw)[names(df_raw) == "Randomization_Check"] <- "Arm_balance_marginal"
  names(df_raw)[names(df_raw) == "Stratum_Balance"] <- "Arm_balance_in_stratum"

  v_conclusion <- df_raw$Conclusion
  v_is_ci <- df_raw$Section != SUPP_S3_OMITTED_SECTION
  v_conclusion[v_is_ci] <- supp_ci_conclusion(df_raw$Conclusion[v_is_ci])
  df_table <- data.frame(
    Section = supp_label(df_raw$Section),
    Item = supp_label(df_raw$Item), Test = supp_label(df_raw$Test),
    N = ifelse(is.na(df_raw$N), "--", as.character(df_raw$N)),
    `Effect (95% CI)` = supp_label(df_raw$Effect),
    P = fmt_p(df_raw$p_value), Conclusion = v_conclusion,
    stringsAsFactors = FALSE, check.names = FALSE
  )
  df_raw$Conclusion <- v_conclusion  # the raw file carries the wording of the table

  supp_check_test_families(df_raw$Test[!is.na(df_raw$p_value)],
                           SUPP_S3_TEST_FAMILIES, "Table S3")

  n_ci <- nrow(cov)
  n_def <- sum(cov$definitional)
  v_balance <- supp_label(cov$statement[cov$arm_balance])
  v_stratum <- supp_label(cov$statement[cov$stratum_arm_balance])
  n_main <- nrow(dd$data_main)
  n_tlr <- nrow(dd$data_tlr)

  # Week-4 CRP arm imbalance in the two cohorts it is reported for.
  df_crp_all <- ci_results[ci_results$Item == "CRP _||_ T", ]
  df_crp_main <- run_fisher_test(dd$data_main, "crp", "T_num", "CRP _||_ T",
                                 "Fisher's exact")
  if (nrow(df_crp_all) != 1L) stop("Table S3: CRP _||_ T is missing from the CI tests.")

  v_footnotes <- c(
    "Data from the METIMMOX trial (NCT03388190).",
    "CRP = C-reactive protein (week 4, before the first nivolumab dose); TLR = target lesion reduction, a reduction of at least 10% in the sum of target-lesion diameters at the first on-treatment CT; TMB = tumor mutational burden; T = treatment; T x CRP and T x TMB = treatment-by-biomarker interaction.",
    sprintf("Tests: %s. Firth-corrected models use profile-likelihood 95%% confidence intervals and p-values. Effects are odds ratios for logistic, Fisher and Mantel-Haenszel tests. N is the number of patients in each test. No multiplicity correction applied.",
            paste(SUPP_S3_TEST_FAMILIES, collapse = "; ")),
    sprintf("The DAG implies %d conditional independencies. %d hold by construction and are not tested: one variable is an interaction node whose parents (T and the biomarker) are both in the conditioning set, so the node is constant within every stratum. The other %d are tested.",
            n_ci, n_def, n_ci - n_def),
    sprintf("Arm-balance checks: %s (marginal) and %s (within biomarker-positive patients). Within biomarker-positive patients an interaction node equals T (and is 0 otherwise), so these statements compare the arms within that stratum.",
            paste(v_balance, collapse = ", "), paste(v_stratum, collapse = ", ")),
    sprintf("CRP is the week-4 value (before the first nivolumab dose). CRP %s T therefore tests arm balance at the immunotherapy decision point, not randomization balance, and CRP %s T x TMB | {TMB_BRAF} is the same week-4 CRP comparison within TMB/BRAF-positive patients. CRP %s T is tested in all %d patients with a week-4 CRP value (OR %.2f, p = %s, as in the table); in the %d-patient complete-case cohort of the analysis the same test gives OR %.2f, p = %s.",
            SUPP_PERP, SUPP_PERP, SUPP_PERP, df_crp_all$N, df_crp_all$Estimate,
            fmt_p(df_crp_all$p_value), df_crp_main$N, df_crp_main$Estimate,
            fmt_p(df_crp_main$p_value)),
    sprintf("Omitted edges: T x CRP -> TLR and T x TMB -> TLR are tested as effect modification of T -> TLR (Firth logistic interaction) and are not counted among the conditional independencies. Tests into TLR use every patient with an observed TLR (n = %d), or fewer when the biomarker is missing (the %d TLR-classified patients with a TMB/BRAF result for the tests that include TMB/BRAF).",
            n_tlr, nrow(dd$data_tlr_landmark))
  )

  list(raw = df_raw, table = df_table, footnotes = v_footnotes)
}
