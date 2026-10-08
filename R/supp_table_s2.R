# ===============================================================================
# SUPPLEMENTARY TABLE S2 (clinical effectiveness): BIVARIATE ASSOCIATION TESTS
# ===============================================================================
# The 18 testable direct edges of the canonical DAG, in the manuscript's order.
# The statistics come from run_dag_edge_tests() (dag_association_tests.R); this
# file only orders, labels and annotates them (issue #192). Needs
# assoc_tests.R, supp_table_labels.R and dag_association_tests.R.

SUPP_S2_ORDER <- c(
  "Age -> OS", "Age -> TMB_BRAF", "Sex -> TMB_BRAF",
  "CRP -> PFS", "CRP -> OS", "CRP -> TLR",
  "TMB_BRAF -> PFS", "TMB_BRAF -> OS", "TMB_BRAF -> TLR",
  "T -> PFS", "T -> OS", "T -> TLR",
  "TLR -> PFS", "PFS -> OS",
  "TxCRP -> PFS", "TxCRP -> OS", "TxTMB -> PFS", "TxTMB -> OS"
)

SUPP_S2_TEST_FAMILIES <- c(
  "Wilcoxon rank-sum" = "Wilcoxon rank-sum test",
  "Fisher's exact" = "Fisher's exact test",
  "Firth Cox, time-dependent progression" =
    "Firth Cox model with a time-dependent progression indicator (PFS -> OS)",
  "Firth Cox" = "Firth Cox regression",
  "Firth logistic" = "Firth logistic regression"
)

#' Build Supplementary Table S2 (direct edges) and its footnotes
#'
#' @param edge_results Output of `run_dag_edge_tests()`.
#' @param adj_results Output of `run_adjusted_interaction_models()`.
#' @param dd Output of `prepare_dag_data()`.
#' @param pfs_window Death window of the PFS endpoint in weeks.
#' @return List with `raw` (ASCII item labels, unformatted estimates, intervals
#'   and p-values), `table` (the displayed table: Item, Test, N,
#'   Effect (95% CI), P, Conclusion) and `footnotes` (character vector).
build_supp_table_s2 <- function(edge_results, adj_results, dd,
                                pfs_window = get0("pfs_death_window_weeks")) {
  if (!setequal(edge_results$Item, SUPP_S2_ORDER) ||
      nrow(edge_results) != length(SUPP_S2_ORDER)) {
    stop("Table S2 rows and the DAG edge tests disagree.\n  Missing from the tests: ",
         paste(setdiff(SUPP_S2_ORDER, edge_results$Item), collapse = "; "),
         "\n  Not in the table: ",
         paste(setdiff(edge_results$Item, SUPP_S2_ORDER), collapse = "; "))
  }
  df_e <- edge_results[match(SUPP_S2_ORDER, edge_results$Item), ]

  df_raw <- data.frame(
    Item = df_e$Item, Test = df_e$Test, N = df_e$N,
    Estimate = df_e$Estimate, CI_Lower = df_e$CI_Lower, CI_Upper = df_e$CI_Upper,
    p_value = df_e$p_value, Effect = df_e$Effect, Conclusion = df_e$Conclusion,
    stringsAsFactors = FALSE
  )
  df_table <- data.frame(
    Item = supp_label(df_e$Item), Test = supp_label(df_e$Test),
    N = as.character(df_e$N), `Effect (95% CI)` = df_e$Effect,
    P = fmt_p(df_e$p_value), Conclusion = df_e$Conclusion,
    stringsAsFactors = FALSE, check.names = FALSE
  )

  supp_check_test_families(df_e$Test, SUPP_S2_TEST_FAMILIES, "Table S2")

  # "0.26 (95% CI 0.06 to 1.34)" -> "0.26 (95% CI 0.06 to 1.34; PLRT p = 0.102)"
  with_p <- function(effect, p, label) {
    sub("\\)$", sprintf("; %s = %s)", label, fmt_p(p)), effect)
  }
  adj_text <- function(edge) {
    r <- adj_results[adj_results$Edge == edge, ]
    with_p(r$Effect, r$p_value, "PLRT p")
  }
  marg_text <- function(edge) {
    r <- df_e[df_e$Item == edge, ]
    with_p(r$Effect, r$p_value, "p")
  }
  n_main <- nrow(dd$data_main)
  n_tlr <- nrow(dd$data_tlr)
  n_tlr_marker <- df_e$N[df_e$Item == "TMB_BRAF -> TLR"]
  n_tmb <- df_e$N[df_e$Item == "TMB_BRAF -> PFS"]
  n_tmb_known <- df_e$N[df_e$Item == "Age -> TMB_BRAF"]
  lm_label <- tolower(dd$landmark$spec$short_label)

  v_footnotes <- c(
    "Data from the METIMMOX trial (NCT03388190).",
    "CRP = C-reactive protein (week 4, before the first nivolumab dose); TLR = target lesion reduction, a reduction of at least 10% in the sum of target-lesion diameters at the first on-treatment CT; TMB = tumor mutational burden; T = treatment; T x CRP and T x TMB = treatment-by-biomarker interaction; PFS = progression-free survival; OS = overall survival; PLRT = penalized likelihood ratio test.",
    sprintf("Tests: %s. Firth-corrected models use profile-likelihood 95%% confidence intervals and p-values. N is the number of patients in each test. No multiplicity correction applied.",
            paste(SUPP_S2_TEST_FAMILIES, collapse = "; ")),
    sprintf("Effects are hazard ratios for Cox models and odds ratios for logistic models and Fisher tests. PFS is progression or death, where a death without recorded progression is an event only within %s weeks of the last assessment (the death-window rule); later deaths are censored.",
            pfs_window),
    sprintf("These are marginal tests, not the adjusted estimates of the main text: each model contains only the parent (for the interaction edges, T, the biomarker and T x biomarker), without the other covariates. They check that the DAG's arrows point where associations are. The adjusted primary model (Firth Cox: age, sex, treatment, CRP, TMB/BRAF, CRP x treatment and TMB/BRAF x treatment; n = %d) gives CRP x treatment HR %s for PFS and %s for OS, and TMB/BRAF x treatment HR %s for PFS and %s for OS; the corresponding marginal tests give HR %s (T x CRP -> PFS), %s (T x CRP -> OS), %s (T x TMB -> PFS) and %s (T x TMB -> OS).",
            n_main,
            adj_text("TxCRP -> PFS"), adj_text("TxCRP -> OS"),
            adj_text("TxTMB -> PFS"), adj_text("TxTMB -> OS"),
            marg_text("TxCRP -> PFS"), marg_text("TxCRP -> OS"),
            marg_text("TxTMB -> PFS"), marg_text("TxTMB -> OS")),
    sprintf("Cohorts: tests on survival endpoints use the %d-patient complete-case cohort; Age and Sex on TMB/BRAF use the %d patients with a TMB/BRAF result. Tests into TLR use every patient with an observed TLR (n = %d), or the %d of them with a TMB/BRAF result for TMB_BRAF -> TLR. The TLR landmark cohorts of Figure 1 and Table 4 start from the TLR-classified patients within the complete-case cohort (n = %d).",
            n_main, n_tmb_known, n_tlr, n_tlr_marker, nrow(dd$data_tlr_landmark)),
    sprintf("PFS -> OS: Firth Cox for death with a time-dependent progression indicator (counting-process time), not PFS time as a baseline covariate. TLR -> PFS: Firth Cox on the %s landmark PFS cohort (n = %d of the %d TLR-classified patients in the %d-patient complete-case cohort, progression-free after the landmark, time from the landmark), the cohort used in Figure 1.",
            lm_label, nrow(dd$data_lm_pfs), nrow(dd$data_tlr_landmark), n_main),
    "Edges from the latent node U and into the interaction nodes are not tested."
  )
  if (!identical(n_tmb, n_main)) {
    stop("Table S2: expected the survival-endpoint tests on the complete-case cohort (n = ",
         n_main, "), found N = ", n_tmb, " for TMB_BRAF -> PFS.")
  }

  list(raw = df_raw, table = df_table, footnotes = v_footnotes)
}
