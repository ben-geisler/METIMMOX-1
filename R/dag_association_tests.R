# ===============================================================================
# DAG ASSOCIATION TESTS DERIVED FROM THE CANONICAL DAG (issue #155)
# ===============================================================================
# dag_associations.qmd and clin_effect_table_s2.qmd previously typed their own
# edge and conditional-independence (CI) lists, and the typed DAG lacked the
# TxCRP node, so the "empirical companion" validated a different graph from
# dag.qmd. This file derives both lists from the `dag` object in dag_helpers.R:
#
#   dag_edge_coverage(dag)  - every edge, split into testable / latent /
#                             definitional (edges INTO a constructed interaction
#                             node are definitional: TxTMB = T * TMB_BRAF)
#   dag_ci_coverage(dag)    - every implied CI, split into testable /
#                             definitional (two interaction nodes given T are
#                             constant within every stratum) with the arm-balance
#                             flag derived by rule (marginal statement involving T)
#
# run_dag_edge_tests() and run_dag_ci_tests() map each derived statement to a
# prespecified test and STOP if the DAG implies a statement without a test or a
# test targets a statement the DAG no longer implies, so a DAG edit cannot
# silently leave the reports validating a stale graph.
#
# Two edge tests that were tautological (issue #155) are replaced:
#   PFS -> OS  : time-dependent progression indicator (survival::tmerge on the
#                progression-only event ProgressionExit / TTPwk) in a Firth Cox
#                model with counting-process time; entering PFS time as a
#                fixed covariate guaranteed HR < 1 because PFSwk <= OSwk.
#   TLR -> PFS : Firth Cox on the week-9 landmark PFS cohort (tlr_landmark.R),
#                the same cohort and time origin as clinical_effectiveness.qmd
#                and clin_effect_figure1.qmd; from randomisation, TLR status
#                requires surviving progression-free to the first scan.
#
# Requires dag_helpers.R (format_ci_statement), assoc_tests.R (run_* helpers)
# and tlr_landmark.R (build_tlr_landmark_cohorts) to be sourced first.
# ===============================================================================

DAG_LATENT_NODES <- "U"
DAG_INTERACTION_NODES <- c(TxTMB = "TMB_BRAF", TxCRP = "CRP")
DAG_TREATMENT_NODE <- "T"

#' Edge labels "X -> Y" of a dagitty object, in dagitty order
dag_edge_labels <- function(dag_obj) {
  e <- dagitty::edges(dag_obj)
  e <- e[e$e == "->", , drop = FALSE]
  paste(as.character(e$v), "->", as.character(e$w))
}

#' Split the edges of a DAG into testable, latent and definitional sets
#'
#' @param dag_obj A dagitty object using the canonical node names.
#' @return List with character vectors `testable`, `latent` and `definitional`.
dag_edge_coverage <- function(dag_obj) {
  e <- dagitty::edges(dag_obj)
  e <- e[e$e == "->", , drop = FALSE]
  labels <- paste(as.character(e$v), "->", as.character(e$w))
  latent <- as.character(e$v) %in% DAG_LATENT_NODES
  definitional <- as.character(e$w) %in% names(DAG_INTERACTION_NODES)
  list(
    all = labels,
    testable = labels[!latent & !definitional],
    latent = labels[latent],
    definitional = labels[definitional]
  )
}

#' Implied conditional independencies of a DAG as formatted statements
#'
#' @return Data frame with `statement`, `X`, `Y`, `Z` (list column),
#'   `arm_balance` (marginal statement involving T) and `definitional` (both
#'   variables are interaction nodes and T is in the conditioning set, so at
#'   least one of them is constant within every stratum).
dag_ci_coverage <- function(dag_obj) {
  ci <- dagitty::impliedConditionalIndependencies(dag_obj)
  if (length(ci) == 0) {
    return(data.frame(statement = character(0), X = character(0),
                      Y = character(0), arm_balance = logical(0),
                      definitional = logical(0), stringsAsFactors = FALSE))
  }
  out <- data.frame(
    statement = vapply(ci, format_ci_statement, character(1)),
    X = vapply(ci, function(s) s$X, character(1)),
    Y = vapply(ci, function(s) s$Y, character(1)),
    stringsAsFactors = FALSE
  )
  out$Z <- lapply(ci, function(s) sort(as.character(s$Z)))
  out$arm_balance <- vapply(seq_len(nrow(out)), function(i) {
    length(out$Z[[i]]) == 0 && DAG_TREATMENT_NODE %in% c(out$X[i], out$Y[i])
  }, logical(1))
  out$definitional <- vapply(seq_len(nrow(out)), function(i) {
    all(c(out$X[i], out$Y[i]) %in% names(DAG_INTERACTION_NODES)) &&
      DAG_TREATMENT_NODE %in% out$Z[[i]]
  }, logical(1))
  out
}

#' Prepare the analysis datasets used by the DAG association tests
#'
#' @param data The data frame produced by scripts 02 and 03.
#' @param landmark Landmark id for the TLR -> PFS test (default "week9").
#' @return List with `data` (all rows, coded), `data_main` (complete cases for
#'   the survival and interaction models), `data_tlr` (observed TLR),
#'   `data_lm_pfs` (landmark PFS cohort), `landmark` (the cohort object) and
#'   `data_td` (counting-process data with a time-dependent progression flag).
prepare_dag_data <- function(data, landmark = "week9") {
  data$sex_num <- as.numeric(as.character(data$sex))
  data$T_num <- as.numeric(data$Rx == "Experimental arm")
  data$TxTMB <- data$T_num * data$tmb_braf
  data$TxCRP <- data$T_num * data$crp
  data$tmb_braf_fac <- factor(data$tmb_braf)
  data$crp_fac <- factor(data$crp)
  data$TxTMB_fac <- factor(data$TxTMB)
  data$TxCRP_fac <- factor(data$TxCRP)
  data$T_fac <- factor(data$T_num)
  data$sex_fac <- factor(data$sex_num)

  data_tlr <- data[!is.na(data$tlr), , drop = FALSE]
  main_vars <- c("Age", "sex", "Rx", "crp", "tmb_braf", "T_num", "TxTMB", "TxCRP",
                 "OSwk", "Death", "PFSwk", "Progression")
  data_main <- data[stats::complete.cases(data[, main_vars]), , drop = FALSE]

  lm <- build_tlr_landmark_cohorts(data_tlr, landmark = landmark)

  list(
    data = data,
    data_main = data_main,
    data_tlr = data_tlr,
    landmark = lm,
    data_lm_pfs = lm$pfs,
    data_td = build_time_dependent_progression(data_main)
  )
}

#' Counting-process data with a time-dependent progression indicator
#'
#' Uses the progression-only trial event (`ProgressionExit`, `TTPwk`), not the
#' composite PFS endpoint, because death is the outcome of the OS model.
#'
#' @param df Data frame with `ID`, `OSwk`, `Death`, `TTPwk`, `ProgressionExit`.
build_time_dependent_progression <- function(df) {
  required <- c("ID", "OSwk", "Death", "TTPwk", "ProgressionExit")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop("build_time_dependent_progression(): missing column(s): ",
         paste(missing, collapse = ", "),
         " (ProgressionExit/TTPwk are kept by 03_biomarker_strategies.R).")
  }
  base <- df[, required, drop = FALSE]
  base$ID <- as.character(base$ID)
  prog <- base[base$ProgressionExit == 1, c("ID", "TTPwk"), drop = FALSE]
  td <- survival::tmerge(base, base, id = ID,
                         death = event(OSwk, Death))
  td <- survival::tmerge(td, prog, id = ID,
                         progressed = tdc(TTPwk))
  if (any(td$tstart >= td$tstop)) {
    stop("build_time_dependent_progression(): zero-length interval produced.")
  }
  td
}

#' Firth Cox test of a time-dependent progression indicator on OS
run_time_dependent_cox_test <- function(td, item, test_label) {
  fit <- coxphf::coxphf(
    survival::Surv(tstart, tstop, death) ~ progressed, data = td,
    firth = TRUE, pl = TRUE, maxit = 100, maxstep = 0.1
  )
  tibble::tibble(
    Item = item,
    Test = test_label,
    N = length(unique(td$ID)),
    Effect = fmt_hr(
      exp(unname(fit$coefficients["progressed"])),
      unname(fit$ci.lower["progressed"]),
      unname(fit$ci.upper["progressed"])
    ),
    p_value = unname(fit$prob["progressed"])
  )
}

#' Row for a statement that holds by construction
definitional_ci_row <- function(item, n) {
  tibble::tibble(
    Item = item,
    Test = "Definitional (constant within every stratum)",
    N = n,
    Effect = "--",
    p_value = NA_real_
  )
}

assert_spec_coverage <- function(spec_names, derived, what) {
  missing <- setdiff(derived, spec_names)
  extra <- setdiff(spec_names, derived)
  if (length(missing) > 0 || length(extra) > 0) {
    stop(
      "DAG ", what, " and the prespecified tests disagree.",
      if (length(missing) > 0) paste0("\n  Implied by the DAG but untested: ",
                                      paste(missing, collapse = "; ")),
      if (length(extra) > 0) paste0("\n  Tested but not implied by the DAG: ",
                                    paste(extra, collapse = "; "))
    )
  }
  invisible(TRUE)
}

#' Run every testable direct-edge test implied by the DAG
#'
#' @param dag_obj The canonical dagitty object.
#' @param dd Output of `prepare_dag_data()`.
#' @return Tibble (Item, Test, N, Effect, p_value, Conclusion) in DAG edge
#'   order, with attribute `coverage` (the `dag_edge_coverage()` list).
run_dag_edge_tests <- function(dag_obj, dd) {
  cov <- dag_edge_coverage(dag_obj)
  lm_label <- dd$landmark$spec$short_label
  tlr_pfs_label <- sprintf("Firth Cox (%s landmark PFS cohort)", tolower(lm_label))

  specs <- list(
    "Age -> TMB_BRAF" = function()
      run_wilcox_test(dd$data, "tmb_braf", "Age", "Age -> TMB_BRAF",
                      "Wilcoxon rank-sum", c("TMB/BRAF=0", "TMB/BRAF=1")),
    "Age -> OS" = function()
      run_firth_cox_test(dd$data_main, Surv(OSwk, Death) ~ Age, "Age",
                         "Age -> OS", "Firth Cox"),
    "Sex -> TMB_BRAF" = function()
      run_fisher_test(dd$data, "sex_num", "tmb_braf", "Sex -> TMB_BRAF",
                      "Fisher's exact"),
    "TMB_BRAF -> PFS" = function()
      run_firth_cox_test(dd$data_main, Surv(PFSwk, Progression) ~ tmb_braf,
                         "tmb_braf", "TMB_BRAF -> PFS", "Firth Cox"),
    "TMB_BRAF -> OS" = function()
      run_firth_cox_test(dd$data_main, Surv(OSwk, Death) ~ tmb_braf,
                         "tmb_braf", "TMB_BRAF -> OS", "Firth Cox"),
    "TMB_BRAF -> TLR" = function()
      run_firth_logistic_test(dd$data_tlr, tlr ~ tmb_braf, "tmb_braf",
                              "TMB_BRAF -> TLR", "Firth logistic"),
    "CRP -> PFS" = function()
      run_firth_cox_test(dd$data_main, Surv(PFSwk, Progression) ~ crp, "crp",
                         "CRP -> PFS", "Firth Cox"),
    "CRP -> OS" = function()
      run_firth_cox_test(dd$data_main, Surv(OSwk, Death) ~ crp, "crp",
                         "CRP -> OS", "Firth Cox"),
    "CRP -> TLR" = function()
      run_firth_logistic_test(dd$data_tlr, tlr ~ crp, "crp", "CRP -> TLR",
                              "Firth logistic"),
    "T -> PFS" = function()
      run_firth_cox_test(dd$data_main, Surv(PFSwk, Progression) ~ T_num, "T_num",
                         "T -> PFS", "Firth Cox"),
    "T -> OS" = function()
      run_firth_cox_test(dd$data_main, Surv(OSwk, Death) ~ T_num, "T_num",
                         "T -> OS", "Firth Cox"),
    "T -> TLR" = function()
      run_firth_logistic_test(dd$data_tlr, tlr ~ T_num, "T_num", "T -> TLR",
                              "Firth logistic"),
    "TxTMB -> PFS" = function()
      run_firth_cox_test(dd$data_main,
                         Surv(PFSwk, Progression) ~ T_num + tmb_braf + T_num:tmb_braf,
                         "T_num:tmb_braf", "TxTMB -> PFS", "Firth Cox interaction"),
    "TxTMB -> OS" = function()
      run_firth_cox_test(dd$data_main,
                         Surv(OSwk, Death) ~ T_num + tmb_braf + T_num:tmb_braf,
                         "T_num:tmb_braf", "TxTMB -> OS", "Firth Cox interaction"),
    "TxCRP -> PFS" = function()
      run_firth_cox_test(dd$data_main,
                         Surv(PFSwk, Progression) ~ T_num + crp + T_num:crp,
                         "T_num:crp", "TxCRP -> PFS", "Firth Cox interaction"),
    "TxCRP -> OS" = function()
      run_firth_cox_test(dd$data_main,
                         Surv(OSwk, Death) ~ T_num + crp + T_num:crp,
                         "T_num:crp", "TxCRP -> OS", "Firth Cox interaction"),
    "TLR -> PFS" = function()
      run_firth_cox_test(dd$data_lm_pfs, Surv(PFSwk_lm, Progression) ~ tlr, "tlr",
                         "TLR -> PFS", tlr_pfs_label),
    "PFS -> OS" = function()
      run_time_dependent_cox_test(dd$data_td, "PFS -> OS",
                                  "Firth Cox, time-dependent progression")
  )

  assert_spec_coverage(names(specs), cov$testable, "edges")
  res <- dplyr::bind_rows(lapply(cov$testable, function(e) specs[[e]]()))
  res$Conclusion <- edge_support_label(res$p_value)
  attr(res, "coverage") <- cov
  res
}

#' Run every conditional-independence test implied by the DAG
#'
#' @inheritParams run_dag_edge_tests
#' @param seed Seed for the conditional permutation test.
#' @return Tibble (Item, Test, N, Effect, p_value, Randomization_Check,
#'   Definitional, Conclusion) in dagitty order, with attribute `coverage`.
run_dag_ci_tests <- function(dag_obj, dd, seed = 123L) {
  cov <- dag_ci_coverage(dag_obj)
  d <- dd$data
  n_all <- nrow(d)

  specs <- list(
    "Age _||_ CRP" = function()
      run_wilcox_test(d, "crp", "Age", "Age _||_ CRP", "Wilcoxon rank-sum",
                      c("CRP=0", "CRP=1")),
    "Age _||_ Sex" = function()
      run_wilcox_test(d, "sex_num", "Age", "Age _||_ Sex", "Wilcoxon rank-sum",
                      c("Sex=0", "Sex=1")),
    "Age _||_ T" = function()
      run_wilcox_test(d, "T_num", "Age", "Age _||_ T",
                      "Wilcoxon rank-sum (randomization check)", c("T=0", "T=1")),
    "Age _||_ TLR | {CRP, TMB_BRAF}" = function()
      run_firth_logistic_test(dd$data_tlr, tlr ~ crp + tmb_braf + Age, "Age",
                              "Age _||_ TLR | {CRP, TMB_BRAF}", "Firth logistic"),
    "Age _||_ TxCRP" = function()
      run_wilcox_test(d, "TxCRP", "Age", "Age _||_ TxCRP", "Wilcoxon rank-sum",
                      c("TxCRP=0", "TxCRP=1")),
    "Age _||_ TxTMB | {TMB_BRAF}" = function()
      run_permutation_test(d, Age ~ TxTMB_fac | tmb_braf_fac,
                           "Age _||_ TxTMB | {TMB_BRAF}",
                           "Conditional permutation test", seed = seed),
    "CRP _||_ Sex" = function()
      run_fisher_test(d, "crp", "sex_num", "CRP _||_ Sex", "Fisher's exact"),
    "CRP _||_ T" = function()
      run_fisher_test(d, "crp", "T_num", "CRP _||_ T",
                      "Fisher's exact (arm balance at week 4)"),
    "CRP _||_ TxTMB | {TMB_BRAF}" = function()
      run_mh_test(d, "crp", "TxTMB", "tmb_braf", "CRP _||_ TxTMB | {TMB_BRAF}",
                  "Mantel-Haenszel"),
    "Sex _||_ T" = function()
      run_fisher_test(d, "sex_num", "T_num", "Sex _||_ T",
                      "Fisher's exact (randomization check)"),
    "Sex _||_ TLR | {CRP, TMB_BRAF}" = function()
      run_firth_logistic_test(dd$data_tlr, tlr ~ crp + tmb_braf + sex_num, "sex_num",
                              "Sex _||_ TLR | {CRP, TMB_BRAF}", "Firth logistic"),
    "Sex _||_ TxCRP" = function()
      run_fisher_test(d, "sex_num", "TxCRP", "Sex _||_ TxCRP", "Fisher's exact"),
    "Sex _||_ TxTMB | {TMB_BRAF}" = function()
      run_mh_test(d, "sex_num", "TxTMB", "tmb_braf", "Sex _||_ TxTMB | {TMB_BRAF}",
                  "Mantel-Haenszel"),
    "T _||_ TMB_BRAF" = function()
      run_fisher_test(d, "T_num", "tmb_braf", "T _||_ TMB_BRAF",
                      "Fisher's exact (randomization check)"),
    "TLR _||_ TxCRP | {CRP, T}" = function()
      run_firth_logistic_test(dd$data_tlr, tlr ~ T_num + crp + TxCRP, "TxCRP",
                              "TLR _||_ TxCRP | {CRP, T}", "Firth logistic"),
    "TLR _||_ TxTMB | {T, TMB_BRAF}" = function()
      run_firth_logistic_test(dd$data_tlr, tlr ~ T_num + tmb_braf + TxTMB, "TxTMB",
                              "TLR _||_ TxTMB | {T, TMB_BRAF}", "Firth logistic"),
    "TMB_BRAF _||_ TxCRP | {CRP}" = function()
      run_mh_test(d, "tmb_braf", "TxCRP", "crp", "TMB_BRAF _||_ TxCRP | {CRP}",
                  "Mantel-Haenszel")
  )

  # Definitional statements get a row without a test.
  for (s in cov$statement[cov$definitional]) {
    specs[[s]] <- local({
      stmt <- s
      function() definitional_ci_row(stmt, n_all)
    })
  }

  assert_spec_coverage(names(specs), cov$statement, "conditional independencies")
  res <- dplyr::bind_rows(lapply(cov$statement, function(s) specs[[s]]()))
  res$Randomization_Check <- cov$arm_balance[match(res$Item, cov$statement)]
  res$Definitional <- cov$definitional[match(res$Item, cov$statement)]
  res$Conclusion <- ifelse(
    res$Definitional, "Holds by construction",
    ci_support_label(res$p_value, res$Randomization_Check)
  )
  attr(res, "coverage") <- cov
  res
}

message("DAG association-test helpers loaded from dag_association_tests.R")
