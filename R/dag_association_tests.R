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
#                             definitional. A statement is definitional when
#                             either tested variable is an interaction node
#                             whose parents (T and its biomarker) are all in
#                             the conditioning set: the node is then constant
#                             within every stratum and the statement holds by
#                             construction (issue #170). Two arm-balance flags
#                             are derived by rule: marginal statements involving
#                             T, and statements relating a variable to an
#                             interaction node given its biomarker but not T
#                             (within biomarker-positive patients the node
#                             equals T, so these compare the arms there).
#   dag_effect_modification_coverage(dag)
#                           - omitted edges from an interaction node into a
#                             child of T (TxCRP -> TLR, TxTMB -> TLR). Their
#                             assumption is "no effect modification of T -> C",
#                             which is checked separately and is NOT counted
#                             among the conditional independencies.
#
# run_dag_edge_tests(), run_dag_ci_tests() and
# run_dag_effect_modification_tests() map each derived statement to a
# prespecified test and STOP if the DAG implies a statement without a test or a
# test targets a statement the DAG no longer implies, so a DAG edit cannot
# silently leave the reports validating a stale graph.
#
# Every direct-edge test is MARGINAL: it contains the parent alone (for the
# interaction edges: T, the biomarker and their product) and no other
# covariate. run_adjusted_interaction_models() fits the adjusted primary model
# of clinical_effectiveness.qmd so that reports can reconcile the marginal
# interaction edges with it (issue #184).
#
# Two edge tests that were tautological (issue #155) are replaced:
#   PFS -> OS  : time-dependent progression indicator (survival::tmerge on the
#                progression-only event ProgressionExit / TTPwk) in a Firth Cox
#                model with counting-process time; entering PFS time as a
#                fixed covariate guaranteed HR < 1 because PFSwk <= OSwk.
#   TLR -> PFS : Firth Cox on the week-9 landmark PFS cohort (tlr_landmark.R)
#                built from tlr_landmark_base_cohort(), the same cohort and
#                time origin as clinical_effectiveness.qmd and
#                clin_effect_figure1.qmd; from randomisation, TLR status
#                requires surviving progression-free to the first scan.
#
# Requires dag_helpers.R (format_ci_statement), assoc_tests.R (run_* helpers),
# tlr_landmark.R (tlr_landmark_base_cohort, build_tlr_landmark_cohorts) and,
# for run_adjusted_interaction_models(), cox_extract.R (fit_firth_cox) to be
# sourced first.
# ===============================================================================

DAG_LATENT_NODES <- "U"
DAG_INTERACTION_NODES <- c(TxTMB = "TMB_BRAF", TxCRP = "CRP")
DAG_TREATMENT_NODE <- "T"
DAG_NODE_LABELS <- c(TMB_BRAF = "TMB/BRAF", CRP = "CRP")

#' Interaction nodes fixed by a conditioning set
#'
#' An interaction node (TxTMB = T * TMB_BRAF, TxCRP = T * CRP) is a
#' deterministic function of its parents, so it is constant within every
#' stratum of `z` when T and its biomarker are both in `z`.
#'
#' @param z Character vector of conditioning-set node names.
#' @return Names of the interaction nodes determined by `z`.
dag_determined_nodes <- function(z) {
  nodes <- names(DAG_INTERACTION_NODES)
  nodes[vapply(nodes, function(node) {
    all(c(DAG_TREATMENT_NODE, DAG_INTERACTION_NODES[[node]]) %in% z)
  }, logical(1))]
}

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
#'   `arm_balance` (marginal statement involving T), `stratum_arm_balance`
#'   (one variable is an interaction node whose biomarker is in Z while T is
#'   not; within biomarker-positive patients the node equals T and it is 0
#'   otherwise, so the statement compares the arms within that stratum),
#'   `balance_stratum` (that biomarker, else NA), `determined` (interaction
#'   node(s) among X and Y fixed by Z, comma-separated, else "") and
#'   `definitional` (`determined` is non-empty: the statement holds by
#'   construction; issue #170).
dag_ci_coverage <- function(dag_obj) {
  ci <- dagitty::impliedConditionalIndependencies(dag_obj)
  if (length(ci) == 0) {
    out <- data.frame(statement = character(0), X = character(0),
                      Y = character(0), arm_balance = logical(0),
                      stratum_arm_balance = logical(0),
                      balance_stratum = character(0), determined = character(0),
                      definitional = logical(0), stringsAsFactors = FALSE)
    out$Z <- list()
    return(out)
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
  out$balance_stratum <- vapply(seq_len(nrow(out)), function(i) {
    xy <- c(out$X[i], out$Y[i])
    z <- out$Z[[i]]
    nodes <- intersect(xy, names(DAG_INTERACTION_NODES))
    other <- setdiff(xy, nodes)
    if (length(nodes) != 1 || length(other) != 1 ||
        other == DAG_TREATMENT_NODE || DAG_TREATMENT_NODE %in% z ||
        !DAG_INTERACTION_NODES[[nodes]] %in% z) {
      return(NA_character_)
    }
    unname(DAG_INTERACTION_NODES[[nodes]])
  }, character(1))
  out$stratum_arm_balance <- !is.na(out$balance_stratum)
  out$determined <- vapply(seq_len(nrow(out)), function(i) {
    paste(intersect(c(out$X[i], out$Y[i]), dag_determined_nodes(out$Z[[i]])),
          collapse = ", ")
  }, character(1))
  out$definitional <- nzchar(out$determined)
  out
}

#' Omitted edges from an interaction node into a child of T
#'
#' For every interaction node I (T x biomarker) and every other child C of T,
#' the DAG either draws I -> C (tested as an edge) or omits it. An omitted
#' I -> C asserts that the biomarker does not modify the effect of T on C. The
#' implied statement "C _||_ I | {T, biomarker, ...}" cannot test that
#' assumption because I is constant within those strata (issue #170), so the
#' assumption is checked here as an interaction test and kept outside the
#' conditional-independence counts.
#'
#' @param dag_obj A dagitty object using the canonical node names.
#' @return Data frame with `item` ("I -> C"), `interaction`, `biomarker` and
#'   `child`, one row per omitted edge.
dag_effect_modification_coverage <- function(dag_obj) {
  e <- dagitty::edges(dag_obj)
  e <- e[e$e == "->", , drop = FALSE]
  from <- as.character(e$v)
  to <- as.character(e$w)
  t_children <- setdiff(to[from == DAG_TREATMENT_NODE], names(DAG_INTERACTION_NODES))
  rows <- lapply(names(DAG_INTERACTION_NODES), function(node) {
    if (!node %in% c(from, to)) return(NULL)
    omitted <- setdiff(t_children, to[from == node])
    if (length(omitted) == 0) return(NULL)
    data.frame(item = paste(node, "->", omitted), interaction = node,
               biomarker = unname(DAG_INTERACTION_NODES[[node]]), child = omitted,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) {
    out <- data.frame(item = character(0), interaction = character(0),
                      biomarker = character(0), child = character(0),
                      stringsAsFactors = FALSE)
  }
  out
}

#' Prepare the analysis datasets used by the DAG association tests
#'
#' @param data The data frame produced by scripts 02 and 03.
#' @param landmark Landmark id for the TLR -> PFS test (default "week9").
#' @return List with `data` (all rows, coded), `data_main` (complete cases for
#'   the survival and interaction models), `data_tlr` (every patient with an
#'   observed TLR; used by the logistic tests into TLR), `data_tlr_landmark`
#'   (TLR-classified patients within the complete-case cohort, from
#'   `tlr_landmark_base_cohort()`; the base of every TLR landmark cohort),
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

  # TLR landmark cohorts start from the TLR-classified patients within the
  # complete-case cohort (issue #184), as in clinical_effectiveness.qmd and
  # clin_effect_figure1.qmd; data_tlr keeps every observed TLR for the
  # non-landmark logistic tests.
  data_tlr_landmark <- tlr_landmark_base_cohort(data)
  lm <- build_tlr_landmark_cohorts(data_tlr_landmark, landmark = landmark)

  list(
    data = data,
    data_main = data_main,
    data_tlr = data_tlr,
    data_tlr_landmark = data_tlr_landmark,
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
#'
#' @param item The formatted statement.
#' @param determined Interaction node(s) fixed by the conditioning set, as in
#'   `dag_ci_coverage()$determined`.
definitional_ci_row <- function(item, determined) {
  nodes <- strsplit(determined, ", ", fixed = TRUE)[[1]]
  fixed_by <- vapply(nodes, function(node) {
    sprintf("%s is fixed by {%s}", node,
            paste(sort(c(DAG_TREATMENT_NODE, DAG_INTERACTION_NODES[[node]])),
                  collapse = ", "))
  }, character(1))
  tibble::tibble(
    Item = item,
    Test = paste0("Definitional (", paste(fixed_by, collapse = "; "), ")"),
    N = NA_integer_,
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
  # Every edge test is marginal (issue #184): the parent alone, or T, the
  # biomarker and their product for the interaction edges.
  tlr_pfs_label <- sprintf("Firth Cox, marginal (%s landmark PFS cohort)",
                           tolower(lm_label))
  cox_m <- "Firth Cox (marginal)"
  logit_m <- "Firth logistic (marginal)"
  interaction_label <- function(biomarker) {
    lab <- DAG_NODE_LABELS[[biomarker]]
    sprintf("Firth Cox interaction (marginal: T, %s, T x %s)", lab, lab)
  }

  specs <- list(
    "Age -> TMB_BRAF" = function()
      run_wilcox_test(dd$data, "tmb_braf", "Age", "Age -> TMB_BRAF",
                      "Wilcoxon rank-sum", c("TMB/BRAF=0", "TMB/BRAF=1")),
    "Age -> OS" = function()
      run_firth_cox_test(dd$data_main, Surv(OSwk, Death) ~ Age, "Age",
                         "Age -> OS", cox_m),
    "Sex -> TMB_BRAF" = function()
      run_fisher_test(dd$data, "sex_num", "tmb_braf", "Sex -> TMB_BRAF",
                      "Fisher's exact"),
    "TMB_BRAF -> PFS" = function()
      run_firth_cox_test(dd$data_main, Surv(PFSwk, Progression) ~ tmb_braf,
                         "tmb_braf", "TMB_BRAF -> PFS", cox_m),
    "TMB_BRAF -> OS" = function()
      run_firth_cox_test(dd$data_main, Surv(OSwk, Death) ~ tmb_braf,
                         "tmb_braf", "TMB_BRAF -> OS", cox_m),
    "TMB_BRAF -> TLR" = function()
      run_firth_logistic_test(dd$data_tlr, tlr ~ tmb_braf, "tmb_braf",
                              "TMB_BRAF -> TLR", logit_m),
    "CRP -> PFS" = function()
      run_firth_cox_test(dd$data_main, Surv(PFSwk, Progression) ~ crp, "crp",
                         "CRP -> PFS", cox_m),
    "CRP -> OS" = function()
      run_firth_cox_test(dd$data_main, Surv(OSwk, Death) ~ crp, "crp",
                         "CRP -> OS", cox_m),
    "CRP -> TLR" = function()
      run_firth_logistic_test(dd$data_tlr, tlr ~ crp, "crp", "CRP -> TLR",
                              logit_m),
    "T -> PFS" = function()
      run_firth_cox_test(dd$data_main, Surv(PFSwk, Progression) ~ T_num, "T_num",
                         "T -> PFS", cox_m),
    "T -> OS" = function()
      run_firth_cox_test(dd$data_main, Surv(OSwk, Death) ~ T_num, "T_num",
                         "T -> OS", cox_m),
    "T -> TLR" = function()
      run_firth_logistic_test(dd$data_tlr, tlr ~ T_num, "T_num", "T -> TLR",
                              logit_m),
    "TxTMB -> PFS" = function()
      run_firth_cox_test(dd$data_main,
                         Surv(PFSwk, Progression) ~ T_num + tmb_braf + T_num:tmb_braf,
                         "T_num:tmb_braf", "TxTMB -> PFS", interaction_label("TMB_BRAF")),
    "TxTMB -> OS" = function()
      run_firth_cox_test(dd$data_main,
                         Surv(OSwk, Death) ~ T_num + tmb_braf + T_num:tmb_braf,
                         "T_num:tmb_braf", "TxTMB -> OS", interaction_label("TMB_BRAF")),
    "TxCRP -> PFS" = function()
      run_firth_cox_test(dd$data_main,
                         Surv(PFSwk, Progression) ~ T_num + crp + T_num:crp,
                         "T_num:crp", "TxCRP -> PFS", interaction_label("CRP")),
    "TxCRP -> OS" = function()
      run_firth_cox_test(dd$data_main,
                         Surv(OSwk, Death) ~ T_num + crp + T_num:crp,
                         "T_num:crp", "TxCRP -> OS", interaction_label("CRP")),
    "TLR -> PFS" = function()
      run_firth_cox_test(dd$data_lm_pfs, Surv(PFSwk_lm, Progression) ~ tlr, "tlr",
                         "TLR -> PFS", tlr_pfs_label),
    "PFS -> OS" = function()
      run_time_dependent_cox_test(dd$data_td, "PFS -> OS",
                                  "Firth Cox, time-dependent progression (marginal)")
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
#'   Stratum_Balance, Balance_Stratum, Definitional, Conclusion) in dagitty
#'   order, with attribute `coverage`. Definitional rows have no test (N,
#'   p_value NA).
run_dag_ci_tests <- function(dag_obj, dd, seed = 123L) {
  cov <- dag_ci_coverage(dag_obj)
  d <- dd$data

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
    "TMB_BRAF _||_ TxCRP | {CRP}" = function()
      run_mh_test(d, "tmb_braf", "TxCRP", "crp", "TMB_BRAF _||_ TxCRP | {CRP}",
                  "Mantel-Haenszel")
  )

  # Definitional statements get a row without a test. A prespecified test for
  # a definitional statement is an error: it would test something else (the
  # TLR x interaction-node statements tested effect modification of T -> TLR,
  # now run_dag_effect_modification_tests(); issue #170).
  defn <- cov$statement[cov$definitional]
  tested_defn <- intersect(names(specs), defn)
  if (length(tested_defn) > 0) {
    stop("Definitional conditional independencies have an empirical test: ",
         paste(tested_defn, collapse = "; "))
  }
  for (i in which(cov$definitional)) {
    specs[[cov$statement[i]]] <- local({
      stmt <- cov$statement[i]
      det <- cov$determined[i]
      function() definitional_ci_row(stmt, det)
    })
  }

  assert_spec_coverage(names(specs), cov$statement, "conditional independencies")
  res <- dplyr::bind_rows(lapply(cov$statement, function(s) specs[[s]]()))
  idx <- match(res$Item, cov$statement)
  res$Randomization_Check <- cov$arm_balance[idx]
  res$Stratum_Balance <- cov$stratum_arm_balance[idx]
  res$Balance_Stratum <- cov$balance_stratum[idx]
  res$Definitional <- cov$definitional[idx]
  stratum_note <- ifelse(
    res$Stratum_Balance,
    paste0(" (within ", DAG_NODE_LABELS[res$Balance_Stratum], "-positive patients)"),
    ""
  )
  res$Conclusion <- ifelse(
    res$Definitional, "Holds by construction",
    paste0(ci_support_label(res$p_value,
                            res$Randomization_Check | res$Stratum_Balance),
           stratum_note)
  )
  attr(res, "coverage") <- cov
  res
}

#' Test the omitted interaction-to-child edges (effect modification of T -> C)
#'
#' Each omitted edge I -> C of `dag_effect_modification_coverage()` is tested
#' as the product term of a Firth logistic (binary C) model with T, the
#' biomarker and T x biomarker, in every patient with the child observed.
#' These rows are not conditional-independence tests and are not counted as
#' such (issue #170).
#'
#' @inheritParams run_dag_edge_tests
#' @return Tibble (Item, Test, N, Effect, p_value, Conclusion) in coverage
#'   order, with attribute `coverage`.
run_dag_effect_modification_tests <- function(dag_obj, dd) {
  cov <- dag_effect_modification_coverage(dag_obj)
  label <- function(biomarker) {
    lab <- DAG_NODE_LABELS[[biomarker]]
    sprintf("Firth logistic interaction (T, %s, T x %s)", lab, lab)
  }
  specs <- list(
    "TxCRP -> TLR" = function()
      run_firth_logistic_test(dd$data_tlr, tlr ~ T_num + crp + T_num:crp, "T_num:crp",
                              "TxCRP -> TLR", label("CRP")),
    "TxTMB -> TLR" = function()
      run_firth_logistic_test(dd$data_tlr, tlr ~ T_num + tmb_braf + T_num:tmb_braf,
                              "T_num:tmb_braf", "TxTMB -> TLR", label("TMB_BRAF"))
  )
  assert_spec_coverage(names(specs), cov$item, "omitted interaction edges")
  if (nrow(cov) == 0) {
    res <- tibble::tibble(Item = character(0), Test = character(0),
                          N = integer(0), Effect = character(0),
                          p_value = numeric(0), Conclusion = character(0))
  } else {
    res <- dplyr::bind_rows(lapply(cov$item, function(s) specs[[s]]()))
    res$Conclusion <- ifelse(
      is.na(res$p_value), "Not estimable",
      ifelse(res$p_value < 0.05,
             "Effect modification detected (p < 0.05)",
             "No effect modification detected (p >= 0.05)")
    )
  }
  attr(res, "coverage") <- cov
  res
}

#' Adjusted primary-model interaction estimates (reconciliation, issue #184)
#'
#' Fits the primary clinical model of clinical_effectiveness.qmd, a Firth Cox
#' model `Age + sex + Rx + crp + tmb_braf + crp:Rx + tmb_braf:Rx`, for OS and
#' PFS on the complete-case cohort `dd$data_main` with `fit_firth_cox()`
#' (cox_extract.R), so that reports can set the marginal interaction edges
#' against the adjusted estimates of the main text.
#'
#' @inheritParams run_dag_edge_tests
#' @return Tibble with `Edge` (e.g. "TxCRP -> PFS"), `Endpoint`, `Term`
#'   ("CRP x treatment", "TMB/BRAF x treatment"), `N`, `HR`, `CI_Lower`,
#'   `CI_Upper`, `p_value` (per-coefficient penalized likelihood ratio test)
#'   and `Effect` (formatted HR and 95% CI).
run_adjusted_interaction_models <- function(dd) {
  if (!exists("fit_firth_cox", mode = "function")) {
    stop("run_adjusted_interaction_models() needs fit_firth_cox(); source R/cox_extract.R.")
  }
  d <- dd$data_main
  formulas <- list(
    OS = Surv(OSwk, Death) ~ Age + sex + Rx + crp + tmb_braf + crp:Rx + tmb_braf:Rx,
    PFS = Surv(PFSwk, Progression) ~ Age + sex + Rx + crp + tmb_braf +
      crp:Rx + tmb_braf:Rx
  )
  terms <- c(TxCRP = "crp", TxTMB = "tmb_braf")
  rows <- list()
  for (endpoint in names(formulas)) {
    fit <- fit_firth_cox(formulas[[endpoint]], d,
                         paste("adjusted primary model", endpoint))
    if (is.null(fit)) stop("Adjusted primary model failed for ", endpoint)
    coefs <- names(fit$coefficients)
    for (node in names(terms)) {
      idx <- which(grepl(":", coefs, fixed = TRUE) &
                     grepl(paste0("(^|:)", terms[[node]], "($|:)"), coefs))
      if (length(idx) != 1) {
        stop("Expected one ", terms[[node]], " x Rx coefficient in the adjusted ",
             endpoint, " model; found ", length(idx))
      }
      hr <- exp(unname(fit$coefficients[idx]))
      lo <- unname(fit$ci.lower[idx])
      hi <- unname(fit$ci.upper[idx])
      rows[[length(rows) + 1]] <- tibble::tibble(
        Edge = paste(node, "->", endpoint),
        Endpoint = endpoint,
        Term = paste(DAG_NODE_LABELS[[DAG_INTERACTION_NODES[[node]]]], "x treatment"),
        N = nrow(d),
        HR = hr, CI_Lower = lo, CI_Upper = hi,
        p_value = unname(fit$prob[idx]),
        Effect = fmt_hr(hr, lo, hi)
      )
    }
  }
  dplyr::bind_rows(rows)
}

message("DAG association-test helpers loaded from dag_association_tests.R")
