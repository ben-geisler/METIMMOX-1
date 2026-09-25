# Contracts for the DAG association tests and the TLR landmark cohorts
# (issue #155). Requires the confidential trial data (data/tidy/METIMMOX.rds);
# the structural checks that need no data run first, the data-dependent checks
# are skipped with a message when the file is absent.
#
# Checked:
#   1. The edge and conditional-independence (CI) lists tested by the DAG
#      reports are derived from the canonical `dag` object: every implied
#      statement has a test or is classified latent/definitional, and both
#      TxCRP edges and all TxCRP independencies are covered.
#   2. A DAG edit that adds an untested edge makes run_dag_edge_tests() stop.
#   3. Landmark cohorts contain only patients with the endpoint time strictly
#      after the landmark; the per-patient landmark retains no first-scan
#      progressor and no patient whose TLR was read after the landmark; the
#      week-9 landmark reports how many of each it retains.
#   4. The PFS -> OS edge is tested with a time-dependent progression
#      indicator (counting-process rows), not with PFS time as a covariate, and
#      the TLR -> PFS edge is tested on the landmark cohort.
#   5. (issue #170) A CI statement is definitional whenever either tested
#      variable is an interaction node whose parents are all in the
#      conditioning set; the synthetic case of the issue (TxCRP constant in
#      every {CRP, T} stratum, yet the old interaction test rejects) is
#      classified as holding by construction, and the TLR interaction tests
#      are reported as effect-modification checks outside the CI counts.
#   6. (issue #184) Every TLR landmark cohort starts from the TLR-classified
#      patients within the complete-case cohort (tlr_landmark_base_cohort()),
#      the edge tests are labelled marginal, statements relating a variable to
#      an interaction node given its biomarker are arm-balance checks within
#      that stratum, and the adjusted primary-model interaction terms are
#      available for reconciliation.
#   7. (issues #176, #182) Clinical-report helpers: the exit reasons of the
#      patients without a TLR value (tlr_missing_exit_summary()) and the
#      coxphf convergence diagnostics (coxphf_convergence()); the week-12
#      landmark still retains a patient scanned after week 12.

suppressPackageStartupMessages({
  library(here); library(dplyr); library(survival); library(dagitty)
})

source(here("R/pfs_endpoint.R"))
source(here("R/dag_helpers.R"))
source(here("R/assoc_tests.R"))
source(here("R/cox_extract.R"))
source(here("R/tlr_landmark.R"))
source(here("R/dag_association_tests.R"))

# ---------------------------------------------------------------------------
# 1. Coverage derived from the DAG (no data needed)
# ---------------------------------------------------------------------------
edge_cov <- dag_edge_coverage(dag)
ci_cov <- dag_ci_coverage(dag)

stopifnot(
  length(edge_cov$all) == length(edge_cov$testable) + length(edge_cov$latent) +
    length(edge_cov$definitional),
  all(c("TxCRP -> PFS", "TxCRP -> OS", "TxTMB -> PFS", "TxTMB -> OS",
        "TLR -> PFS", "PFS -> OS") %in% edge_cov$testable),
  all(c("T -> TxCRP", "CRP -> TxCRP", "T -> TxTMB", "TMB_BRAF -> TxTMB") %in%
        edge_cov$definitional),
  all(grepl("^U ->", edge_cov$latent)),
  nrow(ci_cov) == length(impliedConditionalIndependencies(dag)),
  # Every implied statement mentioning TxCRP is in the derived list.
  sum(grepl("TxCRP", ci_cov$statement)) == 6,
  # Issue #170: a statement holds by construction whenever an interaction
  # node among X and Y has both parents (T and its biomarker) in Z. This
  # covers the two statements relating the interaction nodes given T and the
  # two TLR statements; nothing else is definitional.
  setequal(ci_cov$statement[ci_cov$definitional],
           c("TLR _||_ TxCRP | {CRP, T}", "TLR _||_ TxTMB | {T, TMB_BRAF}",
             "TxCRP _||_ TxTMB | {T, TMB_BRAF}", "TxCRP _||_ TxTMB | {CRP, T}")),
  identical(ci_cov$determined[ci_cov$statement == "TLR _||_ TxCRP | {CRP, T}"], "TxCRP"),
  identical(ci_cov$determined[ci_cov$statement == "TLR _||_ TxTMB | {T, TMB_BRAF}"], "TxTMB"),
  sum(ci_cov$definitional) == 4, sum(!ci_cov$definitional) == 15,
  # Arm-balance checks are exactly the marginal statements involving T.
  setequal(ci_cov$statement[ci_cov$arm_balance],
           c("Age _||_ T", "CRP _||_ T", "Sex _||_ T", "T _||_ TMB_BRAF")),
  # Issue #184: a variable against an interaction node given its biomarker but
  # not T compares the arms within biomarker-positive patients (the node equals
  # T there and is 0 otherwise); CRP _||_ TxTMB | {TMB_BRAF} is one of them.
  setequal(ci_cov$statement[ci_cov$stratum_arm_balance],
           c("Age _||_ TxTMB | {TMB_BRAF}", "CRP _||_ TxTMB | {TMB_BRAF}",
             "Sex _||_ TxTMB | {TMB_BRAF}", "TMB_BRAF _||_ TxCRP | {CRP}")),
  identical(ci_cov$balance_stratum[ci_cov$statement == "CRP _||_ TxTMB | {TMB_BRAF}"],
            "TMB_BRAF"),
  !any(ci_cov$stratum_arm_balance & (ci_cov$definitional | ci_cov$arm_balance)),
  # The omitted interaction edges into TLR are effect-modification checks.
  setequal(dag_effect_modification_coverage(dag)$item, c("TxCRP -> TLR", "TxTMB -> TLR"))
)

# The simplified DAG (issue #184) draws the same parent-child relations among
# its nodes as the full DAG, with PFS and OS merged into Survival: in
# particular Age -> TMB_BRAF, Age -> Survival, Sex -> TMB_BRAF and no
# Sex -> Survival.
children_of <- function(g, v) {
  e <- dagitty::edges(g)
  sort(unique(as.character(e$w[e$e == "->" & as.character(e$v) == v])))
}
simple_nodes <- setdiff(names(dag_simple), "Survival")
for (v in simple_nodes) {
  full_ch <- intersect(children_of(dag, v), c(simple_nodes, "PFS", "OS"))
  full_ch <- sort(unique(ifelse(full_ch %in% c("PFS", "OS"), "Survival", full_ch)))
  if (!identical(children_of(dag_simple, v), full_ch)) {
    stop("Simplified DAG disagrees with the full DAG for the children of ", v, ": ",
         paste(children_of(dag_simple, v), collapse = ", "), " versus ",
         paste(full_ch, collapse = ", "))
  }
}
stopifnot(!"Demo" %in% names(dag_simple),
          identical(children_of(dag_simple, "Sex"), "TMB_BRAF"),
          identical(children_of(dag_simple, "Age"), c("Survival", "TMB_BRAF")))

# ---------------------------------------------------------------------------
# 1b. Synthetic case from issue #170 (no trial data needed)
# ---------------------------------------------------------------------------
# A minimal DAG in which TxCRP = T * CRP is a child of T and CRP only. Given
# {CRP, T}, TxCRP has exactly one value in every stratum, so
# "TLR _||_ TxCRP | {CRP, T}" holds by construction. The former test (Firth
# logistic, TLR ~ T + CRP + TxCRP) nevertheless rejects it on data with a
# strong T x CRP effect on TLR: it tests effect modification of T -> TLR, not
# the stated independence.
mini_dag <- dagitty("dag { T -> TLR\n CRP -> TLR\n T -> TxCRP\n CRP -> TxCRP }")
mini_cov <- dag_ci_coverage(mini_dag)
stopifnot(
  setequal(mini_cov$statement, c("CRP _||_ T", "TLR _||_ TxCRP | {CRP, T}")),
  isTRUE(mini_cov$definitional[mini_cov$statement == "TLR _||_ TxCRP | {CRP, T}"]),
  !mini_cov$definitional[mini_cov$statement == "CRP _||_ T"],
  identical(dag_effect_modification_coverage(mini_dag)$item, "TxCRP -> TLR"),
  # One interaction node is fixed only when BOTH parents are conditioned on.
  identical(dag_determined_nodes(c("CRP", "T")), "TxCRP"),
  length(dag_determined_nodes("T")) == 0, length(dag_determined_nodes("CRP")) == 0
)
toy <- expand.grid(T_num = 0:1, crp = 0:1, replicate = seq_len(100))
toy$TxCRP <- toy$T_num * toy$crp
toy$tlr <- as.numeric(toy$replicate <= ifelse(toy$TxCRP == 1, 90, 10))
stopifnot(all(tapply(toy$TxCRP, interaction(toy$T_num, toy$crp),
                     function(x) length(unique(x))) == 1))
toy_test <- run_firth_logistic_test(toy, tlr ~ T_num + crp + T_num:crp, "T_num:crp",
                                    "TxCRP -> TLR", "Firth logistic interaction")
stopifnot(toy_test$p_value < 1e-6)
defn_row <- definitional_ci_row("TLR _||_ TxCRP | {CRP, T}", "TxCRP")
stopifnot(is.na(defn_row$p_value), is.na(defn_row$N),
          grepl("TxCRP is fixed by {CRP, T}", defn_row$Test, fixed = TRUE))

# ---------------------------------------------------------------------------
# 1c. Clinical-report helpers (issues #176, #182; no trial data needed)
# ---------------------------------------------------------------------------
# Exit reasons of patients without TLR: progression exit first (death after
# progression or progression), then the raw AE exit flag, else no recorded
# reason. Patients with a TLR value are ignored.
toy_cohort <- data.frame(
  ID = paste0("p", 1:6), Rx = factor(c("Control arm", "Control arm", "Control arm",
                                        "Experimental arm", "Control arm", "Control arm")),
  Death = c(1, 1, 0, 1, 1, 0), ProgressionExit = c(1, 0, 0, 0, 1, 0),
  tlr = c(NA, NA, NA, NA, NA, 1))
toy_raw <- data.frame(ID = paste0("p", 6:1), `AE exit` = c(0, 1, 0, 0, 1, 1),
                      check.names = FALSE)
toy_exit <- tlr_missing_exit_summary(toy_cohort, raw = toy_raw)
stopifnot(
  toy_exit$n == 5,
  identical(toy_exit$reasons, c("death after progression" = 2L, "adverse event" = 1L,
                                "no recorded reason" = 2L)),
  identical(toy_exit$reason_text,
            "2 deaths after progression, 1 adverse event, 2 with no recorded reason"),
  identical(toy_exit$arm_text, "4 in the control arm and 1 in the experimental arm"),
  inherits(tryCatch(tlr_missing_exit_summary(toy_cohort, raw = toy_raw[-2, ]),
                    error = function(e) e), "error")
)
# coxphf convergence diagnostics: a converged fit, and a fit that reaches the
# iteration limit is flagged.
set.seed(1)
toy_surv <- data.frame(time = rexp(40), status = rbinom(40, 1, 0.8), x = rbinom(40, 1, 0.5))
toy_fit <- fit_firth_cox(Surv(time, status) ~ x, toy_surv, "toy")
toy_conv <- coxphf_convergence(toy_fit)
stopifnot(toy_conv$Status == "Converged", toy_conv$Maxit == FIRTH_COX_MAXIT,
          toy_conv$N == 40, toy_conv$Events == sum(toy_surv$status),
          toy_conv$Iterations < FIRTH_COX_MAXIT,
          coxphf_convergence(toy_fit, maxit = toy_conv$Iterations)$Status == "Reached maxit",
          coxphf_convergence(NULL)$Status == "Fit failed")

# ---------------------------------------------------------------------------
# 2. A DAG edit that adds an untested edge must stop the report
# ---------------------------------------------------------------------------
dag_extra <- dagitty(sub("PFS -> OS\n}", "PFS -> OS\nSex -> OS\n}", DAG_SPEC, fixed = TRUE))
stopifnot("Sex -> OS" %in% dag_edge_coverage(dag_extra)$testable)
fake_dd <- list(data = data.frame(), data_main = data.frame(), data_tlr = data.frame(),
                data_lm_pfs = data.frame(), data_td = data.frame(),
                landmark = list(spec = tlr_landmark_spec("week9")))
err <- tryCatch(run_dag_edge_tests(dag_extra, fake_dd), error = function(e) conditionMessage(e))
stopifnot(is.character(err), grepl("Sex -> OS", err, fixed = TRUE))

# ---------------------------------------------------------------------------
# Data-dependent checks
# ---------------------------------------------------------------------------
rds_path <- here("data", "tidy", "METIMMOX.rds")
if (!file.exists(rds_path)) {
  cat("Trial data not found; structural checks passed, data checks skipped.\n")
  quit(save = "no", status = 0)
}

invisible(capture.output(suppressMessages({
  source(here("analysis/02_setup_and_global_variables.R"))
  source(here("analysis/03_biomarker_strategies.R"))
})))

stopifnot(
  all(c("CT1wk", "ProgressionExit", "TTPwk", "LastEvalwk", "PFS_rule") %in% names(data)),
  identical(is.na(data$CT1wk), is.na(data$tlr))
)

dd <- prepare_dag_data(data, landmark = "week9")

# 3. Landmark cohorts -------------------------------------------------------
# Issue #184: one base cohort for every TLR landmark analysis, the
# TLR-classified patients within the complete-case cohort. The logistic tests
# into TLR keep every patient with an observed TLR.
base <- tlr_landmark_base_cohort(data)
stopifnot(
  identical(base$ID, dd$data_tlr_landmark$ID),
  setequal(base$ID, data_complete$ID[!is.na(data_complete$tlr)]),
  !anyNA(base$tlr), all(base$ID %in% data_complete$ID),
  nrow(dd$data_tlr) >= nrow(base),
  dd$landmark$diagnostics$n_total == nrow(base),
  all(dd$data_lm_pfs$ID %in% base$ID),
  inherits(tryCatch(tlr_landmark_base_cohort(data[, setdiff(names(data), "tmb_braf")]),
                    error = function(e) e), "error")
)
all_lm <- build_all_tlr_landmark_cohorts(base)
for (k in names(all_lm)) {
  lm <- all_lm[[k]]
  stopifnot(
    all(lm$os$OSwk > lm$os$lm_wk),
    all(lm$pfs$PFSwk > lm$pfs$lm_wk),
    all(lm$os$OSwk_lm > 0), all(lm$pfs$PFSwk_lm > 0),
    lm$diagnostics$n_pfs <= lm$diagnostics$n_os,
    lm$diagnostics$n_pfs_excl_progression_exit + lm$diagnostics$n_pfs_excl_death_pf +
      lm$diagnostics$n_pfs_excl_censored == lm$diagnostics$n_pfs_excluded,
    lm$diagnostics$n_pfs_excl_progression_exit_tlr_neg == lm$diagnostics$n_pfs_excl_progression_exit
  )
}
scan_diag <- all_lm$scan$diagnostics
week9_diag <- all_lm$week9$diagnostics
stopifnot(
  scan_diag$n_first_scan_progressors_kept == 0,
  scan_diag$n_scan_after_landmark_kept == 0,
  scan_diag$n_pfs_excluded == scan_diag$n_first_scan_progressors +
    scan_diag$n_pfs_excl_censored + scan_diag$n_pfs_excl_death_pf,
  week9_diag$n_pfs_excl_censored > 0,
  all_lm$week12$diagnostics$n_first_scan_progressors_kept == 0,
  # The primary week-9 cut is known to retain first-scan progressors and
  # patients scanned after week 9; the reports state these counts.
  week9_diag$n_first_scan_progressors_kept > 0,
  week9_diag$n_scan_after_landmark_kept > 0,
  week9_diag$n_first_scan_progressors_kept < week9_diag$n_first_scan_progressors
)
old_data <- derive_pfs_endpoint(data, Inf, pfs_last_assessment, verbose = FALSE)
old_dd <- prepare_dag_data(old_data, landmark = "week9")
stopifnot(old_dd$landmark$diagnostics$n_pfs_excl_censored == 0,
          nrow(data_complete) == 68,
          nrow(economic_prediction_population(old_data)) == nrow(data_complete),
          sum(data$PFS_rule == "death_censored_at_assessment") > 0)
# Issue #182: the complete-case patients without TLR left the study before the
# first on-treatment CT for three different recorded reasons; they are not all
# early deaths. Issue #176: the week-12 landmark is not after every first scan.
live_exit <- tlr_missing_exit_summary(data_complete)
stopifnot(live_exit$n == nrow(data_complete) - nrow(base),
          identical(live_exit$reason_text,
                    "1 death after progression, 1 adverse event, 1 with no recorded reason"),
          identical(live_exit$arm_text, "all in the control arm"),
          all_lm$week12$diagnostics$n_scan_after_landmark_kept ==
            sum(base$CT1wk > 12 & base$PFSwk > 12),
          max(base$CT1wk) > 12)
cat(sprintf("Landmark cohorts (base n=%d): week9 OS n=%d, PFS n=%d (first-scan progressors kept %d, scanned after %d); scan PFS n=%d; week12 PFS n=%d\n",
            nrow(base), week9_diag$n_os,
            week9_diag$n_pfs, week9_diag$n_first_scan_progressors_kept,
            week9_diag$n_scan_after_landmark_kept, scan_diag$n_pfs,
            all_lm$week12$diagnostics$n_pfs))

# 4. Edge and CI tests ------------------------------------------------------
edges <- run_dag_edge_tests(dag, dd)
cis <- run_dag_ci_tests(dag, dd, seed = 123L)
em <- run_dag_effect_modification_tests(dag, dd)
adj <- run_adjusted_interaction_models(dd)

stopifnot(
  setequal(edges$Item, edge_cov$testable),
  setequal(cis$Item, ci_cov$statement),
  all(is.finite(edges$p_value)),
  all(is.finite(cis$p_value[!cis$Definitional])),
  all(is.na(cis$p_value[cis$Definitional])),
  all(is.na(cis$N[cis$Definitional])),
  all(cis$Conclusion[cis$Definitional] == "Holds by construction"),
  # The TLR interaction statements hold by construction and are no longer
  # tested as CIs (issue #170).
  all(cis$Definitional[cis$Item %in% c("TLR _||_ TxCRP | {CRP, T}",
                                        "TLR _||_ TxTMB | {T, TMB_BRAF}")]),
  # ... their former tests are effect-modification checks of T -> TLR.
  setequal(em$Item, c("TxCRP -> TLR", "TxTMB -> TLR")),
  all(is.finite(em$p_value)), all(grepl("^Firth logistic interaction", em$Test)),
  all(em$N <= nrow(dd$data_tlr)),
  # Issue #184 labels: marginal edge tests, within-stratum arm balance.
  all(grepl("marginal: T, ", edges$Test[grepl("^Tx", edges$Item)], fixed = TRUE)),
  all(grepl("marginal", edges$Test[grepl("Firth", edges$Test)], fixed = TRUE)),
  grepl("within TMB/BRAF-positive patients",
        cis$Conclusion[cis$Item == "CRP _||_ TxTMB | {TMB_BRAF}"], fixed = TRUE),
  # Same label family as CRP _||_ T: "Arm imbalance in realised sample" when
  # p < 0.05, "Compatible with arm balance" otherwise; never "CI contradicted".
  identical(startsWith(cis$Conclusion[cis$Item == "CRP _||_ TxTMB | {TMB_BRAF}"],
                       "Arm imbalance in realised sample"),
            cis$p_value[cis$Item == "CRP _||_ TxTMB | {TMB_BRAF}"] < 0.05),
  !any(grepl("contradicted", cis$Conclusion[cis$Randomization_Check | cis$Stratum_Balance])),
  # Adjusted primary model: both interaction terms for both endpoints on the
  # complete-case cohort.
  nrow(adj) == 4,
  setequal(adj$Edge, c("TxCRP -> OS", "TxCRP -> PFS", "TxTMB -> OS", "TxTMB -> PFS")),
  all(adj$N == nrow(dd$data_main)),
  all(is.finite(adj$HR)), all(adj$CI_Lower < adj$HR & adj$HR < adj$CI_Upper),
  all(is.finite(adj$p_value)),
  grepl("time-dependent", edges$Test[edges$Item == "PFS -> OS"]),
  grepl("landmark", edges$Test[edges$Item == "TLR -> PFS"]),
  edges$N[edges$Item == "TLR -> PFS"] == nrow(dd$data_lm_pfs),
  # Counting-process data: more rows than patients, one row per interval,
  # progression flag switches on only after a recorded progression.
  nrow(dd$data_td) > nrow(dd$data_main),
  all(dd$data_td$tstart < dd$data_td$tstop),
  all(dd$data_td$progressed[dd$data_td$tstart == 0] == 0)
)

# The time-dependent HR must not be the tautological fixed-covariate result
# (HR < 1 per week of PFS time).
hr_td <- as.numeric(sub(" .*", "", edges$Effect[edges$Item == "PFS -> OS"]))
stopifnot(is.finite(hr_td), hr_td > 1)

cat("test_dag_landmark_contracts.R: all checks passed\n")
