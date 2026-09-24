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

suppressPackageStartupMessages({
  library(here); library(dplyr); library(survival); library(dagitty)
})

source(here("R/pfs_endpoint.R"))
source(here("R/dag_helpers.R"))
source(here("R/assoc_tests.R"))
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
  # The two statements relating the interaction nodes given T hold by
  # construction; nothing else is definitional.
  identical(ci_cov$statement[ci_cov$definitional],
            c("TxCRP _||_ TxTMB | {T, TMB_BRAF}", "TxCRP _||_ TxTMB | {CRP, T}")),
  # Arm-balance checks are exactly the marginal statements involving T.
  setequal(ci_cov$statement[ci_cov$arm_balance],
           c("Age _||_ T", "CRP _||_ T", "Sex _||_ T", "T _||_ TMB_BRAF"))
)

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
all_lm <- build_all_tlr_landmark_cohorts(dd$data_tlr)
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
stopifnot(build_tlr_landmark_cohorts(old_dd$data_tlr)$diagnostics$n_pfs_excl_censored == 0,
          nrow(data_complete) == 68,
          nrow(economic_prediction_population(old_data)) == nrow(data_complete),
          sum(data$PFS_rule == "death_censored_at_assessment") > 0)
cat(sprintf("Landmark cohorts: week9 PFS n=%d (first-scan progressors kept %d, scanned after %d); scan PFS n=%d; week12 PFS n=%d\n",
            week9_diag$n_pfs, week9_diag$n_first_scan_progressors_kept,
            week9_diag$n_scan_after_landmark_kept, scan_diag$n_pfs,
            all_lm$week12$diagnostics$n_pfs))

# 4. Edge and CI tests ------------------------------------------------------
edges <- run_dag_edge_tests(dag, dd)
cis <- run_dag_ci_tests(dag, dd, seed = 123L)

stopifnot(
  setequal(edges$Item, edge_cov$testable),
  setequal(cis$Item, ci_cov$statement),
  all(is.finite(edges$p_value)),
  all(is.finite(cis$p_value[!cis$Definitional])),
  all(is.na(cis$p_value[cis$Definitional])),
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
