# Output reconciliation for the clinical effectiveness outputs (issue #192).
#
# Reads the GENERATED outputs and fails loudly (exit status 1, every mismatch
# listed) when two of them state the same quantity differently:
#   outputs/reports/clinical_effectiveness.md, outputs/reports/dag_associations.md,
#   outputs/tables/clin_effect_table_s2*.csv, clin_effect_table_s3*.csv.
# Cohort counts and test results are also recomputed from the trial data, so a
# stale render is caught as well. It is not part of tests/run_all.R (it needs
# rendered reports); render the reports and Tables S2/S3 first, then run from
# the repository root:
#   Rscript tests/check_output_consistency.R
#
# Sections: 1 cohort counts, 2 adjusted interaction estimates, 3 marginal edge
# estimates and the S2/S3 tables against dag_associations, 4 week-4 CRP arm
# imbalance (two cohorts, each labelled with its N), 5 N of the tests into TLR.

suppressPackageStartupMessages({
  pacman::p_load(here, dplyr, survival, coxphf, logistf, coin, dagitty, tibble)
})
setwd(here::here())
suppressMessages({
  source("analysis/02_setup_and_global_variables.R")
  source("analysis/03_biomarker_strategies.R")
  for (f in c("assoc_tests", "dag_helpers", "cox_extract", "tlr_landmark",
              "dag_association_tests", "supp_table_labels")) {
    source(file.path("R", paste0(f, ".R")))
  }
})

# Counts the manuscript states (AGENTS.md); a change must be made deliberately.
EXPECTED <- list(n_main = 68L, deaths = 59L, pfs_events = 50L,
                 n_tlr_base = 65L, n_lm_pfs = 56L,
                 crp_all = list(n = 71L, or = "3.51", p = "0.023"),
                 crp_main = list(n = 68L, or = "3.80", p = "0.020"))

v_problems <- character()
n_checks <- 0L
check <- function(ok, ...) {
  n_checks <<- n_checks + 1L
  if (!isTRUE(ok)) v_problems <<- c(v_problems, paste0(...))
  invisible(ok)
}

# ---- Markdown readers -------------------------------------------------------
PIPE <- "\u0001"
read_md <- function(path) {
  if (!file.exists(path)) stop("Missing output (render it first): ", path)
  readLines(path, warn = FALSE, encoding = "UTF-8")
}
md_cells <- function(line) {
  line <- gsub("\\\\\\|", PIPE, line)
  v_cells <- trimws(strsplit(sub("\\|\\s*$", "", sub("^\\s*\\|", "", line)), "|", fixed = TRUE)[[1]])
  v_cells <- gsub(PIPE, "|", v_cells, fixed = TRUE)
  v_cells <- gsub("\\\\([<>])", "\\1", v_cells)
  v_cells <- gsub("*||*", "_||_", v_cells, fixed = TRUE)
  gsub("’", "'", v_cells, fixed = TRUE)
}
# All pipe tables of a file: list of (line, header, rows).
md_tables <- function(lines) {
  v_start <- which(grepl("^\\s*\\|", lines) & !grepl("^\\s*\\|", c("", lines[-length(lines)])))
  lapply(v_start, function(s) {
    e <- s
    while (e + 1 <= length(lines) && grepl("^\\s*\\|", lines[e + 1])) e <- e + 1
    v_header <- md_cells(lines[s])
    v_body <- lines[(s + 2):e]
    v_body <- v_body[nzchar(v_body)]
    m <- do.call(rbind, lapply(v_body, function(l) {
      v <- md_cells(l); length(v) <- length(v_header); v
    }))
    df <- as.data.frame(m, stringsAsFactors = FALSE)
    names(df) <- v_header
    list(line = s, header = v_header, rows = df)
  })
}
# k-th table after the first heading matching `heading`, whose header's first
# cell (or any cell) matches `header_regex`, before the next heading matching `stop`.
md_table_after <- function(lines, tables, heading, header_regex, k = 1, stop = NULL, from = 1) {
  v_h <- grep(heading, lines)
  v_h <- v_h[v_h >= from]
  if (length(v_h) == 0) stop("Heading not found in output: ", heading)
  h <- v_h[1]
  e <- length(lines) + 1
  if (!is.null(stop)) {
    v_s <- grep(stop, lines)
    v_s <- v_s[v_s > h]
    if (length(v_s) > 0) e <- v_s[1]
  }
  v_ok <- Filter(function(t) t$line > h && t$line < e &&
                   any(grepl(header_regex, t$header)), tables)
  if (length(v_ok) < k) stop("Table '", header_regex, "' (#", k, ") not found after '", heading, "'")
  v_ok[[k]]$rows
}
# "0.65 (0.19-2.60)" or "0.65 (95% CI 0.19 to 2.60)" -> "0.65 (0.19, 2.60)"
hr_key <- function(x) {
  m <- regmatches(x, regexec("^\\s*([0-9.]+) \\((?:95% CI )?([0-9.]+)(?: to |-)([0-9.]+)\\)", x))
  vapply(m, function(v) if (length(v) == 4) sprintf("%s (%s, %s)", v[2], v[3], v[4]) else NA_character_,
         character(1))
}
row_of <- function(df, term_regex, col = 1) {
  i <- which(grepl(term_regex, df[[col]]))
  if (length(i) != 1) stop("Expected exactly one row matching '", term_regex, "', found ", length(i))
  df[i, , drop = FALSE]
}
same <- function(a, b, what) {
  check(!is.na(a) && !is.na(b) && identical(as.character(a), as.character(b)),
        what, ": '", a, "' vs '", b, "'")
}

lines_clin <- read_md("outputs/reports/clinical_effectiveness.md")
lines_dag <- read_md("outputs/reports/dag_associations.md")
tab_clin <- md_tables(lines_clin)
tab_dag <- md_tables(lines_dag)
rd <- function(f, ...) read.csv(f, check.names = FALSE, stringsAsFactors = FALSE, encoding = "UTF-8", ...)
t2 <- rd("outputs/tables/clin_effect_table_s2.csv", colClasses = "character")
t3 <- rd("outputs/tables/clin_effect_table_s3.csv", colClasses = "character")
r2 <- rd("outputs/tables/clin_effect_table_s2_raw.csv")
r3 <- rd("outputs/tables/clin_effect_table_s3_raw.csv")

dd <- prepare_dag_data(data, landmark = "week9")
d_main <- dd$data_main

# ---- 1. Cohort counts -------------------------------------------------------
cat("1. Cohort counts\n")
cnt <- function(d) c(N = nrow(d), Deaths = sum(d$Death == 1), PFS_events = sum(d$Progression == 1))
cohorts <- list(
  "Source dataset after scripts 02 and 03" = cnt(dd$data),
  "Main complete-case dataset for survival and interaction models" = cnt(d_main),
  "TLR subset with observed tlr (logistic tests into TLR)" = cnt(dd$data_tlr),
  "TLR landmark base cohort (TLR-classified, complete-case)" = cnt(dd$data_tlr_landmark),
  "TLR landmark PFS cohort (progression-free after week 9)" =
    c(N = nrow(dd$data_lm_pfs), Deaths = sum(dd$data_lm_pfs$Death == 1),
      PFS_events = sum(dd$data_lm_pfs$Progression == 1))
)
check(cohorts[[2]][["N"]] == EXPECTED$n_main, "Complete-case N is ", cohorts[[2]][["N"]], ", expected ", EXPECTED$n_main)
check(cohorts[[2]][["Deaths"]] == EXPECTED$deaths, "Complete-case deaths ", cohorts[[2]][["Deaths"]], ", expected ", EXPECTED$deaths)
check(cohorts[[2]][["PFS_events"]] == EXPECTED$pfs_events, "Complete-case PFS events ", cohorts[[2]][["PFS_events"]], ", expected ", EXPECTED$pfs_events)
check(cohorts[[4]][["N"]] == EXPECTED$n_tlr_base, "TLR-classified N is ", cohorts[[4]][["N"]], ", expected ", EXPECTED$n_tlr_base)
check(cohorts[[5]][["N"]] == EXPECTED$n_lm_pfs, "Week-9 PFS landmark N is ", cohorts[[5]][["N"]], ", expected ", EXPECTED$n_lm_pfs)

# dag_associations: dataset table.
df_ds <- md_table_after(lines_dag, tab_dag, "^# Methods|^## Data and Variables", "^Dataset$")
for (nm in names(cohorts)) {
  r <- row_of(df_ds, paste0("^", gsub("([()])", "\\\\\\1", nm), "$"))
  for (k in c("N", "Deaths", "PFS_events")) {
    same(r[[k]], cohorts[[nm]][[k]], paste0("dag_associations dataset table, ", nm, ", ", k))
  }
}
# clinical_effectiveness: convergence table (N and events of every Firth fit).
df_conv <- md_table_after(lines_clin, tab_clin, "^# Firth Model Convergence", "^Model$")
for (spec in list(c("Primary model, OS", cohorts[[2]][["N"]], cohorts[[2]][["Deaths"]]),
                  c("Primary model, PFS", cohorts[[2]][["N"]], cohorts[[2]][["PFS_events"]]),
                  c("Landmark week 9, OS", cohorts[[4]][["N"]], cohorts[[4]][["Deaths"]]),
                  c("Landmark week 9, PFS", cohorts[[5]][["N"]], cohorts[[5]][["PFS_events"]]))) {
  r <- row_of(df_conv, paste0("^", spec[1], "$"))
  same(r$N, spec[2], paste0("clinical convergence table, ", spec[1], ", N"))
  same(r$Events, spec[3], paste0("clinical convergence table, ", spec[1], ", events"))
}
# Event counts quoted in prose.
txt_clin <- gsub("\\s+", " ", paste(lines_clin, collapse = " "))
for (m in regmatches(txt_clin, gregexpr("[0-9]+ PFS events", txt_clin))[[1]]) {
  check(sub(" PFS events", "", m) %in% as.character(c(EXPECTED$pfs_events, 63, 51)),
        "Unexpected PFS event count in the clinical report text: ", m)
}
for (m in regmatches(txt_clin, gregexpr("[0-9]+ OS events", txt_clin))[[1]]) {
  check(sub(" OS events", "", m) == EXPECTED$deaths, "Unexpected OS event count in the clinical report text: ", m)
}
check(grepl(sprintf("n = %d, %d PFS events", EXPECTED$n_main, EXPECTED$pfs_events), txt_clin) &&
        grepl(sprintf("%d deaths", EXPECTED$deaths), txt_clin),
      "The clinical report no longer states the complete-case cohort (n, PFS events, deaths)")

# ---- 2. Adjusted interaction estimates --------------------------------------
cat("2. Adjusted interaction estimates (clinical report, dag_associations, live refit)\n")
df_uni_os <- md_table_after(lines_clin, tab_clin, "^## Unified Model Results", "^Term$", k = 1,
                            stop = "^## Forest Plot")
df_uni_pfs <- md_table_after(lines_clin, tab_clin, "^## Unified Model Results", "^Term$", k = 2,
                             stop = "^## Forest Plot")
df_dag_adj <- md_table_after(lines_dag, tab_dag, "^## Interaction edges", "^Edge$")
live <- run_adjusted_interaction_models(dd)

edge_info <- list(
  list(edge = "TxCRP -> OS", term = "^CRP x Rx$", tab = df_uni_os),
  list(edge = "TxCRP -> PFS", term = "^CRP x Rx$", tab = df_uni_pfs),
  list(edge = "TxTMB -> OS", term = "^TMB/BRAF x Rx$", tab = df_uni_os),
  list(edge = "TxTMB -> PFS", term = "^TMB/BRAF x Rx$", tab = df_uni_pfs)
)
unified <- list()
for (e in edge_info) {
  rc <- row_of(e$tab, e$term)
  rd_ <- row_of(df_dag_adj, paste0("^", e$edge, "$"))
  rl <- live[live$Edge == e$edge, ]
  k_clin <- hr_key(rc[["HR (95% CI)"]])
  k_dag <- hr_key(rd_[["Adjusted HR (95% CI)"]])
  k_live <- sprintf("%.2f (%.2f, %.2f)", rl$HR, rl$CI_Lower, rl$CI_Upper)
  same(k_clin, k_dag, paste0(e$edge, ": adjusted HR and CI, clinical report vs dag_associations"))
  same(k_clin, k_live, paste0(e$edge, ": adjusted HR and CI, clinical report vs live refit"))
  same(rc[["PLRT p-value"]], rd_[["Adjusted PLRT p"]], paste0(e$edge, ": PLRT p, clinical report vs dag_associations"))
  same(rc[["PLRT p-value"]], fmt_p(rl$p_value), paste0(e$edge, ": PLRT p, clinical report vs live refit"))
  unified[[e$edge]] <- list(key = k_clin, p = rc[["PLRT p-value"]], term = e$term,
                            endpoint = sub(".* -> ", "", e$edge))
}
# The same estimates in the other clinical tables (Firth columns of Tables S7/S8,
# the Firth column of Table 3, and Table S11's primary rule).
df_sf <- list(
  OS = md_table_after(lines_clin, tab_clin, "^## Standard Cox vs Firth", "Standard Cox HR", k = 1),
  PFS = md_table_after(lines_clin, tab_clin, "^## Standard Cox vs Firth", "Standard Cox HR", k = 2))
df_ridge <- list(
  OS = md_table_after(lines_clin, tab_clin, "^# Sensitivity Analysis: Ridge Regression", "^Interaction$", k = 1),
  PFS = md_table_after(lines_clin, tab_clin, "^# Sensitivity Analysis: Ridge Regression", "^Interaction$", k = 2))
df_s11 <- md_table_after(lines_clin, tab_clin, "^## PFS Endpoint Sensitivity", "^Rule$")
r_primary <- row_of(df_s11, "^primary$")
for (nm in names(unified)) {
  u <- unified[[nm]]
  same(hr_key(row_of(df_sf[[u$endpoint]], u$term)[["Firth HR (95% CI)"]]), u$key,
       paste0(nm, ": Firth column of the standard-vs-Firth table (S7/S8) vs the primary table (Table 2)"))
  same(hr_key(row_of(df_ridge[[u$endpoint]], u$term)[["Firth HR (95% CI)"]]), u$key,
       paste0(nm, ": Firth column of the ridge table (Table 3) vs the primary table (Table 2)"))
  if (u$endpoint == "PFS") {
    is_crp <- grepl("CRP", u$term)
    same(hr_key(r_primary[[if (is_crp) "CRP x Rx HR (95% CI)" else "TMB/BRAF x Rx HR (95% CI)"]]), u$key,
         paste0(nm, ": primary rule of the death-censoring table (S11) vs the primary table (Table 2)"))
    same(r_primary[[if (is_crp) 7 else 9]], u$p, paste0(nm, ": PLRT p, S11 vs Table 2"))
  }
}
same(r_primary[["PFS events"]], cohorts[[2]][["PFS_events"]], "Table S11 primary rule: PFS events vs the complete-case cohort")

# Week-9 landmark model (Table 4, S9/S10) against the sensitivity table.
df_lm_os <- md_table_after(lines_clin, tab_clin, "^## Firth-Corrected Cox Model \\(Landmark\\)", "^Term$", k = 1,
                           stop = "^### Progression-Free Survival")
df_lm_pfs <- md_table_after(lines_clin, tab_clin, "^### Progression-Free Survival", "^Term$", k = 1,
                            from = grep("^## Firth-Corrected Cox Model \\(Landmark\\)", lines_clin)[1])
df_lm_sf <- list(
  OS = md_table_after(lines_clin, tab_clin, "^## Standard Cox vs Firth \\(Landmark\\)", "Standard Cox HR", k = 1),
  PFS = md_table_after(lines_clin, tab_clin, "^## Standard Cox vs Firth \\(Landmark\\)", "Standard Cox HR", k = 2))
df_lm_sens <- md_table_after(lines_clin, tab_clin, "^### Landmark Sensitivity", "^Landmark$")
r_w9 <- row_of(df_lm_sens, "^Fixed week 9 \\(primary\\)$")
for (ep in c("OS", "PFS")) {
  r_t4 <- row_of(if (ep == "OS") df_lm_os else df_lm_pfs, "^TLR x Rx$")
  k4 <- hr_key(r_t4[["HR (95% CI)"]])
  same(hr_key(row_of(df_lm_sf[[ep]], "^TLR x Rx$")[["Firth HR (95% CI)"]]), k4,
       paste0("Week-9 TLR x Rx ", ep, ": standard-vs-Firth table (S9/S10) vs the landmark model (Table 4)"))
  same(hr_key(r_w9[[paste0(ep, " TLR x Rx HR (95% CI)")]]), k4,
       paste0("Week-9 TLR x Rx ", ep, ": landmark sensitivity table vs the landmark model (Table 4)"))
  same(r_w9[[if (ep == "OS") 9 else 11]], r_t4[["PLRT p-value"]],
       paste0("Week-9 TLR x Rx ", ep, ": PLRT p, landmark sensitivity table vs Table 4"))
}
same(r_w9[["N (OS)"]], cohorts[[4]][["N"]], "Landmark sensitivity table, week 9: N (OS)")
same(r_w9[["N (PFS)"]], cohorts[[5]][["N"]], "Landmark sensitivity table, week 9: N (PFS)")

# ---- 3. Marginal edges and the S2/S3 tables against dag_associations ---------
cat("3. Tables S2 and S3 against dag_associations and the interaction table\n")
for (edge in names(unified)) {
  rr <- r2[r2$Item == edge, ]
  rd_ <- row_of(df_dag_adj, paste0("^", edge, "$"))
  same(hr_key(rd_[["Marginal HR (95% CI)"]]), hr_key(rr$Effect), paste0(edge, ": marginal HR, dag_associations vs Table S2"))
  same(rd_[["Marginal p"]], fmt_p(rr$p_value), paste0(edge, ": marginal p, dag_associations vs Table S2"))
}
# Every S2/S3 row against the dag_associations tables (N, effect, p).
dag_rows <- list()
for (t in tab_dag) {
  if (!any(t$header == "Effect")) next  # only the tables with N, Effect and p columns
  for (i in seq_len(nrow(t$rows))) {
    v <- unlist(t$rows[i, ], use.names = FALSE)
    j <- which(grepl("^(\\d+|--)$", v[-1]))[1] + 1
    if (!is.na(j) && j + 2 <= length(v)) {
      dag_rows[[length(dag_rows) + 1]] <- list(item = v[1], n = v[j], effect = v[j + 1], p = v[j + 2])
    }
  }
}
dag_item <- vapply(dag_rows, `[[`, "", "item")
for (src in list(list(r2, "S2"), list(r3, "S3"))) {
  r <- src[[1]]
  for (i in seq_len(nrow(r))) {
    k <- which(dag_item == r$Item[i])
    # Definitional statements have no test row; dag_associations names them in prose.
    if (isTRUE(r$Definitional[i])) {
      check(any(grepl(r$Item[i], gsub("\\\\([<>])", "\\1", lines_dag), fixed = TRUE)),
            "Definitional statement '", r$Item[i], "' is not named in dag_associations")
      next
    }
    if (length(k) == 0) { check(FALSE, "Table ", src[[2]], " row '", r$Item[i], "' is not in dag_associations"); next }
    for (kk in k) {
      same(dag_rows[[kk]]$n, ifelse(is.na(r$N[i]), "--", r$N[i]), paste0("Table ", src[[2]], " N, '", r$Item[i], "'"))
      same(dag_rows[[kk]]$p, fmt_p(r$p_value[i]), paste0("Table ", src[[2]], " p, '", r$Item[i], "'"))
      if (grepl("^[0-9.]+ \\(95% CI", r$Effect[i])) {
        same(hr_key(dag_rows[[kk]]$effect), hr_key(r$Effect[i]), paste0("Table ", src[[2]], " effect, '", r$Item[i], "'"))
      }
    }
  }
}
# The formatted tables follow from the raw files.
check(identical(t2$N, as.character(r2$N)) && identical(t2$P, fmt_p(r2$p_value)),
      "Table S2 formatted N/P do not follow from the raw file")
check(identical(t3$P, fmt_p(r3$p_value)), "Table S3 formatted P does not follow from the raw file")
check(nrow(r2) == 18, "Table S2 has ", nrow(r2), " rows, expected 18 direct edges")
check(sum(!r3$Definitional & r3$Section == "Conditional independence") == 15 &&
        sum(r3$Definitional) == 4 && sum(r3$Section != "Conditional independence") == 2,
      "Table S3 does not have 15 testable, 4 definitional and 2 omitted-edge rows")

# ---- 4. Week-4 CRP arm imbalance: two cohorts, each labelled with N ----------
cat("4. Week-4 CRP arm imbalance (n = 71 and n = 68)\n")
f71 <- run_fisher_test(dd$data, "crp", "T_num", "CRP _||_ T", "Fisher's exact")
f68 <- run_fisher_test(d_main, "crp", "T_num", "CRP _||_ T", "Fisher's exact")
for (v in list(list(f71, EXPECTED$crp_all), list(f68, EXPECTED$crp_main))) {
  same(v[[1]]$N, v[[2]]$n, "CRP arm-imbalance test N")
  same(sprintf("%.2f", v[[1]]$Estimate), v[[2]]$or, paste0("CRP arm-imbalance OR, n = ", v[[2]]$n))
  same(fmt_p(v[[1]]$p_value), v[[2]]$p, paste0("CRP arm-imbalance p, n = ", v[[2]]$n))
}
r_crp <- r3[r3$Item == "CRP _||_ T", ]
same(r_crp$N, EXPECTED$crp_all$n, "Table S3 CRP _||_ T: N")
same(sprintf("%.2f", r_crp$Estimate), EXPECTED$crp_all$or, "Table S3 CRP _||_ T: OR")
same(fmt_p(r_crp$p_value), EXPECTED$crp_all$p, "Table S3 CRP _||_ T: p")

# Every statement of either result must carry its N in the same paragraph
# (prose) or row (tables), and both versions must be stated in both reports.
block_of <- function(lines, i) {
  if (grepl("^\\s*\\|", lines[i])) return(lines[i])
  s <- i; while (s > 1 && nzchar(trimws(lines[s - 1]))) s <- s - 1
  e <- i; while (e < length(lines) && nzchar(trimws(lines[e + 1]))) e <- e + 1
  lines[s:e]
}
label_check <- function(lines, report) {
  for (v in list(EXPECTED$crp_all, EXPECTED$crp_main)) {
    pat <- paste0("p = ", v$p, "\\b|\\| ", v$p, " \\||\\b", gsub(".", "\\.", v$or, fixed = TRUE), "\\b")
    n_stated <- 0L
    for (i in grep(pat, lines)) {
      txt <- gsub("\\s+", " ", paste(block_of(lines, i), collapse = " "))
      if (!grepl("CRP", txt)) next
      n_stated <- n_stated + 1L
      check(grepl(sprintf("n = %d\\b|\\| %d \\|", v$n, v$n), txt),
            report, ", line ", i, ": the week-4 CRP arm-imbalance result (p = ", v$p,
            ", OR ", v$or, ") is quoted without its N (n = ", v$n, ")")
    }
    check(n_stated > 0, report, " does not state the week-4 CRP arm-imbalance result for n = ",
          v$n, " (p = ", v$p, ")")
  }
}
label_check(lines_clin, "clinical_effectiveness.md")
label_check(lines_dag, "dag_associations.md")

# ---- 5. N of the tests into TLR ---------------------------------------------
cat("5. Cohorts behind the tests into TLR\n")
d <- dd$data
n_tlr_obs <- sum(!is.na(d$tlr))
n_tlr_tmb <- sum(!is.na(d$tlr) & !is.na(d$tmb_braf))
n_tlr_crp_tmb <- sum(!is.na(d$tlr) & !is.na(d$tmb_braf) & !is.na(d$crp))
cohort_rules <- tibble::tribble(
  ~Item, ~Cohort, ~Expected,
  "CRP -> TLR", "observed TLR and CRP", sum(!is.na(d$tlr) & !is.na(d$crp)),
  "T -> TLR", "observed TLR", n_tlr_obs,
  "TxCRP -> TLR", "observed TLR and CRP", sum(!is.na(d$tlr) & !is.na(d$crp)),
  "TMB_BRAF -> TLR", "observed TLR and TMB/BRAF", n_tlr_tmb,
  "TxTMB -> TLR", "observed TLR and TMB/BRAF", n_tlr_tmb,
  "Age _||_ TLR | {CRP, TMB_BRAF}", "observed TLR, CRP and TMB/BRAF", n_tlr_crp_tmb,
  "Sex _||_ TLR | {CRP, TMB_BRAF}", "observed TLR, CRP and TMB/BRAF", n_tlr_crp_tmb,
  "TLR -> PFS", "week-9 landmark PFS cohort", nrow(dd$data_lm_pfs)
)
all_raw <- rbind(r2[, c("Item", "N")], r3[, c("Item", "N")])
cohort_rules$Table_N <- all_raw$N[match(cohort_rules$Item, all_raw$Item)]
for (i in seq_len(nrow(cohort_rules))) {
  same(cohort_rules$Table_N[i], cohort_rules$Expected[i], paste0("N of '", cohort_rules$Item[i], "' (", cohort_rules$Cohort[i], ")"))
}
# Survival-type tests all use the complete-case cohort.
v_surv <- c(r2$Item[grepl("-> (OS|PFS)$", r2$Item) & r2$Item != "TLR -> PFS"])
for (it in v_surv) same(r2$N[r2$Item == it], EXPECTED$n_main, paste0("N of '", it, "' (complete-case cohort)"))
check(n_tlr_obs == 68L && n_tlr_tmb == 65L && n_tlr_crp_tmb == 65L, "Observed-TLR counts changed (68/65/65): ",
      n_tlr_obs, "/", n_tlr_tmb, "/", n_tlr_crp_tmb)

cat("\nCohort definitions behind the tests into TLR (not a failure; a decision is pending):\n")
print(as.data.frame(cohort_rules), row.names = FALSE)
n_out <- sum(!is.na(d$tlr) & !(d$ID %in% d_main$ID))
n_out_tmb <- sum(!is.na(d$tlr) & !(d$ID %in% d_main$ID) & is.na(d$tmb_braf))
cat(sprintf(paste0("NOTE: %d patients have an observed TLR but are outside the %d-patient complete-case cohort; %d of them lack a TMB/BRAF result.\n",
                   "      The logistic tests use every patient with the variables of the test (%d for CRP and T, %d once TMB/BRAF is needed);\n",
                   "      the clinical report, Figure 1 and the landmark analyses use the %d TLR-classified patients of the complete-case cohort.\n"),
            n_out, nrow(d_main), n_out_tmb, n_tlr_obs, n_tlr_tmb, nrow(dd$data_tlr_landmark)))

# ---- Verdict ------------------------------------------------------------------
cat(sprintf("\n%d checks, %d failed.\n", n_checks, length(v_problems)))
if (length(v_problems) > 0) {
  cat("\nFAIL:\n", paste0(" - ", v_problems, collapse = "\n"), "\n", sep = "")
  quit(status = 1)
}
cat("PASS: the clinical effectiveness outputs agree.\n")
