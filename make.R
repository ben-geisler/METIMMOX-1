#!/usr/bin/env Rscript
# ===============================================================================
# RUN THE WHOLE METIMMOX-1 PIPELINE
# ===============================================================================
# 1. analysis: the numbered scripts 01-12 (08b included). 01 runs in its own R
#    process because it clears the workspace; 02-12 run together in a second
#    process, as the reports expect. 13_save_snapshot.R is the interactive
#    snapshot utility and is not part of the pipeline.
# 2. vignettes: the publication vignettes in reports/vignettes/, one at a time
#    (they share outputs/manifest.csv); they write outputs/figs/ and outputs/tables/.
# 3. reports: every report in reports/ and reports/technical/, PDF then
#    Markdown, moved to outputs/reports/ (R/render.R).
# 4. tests: tests/run_all.R.
#
# Needs the confidential trial export in data/sensitive/ (see README, "Data").
# Scripts 06 and 10 reuse valid caches; 11 and 12 always recompute. From an
# empty cache the analysis takes about 1.5 hours (06 several minutes, 10 about
# 20 minutes, 12 about 40 minutes); the vignettes and reports take about 1 hour.
#
# Usage (from the repository root):
#   Rscript make.R                         # all four steps
#   Rscript make.R --from 09               # analysis from 09 on (02-06 re-sourced as setup)
#   Rscript make.R --no-analysis           # render and test only
#   Rscript make.R --no-vignettes --no-reports --no-tests
#   Rscript make.R --no-analysis --no-vignettes --no-tests --reports CEA,OWSA
#   Rscript make.R --no-analysis --no-reports --no-tests --vignettes figure1,table_5
#   Rscript make.R --dry-run               # list the steps without running them
# ===============================================================================

args <- commandArgs(trailingOnly = TRUE)
flag <- function(x) x %in% args
opt <- function(x) {
  i <- match(x, args)
  if (is.na(i)) return(NULL)
  if (i == length(args)) stop(x, " needs a value.")
  args[i + 1]
}
split_names <- function(x) if (is.null(x)) NULL else trimws(strsplit(x, ",", fixed = TRUE)[[1]])

root <- here::here()
source(file.path(root, "R", "render.R"))
rscript <- file.path(R.home("bin"), "Rscript")

ANALYSIS_SCRIPTS <- c("01_data_prep", "02_setup_and_global_variables",
                      "03_biomarker_strategies", "04_parametric_survival_analysis",
                      "05_basecase_input_parameters", "06_sampling", "07_traces",
                      "08_basecase_analysis", "08b_enriched_population_analysis",
                      "09_DSA", "10_PSA", "11_EVPPIs", "12_scenario_EVPPIs")
SETUP_SCRIPTS <- ANALYSIS_SCRIPTS[2:5]
PIPELINE_PACKAGES <- c("devtools", "readxl", "dplyr", "tableone", "ggplot2", "flexsurv",
                       "survival", "survminer", "gems", "mstate", "tidyverse", "xtable",
                       "darthtools", "dampack", "mvtnorm", "Matrix", "here", "voi")

#' Analysis scripts to run for a --from value
#'
#' Later scripts expect the session that earlier ones build: 09 and 12 need the
#' parameter distributions of 06, and 11 needs the PSA object of 10. A partial
#' run therefore re-sources 02-06 (06 loads the sampling cache) and, from 11 on,
#' 10 (which loads the PSA cache) before its own steps.
#'
#' @param from Script number prefix ("01" to "12", "08b"), or NULL for all.
#' @return List with `run` (scripts run as steps) and `setup` (earlier scripts
#'   sourced first only to rebuild the session).
analysis_plan <- function(from = NULL) {
  if (is.null(from)) return(list(run = ANALYSIS_SCRIPTS, setup = character(0)))
  v_prefix <- sub("_.*$", "", ANALYSIS_SCRIPTS)
  start <- match(from, v_prefix)
  if (is.na(start)) stop("--from must be one of: ", paste(v_prefix, collapse = ", "))
  v_before <- ANALYSIS_SCRIPTS[seq_len(start - 1)]
  v_needed <- c(SETUP_SCRIPTS, "06_sampling",
                if (start > match("10_PSA", ANALYSIS_SCRIPTS)) "10_PSA")
  list(run = ANALYSIS_SCRIPTS[start:length(ANALYSIS_SCRIPTS)],
       setup = intersect(v_before, v_needed))
}

#' Run R code in a fresh Rscript process from the repository root
#'
#' @param code Character vector of R statements.
#' @param step Label used in the error message.
#' @return Invisible NULL; stops on a non-zero exit status.
run_r <- function(code, step) {
  script <- tempfile(fileext = ".R")
  writeLines(c(sprintf("setwd(%s)", deparse(root)), code), script)
  status <- system2(rscript, shQuote(script))
  if (!identical(status, 0L)) stop(step, " failed (exit status ", status, ")")
  invisible(NULL)
}

#' Source a sequence of analysis scripts in one R process
#'
#' @param v_scripts Script names to run as steps, in order.
#' @param v_setup Scripts sourced first, only to rebuild the session.
#' @return Invisible NULL.
run_analysis <- function(v_scripts, v_setup) {
  if ("01_data_prep" %in% v_scripts) {
    run_r("source('analysis/01_data_prep.R')", "01_data_prep")
    v_scripts <- setdiff(v_scripts, "01_data_prep")
  }
  if (!length(v_scripts)) return(invisible(NULL))
  quoted <- function(x) paste0("'", x, "'", collapse = ", ")
  run_r(c(
    sprintf("pacman::p_load(%s)", paste(PIPELINE_PACKAGES, collapse = ", ")),
    # 06 and 07 expect the model functions; scripts 08-12 source them themselves.
    "for (f in c('model_fun', 'calculate_outcomes', 'prediction_functions')) source(file.path('R', paste0(f, '.R')))",
    if (length(v_setup)) sprintf("for (s in c(%s)) source(file.path('analysis', paste0(s, '.R')))", quoted(v_setup)),
    sprintf("for (s in c(%s)) {", quoted(v_scripts)),
    "  cat('\\n== analysis/', s, '.R started ', format(Sys.time()), '\\n', sep = '')",
    "  source(file.path('analysis', paste0(s, '.R')))",
    "  cat('== analysis/', s, '.R finished ', format(Sys.time()), '\\n', sep = '')",
    "}"
  ), "Analysis")
}

plan <- analysis_plan(opt("--from"))
v_vignettes <- vignette_sources(split_names(opt("--vignettes")))
v_reports <- report_sources(split_names(opt("--reports")))
steps <- list(
  analysis = !flag("--no-analysis"),
  vignettes = !flag("--no-vignettes"),
  reports = !flag("--no-reports"),
  tests = !flag("--no-tests")
)

if (flag("--dry-run")) {
  if (steps$analysis) {
    if (length(plan$setup)) cat("Setup (sourced, not rerun as steps):", paste(plan$setup, collapse = ", "), "\n")
    cat("Analysis:", paste(plan$run, collapse = ", "), "\n")
  }
  if (steps$vignettes) cat("Vignettes (", length(v_vignettes), "):\n  ", paste(v_vignettes, collapse = "\n  "), "\n", sep = "")
  if (steps$reports) cat("Reports (", length(v_reports), "):\n  ", paste(v_reports, collapse = "\n  "), "\n", sep = "")
  if (steps$tests) cat("Tests: tests/run_all.R\n")
  quit(save = "no", status = 0)
}

started <- Sys.time()
if (steps$analysis) run_analysis(plan$run, plan$setup)
if (steps$vignettes) for (v in v_vignettes) render_vignette(v)
if (steps$reports) for (r in v_reports) render_report(r)
if (steps$tests) {
  status <- system2(rscript, shQuote(file.path(root, "tests", "run_all.R")))
  if (!identical(status, 0L)) stop("tests/run_all.R failed (exit status ", status, ")")
}
cat(sprintf("make.R finished in %.1f minutes.\n", as.numeric(difftime(Sys.time(), started, units = "mins"))))
