# ===============================================================================
# TEST: cache contracts and report setup-chunk wiring
# ===============================================================================
# Closes the blind spot that let two defects reach rendered output despite a
# fully passing unit suite. Neither needed a render to detect -- both were
# contract failures:
#
#   1. extract_interaction_coefficients() returned 0% for every draw, behind a
#      warning(). The EVPPI table looked complete but silently omitted all four
#      treatment-by-biomarker rows.  -> Section 1
#   2. scenario_effect.qmd compared failed_draw_policy against a stale literal
#      ("drop_unreplaced_v1") instead of calling psa_failed_draw_policy(), so it
#      rejected every freshly generated cache and could not render at all.
#      -> Section 2
#   3. Fixing (2) required psa_functions in setup_report(funs=). A report that
#      calls a helper it never sources fails only at render time. -> Section 4
#
# Deliberately does NOT render. Rendering all 41 documents takes ~1.5h; these
# checks run in seconds and catch the same failures.
#
# Cache-dependent sections SKIP (not fail) when a cache is absent, so the test is
# usable on a fresh clone before the pipeline has been run.
#
# Run: "C:\Program Files\R\R-4.3.2\bin\x64\Rscript.exe" \
#        tests/test_report_contracts.R
# ===============================================================================

suppressMessages({
  source(here::here("R", "model_configs.R"))
  source(here::here("R", "cache_paths.R"))
  source(here::here("R", "psa_functions.R"))
})

utility_source_label <- "correct"   # cache_paths resolves label from this
failures <- character(0)
skips <- character(0)
checks <- 0L

check <- function(ok, msg) {
  checks <<- checks + 1L
  if (!isTRUE(ok)) failures <<- c(failures, msg)
}
skip <- function(msg) skips <<- c(skips, msg)

# ---------------------------------------------------------------------------
# 1. EVPPI cache carries the treatment-by-biomarker interaction rows
# ---------------------------------------------------------------------------
evppi_file <- evppi_path()
if (!file.exists(evppi_file)) {
  skip(paste("EVPPI cache absent:", basename(evppi_file)))
} else {
  e <- new.env()
  load(evppi_file, envir = e)
  check(!is.null(e$evppi_results), "EVPPI cache has no evppi_results object")

  if (!is.null(e$evppi_results) && nrow(e$evppi_results) == 0 &&
      isTRUE(e$evpi_manual <= 0.01)) {
    # A zero EVPI (control optimal in every PSA draw) makes every EVPPI zero
    # and run_evppi_analysis() returns no rows by design; the interaction-row
    # and group checks below have nothing to test. This is the state of the
    # cache after the issue #156 sampling change at WTP 51,000.
    skip(sprintf("EVPPI cache has no rows because EVPI = %s; row checks not applicable",
                 format(e$evpi_manual)))
  } else if (!is.null(e$evppi_results)) {
    params <- e$evppi_results$parameter
    expected <- unlist(lapply(get_biomarkers(), function(b) {
      paste0("b_", b, "_rx_", c("os", "pfs"))
    }))
    missing <- setdiff(expected, params)
    check(
      length(missing) == 0,
      paste0("EVPPI cache is missing interaction rows: ",
             paste(missing, collapse = ", "),
             ". extract_interaction_coefficients() is matching nothing -- ",
             "check the coefficient-name pattern against a real fit.")
    )
    # Since issue #152 a zero EVPPI is a legitimate result for a parameter that
    # never changes the optimal decision, so a zero value is not a failure
    # signature. An extraction failure now surfaces as an NA estimate with a
    # populated error column and no standard error.
    present <- intersect(expected, params)
    if (length(present) > 0) {
      rows <- e$evppi_results[match(present, params), , drop = FALSE]
      check(
        all(is.finite(rows$evppi)) &&
          (!"error" %in% names(rows) || all(!nzchar(rows$error))),
        "Interaction EVPPI rows are NA or carry an error; extraction likely failed"
      )
      check(
        "evppi_se" %in% names(rows) && all(is.finite(rows$evppi_se)),
        "Interaction EVPPI rows lack a finite Monte Carlo standard error"
      )
    }

    # Every group row must be at least its largest member within Monte Carlo
    # tolerance (issue #152). Rebuild the group definitions the pipeline uses.
    suppressMessages({
      source(here::here("R", "parameter_distributions.R"))
      source(here::here("R", "evppi_functions.R"))
    })
    groups <- add_interaction_param_groups(create_parameter_groups())
    consistency <- check_evppi_group_consistency(e$evppi_results, groups)
    check(nrow(consistency) > 0, "EVPPI cache has no [GROUP] rows to check")
    check(
      !any(consistency$violation),
      paste0("Group EVPPI below its largest member beyond tolerance: ",
             paste(consistency$group[consistency$violation], collapse = ", "))
    )
    check(
      "[GROUP] interaction_all" %in% params,
      "EVPPI cache lacks the joint interaction group used by Figure 5 / Table S7"
    )
  }
}

# ---------------------------------------------------------------------------
# 2. Cache guards that reports enforce must pass against the current caches
# ---------------------------------------------------------------------------
scenario_file <- scenario_evppi_path()
if (!file.exists(scenario_file)) {
  skip(paste("Scenario cache absent:", basename(scenario_file)))
} else {
  sc <- readRDS(scenario_file)
  check(
    identical(sc$failed_draw_policy, psa_failed_draw_policy()),
    paste0("Scenario cache failed_draw_policy is '",
           if (is.null(sc$failed_draw_policy)) "none" else sc$failed_draw_policy,
           "' but psa_failed_draw_policy() is '", psa_failed_draw_policy(),
           "'. scenario_effect.qmd will refuse to render.")
  )
}

psa_file <- psa_obj_path()
if (!file.exists(psa_file)) {
  skip(paste("PSA cache absent:", basename(psa_file)))
} else {
  po <- readRDS(psa_file)
  check(
    identical(po$failed_draw_policy, psa_failed_draw_policy()),
    "PSA cache failed_draw_policy does not match psa_failed_draw_policy()"
  )
  check(!anyNA(po$cost) && !anyNA(po$effectiveness),
        "PSA cache contains NA cost/effectiveness entries")
  # Issue #156: per-row model index and input fingerprint.
  check(!is.null(po$model_idx) && length(po$model_idx) == po$n_sim,
        "PSA cache lacks a model_idx vector aligned with its rows (issue #156)")
  check(is.character(po$fingerprint) && nchar(po$fingerprint) > 0,
        "PSA cache lacks an input fingerprint (issue #156)")
  params_file <- psa_params_path()
  if (file.exists(params_file)) {
    pp <- readRDS(params_file)
    check("model_idx" %in% names(pp) && identical(pp$model_idx, po$model_idx),
          "psa_params model_idx column missing or not aligned with the PSA object")
  }
  evppi_file <- evppi_path()
  if (file.exists(evppi_file)) {
    ev_env <- new.env()
    load(evppi_file, envir = ev_env)
    check(identical(ev_env$evppi_psa_fingerprint, po$fingerprint),
          "EVPPI cache was built from a PSA with a different fingerprint (issue #156)")
  }
}

sampling_file <- sampling_cache_path(5000)
if (!file.exists(sampling_file)) {
  skip(paste("Sampling cache absent:", basename(sampling_file)))
} else {
  sm <- readRDS(sampling_file)
  check(!is.null(sm$joint) && is.null(sm$crp),
        "Sampling cache should hold one joint component, not per-biomarker copies (issue #156)")
  check(is.character(sm$fingerprint) && nchar(sm$fingerprint) > 0,
        "Sampling cache lacks an input fingerprint (issue #156)")
  check(identical(sm$joint$method, "mvn_v1") && !is.null(sm$joint$draws$os),
        "Sampling cache is not a multivariate-normal coefficient-draw cache (issue #156)")
  if (exists("po") && !is.null(po$fingerprint_inputs)) {
    check(identical(po$fingerprint_inputs$sampling_fingerprint, sm$fingerprint),
          "PSA cache was built from a different sampling cache (issue #156)")
  }
  if (exists("sc") && !is.null(sc$sampling_fingerprint)) {
    check(identical(sc$sampling_fingerprint, sm$fingerprint),
          "Scenario cache was built from a different sampling cache (issue #156)")
  }
}

# ---------------------------------------------------------------------------
# 3. Every setup_report(sources = ...) prefix resolves to an analysis script
#    Guards against the 06->04 style renumbering silently repointing a report.
# ---------------------------------------------------------------------------
qmd_files <- c(
  list.files(here::here("reports"), pattern = "\\.qmd$",
             recursive = TRUE, full.names = TRUE),
  list.files(here::here("outputs", "vignettes"), pattern = "\\.qmd$",
             full.names = TRUE)
)
analysis_dir <- here::here("analysis")
fun_dir <- here::here("R")

setup_chunk <- function(path) {
  lines <- readLines(path, warn = FALSE)
  starts <- grep("^```\\{r[ ,].*setup", lines)
  if (length(starts) == 0) return(character(0))
  start <- starts[1]
  ends <- grep("^```\\s*$", lines)
  ends <- ends[ends > start]
  if (length(ends) == 0) return(character(0))
  lines[(start + 1):(ends[1] - 1)]
}

arg_values <- function(chunk, arg) {
  txt <- paste(chunk, collapse = " ")
  m <- regmatches(txt, regexpr(paste0(arg, "\\s*=\\s*c\\([^)]*\\)"), txt))
  if (length(m) == 0) {
    m <- regmatches(txt, regexpr(paste0(arg, '\\s*=\\s*"[^"]*"'), txt))
  }
  if (length(m) == 0) return(character(0))
  unlist(regmatches(m, gregexpr('"[^"]*"', m))) |> gsub(pattern = '"', replacement = "")
}

for (qmd in qmd_files) {
  chunk <- setup_chunk(qmd)
  if (length(chunk) == 0) next
  rel <- sub(".*/(reports|outputs)/", "", gsub("\\\\", "/", qmd))

  for (pfx in arg_values(chunk, "sources")) {
    hits <- list.files(analysis_dir, pattern = paste0("^", pfx, "_.*\\.R$"))
    check(length(hits) > 0,
          paste0(rel, ": setup_report(sources=) prefix '", pfx,
                 "' matches no script in analysis/"))
  }
  for (fn in arg_values(chunk, "funs")) {
    check(file.exists(file.path(fun_dir, paste0(fn, ".R"))),
          paste0(rel, ": setup_report(funs=) names '", fn,
                 "' but R/", fn, ".R does not exist"))
  }
}

# ---------------------------------------------------------------------------
# 4. Helpers called in a setup chunk must actually be reachable from it
#    This is the class that made the scenario_effect fix fail on first attempt:
#    the report called psa_failed_draw_policy() without sourcing psa_functions.
# ---------------------------------------------------------------------------
fun_files <- list.files(fun_dir, pattern = "\\.R$", full.names = TRUE)
defines <- list()
for (f in fun_files) {
  src <- readLines(f, warn = FALSE)
  nm <- unique(sub("^\\s*([A-Za-z._][A-Za-z0-9._]*)\\s*(<-|=)\\s*function.*$", "\\1",
                   grep("^\\s*[A-Za-z._][A-Za-z0-9._]*\\s*(<-|=)\\s*function", src, value = TRUE)))
  for (n in nm) defines[[n]] <- c(defines[[n]], tools::file_path_sans_ext(basename(f)))
}
# Sourced unconditionally by 02_setup_and_global_variables.R, so always in scope.
always <- c("model_configs", "cache_paths", "report_setup", "report_tables")

# A declared analysis script brings its own source() calls with it. e.g. 04
# sources prediction_functions.R, so a report declaring sources = "04" may call
# generate_population_averaged_predictions() without naming it in funs=.
# Resolving this transitively is what separates a real wiring gap from a
# false positive.
sourced_by_analysis <- function(prefix) {
  hits <- list.files(analysis_dir, pattern = paste0("^", prefix, "_.*\\.R$"),
                     full.names = TRUE)
  if (length(hits) == 0) return(character(0))
  src <- readLines(hits[1], warn = FALSE)
  refs <- unlist(regmatches(
    src, gregexpr('"?R"?[/", ]+[A-Za-z0-9._]+\\.R', src)
  ))
  unique(tools::file_path_sans_ext(basename(gsub('[", ]+', "/", refs))))
}

for (qmd in qmd_files) {
  chunk <- setup_chunk(qmd)
  if (length(chunk) == 0) next
  rel <- sub(".*/(reports|outputs)/", "", gsub("\\\\", "/", qmd))
  available <- c(always, arg_values(chunk, "funs"))
  for (pfx in arg_values(chunk, "sources")) {
    available <- c(available, sourced_by_analysis(pfx))
  }
  available <- unique(available)

  called <- unique(unlist(regmatches(
    chunk, gregexpr("[A-Za-z._][A-Za-z0-9._]*(?=\\s*\\()", chunk, perl = TRUE)
  )))
  for (fname in called) {
    owners <- defines[[fname]]
    if (is.null(owners)) next          # not one of ours; base/pkg function
    if (!any(owners %in% available)) {
      check(FALSE,
            paste0(rel, ": setup chunk calls ", fname, "() from ",
                   paste(unique(owners), collapse = "/"),
                   " but does not source it via setup_report(funs=)"))
    }
  }
}

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
for (s in skips) cat("SKIP:", s, "\n")
if (length(failures) > 0) {
  cat("\n", length(failures), " CONTRACT FAILURE(S):\n", sep = "")
  for (f in failures) cat("  - ", f, "\n", sep = "")
  stop("Report/cache contract check failed.")
}
cat("PASS: ", checks, " report/cache contract checks (",
    length(skips), " skipped).\n", sep = "")
