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
v_failures <- character(0)
v_skips <- character(0)
checks <- 0L

#' Count a contract check and record `msg` as a failure unless `ok` is TRUE.
check <- function(ok, msg) {
  checks <<- checks + 1L
  if (!isTRUE(ok)) v_failures <<- c(v_failures, msg)
}
#' Record `msg` as a skipped (not applicable) check.
skip <- function(msg) v_skips <<- c(v_skips, msg)

# ---------------------------------------------------------------------------
# 1. EVPPI cache carries the treatment-by-biomarker interaction rows
# ---------------------------------------------------------------------------
evppi_file <- evppi_path()
if (!file.exists(evppi_file)) {
  skip(paste("EVPPI cache absent:", basename(evppi_file)))
} else {
  e <- new.env()
  load(evppi_file, envir = e)
  check(!is.null(e$evppi_fingerprint) && !is.null(e$evppi_fingerprint_inputs),
        "EVPPI cache predates independent input identity (#167); regenerate script 11")
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
    v_params <- e$evppi_results$parameter
    v_expected <- unlist(lapply(get_biomarkers(), function(b) {
      paste0("b_", b, "_rx_", c("os", "pfs"))
    }))
    v_missing <- setdiff(v_expected, v_params)
    check(
      length(v_missing) == 0,
      paste0("EVPPI cache is missing interaction rows: ",
             paste(v_missing, collapse = ", "),
             ". extract_interaction_coefficients() is matching nothing -- ",
             "check the coefficient-name pattern against a real fit.")
    )
    # Since issue #152 a zero EVPPI is a legitimate result for a parameter that
    # never changes the optimal decision, so a zero value is not a failure
    # signature. An extraction failure now surfaces as an NA estimate with a
    # populated error column and no standard error.
    v_present <- intersect(v_expected, v_params)
    if (length(v_present) > 0) {
      df_rows <- e$evppi_results[match(v_present, v_params), , drop = FALSE]
      check(
        all(is.finite(df_rows$evppi)) &&
          (!"error" %in% names(df_rows) || all(!nzchar(df_rows$error))),
        "Interaction EVPPI rows are NA or carry an error; extraction likely failed"
      )
      check(
        "evppi_se" %in% names(df_rows) && all(is.finite(df_rows$evppi_se)),
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
    df_consistency <- check_evppi_group_consistency(e$evppi_results, groups)
    check(nrow(df_consistency) > 0, "EVPPI cache has no [GROUP] rows to check")
    check(
      !any(df_consistency$violation),
      paste0("Group EVPPI below its largest member beyond tolerance: ",
             paste(df_consistency$group[df_consistency$violation], collapse = ", "))
    )
    check(
      "[GROUP] interaction_all" %in% v_params,
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
  check(!is.null(sc$fingerprint), "Scenario cache lacks its own input identity (#167)")
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
    pair_valid <- tryCatch({ validate_psa_pair(po, pp); TRUE }, error = function(e) FALSE)
    check(pair_valid, "PSA files lack valid generation/content binding (#167); regenerate script 10")
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
  check(identical(sm$fingerprint_inputs$schema, "sampling_v2") &&
          !is.null(sm$canonical_data) && !is.null(sm$runtime),
        "Sampling cache lacks canonical inputs/runtime diagnostics (#169); regenerate script 06")
  check(!is.null(sm$joint) && is.null(sm$crp),
        "Sampling cache should hold one joint component, not per-biomarker copies (issue #156)")
  check(is.character(sm$fingerprint) && nchar(sm$fingerprint) > 0,
        "Sampling cache lacks an input fingerprint (issue #156)")
  check(identical(sm$joint$method, "mvn_joint_v2") && !is.null(sm$joint$draws$os) &&
          !is.null(sm$joint$joint_covariance$covariance),
        "Sampling cache lacks joint OS/PFS coefficient covariance (issue #159)")
  if (!is.null(sm$joint$joint_covariance$covariance)) {
    joint <- sm$joint
    m_covariance <- joint$joint_covariance$covariance
    v_os_idx <- seq_along(joint$original_os$opt$par)
    v_pfs_idx <- length(v_os_idx) + seq_along(joint$original_pfs$opt$par)
    check(max(abs(m_covariance[v_os_idx, v_os_idx] - joint$original_os$cov)) < 1e-10 &&
            max(abs(m_covariance[v_pfs_idx, v_pfs_idx] - joint$original_pfs$cov)) < 1e-10,
          "Joint sampling changed the fitted endpoint marginal covariances (#159)")
    m_realised <- cov(cbind(joint$draws$os[, joint$original_os$optpars],
                          joint$draws$pfs[, joint$original_pfs$optpars]))
    m_mcse <- sqrt((outer(diag(m_covariance), diag(m_covariance)) + m_covariance^2) /
                   (joint$n_samples - 1))
    check(all(abs(m_realised - m_covariance) <= 5 * m_mcse),
          "Realised coefficient covariance differs from joint target by >5 MCSE (#159)")
    check(joint$joint_covariance$n_successful + joint$joint_covariance$n_failed ==
            joint$joint_covariance$n_attempted,
          "Covariance bootstrap convergence counts do not reconcile (#159)")
  }
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
# Reports, technical reports and publication vignettes (reports/vignettes/).
v_qmd_files <- list.files(here::here("reports"), pattern = "\\.qmd$",
                          recursive = TRUE, full.names = TRUE)
analysis_dir <- here::here("analysis")
fun_dir <- here::here("R")

#' Lines of the first setup chunk of the Quarto file at `path`, without its
#' fences; character(0) when there is none.
setup_chunk <- function(path) {
  v_lines <- readLines(path, warn = FALSE)
  v_starts <- grep("^```\\{r[ ,].*setup", v_lines)
  if (length(v_starts) == 0) return(character(0))
  start <- v_starts[1]
  v_ends <- grep("^```\\s*$", v_lines)
  v_ends <- v_ends[v_ends > start]
  if (length(v_ends) == 0) return(character(0))
  v_lines[(start + 1):(v_ends[1] - 1)]
}

#' String values passed to argument `arg` (a c(...) vector or a single string)
#' in a setup chunk; returns a character vector.
arg_values <- function(chunk, arg) {
  txt <- paste(chunk, collapse = " ")
  m <- regmatches(txt, regexpr(paste0(arg, "\\s*=\\s*c\\([^)]*\\)"), txt))
  if (length(m) == 0) {
    m <- regmatches(txt, regexpr(paste0(arg, '\\s*=\\s*"[^"]*"'), txt))
  }
  if (length(m) == 0) return(character(0))
  unlist(regmatches(m, gregexpr('"[^"]*"', m))) |> gsub(pattern = '"', replacement = "")
}

for (qmd in v_qmd_files) {
  v_chunk <- setup_chunk(qmd)
  if (length(v_chunk) == 0) next
  rel <- sub(".*/(reports|outputs)/", "", gsub("\\\\", "/", qmd))

  for (pfx in arg_values(v_chunk, "sources")) {
    v_hits <- list.files(analysis_dir, pattern = paste0("^", pfx, "_.*\\.R$"))
    check(length(v_hits) > 0,
          paste0(rel, ": setup_report(sources=) prefix '", pfx,
                 "' matches no script in analysis/"))
  }
  for (fn in arg_values(v_chunk, "funs")) {
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
v_fun_files <- list.files(fun_dir, pattern = "\\.R$", full.names = TRUE)
defines <- list()
for (f in v_fun_files) {
  v_src <- readLines(f, warn = FALSE)
  nm <- unique(sub("^\\s*([A-Za-z._][A-Za-z0-9._]*)\\s*(<-|=)\\s*function.*$", "\\1",
                   grep("^\\s*[A-Za-z._][A-Za-z0-9._]*\\s*(<-|=)\\s*function", v_src, value = TRUE)))
  for (n in nm) defines[[n]] <- c(defines[[n]], tools::file_path_sans_ext(basename(f)))
}
# Sourced unconditionally by 02_setup_and_global_variables.R, so always in scope.
v_always <- c("model_configs", "cache_paths", "cache_provenance", "report_setup", "report_tables")

# A declared analysis script brings its own source() calls with it. e.g. 04
# sources prediction_functions.R, so a report declaring sources = "04" may call
# generate_population_averaged_predictions() without naming it in funs=.
# Resolving this transitively is what separates a real wiring gap from a
# false positive.
#' Names (without .R) of the R/ files sourced by the analysis script whose name
#' starts with `prefix`.
sourced_by_analysis <- function(prefix) {
  v_hits <- list.files(analysis_dir, pattern = paste0("^", prefix, "_.*\\.R$"),
                     full.names = TRUE)
  if (length(v_hits) == 0) return(character(0))
  v_src <- readLines(v_hits[1], warn = FALSE)
  v_refs <- unlist(regmatches(
    v_src, gregexpr('"?R"?[/", ]+[A-Za-z0-9._]+\\.R', v_src)
  ))
  unique(tools::file_path_sans_ext(basename(gsub('[", ]+', "/", v_refs))))
}

for (qmd in v_qmd_files) {
  v_chunk <- setup_chunk(qmd)
  if (length(v_chunk) == 0) next
  rel <- sub(".*/(reports|outputs)/", "", gsub("\\\\", "/", qmd))
  v_available <- c(v_always, arg_values(v_chunk, "funs"))
  if (any(grepl("R/publication_artifacts.R", v_chunk, fixed = TRUE)))
    v_available <- c(v_available, "publication_artifacts")
  for (pfx in arg_values(v_chunk, "sources")) {
    v_available <- c(v_available, sourced_by_analysis(pfx))
  }
  v_available <- unique(v_available)

  v_called <- unique(unlist(regmatches(
    v_chunk, gregexpr("[A-Za-z._][A-Za-z0-9._]*(?=\\s*\\()", v_chunk, perl = TRUE)
  )))
  for (fname in v_called) {
    v_owners <- defines[[fname]]
    if (is.null(v_owners)) next          # not one of ours; base/pkg function
    if (!any(v_owners %in% v_available)) {
      check(FALSE,
            paste0(rel, ": setup chunk calls ", fname, "() from ",
                   paste(unique(v_owners), collapse = "/"),
                   " but does not source it via setup_report(funs=)"))
    }
  }
}

# Shared cache loaders are required even in load-data chunks outside setup.
for (file in c("figure2", "figure3", "figure5", "table_5", "table_s7")) {
  src <- paste(readLines(here::here("reports", "vignettes", paste0(file, ".qmd")), warn = FALSE), collapse = "\n")
  check(grepl("load_current_psa_cache|load_scenario_cache", src), paste(file, "bypasses validated cache loading"))
  check(!grepl("readRDS", src, fixed = TRUE), paste(file, "still reads a cache directly"))
}

# ---------------------------------------------------------------------------
# Publication outputs (#168): cached PSA -> Table 5, and manifest -> artifacts.
# ---------------------------------------------------------------------------
source(here::here("R/publication_artifacts.R"))
source(here::here("R/report_format.R"))
table5_path <- here::here("outputs/tables/table_5.csv")
if (exists("po") && file.exists(table5_path)) {
  df_table5 <- read.csv(table5_path, colClasses = "character", check.names = FALSE)
  wtp <- 51000
  m_nmb <- as.matrix(po$effect) * wtp - as.matrix(po$cost)
  # Independent calculation, including dampack's equal sharing of ties.
  m_winners <- m_nmb == apply(m_nmb, 1, max)
  v_probabilities <- colMeans(m_winners / rowSums(m_winners))
  check(identical(df_table5$Strategy, as.character(po$strategies)),
        "Table 5 strategy rows differ from the PSA cache")
  control <- match(get_control_strategy(), po$strategies)
  m_cost <- as.matrix(po$cost)
  m_effect <- as.matrix(po$effect)
  m_dc <- sweep(m_cost, 1, m_cost[, control])
  m_de <- sweep(m_effect, 1, m_effect[, control])
  expected <- list(Mean_Cost = colMeans(m_cost), Mean_QALY = colMeans(m_effect),
    Mean_Inc_Cost = colMeans(m_dc), Mean_Inc_QALY = colMeans(m_de), Prob_CE = v_probabilities)
  for (entry in list(list("Cost", m_cost), list("QALY", m_effect),
                     list("Inc_Cost", m_dc), list("Inc_QALY", m_de))) {
    expected[[paste0(entry[[1]], "_Lower")]] <- apply(entry[[2]], 2, quantile, 0.025)
    expected[[paste0(entry[[1]], "_Upper")]] <- apply(entry[[2]], 2, quantile, 0.975)
  }
  expected$ICER <- colMeans(m_dc) / colMeans(m_de)
  expected$ICER[control] <- NA_real_
  for (column in names(expected)) {
    check(isTRUE(all.equal(as.numeric(df_table5[[column]]), unname(expected[[column]]),
                           tolerance = 1e-10)),
          paste("Table 5 differs from independent PSA calculation:", column))
  }
  v_expected_status <- ifelse(colMeans(m_dc) > 0 & colMeans(m_de) <= 0,
    "Dominated by SoC", "Pairwise ICER vs SoC")
  v_expected_status[control] <- "Reference"
  check(identical(df_table5$Status, unname(v_expected_status)),
        "Table 5 pairwise status differs from PSA means")
} else skip("Table 5 versus PSA comparison requires the table and PSA cache")

registry <- publication_outputs()
# Every declared vignette exists and runs the artifact lifecycle (a source
# check, so it also runs in a tree without rendered outputs).
for (vignette in names(registry)) {
  vignette_path <- here::here("reports", "vignettes", paste0(vignette, ".qmd"))
  check(file.exists(vignette_path), paste("Missing vignette for a declared output:", vignette))
  if (file.exists(vignette_path)) {
    src <- paste(readLines(vignette_path, warn = FALSE), collapse = "\n")
    check(grepl(paste0('begin_artifact_render("', vignette, '")'), src, fixed = TRUE) &&
            grepl("finish_artifact_render()", src, fixed = TRUE),
          paste("Vignette lacks an artifact lifecycle:", vignette))
  }
}

# The manifest and the artifacts are render products. A tree without any (a
# fresh clone of the public repository) skips these checks; a tree with
# publication CSVs but no manifest fails.
manifest_path <- here::here("outputs/manifest.csv")
v_output_csvs <- file.path("outputs/tables",
                           list.files(here::here("outputs/tables"), pattern = "\\.csv$"))
if (!file.exists(manifest_path) && length(v_output_csvs) == 0L) {
  skip("No rendered publication outputs in this tree (manifest and artifact checks not applicable)")
} else check(file.exists(manifest_path), "Publication artifact manifest is absent")
if (file.exists(manifest_path)) {
  df_manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
  v_owners <- setNames(rep(paste0("reports/vignettes/", names(registry), ".qmd"), lengths(registry)),
                     unlist(registry, use.names = FALSE))
  check(!anyDuplicated(df_manifest$file), "Manifest has duplicate artifact rows")
  check(setequal(df_manifest$file, unlist(registry)), "Manifest does not cover every declared output")
  check(all(v_output_csvs %in% df_manifest$file), "A publication CSV is missing from the manifest")
  for (i in seq_len(nrow(df_manifest))) {
    df_row <- df_manifest[i, ]
    path <- here::here(df_row$file)
    check(file.exists(here::here(df_row$generating_vignette)), paste("Missing generator:", df_row$file))
    check(identical(df_row$generating_vignette, unname(v_owners[df_row$file])),
          paste("Manifest names the wrong generator:", df_row$file))
    check(nzchar(df_row$render_time), paste("Missing render time:", df_row$file))
    if (df_row$status == "deleted_empty") {
      check(!file.exists(path), paste("Obsolete empty-state predecessor survives:", df_row$file))
    } else {
      check(file.exists(path), paste("Missing publication artifact:", df_row$file))
      check(file.exists(path) && identical(unname(tools::md5sum(path)), df_row$artifact_md5),
            paste("Artifact differs from its manifest:", df_row$file))
      if (grepl("\\.csv$", df_row$file))
        check(df_row$status == "csv_verified", paste("CSV not checked against its generating object:", df_row$file))
    }
  }
  df_table5_manifest <- subset(df_manifest, file == "outputs/tables/table_5.csv")
  if (exists("po")) check(nrow(df_table5_manifest) == 1L &&
    grepl(paste0("psa_obj=", po$fingerprint), df_table5_manifest$source_cache_fingerprints, fixed = TRUE),
    "Table 5 manifest refers to a different PSA input fingerprint")
  if (exists("sc")) {
    for (file in c("outputs/figs/figure3.png", "outputs/figs/figure_s4.png",
                   "outputs/figs/figure5.png", "outputs/tables/table_s7.csv")) {
      df_entry <- df_manifest[df_manifest$file == file, ]
      check(nrow(df_entry) == 1L && grepl(sc$fingerprint, df_entry$source_cache_fingerprints, fixed = TRUE),
            paste("Manifest refers to a different scenario input fingerprint:", file))
    }
  }
}
if (exists("sc") && file.exists(here::here("outputs/tables/table_s7.csv"))) {
  df_s7 <- read.csv(here::here("outputs/tables/table_s7.csv"))
  df_totals <- subset(df_s7, group == "evpi")
  check(setequal(df_totals$scenario_id, sc$scenarios$scenario_id), "Table S7 drops a scenario EVPI row")
  for (result in sc$all_scenario_results) {
    id <- result$scenario_info$scenario_id
    m_nmb <- as.matrix(result$psa_obj$effect) * result$scenario_info$wtp - as.matrix(result$psa_obj$cost)
    evpi <- mean(apply(m_nmb, 1, max)) - max(colMeans(m_nmb))
    check(isTRUE(all.equal(df_totals$evpi[df_totals$scenario_id == id], evpi, tolerance = 1e-10)),
          paste("Table S7 total EVPI differs from scenario PSA:", id))
    if (evpi == 0) {
      df_rows <- subset(df_s7, scenario_id == id & group != "evpi")
      check(nrow(df_rows) > 0 && all(df_rows$evppi == 0) && all(df_rows$evppi_se == 0) &&
        all(is.na(df_rows$evppi_percent_of_evpi)), paste("Table S7 loses explicit zero groups:", id))
    }
  }
}

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
for (s in v_skips) cat("SKIP:", s, "\n")
if (length(v_failures) > 0) {
  cat("\n", length(v_failures), " CONTRACT FAILURE(S):\n", sep = "")
  for (f in v_failures) cat("  - ", f, "\n", sep = "")
  stop("Report/cache contract check failed.")
}
cat("PASS: ", checks, " report/cache contract checks (",
    length(v_skips), " skipped).\n", sep = "")
