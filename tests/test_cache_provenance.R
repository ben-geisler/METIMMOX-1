# Synthetic adversarial cache contracts: no production cache is read or written.
source("R/cache_paths.R")
source("R/cea_helpers.R")

#' Run every synthetic cache-provenance contract in a temporary cache directory
#' (restored on exit); stops at the first failed assertion, returns nothing.
run_cache_provenance_tests <- function() {
  directory <- tempfile("cache-contracts-")
  dir.create(directory)
  old_options <- options(metimmox.cache_dir = directory,
                         metimmox.sampling_allow_regenerate = FALSE)
  on.exit(options(old_options), add = TRUE)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  fails <- function(expr, pattern) {
    error <- tryCatch({ force(expr); NULL }, error = identity)
    stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
  }

  # Stable value hashing: lazy/materialised integers and numbers, UTF-8 strings,
  # column order and irrelevant columns. Factor coding still changes the fit.
  d <- data.frame(ID = 1:10, time = seq_len(10), event = rep(1, 10),
                  x = factor(rep(c("a", "b"), 5)), unused = 11:20)
  f <- list(os = survival::Surv(time, event) ~ x,
            pfs = survival::Surv(time, event) ~ x)
  sf <- function(data = d, formulas = f) sampling_cache_fingerprint(
    formulas, data, list(os = "gamma", pfs = "gamma"), 10L, 123L, "mvn_v1")
  expected <- sf()
  df_changed <- d
  df_changed$unused <- 0
  df_changed$time <- as.double(df_changed$time)
  stopifnot(identical(sf(df_changed)$fingerprint, expected$fingerprint),
            identical(sf(d[, rev(names(d))])$fingerprint, expected$fingerprint),
            identical(cache_fingerprint(1:100), cache_fingerprint(c(1:100))),
            identical(names(expected$canonical_data), c("ID", "event", "time", "x")),
            !is.null(expected$runtime$libraries))
  df_changed$time[1] <- 99
  stopifnot(sf(df_changed)$fingerprint != expected$fingerprint)
  df_changed <- d
  df_changed$x <- stats::relevel(df_changed$x, "b")
  stopifnot(sf(df_changed)$fingerprint != expected$fingerprint)
  v_diagnostics <- capture.output(report_fingerprint_differences(expected$inputs, sf(df_changed)$inputs))
  stopifnot(any(grepl("x:", v_diagnostics, fixed = TRUE)))
  # Fresh R processes must agree; serialize v2 removes ALTREP representation
  # differences without depending on the review session's unknown environment.
  hash_script <- file.path(directory, "hash-session.R")
  writeLines(c('source("R/cache_paths.R")', 'source("R/model_configs.R")',
    paste0("d <- ", paste(capture.output(dput(d)), collapse = "\n")),
    paste0("f <- ", paste(capture.output(dput(f)), collapse = "\n")),
    'cat("HASH:", sampling_cache_fingerprint(f, d, list(os="gamma", pfs="gamma"), 10L, 123L, "mvn_v1")$fingerprint, "\\n")'),
    hash_script)
  rscript <- file.path(R.home("bin"), "Rscript.exe")
  v_hashes <- replicate(2, {
    output <- system2(rscript, shQuote(hash_script), stdout = TRUE, stderr = TRUE)
    sub("^HASH: ", "", trimws(grep("^HASH:", output, value = TRUE)))
  })
  stopifnot(length(v_hashes) == 2L, all(v_hashes == expected$fingerprint))
  cat("PASS: sampling identity hashes relevant values and records column diagnostics\n")

  # Check the real regeneration branch from script 06, not just the guard in
  # isolation; the sentinel must survive missing and stale sampling inputs.
  exprs <- parse("analysis/06_sampling.R", keep.source = FALSE)
  blocks <- Filter(function(e) is.call(e) && identical(e[[1]], as.name("if")) &&
    identical(paste(deparse(e[[2]]), collapse = ""), "is.null(sampling_models)"), as.list(exprs))
  stopifnot(length(blocks) == 1L)
  e <- new.env(parent = environment())
  e$sampling_models <- NULL
  e$cache_file <- file.path(directory, "sentinel.rds")
  saveRDS(list(untouched = TRUE), e$cache_file)
  before <- tools::md5sum(e$cache_file)
  fails(eval(blocks[[1]], e), "regeneration is disabled")
  stopifnot(identical(before, tools::md5sum(e$cache_file)))
  unlink(e$cache_file)
  fails(eval(blocks[[1]], e), "regeneration is disabled")
  stopifnot(!file.exists(e$cache_file))

  # An interactive calculation edit must invalidate the PSA even if the input
  # parameters and source file have not changed.
  fingerprint <- function() psa_cache_fingerprint("sampling", list(cost = 1),
    list(x = list(dist = "norm", mean = 0, sd = 1)), "control", 3L, 1L, 52, 1/52)
  original <- get0("calculate_outcomes", envir = .GlobalEnv, inherits = FALSE)
  on.exit(if (is.null(original)) {
    if (exists("calculate_outcomes", envir = .GlobalEnv, inherits = FALSE))
      rm("calculate_outcomes", envir = .GlobalEnv)
  } else
    assign("calculate_outcomes", original, .GlobalEnv), add = TRUE)
  baseline <- fingerprint()
  assign("calculate_outcomes", function(...) 2, envir = .GlobalEnv)
  stopifnot(fingerprint()$fingerprint != baseline$fingerprint)
  if (is.null(original)) rm("calculate_outcomes", envir = .GlobalEnv) else
    assign("calculate_outcomes", original, .GlobalEnv)
  cat("PASS: readers cannot regenerate sampling; interactive model edits invalidate PSA\n")

  po <- list(cost = data.frame(control = c(1, 2, 3), guided = c(3, 2, 1)),
    effect = data.frame(control = c(1, 1, 1), guided = c(1, 2, 3)),
    strategies = c("control", "guided"), n_sim = 3L, model_idx = c(1L, 3L, 3L),
    fingerprint = baseline$fingerprint)
  df_pp <- data.frame(sim = 1:3, model_idx = po$model_idx, x = c(1, 2, 3))
  attr(df_pp, "seed") <- 1L
  pair <- bind_psa_pair(po, df_pp)
  validate_psa_pair(pair$psa_obj, pair$psa_params)
  saveRDS(pair$psa_obj, psa_obj_path())
  saveRDS(pair$psa_params, psa_params_path())
  load_psa_params_cache(required = TRUE, seed = 1L, n_sim = 3L, verbose = FALSE)
  load_psa_cache(required = TRUE, verbose = FALSE)
  unlink(psa_params_path())
  fails(load_psa_cache(required = TRUE), "parameter cache missing")
  saveRDS(pair$psa_params, psa_params_path())
  unlink(psa_obj_path())
  fails(load_psa_params_cache(required = TRUE), "outcome cache missing")
  saveRDS(pair$psa_obj, psa_obj_path())
  df_reversed <- pair$psa_params[3:1, ]
  saveRDS(df_reversed, psa_params_path())
  fails(load_psa_params_cache(required = TRUE), "content hash mismatch")
  fails(load_psa_cache(required = TRUE), "content hash mismatch")
  df_pp$x[1] <- 10
  other <- bind_psa_pair(po, df_pp)
  fails(validate_psa_pair(pair$psa_obj, other$psa_params), "mixed generation")
  corrupted <- pair$psa_obj
  corrupted$cost[1, 1] <- 999
  fails(validate_psa_pair(corrupted, pair$psa_params), "content hash mismatch")
  fails(validate_psa_pair(po, df_pp), "missing or mixed")
  df_duplicate <- df_pp
  df_duplicate$sim <- c(1L, 1L, 3L)
  fails(bind_psa_pair(po, df_duplicate), "indices are not aligned")
  cat("PASS: shuffled parameters, mixed saves, changed outcomes and duplicate IDs rejected\n")

  config <- list(params = "x", groups = list(group = "x"))
  pop <- list(annual_incidence = 1500, research_horizon = 10, discount_rate = 0.035)
  ef <- function(wtp = 1, cfg = config, seed = 123L, population = pop,
                 settings = list(method = "gam", B = 1000L, se = TRUE))
    evppi_cache_fingerprint(pair$psa_obj, wtp, cfg, seed, population, settings)
  base <- ef()
  stopifnot(ef(2)$fingerprint != base$fingerprint,
    ef(cfg = list(params = "y"))$fingerprint != base$fingerprint,
    ef(seed = 124L)$fingerprint != base$fingerprint,
    ef(population = modifyList(pop, list(annual_incidence = 100)))$fingerprint != base$fingerprint,
    ef(settings = list(B = 2L))$fingerprint != base$fingerprint)
  evppi_results <- data.frame()
  m_nmb <- as.matrix(pair$psa_obj$effect) - as.matrix(pair$psa_obj$cost)
  evpi_manual <- mean(apply(m_nmb, 1, max)) - max(colMeans(m_nmb))
  evppi_fingerprint <- base$fingerprint
  save(evppi_results, evpi_manual, evppi_fingerprint, file = evppi_path())
  load_evppi_cache(pair$psa_obj, base, 1)
  fails(load_evppi_cache(pair$psa_obj, ef(2), 2), "fingerprint mismatch")
  evpi_manual <- evpi_manual + 10
  save(evppi_results, evpi_manual, evppi_fingerprint, file = evppi_path())
  fails(load_evppi_cache(pair$psa_obj, base, 1), "Cached EVPI disagrees")

  df_scenarios <- data.frame(scenario_id = "base", wtp = 1, c_drug_nivo = 100)
  scf <- function(psa = baseline$fingerprint, sc = df_scenarios) scenario_cache_fingerprint(
    psa, sc, config, 123L, pop)
  sid <- scf()
  cache <- list(fingerprint = sid$fingerprint, scenarios = df_scenarios,
                all_scenario_results = list(base = pair))
  cache$result_hash <- scenario_result_hash(cache)
  saveRDS(cache, scenario_evppi_path())
  load_scenario_cache(sid)
  fails(load_scenario_cache(scf(psa = "new code or parameters")), "fingerprint mismatch")
  df_scenarios$wtp <- 2
  fails(load_scenario_cache(scf(sc = df_scenarios)), "fingerprint mismatch")
  df_scenarios$c_drug_nivo <- 50
  fails(load_scenario_cache(scf(sc = df_scenarios)), "fingerprint mismatch")
  cat("PASS: EVPPI/scenario identities track WTP, groups, estimator, seed, population and PSA\n")

  # A forced failed assertion must be visible to an executable-test runner.
  v_index_script <- readLines("tests/test_sim_idx_validation.R", warn = FALSE)
  v_index_script <- sub('^results\\["valid"\\] <-.*', 'results["valid"] <- FALSE', v_index_script)
  failed_script <- file.path(directory, "forced-index-failure.R")
  writeLines(v_index_script, failed_script)
  status <- suppressWarnings(system2(rscript, shQuote(failed_script),
    stdout = file.path(directory, "forced.out"), stderr = file.path(directory, "forced.err")))
  stopifnot(status != 0L)
  cat("PASS: forced index assertion failure returns a failing process status\n")
}

run_cache_provenance_tests()
cat("All cache provenance tests passed.\n")
