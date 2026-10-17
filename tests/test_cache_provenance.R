# Synthetic adversarial cache contracts: no production cache is read or written.
source("R/cache_paths.R")
source("R/cea_helpers.R")

#' The calculation identity of the four model categories as computed before
#' issue #205 (file and package lists written out here, files read through
#' here::here()), against which calculation_identity() must stay identical.
legacy_calculation_identity <- function(stage) {
  v_files <- switch(stage,
    sampling = c("analysis/06_sampling.R", "R/model_configs.R", "R/joint_survival_sampling.R"),
    psa = c("R/model_fun.R", "R/calculate_outcomes.R", "R/prediction_functions.R",
            "R/psa_functions.R", "R/parameter_distributions.R", "R/model_configs.R",
            "R/prediction_population.R",
            "analysis/06_sampling.R", "R/joint_survival_sampling.R"),
    evppi = "R/evppi_functions.R",
    scenario = c("R/scenario_analysis.R", "R/evppi_functions.R"))
  code <- lapply(v_files, function(file) {
    exprs <- parse(here::here(file), keep.source = FALSE)
    definitions <- list()
    for (expr in exprs) {
      if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
          is.symbol(expr[[2]]) && is.call(expr[[3]]) &&
          identical(expr[[3]][[1]], as.name("function"))) {
        name <- as.character(expr[[2]])
        fn <- get0(name, envir = .GlobalEnv, mode = "function", inherits = FALSE)
        if (is.null(fn)) fn <- eval(expr[[3]], envir = baseenv())
        definitions[[name]] <- list(formals = paste(deparse(formals(fn)), collapse = "\n"),
                                    body = paste(deparse(body(fn)), collapse = "\n"))
      }
    }
    list(source = paste(deparse(exprs), collapse = "\n"), definitions = definitions)
  })
  names(code) <- v_files
  v_packages <- switch(stage, sampling = c("flexsurv", "survival", "mvtnorm"),
    psa = c("flexsurv", "survival", "dampack", "mvtnorm"),
    evppi = c("voi", "mgcv", "dampack"), scenario = c("voi", "mgcv", "dampack"))
  list(code = vapply(code, cache_fingerprint, character(1)),
       packages = setNames(vapply(v_packages, function(p)
         as.character(utils::packageVersion(p)), character(1)), v_packages),
       R = paste(R.version$major, R.version$minor, sep = "."),
       contrasts = getOption("contrasts"), rng = RNGkind())
}

#' Contracts of the value-of-information caches (issue #205) on a bound
#' synthetic PSA pair and scenario cache: the `voi` category leaves the model
#' categories unchanged; every VoI fingerprint changes when a paper 3 source is
#' edited (in a temporary copy of the sources) or an upstream cache changes,
#' and not when only comments change; the loaders name the differing input.
run_voi_identity_tests <- function(pair, scenario_cache, config, pop, directory, fails) {
  v_model_stages <- c("psa", "sampling", "evppi", "scenario")
  for (stage in v_model_stages) {
    stopifnot(identical(calculation_identity(stage), legacy_calculation_identity(stage)))
  }
  v_voi_files <- calculation_identity_sources("voi")$files
  stopifnot(identical(v_voi_files, c(
    "R/voi_price_axis.R", "R/voi_survival_parameters.R", "R/voi_policies.R", "R/voi_incidence.R",
    "R/design_information.R", "R/evidence_update.R", "R/evsi_functions.R", "R/voi_cache.R",
    "R/voi_designs.R", "R/voi_research_value.R", "R/voi_evidence_psa.R", "analysis/14_VOI.R")),
    all(file.exists(here::here(v_voi_files))),
    !length(intersect(v_voi_files, unlist(lapply(v_model_stages, function(s)
      calculation_identity_sources(s)$files)))),
    all(c("voi", "mgcv", "flexsurv") %in% calculation_identity_sources("voi")$packages))
  cat("PASS: psa/sampling/evppi/scenario identities equal the pre-#205 algorithm; voi lists its own sources\n")

  for (f in c("voi_cache", "voi_price_axis", "voi_research_value")) source(here::here("R", paste0(f, ".R")))
  sampling <- list(fingerprint = "sampling-A")
  df_threshold <- data.frame(strategy = c("control", "guided"), price = c(1, 2),
                             date = as.Date(c("2026-10-01", "2026-10-02")))
  fingerprints <- function(sampling_models = sampling, psa_obj = pair$psa_obj, scenario = scenario_cache,
                           base = "sampling-A", evidence_states = NULL) list(
    step1a = voi_step1a_identity(scenario, sampling_models, df_threshold, NULL, config,
                                 list(wtp = 1, scale = function(x) 2 * x), pop),
    research = voi_research_fingerprint(list(
      upstream = voi_upstream_fingerprints(sampling_models, psa_obj, scenario), settings = list(wtp = 1))),
    evidence_state = evidence_state_fingerprint(base, "E2", NULL, NA_character_, list(), NULL,
                                                list(shared = 0), data.frame(Age = 60), NULL, 10L, 1L, 1e-6, 1),
    evsi = evsi_fingerprint(voi_upstream_fingerprints(psa_obj = psa_obj, evidence_states = evidence_states),
                            list(), c(10, 20), 1, 1, 1L, 10L, "gam_regression"))
  hashes <- function(...) vapply(fingerprints(...), function(x) x$fingerprint, character(1))

  # Edited copies of the sources: the option points calculation_identity() at them.
  root <- file.path(directory, "sources")
  dir.create(root)
  file.copy(here::here(c("R", "analysis")), root, recursive = TRUE)
  old_root <- options(metimmox.calculation_root = root)
  on.exit(options(old_root), add = TRUE)
  base <- hashes()
  v_model_identity <- lapply(setNames(nm = v_model_stages), calculation_identity)
  options(old_root)
  stopifnot(identical(base, hashes()))
  options(metimmox.calculation_root = root)
  #' Hashes after replacing one source file of the copy by `lines`; restores it.
  edited_hashes <- function(file, lines) {
    path <- file.path(root, file)
    original <- readLines(path, warn = FALSE)
    on.exit(writeLines(original, path))
    writeLines(lines, path)
    list(voi = hashes(), model = lapply(setNames(nm = v_model_stages), calculation_identity))
  }
  for (file in v_voi_files) {
    original <- readLines(file.path(root, file), warn = FALSE)
    changed <- edited_hashes(file, c(original, "metimmox_issue205_probe <- 1L"))
    stopifnot(all(changed$voi != base), identical(changed$model, v_model_identity))
    commented <- edited_hashes(file, c("# a comment-only edit (#205)", original, "  # trailing comment"))
    stopifnot(identical(commented$voi, base), identical(commented$model, v_model_identity))
  }
  # An edit inside an existing definition (its formals) is caught as well.
  original <- readLines(file.path(root, "R/voi_policies.R"), warn = FALSE)
  v_first <- grep("<- function(", original, fixed = TRUE)[1]
  edited <- original
  edited[v_first] <- sub("<- function(", "<- function(probe_205 = NULL, ", edited[v_first], fixed = TRUE)
  stopifnot(all(edited_hashes("R/voi_policies.R", edited)$voi != base))
  # Step 1a and steps 5-6 also depend on the EVPPI estimator.
  original <- readLines(file.path(root, "R/evppi_functions.R"), warn = FALSE)
  changed <- edited_hashes("R/evppi_functions.R", c(original, "metimmox_issue205_probe <- 1L"))$voi
  stopifnot(changed[c("step1a", "research")] != base[c("step1a", "research")],
            identical(changed[c("evidence_state", "evsi")], base[c("evidence_state", "evsi")]))
  options(old_root)
  cat("PASS: editing any paper 3 source changes every VoI fingerprint and no model identity;",
      "comment-only edits change nothing\n")

  # Upstream caches: each recorded fingerprint invalidates the results.
  other_pair <- pair
  other_pair$psa_obj$pair_metadata$generation_id <- "another PSA generation"
  other_psa <- pair$psa_obj
  other_psa$fingerprint <- "another PSA input fingerprint"
  other_scenario <- scenario_cache
  other_scenario$fingerprint <- "another scenario fingerprint"
  regenerated_scenario <- scenario_cache
  regenerated_scenario$all_scenario_results$base <- other_pair
  stopifnot(
    all(hashes(sampling_models = list(fingerprint = "sampling-B"))[c("step1a", "research")] !=
          base[c("step1a", "research")]),
    hashes(scenario = other_scenario)[c("step1a", "research")] != base[c("step1a", "research")],
    hashes(scenario = regenerated_scenario)[["step1a"]] != base[["step1a"]],
    all(hashes(psa_obj = other_pair$psa_obj)[c("research", "evsi")] != base[c("research", "evsi")]),
    all(hashes(psa_obj = other_psa)[c("research", "evsi")] != base[c("research", "evsi")]),
    hashes(base = "sampling-B")[["evidence_state"]] != base[["evidence_state"]],
    hashes(evidence_states = list(E2 = list(fingerprint = "E2-state")))[["evsi"]] != base[["evsi"]])
  unbound <- pair$psa_obj
  unbound$pair_metadata <- NULL
  fails(voi_upstream_fingerprints(psa_obj = unbound), "PSA pair generation ID is missing")
  fails(voi_upstream_fingerprints(sampling_models = list()), "sampling-cache fingerprint is missing")
  fails(voi_research_fingerprint(list(settings = 1)), "needs inputs\\$upstream")
  fails(evsi_fingerprint(list(sampling = "s"), list(), 1, 1, 1, 1L, 1L, "m"), "needs upstream\\$psa")
  cat("PASS: VoI fingerprints record the sampling, PSA, scenario and evidence-state fingerprints\n")

  # Loaders: unchanged inputs read back (after an RDS round trip); a stale cache
  # names the differing input, stops for readers and returns NULL for producers.
  identities <- fingerprints()
  path <- evsi_cache_path("research_contract", directory)
  save_voi_results(list(value = 42), identities$research, path)
  stopifnot(identical(load_voi_results(path, identities$research), list(value = 42)))
  stale <- fingerprints(psa_obj = other_pair$psa_obj)$research
  fails(load_voi_results(path, stale), "upstream\\$psa\\$generation_id")
  stopifnot(is.null(suppressMessages(load_voi_results(path, stale, required = FALSE))))
  fails(load_voi_results(file.path(directory, "evsi_results_absent.rds"), stale), "Missing")
  stopifnot(is.null(load_voi_results(file.path(directory, "evsi_results_absent.rds"), stale, required = FALSE)))
  options(metimmox.calculation_root = root)
  original <- readLines(file.path(root, "R/voi_policies.R"), warn = FALSE)
  writeLines(c(original, "metimmox_issue205_probe <- 1L"), file.path(root, "R/voi_policies.R"))
  edited_identity <- fingerprints()$research
  writeLines(original, file.path(root, "R/voi_policies.R"))
  options(old_root)
  fails(load_voi_results(path, edited_identity), "implementation\\$voi\\$code\\$R/voi_policies.R")
  tampered <- readRDS(path)
  tampered$inputs$settings$wtp <- 2
  saveRDS(tampered, path)
  fails(load_voi_results(path, identities$research), "do not match its saved fingerprint")
  fails(save_voi_results(1, list(fingerprint = "x", inputs = list()), path), "upstream fingerprints")

  state <- list(joint = list(method = "mvn_joint_v2", n_samples = 10L), evidence_state = "E2",
                fingerprint = identities$evidence_state$fingerprint,
                fingerprint_inputs = identities$evidence_state$inputs)
  state_path <- save_evidence_state(state, evidence_state_cache_path("E2", 10L, directory))
  stopifnot(identical(load_evidence_state(state_path, identities$evidence_state), state),
            identical(load_evidence_state(state_path, identities$evidence_state$fingerprint), state))
  stale_state <- fingerprints(base = "sampling-B")$evidence_state
  fails(load_evidence_state(state_path, stale_state, required = TRUE), "base_sampling_fingerprint")
  stopifnot(is.null(suppressMessages(load_evidence_state(state_path, stale_state))))
  fails(load_evidence_state(file.path(directory, "absent.rds"), stale_state, required = TRUE), "Missing")

  fails(check_voi_identity(identities$step1a$fingerprint, identities$step1a$inputs,
                           fingerprints(scenario = other_scenario)$step1a, "step 1a", "Rerun."),
        "upstream\\$scenario\\$fingerprint")
  fails(check_voi_identity(identities$step1a$fingerprint, identities$step1a$inputs,
                           fingerprints(scenario = regenerated_scenario)$step1a, "step 1a", "Rerun."),
        "upstream\\$scenario\\$psa_generations\\$base")
  # load_or_run_voi_step1a() reuses a current result without recomputing.
  step1a_file <- file.path(directory, "voi", "voi_step1a_results.rds")
  dir.create(dirname(step1a_file))
  saveRDS(list(price_axis = data.frame(c_drug_nivo = 1), evppi = data.frame(group = "g"),
               identity = identities$step1a), step1a_file)
  reused <- load_or_run_voi_step1a(scenario_cache, sampling, df_threshold, NULL, config, pop,
                                   list(wtp = 1, scale = function(x) 2 * x), cache_file = step1a_file)
  stopifnot(identical(reused$computed, FALSE),
            file.exists(file.path(dirname(step1a_file), c("voi_price_axis.csv", "evppi_extended.csv"))))
  cat("PASS: VoI loaders validate the upstream fingerprints and name the differing input\n")
}

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
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
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
    ly_disc = data.frame(control = c(2, 2, 2), guided = c(2, 3, 4)),
    ly_undisc = data.frame(control = c(3, 3, 3), guided = c(3, 4, 5)),
    strategies = c("control", "guided"), n_sim = 3L, model_idx = c(1L, 3L, 3L),
    fingerprint = baseline$fingerprint)
  df_pp <- data.frame(sim = 1:3, model_idx = po$model_idx, x = c(1, 2, 3))
  attr(df_pp, "seed") <- 1L
  # A cost/QALY-only object (the psa_pair_v1 content) is refused, never stamped
  # (issues #197, #200).
  po_cost_qaly <- po[setdiff(names(po), c("ly_disc", "ly_undisc"))]
  fails(bind_psa_pair(po_cost_qaly, df_pp), "Refusing to bind a legacy PSA object: no ly_disc")
  # With life-years but without cells (the content of the #197 branch before
  # the merge with #200) the object is refused too.
  fails(bind_psa_pair(po, df_pp), "no cell-level component")
  v_cells <- c("00", "01", "10", "11")
  v_arms <- c("control", "experimental")
  v_outcomes <- c("cost", "qalys", "ly_disc", "ly_undisc", "nivo_admin")
  po$cells <- list(schema = "psa_cells_v1", cells = v_cells, arms = v_arms, outcomes = v_outcomes,
    values = array(1, c(3, 4, 2, 5), dimnames = list(NULL, v_cells, v_arms, v_outcomes)),
    weights = matrix(0.25, 3, 4, dimnames = list(NULL, v_cells)),
    violations = array(0L, c(3, 4, 2)), nivolumab_price = rep(1, 3))
  pair <- bind_psa_pair(po, df_pp)
  stopifnot(identical(pair$psa_obj$pair_metadata$schema, "psa_pair_v2"))
  legacy <- pair
  legacy$psa_obj$pair_metadata$schema <- "psa_pair_v1"
  attr(legacy$psa_params, "pair_metadata") <- legacy$psa_obj$pair_metadata
  fails(validate_psa_pair(legacy$psa_obj, legacy$psa_params), "legacy pair \\(schema psa_pair_v1\\)")
  corrupted_cells <- pair$psa_obj
  corrupted_cells$cells$values[1, 1, 1, 1] <- 2
  fails(validate_psa_pair(corrupted_cells, pair$psa_params), "content hash mismatch")
  cat("PASS: legacy cost/QALY-only PSA objects and v1 pairs are rejected; cells are hashed\n")
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
  # Life-years are bound outcomes (#197): changed, missing or misaligned
  # life-year draws invalidate the pair.
  corrupted <- pair$psa_obj
  corrupted$ly_undisc[2, 2] <- 99
  fails(validate_psa_pair(corrupted, pair$psa_params), "content hash mismatch")
  corrupted <- pair$psa_obj
  corrupted$ly_disc <- NULL
  fails(validate_psa_pair(corrupted, pair$psa_params), "life-year draws missing")
  corrupted <- pair$psa_obj
  corrupted$ly_undisc <- corrupted$ly_undisc[1:2, ]
  fails(validate_psa_pair(corrupted, pair$psa_params), "life-year draws missing or not aligned")
  corrupted <- pair$psa_obj
  corrupted$ly_disc <- corrupted$ly_disc[, c("guided", "control")]
  fails(validate_psa_pair(corrupted, pair$psa_params), "not aligned")
  po_no_ly <- po
  po_no_ly$ly_disc <- NULL
  fails(bind_psa_pair(po_no_ly, df_pp), "life-year draws missing")
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

  run_voi_identity_tests(pair, cache, config, pop, directory, fails)

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
