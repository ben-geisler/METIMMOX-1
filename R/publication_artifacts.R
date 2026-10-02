# Publication output lifecycle (issue #168). Vignettes declare ownership before
# setup, remove predecessors, verify CSVs against their generating objects, and
# finish with an inventory. Render vignettes sequentially (one manifest writer).
publication_outputs <- function() {
  v_figs <- c("figure1", "figure2", "figure3", "figure4", "figure5", "figure_s1",
    "figure_s2", "figure_s5", "clin_effect_figure1", "clin_effect_figure_s1",
    "clin_effect_figure_sensitivity_dag", "clin_effect_figure_simplified_dag")
  v_tabs <- c("table_1", "table_4", "table_5", "table_6", paste0("table_s", 1:4),
    "table_s6", "table_s7", "table_s8", "clin_effect_table_s2")
  registry <- c(setNames(lapply(v_figs, function(x) paste0("outputs/figs/", x, ".png")), v_figs),
    setNames(lapply(v_tabs, function(x) paste0("outputs/tables/", x, ".csv")), v_tabs))
  registry$figure3 <- c(registry$figure3, "outputs/figs/figure_s3.png")
  registry$figure4 <- c(registry$figure4, "outputs/figs/figure_s4.png")
  registry$clin_effect_figure1 <- c(registry$clin_effect_figure1, "outputs/figs/clin_effect_figure1.eps")
  registry$table_4 <- c(registry$table_4, "outputs/tables/table_s5.csv")
  registry$table_s3 <- c(registry$table_s3, "outputs/tables/table_s3_pairs.csv")
  registry$figure_pfs_os_curves <- paste0("outputs/figs/figure_pfs_os_curves_",
    c("soc", "crp", "tmb_braf"), ".png")
  registry$poster_SMDM <- c(paste0("outputs/figs/poster_", c("forest", "ceac"), ".png"),
    paste0("outputs/tables/poster_", c("cea", "forest_data", "ceac_data"), ".csv"))
  registry
}

artifact_manifest_path <- function(root) file.path(root, "outputs/manifest.csv")

begin_artifact_render <- function(vignette, root = here::here()) {
  v_files <- publication_outputs()[[vignette]]
  if (is.null(v_files)) stop("No publication output contract for: ", vignette)
  ctx <- new.env(parent = emptyenv())
  ctx$root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  ctx$vignette <- paste0("outputs/vignettes/", vignette, ".qmd")
  ctx$files <- v_files
  ctx$written <- list()
  ctx$inputs <- character()
  ctx$previous <- lapply(file.path(ctx$root, v_files), function(f)
    if (file.exists(f)) unname(tools::md5sum(f)) else NA_character_)
  names(ctx$previous) <- v_files
  # Every deletion is an exact, registered file beneath this repository root.
  for (f in file.path(ctx$root, v_files)) {
    if (file.exists(f) && !file.remove(f)) stop("Cannot remove predecessor: ", f)
  }
  manifest <- artifact_manifest_path(ctx$root)
  if (file.exists(manifest)) {
    df_old <- read.csv(manifest, stringsAsFactors = FALSE)
    utils::write.csv(df_old[!df_old$file %in% v_files, , drop = FALSE], manifest, row.names = FALSE)
  }
  options(metimmox.artifact_context = ctx)
  invisible(ctx)
}

#' Record an input file other than data/tidy/METIMMOX.rds that the active
#' render reads (issue #161), so the manifest carries its md5 digest.
#'
#' @param file Path relative to the repository root.
#' @return The path, invisibly.
register_artifact_input <- function(file) {
  ctx <- getOption("metimmox.artifact_context")
  if (is.null(ctx)) stop("Call begin_artifact_render() before registering an input")
  path <- file.path(ctx$root, file)
  if (!file.exists(path)) stop("Registered artifact input does not exist: ", file)
  ctx$inputs[file] <- unname(tools::md5sum(path))
  invisible(file)
}

artifact_target <- function(file) {
  ctx <- getOption("metimmox.artifact_context")
  if (is.null(ctx)) stop("Call begin_artifact_render() before writing an artifact")
  path <- normalizePath(file, winslash = "/", mustWork = FALSE)
  v_targets <- file.path(ctx$root, ctx$files)
  idx <- match(path, v_targets)
  if (is.na(idx)) stop("Undeclared publication output: ", file)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  list(ctx = ctx, path = path, relative = ctx$files[idx])
}

write_artifact_csv <- function(x, file, ...) {
  target <- artifact_target(file)
  # Serialize the generating object independently in memory, then compare the
  # actual on-disk CSV including column names, row names, NA and numeric text.
  expected <- character()
  con <- textConnection("expected", "w", local = TRUE)
  tryCatch(utils::write.csv(x, con, ...), finally = close(con))
  utils::write.csv(x, target$path, ...)
  if (!identical(enc2utf8(readLines(target$path, warn = FALSE)), enc2utf8(expected)))
    stop("CSV differs from its generating object: ", file)
  target$ctx$written[[target$relative]] <- "csv_verified"
  invisible(file)
}

save_artifact_plot <- function(filename, plot, ...) {
  target <- artifact_target(filename)
  ggplot2::ggsave(filename = target$path, plot = plot, ...)
  target$ctx$written[[target$relative]] <- "written"
  invisible(filename)
}

finish_artifact_render <- function(envir = knitr::knit_global()) {
  ctx <- getOption("metimmox.artifact_context")
  if (is.null(ctx)) stop("No active publication render")
  # Record only cache objects actually brought into the vignette environment.
  v_fingerprints <- character()
  for (name in c("sampling_models", "psa_obj", "base_psa", "pipeline_psa",
                 "scenario_cache", "evppi_cache")) {
    obj <- get0(name, envir = envir, inherits = FALSE)
    if (is.list(obj) && length(obj$fingerprint) == 1L)
      v_fingerprints[name] <- obj$fingerprint
  }
  v_source_files <- c(list.files(file.path(ctx$root, "R"), "\\.R$", full.names = TRUE),
    list.files(file.path(ctx$root, "analysis"), "\\.R$", full.names = TRUE),
    file.path(ctx$root, ctx$vignette))
  v_source_files <- v_source_files[file.exists(v_source_files)]
  code_hash <- digest::digest(unname(tools::md5sum(sort(v_source_files, method = "radix"))), algo = "sha256")
  trial <- file.path(ctx$root, "data/tidy/METIMMOX.rds")
  df_rows <- lapply(ctx$files, function(f) {
    path <- file.path(ctx$root, f)
    exists <- file.exists(path)
    status <- ctx$written[[f]]
    if (exists && is.null(status)) stop("Output bypassed artifact writer: ", f)
    hash <- if (exists) unname(tools::md5sum(path)) else NA_character_
    data.frame(file = f, generating_vignette = ctx$vignette,
      source_cache_fingerprints = if (length(v_fingerprints))
        paste(names(v_fingerprints), v_fingerprints, sep = "=", collapse = ";") else "",
      render_time = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
      status = if (exists) status else "deleted_empty",
      artifact_md5 = hash, source_code_sha256 = code_hash,
      trial_data_md5 = if (file.exists(trial)) unname(tools::md5sum(trial)) else NA_character_,
      additional_inputs_md5 = if (length(ctx$inputs)) paste(names(ctx$inputs), ctx$inputs,
        sep = "=", collapse = ";") else "",
      comparison_with_predecessor = if (!exists) "removed" else if (is.na(ctx$previous[[f]]))
        "created" else if (identical(hash, ctx$previous[[f]])) "unchanged" else "changed",
      stringsAsFactors = FALSE)
  })
  df_rows <- do.call(rbind, df_rows)
  manifest <- artifact_manifest_path(ctx$root)
  if (file.exists(manifest)) {
    df_old <- read.csv(manifest, stringsAsFactors = FALSE)
    # Manifests written before issue #161 lack the additional-inputs column.
    if (!"additional_inputs_md5" %in% names(df_old)) df_old$additional_inputs_md5 <- ""
    df_old <- df_old[, names(df_rows), drop = FALSE]
    df_rows <- rbind(df_old[!df_old$file %in% ctx$files, , drop = FALSE], df_rows)
  }
  # Locale-independent row order, so the manifest diff is stable (#185).
  df_rows <- df_rows[order(df_rows$file, method = "radix"), ]
  utils::write.csv(df_rows, manifest, row.names = FALSE, na = "")
  options(metimmox.artifact_context = NULL)
  invisible(df_rows)
}
