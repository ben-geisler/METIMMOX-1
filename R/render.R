# ===============================================================================
# RENDERING: REPORTS AND PUBLICATION VIGNETTES
# ===============================================================================
# Sources live in reports/: the reports at the top level and in technical/, the
# publication vignettes in vignettes/. Everything a render produces lives in
# outputs/: rendered reports (.pdf, .md and <name>_files/) in outputs/reports/
# with the same subfolders, and the vignettes' HTML in outputs/vignettes/. The
# vignettes write their own figures and tables (outputs/figs/, outputs/tables/
# and outputs/manifest.csv; see R/publication_artifacts.R).
#
# Render order matters: the PDF pass deletes <name>_files/, so every report is
# rendered --to pdf first and --to gfm second, and the products are moved to
# outputs/ only afterwards. Vignettes are rendered one at a time because they
# share outputs/manifest.csv. Rendering needs the confidential trial data and
# the local caches; make.R runs it after the analysis scripts.
# ===============================================================================

#' Report sources under reports/ (excluding the vignettes)
#'
#' @param only Optional character vector of report basenames (without .qmd).
#' @return Repository-relative paths of the report .qmd files.
report_sources <- function(only = NULL) {
  v_files <- list.files(here::here("reports"), pattern = "\\.qmd$", recursive = TRUE)
  v_files <- file.path("reports", v_files[!startsWith(v_files, "vignettes/")])
  select_sources(v_files, only)
}

#' Publication vignette sources under reports/vignettes/
#'
#' @param only Optional character vector of vignette basenames (without .qmd).
#' @return Repository-relative paths of the vignette .qmd files.
vignette_sources <- function(only = NULL) {
  v_files <- file.path("reports", "vignettes",
                       list.files(here::here("reports", "vignettes"), pattern = "\\.qmd$"))
  select_sources(v_files, only)
}

#' Keep the sources named in `only`
#'
#' @param v_files Repository-relative .qmd paths.
#' @param only Optional character vector of basenames; NULL keeps every file.
#' @return The selected paths; stops if a requested name does not exist.
select_sources <- function(v_files, only) {
  if (is.null(only)) return(v_files)
  v_names <- tools::file_path_sans_ext(basename(v_files))
  unknown <- setdiff(only, v_names)
  if (length(unknown)) stop("No such source: ", paste(unknown, collapse = ", "))
  v_files[v_names %in% only]
}

#' Output folder for the products of one source file
#'
#' @param qmd Repository-relative path of a report or vignette .qmd.
#' @return Repository-relative output folder: outputs/vignettes/ for vignettes,
#'   outputs/reports/ (with the same subfolder) for reports.
render_output_dir <- function(qmd) {
  source_dir <- dirname(qmd)
  if (source_dir == "reports/vignettes") return("outputs/vignettes")
  sub("^reports", "outputs/reports", source_dir)
}

#' Run quarto render on one file
#'
#' @param qmd Repository-relative path of the .qmd file.
#' @param to Quarto output format ("pdf", "gfm" or "html").
#' @return Invisible NULL; stops if quarto exits with a non-zero status.
run_quarto <- function(qmd, to) {
  cat(sprintf("quarto render %s --to %s\n", qmd, to))
  status <- system2("quarto", c("render", shQuote(here::here(qmd)), "--to", to))
  if (!identical(status, 0L)) stop("Render failed: ", qmd, " (", to, ")")
  invisible(NULL)
}

#' Move the render products of one source into outputs/
#'
#' Replaces any earlier products of the same source in the output folder, so a
#' removed image or format never survives from an older render.
#'
#' @param qmd Repository-relative path of the .qmd file.
#' @param v_suffixes File suffixes to move (e.g. ".pdf", ".md").
#' @return Repository-relative paths of the moved files and folders.
move_render_products <- function(qmd, v_suffixes) {
  name <- tools::file_path_sans_ext(basename(qmd))
  source_dir <- here::here(dirname(qmd))
  out_dir <- render_output_dir(qmd)
  dir.create(here::here(out_dir), recursive = TRUE, showWarnings = FALSE)
  v_products <- c(paste0(name, v_suffixes), paste0(name, "_files"))
  v_products <- v_products[file.exists(file.path(source_dir, v_products))]
  for (product in v_products) {
    target <- here::here(out_dir, product)
    if (file.exists(target)) unlink(target, recursive = TRUE)
    if (!file.rename(file.path(source_dir, product), target))
      stop("Could not move ", file.path(dirname(qmd), product), " to ", out_dir)
  }
  file.path(out_dir, v_products)
}

#' Run git in the repository and return its output
#'
#' @param args Character vector of git arguments.
#' @return Trimmed standard output; stops on a non-zero exit status.
git_output <- function(args) {
  result <- suppressWarnings(system2("git", c("-C", shQuote(here::here()), args),
                                     stdout = TRUE, stderr = TRUE))
  status <- attr(result, "status")
  if (!is.null(status) && status != 0L)
    stop("git ", paste(args, collapse = " "), " failed: ", paste(result, collapse = "\n"))
  trimws(paste(result, collapse = "\n"))
}

#' Write the provenance record of a rendered PDF
#'
#' @param qmd Repository-relative path of the report source.
#' @param pdf Absolute path of the freshly rendered PDF.
#' @return Invisible path of the written <name>.pdf.provenance.json sidecar.
write_render_provenance <- function(qmd, pdf) {
  v_caches <- list.files(here::here("data", "tidy"),
    pattern = "^(sampling_models_n.*_full|psa_obj_.*|psa_params_.*|evppi_results_.*|scenario_evppi_results_.*)\\.(rds|RData)$",
    full.names = TRUE)
  provenance <- list(
    source_commit = git_output(c("rev-parse", "HEAD")),
    source_md5 = unname(tools::md5sum(here::here(qmd))),
    source_worktree_dirty = nzchar(git_output(c("status", "--porcelain", "--", qmd, "R", "analysis"))),
    render_time = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    cache_fingerprints = as.list(stats::setNames(unname(tools::md5sum(v_caches)), basename(v_caches))),
    pdf_md5 = unname(tools::md5sum(pdf)))
  sidecar <- paste0(pdf, ".provenance.json")
  jsonlite::write_json(provenance, sidecar, auto_unbox = TRUE, pretty = TRUE)
  invisible(sidecar)
}

#' Render one report to PDF and GitHub-flavoured Markdown
#'
#' @param qmd Repository-relative path of the report .qmd.
#' @return Repository-relative paths of the products in outputs/reports/.
render_report <- function(qmd) {
  run_quarto(qmd, "pdf")
  pdf <- here::here(sub("\\.qmd$", ".pdf", qmd))
  if (!file.exists(pdf)) stop("Render returned without PDF: ", qmd)
  write_render_provenance(qmd, pdf)
  run_quarto(qmd, "gfm")
  move_render_products(qmd, c(".pdf", ".pdf.provenance.json", ".md"))
}

#' Render one publication vignette
#'
#' The vignette writes its declared figures and tables itself; its HTML page is
#' moved to outputs/vignettes/.
#'
#' @param qmd Repository-relative path of the vignette .qmd.
#' @return Repository-relative path of the moved HTML page.
render_vignette <- function(qmd) {
  run_quarto(qmd, "html")
  move_render_products(qmd, ".html")
}

#' Render every report (or a subset), one after the other
#'
#' @param only Optional character vector of report basenames.
#' @return Invisible list of moved products per report.
render_reports <- function(only = NULL) {
  invisible(lapply(stats::setNames(nm = report_sources(only)), render_report))
}

#' Render every publication vignette (or a subset), one after the other
#'
#' @param only Optional character vector of vignette basenames.
#' @return Invisible list of moved products per vignette.
render_vignettes <- function(only = NULL) {
  invisible(lapply(stats::setNames(nm = vignette_sources(only)), render_vignette))
}
