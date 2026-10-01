#!/usr/bin/env Rscript
# ===============================================================================
# RENDER AND PUBLISH THE REPORTS
# ===============================================================================
# Every report under reports/ (and reports/technical/) renders to PDF and to
# GitHub-flavoured Markdown. The .md renders and their <name>_files/ assets are
# tracked in git (LLM-readable); the PDFs are gitignored and published from a
# generated _site/ directory to the gh-pages branch (GitHub Pages) and/or as a
# GitHub release. Rendering needs the confidential trial data and the local
# caches, so this runs locally, never in CI.
#
# Usage (from the repository root):
#   Rscript publish/publish_reports.R                # render all reports, build _site/
#   Rscript publish/publish_reports.R --vignettes    # also render outputs/vignettes/*.qmd
#   Rscript publish/publish_reports.R --no-render    # rebuild _site/ from existing PDFs
#   Rscript publish/publish_reports.R --push         # ... and push _site/ to gh-pages
#   Rscript publish/publish_reports.R --release      # ... and attach the PDFs to a release
#   Rscript publish/publish_reports.R --only CEA,OWSA   # subset by report basename
#
# Render order matters: the PDF pass deletes <name>_files/, so each report is
# rendered --to pdf first and --to gfm second.
# ===============================================================================

args <- commandArgs(trailingOnly = TRUE)
flag <- function(x) x %in% args
opt  <- function(x) { i <- match(x, args); if (is.na(i) || i == length(args)) NULL else args[i + 1] }

root <- here::here()
setwd(root)
site_dir <- file.path(root, "_site")

report_files <- list.files("reports", pattern = "\\.qmd$", recursive = TRUE, full.names = TRUE)
if (!is.null(opt("--only"))) {
  keep <- strsplit(opt("--only"), ",")[[1]]
  report_files <- report_files[tools::file_path_sans_ext(basename(report_files)) %in% keep]
}
vignette_files <- if (flag("--vignettes")) list.files("outputs/vignettes", pattern = "\\.qmd$", full.names = TRUE) else character(0)

render <- function(file, to) {
  cat(sprintf("quarto render %s --to %s\n", file, to))
  status <- system2("quarto", c("render", shQuote(file), "--to", to))
  if (status != 0) stop("Render failed: ", file, " (", to, ")")
}

run_checked <- function(command, args, step) {
  status <- system2(command, args)
  if (!identical(status, 0L)) stop(step, " failed (exit status ", status, ")")
  invisible(NULL)
}

git_output <- function(args) {
  result <- suppressWarnings(system2("git", args, stdout = TRUE, stderr = TRUE))
  status <- attr(result, "status")
  if (!is.null(status) && status != 0L) stop("git ", paste(args, collapse = " "),
                                            " failed: ", paste(result, collapse = "\n"))
  trimws(paste(result, collapse = "\n"))
}

commit <- git_output(c("rev-parse", "HEAD"))
render_time <- function(path) format(file.info(path)$mtime, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
cache_fingerprints <- function() {
  paths <- list.files("data/tidy", pattern = "^(sampling_models_n.*_full|psa_obj_.*|psa_params_.*|evppi_results_.*|scenario_evppi_results_.*)\\.(rds|RData)$",
                      full.names = TRUE)
  if (!length(paths)) return(list())
  as.list(stats::setNames(unname(tools::md5sum(paths)), basename(paths)))
}
provenance_path <- function(pdf) paste0(pdf, ".provenance.json")
read_provenance <- function(pdf) {
  sidecar <- provenance_path(pdf)
  if (file.exists(sidecar)) {
    value <- jsonlite::read_json(sidecar, simplifyVector = TRUE)
    if (identical(value$pdf_md5, unname(tools::md5sum(pdf)))) return(value)
    stop("PDF changed since its render provenance was recorded: ", pdf)
  }
  list(source_commit = NA_character_, render_time = render_time(pdf),
       time_basis = "PDF file modification time (unverified legacy render)",
       cache_fingerprints = list(), pdf_md5 = unname(tools::md5sum(pdf)))
}

if (!flag("--no-render")) {
  for (f in report_files) {
    render(f, "pdf")
    pdf <- sub("\\.qmd$", ".pdf", f)
    if (!file.exists(pdf)) stop("Render returned without PDF: ", f)
    provenance <- list(source_commit = commit,
      source_md5 = unname(tools::md5sum(f)),
      source_worktree_dirty = nzchar(git_output(c("status", "--porcelain", "--", f, "R", "analysis"))),
      render_time = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
      time_basis = "recorded at render", cache_fingerprints = cache_fingerprints(),
      pdf_md5 = unname(tools::md5sum(pdf)))
    jsonlite::write_json(provenance, provenance_path(pdf), auto_unbox = TRUE, pretty = TRUE)
    render(f, "gfm")
  }
  for (f in vignette_files) render(f, "html")
}

# ---- Assemble _site/ ----------------------------------------------------------
unlink(site_dir, recursive = TRUE)
if (dir.exists(site_dir)) stop("Could not clear _site/")
if (!dir.create(file.path(site_dir, "technical"), recursive = TRUE, showWarnings = FALSE))
  stop("Could not create _site/technical/")
pdfs <- sub("\\.qmd$", ".pdf", report_files)
missing_pdfs <- pdfs[!file.exists(pdfs)]
if (length(missing_pdfs)) stop("Missing rendered PDFs: ", paste(missing_pdfs, collapse = ", "))
for (p in pdfs) {
  if (!file.copy(p, file.path(site_dir, sub("^reports/", "", p)), overwrite = TRUE))
    stop("Could not copy PDF: ", p)
}
records <- lapply(pdfs, function(p) c(list(file = sub("^reports/", "", p)), read_provenance(p)))
jsonlite::write_json(list(schema = "report_pdf_manifest_v1", reports = records),
                     file.path(site_dir, "manifest.json"), auto_unbox = TRUE, pretty = TRUE,
                     null = "null")

row <- function(p) {
  rel <- sub("^reports/", "", p)
  record <- records[[match(p, pdfs)]]
  source_commit <- if (is.na(record$source_commit)) "source commit unknown" else
    paste0("source ", substr(record$source_commit, 1, 7))
  if (isTRUE(record$source_worktree_dirty))
    source_commit <- paste0(source_commit, " (working tree had edits)")
  sprintf('<li><a href="%s">%s</a> <small>(%s; %s; %s KB)</small></li>', rel,
          tools::file_path_sans_ext(basename(p)),
          record$render_time, source_commit,
          format(round(file.info(p)$size / 1024), big.mark = ","))
}
index <- c(
  "<!doctype html><html><head><meta charset='utf-8'><title>METIMMOX-1 reports</title>",
  "<style>body{font-family:system-ui,sans-serif;max-width:48rem;margin:3rem auto;line-height:1.5}</style></head><body>",
  "<h1>METIMMOX-1 economic evaluation: rendered reports</h1>",
  "<p>Dates and source commits are recorded per PDF. Legacy PDFs without a provenance record show their file date and an unknown source commit. Cache fingerprints and PDF hashes are in <a href='manifest.json'>manifest.json</a>.</p>",
  "<h2>Reports</h2><ul>", vapply(pdfs[!grepl("/technical/", pdfs)], row, ""), "</ul>",
  "<h2>Technical documentation</h2><ul>", vapply(pdfs[grepl("/technical/", pdfs)], row, ""), "</ul>",
  "</body></html>"
)
writeLines(index, file.path(site_dir, "index.html"))
writeLines("", file.path(site_dir, ".nojekyll"))
cat(sprintf("_site/ assembled: %d PDFs + index.html\n", length(pdfs)))

# ---- Publish ------------------------------------------------------------------
if (flag("--push")) {
  publish_pages <- function() {
    wt <- file.path(tempdir(), "gh-pages-worktree")
    if (dir.exists(wt)) stop("Temporary worktree already exists: ", wt)
    local_branch <- system2("git", c("show-ref", "--quiet", "refs/heads/gh-pages")) == 0L
    remote_status <- system2("git", c("ls-remote", "--exit-code", "--heads", "origin", "gh-pages"),
                             stdout = FALSE, stderr = FALSE)
    if (!remote_status %in% c(0L, 2L)) stop("Could not check origin/gh-pages")
    remote_branch <- remote_status == 0L
    if (local_branch && remote_branch) {
      run_checked("git", c("fetch", "origin", "gh-pages"), "Fetch remote gh-pages")
      local_sha <- git_output(c("rev-parse", "gh-pages"))
      remote_sha <- git_output(c("rev-parse", "FETCH_HEAD"))
      if (!identical(local_sha, remote_sha))
        stop("Local gh-pages differs from origin/gh-pages; reconcile it before publishing")
    } else if (!local_branch && remote_branch) {
      run_checked("git", c("fetch", "origin", "gh-pages"), "Fetch remote gh-pages")
    }
    if (local_branch) {
      run_checked("git", c("worktree", "add", shQuote(wt), "gh-pages"), "Add gh-pages worktree")
    } else if (remote_branch) {
      run_checked("git", c("worktree", "add", "-b", "gh-pages", shQuote(wt), "FETCH_HEAD"),
                  "Add remote gh-pages worktree")
    } else {
      run_checked("git", c("worktree", "add", "--detach", shQuote(wt)), "Add detached worktree")
      run_checked("git", c("-C", shQuote(wt), "checkout", "--orphan", "gh-pages"),
                  "Create gh-pages orphan branch")
      tracked <- git_output(c("-C", shQuote(wt), "ls-files"))
      if (nzchar(tracked)) run_checked("git", c("-C", shQuote(wt), "rm", "-rfq", "."),
                                      "Clear orphan worktree index")
    }
    on.exit(run_checked("git", c("worktree", "remove", "--force", shQuote(wt)),
                        "Remove gh-pages worktree"), add = TRUE)
    old <- list.files(wt, all.files = TRUE, no.. = TRUE, full.names = TRUE)
    old <- old[basename(old) != ".git"]
    if (length(old)) unlink(old, recursive = TRUE)
    if (length(list.files(wt, all.files = TRUE, no.. = TRUE)) != 1L)
      stop("Could not clear stale gh-pages content")
    files <- list.files(site_dir, full.names = TRUE, all.files = TRUE, no.. = TRUE)
    if (!all(file.copy(files, wt, recursive = TRUE, overwrite = TRUE)))
      stop("Could not copy all site files to gh-pages worktree")
    run_checked("git", c("-C", shQuote(wt), "add", "-A"), "Stage gh-pages content")
    changed <- system2("git", c("-C", shQuote(wt), "diff", "--cached", "--quiet"))
    if (changed == 1L) {
      run_checked("git", c("-C", shQuote(wt), "commit", "-m", shQuote("Publish reports")),
                  "Commit gh-pages content")
    } else if (changed != 0L) stop("Could not inspect staged gh-pages changes")
    run_checked("git", c("-C", shQuote(wt), "push", "origin", "gh-pages"),
                "Push gh-pages")
  }
  publish_pages()
  cat("Pushed _site/ to gh-pages. Enable Pages on that branch in the repository settings.\n")
}

if (flag("--release")) {
  tag <- paste0("reports-", format(Sys.Date(), "%Y%m%d"))
  release_files <- c(list.files(site_dir, "\\.pdf$", recursive = TRUE, full.names = TRUE),
                     file.path(site_dir, "manifest.json"))
  status <- system2("gh", c("release", "create", tag, shQuote(release_files),
                            "--title", shQuote(paste("Reports", format(Sys.Date()))),
                            "--notes", shQuote("Per-PDF source commits, render times and cache fingerprints are recorded in manifest.json")))
  if (status != 0) stop("gh release create failed")
  cat("Release", tag, "created with", length(pdfs), "PDFs.\n")
}
