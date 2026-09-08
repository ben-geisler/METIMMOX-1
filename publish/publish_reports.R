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

if (!flag("--no-render")) {
  for (f in report_files) { render(f, "pdf"); render(f, "gfm") }
  for (f in vignette_files) render(f, "html")
}

# ---- Assemble _site/ ----------------------------------------------------------
unlink(site_dir, recursive = TRUE)
dir.create(file.path(site_dir, "technical"), recursive = TRUE, showWarnings = FALSE)
pdfs <- sub("\\.qmd$", ".pdf", report_files)
pdfs <- pdfs[file.exists(pdfs)]
for (p in pdfs) file.copy(p, file.path(site_dir, sub("^reports/", "", p)), overwrite = TRUE)

commit <- tryCatch(trimws(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE)), error = function(e) "unknown")
row <- function(p) {
  rel <- sub("^reports/", "", p)
  sprintf('<li><a href="%s">%s</a> <small>(%s, %s)</small></li>', rel,
          tools::file_path_sans_ext(basename(p)),
          format(file.info(p)$mtime, "%Y-%m-%d"), format(round(file.info(p)$size / 1024), big.mark = ","))
}
index <- c(
  "<!doctype html><html><head><meta charset='utf-8'><title>METIMMOX-1 reports</title>",
  "<style>body{font-family:system-ui,sans-serif;max-width:48rem;margin:3rem auto;line-height:1.5}</style></head><body>",
  "<h1>METIMMOX-1 economic evaluation: rendered reports</h1>",
  sprintf("<p>Rendered %s from commit <code>%s</code>. Sizes in KB.</p>", format(Sys.Date()), commit),
  "<h2>Reports</h2><ul>", vapply(pdfs[!grepl("/technical/", pdfs)], row, ""), "</ul>",
  "<h2>Technical documentation</h2><ul>", vapply(pdfs[grepl("/technical/", pdfs)], row, ""), "</ul>",
  "</body></html>"
)
writeLines(index, file.path(site_dir, "index.html"))
writeLines("", file.path(site_dir, ".nojekyll"))
cat(sprintf("_site/ assembled: %d PDFs + index.html\n", length(pdfs)))

# ---- Publish ------------------------------------------------------------------
if (flag("--push")) {
  wt <- file.path(tempdir(), "gh-pages-worktree")
  unlink(wt, recursive = TRUE)
  has_branch <- system2("git", c("show-ref", "--quiet", "refs/heads/gh-pages")) == 0
  if (has_branch) {
    system2("git", c("worktree", "add", shQuote(wt), "gh-pages"))
  } else {
    system2("git", c("worktree", "add", "--detach", shQuote(wt)))
    system2("git", c("-C", shQuote(wt), "checkout", "--orphan", "gh-pages"))
    system2("git", c("-C", shQuote(wt), "rm", "-rfq", "."))
  }
  file.copy(list.files(site_dir, full.names = TRUE, all.files = TRUE, no.. = TRUE), wt, recursive = TRUE, overwrite = TRUE)
  system2("git", c("-C", shQuote(wt), "add", "-A"))
  system2("git", c("-C", shQuote(wt), "commit", "-qm", shQuote(paste("Publish reports from", commit))))
  status <- system2("git", c("-C", shQuote(wt), "push", "origin", "gh-pages"))
  system2("git", c("worktree", "remove", "--force", shQuote(wt)))
  if (status != 0) stop("git push to gh-pages failed")
  cat("Pushed _site/ to gh-pages. Enable Pages on that branch in the repository settings.\n")
}

if (flag("--release")) {
  tag <- paste0("reports-", format(Sys.Date(), "%Y%m%d"))
  status <- system2("gh", c("release", "create", tag, shQuote(list.files(site_dir, "\\.pdf$", recursive = TRUE, full.names = TRUE)),
                            "--title", shQuote(paste("Reports", format(Sys.Date()))),
                            "--notes", shQuote(paste("Rendered from commit", commit))))
  if (status != 0) stop("gh release create failed")
  cat("Release", tag, "created with", length(pdfs), "PDFs.\n")
}
