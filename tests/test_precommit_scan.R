# Pre-commit public-safety scan (tools/hooks/pre_commit_scan.R; private working
# repository only). On a temporary git repository: a clean staged public file
# passes; a staged public file with a trial ID is refused, also when the working
# copy has been cleaned after staging; a staged file of a forbidden type under a
# public folder is refused; a violation outside the public paths is ignored; a
# commit that stages nothing public passes. Skips when the hook script is absent
# (the public mirror). Needs no trial data. Run from the repository root.
hook <- file.path("tools", "hooks", "pre_commit_scan.R")
if (!file.exists(hook)) {
  cat("SKIP: tools/hooks/pre_commit_scan.R is not in this tree (public mirror)\n")
  quit(save = "no", status = 0L)
}
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")

repo <- tempfile("precommit_")
dir.create(repo)
#' Run git in the temporary repository and stop on a nonzero status.
git <- function(...) {
  out <- suppressWarnings(system2("git", c("-C", shQuote(repo), ...), stdout = TRUE, stderr = TRUE))
  stopifnot(is.null(attr(out, "status")) || attr(out, "status") == 0L)
  invisible(out)
}
#' Write a file into the temporary repository, creating its folder.
put <- function(rel, lines) {
  dir.create(dirname(file.path(repo, rel)), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, file.path(repo, rel))
}
#' Run the hook script with the temporary repository as working directory;
#' returns its exit status with the output as attribute "output".
run_scan <- function() {
  old <- setwd(repo)
  on.exit(setwd(old), add = TRUE)
  out <- suppressWarnings(system2(rscript, shQuote(file.path(repo, hook)), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  structure(if (is.null(status)) 0L else as.integer(status), output = paste(out, collapse = "\n"))
}

git("init", "-q", "-b", "main")
git("config", "user.name", "hook test")
git("config", "user.email", "hook@users.noreply.github.com")
put("R/public_safety.R", readLines(file.path("R", "public_safety.R"), warn = FALSE))
put(hook, readLines(hook, warn = FALSE))
git("add", "-A")
git("commit", "-q", "-m", "init")

# Forbidden strings are assembled at run time so that this file passes the scan.
trial_id <- paste0("1-0", "42")

# 1. A clean staged public file passes.
put("R/ok.R", "x <- 1")
git("add", "R/ok.R")
stopifnot(run_scan() == 0L)
cat("PASS: a clean staged public file passes\n")

# 2. A trial ID in a staged public file is refused, with the rule named ...
put("R/bad.R", c("# excluded:", trial_id))
git("add", "R/bad.R")
res <- run_scan()
stopifnot(res == 1L, grepl("trial_id", attr(res, "output"), fixed = TRUE),
          grepl("R/bad.R", attr(res, "output"), fixed = TRUE))
# ... also when the working copy has been cleaned after staging (the staged
# content is what the commit would contain).
put("R/bad.R", "x <- 2")
stopifnot(run_scan() == 1L)
cat("PASS: a staged trial ID is refused, from the staged content\n")
git("reset", "-q", "HEAD", "--", "R/bad.R")
unlink(file.path(repo, "R", "bad.R"))

# 3. A forbidden file type under a public folder is refused.
put("tests/fixture.RData", "not really a workspace")
git("add", "-f", "tests/fixture.RData")
res <- run_scan()
stopifnot(res == 1L, grepl("file_type", attr(res, "output"), fixed = TRUE))
cat("PASS: a forbidden file type under a public folder is refused\n")
git("reset", "-q", "HEAD", "--", "tests/fixture.RData")

# 4. The same violation outside the public paths is not the hook's business.
put("docs/notes.md", trial_id)
git("add", "docs/notes.md")
stopifnot(run_scan() == 0L)
cat("PASS: a violation outside the public paths is ignored\n")
git("reset", "-q", "HEAD", "--", "docs/notes.md")

# 5. Nothing public staged: passes without scanning.
git("reset", "-q")
stopifnot(run_scan() == 0L)
cat("PASS: a commit that stages nothing public passes\n")

unlink(repo, recursive = TRUE)
