# Public-mirror safety scan (issue #189). Synthetic fixtures check every rule;
# then the files this tree would publish (R/, analysis/, tests/, reports/ and
# the root files of R/public_safety.R) are scanned for real, in the private
# working repository and in the public mirror alike. Needs no trial data.
# Run from the repository root: Rscript tests/test_public_safety.R
source("R/public_safety.R")

#' Write a synthetic tree and return its root
#'
#' @param files Named list: relative path -> character vector of lines.
#' @return Path of the temporary root directory.
make_tree <- function(files) {
  root <- tempfile("public_tree_")
  for (path in names(files)) {
    dir.create(dirname(file.path(root, path)), recursive = TRUE, showWarnings = FALSE)
    writeLines(files[[path]], file.path(root, path))
  }
  root
}

clean <- list(
  "R/model.R" = c("# Model code", "f <- function(x) x + 1"),
  "analysis/01_data_prep.R" = c("read.csv('data/sensitive/excluded_ids.csv')"),
  "reports/CEA.qmd" = c("---", "title: CEA", "---", "Dated 2026-10-08; R 4.3.2-1."),
  "DESCRIPTION" = c("Authors@R: person(email = \"8674152+ben-geisler@users.noreply.github.com\")"),
  "README.md" = c("Co-Authored-By: Claude <noreply@anthropic.com>")
)
root <- make_tree(clean)
stopifnot(nrow(scan_public_tree(root, names(clean))) == 0)
cat("PASS: a clean tree has no violations (no-reply addresses allowed)\n")

#' Rules violated by one extra file added to the clean tree
rules_for <- function(path, lines = "x") {
  files <- c(clean, setNames(list(lines), path))
  root <- make_tree(files)
  unique(scan_public_tree(root, names(files))$rule)
}
# The forbidden strings are assembled at run time, so that this file itself
# passes the scan of the real tree below.
stopifnot(
  identical(rules_for("R/ids.R", paste0("exclude <- c('", "9", "-0", "42')")), "trial_id"),
  identical(rules_for("tests/t.R", paste0("# C:", "\\", "Users", "\\", "someone", "\\", "data.rds")), "machine_path"),
  identical(rules_for("README2.md", paste0("Contact: someone", "@", "hospital.no")), "email"),
  identical(rules_for("reports/x.qmd", paste0("Two patients were de-", "randomized.")), "protocol_deviation"),
  identical(rules_for("analysis/data.rds"), "file_type"),
  identical(rules_for("reports/figs/plot.png"), "file_type"),
  identical(rules_for("reports/table.csv", "a,b"), "file_type"),
  identical(rules_for("AGENTS.md", "notes"), "private_file"),
  identical(rules_for("data/tidy/notes.md", "x"), "private_folder"),
  identical(rules_for("outputs/reports/CEA.md", "x"), "private_folder")
)
cat("PASS: every rule catches its synthetic violation\n")

# Which files belong in the mirror: the public folders and root files only.
root <- make_tree(c(clean, list("docs/manuscript.md" = "x", "AGENTS.md" = "x",
                                "outputs/tables/t.csv" = "a", "make.R" = "x")))
stopifnot(identical(public_mirror_files(root), sort(c(names(clean), "make.R"), method = "radix")))
cat("PASS: public_mirror_files() keeps only the mirrored folders and root files\n")

# The real tree: what this repository would publish must pass the scan.
v_public <- public_mirror_files(".")
stopifnot(length(v_public) > 100, all(c("README.md", "make.R") %in% v_public),
          !any(grepl("^(data|outputs|docs|tools)/", v_public)))
df_violations <- scan_public_tree(".", v_public)
if (nrow(df_violations)) {
  print(df_violations, row.names = FALSE)
  stop(nrow(df_violations), " violation(s) in the files of the public mirror")
}
cat("PASS: the", length(v_public), "files of the public mirror pass the safety scan\n")
