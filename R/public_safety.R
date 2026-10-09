# ===============================================================================
# PUBLIC MIRROR: WHAT IS PUBLISHED AND THE SAFETY SCAN
# ===============================================================================
# The public repository METIMMOX-1 mirrors the folders and root files listed
# below from a private working repository, which also holds the confidential
# trial data, the rendered outputs and the manuscripts (tools/sync_public.R
# there). Every mirrored file is scanned before it is committed: only plain-text
# source files are allowed, and none may contain a trial patient ID, a local
# machine path, a personal e-mail address or the private agent notes. The scan
# cannot recognise patient-level facts written as prose, so those must never
# be typed into these sources in the first place.
# ===============================================================================

# .github holds the GitHub Actions workflow of the public repository (issue
# #194); it is mirrored and scanned like the code.
PUBLIC_MIRROR_DIRS <- c("R", "analysis", "tests", "reports", ".github")
PUBLIC_MIRROR_FILES <- c("README.md", "LICENSE", "DESCRIPTION", "renv.lock",
                         "METIMMOX-1.Rproj", "make.R", "CITATION.cff")
# File types allowed in the mirror (lower case; "" for DESCRIPTION and LICENSE;
# yml and yaml for the workflow files under .github).
PUBLIC_ALLOWED_EXTENSIONS <- c("r", "qmd", "rmd", "md", "lock", "rproj", "cff",
                               "gitignore", "yml", "yaml", "")
# Text rules: a trial ID is the site digit, a hyphen and a three-digit patient
# number starting with 0.
PUBLIC_TEXT_RULES <- c(
  trial_id = "\\b[0-9]-0[0-9]{2}\\b",
  machine_path = "[A-Za-z]:[\\\\/]Users[\\\\/]|/Users/[A-Za-z]|/home/[A-Za-z]",
  email = "[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\\.[A-Za-z]{2,}",
  protocol_deviation = "de-?randomi[sz]"
)
# E-mail addresses that may appear (GitHub and tool no-reply addresses).
PUBLIC_ALLOWED_EMAIL <- "noreply"
# Files exempt from one text rule: renv.lock records the CRAN maintainer
# addresses of every package.
PUBLIC_RULE_EXEMPTIONS <- list(email = "renv.lock")
PUBLIC_FORBIDDEN_NAMES <- c("AGENTS.md", "CLAUDE.md", "CLAUDE.local.md", ".Rhistory", ".RData")

#' Files of a tree that belong in the public mirror
#'
#' @param root Repository root (the private working repository or the mirror).
#' @return Sorted repository-relative paths of the tracked files in
#'   PUBLIC_MIRROR_DIRS and PUBLIC_MIRROR_FILES (all files there when `root` is
#'   not a git repository).
public_mirror_files <- function(root) {
  # core.quotePath=false: one unquoted path per line, also for non-ASCII names.
  v_files <- suppressWarnings(system2("git", c("-C", shQuote(root), "-c", "core.quotePath=false",
                                               "ls-files"), stdout = TRUE, stderr = FALSE))
  if (!is.null(attr(v_files, "status")) || !length(v_files)) {
    v_files <- list.files(root, recursive = TRUE, all.files = TRUE)
  }
  v_top <- sub("/.*$", "", v_files)
  v_keep <- (v_top %in% PUBLIC_MIRROR_DIRS & grepl("/", v_files)) | v_files %in% PUBLIC_MIRROR_FILES
  sort(unique(v_files[v_keep]), method = "radix")
}

#' Scan files for content that must never be public
#'
#' @param root Root directory the paths are relative to.
#' @param v_files Repository-relative paths to scan.
#' @return Data frame of violations (file, line, rule, text); zero rows when
#'   the files are safe to publish.
scan_public_tree <- function(root, v_files) {
  violation <- function(file, line, rule, text)
    data.frame(file = file, line = line, rule = rule, text = substr(text, 1, 120),
               stringsAsFactors = FALSE)
  v_out <- list()
  for (file in v_files) {
    name <- basename(file)
    ext <- tolower(sub("^.*\\.", "", name))
    if (!grepl(".", name, fixed = TRUE)) ext <- ""
    if (name == ".gitignore") ext <- "gitignore"
    if (name %in% PUBLIC_FORBIDDEN_NAMES) {
      v_out[[length(v_out) + 1]] <- violation(file, NA_integer_, "private_file", name)
      next
    }
    if (sub("/.*$", "", file) %in% c("data", "outputs", "docs", "tools")) {
      v_out[[length(v_out) + 1]] <- violation(file, NA_integer_, "private_folder", file)
      next
    }
    if (!ext %in% PUBLIC_ALLOWED_EXTENSIONS) {
      v_out[[length(v_out) + 1]] <- violation(file, NA_integer_, "file_type", ext)
      next
    }
    v_lines <- readLines(file.path(root, file), warn = FALSE, encoding = "UTF-8")
    for (rule in names(PUBLIC_TEXT_RULES)) {
      if (file %in% PUBLIC_RULE_EXEMPTIONS[[rule]]) next
      v_hits <- grep(PUBLIC_TEXT_RULES[[rule]], v_lines, perl = TRUE)
      if (rule == "email") {
        v_hits <- v_hits[vapply(v_hits, function(i) {
          v_addresses <- regmatches(v_lines[i], gregexpr(PUBLIC_TEXT_RULES[["email"]], v_lines[i], perl = TRUE))[[1]]
          any(!grepl(PUBLIC_ALLOWED_EMAIL, v_addresses, fixed = TRUE))
        }, logical(1))]
      }
      for (i in v_hits) v_out[[length(v_out) + 1]] <- violation(file, i, rule, v_lines[i])
    }
  }
  if (!length(v_out)) return(violation(character(0), integer(0), character(0), character(0)))
  do.call(rbind, v_out)
}
