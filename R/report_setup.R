# ===============================================================================
# REPORT SETUP HELPER
# ===============================================================================
# One call to set up a Quarto report environment: sets the common knitr chunk and
# root-directory options, loads the shared package set, applies the shared ggplot
# theme, and sources the requested analysis scripts (by number prefix, e.g. "02"
# or "08b") and function files. Replaces the boilerplate setup chunk that was
# duplicated across the economic reports.
#
# Usage (in a report setup chunk):
#   source(here::here("R/report_setup.R"))
#   setup_report(sources = c("02", "03", "04", "05"),
#                funs    = c("model_fun", "calculate_outcomes", "cea_helpers"))
# ===============================================================================

# Format-aware table wrappers (tbl_style(), tbl_footnote(), ...): kableExtra under
# LaTeX only, so every report also renders to gfm.
source(here::here("R", "report_tables.R"))

#' Default package set for reports (superset of what the economic reports use)
REPORT_PACKAGES <- c(
  "here", "knitr", "kableExtra", "ggplot2", "dplyr", "tidyr", "scales",
  "survival", "flexsurv", "dampack"
)

#' Set up a Quarto report environment
#'
#' @param sources Character vector of analysis-script number prefixes to source,
#'   in order (e.g. c("02", "03", "04", "05"), "08b", "09"). Each is resolved to
#'   the matching analysis/<prefix>_*.R file.
#' @param funs Character vector of function-file basenames (without ".R") to
#'   source from R/.
#' @param packages Character vector of packages to attach (default REPORT_PACKAGES).
#' @param set_knitr Logical; set the common knitr chunk and root-dir options.
#' @param set_theme Logical; apply the shared minimal ggplot theme.
#' @param quiet_sources Logical; capture and discard console output emitted while
#'   sourcing analysis scripts and function files. Messages and warnings are
#'   suppressed regardless.
#' @return Invisibly NULL.
#' @export
setup_report <- function(sources = character(0),
                         funs = character(0),
                         packages = REPORT_PACKAGES,
                         set_knitr = TRUE,
                         set_theme = TRUE,
                         quiet_sources = FALSE) {

  old_sampling_option <- options(metimmox.sampling_allow_regenerate = FALSE)
  on.exit(options(old_sampling_option), add = TRUE)

  if (!is.logical(quiet_sources) || length(quiet_sources) != 1L ||
      is.na(quiet_sources)) {
    stop("quiet_sources must be TRUE or FALSE.")
  }

  if (!require("pacman")) install.packages("pacman")
  # kableExtra sets knitr.table.format to "html" on load whenever the output is
  # not LaTeX, which turns every knitr::kable() into an HTML table in the gfm
  # render. Disable that and pin the default to pipe tables outside LaTeX.
  options(kableExtra.auto_format = FALSE)
  pacman::p_load(char = packages)
  options(knitr.table.format = report_table_format())

  if (isTRUE(set_knitr)) {
    knitr::opts_chunk$set(
      echo = FALSE, warning = FALSE, message = FALSE,
      fig.align = "center", fig.width = 8, fig.height = 6,
      dpi = 300, out.width = "100%"
    )
    knitr::opts_knit$set(root.dir = here::here())
  }

  if (isTRUE(set_theme)) {
    ggplot2::theme_set(
      ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 14, face = "bold"),
          plot.subtitle = ggplot2::element_text(size = 12),
          legend.position = "bottom"
        )
    )
  }

  analysis_dir <- here::here("analysis")
  fun_dir <- here::here("R")

  source_report_file <- function(path) {
    if (quiet_sources) {
      invisible(utils::capture.output(source(path)))
    } else {
      source(path)
    }
  }

  suppressMessages(suppressWarnings({
    for (num in sources) {
      matches <- list.files(analysis_dir, pattern = paste0("^", num, "_.*\\.R$"),
                            full.names = TRUE)
      if (length(matches) == 0) {
        stop("No analysis script found for prefix '", num, "' in ", analysis_dir)
      }
      source_report_file(matches[1])
    }
    for (fn in funs) {
      fn_path <- file.path(fun_dir, paste0(fn, ".R"))
      if (!file.exists(fn_path)) {
        stop("Function file not found: ", fn_path)
      }
      source_report_file(fn_path)
    }
  }))

  invisible(NULL)
}
