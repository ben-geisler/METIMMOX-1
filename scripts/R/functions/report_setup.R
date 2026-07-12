# ===============================================================================
# REPORT SETUP HELPER
# ===============================================================================
# One call to set up a Quarto report environment: sets the common knitr chunk and
# root-directory options, loads the shared package set, applies the shared ggplot
# theme, and sources the requested analysis scripts (by number prefix, e.g. "02"
# or "10b") and function files. Replaces the boilerplate setup chunk that was
# duplicated across the economic reports.
#
# Usage (in a report setup chunk):
#   source(here::here("scripts/R/functions/report_setup.R"))
#   setup_report(sources = c("02", "03", "06", "07"),
#                funs    = c("model_fun", "calculate_outcomes", "cea_helpers"))
# ===============================================================================

#' Default package set for reports (superset of what the economic reports use)
REPORT_PACKAGES <- c(
  "here", "knitr", "kableExtra", "ggplot2", "dplyr", "tidyr", "scales",
  "survival", "flexsurv", "dampack"
)

#' Set up a Quarto report environment
#'
#' @param sources Character vector of analysis-script number prefixes to source,
#'   in order (e.g. c("02", "03", "06", "07"), "10b", "11"). Each is resolved to
#'   the matching scripts/R/analysis/<prefix>_*.R file.
#' @param funs Character vector of function-file basenames (without ".R") to
#'   source from scripts/R/functions/.
#' @param packages Character vector of packages to attach (default REPORT_PACKAGES).
#' @param set_knitr Logical; set the common knitr chunk and root-dir options.
#' @param set_theme Logical; apply the shared minimal ggplot theme.
#' @return Invisibly NULL.
#' @export
setup_report <- function(sources = character(0),
                         funs = character(0),
                         packages = REPORT_PACKAGES,
                         set_knitr = TRUE,
                         set_theme = TRUE) {

  if (!require("pacman")) install.packages("pacman")
  pacman::p_load(char = packages)

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

  analysis_dir <- here::here("scripts", "R", "analysis")
  fun_dir <- here::here("scripts", "R", "functions")

  suppressMessages(suppressWarnings({
    for (num in sources) {
      matches <- list.files(analysis_dir, pattern = paste0("^", num, "_.*\\.R$"),
                            full.names = TRUE)
      if (length(matches) == 0) {
        stop("No analysis script found for prefix '", num, "' in ", analysis_dir)
      }
      source(matches[1])
    }
    for (fn in funs) {
      fn_path <- file.path(fun_dir, paste0(fn, ".R"))
      if (!file.exists(fn_path)) {
        stop("Function file not found: ", fn_path)
      }
      source(fn_path)
    }
  }))

  invisible(NULL)
}
