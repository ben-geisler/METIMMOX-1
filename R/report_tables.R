# ===============================================================================
# FORMAT-AWARE TABLE HELPERS
# ===============================================================================
# Every report renders to PDF (LaTeX) and to GitHub-flavoured Markdown (gfm).
# kableExtra's styling functions force HTML output, which aborts the gfm pass
# ("Functions that produce HTML output found in document targeting commonmark
# output"). These wrappers apply the kableExtra call only under LaTeX and pass
# the plain knitr::kable() (a pipe table under gfm) through unchanged otherwise.
# Sourced by report_setup.R; use them in place of the kableExtra functions.
#
#   knitr::kable(df, booktabs = TRUE, caption = "...") %>%
#     tbl_style(latex_options = "HOLD_position") %>%
#     tbl_column_spec(2, width = "10cm") %>%
#     tbl_footnote(general = "Note text")
#
# Under gfm: tbl_style/tbl_column_spec/tbl_row_spec/tbl_header_above/
# tbl_pack_rows/tbl_landscape return the table unchanged (group headers and
# spanners are lost); tbl_footnote appends the note as a paragraph below the
# table so the text survives in the Markdown render.
# ===============================================================================

#' Table format for knitr::kable(): "latex" under PDF, "pipe" otherwise
report_table_format <- function() {
  if (knitr::is_latex_output()) "latex" else "pipe"
}

.tbl_latex_only <- function(fun) {
  function(kable_input, ...) {
    if (knitr::is_latex_output()) fun(kable_input, ...) else kable_input
  }
}

tbl_style        <- .tbl_latex_only(kableExtra::kable_styling)
tbl_column_spec  <- .tbl_latex_only(kableExtra::column_spec)
tbl_row_spec     <- .tbl_latex_only(kableExtra::row_spec)
tbl_header_above <- .tbl_latex_only(kableExtra::add_header_above)
tbl_pack_rows    <- .tbl_latex_only(kableExtra::pack_rows)
tbl_landscape    <- .tbl_latex_only(kableExtra::landscape)

#' Footnote that survives the Markdown render
tbl_footnote <- function(kable_input, general = NULL, general_title = "Note: ", ...) {
  if (knitr::is_latex_output()) {
    return(kableExtra::footnote(kable_input, general = general,
                                general_title = general_title, ...))
  }
  if (is.null(general)) return(kable_input)
  knitr::asis_output(paste0(
    paste(kable_input, collapse = "\n"), "\n\n*",
    trimws(general_title), "* ", paste(general, collapse = " "), "\n"
  ))
}
