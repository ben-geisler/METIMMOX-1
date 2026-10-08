# ===============================================================================
# DISPLAY LABELS SHARED BY SUPPLEMENTARY TABLES S2 AND S3 (clinical effectiveness)
# ===============================================================================
# The test functions of dag_association_tests.R keep ASCII item labels
# ("TxCRP -> PFS", "Age _||_ CRP") because dag_associations.qmd, the DAG coverage
# checks and check_output_consistency.R join on them. The supplementary tables
# show the typeset forms. Needs assoc_tests.R (fmt_p).

SUPP_ARROW <- "→"
SUPP_PERP <- "⊥"

#' Typeset a DAG item or test label for the supplementary tables
#'
#' Replaces `TxCRP`/`TxTMB` by "T x CRP"/"T x TMB", `_||_` by the independence
#' sign and `->` by an arrow.
#'
#' @param v_x Character vector of item or test labels.
#' @return Character vector of the same length.
supp_label <- function(v_x) {
  v_x <- gsub("TxCRP", "T x CRP", v_x, fixed = TRUE)
  v_x <- gsub("TxTMB", "T x TMB", v_x, fixed = TRUE)
  v_x <- gsub("_||_", SUPP_PERP, v_x, fixed = TRUE)
  gsub("->", SUPP_ARROW, v_x, fixed = TRUE)
}

#' Conclusion labels of the conditional-independence table (Table S3)
#'
#' Maps the labels of `ci_support_label()` to the manuscript wording and
#' shortens the stratum note "(within X-positive patients)" to "(X-positive)".
#' An unmapped label stops, so a new conclusion cannot reach the table unlabelled.
#'
#' @param v_x Character vector of `Conclusion` values from `run_dag_ci_tests()`.
#' @return Character vector of the same length.
supp_ci_conclusion <- function(v_x) {
  v_map <- c(
    "Compatible with DAG-implied CI" = "Compatible with assumed structure",
    "Compatible with arm balance" = "Compatible with arm balance",
    "Arm imbalance in realised sample" = "Arm imbalance in realized sample",
    "CI contradicted by data" = "Incompatible with assumed structure",
    "Holds by construction" = "Holds by construction",
    "Not estimable" = "Not estimable"
  )
  v_stratum <- ifelse(grepl("\\(within .+-positive patients\\)$", v_x),
                      sub("^.*\\(within (.+)-positive patients\\)$", " (\\1-positive)", v_x),
                      "")
  v_base <- sub(" \\(within .+-positive patients\\)$", "", v_x)
  v_unmapped <- setdiff(unique(v_base), names(v_map))
  if (length(v_unmapped) > 0) {
    stop("Conclusion label without a Table S3 wording: ",
         paste(v_unmapped, collapse = "; "))
  }
  paste0(unname(v_map[v_base]), v_stratum)
}

#' Check that every test family used in a table is named in its footnote
#'
#' @param v_tests The `Test` column of the table, tests without a p-value
#'   (definitional rows) excluded by the caller.
#' @param v_families Named character vector: names are matched against the
#'   start of each test label, values are the footnote wording.
#' @param table_name Table label for the error message.
#' @return `v_families` (the wording), invisibly, after both directions pass.
supp_check_test_families <- function(v_tests, v_families, table_name) {
  v_known <- vapply(v_tests, function(t) any(startsWith(t, names(v_families))), logical(1))
  if (!all(v_known)) {
    stop("Test type missing from the ", table_name, " footnote: ",
         paste(unique(v_tests[!v_known]), collapse = "; "))
  }
  v_used <- vapply(names(v_families), function(f) any(startsWith(v_tests, f)), logical(1))
  if (!all(v_used)) {
    stop("The ", table_name, " footnote lists a test type the table does not use: ",
         paste(names(v_families)[!v_used], collapse = "; "))
  }
  invisible(v_families)
}
