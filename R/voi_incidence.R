# ===============================================================================
# ELIGIBLE POPULATION FOR POPULATION-LEVEL VALUE OF INFORMATION (issues #200, #201, #212)
# ===============================================================================
# data/external/voi_incidence.csv holds, per region, the factors of the
# eligible-population chain (CRC incidence x share metastatic x share MSS/pMMR x
# share fit for first-line doublet chemotherapy) with a source per row, and the
# derived `eligible_per_year` row with extreme low/high scenarios. The
# derivation is documented in the private repository
# (docs/manuscripts/3_value_of_information/inputs/incidence_derivation.md).
# Only aggregate published figures; no patient data.
# ===============================================================================

# Repository-relative path of the eligible-population table.
VOI_INCIDENCE_PATH <- "data/external/voi_incidence.csv"

# Columns the table must carry.
VOI_INCIDENCE_COLUMNS <- c("region", "factor", "value", "lower", "upper", "source", "doi_or_url")

#' Read the eligible patients per year by region and check the contract
#'
#' Contract: the required columns exist; the `eligible_per_year` factor has
#' exactly one row per region, every region of the table has one, and its
#' value and bounds are finite, positive and ordered (lower <= value <=
#' upper); every row names a source; every row names a DOI or URL in
#' `doi_or_url` ("derived" for the derived rows) unless its source is a
#' personal communication, which has none (issue #212).
#'
#' @param path Path of the CSV; default VOI_INCIDENCE_PATH under the
#'   repository root.
#' @param required_regions Regions that must be present.
#' @return Data frame with `region`, `value`, `lower`, `upper`, `source`,
#'   `doi_or_url` (one row per region), with attribute `path`.
read_voi_incidence <- function(path = here::here(VOI_INCIDENCE_PATH),
                               required_regions = c("Norway", "Nordic", "Europe")) {
  if (!file.exists(path)) {
    stop("File not found: ", path, ". The eligible-population table (issue #201) ",
         "is required for population-level value of information.", call. = FALSE)
  }
  df_raw <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  v_problems <- character(0)
  v_missing <- setdiff(VOI_INCIDENCE_COLUMNS, names(df_raw))
  if (length(v_missing)) {
    stop("voi_incidence.csv lacks columns: ", paste(v_missing, collapse = ", "), call. = FALSE)
  }
  v_no_source <- which(is.na(df_raw$source) | !nzchar(trimws(df_raw$source)))
  if (length(v_no_source)) v_problems <- c(v_problems, paste0("rows without a source: ",
                                                               paste(v_no_source, collapse = ", ")))
  v_personal <- grepl("personal communication", df_raw$source, ignore.case = TRUE)
  v_no_doi <- which((is.na(df_raw$doi_or_url) | !nzchar(trimws(df_raw$doi_or_url))) & !v_personal)
  if (length(v_no_doi)) v_problems <- c(v_problems, paste0(
    "rows without a DOI or URL that are not a personal communication: ",
    paste(v_no_doi, collapse = ", ")))
  df_eligible <- df_raw[df_raw$factor == "eligible_per_year", , drop = FALSE]
  v_counts <- table(factor(df_eligible$region, levels = unique(df_raw$region)))
  if (any(v_counts != 1L)) {
    v_problems <- c(v_problems, paste0("eligible_per_year must have one row per region; counts: ",
      paste(names(v_counts), v_counts, sep = " = ", collapse = ", ")))
  }
  v_absent <- setdiff(required_regions, df_eligible$region)
  if (length(v_absent)) v_problems <- c(v_problems, paste0("regions without eligible_per_year: ",
                                                           paste(v_absent, collapse = ", ")))
  for (column in c("value", "lower", "upper")) {
    v_x <- suppressWarnings(as.numeric(df_eligible[[column]]))
    v_bad <- df_eligible$region[!is.finite(v_x) | v_x <= 0]
    if (length(v_bad)) v_problems <- c(v_problems, paste0(column, " not finite and positive for: ",
                                                          paste(v_bad, collapse = ", ")))
    df_eligible[[column]] <- v_x
  }
  v_unordered <- df_eligible$region[is.finite(df_eligible$value) &
    !(df_eligible$lower <= df_eligible$value & df_eligible$value <= df_eligible$upper)]
  v_unordered <- v_unordered[!is.na(v_unordered)]
  if (length(v_unordered)) v_problems <- c(v_problems, paste0("bounds not ordered for: ",
                                                              paste(v_unordered, collapse = ", ")))
  if (length(v_problems)) {
    stop("voi_incidence.csv violates its contract:\n",
         paste0("  - ", v_problems, collapse = "\n"), call. = FALSE)
  }
  df_out <- df_eligible[c("region", "value", "lower", "upper", "source", "doi_or_url")]
  rownames(df_out) <- df_out$region
  attr(df_out, "path") <- path
  df_out
}

#' Eligible patients per year of one region
#'
#' @param df_incidence Output of read_voi_incidence().
#' @param region Region name.
#' @param bound "value", "lower" or "upper".
#' @return Numeric scalar.
voi_eligible_per_year <- function(df_incidence, region, bound = c("value", "lower", "upper")) {
  bound <- match.arg(bound)
  if (!region %in% df_incidence$region) stop("No eligible_per_year row for region ", region)
  df_incidence[[bound]][df_incidence$region == region]
}
