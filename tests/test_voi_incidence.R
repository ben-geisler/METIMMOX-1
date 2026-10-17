# ===============================================================================
# TEST: eligible-population reader for population VOI (issues #200, #201, #212)
# ===============================================================================
# 1. A synthetic table that meets the contract is read into one row per region
#    with value and bounds; voi_eligible_per_year() returns them.
# 2. Every contract violation stops with a message naming it: missing file,
#    missing column, two eligible_per_year rows for a region, a region without
#    one, a non-positive or non-finite value, unordered bounds, a row without a
#    source, a row without a DOI or URL that is not a personal communication
#    (a personal communication without one passes; #212).
# 3. Live (when data/external/voi_incidence.csv exists; it is not in the
#    public mirror): the sourced table meets the contract, Norway's value lies
#    within its bounds and below the former placeholder of 1,500, the Nordic
#    and European values exceed Norway's, and every eligible_per_year value
#    and bound equals the product of its factor rows within 0.5%.
#
# Run: "C:\Program Files\R\R-4.3.2\bin\x64\Rscript.exe" tests/test_voi_incidence.R
# ===============================================================================

source(here::here("R", "voi_incidence.R"))

#' Stop unless `expr` raises an error whose message matches `pattern`.
expect_error_matching <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  if (!inherits(error, "error") || !grepl(pattern, conditionMessage(error))) {
    stop("Expected an error matching '", pattern, "', got: ",
         if (inherits(error, "error")) conditionMessage(error) else "no error")
  }
  invisible(TRUE)
}

df_valid <- data.frame(
  region = c("Norway", "Norway", "Nordic", "Europe"),
  factor = c("crc_incidence", "eligible_per_year", "eligible_per_year", "eligible_per_year"),
  value = c(5000, 500, 2000, 40000), lower = c(4800, 400, 1500, 30000),
  upper = c(5100, 800, 3500, 60000), source = c("Registry", "Derived", "Derived", "Derived"),
  doi_or_url = c("https://example.org", "derived", "derived", "derived"),
  stringsAsFactors = FALSE)
path <- tempfile(fileext = ".csv")
#' Write a variant of the synthetic table and return its path.
write_variant <- function(df) {
  utils::write.csv(df, path, row.names = FALSE)
  path
}

df_read <- read_voi_incidence(write_variant(df_valid))
stopifnot(identical(df_read$region, c("Norway", "Nordic", "Europe")),
          voi_eligible_per_year(df_read, "Norway") == 500,
          voi_eligible_per_year(df_read, "Norway", "lower") == 400,
          voi_eligible_per_year(df_read, "Europe", "upper") == 60000)
expect_error_matching(voi_eligible_per_year(df_read, "Mars"), "No eligible_per_year row")
cat("PASS: a valid table gives one eligible_per_year row per region\n")

expect_error_matching(read_voi_incidence(file.path(tempdir(), "absent.csv")), "File not found")
expect_error_matching(read_voi_incidence(write_variant(df_valid[, names(df_valid) != "upper"])),
                      "lacks columns: upper")
expect_error_matching(read_voi_incidence(write_variant(rbind(df_valid, df_valid[2, ]))),
                      "one row per region")
expect_error_matching(read_voi_incidence(write_variant(df_valid[-4, ])),
                      "regions without eligible_per_year: Europe")
df_bad <- df_valid; df_bad$value[3] <- 0
expect_error_matching(read_voi_incidence(write_variant(df_bad)), "value not finite and positive for: Nordic")
df_bad <- df_valid; df_bad$lower[2] <- NA
expect_error_matching(read_voi_incidence(write_variant(df_bad)), "lower not finite and positive for: Norway")
df_bad <- df_valid; df_bad$upper[4] <- 35000
expect_error_matching(read_voi_incidence(write_variant(df_bad)), "bounds not ordered for: Europe")
df_bad <- df_valid; df_bad$source[1] <- ""
expect_error_matching(read_voi_incidence(write_variant(df_bad)), "rows without a source: 1")
df_bad <- df_valid; df_bad$doi_or_url[1] <- ""
expect_error_matching(read_voi_incidence(write_variant(df_bad)),
                      "rows without a DOI or URL that are not a personal communication: 1")
df_personal <- df_valid; df_personal$doi_or_url[1] <- ""
df_personal$source[1] <- "Personal communication, A. H. Ree, 15 October 2026"
stopifnot(voi_eligible_per_year(read_voi_incidence(write_variant(df_personal)), "Norway") == 500)
unlink(path)
cat("PASS: every contract violation stops with a message naming it\n")

live_path <- here::here(VOI_INCIDENCE_PATH)
if (!file.exists(live_path)) {
  cat("SKIP: live eligible-population table absent (", VOI_INCIDENCE_PATH, ")\n", sep = "")
} else {
  df_live <- read_voi_incidence(live_path)
  v_norway <- unlist(df_live["Norway", c("lower", "value", "upper")])
  stopifnot(all(c("Norway", "Nordic", "Europe") %in% df_live$region),
            v_norway[["lower"]] <= v_norway[["value"]], v_norway[["value"]] <= v_norway[["upper"]],
            v_norway[["value"]] < 1500,
            voi_eligible_per_year(df_live, "Nordic") > v_norway[["value"]],
            voi_eligible_per_year(df_live, "Europe") > voi_eligible_per_year(df_live, "Nordic"))
  # Each eligible_per_year value and bound is the product of its region's
  # factor rows (bounds with bounds; a missing bound falls back to the value).
  # Tolerance 0.5%: the rounded rows, and those derived from unrounded shares
  # before #212.
  df_factors <- utils::read.csv(live_path, stringsAsFactors = FALSE)
  v_chain <- c("crc_incidence", "share_synchronous_metastatic", "share_metachronous",
               "share_mss_pmmr", "share_fit_first_line_doublet")
  for (region in df_live$region) for (bound in c("value", "lower", "upper")) {
    df_region <- df_factors[df_factors$region == region, ]
    v_x <- vapply(v_chain, function(k) {
      x <- suppressWarnings(as.numeric(df_region[df_region$factor == k, bound]))
      if (length(x) != 1L) stop("Factor ", k, " missing or duplicated for ", region)
      if (is.na(x)) x <- as.numeric(df_region[df_region$factor == k, "value"])
      x
    }, numeric(1))
    product <- v_x[[1]] * (v_x[[2]] + v_x[[3]]) * v_x[[4]] * v_x[[5]]
    stated <- voi_eligible_per_year(df_live, region, bound)
    if (abs(stated / product - 1) > 0.005) {
      stop(sprintf("eligible_per_year %s %s = %g, but the factor chain gives %.1f",
                   region, bound, stated, product))
    }
  }
  cat("PASS: live table: every eligible_per_year row equals its factor chain\n")
  cat(sprintf("PASS: live table: Norway %g (%g to %g) eligible patients per year\n",
              v_norway[["value"]], v_norway[["lower"]], v_norway[["upper"]]))
}

cat("PASS: test_voi_incidence.R\n")
