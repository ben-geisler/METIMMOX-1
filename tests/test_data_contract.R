# Structural test of the raw export contract on a synthetic export (issue #174).
# Needs no trial data.
source("R/data_contract.R")
contract <- export_contract()
#' Assert that evaluating `expr` raises an error whose message contains the
#' literal text `pattern`.
fails <- function(expr, pattern) {
  e <- tryCatch({force(expr); NULL}, error = identity)
  stopifnot(inherits(e, "error"), grepl(pattern, conditionMessage(e), fixed = TRUE))
}

# Synthetic export: 541 columns, contract headers at their positions, filler
# elsewhere; 74 rows (36 control, 38 experimental) with valid values.
n <- contract$n_rows
width <- max(as.integer(c(names(contract$headers), names(contract$top_labels))))
v_hdr <- sprintf("Filler...%d", seq_len(width))
v_hdr[as.integer(names(contract$headers))] <- contract$headers
v_top <- rep("", width)
v_top[as.integer(names(contract$top_labels))] <- contract$top_labels
df_syn <- as.data.frame(setNames(replicate(width, rep(NA_real_, n), simplify = FALSE), v_hdr),
                     check.names = FALSE)
#' Convert a Date vector to the POSIXct (UTC) type that readxl returns for dates.
posix <- function(d) as.POSIXct(format(d), tz = "UTC")
v_inclusion <- as.Date("2019-01-01") + seq_len(n) * 7
v_no_ct1 <- c(3, 40)
df_syn$ID <- sprintf("1-%03d", seq_len(n))
df_syn$`Study1/control0` <- rep(0:1, c(36, 38))
df_syn$`Date of inclusion` <- posix(v_inclusion)
df_syn$`Days until last evaluation` <- 100
df_syn$`Days until progression` <- 120
df_syn$`Progression exit` <- rep(0:1, length.out = n)
df_syn$`Days until death/last follow up` <- 300
df_syn$Death <- rep(1:0, length.out = n)
df_syn$Mutation <- rep(c("BRAF", "KRAS", "NRAS", "wt"), length.out = n)
df_syn$TMB <- ifelse(seq_len(n) %% 10 == 0, "Missing", "5.5")
df_syn$`Sex 0female` <- 0
df_syn$Age <- 60
df_syn$`Date...114` <- posix(v_inclusion - 7)
df_syn$`TL LD...121` <- 50
df_syn$`Date...122` <- posix(v_inclusion + 56)
df_syn$`Date...122`[v_no_ct1] <- NA
df_syn$`TL LD...130` <- ifelse(seq_len(n) %in% v_no_ct1, "NA", "40")
df_syn$`Date...468` <- posix(v_inclusion + 5)
df_syn$`Date...530` <- posix(v_inclusion + 35)
df_syn$`CRP...478` <- 3
df_syn$`CRP...540` <- c(NA, rep(8, n - 1))

stopifnot(identical(check_export_contract(df_syn, v_top), df_syn))
stopifnot(identical(parse_export_numeric(c("1.5", "NA", NA, " 2 "), "x", "NA"),
                    c(1.5, NA, NA, 2)))
fails(parse_export_numeric(c("1", "NE"), "x", "NA"), "x has 1 non-numeric entry")

#' Apply `f` to a copy of the synthetic export and return the modified copy.
modify <- function(f) { x <- df_syn; f(x) }
# Layout: an inserted column shifts every positional name.
df_shifted <- df_syn[, c(1:99, 1, 100:width)]
v_tail_hdr <- v_hdr[100:width]
v_repaired <- grepl("[.]{3}[0-9]+$", v_tail_hdr)  # readxl renumbers these
v_tail_hdr[v_repaired] <- paste0(sub("[.]{3}[0-9]+$", "", v_tail_hdr[v_repaired]), "...",
                             (101:(width + 1))[v_repaired])
names(df_shifted) <- c(v_hdr[1:99], "Inserted", v_tail_hdr)
fails(check_export_contract(df_shifted, c(v_top[1:99], "", v_top[100:width])),
      "column 114 is 'Filler...114', expected 'Date...114'")
# A whole lab block inserted keeps "CRP...540" but moves the top labels.
fails(check_export_contract(df_syn, replace(v_top, 496, "SEQ1 V2")),
      "top header at column 496 is 'SEQ1 V2', expected 'SEQ1 V3'")
# ... and the lab date of the block no longer matches cycle 3 day 1.
fails(check_export_contract(modify(function(x) { x$`Date...530` <- posix(v_inclusion + 21); x }), v_top),
      "Date...530: median 21 days from inclusion")
fails(check_export_contract(modify(function(x) { x$`Date...530` <- posix(v_inclusion + 52); x }), v_top),
      "Date...530: median 52 days from inclusion")
fails(check_export_contract(modify(function(x) { x$`Date...468` <- posix(v_inclusion); x }), v_top),
      "Date...468: median 0 days from inclusion")
# Identifiers and row count.
fails(check_export_contract(modify(function(x) { x$ID[2] <- x$ID[1]; x }), v_top),
      "ID must be non-missing and unique (0 missing, 1 duplicated)")
fails(check_export_contract(df_syn[-1, ], v_top), "73 rows after filtering, expected 74")
# Types and ranges.
fails(check_export_contract(modify(function(x) { x$`TL LD...130`[1] <- "NE"; x }), v_top),
      "TL LD...130 has 1 non-numeric entry")
fails(check_export_contract(modify(function(x) { x$`TL LD...121`[1] <- 0; x }), v_top),
      "TL LD...121 (baseline target-lesion sum) must be positive")
fails(check_export_contract(modify(function(x) { x$`CRP...540`[2] <- -1; x }), v_top),
      "CRP...540 must be nonnegative numeric or missing")
fails(check_export_contract(modify(function(x) { x$`CRP...478` <- "3"; x }), v_top),
      "CRP...478 must be nonnegative numeric or missing")
fails(check_export_contract(modify(function(x) { x$TMB[1] <- "unknown"; x }), v_top),
      "TMB has 1 non-numeric entry")
fails(check_export_contract(modify(function(x) { x$Mutation[1] <- "HRAS"; x }), v_top),
      "Mutation must be one of")
fails(check_export_contract(modify(function(x) { x$Death[1] <- 2; x }), v_top),
      "Death: expected 0/1")
fails(check_export_contract(modify(function(x) { x$`Days until progression`[1] <- NA; x }), v_top),
      "Days until progression: expected nonnegative finite numeric values")
fails(check_export_contract(modify(function(x) { x$`Days until last evaluation`[1] <- 400; x }), v_top),
      "last evaluation exceeds death/last follow-up time")
# Dates.
fails(check_export_contract(modify(function(x) { x$`Date...122` <- as.character(x$`Date...122`); x }), v_top),
      "Date...122 is not a date column")
fails(check_export_contract(modify(function(x) { x$`Date...122`[5] <- posix(as.Date("2030-01-01")); x }), v_top),
      "Date...122: 1 date(s) outside 2018-01-01 to 2023-12-31")
fails(check_export_contract(modify(function(x) { x$`Date...122`[5] <- posix(v_inclusion[5] + 7); x }), v_top),
      "first on-treatment CT (Date...122) outside (4, 20] weeks")
fails(check_export_contract(modify(function(x) { x$`Date...114`[5] <- posix(v_inclusion[5] + 100); x }), v_top),
      "baseline scan (Date...114) missing or outside")
fails(check_export_contract(modify(function(x) { x$`Date...122`[1] <- NA; x }), v_top),
      "TL LD...130 and Date...122 (first on-treatment CT) are not missing together")

# Biomarker counts on a synthetic derived data set with the known counts.
#' Binary factor of length `total`: `pos` ones, `n - pos` zeros, then NA for
#' the `total - n` unclassified patients.
cat_vec <- function(pos, n, total) factor(c(rep(1, pos), rep(0, n - pos), rep(NA, total - n)), levels = 0:1)
df_derived <- data.frame(
  Rx = factor(rep(c("Control arm", "Experimental arm"), c(36, 38)),
              levels = c("Control arm", "Experimental arm")),
  CRP1cat = c(cat_vec(7, 35, 36), cat_vec(17, 36, 38)),
  CRP0cat = c(cat_vec(7, 35, 36), cat_vec(8, 38, 38)),
  TMBcat = cat_vec(31, 69, 74),
  Mutation = "KRAS",
  TLRcat = cat_vec(44, 68, 74))
stopifnot(identical(check_export_counts(df_derived), df_derived))
fails(check_export_counts(transform(df_derived, CRP1cat = c(cat_vec(7, 35, 36), cat_vec(16, 36, 38)))),
      "crp_week4 positives 7/35 vs 16/36, expected 7/35 vs 17/36")
fails(check_export_counts(transform(df_derived, Mutation = replace(Mutation, 74, "BRAF"))),
      "tmb_braf positives 32/70, expected 31/69")
fails(check_export_counts(transform(df_derived, TLRcat = cat_vec(44, 67, 74))),
      "tlr positives 44/67, expected 44/68")
cat("PASS: export layout, top-row block labels, identifiers, row count, types, ranges, dates and biomarker counts\n")
