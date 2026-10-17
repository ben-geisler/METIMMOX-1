# Structural test of the raw export contract on a synthetic export (issue #174),
# including the treatment-administration block and the "Dates" sheet (issue
# #208). Needs no trial data.
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

# Treatment-administration block (issue #208): sequence 1 for everyone,
# sequence 2 for the first 16 control and 20 experimental patients, sequence 3
# for 5 and 13 of them, sequence 4 for 4 experimental patients; FLOX every two
# weeks in the control arm, 2 FLOX + 2 nivolumab alternating otherwise.
# Sequences 1-2 are date columns, 3-4 text columns of Excel serial numbers
# with a "Pause" token, as readxl returns them.
v_arm <- df_syn$`Study1/control0`
v_rank <- ave(seq_len(n), v_arm, FUN = seq_along)
v_seq_n <- list(c(36, 38), c(16, 20), c(5, 13), c(0, 4))
v_seq_start <- c(0, 24, 52, 80)
#' Patients of the synthetic export who have sequence k.
seq_has <- function(k) v_rank <= v_seq_n[[k]][v_arm + 1]
#' Excel serial number (days since 1899-12-30) as text.
serial <- function(d) as.character(as.numeric(d - as.Date("1899-12-30")))
for (k in 1:4) for (v in 1:8) {
  j <- sequence_block_column(k, v)
  v_given <- seq_has(k)
  v_date <- v_inclusion + 7 * (v_seq_start[k] + 2 * (v - 1))
  if (k <= 2) {
    df_syn[[j]] <- posix(replace(v_date, !v_given, NA))
  } else {
    df_syn[[j]] <- ifelse(v_given, serial(v_date), NA)
    if (v == 1) df_syn[[j]][which(!v_given)[1]] <- "Pause"
  }
  v_nivo_visit <- v %in% c(3, 4, 7, 8) & v_arm == 1
  df_syn[[j + 1]] <- ifelse(v_given, ifelse(v_nivo_visit, "Nivolumab", "Oxaliplatin"), NA)
  df_syn[[j + 3]] <- ifelse(v_given & !v_nivo_visit, "5-FU", NA)
}
# A spelling variant counts as an administration of the block but not in the
# export's totals; one further patient's total misses one nivolumab dose.
df_syn[[sequence_block_column(1, 1) + 1]][1] <- "Oxalliplatin"
df_syn[[sequence_block_column(1, 3) + 1]][37] <- "Nivoliumab"
v_flox <- rowSums(sapply(1:4, function(k) seq_has(k) * ifelse(v_arm == 1, 4, 8)))
v_nivo <- rowSums(sapply(1:4, function(k) seq_has(k) * ifelse(v_arm == 1, 4, 0)))
df_syn$`Cycles FLOX` <- v_flox - (seq_len(n) == 1)
df_syn$`Cycles Nivolumab` <- v_nivo - (seq_len(n) %in% c(37, 38))

stopifnot(identical(check_export_contract(df_syn, v_top), df_syn))
df_seq <- derive_treatment_sequences(df_syn)
stopifnot(identical(names(df_seq), c(paste0("Seq", rep(2:4, each = 4),
                                            c("", "wk", "_FLOX", "_nivo")), "AdmFLOX", "AdmNivo")),
          all(sapply(2:4, function(k) all(tapply(df_seq[[paste0("Seq", k)]], v_arm, sum) == v_seq_n[[k]]))),
          all(df_seq$Seq3wk[df_seq$Seq3 == 1] == 52), all(is.na(df_seq$Seq3wk[df_seq$Seq3 == 0])),
          all(df_seq$Seq4_nivo[df_seq$Seq4 == 1] == 4),
          all(df_seq$Seq3_FLOX[df_seq$Seq3 == 1 & v_arm == 0] == 8),
          all(df_seq$AdmFLOX == v_flox), all(df_seq$AdmNivo == v_nivo))
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

# Treatment-administration block (issue #208).
# A whole visit block inserted moves the block's first-row labels.
fails(check_export_contract(df_syn, replace(v_top, sequence_block_column(3, 1) - 34L, "S2V8")),
      "top header at column 2686 is 'S2V8', expected 'S3V1'")
fails(check_export_contract(modify(function(x) {
  x[[sequence_block_column(2, 1) + 1]][1] <- "Irinotecan"; x }), v_top),
  "Treatment...2681 has 1 unknown treatment entry")
fails(check_export_contract(modify(function(x) {
  x[[sequence_block_column(3, 1)]][2] <- "later"; x }), v_top),
  "Date...2720 has 1 entry that are neither dates nor")
fails(check_export_contract(modify(function(x) {
  x[[sequence_block_column(2, 2)]][1] <- NA; x }), v_top),
  "1 visit(s) with a drug but no date")
fails(check_export_contract(modify(function(x) {
  for (v in 1:8) { j <- sequence_block_column(2, v); x[[j]][1] <- NA; x[[j + 1]][1] <- NA }; x }), v_top),
  "sequence 3 started without sequence 2")
fails(check_export_contract(modify(function(x) {
  x[[sequence_block_column(3, 1)]][1] <- serial(v_inclusion[1] + 7 * 30); x }), v_top),
  "sequence 3 starts before sequence 2 ends")
fails(check_export_contract(modify(function(x) { x$`Cycles FLOX`[3] <- NA; x }), v_top),
      "Cycles FLOX must be complete nonnegative whole numbers")

# "Dates" sheet: formulas reading the block. Entries 74/37/22/17; a start cell
# referring to an empty block cell reads as Excel day zero.
df_dates <- as.data.frame(setNames(replicate(65, rep(NA_character_, n), simplify = FALSE),
                                   sprintf("Col...%d", 1:65)), check.names = FALSE)
names(df_dates)[as.integer(names(contract$dates_sheet$headers))] <- contract$dates_sheet$headers
#' Block date of the first visit of sequence k.
block_start <- function(k) parse_export_date(df_syn[[sequence_block_column(k, 1)]], "x")
#' Dates-sheet start cells: n_real block dates and n_zero day-zero entries.
start_cells <- function(k, n_real, n_zero) {
  v <- rep(NA_character_, n)
  v_real <- which(seq_has(k))[seq_len(n_real)]
  v[v_real] <- serial(block_start(k)[v_real])
  v[setdiff(seq_len(n), which(seq_has(k)))[seq_len(n_zero)]] <- "0"
  v
}
df_dates$ID <- df_syn$ID
df_dates$`Study1/control0` <- v_arm
df_dates$`Date of inclusion` <- posix(v_inclusion)
df_dates$`Start first cycle` <- posix(block_start(1))
df_dates$`Start second cycle` <- start_cells(2, 22, 15)
df_dates$`Start third cycle` <- start_cells(3, 7, 15)
df_dates$`Start fourth cycle` <- start_cells(4, 0, 17)
df_dates$`Start fourth cycle`[which(!is.na(df_dates$`Start fourth cycle`))[1]] <- "Pause"
df_dates <- df_dates[rev(seq_len(n)), ]  # row order need not match the Datasheet
df_dates_summary <- check_dates_sheet(df_dates, df_syn)
stopifnot(identical(df_dates_summary$entries, c(74, 37, 22, 17)),
          identical(df_dates_summary$real_dates, c(74, 22, 7, 0)),
          identical(df_dates_summary$with_administration, c(74, 22, 7, 0)),
          identical(df_dates_summary$administered, c(74, 36, 18, 4)))
fails(check_dates_sheet(within(df_dates, `Start third cycle` <- `Start fourth cycle`),
                        df_syn),
      "Dates sheet start entries 74/37/17/17, expected 74/37/22/17")
i_real <- which(!is.na(df_dates$`Start second cycle`) & df_dates$`Start second cycle` != "0")[1]
fails(check_dates_sheet(replace(df_dates, "Start second cycle",
                                list(replace(df_dates$`Start second cycle`, i_real,
                                             serial(as.Date("2019-06-01"))))), df_syn),
      "Start second cycle: 1 date(s) differ from the first visit of sequence 2")
fails(check_dates_sheet(setNames(df_dates, replace(names(df_dates), 62, "Start 3rd cycle")), df_syn),
      "Dates sheet column 62 is 'Start 3rd cycle', expected 'Start third cycle'")
fails(check_dates_sheet(df_dates[-1, ], df_syn), "does not hold exactly the analysed patients")

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
df_derived <- cbind(df_derived, df_seq[c("Seq2", "Seq3", "Seq4", "AdmFLOX", "AdmNivo")],
                    df_syn[c("Cycles FLOX", "Cycles Nivolumab")])
stopifnot(identical(check_export_counts(df_derived), df_derived))
fails(check_export_counts(transform(df_derived, Seq4 = replace(Seq4, 1, 1))),
      "sequence_4 positives 1/36 vs 4/38, expected 0/36 vs 4/38")
fails(check_export_counts(replace(df_derived, "Cycles FLOX", list(df_derived$AdmFLOX))),
      "agrees with 'Cycles FLOX'/'Cycles Nivolumab' for 74/72 patients, expected 73/72")
fails(check_export_counts(transform(df_derived, CRP1cat = c(cat_vec(7, 35, 36), cat_vec(16, 36, 38)))),
      "crp_week4 positives 7/35 vs 16/36, expected 7/35 vs 17/36")
fails(check_export_counts(transform(df_derived, Mutation = replace(Mutation, 74, "BRAF"))),
      "tmb_braf positives 32/70, expected 31/69")
fails(check_export_counts(transform(df_derived, TLRcat = cat_vec(44, 67, 74))),
      "tlr positives 44/67, expected 44/68")
cat("PASS: export layout, top-row block labels, identifiers, row count, types, ranges, dates, biomarker counts, administration block, treatment sequences and Dates sheet\n")
