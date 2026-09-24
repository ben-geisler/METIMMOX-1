source("R/pfs_endpoint.R")
derive <- function(x, ...) derive_pfs_endpoint(x, ..., verbose = FALSE)
fails <- function(expr) {
  e <- tryCatch({force(expr); NULL}, error = identity)
  stopifnot(inherits(e, "error"), grepl("derive_pfs_endpoint()", conditionMessage(e), fixed = TRUE))
}
d <- data.frame(Progression = c(1, 0, 0, 0, 0),
                PFSwk = c(5, 10, 10, 20, 10), Death = c(1, 1, 1, 0, 1),
                OSwk = c(30, 15, 40, 30, 26), LastEvalwk = c(6, 12, 25, 22, 12))
x <- derive(d)
stopifnot(identical(x$Progression, c(1, 1, 0, 0, 1)),
          identical(x$PFSwk, c(5, 15, 10, 20, 26)),
          identical(as.character(x$PFS_rule), c("progression", "death_within_window",
            "death_censored_at_assessment", "censored_alive", "death_within_window")),
          identical(x, derive(x)), identical(x$ProgressionExit, d$Progression),
          identical(x$TTPwk, d$PFSwk), sum(attr(x, "pfs_endpoint")$counts) == nrow(d))
old <- derive(x, Inf)
stopifnot(identical(old$Progression, as.numeric(d$Progression == 1 | d$Death == 1)),
          identical(old$PFSwk, ifelse(d$Progression == 1, d$PFSwk,
                                    ifelse(d$Death == 1, d$OSwk, d$PFSwk))),
          identical(derive(old), x), identical(derive(x, 0), derive(d, 0)))
alt <- derive(x, 8, "LastEvalwk")
stopifnot(alt$PFSwk[3] == 25, alt$Progression[3] == 0, alt$PFSwk[4] == 22,
          identical(derive(alt), x))
fails(derive(d[, setdiff(names(d), "LastEvalwk")], 16, "LastEvalwk"))
fails(derive(transform(d, PFSwk = OSwk + 1)))
fails(derive(transform(d, LastEvalwk = OSwk + 1), 16, "LastEvalwk"))
fails(derive(transform(d, Death = NA_real_)))
fails(derive(transform(d, Progression = NA_real_)))
fails(derive(transform(d, Death = 2)))
for (window in list(-1, NA_real_, NaN, -Inf, numeric(), c(8, 16), "16")) fails(derive(d, window))
fails(derive(d, 16, "bad_anchor"))
missing <- d
missing$PFSwk[2] <- NA_real_
for (window in c(16, Inf)) {
  out <- derive(missing, window)
  stopifnot(is.na(out$Progression[2]), is.na(out$PFSwk[2]),
            out$PFS_rule[2] == "anchor_missing")
}
stopifnot(nrow(derive(d[FALSE, ])) == 0)
cat("PASS: PFS window, anchors, boundary, raw preservation, re-derivation and validation\n")
