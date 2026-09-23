# Independent arithmetic examples for primary probabilistic reporting (#180).
source("R/cea_helpers.R")
source("R/report_format.R")
make_psa <- function(cost, effect, strategies) {
  dampack::make_psa_obj(cost = as.data.frame(cost),
    effectiveness = as.data.frame(effect), strategies = strategies)
}
# Put control second to exercise lookup; strongly correlated levels make
# subtracting marginal interval endpoints visibly wrong for paired increments.
cost <- cbind(c(110, 130, 150, 170), c(100, 110, 120, 130), c(120, 130, 140, 150))
effect <- cbind(c(2, 4, 6, 8), c(1, 2, 3, 4), c(0.5, 1.5, 2.5, 3.5))
p <- make_psa(cost, effect, c("crp", "control", "tmb_braf"))
x <- create_psa_summary_table(p, wtp = 100)
stopifnot(x$Mean_Cost[1] == 140, x$Mean_QALY[1] == 5,
  x$Mean_Inc_Cost[1] == 25, x$Mean_Inc_QALY[1] == 2.5,
  x$ICER[1] == 10, abs(x$ICER[1] - mean(c(10, 20, 30, 40) / c(1, 2, 3, 4))) < 1e-10,
  x$Cost_Lower[1] == 111.5, x$Cost_Upper[1] == 168.5,
  x$QALY_Lower[1] == 2.15, x$QALY_Upper[1] == 7.85,
  x$Inc_Cost_Lower[1] == 10.75, x$Inc_Cost_Upper[1] == 39.25,
  x$Inc_QALY_Lower[1] == 1.075, x$Inc_QALY_Upper[1] == 3.925,
  x$Mean_Inc_Cost[2] == 0, x$Inc_Cost_Upper[2] == 0,
  is.na(x$ICER[2]), x$Status[2] == "Reference",
  x$Status[3] == "Dominated by SoC", format_psa_summary(x)$ICER[3] == "Dominated")
# Nonconstant per-draw ratios: ratio of means 10, mean of ratios 25/3.
p$cost[, 1] <- p$cost[, 2] + c(0, 20, 40, 40)
x <- create_psa_summary_table(p, 100)
stopifnot(x$ICER[1] == 10,
          abs(x$ICER[1] - mean(c(0, 20, 40, 40) / c(1, 2, 3, 4))) > 1)
# Zero incremental effects never produce an infinite ICER; dominance still shows.
p$effect[, 1] <- p$effect[, 2]
x <- create_psa_summary_table(p, 100)
stopifnot(is.na(x$ICER[1]), x$Status[1] == "Dominated by SoC")
p$cost[, 1] <- p$cost[, 2] - 10
x <- create_psa_summary_table(p, 100)
stopifnot(x$Status[1] == "Cost-saving vs SoC")
stopifnot(identical(names(x), names(create_psa_summary_table(NULL))),
          nrow(create_psa_summary_table(NULL)) == 0)
p$cost[1, 1] <- NA_real_
stopifnot(inherits(try(create_psa_summary_table(p), silent = TRUE), "try-error"))
cat("PASS: paired percentile intervals, ratio-of-means ICER, reference and dominance semantics.\n")
