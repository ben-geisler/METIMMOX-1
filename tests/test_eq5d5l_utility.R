# Regression tests for the EQ-5D-5L value-set functions extracted from the
# archived 05_QALYs.Rmd notebook. Run from the repository root with Rscript.

source("R/eq5d5l_utility.R")

stopifnot(
  identical(eq5d5l_utility_dk(1, 1, 1, 1, 1), 1),
  isTRUE(all.equal(eq5d5l_utility_dk(5, 5, 5, 5, 5), -0.758)),
  identical(eq5d5l_utility_uk(1, 1, 1, 1, 1), 1),
  isTRUE(all.equal(eq5d5l_utility_uk(5, 5, 5, 5, 5), -0.285)),
  identical(
    eq5d5l_utility_dk(c(1, 2), 1, 1, 1, 1),
    c(1, 0.959)
  ),
  is.na(eq5d5l_utility_uk(NA_real_, 1, 1, 1, 1))
)

cat("EQ-5D-5L utility tests passed.\n")
