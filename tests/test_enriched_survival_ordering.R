# ===============================================================================
# TEST: partitioned_survival_states() enforces OS >= PFS loudly, not silently
# ===============================================================================
# Regression guard for the enriched-population path (08b). That path previously
# passed enforce_order = TRUE, which applied pmin(pfs, os) with no warning and no
# error, so an ordering violation in the enriched curves could not be observed.
# model_fun() stop()s on the same condition in deterministic mode; the enriched
# analysis is deterministic and must behave identically.
#
# Run: "C:\Program Files\R\R-4.3.2\bin\x64\Rscript.exe" \
#        tests/test_enriched_survival_ordering.R
# ===============================================================================

source(here::here("R", "calculate_outcomes.R"))

v_os <- c(1.00, 0.90, 0.80, 0.70)
v_pfs_ok <- c(1.00, 0.85, 0.70, 0.60)
v_pfs_bad <- c(1.00, 0.95, 0.70, 0.60)  # violates at index 2

# ---------------------------------------------------------------------------
# 1. Unlabelled calls stay permissive (model_fun has already enforced ordering)
# ---------------------------------------------------------------------------
states <- partitioned_survival_states(v_os, v_pfs_ok)
stopifnot(
  identical(states$p_pf[1], 1),
  identical(states$p_p[1], 0),
  identical(states$p_d[1], 0),
  all(states$p_p >= 0),
  isTRUE(all.equal(states$p_pf + states$p_p + states$p_d, rep(1, length(v_os))))
)

# An unlabelled call must NOT error even on a violating pair, because model_fun
# clamps upstream and relies on this function staying a pure constructor.
invisible(partitioned_survival_states(v_os, v_pfs_bad))

# ---------------------------------------------------------------------------
# 2. A labelled call accepts an ordered pair
# ---------------------------------------------------------------------------
labelled_ok <- partitioned_survival_states(v_os, v_pfs_ok, "enriched crp+ control")
stopifnot(isTRUE(all.equal(labelled_ok$p_pf, states$p_pf)))

# ---------------------------------------------------------------------------
# 3. A labelled call REFUSES a violating pair, naming the curve and the excess
# ---------------------------------------------------------------------------
violation <- tryCatch({
  partitioned_survival_states(v_os, v_pfs_bad, "enriched crp+ control")
  NULL
}, error = function(e) conditionMessage(e))

stopifnot(
  !is.null(violation),
  grepl("enriched crp+ control", violation, fixed = TRUE),
  grepl("1 time points", violation, fixed = TRUE),
  grepl("max excess=0.05", violation, fixed = TRUE)
)

# ---------------------------------------------------------------------------
# 4. The enriched analysis script uses labelled calls (no silent clamp left)
# ---------------------------------------------------------------------------
v_enriched_src <- readLines(
  here::here("analysis", "08b_enriched_population_analysis.R"),
  warn = FALSE
)
v_calls <- grep("partitioned_survival_states", v_enriched_src, value = TRUE)
stopifnot(
  length(v_calls) >= 2,
  !any(grepl("TRUE", v_calls, fixed = TRUE))
)

cat("PASS: enriched-population ordering invariant is enforced loudly.\n")
