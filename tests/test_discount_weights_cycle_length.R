# ===============================================================================
# TEST: discount_weights() honours the model's cycle length (issue #82)
# ===============================================================================
# discount_weights() previously hard-coded a 52-week year:
#   years <- (seq_len(n_cycles) - 1) / 52
# That silently assumes weekly cycles. If cl ever changed (monthly, or a
# different horizon granularity), costs and QALYs would be discounted over the
# wrong number of years while every other part of the model used cl, and nothing
# would error -- the results would just be quietly wrong.
#
# The weekly result must be unchanged within the 1e-10 tolerance the test
# protocol (docs/validation/README.md of the private repository)
# documents for numeric comparisons, so this fix does not invalidate the
# regenerated caches or snapshots.
#
# Run: "C:\Program Files\R\R-4.3.2\bin\x64\Rscript.exe" \
#        tests/test_discount_weights_cycle_length.R
# ===============================================================================

source(here::here("R", "calculate_outcomes.R"))

TOL <- 1e-10
params <- list(dr_costs = 0.04, dr_effects = 0.04, cl = 1 / 52)
n_cycles <- 521L   # time_horizon 520 + 1

# ---------------------------------------------------------------------------
# 1. Weekly weights match the legacy /52 form within the documented tolerance
# ---------------------------------------------------------------------------
v_legacy_years <- (seq_len(n_cycles) - 1) / 52
v_legacy_cost <- 1 / (1 + params$dr_costs)^v_legacy_years
v_legacy_effect <- 1 / (1 + params$dr_effects)^v_legacy_years

w <- discount_weights(params, n_cycles)
stopifnot(
  length(w$cost) == n_cycles,
  length(w$effect) == n_cycles,
  max(abs(w$cost - v_legacy_cost)) < TOL,
  max(abs(w$effect - v_legacy_effect)) < TOL
)

# ---------------------------------------------------------------------------
# 2. No discounting at t = 0, and weights decrease monotonically thereafter
# ---------------------------------------------------------------------------
stopifnot(
  w$cost[1] == 1,
  w$effect[1] == 1,
  all(diff(w$cost) < 0),
  all(w$cost > 0)
)

# ---------------------------------------------------------------------------
# 3. cl is actually honoured -- the point of the fix
#    A monthly cycle spans ~4.33x more calendar time per cycle than a weekly
#    one, so by the final cycle the weights must differ substantially. Under the
#    old hard-coded /52 this assertion fails: cl was ignored entirely.
# ---------------------------------------------------------------------------
w_monthly <- discount_weights(params, n_cycles, cl = 1 / 12)
stopifnot(
  abs(w_monthly$cost[n_cycles] - w$cost[n_cycles]) > 0.1,
  w_monthly$cost[n_cycles] < w$cost[n_cycles]   # longer horizon, more discounting
)

# An explicit cl overrides params$cl rather than being ignored.
w_override <- discount_weights(params, n_cycles, cl = 1 / 12)
stopifnot(isTRUE(all.equal(w_override$cost, w_monthly$cost)))

# ---------------------------------------------------------------------------
# 4. A missing or invalid cl errors rather than silently guessing
# ---------------------------------------------------------------------------
for (bad in list(NULL, 0, -1 / 52, NA_real_, c(1 / 52, 1 / 12), "1/52")) {
  bad_params <- list(dr_costs = 0.04, dr_effects = 0.04, cl = bad)
  err <- tryCatch({
    discount_weights(bad_params, 10L)
    NULL
  }, error = function(e) conditionMessage(e))
  stopifnot(!is.null(err), grepl("cycle length", err, fixed = TRUE))
}

# ---------------------------------------------------------------------------
# 5. Zero discount rates leave every weight at 1, whatever the cycle length
# ---------------------------------------------------------------------------
w0 <- discount_weights(list(dr_costs = 0, dr_effects = 0, cl = 1 / 52), n_cycles)
stopifnot(all(w0$cost == 1), all(w0$effect == 1))

cat("PASS: discount_weights() honours cycle length; weekly result unchanged.\n")
