# Horizon conservation and separation of occupancy rates from scheduled events (#163).
source("R/calculate_outcomes.R")
source("R/model_configs.R")
source("R/report_format.R")
source("R/pfs_os_violation_diagnostics.R")

params <- list(u_np = 1, u_p = 1, c_drug_nivo = 0, c_drug_FLOX = 0,
  c_test_CT = 0, c_test_blood = 0, c_test_CRP = 0, c_test_NGS = 0,
  c_other_visit = 0, c_other_baseline = 0, c_other_follow = 0,
  c_other_pp = 0, c_other_last = 0,
  l_nivo = 0, l_FLOX_exp = 0, l_FLOX_control = 0,
  l_CT = 0, l_blood = 0, l_visit = 0)
#' calculate_outcomes() for an experimental-arm, CRP-guided trace with the given
#' state occupancies, discount weights and cycle length; returns its outcome list.
run <- function(p, pf, progressed = 1 - pf, dead = rep(0, length(pf)),
                dw = rep(1, length(pf)), cl = 1 / 52) {
  calculate_outcomes(p, pf, progressed, dead, "experimental", "crp", dw, dw, cl)
}
#' Assert that two numeric vectors agree within 1e-10.
equal <- function(x, y) stopifnot(max(abs(x - y)) < 1e-10)
for (years in c(1, 10, 20)) {
  n <- years * 52 + 1
  v_pf <- rep(1, n)
  equal(run(params, v_pf)$qalys_total, years)
  equal(run(params, v_pf)$qalys_total, restricted_mean_survival(v_pf, 1 / 52))
  p <- params
  p$c_other_follow <- 33
  p$c_other_pp <- 5000
  p$c_test_CT <- 386
  out <- run(p, rep(0, n))
  equal(sum(out$follow_up_costs), 33 * 4 * years)
  equal(sum(out$post_progression_costs), 5000 * 4 * years)
  # Quarterly CT in the progressed state (issue #164).
  equal(sum(out$progressed_imaging_costs), 386 * 4 * years)
}
equal(run(params, 1)$qalys_total, 0)
equal(run(params, rep(1, 5), cl = 1 / 4)$qalys_total, 1)

# Discounted, nonconstant occupancy matches an independent interval sum.
v_pf <- c(1, 0.7, 0.4, 0)
v_pp <- c(0, 0.2, 0.4, 0.5)
v_dw <- 1 / 1.04^(0:3 / 52)
p <- params
p$u_np <- 0.73
p$u_p <- 0.59
p$c_other_follow <- 33
p$c_other_pp <- 5000
#' Trapezoidal integral of a four-point weekly curve over weeks 0-3, in years.
area <- function(y) sum(diff(0:3 / 52) * (head(y, -1) + tail(y, -1)) / 2)
out <- run(p, v_pf, v_pp, 1 - v_pf - v_pp, v_dw)
equal(out$qalys_total, area((v_pf * p$u_np + v_pp * p$u_p) * v_dw))
equal(out$costs_total, 4 * (33 + 5000) * area(v_pp * v_dw))
# With a CT price and no l_CT scan, progressed imaging is the only CT charge.
p$c_test_CT <- 386
out <- run(p, v_pf, v_pp, 1 - v_pf - v_pp, v_dw)
equal(out$costs_total, 4 * (33 + 386 + 5000) * area(v_pp * v_dw))

# Full scheduled costs at BOTH horizon endpoints, even for a one-point grid.
p <- params
p$c_drug_nivo <- 100
p$c_drug_FLOX <- 20
p$c_test_CT <- 30
p$c_test_blood <- 10
p$c_other_visit <- 40
p$c_other_baseline <- 50
p$c_test_CRP <- 60
for (nm in c("l_nivo", "l_FLOX_exp", "l_CT", "l_blood", "l_visit")) p[[nm]] <- c(1, 0, 1)
equal(run(p, rep(1, 3), dw = c(1, 0.9, 0.8))$costs_total, 200 * 1.8 + 110)
equal(run(p, 1)$costs_total, 310)
p <- params
p$c_other_last <- 1000
equal(run(p, c(1, 0), c(0, 0), c(0, 1), c(1, 0.8))$costs_total, 800)

# Diagnostic QALYs use exactly the same endpoint weights as the economic model.
v_os <- c(1, 0.6, 0.3, 0.1)
v_pfs <- c(1, 0.7, 0.2, 0.2)
p <- params
p$u_np <- 0.73
p$u_p <- 0.59
metrics <- pfs_os_violation_metrics(v_os, v_pfs, 0:3, 1 / 52, v_dw, p$u_np, p$u_p)
equal(metrics$qaly_clamped, run(p, pmin(v_pfs, v_os), pmax(v_os - v_pfs, 0), 1 - v_os, v_dw)$qalys_total)
equal(metrics$qaly_raw_minus_clamped, area(pmax(v_pfs - v_os, 0) * (p$u_np - p$u_p) * v_dw))
cat("Interval integration and scheduled-event tests passed.\n")
