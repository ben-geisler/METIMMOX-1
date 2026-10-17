# Horizon conservation and separation of occupancy rates from scheduled events (#163);
# life-years integrated like the QALYs (#197); progressed-state follow-up and CT
# every 2.5 months (4.8 per year) and two contacts per FLOX administration (#207).
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
  # Undiscounted life-years with no mortality equal the horizon (#197).
  equal(run(params, v_pf)$ly_undisc_total, years)
  p <- params
  p$c_other_follow <- 33
  p$c_other_pp <- 5000
  p$c_test_CT <- 386
  out <- run(p, rep(0, n))
  # Follow-up visit and CT every 2.5 months in the progressed state (issues
  # #164, #207); the post-progression treatment cost stays quarterly.
  equal(sum(out$follow_up_costs), 33 * 4.8 * years)
  equal(sum(out$post_progression_costs), 5000 * 4 * years)
  equal(sum(out$progressed_imaging_costs), 386 * 4.8 * years)
}
equal(run(params, 1)$qalys_total, 0)
equal(run(params, 1)$ly_undisc_total, 0)
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
equal(out$costs_total, (4.8 * 33 + 4 * 5000) * area(v_pp * v_dw))
# Life-years (#197): the alive occupancy with the QALY interval weights and the
# effect discount factors. With both utilities at 1 (params), discounted QALYs
# equal discounted life-years; with other utilities the life-years are unchanged.
out_u1 <- run(params, v_pf, v_pp, 1 - v_pf - v_pp, v_dw)
equal(out_u1$ly_disc_total, out_u1$qalys_total)
equal(out_u1$ly_disc_total, area((v_pf + v_pp) * v_dw))
equal(out_u1$ly_undisc_total, area(v_pf + v_pp))
equal(out$ly_disc_total, out_u1$ly_disc_total)
stopifnot(out_u1$ly_disc_total < out_u1$ly_undisc_total)
# Undiscounted life-years are the restricted mean survival of the OS curve.
v_os_ly <- c(1, 0.9, 0.6, 0.2)
equal(run(params, c(1, 0.5, 0.3, 0.1), v_os_ly - c(1, 0.5, 0.3, 0.1), 1 - v_os_ly)$ly_undisc_total,
      restricted_mean_survival(v_os_ly, 1 / 52))
# Life-years come from OS (1 - p_d), not from p_pf + p_p: an uncapped PFS that
# crosses OS (p_p = pmax(os - pfs, 0) = 0 there) would otherwise give
# max(pfs, os) and overstate them. With the crossing the QALYs at utility 1
# exceed the life-years, so the QALY = life-year identity holds only without one.
v_pfs_cross <- c(1, 0.95, 0.7, 0.4)
out_cross <- run(params, v_pfs_cross, pmax(v_os_ly - v_pfs_cross, 0), 1 - v_os_ly, v_dw)
equal(out_cross$ly_undisc_total, area(v_os_ly))
equal(out_cross$ly_disc_total, area(v_os_ly * v_dw))
stopifnot(out_cross$ly_disc_total < out_cross$qalys_total,
          area(pmax(v_pfs_cross, v_os_ly)) - out_cross$ly_undisc_total > 1e-4)

# With a CT price and no l_CT scan, progressed imaging is the only CT charge.
p$c_test_CT <- 386
out <- run(p, v_pf, v_pp, 1 - v_pf - v_pp, v_dw)
equal(out$costs_total, (4.8 * (33 + 386) + 4 * 5000) * area(v_pp * v_dw))

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
# Per scheduled point: nivolumab 100, FLOX 20, CT 30, blood 10 and two visits of
# 40 (the FLOX administration is two contacts, issue #207) = 240.
equal(run(p, rep(1, 3), dw = c(1, 0.9, 0.8))$costs_total, 240 * 1.8 + 110)
equal(run(p, 1)$costs_total, 350)
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
