# PFS \> OS Ordering Violations in the PSA
Ben Geisler
2026-09-30

- [Overview](#overview)
- [Methods](#methods)
  - [Two sources of variation in a PSA
    draw](#two-sources-of-variation-in-a-psa-draw)
  - [Curves examined](#curves-examined)
  - [Metrics](#metrics)
  - [Mapping draws to PSA rows](#mapping-draws-to-psa-rows)
  - [Cache provenance](#cache-provenance)
- [Results](#results)
  - [How many draws are affected](#how-many-draws-are-affected)
  - [Where and how far the curves
    cross](#where-and-how-far-the-curves-cross)
  - [PFS-over-OS excess area before and after the
    clamp](#pfs-over-os-excess-area-before-and-after-the-clamp)
  - [QALY consequence of the clamp](#qaly-consequence-of-the-clamp)
  - [Is the clamp why the PSA increments differ from the base
    case?](#is-the-clamp-why-the-psa-increments-differ-from-the-base-case)
- [Interpretation](#interpretation)
- [Notes](#notes)

# Overview

The probabilistic sensitivity analysis (PSA) draws the coefficient
vectors of the overall survival (OS) and progression-free survival (PFS)
gamma models jointly (issue \#159). Their cross-endpoint dependence is
estimated by fitting both endpoints to the same patient bootstrap
resamples, while preserving the fitted marginal covariance matrices.
This replaces the independent draws of issue \#156. Correlation does not
enforce survival ordering: a sampled pair can still predict PFS above OS
on part of the 10-year weekly grid. A partitioned survival model cannot
represent a negative progressed-state occupancy OS - PFS.

The model handles the problem with a clamp. In PSA mode `model_fun()`
replaces the PFS curve by min(PFS, OS) for any curve that crosses and
emits a `survival_ordering_warning`; the PSA script counts those
warnings per subgroup but does not retain the curves. In the
deterministic base case the ordering is guaranteed by the
ordering-constrained distribution selection and any crossing stops the
model.

This report quantifies the crossings behind those warnings. For every
one of the 5,000 retained PSA rows it re-predicts the five
population-averaged subgroup curves that the PSA uses and reports:

1.  the number and share of draws, and of retained PSA rows, with at
    least one PFS \> OS violation, overall and per curve;
2.  where in the horizon the curves cross and by how much;
3.  the **PFS-over-OS excess area**, the area between the two curves on
    the stretch where PFS lies above OS (the trough in the figures),
    before the clamp (after the clamp it is zero by construction),
    together with the progression-free and progressed-state areas before
    and after the clamp;
4.  the discounted-QALY consequence of the clamp per curve and per
    strategy.

# Methods

## Two sources of variation in a PSA draw

Each PSA iteration draws one combined OS/PFS coefficient vector (log
shapes, log baseline rates, and covariate effects on log rate) from a
joint multivariate normal. Each endpoint retains the mean and marginal
covariance of its own fit. Within the iteration the OS vector is shared
by every patient and every curve, and so is the PFS vector. Patient
heterogeneity is integrated out by averaging the patient curves with the
shared weights from that row’s joint biomarker distribution (issue
\#166). The paired-bootstrap cross-covariance links the two endpoints
but does not impose OS \>= PFS. The [joint sampling
specification](joint_survival_sampling.md) reports bootstrap convergence
and coefficient dependence.

## Curves examined

The PSA evaluates five curves per draw, each a population average over
the complete-case patients of the subgroup with `Rx` set as the strategy
dictates (issue \#151). They are the curves that `model_fun()` checks
and clamps one by one. The short labels in the first column are used in
the tables that follow.

| Label     | Curve                            | Patients | Rx set to    | Strategy |
|:----------|:---------------------------------|---------:|:-------------|:---------|
| SoC       | Standard of care (control)       |       68 | Control      | Control  |
| CRP+      | CRP positive (experimental)      |       23 | Experimental | CRP      |
| CRP-      | CRP negative (control)           |       45 | Control      | CRP      |
| TMB/BRAF+ | TMB/BRAF positive (experimental) |       30 | Experimental | TMB/BRAF |
| TMB/BRAF- | TMB/BRAF negative (control)      |       38 | Control      | TMB/BRAF |

Subgroup curves evaluated per draw

The curves are computed in closed form from the gamma
accelerated-failure-time parameters of each draw (shape and
covariate-specific rate) rather than through `flexsurv::predict()`,
which makes the full pass about twenty times faster. The closed form was
checked against the pipeline’s own prediction helper
(`generate_psa_population_averaged_predictions()`) on 3 draws and all
five curves for both outcomes: the largest absolute difference at any of
the 521 weekly points is 3.3e-16. The clamped QALYs computed here
reproduce the strategy QALYs that `model_fun()` returns for the same
draw.

## Metrics

For each draw and curve, with $t$ running over the weekly grid (0 to 520
weeks) and $\Delta_t = \max(\mathrm{PFS}_t - \mathrm{OS}_t, 0)$ the
excess of PFS over OS:

- **Violation**: $\Delta_t > 0$ at one or more time points. The number
  of violating time points, the first and last violating week, and the
  largest excess are recorded.
- **PFS-over-OS excess area (before the clamp)**: $\int \Delta_t \, dt$
  by the trapezoidal rule on the weekly grid, in patient-years. It is
  the area between the PFS and OS curves on the stretch where PFS lies
  on top. After the clamp, PFS is replaced by min(PFS, OS) and the area
  is zero.
- **State areas**: the progression-free area $\int \mathrm{PFS}_t\,dt$
  and the progressed-state area
  $\int (\mathrm{OS}_t - \mathrm{PFS}_t)\,dt$ before the clamp (the
  latter negative on the violating stretch), and
  $\int \min(\mathrm{PFS}_t, \mathrm{OS}_t)\,dt$ and
  $\int \max(\mathrm{OS}_t - \mathrm{PFS}_t, 0)\,dt$ after it, in
  patient-years. The clamp moves exactly the excess area from the
  progression-free state to the progressed state; OS is untouched.
- **QALY consequence**: discounted QALYs the raw curves would have
  produced minus the clamped QALYs. With $p_{PF} = \mathrm{PFS}$ and
  $p_P = \mathrm{OS} - \mathrm{PFS}$, the difference is
  $\sum_t h_t \Delta_t \,(u_{np} - u_p)\, cl\, w_t$, where $h_t$ is 0.5
  at the two horizon endpoints and 1 elsewhere. This uses each PSA row’s
  sampled utilities, $cl = 1/52$ and effect discount weights $w_t$ at 4%
  per year. It is the amount by which the raw curves overstate QALYs by
  counting the excess as progression-free time at the higher utility.

Areas and discounted QALYs both use the trapezoidal rule, matching
`restricted_mean_survival()` and `model_fun()` (issue \#163).

## Mapping draws to PSA rows

The PSA object records, for every retained row (`sim`), the
sampling-cache draw that produced it (`model_idx`; issue \#156).
Diagnostics retain the row’s sampled joint-population weights and
utilities and use its `model_idx` for the survival coefficients. Rows
whose original draw failed therefore keep their economic parameters
while using the replacement model. In the current cache 0 of 5,000 rows
were replaced. Counts below refer to these retained PSA rows.

## Cache provenance

The diagnostics are cached with a fingerprint of the sampling-cache
fingerprint, the explicit prediction population, the PSA parameter table
(including joint weights, utilities and model indices), the base
parameters, the discount rate, the horizon and the cycle length. Changes
to these inputs force regeneration. The PSA cache stores the same
sampling fingerprint. Before rendering, this report stops unless all
three caches carry the same sampling fingerprint and the diagnostics
file is younger than the sampling-cache file.

| Cache | Modified | Sampling fingerprint |
|:---|:---|:---|
| Sampling draws | 2026-09-30 15:03 | 09df4ea65ee297a56ef6fdef1e6207c2d052da3f908229ed359a84fd16b9ced8 |
| PSA results | 2026-09-30 16:51 | 09df4ea65ee297a56ef6fdef1e6207c2d052da3f908229ed359a84fd16b9ced8 |
| Violation diagnostics | 2026-09-30 19:00 | 09df4ea65ee297a56ef6fdef1e6207c2d052da3f908229ed359a84fd16b9ced8 |

Caches read by this report

*Note:* Files in data/tidy: Sampling draws:
sampling_models_n5000_full.rds; PSA results: psa_obj_correct.rds;
Violation diagnostics: pfs_os_violations_n5000.rds.

# Results

## How many draws are affected

| Quantity                                       |        Value |
|:-----------------------------------------------|-------------:|
| Sampling-cache draws evaluated                 |        5,000 |
| Draws with PFS \> OS on at least one curve     |        1,160 |
| Share of draws affected (Monte Carlo SE)       | 23.2% (0.6%) |
| Retained PSA rows                              |        5,000 |
| PSA rows whose draw has at least one violation |        1,160 |
| Share of PSA rows affected                     |        23.2% |

Draws and PSA rows with a PFS \> OS violation

**23.2% of the 5,000 coefficient draws** put PFS above OS on at least
one of the five curves, and 23.2% of the retained PSA rows come from
such a draw. The clamp is therefore not a rare safety net but acts in
roughly a quarter to a half of the PSA.

| Curve     | Draws violating | Share of draws | MC SE |
|:----------|----------------:|---------------:|------:|
| SoC       |             155 |           3.1% | 0.25% |
| CRP+      |             653 |          13.1% | 0.48% |
| CRP-      |              99 |           2.0% | 0.20% |
| TMB/BRAF+ |             923 |          18.5% | 0.55% |
| TMB/BRAF- |              60 |           1.2% | 0.15% |

Draws with PFS \> OS, by curve

*Note:* Each draw contributes one OS and one PFS coefficient vector
shared by all five curves, so the same draw can violate on several
curves. MC SE: binomial Monte Carlo standard error of the share over
5,000 draws.

| Curves violating in the draw | Draws | Share |
|-----------------------------:|------:|------:|
|                            0 | 3,840 | 76.8% |
|                            1 |   585 | 11.7% |
|                            2 |   474 |  9.5% |
|                            3 |    58 |  1.2% |
|                            4 |    32 |  0.6% |
|                            5 |    11 |  0.2% |

Number of curves (out of five) with PFS \> OS per draw

## Where and how far the curves cross

| Curve | Violating draws | Points, median \[IQR\] | Points, max | First week | Last week |
|:---|---:|---:|---:|---:|---:|
| SoC | 155 | 166 \[66, 244\] | 372 | 322 (1) | 520 (520) |
| CRP+ | 653 | 106 \[25, 183\] | 474 | 355 (1) | 520 (520) |
| CRP- | 99 | 34 \[2, 145\] | 350 | 307 (1) | 520 (520) |
| TMB/BRAF+ | 923 | 145 \[52, 222\] | 520 | 342 (1) | 520 (520) |
| TMB/BRAF- | 60 | 147 \[67, 233\] | 372 | 306 (1) | 520 (520) |

Where the curves cross, among violating draws

*Note:* Points: number of the 521 weekly grid points (weeks 0 to 520) at
which PFS exceeds OS. First week: median (minimum) of the first
violating week; last week: median (maximum) of the last violating week.

| Curve     | Violating draws | Median | 95th percentile | Maximum |
|:----------|----------------:|-------:|----------------:|--------:|
| SoC       |             155 | 0.0016 |          0.0226 |  0.0648 |
| CRP+      |             653 | 0.0012 |          0.0170 |  0.1053 |
| CRP-      |              99 | 0.0001 |          0.0021 |  0.0108 |
| TMB/BRAF+ |             923 | 0.0011 |          0.0131 |  0.0673 |
| TMB/BRAF- |              60 | 0.0004 |          0.0175 |  0.0410 |

Largest PFS - OS within the draw (survival-probability units), among
violating draws

<img src="pfs_os_violations_files/figure-commonmark/timing-plot-1.png"
style="width:100.0%" data-fig-align="center"
alt="Share of the 5,000 draws in which PFS exceeds OS at each week of the horizon, by curve. The share rises through the horizon and is highest at its end (week 520): SoC 3%; CRP+ 10%; CRP- 1%; TMB/BRAF+ 17%; TMB/BRAF- 1%." />

<img
src="pfs_os_violations_files/figure-commonmark/mean-excess-plot-1.png"
style="width:100.0%" data-fig-align="center"
alt="Mean excess of PFS over OS at each week, averaged over all draws (draws without a violation at that week contribute zero). The area under each line, in weeks, is 52 times the all-draw mean excess area of that curve in the state-areas table." />

## PFS-over-OS excess area before and after the clamp

| Curve     |   Mean | Median | 95th percentile | Maximum | After clamp |
|:----------|-------:|-------:|----------------:|--------:|------------:|
| SoC       | 0.0206 | 0.0029 |          0.0972 |  0.3557 |           0 |
| CRP+      | 0.0130 | 0.0014 |          0.0601 |  0.5647 |           0 |
| CRP-      | 0.0012 | 0.0000 |          0.0044 |  0.0389 |           0 |
| TMB/BRAF+ | 0.0112 | 0.0022 |          0.0486 |  0.3661 |           0 |
| TMB/BRAF- | 0.0109 | 0.0005 |          0.0547 |  0.1843 |           0 |

PFS-over-OS excess area (patient-years) among violating draws, before
and after the clamp

*Note:* Area = integral of max(PFS - OS, 0) over the 10-year horizon,
trapezoidal rule on the weekly grid, in patient-years per patient. After
the clamp PFS = min(PFS, OS) and the area is zero for every draw.

| State            | Curve     | Before clamp | After clamp |  Change |
|:-----------------|:----------|-------------:|------------:|--------:|
| Progression-free | CRP+      |       1.8703 |      1.8686 | -0.0017 |
| Progression-free | CRP-      |       0.9677 |      0.9677 | -0.0000 |
| Progression-free | SoC       |       1.0029 |      1.0022 | -0.0006 |
| Progression-free | TMB/BRAF+ |       1.4970 |      1.4949 | -0.0021 |
| Progression-free | TMB/BRAF- |       0.8694 |      0.8693 | -0.0001 |
| Progressed       | CRP+      |       0.9901 |      0.9918 |  0.0017 |
| Progressed       | CRP-      |       0.9423 |      0.9423 |  0.0000 |
| Progressed       | SoC       |       1.2528 |      1.2535 |  0.0006 |
| Progressed       | TMB/BRAF+ |       0.7797 |      0.7818 |  0.0021 |
| Progressed       | TMB/BRAF- |       1.2810 |      1.2811 |  0.0001 |

Mean state areas over all draws (patient-years), before and after the
clamp

*Note:* Means over all 5,000 draws (non-violating draws are unchanged by
the clamp). Progression-free area = integral of PFS (min(PFS, OS) after
the clamp); progressed area = integral of OS - PFS (negative on the
violating stretch before the clamp; max(OS - PFS, 0) after). The change
is the all-draw mean excess area, moved from the progression-free to the
progressed state; the OS area (SoC 2.256, CRP+ 2.860, CRP- 1.910,
TMB/BRAF+ 2.277, TMB/BRAF- 2.150 patient-years) is unchanged.

Averaged over all draws, the clamp moves 0.0006 patient-years from the
progression-free to the progressed state on the standard-of-care curve,
against a mean progressed-state area after the clamp of 1.253
patient-years. The largest excess area on any curve is 0.565
patient-years (CRP positive (experimental)).

<img
src="pfs_os_violations_files/figure-commonmark/area-distribution-plot-1.png"
style="width:100.0%" data-fig-align="center"
alt="Distribution of the PFS-over-OS excess area among violating draws, by curve (log scale)." />

<img src="pfs_os_violations_files/figure-commonmark/example-plot-1.png"
style="width:100.0%" data-fig-align="center"
alt="The draw with the largest PFS-over-OS excess area for each curve (draws 1367, 4511, 1256, 3722, 3528). Solid lines are the raw sampled OS and PFS; the shaded region is the excess area, the trough where PFS lies above OS; the dashed line is the clamped PFS = min(PFS, OS) that the model uses." />

## QALY consequence of the clamp

| Curve     | All draws | Violating draws | Maximum | Clamped QALYs | Relative |
|:----------|----------:|----------------:|--------:|--------------:|---------:|
| SoC       |    0.0001 |          0.0022 |  0.0366 |         1.378 |    0.00% |
| CRP+      |    0.0002 |          0.0014 |  0.0628 |         1.800 |    0.01% |
| CRP-      |    0.0000 |          0.0001 |  0.0040 |         1.197 |    0.00% |
| TMB/BRAF+ |    0.0002 |          0.0012 |  0.0380 |         1.449 |    0.02% |
| TMB/BRAF- |    0.0000 |          0.0012 |  0.0198 |         1.307 |    0.00% |

Discounted QALYs the raw curves would add relative to the clamped
curves, by curve

*Note:* Raw minus clamped discounted QALYs per patient = sum over cycles
of max(PFS - OS, 0) x (u_np - u_p) x cycle length x discount weight x
trapezoidal endpoint weight; mean over all draws and over violating
draws, and the maximum. Clamped QALYs: mean over all draws. Relative:
all-draw mean difference as a share of the mean clamped QALYs of the
curve.

| Strategy | Mean overstatement | Maximum | Mean PSA QALYs |
|:---------|-------------------:|--------:|---------------:|
| Control  |             0.0001 |  0.0366 |          1.378 |
| CRP      |             0.0001 |  0.0162 |          1.401 |
| TMB/BRAF |             0.0001 |  0.0176 |          1.370 |

QALY overstatement of the raw curves by strategy (per patient,
discounted)

*Note:* Strategy values use each PSA draw’s joint-population weights,
marginal prevalences and utilities. Mean PSA QALYs are the cached PSA
means (clamped curves).

| Strategy | Mean shift | Largest downward | Largest upward |
|:---------|-----------:|-----------------:|---------------:|
| CRP      |     0.0000 |          -0.0162 |         0.0360 |
| TMB/BRAF |    -0.0000 |          -0.0176 |         0.0317 |

Shift of the incremental QALYs versus standard of care caused by the
clamp (clamped minus raw, per patient)

*Note:* Per draw, the clamp’s effect on the strategy’s QALYs minus its
effect on standard of care. A negative value means the clamp makes the
guided strategy look worse relative to standard of care than the raw
curves would.

The clamp lowers the QALYs of every strategy, because the excess is
progression-free time that becomes progressed time at the lower sampled
utility. The mean overstatement the raw curves would have produced is
0.0001 QALYs for standard of care, 0.0001 for the CRP-guided strategy
and 0.0001 for the TMB/BRAF-guided strategy. What matters for the
decision is the effect on the increment: the clamp shifts the mean
incremental QALYs of the CRP-guided strategy versus standard of care by
0.0000 and of the TMB/BRAF-guided strategy by -0.0000, with a per-draw
range from -0.0162 to 0.0360 for CRP.

## Is the clamp why the PSA increments differ from the base case?

| Strategy | Base case | PSA, clamped | PSA, raw |
|:---------|----------:|-------------:|---------:|
| Control  |    1.3400 |       1.3775 |   1.3776 |
| CRP      |    1.3775 |       1.4007 |   1.4008 |
| TMB/BRAF |    1.3389 |       1.3696 |   1.3697 |

Discounted QALYs per patient: deterministic base case, PSA mean with the
clamp (the cached PSA), and PSA mean the raw curves would have given

*Note:* PSA, raw = cached PSA mean + mean overstatement of the raw
curves for the strategy.

| Strategy | Base case | PSA, clamped | PSA, raw |
|:---------|----------:|-------------:|---------:|
| CRP      |    0.0375 |       0.0232 |   0.0232 |
| TMB/BRAF |   -0.0011 |      -0.0079 |  -0.0079 |

Incremental discounted QALYs versus standard of care: base case, PSA
mean with the clamp, and PSA mean without it

The PSA incremental QALYs of the CRP-guided strategy are 0.0232,
compared with a deterministic increment of 0.0375. Removing the clamp
would give 0.0232: the clamp narrows the gap by an amount equal to 0.0%
of the reported gap. The remaining difference reflects nonlinear
averaging of survival curves over coefficient uncertainty and Monte
Carlo error. Cross-endpoint dependence cannot change the expected
raw-curve QALYs when the endpoint marginal distributions and independent
utility/population distributions are held fixed. Earlier wording
attributed the whole gap to endpoint independence; that causal
attribution was too strong.

# Interpretation

- **The crossings are frequent but shallow.** About 23% of draws cross
  on at least one curve, yet the median violating draw crosses by at
  most 0.0010 in survival probability and carries an excess area of
  0.0015 patient-years. The crossings sit in the extrapolated tail where
  both curves are close to their floor, so the area involved is small
  relative to the restricted OS area of about 2.29 patient-years.
- **The clamp is a one-sided correction.** min(PFS, OS) accepts the OS
  draw as given and truncates PFS; it never raises OS. It lowers QALYs
  when the progression-free utility exceeds the progressed utility,
  without targeting the deterministic base case.
- **The consequence is differential but small.** The five curves cross
  at different rates (SoC 3%; CRP+ 13%; CRP- 2%; TMB/BRAF+ 18%;
  TMB/BRAF- 1%), so the clamp moves the increments of the guided
  strategies versus standard of care, not only their levels. The mean
  shift on the CRP increment is 0.0000 QALYs against a base-case
  incremental gain of 0.0375. Individual draws can move by up to 0.036
  QALYs, substantially more than the mean shift.
- **Report both estimands.** The deterministic CRP increment is 0.0375
  QALYs; the PSA mean is 0.0232 with the clamp and 0.0232 without it. A
  nonlinear model can produce this difference under correctly centred
  coefficient uncertainty. The original five-SE diagnostic is retained;
  the separate zero-uncertainty regression verifies equality when all
  uncertainty is removed.
- **Joint covariance is implemented.** Issue \#159 uses the first
  proposed fix: a paired patient bootstrap estimates dependence for
  joint normal coefficient draws. It does not guarantee ordered curves.
  A survival model parameterised to enforce ordering would be a separate
  structural change; correlation or conditioning alone is insufficient
  without that constraint.

# Notes

- The diagnostics are cached in
  `data/tidy/pfs_os_violations_n{n_samples}.rds`;
  `run_pfs_os_violation_diagnostics()` regenerates the cache on any
  fingerprint mismatch. The current cache was built from sampling cache
  09df4ea65ee297a56ef6fdef1e6207c2d052da3f908229ed359a84fd16b9ced8
  (method mvn_joint_v2) in 852 seconds on 2026-09-30 19:00.
- Curves are computed in closed form from the gamma parameters of each
  draw; `validate_direct_curves()` checks them against the pipeline’s
  `predict()`-based helper before every regeneration and stops if they
  differ by more than 1e-8.
- All values are per patient over the 10-year horizon; areas and QALY
  differences use the trapezoidal rule.
- “PFS-over-OS excess area” is this report’s term for the area between
  the two curves on the stretch where they are inverted; it is not a
  standard quantity. The figures call the same region the trough.

------------------------------------------------------------------------

**Report completed on:** 2026-09-30  
**Repository:** ben-geisler/METIMMOX-1  
**Report version:** 2.1
