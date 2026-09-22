# Combined fixes: #173, #171 and #163

The baseline snapshot was taken before editing the calculation sources. Both
snapshots use issue number 173 and precede the batch commit, so their filenames
carry the original HEAD, `1b50e8d`. The
[impact report](../../reports/technical/bug_fix_impact.md) identifies this as a
combined batch rather than attributing all outcome changes to diagnostic prices.

## Changes

- **#173:** Diagnostic charges resolve the canonical `c_test_CRP` or
  `c_test_NGS` scalar at calculation time. The model no longer requires the
  duplicate lookup, and direct, base-case, enriched, DSA and PSA callers need
  no synchronization step. A legacy lookup cannot override the scalar.
- **#171:** Initial and replacement PSA calls use the same completeness check:
  exactly one finite cost and effect for each strategy. Invalid outcomes count
  toward the 2% initial-failure threshold and receive replacement attempts.
  PSA requires a draw index; supplied-curve evaluation is explicitly `curves`.
  The SMDM poster uses that mode.
- **#163:** QALYs and ongoing progressed-state cost rates use trapezoidal
  integration. Quarterly rates accrue as four quarters per patient-year.
  Scheduled administrations, tests and visits, baseline/screening charges and
  incremental-death costs retain their full time-point charges. Ordering
  diagnostics use the same QALY integration convention.

## Cache recomputation

[recompute.R](recompute.R) uses the unchanged 5,000 paired survival coefficient
draws and saved economic/population draws. It verifies that those parameter
values equal the current seeded draws and disables sampling regeneration.
The standard prediction path estimated approximately one hour per PSA run on
this machine. The existing closed-form gamma helper supplies identical
population-averaged curves to the explicit supplied-curve mode more quickly.

Validation of this route includes all five curves and both endpoints for 12
draws spread across the cache, comparison of corrected outcomes against normal
PSA mode for both price scenarios on those draws, and reproduction of **every
pre-fix PSA row** with the original calculation functions from `1b50e8d`.
The script stops unless cost differences are below EUR 1e-7 and QALY
differences below 1e-11; curve comparisons use 1e-12. It checks each new outcome
for completeness, rebuilds the PSA pairs and scenario results, runs the
canonical EVPPI functions, refreshes ordering diagnostics, and verifies that
the sampling-cache MD5 is unchanged. Recalculation provenance is stored in
each PSA object's `recomputation` field.
The numerical checks are also saved in [recomputation_checks.json](recomputation_checks.json).

Observed maximum discrepancies were 5.55e-16 per survival probability,
EUR 2.91e-11 and 4.44e-16 QALYs across all pre-fix PSA rows, and
EUR 1.46e-11 and 2.22e-16 QALYs against standard corrected PSA mode.
All 5,000 draws produced complete finite outcomes for both price scenarios.

The completed ordering diagnostics reproduce every PSA strategy's QALYs across
all 5,000 draws (maximum difference 8.88e-16). Run
`Rscript validation/issue173_batch_2026-09-21/check_outputs.R` to check this
agreement and the snapshot/cache hashes. The calculation driver encountered a
trailing parse error after saving all caches and verifying the sampling hash,
because its final output statement was edited during execution. The saved
script parses cleanly; the separate output checks passed.
The ordering report uses `l_params_base$prediction_population`, asserting that
its rows and retained columns equal the broader complete-case data. This
avoids redundant regeneration caused by unused clinical columns in the cache
fingerprint.

## Deterministic impact

| Strategy | Cost change (EUR) | QALY change |
|---|---:|---:|
| Standard of care | -0.02251549 | -0.007025276 |
| CRP-guided | -0.01891231 | -0.007026462 |
| TMB/BRAF-guided | -0.02058024 | -0.007026925 |

Standard of care remains preferred at EUR 51,000/QALY. The CRP frontier ICER
is EUR 1,578,614/QALY; TMB/BRAF remains dominated. Changes to default outcomes
come from interval integration; the other fixes protect changed-price and
invalid-call/draw paths.

## Regression coverage

Focused regressions cover exact 1-, 10- and 20-year horizons; discounted
nonconstant occupancy; zero-duration grids; full scheduled and event costs;
diagnostic-price changes through direct and PSA callers; non-finite, duplicate
and incomplete outcomes at 3% (failure) and 2% (replacement); invalid indices
and modes; and explicit supplied-curve ordering behavior.

The survival-ordering, zero-uncertainty PSA, canonical-prevalence, model-config,
MVN sampling, cache-provenance and sampling-specification checks also pass.
The historical five-Monte-Carlo-SE comparison of uncertain PSA means with
fitted-parameter results remains a separate diagnostic; this batch does not
change its threshold or recenter the coefficient draws.

Rerunning that historical diagnostic returns failure, as before: four of six
level comparisons exceed five SE. The QALY mean differences are +0.03680,
+0.02255 and +0.03034 (8.77, 5.94 and 7.80 SE) for control, CRP and TMB/BRAF;
TMB/BRAF cost differs by -6.11 SE. The structural control-curve check passes.

All **301 report/cache contract checks pass**, with one inapplicable EVPPI-row
block skipped because base-case EVPI is zero. Both #173 snapshots and their
comparison have been saved. Their PSA files have different hashes, while their
sampling-cache provenance is unchanged.
