# Issue #181: off-study deaths and PFS ascertainment

Implementation date: 24 September 2026. Only aggregate trial results are recorded.

The attached verification plan reports 14 deaths without recorded progression
in the 68-patient cohort (12 after adverse-event exit; 10 control, 4 experimental),
a median assessment-to-death gap of 64 weeks (range 5.7-119.7), and 13 gaps of at
least 24 weeks. The independent SAM cross-check below reproduces the plan's
endpoint agreement counts. The AE-exit and gap summaries are carried from the plan.

## SAM endpoint cross-check (25 September 2026)

The user clarified that SAM OS and SAM PFS are Excel-export tuples of
`(days until endpoint, event reached)`, not additional biomarkers. The supplied
example reproduces both parsed tuples and agrees with our derived OS endpoint.
Its SAM PFS date differs from both Days until progression and Days until last
evaluation: it is one of the nine progressor-date discrepancies. The analysis
retains the recorded Days until progression time, as specified by the plan.
`check_sam_endpoints.R` checks headers and tuple syntax, applies the same two
exclusions as script 01 (74 eligible export rows), and matches the 68 analysis
patients by ID. No patient-level records are written or printed.

| Comparison | Exact time-and-event agreement | Time only differs | Event only differs | Both differ |
|---|---:|---:|---:|---:|
| SAM OS vs export OS (74 eligible rows) | 71 | 0 | 2 | 1 |
| SAM OS vs analysis OS (68 patients) | 65 | 0 | 2 | 1 |
| SAM PFS vs primary 16-week/TTPwk rule | 46 | 17 | 3 | 2 |
| SAM PFS vs former all-deaths rule | 42 | 14 | 1 | 11 |
| SAM PFS vs 16-week/LastEvalwk alternative | 51 | 12 | 2 | 3 |

Among the 14 deaths without recorded progression, SAM PFS censors 11 and our
primary rule censors 13: both censor 10, SAM counts three as events that the
primary rule censors, and SAM censors one that the primary rule counts as an
event. Nine recorded progressors have different SAM dates, six equal to the
last-evaluation date. Across the 74 eligible export rows, 14 of 15 SAM-censored
deaths use the last-evaluation date. Machine-readable aggregate results are in
`sam_endpoint_comparison.csv` and `sam_reconciliation_counts.csv`.

The SAM OS audit corrects one detail in the attached plan: all three
discrepancies censor recorded deaths, but only one SAM date is earlier; two
have the same time and a different event flag.

This is an independent consistency check, not complete agreement or proof that
one anchor is correct. The primary rule remains as specified in the approved
plan; the example check asserts the supplied SAM tuples and independently
asserts that the derived progression endpoint uses the raw progression time; neither SAM endpoint is silently substituted for the analysis endpoint.

The implemented primary rule uses a 16-week inclusive death window (two scheduled
8-week CT intervals), anchored at TTPwk (Days until progression). LastEvalwk
(Days until last evaluation) is selectable. Progression-exit times are preserved.
An infinite window recovers the former rule; a missing non-progressor assessment
yields missing PFS time and event status. The clinical report compares both anchors,
8/16/24-week windows, and the former all-deaths rule on identical patient IDs.

The plan corrects the issue's preview attribution: the approximate CRP interaction
HR 0.26 comes from censoring at Days until progression; Days until last evaluation
gives approximately 0.36. The clinical report computes these values afresh.

## Questions for the data provider (not sent)

- Which source fields and censoring rules generated the SAM tuples? Their
  (days, event) encoding is now known, but the 26 differences from the former
  derived PFS and 22 differences from the new primary PFS still need explanation.
- Why do nine recorded progressors have different SAM PFS dates, six equal to
  Days until last evaluation?
- Does Days until last evaluation represent the last imaging assessment or the
  last clinical evaluation? Which field identifies the last adequate assessment?
- How should five post-exit PD reads with Progression exit = 0 be interpreted?
- Why are the two non-AE deaths without progression SAM PFS events at the
  Days until progression date? Is a clinical-failure flag missing from the export?
- Why does SAM OS censor three recorded deaths (two at the same exported OS
  time and one at an earlier time)?

## Execution record

- Baseline saved before source edits: `snapshot_181_baseline_2a48ae8.rds` and
  `psa_181_baseline_2a48ae8.rds`.
- Synthetic endpoint test passes. Data preparation re-created the ignored RDS
  with header/type/range checks and LastEvalwk.
- Live distribution audit (`check_selection.R`): 68 complete cases, 50 PFS events;
  infinite window restores 63 events on the same cohort. Selected OS/PFS remains
  gamma/gamma, combined-AIC ordering penalty 4.431952.
- Landmark contract checks pass: DAG TLR subset week-9 PFS n = 59, scan-date n = 55,
  week-12 n = 54; week 9 retains 3 first-scan progressors and 19 later scans.
- Clinical report initial PDF/GFM render passed, including the asserted primary
  model/sensitivity fit agreement. CRP x Rx PFS: primary 0.26 (0.06-1.34),
  p = 0.102; all deaths 0.42 (0.11-1.68), p = 0.213; alternative anchor
  0.36 (0.09-1.78), p = 0.197. TMB/BRAF primary 0.65 (0.18-2.41), p = 0.516.
- Clinical TLR PFS cohorts: week 9 n = 56 (4 censored before landmark),
  scan date n = 52 (5 censored), week 12 n = 51 (6 censored). All retain the
  same 65-patient OS cohort. PFS TLR interaction HRs are 0.88, 0.75 and 0.85.
- The report review also corrected two inherited interpretations exposed by the
  regenerated results: per-patient scan-date resetting can change the OS risk
  sets even without early exclusions; the current Schoenfeld tests do not
  identify biomarker main effects as the sole cause of the global PFS departure.
- `test_survival_ordering.R`, `test_psa_zero_uncertainty.R`,
  `test_publication_artifacts.R`, `test_cache_provenance.R`, and
  `test_prediction_population_live.R` pass. The live population test used its
  own temporary sampling cache and checked weighted diagnostic curves against
  the regular prediction path and PSA QALYs with repeated model indices.
- Source implementation commit: `791443e`; regeneration started after that commit.
- Sampling, PSA, EVPPI and scenario regeneration completed. Of 1,000 paired
  bootstrap attempts, 996 succeeded and four failed. Base-case EVPI remains
  EUR 0 per patient at WTP EUR 51,000; the base EVPPI cache has no rows.
- Fixed snapshot: `snapshot_181_fixed_791443e.rds` and
  `psa_181_fixed_791443e.rds`; snapshot comparison passes and confirms different
  PSA fingerprints. `results.txt` records fingerprints and aggregate outcomes.
- Both ordering-diagnostic and extrapolation reports rendered PDF then GFM.
- Full report rendering exposed an inherited EVPPI table schema mismatch after
  the shared PSA summary acquired more columns. The report now explicitly
  selects its four displayed columns (version 5.3); PDF/GFM renders pass.
- Completed 25 September 2026: all 20 reports rendered PDF then GFM, and all
  26 publication vignettes rendered sequentially. The local `_site/` contains
  20 PDFs; nothing was pushed or published remotely.
- Final report/cache contracts: PASS, 576 checks, one inapplicable empty-base-
  EVPPI row block skipped because EVPI is zero. All eight tests specified in
  the plan passed, as did the SAM tuple audit and snapshot comparison.
- Refreshed ordering diagnostics: 1,160/5,000 draws (23.2%) cross on at least
  one subgroup curve; the direct gamma predictions agree with the regular
  prediction path to 3.3e-16.
- The final commit includes the source/report corrections, regenerated tracked
  outputs, fixed snapshots, and aggregate SAM audit, with `Fixes #181` in its
  message. The fixed snapshots retain the pre-regeneration implementation
  commit `791443e` for provenance.

The ESMO poster is out of scope and retains clinical report v3.7 results. The
issue #160 note, remote push, and sending questions remain with the user.
