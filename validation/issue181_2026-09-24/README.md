# Issue #181: off-study deaths and PFS ascertainment

Implementation date: 24 September 2026. Only aggregate trial results are recorded.

The attached verification plan reports 14 deaths without recorded progression
in the 68-patient cohort (12 after adverse-event exit; 10 control, 4 experimental).
It reports a median assessment-to-death gap of 64 weeks (range 5.7–119.7), with
13 at least 24 weeks. Its export comparison found SAM OS agreement in 71/74 rows,
SAM PFS agreement in 42/68, censoring of 11/14 deaths without progression,
14/15 censored deaths at Days until last evaluation, and different dates for
nine progressors. These export reconciliation findings are carried from the
plan; they are not a decision to adopt the undocumented SAM columns.

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

- What are the provenance and definitions of SAM OS and SAM PFS, and how do the
  26 differences from the former derived PFS arise?
- Why do nine recorded progressors have different SAM PFS dates, six equal to
  Days until last evaluation?
- Does Days until last evaluation represent the last imaging assessment or the
  last clinical evaluation? Which field identifies the last adequate assessment?
- How should five post-exit PD reads with Progression exit = 0 be interpreted?
- Why are the two non-AE deaths without progression SAM PFS events at the
  Days until progression date? Is a clinical-failure flag missing from the export?
- Why does SAM OS censor three deaths at earlier dates?

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
- Full regeneration, fixed snapshot, renders and remaining checks: pending.

The ESMO poster is out of scope and retains clinical report v3.7 results. The
issue #160 note, remote push, and sending questions remain with the user.
