# Issues #36, #164 and #56 (batch, 30 September 2026)

## Decisions (user, 30 September 2026)

- **#36 visit costs:** Norwegian DRG (ISF) weight 0.047 for every outpatient
  visit type times the 2023 somatic-care unit price NOK 49,484 (EUR 4,541.60)
  = NOK 2,325.75 = EUR 213.46. `c_other_visit`, `c_other_baseline` and
  `c_other_follow` all take this value (previously EUR 33, 530 and 33 without
  a source) and stay separate parameters.
- **#164 monitoring:** a visit with every CT scan (`visit_schedule()`), which
  adds a surveillance visit every 12 weeks after the last administration
  (positions 49, 61, ..., 517) for progression-free patients; progressed
  patients receive one CT per quarter (`progressed_imaging_costs`). Point 1
  of the issue (12- vs 13-week cadence) had already gone with #163.
- **#56 CVs:** CV 0.20 (resource-use costs) and 0.15 (`u_np`, `u_decrement`)
  are documented as assumptions without a source; no new analysis.

## Sequence

1. Baseline snapshot `snapshot_36_baseline_08ec4a9` before any source edit.
   The sampling cache was stale on entry (`data_hash`/`data_columns`, from the
   #185 C-order fingerprints), so this run regenerated sampling, PSA and EVPPI
   under the unchanged model.
2. Fix committed as `e5ff9c1`; synthetic tests pass
   (`test_interval_integration`, `test_model_input_contract`,
   `test_survival_ordering`, `test_biomarker_test_cost_mapping`,
   `test_canonical_prevalence`, `test_sampling_rework`,
   `test_sim_idx_validation`, `test_model_config_coupling`).
3. `run_all.ps1`: scripts 06-12 (`regenerate.R`; sampling reused), all
   reports (`publish/publish_reports.R`), all publication vignettes,
   `tests/test_report_contracts.R` (PASS: 593 checks, 1 skipped: empty
   base-EVPPI rows because EVPI = 0), fixed snapshot
   `snapshot_36_fixed_e5ff9c1` (reused the regenerated PSA), 
   `tests/compare_snapshots.R 36`, `bug_fix_impact.qmd`.

## Results

| | Control | CRP-guided | TMB/BRAF-guided |
|---|---:|---:|---:|
| Deterministic cost, before | 20,967 | 53,686 | 62,382 |
| Deterministic cost, after | 25,764 | 58,203 | 66,916 |
| PSA mean cost, before | 20,931.85 | 53,228.19 | 61,751.88 |
| PSA mean cost, after | 25,732.03 | 57,739.78 | 66,302.48 |

QALYs are identical before and after (deterministic and every PSA draw means:
1.377530 / 1.400750 / 1.369621). CRP pairwise ICER: deterministic EUR 871,972
-> 864,492; probabilistic ratio of mean increments EUR 1,390,876 -> 1,378,448.
TMB/BRAF remains dominated by control, control remains optimal, EVPI remains
EUR 0 at WTP EUR 51,000. `compare_snapshots.R` reports "No parameter changes
detected" because its metadata comparison does not track the unit costs; the
cost changes are visible in the result tables.
