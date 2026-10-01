# Bug Fix Impact Analysis
Ben Geisler
2026-10-01

- [Snapshot Inventory](#snapshot-inventory)
- [Impact Comparisons](#impact-comparisons)
  - [Issue 36: baseline vs fixed](#issue-36-baseline-vs-fixed)
  - [Issue 145: baseline vs fixed](#issue-145-baseline-vs-fixed)
  - [Issue 146: baseline vs fixed](#issue-146-baseline-vs-fixed)
  - [Issue 147: baseline vs fixed](#issue-147-baseline-vs-fixed)
  - [Issue 149: baseline vs fixed](#issue-149-baseline-vs-fixed)
  - [Issue 151: baseline vs fixed](#issue-151-baseline-vs-fixed)
  - [Issue 152: baseline vs fixed](#issue-152-baseline-vs-fixed)
  - [Issue 153: baseline vs fixed](#issue-153-baseline-vs-fixed)
  - [Issue 154: baseline vs fixed](#issue-154-baseline-vs-fixed)
  - [Issue 155: baseline vs fixed](#issue-155-baseline-vs-fixed)
  - [Issue 156: baseline vs fixed](#issue-156-baseline-vs-fixed)
  - [Issue 157: baseline vs fixed](#issue-157-baseline-vs-fixed)
  - [Issue 159: baseline vs fixed](#issue-159-baseline-vs-fixed)
  - [Issue 165: baseline vs fixed](#issue-165-baseline-vs-fixed)
  - [Issue 166: baseline vs fixed](#issue-166-baseline-vs-fixed)
  - [Issue 168: baseline vs fixed](#issue-168-baseline-vs-fixed)
  - [Issue 172: baseline vs fixed](#issue-172-baseline-vs-fixed)
  - [Batch \#173 / \#171 / \#163: baseline vs
    fixed](#batch-173--171--163-baseline-vs-fixed)
  - [Issue 181: baseline vs fixed](#issue-181-baseline-vs-fixed)
  - [Cumulative: issue 145 baseline vs issue 165
    fixed](#cumulative-issue-145-baseline-vs-issue-165-fixed)
- [Interpretation](#interpretation)
  - [Validation warnings](#validation-warnings)

Issue \#181 changes PFS ascertainment for deaths without recorded
progression: deaths more than 16 weeks after the last assessment are
censored at that assessment (primary proxy: `TTPwk`). Recorded
progression times and OS are unchanged. The clinical report compares the
former all-deaths rule and the alternative `LastEvalwk` anchor. The
baseline/fixed snapshots below capture the resulting economic changes
after refitting and regenerating every downstream cache.

# Snapshot Inventory

| filename                          | issue | status   | commit  | mtime               |
|:----------------------------------|:------|:---------|:--------|:--------------------|
| snapshot_145_baseline_ce8ec08.rds | 145   | baseline | ce8ec08 | 2026-07-12 11:55:14 |
| snapshot_145_fixed_ea4431f.rds    | 145   | fixed    | ea4431f | 2026-07-12 19:54:11 |
| snapshot_146_baseline_a3befd2.rds | 146   | baseline | a3befd2 | 2026-07-27 13:57:12 |
| snapshot_146_fixed_7f76f9d.rds    | 146   | fixed    | 7f76f9d | 2026-07-27 13:57:12 |
| snapshot_147_baseline_7f76f9d.rds | 147   | baseline | 7f76f9d | 2026-07-27 13:57:12 |
| snapshot_147_fixed_d1f203f.rds    | 147   | fixed    | d1f203f | 2026-07-27 13:57:12 |
| snapshot_149_baseline_8f21e24.rds | 149   | baseline | 8f21e24 | 2026-09-08 08:51:07 |
| snapshot_149_fixed_4140474.rds    | 149   | fixed    | 4140474 | 2026-09-08 08:51:07 |
| snapshot_151_baseline_89d7a28.rds | 151   | baseline | 89d7a28 | 2026-09-08 08:51:07 |
| snapshot_151_fixed_70d672d.rds    | 151   | fixed    | 70d672d | 2026-09-08 08:51:07 |
| snapshot_152_baseline_e27b60f.rds | 152   | baseline | e27b60f | 2026-09-08 08:51:07 |
| snapshot_152_fixed_9cb460b.rds    | 152   | fixed    | 9cb460b | 2026-09-08 08:51:07 |
| snapshot_153_baseline_66a995c.rds | 153   | baseline | 66a995c | 2026-09-08 08:51:07 |
| snapshot_153_fixed_3b291ff.rds    | 153   | fixed    | 3b291ff | 2026-09-08 08:51:07 |
| snapshot_154_baseline_3d5938d.rds | 154   | baseline | 3d5938d | 2026-09-08 08:51:07 |
| snapshot_154_fixed_1f7202c.rds    | 154   | fixed    | 1f7202c | 2026-09-08 08:51:07 |
| snapshot_155_baseline_331a1a2.rds | 155   | baseline | 331a1a2 | 2026-09-08 08:51:07 |
| snapshot_155_fixed_9e29ea6.rds    | 155   | fixed    | 9e29ea6 | 2026-09-08 08:51:07 |
| snapshot_156_baseline_2806457.rds | 156   | baseline | 2806457 | 2026-09-08 08:51:07 |
| snapshot_156_fixed_8a89056.rds    | 156   | fixed    | 8a89056 | 2026-09-08 08:51:07 |
| snapshot_157_baseline_efcd775.rds | 157   | baseline | efcd775 | 2026-09-08 08:51:07 |
| snapshot_157_fixed_cdf852f.rds    | 157   | fixed    | cdf852f | 2026-09-08 08:51:07 |
| snapshot_159_baseline_6ca19ab.rds | 159   | baseline | 6ca19ab | 2026-09-18 16:40:08 |
| snapshot_159_fixed_6ca19ab.rds    | 159   | fixed    | 6ca19ab | 2026-09-21 09:46:55 |
| snapshot_165_baseline_6499092.rds | 165   | baseline | 6499092 | 2026-10-01 12:38:19 |
| snapshot_165_fixed_6499092.rds    | 165   | fixed    | 6499092 | 2026-10-01 14:48:51 |
| snapshot_166_baseline_96d0710.rds | 166   | baseline | 96d0710 | 2026-09-18 09:55:46 |
| snapshot_166_fixed_96d0710.rds    | 166   | fixed    | 96d0710 | 2026-09-18 09:57:44 |
| snapshot_168_baseline_b8e1b02.rds | 168   | baseline | b8e1b02 | 2026-09-22 11:25:26 |
| snapshot_168_fixed_b8e1b02.rds    | 168   | fixed    | b8e1b02 | 2026-09-22 12:04:46 |
| snapshot_172_baseline_15999f9.rds | 172   | baseline | 15999f9 | 2026-09-29 12:50:41 |
| snapshot_172_fixed_15999f9.rds    | 172   | fixed    | 15999f9 | 2026-09-29 15:29:53 |
| snapshot_173_baseline_1b50e8d.rds | 173   | baseline | 1b50e8d | 2026-09-21 15:09:17 |
| snapshot_173_fixed_1b50e8d.rds    | 173   | fixed    | 1b50e8d | 2026-09-21 16:04:39 |
| snapshot_181_baseline_2a48ae8.rds | 181   | baseline | 2a48ae8 | 2026-09-24 17:07:25 |
| snapshot_181_fixed_791443e.rds    | 181   | fixed    | 791443e | 2026-09-24 19:19:42 |
| snapshot_36_baseline_08ec4a9.rds  | 36    | baseline | 08ec4a9 | 2026-09-30 15:50:08 |
| snapshot_36_fixed_e5ff9c1.rds     | 36    | fixed    | e5ff9c1 | 2026-09-30 20:03:23 |

Available single-model snapshots

# Impact Comparisons

The issue \#173 snapshot pair covers the combined fixes for **\#173,
\#171 and \#163**. Diagnostic prices now use the scalar inputs on every
call; invalid PSA outcomes enter the failure threshold and replacement
process, and supplied curves have an explicit evaluation mode. QALYs and
ongoing progressed-state costs now use trapezoidal interval integration,
while scheduled and event costs retain their full time-point charges.
These changes require new economic results but leave the survival
coefficient draws unchanged. The baseline was saved before editing the
calculation sources; both snapshots precede the batch commit, so their
filenames carry the same pre-fix HEAD identifier.

| Comparison | Before | After |
|:---|:---|:---|
| Issue 36: baseline vs fixed | snapshot_36_baseline_08ec4a9.rds | snapshot_36_fixed_e5ff9c1.rds |
| Issue 145: baseline vs fixed | snapshot_145_baseline_ce8ec08.rds | snapshot_145_fixed_ea4431f.rds |
| Issue 146: baseline vs fixed | snapshot_146_baseline_a3befd2.rds | snapshot_146_fixed_7f76f9d.rds |
| Issue 147: baseline vs fixed | snapshot_147_baseline_7f76f9d.rds | snapshot_147_fixed_d1f203f.rds |
| Issue 149: baseline vs fixed | snapshot_149_baseline_8f21e24.rds | snapshot_149_fixed_4140474.rds |
| Issue 151: baseline vs fixed | snapshot_151_baseline_89d7a28.rds | snapshot_151_fixed_70d672d.rds |
| Issue 152: baseline vs fixed | snapshot_152_baseline_e27b60f.rds | snapshot_152_fixed_9cb460b.rds |
| Issue 153: baseline vs fixed | snapshot_153_baseline_66a995c.rds | snapshot_153_fixed_3b291ff.rds |
| Issue 154: baseline vs fixed | snapshot_154_baseline_3d5938d.rds | snapshot_154_fixed_1f7202c.rds |
| Issue 155: baseline vs fixed | snapshot_155_baseline_331a1a2.rds | snapshot_155_fixed_9e29ea6.rds |
| Issue 156: baseline vs fixed | snapshot_156_baseline_2806457.rds | snapshot_156_fixed_8a89056.rds |
| Issue 157: baseline vs fixed | snapshot_157_baseline_efcd775.rds | snapshot_157_fixed_cdf852f.rds |
| Issue 159: baseline vs fixed | snapshot_159_baseline_6ca19ab.rds | snapshot_159_fixed_6ca19ab.rds |
| Issue 165: baseline vs fixed | snapshot_165_baseline_6499092.rds | snapshot_165_fixed_6499092.rds |
| Issue 166: baseline vs fixed | snapshot_166_baseline_96d0710.rds | snapshot_166_fixed_96d0710.rds |
| Issue 168: baseline vs fixed | snapshot_168_baseline_b8e1b02.rds | snapshot_168_fixed_b8e1b02.rds |
| Issue 172: baseline vs fixed | snapshot_172_baseline_15999f9.rds | snapshot_172_fixed_15999f9.rds |
| Batch \#173 / \#171 / \#163: baseline vs fixed | snapshot_173_baseline_1b50e8d.rds | snapshot_173_fixed_1b50e8d.rds |
| Issue 181: baseline vs fixed | snapshot_181_baseline_2a48ae8.rds | snapshot_181_fixed_791443e.rds |
| Cumulative: issue 145 baseline vs issue 165 fixed | snapshot_145_baseline_ce8ec08.rds | snapshot_165_fixed_6499092.rds |

Impact comparisons included in this report

Each comparison below reports the change in base case costs and QALYs,
and in net monetary benefit at the willingness-to-pay threshold, from
the *before* snapshot to the *after* snapshot. Where both snapshots
carry the EVPPI cache (snapshots taken from issue \#152 onwards), the
per-patient EVPI and the EVPPI table are compared as well; the *after*
column includes the Monte Carlo standard error of the regression
estimator. From issue \#156 onwards each snapshot records the md5 and
modification time of the PSA cache it copied, and a *fixed* snapshot
regenerates the PSA and EVPPI caches when they predate the fix commit;
the provenance line under each comparison states whether the two PSA
files are the same cache.

## Issue 36: baseline vs fixed

| Strategy | Cost Before | QALYs Before | Cost After | QALYs After | Cost Change | QALY Change |
|:---|---:|---:|---:|---:|---:|---:|
| Standard of Care | EUR 20,967 | 1.3400 | EUR 25,764 | 1.3400 | EUR 4,797 | 0.0000 |
| CRP-guided | EUR 53,686 | 1.3775 | EUR 58,203 | 1.3775 | EUR 4,517 | 0.0000 |
| TMB/BRAF-guided | EUR 62,382 | 1.3389 | EUR 66,916 | 1.3389 | EUR 4,534 | 0.0000 |

Base case impact – Issue 36: baseline vs fixed

| Strategy         | NMB Before |  NMB After | NMB Change |
|:-----------------|-----------:|-----------:|-----------:|
| Standard of Care | EUR 47,371 | EUR 42,574 | -EUR 4,797 |
| CRP-guided       | EUR 16,565 | EUR 12,048 | -EUR 4,517 |
| TMB/BRAF-guided  |  EUR 5,901 |  EUR 1,368 | -EUR 4,534 |

Net monetary benefit impact – Issue 36: baseline vs fixed

PSA cache md5 before 5226285758db337ba15d5ad131d00649 (modified
2026-09-30 15:50), after 7a04965187ff820c11f6d3ad0ccab9e8 (modified
2026-09-30 16:51).

Per-patient EVPI: EUR 0.00 before, EUR 0.00 after.

## Issue 145: baseline vs fixed

| Strategy | Cost Before | QALYs Before | Cost After | QALYs After | Cost Change | QALY Change |
|:---|---:|---:|---:|---:|---:|---:|
| Standard of Care | EUR 20,867 | 1.4141 | EUR 21,056 | 1.3498 | EUR 189 | -0.0644 |
| CRP-guided | EUR 55,026 | 1.3842 | EUR 53,791 | 1.3886 | -EUR 1,235 | 0.0043 |
| TMB/BRAF-guided | EUR 60,798 | 1.3530 | EUR 63,235 | 1.3514 | EUR 2,438 | -0.0017 |

Base case impact – Issue 145: baseline vs fixed

| Strategy         | NMB Before |  NMB After | NMB Change |
|:-----------------|-----------:|-----------:|-----------:|
| Standard of Care | EUR 51,254 | EUR 47,783 | -EUR 3,471 |
| CRP-guided       | EUR 15,570 | EUR 17,026 |  EUR 1,456 |
| TMB/BRAF-guided  |  EUR 8,207 |  EUR 5,684 | -EUR 2,523 |

Net monetary benefit impact – Issue 145: baseline vs fixed

PSA provenance: before snapshot predates issue \#156 and did not record
the PSA cache md5; after snapshot predates issue \#156 and did not
record the PSA cache md5. Identity of the two PSA files cannot be
established from metadata.

## Issue 146: baseline vs fixed

| Strategy | Cost Before | QALYs Before | Cost After | QALYs After | Cost Change | QALY Change |
|:---|---:|---:|---:|---:|---:|---:|
| Standard of Care | EUR 21,056 | 1.3498 | EUR 21,056 | 1.3498 | EUR 0 | 0.0000 |
| CRP-guided | EUR 53,791 | 1.3886 | EUR 53,791 | 1.3886 | EUR 0 | 0.0000 |
| TMB/BRAF-guided | EUR 63,235 | 1.3514 | EUR 63,235 | 1.3514 | EUR 0 | 0.0000 |

Base case impact – Issue 146: baseline vs fixed

| Strategy         | NMB Before |  NMB After | NMB Change |
|:-----------------|-----------:|-----------:|-----------:|
| Standard of Care | EUR 47,783 | EUR 47,783 |      EUR 0 |
| CRP-guided       | EUR 17,026 | EUR 17,026 |      EUR 0 |
| TMB/BRAF-guided  |  EUR 5,684 |  EUR 5,684 |      EUR 0 |

Net monetary benefit impact – Issue 146: baseline vs fixed

PSA provenance: before snapshot predates issue \#156 and did not record
the PSA cache md5; after snapshot predates issue \#156 and did not
record the PSA cache md5. Identity of the two PSA files cannot be
established from metadata.

## Issue 147: baseline vs fixed

| Strategy | Cost Before | QALYs Before | Cost After | QALYs After | Cost Change | QALY Change |
|:---|---:|---:|---:|---:|---:|---:|
| Standard of Care | EUR 21,056 | 1.3498 | EUR 21,056 | 1.3498 | EUR 0 | 0.0000 |
| CRP-guided | EUR 53,791 | 1.3886 | EUR 53,791 | 1.3886 | EUR 0 | -0.0000 |
| TMB/BRAF-guided | EUR 63,235 | 1.3514 | EUR 63,235 | 1.3514 | EUR 0 | -0.0000 |

Base case impact – Issue 147: baseline vs fixed

| Strategy         | NMB Before |  NMB After | NMB Change |
|:-----------------|-----------:|-----------:|-----------:|
| Standard of Care | EUR 47,783 | EUR 47,783 |      EUR 0 |
| CRP-guided       | EUR 17,026 | EUR 17,026 |      EUR 0 |
| TMB/BRAF-guided  |  EUR 5,684 |  EUR 5,684 |      EUR 0 |

Net monetary benefit impact – Issue 147: baseline vs fixed

PSA provenance: before snapshot predates issue \#156 and did not record
the PSA cache md5; after snapshot predates issue \#156 and did not
record the PSA cache md5. Identity of the two PSA files cannot be
established from metadata.

## Issue 149: baseline vs fixed

| Strategy | Cost Before | QALYs Before | Cost After | QALYs After | Cost Change | QALY Change |
|:---|---:|---:|---:|---:|---:|---:|
| Standard of Care | EUR 21,056 | 1.3498 | EUR 21,523 | 1.3714 | EUR 468 | 0.0216 |
| CRP-guided | EUR 53,791 | 1.3886 | EUR 53,928 | 1.3918 | EUR 137 | 0.0033 |
| TMB/BRAF-guided | EUR 63,235 | 1.3514 | EUR 63,399 | 1.3609 | EUR 164 | 0.0095 |

Base case impact – Issue 149: baseline vs fixed

| Strategy         | NMB Before |  NMB After | NMB Change |
|:-----------------|-----------:|-----------:|-----------:|
| Standard of Care | EUR 47,783 | EUR 48,419 |    EUR 636 |
| CRP-guided       | EUR 17,026 | EUR 17,056 |     EUR 30 |
| TMB/BRAF-guided  |  EUR 5,684 |  EUR 6,005 |    EUR 320 |

Net monetary benefit impact – Issue 149: baseline vs fixed

PSA provenance: before snapshot predates issue \#156 and did not record
the PSA cache md5; after snapshot predates issue \#156 and did not
record the PSA cache md5. Identity of the two PSA files cannot be
established from metadata.

## Issue 151: baseline vs fixed

| Strategy | Cost Before | QALYs Before | Cost After | QALYs After | Cost Change | QALY Change |
|:---|---:|---:|---:|---:|---:|---:|
| Standard of Care | EUR 21,523 | 1.3714 | EUR 21,523 | 1.3714 | EUR 0 | 0.0000 |
| CRP-guided | EUR 53,928 | 1.3918 | EUR 53,928 | 1.3918 | EUR 0 | 0.0000 |
| TMB/BRAF-guided | EUR 63,399 | 1.3609 | EUR 63,399 | 1.3609 | EUR 0 | 0.0000 |

Base case impact – Issue 151: baseline vs fixed

| Strategy         | NMB Before |  NMB After | NMB Change |
|:-----------------|-----------:|-----------:|-----------:|
| Standard of Care | EUR 48,419 | EUR 48,419 |      EUR 0 |
| CRP-guided       | EUR 17,056 | EUR 17,056 |      EUR 0 |
| TMB/BRAF-guided  |  EUR 6,005 |  EUR 6,005 |      EUR 0 |

Net monetary benefit impact – Issue 151: baseline vs fixed

PSA provenance: before snapshot predates issue \#156 and did not record
the PSA cache md5; after snapshot predates issue \#156 and did not
record the PSA cache md5. Identity of the two PSA files cannot be
established from metadata.

## Issue 152: baseline vs fixed

| Strategy | Cost Before | QALYs Before | Cost After | QALYs After | Cost Change | QALY Change |
|:---|---:|---:|---:|---:|---:|---:|
| Standard of Care | EUR 21,523 | 1.3714 | EUR 21,523 | 1.3714 | EUR 0 | 0.0000 |
| CRP-guided | EUR 53,928 | 1.3918 | EUR 53,928 | 1.3918 | EUR 0 | 0.0000 |
| TMB/BRAF-guided | EUR 63,399 | 1.3609 | EUR 63,399 | 1.3609 | EUR 0 | 0.0000 |

Base case impact – Issue 152: baseline vs fixed

| Strategy         | NMB Before |  NMB After | NMB Change |
|:-----------------|-----------:|-----------:|-----------:|
| Standard of Care | EUR 48,419 | EUR 48,419 |      EUR 0 |
| CRP-guided       | EUR 17,056 | EUR 17,056 |      EUR 0 |
| TMB/BRAF-guided  |  EUR 6,005 |  EUR 6,005 |      EUR 0 |

Net monetary benefit impact – Issue 152: baseline vs fixed

PSA provenance: before snapshot predates issue \#156 and did not record
the PSA cache md5; after snapshot predates issue \#156 and did not
record the PSA cache md5. Identity of the two PSA files cannot be
established from metadata.

Per-patient EVPI: EUR 45.20 before, EUR 45.20 after.

| Parameter / group | EVPPI Before | % EVPI Before | EVPPI After | SE After | % EVPI After |
|:---|---:|---:|---:|---:|---:|
| All costs (group) | EUR 45.20 | 100.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| Drug costs (group) | EUR 45.20 | 100.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| TMB/BRAF-treatment interaction (group) | EUR 45.20 | 100.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| Resource-use costs (group) | EUR 45.20 | 100.0% | EUR 0.00 | EUR 0.04 | 0.0% |
| b_tmb_braf_rx_os | EUR 45.20 | 100.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| b_tmb_braf_rx_pfs | EUR 45.20 | 100.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| c_drug_FLOX | EUR 45.20 | 100.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| c_drug_nivo | EUR 45.20 | 100.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| u_np | EUR 45.20 | 100.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| c_test_CT | EUR 43.76 | 96.8% | EUR 0.00 | EUR 0.00 | 0.0% |
| c_other_last | EUR 41.57 | 92.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| p_tmb_braf | EUR 40.09 | 88.7% | EUR 0.00 | EUR 0.00 | 0.0% |
| CRP-treatment interaction (group) | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.74 | 0.0% |
| Biomarker prevalence (group) | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| Test costs (group) | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| Utilities (group) | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| b_crp_rx_os | EUR 0.00 | 0.0% | EUR 0.10 | EUR 1.60 | 0.2% |
| b_crp_rx_pfs | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.12 | 0.0% |
| c_other_baseline | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| c_other_follow | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| c_other_visit | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| c_test_CRP | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| c_test_NGS | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| c_test_blood | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| p_crp | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| u_p | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| Biomarker-treatment interaction (group) | – | – | EUR 0.00 | EUR 1.61 | 0.0% |

EVPPI impact – Issue 152: baseline vs fixed

## Issue 153: baseline vs fixed

| Strategy | Cost Before | QALYs Before | Cost After | QALYs After | Cost Change | QALY Change |
|:---|---:|---:|---:|---:|---:|---:|
| Standard of Care | EUR 21,523 | 1.3714 | EUR 21,523 | 1.3714 | EUR 0 | 0.0000 |
| CRP-guided | EUR 53,928 | 1.3918 | EUR 53,928 | 1.3918 | EUR 0 | 0.0000 |
| TMB/BRAF-guided | EUR 63,399 | 1.3609 | EUR 63,399 | 1.3609 | EUR 0 | 0.0000 |

Base case impact – Issue 153: baseline vs fixed

| Strategy         | NMB Before |  NMB After | NMB Change |
|:-----------------|-----------:|-----------:|-----------:|
| Standard of Care | EUR 48,419 | EUR 48,419 |      EUR 0 |
| CRP-guided       | EUR 17,056 | EUR 17,056 |      EUR 0 |
| TMB/BRAF-guided  |  EUR 6,005 |  EUR 6,005 |      EUR 0 |

Net monetary benefit impact – Issue 153: baseline vs fixed

PSA provenance: before snapshot predates issue \#156 and did not record
the PSA cache md5; after snapshot predates issue \#156 and did not
record the PSA cache md5. Identity of the two PSA files cannot be
established from metadata.

Per-patient EVPI: EUR 45.20 before, EUR 45.20 after.

| Parameter / group | EVPPI Before | % EVPI Before | EVPPI After | SE After | % EVPI After |
|:---|---:|---:|---:|---:|---:|
| b_crp_rx_os | EUR 0.10 | 0.2% | EUR 0.10 | EUR 1.60 | 0.2% |
| All costs (group) | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| Drug costs (group) | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| Biomarker-treatment interaction (group) | EUR 0.00 | 0.0% | EUR 0.00 | EUR 1.61 | 0.0% |
| CRP-treatment interaction (group) | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.74 | 0.0% |
| TMB/BRAF-treatment interaction (group) | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| Resource-use costs (group) | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.04 | 0.0% |
| Biomarker prevalence (group) | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| Test costs (group) | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| Utilities (group) | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| b_crp_rx_pfs | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.12 | 0.0% |
| b_tmb_braf_rx_os | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| b_tmb_braf_rx_pfs | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| c_drug_FLOX | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| c_drug_nivo | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| c_other_baseline | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| c_other_follow | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| c_other_last | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| c_other_visit | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| c_test_CRP | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| c_test_CT | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| c_test_NGS | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| c_test_blood | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| p_crp | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| p_tmb_braf | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| u_np | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| u_p | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |

EVPPI impact – Issue 153: baseline vs fixed

## Issue 154: baseline vs fixed

| Strategy | Cost Before | QALYs Before | Cost After | QALYs After | Cost Change | QALY Change |
|:---|---:|---:|---:|---:|---:|---:|
| Standard of Care | EUR 21,523 | 1.3714 | EUR 21,523 | 1.3714 | EUR 0 | 0.0000 |
| CRP-guided | EUR 53,928 | 1.3918 | EUR 53,928 | 1.3918 | EUR 0 | 0.0000 |
| TMB/BRAF-guided | EUR 63,399 | 1.3609 | EUR 63,399 | 1.3609 | EUR 0 | 0.0000 |

Base case impact – Issue 154: baseline vs fixed

| Strategy         | NMB Before |  NMB After | NMB Change |
|:-----------------|-----------:|-----------:|-----------:|
| Standard of Care | EUR 48,419 | EUR 48,419 |      EUR 0 |
| CRP-guided       | EUR 17,056 | EUR 17,056 |      EUR 0 |
| TMB/BRAF-guided  |  EUR 6,005 |  EUR 6,005 |      EUR 0 |

Net monetary benefit impact – Issue 154: baseline vs fixed

PSA provenance: before snapshot predates issue \#156 and did not record
the PSA cache md5; after snapshot predates issue \#156 and did not
record the PSA cache md5. Identity of the two PSA files cannot be
established from metadata.

Per-patient EVPI: EUR 45.20 before, EUR 22.13 after.

| Parameter / group | EVPPI Before | % EVPI Before | EVPPI After | SE After | % EVPI After |
|:---|---:|---:|---:|---:|---:|
| b_crp_rx_os | EUR 0.10 | 0.2% | EUR 0.09 | EUR 0.95 | 0.4% |
| All costs (group) | EUR 0.00 | 0.0% | – | – | – |
| Drug costs (group) | EUR 0.00 | 0.0% | – | – | – |
| Biomarker-treatment interaction (group) | EUR 0.00 | 0.0% | EUR 1.18 | EUR 1.94 | 5.3% |
| CRP-treatment interaction (group) | EUR 0.00 | 0.0% | EUR 0.04 | EUR 0.41 | 0.2% |
| TMB/BRAF-treatment interaction (group) | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| Resource-use costs (group) | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| Biomarker prevalence (group) | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| Test costs (group) | EUR 0.00 | 0.0% | – | – | – |
| Utilities (group) | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| b_crp_rx_pfs | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.08 | 0.0% |
| b_tmb_braf_rx_os | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| b_tmb_braf_rx_pfs | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| c_drug_FLOX | EUR 0.00 | 0.0% | – | – | – |
| c_drug_nivo | EUR 0.00 | 0.0% | – | – | – |
| c_other_baseline | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| c_other_follow | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| c_other_last | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| c_other_visit | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| c_test_CRP | EUR 0.00 | 0.0% | – | – | – |
| c_test_CT | EUR 0.00 | 0.0% | – | – | – |
| c_test_NGS | EUR 0.00 | 0.0% | – | – | – |
| c_test_blood | EUR 0.00 | 0.0% | – | – | – |
| p_crp | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| p_tmb_braf | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| u_np | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| u_p | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| u_decrement | – | – | EUR 0.00 | EUR 0.00 | 0.0% |

EVPPI impact – Issue 154: baseline vs fixed

## Issue 155: baseline vs fixed

| Strategy | Cost Before | QALYs Before | Cost After | QALYs After | Cost Change | QALY Change |
|:---|---:|---:|---:|---:|---:|---:|
| Standard of Care | EUR 21,523 | 1.3714 | EUR 21,523 | 1.3714 | EUR 0 | 0.0000 |
| CRP-guided | EUR 53,928 | 1.3918 | EUR 53,928 | 1.3918 | EUR 0 | 0.0000 |
| TMB/BRAF-guided | EUR 63,399 | 1.3609 | EUR 63,399 | 1.3609 | EUR 0 | 0.0000 |

Base case impact – Issue 155: baseline vs fixed

| Strategy         | NMB Before |  NMB After | NMB Change |
|:-----------------|-----------:|-----------:|-----------:|
| Standard of Care | EUR 48,419 | EUR 48,419 |      EUR 0 |
| CRP-guided       | EUR 17,056 | EUR 17,056 |      EUR 0 |
| TMB/BRAF-guided  |  EUR 6,005 |  EUR 6,005 |      EUR 0 |

Net monetary benefit impact – Issue 155: baseline vs fixed

PSA provenance: before snapshot predates issue \#156 and did not record
the PSA cache md5; after snapshot predates issue \#156 and did not
record the PSA cache md5. Identity of the two PSA files cannot be
established from metadata.

Per-patient EVPI: EUR 22.13 before, EUR 22.13 after.

| Parameter / group | EVPPI Before | % EVPI Before | EVPPI After | SE After | % EVPI After |
|:---|---:|---:|---:|---:|---:|
| Biomarker-treatment interaction (group) | EUR 1.18 | 5.3% | EUR 1.18 | EUR 1.94 | 5.3% |
| b_crp_rx_os | EUR 0.09 | 0.4% | EUR 0.09 | EUR 0.95 | 0.4% |
| CRP-treatment interaction (group) | EUR 0.04 | 0.2% | EUR 0.04 | EUR 0.41 | 0.2% |
| TMB/BRAF-treatment interaction (group) | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| Resource-use costs (group) | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| Biomarker prevalence (group) | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| Utilities (group) | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| b_crp_rx_pfs | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.08 | 0.0% |
| b_tmb_braf_rx_os | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| b_tmb_braf_rx_pfs | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| c_other_baseline | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| c_other_follow | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| c_other_last | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| c_other_visit | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| p_crp | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| p_tmb_braf | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| u_decrement | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| u_np | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |
| u_p | EUR 0.00 | 0.0% | EUR 0.00 | EUR 0.00 | 0.0% |

EVPPI impact – Issue 155: baseline vs fixed

## Issue 156: baseline vs fixed

| Strategy | Cost Before | QALYs Before | Cost After | QALYs After | Cost Change | QALY Change |
|:---|---:|---:|---:|---:|---:|---:|
| Standard of Care | EUR 21,523 | 1.3714 | EUR 21,523 | 1.3714 | EUR 0 | 0.0000 |
| CRP-guided | EUR 53,928 | 1.3918 | EUR 53,928 | 1.3918 | EUR 0 | 0.0000 |
| TMB/BRAF-guided | EUR 63,399 | 1.3609 | EUR 63,399 | 1.3609 | EUR 0 | 0.0000 |

Base case impact – Issue 156: baseline vs fixed

| Strategy         | NMB Before |  NMB After | NMB Change |
|:-----------------|-----------:|-----------:|-----------:|
| Standard of Care | EUR 48,419 | EUR 48,419 |      EUR 0 |
| CRP-guided       | EUR 17,056 | EUR 17,056 |      EUR 0 |
| TMB/BRAF-guided  |  EUR 6,005 |  EUR 6,005 |      EUR 0 |

Net monetary benefit impact – Issue 156: baseline vs fixed

PSA provenance: before snapshot predates issue \#156 and did not record
the PSA cache md5; after PSA cache md5 fa17a8dab004ee9fc4c850ff9897b656
(modified 2026-09-06 23:39). Identity of the two PSA files cannot be
established from metadata.

Per-patient EVPI: EUR 22.13 before, EUR 0.00 after.

## Issue 157: baseline vs fixed

| Strategy | Cost Before | QALYs Before | Cost After | QALYs After | Cost Change | QALY Change |
|:---|---:|---:|---:|---:|---:|---:|
| Standard of Care | EUR 21,523 | 1.3714 | EUR 21,523 | 1.3714 | EUR 0 | 0.0000 |
| CRP-guided | EUR 53,928 | 1.3918 | EUR 53,928 | 1.3918 | EUR 0 | 0.0000 |
| TMB/BRAF-guided | EUR 63,399 | 1.3609 | EUR 63,399 | 1.3609 | EUR 0 | 0.0000 |

Base case impact – Issue 157: baseline vs fixed

| Strategy         | NMB Before |  NMB After | NMB Change |
|:-----------------|-----------:|-----------:|-----------:|
| Standard of Care | EUR 48,419 | EUR 48,419 |      EUR 0 |
| CRP-guided       | EUR 17,056 | EUR 17,056 |      EUR 0 |
| TMB/BRAF-guided  |  EUR 6,005 |  EUR 6,005 |      EUR 0 |

Net monetary benefit impact – Issue 157: baseline vs fixed

PSA cache md5 before fa17a8dab004ee9fc4c850ff9897b656 (modified
2026-09-06 23:39), after a057d7d7856e7c3bc7602058118d9a92 (modified
2026-09-07 09:28). The after PSA was regenerated inside the snapshot run
because the cache predated the fix commit.

Per-patient EVPI: EUR 0.00 before, EUR 0.00 after.

## Issue 159: baseline vs fixed

| Strategy | Cost Before | QALYs Before | Cost After | QALYs After | Cost Change | QALY Change |
|:---|---:|---:|---:|---:|---:|---:|
| Standard of Care | EUR 21,523 | 1.3714 | EUR 21,523 | 1.3714 | EUR 0 | 0.0000 |
| CRP-guided | EUR 53,948 | 1.3920 | EUR 53,948 | 1.3920 | EUR 0 | 0.0000 |
| TMB/BRAF-guided | EUR 62,685 | 1.3599 | EUR 62,685 | 1.3599 | EUR 0 | 0.0000 |

Base case impact – Issue 159: baseline vs fixed

| Strategy         | NMB Before |  NMB After | NMB Change |
|:-----------------|-----------:|-----------:|-----------:|
| Standard of Care | EUR 48,419 | EUR 48,419 |      EUR 0 |
| CRP-guided       | EUR 17,042 | EUR 17,042 |      EUR 0 |
| TMB/BRAF-guided  |  EUR 6,670 |  EUR 6,670 |      EUR 0 |

Net monetary benefit impact – Issue 159: baseline vs fixed

PSA cache md5 before 39f39bd3f4f56a8fcf25acd341b6e560 (modified
2026-09-17 16:48), after 90645792c191da22ce9a8c81ed323ff9 (modified
2026-09-18 17:20).

Per-patient EVPI: EUR 0.22 before, EUR 0.00 after.

## Issue 165: baseline vs fixed

| Strategy | Cost Before | QALYs Before | Cost After | QALYs After | Cost Change | QALY Change |
|:---|---:|---:|---:|---:|---:|---:|
| Standard of Care | EUR 25,764 | 1.3400 | EUR 25,764 | 1.3400 | EUR 0 | 0.0000 |
| CRP-guided | EUR 58,203 | 1.3775 | EUR 58,203 | 1.3775 | EUR 0 | 0.0000 |
| TMB/BRAF-guided | EUR 66,916 | 1.3389 | EUR 66,916 | 1.3389 | EUR 0 | 0.0000 |

Base case impact – Issue 165: baseline vs fixed

| Strategy         | NMB Before |  NMB After | NMB Change |
|:-----------------|-----------:|-----------:|-----------:|
| Standard of Care | EUR 42,574 | EUR 42,574 |      EUR 0 |
| CRP-guided       | EUR 12,048 | EUR 12,048 |      EUR 0 |
| TMB/BRAF-guided  |  EUR 1,368 |  EUR 1,368 |      EUR 0 |

Net monetary benefit impact – Issue 165: baseline vs fixed

PSA cache md5 before 7a04965187ff820c11f6d3ad0ccab9e8 (modified
2026-09-30 16:51), after 240069919ae165e572113ebffdd9dc9d (modified
2026-10-01 13:31).

Per-patient EVPI: EUR 0.00 before, EUR 0.00 after.

## Issue 166: baseline vs fixed

| Strategy | Cost Before | QALYs Before | Cost After | QALYs After | Cost Change | QALY Change |
|:---|---:|---:|---:|---:|---:|---:|
| Standard of Care | EUR 21,523 | 1.3714 | EUR 21,523 | 1.3714 | EUR 0 | 0.0000 |
| CRP-guided | EUR 53,928 | 1.3918 | EUR 53,948 | 1.3920 | EUR 20 | 0.0001 |
| TMB/BRAF-guided | EUR 63,399 | 1.3609 | EUR 62,685 | 1.3599 | -EUR 715 | -0.0010 |

Base case impact – Issue 166: baseline vs fixed

| Strategy         | NMB Before |  NMB After | NMB Change |
|:-----------------|-----------:|-----------:|-----------:|
| Standard of Care | EUR 48,419 | EUR 48,419 |      EUR 0 |
| CRP-guided       | EUR 17,056 | EUR 17,042 |    -EUR 14 |
| TMB/BRAF-guided  |  EUR 6,005 |  EUR 6,670 |    EUR 665 |

Net monetary benefit impact – Issue 166: baseline vs fixed

PSA cache md5 before 8f061b8959fee3e221a34f65ce0d9dc1 (modified
2026-09-18 09:55), after 39f39bd3f4f56a8fcf25acd341b6e560 (modified
2026-09-17 16:48).

Per-patient EVPI: EUR 0.00 before, EUR 0.22 after.

## Issue 168: baseline vs fixed

| Strategy | Cost Before | QALYs Before | Cost After | QALYs After | Cost Change | QALY Change |
|:---|---:|---:|---:|---:|---:|---:|
| Standard of Care | EUR 21,523 | 1.3644 | EUR 21,523 | 1.3644 | EUR 0 | 0.0000 |
| CRP-guided | EUR 53,948 | 1.3849 | EUR 53,948 | 1.3849 | EUR 0 | 0.0000 |
| TMB/BRAF-guided | EUR 62,685 | 1.3529 | EUR 62,685 | 1.3529 | EUR 0 | 0.0000 |

Base case impact – Issue 168: baseline vs fixed

| Strategy         | NMB Before |  NMB After | NMB Change |
|:-----------------|-----------:|-----------:|-----------:|
| Standard of Care | EUR 48,061 | EUR 48,061 |      EUR 0 |
| CRP-guided       | EUR 16,684 | EUR 16,684 |      EUR 0 |
| TMB/BRAF-guided  |  EUR 6,312 |  EUR 6,312 |      EUR 0 |

Net monetary benefit impact – Issue 168: baseline vs fixed

PSA cache md5 before 943be17d86d24aa1a30a17ec1d4238fa (modified
2026-09-21 15:41), after ca7e7d9091880b1f5879fe9cd8a7f967 (modified
2026-09-22 12:04). The after PSA was regenerated inside the snapshot run
because the cache predated the fix commit.

Per-patient EVPI: EUR 0.00 before, EUR 0.00 after.

## Issue 172: baseline vs fixed

| Strategy | Cost Before | QALYs Before | Cost After | QALYs After | Cost Change | QALY Change |
|:---|---:|---:|---:|---:|---:|---:|
| Standard of Care | EUR 20,967 | 1.3400 | EUR 20,967 | 1.3400 | EUR 0 | 0.0000 |
| CRP-guided | EUR 53,686 | 1.3775 | EUR 53,686 | 1.3775 | EUR 0 | 0.0000 |
| TMB/BRAF-guided | EUR 62,382 | 1.3389 | EUR 62,382 | 1.3389 | EUR 0 | 0.0000 |

Base case impact – Issue 172: baseline vs fixed

| Strategy         | NMB Before |  NMB After | NMB Change |
|:-----------------|-----------:|-----------:|-----------:|
| Standard of Care | EUR 47,371 | EUR 47,371 |      EUR 0 |
| CRP-guided       | EUR 16,565 | EUR 16,565 |      EUR 0 |
| TMB/BRAF-guided  |  EUR 5,901 |  EUR 5,901 |      EUR 0 |

Net monetary benefit impact – Issue 172: baseline vs fixed

PSA cache md5 before 05d5d3170708cebb2432ec90648ff1f8 (modified
2026-09-29 12:50), after cfeb0a0395569b646c13e0ce3d026c12 (modified
2026-09-29 13:57).

Per-patient EVPI: EUR 0.00 before, EUR 0.00 after.

## Batch \#173 / \#171 / \#163: baseline vs fixed

| Strategy | Cost Before | QALYs Before | Cost After | QALYs After | Cost Change | QALY Change |
|:---|---:|---:|---:|---:|---:|---:|
| Standard of Care | EUR 21,523 | 1.3714 | EUR 21,523 | 1.3644 | EUR 0 | -0.0070 |
| CRP-guided | EUR 53,948 | 1.3920 | EUR 53,948 | 1.3849 | EUR 0 | -0.0070 |
| TMB/BRAF-guided | EUR 62,685 | 1.3599 | EUR 62,685 | 1.3529 | EUR 0 | -0.0070 |

Base case impact – Batch \#173 / \#171 / \#163: baseline vs fixed

| Strategy         | NMB Before |  NMB After | NMB Change |
|:-----------------|-----------:|-----------:|-----------:|
| Standard of Care | EUR 48,419 | EUR 48,061 |   -EUR 358 |
| CRP-guided       | EUR 17,042 | EUR 16,684 |   -EUR 358 |
| TMB/BRAF-guided  |  EUR 6,670 |  EUR 6,312 |   -EUR 358 |

Net monetary benefit impact – Batch \#173 / \#171 / \#163: baseline vs
fixed

PSA cache md5 before 90645792c191da22ce9a8c81ed323ff9 (modified
2026-09-18 17:20), after 943be17d86d24aa1a30a17ec1d4238fa (modified
2026-09-21 15:41).

Per-patient EVPI: EUR 0.00 before, EUR 0.00 after.

## Issue 181: baseline vs fixed

| Strategy | Cost Before | QALYs Before | Cost After | QALYs After | Cost Change | QALY Change |
|:---|---:|---:|---:|---:|---:|---:|
| Standard of Care | EUR 21,523 | 1.3644 | EUR 20,967 | 1.3400 | -EUR 557 | -0.0244 |
| CRP-guided | EUR 53,948 | 1.3849 | EUR 53,686 | 1.3775 | -EUR 262 | -0.0075 |
| TMB/BRAF-guided | EUR 62,685 | 1.3529 | EUR 62,382 | 1.3389 | -EUR 302 | -0.0140 |

Base case impact – Issue 181: baseline vs fixed

| Strategy         | NMB Before |  NMB After | NMB Change |
|:-----------------|-----------:|-----------:|-----------:|
| Standard of Care | EUR 48,061 | EUR 47,371 |   -EUR 690 |
| CRP-guided       | EUR 16,684 | EUR 16,565 |   -EUR 119 |
| TMB/BRAF-guided  |  EUR 6,312 |  EUR 5,901 |   -EUR 410 |

Net monetary benefit impact – Issue 181: baseline vs fixed

PSA cache md5 before ca7e7d9091880b1f5879fe9cd8a7f967 (modified
2026-09-22 12:04), after ace97153eb4660b2df5642cbc4348c60 (modified
2026-09-24 17:54).

Per-patient EVPI: EUR 0.00 before, EUR 0.00 after.

## Cumulative: issue 145 baseline vs issue 165 fixed

| Strategy | Cost Before | QALYs Before | Cost After | QALYs After | Cost Change | QALY Change |
|:---|---:|---:|---:|---:|---:|---:|
| Standard of Care | EUR 20,867 | 1.4141 | EUR 25,764 | 1.3400 | EUR 4,897 | -0.0742 |
| CRP-guided | EUR 55,026 | 1.3842 | EUR 58,203 | 1.3775 | EUR 3,177 | -0.0068 |
| TMB/BRAF-guided | EUR 60,798 | 1.3530 | EUR 66,916 | 1.3389 | EUR 6,118 | -0.0141 |

Base case impact – Cumulative: issue 145 baseline vs issue 165 fixed

| Strategy         | NMB Before |  NMB After | NMB Change |
|:-----------------|-----------:|-----------:|-----------:|
| Standard of Care | EUR 51,254 | EUR 42,574 | -EUR 8,680 |
| CRP-guided       | EUR 15,570 | EUR 12,048 | -EUR 3,522 |
| TMB/BRAF-guided  |  EUR 8,207 |  EUR 1,368 | -EUR 6,839 |

Net monetary benefit impact – Cumulative: issue 145 baseline vs issue
165 fixed

PSA provenance: before snapshot predates issue \#156 and did not record
the PSA cache md5; after PSA cache md5 240069919ae165e572113ebffdd9dc9d
(modified 2026-10-01 13:31). Identity of the two PSA files cannot be
established from metadata.

# Interpretation

This report compares baseline and fixed single-model snapshots for each
bug-fix issue that has both, plus a cumulative comparison from the
earliest baseline to the latest fixed snapshot. Snapshots are produced
by the `13_save_snapshot.R` analysis script (in `analysis/`). If older
snapshots contain obsolete multi-model components, `load_snapshot()`
warns during render and the tables use the top-level single-model
results.

Issue \#166 aligns every economic strategy to the same 68-patient
complete-case target. Prevalence uncertainty now samples joint biomarker
cell masses and reweights control and guided strategies together. Both
snapshots precede the requested fix commit and therefore carry HEAD
`96d0710`; their distinct PSA fingerprints and result hashes identify
the original and repaired calculations. Sampling, PSA, EVPPI and
scenario caches were regenerated for the fix. The identical-treatment
population tests and all 296 report/cache contracts pass; the existing
PSA-versus-base-case numerical alignment criterion remains a failure,
documented in `AGENTS.md`.

Issue \#159 implements the first proposed fix: paired patient bootstrap
fits estimate OS/PFS coefficient dependence for joint normal draws,
retaining both fitted marginal covariance matrices. The deterministic
results are unchanged. Crossings fell from 2,422/5,000 draws (48.44%) to
1,198/5,000 (23.96%), and EVPI at EUR 51,000 fell from EUR 0.22 to zero.
The CRP incremental PSA mean changed from 0.00175 to 0.00628 QALYs,
versus a deterministic increment of 0.02054. Removing all uncertainty
reproduces every deterministic result, but the original five-SE
numerical alignment criterion still fails for four of six level
comparisons. No draws were recentered and no threshold was relaxed;
nonlinear outcome averaging remains distinct from evaluation at fitted
coefficients. See [the joint sampling
specification](joint_survival_sampling.md).

The \#159 baseline and fixed snapshots both precede the requested commit
and carry HEAD `6ca19ab`, with distinct PSA fingerprints and result
hashes. All active generated caches were cleared and rebuilt, and
reports and publication outputs were rendered before the fixed snapshot.
The cumulative comparison uses saved snapshot timestamps, because \#159
was implemented after \#166.

Issue \#168 repairs publication artifacts and their provenance. Table S7
and Figure 5 now retain the base case with EVPI zero, distinguish failed
estimates from zeros, and use scenario PSA objects to establish scenario
totals. All 26 publication vignettes were rendered, all 17 tracked CSVs
were checked against their generating objects, and the 37 vignette-owned
artifacts are recorded in `outputs/manifest.csv`. Table 5 already agreed
with the current PSA at published precision; its new regression checks
guard against future drift.

The \#168 baseline and fixed snapshots were saved with the usual
`baseline` and `fixed` arguments before the fix commit, so both carry
HEAD `b8e1b02`. The fixed snapshot regenerated all 5,000 PSA draws and
EVPPI because the previous cache predated HEAD. Deterministic costs and
QALYs are identical; maximum absolute PSA differences are EUR 2.91e-11
and 8.88e-16 QALYs, and EVPI remains zero. Scenario compilation was
rebuilt from validated existing scenario results after verifying that
the estimation functions and inputs were unchanged. The existing five-SE
diagnostic still fails for four of six level comparisons; its criterion
is unchanged. See [the validation
record](../../validation/issue168_2026-09-22/README.md).

------------------------------------------------------------------------

**Report completed on:** 2026-10-01  
**Repository:** ben-geisler/METIMMOX-1  
**Report version:** 4.7

## Validation warnings

Warnings recorded during rendering: 1 (1 distinct).

- report setup: package ‘survival’ was built under R version 4.3.3
